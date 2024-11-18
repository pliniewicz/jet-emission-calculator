#include <gsl/gsl_math.h>
#include <gsl/gsl_pow_int.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>

#include <gsl/gsl_mode.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_mksa.h>

#include "rad_functions.h"
#include "profiles.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>


double synchrotron_integrand(double gamma, double frequency, double magfield, double gmin, double gmax, double angle)
{
  // struct synchrotron_kernel_parameters *p = (struct synchrotron_kernel_parameters *)parameters;
  //
  // double frequency = p->frequency;
  // double magnetic_field = p->magnetic_field;
  // double minimal_gamma = p->minimal_gamma;
  // double maximal_gamma = p->maximal_gamma;
  // double angle = p->angle;

  double function_argument = (frequency / nu_critical(gamma, magfield));
  if (function_argument <= 100){
    double k43 = gsl_sf_bessel_Knu(4./3., 0.5 * function_argument);
    double k13 = gsl_sf_bessel_Knu(1./3., 0.5 * function_argument);

    double x22 = 0.5*gsl_pow_2(function_argument);
    double x33 = 0.15*gsl_pow_3(function_argument);

    double Rx = x22*k43*k13 - x33*(k43*k43 - k13*k13);
    double ne = exp_cutoff_power_law_2(gamma, gmax, gmin);
    return Rx*ne;
  } else

    return 0; 
  } 

struct synchro_integrand_params {double frequency; double gmin; double gmax; double angle;};

double g (double *k, size_t dim, void *params){
  (void)(dim);
  struct synchro_integrand_params * fp = (struct synchro_integrand_params *)params;
  double frequency = fp->frequency;
  double gmin = fp->gmin;
  double gmax = fp->gmax;
  double angle = fp->angle;
  double x = k[1];
  double gamma = k[0];
  return x*b1(x)*gsl_pow_2(doppler_profile(x, angle))*f1(x)*exp_cutoff_power_law_2(gamma, gmin, gmax)*synchrotron_integrand(gamma, frequency, b1(x), gmin, gmax, angle);
}


double synchroA(double freq, double gmin, double gmax, double angle){
  double result = 0;
  double error = 0;
  double xl[2] = {gmin, 0};
  double xu[2] = {gmax, 1};

  struct synchro_integrand_params parameters = {freq, gmin, gmax, angle};

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G = {&g, 2, 0};

  G.params = &parameters;
  G.dim = 2;

  size_t calls = 500000;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  {
    gsl_monte_miser_state *s = gsl_monte_miser_alloc (2);
    // gsl_monte_miser_params_set(s, &parameters);
    gsl_monte_miser_integrate(&G, xl, xu, 2, calls, r, s, &result, &error);

    gsl_monte_miser_free(s);
  }
  return result;
}

