#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_pow_int.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_mode.h>
#include <gsl/gsl_sys.h>
#include <math.h>
#include "logspace.h"

struct synchrotron_kernel_parameters {double frequency; double magnetic_field; double minimal_gamma; double maximal_gamma; double angle;};
struct synchrotron_kernel_ball_parameters {double frequency; double bulk_gamma; double magnetic_field; double minimal_gamma; double maximal_gamma; double breaking_gamma; double k1; double sL; double sH; double angle;};
struct IC_kernel_parameters {double target_energy; double scattered_energy; double minimal_gamma; double maximal_gamma;};
struct IC_kernel_ball_parameters {double target_energy; double scattered_energy; double minimal_gamma; double maximal_gamma; double breaking_gamma; double k1; double sL; double sH;};

double internal_e=GSL_CONST_CGSM_ELECTRON_CHARGE * GSL_CONST_CGSM_SPEED_OF_LIGHT;
double internal_me=GSL_CONST_CGSM_MASS_ELECTRON;
double internal_c=GSL_CONST_CGSM_SPEED_OF_LIGHT;
double internal_pc=GSL_CONST_CGSM_PARSEC;
double internal_h = GSL_CONST_CGSM_PLANCKS_CONSTANT_H;
double internal_sigmaT = GSL_CONST_CGSM_THOMSON_CROSS_SECTION;
double internal_kB = GSL_CONST_CGSM_BOLTZMANN;


// double[] logspace(double first, double last, int num)
// {
//   double logspace[num];
//   logspace[0] = first;
//   double step = pow(last / first , 1./(num - 1.));
//   for (int i = 0 ; i<=num ; i++) {
//     logspace[i+1] = logspace[i]*step;
//   }
//   return logspace;
// }

double doppler(double gamma, double angle)
{
  double cosine_degree = gsl_sf_cos(angle * M_PI / 180.);
  double beta = sqrt(1. - 1. / (gsl_pow_2(gamma)));
  return 1. / (gamma * (1. - beta*cosine_degree));
}

double nu_critical(double gamma, double magnetic_field)
{
  return ( 3. * internal_e * magnetic_field * gsl_pow_2(gamma) / (4. * M_PI * internal_me * internal_c) );
}

// struct SimplePowerLawParams{
//   double gmin;
//   double gmax;
//   double index;
// };

double simple_power_law(double x, void *params)
{
  struct normalizationParametersSimplePowerLaw *p = (struct normalizationParametersSimplePowerLaw *)params;
  double gmin = p->gmin;
  double gmax = p->gmax;
  double index = p->index;
  if (x <= gmax || x >= gmin) {
    return x*pow(x, -1.*index);
  } else {
    return 0;
  }
  // if (x < gmin ) {
  //   return 0;
  // } else if (x > gmax){
  //   return 0;
  // } else {
  // return x*pow(x, -1.*index);
  // }
}

// double simple_power_law(double gamma, double minimal_gamma, double maximal_gamma, double power_law_index)
// {
//   if (gamma <= minimal_gamma ) {
//     return 0;
//   } else if (gamma >= maximal_gamma){
//     return 0;
//   } else
//   return pow(gamma, -1.*power_law_index);
// }

double simple_power_law_2(double gamma, double minimal_gamma, double maximal_gamma)
{
  if (gamma <= minimal_gamma ) {
    return 0;
  } else if (gamma >= maximal_gamma){
    return 0;
  } else
  return 1. / gsl_pow_2(gamma);
}

double broken_power_law(double gamma, double minimal_gamma, double maximal_gamma, double breaking_gamma, double k1, double sL, double sH)
{
  double k2 = k1 * (gsl_pow_int(breaking_gamma, -1.*sL) / gsl_pow_int(breaking_gamma, -1.*sH));
  if (gamma <= minimal_gamma) {
    return 0;
  } else if (gamma >= maximal_gamma){
    return 0;
  } else if (gamma > minimal_gamma && gamma < breaking_gamma) {
    return k1*pow(gamma, -1.*sL);
  } else {
    return k2*pow(gamma , -1.*sH);
    // double Heaviside = 1. / (1. + gsl_sf_exp(gamma - breaking_gamma));
    // return k1*pow(gamma, -1.*sL)*Heaviside + k2*pow(gamma, -1.*sH )*(1. - Heaviside);
}
}

double exp_cutoff_power_law_2(double gamma, double minimal_gamma, double maximal_gamma)
{
  return gsl_sf_exp(-gamma/maximal_gamma)*gsl_sf_exp(-minimal_gamma/gamma)/gsl_pow_2(gamma);
}
double exp_cutoff_power_law(double gamma, double minimal_gamma, double maximal_gamma, double power_law_index)
{
  return pow(gamma, -1.*power_law_index)*gsl_sf_exp(-gamma/maximal_gamma)*gsl_sf_exp(-minimal_gamma/gamma);
}

double broken_exp_cutoff_power_law(double gamma); // %TODO

double synchrotron_kernel_ball(double gamma, void *parameters)
{
  struct synchrotron_kernel_ball_parameters *p = (struct synchrotron_kernel_ball_parameters *)parameters;

  double frequency = p->frequency;
  double bulk_gamma = p->bulk_gamma;
  double magnetic_field = p->magnetic_field;
  double minimal_gamma = p->minimal_gamma;
  double maximal_gamma = p->maximal_gamma;
  double breaking_gamma = p->breaking_gamma;
  double k1 = p->k1;
  double sL = p->sL;
  double sH = p->sH;
  double angle = p->angle;

  double function_argument = (frequency / nu_critical(gamma, magnetic_field) );
  if (function_argument <= 500){
    double k43 = gsl_sf_bessel_Knu(4./3., 0.5 * function_argument);
    double k13 = gsl_sf_bessel_Knu(1./3., 0.5 * function_argument);

    double x22 = 0.5*gsl_pow_2(function_argument);
    double x33 = 0.15*gsl_pow_3(function_argument);

    double Rx = x22*k43*k13 - x33*(k43*k43 - k13*k13);
    double ne = broken_power_law(gamma, minimal_gamma, maximal_gamma, breaking_gamma, k1, sL, sH);
    return Rx*ne;
  } else

    return 0; 
  } 


double synchrotron_kernel(double gamma, void *parameters)
{
  struct synchrotron_kernel_parameters *p = (struct synchrotron_kernel_parameters *)parameters;

  double frequency = p->frequency;
  double magnetic_field = p->magnetic_field;
  double minimal_gamma = p->minimal_gamma;
  double maximal_gamma = p->maximal_gamma;
  double angle = p->angle;

  double function_argument = (frequency / nu_critical(gamma, magnetic_field));
  if (function_argument <= 1003){
    double k43 = gsl_sf_bessel_Knu(4./3., 0.5 * function_argument);
    double k13 = gsl_sf_bessel_Knu(1./3., 0.5 * function_argument);

    double x22 = 0.5*gsl_pow_2(function_argument);
    double x33 = 0.15*gsl_pow_3(function_argument);

    double Rx = x22*k43*k13 - x33*(k43*k43 - k13*k13);
    double ne = exp_cutoff_power_law_2(gamma, minimal_gamma, maximal_gamma);
    return Rx*ne;
  } else

    return 0; 
  } 

double qu(double scattered_energy, double target_energy, double gamma){return scattered_energy / (4*target_energy*gamma*(gamma - scattered_energy));}
double Qu(double target_energy, double gamma){return 4*target_energy*gamma;}

double f_IC(double scattered_energy, double target_energy, double gamma)
{
  double qq = qu(scattered_energy, target_energy, gamma);
  double Qq = Qu(target_energy, gamma);

  // double f;
  // if (qq <= 1 && qq >= 1 / (4*gamma*gamma) ) {
  //   f = 2 * qq * gsl_sf_log(qq) + qq + 1. - 2*gsl_pow_2(qq) + ( (gsl_pow_2(qq*Qq))*(1. - qq) / (2 * (1 + Qq * qq)) );
  //   return f;
  // } else {
  //   return 0;
  // }
  
  double f;
  if (qq >= 1) {
   return 0;
  } else if (qq <= 1 / (4*gamma*gamma)) {
    return 0;
  } else {
    f = 2 * qq * gsl_sf_log(qq) + qq + 1. - 2*gsl_pow_2(qq) + ( (gsl_pow_2(qq*Qq))*(1. - qq) / (2 * (1 + Qq * qq)) );
    return f;
  }
}

double inverse_compton_kernel(double gamma, void *parameters)
{
  struct IC_kernel_parameters * p = (struct IC_kernel_parameters *)parameters;
  double target_energy = p->target_energy;
  double scattered_energy = p->scattered_energy;
  double minimal_gamma = p->minimal_gamma;
  double maximal_gamma = p->maximal_gamma;

  double f = f_IC(target_energy, scattered_energy, gamma);
  double ne = exp_cutoff_power_law_2(gamma, minimal_gamma, maximal_gamma);

  return ne * f / gsl_pow_2(gamma);
}

double inverse_compton_ball_kernel(double gamma, void *parameters)
{
  struct IC_kernel_ball_parameters * p = (struct IC_kernel_ball_parameters *)parameters;
  double target_energy = p->target_energy;
  double scattered_energy = p->scattered_energy;
  double minimal_gamma = p->minimal_gamma;
  double maximal_gamma = p->maximal_gamma;
  double breaking_gamma = p->breaking_gamma;
  double k1 = p->k1;
  double sL = p->sL;
  double sH = p->sH;

  double f = f_IC(target_energy, scattered_energy, gamma);

  double ne = broken_power_law(gamma, minimal_gamma, maximal_gamma, breaking_gamma, k1, sL, sH);

  return ne * f / (target_energy * gamma * gamma);
}
