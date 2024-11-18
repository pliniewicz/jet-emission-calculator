#include "profiles.h"
#include <gsl/gsl_pow_int.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <math.h>
#include <gsl/gsl_integration.h>

double Gamma(double x)
{
  double G;
  G = 1 + 9*(1 - gsl_pow_2(x));
  return G;
}

double b1(double x)
{
  double b;
  b = 81.*gsl_pow_3(x)/(1. + 80.*gsl_pow_3(x));
  return b;
}

double b2(double x)
{
  double b;
  b = 81.*x/(1. + 80.*gsl_pow_4(x));
  return b;
}
double b3(double x)
{
  if (x<= 0) {
    return 0;
  } else {
  return 2*gsl_sf_sin(4.*M_PI*x)/x + gsl_pow_3(x);
  }
}

double f1(double x){
  return gsl_pow_2(b1(x))/gsl_pow_2(Gamma(x));
}
double f2(double x){
  return gsl_pow_2(b2(x))/gsl_pow_2(Gamma(x));
}
double f3(double x){
  return gsl_pow_2(b3(x))/gsl_pow_2(Gamma(x));
}

double p1_integrand(double s, void * parameters){return f1(s)/s;};
double p2_integrand(double s, void * parameters){return f2(s)/s;};
double p3_integrand(double s, void * parameters){return f3(s)/s;};


double p1(double x, double q){
 gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(100);
 gsl_function F;
 F.function = &p1_integrand;
 double p1;
 gsl_integration_cquad(&F, x, 1, 0, 1e-12, w, &p1, NULL, NULL);
 gsl_integration_cquad_workspace_free(w);
 p1 = 1 + q - q*f1(x) + 2*q*p1;
 return p1;
}

double p2(double x, double q){
 gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(100);
 gsl_function F;
 F.function = &p2_integrand;
 double p;
 gsl_integration_cquad(&F, x, 1, 0, 1e-12, w, &p, NULL, NULL);
 gsl_integration_cquad_workspace_free(w);
 p = 1 + q - q*f1(x) + 2*q*p;
 return p;
}

double p3(double x, double q){
 gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(100);
 gsl_function F;
 F.function = &p3_integrand;
 double p;
 gsl_integration_cquad(&F, x, 1, 0, 1e-12, w, &p, NULL, NULL);
 gsl_integration_cquad_workspace_free(w);
 p = 1 + q - q*f1(x) + 2*q*p;
 return p;
}

double relativistic_beta(double x){
 return sqrt( 1. - 1./gsl_pow_2( Gamma(x) ) );
}

double doppler_profile(double x, double theta){
  double cosine_degree = gsl_sf_cos(theta * M_PI / 180.);
  double d = 1. / ( Gamma(x) * (1. - relativistic_beta(x)*cosine_degree) );
  return d;
}

double pBint1(double x, void * params){
  double q = *(double *) params;
  return x*relativistic_beta(x)*gsl_pow_2( Gamma(x) )*p1(x, q);
}
double pBint2(double x, void * params){
  double q = *(double *) params;
  return x*relativistic_beta(x)*gsl_pow_2( Gamma(x) )*p2(x, q);
}
double pBint3(double x, void * params){
  double q = *(double *) params;
  return x*relativistic_beta(x)*gsl_pow_2( Gamma(x) )*p3(x, q);
}

