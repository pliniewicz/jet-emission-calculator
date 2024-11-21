#ifndef ELECTRON_DISTRIBUTION
#define ELECTRON_DISTRIBUTION

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_exp.h>
#include <math.h>
#include "logspace.h"

struct normalizationParametersSimplePowerLaw {
  double gmin;
  double gmax;
  int index;
};

double simple_power_law(double x, void *params);


#endif // !ELECTRON_DISTRIBUTION

