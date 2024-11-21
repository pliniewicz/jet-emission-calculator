#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_exp.h>
#include <math.h>
#include "logspace.h"

struct normalizationParametersSimplePowerLaw {double gmin; double gmax; int index;};

double simple_power_law(double x, void *params)
{
  struct normalizationParametersSimplePowerLaw *p = (struct normalizationParametersSimplePowerLaw *)params;
  double gmin = p->gmin;
  double gmax = p->gmax;
  int index = p->index;
  if (x <= gmax || x >= gmin) {
    return 1. / gsl_sf_pow_int(x, index);
  } else {
    return 0;
  }
}
