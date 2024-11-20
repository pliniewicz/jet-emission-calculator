#ifndef LOGSPACE
#define LOGSPACE

#include "math.h"
#include "profiles.h"
#include <gsl/gsl_integration.h>
#include "config.h"

void logspaced(double first, double last, int N, double space[]);

typedef double (*F_ptr)(double, void *);

double find_Bfield_normalization(F_ptr F, void *params);

struct normalizationParametersSimplePowerLaw {double gmin; double gmax; double index;};

// double normalize_distribution(F_ptr F, void *params, int profile);
double normalize_distribution(int profile, void *params);

#endif // !LOGSPACE
