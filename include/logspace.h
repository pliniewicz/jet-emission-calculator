#ifndef LOGSPACE
#define LOGSPACE

#include "math.h"
#include "profiles.h"
#include <gsl/gsl_integration.h>
#include "config.h"

void logspaced(double first, double last, int N, double space[]);

typedef double (*F_ptr)(double, void *);

double find_Bfield_normalization(F_ptr F, void *params);

#endif // !LOGSPACE
