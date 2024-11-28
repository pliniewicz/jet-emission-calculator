#include <math.h>
#include "profiles.h"
#include "rad_functions.h"
#include <gsl/gsl_integration.h>
#include <stdio.h>
#include "config.h"
#include "electron_distributions.h"

void logspaced(double first, double last, int N, double space[]) {

  space[0] = first;

  double step = pow(last / first, 1. / (N - 1.));

  for (int i = 0; i < N; i++) {
    space[i+1] = space[i]*step;
  };

}

typedef double (*F_ptr)(double, void *);

double find_Bfield_normalization(F_ptr F, void *params){

  gsl_integration_cquad_workspace *BfieldNormalization = gsl_integration_cquad_workspace_alloc(5000);
  double pB1 = 0;
  double q = *(double *)params;
  // gsl_function pBintegration;
  // pBintegration.function = profile_function;
  // pBintegration.params = &q;
  gsl_function Func;
  Func.function = F;
  Func.params = params;

  gsl_integration_cquad(&Func, 0, 1, 0, 1e-15, BfieldNormalization, &pB1, NULL, NULL);
  gsl_integration_cquad_workspace_free(BfieldNormalization);

  pB1 = q * LUMJET / (8 * M_PI * SPEEDC * gsl_pow_2(RJET)*pB1);

  return pB1;

}

double normalize_distribution(int profile, void *params){
  gsl_integration_cquad_workspace *distributionNormalization = gsl_integration_cquad_workspace_alloc(200);

  // W ogóle nie następuje przekazanie parametrów ze struktury do power law

  double Ke = 0;

  switch (profile) {
    default:
      printf("Normalizing for the simple power law\n");
      struct normalizationParametersSimplePowerLaw *p = (struct normalizationParametersSimplePowerLaw *)params;
      double gmin = p->gmin;
      double gmax = p->gmax;

      gsl_function F;
      F.function = &simple_power_law;
      F.params = &p;
      
      gsl_integration_cquad(&F, gmin, gmax, 1e-15, 1e-15, distributionNormalization, &Ke, NULL, NULL);

      gsl_integration_cquad_workspace_free(distributionNormalization);
      return Ke;
      break;
  }

}
