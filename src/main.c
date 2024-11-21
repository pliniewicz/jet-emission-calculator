#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_integration.h>

#include "config.h"
#include "profiles.h"
#include "logspace.h"
#include "electron_distributions.h"
// #include "rad_functions.h"


int main(int argc, char *argv[])
{

int N = 250;


  double frequency[N];
  double first = 1e10;
  double last = 1e30;


  logspaced(first, last, N, frequency);

  double pB1 = 0, q = EQUIPART;
  pB1 = find_Bfield_normalization(pBint1, &q);

  printf("Normalized magnetic field pressure P_B: %e\n", pB1);

  double gmin = 1e2;
  double gmax = 1e6;
  struct normalizationParametersSimplePowerLaw params = {gmin, gmax, 1};

  double Ke = normalize_distribution(1, &params);


  printf("Normalization constant for the electron distribution: %e\n", Ke);

  return EXIT_SUCCESS;
}

