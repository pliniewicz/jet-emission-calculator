#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_integration.h>

#include "config.h"
#include "profiles.h"
#include "logspace.h"
#include "rad_functions.h"


int main(int argc, char *argv[])
{

int N = 250;


  double frequency[N];
  double first = 1e10;
  double last = 1e30;


  logspaced(first, last, N, frequency);

  double pB1 = 0, q = EQUIPART;
  pB1 = find_Bfield_normalization(pBint2, &q);

  printf("%e\n", pB1);

  return EXIT_SUCCESS;
}

