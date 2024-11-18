#ifndef PROFILES
#define PROFILES

double Gamma(double x);
double b1(double x);
double b2(double x);
double b3(double x);
double f1(double x);
double f2(double x);
double f3(double x);
double p1_integrand(double s, void * parameters);
double p2_integrand(double s, void * parameters);
double p3_integrand(double s, void * parameters);
double p1(double x, double q);
double p2(double x, double q);
double p3(double x, double q);
double relativistic_beta(double x);
double doppler_profile(double x, double theta);
double pBint1(double x, void * params);
double pBint2(double x, void * params);
double pBint3(double x, void * params);

#endif // !PROFILES
