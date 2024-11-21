#ifndef RAD_FUNCTIONS
#define RAD_FUNCTIONS

struct synchrotron_kernel_parameters {double frequency; double bulk_gamma; double magnetic_field; double minimal_gamma; double maximal_gamma; double angle;};
struct synchrotron_kernel_ball_parameters {double frequency; double bulk_gamma; double magnetic_field; double minimal_gamma; double maximal_gamma; double breaking_gamma; double k1; double k2; double sL; double sH; double angle;};
struct IC_kernel_parameters {double target_energy; double scattered_energy; double minimal_gamma; double maximal_gamma;};
struct IC_kernel_ball_parameters {double target_energy; double scattered_energy; double minimal_gamma; double maximal_gamma; double breaking_gamma; double k1; double sL; double sH;};

double doppler(double gamma, double theta);
double nu_critical(double gamma, double magnetic_field);
// struct SimplePowerLawParameters {double gmin; double gmax; double index;};
// double simple_power_law(double gamma, void *params);
double simple_power_law_2(double gamma, double minimal_gamma, double maximal_gamma);
double broken_power_law(double gamma, double minimal_gamma, double maximal_gamma, double breaking_gamma, double k1, double sL, double sH);
double exp_cutoff_power_law_2(double gamma, double minimal_gamma, double maximal_gamma);
double exp_cutoff_power_law(double gamma, double minimal_gamma, double maximal_gamma, double power_law_index);
double broken_exp_cutoff_power_law(double gamma);

double synchrotron_kernel(double gamma, void *parameters);
double synchrotron_kernel_ball(double gamma, void *parameters);

double qu(double scattered_energy, double target_energy, double gamma);
double Qu(double target_energy, double gamma);

double f_IC(double target_energy, double scattered_energy, double gamma);

double inverse_compton_kernel(double gamma, void *parameters);
double inverse_compton_ball_kernel(double gamma, void *parameters);

#endif // !RAD_FUNCTIONS
