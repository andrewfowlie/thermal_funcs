/*

*/

#ifndef _THERMAL_FUNCS_H_
#define _THERMAL_FUNCS_H_


#ifdef __cplusplus
extern "C" {
#endif

double J_B_quad(double y_squared,
                double abs_error = 1E-7,
                double rel_error = 1E-7,
                int max_n = 10000);
double J_F_quad(double y_squared,
                double abs_error = 1E-7,
                double rel_error = 1E-7,
                int max_n = 10000);
double J_B_taylor(double y_squared,
                  double abs_error = 1E-7,
                  double rel_error = 1E-7,
                  int max_n = 10000);
double J_F_taylor(double y_squared,
                  double abs_error = 1E-7,
                  double rel_error = 1E-7,
                  int max_n = 10000);
double J_F_bessel(double y_squared,
                  double abs_error = 1E-7,
                  double rel_error = 1E-7,
                  int max_n = 10000,
                  bool fast = false);
double J_B_bessel(double y_squared,
                  double abs_error = 1E-7,
                  double rel_error = 1E-7,
                  int max_n = 10000,
                  bool fast = false);
double J_B_lim(double y_squared, bool upper = true);
double J_F_lim(double y_squared, bool upper = true);
double J_B_approx(double y_squared);
double J_F_approx(double y_squared);
double J_B_zeta(double y_squared, int n = 25);
double J_F_zeta(double y_squared, int n = 25);

#ifdef __cplusplus
}
#endif

#endif  // _THERMAL_FUNCS_H_
