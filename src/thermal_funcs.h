/*

*/

#ifndef _THERMAL_FUNCS_H_
#define _THERMAL_FUNCS_H_

extern "C" double J_B_quad(double y_squared,
                           double abs_error = 1E-8,
                           double rel_error = 1E-8,
                           int max_n = 10000);
extern "C" double J_F_quad(double y_squared,
                           double abs_error = 1E-8,
                           double rel_error = 1E-8,
                           int max_n = 10000);
extern "C" double J_B_taylor(double y_squared,
                             double abs_error = 1E-8,
                             double rel_error = 1E-8,
                             int max_n = 10000);
extern "C" double J_F_taylor(double y_squared,
                             double abs_error = 1E-8,
                             double rel_error = 1E-8,
                             int max_n = 10000);
extern "C" double J_F_bessel(double y_squared,
                             double abs_error = 1E-8,
                             double rel_error = 1E-8,
                             int max_n = 10000,
                             bool fast = false);
extern "C" double J_B_bessel(double y_squared,
                             double abs_error = 1E-8,
                             double rel_error = 1E-8,
                             int max_n = 10000,
                             bool fast = false);
extern "C" double J_B_lim(double y_squared);
extern "C" double J_F_lim(double y_squared);
extern "C" double J_B_approx(double y_squared);
extern "C" double J_F_approx(double y_squared);
#endif  // _THERMAL_FUNCS_H_
