/*

*/

#ifndef _THERMAL_FUNCS_H_
#define _THERMAL_FUNCS_H_


#ifdef __cplusplus
#define DEFAULT(x) = x
#else
#include <stdbool.h>
#define DEFAULT(x)
#endif


#ifdef __cplusplus
extern "C" {
#endif

double J_B_quad(double y_squared,
                double abs_error DEFAULT(1E-7),
                double rel_error DEFAULT(1E-7),
                int max_n DEFAULT(10000));
double J_F_quad(double y_squared,
                double abs_error DEFAULT(1E-7),
                double rel_error DEFAULT(1E-7),
                int max_n DEFAULT(10000));
double J_B_taylor(double y_squared,
                  double abs_error DEFAULT(1E-7),
                  double rel_error DEFAULT(1E-7),
                  int max_n DEFAULT(10000));
double J_F_taylor(double y_squared,
                  double abs_error DEFAULT(1E-7),
                  double rel_error DEFAULT(1E-7),
                  int max_n DEFAULT(10000));
double J_F_bessel(double y_squared,
                  double abs_error DEFAULT(1E-7),
                  double rel_error DEFAULT(1E-7),
                  int max_n DEFAULT(10000),
                  bool fast DEFAULT(false));
double J_B_bessel(double y_squared,
                  double abs_error DEFAULT(1E-7),
                  double rel_error DEFAULT(1E-7),
                  int max_n DEFAULT(10000),
                  bool fast DEFAULT(false));
double J_B_lim(double y_squared, bool upper DEFAULT(true));
double J_F_lim(double y_squared, bool upper DEFAULT(true));
double J_B_approx(double y_squared);
double J_F_approx(double y_squared);
double J_B_zeta(double y_squared, int n DEFAULT(25));
double J_F_zeta(double y_squared, int n DEFAULT(25));

#ifdef __cplusplus
}
#endif

#endif  // _THERMAL_FUNCS_H_
