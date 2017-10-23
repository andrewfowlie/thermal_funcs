/*

*/

#ifndef _DERIVATIVES_H_
#define _DERIVATIVES_H_


#ifdef __cplusplus
#define DEFAULT(x) = x
#else
#define DEFAULT(x)
#endif


#ifdef __cplusplus
extern "C" {
#endif


double D1_J_F_bessel(double y_squared,
                     double abs_error DEFAULT(1E-7),
                     double rel_error DEFAULT(1E-7),
                     int max_n DEFAULT(10000));
double D1_J_B_bessel(double y_squared,
                     double abs_error DEFAULT(1E-7),
                     double rel_error DEFAULT(1E-7),
                     int max_n DEFAULT(10000));

double D2_J_F_bessel(double y_squared,
                     double abs_error DEFAULT(1E-7),
                     double rel_error DEFAULT(1E-7),
                     int max_n DEFAULT(10000));
double D2_J_B_bessel(double y_squared,
                     double abs_error DEFAULT(1E-7),
                     double rel_error DEFAULT(1E-7),
                     int max_n DEFAULT(10000));
                     
double D1_J_F_approx(double y_squared, double step DEFAULT(1E-2));
double D1_J_B_approx(double y_squared, double step DEFAULT(1E-2));

double D2_J_F_approx(double y_squared, double step DEFAULT(1E-1));
double D2_J_B_approx(double y_squared, double step DEFAULT(1E-1));

#ifdef __cplusplus
}
#endif

#endif  // _DERIVATIVES_H_
