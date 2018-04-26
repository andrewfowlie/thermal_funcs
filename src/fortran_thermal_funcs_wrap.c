/*
    Fortran wrapper for thermal_funcs
*/

#include <thermal_funcs.h>

double j_b(double *y_squared) {
    return J_B_bessel(*y_squared, 1E-7, 1E-7, 10000, false);
}
