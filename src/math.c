#include <stdbool.h>
#include <wstp.h>
#include <thermal_funcs.h>

double J_F(double y_squared) {
   return J_F_bessel(y_squared, 1E-7, 1E-7, 10000, false);
}

double J_B(double y_squared) {
   return J_B_bessel(y_squared, 1E-7, 1E-7, 10000, false);
}

int main(int argc, char *argv[]) {
   return WSMain(argc, argv);
}
