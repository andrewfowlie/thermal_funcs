#include "wstp.h"

double J_F_quad(double y_squared, double abs_error, double rel_error, int max_n);
double J_B_quad(double y_squared, double abs_error, double rel_error, int max_n);

int J_F(double y_squared) {
   return J_F_quad(y_squared, 1E-8, 1E-8, 10000);
}

int J_B(double y_squared) {
   return J_B_quad(y_squared, 1E-8, 1E-8, 10000);
}

int main(int argc, char *argv[]) {
   return WSMain(argc, argv);
}
