/*
    Example program that links to thermal_funcs library. To build,
    
    make example
    
    builds an executable ./bin/example
*/

#include <stdio.h>
#include <thermal_funcs.h>
#include <derivatives.h>


int main() {
    printf("J_B = %e\n", J_B_bessel(100.));
    printf("D1_J_B_bessel = %e\n", D1_J_B_bessel(100.));
    printf("D2_J_B_bessel = %e\n", D2_J_B_bessel(100.));
    return 0;
}
