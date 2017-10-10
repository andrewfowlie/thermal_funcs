/*
    Example program that links to thermal_funcs library. To build,
    
    make example
    
    builds an executable ./bin/example
*/

#include  <stdio.h>
#include <thermal_funcs.h>


int main() {
    printf("J_B = %e\n", J_B_bessel(100., 1E-7, 1E-7, 10000, false));
    return 0;
}
