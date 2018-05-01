/**
    @file
    @brief Example program that links to the `thermal_funcs` library  in `C`.
    
    To build, `make example` builds an executable `./bin/example`.
    
    @example example.cpp
*/

#include <stdio.h>
#include <thermal_funcs.h>
#include <derivatives.h>


int main() {
  /**
      This is an example of a call to the `thermal_funcs` library  in `C`.
  */
  printf("J_B = %e\n", J_B_bessel(100.));
  printf("D1_J_B_bessel = %e\n", D1_J_B_bessel(100.));
  printf("D2_J_B_bessel = %e\n", D2_J_B_bessel(100.));
  return 0;
}
