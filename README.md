# Thermal functions

We provide a C++ library and Python interface to thermal functions, defined

    J_{B/F}(y^2)=\Re\int_0^{\infty}dx\,x^2 \ln(1\mp\exp(-\sqrt{x^2+y^2}))
    
<p align="center">
  <img src="https://latex.codecogs.com/png.latex?J_{B/F}(y^2)=\Re\int_0^{\infty}dx\,x^2&space;\ln(1\mp\exp(-\sqrt{x^2&plus;y^2}))"/>
</p>    

  
We offer Taylor expansion, numerical integration (quadrature), a Bessel function 
representation, an approximation, a Hurwitz zeta function representation, and an upper bound for the integrals.

# Dependencies

The C++ requires `gsl` and `gslcblas`.

# Build

Build the library via 

    make
    
This should build `./lib/thermal_funcs.so`. The header file is `./src/thermal_funcs.h`. 

# Python interface

The Python interface 

    from thermal_funcs import J_B, J_F
    J_F(100., method='quad')
    
is compatible with Python 2 and 3, though must be built for a specific version. It has no dependencies. By default,
SWIG will build for your `python --version`. To alter this, change the `PYTHON` variable in the makefile to compile with
your chosen `Python.h` header.

# Performance

The functions are fast and accurate even for `y^2 << 0`:

<p align="center">
  <img src="https://user-images.githubusercontent.com/3758193/27900262-fad8eaf8-6270-11e7-8324-4e745fd04301.png"/>
</p>   

The time per evaluation is typically about `1E-5` seconds.

