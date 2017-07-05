# Thermal functions

We provide a C++ library and Python interface to thermal functions, defined

    J_{B/F}(y^2) = Real{\int_0^\infnty dx x^2 \log(1 \mp \exp(-\sqrt{x^2 + y^2}))}
  
We offer Taylor expansion, numerical integration (quadrature), a Bessel function 
representation, an approximation and an upper bound for the integrals.

# Dependencies

The C++ requires gsl and gslcblas.

# Build

Build the library via 

    make
    
This should build `./lib/thermal_funcs.so`. The header file is `./src/thermal_funcs.h`. 

# Python interface

The Python interface 

    from thermal_funcs import J_B, J_F
    J_F(100., method='quad')
    
is compatible with Python 2 and 3, though must be built for a specific version. It has no dependencies. By default,
SWIG will build for your `python --version`. To alter this, change the PYTHON variable in the makefile to compile with
your chosen `Python.h` header.

# Results

The functions are fast and accurate even for `y^2 << 0`. See e.g. ![this plot](https://github.com/andrewfowlie/thermal_funcs/blob/master/test/J_B_neg.pdf). The time per evaluation is typically 1E-5 seconds or so.

