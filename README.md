# Thermal functions

We provide a C++ library and Python interface to thermal functions, defined

    J_{B/F}(y^2) = \int_0^\infnty dx x^2 \log(1 \mp \exp(-\sqrt{x^2 + y^2}))
  
We offer Taylor expansion, numerical integration (quadrature), and a Bessel function 
representation of the integrals.

# Dependencies

The C++ requires gsl and gslcblas.

# Build

Build the library via 

    make
    
THe *.so should be built in the ./lib/ directory.

# Python interface

The Python interface 

    from thermal_funcs import J_B, J_F
    J_F(100., method='quad')
    
is compatible with Python 2 and 3. It requires CppHeaderParser, future, 
and functools32 for Python 2, all of which may be installed via pip.
