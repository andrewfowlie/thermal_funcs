# Thermal functions

We provide a C++ library and Python and Mathematica interfaces to thermal functions, defined

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
    
This should build `./lib/thermal_funcs.so`. The header file is `./src/thermal_funcs.h`. It should also build a Python interface, but not a Mathematica interface. 
If you only want the C++ library,

    make thermal_funcs.so 

# Python interface

The Python interface 

    from thermal_funcs import J_B, J_F
    J_F(100., method='quad')
    
is compatible with Python 2 and 3, though must be built for a specific version. It has no dependencies. By default,
SWIG will build for your `python --version`. To alter this, change the `PYTHON` variable in the makefile to compile with
your chosen `Python.h` header.

# Mathematica interface

This is slightly more involved. This may work:

    make math.exe
    
But you may have to tweak the `./src/makefile` variables `MATH` and `MATH_INC` for the locations of your `wscc` linker and `wstp.h` header file. Try `locate wscc` and `locate wstp.h`. They should be somewhere within `GetEnvironment[{"MATHEMATICA_HOME"}]`.

Then within Mathematica,

    Install["./src/math.exe"];
    Plot[{JB[ysq], JF[ysq]}, {ysq, -100, 100}]
    
Note well that you should use the correct (relative or absolute) path to `./src/math.exe` in the command `Install["./src/math.exe"]`.

## Debugging

If the executable `./src/math.exe` was built but `Install` fails, try installing step by step to find debugging information. First, run the created executable,

    ./src/math.exe
    
This should prompt you to `Create link:`. Enter e.g. `foo`. Don't exit that session. In Mathematica, try

    $VersionNumber
    link = LinkConnect["foo"]
    Install[link]
    JB[100]
    
to find the step that fails.

You can also try one of the pre-built examples provided by Mathematica, e.g.,

    Install["/usr/local/Wolfram/Mathematica/11.1/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64/PrebuiltExamples/addtwo"]
    AddTwo[2, 2]
    
and re-building it locally,
    
    MATH=/usr/local/Wolfram/Mathematica/11.1/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64/
    mkdir ~/addtwo
    cd ~/addtwo
    cp $MATH/WSTPExamples/addtwo* ./
    $MATH/CompilerAdditions/wscc addtwo.tm addtwo.c -o addtwo
    
then in Mathematica,

    Install["~/addtwo/addtwo"]
    AddTwo[2, 2]
    
This may help find the origin of any problems. You must, of course, replace the paths to the ones on your machine.

# Performance

The functions are fast and accurate even for `y^2 << 0`:

<p align="center">
  <img src="https://user-images.githubusercontent.com/3758193/27900262-fad8eaf8-6270-11e7-8324-4e745fd04301.png"/>
</p>

The time per evaluation is typically about `1E-5` seconds.

# Acknowledgements

This [Stack Exchange answer](https://mathematica.stackexchange.com/a/154643/38645) was helpful for removing linker warnings from `wscc`.
