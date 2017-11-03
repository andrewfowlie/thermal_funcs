# Thermal functions

We provide a C++ library and Python and Mathematica interfaces to thermal functions, defined

    J_{B/F}(y^2)=\Re\int_0^{\infty}dx\,x^2 \ln(1\mp\exp(-\sqrt{x^2+y^2}))
    
<p align="center">
  <img src="https://latex.codecogs.com/png.latex?J_{B/F}(y^2)=\Re\int_0^{\infty}dx\,x^2&space;\ln(1\mp\exp(-\sqrt{x^2&plus;y^2}))"/>
</p>    

  
We offer Taylor expansion, numerical integration (quadrature), a Bessel function 
representation, an approximation, a Hurwitz zeta function representation, and an upper bound for the integrals. First and second
derivatives are also implemented.

# Dependencies

The C++ requires `gsl` and `gslcblas`. The Python interface requires Python 2 or 3, SWIG and a Python.h header file (which is part of python-dev in Ubuntu). The Mathematica interface was tested for Mathematica 11.

# Build

Build the library via 

    make lib
    
This should build `./lib/thermal_funcs.so`. The header file is `./src/thermal_funcs.h`. 

# Example

There is a C example ./src/example.cpp, built by

    make example
    
This should build a program `./bin/example`, which when executed prints the result of evaluating a thermal function.

# Python interface

Build the interface via 

    make python

The interface 

    from thermal_funcs import J_B, J_F
    J_F(100., method='quad')
    
is compatible with Python 2 and 3, though must be built for a specific version. It has no module dependencies. By default,
SWIG will build for your `python --version`. To alter this, change the `PYTHON` variable in the makefile to compile with
your chosen `Python.h` header. The derivatives are called by a keyword argument e.g., `J_F(100., derivative=1)`.

# Mathematica interface

This is slightly more involved. This may work in Linux if `math` is in your `PATH`:

    make mathematica
    
But otherwise you may have to tweak the `./src/makefile` variable `MATH_INC` for the locations of your `wscc` linker and `wstp.h` header file. You can find this on any platform in Mathematica from `FileNameJoin[{$InstallationDirectory, "SystemFiles", "Links", "WSTP", "DeveloperKit", $SystemID, "CompilerAdditions"}]`. 

Then within Mathematica,

    Install["./src/math.exe"];
    Plot[{JB[ysq], JF[ysq]}, {ysq, -100, 100}]
    
Note well that you should use the correct (relative or absolute) path to `./src/math.exe` in the command `Install["./src/math.exe"]`. The interface was built and tested with Mathematica 11.1.1. The derivatives are called by a keyword argument e.g., `JB[100., derivative->1]`.

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

# Files
<pre>
thermal_funcs/
|-- bin : make example builds example program here (intentionally empty)
|-- lib : make lib builds thermal_funcs.so here (intentionally empty)
|-- src : Source code
|   |-- bernoulli.h : C++ header containing Bernoulli numbers 
|   |-- derivatives.cpp : C++ source for derivatives of thermal functions
|   |-- derivatives.h : C/C++ header for derivatives of thermal functions
|   |-- example.cpp : C++ example program
|   |-- makefile : Makefile for building source code
|   |-- math.c : C source for Mathematica interface
|   |-- math.tm : Mathematica template for interface
|   |-- thermal_funcs.cpp : C++ source for thermal functions
|   |-- thermal_funcs.h : C/C++ header for thermal functions
|   |-- thermal_funcs.i : SWIG Python interface for thermal functions
|   |-- zeta.cpp : C++ source for Hurwitz zeta function
|   `-- zeta.h : C++ header for Hurwitz zeta function
|-- thermal_funcs : Python module
|   |-- __init__.py : Python module file
|   `-- thermal_funcs.py : Python source for wrapper for thermal functions
|-- LICENSE : BSD 3-Clause License
|-- makefile : Wrapper for src/makefile
`-- README.md : This README file
</pre>
# Acknowledgements

This [Stack Exchange answer](https://mathematica.stackexchange.com/a/154643/38645) was helpful for removing linker warnings from `wscc`, and [this one](https://mathematica.stackexchange.com/a/154664/38645) was helpful for automatically locating Mathematica header files. 
