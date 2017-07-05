 %module thermal_funcs
 %{
 #define SWIG_FILE_WITH_INIT
 #include <thermal_funcs.h>
 %}
 
%include "numpy.i"

%init %{
import_array();
%}
 
%feature("kwargs") J_B_quad; 
extern "C" double J_B_quad(double y_squared,
                           double abs_error = 1E-8,
                           double rel_error = 1E-8,
                           int max_n = 10000);
%feature("kwargs") J_F_quad; 
extern "C" double J_F_quad(double y_squared,
                           double abs_error = 1E-8,
                           double rel_error = 1E-8,
                           int max_n = 10000);
%feature("kwargs") J_B_taylor; 
extern "C" double J_B_taylor(double y_squared,
                             double abs_error = 1E-8,
                             double rel_error = 1E-8,
                             int max_n = 10000);
%feature("kwargs") J_F_taylor;
extern "C" double J_F_taylor(double y_squared,
                             double abs_error = 1E-8,
                             double rel_error = 1E-8,
                             int max_n = 10000);
%feature("kwargs") J_F_bessel;
extern "C" double J_F_bessel(double y_squared,
                             double abs_error = 1E-8,
                             double rel_error = 1E-8,
                             int max_n = 10000,
                             bool fast = false);
%feature("kwargs") J_B_bessel;
extern "C" double J_B_bessel(double y_squared,
                             double abs_error = 1E-8,
                             double rel_error = 1E-8,
                             int max_n = 10000,
                             bool fast = false);
extern "C" double J_B_lim(double y_squared);
extern "C" double J_F_lim(double y_squared);
extern "C" double J_B_approx(double y_squared);
extern "C" double J_F_approx(double y_squared);
