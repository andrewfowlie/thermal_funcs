 %module thermal_funcs
 %{
  #include <thermal_funcs.h>
  #include <derivatives.h>
 %}
  
%feature("kwargs") J_B_quad; 
%feature("kwargs") J_F_quad; 
%feature("kwargs") J_B_taylor; 
%feature("kwargs") J_F_taylor;
%feature("kwargs") J_F_bessel;
%feature("kwargs") J_B_bessel;
%feature("kwargs") J_F_zeta;
%feature("kwargs") J_B_zeta;
%feature("kwargs") J_F_lim;
%feature("kwargs") J_B_lim;
%feature("kwargs") D1_J_B_bessel;
%feature("kwargs") D1_J_F_bessel;
%feature("kwargs") D2_J_B_bessel;
%feature("kwargs") D2_J_F_bessel;
%feature("kwargs") D1_J_B_approx;
%feature("kwargs") D1_J_F_approx;
%feature("kwargs") D2_J_B_approx;
%feature("kwargs") D2_J_F_approx;
%include <thermal_funcs.h>;
%include <derivatives.h>;
