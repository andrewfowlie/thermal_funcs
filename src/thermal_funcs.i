 %module thermal_funcs
 %{
  #include <thermal_funcs.h>
 %}
  
%feature("kwargs") J_B_quad; 
%feature("kwargs") J_F_quad; 
%feature("kwargs") J_B_taylor; 
%feature("kwargs") J_F_taylor;
%feature("kwargs") J_F_bessel;
%feature("kwargs") J_B_bessel;
%feature("kwargs") J_F_zeta;
%feature("kwargs") J_B_zeta;
%include <thermal_funcs.h>;
