! @file
! @example fortran_example.f
! @brief Fortran interface to `thermal_funcs`
! Build with `make fortran` to create `./bin/fortran_example`

program thermal_funcs

  use iso_c_binding

  implicit none
  
  interface
    function j_b(y_squared) bind (c)
      use iso_c_binding
      real(c_double) :: y_squared
      real(c_double) :: j_b
    end function
  end interface
  
  real(c_double) :: y_squared = 10. 
  real(c_double) :: r 
  
  r = j_b(y_squared)  ! This calls the c interface
  write(*, *) "j_b(", y_squared, ") = ", r

end program
