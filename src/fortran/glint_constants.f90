module glint_constants

  use glimmer_global

  implicit none

  ! ------------------------------------------------------------
  ! global parameters/constants
  ! ------------------------------------------------------------

  integer, parameter :: days_in_year=360                   !*FD The number of days in a year  
  real(rk),parameter :: pi=3.141592654                     !*FD The value of pi
  real(rk),parameter :: days2hours=24.0
  real(rk),parameter :: years2hours=days2hours*days_in_year !*FD Years to hours conversion factor
  real(rk),parameter :: hours2years=1/years2hours          !*FD Hours to years conversion factor
  real(rk),parameter :: hours2seconds=3600.0               !*FD Hours to seconds conversion factor

end module glint_constants
