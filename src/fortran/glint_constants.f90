module glint_constants

  use glimmer_global

  implicit none

  ! ------------------------------------------------------------
  ! global parameters
  ! ------------------------------------------------------------

  real(rk),parameter :: pi=3.141592654          !*FD The value of pi
  real(rk),parameter :: days2hours=24.0
  real(rk),parameter :: hours2seconds=3600.0    !*FD Hours to seconds conversion factor

  integer, parameter :: default_diy=360                    !*FD Default number of days in year
  integer, parameter :: default_y2h=days2hours*default_diy !*FD Default years to hours conversion

  ! Constants set at run-time

  integer  :: days_in_year=default_diy        !*FD The number of days in a year  
  real(rk) :: years2hours =default_y2h        !*FD Years to hours conversion factor
  real(rk) :: hours2years =1.0_rk/default_y2h !*FD Hours to years conversion factor

  private :: default_diy,default_y2h

contains

  subroutine glint_set_year_length(daysinyear)

    integer, intent(in) :: daysinyear

    days_in_year=daysinyear
    years2hours=days2hours*days_in_year 
    hours2years=1.0_rk/years2hours      

  end subroutine glint_set_year_length

end module glint_constants
