module smb_const

  use glimmer_global

  implicit none

  ! Fundamental constants:
  real(rk), parameter :: pi     = 3.1415926535897
  real(rk), parameter :: grav   = 9.81             ! m s-2
  real(rk), parameter :: gascda = 287.0            ! gas const for dry air J kg-1 K-1
  real(rk), parameter :: gascwv = 461.89           ! gas const for water vapour J kg-1 K-1
  real(rk), parameter :: tkel   = 273.15           ! to get temperature in Kelvin
  real(rk), parameter :: sb     = 5.67e-08         ! stefan-boltzmann constant, W m-2 K-4
  real(rk), parameter :: cp     = 2009.0           ! spec heat capacity ICE J kg-1 K-1
  real(rk), parameter :: cpwat  = 4200.0           ! spec heat capacity WATER J kg-1 K-1
  real(rk), parameter :: fus    = 0.335e06         ! latent heat fusion J kg-1
  real(rk), parameter :: vap    = 2.5e06           ! latent heat of sublimation J kg-1
  real(rk), parameter :: densn  = 300.0            ! snow density kg m-3
  real(rk), parameter :: denice = 910.0            ! ice density  kg m-3
  real(rk), parameter :: denwat = 1000.0           ! water density kg m-3
  real(rk), parameter :: cpair  = 1005.0           ! spec. heat capacity AIR J kg-1 K-1
  real(rk), parameter :: daymin = 24.0*60.0        ! minutes in one day  
  ! Other constants:
  integer,  parameter :: ndyear = 365              ! 365.25 for glimmer, 360 for IGCM ??!!
  real(rk), parameter :: epsden = 1.0              ! epsilon for density difference
  real(rk), parameter :: tav    = 10.0             !
  real(rk), parameter :: accden = 1000.0           !
  real(rk), parameter :: denz0  = 600.0            ! parameterisation of roughness lengths of snow
  real(rk), parameter :: denlow = 300.0            ! parameterisation of roughness lengths of snow  
  real(rk), parameter :: blc    = 0.5              ! to calculate and maintain the vertical grid
  real(rk), parameter :: bhc    = 2.0              ! to calculate and maintain the vertical grid

end module smb_const
