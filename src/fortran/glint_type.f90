! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_type.f90 - part of the GLIMMER ice model           + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004 GLIMMER contributors - see COPYRIGHT file 
! for list of contributors.
!
! This program is free software; you can redistribute it and/or 
! modify it under the terms of the GNU General Public License as 
! published by the Free Software Foundation; either version 2 of 
! the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License 
! along with this program; if not, write to the Free Software 
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 
! 02111-1307 USA
!
! GLIMMER is maintained by:
!
! Ian Rutt
! School of Geographical Sciences
! University of Bristol
! University Road
! Bristol
! BS8 1SS
! UK
!
! email: <i.c.rutt@bristol.ac.uk> or <ian.rutt@physics.org>
!
! GLIMMER is hosted on NeSCForge:
!
! http://forge.nesc.ac.uk/projects/glimmer/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module glint_type
  !*FD contains type definitions for GLINT
  use glimmer_global
  use glint_proj
  use glint_interp
  use glimmer_degd
  use glide_types

  implicit none

  type glint_forcdata
     !*FD Holds information relating to time-series forcing.
     character(fname_length) :: forcfile = ''             !*FD Temperature forcing file
     real(sp),dimension(:,:),pointer :: forcing => null() !*FD Forcing data(?)
     integer  :: flines = 0
     real(sp) :: trun   = 10000
  end type glint_forcdata


  type glint_climate
     !*FD Derived type holding climate variables for GLINT

     integer :: whichartm = 3
     !*FD Method of calculation of surface air temperature:
     !*FD \begin{description} 
     !*FD \item[0] Linear function of surface elevation 
     !*FD \item[1] Cubic function of distance from domain centre
     !*FD \item[2] Linear function of distance from domain centre
     !*FD \item[3] Greenland conditions (function of surface elevation 
     !*FD and latitude) including forcing
     !*FD \item[4] Antarctic conditions (sea-level air temperature 
     !*FD -- function of position)
     !*FD \item[5] Uniform temperature, zero range (temperature set in
     !*FD \texttt{cons} namelist) 
     !*FD \item[6] Uniform temperature, corrected for height, zero range.
     !*FD \item[7] Use large-scale temperature and range. 
     !*FD \end{description}

     integer :: whichprecip = 0
     !*FD Source of precipitation:
     !*FD \begin{description}
     !*FD \item[0] Uniform precipitation rate (set internally 
     !*FD at present) 
     !*FD \item[1] Use large-scale precipitation rate 
     !*FD \item[2] Use parameterization of \emph{Roe and Lindzen} 
     !*FD \end{description}         

     integer :: whichacab = 2
     !*FD Net accumulation: 
     !*FD \begin{description}
     !*FD \item[0] EISMINT moving margin 
     !*FD \item[1] PDD mass-balance model [recommended] 
     !*FD \item[2] Accumulation only 
     !*FD \end{description}

     real(sp),dimension(:,:),pointer :: arng     => null() !*FD Annual air temperature range
     real(sp),dimension(:,:),pointer :: ablt     => null() !*FD Annual ablation.
     real(sp),dimension(:,:),pointer :: g_arng   => null() !*FD Global annual air temperature range
     real(sp),dimension(:,:),pointer :: g_artm   => null() !*FD Global annual mean air temperature
     real(sp),dimension(:,:),pointer :: prcp     => null() !*FD Annual precipitation.
     real(sp),dimension(:,:),pointer :: presprcp => null() !*FD Present-day annual precip (mm)
     real(sp),dimension(:,:),pointer :: presartm => null() !*FD Present-day mean air-temp ($^{\circ}$C)
     real(dp),dimension(:,:),pointer :: presusrf => null() !*FD Present-day upper surface (km)
     integer, dimension(:,:),pointer :: out_mask => null() !*FD Array indicating whether a point 
                                                           !*FD should be considered or ignored 
                                                           !*FD when upscaling data for output. 
                                                           !*FD 1 means use, 0 means ignore.
     real(sp) :: uprecip_rate =   0.5 !*FD Uniform precipitaion rate in m/a
     real(sp) :: usurftemp    = -20.0 !*FD Uniform surface temperature in $^{\circ}$C.
     real(sp) :: ice_albedo   =   0.4 !*FD Ice albedo. (fraction)
     real(sp) :: ulapse_rate  =  -8.0 !*FD Uniform lapse rate in deg C/km 
                                      !*FD (N.B. This should be \emph{negative}!)
     real(sp),dimension(2) :: airt = (/ -3.150, -1.0e-2 /)       ! K, K km^{-3}
     real(sp),dimension(3) :: nmsb = (/ 0.5, 1.05e-5, 450.0e3 /) ! m yr^{-1}, yr^{-1}, m                
  end type glint_climate

  type glint_instance
     !*FD Derived type holding information about ice model instance. 
     type(projection)                 :: proj               !*FD The projection definition of the instance.
     type(downscale)                  :: downs              !*FD Downscaling parameters.
     type(upscale)                    :: ups                !*FD Upscaling parameters
     type(upscale)                    :: ups_orog           !*FD Upscaling parameters for orography (to cope
     !*FD with need to convert to spectral form).
     type(glide_global_type)          :: model              !*FD The instance and all its arrays.
     type(glint_forcdata)             :: forcdata           !*FD temperature forcing
     type(glint_climate)              :: climate            !*FD climate data
     type(glimmer_pddcalc)            :: pddcalc            !*FD positive degree-day data
     character(fname_length)          :: paramfile          !*FD The name of the configuration file.
     logical                          :: newtemps           !*FD Flag to say we have new temperatures.
     logical                          :: first     = .true. !*FD Is this the first timestep?

     ! Arrays to hold downscaled versions of input data --------------------------
     real(rk), dimension(:,:),pointer :: xwind         => null() 
     !*FD $x$-component of surface winds on local grid.

     real(rk), dimension(:,:),pointer :: ywind         => null() 
     !*FD $y$-component of surface winds on local grid.

     real(dp), dimension(:,:),pointer :: global_orog   => null() 
     !*FD Global orography on local coordinates.

     real(rk), dimension(:,:),pointer :: local_orog    => null() 
     !*FD Local orography on local coordinates.

     ! Fractional coverage information ------------------------------------------- 
     real(rk) ,dimension(:,:),pointer :: frac_coverage => null() 
     !*FD Fractional coverage of each global gridbox by the projected grid.

     real(rk) ,dimension(:,:),pointer :: frac_cov_orog => null() 
     !*FD Fractional coverage of each global gridbox by the projected grid (orography).
  end type glint_instance

  type output_flags
     !*FD A derived type used internally to communicate the outputs which need
     !*FD to be upscaled, thus avoiding unnecessary calculation

     logical :: orog         !*FD Set if we need to upscale the orography
     logical :: albedo       !*FD Set if we need to upscale the albedo
     logical :: ice_frac     !*FD Set if we need to upscale the ice fraction
     logical :: water_in     !*FD Set if we need to upscale the input water flux
     logical :: water_out    !*FD Set if we need to upscale the output water flux
     logical :: total_win    !*FD Set if we need to sum the total water taken up by ice sheet
     logical :: total_wout   !*FD Set if we need to sum the total ablation by the ice sheet
     logical :: ice_vol      !*FD Set if we need to calculate the total ice volume
  end type output_flags

contains
    subroutine glint_i_allocate(instance,nxg,nyg)

    !*FD Allocate top-level arrays in
    !*FD the model instance, and ice model arrays.

    implicit none

    type(glint_instance),intent(inout) :: instance    !*FD Instance whose elements are to be allocated.
    integer,               intent(in)    :: nxg         !*FD Longitudinal size of global grid (grid-points).
    integer,               intent(in)    :: nyg         !*FD Latitudinal size of global grid (grid-points).

    integer ewn,nsn

    ewn=instance%model%general%ewn
    nsn=instance%model%general%nsn

    ! First deallocate if necessary

    if (associated(instance%xwind))         deallocate(instance%xwind)
    if (associated(instance%ywind))         deallocate(instance%ywind)
    if (associated(instance%global_orog))   deallocate(instance%global_orog) 
    if (associated(instance%local_orog))    deallocate(instance%local_orog)   
    if (associated(instance%frac_coverage)) deallocate(instance%frac_coverage)

    ! Then reallocate...
    ! Wind field arrays

    allocate(instance%xwind(ewn,nsn)) 
    allocate(instance%ywind(ewn,nsn))

    ! Local-global orog

    allocate(instance%global_orog(ewn,nsn))

    ! Local-local orog

    allocate(instance%local_orog(ewn,nsn))

    ! Global box indices and number of points contained therein

    ! Fractional coverage map

    allocate(instance%frac_coverage(nxg,nyg)) ! allocate fractional coverage map
    allocate(instance%frac_cov_orog(nxg,nyg)) ! allocate fractional coverage map (orog)

    ! allocate climate data
    allocate(instance%climate%arng(ewn,nsn));            instance%climate%arng = 0.0
    allocate(instance%climate%g_artm(ewn,nsn));          instance%climate%g_artm = 0.0
    allocate(instance%climate%g_arng(ewn,nsn));          instance%climate%g_arng = 0.0
    allocate(instance%climate%out_mask(ewn,nsn));        instance%climate%out_mask = 1.0
    allocate(instance%climate%ablt(ewn,nsn));            instance%climate%ablt = 0.0
    allocate(instance%climate%prcp(ewn,nsn));            instance%climate%prcp = 0.0
    allocate(instance%climate%presprcp(ewn,nsn));        instance%climate%presprcp = 0.0
    allocate(instance%climate%presartm(ewn,nsn));        instance%climate%presartm = 0.0
    allocate(instance%climate%presusrf(ewn,nsn));        instance%climate%presusrf = 0.0    

  end subroutine glint_i_allocate
end module glint_type
