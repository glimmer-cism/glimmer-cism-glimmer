! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_example.f90 - part of the GLIMMER ice model      + 
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

program glint_example

  !*FD This program demonstrates the use of GLIMMER. It loads in
  !*FD some example global fields and associated grid data,
  !*FD Initialises the model, and then runs it for 1000 years.

  use glint_main
  use glimmer_log
  use glint_global_interp
  use glint_example_clim
  implicit none

  ! Program variables -------------------------------------------------------------------

  type(glint_params)      :: ice_sheet   ! This is the derived type variable that holds all 
                                         ! domains of the ice model
  type(glex_climate)      :: climate     ! Climate parameters and fields
  character(fname_length) :: paramfile   ! Name of the top-level configuration file
  character(fname_length) :: climatefile ! Name of climate configuration file

  ! Arrays which hold the global fields used as input to GLIMMER ------------------------

  real(rk),dimension(:,:),allocatable :: temp      ! Temperature     (degC)
  real(rk),dimension(:,:),allocatable :: precip    ! Precipitation   (mm/s)
  real(rk),dimension(:,:),allocatable :: zonwind   ! Zonal wind      (m/s)
  real(rk),dimension(:,:),allocatable :: merwind   ! Meridional wind (m/s)
  real(rk),dimension(:,:),allocatable :: orog      ! Orography       (m)

  ! Arrays which hold information about the ice model instances -------------------------

  real(rk),dimension(:,:),allocatable :: coverage ! Coverage map for normal global grid
  real(rk),dimension(:,:),allocatable :: cov_orog ! Coverage map for orography grid

  ! Arrays which hold output from the model ---------------------------------------------
  ! These are all on the normal global grid, except for the orography

  real(rk),dimension(:,:),allocatable :: albedo   ! Fractional albedo
  real(rk),dimension(:,:),allocatable :: orog_out ! Output orography (m)
  real(rk),dimension(:,:),allocatable :: ice_frac ! Ice coverage fraction
  real(rk),dimension(:,:),allocatable :: fw       ! Freshwater output flux (mm/s)
  real(rk),dimension(:,:),allocatable :: fw_in    ! Freshwater input flux (mm/s)

  ! Arrays which hold information about the global grid ---------------------------------

  real(rk),dimension(:),  allocatable :: lats_orog      ! Latitudes of global orography gridpoints
  real(rk),dimension(:),  allocatable :: lons_orog      ! Longitudes of global oropraphy gridpoints

  ! Scalars which hold information about the global grid --------------------------------

  integer :: nx,ny   ! Size of normal global grid
  integer :: nxo,nyo ! Size of global orography grid

  ! Scalar model outputs ----------------------------------------------------------------

  real(rk) :: twin     ! Timestep-integrated input water flux (kg)
  real(rk) :: twout    ! Timestep-integrated output water flux (kg)
  real(rk) :: ice_vol  ! Total ice volume (m^3)

  ! Other variables ---------------------------------------------------------------------

  logical :: out    ! Outputs set flag
  integer :: i,j    ! Array index counters
  integer :: time   ! Current time (hours)

  ! -------------------------------------------------------------------------------------
  ! Executable code starts here - Basic initialisation
  ! -------------------------------------------------------------------------------------

#ifdef GLEX_COM_LINE

  ! Non-f95-standard command-line interface using Intel Compiler
  ! features - possibly portable to other compilers, but untested

  if (nargs().eq.3) then
     call getarg(1,climatefile)
     call getarg(2,paramfile)
     Print*,'Using climate configuration: ',climatefile
     Print*,'Using ice-model configuration: ',paramfile
  else

#endif

     ! These are the default inputs
     Print*,'Enter name of climate configuration file:'
     read*,climatefile
     Print*,'Enter name of ice model configuration file:'
     read*,paramfile

#ifdef GLEX_COM_LINE
  endif
#endif

  ! Initialise climate

  call glex_clim_init(climate,climatefile)

  ! Set dimensions of global grids

  call get_grid_dims(climate%all_grid,nx,ny) ! Normal global grid
  nxo=200 ; nyo=100                          ! Example grid used for orographic output

  ! start logging
  call open_log(unit=101)  

  ! Allocate arrays appropriately

  allocate(temp(nx,ny),precip(nx,ny),zonwind(nx,ny))
  allocate(merwind(nx,ny),orog(nx,ny))
  allocate(coverage(nx,ny),orog_out(nxo,nyo),albedo(nx,ny),ice_frac(nx,ny),fw(nx,ny))
  allocate(lats_orog(nyo),lons_orog(nxo),cov_orog(nxo,nyo),fw_in(nx,ny))

  ! Initialise array contents

  temp=0.0
  precip=0.0
  zonwind=0.0
  merwind=0.0
  albedo=0.0
  orog_out=0.0
  orog=real(climate%orog_clim)                    ! Put orography where it belongs

  ! Set up global grids ----------------------------------------------------------------

  ! Calculate example orographic latitudes

  do j=1,nyo
     lats_orog(j)=-(180.0/nyo)*j+90.0+(90.0/nyo)
  enddo

  ! Calculate example orographic longitudes

  do i=1,nxo
     lons_orog(i)=(360.0/nxo)*i-(180.0/nxo)
  enddo

  ! Set the message level (6 is the default - all messages on)

  call glimmer_set_msg_level(6)

  ! Initialise the ice model

  call initialise_glint(ice_sheet, &
       climate%all_grid%lats, &
       climate%all_grid%lons, &
       (/paramfile/), &
       orog=orog_out, &
       ice_frac=ice_frac, &
       albedo=albedo, &
       orog_longs=lons_orog, &
       orog_lats=lats_orog, &
       daysinyear=365)

  ! Get coverage maps for the ice model instances

  if (glint_coverage_map(ice_sheet,coverage,cov_orog).ne.0) then
     call write_log('Unable to get coverage maps',GM_FATAL,__FILE__,__LINE__)
     stop
  endif

  ! Do initial timesteps ---------------------------------------------------------------------------

  time=0

  do
     call example_climate(climate,precip,temp,real(time,rk))
     call glint(ice_sheet,time,temp,precip,zonwind,merwind,orog, &
          orog_out=orog_out,   albedo=albedo,         output_flag=out, &
          ice_frac=ice_frac,   water_out=fw,          water_in=fw_in, &
          total_water_in=twin, total_water_out=twout, ice_volume=ice_vol) 
     time=time+climate%climate_tstep
     if (time>climate%initial_years*climate%hours_in_year) exit
  end do

  ! Do main loop

  do
     do i=1,climate%hours_in_year,climate%climate_tstep
        call example_climate(climate,precip,temp,real(time,rk))
        call glint(ice_sheet,time,temp,precip,zonwind,merwind,orog, &
             orog_out=orog_out,   albedo=albedo,         output_flag=out, &
             ice_frac=ice_frac,   water_out=fw,          water_in=fw_in, &
             total_water_in=twin, total_water_out=twout, ice_volume=ice_vol)
        time=time+climate%climate_tstep
        if (time>climate%total_years*climate%hours_in_year) exit
     end do

     if (time>climate%total_years*climate%hours_in_year) exit
     time=time+climate%hours_in_year-climate%climate_tstep

     do i=1,climate%years_ratio-1
        call glint(ice_sheet,time,temp,precip,zonwind,merwind,orog, &
             orog_out=orog_out,   albedo=albedo,         output_flag=out, &
             ice_frac=ice_frac,   water_out=fw,          water_in=fw_in, &
             total_water_in=twin, total_water_out=twout, ice_volume=ice_vol, &
             skip_mbal=.true.)
        time=time+climate%hours_in_year
     end do

     time=time-climate%hours_in_year+climate%climate_tstep
     if (time>climate%total_years*climate%hours_in_year) exit
  enddo

  ! Finalise/tidy up everything ------------------------------------------------------------

  call end_glint(ice_sheet)

100 format(f9.5)
101 format(e12.5)

end program glint_example
