
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

program glimmer_example

!*FD This program demonstrates the use of GLIMMER. It loads in
!*FD some example global fields and associated grid data,
!*FD Initialises the model, and then runs it for 1000 years.

  use glint_main
  use glimmer_log
  implicit none

  ! Program variables -------------------------------------------------------------------

  type(glint_params) :: ice_sheet    ! This is the derived type variable that holds all 
                                       ! domains of the ice model

  character(fname_length) :: paramfile='g_land.config'  ! The top-level configuration file

  integer :: total_years ! Length of run in years

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

  real(rk),dimension(:),  allocatable :: lats      ! Latitudes of normal global gridpoints
  real(rk),dimension(:),  allocatable :: lons      ! Longitudes of normal global gridpoints
  real(rk),dimension(:),  allocatable :: lats_orog ! Latitudes of global orography gridpoints
  real(rk),dimension(:),  allocatable :: lons_orog ! Longitudes of global oropraphy gridpoints

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

  ! Set dimensions of global grids

  nx=64   ; ny=32           ! Normal global grid
  nxo=200 ; nyo=100         ! Grid used for orographic output

  ! Allocate arrays appropriately

  ! start logging
  call open_log(unit=101)  

  allocate(temp(nx,ny),precip(nx,ny),zonwind(nx,ny))
  allocate(merwind(nx,ny),lats(ny),lons(nx),orog(nx,ny))
  allocate(coverage(nx,ny),orog_out(nxo,nyo),albedo(nx,ny),ice_frac(nx,ny),fw(nx,ny))
  allocate(lats_orog(nyo),lons_orog(nxo),cov_orog(nxo,nyo),fw_in(nx,ny))

  ! Set array contents to zero

  temp=0.0
  precip=0.0
  zonwind=0.0
  merwind=0.0
  lats=0.0
  lons=0.0
  orog=0.0
  albedo=0.0
  orog_out=0.0

  ! Set up global grids ----------------------------------------------------------------

  ! Calculate example orographic latitudes

  do j=1,nyo
    lats_orog(j)=-(180.0/nyo)*j+90.0+(90.0/nyo)
  enddo

  ! Calculate example orographic longitudes

  do i=1,nxo
    lons_orog(i)=(360.0/nxo)*i-(180.0/nxo)
  enddo

  ! Load in latitudes for normal global grid

  open(20,file='latitudes.txt')
  do j=1,ny
    read(20,100)lats(j)
  enddo
  close(20)

  ! Load in longitudes for normal global grid

  open(20,file='longitudes.txt')
  do i=1,nx
    read(20,100)lons(i)
  enddo
  close(20)
 
  ! Load in example global fields ------------------------------------------------------

  ! Load in sample temperature field

  open(20,file='tempfield.dat')
  do i=1,nx
    do j=1,ny
      read(20,101)temp(i,j)
    enddo
  enddo
  close(20)

  ! Load in sample precip field

  open(20,file='precipfield.dat')
  do i=1,nx
    do j=1,ny
      read(20,101)precip(i,j)
    enddo
  enddo
  close(20)

  ! Load in sample zonal wind field

  open(20,file='zonwindfield.dat')
  do i=1,nx
    do j=1,ny
      read(20,101)zonwind(i,j)
    enddo
  enddo
  close(20)

  ! Load in sample meridional wind field

  open(20,file='merwindfield.dat')
  do i=1,nx
    do j=1,ny
      read(20,101)merwind(i,j)
    enddo
  enddo
  close(20)

  ! Load global (large-scale) orography

  open(20,file='largescale.orog')
  do i=1,nx
    do j=1,ny
      read(20,101)orog(i,j)
    enddo
  enddo
  close(20)

  ! Do various bits of initialisation -------------------------------------------------------

  precip=precip/3600.0   ! Convert background precip rate from mm/h to mm/s

  ! Set the message level (6 is the default - all messages on)

  call glimmer_set_msg_level(6)

  ! Initialise the ice model

  call initialise_glint(ice_sheet,lats,lons,paramfile,orog=orog,ice_frac=ice_frac, &
       albedo=albedo,orog_longs=lons_orog,orog_lats=lats_orog)

  ! Get coverage maps for the ice model instances

  if (glint_coverage_map(ice_sheet,coverage,cov_orog).ne.0) then
    call write_log('Unable to get coverage maps',GM_FATAL,__FILE__,__LINE__)
    stop
  endif

  ! Get run length -------------------------------------------------------------------------

  Print*,'* Enter length of run in years:'
  Read*,total_years

  ! Do timesteps ---------------------------------------------------------------------------

  do time=0*24*360,total_years*24*360,24*360
    call glint(ice_sheet,time,temp,precip,zonwind,merwind,orog, &
                 orog_out=orog_out,   albedo=albedo,         output_flag=out, &
                 ice_frac=ice_frac,   water_out=fw,          water_in=fw_in, &
                 total_water_in=twin, total_water_out=twout, ice_volume=ice_vol) 
    call write_log_div ! Print a row of stars
  enddo

  ! Finalise/tidy up everything ------------------------------------------------------------

  call end_glint(ice_sheet)

100 format(f9.5)
101 format(e12.5)

end program glimmer_example
