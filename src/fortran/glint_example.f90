
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
  implicit none

  ! Program variables -------------------------------------------------------------------

  type(glint_params) :: ice_sheet    ! This is the derived type variable that holds all 
                                       ! domains of the ice model

  character(fname_length) :: paramfile  ! The top-level configuration file

  integer :: total_years ! Length of run in years

  ! Pointer arrays to the global climate data -------------------------------------------

  real(rk),dimension(:,:,:),pointer :: orog_clim     => null()  ! Orography
  real(rk),dimension(:,:,:),pointer :: precip_clim   => null()  ! Precip
  real(rk),dimension(:,:,:),pointer :: precip_clim2  => null()  ! Precip
  real(rk),dimension(:,:,:),pointer :: surftemp_clim => null()  ! Surface temperature

  type(global_grid),pointer :: orog_grid   => null()
  type(global_grid),pointer :: precip_grid => null()
  type(global_grid),pointer :: temp_grid   => null()

  ! arrays for interpolation ------------------------------------------------------------

  logical,dimension(:,:),allocatable :: maskin,maskout
  real,dimension(:,:),allocatable :: intin,intout

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

  real(rk),dimension(:),  pointer :: lats => null()     ! Latitudes of normal global gridpoints
  real(rk),dimension(:),  pointer :: lons => null()     ! Longitudes of normal global gridpoints
  real(rk),dimension(:),  allocatable :: lats_orog ! Latitudes of global orography gridpoints
  real(rk),dimension(:),  allocatable :: lons_orog ! Longitudes of global oropraphy gridpoints

  ! Scalars which hold information about the global grid --------------------------------

  integer :: nx,ny   ! Size of normal global grid
  integer :: nxo,nyo ! Size of global orography grid
  integer :: nxp,nyp ! Precip grid

  ! Scalar model outputs ----------------------------------------------------------------

  real(rk) :: twin     ! Timestep-integrated input water flux (kg)
  real(rk) :: twout    ! Timestep-integrated output water flux (kg)
  real(rk) :: ice_vol  ! Total ice volume (m^3)
  
  ! Other variables ---------------------------------------------------------------------

  logical :: out    ! Outputs set flag
  integer :: i,j    ! Array index counters
  integer :: time   ! Current time (hours)
  integer :: ierr

  ! -------------------------------------------------------------------------------------
  ! Executable code starts here - Basic initialisation
  ! -------------------------------------------------------------------------------------

  Print*,'Enter name of configuration file:'
  read*,paramfile

  ! Read in climate data

  call read_ncdf_3d('monthly_precip_mean_1974-2003.nc','prate',precip_clim2,precip_grid)
  call read_ncdf_3d('surf_temp_6h_1974-2003.nc','air_temperature',surftemp_clim,temp_grid)
  call read_ncdf_3d('global_orog.nc','hgt',orog_clim,orog_grid)

  ! Fix up a few things

  surftemp_clim=surftemp_clim-273.15       ! Convert temps to degreesC

  ! Set dimensions of global grids

  nx=size(surftemp_clim,1) ; ny=size(surftemp_clim,2) ! Normal global grid
  nxp=size(precip_clim2,1) ; nyp=size(precip_clim2,2) ! Precip input grid
  nxo=200 ; nyo=100                               ! Grid used for orographic output

  ! Do precip interpolation

  allocate(precip_clim(nx,ny,12))
  allocate(intin(nxp,nyp),maskin(nxp,nyp))
  allocate(intout(nx,ny),maskout(nx,ny))

  do i=1,12
     intin=precip_clim2(:,:,i)
     maskin=.true.
     maskout=.true.
     call global_interp(precip_grid,intin,temp_grid,intout,in_mask=maskin,out_mask=maskout,error=ierr)
     precip_clim(:,:,i)=intout
  end do

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
  orog=real(orog_clim(:,:,1))                    ! Put orography where it belongs

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

  print*,size(lats),size(lons),size(orog_out),size(ice_frac),size(albedo),size(lons_orog),size(lats_orog)

  call initialise_glint(ice_sheet,temp_grid%lats,temp_grid%lons,paramfile,orog=orog_out,ice_frac=ice_frac, &
       albedo=albedo,orog_longs=lons_orog,orog_lats=lats_orog,daysinyear=365)

  ! Get coverage maps for the ice model instances

  if (glint_coverage_map(ice_sheet,coverage,cov_orog).ne.0) then
    call write_log('Unable to get coverage maps',GM_FATAL,__FILE__,__LINE__)
    stop
  endif

  ! Get run length -------------------------------------------------------------------------

  Print*,'* Enter length of run in years:'
  Read*,total_years

  ! Do timesteps ---------------------------------------------------------------------------

  do time=0,total_years*24*365+6,6
     call example_climate(precip_clim,surftemp_clim,precip,temp,real(time,rk))
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

contains

  subroutine read_ncdf_3d(filename,varname,array,grid)

    use netcdf

    character(*)                      :: filename,varname
    real(rk),dimension(:,:,:),pointer :: array
    type(global_grid),        pointer :: grid

    real(rk),dimension(:),allocatable :: dim1,dim2,dim3
    real(rk),dimension(:),allocatable :: lonbound,latbound

    integer  :: ncerr     ! NetCDF error 
    integer  :: ncid      ! NetCDF file id
    integer  :: varid     ! NetCDF variable id
    integer  :: ndims     ! Number of dimensions
    integer  :: i,args
    real(rk) :: offset=0.0,scale=1.0
    integer,      dimension(3) :: dimids,dimlens
    character(20),dimension(3) :: dimnames
    real(dp),     dimension(:,:),allocatable :: lnb,ltb
    logical :: lonb_present,latb_present

    if (associated(array)) deallocate(array)
    if (associated(grid))  deallocate(grid)

    ! Open file ----------------------------------------

    ncerr=nf90_open(filename,0,ncid)
    call handle_err(ncerr,__LINE__)

    ! Find out the id of variable and its dimensions ---

    ncerr=nf90_inq_varid(ncid,varname,varid)
    call handle_err(ncerr,__LINE__)
    ncerr=nf90_inquire_variable(ncid, varid, ndims=ndims)
    call handle_err(ncerr,__LINE__)

    ! If not a 3d variable, flag and error and exit ----

    if (ndims/=3) then
       print*,'NetCDF: Requested variable only has ',ndims,' dimensions'
       stop
    end if

    ! Get dimensions ids -------------------------------

    ncerr=nf90_inquire_variable(ncid, varid, dimids=dimids)
    call handle_err(ncerr,__LINE__)

    ! Retrieve dimension names -------------------------

    do i=1,3
       ncerr=nf90_inquire_dimension(ncid, dimids(i), &
            name=dimnames(i),len=dimlens(i))
       call handle_err(ncerr,__LINE__)
    end do

    ! Allocate output and dimension arrays -------------

    allocate(array(dimlens(1),dimlens(2),dimlens(3)))
    allocate(dim1(dimlens(1)))
    allocate(dim2(dimlens(2)))
    allocate(dim3(dimlens(3)))
    
    ! Retrieve variable contents -----------------------

    ncerr=nf90_get_var(ncid, varid, array)
    call handle_err(ncerr,__LINE__)

    ! Get scaling and offset, if present, and apply ----

    ncerr=nf90_get_att(ncid, varid, 'add_offset', offset)
    if (ncerr/=NF90_NOERR) then
       offset=0.0
       ncerr=NF90_NOERR
    end if

    ncerr=nf90_get_att(ncid, varid, 'scale_factor', scale)
    if (ncerr/=NF90_NOERR) then
       scale=1.0
       ncerr=NF90_NOERR
    end if

    array=offset+(array*scale)

    ! Get dimension variables --------------------------

    ncerr=nf90_inq_varid(ncid,dimnames(1),varid)
    call handle_err(ncerr,__LINE__)
    ncerr=nf90_get_var(ncid, varid, dim1)
    call handle_err(ncerr,__LINE__)

    ncerr=nf90_inq_varid(ncid,dimnames(2),varid)
    call handle_err(ncerr,__LINE__)
    ncerr=nf90_get_var(ncid, varid, dim2)
    call handle_err(ncerr,__LINE__)

    ncerr=nf90_inq_varid(ncid,dimnames(3),varid)
    call handle_err(ncerr,__LINE__)
    ncerr=nf90_get_var(ncid, varid, dim3)
    call handle_err(ncerr,__LINE__)

    ! Get boundary arrays, if present ------------------

    ncerr=nf90_inq_varid(ncid,'bounds_lon',varid)
    if (ncerr/=NF90_NOERR) then
       lonb_present=.false.
       ncerr=NF90_NOERR
    else
       ncerr=nf90_inquire_variable(ncid, varid, dimids=dimids)
       call handle_err(ncerr,__LINE__)
       do i=1,2
          ncerr=nf90_inquire_dimension(ncid, dimids(i), name=dimnames(i),len=dimlens(i))
          call handle_err(ncerr,__LINE__)
       end do
       allocate(lnb(dimlens(1),dimlens(2)),lonbound(dimlens(2)+1))
       ncerr=nf90_get_var(ncid, varid,lnb)
       call handle_err(ncerr,__LINE__)
       do i=1,dimlens(2)
          lonbound(i)=lnb(1,i)
          lonbound(i+1)=lnb(2,i)
       end do
       deallocate(lnb)
       lonb_present=.true.
    end if

    ncerr=nf90_inq_varid(ncid,'bounds_lat',varid)
    if (ncerr/=NF90_NOERR) then
       latb_present=.false.
       ncerr=NF90_NOERR
    else
       ncerr=nf90_inquire_variable(ncid, varid, dimids=dimids)
       call handle_err(ncerr,__LINE__)
       do i=1,2
          ncerr=nf90_inquire_dimension(ncid, dimids(i), name=dimnames(i),len=dimlens(i))
          call handle_err(ncerr,__LINE__)
       end do
       allocate(ltb(dimlens(1),dimlens(2)),latbound(dimlens(2)+1))
       ncerr=nf90_get_var(ncid, varid,ltb)
       call handle_err(ncerr,__LINE__)
       do i=1,dimlens(2)
          latbound(i)=ltb(1,i)
          latbound(i+1)=ltb(2,i)
       end do
       deallocate(ltb)
       latb_present=.true.
    end if

    ! Construct grid type ---------------------

    allocate(grid)
    args=0
    if (lonb_present) args=args+1
    if (latb_present) args=args+2

    select case(args)
    case(0)
       call new_global_grid(grid,dim1,dim2)
    case(1)
       call new_global_grid(grid,dim1,dim2,lonb=lonbound)
    case(2)
       call new_global_grid(grid,dim1,dim2,latb=latbound)
    case(3)
       call new_global_grid(grid,dim1,dim2,lonb=lonbound,latb=latbound)
    end select

    ! Tidy up ---------------------------------

    deallocate(dim1,dim2,dim3)
    if (allocated(latbound)) deallocate(latbound)
    if (allocated(lonbound)) deallocate(lonbound)

  end subroutine read_ncdf_3d

  subroutine handle_err(status,line)

    use netcdf

    integer, intent (in) :: status
    integer, intent (in) :: line
    
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       print *, 'Line:',line
       stop "Stopped"
    end if
  end subroutine handle_err

  subroutine example_climate(precip_clim,temp_clim,precip,temp,time)

    real(rk),dimension(:,:,:),intent(in) :: precip_clim,temp_clim
    real(rk),dimension(:,:),intent(out)  :: precip,temp
    real(rk),intent(in) :: time

    integer :: ntemp,nprecip
    real(rk) :: tsp,tst
    real(rk) :: hoursinyear=365.0*24.0
    real(rk) :: pos
    integer :: lower,upper

    ntemp=size(temp_clim,3) ; nprecip=size(precip_clim,3)
    tst=hoursinyear/ntemp ; tsp=hoursinyear/nprecip
    
    ! Temperature first
    
    lower=int(time/tst)
    upper=lower+1
    pos=mod(time,tst)/tst
    call fixbounds(lower,1,ntemp)
    call fixbounds(upper,1,ntemp)
    temp=linear_interp(temp_clim(:,:,lower),temp_clim(:,:,upper),pos)

    ! precip

    lower=int(time/tsp)
    upper=lower+1
    pos=mod(time,tsp)/tsp
    call fixbounds(lower,1,nprecip)
    call fixbounds(upper,1,nprecip)
    precip=linear_interp(precip_clim(:,:,lower),precip_clim(:,:,upper),pos)

  end subroutine example_climate

  function linear_interp(a,b,pos)

    real(rk),dimension(:,:),intent(in) :: a,b
    real(rk),dimension(size(a,1),size(a,2)) :: linear_interp
    real(rk),               intent(in) :: pos

    linear_interp=a*(1.0-pos)+b*pos

  end function linear_interp

  subroutine fixbounds(in,bottom,top)

    integer :: in,top,bottom

    do
       if (in<=top) exit
       in=in-(top-bottom+1)
    end do

    do
       if (in>=bottom) exit
       in=in+(top-bottom+1)
    end do

  end subroutine fixbounds

end program glint_example
