
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer.f90 - part of the GLIMMER ice model              + 
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

module glimmer_main

  !*FD  This is the main glimmer module, which contains the top-level 
  !*FD  subroutines and derived types comprising the glimmer ice model.

  use glimmer_global
  use glimmer_object

! Derived type definitions

! the glimmer_params type contains parameters relevant to all instances of
! the model - i.e. those parameters which pertain to the global model 

! quantities which vary between instances are contained in the glimmer_instance object

  type glimmer_params !*FD derived type containing parameters 
                      !*FD relevant to all instances of the model - i.e. those 
                      !*FD parameters which pertain to the global model. 

    integer :: nxg,nyg                                     !*FD Dimensions of the global input fields
                                                           !*FD (grid-points)
    real(rk),pointer,dimension(:) :: glat,glon             !*FD Locations of data-points in global fields (degrees)
    real(rk),pointer,dimension(:) :: glat_bound,glon_bound !*FD Boundaries of grid boxes latitudinally
                                                           !*FD and longitudinally (degrees)
    integer :: ninstances                                  !*FD Number of ice model instances
    character(fname_length) :: paramfile                   !*FD Name of global parameter file
    type(glimmer_instance),pointer,dimension(:) :: instances  !*FD Array of glimmer\_instances
    real(rk) :: radea                                      !*FD Radius of the earth (m)

    real(rk) :: tstep_main                                 !*FD Main timestep (years)

    real(rk) :: av_start_time                              !*FD Holds the value of time from the last occasion
                                                           !*FD averaging was restarted.
    integer :: av_steps                                    !*FD Holds the number of times glimmer has been
                                                           !*FD called in current round of averaging

    real(rk),pointer,dimension(:,:) :: g_av_precip         !*FD globally averaged precip
    real(rk),pointer,dimension(:,:) :: g_av_temp           !*FD globally averaged temperature 
    real(rk),pointer,dimension(:,:) :: g_max_temp          !*FD global maximum temperature
    real(rk),pointer,dimension(:,:) :: g_min_temp          !*FD global minimum temperature
    real(rk),pointer,dimension(:,:) :: g_temp_range        !*FD global temperature range
    real(rk),pointer,dimension(:,:) :: g_av_zonwind        !*FD globally averaged zonal wind 
    real(rk),pointer,dimension(:,:) :: g_av_merwind        !*FD globally averaged meridional wind 

    real(rk),pointer,dimension(:,:) :: total_coverage      !*FD fractional coverage by all ice model instances
    logical :: coverage_calculated                         !*FD have we calculated this yet?
    real(rk),pointer,dimension(:,:) :: cov_normalise       !*FD normalisation values for coverage calculation

    integer :: file_unit                                   !*FD File unit to use for all ouput operations
    integer :: logf_unit                                   !*FD File unit for logging. Could be confusing
                                                           !*FD as the same unit is used by all model instances...

  end type glimmer_params

! global parameters/constants

  integer,parameter  :: days_in_year=360                   !*FD The number of days in a year  
  real(rk),parameter :: pi=3.141592654                     !*FD The value of pi
  real(rk),parameter :: years2hours=24.0*days_in_year      !*FD Years to hours conversion factor
  real(rk),parameter :: hours2years=1/years2hours          !*FD Hours to years conversion factor

! Private names

  private glimmer_allocate_arrays,glimmer_init_common
  private get_globals,get_fnamelist,calc_bounds,pi

contains

  subroutine initialise_glimmer(params,lats,longs,paramfile)

    !*FD Initialises the model

    use glimmer_project

    implicit none

    ! Subroutine argument declarations

    type(glimmer_params),   intent(inout) :: params      !*FD parameters to be set
    real(rk),dimension(:),  intent(in)    :: lats,longs  !*FD location of gridpoints in global data
    character(fname_length),intent(in)    :: paramfile   !*FD name of file containing parameters for all 
                                                         !*FD required instances. Eventually, this will be
                                                         !*FD defined in XML, but for the moment, it's a namelist

    ! Internal variables

    character(fname_length),dimension(:),allocatable :: fnamelist   ! The parameter filenames for each instance
    integer :: i,j
    real(rk) :: dlon

    ! ---------------------------------------------------------------
    ! Basic initialisation beginning with common initialisation
    ! ---------------------------------------------------------------

    call glimmer_init_common(params%logf_unit,params)

    ! Set up constants and arrays in the parameter structure

    params%paramfile = paramfile
    params%nxg       = size(longs)
    params%nyg       = size(lats)

    ! Allocate arrays

    call glimmer_allocate_arrays(params)

    ! Initialise arrays

    params%g_av_precip  = 0.0
    params%g_av_temp    = 0.0
    params%g_max_temp   = -1000.0
    params%g_min_temp   = 1000.0
    params%g_temp_range = 0.0
    params%g_av_zonwind = 0.0
    params%g_av_merwind = 0.0

    ! Copy lats and lons over

    params%glat=lats
    params%glon=longs

    ! ---------------------------------------------------------------
    ! Calculate box boundaries
    ! ---------------------------------------------------------------

    call calc_bounds(params%glon,params%glat,params%glon_bound,params%glat_bound)
   
    ! ---------------------------------------------------------------
    ! Open the namelist file and read in the parameters for each run,
    ! allocating appropriately
    ! ---------------------------------------------------------------

    open (params%file_unit,file=paramfile)

    ! Read in global parameters from this file

    call get_globals(params)

    ! Allocate array of glimmer instances
    ! and the list of associated namelist files

    allocate(fnamelist(params%ninstances),params%instances(params%ninstances))

    ! Get the list of namelist files

    call get_fnamelist(params%file_unit,fnamelist)

    ! Close the glimmer parameter file

    close(params%file_unit)

    ! ---------------------------------------------------------------
    ! Read namelist file for each glimmer instance, and
    ! initialise accordingly.
    ! Also add up coverage map - this isn't too sophisticated 
    ! at the moment.
    ! ---------------------------------------------------------------

    params%total_coverage=0.0

    do i=1,params%ninstances
      call glimmer_i_initialise(params%file_unit,    &
                             fnamelist(i),        &
                             params%instances(i), &
                             params%radea,        &
                             params%glon,         &
                             params%glat,         &
                             params%glon_bound,   &
                             params%glat_bound,   &
                             params%tstep_main)
      params%total_coverage = params%total_coverage &
                            + params%instances(i)%frac_coverage
    enddo

    ! Normalisation array set to the sum of all coverage

    params%cov_normalise=params%total_coverage

    ! Check we don't have coverage greater than one at any point.

    where (params%total_coverage>1.0) params%total_coverage=1.0
    params%coverage_calculated=.true.

  end subroutine initialise_glimmer

!================================================================================

  subroutine glimmer(params,time,temp,precip,zonwind,merwind,orog, &
                     output_flag,orog_out,albedo,ice_frac,fw_flux)

    !*FD Main Glimmer subroutine.
    !*FD
    !*FD This should be called daily or hourly, depending on
    !*FD the mass-balance scheme being used. It does all necessary 
    !*FD spatial and temporal averaging, and calls the dynamics 
    !*FD part of the model when required. 
    !*FD
    !*FD Input fields should 
    !*FD be taken as means over the period since the last call, or, 
    !*FD in the case of precip, the total accumulation since the last call. See
    !*FD the user documentation for more information.

    use glimmer_utils
    use glimmer_mbal
    use glimmer_interp

    implicit none

    ! Subroutine argument declarations

    type(glimmer_params),            intent(inout) :: params          !*FD parameters for this run
    real(rk),                        intent(in)    :: time            !*FD Current model time        (hours)
    real(rk),dimension(:,:),         intent(in)    :: temp            !*FD Surface temperature field (celcius)
    real(rk),dimension(:,:),         intent(in)    :: precip          !*FD Precipitation fields  (mm accumulated)
    real(rk),dimension(:,:),         intent(in)    :: zonwind,merwind !*FD Zonal and meridional components 
                                                                      !*FD of the wind field (m/s)
    real(rk),dimension(:,:),         intent(inout) :: orog            !*FD The large-scale orography (m)
    logical,                optional,intent(out)   :: output_flag     !*FD Set if outputs set
    real(rk),dimension(:,:),optional,intent(inout) :: orog_out        !*FD The fed-back, output orography
    real(rk),dimension(:,:),optional,intent(inout) :: albedo          !*FD surface albedo
    real(rk),dimension(:,:),optional,intent(inout) :: ice_frac        !*FD grid-box ice-fraction
    real(rk),dimension(:,:),optional,intent(inout) :: fw_flux         !*FD Freshwater flux (mm)
    
    ! Internal variables

    integer :: i
    real(rk),dimension(size(orog,1),size(orog,2)) :: orog_out_temp,albedo_temp,if_temp,fw_temp
    type(output_flags) :: out_f

    ! Reset output flag

    if (present(output_flag)) output_flag=.false.

    ! ---------------------------------------------------------
    ! Do averaging and so on...
    ! Averages
    ! ---------------------------------------------------------

    params%g_av_temp    = params%g_av_temp    + temp
    params%g_av_precip  = params%g_av_precip  + precip
    params%g_av_zonwind = params%g_av_zonwind + zonwind
    params%g_av_merwind = params%g_av_merwind + merwind

    ! Ranges of temperature

    where (temp > params%g_max_temp) params%g_max_temp=temp
    where (temp < params%g_min_temp) params%g_min_temp=temp

    ! Increment step counter

    params%av_steps=params%av_steps+1

    ! ---------------------------------------------------------
    ! If this is a timestep, prepare global fields, and do a timestep
    ! for each model instance
    ! ---------------------------------------------------------

    if (time-params%av_start_time.ge.params%tstep_main*years2hours) then

      ! Set output_flag

      if (present(output_flag)) output_flag=.true.

      ! Reset output fields

      if (present(orog_out)) then
        orog_out  = 0.0
        out_f%orog=.true.
      else
        out_f%orog=.false.
      endif

      if (present(albedo)) then
        albedo    = 0.0
        out_f%albedo=.true.
      else
        out_f%albedo=.false.
      endif

      if (present(ice_frac)) then
        ice_frac  = 0.0
        out_f%ice_frac=.true.
      else
        out_f%ice_frac=.false.
      endif

      if (present(fw_flux)) then
        fw_flux = 0.0
        out_f%fresh_water=.true.
      else
        out_f%fresh_water=.false.
      endif

      ! Calculate averages by dividing by elapsed time,
      ! except for precip, where we want an accumulated total (mm)

      params%g_av_temp    = params%g_av_temp   /params%av_steps
      params%g_av_zonwind = params%g_av_zonwind/params%av_steps
      params%g_av_merwind = params%g_av_merwind/params%av_steps

      ! Calculate temperature range

      params%g_temp_range=params%g_max_temp-params%g_min_temp

      ! Do a timestep for each instance

      do i=1,params%ninstances
        call glimmer_i_tstep(params%file_unit,          &
                          params%logf_unit,             &
                          time*hours2years,             &
                          params%instances(i),          &
                          params%glat,                  &
                          params%glon,                  &
                          params%g_av_temp,             &
                          params%g_temp_range,          &
                          params%g_av_precip,           &
                          params%g_av_zonwind,          &
                          params%g_av_merwind,          &
                          orog,                         &
                          orog_out_temp,                &
                          albedo_temp,                  &
                          if_temp,                      &
                          fw_temp,                      &
                          out_f)

        ! Add this contribution to the output orography
  
        if (present(orog_out)) then
           where (params%cov_normalise==0.0)
              orog_out=0.0
           elsewhere
              orog_out=orog_out+(orog_out_temp* &
                   params%instances(i)%frac_coverage/ &
                   params%cov_normalise)
           end where
        endif

        if (present(albedo)) then
           where (params%cov_normalise==0.0) 
              albedo=0.0
           elsewhere
              albedo=albedo+(albedo_temp* &
                   params%instances(i)%frac_coverage/ &
                   params%cov_normalise)
           end where
        endif

        if (present(ice_frac)) then
           where (params%cov_normalise==0.0) 
              ice_frac=0.0
           elsewhere
              ice_frac=ice_frac+(if_temp* &
                   params%instances(i)%frac_coverage/ &
                   params%cov_normalise)
           end where
        endif

        if (present(fw_flux)) then
           where (params%cov_normalise==0.0) 
              fw_flux=0.0
           elsewhere
              fw_flux=fw_flux+(fw_temp* &
                   params%instances(i)%frac_coverage/ &
                   params%cov_normalise)
           end where
        endif
      enddo

      ! ---------------------------------------------------------
      ! Reset averaging fields, flags and counters
      ! ---------------------------------------------------------

      params%g_av_temp    = 0.0
      params%g_av_precip  = 0.0
      params%g_av_zonwind = 0.0
      params%g_av_merwind = 0.0
      params%g_temp_range = 0.0
      params%g_max_temp   = -1000.0
      params%g_min_temp   = 1000.0

      params%av_steps     = 0
      params%av_start_time = time

    endif

  end subroutine glimmer

!===================================================================

  subroutine end_glimmer(params)

    !*FD perform tidying-up operations for glimmer

    type(glimmer_params),intent(inout) :: params          !*FD parameters for this run
       
    ! end individual instances

    do i=1,params%ninstances
      call glimmer_i_end(params%instances(i)%model,params%file_unit)
    enddo

    ! close log file

    close(params%logf_unit)

  end subroutine end_glimmer

!=====================================================

  integer function glimmer_coverage_map(params,coverage)

    !*FD Retrieve ice model fractional 
    !*FD coverage map. This function is provided so that glimmer may
    !*FD be restructured without altering the interface.
    !*FDRV Three return values are possible:
    !*FDRV \begin{description}
    !*FDRV \item[0 ---] Successful return
    !*FDRV \item[1 ---] Coverage map not calculated yet (fail)
    !*FDRV \item[2 ---] Coverage array is the wrong size (fail)
    !*FDRV \end{description}

    type(glimmer_params),intent(in) :: params       !*FD ice model parameters
    real(rk),dimension(:,:),intent(out) :: coverage !*FD array to hold coverage map
    
    if (.not.params%coverage_calculated) then
      glimmer_coverage_map=1
      return
    endif

    if (size(coverage,1).ne.params%nxg.or. &
        size(coverage,2).ne.params%nyg) then
      glimmer_coverage_map=2
      return
    endif

    glimmer_coverage_map=0
    coverage=params%total_coverage

  end function glimmer_coverage_map

!=======================================================

  integer function glimmer_main_funit(params)

!*FD Returns the value of the 
!*FD main logical file unit used by glimmer
!*FDRV The value of the 
!*FDRV main logical file unit used by glimmer.

    type(glimmer_params),intent(in) :: params !*FD Ice model parameters
   
    glimmer_main_funit=params%file_unit

  end function glimmer_main_funit

!=====================================================

  subroutine glimmer_set_main_funit(params,unit)

  !*FD sets the main logical 
  !*FD file unit used by glimmer.

    type(glimmer_params),intent(inout) :: params !*FD ice model parameters
    integer,intent(in) :: unit !*FD number of logical file unit to be used

    params%file_unit=unit

  end subroutine glimmer_set_main_funit

!=====================================================

  subroutine glimmer_write_restart(params,unit,filename)

    !*FD write the restart
    !*FD file for the whole ice model. Note that the logical file unit 
    !*FD is not open at entry or exit to this subroutine.

    use glimmer_restart

    type(glimmer_params),intent(in) :: params   !*FD ice model parameters
    integer,             intent(in) :: unit     !*FD logical file unit to use
    character(*),        intent(in) :: filename !*FD filename to write

    integer :: i

    open(unit,file=filename,form='unformatted')
    
    write(unit) params% nxg
    write(unit) params% nyg
    write(unit) params% glat
    write(unit) params% glon
    write(unit) params% glat_bound
    write(unit) params% glon_bound
    write(unit) params% ninstances
    write(unit) params% paramfile
    write(unit) params% radea
    write(unit) params% tstep_main
    write(unit) params% av_start_time
    write(unit) params% av_steps
    write(unit) params% g_av_precip
    write(unit) params% g_av_temp
    write(unit) params% g_max_temp
    write(unit) params% g_min_temp
    write(unit) params% g_temp_range
    write(unit) params% g_av_zonwind
    write(unit) params% g_av_merwind
    write(unit) params% total_coverage 
    write(unit) params% coverage_calculated
    write(unit) params% cov_normalise
    write(unit) params% file_unit
    write(unit) params% logf_unit

    do i=1,params%ninstances
      call glimmer_i_write_restart(params%instances(i),unit)
    enddo

    close(unit)

  end subroutine glimmer_write_restart

!================================================

  subroutine glimmer_read_restart(params,unit,filename)

    !*FD Read the restart file
    !*FD for the ice model, and allocate necessary
    !*FD arrays.  Note that the file unit is closed on entry and exit.
    !*FD Also, this subroutine must be run `from cold' --- i.e.
    !*FD without any arrays having been allocated. The
    !*FD model cannot thus be restarted mid-run.

    use glimmer_restart

    implicit none

    type(glimmer_params),intent(out) :: params !*FD ice model parameters
    integer,intent(in) :: unit                 !*FD logical file unit to use
    character(*),intent(in) :: filename        !*FD filename to read

    integer :: i

    if (associated(params%glat))           deallocate(params%glat)
    if (associated(params%glon))           deallocate(params%glon)
    if (associated(params%glat_bound))     deallocate(params%glat_bound)
    if (associated(params%glon_bound))     deallocate(params%glon_bound)
    if (associated(params%instances))      deallocate(params%instances)
    if (associated(params%g_av_precip))    deallocate(params%g_av_precip)
    if (associated(params%g_av_temp))      deallocate(params%g_av_temp)
    if (associated(params%g_max_temp))     deallocate(params%g_max_temp)
    if (associated(params%g_min_temp))     deallocate(params%g_min_temp)
    if (associated(params%g_temp_range))   deallocate(params%g_temp_range)
    if (associated(params%g_av_zonwind))   deallocate(params%g_av_zonwind)
    if (associated(params%g_av_merwind))   deallocate(params%g_av_merwind)
    if (associated(params%total_coverage)) deallocate(params%total_coverage)
    if (associated(params%cov_normalise))  deallocate(params%cov_normalise)

    call glimmer_init_common(params%logf_unit,params)

    open(unit,file=filename,form='unformatted')

    read(unit) params% nxg
    read(unit) params% nyg

    call glimmer_allocate_arrays(params)

    read(unit) params% glat
    read(unit) params% glon
    read(unit) params% glat_bound
    read(unit) params% glon_bound
    read(unit) params% ninstances
    read(unit) params% paramfile
    read(unit) params% radea
    read(unit) params% tstep_main
    read(unit) params% av_start_time
    read(unit) params% av_steps
    read(unit) params% g_av_precip
    read(unit) params% g_av_temp
    read(unit) params% g_max_temp
    read(unit) params% g_min_temp
    read(unit) params% g_temp_range
    read(unit) params% g_av_zonwind
    read(unit) params% g_av_merwind
    read(unit) params% total_coverage 
    read(unit) params% coverage_calculated
    read(unit) params% cov_normalise
    read(unit) params% file_unit
    read(unit) params% logf_unit

    allocate(params%instances(params%ninstances))

    do i=1,params%ninstances
      call glimmer_i_read_restart(params%instances(i),unit,params%nxg,params%nyg)
    enddo  

    close(unit)

  end subroutine glimmer_read_restart

!----------------------------------------------------------------------
! PRIVATE INTERNAL GLIMMER SUBROUTINES FOLLOW.............
!----------------------------------------------------------------------

  subroutine glimmer_allocate_arrays(params)

  !*FD allocates glimmer arrays

    type(glimmer_params),intent(inout) :: params !*FD ice model parameters

    allocate(params%glat        (params%nyg))
    allocate(params%glon        (params%nxg))
    allocate(params%glat_bound  (params%nyg+1))
    allocate(params%glon_bound  (params%nxg+1))
    allocate(params%g_av_precip (params%nxg,params%nyg))
    allocate(params%g_av_temp   (params%nxg,params%nyg))
    allocate(params%g_max_temp  (params%nxg,params%nyg))
    allocate(params%g_min_temp  (params%nxg,params%nyg))
    allocate(params%g_temp_range(params%nxg,params%nyg))
    allocate(params%g_av_zonwind(params%nxg,params%nyg))
    allocate(params%g_av_merwind(params%nxg,params%nyg))
    allocate(params%total_coverage(params%nxg,params%nyg))
    allocate(params%cov_normalise(params%nxg,params%nyg))

  end subroutine glimmer_allocate_arrays

!====================================================

  subroutine glimmer_init_common(unit,params)

    !*FD Initialises common
    !*FD variables that need to be initialised, regardless
    !*FD of whether a restart is being performed. Also opens the log file

    use paramets, only: tau0,vel0,vis0,len0
    use physcon,  only: gn

    integer :: unit !*FD file unit to use for the log file
    type(glimmer_params) :: params !*FD ice model parameters

    ! THIS IS WHERE INITIALISATION THAT IS NEEDED REGARDLESS
    ! OF WHETHER WE ARE RESTARTING BELONGS!!!

    params%ninstances=1
    params%radea=6.37e6     
    params%tstep_main=1.0
    params%av_start_time=0.0
    params%av_steps=0
    params%coverage_calculated=.false.
    params%file_unit=20
    params%logf_unit=21

    ! Set up tau0 - this a moved from a the glimmer init subroutine

    tau0 = (vel0/(vis0*len0))**(1.0/gn) 

    ! Open log file - this is where all sorts of junk gets written, apparently.

    open(unit,file='glimmer.gll')

  end subroutine glimmer_init_common

!======================================================

  subroutine get_globals(params)

    !*FD reads namelist of global
    !*FD parameters from already-open file

    implicit none

    type(glimmer_params),intent(inout) :: params !*FD ice model parameters
    integer :: ninst
    real(rk) :: tinc

    namelist /file_paras/ ninst
    namelist /timesteps/ tinc

    read(params%file_unit,nml=timesteps)
    read(params%file_unit,nml=file_paras)

    params%ninstances=ninst
    params%tstep_main=tinc
 
   end subroutine get_globals

!=====================================================

  subroutine get_fnamelist(unit,fnamelist)

    !*FD reads namelist of parameter files 
    !*FD for each instance from already open file

    implicit none

    integer,intent(in) :: unit !*FD logical file unit to use
    character(*),dimension(:),intent(out) :: fnamelist !*FD array of namelist filenames
    integer :: i
    
    do i=1,size(fnamelist)
      read(unit,100)fnamelist(i)
    enddo

100 format(A)

  end subroutine get_fnamelist

!========================================================

  subroutine calc_bounds(lon,lat,lonb,latb)

  !*FD Calculates the boundaries between
  !*FD global grid-boxes. Note that we assume that the boundaries lie 
  !*FD half-way between the 
  !*FD points, both latitudinally and longitudinally, although 
  !*FD this isn't strictly true for a Gaussian grid.

    implicit none

    real(rk),dimension(:),intent(in) :: lon,lat !*FD locations of global grid-points (degrees)
    real(rk),dimension(:),intent(out) :: lonb,latb !*FD boundaries of grid-boxes (degrees)

    real(rk) :: dlon

    integer :: nxg,nyg,i,j

    nxg=size(lon) ; nyg=size(lat)

    ! Latitudes first - we assume the boundaries of the first and 
    ! last boxes coincide with the poles. Not sure how to
    ! handle it if they don't...

    latb(1)=90.0
    latb(nyg+1)=-90.0

    do j=2,nyg
      latb(j)=lat(j-1)-(lat(j-1)-lat(j))/2.0
    enddo

    ! Longitudes
    
    if (lon(1)<lon(nxg)) then
      dlon=lon(1)-lon(nxg)+360.0
    else
      dlon=lon(1)-lon(nxg)
    endif
    lonb(1)=lon(nxg)+dlon/2
    lonb(1)=loncorrect(lonb(1))      

    lonb(nxg+1)=lonb(1)
 
    do i=2,nxg
      if (lon(i)<lon(i-1)) then
        dlon=lon(i)-lon(i-1)+360.0
      else
        dlon=lon(i)-lon(i-1)
      endif
      lonb(i)=lon(i-1)+dlon/2
      lonb(i)=loncorrect(lonb(i))      
    enddo

  end subroutine calc_bounds

end module glimmer_main
