
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
  use glimmer_global_grid

  ! ------------------------------------------------------------
  ! GLIMMER_PARAMS derived type definition
  ! This is where default values are set.
  ! ------------------------------------------------------------

  type glimmer_params 
  
    !*FD Derived type containing parameters relevant to all instances of 
    !*FD the model - i.e. those parameters which pertain to the global model. 
  
    ! Global grids used ----------------------------------------

    type(global_grid) :: g_grid      !*FD The main global grid, used for 
                                     !*FD input and most outputs
    type(global_grid) :: g_grid_orog !*FD Global grid used for orography output.

    ! Ice model instances --------------------------------------

    integer                                     :: ninstances = 1       !*FD Number of ice model instances
    type(glimmer_instance),pointer,dimension(:) :: instances  => null() !*FD Array of glimmer\_instances

    ! Global model parameters ----------------------------------

    real(rk) :: radea      = 6.37e6 !*FD Radius of the earth (m)
    real(rk) :: tstep_main = 1.0    !*FD Main timestep (years)

    ! Averaging parameters -------------------------------------

    real(rk) :: av_start_time = 0.0 !*FD Holds the value of time from 
                                    !*FD the last occasion averaging was restarted (hours)
    integer  :: av_steps      = 0   !*FD Holds the number of times glimmer has 
                                    !*FD been called in current round of averaging.
    ! Averaging arrays -----------------------------------------

    real(rk),pointer,dimension(:,:) :: g_av_precip  => null() !*FD globally averaged precip
    real(rk),pointer,dimension(:,:) :: g_av_temp    => null() !*FD globally averaged temperature 
    real(rk),pointer,dimension(:,:) :: g_max_temp   => null() !*FD global maximum temperature
    real(rk),pointer,dimension(:,:) :: g_min_temp   => null() !*FD global minimum temperature
    real(rk),pointer,dimension(:,:) :: g_temp_range => null() !*FD global temperature range
    real(rk),pointer,dimension(:,:) :: g_av_zonwind => null() !*FD globally averaged zonal wind 
    real(rk),pointer,dimension(:,:) :: g_av_merwind => null() !*FD globally averaged meridional wind 

    ! Fractional coverage information --------------------------

    real(rk),pointer,dimension(:,:) :: total_coverage  => null()     !*FD Fractional coverage by 
                                                                     !*FD all ice model instances.
    real(rk),pointer,dimension(:,:) :: cov_normalise   => null()     !*FD Normalisation values 
                                                                     !*FD for coverage calculation.
    real(rk),pointer,dimension(:,:) :: total_cov_orog  => null()     !*FD Fractional coverage by 
                                                                     !*FD all ice model instances (orog).
    real(rk),pointer,dimension(:,:) :: cov_norm_orog   => null()     !*FD Normalisation values 
                                                                     !*FD for coverage calculation (orog).
    logical                         :: coverage_calculated = .false. !*FD Have we calculated the
                                                                     !*FD coverage map yet?
    ! File information -----------------------------------------

    character(fname_length) :: paramfile      !*FD Name of global parameter file
    integer                 :: file_unit = 20 !*FD File unit to use for all ouput operations
    integer                 :: logf_unit = 21 !*FD File unit for logging. Could be confusing
                                              !*FD as the same unit is used by all model instances...
  end type glimmer_params

  ! ------------------------------------------------------------
  ! global parameters/constants
  ! ------------------------------------------------------------

  integer, parameter :: days_in_year=360                   !*FD The number of days in a year  
  real(rk),parameter :: pi=3.141592654                     !*FD The value of pi
  real(rk),parameter :: years2hours=24.0*days_in_year      !*FD Years to hours conversion factor
  real(rk),parameter :: hours2years=1/years2hours          !*FD Hours to years conversion factor
  real(rk),parameter :: hours2seconds=3600.0               !*FD Hours to seconds conversion factor

  ! Private names -----------------------------------------------

  private glimmer_allocate_arrays,glimmer_init_common
  private get_globals,get_fnamelist,calc_bounds,pi

contains

  subroutine initialise_glimmer(params,lats,longs,paramfile,latb,lonb)

    !*FD Initialises the model

    use glimmer_project

    implicit none

    ! Subroutine argument declarations

    type(glimmer_params),          intent(inout) :: params      !*FD parameters to be set
    real(rk),dimension(:),         intent(in)    :: lats,longs  !*FD location of gridpoints 
                                                                !*FD in global data.
    character(fname_length),       intent(in)    :: paramfile   !*FD name of file containing 
                                                                !*FD parameters for all required 
                                                                !*FD instances. Eventually, this
                                                                !*FD will be defined in XML, but
                                                                !*FD for the moment, it's a namelist.
    real(rk),dimension(:),optional,intent(in)    :: latb        !*FD Locations of the latitudinal 
                                                                !*FD boundaries of the grid-boxes.
    real(rk),dimension(:),optional,intent(in)    :: lonb        !*FD Locations of the longitudinal
                                                                !*FD boundaries of the grid-boxes.

    ! Internal variables

    character(fname_length),dimension(:),allocatable :: fnamelist   ! The parameter filenames for each instance
    integer :: i,args

    ! ---------------------------------------------------------------
    ! Basic initialisation beginning with common initialisation
    ! ---------------------------------------------------------------

    call glimmer_init_common(params%logf_unit,params)

    ! Set up constants and arrays in the parameter structure --------

    params%paramfile = paramfile
 
    ! Initialise main global grid -----------------------------------

    args=0

    if (present(lonb)) args=args+1
    if (present(latb)) args=args+2

    select case(args)
    case(0)
      call new_global_grid(params%g_grid,longs,lats)
    case(1)
      call new_global_grid(params%g_grid,longs,lats,lonb=lonb)
    case(2)
      call new_global_grid(params%g_grid,longs,lats,latb=latb)
    case(3)
      call new_global_grid(params%g_grid,longs,lats,lonb=lonb,latb=latb)
    end select

    ! Initialise orography grid identically at the moment -----------

    call copy_global_grid(params%g_grid,params%g_grid_orog)

    ! Allocate arrays -----------------------------------------------

    call glimmer_allocate_arrays(params)

    ! Initialise arrays ---------------------------------------------

    params%g_av_precip  = 0.0
    params%g_av_temp    = 0.0
    params%g_max_temp   = -1000.0
    params%g_min_temp   = 1000.0
    params%g_temp_range = 0.0
    params%g_av_zonwind = 0.0
    params%g_av_merwind = 0.0

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
    params%total_cov_orog=0.0

    do i=1,params%ninstances
      call glimmer_i_initialise(params%file_unit,    &
                                fnamelist(i),        &
                                params%instances(i), &
                                params%radea,        &
                                params%g_grid,       &
                                params%tstep_main)

      params%total_coverage = params%total_coverage &
                            + params%instances(i)%frac_coverage
      params%total_cov_orog = params%total_cov_orog &
                            + params%instances(i)%frac_cov_orog
    enddo

    ! Normalisation array set to the sum of all coverage

    params%cov_normalise=params%total_coverage
    params%cov_norm_orog=params%total_cov_orog

    ! Check we don't have coverage greater than one at any point.

    where (params%total_coverage>1.0) params%total_coverage=1.0
    where (params%total_cov_orog>1.0) params%total_cov_orog=1.0
    params%coverage_calculated=.true.

  end subroutine initialise_glimmer

!================================================================================

  subroutine glimmer(params,time,temp,precip,zonwind,merwind,orog, &
                     output_flag,orog_out,albedo,ice_frac,water_in, &
                     water_out,elapsed_time,total_water_in,total_water_out, &
                     ice_volume)

    !*FD Main Glimmer subroutine.
    !*FD
    !*FD This should be called daily or hourly, depending on
    !*FD the mass-balance scheme being used. It does all necessary 
    !*FD spatial and temporal averaging, and calls the dynamics 
    !*FD part of the model when required. 
    !*FD
    !*FD Input fields should be taken as means over the period since the last call.
    !*FD See the user documentation for more information.

    use glimmer_utils
    use glimmer_mbal
    use glimmer_interp

    implicit none

    ! Subroutine argument declarations

    type(glimmer_params),            intent(inout) :: params          !*FD parameters for this run
    real(rk),                        intent(in)    :: time            !*FD Current model time        (hours)
    real(rk),dimension(:,:),         intent(in)    :: temp            !*FD Surface temperature field (celcius)
    real(rk),dimension(:,:),         intent(in)    :: precip          !*FD Precipitation rate        (mm/s)
    real(rk),dimension(:,:),         intent(in)    :: zonwind,merwind !*FD Zonal and meridional components 
                                                                      !*FD of the wind field         (m/s)
    real(rk),dimension(:,:),         intent(inout) :: orog            !*FD The large-scale orography (m)
    logical,                optional,intent(out)   :: output_flag     !*FD Set true if outputs set
    real(rk),dimension(:,:),optional,intent(inout) :: orog_out        !*FD The fed-back, output orography (m)
    real(rk),dimension(:,:),optional,intent(inout) :: albedo          !*FD surface albedo
    real(rk),dimension(:,:),optional,intent(inout) :: ice_frac        !*FD grid-box ice-fraction
    real(rk),dimension(:,:),optional,intent(inout) :: water_in        !*FD Input water flux          (mm/s)
    real(rk),dimension(:,:),optional,intent(inout) :: water_out       !*FD Output water flux         (mm/s)
    real(rk),               optional,intent(inout) :: elapsed_time    !*FD Time over which output is valid (hours)
    real(rk),               optional,intent(inout) :: total_water_in  !*FD Area-integrated water flux in (kg)
    real(rk),               optional,intent(inout) :: total_water_out !*FD Area-integrated water flux out (kg)
    real(rk),               optional,intent(inout) :: ice_volume      !*FD Total ice volume (m$^3$)
    
    ! Internal variables

    integer :: i
    real(rk),dimension(:,:),allocatable :: albedo_temp,if_temp,wout_temp,orog_out_temp,win_temp
    real(rk) :: twin_temp,twout_temp,icevol_temp
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
        allocate(orog_out_temp(size(orog_out,1),size(orog_out,2)))
        out_f%orog=.true.
      else
        out_f%orog=.false.
      endif

      if (present(albedo)) then
        albedo    = 0.0
        allocate(albedo_temp(size(orog,1),size(orog,2)))
        out_f%albedo=.true.
      else
        out_f%albedo=.false.
      endif

      if (present(ice_frac)) then
        ice_frac  = 0.0
        allocate(if_temp(size(orog,1),size(orog,2)))
        out_f%ice_frac=.true.
      else
        out_f%ice_frac=.false.
      endif

      if (present(water_out)) then
        water_out = 0.0
        allocate(wout_temp(size(orog,1),size(orog,2)))
        out_f%water_out=.true.
      else
        out_f%water_out=.false.
      endif

      if (present(water_in)) then
        water_in = 0.0
        allocate(win_temp(size(orog,1),size(orog,2)))
        out_f%water_in=.true.
      else
        out_f%water_in=.false.
      endif

      ! Reset output total variables and set flags

      if (present(total_water_in))  then
        total_water_in  = 0.0
        out_f%total_win = .true.
      else
        out_f%total_win = .false.
      endif

      if (present(total_water_out)) then
        total_water_out = 0.0
        out_f%total_wout = .true.
      else
        out_f%total_wout = .false.
      endif

      if (present(ice_volume)) then
        ice_volume = 0.0
        out_f%ice_vol = .true.
      else
        out_f%ice_vol = .false.
      endif

      ! Calculate averages by dividing by number of steps elapsed
      ! since last model timestep.

      params%g_av_temp    = params%g_av_temp   /params%av_steps
      params%g_av_zonwind = params%g_av_zonwind/params%av_steps
      params%g_av_merwind = params%g_av_merwind/params%av_steps
      params%g_av_precip  = params%g_av_precip /params%av_steps

      ! Calculate total accumulated precipitation - multiply
      ! by time since last model timestep

      params%g_av_precip = params%g_av_precip*(time-params%av_start_time)*hours2seconds

      ! Calculate temperature range

      params%g_temp_range=params%g_max_temp-params%g_min_temp

      ! Do a timestep for each instance

      do i=1,params%ninstances
        call glimmer_i_tstep(params%file_unit,          &
                          params%logf_unit,             &
                          time*hours2years,             &
                          params%instances(i),          &
                          params%g_grid%lats,           &
                          params%g_grid%lons,           &
                          params%g_av_temp,             &
                          params%g_temp_range,          &
                          params%g_av_precip,           &
                          params%g_av_zonwind,          &
                          params%g_av_merwind,          &
                          orog,                         &
                          orog_out_temp,                &
                          albedo_temp,                  &
                          if_temp,                      &
                          win_temp,                     &
                          wout_temp,                    &
                          twin_temp,                    &
                          twout_temp,                   &
                          icevol_temp,                  &
                          out_f)

        ! Add this contribution to the output orography
  
        if (present(orog_out)) then
           where (params%cov_norm_orog==0.0)
              orog_out=0.0
           elsewhere
              orog_out=orog_out+(orog_out_temp* &
                   params%instances(i)%frac_cov_orog/ &
                   params%cov_norm_orog)
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

        if (present(water_in)) then
           where (params%cov_normalise==0.0) 
              water_in=0.0
           elsewhere
              water_in=water_in+(win_temp* &
                   params%instances(i)%frac_coverage/ &
                   params%cov_normalise)
           end where
        endif

        if (present(water_out)) then
           where (params%cov_normalise==0.0) 
              water_out=0.0
           elsewhere
              water_out=water_out+(wout_temp* &
                   params%instances(i)%frac_coverage/ &
                   params%cov_normalise)
           end where
        endif

        ! Add total water variables to running totals

        if (present(total_water_in))  total_water_in =total_water_in +twin_temp
        if (present(total_water_out)) total_water_out=total_water_out+twout_temp
        if (present(ice_volume))      ice_volume=ice_volume+icevol_temp

      enddo

      ! Scale output water fluxes to be in mm/s

      if (present(water_in)) water_in=water_in/ &
                        ((time-params%av_start_time)*hours2seconds)

      if (present(water_out)) water_out=water_out/ &
                        ((time-params%av_start_time)*hours2seconds)

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

    ! Deallocate temporary arrays if necessary

    if (allocated(albedo_temp))   deallocate(albedo_temp)
    if (allocated(if_temp))       deallocate(if_temp)
    if (allocated(wout_temp))     deallocate(wout_temp)
    if (allocated(win_temp))      deallocate(win_temp)
    if (allocated(orog_out_temp)) deallocate(orog_out_temp)

  end subroutine glimmer

!===================================================================

  subroutine end_glimmer(params)

    !*FD perform tidying-up operations for glimmer
    implicit none

    type(glimmer_params),intent(inout) :: params          !*FD parameters for this run

    integer i
    ! end individual instances

    do i=1,params%ninstances
      call glimmer_i_end(params%instances(i)%model,params%file_unit)
    enddo

    ! close log file

    close(params%logf_unit)

  end subroutine end_glimmer

!=====================================================

  integer function glimmer_coverage_map(params,coverage,cov_orog)

    !*FD Retrieve ice model fractional 
    !*FD coverage map. This function is provided so that glimmer may
    !*FD be restructured without altering the interface.
    !*RV Three return values are possible:
    !*RV \begin{description}
    !*RV \item[0 ---] Successful return
    !*RV \item[1 ---] Coverage map not calculated yet (fail)
    !*RV \item[2 ---] Coverage array is the wrong size (fail)
    !*RV \end{description}

    implicit none
    
    type(glimmer_params),intent(in) :: params       !*FD ice model parameters
    real(rk),dimension(:,:),intent(out) :: coverage !*FD array to hold coverage map
    real(rk),dimension(:,:),intent(out) :: cov_orog !*FD Orography coverage

    if (.not.params%coverage_calculated) then
      glimmer_coverage_map=1
      return
    endif

    if (size(coverage,1).ne.params%g_grid%nx.or. &
        size(coverage,2).ne.params%g_grid%ny) then
      glimmer_coverage_map=2
      return
    endif

    glimmer_coverage_map=0
    coverage=params%total_coverage
    cov_orog=params%total_cov_orog

  end function glimmer_coverage_map

!=======================================================

  integer function glimmer_main_funit(params)

    !*FD Returns the value of the 
    !*FD main logical file unit used by glimmer
    !*RV The value of the 
    !*RV main logical file unit used by glimmer.

    implicit none

    type(glimmer_params),intent(in) :: params !*FD Ice model parameters
   
    glimmer_main_funit=params%file_unit

  end function glimmer_main_funit

!=====================================================

  subroutine glimmer_set_main_funit(params,unit)

    !*FD sets the main logical 
    !*FD file unit used by glimmer.

    implicit none

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

    implicit none

    type(glimmer_params),intent(in) :: params   !*FD ice model parameters
    integer,             intent(in) :: unit     !*FD logical file unit to use
    character(*),        intent(in) :: filename !*FD filename to write

    integer :: i

    open(unit,file=filename,form='unformatted')
    
!    write(unit) params% nxg
!    write(unit) params% nyg
!    write(unit) params% glat
!    write(unit) params% glon
!    write(unit) params% glat_bound
!    write(unit) params% glon_bound
    write(unit) params% ninstances
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
    write(unit) params% cov_normalise
    write(unit) params% coverage_calculated
    write(unit) params% paramfile
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

!    if (associated(params%glat))           deallocate(params%glat)
!    if (associated(params%glon))           deallocate(params%glon)
!    if (associated(params%glat_bound))     deallocate(params%glat_bound)
!    if (associated(params%glon_bound))     deallocate(params%glon_bound)
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

!    read(unit) params% nxg
!    read(unit) params% nyg

    call glimmer_allocate_arrays(params)

!    read(unit) params% glat
!    read(unit) params% glon
!    read(unit) params% glat_bound
!    read(unit) params% glon_bound
    read(unit) params% ninstances
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
    read(unit) params% cov_normalise
    read(unit) params% coverage_calculated
    read(unit) params% paramfile
    read(unit) params% file_unit
    read(unit) params% logf_unit

    allocate(params%instances(params%ninstances))

    do i=1,params%ninstances
!      call glimmer_i_read_restart(params%instances(i),unit,params%nxg,params%nyg)
    enddo  

    close(unit)

  end subroutine glimmer_read_restart

!=============================================================================

  subroutine glimmer_set_orog_res(params,lons,lats,lonb,latb)

    !*FD Sets the output resolution of the upscaled orography, which needs
    !*FD to be different when using a spectral-transform atmosphere. If present, the boundary
    !*FD arrays \texttt{lonb} and \texttt{latb} have one more element the corresponding
    !*FD grid-point location array (\texttt{lons} and \texttt{lats}), and are indexed from 1.
    !*FD The elements are arranged such that the boundarys of \texttt{lons(i)} are found in
    !*FD \texttt{lonb(i)} and \texttt{lonb(i+1)}.

    implicit none

    type(glimmer_params),          intent(inout) :: params !*FD Ice model parameters.
    real(rk),dimension(:),         intent(in)    :: lons   !*FD Global grid longitude locations (degrees)
    real(rk),dimension(:),         intent(in)    :: lats   !*FD Global grid latitude locations (degrees)
    real(rk),dimension(:),optional,intent(in)    :: lonb   !*FD Global grid-box boundaries in longitude (degrees)
    real(rk),dimension(:),optional,intent(in)    :: latb   !*FD Global grid-box boundaries in latitude (degrees)

    integer :: i, args

    ! Reset grid variable - the call varies according to the optional arguments

    args=0

    if (present(lonb)) args=args+1
    if (present(latb)) args=args+2

    select case(args)
    case(0)
      call new_global_grid(params%g_grid_orog,lons,lats)
    case(1)
      call new_global_grid(params%g_grid_orog,lons,lats,lonb=lonb)
    case(2)
      call new_global_grid(params%g_grid_orog,lons,lats,latb=latb)
    case(3)
      call new_global_grid(params%g_grid_orog,lons,lats,lonb=lonb,latb=latb)
    end select

    ! Deallocate and reallocate coverage maps

    if (associated(params%total_cov_orog)) deallocate(params%total_cov_orog)
    if (associated(params%cov_norm_orog))  deallocate(params%cov_norm_orog)

    allocate(params%total_cov_orog(params%g_grid_orog%nx,params%g_grid_orog%ny))
    allocate(params%cov_norm_orog (params%g_grid_orog%nx,params%g_grid_orog%ny))

    ! Set total coverage to zero

    params%total_cov_orog=0.0

    ! Loop over instances

    do i=1,params%ninstances

      ! Initialise upscaling

      call new_upscale(params%instances(i)%ups_orog, &
                       params%g_grid_orog, &
                       params%instances(i)%proj, &
                       params%instances(i)%model%climate%out_mask)

      ! Deallocate fractional coverage if necessary, and reallocate

      if (associated(params%instances(i)%frac_cov_orog)) &
                             deallocate(params%instances(i)%frac_cov_orog)
      allocate(params%instances(i)%frac_cov_orog(params%g_grid_orog%nx,params%g_grid_orog%ny))

      ! Calculate fractional coverage

      call calc_coverage(params%instances(i)%proj, &
                         params%instances(i)%ups_orog,&             ! Calculate coverage map
                         params%g_grid_orog, &
                         params%radea, &
                         params%instances(i)%model%climate%out_mask, &
                         params%instances(i)%frac_cov_orog)

      ! Add to total

      params%total_cov_orog = params%total_cov_orog &
                   + params%instances(i)%frac_cov_orog

    enddo

    ! Copy to normalisation, and check not greater than 1

    params%cov_norm_orog=params%total_cov_orog
    where (params%total_cov_orog>1.0) params%total_cov_orog=1.0

  end subroutine glimmer_set_orog_res

!----------------------------------------------------------------------
! PRIVATE INTERNAL GLIMMER SUBROUTINES FOLLOW.............
!----------------------------------------------------------------------

  subroutine glimmer_allocate_arrays(params)

  !*FD allocates glimmer arrays

    implicit none

    type(glimmer_params),intent(inout) :: params !*FD ice model parameters

    allocate(params%g_av_precip (params%g_grid%nx,params%g_grid%ny))
    allocate(params%g_av_temp   (params%g_grid%nx,params%g_grid%ny))
    allocate(params%g_max_temp  (params%g_grid%nx,params%g_grid%ny))
    allocate(params%g_min_temp  (params%g_grid%nx,params%g_grid%ny))
    allocate(params%g_temp_range(params%g_grid%nx,params%g_grid%ny))
    allocate(params%g_av_zonwind(params%g_grid%nx,params%g_grid%ny))
    allocate(params%g_av_merwind(params%g_grid%nx,params%g_grid%ny))

    allocate(params%total_coverage(params%g_grid%nx,params%g_grid%ny))
    allocate(params%cov_normalise (params%g_grid%nx,params%g_grid%ny))

    allocate(params%total_cov_orog(params%g_grid_orog%nx,params%g_grid_orog%ny))
    allocate(params%cov_norm_orog (params%g_grid_orog%nx,params%g_grid_orog%ny))

  end subroutine glimmer_allocate_arrays

!====================================================

  subroutine glimmer_init_common(unit,params)

    !*FD Initialises common
    !*FD variables that need to be initialised, regardless
    !*FD of whether a restart is being performed. Also opens the log file

    use paramets, only: tau0,vel0,vis0,len0
    use physcon,  only: gn

    implicit none

    integer :: unit !*FD file unit to use for the log file
    type(glimmer_params) :: params !*FD ice model parameters

    ! THIS IS WHERE INITIALISATION THAT IS NEEDED REGARDLESS
    ! OF WHETHER WE ARE RESTARTING BELONGS!!!

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

    ninst=params%ninstances
    tinc=params%tstep_main

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
