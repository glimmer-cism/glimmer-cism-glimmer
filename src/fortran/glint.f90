

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

module glint_main

  !*FD  This is the main glimmer module, which contains the top-level 
  !*FD  subroutines and derived types comprising the glimmer ice model.

  use glimmer_global
  use glint_type
  use glint_global_grid
  use glint_constants

  ! ------------------------------------------------------------
  ! GLIMMER_PARAMS derived type definition
  ! This is where default values are set.
  ! ------------------------------------------------------------

  type glint_params 

     !*FD Derived type containing parameters relevant to all instances of 
     !*FD the model - i.e. those parameters which pertain to the global model. 

     ! Global grids used ----------------------------------------

     type(global_grid) :: g_grid      !*FD The main global grid, used for 
                                      !*FD input and most outputs
     type(global_grid) :: g_grid_orog !*FD Global grid used for orography output.

     ! Ice model instances --------------------------------------

     integer                                   :: ninstances = 1       !*FD Number of ice model instances
     type(glint_instance),pointer,dimension(:) :: instances  => null() !*FD Array of glimmer\_instances

     ! Global model parameters ----------------------------------

     real(rk) :: tstep_mbal = 1.0    !*FD Mass-balance timestep (hours)

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
     real(rk),pointer,dimension(:,:) :: g_av_humid   => null() !*FD globally averaged humidity (%)
     real(rk),pointer,dimension(:,:) :: g_av_lwdown  => null() !*FD globally averaged downwelling longwave (W/m^2)
     real(rk),pointer,dimension(:,:) :: g_av_swdown  => null() !*FD globally averaged downwelling shortwave (W/m^2)
     real(rk),pointer,dimension(:,:) :: g_av_airpress => null() !*FD globally averaged surface air pressure (Pa)

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

     ! Start flag -----------------------------------------------

     logical :: first = .true. !*FD Set if this is the first call to glimmer - make sure we set up
                               !*FD start times correctly.

     ! Accumulation/averaging flags -----------------------------

     logical :: need_winds=.false. !*FD Set if we need the winds to be accumulated/downscaled
     logical :: enmabal=.false.    !*FD Set if we're using the energy balance mass balance model anywhere

  end type glint_params

  ! Private names -----------------------------------------------

  private glint_allocate_arrays
  private glint_readconfig, calc_bounds

contains

  subroutine initialise_glint(params,lats,longs,paramfile,latb,lonb,orog,albedo, &
       ice_frac,orog_lats,orog_longs,orog_latb,orog_lonb,output_flag,daysinyear)

    !*FD Initialises the model

    use glint_proj
    use glimmer_config
    use glint_initialise
    use glimmer_log
    implicit none

    ! Subroutine argument declarations --------------------------------------------------------

    type(glint_params),              intent(inout) :: params      !*FD parameters to be set
    real(rk),dimension(:),           intent(in)    :: lats,longs  !*FD location of gridpoints 
                                                                  !*FD in global data.
    character(fname_length),         intent(in)    :: paramfile   !*FD name of file containing 
                                                                  !*FD parameters for all required 
                                                                  !*FD instances. 
    real(rk),dimension(:),  optional,intent(in)    :: latb        !*FD Locations of the latitudinal 
                                                                  !*FD boundaries of the grid-boxes.
    real(rk),dimension(:),  optional,intent(in)    :: lonb        !*FD Locations of the longitudinal
                                                                  !*FD boundaries of the grid-boxes.
    real(rk),dimension(:,:),optional,intent(out)   :: orog        !*FD Initial global orography
    real(rk),dimension(:,:),optional,intent(out)   :: albedo      !*FD Initial albedo
    real(rk),dimension(:,:),optional,intent(out)   :: ice_frac    !*FD Initial ice fraction 
    real(rk),dimension(:),  optional,intent(in)    :: orog_lats   !*FD Latitudinal location of gridpoints 
                                                                  !*FD for global orography output.
    real(rk),dimension(:),  optional,intent(in)    :: orog_longs  !*FD Longitudinal location of gridpoints 
                                                                  !*FD for global orography output.
    real(rk),dimension(:),  optional,intent(in)    :: orog_latb   !*FD Locations of the latitudinal 
                                                                  !*FD boundaries of the grid-boxes (orography).
    real(rk),dimension(:),  optional,intent(in)    :: orog_lonb   !*FD Locations of the longitudinal
                                                                  !*FD boundaries of the grid-boxes (orography).
    logical,                optional,intent(out)   :: output_flag !*FD Flag to show output set (provided for
                                                                  !*FD consistency)
    integer,                optional,intent(in)    :: daysinyear  !*FD Number of days in the year

    ! Internal variables -----------------------------------------------------------------------

    type(ConfigSection), pointer :: global_config, instance_config, section  ! configuration stuff
    character(len=100) :: message                 ! For log-writing
    character(fname_length):: instance_fname      ! name of instance specific configuration file
    integer :: i,args,o_args
    real(rk),dimension(:,:),allocatable :: orog_temp,if_temp,alb_temp ! Temporary output arrays
    integer,dimension(:),allocatable :: mbts ! Array of mass-balance timesteps

    ! Initialise year-length -------------------------------------------------------------------

    if (present(daysinyear)) then
       call glint_set_year_length(daysinyear)
    end if

    ! Initialise main global grid --------------------------------------------------------------

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

    ! Initialise orography grid ------------------------- -----------

    o_args=0

    if (present(orog_longs)) o_args=o_args+1
    if (present(orog_lats)) o_args=o_args+2
    if (present(orog_lonb)) o_args=o_args+4
    if (present(orog_latb)) o_args=o_args+8

    select case(o_args)
    case(0)
       ! Copy the main grid
       call copy_global_grid(params%g_grid,params%g_grid_orog)
    case(3)
       ! Initialise with grid points only
       call new_global_grid(params%g_grid_orog,orog_longs,orog_lats)
    case(7)
       ! Initialise with grid points and longitudinal boundaries
       call new_global_grid(params%g_grid_orog,orog_longs,orog_lats,lonb=orog_lonb)
    case(11)
       ! Initialise with grid points and latitudinal boundaries
       call new_global_grid(params%g_grid_orog,orog_longs,orog_lats,latb=orog_latb)
    case(15)
       ! Initialise with grid points and both sets of boundaries
       call new_global_grid(params%g_grid_orog,orog_longs,orog_lats,lonb=orog_lonb,latb=orog_latb)
    case default
       ! Unexpected combination
       call write_log('Unexpected combination of arguments to initialise_glint',GM_FATAL,__FILE__,__LINE__)
    end select

    ! Allocate arrays -----------------------------------------------

    call glint_allocate_arrays(params)

    ! Initialise arrays ---------------------------------------------

    params%g_av_precip  = 0.0
    params%g_av_temp    = 0.0
    params%g_max_temp   = -1000.0
    params%g_min_temp   = 1000.0
    params%g_temp_range = 0.0
    params%g_av_zonwind = 0.0
    params%g_av_merwind = 0.0
    params%g_av_humid   = 0.0
    params%g_av_lwdown  = 0.0
    params%g_av_swdown  = 0.0
    params%g_av_airpress = 0.0

    ! ---------------------------------------------------------------
    ! Open the global configuration file and read in the parameters
    ! ---------------------------------------------------------------

    call ConfigRead(paramfile,global_config)    ! Load the configuration file into the linked list
    call glint_readconfig(params,global_config) ! Parse the list

    ! Allocate array of glimmer instances
    ! and the array of mass-balance timesteps

    allocate(params%instances(params%ninstances))
    allocate(mbts(params%ninstances))

    ! ---------------------------------------------------------------
    ! Zero coverage maps and normalisation fields for main grid and
    ! orography grid
    ! ---------------------------------------------------------------

    params%total_coverage=0.0
    params%total_cov_orog=0.0

    params%cov_normalise=0.0
    params%cov_norm_orog=0.0

    ! ---------------------------------------------------------------
    ! Read config files, and initialise instances accordingly
    ! ---------------------------------------------------------------

    call write_log('Reading instance configurations')
    call write_log('-------------------------------')

    ! First, see if there are multiple instances --------------------

    call GetSection(global_config,section,'GLINT instance')

    if (.not.associated(section)) then

       ! If there aren't any GLINT instance sections, then we must
       ! only have one instance. Otherwise, flag an error.

       if (params%ninstances.gt.1) then
          write(message,*) 'Must specify ',params%ninstances,' instance config files'
          call write_log(message,GM_FATAL,__FILE__,__LINE__)
       end if

       ! In this situation, we write the name of the parameter file, and
       ! initialise the single instance

       call write_log(trim(paramfile))
       call glint_i_initialise(global_config,params%instances(1),params%g_grid,params%g_grid_orog, &
            mbts(1),params%need_winds,params%enmabal)
       call write_log('')

       ! Update the coverage and normalisation fields


       params%total_coverage = params%instances(1)%frac_coverage
       params%total_cov_orog = params%instances(1)%frac_cov_orog

       where (params%total_coverage>0.0) params%cov_normalise=1.0
       where (params%total_cov_orog>0.0) params%cov_norm_orog=1.0
 
    else

       ! Otherwise, loop through instances

       do i=1,params%ninstances

          ! The section contains one element, 'name', which
          ! points to the config file for an individual instance.
          ! This is read in to instance_config, and passed to the initialisation
          ! for the instance

          instance_fname = ''
          call GetValue(section,'name',instance_fname)
          call write_log(trim(instance_fname))
          call ConfigRead(instance_fname,instance_config)
          call glint_i_initialise(instance_config,params%instances(i),params%g_grid,params%g_grid_orog, &
               mbts(i),params%need_winds,params%enmabal)

          params%total_coverage = params%total_coverage + params%instances(i)%frac_coverage
          params%total_cov_orog = params%total_cov_orog + params%instances(i)%frac_cov_orog

          where (params%total_coverage>0.0) params%cov_normalise=params%cov_normalise+1.0
          where (params%total_cov_orog>0.0) params%cov_norm_orog=params%cov_norm_orog+1.0

          ! Get the next GLINT instance section, and if not present, flag an error

          !***** I (IR) think there's a bug here - won't the code always try and find one
          !***** more section than there needs to be? Need to test and possibly
          !***** think of a way of fixing this.

          call GetSection(section%next,section,'GLINT instance')
          if (.not.associated(section)) then
             write(message,*) 'Must specify ',params%ninstances,' instance config files'
             call write_log(message,GM_FATAL,__FILE__,__LINE__)
          end if
          ! check if we used all sections
          call CheckSections(instance_config)
       end do
    end if
    ! check if we used all sections
    call CheckSections(global_config)

    ! Check that all mass-balance time-steps are the same length and 
    ! assign that value to the top-level variable

    params%tstep_mbal=check_mbts(mbts)

    ! Check we don't have coverage greater than one at any point.

    where (params%total_coverage>1.0) params%total_coverage=1.0
    where (params%total_cov_orog>1.0) params%total_cov_orog=1.0
    params%coverage_calculated=.true.

    ! Zero optional outputs, if present

    if (present(orog))     orog=0.0
    if (present(albedo))   albedo=0.0
    if (present(ice_frac)) ice_frac=0.0

    ! Get initial fields from instances, splice together and return
    ! First, orography

    if (present(orog)) then
       allocate(orog_temp(size(orog,1),size(orog,2)))
       do i=1,params%ninstances
          call get_i_upscaled_fields(params%instances(i),orog=orog_temp)
          orog=splice_field(orog,orog_temp,params%instances(i)%frac_cov_orog, &
               params%cov_norm_orog)
       enddo
       deallocate(orog_temp)
    endif

    ! Then albedo and ice fraction

    if (present(albedo).or.present(ice_frac)) then
       ! Allocate temporary upscaled arrays for passing to instances
       if (present(albedo)) then
          allocate(if_temp(size(albedo,1),size(albedo,2)))
          allocate(alb_temp(size(albedo,1),size(albedo,2)))
       else
          allocate(if_temp(size(ice_frac,1),size(ice_frac,2)))
          allocate(alb_temp(size(ice_frac,1),size(ice_frac,2)))
       endif
       ! Loop through instances in turn and get fields
       do i=1,params%ninstances
          call get_i_upscaled_fields(params%instances(i),ice_frac=if_temp,albedo=alb_temp)
          if (present(albedo)) then
             albedo=splice_field(albedo,alb_temp,params%instances(i)%frac_coverage, &
                  params%cov_normalise)
          endif
          if (present(ice_frac)) then
             ice_frac=splice_field(ice_frac,if_temp,params%instances(i)%frac_coverage, &
                  params%cov_normalise)
          endif
       enddo
       deallocate(if_temp,alb_temp)
    endif

    if (present(output_flag)) output_flag=.true.

  end subroutine initialise_glint

  !================================================================================

  subroutine glint(params,time,temp,precip,zonwind,merwind,orog, &
       humid,lwdown,swdown,airpress, &
       output_flag,orog_out,albedo,ice_frac,water_in, &
       water_out,total_water_in,total_water_out, &
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
    !*FD
    !*FD Note that the total ice volume returned is the total at the end of the time-step;
    !*FD the water fluxes are valid over the duration of the timestep. Thus the difference
    !*FD between \texttt{total\_water\_in} and \texttt{total\_water\_out} should be equal
    !*FD to the change in \texttt{ice\_volume}, after conversion between m$^3$ and kg.

    use glimmer_utils
    use glint_interp
    use glint_timestep
    use glimmer_log
    implicit none

    ! Subroutine argument declarations -------------------------------------------------------------

    type(glint_params),              intent(inout) :: params          !*FD parameters for this run
    integer,                         intent(in)    :: time            !*FD Current model time        (hours)
    real(rk),dimension(:,:),         intent(in)    :: temp            !*FD Surface temperature field (celcius)
    real(rk),dimension(:,:),         intent(in)    :: precip          !*FD Precipitation rate        (mm/s)
    real(rk),dimension(:,:),         intent(in)    :: zonwind,merwind !*FD Zonal and meridional components 
    !*FD of the wind field         (m/s)
    real(rk),dimension(:,:),         intent(inout) :: orog            !*FD The large-scale orography (m)
    real(rk),dimension(:,:),optional,intent(in)    :: humid           !*FD Surface humidity (%)
    real(rk),dimension(:,:),optional,intent(in)    :: lwdown          !*FD Downwelling longwave (W/m^2)
    real(rk),dimension(:,:),optional,intent(in)    :: swdown          !*FD Downwelling shortwave (W/m^2)
    real(rk),dimension(:,:),optional,intent(in)    :: airpress        !*FD surface air pressure (Pa)
    logical,                optional,intent(out)   :: output_flag     !*FD Set true if outputs set
    real(rk),dimension(:,:),optional,intent(inout) :: orog_out        !*FD The fed-back, output orography (m)
    real(rk),dimension(:,:),optional,intent(inout) :: albedo          !*FD surface albedo
    real(rk),dimension(:,:),optional,intent(inout) :: ice_frac        !*FD grid-box ice-fraction
    real(rk),dimension(:,:),optional,intent(inout) :: water_in        !*FD Input water flux          (mm)
    real(rk),dimension(:,:),optional,intent(inout) :: water_out       !*FD Output water flux         (mm)
    real(rk),               optional,intent(inout) :: total_water_in  !*FD Area-integrated water flux in (kg)
    real(rk),               optional,intent(inout) :: total_water_out !*FD Area-integrated water flux out (kg)
    real(rk),               optional,intent(inout) :: ice_volume      !*FD Total ice volume (m$^3$)

    ! Internal variables ----------------------------------------------------------------------------

    integer :: i
    real(rk),dimension(:,:),allocatable :: albedo_temp,if_temp,wout_temp,orog_out_temp,win_temp
    real(rk) :: twin_temp,twout_temp,icevol_temp
    type(output_flags) :: out_f

    ! Check we have necessary input fields ----------------------------------------------------------

    if (params%enmabal) then
       if (.not.(present(humid).and.present(lwdown).and. &
            present(swdown).and.present(airpress))) &
            call write_log('Necessary fields not supplied for Energy Balance Mass Balance model',GM_FATAL,__FILE__,__LINE__)
    end if

    ! Set averaging start if necessary and return if this is not a mass-balance timestep
    ! Still not sure if this is the correct solution, but should prevent averaging of one-
    ! too-many steps at first

    if (params%first) then
       params%av_start_time=time
       params%first=.false.
       return
    endif

    ! Reset output flag

    if (present(output_flag)) output_flag=.false.

    ! ---------------------------------------------------------
    ! Do averaging and so on...
    ! Averages
    ! ---------------------------------------------------------

    params%g_av_temp    = params%g_av_temp    + temp
    params%g_av_precip  = params%g_av_precip  + precip

    if (params%need_winds) params%g_av_zonwind = params%g_av_zonwind + zonwind
    if (params%need_winds) params%g_av_merwind = params%g_av_merwind + merwind

    if (params%enmabal) then
       params%g_av_humid    = params%g_av_humid    + humid
       params%g_av_lwdown   = params%g_av_lwdown   + lwdown
       params%g_av_swdown   = params%g_av_swdown   + swdown
       params%g_av_airpress = params%g_av_airpress + airpress
    endif

    ! Ranges of temperature

    where (temp > params%g_max_temp) params%g_max_temp=temp
    where (temp < params%g_min_temp) params%g_min_temp=temp

    ! Increment step counter

    params%av_steps=params%av_steps+1

    ! ---------------------------------------------------------
    ! If this is a timestep, prepare global fields, and do a timestep
    ! for each model instance
    ! ---------------------------------------------------------

    if (time-params%av_start_time.ge.params%tstep_mbal) then

       ! Set output_flag

       ! At present, outputs are done for each mass-balance timestep, since
       ! that involved least change to the code. However, it might be good
       ! to change the output to occur with user-specified frequency.

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

       params%g_av_temp    = params%g_av_temp   /real(params%av_steps)
       params%g_av_precip  = params%g_av_precip /real(params%av_steps)
       if (params%need_winds) params%g_av_zonwind = params%g_av_zonwind/real(params%av_steps)
       if (params%need_winds) params%g_av_merwind = params%g_av_merwind/real(params%av_steps)
       if (params%enmabal) then
          params%g_av_humid    = params%g_av_humid   /real(params%av_steps)
          params%g_av_lwdown   = params%g_av_lwdown  /real(params%av_steps)
          params%g_av_swdown   = params%g_av_swdown  /real(params%av_steps)
          params%g_av_airpress = params%g_av_airpress/real(params%av_steps)
       endif

       ! Calculate total accumulated precipitation - multiply
       ! by time since last model timestep

       params%g_av_precip = params%g_av_precip*(time-params%av_start_time)*hours2seconds

       ! Calculate temperature half-range

       params%g_temp_range=(params%g_max_temp-params%g_min_temp)/2.0

       ! Do a timestep for each instance

       do i=1,params%ninstances
          call glint_i_tstep(time,&
               params%instances(i),          &
               params%g_av_temp,             &
               params%g_temp_range,          &
               params%g_av_precip,           &
               params%g_av_zonwind,          &
               params%g_av_merwind,          &
               params%g_av_humid,            &
               params%g_av_lwdown,           &
               params%g_av_swdown,           &
               params%g_av_airpress,         &
               orog,                         &
               orog_out_temp,                &
               albedo_temp,                  &
               if_temp,                      &
               win_temp,                     &
               wout_temp,                    &
               twin_temp,                    &
               twout_temp,                   &
               icevol_temp,                  &
               out_f,                        &
               .true.)

          ! Add this contribution to the output orography

          if (present(orog_out)) &
               orog_out=splice_field(orog_out, &
               orog_out_temp, &
               params%instances(i)%frac_cov_orog, &
               params%cov_norm_orog)

          if (present(albedo)) &
               albedo=splice_field(albedo,&
               albedo_temp, &
               params%instances(i)%frac_coverage, &
               params%cov_normalise)

          if (present(ice_frac)) &
               ice_frac=splice_field(ice_frac, &
               if_temp, &
               params%instances(i)%frac_coverage, &
               params%cov_normalise)

          if (present(water_in)) &
               water_in=splice_field(water_in, &
               win_temp, &
               params%instances(i)%frac_coverage, &
               params%cov_normalise)

          if (present(water_out)) &
               water_out=splice_field(water_out, &
               wout_temp, &
               params%instances(i)%frac_coverage, &
               params%cov_normalise)

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
       params%g_av_humid   = 0.0
       params%g_av_lwdown  = 0.0
       params%g_av_swdown  = 0.0
       params%g_av_airpress = 0.0
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

  end subroutine glint

  !===================================================================

  subroutine end_glint(params)

    !*FD perform tidying-up operations for glimmer
    use glint_initialise
    implicit none

    type(glint_params),intent(inout) :: params          !*FD parameters for this run

    integer i
    ! end individual instances

    do i=1,params%ninstances
       call glint_i_end(params%instances(i))
    enddo

  end subroutine end_glint

  !=====================================================

  integer function glint_coverage_map(params,coverage,cov_orog)

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

    type(glint_params),intent(in) :: params         !*FD ice model parameters
    real(rk),dimension(:,:),intent(out) :: coverage !*FD array to hold coverage map
    real(rk),dimension(:,:),intent(out) :: cov_orog !*FD Orography coverage

    if (.not.params%coverage_calculated) then
       glint_coverage_map=1
       return
    endif

    if (size(coverage,1).ne.params%g_grid%nx.or. &
         size(coverage,2).ne.params%g_grid%ny) then
       glint_coverage_map=2
       return
    endif

    glint_coverage_map=0
    coverage=params%total_coverage
    cov_orog=params%total_cov_orog

  end function glint_coverage_map

  !----------------------------------------------------------------------
  ! PRIVATE INTERNAL GLIMMER SUBROUTINES FOLLOW.............
  !----------------------------------------------------------------------

  subroutine glint_allocate_arrays(params)

    !*FD allocates glimmer arrays

    implicit none

    type(glint_params),intent(inout) :: params !*FD ice model parameters

    allocate(params%g_av_precip (params%g_grid%nx,params%g_grid%ny))
    allocate(params%g_av_temp   (params%g_grid%nx,params%g_grid%ny))
    allocate(params%g_max_temp  (params%g_grid%nx,params%g_grid%ny))
    allocate(params%g_min_temp  (params%g_grid%nx,params%g_grid%ny))
    allocate(params%g_temp_range(params%g_grid%nx,params%g_grid%ny))
    allocate(params%g_av_zonwind(params%g_grid%nx,params%g_grid%ny))
    allocate(params%g_av_merwind(params%g_grid%nx,params%g_grid%ny))
    allocate(params%g_av_humid  (params%g_grid%nx,params%g_grid%ny))
    allocate(params%g_av_lwdown (params%g_grid%nx,params%g_grid%ny))
    allocate(params%g_av_swdown (params%g_grid%nx,params%g_grid%ny))
    allocate(params%g_av_airpress(params%g_grid%nx,params%g_grid%ny))

    allocate(params%total_coverage(params%g_grid%nx,params%g_grid%ny))
    allocate(params%cov_normalise (params%g_grid%nx,params%g_grid%ny))

    allocate(params%total_cov_orog(params%g_grid_orog%nx,params%g_grid_orog%ny))
    allocate(params%cov_norm_orog (params%g_grid_orog%nx,params%g_grid_orog%ny))

  end subroutine glint_allocate_arrays

  !========================================================

  function splice_field(global,local,coverage,normalise)

    !*FD Splices an upscaled field into a global field

    real(rk),dimension(:,:),intent(in) :: global  !*FD Field to receive the splice
    real(rk),dimension(:,:),intent(in) :: local    !*FD The field to be spliced in
    real(rk),dimension(:,:),intent(in) :: coverage  !*FD The coverage fraction
    real(rk),dimension(:,:),intent(in) :: normalise  !*FD The normalisation field

    real(rk),dimension(size(global,1),size(global,2)) :: splice_field

    where (coverage==0.0)
       splice_field=global
    elsewhere
       splice_field=(global*(1-coverage/normalise))+(local*coverage/normalise)
    end where

  end function splice_field

  !========================================================

  subroutine glint_readconfig(params,config)

    !*FD read global parameters for GLINT model

    use glimmer_config
    use glimmer_log
    implicit none

    ! Arguments -------------------------------------------------------------

    type(glint_params),intent(inout) :: params !*FD ice model parameters
    type(ConfigSection), pointer :: config     !*FD structure holding sections of configuration file

    ! Internal variables ----------------------------------------------------

    type(ConfigSection), pointer :: section
    character(len=100) :: message

    ! -----------------------------------------------------------------------
    ! If there's a section called 'GLINT' in the config file,
    ! then we use that information to overwrite the defaults.
    ! Otherwise, it's one instance with a timestep of one year.

    call GetSection(config,section,'GLINT')
    if (associated(section)) then
       call GetValue(section,'n_instance',params%ninstances)
    end if

    ! Print some configuration information

    call write_log('GLINT global')
    call write_log('------------')
    write(message,*) 'number of instances :',params%ninstances
    call write_log(message)
    call write_log('')

  end subroutine glint_readconfig

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


  !========================================================

  real(rk) function check_mbts(timesteps)

    !*FD Checks to see that all mass-balance time-steps are
    !*FD the same. Flags a fatal error if not, else assigns that
    !*FD value to the output

    use glimmer_log

    implicit none

    integer,dimension(:) :: timesteps !*FD Array of mass-balance timsteps

    integer :: n,i

    n=size(timesteps)
    if (n==0) then
       check_mbts=0
       return
    endif

    check_mbts=timesteps(1)

    do i=2,n
       if (timesteps(i)/=check_mbts) then
          call write_log('All mass-balance schemes must have the same timestep', &
               GM_FATAL,__FILE__,__LINE__)
       endif
    enddo

  end function check_mbts


end module glint_main

