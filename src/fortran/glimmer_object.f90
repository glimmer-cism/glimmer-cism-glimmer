
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_object.f90 - part of the GLIMMER ice model       + 
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

module glimmer_object

  !*FD Holds the dervied type for the ice model
  !*FD instances, and the necessary code for calling them.

  use glimmer_global
  use glimmer_project
  use glimmer_interp
  use glimmer_types

  implicit none

  type glimmer_instance

    !*FD Derived type holding information about ice model instance. 

    type(projection)                 :: proj               !*FD The projection definition of the instance.
    type(downscale)                  :: downs              !*FD Downscaling parameters.
    type(upscale)                    :: ups                !*FD Upscaling parameters
    type(upscale)                    :: ups_orog           !*FD Upscaling parameters for orography (to cope
                                                           !*FD with need to convert to spectral form).
    type(glimmer_global_type)        :: model              !*FD The instance and all its arrays.
    character(fname_length)          :: paramfile          !*FD The name list file of parameters.
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

    ! Accumulation information --------------------------------------------------

    real(rk) :: accum_start = 0.0 ! Time when mass-balance accumulation started
    logical  :: first_accum = .true.  ! First time accumulation (difference from first, since
                                      ! that relates to first dynamics step)

  end type glimmer_instance

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

  subroutine glimmer_i_initialise(unit,nmlfile,instance,radea,grid,time_step,tstep_mbal,start_time)

    !*FD Initialise an ice model (glimmer) instance

    use glimmer_setup
    use glimmer_temp
    use glimmer_velo
    use glimmer_outp
    use glimmer_global_grid
    use glimmer_ncfile

    ! Arguments

    integer,               intent(in)    :: unit        !*FD Filename unit to use when opening namelist file.
    character(*),          intent(in)    :: nmlfile     !*FD Name of namelist file.
    type(glimmer_instance),intent(inout) :: instance    !*FD The instance being initialised.
    real(rk),              intent(in)    :: radea       !*FD Radius of the earth (m).
    type(global_grid),     intent(in)    :: grid        !*FD Global grid to use
    real(rk),              intent(in)    :: time_step   !*FD Model time-step (years).
    real(rk),              intent(out)   :: tstep_mbal  !*FD Mass-balance time-step (years).
    real(rk),optional,     intent(in)    :: start_time  !*FD Start time of model (years).

    ! Internal variables

    instance%model%numerics%tinc=time_step            ! Initialise the model time step

    call read_config_file(unit,nmlfile,instance%model,instance%proj) ! Read config file and initialise variables
    call new_proj(instance%proj,radea)                          ! Initialise the projection
    call new_downscale(instance%downs,instance%proj,grid)       ! Initialise the downscaling
    call glimmer_i_allocate(instance,grid%nx,grid%ny)           ! Allocate arrays appropriately
    call glimmer_load_sigma(instance%model,unit)                ! Load the sigma file
    call calc_lats(instance%proj,instance%model%climate%lati)   ! Initialise the local latitude array. 
                                                                ! This may be redundant, though.
    call testinisthk(instance%model,unit,instance%first, &
	                 instance%newtemps,real(instance%global_orog,rk)) ! Load in initial surfaces, and other 2d fields
    call new_upscale(instance%ups,grid,instance%proj, &
                     instance%model%climate%out_mask)           ! Initialise upscaling parameters
    call calc_coverage(instance%proj, &                         ! Calculate coverage map
                       instance%ups,  &             
                       grid,          &
                       radea,         &
                       instance%model%climate%out_mask, &
                       instance%frac_coverage)
    call copy_upscale(instance%ups,instance%ups_orog)    ! Set upscaling for orog output to same as for 
                                                         ! other fields.
    instance%frac_cov_orog=instance%frac_coverage        ! Set fractional coverage for orog to be same as
                                                         ! for other fields.

    if (present(start_time)) then
      instance%model%numerics%time = start_time       ! Initialise the counter.
    else                                              ! Despite being in the GLIMMER framework,
      instance%model%numerics%time = 0.0              ! each instance has a copy of the counter
    endif                                             ! for simplicity.

    call openall_out(instance%model)                  ! Initialise output files
    call writeall(instance%model)

    tstep_mbal=instance%model%numerics%tinc_mbal      ! Initialise the mass-balance timestep

  end subroutine glimmer_i_initialise

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glimmer_i_allocate(instance,nxg,nyg)

    !*FD Allocate top-level arrays in
    !*FD the model instance, and ice model arrays.

    use glimmer_setup

    type(glimmer_instance),intent(inout) :: instance    !*FD Instance whose elements are to be allocated.
    integer,               intent(in)    :: nxg         !*FD Longitudinal size of global grid (grid-points).
    integer,               intent(in)    :: nyg         !*FD Latitudinal size of global grid (grid-points).

    ! First deallocate if necessary

    if (associated(instance%xwind))         deallocate(instance%xwind)
    if (associated(instance%ywind))         deallocate(instance%ywind)
    if (associated(instance%global_orog))   deallocate(instance%global_orog) 
    if (associated(instance%local_orog))    deallocate(instance%local_orog)   
    if (associated(instance%frac_coverage)) deallocate(instance%frac_coverage)

    ! Then reallocate...
    ! Wind field arrays

    allocate(instance%xwind(instance%model%general%ewn,instance%model%general%nsn)) 
    allocate(instance%ywind(instance%model%general%ewn,instance%model%general%nsn))

    ! Local-global orog

    allocate(instance%global_orog(instance%model%general%ewn,instance%model%general%nsn))

    ! Local-local orog

    allocate(instance%local_orog(instance%model%general%ewn,instance%model%general%nsn))

    ! Global box indices and number of points contained therein

    ! Fractional coverage map

    allocate(instance%frac_coverage(nxg,nyg)) ! allocate fractional coverage map
    allocate(instance%frac_cov_orog(nxg,nyg)) ! allocate fractional coverage map (orog)

    ! Allocate the model arrays - use original ice model code

    call allocarr(instance%model)

  end subroutine glimmer_i_allocate

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glimmer_i_tstep(unit,logunit,time,instance,g_temp,g_temp_range, &
                          g_precip,g_zonwind,g_merwind,g_orog,g_orog_out,g_albedo,g_ice_frac,&
                          g_water_in,g_water_out,t_win,t_wout,ice_vol,out_f)

    !*FD Performs time-step of an ice model instance. Note that this 
    !*FD code will need to be altered to take account of the 
    !*FD energy-balance mass-balance model when it is completed.
    !*FD
    !*FD Note also that input quantities here are accumulated totals since the
    !*FD last call.

    use glimmer_isot
    use glimmer_thck
    use glimmer_velo
    use glimmer_temp
    use glimmer_setup
    use glimmer_ncfile
    use glimmer_interp
    use glimmer_mbal
    use glide_messages
    use paramets

    ! ------------------------------------------------------------------------  
    ! Arguments
    ! ------------------------------------------------------------------------  

    integer,                intent(in)   :: unit         !*FD Logical file unit to be used for write operations
    integer,                intent(in)   :: logunit      !*FD Unit for log file
    real(rk),               intent(in)   :: time         !*FD Current time in years
    type(glimmer_instance), intent(inout):: instance     !*FD Model instance
    real(rk),dimension(:,:),intent(in)   :: g_temp       !*FD Global mean surface temperature field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_temp_range !*FD Global surface temperature half-range field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_precip     !*FD Global precip field total (mm)
    real(rk),dimension(:,:),intent(in)   :: g_zonwind    !*FD Global mean surface zonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_merwind    !*FD Global mean surface meridonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_orog       !*FD Input global orography (m)
    real(rk),dimension(:,:),intent(out)  :: g_orog_out   !*FD Output orography (m)
    real(rk),dimension(:,:),intent(out)  :: g_albedo     !*FD Output surface albedo 
    real(rk),dimension(:,:),intent(out)  :: g_ice_frac   !*FD Output ice fraction
    real(rk),dimension(:,:),intent(out)  :: g_water_in   !*FD Input water flux (mm)
    real(rk),dimension(:,:),intent(out)  :: g_water_out  !*FD Output water flux (mm)
    real(rk),               intent(out)  :: t_win        !*FD Total water input (kg)
    real(rk),               intent(out)  :: t_wout       !*FD Total water output (kg)
    real(rk),               intent(out)  :: ice_vol      !*FD Output ice volume (m$^3$)
    type(output_flags),     intent(in)   :: out_f        !*FD Flags to tell us whether to do output       

    ! ------------------------------------------------------------------------  
    ! Internal variables
    ! ------------------------------------------------------------------------  

    real(rk),dimension(:,:),allocatable :: upscale_temp  ! temporary array for upscaling
    real(rk),dimension(:,:),allocatable :: accum_temp    ! temporary array for accumulation
    real(rk),dimension(:,:),allocatable :: ablat_temp    ! temporary array for ablation
    integer, dimension(:,:),allocatable :: fudge_mask    ! temporary array for fudging
    real(rk) :: f1 ! Scaling factor for converting precip and run-off amounts.
    character(40) :: timetxt
    real(rk) :: start_volume,end_volume,flux_fudge

    ! Set scaling factor for water flux calculations -------------------------

    f1 = scyr * thk0 / tim0

    ! ------------------------------------------------------------------------  
    ! Update internal clock
    ! ------------------------------------------------------------------------  

    instance%model%numerics%time=time  

    ! Set mass-balance accumulation start if necessary -----------------------

    if (instance%first_accum) then
      instance%accum_start=nint(time-instance%model%numerics%tinc_mbal)
      instance%first_accum=.false.
    endif

    ! ------------------------------------------------------------------------  
    ! Downscale input fields, but only if required by options selected
    ! ------------------------------------------------------------------------  

    ! Temperature downscaling ------------------------------------------------

    select case(instance%model%options%whichartm)
    case(7)
      call interp_to_local(instance%proj,  &
                           g_temp,         &
                           instance%downs, &
                           localsp=instance%model%climate%g_artm)

      call interp_to_local(instance%proj,  &
                           g_temp_range,   &
                           instance%downs, &
                           localsp=instance%model%climate%g_arng)
    end select

    ! Precip downscaling -----------------------------------------------------

    select case(instance%model%options%whichprecip)
    case(1,2)
      call interp_to_local(instance%proj,  &
                           g_precip,       &
                           instance%downs, &
                           localsp=instance%model%climate%prcp)
    end select

    ! Orography downscaling --------------------------------------------------

    call interp_to_local(instance%proj,    &
                         g_orog,           &
                         instance%downs,   &
                         localdp=instance%global_orog)

    ! Wind downscaling -------------------------------------------------------

    select case(instance%model%options%whichprecip)
    case(2)
      call interp_wind_to_local(instance%proj,    &
                                g_zonwind,        &
                                g_merwind,        &
                                instance%downs,   &
                                instance%xwind,instance%ywind)
    end select

    ! ------------------------------------------------------------------------  
    ! Sort out some local orography - scale orography by correct amount
    ! ------------------------------------------------------------------------  

    instance%local_orog=instance%model%geometry%usrf*thk0

    ! Remove bathymetry. This relies on the point 1,1 being underwater.
    ! However, it's a better method than just setting all points < 0.0 to zero

    call glimmer_remove_bath(instance%local_orog,1,1)

    ! ------------------------------------------------------------------------  
    ! Calculate or process the surface temperature and half-range
    ! ------------------------------------------------------------------------  

    call calcartm(instance%model,                       &
                  instance%model%options%whichartm,     &
                  instance%model%geometry%usrf,         &
                  instance%model%climate%lati,          &
                  instance%model%climate%artm,          &
                  instance%model%climate%arng,          &
                  g_orog=instance%global_orog,          &
                  g_artm=instance%model%climate%g_artm, &
                  g_arng=instance%model%climate%g_arng)

    ! ------------------------------------------------------------------------  
    ! Calculate or process the precipitation field 
    ! ------------------------------------------------------------------------  

    call calcprcp(instance%model%options%whichprecip, &
                  instance%model%climate%prcp, &
                  instance%model%climate%uprecip_rate, &
                  f1, &
                  instance%xwind, &
                  instance%ywind, &
                  instance%model%climate%artm, &
                  instance%local_orog, &
                  instance%proj%dx, &
                  instance%proj%dy, &
                  instance%model%climate%presprcp, &
                  instance%model%climate%presartm, &
                  instance%model%climate%pfac)

    ! ------------------------------------------------------------------------ 
    ! Calculate ablation, and thus mass-balance
    ! ------------------------------------------------------------------------ 

    call calcacab(instance%model%numerics, &
                  instance%model%paramets, &
                  instance%model%pddcalc,  &
                  instance%model%options%  whichacab, &
                  instance%model%geometry% usrf,      &
                  instance%model%climate%  artm,      &
                  instance%model%climate%  arng,      &
                  instance%model%climate%  prcp,      &
                  instance%model%climate%  ablt,      &
                  instance%model%climate%  lati,      &
                  instance%model%climate%  acab,      &
                  instance%model%geometry% thck)

    ! Accumulate mass-balance if necessary

    instance%model%climate%prcp_save = instance%model%climate%prcp_save + instance%model%climate%prcp
    instance%model%climate%ablt_save = instance%model%climate%ablt_save + instance%model%climate%ablt
    instance%model%climate%acab_save = instance%model%climate%acab_save + instance%model%climate%acab

    ! If it's not time for a dynamics/temp/velocity timestep, return

    if (time-instance%accum_start.lt.instance%model%numerics%tinc) return

    ! ------------------------------------------------------------------------  
    ! ICE TIMESTEP begins HERE ***********************************************
    ! ------------------------------------------------------------------------  
 
    ! Copy and zero accumulated totals ---------------------------------------

    instance%model%climate%prcp = instance%model%climate%prcp_save
    instance%model%climate%ablt = instance%model%climate%ablt_save
    instance%model%climate%acab = instance%model%climate%acab_save

    instance%model%climate%prcp_save = 0.0
    instance%model%climate%ablt_save = 0.0
    instance%model%climate%acab_save = 0.0

    instance%accum_start=time

    ! Calculate the initial ice volume (scaled) ------------------------------

    start_volume=sum(instance%model%geometry%thck)

    ! Allocate temporary upscaling array -------------------------------------

    allocate(upscale_temp(instance%model%general%ewn,instance%model%general%nsn))
    allocate(accum_temp(instance%model%general%ewn,instance%model%general%nsn))
    allocate(ablat_temp(instance%model%general%ewn,instance%model%general%nsn))
    allocate(fudge_mask(instance%model%general%ewn,instance%model%general%nsn))

    accum_temp=0.0 ; ablat_temp=0.0

    ! Do water budget accounting ---------------------------------------------

    if (out_f%water_out.or.out_f%total_wout.or.out_f%water_in .or.out_f%total_win) then
      accum_temp=instance%model%climate%prcp*instance%model%numerics%dt
      ablat_temp=instance%model%climate%ablt*instance%model%numerics%dt
    endif

    ! ------------------------------------------------------------------------  
    ! If first step, use as seed...
    ! ------------------------------------------------------------------------  

    if (instance%first) then

	    call timeevoltemp(instance%model,0,instance%model%climate%artm)
	    instance%newtemps = .true.          ! we have new temperatures
      
      where (instance%local_orog>0.0)
        instance%model%geometry%thck = max(0.0d0, &
                                         instance%model%climate%acab* &
                                         instance%model%numerics%dt)
      endwhere

      call calclsrf(instance%model%geometry%thck, &
                    instance%model%geometry%topg, &
                    instance%model%geometry%lsrf)

      instance%model%geometry%usrf = instance%model%geometry%thck + &
                                     instance%model%geometry%lsrf

      ! Reset mass-balance variables to zero, as we've already used these

      instance%model%climate%prcp = 0.0
      instance%model%climate%ablt = 0.0
      instance%model%climate%acab = 0.0

    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate various derivatives...
    ! ------------------------------------------------------------------------ 

    call stagvarb(instance%model%geometry% thck, &
                  instance%model%geomderv% stagthck,&
                  instance%model%general%  ewn, &
                  instance%model%general%  nsn)

    call geomders(instance%model%numerics, &
                  instance%model%geometry% usrf, &
                  instance%model%geomderv% stagthck, &
                  instance%model%geomderv% dusrfdew, &
                  instance%model%geomderv% dusrfdns)

    call geomders(instance%model%numerics, &
                  instance%model%geometry% thck, &
                  instance%model%geomderv% stagthck, &
                  instance%model%geomderv% dthckdew, &
                  instance%model%geomderv% dthckdns)

    ! ------------------------------------------------------------------------ 
    ! Do velocity calculation if necessary
    ! ------------------------------------------------------------------------ 

    if (instance%model%numerics%tinc > &
        mod(instance%model%numerics%time,instance%model%numerics%nvel) .or. &
        instance%model%numerics%time == instance%model%numerics%tinc ) then
           
       if (instance%model%options%whichevol/=2) then

          call slipvelo(instance%model%numerics, &
                        instance%model%velowk,   &
                        instance%model%paramets, &
                        instance%model%geomderv, &
                        (/instance%model%options%whichslip,&
                        instance%model%options%whichbtrc/), &
                        instance%model%temper%   bwat,     &
                        instance%model%velocity% btrc,     &
                        instance%model%geometry% relx,     &
                        instance%model%velocity% ubas,     &
                        instance%model%velocity% vbas)

          call zerovelo(instance%model%velowk,             &
                        instance%model%numerics%sigma,     &
                        0,                                 &
                        instance%model%geomderv% stagthck, &
                        instance%model%geomderv% dusrfdew, &
                        instance%model%geomderv% dusrfdns, &
                        instance%model%temper%   flwa,     &
                        instance%model%velocity% ubas,     &
                        instance%model%velocity% vbas,     &
                        instance%model%velocity% uvel,     &
                        instance%model%velocity% vvel,     &
                        instance%model%velocity% uflx,     &
                        instance%model%velocity% vflx,     &
                        instance%model%velocity% diffu)

       end if

    end if

    ! ------------------------------------------------------------------------ 
    ! Do thickness masking
    ! ------------------------------------------------------------------------ 

    call maskthck(instance%model%options%  whichthck, &
                  instance%model%geometry% thck,      &
                  instance%model%climate%  acab,      &
                  instance%model%geometry% dom,       &
                  instance%model%geometry% mask,      &
                  instance%model%geometry% totpts,    &
                  instance%model%geometry% empty)

    ! ------------------------------------------------------------------------ 
    ! Calculate temperature evolution and Glenn's A, if necessary
    ! ------------------------------------------------------------------------ 

    if ( instance%model%numerics%tinc >  &
         mod(instance%model%numerics%time,instance%model%numerics%ntem) ) then
      call timeevoltemp(instance%model, &
                        instance%model%options%whichtemp, &
                        instance%model%climate%artm)
      instance%newtemps = .true.
    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate flow evolution by various different methods
    ! ------------------------------------------------------------------------ 

    select case(instance%model%options%whichevol)
    case(0) ! Use precalculated uflx, vflx -----------------------------------

      call timeevolthck(instance%model, &
                        instance%model%options%  whichthck, &
                        instance%model%geometry% usrf,      &
                        instance%model%geometry% thck,      &
                        instance%model%geometry% lsrf,      &
                        instance%model%climate%  acab,      &
                        instance%model%geometry% mask,      &
                        instance%model%velocity% uflx,      &
                        instance%model%velocity% vflx,      &
                        instance%model%geomderv% dusrfdew,  &
                        instance%model%geomderv% dusrfdns,  &
                        instance%model%geometry% totpts,    &
                        logunit)

    case(1) ! Use explicit leap frog method with uflx,vflx -------------------

      call stagleapthck(instance%model, &
                        instance%model%velocity% uflx, &
                        instance%model%velocity% vflx, &
                        instance%model%geometry% thck, &
                        instance%model%geometry% usrf, &
                        instance%model%geometry% lsrf, &
                        instance%model%climate%  acab)

    case(2) ! Use non-linear calculation that incorporates velocity calc -----

      call nonlevolthck(instance%model, &
                        instance%model%options%whichthck, &
                        instance%newtemps,logunit)

    case default

       call glide_msg(GM_FATAL,__FILE__,__LINE__,'Unrecognised value of whichevol')

    end select

    ! ------------------------------------------------------------------------ 
    ! Remove ice which is either floating, or is present below prescribed
    ! depth, depending on value of whichmarn
    ! ------------------------------------------------------------------------ 

    call marinlim(instance%model%options%  whichmarn, &
                  instance%model%options%  whichthck, &
                  instance%model%geometry% thck,      &
                  instance%model%geometry% usrf,      &
                  instance%model%geometry% relx,      &
                  instance%model%geometry% topg,      &
                  instance%model%climate%  lati,      &
                  instance%model%numerics%mlimit,     &
                  ablat_temp)

    ! ------------------------------------------------------------------------ 
    ! Calculate isostasy
    ! ------------------------------------------------------------------------ 

    call isosevol(instance%model%numerics,            & 
                  instance%model%paramets,            &
                  instance%model%isotwk,              &                                      
                  instance%model%options%  whichisot, &
                  instance%model%geometry% thck,      &
                  instance%model%geometry% topg,      &
                  instance%model%geometry% relx)

    call calclsrf(instance%model%geometry%thck, &
                  instance%model%geometry%topg, &
                  instance%model%geometry%lsrf)

    instance%model%geometry%usrf = instance%model%geometry%thck + &
                                   instance%model%geometry%lsrf

    ! ------------------------------------------------------------------------ 
    ! Do outputs if necessary
    ! ------------------------------------------------------------------------ 

    call writeall(instance%model)

    instance%newtemps = .false.

    write(timetxt,*)instance%model%numerics%time
    call glide_msg(GM_TIMESTEP,__FILE__,__LINE__,"completed time "//trim(adjustl(timetxt)))

    ! Calculate flux fudge factor --------------------------------------------

    if (out_f%water_out.or.out_f%total_wout.or.out_f%water_in .or.out_f%total_win) then

      end_volume=sum(instance%model%geometry%thck)

      where (instance%model%geometry%thck>0.0)
        fudge_mask=1
      elsewhere
        fudge_mask=0
      endwhere

      flux_fudge=(start_volume+sum(accum_temp)-sum(ablat_temp)-end_volume)/sum(fudge_mask)

      ! Apply fudge_factor

      where(instance%model%geometry%thck>0.0)
        ablat_temp=ablat_temp+flux_fudge
      endwhere

    endif

    ! ------------------------------------------------------------------------ 
    ! Upscaling of output
    ! ------------------------------------------------------------------------ 

    ! Upscale the output orography field, and re-dimensionalise --------------

    if (out_f%orog) call get_i_upscaled_fields(instance,orog=g_orog_out)

    ! Use thickness to calculate albedo and ice fraction ---------------------

    if (out_f%albedo.or.out_f%ice_frac.or.out_f%water_out) then
      call get_i_upscaled_fields(instance,albedo=g_albedo,ice_frac=g_ice_frac)
    endif

    ! Calculate ice volume ---------------------------------------------------

    if (out_f%ice_vol) then
      ice_vol=f1*1000.0*sum(instance%model%geometry%thck)*instance%proj%dx*instance%proj%dy
    endif

    ! Upscale water flux fields ----------------------------------------------
    ! First water input (i.e. mass balance + ablation)

    if (out_f%water_in) then

      where (instance%model%geometry%thck>0.0)
        upscale_temp=f1*1000.0*accum_temp
      elsewhere
        upscale_temp=0.0
      endwhere

      call mean_to_global(instance%proj,  &
                          instance%ups,   &
                          upscale_temp,   &
                          g_water_in,     &
                          instance%model%climate%out_mask)

    endif

    ! Now water output (i.e. ablation)

    if (out_f%water_out) then

      where (instance%model%geometry%thck>0.0)
        upscale_temp=f1*1000.0*ablat_temp
      elsewhere
        upscale_temp=0.0
      endwhere

      call mean_to_global(instance%proj,  &
                          instance%ups,   &
                          upscale_temp,   &
                          g_water_out,    &
                          instance%model%climate%out_mask)

    endif

    ! Sum water fluxes and convert if necessary ------------------------------

    if (out_f%total_win) then
      t_win  = sum(accum_temp)*f1*1000.0*instance%proj%dx*instance%proj%dy
    endif

    if (out_f%total_wout) then
      t_wout = sum(ablat_temp)*f1*1000.0*instance%proj%dx*instance%proj%dy
    endif

    ! Tidy up ----------------------------------------------------------------

    deallocate(upscale_temp,accum_temp,ablat_temp,fudge_mask)
    instance%first=.false.
 
  end subroutine glimmer_i_tstep

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glimmer_i_end(model,unit)

    !*FD Performs tidying-up for an ice model. Logically, this 
    !*FD should take type {\tt glimmer\_instance} as
    !*FD input, rather than {\tt glimmer\_global\_type}, but it doesn't!

    use glimmer_setup
    use glimmer_ncfile
    use glimmer_ncinfile

    ! Arguments

    type(glimmer_global_type),intent(inout) :: model !*FD Model to be tidyed-up.
    integer,                  intent(in)    :: unit  !*FD Logical file unit to use for writing.

    ! Beginning of code

    call writeall(model,.true.)
    call closeall_out(model)
    call closeall_in(model)

  end subroutine glimmer_i_end

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_instance_params(unit,instance,filename)

    !*FD Read namelist of parameters 
    !*FD for each instance.

    implicit none

    ! Arguments

    integer,               intent(in)    :: unit     !*FD Logical file unit to use
    type(glimmer_instance),intent(inout) :: instance !*FD Instance being initialised
    character(*),          intent(in)    :: filename !*FD Filename of namelist file

    ! Internal variables

    integer  :: nx,ny
    real(rk) :: dx,dy,cpx,cpy,latc,lonc

    ! Namelist definition

    namelist /project/ nx,ny,dx,dy,cpx,cpy,latc,lonc

    ! Beginning of code

    instance%paramfile=filename

    open(unit,file=filename)
    read(unit,nml=project)
    close(unit)

    instance%proj%nx=nx
    instance%proj%ny=ny
    instance%proj%dx=dx
    instance%proj%dy=dy
    instance%proj%cpx=cpx
    instance%proj%cpy=cpy
    instance%proj%latc=latc
    instance%proj%lonc=lonc

  end subroutine get_instance_params

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_i_upscaled_fields(instance,orog,albedo,ice_frac)

    !*FD Upscales and returns certain fields, according to the 
    !*FD arguments supplied

    use paramets

    type(glimmer_instance),          intent(in)  :: instance
    real(rk),dimension(:,:),optional,intent(out) :: orog
    real(rk),dimension(:,:),optional,intent(out) :: albedo
    real(rk),dimension(:,:),optional,intent(out) :: ice_frac

    real(rk),dimension(:,:),allocatable :: if_temp,upscale_temp

	  ! Calculate orography

    if (present(orog)) then
      call mean_to_global(instance%proj, &
                          instance%ups_orog, &
                          instance%model%geometry%usrf, &
                          orog,    &
                          instance%model%climate%out_mask)
      orog=thk0*orog
    endif

    if (present(albedo).or.present(ice_frac)) then

      if (present(albedo)) then
        allocate(if_temp(size(albedo,1),size(albedo,2)))
      else
        allocate(if_temp(size(ice_frac,1),size(ice_frac,2)))
      endif
      allocate(upscale_temp(instance%model%general%ewn,instance%model%general%nsn))

      ! First, ice coverage on local grid 
  
      where (instance%model%geometry%thck>0.0)
        upscale_temp=1.0
      elsewhere
        upscale_temp=0.0
      endwhere

      ! Upscale it...

      call mean_to_global(instance%proj, &
                          instance%ups, &
                          upscale_temp, &
                          if_temp,    &
                          instance%model%climate%out_mask)

      if (present(ice_frac)) ice_frac=if_temp

    endif

    ! Calculate albedo -------------------------------------------------------

    if (present(albedo)) then 
      where (if_temp>0.0)
        albedo=instance%model%climate%ice_albedo
      elsewhere
        albedo=0.0
      endwhere
    endif

    if (allocated(if_temp)) deallocate(if_temp)
    if (allocated(upscale_temp)) deallocate(upscale_temp)

  end subroutine get_i_upscaled_fields

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine calc_coverage(proj,ups,grid,radea,mask,frac_coverage)

    !*FD Calculates the fractional
    !*FD coverage of the global grid-boxes by the ice model
    !*FD domain.

    use glimmer_project
    use glimmer_global_grid
    use gmt, only: D2R

    ! Arguments

    type(projection),       intent(in)  :: proj          !*FD Projection to be used
    type(upscale),          intent(in)  :: ups           !*FD Upscaling used
    type(global_grid),      intent(in)  :: grid          !*FD Global grid used
    real(rk),               intent(in)  :: radea         !*FD Radius of the earth (m)
    integer, dimension(:,:),intent(in)  :: mask          !*FD Mask of points for upscaling
    real(rk),dimension(:,:),intent(out) :: frac_coverage !*FD Map of fractional 
                                                         !*FD coverage of global by local grid-boxes.
    ! Internal variables

    integer,dimension(grid%nx,grid%ny) :: tempcount
    integer :: i,j

    ! Beginning of code

    tempcount=0

    do i=1,proj%nx
      do j=1,proj%ny
        tempcount(ups%gboxx(i,j),ups%gboxy(i,j))=tempcount(ups%gboxx(i,j),ups%gboxy(i,j))+mask(i,j)
      enddo
    enddo

    do i=1,grid%nx
      do j=1,grid%ny
        if (tempcount(i,j)==0) then
          frac_coverage(i,j)=0.0
        else
          frac_coverage(i,j)=(tempcount(i,j)*proj%dx*proj%dy)/ &
                 (lon_diff(grid%lon_bound(i+1),grid%lon_bound(i))*D2R*radea**2*    &
                 (sin(grid%lat_bound(j)*D2R)-sin(grid%lat_bound(j+1)*D2R)))
        endif
      enddo
    enddo

    ! Fix points that should be 1.0 by checking their surroundings

    ! Interior points first

    do i=2,grid%nx-1
      do j=2,grid%ny-1
        if ((frac_coverage(i,j).ne.0).and. &
            (frac_coverage(i+1,j).ne.0).and. &
            (frac_coverage(i,j+1).ne.0).and. &
            (frac_coverage(i-1,j).ne.0).and. &
            (frac_coverage(i,j-1).ne.0)) &
                        frac_coverage(i,j)=1.0
      enddo
    enddo

    ! top and bottom edges

    do i=2,grid%nx/2
      if ((frac_coverage(i,1).ne.0).and. &
          (frac_coverage(i+1,1).ne.0).and. &
          (frac_coverage(i,2).ne.0).and. &
          (frac_coverage(i-1,1).ne.0).and. &
          (frac_coverage(i+grid%nx/2,1).ne.0)) &
                      frac_coverage(i,1)=1.0
    enddo

    do i=grid%nx/2+1,grid%nx-1
      if ((frac_coverage(i,1).ne.0).and. &
          (frac_coverage(i+1,1).ne.0).and. &
          (frac_coverage(i,2).ne.0).and. &
          (frac_coverage(i-1,1).ne.0).and. &
          (frac_coverage(i-grid%nx/2,1).ne.0)) &
                      frac_coverage(i,1)=1.0
    enddo

    do i=2,grid%nx/2
      if ((frac_coverage(i,grid%ny).ne.0).and. &
          (frac_coverage(i+1,grid%ny).ne.0).and. &
          (frac_coverage(i+grid%nx/2,grid%ny).ne.0).and. &
          (frac_coverage(i-1,grid%ny).ne.0).and. &
          (frac_coverage(i,grid%ny-1).ne.0)) &
                      frac_coverage(i,grid%ny)=1.0
    enddo

    do i=grid%nx/2+1,grid%nx-1
      if ((frac_coverage(i,grid%ny).ne.0).and. &
          (frac_coverage(i+1,grid%ny).ne.0).and. &
          (frac_coverage(i-grid%nx/2,grid%ny).ne.0).and. &
          (frac_coverage(i-1,grid%ny).ne.0).and. &
          (frac_coverage(i,grid%ny-1).ne.0)) &
                      frac_coverage(i,grid%ny)=1.0
    enddo
 
    ! left and right edges

    do j=2,grid%ny-1
      if ((frac_coverage(1,j).ne.0).and. &
          (frac_coverage(2,j).ne.0).and. &
          (frac_coverage(1,j+1).ne.0).and. &
          (frac_coverage(grid%nx,j).ne.0).and. &
          (frac_coverage(1,j-1).ne.0)) &
                      frac_coverage(1,j)=1.0
      if ((frac_coverage(grid%nx,j).ne.0).and. &
          (frac_coverage(1,j).ne.0).and. &
          (frac_coverage(grid%nx,j+1).ne.0).and. &
          (frac_coverage(grid%nx-1,j).ne.0).and. &
          (frac_coverage(grid%nx,j-1).ne.0)) &
                      frac_coverage(grid%nx,j)=1.0
    enddo

    ! corners

    if ((frac_coverage(1,1).ne.0).and. &
        (frac_coverage(2,1).ne.0).and. &
        (frac_coverage(1,2).ne.0).and. &
        (frac_coverage(grid%nx,1).ne.0).and. &
        (frac_coverage(grid%nx/2+1,1).ne.0)) &
                    frac_coverage(1,1)=1.0

    if ((frac_coverage(1,grid%ny).ne.0).and. &
        (frac_coverage(2,grid%ny).ne.0).and. &
        (frac_coverage(grid%nx/2+1,grid%ny).ne.0).and. &
        (frac_coverage(grid%nx,grid%ny).ne.0).and. &
        (frac_coverage(1,grid%ny-1).ne.0)) &
                    frac_coverage(1,grid%ny)=1.0

    if ((frac_coverage(grid%nx,1).ne.0).and. &
        (frac_coverage(1,1).ne.0).and. &
        (frac_coverage(grid%nx,2).ne.0).and. &
        (frac_coverage(grid%nx-1,1).ne.0).and. &
        (frac_coverage(grid%nx/2,1).ne.0)) &
                   frac_coverage(grid%nx,1)=1.0

    if ((frac_coverage(grid%nx,grid%ny).ne.0).and. &
        (frac_coverage(1,grid%ny).ne.0).and. &
        (frac_coverage(grid%nx/2,grid%ny).ne.0).and. &
        (frac_coverage(grid%nx-1,grid%ny).ne.0).and. &
        (frac_coverage(grid%nx,grid%ny-1).ne.0)) &
                   frac_coverage(grid%nx,grid%ny)=1.0

    ! Finally fix any rogue points > 1.0

    where (frac_coverage>1.0) frac_coverage=1.0

  end subroutine calc_coverage

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glimmer_load_sigma(model,unit)

    use glide_messages

    !*FD Loads a file containing
    !*FD sigma vertical coordinates.

    ! Arguments

    type(glimmer_global_type),intent(inout) :: model !*FD Ice model to use
    integer,               intent(in)    :: unit  !*FD Logical file unit to use. 
                                                  !*FD The logical file unit specified 
                                                  !*FD must not already be in use

    ! Internal variables

    integer :: up,upn
    logical :: there

    ! Beginning of code

    upn=model%general%upn

    inquire (exist=there,file=model%funits%sigfile)
  
    if (there) then
      open(unit,file=model%funits%sigfile)
      read(unit,'(f5.2)',err=10,end=10) (model%numerics%sigma(up), up=1,upn)
      close(unit)
      return
    end if

10  call glide_msg(GM_FATAL,__FILE__,__LINE__,'something wrong with sigma coord file')
 
   end subroutine glimmer_load_sigma

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glimmer_remove_bath(orog,x,y)

    !*FD Sets ocean areas to zero height, working recursively from
    !*FD a known ocean point.

    real(rk),dimension(:,:),intent(inout) :: orog !*FD Orography --- used for input and output
    integer,                intent(in)    :: x,y  !*FD Location of starting point (index)

    integer :: nx,ny

    nx=size(orog,1) ; ny=size(orog,2)

    if (orog(x,y).lt.0.0) orog(x,y)=0.0
    call glimmer_find_bath(orog,x,y,nx,ny)

  end subroutine glimmer_remove_bath

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  recursive subroutine glimmer_find_bath(orog,x,y,nx,ny)

    !*FD Recursive subroutine called by {\tt glimmer\_remove\_bath}.

    real(rk),dimension(:,:),intent(inout) :: orog  !*FD Orography --- used for input and output
    integer,                intent(in)    :: x,y   !*FD Starting point
    integer,                intent(in)    :: nx,ny !*FD Size of array {\tt orography}

    integer,dimension(4) :: xi=(/ -1,1,0,0 /)
    integer,dimension(4) :: yi=(/ 0,0,-1,1 /)
    integer :: ns=4,i

    do i=1,ns
      if (x+xi(i).le.nx.and.x+xi(i).gt.0.and. &
          y+yi(i).le.ny.and.y+yi(i).gt.0) then
        if (orog(x+xi(i),y+yi(i)).lt.0.0) then
          orog(x+xi(i),y+yi(i))=0.0
          call glimmer_find_bath(orog,x+xi(i),y+yi(i),nx,ny)
        endif
      endif
    enddo

  end subroutine glimmer_find_bath

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(rk) function lon_diff(a,b)

    implicit none

    real(rk),intent(in) :: a,b
    real(rk) :: aa,bb

    aa=a ; bb=b

      do
        if (aa>bb) exit
        aa=aa+360.0
      enddo

    lon_diff=aa-bb

  end function lon_diff

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine calcprcp(which,prcp,up_rate,f1,xwind,ywind,artm,orog,dx,dy,presprcp, &
                      presartm,pfac)

    use glimmer_mbal

    integer,                intent(in)    :: which    !*FD Precipitation option
    real(sp),dimension(:,:),intent(inout) :: prcp     !*FD Large-scale precip on input,
                                                      !*FD desired precip field on output
    real(sp),               intent(in)    :: up_rate  !*FD user-specified uniform precip rate
    real(rk),               intent(in)    :: f1       !*FD scaling parameter
    real(rk),dimension(:,:),intent(in)    :: xwind    !*FD x-wind
    real(rk),dimension(:,:),intent(in)    :: ywind    !*FD y-wind
    real(sp),dimension(:,:),intent(in)    :: artm     !*FD air temperature
    real(rk),dimension(:,:),intent(in)    :: orog     !*FD local orography
    real(rk),               intent(in)    :: dx       !*FD x-spacing
    real(rk),               intent(in)    :: dy       !*FD y-spacing
    real(sp),dimension(:,:),intent(in)    :: presprcp !*FD present-day precip
    real(sp),dimension(:,:),intent(in)    :: presartm !*FD present-day air temp
    real(sp),               intent(in)    :: pfac     !*FD precip parameterisation factor

   select case(which)

    case(0)   ! Uniform precipitation ----------------------------------------

      prcp=up_rate/f1

    case(1)   ! Large scale precip -------------------------------------------
              ! Note that we / by 1000 to convert to meters, and then by
              ! f1 for scaling in the model.

      prcp=prcp/(1000.0*f1)

    case(2)   ! Precip parameterization --------------------------------------

      call glimmer_precip(prcp,xwind,ywind,artm,orog,dx,dy,fixed_a=.true.)
      prcp=prcp/(1000.0*f1)

    case(3)   ! Prescribed small-scale precip from file, ---------------------
              ! adjusted for temperature

      prcp = presprcp * pfac ** (artm - presartm)

    end select

  end subroutine calcprcp

end module glimmer_object
