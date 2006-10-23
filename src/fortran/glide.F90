! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide.f90 - part of the GLIMMER ice model                + 
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

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

module glide
  !*FD the top-level GLIDE module

  use glide_types
  use glide_stop
  use glide_nc_custom
  use glide_io
  use glide_lithot
  use glide_profile
  integer, private, parameter :: dummyunit=99

contains

  subroutine glide_config(model,config)
    !*FD read glide configuration from file and print it to the log
    use glide_setup
    use isostasy
    use glimmer_ncparams
    use glimmer_config
    use glimmer_map_init
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file

    type(ConfigSection), pointer :: ncconfig
   
    ! read configuration file
    call glide_readconfig(model,config)
    call glide_printconfig(model)
    ! Read alternate sigma levels from config file, if necessary
    call glide_read_sigma(model,config)
    ! read isostasy configuration file
    call isos_readconfig(model%isos,config)
    call isos_printconfig(model%isos)
    ! read mapping from config file
    ! **** Use of dew and dns here is an ugly fudge that
    ! **** allows the use of old [GLINT projection] config section
    ! **** for backwards compatibility. It will be deleted soon.
    ! **** (You have been warned!)
    ! **** N.B. Here, dew and dns are unscaled - i.e. real distances in m
    call glimmap_readconfig(model%projection,config, &
         model%numerics%dew, &
         model%numerics%dns)

    ! netCDF I/O
    if (trim(model%funits%ncfile).eq.'') then
       ncconfig => config
    else
       call ConfigRead(model%funits%ncfile,ncconfig)
    end if
    call glimmer_nc_readparams(model,ncconfig)
  end subroutine glide_config

  subroutine glide_initialise(model)
    !*FD initialise GLIDE model instance
    use glide_setup
    use glimmer_ncio
    use glide_velo
    use glide_thck
    use glide_temp
    use glimmer_log
    use glide_mask
    use glimmer_scales
    use glide_mask
    use isostasy
    use glimmer_map_init
    use glide_remap   !lipscomb - for remapping scheme

    implicit none
    type(glide_global_type) :: model        !*FD model instance

    call write_log(glimmer_version)

    ! initialise scales
    call glimmer_init_scales

    ! scale parameters
    call glide_scale_params(model)
    ! set up coordinate systems
    model%general%ice_grid = coordsystem_new(0.d0, 0.d0, &
         model%numerics%dew, model%numerics%dns, &
         model%general%ewn, model%general%nsn)
    model%general%velo_grid = coordsystem_new(model%numerics%dew/2.,model%numerics%dns/2., &
         model%numerics%dew,model%numerics%dns, &
         model%general%ewn-1,model%general%nsn-1)

    ! allocate arrays
    call glide_allocarr(model)

    ! initialise bed softness to uniform parameter
    model%velocity%bed_softness = model%velowk%btrac_const
    
    ! load sigma file
    call glide_load_sigma(model,dummyunit)

    ! open all input files
    call openall_in(model)
    ! and read first time slice
    call glide_io_readall(model,model)
    ! Write projection info to log
    call glimmap_printproj(model%projection)

    ! handle relaxed/equilibrium topo
    ! Initialise isostasy first
    call init_isostasy(model)
    select case(model%options%whichrelaxed)
    case(1) ! Supplied topography is relaxed
       model%isos%relx = model%geometry%topg
    case(2) ! Supplied topography is in equilibrium
       call isos_relaxed(model)
    end select

    ! set uniform basal heat flux
    model%temper%bheatflx = model%paramets%geot

    ! open all output files
    call openall_out(model)
    ! create glide variables
    call glide_io_createall(model)

    ! initialise glide components
    call init_velo(model)
    call init_temp(model)
    call init_thck(model)
    if (model%options%gthf.gt.0) then
       call init_lithot(model)
    end if

!lipscomb - new initialisations for remapping, ice age 
    ! initialise grid-related arrays for remap transport 
    if (model%options%whichevol==3 .or. model%options%whichevol==4) then 
       call init_remap(model%general%ewn,    & 
                       model%general%nsn,    & 
                       model%numerics%dew,   & 
                       model%numerics%dns,   & 
                       model%gridwk) 
    endif 
 
    !initialise ice age 
    model%geometry%age(:,:,:) = 0._dp 

    if (model%options%hotstart.ne.1) then
       ! initialise Glen's flow parameter A using an isothermal temperature distribution
       call timeevoltemp(model,0)
    end if

    ! calculate mask
    call glide_set_mask(model)

    ! and calculate lower and upper ice surface
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus,model%geometry%lsrf)
    model%geometry%usrf = model%geometry%thck + model%geometry%lsrf

    ! initialise profile
#ifdef PROFILING
    call glide_prof_init(model)
#endif
  end subroutine glide_initialise
  
  subroutine glide_tstep_p1(model,time)
    !*FD Performs first part of time-step of an ice model instance.
    !*FD calculate velocity and temperature
    use glimmer_global, only : rk
    use glide_thck
    use glide_velo
    use glide_setup
    use glide_temp
    use glide_mask
    implicit none

    type(glide_global_type) :: model        !*FD model instance
    real(rk),  intent(in)   :: time         !*FD Current time in years

    ! Update internal clock
    model%numerics%time=time  
    model%temper%newtemps = .false.

    ! ------------------------------------------------------------------------ 
    ! Calculate various derivatives...
    ! ------------------------------------------------------------------------     
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%geomderv)
#endif
    call stagvarb(model%geometry% thck, &
         model%geomderv% stagthck,&
         model%general%  ewn, &
         model%general%  nsn)

    call geomders(model%numerics, &
         model%geometry% usrf, &
         model%geomderv% stagthck,&
         model%geomderv% dusrfdew, &
         model%geomderv% dusrfdns)

    call geomders(model%numerics, &
         model%geometry% thck, &
         model%geomderv% stagthck,&
         model%geomderv% dthckdew, &
         model%geomderv% dthckdns)
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%geomderv)
#endif

#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%ice_mask1)
#endif
    call glide_maskthck(0, &                                    !magi a hack, someone explain what whichthck=5 does
         model%geometry% thck,      &
         model%climate%  acab,      &
         model%geometry% dom,       &
         model%geometry% mask,      &
         model%geometry% totpts,    &
         model%geometry% empty)
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%ice_mask1)
#endif

    ! ------------------------------------------------------------------------ 
    ! calculate geothermal heat flux
    ! ------------------------------------------------------------------------ 
    if (model%options%gthf.gt.0) then
       call calc_lithot(model)
    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate temperature evolution and Glenn's A, if necessary
    ! ------------------------------------------------------------------------ 
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%temperature)
#endif
    if ( model%numerics%tinc >  mod(model%numerics%time,model%numerics%ntem)) then
       call timeevoltemp(model, model%options%whichtemp)
       model%temper%newtemps = .true.
    end if
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%temperature)
#endif

    ! ------------------------------------------------------------------------ 
    ! Calculate basal traction factor
    ! ------------------------------------------------------------------------ 
    call calc_btrc(model,model%options%whichbtrc,model%velocity%btrc)

  end subroutine glide_tstep_p1


  subroutine glide_tstep_p2(model)
    !*FD Performs second part of time-step of an ice model instance.
    !*FD write data and move ice
    use glide_thck
    use glide_velo
    use glide_setup
    use glide_temp
    use glide_mask
    use isostasy
    use glide_remap    !lipscomb - for remapping scheme

    implicit none

    type(glide_global_type) :: model        !*FD model instance

    ! ------------------------------------------------------------------------ 
    ! write to netCDF file
    ! ------------------------------------------------------------------------ 
    call glide_io_writeall(model,model)

    ! ------------------------------------------------------------------------ 
    ! Calculate flow evolution by various different methods
    ! ------------------------------------------------------------------------ 
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%ice_evo)
#endif
    select case(model%options%whichevol)
    case(0) ! Use precalculated uflx, vflx -----------------------------------

       call thck_lin_evolve(model,model%temper%newtemps, 6)

    case(1) ! Use explicit leap frog method with uflx,vflx -------------------

       call stagleapthck(model,model%temper%newtemps, 6)

    case(2) ! Use non-linear calculation that incorporates velocity calc -----

       call thck_nonlin_evolve(model,model%temper%newtemps, 6)

!lipscomb - added incremental remapping options 3 and 4
 
    case(3) ! Use incremental remapping scheme for advecting ice thickness ---
            ! (Temperature is advected by glide_temp)
 
       call thck_remap_evolve(model, model%temper%newtemps, 6, .false.)
 
    case(4) ! Use incremental remapping scheme for advecting thickness
            ! and temperature, as well as tracers such as ice age.
 
       call thck_remap_evolve(model, model%temper%newtemps, 6, .true.)

    end select
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%ice_evo)
#endif

    ! ------------------------------------------------------------------------
    ! get new mask
    ! ------------------------------------------------------------------------
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%ice_mask2)
#endif
    call glide_set_mask(model)
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%ice_mask2)
#endif

    ! ------------------------------------------------------------------------ 
    ! Remove ice which is either floating, or is present below prescribed
    ! depth, depending on value of whichmarn
    ! ------------------------------------------------------------------------ 
    call glide_marinlim(model%options%  whichmarn, &
         model%geometry% thck,      &
         model%isos% relx,      &
         model%geometry%topg,   &
         model%geometry%thkmask,    &
         model%numerics%mlimit,     &
         model%numerics%calving_fraction, &
         model%climate%eus,         &
         model%climate%calving)

    ! ------------------------------------------------------------------------
    ! update ice/water load if necessary
    ! ------------------------------------------------------------------------
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%isos_water)
#endif
    if (model%isos%do_isos) then
       if (model%numerics%time.ge.model%isos%next_calc) then
          model%isos%next_calc = model%isos%next_calc + model%isos%period
          call isos_icewaterload(model)
          model%isos%new_load = .true.
       end if
    end if
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%isos_water)
#endif
    
    ! basal shear stress calculations
    call calc_basal_shear(model)
  end subroutine glide_tstep_p2

  subroutine glide_tstep_p3(model)
    !*FD Performs third part of time-step of an ice model instance:
    !*FD calculate isostatic adjustment and upper and lower ice surface
    use isostasy
    use glide_setup
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    
    ! ------------------------------------------------------------------------ 
    ! Calculate isostasy
    ! ------------------------------------------------------------------------ 
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%isos)
#endif
    if (model%isos%do_isos) then
       call isos_isostasy(model)
    end if
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%isos)
#endif

    ! ------------------------------------------------------------------------
    ! calculate upper and lower ice surface
    ! ------------------------------------------------------------------------
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0,model%geometry%thck + model%geometry%lsrf)

    ! increment time counter
    model%numerics%timecounter = model%numerics%timecounter + 1

  end subroutine glide_tstep_p3

!------------------------------------------------------------------------
!lipscomb - added a subroutine to print glide diagnostic output

  subroutine glide_print_diag (model, time, idiag, jdiag)
    !*FD Writes diagnostic output
 
    use glimmer_log
    use paramets, only: thk0, len0, vel0, tim0
    use physcon, only: scyr
 
    implicit none
 
    type(glide_global_type), intent(in) :: model    ! model instance
    real(rk),  intent(in)   :: time                 ! current time in years
    integer, intent(in), optional :: idiag, jdiag
 
    real(dp) ::          &
         tot_area,       &    ! total ice area (km^2)
         tot_volume,     &    ! total ice volume (km^3)
         tot_energy,     &    ! total ice energy (J)
         tot_age,        &    ! sum of volume*age
         max_thck,       &    ! max ice thickness (m)
         min_temp,       &    ! min ice temperature (deg C)
         mean_thck,      &    ! mean ice thickness (m)
         mean_temp,      &    ! mean ice temperature (deg C)
         mean_age,       &    ! mean ice age (yr)
         max_spd_sfc,    &    ! max surface ice speed (m/yr)
         max_spd_bas,    &    ! max basal ice speed (m/yr)
         thck,           &    ! thickness
         spd,            &    ! speed
         age                  ! age
 
    integer :: i, j, k,            &
               imin, jmin, kmin,   &
               imax, jmax, kmax,   &
               ewn, nsn, upn       ! model%numerics%ewn, etc.
 
    character(len=100) :: message
 
    real(dp), parameter ::   &
         eps = 1.0e-11_dp     ! small number
 
 
    ewn = model%general%ewn
    nsn = model%general%nsn
    upn = model%general%upn
 
    !-----------------------------------------------------------------
    ! Compute and write global diagnostics
    !-----------------------------------------------------------------
 
    print*, ' '
    print*, 'Writing diagnostics to log file, time =', time
 
    call write_log(' ')
    write(message,'(a32,f10.2)') 'Global diagnostic output, time =', time
    call write_log(trim(message), type = GM_DIAGNOSTIC)
    call write_log(' ')
 
    ! total ice area
 
    tot_area = 0._dp
    do j = 1, nsn
       do i = 1, ewn
          if (model%geometry%thck(i,j) > model%numerics%thklim) then
             tot_area = tot_area + 1.0_dp
          endif
       enddo
    enddo
    tot_area = tot_area * model%numerics%dew * model%numerics%dns * len0**2
    write(message,'(a25,e20.8)') 'Total ice area (km^2)     ',   &
                                   tot_area*1.0e-6_dp   ! convert to km^2
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! total ice volume
 
    tot_volume = 0._dp
    do j = 1, nsn
       do i = 1, ewn
          if (model%geometry%thck(i,j) > model%numerics%thklim) then
             tot_volume = tot_volume + model%geometry%thck(i,j)
          endif
       enddo
    enddo
    tot_volume = tot_volume * model%numerics%dew * model%numerics%dns  &
               * thk0 * len0**2
    write(message,'(a25,e20.8)') 'Total ice volume (km^3)  ',   &
                                   tot_volume*1.0e-9_dp   ! convert to km^3
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! total ice energy relative to T = 0 deg C
 
    tot_energy = 0._dp
    do j = 1, nsn
       do i = 1, ewn
          if (model%geometry%thck(i,j) > model%numerics%thklim) then
             thck = model%geometry%thck(i,j)
             do k = 1, upn-1
                tot_energy = tot_energy  &
                     + thck * model%velowk%dups(k) * model%temper%temp(k,i,j)
             enddo
          endif
       enddo
    enddo
    tot_energy = tot_energy * model%numerics%dew * model%numerics%dns  &
               * thk0 * len0**2
    write(message,'(a25,e20.8)') 'Total ice energy (J)     ', tot_energy
    call write_log(trim(message), type = GM_DIAGNOSTIC)
    
    ! total sum of volume * age
 
    tot_age = 0._dp
    do j = 1, nsn
       do i = 1, ewn
          if (model%geometry%thck(i,j) > model%numerics%thklim) then
             thck = model%geometry%thck(i,j)
             do k = 1, upn-1
                tot_age = tot_age  &
                     + thck * model%velowk%dups(k) * model%geometry%age(k,i,j)
             enddo
          endif
       enddo
    enddo
    tot_age = tot_age * model%numerics%dew * model%numerics%dns  &
               * thk0 * len0**2 * tim0/scyr
 
    ! mean thickness
 
    if (tot_area > eps) then
       mean_thck = tot_volume/tot_area
    else
       mean_thck = 0._dp
    endif
    write(message,'(a25,f20.12)') 'Mean thickness (m)       ', mean_thck
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! max thickness
 
    imax = 0
    jmax = 0
    max_thck = 0._dp
    do j = 1, nsn
       do i = 1, ewn
          if (model%geometry%thck(i,j) > max_thck) then
             max_thck = model%geometry%thck(i,j)
             imax = i
             jmax = j
          endif
       enddo
    enddo
    write(message,'(a25,f20.12,2i4)') 'Max thickness (m), i, j  ',   &
                                       max_thck*thk0, imax, jmax
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! mean temperature
 
    if (tot_volume > eps) then
       mean_temp = tot_energy/tot_volume
    else
       mean_temp = 0._dp
    endif
    write(message,'(a25,f20.12)') 'Mean temperature (C)     ', mean_temp
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! min temperature
 
    min_temp =  999._dp
    do j = 1, nsn
       do i = 1, ewn
          do k = 1, upn-1
             if (model%temper%temp(k,i,j) < min_temp) then
                min_temp = model%temper%temp(k,i,j)
                imin = i
                jmin = j
                kmin = k
             endif
          enddo
       enddo
    enddo
    write(message,'(a25,f20.12,3i4)') 'Min temperature, i, j, k ',   &
                                       min_temp, imin, jmin, kmin
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! mean age
 
    if (tot_age > eps) then
       mean_age = tot_age/tot_volume
    else
       mean_age = 0._dp
    endif
    write(message,'(a25,f20.12)') 'Mean ice age (yr)        ', mean_age
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! max surface speed
    ! Note: uvel and vvel not defined at i = ewn, j = nsn
    imax = 0
    jmax = 0
    max_spd_sfc = 0._dp
    k = 1
    do j = 1, nsn-1
       do i = 1, ewn-1
          spd = sqrt(model%velocity%uvel(k,i,j)**2   &
                   + model%velocity%vvel(k,i,j)**2)
          if (model%geometry%thck(i,j) > eps .and. spd > max_spd_sfc) then
             max_spd_sfc = spd
             imax = i
             jmax = j
          endif
       enddo
    enddo
    write(message,'(a25,f20.12,2i4)') 'Max sfc spd (m/yr), i, j ',   &
                                       max_spd_sfc*vel0*scyr, imax, jmax
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! max basal speed
 
    imax = 0
    jmax = 0
    max_spd_bas = 0._dp
    do j = 1, nsn-1
       do i = 1, ewn-1
          spd = sqrt(model%velocity%ubas(i,j)**2   &
                   + model%velocity%vbas(i,j)**2)
          if (model%geometry%thck(i,j) > eps .and. spd > max_spd_bas) then
             max_spd_bas = spd
             imax = i
             jmax = j
          endif
       enddo
    enddo
    write(message,'(a25,f20.12,2i4)') 'Max base spd (m/yr), i, j',   &
                                       max_spd_bas*vel0*scyr, imax, jmax
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
 
    ! local diagnostics
 
    if (present(idiag) .and. present(jdiag)) then
       if (idiag >= 1 .and. idiag <=model%general%ewn  &
                      .and.                            &
           jdiag >= 1 .and. jdiag <=model%general%nsn) then
          call write_log(' ')
          write(message,'(a30,2i4)')  &
               'Grid point diagnostics: i, j =', idiag, jdiag
          call write_log(trim(message), type = GM_DIAGNOSTIC)
          call write_log(' ')
 
          i = idiag
          j = jdiag
          write(message,'(a25,f20.12)')   &
               'Thickness (m)            ',model%geometry%thck(i,j)*thk0
          call write_log(trim(message), type = GM_DIAGNOSTIC)
 
          write(message,'(a25,f20.12)')   &
               'Sfc air temperature (C)  ',model%climate%artm(i,j)
          call write_log(trim(message), type = GM_DIAGNOSTIC)
 
          write(message,'(a25,f20.12)')   &
               'Basal temperature (C)    ',model%temper%temp(upn,i,j)
          call write_log(trim(message), type = GM_DIAGNOSTIC)
 
          spd = sqrt(model%velocity%ubas(i,j)**2   &
                   + model%velocity%vbas(i,j)**2) * vel0*scyr
          write(message,'(a25,f20.12)')   &
               'Basal speed (m/yr)       ', spd
          call write_log(trim(message), type = GM_DIAGNOSTIC)
 
          call write_log(' ')
          write(message,'(a52)')   &
               'Layer    Speed (m/yr)  Temperature (C)     Age (yr) '
          call write_log(trim(message), type = GM_DIAGNOSTIC)
 
          do k = 1, upn-1
             spd = sqrt(model%velocity%uvel(k,i,j)**2   &
                      + model%velocity%vvel(k,i,j)**2) * vel0*scyr
             age = model%geometry%age(k,i,j)*tim0/scyr
             write (message,'(i4,f16.6,f16.8,f16.6)')  &
                k, spd, model%temper%temp(k,i,j), age
             call write_log(trim(message), type = GM_DIAGNOSTIC)
          enddo
 
       endif  ! idiag and jdiag in bounds
    endif     ! idiag and jdiag present
 
  end subroutine glide_print_diag


end module glide
