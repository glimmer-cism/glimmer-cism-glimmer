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

module glide
  !*FD the top-level GLIDE module

  use glide_types

  integer, private, parameter :: dummyunit=99

contains
  subroutine glide_initialise(model,config)
    !*FD initialise GLIDE model instance
    use glide_setup
    use glimmer_ncparams
    use glimmer_ncfile
    use glimmer_ncinfile
    use glide_isot
    use glide_velo
    use glide_thck
    use glide_temp
    use glimmer_log
    use glimmer_config
    use glide_mask
    use glimmer_scales
    use glide_mask
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file
   
    type(ConfigSection), pointer :: ncconfig

    ! initialise scales
    call glimmer_init_scales
    ! read configuration file
    call glide_readconfig(model,config)
    call glide_printconfig(model)
    ! scale parameters
    call glide_scale_params(model)
    ! allocate arrays
    call glide_allocarr(model)
    
    ! load sigma file
    call glide_load_sigma(model,dummyunit)

    ! netCDF I/O
    if (trim(model%funits%ncfile).eq.'') then
       ncconfig => config
    else
       call ConfigRead(model%funits%ncfile,ncconfig)
    end if
    call ReadNCParams(model,ncconfig)
    ! open all input files
    call openall_in(model)
    ! and read first time slice
    call readall(model)

    ! handle relaxed topo
    if (model%options%whichrelaxed.eq.1) then
       model%geometry%relx = model%geometry%topg
    end if

    ! open all output files
    call openall_out(model)

    ! initialise glide components
    call init_isostasy(model)
    call init_velo(model)
    call init_temp(model)
    call init_thck(model)

    if (model%options%hotstart.ne.1) then
       ! initialise Glen's flow parameter A using an isothermal temperature distribution
       call timeevoltemp(model,0)
    end if

    ! calculate mask
    call glide_set_mask(model)

    ! and calculate lower and upper ice surface
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus,model%geometry%lsrf)
    model%geometry%usrf = model%geometry%thck + model%geometry%lsrf

  end subroutine glide_initialise
  
  subroutine glide_tstep(model,time)
    !*FD Performs time-step of an ice model instance.
    use glimmer_global, only : rk
    use glimmer_ncfile
    use glide_isot
    use glide_thck
    use glide_velo
    use glide_setup
    use glide_temp
    use glide_mask
    implicit none

    type(glide_global_type) :: model        !*FD model instance
    real(rk),  intent(in)   :: time         !*FD Current time in years

    ! local variables
    logical :: newtemps

    ! Update internal clock
    model%numerics%time=time  
    ! and possibly write to file
    call writeall(model)
    newtemps = .false.

    ! ------------------------------------------------------------------------ 
    ! Calculate various derivatives...
    ! ------------------------------------------------------------------------     
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

    ! ------------------------------------------------------------------------ 
    ! Do velocity calculation if necessary
    ! ------------------------------------------------------------------------ 
    if (model%numerics%tinc > mod(model%numerics%time,model%numerics%nvel) .or. &
         model%numerics%time == model%numerics%tinc ) then

       call slipvelo(model%numerics, &
            model%velowk,   &
            model%geomderv, &
            (/model%options%whichslip,model%options%whichbtrc/), &
            model%temper%   bwat,     &
            model%velocity% btrc,     &
            model%geometry% relx,     &
            model%velocity% ubas,     &
            model%velocity% vbas)

       call zerovelo(model%velowk,             &
            model%numerics%sigma,     &
            0,                                 &
            model%geomderv% stagthck, &
            model%geomderv% dusrfdew, &
            model%geomderv% dusrfdns, &
            model%temper%   flwa,     &
            model%velocity% ubas,     &
            model%velocity% vbas,     &
            model%velocity% uvel,     &
            model%velocity% vvel,     &
            model%velocity% uflx,     &
            model%velocity% vflx,     &
            model%velocity% diffu)

    end if

    call glide_maskthck(0, &                                    !magi a hack, someone explain what whichthck=5 does
                  model%geometry% thck,      &
                  model%climate%  acab,      &
                  model%geometry% dom,       &
                  model%geometry% mask,      &
                  model%geometry% totpts,    &
                  model%geometry% empty)

    ! ------------------------------------------------------------------------ 
    ! Calculate temperature evolution and Glenn's A, if necessary
    ! ------------------------------------------------------------------------ 
    if ( model%numerics%tinc >  mod(model%numerics%time,model%numerics%ntem)) then
       call timeevoltemp(model, model%options%whichtemp)
       newtemps = .true.
    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate flow evolution by various different methods
    ! ------------------------------------------------------------------------ 
    select case(model%options%whichevol)
    case(0) ! Use precalculated uflx, vflx -----------------------------------

       call timeevolthck(model, &
            0, &                                        !magi a hack, someone explain what whichthck=5 does
            model%geometry% usrf,      &
            model%geometry% thck,      &
            model%climate%  acab,      &
            model%geometry% mask,      &
            model%velocity% uflx,      &
            model%velocity% vflx,      &
            model%geomderv% dusrfdew,  &
            model%geomderv% dusrfdns,  &
            model%geometry% totpts,    &
            6)

    case(1) ! Use explicit leap frog method with uflx,vflx -------------------

       call stagleapthck(model, &
            model%velocity% uflx, &
            model%velocity% vflx, &
            model%geometry% thck)

    case(2) ! Use non-linear calculation that incorporates velocity calc -----

       call nonlevolthck(model, 0, newtemps, 6) !magi a hack, someone explain what whichthck=5 does

    end select

    ! ------------------------------------------------------------------------
    ! get new mask
    ! ------------------------------------------------------------------------
    call glide_set_mask(model)

    ! ------------------------------------------------------------------------ 
    ! Remove ice which is either floating, or is present below prescribed
    ! depth, depending on value of whichmarn
    ! ------------------------------------------------------------------------ 
    call glide_marinlim(model%options%  whichmarn, &
         0, &                                        !magi a hack, someone explain what whichthck=6 does
         model%geometry% thck,      &
         model%geometry% relx,      &
         model%geometry% topg,      &
         model%climate%  lati,      &
         model%geometry%thkmask,    &
         model%numerics%mlimit,     &
         model%climate%eus)

    ! ------------------------------------------------------------------------ 
    ! Calculate isostasy
    ! ------------------------------------------------------------------------ 
    call isosevol(model%numerics,   & 
         model%isotwk,              &                                      
         model%options%  whichisot, &
         model%geometry% thck,      &
         model%geometry% topg,      &
         model%geometry% relx)

    ! ------------------------------------------------------------------------
    ! calculate upper and lower ice surface
    ! ------------------------------------------------------------------------
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0,model%geometry%thck + model%geometry%lsrf)

  end subroutine glide_tstep

  subroutine glide_finalise(model)
    !*FD finalise GLIDE model instance
    use glimmer_ncfile
    use glimmer_ncinfile
    use glimmer_log
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    
    call closeall_in(model)
    call closeall_out(model)
    
    call glide_deallocarr(model)

    call close_log
  end subroutine glide_finalise
end module glide
