! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_mbal.f90 - part of the GLIMMER ice model           + 
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

module glint_mbal

  use glimmer_pdd
  use glimmer_enmabal
  use glimmer_global

  implicit none

  !*FD Unified wrapper for different mass-balance codes

  type glint_mbal_params
    type(glimmer_pdd_params),pointer :: annual_pdd => null() !*FD Pointer to annual PDD params
    integer :: which !*FD Flag for chosen mass-balance type
    integer :: tstep !*FD Timestep of mass-balance scheme in hours
  end type glint_mbal_params

contains

  subroutine glint_mbal_init(params,config,which)

    use glimmer_config
    use glimmer_log
    use glint_constants

    !*FD Initialise mass-balance schemes

    type(glint_mbal_params)      :: params !*FD parameters to be initialised
    type(ConfigSection), pointer :: config !*FD structure holding sections of configuration file
    integer,intent(in)           :: which  !*FD selector for pdd type

    ! Copy selector

    params%which=which

    ! Deallocate if necessary
  
    if (associated(params%annual_pdd)) deallocate(params%annual_pdd)

    ! Allocate desired type and initialise
    ! Also check we have a valid value of which

    select case(which)
    case(1)
      allocate(params%annual_pdd)
      call glimmer_pdd_init(params%annual_pdd,config)
      params%tstep=years2hours
    case(2)
      params%tstep=years2hours
    case(3)
      ! The energy-balance model will go here...
      call enmabal_init
      call write_log('Energy-balance mass-balance model not implemented yet',GM_FATAL,__FILE__,__LINE__)
      params%tstep=1
    case default
      call write_log('Invalid value of whichacab',GM_FATAL,__FILE__,__LINE__)
    end select

  end subroutine glint_mbal_init

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_mbal_calc(params,artm,arng,prcp,ablt,acab)

    use glimmer_log

    type(glint_mbal_params)      :: params !*FD parameters to be initialised
    real(sp), dimension(:,:), intent(in)    :: artm    !*FD Annual mean air-temperature 
                                                       !*FD ($^{\circ}$C)
    real(sp), dimension(:,:), intent(in)    :: arng    !*FD Annual temperature half-range ($^{\circ}$C)
    real(sp), dimension(:,:), intent(in)    :: prcp    !*FD Annual accumulated precipitation 
                                                       !*FD (mm water equivalent)
    real(sp), dimension(:,:), intent(out)   :: ablt    !*FD Annual ablation (mm water equivalent)
    real(sp), dimension(:,:), intent(out)   :: acab    !*FD Annual mass-balance (mm water equivalent)

    select case(params%which)
    case(1)
       call glimmer_pdd_mbal(params%annual_pdd,artm,arng,prcp,ablt,acab) 
    case(2) 
       acab = prcp
    case(3)
       ! The energy-balance model will go here...
       call enmabal
       call write_log('Energy-balance mass-balance model not implemented yet',GM_FATAL,__FILE__,__LINE__)
    end select

  end subroutine glint_mbal_calc

end module glint_mbal
