! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  eis_forcing.f90 - part of the GLIMMER ice model          + 
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

module eis_forcing
  !*FD climate forcing similar to the old Edinburgh Ice Sheet model
  !*FD Magnus Hagdorn, June 2004
  
  use eis_ela
  use eis_temp

  type eis_climate
     !*FD holds EIS climate
     type(eis_ela_type)  :: ela  !*FD ELA forcing
     type(eis_temp_type) :: temp !*FD temperature forcing
  end type eis_climate
  
contains
  subroutine eis_initialise(climate,fname,model)
    !*FD initialise EIS climate
    use glide_types
    implicit none
    type(eis_climate) :: climate         !*FD structure holding EIS climate
    character(len=*),intent(in) :: fname !*FD name of file containing configuration
    type(glide_global_type) :: model     !*FD model instance

    ! read config
    call eis_readconfig(climate,fname)
    ! print config
    call eis_printconfig(climate)

    ! initialise subsystems
    call eis_init_ela(climate%ela,model)
    call eis_init_temp(climate%temp,model)
  end subroutine eis_initialise

  subroutine eis_readconfig(climate,fname)
    !*FD read EIS configuration
    use glimmer_config
    implicit none
    type(eis_climate) :: climate         !*FD structure holding EIS climate
    character(len=*),intent(in) :: fname !*FD name of file containing configuration
    ! local variables
    type(ConfigSection), pointer :: config
    type(ConfigSection), pointer :: section
 
    ! read configuration
    call ConfigRead(fname,config)

    call GetSection(config,section,'EIS ELA')
    if (associated(section)) then
       call eis_ela_config(section,climate%ela)
    end if
    call GetSection(config,section,'EIS Temperature')
    if (associated(section)) then
       call eis_temp_config(section,climate%temp)
    end if
  end subroutine eis_readconfig

  subroutine eis_printconfig(climate)
    !*FD print EIS configuration
    use glimmer_log
    implicit none
    type(eis_climate) :: climate         !*FD structure holding EIS climate
    
    call write_log_div
    call write_log('Edinburgh Ice Model')
    call eis_ela_printconfig(climate%ela)
    call eis_temp_printconfig(climate%temp)
  end subroutine eis_printconfig
end module eis_forcing
