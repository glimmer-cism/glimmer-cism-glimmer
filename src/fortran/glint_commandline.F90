! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_commandline.f90 - part of the GLIMMER ice model  + 
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

#include "config.inc"

#ifdef HAVE_2003ARGS
#define NARGS   command_argument_count
#define GETARG  get_command_argument
#else
#define NARGS   iargc
#define GETARG  getarg
#endif


module glint_commandline

  use glimmer_global, only:fname_length

  character(len=5000)         :: commandline_history     !< complete command line
  character(len=fname_length) :: commandline_configname  !< name of the configuration file
  character(len=fname_length) :: commandline_climatename !< name of climate configuration file

contains

  !> get the command line and parse it
  !!
  !! \author Magnus Hagdorn
  !! \date April 2009
  subroutine glint_GetCommandline()
    implicit none

    integer numargs
    integer :: i
    integer, external :: iargc
    character(len=100) :: argument
    
    ! get number of arguments and file names
    numargs = NARGS()
    ! reconstruct command line to store commandline_history
    call GETARG(0,commandline_history)
    do i=1,numargs
       call GETARG(i,argument)
       commandline_history = trim(commandline_history)//" "//trim(argument)
    end do

    if (numargs.gt.1) then
       call GETARG(1,commandline_configname)
       call GETARG(2,commandline_climatename)
    else
       write(*,*) 'Enter name of climate configuration file'
       read(*,'(a)') commandline_climatename
       write(*,*) 'Enter name of GLIDE configuration file to be read'
       read(*,'(a)') commandline_configname
    end if
  end subroutine glint_GetCommandline

  !> print out command line
  !!
  !! \author Magnus Hagdorn
  !! \date April 2009
  subroutine glint_PrintCommandline()
    implicit none

    write(*,*) 'Entire commandline'
    write(*,*) trim(commandline_history)
    write(*,*)
    write(*,*) 'commandline_climatename: ',trim(commandline_climatename)
    write(*,*) 'commandline_configname:  ', trim(commandline_configname)
  end subroutine glint_PrintCommandline
end module glint_commandline
