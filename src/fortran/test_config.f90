! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  testconfig.f90 - part of the GLIMMER ice model           + 
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

program testconfig
  !*FD testing config module
  !*FD written by Magnus Hagdorn, May 2004
  use glimmer_config
  implicit none

  character(len=30) fname
  type(ConfigSection), pointer :: config,section

  character(len=100) :: charval
  integer :: intval
  real :: realval
  real, dimension(:), pointer :: realarray

  write(*,*) 'Enter name of configuration file'
  read(*,*) fname

  call ConfigRead(fname,config)

  call PrintConfig(config)

  call GetSection(config,section,'a section')

  write(*,*) associated(section)
  if (.not.associated(section)) then
     write(*,*) 'Huh?1'
     stop
  end if

  call GetValue(section,'an_int',intval)
  call GetValue(section,'a_float',realval)
  call GetValue(section,'a_char',charval)
  call GetValue(section,'an_array',realarray)
  write(*,*) intval,realval,trim(charval)
  write(*,*) realarray

  intval = -1
  realval = -10.
  call GetValue(section,'an_int',realval)
  call GetValue(section,'a_float',charval)
  call GetValue(section,'a_char',intval)
  write(*,*) intval,realval,trim(charval)
end program testconfig
