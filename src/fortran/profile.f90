! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  profile.f90 - part of the GLIMMER ice model              + 
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

module profile
  !*FD Magnus Hagdorn
  !*FD January 2005
  !*FD module for profiling programs

  integer, private :: current_unit = 200
  integer, private,parameter :: max_prof = 100

  type profile_type
     integer :: profile_unit    !*FD file unit to be written to
     real :: start_time         !*FD CPU time at start of log

     real, dimension(max_prof) :: pstart,ptotal  !*FD for each log store start and totals
  end type profile_type

contains
  
  subroutine profile_init(prof,name)
    !*FD initialise a profile
    implicit none
    type(profile_type), intent(out) :: prof !*FD structure storing profile definitions
    character(len=*), intent(in) :: name    !*FD name of file
    ! local variables
    character(len=8)  :: date
    character(len=10) :: time    

    prof%profile_unit = current_unit
    current_unit = current_unit + 1
    call cpu_time(prof%start_time)
    call date_and_time (date, time)
    open(unit=prof%profile_unit,file=name,status='unknown')
    write(unit=prof%profile_unit,fmt="(a,a4,'-',a2,'-',a2,' ',a2,':',a2,':',a6)") '# Started profile on ',&
         date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:10)
  end subroutine profile_init
  
  subroutine profile_start(prof,profn)
    !*FD start profiling
    implicit none
    type(profile_type) :: prof !*FD structure storing profile definitions
    integer, intent(in) :: profn
    
    call cpu_time(prof%pstart(profn))
  end subroutine profile_start

  subroutine profile_stop(prof,profn)
    !*FD stop profiling
    implicit none
    type(profile_type)  :: prof !*FD structure storing profile definitions
    integer, intent(in) :: profn
    
    real t
    call cpu_time(t)
    prof%ptotal(profn) = prof%ptotal(profn) + t-prof%pstart(profn)
  end subroutine profile_stop

  subroutine profile_log(prof,profn,msg)
    !*FD log a message to profile
    implicit none
    type(profile_type)           :: prof !*FD structure storing profile definitions
    integer, intent(in) :: profn
    character(len=*), intent(in) :: msg    !*FD name of file

    real t

    call cpu_time(t)
    write(prof%profile_unit,*) t-prof%start_time,prof%ptotal(profn),profn,msg
    prof%ptotal(profn) = 0.
    prof%pstart(profn) = 0.
  end subroutine profile_log

  subroutine profile_close(prof)
    !*FD close profile
    implicit none
    type(profile_type), intent(in) :: prof !*FD structure storing profile definitions
    ! local variables
    character(len=8)  :: date
    character(len=10) :: time    
    real t

    call cpu_time(t)
    call date_and_time (date, time)
    write(prof%profile_unit,*) '# total elapse cpu time: ',t-prof%start_time
    write(unit=prof%profile_unit,fmt="(a,a4,'-',a2,'-',a2,' ',a2,':',a2,':',a6)") '# Finished profile on ',&
         date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:10)
    close(prof%profile_unit)
  end subroutine profile_close
end module profile
