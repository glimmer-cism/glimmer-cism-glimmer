! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_log.f90 - part of the GLIMMER ice model          + 
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

module glimmer_log
  !*FD module providing file logging

  use glimmer_global, only : fname_length

  character(len=fname_length),private :: glimmer_logname !*FD name of log file
  integer,private :: glimmer_unit=6                      !*FD log unit

contains  
  subroutine open_log(unit,fname)
    !*FD opens log file
    implicit none
    integer, optional          :: unit   !*FD file unit to use
    character(len=*), optional :: fname  !*FD name of log file

    ! local variables
    character(len=8) :: date
    character(len=10) :: time

    if (present(unit)) then
       glimmer_unit = unit
    end if
    if (present(fname)) then
       glimmer_logname = adjustl(trim(fname))
    else
       glimmer_logname = 'glide.log'
    end if

    if (glimmer_unit.ne.6) then
       open(unit=glimmer_unit,file=glimmer_logname,status='unknown')
    end if

    call date_and_time(date,time)
    call write_log_div
    write(unit=glimmer_unit,fmt="(a,a4,'-',a2,'-',a2,' ',a2,':',a2,':',a6)") ' Started logging at ',&
         date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:10)
    call write_log_div
  end subroutine open_log

  subroutine write_log(message,file,line)
    !*FD write to log
    implicit none
    character(len=*),intent(in)          :: message !*FD message to be written
    character(len=*),intent(in),optional :: file    !*FD the name of the file which triggered the message
    integer,intent(in),optional          :: line    !*FD the line number at the which the message was triggered

    if (present(file) .and. present(line)) then
       write(glimmer_unit,*) '(',file,line,') ',message
    else
       write(glimmer_unit,*) trim(message)
    end if
  end subroutine write_log
    
  subroutine error_log(message,file,line)
    !*FD write an error to log and shut down
    implicit none
    character(len=*),intent(in)          :: message !*FD message to be written
    character(len=*),intent(in),optional :: file    !*FD the name of the file which triggered the message
    integer,intent(in),optional          :: line    !*FD the line number at the which the message was triggered

    call write_log(message,file,line)
    call close_log
  end subroutine error_log

  subroutine write_log_div
    !*FD start a new section
    implicit none
    write(glimmer_unit,*) '*******************************************************************************'
  end subroutine write_log_div

  subroutine close_log
    !*FD close log file
    implicit none
    ! local variables
    character(len=8) :: date
    character(len=10) :: time

    call date_and_time(date,time)
    call write_log_div
    write(unit=glimmer_unit,fmt="(a,a4,'-',a2,'-',a2,' ',a2,':',a2,':',a6)") ' Finished logging at ',&
         date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:10)
    call write_log_div
    
    close(glimmer_unit)
  end subroutine close_log

  subroutine sync_log
    !*FD synchronise log to disk
    implicit none
    close(glimmer_unit)
    open(unit=glimmer_unit,file=glimmer_logname, position="append", status='old')
  end subroutine sync_log
end module glimmer_log
