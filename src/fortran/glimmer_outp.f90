
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_outp.f90 - part of the GLIMMER ice model         + 
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

module glimmer_outp

  use glimmer_types

contains

  subroutine redtsout(fname,unit,forcdata)

    !*FD Sets up forcing data from time-series file (e.g. GRIP data)

    use glimmer_log

    implicit none

    character(*), intent(in)          :: fname     !*FD filename to use
    integer,      intent(in)          :: unit      !*FD File unit to use
    type(glimmer_forcdata),intent(inout) :: forcdata  !*FD Parameters to be set
  
    ! Internal variables

    integer :: count, ios = 0
    logical :: there

    ! ---------------------------------------------------------------
    ! Check to see whether file exists
    ! ---------------------------------------------------------------
 
    inquire(file=fname,exist=there)    

    if ( .not. there ) then
      call glide_msg(GM_FATAL,__FILE__,__LINE__,'Time series file not found')
    endif

    ! ---------------------------------------------------------------
    ! Read in the whole file so we know how many lines there are
    ! ---------------------------------------------------------------

    open(unit,file=fname,form='formatted')

    forcdata%flines = 0

    do while (ios == 0)
      forcdata%flines = forcdata%flines + 1
      read(unit,*,iostat=ios)  
    end do

    forcdata%flines = forcdata%flines - 1

    ! ---------------------------------------------------------------
    ! Allocate array appropriately, then read in data
    ! ---------------------------------------------------------------

    allocate(forcdata%forcing(forcdata%flines,2))

    rewind(unit)

    do count = 1, forcdata%flines
      read(unit,*) forcdata%forcing(count,:)
    end do

    close(unit)

  end subroutine redtsout

end module glimmer_outp
