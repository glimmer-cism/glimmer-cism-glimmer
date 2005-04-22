! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_ts.f90 - part of the GLIMMER ice model           + 
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

module glimmer_ts
  !*FD Handling time series

  type glimmer_tseries
     !*FD time series type
     integer :: numt=0                    !*FD number of times in time series
     integer :: numv=1                    !*FD number of values per time
     integer :: current=1                 !*FD current position in ts
     real, dimension(:), pointer :: times=>NULL() !*FD array holding times
     real, dimension(:,:), pointer :: values=>NULL()!*FD array holding values
  end type glimmer_tseries

  interface glimmer_ts_step
     module procedure glimmer_ts_step_array, glimmer_ts_step_scalar
  end interface

  interface glimmer_ts_linear
     module procedure glimmer_ts_linear_array,glimmer_ts_linear_scalar
  end interface

  private :: get_i

contains
  
  subroutine glimmer_read_ts(ts,fname,numv)
    !*FD read time series from file
    use glimmer_log
    implicit none
    type(glimmer_tseries) :: ts           !*FD time series data
    character(len=*), intent(in) :: fname !*FD read from this file
    integer, intent(in),optional :: numv  !*FD number of values per time

    ! local variables
    real dummy
    integer i,j,ios

    if (present(numv)) then
       ts%numv = numv
    else
       ts%numv = 1
    end if

    open(99,file=trim(fname),status='old',iostat=ios)
    
    if (ios.ne.0) then
       call write_log('Error opening file: '//trim(fname),type=GM_FATAL)
    end if

    ! find number of times
    ios=0
    do
       read(99,*,iostat=ios) dummy
       if (ios.ne.0) then
          exit
       end if
       ts%numt = ts%numt + 1
    end do
    rewind(99)

    allocate(ts%times(ts%numt))
    allocate(ts%values(ts%numv,ts%numt))
    ! read data
    do i=1,ts%numt
       read(99,*) ts%times(i),(ts%values(j,i),j=1,ts%numv)
    end do
    close(99)
  end subroutine glimmer_read_ts

  subroutine glimmer_ts_step_array(ts,time,value)
    !*FD interpolate time series by stepping
    use glimmer_log
    implicit none
    type(glimmer_tseries) :: ts     !*FD time series data
    real, intent(in)      :: time   !*FD time value to get
    real, dimension(:)    :: value  !*FD interpolated value
       
    if (size(value).ne.ts%numv) then
       call write_log('Error, wrong number of values',GM_FATAL,__FILE__,__LINE__)
    end if

    value = ts%values(:,get_i(ts,time))
  end subroutine glimmer_ts_step_array

  subroutine glimmer_ts_step_scalar(ts,time,value)
    !*FD interpolate time series by stepping
    use glimmer_log
    implicit none
    type(glimmer_tseries) :: ts     !*FD time series data
    real, intent(in)      :: time   !*FD time value to get
    real                  :: value  !*FD interpolated value
       
    value = ts%values(1,get_i(ts,time))
  end subroutine glimmer_ts_step_scalar
  
  subroutine glimmer_ts_linear_array(ts,time,value)
    !*FD linear interpolate time series
    use glimmer_log
    implicit none
    type(glimmer_tseries) :: ts     !*FD time series data
    real, intent(in)      :: time   !*FD time value to get
    real, dimension(:)    :: value  !*FD interpolated value
       
    integer i
    real,dimension(size(value)) :: slope

    if (size(value).ne.ts%numv) then
       call write_log('Error, wrong number of values',GM_FATAL,__FILE__,__LINE__)
    end if
    
    i = get_i(ts,time)
    slope(:) = (ts%values(:,i+1)-ts%values(:,i))/(ts%times(i+1)-ts%times(i))
    value(:) = ts%values(:,i) + slope(:)*(time-ts%times(i))
  end subroutine glimmer_ts_linear_array

  subroutine glimmer_ts_linear_scalar(ts,time,value)
    !*FD linear interpolate time series
    use glimmer_log
    implicit none
    type(glimmer_tseries) :: ts     !*FD time series data
    real, intent(in)      :: time   !*FD time value to get
    real                  :: value  !*FD interpolated value
       
    integer i
    real :: slope

    i = get_i(ts,time)
    slope = (ts%values(1,i+1)-ts%values(1,i))/(ts%times(i+1)-ts%times(i))
    value = ts%values(1,i) + slope*(time-ts%times(i))
  end subroutine glimmer_ts_linear_scalar
  

  function get_i(ts,time)
    !*FD get index
    implicit none
    type(glimmer_tseries) :: ts     !*FD time series data
    real, intent(in)      :: time   !*FD time value to get
    integer get_i
    integer upper,lower

    ! BC
    if (time.le.ts%times(1)) then
       get_i = 1
       return
    end if
    if (time.ge.ts%times(ts%numt)) then
       get_i = ts%numt
       return
    end if
    ! first try if the interpolated value is around the last value
    ts%current=min(ts%current,ts%numt-1)
    if (time.ge.ts%times(ts%current) .and. time.lt.ts%times(ts%current+1)) then
       get_i = ts%current
       return
    end if
    ! this didn't work, let's try the next interval
    ts%current=ts%current+1
    if (time.ge.ts%times(ts%current) .and. time.lt.ts%times(ts%current+1)) then
       get_i = ts%current
       return
    end if
    ! nope, let's do a Newton search
    lower = 1
    upper = ts%numt
    do
       ts%current = lower+int((upper-lower)/2.)
       if (time.ge.ts%times(ts%current) .and. time.lt.ts%times(ts%current+1)) then
          get_i = ts%current
          return
       end if
       if (time.gt.ts%times(ts%current)) then
          lower = ts%current
       else
          upper = ts%current
       end if
    end do
  end function get_i
end module glimmer_ts
