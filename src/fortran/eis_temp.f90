! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  eis_temp.f90 - part of the GLIMMER ice model             + 
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

module eis_temp
  !*FD climate forcing similar to the old Edinburgh Ice Sheet model
  !*FD Magnus Hagdorn, June 2004
  !*FD temperature forcing

  use glimmer_ts
  use glimmer_global, only : fname_length

  type eis_temp_type
     real :: lapse_rate  = -0.008                !*FD lapse rate
     character(len=fname_length) :: fname=''     !*FD name of file containing temperature ts
     integer :: torder = 2                       !*FD order of temperature polynomial
     type(glimmer_tseries) :: temp_ts            !*FD temperature time series 
     real, dimension(:),pointer :: tvalue        !*FD temperature value
  end type eis_temp_type

contains
   subroutine eis_temp_config(section,temp)
    !*FD get temperature configuration from config file
    use glimmer_config
    implicit none
    type(ConfigSection), pointer :: section !*FD section containing config values
    type(eis_temp_type)           :: temp     !*FD ela data

    call GetValue(section,'temp_file',temp%fname)
    call GetValue(section,'order',temp%torder)
    call GetValue(section,'lapse_rate',temp%lapse_rate)
  end subroutine eis_temp_config

  subroutine eis_temp_printconfig(temp)
    !*FD print configuration to log
    use glimmer_log
    implicit none
    type(eis_temp_type)      :: temp   !*FD temperature data
    ! local variables
    character(len=100) :: message
    call write_log('EIS Temperature')
    call write_log('---------------')
    write(message,*) 'temperature file  : ',trim(temp%fname)
    call write_log(message)
    write(message,*) 'order of temp poly: ',temp%torder
    call write_log(message)
    write(message,*) 'lapse rate        : ',temp%lapse_rate
    call write_log(message)
    call write_log('')
  end subroutine eis_temp_printconfig

  subroutine eis_init_temp(temp,model)
    !*FD initialise temperature forcing
    use glide_types
    use paramets, only: thk0
    implicit none
    type(eis_temp_type)     :: temp  !*FD ela data
    type(glide_global_type) :: model !*FD model instance

    call glimmer_read_ts(temp%temp_ts,temp%fname,temp%torder)
    allocate(temp%tvalue(temp%torder))
    ! scale parameters
    temp%lapse_rate = temp%lapse_rate*thk0
  end subroutine eis_init_temp

  subroutine eis_surftemp(temp,model,time)
    !*FD calculate surface temperature
    use glide_types
    use glimmer_global, only : rk
    implicit none
    type(eis_temp_type)       :: temp  !*FD temperature data
    type(glide_global_type)   :: model !*FD model instance
    real(kind=rk), intent(in) :: time  !*FD current time

    integer i

    call glimmer_ts_step(temp%temp_ts,real(time),temp%tvalue)
    ! spatial temp distrib
    model%climate%artm = temp%tvalue(1)    
    do i=2,temp%torder
       model%climate%artm(:,:) = model%climate%artm(:,:) + temp%tvalue(i)*model%climate%lati(:,:)**(i-1)
    end do
    ! vertical temp distrib
    model%climate%artm(:,:) = model%climate%artm(:,:) + temp%lapse_rate * & 
         (model%geometry%usrf(:,:) - model%climate%eus)
  end subroutine eis_surftemp
end module eis_temp
