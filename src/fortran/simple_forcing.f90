! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  simple_forcing.f90 - part of the GLIMMER ice model       + 
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

module simple_forcing
  !*FD read configuration and generate simple massbalance and 
  !*FD temperature fields

  use glimmer_global, only : sp

  type simple_climate
     ! holds parameters for the simple climate

     real(kind=sp), dimension(2) :: airt = (/ -3.150, -1.e-2 /)  
     !*FD air temperature parameterisation K, K km$^{-3}$
     real(kind=sp), dimension(3) :: nmsb = (/ 0.5, 1.05e-5, 450.0e3 /)
     !*FD mass balance parameterisation m yr$^{-1}$, yr$^{-1}$, m
  end type simple_climate

contains
  subroutine simple_initialise(climate,fname)
    !*FD initialise simple climate model
    use paramets, only: thk0, acc0, scyr
    implicit none
    type(simple_climate) :: climate         !*FD structure holding climate info
    character(len=*), intent(in) :: fname   !*FD name of paramter file

    call simple_readconfig(climate,fname)
    call simple_printconfig(6,climate)

    ! scale parameters
    climate%airt(2) = climate%airt(2) * thk0
    climate%nmsb(1) = climate%nmsb(1) / (acc0 * scyr)
    climate%nmsb(2) = climate%nmsb(2) / (acc0 * scyr)
  end subroutine simple_initialise

  subroutine simple_readconfig(climate,fname)
    !*FD read configuration
    use glimmer_config
    implicit none
    type(simple_climate) :: climate         !*FD structure holding climate info
    character(len=*), intent(in) :: fname   !*FD name of paramter file

    ! local variables
    type(ConfigSection), pointer :: config
    type(ConfigSection), pointer :: section
    real(kind=sp), dimension(:), pointer :: dummy

    ! read configuration
    call ConfigRead(fname,config)

    call GetSection(config,section,'simple')
    if (associated(section)) then
       dummy=>NULL()
       call GetValue(section,'temperature',dummy,2)
       if (associated(dummy)) then
          climate%airt = dummy
          deallocate(dummy)
          dummy=>NULL()
       end if
       call GetValue(section,'massbalance',dummy,3)
       if (associated(dummy)) then
          climate%nmsb = dummy
          deallocate(dummy)
          dummy=>NULL()
       end if
    end if
  end subroutine simple_readconfig

  subroutine simple_printconfig(unit,climate)
    !*FD print simple climate configuration
    implicit none
    integer, intent(in) :: unit       !*FD unit to print to
    type(simple_climate) :: climate   !*FD structure holding climate info
    write(unit,*) '*******************************************************************************'
    write(unit,*) 'Simple Climate configuration'
    write(unit,*) '----------------------------'
    write(unit,*) 'temperature  : ',climate%airt(1)
    write(unit,*) '               ',climate%airt(2)
    write(unit,*) 'massbalance  : ',climate%nmsb(1)
    write(unit,*) '               ',climate%nmsb(2)
    write(unit,*) '               ',climate%nmsb(3)
    write(unit,*) 
  end subroutine simple_printconfig

  subroutine simple_massbalance(climate,model)
    !*FD calculate simple mass balance
    use glide_types
    use paramets, only : len0
    implicit none
    type(simple_climate) :: climate         !*FD structure holding climate info
    type(glide_global_type) :: model        !*FD model instance

    ! local variables
    integer  :: ns,ew
    real :: dist, ewct, nsct, grid

    ! simple EISMINT mass balance
    ewct = real(model%general%ewn+1) / 2.0
    nsct = real(model%general%nsn+1) / 2.0
    grid = model%numerics%dew * len0

    do ns = 1,model%general%nsn
       do ew = 1,model%general%ewn
          dist = grid * sqrt((real(ew) - ewct)**2 + (real(ns) - nsct)**2)
          model%climate%acab(ew,ns) = min(climate%nmsb(1), climate%nmsb(2) * (climate%nmsb(3) - dist))
       end do
    end do
  end subroutine simple_massbalance

  subroutine simple_surftemp(climate,model)
    !*FD calculate simple air surface temperature
    use glide_types
    implicit none
    type(simple_climate) :: climate         !*FD structure holding climate info
    type(glide_global_type) :: model        !*FD model instance

    model%climate%artm(:,:) = climate%airt(1) - model%geometry%thck(:,:) * climate%airt(2)
  end subroutine simple_surftemp
end module simple_forcing
