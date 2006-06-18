! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  test_lithot.f90 - part of the GLIMMER ice model          + 
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

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

program test_lithot
  !*FD testing the temperature diffusion in lithosphere layer
  use glimmer_global, only:rk
  use glide
  use glimmer_log
  use glimmer_config
  use simple_forcing
  implicit none

  type(glide_global_type) :: model        ! model instance
  type(simple_climate) :: climate         ! climate
  type(ConfigSection), pointer :: config  ! configuration stuff
  character(len=50) :: fname   ! name of paramter file
  real(kind=rk) time

  integer si,sj
  
  write(*,*) 'Enter name of GLIDE configuration file to be read'
  read(*,*) fname

  ! start logging
  call open_log(unit=50)
  
  ! read configuration
  call ConfigRead(fname,config)

  ! initialise GLIDE
  call glide_config(model,config)
  call simple_initialise(climate,config)
  call glide_initialise(model)
  call CheckSections(config)
  ! fill dimension variables
  call glide_nc_fillall(model)
  time = model%numerics%tstart

  call simple_surftemp(climate,model,time)
  
  ! set all to 0
  model%lithot%temp(:,:,:) = 0.
  model%climate%artm(:,:)  = 0.
  model%geometry%thkmask = -1
  call spinup_lithot(model)

  si = model%general%ewn/4
  sj = model%general%nsn/4

  model%climate%artm(si:3*si,sj:3*sj) = 10.

  do while(time.le.model%numerics%tend)
     model%numerics%time=time  
     call glide_io_writeall(model,model)
     call calc_lithot(model)
     call calc_geoth(model)     
     ! override masking stuff for now
     time = time + model%numerics%tinc
  end do

  ! finalise GLIDE
  call glide_finalise(model)
end program test_lithot
