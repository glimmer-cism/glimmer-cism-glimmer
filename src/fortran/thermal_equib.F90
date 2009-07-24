! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  simple_glide.f90 - part of the GLIMMER ice model         + 
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
#include "config.inc"
#endif

program simple_glide
  !*FD This is GLIDE driver that is designed to run a model to thermal equilibrium
  !*FD using higher-order velocity computations.  It differs from a normal run of
  !*FD GLIDE in that the thickness is forced to be held constant and the temperature
  !*FD is hacked into using the higher-order velocities rather than the shallow ice
  !*FD velocities.
  use glimmer_global, only:rk,fname_length
  use glide
  use simple_forcing
  use glimmer_log
  use glimmer_config
  use glide_nonlin
  implicit none

  type(glide_global_type) :: model        ! model instance
  type(simple_climate) :: climate         ! climate
  type(ConfigSection), pointer :: config  ! configuration stuff
  character(len=fname_length) :: fname   ! name of paramter file
  real(kind=rk) time

  real(dp), dimension(:,:,:), pointer :: uvel_sia
  real(dp), dimension(:,:,:), pointer :: vvel_sia
  real(dp), dimension(:,:), pointer :: orig_thck
  real(dp), dimension(:,:,:), pointer :: old_temp
  real(dp) :: err
  integer :: iter

  ! For unstable manifold correction
  real(dp), dimension(:), pointer :: vec_new, vec_old, vec_correction
  integer :: vec_size, lin_start
  logical :: umc_continue

  write(*,*) 'Enter name of GLIDE configuration file to be read'
  read(*,*) fname
  
  ! start logging
  call open_log(unit=50, fname=logname(fname))
  
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

  ! Save the original thickness, as we need to maintain it during the iteration.
  allocate(orig_thck(model%general%ewn, model%general%nsn))
  orig_thck = model%geometry%thck

  call simple_massbalance(climate,model,time)
  call simple_surftemp(climate,model,time)
  call spinup_lithot(model)

  iter = 0

  allocate(old_temp(size(model%temper%temp, 1), size(model%temper%temp, 2), size(model%temper%temp, 3)))
  vec_size = size(model%temper%temp, 1)*size(model%temper%temp, 2)*size(model%temper%temp,3)
  allocate(vec_new(vec_size))
  allocate(vec_old(vec_size))
  allocate(vec_correction(vec_size))
  vec_new = 0
  vec_old = -10
  vec_correction = 0
  old_temp = -10
  model%temper%temp = -10 !Begin the temperature iteration as isothermal

  !Force a non-diagnostic run!
  model%options%diagnostic_run = 0

  write(*,*) shape(old_temp)
  write(*,*) shape(model%temper%temp)
  
  call glide_tstep_p1(model,time)

  do while(iter < 100)
     !Because this is not a "normal" simulation and is instead meant to
     !equilibriate temperatures, I will run thickness first then temperature,
     !and completely ignore isostasy.  This is backwards from what is normally
     !done.
     call glide_tstep_p2(model)
     call reset_thck(model, orig_thck) 
     
     !We want to base temperature calculations off of the higher-order uvel and vvel
     !rather than the SIA uvel and vvel.  I'll hack this in for now by reassigning
     !the pointers in glide_types to trick the module into thinking it's using
     !the SIA guesses...
     uvel_sia => model%velocity%uvel
     vvel_sia => model%velocity%vvel
     model%velocity%uvel => model%velocity_hom%uvel
     model%velocity%vvel => model%velocity_hom%vvel
     
     !We are running this (supposedly) with a frozen bed, so
     !set the basal velocities to 0
     model%velocity%ubas = 0
     model%velocity%vbas = 0

     !Temperature calculation
     call glide_tstep_p1(model,time)
     
     !Put the pointers back the way they were
     model%velocity%uvel => uvel_sia
     model%velocity%vvel => vvel_sia
      
     ! override masking stuff for now
     time = time + model%numerics%tinc
    
     lin_start = 1
     call linearize_3d(vec_new, lin_start, model%temper%temp)
     umc_continue = unstable_manifold_correction(vec_new, vec_old, vec_correction, &
                                                 vec_size, 1d-2, err)

     lin_start = 1
     call delinearize_3d(vec_new, lin_start, model%temper%temp)

     write(*,*) "TEMPERATURE ITERATION: iter = ", iter, ", err = ", err
     if (err < 1e-2) exit
     iter = iter + 1

     call simple_massbalance(climate,model,time)
     call simple_surftemp(climate,model,time)     
  end do

  ! finalise GLIDE
  call glide_finalise(model)
  call close_log
contains
  subroutine reset_thck(model, oldthck)
    type(glide_global_type) :: model
    real(dp), dimension(:,:) :: oldthck

    integer :: i

    model%geometry%thck = oldthck
    do i = 1, model%thckwk%nwhich
        model%thckwk%olds(:,:,2) = oldthck
    end do
  end subroutine
end program simple_glide
