! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_lithot.f90 - part of the GLIMMER ice model         + 
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

! module for temperature calculations in the upper lithosphere

module glide_lithot

contains  
  subroutine init_lithot(model)
    use glide_types
    use glide_setup
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    ! local variables
    integer k

    ! set up vertical grid
    do k=1,model%lithot%nlayer
       model%lithot%deltaz(k) = (1-glide_calc_sigma(real(model%lithot%nlayer-k)/real(model%lithot%nlayer-1),2.)) &
            *model%lithot%rock_base
    end do

    ! set up factors for vertical finite differences
    do k=2,model%lithot%nlayer-1
       model%lithot%zfactors(1,k) = 1./((model%lithot%deltaz(k)-model%lithot%deltaz(k-1)) * &
            (model%lithot%deltaz(k+1)-model%lithot%deltaz(k-1)))
       model%lithot%zfactors(2,k) = 1./((model%lithot%deltaz(k+1)-model%lithot%deltaz(k)) * &
            (model%lithot%deltaz(k)-model%lithot%deltaz(k-1)))
       model%lithot%zfactors(3,k) = 1./((model%lithot%deltaz(k+1)-model%lithot%deltaz(k)) * &
            (model%lithot%deltaz(k+1)-model%lithot%deltaz(k-1)))
    end do

    ! calculate diffusion coefficient
    model%lithot%diffu = model%lithot%con_r/(model%lithot%rho_r*model%lithot%shc_r)

  end subroutine init_lithot


end module glide_lithot
