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

private :: linearise

contains  
  subroutine init_lithot(model)
    use glide_types
    use glide_setup
    use paramets, only: len0,tim0
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    ! local variables
    integer i,j,k,r
    real(kind=dp) :: factor

    ! set up vertical grid
    do k=1,model%lithot%nlayer
       model%lithot%deltaz(k) = (1-glide_calc_sigma(real(model%lithot%nlayer-k)/real(model%lithot%nlayer-1),2.)) &
            *model%lithot%rock_base
    end do

    ! calculate diffusion coefficient
    model%lithot%diffu = model%lithot%con_r/(model%lithot%rho_r*model%lithot%shc_r)

    ! set up factors for finite differences
    model%lithot%xfactor = 0.5*model%lithot%diffu*tim0*model%numerics%dt / (model%numerics%dew*len0)
    model%lithot%yfactor = 0.5*model%lithot%diffu*tim0*model%numerics%dt / (model%numerics%dns*len0)
    do k=2,model%lithot%nlayer-1
       model%lithot%zfactors(1,k) =  model%lithot%diffu*tim0*model%numerics%dt / &
            ((model%lithot%deltaz(k)-model%lithot%deltaz(k-1)) * (model%lithot%deltaz(k+1)-model%lithot%deltaz(k-1)))
       model%lithot%zfactors(2,k) = model%lithot%diffu*tim0*model%numerics%dt / &
            ((model%lithot%deltaz(k+1)-model%lithot%deltaz(k)) * (model%lithot%deltaz(k)-model%lithot%deltaz(k-1)))
       model%lithot%zfactors(3,k) = model%lithot%diffu*tim0*model%numerics%dt / &
            ((model%lithot%deltaz(k+1)-model%lithot%deltaz(k)) * (model%lithot%deltaz(k+1)-model%lithot%deltaz(k-1)))
    end do
    k = model%lithot%nlayer
    model%lithot%zfactors(:,k) = 0.5*model%lithot%diffu*tim0*model%numerics%dt / &
         (model%lithot%deltaz(k)-model%lithot%deltaz(k-1))**2

    if (model%options%hotstart.ne.1) then
       ! set initial temp distribution to thermal gradient
       factor = model%paramets%geot/model%lithot%con_r
       do k=1,model%lithot%nlayer
          model%lithot%temp(:,:,k) = model%climate%artm(:,:)+model%lithot%deltaz(k)*factor
       end do
    end if

    ! calculate finite difference coefficient matrix
    ! interior
    do k=2, model%lithot%nlayer-1
       do j=2,model%general%nsn-1
          do i=2,model%general%ewn-1
             r = linearise(model,i,j,k)
             ! i+1,j,k
             call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i+1,j,k), -model%lithot%xfactor)
             ! i-1,j,k
             call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i-1,j,k), -model%lithot%xfactor)
             ! i,j+1,k
             call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j+1,k), -model%lithot%yfactor)
             ! i,j-1,k
             call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j-1,k), -model%lithot%yfactor)
             ! i,j,k+1
             call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k+1), -model%lithot%zfactors(3,k))
             ! i,j,k-1
             call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k-1), -model%lithot%zfactors(1,k))
             ! i,j,k
             call sparse_insert_val(model%lithot%fd_coeff,r,r, &
                  2.*model%lithot%xfactor + 2.*model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)
          end do
       end do
    end do
    ! boundaries
    ! left and right face
    do k=2, model%lithot%nlayer-1
       do j=2,model%general%nsn-1
          ! left face
          i = 1
          r = linearise(model,i,j,k)
          ! i+1
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i+1,j,k), -model%lithot%xfactor)
          ! i-1
          ! 0
          ! i,j+1,k
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j+1,k), -model%lithot%yfactor)
          ! i,j-1,k
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j-1,k), -model%lithot%yfactor)
          ! i,j,k+1
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k+1), -model%lithot%zfactors(3,k))
          ! i,j,k-1
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k-1), -model%lithot%zfactors(1,k))
          ! i,j,k
          call sparse_insert_val(model%lithot%fd_coeff,r,r, &
               model%lithot%xfactor + 2.*model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)
          
          ! right face
          i = model%general%ewn
          r = linearise(model,i,j,k)
          ! i+1
          ! 0
          ! i-1,j,k
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i-1,j,k), -model%lithot%xfactor)
          ! i,j+1,k
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j+1,k), -model%lithot%yfactor)
          ! i,j-1,k
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j-1,k), -model%lithot%yfactor)
          ! i,j,k+1
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k+1), -model%lithot%zfactors(3,k))
          ! i,j,k-1
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k-1), -model%lithot%zfactors(1,k))
          ! i,j,k
          call sparse_insert_val(model%lithot%fd_coeff,r,r, &
               model%lithot%xfactor + 2.*model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)
       end do
    end do
    ! front and back face
    do k=2, model%lithot%nlayer-1
       do i=2,model%general%ewn-1
          ! front face
          j = 1
          r = linearise(model,i,j,k)
          ! i+1,j,k
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i+1,j,k), -model%lithot%xfactor)
          ! i-1,j,k
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i-1,j,k), -model%lithot%xfactor)
          ! i,j+1,k
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j+1,k), -model%lithot%yfactor)
          ! i,j-1,k
          ! 0          
          ! i,j,k+1
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k+1), -model%lithot%zfactors(3,k))
          ! i,j,k-1
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k-1), -model%lithot%zfactors(1,k))
          ! i,j,k
          call sparse_insert_val(model%lithot%fd_coeff,r,r, &
               2.*model%lithot%xfactor + model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)

          ! back face
          j=model%general%nsn
          r = linearise(model,i,j,k)
          ! i+1,j,k
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i+1,j,k), -model%lithot%xfactor)
          ! i-1,j,k
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i-1,j,k), -model%lithot%xfactor)
          ! i,j+1,k
          ! 0
          ! i,j-1,k
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j-1,k), -model%lithot%yfactor)
          ! i,j,k+1
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k+1), -model%lithot%zfactors(3,k))
          ! i,j,k-1
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k-1), -model%lithot%zfactors(1,k))
          ! i,j,k
          call sparse_insert_val(model%lithot%fd_coeff,r,r, &
               2.*model%lithot%xfactor + model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)
       end do
    end do
    ! vertical edges
    do k=2, model%lithot%nlayer-1
       i = 1
       j = 1
       r = linearise(model,i,j,k)
       ! i+1,j,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i+1,j,k), -model%lithot%xfactor)
       ! i,j+1,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j+1,k), -model%lithot%yfactor)
       ! i,j,k+1
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k+1), -model%lithot%zfactors(3,k))
       ! i,j,k-1
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k-1), -model%lithot%zfactors(1,k))
       ! i,j,k
       call sparse_insert_val(model%lithot%fd_coeff,r,r, &
            model%lithot%xfactor + model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)

       i = 1
       j = model%general%nsn
       r = linearise(model,i,j,k)
       ! i+1,j,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i+1,j,k), -model%lithot%xfactor)
       ! i,j-1,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j-1,k), -model%lithot%yfactor)
       ! i,j,k+1
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k+1), -model%lithot%zfactors(3,k))
       ! i,j,k-1
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k-1), -model%lithot%zfactors(1,k))
       ! i,j,k
       call sparse_insert_val(model%lithot%fd_coeff,r,r, &
            model%lithot%xfactor + model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)

       i = model%general%ewn
       j = model%general%nsn
       r = linearise(model,i,j,k)
       ! i-1,j,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i-1,j,k), -model%lithot%xfactor)
       ! i,j-1,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j-1,k), -model%lithot%yfactor)
       ! i,j,k+1
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k+1), -model%lithot%zfactors(3,k))
       ! i,j,k-1
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k-1), -model%lithot%zfactors(1,k))
       ! i,j,k
       call sparse_insert_val(model%lithot%fd_coeff,r,r, &
            model%lithot%xfactor + model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)

       i = model%general%ewn
       j = 1
       r = linearise(model,i,j,k)
       ! i-1,j,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i-1,j,k), -model%lithot%xfactor)
       ! i,j+1,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j+1,k), -model%lithot%yfactor)
       ! i,j,k+1
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k+1), -model%lithot%zfactors(3,k))
       ! i,j,k-1
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k-1), -model%lithot%zfactors(1,k))
       ! i,j,k
       call sparse_insert_val(model%lithot%fd_coeff,r,r, &
            model%lithot%xfactor + model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)
    end do
    ! bottom face
    k = model%lithot%nlayer
    do j=2,model%general%nsn-1
       do i=2,model%general%ewn-1
          r = linearise(model,i,j,k)
          ! i+1,j,k
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i+1,j,k), -model%lithot%xfactor)
          ! i-1,j,k
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i-1,j,k), -model%lithot%xfactor)
          ! i,j+1,k
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j+1,k), -model%lithot%yfactor)
          ! i,j-1,k
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j-1,k), -model%lithot%yfactor)
          ! i,j,k-1
          call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j-1,k-1), -model%lithot%zfactors(1,k))
          ! i,j,k
          call sparse_insert_val(model%lithot%fd_coeff,r,r, &
                  2.*model%lithot%xfactor + 2.*model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)
       end do
    end do
    ! bottom edges
    do j=2,model%general%nsn-1
       ! left edge
       i = 1
       r = linearise(model,i,j,k)
       ! i+1
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i+1,j,k), -model%lithot%xfactor)
       ! i-1
       ! 0
       ! i,j+1,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j+1,k), -model%lithot%yfactor)
       ! i,j-1,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j-1,k), -model%lithot%yfactor)
       ! i,j,k-1
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j-1,k-1), -model%lithot%zfactors(1,k))
       ! i,j,k
       call sparse_insert_val(model%lithot%fd_coeff,r,r, &
            model%lithot%xfactor + 2.*model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)

       ! right edge
       i = model%general%ewn
       r = linearise(model,i,j,k)
       ! i+1
       ! 0
       ! i-1,j,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i-1,j,k), -model%lithot%xfactor)
       ! i,j+1,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j+1,k), -model%lithot%yfactor)
       ! i,j-1,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j-1,k), -model%lithot%yfactor)
       ! i,j,k-1
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k-1), -model%lithot%zfactors(1,k))
       ! i,j,k
       call sparse_insert_val(model%lithot%fd_coeff,r,r, &
            model%lithot%xfactor + 2.*model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)
    end do
    do i=2,model%general%ewn-1
       ! front edge
       j = 1
       r = linearise(model,i,j,k)
       ! i+1,j,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i+1,j,k), -model%lithot%xfactor)
       ! i-1,j,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i-1,j,k), -model%lithot%xfactor)
       ! i,j+1,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j+1,k), -model%lithot%yfactor)
       ! i,j-1,k
       ! 0          
       ! i,j,k-1
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k-1), -model%lithot%zfactors(1,k))
       ! i,j,k
       call sparse_insert_val(model%lithot%fd_coeff,r,r, &
            2.*model%lithot%xfactor + model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)

       ! back edge
       j=model%general%nsn
       r = linearise(model,i,j,k)
       ! i+1,j,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i+1,j,k), -model%lithot%xfactor)
       ! i-1,j,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i-1,j,k), -model%lithot%xfactor)
       ! i,j+1,k
       ! 0
       ! i,j-1,k
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j-1,k), -model%lithot%yfactor)
       ! i,j,k-1
       call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k-1), -model%lithot%zfactors(1,k))
       ! i,j,k
       call sparse_insert_val(model%lithot%fd_coeff,r,r, &
            2.*model%lithot%xfactor + model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)
    end do
    ! bottom corners
    i = 1
    j = 1
    r = linearise(model,i,j,k)
    ! i+1,j,k
    call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i+1,j,k), -model%lithot%xfactor)
    ! i,j+1,k
    call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j+1,k), -model%lithot%yfactor)
    ! i,j,k-1
    call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k-1), -model%lithot%zfactors(1,k))
    ! i,j,k
    call sparse_insert_val(model%lithot%fd_coeff,r,r, &
         model%lithot%xfactor + model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)

    i = 1
    j = model%general%nsn
    r = linearise(model,i,j,k)
    ! i+1,j,k
    call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i+1,j,k), -model%lithot%xfactor)
    ! i,j-1,k
    call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j-1,k), -model%lithot%yfactor)
    ! i,j,k-1
    call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k-1), -model%lithot%zfactors(1,k))
    ! i,j,k
    call sparse_insert_val(model%lithot%fd_coeff,r,r, &
         model%lithot%xfactor + model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)

    i = model%general%ewn
    j = model%general%nsn
    r = linearise(model,i,j,k)
    ! i-1,j,k
    call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i-1,j,k), -model%lithot%xfactor)
    ! i,j-1,k
    call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j-1,k), -model%lithot%yfactor)
    ! i,j,k-1
    call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k-1), -model%lithot%zfactors(1,k))
    ! i,j,k
    call sparse_insert_val(model%lithot%fd_coeff,r,r, &
         model%lithot%xfactor + model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)

    i = model%general%ewn
    j = 1
    r = linearise(model,i,j,k)
    ! i-1,j,k
    call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i-1,j,k), -model%lithot%xfactor)
    ! i,j+1,k
    call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j+1,k), -model%lithot%yfactor)
    ! i,j,k-1
    call sparse_insert_val(model%lithot%fd_coeff,r,linearise(model,i,j,k-1), -model%lithot%zfactors(1,k))
    ! i,j,k
    call sparse_insert_val(model%lithot%fd_coeff,r,r, &
         model%lithot%xfactor + model%lithot%yfactor + model%lithot%zfactors(2,k) + 1.)
    
    ! top face
    ! simply match air temperature where no ice and basal temperature where ice
    k = 1
    do j=1,model%general%nsn
       do i=1,model%general%ewn
          r = linearise(model,i,j,k)
          call sparse_insert_val(model%lithot%fd_coeff,r,r, 1.d0)
       end do
    end do

    ! convert from SLAP Triad to SLAP Column format
    call copy_sparse_matrix(model%lithot%fd_coeff,model%lithot%fd_coeff_slap)
    call ds2y(model%general%nsn*model%general%ewn*model%lithot%nlayer,model%lithot%fd_coeff_slap%n, &
         model%lithot%fd_coeff_slap%col,model%lithot%fd_coeff_slap%row,model%lithot%fd_coeff_slap%val, 0)

    ! initialise result vector
    do k=1,model%lithot%nlayer
       do j=1,model%general%nsn
          do i=1,model%general%ewn
             model%lithot%answer(linearise(model,i,j,k)) = model%lithot%temp(i,j,k)
          end do
       end do
    end do
  end subroutine init_lithot

  subroutine spinup_lithot(model)
    use glide_types
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    integer i,k,j
    real(dp) :: factor

    if (model%options%hotstart.ne.1) then
       ! set initial temp distribution to thermal gradient
       factor = model%paramets%geot/model%lithot%con_r
       do k=1,model%lithot%nlayer
          model%lithot%temp(:,:,k) = model%climate%artm(:,:)+model%lithot%deltaz(k)*factor
       end do
    end if

    ! initialise result vector
    do k=1,model%lithot%nlayer
       do j=1,model%general%nsn
          do i=1,model%general%ewn
             model%lithot%answer(linearise(model,i,j,k)) = model%lithot%temp(i,j,k)
          end do
       end do
    end do

    call calc_geoth(model)
  end subroutine spinup_lithot

  subroutine calc_lithot(model)
    use glide_types
    use glide_mask
    use glide_stop
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    integer i,j,k,r
    real(kind=dp) :: factor
    integer iter
    real(dp) err
    real(dp), parameter :: tol = 1.0d-12
    integer, parameter :: isym = 0, itol = 2, itmax = 101
    integer :: ierr

    ! calculate RHS
    call sparse_matrix_vec_prod(model%lithot%fd_coeff,model%lithot%answer,model%lithot%rhs)
    model%lithot%rhs = -model%lithot%rhs + 2. * model%lithot%answer
    k = model%lithot%nlayer
    factor = model%paramets%geot/model%lithot%con_r/(model%lithot%deltaz(k)-model%lithot%deltaz(k-1))
    do j=1,model%general%nsn
       do i=1,model%general%ewn
          k = 1
          r = linearise(model,i,j,k)
          if (is_ocean(model%geometry%thkmask(i,j))) then
             model%lithot%rhs(r) = 2.    ! 2degC for bottom of ocean
          else if (is_land(model%geometry%thkmask(i,j))) then
             model%lithot%rhs(r) = model%climate%artm(i,j) ! air temperature outside ice sheet
          else
             model%lithot%rhs(r) = model%temper%temp(model%general%upn,i,j) ! ice basal temperature
          end if

          k = model%lithot%nlayer
          r = linearise(model,i,j,k)
          model%lithot%rhs(r) = model%lithot%rhs(r) + 2.*(model%lithot%xfactor + model%lithot%yfactor)*factor
       end do
    end do

    ! solve matrix equation
    call dslucs(model%general%nsn*model%general%ewn*model%lithot%nlayer, model%lithot%rhs, model%lithot%answer, &
         model%lithot%fd_coeff_slap%n, model%lithot%fd_coeff_slap%col,model%lithot%fd_coeff_slap%row, &
         model%lithot%fd_coeff_slap%val, isym,itol,tol,itmax,iter,err,ierr,0, &
         model%lithot%rwork, model%lithot%mxnelt, model%lithot%iwork, model%lithot%mxnelt)

    if (ierr /= 0) then
      print *, 'pcg error ', ierr, itmax, iter
      write(*,*) model%numerics%time
      call glide_finalise(model,.true.)
      stop
    end if

    ! de-linearise results
    do k=1, model%lithot%nlayer
       do j=1,model%general%nsn
          do i=1,model%general%ewn
             model%lithot%temp(i,j,k) = model%lithot%answer(linearise(model,i,j,k))
          end do
       end do
    end do
      
    call calc_geoth(model)

  end subroutine calc_lithot

  subroutine calc_geoth(model)
    !*FD calculate geothermal heat flux
    use glide_types
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    real(dp) factor

    factor = model%lithot%con_r/(model%lithot%deltaz(2)-model%lithot%deltaz(1))
    model%temper%bheatflx(:,:) = factor*(model%lithot%temp(:,:,2)-model%lithot%temp(:,:,1))
  end subroutine calc_geoth

  function linearise(model,i,j,k)
    use glide_types
    implicit none
    type(glide_global_type),intent(in) :: model   
    integer, intent(in) :: i,j,k
    integer :: linearise
    
    linearise = i + (j-1)*model%general%ewn + (k-1)*model%general%ewn*model%general%nsn
  end function linearise


end module glide_lithot
