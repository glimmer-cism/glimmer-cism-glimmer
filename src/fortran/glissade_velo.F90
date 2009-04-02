! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
! +                                                           + 
! +  glissade_velo.F90                                             + 
! +                                                           + 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!
! This module contains routines for computing the ice velocity
! using the variational approach of Dukowicz, Price and Lipscomb
! (submitted to J. Glac., 2009).
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

  module glissade_velo

    use glimmer_global, only: sp, dp
!    use glide_deriv
!    use glide_grids, only: stagvarb

    implicit none
    private

!whl - to do - Move constants to another module

! constants 

    integer, parameter :: &
       kne = 1,    &       ! indices for the four corners of a grid cell
       knw = 2,    &
       kse = 3,    &
       ksw = 4

    real(dp), parameter ::   &
       c0 = 0.0_dp,     &
       c1 = 1.0_dp,     &
       c2 = 2.0_dp,     &
       c3 = 3.0_dp,     &
       c4 = 4.0_dp,     &
       p25  = c1/c4,    &
       p333 = c1/c3,    &
       p5   = c1/c2,    &
       p667 = c2/c3

    public :: glissade_velo_init, glissade_velo_driver

    contains

!-----------------------------------------------------------------------------
!****************************************************************************

!whl - The next two subroutines are just stubs for now

  subroutine glissade_velo_init(model)

    use glide_types

    type(glide_global_type), intent(inout) :: model


  end subroutine glissade_velo_init

!****************************************************************************

  subroutine glissade_velo_driver(model)

    use glide_types

    type(glide_global_type), intent(inout) :: model

  end subroutine glissade_velo_driver

!****************************************************************************

  subroutine calc_eps2_2d (ewn,       nsn,   &
                           dew,       dns,   &
                           uvel,      vvel,  &
                           thckmask,  eps2)
  !---------------------------------------------------------------------------
  ! Calculate the effective strain rate for the 2D shallow shelf approximation
  !
  ! eps^2 = (du/dx)^2 + (dv/dy)^2 + (du/dx + dv/dy)^2 + (1/2)*(du/dy + dv/dx)^2
  !
  !       =       2 * [(du/dx)^2 + (dv/dy)^2 + (du/dx*dv/dy)] + du/dy*dv/dx
  !         + (1/2) * [(dv/dx)^2 + (du/dy)^2]
  !
  ! Each component in this sum is an integrated average over the grid cell,
  ! given a bilinear reconstruction of u and v as a function of the four
  ! nodal values.
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! in-out variables
  !---------------------------------------------------------------------------

  integer, intent(in) ::  &
      ewn, nsn                  ! horizontal grid dimensions
      
  real(dp), intent(in) ::   &
      dew, dns                  ! E-W and N-S grid lengths

  real(dp), dimension(:,:), intent(in) ::    &
      uvel, vvel                ! E-W and N-S velocity components
 
  integer, dimension(:,:), intent(in) ::     &
      thckmask                  ! = 1 where ice is present, else = 0

  real(dp), dimension(:,:), intent(out) ::   &
       eps2                     ! square of effective strain rate for SSA

  !---------------------------------------------------------------------------
  ! local variables
  !---------------------------------------------------------------------------

  integer ::  i, j           ! horizontal indices

  real(dp) ::    &
     fxx, fyy, fxy,         &! geometric factors
     dun, dus, due, duw,    &! difference in uvel between adjacent nodes
     dvn, dvs, dve, dvw,    &! difference in vvel between adjacent nodes
     dudx2, dvdx2,          &! components of eps2
     dudy2, dvdy2,          &
     dudxdvdy, dvdxdudy

  !---------------------------------------------------------------------------
  ! calculate eps2
  !---------------------------------------------------------------------------

!whl - may be able to speed this up by setting duw(i,j) = due(i-1,j)?

  fxx = c1 / (dew*dew)
  fyy = c1 / (dns*dns)
  fxy = c1 / (dew*dns) 

  do j = 1, nsn
  do i = 1, ewn

     if (thckmask(i,j)==1) then

        dun = uvel(i,j)   - uvel(i-1,j)
        dus = uvel(i,j-1) - uvel(i-1,j-1)
        due = uvel(i,j)   - uvel(i,j-1)
        duw = uvel(i-1,j) - uvel(i-1,j-1)

        dvn = vvel(i,j)   - vvel(i-1,j)
        dvs = vvel(i,j-1) - vvel(i-1,j-1)
        dve = vvel(i,j)   - vvel(i,j-1)
        dvw = vvel(i-1,j) - vvel(i-1,j-1)

        dudx2 = p333 * fxx *  ( dun*dun + dus*dus + dun*dus)
        dvdx2 = p333 * fxx *  ( dvn*dvn + dvs*dvs + dvn*dvs)

        dudy2 = p333 * fyy *  ( due*due + duw*duw + due*duw)
        dvdy2 = p333 * fyy *  ( dve*dve + dvw*dvw + dve*dvw)
  
        dudxdvdy = p25 * fxy *  (due*dvn + due*dvs + duw*dvn + duw*dvs)
        dvdxdudy = p25 * fxy *  (dve*dun + dve*dus + dvw*dun + dvw*dus)

        eps2(i,j) = c2 * (dudx2 + dvdy2 + dudxdvdy)  +  dvdxdudy   &
                  + p5 * (dvdx2 + dudy2)

     else    ! thckmask = 0

        eps2(i,j) = c0

     endif   ! thckmask

  enddo   ! i
  enddo   ! j

  end subroutine calc_eps2_2d

!****************************************************************************

  subroutine calc_Du_2d (ewn,           nsn,          &
                         dew,           dns,          &  
                         uvel,          vvel,         &
                         thckmask,                    &
                         Du_ne,         Du_nw,        &
                         Du_se,         Du_sw,        &
                         Dv_ne,         Dv_nw,        &
                         Dv_se,         Dv_sw )

  !---------------------------------------------------------------------------
  ! Compute the derivative of the cell-integrated average of eps^2 with respect
  ! to each velocity component (u,v) at each cell corner (ne, nw, se, sw).
  !
  ! Do this by computing the derivatives for each of 6 terms in eps^2,
  ! then summing over the terms.
  !
  ! Recall that for the 2D SSA:
  ! eps^2  =      2 * [(du/dx)^2 + (dv/dy)^2 + (du/dx*dv/dy)] + du/dy*dv/dx
  !         + (1/2) * [(dv/dx)^2 + (du/dy)^2]
  !
  !---------------------------------------------------------------------------


  !---------------------------------------------------------------------------
  ! in-out variables
  !---------------------------------------------------------------------------

  integer, intent(in) ::  &
      ewn, nsn                  ! horizontal grid dimensions
      
  real(dp), intent(in) ::   &
      dew, dns                  ! E-W and N-S grid lengths

  real(dp), dimension(:,:), intent(in) ::      &
      uvel, vvel                ! E-W and N-S velocity components
 
  integer, dimension(:,:), intent(in) ::       &
      thckmask                  ! = 1 where ice is present, else = 0

  real(dp), dimension(:,:), intent(out) ::   &
      Du_ne,       &! derivative of eps2 wrt u_ne
      Du_nw,       &! derivative of eps2 wrt u_nw
      Du_se,       &! derivative of eps2 wrt u_se
      Du_sw,       &! derivative of eps2 wrt u_sw
      Dv_ne,       &! derivative of eps2 wrt v_ne
      Dv_nw,       &! derivative of eps2 wrt v_nw
      Dv_se,       &! derivative of eps2 wrt v_se
      Dv_sw         ! derivative of eps2 wrt v_sw
  
  !---------------------------------------------------------------------------
  ! local variables
  !---------------------------------------------------------------------------

  integer ::  i, j           ! horizontal indices

  real(dp) ::    &
     fxx, fyy, fxy,         &! geometric factors
     dun, dus, due, duw,    &! difference in uvel between adjacent nodes
     dvn, dvs, dve, dvw      ! difference in vvel between adjacent nodes

  real(dp) ::   &
      Du_dudx2_ne,     Du_dudx2_nw,     &!  d/du [(du/dx)^2]
      Du_dudx2_se,     Du_dudx2_sw,     &!
      Dv_dvdx2_ne,     Dv_dvdx2_nw,     &!  d/dv [(dv/dx)^2]
      Dv_dvdx2_se,     Dv_dvdx2_sw,     &!
      Du_dudy2_ne,     Du_dudy2_nw,     &!  d/du [(du/dy)^2]
      Du_dudy2_se,     Du_dudy2_sw,     &!
      Dv_dvdy2_ne,     Dv_dvdy2_nw,     &!  d/dv [(dv/dy)^2]
      Dv_dvdy2_se,     Dv_dvdy2_sw,     &!
      Du_dudxdvdy_ne,  Du_dudxdvdy_nw,  &!  d/du [(du/dx)*(dv/dy)]
      Du_dudxdvdy_se,  Du_dudxdvdy_sw,  &!  
      Dv_dudxdvdy_ne,  Dv_dudxdvdy_nw,  &!  d/dv [(du/dx)*(dv/dy)]
      Dv_dudxdvdy_se,  Dv_dudxdvdy_sw,  &!  
      Du_dvdxdudy_ne,  Du_dvdxdudy_nw,  &!  d/du [(dv/dx)*(du/dy)]
      Du_dvdxdudy_se,  Du_dvdxdudy_sw,  &!  
      Dv_dvdxdudy_ne,  Dv_dvdxdudy_nw,  &!  d/dv [(dv/dx)*(du/dy)]
      Dv_dvdxdudy_se,  Dv_dvdxdudy_sw

  !---------------------------------------------------------------------------
  ! For each grid cell, calculate various derivatives of the effective strain rate
  ! with respect to the nodal velocities.
  !---------------------------------------------------------------------------

  fxx = c1 / (dew*dew)
  fyy = c1 / (dns*dns)
  fxy = c1 / (dew*dns) 

  Du_ne(:,:) = c0
  Du_nw(:,:) = c0
  Du_se(:,:) = c0
  Du_sw(:,:) = c0

  Dv_ne(:,:) = c0
  Dv_nw(:,:) = c0
  Dv_se(:,:) = c0
  Dv_sw(:,:) = c0
 
  do j = 1, nsn
  do i = 1, ewn

     if (thckmask(i,j)==1) then

        dun = uvel(i,j)   - uvel(i-1,j)
        dus = uvel(i,j-1) - uvel(i-1,j-1)
        due = uvel(i,j)   - uvel(i,j-1)
        duw = uvel(i-1,j) - uvel(i-1,j-1)

        dvn = vvel(i,j)   - vvel(i-1,j)
        dvs = vvel(i,j-1) - vvel(i-1,j-1)
        dve = vvel(i,j)   - vvel(i,j-1)
        dvw = vvel(i-1,j) - vvel(i-1,j-1)

        ! d/du [(du/dx)^2]
        Du_dudx2_ne = fxx * (p667*dun + p333*dus)
        Du_dudx2_nw = -Du_dudx2_ne
        Du_dudx2_se = fxx * (p333*dun + p667*dus)
        Du_dudx2_sw = -Du_dudx2_se

        ! d/dv [(dv/dx)^2]
        Dv_dvdx2_ne = fxx * (p667*dvn + p333*dvs)
        Dv_dvdx2_nw = -Dv_dvdx2_ne
        Dv_dvdx2_se = fxx * (p333*dvn + p667*dvs)
        Dv_dvdx2_sw = -Dv_dvdx2_se

        ! d/du [(du/dy)^2]
        Du_dudy2_ne = fyy * (p667*due + p333*duw)
        Du_dudy2_nw = fyy * (p333*due + p667*duw)
        Du_dudy2_se = -Du_dudy2_ne
        Du_dudy2_sw = -Du_dudy2_nw

        ! d/dv [(dv/dy)^2]
        Dv_dvdy2_ne = fyy * (p667*dve + p333*dvw)
        Dv_dvdy2_nw = fyy * (p333*dve + p667*dvw)
        Dv_dvdy2_se = -Dv_dvdy2_ne
        Dv_dvdy2_sw = -Dv_dvdy2_nw

        ! d/du [(du/dx)*(dv/dy)]
        Du_dudxdvdy_ne = fxy * p25 * (dve + dvw)
        Du_dudxdvdy_nw = -Du_dudxdvdy_ne
        Du_dudxdvdy_se =  Du_dudxdvdy_ne
        Du_dudxdvdy_sw = -Du_dudxdvdy_ne

        ! d/dv [(du/dx)*(dv/dy)]
        Dv_dudxdvdy_ne = fxy * p25 * (dun + dus)
        Dv_dudxdvdy_nw =  Du_dudxdvdy_ne
        Dv_dudxdvdy_se = -Du_dudxdvdy_ne
        Dv_dudxdvdy_sw = -Du_dudxdvdy_ne

        ! d/du [(dv/dx)*(du/dy)]
        Du_dvdxdudy_ne = fxy * p25 * (dvn + dvs)
        Du_dvdxdudy_nw =  Du_dvdxdudy_ne
        Du_dvdxdudy_se = -Du_dvdxdudy_ne
        Du_dvdxdudy_sw = -Du_dvdxdudy_ne

        ! d/dv [(dv/dx)*(du/dy)]
        Dv_dvdxdudy_ne = fxy * p25 * (due + duw)
        Dv_dvdxdudy_nw = -Du_dvdxdudy_ne
        Dv_dvdxdudy_se =  Du_dvdxdudy_ne
        Dv_dvdxdudy_sw = -Du_dvdxdudy_ne

        Du_ne(i,j) = c2 * Du_dudx2_ne + p5 * Du_dudy2_ne   &
                   + c2 * Du_dudxdvdy_ne + Du_dvdxdudy_ne

        Du_nw(i,j) = c2 * Du_dudx2_nw + p5 * Du_dudy2_nw   &
                   + c2 * Du_dudxdvdy_nw + Du_dvdxdudy_nw

        Du_se(i,j) = c2 * Du_dudx2_se + p5 * Du_dudy2_se   &
                   + c2 * Du_dudxdvdy_se + Du_dvdxdudy_se

        Du_sw(i,j) = c2 * Du_dudx2_sw + p5 * Du_dudy2_sw   &
                   + c2 * Du_dudxdvdy_sw + Du_dvdxdudy_sw

        Dv_ne(i,j) = c2 * Dv_dvdx2_ne + p5 * Dv_dvdy2_ne   &
                   + c2 * Dv_dudxdvdy_ne + Dv_dvdxdudy_ne

        Dv_nw(i,j) = c2 * Dv_dvdx2_nw + p5 * Dv_dvdy2_nw   &
                   + c2 * Dv_dudxdvdy_nw + Dv_dvdxdudy_nw

        Dv_se(i,j) = c2 * Dv_dvdx2_se + p5 * Dv_dvdy2_se   &
                   + c2 * Dv_dudxdvdy_se + Dv_dvdxdudy_se

        Dv_sw(i,j) = c2 * Dv_dvdx2_sw + p5 * Dv_dvdy2_sw   &
                   + c2 * Dv_dudxdvdy_sw + Dv_dvdxdudy_sw

     endif  ! thckmask

  enddo     ! i
  enddo     ! j

  end subroutine calc_Du_2d

!****************************************************************************

  end module glissade_velo
!*************************************************************************
