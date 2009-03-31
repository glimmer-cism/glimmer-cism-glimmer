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
                         Du_dudx2,      Dv_dvdx2,     &
                         Du_dudy2,      Dv_dvdy2,     &
                         Du_dudxdvdy,   Dv_dudxdvdy,  &
                         Du_dvdxdudy,   Dv_dvdxdudy)

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

  real(dp), dimension(:,:,:), intent(out) ::   &
      Du_dudx2,         &!  d/du [(du/dx)^2]
      Dv_dvdx2,         &!  d/dv [(dv/dx)^2]
      Du_dudy2,         &!  d/du [(du/dy)^2]
      Dv_dvdy2,         &!  d/dv [(dv/dy)^2]
      Du_dudxdvdy,      &!  d/du [(du/dx)*(dv/dy)]
      Dv_dudxdvdy,      &!  d/dv [(du/dx)*(dv/dy)]
      Du_dvdxdudy,      &!  d/du [(dv/dx)*(du/dy)]
      Dv_dvdxdudy        !  d/dv [(dv/dx)*(du/dy)]

  !---------------------------------------------------------------------------
  ! local variables
  !---------------------------------------------------------------------------

  integer ::  i, j           ! horizontal indices

  real(dp) ::    &
     fxx, fyy, fxy,         &! geometric factors
     dun, dus, due, duw,    &! difference in uvel between adjacent nodes
     dvn, dvs, dve, dvw      ! difference in vvel between adjacent nodes

  !---------------------------------------------------------------------------
  ! For each grid cell, calculate various derivatives of the effective strain rate
  ! with respect to the nodal velocities.
  !---------------------------------------------------------------------------

  fxx = c1 / (dew*dew)
  fyy = c1 / (dns*dns)
  fxy = c1 / (dew*dns) 

  Du_dudx2(:,:,:) = c0
  Dv_dvdx2(:,:,:) = c0
  Du_dudy2(:,:,:) = c0
  Dv_dvdy2(:,:,:) = c0
  Du_dudxdvdy(:,:,:) = c0  
  Dv_dudxdvdy(:,:,:) = c0  
  Du_dvdxdudy(:,:,:) = c0  
  Dv_dvdxdudy(:,:,:) = c0  
 
!whl - to do - have data strides here with kne as 3rd index
!    - Could define Du_dudx2_ne(i,j) instead of Du_dudx2(i,j,kne)
!    - More arguments to pass, but might be faster 

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
        Du_dudx2(i,j,kne) = fxx * (p667*dun + p333*dus)
        Du_dudx2(i,j,knw) = -Du_dudx2(i,j,kne)
        Du_dudx2(i,j,kse) = fxx * (p333*dun + p667*dus)
        Du_dudx2(i,j,ksw) = -Du_dudx2(i,j,kse)

        ! d/dv [(dv/dx)^2]
        Dv_dvdx2(i,j,kne) = fxx * (p667*dvn + p333*dvs)
        Dv_dvdx2(i,j,knw) = -Dv_dvdx2(i,j,kne)
        Dv_dvdx2(i,j,kse) = fxx * (p333*dvn + p667*dvs)
        Dv_dvdx2(i,j,ksw) = -Dv_dvdx2(i,j,kse)

        ! d/du [(du/dy)^2]
        Du_dudy2(i,j,kne) = fyy * (p667*due + p333*duw)
        Du_dudy2(i,j,knw) = fyy * (p333*due + p667*duw)
        Du_dudy2(i,j,kse) = -Du_dudy2(i,j,kne)
        Du_dudy2(i,j,ksw) = -Du_dudy2(i,j,knw)

        ! d/dv [(dv/dy)^2]
        Dv_dvdy2(i,j,kne) = fyy * (p667*dve + p333*dvw)
        Dv_dvdy2(i,j,knw) = fyy * (p333*dve + p667*dvw)
        Dv_dvdy2(i,j,kse) = -Dv_dvdy2(i,j,kne)
        Dv_dvdy2(i,j,ksw) = -Dv_dvdy2(i,j,knw)

        ! d/du [(du/dx)*(dv/dy)]
        Du_dudxdvdy(i,j,kne) = fxy * p25 * (dve + dvw)
        Du_dudxdvdy(i,j,knw) = -Du_dudxdvdy(i,j,kne)
        Du_dudxdvdy(i,j,kse) =  Du_dudxdvdy(i,j,kne)
        Du_dudxdvdy(i,j,ksw) = -Du_dudxdvdy(i,j,kne)

        ! d/dv [(du/dx)*(dv/dy)]
        Dv_dudxdvdy(i,j,kne) = fxy * p25 * (dun + dus)
        Dv_dudxdvdy(i,j,knw) =  Du_dudxdvdy(i,j,kne)
        Dv_dudxdvdy(i,j,kse) = -Du_dudxdvdy(i,j,kne)
        Dv_dudxdvdy(i,j,ksw) = -Du_dudxdvdy(i,j,kne)

        ! d/du [(dv/dx)*(du/dy)]
        Du_dvdxdudy(i,j,kne) = fxy * p25 * (dvn + dvs)
        Du_dvdxdudy(i,j,knw) =  Du_dvdxdudy(i,j,kne)
        Du_dvdxdudy(i,j,kse) = -Du_dvdxdudy(i,j,kne)
        Du_dvdxdudy(i,j,ksw) = -Du_dvdxdudy(i,j,kne)

        ! d/dv [(dv/dx)*(du/dy)]
        Dv_dvdxdudy(i,j,kne) = fxy * p25 * (due + duw)
        Dv_dvdxdudy(i,j,knw) = -Du_dvdxdudy(i,j,kne)
        Dv_dvdxdudy(i,j,kse) =  Du_dvdxdudy(i,j,kne)
        Dv_dvdxdudy(i,j,ksw) = -Du_dvdxdudy(i,j,kne)

     endif  ! thckmask

  enddo     ! i
  enddo     ! j

  end subroutine calc_Du_2d

!****************************************************************************

  end module glissade_velo
!*************************************************************************
