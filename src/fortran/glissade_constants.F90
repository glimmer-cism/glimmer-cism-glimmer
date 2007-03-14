! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
! +                                                           + 
! +  glissade_constants.F90                                   + 
! +                                                           + 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!
! This module contains constants used by other GLISSADE modules.
!
! Note that most of the required parameters are contained
! in glimmer_physcon and glimmer_params.  The ones defined here
! are mostly some standard constants used in the POP and CICE
! models developed at Los Alamos National Lab.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

  module glissade_constants

  use glimmer_global, only: sp, dp
  implicit none

  real(dp), parameter ::   &
       c0 = 0.0_dp,                  &
       c1 = 1.0_dp,                  &
       c2 = 2.0_dp,                  &
       c3 = 3.0_dp,                  &
       p25= 0.25_dp,                 &
       p333 = c1/c3,                 &
       p5 = 0.5_dp,                  &
       eps = 1.0e-11_dp 

!------------------------------------------------------------------------

  end module glissade_constants

!------------------------------------------------------------------------
