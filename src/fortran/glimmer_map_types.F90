! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_map_types.f90 - part of the GLIMMER ice model    + 
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

module glimmer_map_types

  !*FD New combined map-projection handling code for GLIMMER:
  !*FD This module contains derived types.

  use glimmer_global, only: rk

  implicit none

  type projection
     type(proj_laea),  pointer :: laea  => NULL() !*FD Pointer to Lambert azimuthal equal area type
     type(proj_aea),   pointer :: aea   => NULL() !*FD Pointer to Albers equal area conic type
     type(proj_lcc),   pointer :: lcc   => NULL() !*FD Pointer to Lambert conic conformal type
     type(proj_stere), pointer :: stere => NULL() !*FD Pointer to Stereographic type
     real(rk),dimension(:,:),pointer :: sintheta  => NULL() !*FD sines of grid angle relative to north.
     real(rk),dimension(:,:),pointer :: costheta  => NULL() !*FD coses of grid angle relative to north.
     real(rk),dimension(:,:),pointer :: latitudes => NULL() !*FD The latitude of each grid-point
  end type projection

  !-------------------------------------------------------------

  type proj_stere
     !*FD Stereographic projection derived type 
     logical :: polar                               !*FD Polar projection?
     real :: longitude_of_central_meridian          
     real :: latitude_of_projection_origin          
     real :: scale_factor_at_proj_origin = 0. 
     real :: standard_parallel = 0.                 
     real :: false_easting          
     real :: false_northing
  end type proj_stere

  !-------------------------------------------------------------

  type proj_laea
     !*FD Lambert Azimuthal Equal Area
     real :: longitude_of_central_meridian
     real :: latitude_of_projection_origin
     real :: false_easting
     real :: false_northing
     real :: sinp    !*FD Sine of latitude_of_projection_origin
     real :: cosp    !*FD Cosine of latitude_of_projection_origin
     integer :: pole !&FD Set to 1 for N pole, -1 for S pole, 0 otherise
  end type proj_laea

  !-------------------------------------------------------------

  type proj_aea
     !*FD Albers Equal-Area Conic
     real,dimension(2) :: standard_parallel
     real :: longitude_of_central_meridian
     real :: latitude_of_projection_origin
     real :: false_easting
     real :: false_northing
  end type proj_aea

  !-------------------------------------------------------------

  type proj_lcc
     !*FD Lambert Conic Conformal
     real,dimension(2) :: standard_parallel
     real :: longitude_of_central_meridian
     real :: latitude_of_projection_origin
     real :: false_easting
     real :: false_northing
  end type proj_lcc

  ! Global mapping parameters ----------------------------------

  real(rk),parameter :: pi         = 3.141592654    !*FD The value of $\pi$.
  real(rk),parameter :: M_PI_4     = pi/4           !*FD The value of $\pi/4$.
  real(rk),parameter :: M_PI_2     = pi/2           !*FD The value of $\pi/2$.
  real(rk),parameter :: D2R        = pi/180.0       !*FD Degrees-to-radians conversion factor.
  real(rk),parameter :: R2D        = 180.0/pi       !*FD Radians-to-degrees conversion factor.
  real(rk),parameter :: EQ_RAD     = 6.37e6         !*FD Radius of the earth (m)
  real(rk),parameter :: i_EQ_RAD   = 1.0_rk/EQ_RAD  !*FD Inverse radius of the earth (m^-1)
  real(rk),parameter :: CONV_LIMIT = 1.0e-8         !*FD Convergence limit (a small number).

end module glimmer_map_types
