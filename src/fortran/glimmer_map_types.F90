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

  type glimmap_proj
     logical :: found = .false.
     type(proj_laea),  pointer :: laea  => NULL() !*FD Pointer to Lambert azimuthal equal area type
     type(proj_aea),   pointer :: aea   => NULL() !*FD Pointer to Albers equal area conic type
     type(proj_lcc),   pointer :: lcc   => NULL() !*FD Pointer to Lambert conic conformal type
     type(proj_stere), pointer :: stere => NULL() !*FD Pointer to Stereographic type
  end type glimmap_proj

  !-------------------------------------------------------------

  type proj_laea
     !*FD Lambert Azimuthal Equal Area
     real :: longitude_of_central_meridian
     real :: latitude_of_projection_origin
     real :: false_easting
     real :: false_northing
     real :: sinp    !*FD Sine of latitude_of_projection_origin
     real :: cosp    !*FD Cosine of latitude_of_projection_origin
     integer :: pole !*FD Set to 1 for N pole, -1 for S pole, 0 otherwise
  end type proj_laea

  !-------------------------------------------------------------

  type proj_aea
     !*FD Albers Equal-Area Conic
     real,dimension(2) :: standard_parallel
     real :: longitude_of_central_meridian
     real :: latitude_of_projection_origin
     real :: false_easting
     real :: false_northing
     real :: rho0 !*FD Convenience constant
     real :: rho0_R !*FD Convenience constant (is rho0/EQ_RAD)
     real :: c    !*FD Convenience constant
     real :: n    !*FD Convenience constant
     real :: i_n  !*FD Convenience constant (inverse of n)
  end type proj_aea

  !-------------------------------------------------------------

  type proj_lcc
     !*FD Lambert Conic Conformal
     real,dimension(2) :: standard_parallel
     real :: longitude_of_central_meridian
     real :: latitude_of_projection_origin
     real :: false_easting
     real :: false_northing
     real :: rho0 !*FD Convenience constant
     real :: f    !*FD Convenience constant
     real :: n    !*FD Convenience constant
     real :: i_n  !*FD Convenience constant (inverse of n)
  end type proj_lcc

  !-------------------------------------------------------------

  type proj_stere
     !*FD Stereographic projection derived type 
     real :: longitude_of_central_meridian          
     real :: latitude_of_projection_origin          
     real :: scale_factor_at_proj_origin = 0. 
     real :: standard_parallel = 0.                 
     real :: false_easting          
     real :: false_northing
     integer :: pole  !*FD Set to 1 for N pole, -1 for S pole, 0 otherwise
     logical :: equatorial !*FD Set true if equatorial aspect
     real :: k0 !*FD scale factor or std par converted to scale factor
     real :: sinp,cosp !*FD sin and cos of latitude_of_projection_origin
  end type proj_stere

  ! Global mapping parameters ----------------------------------

  real(rk),parameter :: pi         = 3.141592654    !*FD The value of $\pi$.
  real(rk),parameter :: M_PI_4     = pi/4           !*FD The value of $\pi/4$.
  real(rk),parameter :: M_PI_2     = pi/2           !*FD The value of $\pi/2$.
  real(rk),parameter :: D2R        = pi/180.0       !*FD Degrees-to-radians conversion factor.
  real(rk),parameter :: R2D        = 180.0/pi       !*FD Radians-to-degrees conversion factor.
  real(rk),parameter :: EQ_RAD     = 6.37e6         !*FD Radius of the earth (m)
  real(rk),parameter :: i_EQ_RAD   = 1.0_rk/EQ_RAD  !*FD Inverse radius of the earth (m^-1)
  real(rk),parameter :: CONV_LIMIT = 1.0e-8         !*FD Convergence limit (a small number).

  integer, parameter :: GMAP_LAEA=1
  integer, parameter :: GMAP_AEA=2
  integer, parameter :: GMAP_LCC=3
  integer, parameter :: GMAP_STERE=4

contains

  function glimmap_allocated(proj)

    !*FD return true if structure contains a known projection

    implicit none
    type(glimmap_proj) :: proj
    logical glimmap_allocated

    glimmap_allocated = proj%found
  end function glimmap_allocated

  !-------------------------------------------------------------------------

  subroutine glimmap_proj_define(cfp,ptype, &
       longitude_of_central_meridian, &
       latitude_of_projection_origin, &
       false_easting, &
       false_northing, &
       scale_factor_at_proj_origin, &
       standard_parallel, &
       standard_parallel_2)

    !*FD Defines a projection from scratch, and initialises 
    !*FD the other elements appropriately.

    use glimmer_log

    type(glimmap_proj),intent(inout) :: cfp 
    integer,intent(in) :: ptype
    real,intent(in) :: longitude_of_central_meridian
    real,intent(in) :: latitude_of_projection_origin
    real,intent(in) :: false_easting
    real,intent(in) :: false_northing
    real,optional,intent(in) :: scale_factor_at_proj_origin
    real,optional,intent(in) :: standard_parallel
    real,optional,intent(in) :: standard_parallel_2

    if (associated(cfp%laea))  deallocate(cfp%laea)
    if (associated(cfp%aea))   deallocate(cfp%aea)
    if (associated(cfp%lcc))   deallocate(cfp%lcc)
    if (associated(cfp%stere)) deallocate(cfp%stere)

    cfp%found = .true.
    select case(ptype)
    case(GMAP_LAEA)
       allocate(cfp%laea)
       cfp%laea%longitude_of_central_meridian = longitude_of_central_meridian
       cfp%laea%latitude_of_projection_origin = latitude_of_projection_origin
       cfp%laea%false_easting  = false_easting
       cfp%laea%false_northing = false_northing
       call glimmap_laea_init(cfp%laea)
    case(GMAP_AEA)
       allocate(cfp%aea)
       cfp%aea%longitude_of_central_meridian = longitude_of_central_meridian
       cfp%aea%latitude_of_projection_origin = latitude_of_projection_origin
       cfp%aea%false_easting  = false_easting
       cfp%aea%false_northing = false_northing
       if (present(standard_parallel).and.present(standard_parallel_2)) then
          cfp%aea%standard_parallel = (/ standard_parallel,standard_parallel_2 /)
       else if (present(standard_parallel).and..not.present(standard_parallel_2)) then
          cfp%aea%standard_parallel = (/ standard_parallel,standard_parallel /)
       else
          call write_log('Albers Equal Area: you must supply at least one standard parallel',&
               GM_FATAL,__FILE__,__LINE__)
       end if
       call glimmap_aea_init(cfp%aea)
    case(GMAP_LCC)
       allocate(cfp%lcc)
       cfp%lcc%longitude_of_central_meridian = longitude_of_central_meridian
       cfp%lcc%latitude_of_projection_origin = latitude_of_projection_origin
       cfp%lcc%false_easting  = false_easting
       cfp%lcc%false_northing = false_northing
       if (present(standard_parallel).and.present(standard_parallel_2)) then
          cfp%aea%standard_parallel = (/ standard_parallel,standard_parallel_2 /)
       else if (present(standard_parallel).and..not.present(standard_parallel_2)) then
          cfp%aea%standard_parallel = (/ standard_parallel,standard_parallel /)
       else
          call write_log('Lambert Conformal Conic: you must supply at least one standard parallel',&
               GM_FATAL,__FILE__,__LINE__)
       end if
       call glimmap_lcc_init(cfp%lcc)
    case(GMAP_STERE)
       allocate(cfp%stere)
       cfp%stere%longitude_of_central_meridian = longitude_of_central_meridian
       cfp%stere%latitude_of_projection_origin = latitude_of_projection_origin
       cfp%stere%false_easting  = false_easting
       cfp%stere%false_northing = false_northing
       if(present(scale_factor_at_proj_origin).and.present(standard_parallel)) then
          if (scale_factor_at_proj_origin/=0.0.and.standard_parallel/=0.0) &
               call write_log('Both standard parallel and scale factor specified',GM_FATAL,__FILE__,__LINE__)
       end if
       if(present(scale_factor_at_proj_origin)) &
            cfp%stere%scale_factor_at_proj_origin = scale_factor_at_proj_origin
       if(present(standard_parallel)) &
            cfp%stere%standard_parallel = standard_parallel
       call glimmap_stere_init(cfp%stere)
    case default
       call write_log('Unrecognised projection type',GM_FATAL,__FILE__,__LINE__)
    end select

  end subroutine glimmap_proj_define

end module glimmer_map_types
