! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_map_init.f90 - part of the GLIMMER ice model     + 
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

module glimmer_map_init

  !*FD Contains initialisation routines for 
  !*FD individual projection types.

  use glimmer_map_types

  implicit none

contains

  subroutine glimmap_laea_init(params)

    type(proj_laea),intent(inout) :: params

    params%sinp=sin(params%latitude_of_projection_origin*D2R)
    params%cosp=cos(params%latitude_of_projection_origin*D2R)

    ! Check whether polar

    if (abs(params%latitude_of_projection_origin-90.0)<CONV_LIMIT) then
       params%pole=1
    else if (abs(params%latitude_of_projection_origin+90)<CONV_LIMIT) then
       params%pole=-1
    else
       params%pole=0
    end if

  end subroutine glimmap_laea_init

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glimmap_aea_init(params)

    type(proj_aea),intent(inout) :: params

    params%n = 0.5*(sin(params%standard_parallel(1)*D2R) &
         + sin(params%standard_parallel(2)*D2R))
    params%i_n = 1.0/params%n
    params%c = cos(params%standard_parallel(1)*D2R)**2.0 &
         + 2.0*params%n*sin(params%standard_parallel(1)*D2R)
    params%rho0_R = params%i_n * sqrt(params%c - &
         2.0*params%n*sin(params%latitude_of_projection_origin*D2R))
    params%rho0 = params%rho0_R * EQ_RAD

  end subroutine glimmap_aea_init

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glimmap_lcc_init(params)

    type(proj_lcc),intent(inout) :: params

    if (abs(params%standard_parallel(1)-params%standard_parallel(2))<CONV_LIMIT) then
       params%n = sin(params%standard_parallel(1)*D2R)
    else
       params%n = log(cos(params%standard_parallel(1)*D2R)/cos(params%standard_parallel(2)*D2R))/ &
            log(tan(M_PI_4+params%standard_parallel(2)*D2R/2.0)/tan(M_PI_4+params%standard_parallel(1)*D2R/2.0))
    end if

    params%f = params%i_n*cos(params%standard_parallel(1)*D2R)* &
         (tan(M_PI_4+params%standard_parallel(1)*D2R/2.0))**params%n
    params%rho0 = EQ_RAD * params%f/(tan(M_PI_4+params%latitude_of_projection_origin*D2R/2.0))**params%n

  end subroutine glimmap_lcc_init

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glimmap_stere_init(params)

    use glimmer_log

    type(proj_stere),intent(inout) :: params

    ! Determine polar/equatorial, etc.

    if (abs(params%latitude_of_projection_origin-90.0)<CONV_LIMIT) then
       params%pole=1
    else if (abs(params%latitude_of_projection_origin+90)<CONV_LIMIT) then
       params%pole=-1
    else
       params%pole=0
       if (abs(params%latitude_of_projection_origin)<CONV_LIMIT) then
          params%equatorial = .true.
       else
          params%equatorial = .false.
       end if
    end if

    ! Set up constants accordingly

    if (params%pole==1.or.params%pole==-1) then
       if (params%standard_parallel/=0.0) then
          params%k0 = EQ_RAD * (1 + sin(D2R*params%standard_parallel))/2.0
       else if (params%scale_factor_at_proj_origin/=0.0) then
          params%k0 = EQ_RAD * params%scale_factor_at_proj_origin
       else
          params%k0 = EQ_RAD
       end if
    else
       if (params%scale_factor_at_proj_origin/=0.0) then
          params%k0 = EQ_RAD * params%scale_factor_at_proj_origin
       else
          params%k0 = EQ_RAD
       end if
       if (params%standard_parallel/=0.0) &
            call write_log('Stereographic projection not polar: ignoring standard parallel',GM_WARNING)
       params%sinp = sin(D2R * params%latitude_of_projection_origin)
       params%cosp = cos(D2R * params%latitude_of_projection_origin)
    end if

  end subroutine glimmap_stere_init

end module glimmer_map_init
