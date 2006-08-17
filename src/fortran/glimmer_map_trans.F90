! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_map_trans.f90 - part of the GLIMMER ice model    + 
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

module glimmer_map_trans

  use glimmer_map_types

  implicit none

  private
  public :: ll_to_xy, xy_to_ll

contains

  subroutine ll_to_xy(lon,lat,x,y,proj,grid)

    !*FD Convert lat-long coordinates to grid coordinates. Note that
    !*FD The subroutine returns the x-y coordinates as real values,
    !*FD non-integer values indicating a position between grid-points.

    use glimmer_log
    use glimmer_coordinates

    implicit none

    real(rk),intent(in)  :: lon !*FD The location of the point in lat-lon space (Longitude)
    real(rk),intent(in)  :: lat !*FD The location of the point in lat-lon space (Latitude)
    real(rk),intent(out) :: x   !*FD The location of the point in $x$-$y$ space ($x$ coordinate)
    real(rk),intent(out) :: y   !*FD The location of the point in $x$-$y$ space ($y$ coordinate)
    type(projection),      intent(in) :: proj !*FD The projection being used
    type(coordsystem_type),intent(in) :: grid !*FD the grid definition

    real(rk) :: xx,yy ! These are real-space distances in meters

    if (associated(proj%laea)) then
       call glimmap_laea(lon,lat,xx,yy,proj%laea)
    else if (associated(proj%aea)) then
   !    call glimmap_albers(lon,lat,xx,yy,proj%aea)
    else if (associated(proj%lcc)) then
   !    call glimmap_lambconf(lon,lat,xx,yy,proj%lcc)
    else if (associated(proj%stere)) then
   !    call glimmap_stereo(lon,lat,xx,yy,proj%stere)
    else
       call write_log('No known projection found!',GM_WARNING)
    end if

    ! Now convert the real-space distances to grid-points using the grid type

    call coordsystem_space2grid(xx,yy,x,y,grid)

  end subroutine ll_to_xy

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine xy_to_ll(lon,lat,x,y,proj,grid)

    !*FD Convert grid coordinates to lat-lon coordinates. The 
    !*FD subroutine returns the lat-lon coordinates as real values,
    !*FD non-integer values indicating a position between grid-points.

    use glimmer_log
    use glimmer_coordinates

    implicit none

    real(rk),intent(out) :: lon !*FD The location of the point in lat-lon space (Longitude)
    real(rk),intent(out) :: lat !*FD The location of the point in lat-lon space (Latitude)
    real(rk),intent(in)  :: x   !*FD The location of the point in $x$-$y$ space ($x$ coordinate)
    real(rk),intent(in)  :: y   !*FD The location of the point in $x$-$y$ space ($y$ coordinate)
    type(projection),      intent(in) :: proj !*FD The projection being used
    type(coordsystem_type),intent(in) :: grid !*FD the grid definition

    real(rk) :: xx,yy ! These are real-space distances in meters

    ! First convert grid-point space to real space

    call coordsystem_grid2space(xx,yy,x,y,grid)

    if (associated(proj%laea)) then
       call glimmap_ilaea(lon,lat,xx,yy,proj%laea)
    else if (associated(proj%aea)) then
   !    call glimmap_ialbers(lon,lat,xx,yy,proj%aea)
    else if (associated(proj%lcc)) then
   !    call glimmap_ilambconf(lon,lat,xx,yy,proj%lcc)
    else if (associated(proj%stere)) then
   !    call glimmap_istereo(lon,lat,xx,yy,proj%stere)
    else
       call write_log('No known projection found!',GM_WARNING)
    end if

    lon=loncorrect(lon)

  end subroutine xy_to_ll
  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! PRIVATE subroutines follow
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Lambert azimuthal equal area projection
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glimmap_laea(lon,lat,x,y,params)

    !*FD Forward transformation: lat-lon -> x-y

    use glimmer_log

    real(rk),intent(in)  :: lon
    real(rk),intent(in)  :: lat
    real(rk),intent(out) :: x
    real(rk),intent(out) :: y
    type(proj_laea),intent(in) :: params

    real(rk) :: sin_lat,cos_lat,sin_lon,cos_lon,c,dlon,dlat,tmp,k
    character(80) :: errtxt

    dlon = lon-params%longitude_of_central_meridian

    ! Check domain of longitude

    do while (dlon.lt.-180.0)
       dlon=dlon+360.0
    enddo

    do while (dlon.gt.180.0)
       dlon=dlon-360.0
    enddo

    ! Convert to radians and calculate sine and cos

    dlon = dlon*D2R ; dlat = lat*D2R

    call sincos(dlon,sin_lon,cos_lon);
    call sincos(dlat,sin_lat,cos_lat);
    c = cos_lat * cos_lon

    ! Mapping transformation

    tmp = 1.0 + params%sinp * sin_lat + params%cosp * c

    if (tmp > 0.0) then
       k = EQ_RAD * sqrt (2.0 / tmp)
       x = k * cos_lat * sin_lon
       y = k * (params%cosp * sin_lat - params%sinp * c)
    else
       write(errtxt,*)'LAEA projection error:',lon,lat,params%latitude_of_projection_origin
       call write_log(trim(errtxt),GM_FATAL,__FILE__,__LINE__)
    endif

    ! Apply false eastings and northings

    x = x + params%false_easting
    y = y + params%false_northing

  end subroutine glimmap_laea

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glimmap_ilaea(lon,lat,x,y,params)

    use glimmer_log

    real(rk),intent(out) :: lon
    real(rk),intent(out) :: lat
    real(rk),intent(in)  :: x
    real(rk),intent(in)  :: y
    type(proj_laea),intent(in) :: params

    real(rk) :: rho,c,sin_c,cos_c,xx,yy
    character(80) :: errtxt

    xx=x ; yy=y

    ! Account for false eastings and northings

    xx = xx - params%false_easting
    yy = yy - params%false_northing

    rho=hypot (xx,yy)

    if (abs(rho) < CONV_LIMIT) then
       ! If very near the centre of the map...
       lat = params%latitude_of_projection_origin
       lon = params%longitude_of_central_meridian
    else
       c = 2.0 * asin(0.5 * rho * i_EQ_RAD)
       call sincos (c, sin_c, cos_c)
       lat = asin (cos_c * params%sinp + (yy * sin_c * params%cosp / rho)) * R2D
       select case(params%pole)
       case(1)
          lon = params%longitude_of_central_meridian + R2D * atan2 (xx, -yy)
       case(-1)
          lon = params%longitude_of_central_meridian + R2D * atan2 (xx, yy)
       case(0)
          lon = params%longitude_of_central_meridian + &
               R2D * atan2 (xx * sin_c, (rho * params%cosp * cos_c - yy * params%sinp * sin_c))
       case default
          write(errtxt,*)'Inverse LAEA projection error:',params%pole
          call write_log(trim(errtxt),GM_FATAL,__FILE__,__LINE__)
       end select
    endif

  end subroutine glimmap_ilaea

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Utility routines
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine grid2space(x,y,gx,gy,coordsys)

    use glimmer_coordinates

    implicit none

    real(rk),intent(out) :: x  !*FD x-location in real space
    real(rk),intent(out) :: y  !*FD y-location in real space
    real(rk),intent(in)  :: gx !*FD x-location in grid space
    real(rk),intent(in)  :: gy !*FD y-location in grid space
    type(coordsystem_type), intent(in) :: coordsys  !*FD coordinate system

    x=coordsys%origin%pt(1) + (gx - 1)*coordsys%delta%pt(1)
    y=coordsys%origin%pt(2) + (gy - 1)*coordsys%delta%pt(2)

  end subroutine grid2space

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine space2grid(x,y,gx,gy,coordsys)

    use glimmer_coordinates

    implicit none

    real(rk),intent(in)  :: x  !*FD x-location in real space
    real(rk),intent(in)  :: y  !*FD y-location in real space
    real(rk),intent(out) :: gx !*FD x-location in grid space
    real(rk),intent(out) :: gy !*FD y-location in grid space
    type(coordsystem_type), intent(in) :: coordsys  !*FD coordinate system

    gx = 1 + (x - coordsys%origin%pt(1))/coordsys%delta%pt(1)
    gy = 1 + (y - coordsys%origin%pt(2))/coordsys%delta%pt(2)

  end subroutine space2grid

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sincos(a,s,c)

    !*FD Calculates the sin and cos of an angle.

    implicit none

    real(rk),intent(in)  :: a !*FD Input value (radians).
    real(rk),intent(out) :: s !*FD sin(\texttt{a})
    real(rk),intent(out) :: c !*FD cos(\texttt{a})

    s = sin (a)
    c = cos (a)

  end subroutine sincos

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(rk) function loncorrect(lon)

    !*FD Normalises a value of longitude to the range 0 to 360 degrees.
    !*RV The normalised value of longitude.

    real(rk),intent(in) :: lon !*FD The longitude under consideration (degrees east)

    loncorrect=lon

    do while (loncorrect>360.0)
       loncorrect=loncorrect-360.0
    enddo

    do while (loncorrect<0.0)
       loncorrect=loncorrect+360.0
    enddo

  end function loncorrect

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(rk) function hypot(x,y)

    !*FD This is an implementation of a standard C library function, 
    !*FD returning the value of $\sqrt{x^2+y^2}$.
    !*RV Returns the value of $\sqrt{x^2+y^2}$.

    implicit none

    real(rk),intent(in) :: x !*FD One input value
    real(rk),intent(in) :: y !*FD Another input value

    hypot=sqrt(x*x+y*y)

  end function hypot

end module glimmer_map_trans
