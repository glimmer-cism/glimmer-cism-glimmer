
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This code is taken from the Generic Mapping Tools and was 
! converted into Fortran 90 by Ian Rutt.
!
! Original code (in C) Copyright (c) 1991-2003 by P. Wessel and W. H. F. Smith
!
! Partial translation into Fortran 90 (c) 2004 Ian C. Rutt
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
!
! The Generic Mapping Tools are maintained by Paul Wessel and 
! Walter H. F. Smith. The GMT homepage is:
!
! http://gmt.soest.hawaii.edu/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module gmt

!*FD This code is taken from the Generic Mapping Tools and was 
!*FD converted into Fortran 90 by Ian Rutt.
!*FD 
!*FD Original code (in C) Copyright (c) 1991-2003 by P. Wessel and W. H. F. Smith
!*FD
!*FD Partial translation into Fortran 90 (c) 2004 Ian C. Rutt
!*FD 
!*FD To use this library:
!*FD \begin{itemize}
!*FD \item Include a use gmt statement at the start of the program;
!*FD \item Declare an instance of gmt\_pinf for each instance of the projection you
!*FD    want to use;
!*FD \item initialise the projection by calling \texttt{gmt\_set(latc,lonc,radea,project\_info)}
!*FD    with the centre of the projection (latc,lonc) in degrees, the radius of the earth
!*FD    (radea), in metres (though often it's convenient to set this to 1.0), and the
!*FD    instance of gmt\_pinf.
!*FD \end{itemize}
!*FD
!*FD NB: Longitude and latitude coordinates are in degrees, relative to the global origin,
!*FD     while x and y are relative to the centre of the projection

  use glimmer_global

  implicit none

  type gmt_pinf
    real(dp) :: cosp               !*FD cos ($\phi_p$)
    real(dp) :: sinp               !*FD sin ($\phi_p$)
    real(dp) :: pole               !*FD +90 or -90, depending on hemisphere
    real(dp) :: central_meridian   !*FD central meridian of projection
    real(dp) :: EQ_RAD             !*FD radius of earth
    real(dp) :: i_EQ_RAD           !*FD 1/radius of earth
    real(dp) :: Dx                 !*FD fudge factor for scaling the projection 
    real(dp) :: Dy                 !*FD ditto
    real(dp) :: iDx                !*FD inverse of Dx
    real(dp) :: iDy                !*FD inverse of Dy
    real(dp) :: s_c                !*FD Not sure what this does\ldots
    real(dp) :: s_ic               !*FD Not sure what this does, but it's the inverse of {\tt s\_ic}\ldots
    logical :: n_polar            !*FD not needed - indicates if projection pole is north pole
    logical :: s_polar            !*FD not needed - indicates if projection pole is south pole
    logical :: north_pole         !*FD true if projection in northern hemisphere, false otherwise
    logical :: polar              !*FD Is this a polar projection?
  end type gmt_pinf

  real(dp),parameter :: pi=3.141592654          ! public
  real(dp),parameter :: M_PI_4=pi/4
  real(dp),parameter :: M_PI_2=pi/2
  real(dp),parameter :: D2R=pi/180.0            ! public
  real(dp),parameter :: R2D=180.0/pi            ! public
  real(dp),parameter :: DBL_MAX=huge(pi)        
  real(dp),parameter :: GMT_CONV_LIMIT=1.0e-8
  logical,parameter :: GMT_convert_latitudes=.false.
  real(dp),parameter :: GMT_map_scale_factor = 1.0
  real(dp),parameter :: SMALL=1.0e-4

  private sincos,hypot
  private DBL_MAX,GMT_CONV_LIMIT,GMT_convert_latitudes,M_PI_4

contains

  subroutine gmt_set(p_type,latc,lonc,radea,project_info)

    !*FD Initialise a GMT parameter type

    integer,intent(in) :: p_type  !*FD Type of projection to select
                                  !*FD 
                                  !*FD The projections available are:
                                  !*FD \begin{enumerate}
                                  !*FD \item Lambert Equal Area
                                  !*FD \item Spherical polar
                                  !*FD \item Spherical stereographic (oblique)
                                  !*FD \item Spherical stereographic (equatorial)
                                  !*FD \end{enumerate}
    real(dp),intent(in) :: latc    !*FD Centre of projection (Latitude)
    real(dp),intent(in) :: lonc    !*FD Centre of projection (Longitude)
    real(dp),intent(in) :: radea   !*FD Radius of the Earth (m)
    type(gmt_pinf),intent(inout) :: project_info !*FD GMT parameters to be initialised

    call gmt_type_init(project_info)

    project_info%EQ_RAD=radea                 ! radius of earth
    project_info%i_EQ_RAD=1.0/radea           ! 1/radius of earth

    select case(p_type)
    case(1)   ! Lambert equal area 
      call gmt_lambeq_set(latc,lonc,project_info)
    case(2:4)   ! Spherical/stereographic (polar,oblique equatorial) 
      call gmt_map_init_stereo(latc,lonc,project_info)
    case default
      print*,'* ERROR: ',p_type,' is not a valid projection type.'
      stop 
    end select

  end subroutine gmt_set

  subroutine gmt_type_init(project_info)

    !*FD Initialise the GMT projection parameters to their
    !*FD default values.

    type(gmt_pinf),intent(out) :: project_info !*FD GMT parameters to be initialised

    project_info%cosp = 0.0
    project_info%sinp = 0.0
    project_info%pole = 0.0
    project_info%central_meridian = 0.0
    project_info%EQ_RAD = 0.0
    project_info%i_EQ_RAD = 0.0
    project_info%Dx = 0.0
    project_info%Dy = 0.0
    project_info%iDx = 0.0
    project_info%iDy = 0.0
    project_info%s_c = 0.0
    project_info%s_ic = 0.0
    project_info%n_polar = .false.
    project_info%s_polar = .false.
    project_info%north_pole = .false.
    project_info%polar = .false.

  end subroutine gmt_type_init

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! PROJECTION CONVERSION CODES
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Lambert equal area
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_lambeq_set(latc,lonc,project_info)

    implicit none

    real(dp),intent(in) :: latc,lonc
    type(gmt_pinf),intent(inout) :: project_info

    project_info%cosp=cos(latc*D2R)           ! cos (phi_p)
    project_info%sinp=sin(latc*D2R)           ! sin (phi_p)

    project_info%central_meridian=lonc        ! central meridian of projection
    project_info%pole=latc                    ! centre of map (latitude)

    ! These next four are not used in this incomplete translation to f90

    project_info%Dx=1.0                       ! Fudge factor for scaling the projection 
    project_info%Dy=1.0                       ! ditto
    project_info%iDx=1.0                      ! inverse of Dx
    project_info%iDy=1.0                      ! inverse of Dy

  end subroutine gmt_lambeq_set

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_lambeq(lon,lat,x,y,project_info)

    !  Convert lon/lat to Spherical Lambert Azimuthal Equal-Area x/y

    implicit none

    real(dp),intent(in)  :: lon,lat
    type(gmt_pinf),intent(in) :: project_info
    real(dp),intent(out) :: x,y
    real(dp) :: k,tmp,sin_lat,cos_lat,sin_lon,cos_lon,c,dlon,dlat

    dlon = lon-project_info%central_meridian

    do while (dlon.lt.-180.0)
      dlon=dlon+360.0
    enddo

    do while (dlon.gt.180.0)
      dlon=dlon-360.0
    enddo

    dlon = dlon*D2R
    dlat = lat*D2R

    call sincos(dlat,sin_lat,cos_lat);
    call sincos(dlon,sin_lon,cos_lon);
    c = cos_lat * cos_lon

    tmp = 1.0 + project_info%sinp * sin_lat + project_info%cosp * c

    if (tmp > 0.0) then
      k = project_info%EQ_RAD * sqrt (2.0 / tmp)
      x = k * cos_lat * sin_lon
      y = k * (project_info%cosp * sin_lat - project_info%sinp * c)
      if (GMT_convert_latitudes) then  ! Gotta fudge abit 
        x = x*project_info%Dx
        y = y*project_info%Dy
      endif
    else
      x = -DBL_MAX ; y=-DBL_MAX
    endif    

  end subroutine gmt_lambeq

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_ilambeq(lon,lat,x,y,project_info)

    ! Convert Lambert Azimuthal Equal-Area x/yto lon/lat

    implicit none

    real(dp),intent(in) :: x,y
    type(gmt_pinf),intent(in) :: project_info
    real(dp),intent(out) :: lon,lat
    real(dp) :: rho,c,sin_c,cos_c,xx,yy
  
    xx=x ; yy=y

    if (GMT_convert_latitudes) then  ! Undo effect of fudge factors
      xx = xx*project_info%iDx
      yy = yy*project_info%iDy
    endif

    rho=hypot (xx,yy)
  
    if (abs(rho) < GMT_CONV_LIMIT) then
      lat = project_info%pole
      lon = project_info%central_meridian
    else
       c = 2.0 * asin(0.5 * rho * project_info%i_EQ_RAD)
      call sincos (c, sin_c, cos_c)
      lat = asin (cos_c * project_info%sinp + (yy * sin_c * project_info%cosp / rho)) * R2D
      if (project_info%n_polar) then
        lon = project_info%central_meridian + R2D * atan2 (xx, -yy)
      else if (project_info%s_polar) then
        lon = project_info%central_meridian + R2D * atan2 (xx, yy)
      else
        lon = project_info%central_meridian + &
          R2D * atan2 (xx * sin_c, (rho * project_info%cosp * cos_c - yy * project_info%sinp * sin_c))
      endif
    endif

  end subroutine gmt_ilambeq

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Spherical polar / Polar Stereographic
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_map_init_stereo(latg,rlong,project_info)

    real(dp) :: latg      ! Centre of map projection
    real(dp) :: rlong     ! Centre of map projection
    type(gmt_pinf),intent(inout) :: project_info    

    call gmt_set_polar (latg,project_info)

    ! Equatorial view has a problem with infinite loops.  Untill I find a cure
    !  we set projection center latitude to 0.001 so equatorial works for now 

    if (abs (latg) < SMALL) latg = 0.001

    call gmt_vstereo(rlong,latg,project_info)

  end subroutine gmt_map_init_stereo

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_set_polar (plat,project_info)

    ! Determines if the projection pole is N or S pole 

    implicit none

    real(dp),intent(in) :: plat
    type(gmt_pinf),intent(inout) :: project_info

    if (abs (abs (plat) - 90.0) < GMT_CONV_LIMIT) then
      project_info%polar = .true.
      project_info%north_pole  = (plat > 0.0)
       project_info%n_polar = project_info%north_pole
      project_info%s_polar = (.not.project_info%n_polar)
    endif

  end subroutine gmt_set_polar

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_vstereo (rlong0,plat,project_info)

    ! Set up a Stereographic transformation 

    implicit none 

    real(dp),intent(in) :: rlong0,plat
    type(gmt_pinf),intent(inout) :: project_info
    real(dp) :: clat

    clat = plat

    project_info%central_meridian = rlong0
    project_info%pole = plat          ! This is always geodetic
    project_info%sinp = sin (clat)   ! These may be conformal
    project_info%cosp = cos (clat)
    project_info%north_pole = (plat > 0.0)
    project_info%s_c = 2.0 * project_info%EQ_RAD * GMT_map_scale_factor
    project_info%s_ic = 1.0 / project_info%s_c

  end subroutine gmt_vstereo

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_plrs_sph(lon,lat,x,y,project_info)

  ! Convert lon/lat to x/y using Spherical polar projection

    implicit none

    real(dp),intent(in) :: lon,lat
    real(dp),intent(out) :: x,y
    type(gmt_pinf),intent(in) :: project_info

    real(dp) :: rho, slon, clon,dlon,dlat
      
    dlon = lon-project_info%central_meridian
    dlat = lat

    do 
    if (.not.(dlon < -180.0)) exit
    dlon = dlon+360.0
    enddo 

    do
    if (.not.(dlon > 180.0)) exit
    dlon = dlon-360.0
    enddo  

    dlon = dlon*D2R
    call sincos (lon,slon,clon)

     if (project_info%north_pole) then
      rho = project_info%s_c * tan (M_PI_4 - 0.5 * D2R * dlat)
      y = -rho * clon
      x =  rho * slon
    else
      rho = project_info%s_c * tan (M_PI_4 + 0.5 * D2R * dlat)
      y = rho * clon
       x = rho * slon
    endif

  end subroutine gmt_plrs_sph

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_iplrs_sph (lon,lat,x,y,project_info)

    ! Convert Spherical polar x/y to lon/lat

    implicit none

    real(dp),intent(out) :: lon,lat
    real(dp),intent(in)  :: x,y
    type(gmt_pinf),intent(in) :: project_info

    real(dp) :: c
    real(dp) :: xx,yy

    xx=x ; yy=y
    
    if (x == 0.0.and.y == 0.0) then
      lon = project_info%central_meridian
      lat = project_info%pole
      return
    endif

    c = 2.0 * atan (hypot (xx, yy) * project_info%s_ic)

    if (project_info%north_pole) then
      lon = project_info%central_meridian + d_atan2 (xx, -yy) * R2D  
      lat = d_asin (cos (c)) * R2D
    else
      lon = project_info%central_meridian + d_atan2 (xx, yy) * R2D  
      lat = d_asin (-cos (c)) * R2D
    endif

  end subroutine gmt_iplrs_sph

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_stereo1_sph (lon,lat,x,y,project_info)

    ! Convert lon/lat to x/y using Spherical stereographic projection, oblique view
  
    implicit none

    real(dp),intent(in) :: lon,lat
    real(dp),intent(out) :: x,y
    type(gmt_pinf),intent(in) :: project_info

    real(dp) :: dlon,sin_dlon,cos_dlon,s,c,cc,A,dlat

    dlon = D2R * (lon - project_info%central_meridian)
    dlat = lat*D2R
    call sincos (dlon,sin_dlon,cos_dlon)
    call sincos (dlat,s,c)
    cc=c*cos_dlon
    A = project_info%s_c / (1.0 + project_info%sinp * s + project_info%cosp * cc)
    x = A * c * sin_dlon
    y = A * (project_info%cosp * s - project_info%sinp * cc)

  end subroutine gmt_stereo1_sph

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_istereo_sph (lon,lat,x,y,project_info)

    implicit none

    real(dp),intent(out) :: lon,lat
    real(dp),intent(in)  :: x,y
    type(gmt_pinf),intent(in) :: project_info

    real(dp) :: rho,c,sin_c,cos_c

    if (x == 0.0 .and. y == 0.0) then
      lon = project_info%central_meridian
      lat = project_info%pole
    else 
      rho = hypot(x,y)
      c = 2.0 * atan (rho * project_info%s_ic)
      call sincos (c, sin_c, cos_c)
      lat = d_asin (cos_c * project_info%sinp + (y * sin_c * project_info%cosp / rho)) * R2D
      lon = R2D * atan (x * sin_c / (rho * project_info%cosp * cos_c - y * project_info%sinp * sin_c)) &
             + project_info%central_meridian
    endif

  end subroutine gmt_istereo_sph

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_stereo2_sph (lon,lat,x,y,project_info)

    ! Spherical equatorial view
    ! Convert lon/lat to x/y using stereographic projection, equatorial view 

    implicit none

    real(dp),intent(in) :: lon,lat
    real(dp),intent(out) :: x,y
    type(gmt_pinf),intent(in) :: project_info

    real(dp) :: dlon,dlat,s,c,clon,slon,A
    
    dlon = lon - project_info%central_meridian
    if (abs (dlon - 180.0) < GMT_CONV_LIMIT) then
      x = 0.0
      y = 0.0
    else 
      dlon = dlon*D2R
      dlat = lat*D2R
      call sincos (dlat,s,c)
      call sincos (dlon,slon,clon)
      A = project_info%s_c / (1.0 + c * clon)
      x = A * c * slon
      y = A * s
    endif

  end subroutine gmt_stereo2_sph

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Utility routines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sincos(a,s,c)

    implicit none

    real(dp),intent(in) :: a
    real(dp),intent(out) :: s,c

      s = sin (a)
      c = cos (a)

  end subroutine sincos

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(dp) function hypot(x,y)

  ! This is an implementation of a standard C library function

    implicit none

    real(dp),intent(in) :: x,y

    hypot=sqrt(x*x+y*y)

  end function hypot

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(dp) function d_atan2(x,y)

  ! Macro-implemented function from GMT

    implicit none

    real(dp),intent(in) :: x,y

    if (x==0.0 .and. y==0.0) then
      d_atan2=0.0
    else
      d_atan2=atan2 (y, x)
    endif

  end function d_atan2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(dp) function d_asin(x)

  ! Macro-implemented function from GMT

    implicit none

    real(dp),intent(in) :: x

    if (abs(x) >= 1.0) then
      d_asin=sign(M_PI_2,x)
    else
      d_asin=asin(x)
    endif

  end function d_asin

end module gmt
