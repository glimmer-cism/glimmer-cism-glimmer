
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_proj.f90 - part of the GLIMMER ice model         + 
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

module glimmer_project

  !*FD Holds derived types and subroutines
  !*FD necessary for handling map projections and the associated
  !*FD grids.

  use gmt, only: gmt_pinf,D2R,R2D
  use glimmer_global

  implicit none

  type projection

    !*FD A derived type that holds information
    !*FD relating to the projected grid.
  
    integer        :: p_type                     !*FD Number indicating projection used
                                                 !*FD 
                                                 !*FD The projections available are:
                                                 !*FD \begin{enumerate}
                                                 !*FD \item Lambert Equal Area
                                                 !*FD \item Spherical polar
                                                 !*FD \item Spherical stereographic (oblique)
                                                 !*FD \item Spherical stereographic (equatorial)
                                                 !*FD \end{enumerate}
    type(gmt_pinf) :: gmt_params                 !*FD parameters needed for gmt routines
    integer        :: nx                         !*FD number of grid points in the EW direction
    integer        :: ny                         !*FD number of grid points in the NS direction
    real(rk)       :: dx                         !*FD the nominal $x$ grid-spacing at the centre of the projection
    real(rk)       :: dy                         !*FD the nominal $y$ grid-spacing at the centre of the projection
    real(rk)       :: cpx,cpy                    !*FD The location of the map projection centre within the grid ($x$ and $y$)
    real(rk)       :: latc,lonc                  !*FD The location of the projection centre in lat/lon space (lat and lon)
    real(rk),dimension(:,:),pointer :: sintheta  !*FD sines of grid angle relative to north.
    real(rk),dimension(:,:),pointer :: costheta  !*FD coses of grid angle relative to north.
    real(rk),dimension(:,:),pointer :: latitudes !*FD The latitude of each grid-point
  end type projection

  real(rk),parameter :: pi=3.141592654  !*FD The value of pi

  private pi

contains

  subroutine new_proj(proj,radea,p_type,nx,ny,dx,dy,cpx,cpy,latc,lonc)

    !*FD Initialise new map projection area. The subroutine may be used without the optional
    !*FD arguments to initialise a projection when these parameters
    !*FD have been pre-loaded (as during the glimmer initialisation)

    use gmt

    implicit none

    type(projection),intent(inout) :: proj   !*FD The projection parameters to be initialised
    real(rk),intent(in)            :: radea  !*FD The radius of the Earth (m)   
    integer, intent(in),optional   :: p_type !*FD The type of projection
    integer, intent(in),optional   :: nx     !*FD The number of grid-points in the $x$-direction
    integer, intent(in),optional   :: ny     !*FD The number of grid-points in the $y$-direction
    real(rk),intent(in),optional   :: dx     !*FD The grid-spacing in the $x$-direction (m)
    real(rk),intent(in),optional   :: dy     !*FD The grid-spacing in the $y$-direction (m)
    real(rk),intent(in),optional   :: cpx    !*FD The $x$-location of the projection centre within the grid
    real(rk),intent(in),optional   :: cpy    !*FD The $y$-location of the projection centre within the grid
    real(rk),intent(in),optional   :: latc   !*FD The latitudinal location of the projection centre (degrees north)
    real(rk),intent(in),optional   :: lonc   !*FD The longituidinal location of the projection centre (degrees east)

    if (present(p_type).and. &
        present(nx).and. &
        present(ny).and. &
        present(dx).and. &
        present(dy).and. &
        present(cpx).and. &
        present(cpy).and. &
        present(latc).and. &
        present(lonc)) then

      proj%p_type=p_type

      proj%nx=nx     ; proj%ny=ny
      proj%dx=dx     ; proj%dy=dy
      proj%cpx=cpx   ; proj%cpy=cpy
      proj%latc=latc ; proj%lonc=lonc

    endif

    call proj_allocate(proj)

    proj%sintheta=0.0 ; proj%costheta=0.0

    call gmt_set(proj%p_type,proj%latc,proj%lonc,radea,proj%gmt_params)

    call calc_grid_angle(proj)

  end subroutine new_proj

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine proj_allocate(proj,re_alloc)

    !*FD Allocates the array pointers in the \texttt{projection} type.
    !*FD \texttt{re\_alloc} is necessary because it is not possible to 
    !*FD detect whether the arrays are already allocated in Fortran 90.

    implicit none

    type(projection),intent(inout) :: proj     !*FD The projection being initialised
    logical,intent(in),optional    :: re_alloc !*FD Set if we need to deallocate the arrays first

    ! First, deallocate if necessary

    if (present(re_alloc)) then
      if (re_alloc) then
        if (associated(proj%costheta))  deallocate(proj%costheta)
        if (associated(proj%sintheta))  deallocate(proj%sintheta)
        if (associated(proj%latitudes)) deallocate(proj%latitudes)
      endif
    endif

    allocate(proj%costheta(proj%nx,proj%ny),proj%sintheta(proj%nx,proj%ny))
    allocate(proj%latitudes(proj%nx,proj%ny))
   
  end subroutine proj_allocate

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine ll_to_xy(lon,lat,x,y,proj)

    !*FD Convert lat-long coordinates to grid coordinates. Note that
    !*FD The subroutine returns the x-y coordinates as real values,
    !*FD non-integer values indicating a position between grid-points.

    use gmt

    implicit none

    real(rk),intent(in) :: lat  !*FD The location of the point in lat-lon space (Latitude)
    real(rk),intent(in) :: lon  !*FD The location of the point in lat-lon space (Longitude)
    real(rk),intent(out) :: x   !*FD The location of the point in $x$-$y$ space ($x$ coordinate)
    real(rk),intent(out) :: y   !*FD The location of the point in $x$-$y$ space ($y$ coordinate)
    type(projection),intent(in) :: proj !*FD The projection being used

    real(rk) :: tlat,tlon

    tlat=lat ; tlon=lon

    select case(proj%p_type)
    case(1)
      call gmt_lambeq(tlon,tlat,x,y,proj%gmt_params)
    case(2)
      call gmt_plrs_sph(tlon,tlat,x,y,proj%gmt_params)
    case(3)
      call gmt_stereo1_sph (tlon,tlat,x,y,proj%gmt_params)
    case(4)
      call gmt_stereo2_sph (tlon,tlat,x,y,proj%gmt_params)
    case default
      print*,'ERROR: ',proj%p_type,' is not a valid projection type'
    end select

    x=(x/proj%dx)+proj%cpx
    y=(y/proj%dy)+proj%cpy

  end subroutine ll_to_xy

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine xy_to_ll(lon,lat,x,y,proj)

    !*FD Convert grid coordinates to lat-lon coordinates. The 
    !*FD subroutine returns the lat-lon coordinates as real values,
    !*FD non-integer values indicating a position between grid-points.

    use gmt

    implicit none

    real(rk),intent(in) :: x     !*FD The location of the point in $x$-$y$ space ($x$ coordinate)
    real(rk),intent(in) :: y     !*FD The location of the point in $x$-$y$ space ($y$ coordinate)
    real(rk),intent(out) :: lat  !*FD The location of the point in lat-lon space (Latitude)
    real(rk),intent(out) :: lon  !*FD The location of the point in lat-lon space (Longitude)
    type (projection),intent(in) :: proj !*FD The projection being used

    real(rk) :: xx,yy

    xx=x ; yy=y

    xx=(xx-proj%cpx)*proj%dx
    yy=(yy-proj%cpy)*proj%dy

    select case(proj%p_type)
    case(1)
      call gmt_ilambeq(lon,lat,xx,yy,proj%gmt_params)
    case(2)
      call gmt_iplrs_sph(lon,lat,xx,yy,proj%gmt_params)
    case(3:4)
      call gmt_istereo_sph (lon,lat,xx,yy,proj%gmt_params)
    case default
      print*,'ERROR: ',proj%p_type,' is not a valid projection type'
    end select

    lon=loncorrect(lon)

  end subroutine xy_to_ll

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine calc_lats(proj,array)

    !*FD Calculates the latitude of all points in the
    !*FD projected domain. The output latitudes (\texttt{array})
    !*FD are in degrees north. Note that the array must be the right 
    !*FD size, otherwise an error is generated, and the program halted.

    type(projection),intent(in) :: proj !*FD The projection to be used 
    real(sp),dimension(:,:),intent(out) :: array !*FD The array for output

    integer :: i,j,nx,ny
    real(rk) :: lon,lat,x,y

    nx=size(array,1) ; ny=size(array,2)

    if ((nx.ne.proj%nx).or.(ny.ne.proj%ny)) then
      print*,'Array size wrong in calc_lats'
      stop
    endif

    do i=1,nx
      do j=1,ny
        x=i ; y=j
        call xy_to_ll(lon,lat,x,y,proj)
        array(i,j)=lat
      enddo
    enddo

  end subroutine calc_lats

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine calc_grid_angle(proj)

    !*FD Calculates the angle the projected 
    !*FD grid makes with north at each point and stores the cos 
    !*FD and sin of that angle in the relevant arrays in \texttt{proj}.

    type(projection),intent(inout) :: proj !*FD The projection to be used
 
    integer :: i,j
    real(rk) :: latn,lonn,lats,lons,lat,lon,dlat,dlon,temp

    do i=1,proj%nx

    ! Main, central block

      do j=2,proj%ny-1
        call xy_to_ll(lonn,latn,real(i,rk),real(j+1,rk),proj)
        call xy_to_ll(lon,lat,real(i,rk),real(j,rk),proj)
        call xy_to_ll(lons,lats,real(i,rk),real(j-1,rk),proj)
        dlat=latn-lats
        dlon=lonn-lons
        if (dlon<-90) dlon=dlon+360
        temp=atan(dlon/dlat)
        proj%sintheta(i,j)=sin(temp)
        proj%costheta(i,j)=cos(temp)
        proj%latitudes(i,j)=lat
      enddo

    ! bottom row

      call xy_to_ll(lonn,latn,real(i,rk),real(2,rk),proj)
      call xy_to_ll(lon,lat,real(i,rk),real(1,rk),proj)
      dlat=latn-lat
      dlon=lonn-lon
      if (dlon<-90) dlon=dlon+360
      temp=atan(dlon/dlat)
      proj%sintheta(i,1)=sin(temp)
      proj%costheta(i,1)=cos(temp)
      proj%latitudes(i,1)=lat

    ! top row

      call xy_to_ll(lon,lat,real(i,rk),real(proj%ny,rk),proj)
      call xy_to_ll(lons,lats,real(i,rk),real(proj%ny-1,rk),proj)
      dlat=lat-lats
      dlon=lon-lons
      if (dlon<-90) dlon=dlon+360
      temp=atan(dlon/dlat)
      proj%sintheta(i,proj%ny)=sin(temp)
      proj%costheta(i,proj%ny)=cos(temp)
      proj%latitudes(i,proj%ny)=lat

    enddo

  end subroutine calc_grid_angle

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine print_corners(proj)

    !*FD A diagnostic routine that prints the corners of the 
    !*FD supplied projection in latitude-longitude coordinates 
    !*FD to standard out. Useful for error finding.

    type(projection),intent(in) :: proj !*FD Projection for printing.

    real(rk) :: clat,clon

    call xy_to_ll(clon,clat,real(1,rk),real(1,rk),proj)
    print*,'SW: lon=',clon,' lat=',clat

    call xy_to_ll(clon,clat,real(proj%nx,rk),real(1,rk),proj)
    print*,'SE: lon=',clon,' lat=',clat

    call xy_to_ll(clon,clat,real(1,rk),real(proj%ny,rk),proj)
    print*,'NW: lon=',clon,' lat=',clat

    call xy_to_ll(clon,clat,real(proj%nx,rk),real(proj%ny,rk),proj)
    print*,'NE: lon=',clon,' lat=',clat

  end subroutine print_corners

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine proj_write_restart(proj,unit)

    !*FD Write the part of the
    !*FD restart file relating to the {\tt projection}.
    !*FD Note that the file unit specified must be already
    !*FD open to an unformatted file.

    type(projection),intent(in) :: proj !*FD Projection type to be written
    integer,intent(in) :: unit          !*FD Logical file unit to use

    write(unit) proj%p_type
    write(unit) proj%nx
    write(unit) proj%ny
    write(unit) proj%dx
    write(unit) proj%dy
    write(unit) proj%cpx
    write(unit) proj%cpy
    write(unit) proj%latc
    write(unit) proj%lonc
    write(unit) proj%sintheta
    write(unit) proj%costheta
    write(unit) proj%latitudes

    ! GMT parameters

    write(unit) proj%gmt_params% cosp
    write(unit) proj%gmt_params% sinp
    write(unit) proj%gmt_params% pole
    write(unit) proj%gmt_params% central_meridian
    write(unit) proj%gmt_params% EQ_RAD
    write(unit) proj%gmt_params% i_EQ_RAD
    write(unit) proj%gmt_params% Dx
    write(unit) proj%gmt_params% Dy
    write(unit) proj%gmt_params% iDx
    write(unit) proj%gmt_params% iDy
    write(unit) proj%gmt_params% s_c
    write(unit) proj%gmt_params% s_ic
    write(unit) proj%gmt_params% n_polar
    write(unit) proj%gmt_params% s_polar
    write(unit) proj%gmt_params% north_pole
    write(unit) proj%gmt_params% polar

  end subroutine proj_write_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine proj_read_restart(proj,unit)

    !*FD Read the part of the
    !*FD restart file relating to the projection
    !*FD Note that the file unit specified must be already
    !*FD open to an unformatted file.

    type(projection),intent(inout) :: proj !*FD Projection to be read
    integer,intent(in) :: unit             !*FD Logical file unit to use

    read(unit) proj%p_type
    read(unit) proj%nx
    read(unit) proj%ny
    read(unit) proj%dx
    read(unit) proj%dy
    read(unit) proj%cpx
    read(unit) proj%cpy
    read(unit) proj%latc
    read(unit) proj%lonc

    call proj_allocate(proj)

    read(unit) proj%sintheta
    read(unit) proj%costheta
    read(unit) proj%latitudes

    ! GMT parameters

    read(unit) proj%gmt_params% cosp
    read(unit) proj%gmt_params% sinp
    read(unit) proj%gmt_params% pole
    read(unit) proj%gmt_params% central_meridian
    read(unit) proj%gmt_params% EQ_RAD
    read(unit) proj%gmt_params% i_EQ_RAD
    read(unit) proj%gmt_params% Dx
    read(unit) proj%gmt_params% Dy
    read(unit) proj%gmt_params% iDx
    read(unit) proj%gmt_params% iDy
    read(unit) proj%gmt_params% s_c
    read(unit) proj%gmt_params% s_ic
    read(unit) proj%gmt_params% n_polar
    read(unit) proj%gmt_params% s_polar
    read(unit) proj%gmt_params% north_pole
    read(unit) proj%gmt_params% polar

  end subroutine proj_read_restart

end module glimmer_project
