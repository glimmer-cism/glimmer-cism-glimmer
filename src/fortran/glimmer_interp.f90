! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_interp.f90 - part of the GLIMMER ice model       + 
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

module glimmer_interp

  !*FD Downscaling and upscaling
  !*FD routines for use in GLIMMER.

  use glimmer_global
  use glimmer_project

  implicit none

  type downscale

    !*FD Derived type containing indexing 
    !*FD information for downscaling. This type was 
    !*FD included for speed. Four of the arrays contained in it
    !*FD are arrays of the indices of the corners 
    !*FD of the global grid-boxes within which the given
    !*FD local grid point lies.

    integer, dimension(:,:),pointer :: il    !*FD (see above)
    integer, dimension(:,:),pointer :: ilp   !*FD (see above)
    integer, dimension(:,:),pointer :: jl    !*FD (see above)
    integer, dimension(:,:),pointer :: jlp   !*FD (see above)
    real(rk),dimension(:,:),pointer :: llats !*FD The latitude of each point in x-y space.
    real(rk),dimension(:,:),pointer :: llons !*FD The longitude of each point in x-y space.

  end type downscale

contains

  subroutine new_downscale(downs,proj,lons,lats)

    !*FD Initialises a downscale variable,
    !*FD according to given projected and global grids

    ! Arguments

    type(downscale),intent(out)      :: downs   !*FD Downscaling variable to be set
    type(projection),intent(in)      :: proj    !*FD Projection to use
    real(rk),intent(in),dimension(:) :: lats    !*FD Latitudes of global grid points (degrees).
                                                !*FD As in all cases, latitudes begin at north pole.
    real(rk),intent(in),dimension(:) :: lons    !*FD Longitude of global grid points (degrees)

    ! Internal variables

    real(rk) :: llat,llon
    integer :: i,j

    ! Allocate arrays

    allocate(downs%il   (proj%nx,proj%ny))
    allocate(downs%ilp  (proj%nx,proj%ny))
    allocate(downs%jl   (proj%nx,proj%ny))
    allocate(downs%jlp  (proj%nx,proj%ny))
    allocate(downs%llons(proj%nx,proj%ny))
    allocate(downs%llats(proj%nx,proj%ny))

    ! index local boxes

    call index_local_boxes(downs%il,downs%ilp,downs%jl,downs%jlp,lons,lats,proj)

    ! Find lats and lons

    do i=1,proj%nx
      do j=1,proj%ny
        call xy_to_ll(llon,llat,real(i,rk),real(j,rk),proj)
        downs%llons(i,j)=llon
        downs%llats(i,j)=llat
      end do
    end do

  end subroutine new_downscale

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine interp_wind_to_local(proj,zonwind,merwind,lats,lons,downs,xwind,ywind)

    !*FD Interpolates a global wind field 
    !*FD (or any vector field) onto a given projected grid.

    use glimmer_utils

    ! Argument declarations

    type(projection),       intent(in)  :: proj             !*FD Target map projection
    real(rk),dimension(:,:),intent(in)  :: zonwind          !*FD Zonal component (input)
    real(rk),dimension(:,:),intent(in)  :: merwind          !*FD Meridional components (input)
    real(rk),dimension(:),  intent(in)  :: lats             !*FD Latitudes of global gridpoints 
    real(rk),dimension(:),  intent(in)  :: lons             !*FD Longitudes of global gridpoints 
    type(downscale),        intent(in)  :: downs            !*FD Downscaling parameters
    real(rk),dimension(:,:),intent(out) :: xwind,ywind      !*FD x and y components on the projected grid (output)

    ! Declare two temporary arrays to hold the interpolated zonal and meridional winds

    real(dp),dimension(size(xwind,1),size(xwind,2)) :: tempzw,tempmw

    ! Check input arrays are conformal to one another

    call check_conformal(zonwind,merwind,'interp_wind 1')
    call check_conformal(xwind,ywind,'interp_wind 2')

    ! Interpolate onto the projected grid

    call interp_to_local(proj,zonwind,lats,lons,downs,localdp=tempzw)
    call interp_to_local(proj,merwind,lats,lons,downs,localdp=tempmw)

    ! Apply rotation

    xwind=tempzw*proj%costheta-tempmw*proj%sintheta
    ywind=tempzw*proj%sintheta+tempmw*proj%costheta

  end subroutine interp_wind_to_local

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine interp_to_local(proj,global,lats,lons,downs,localsp,localdp)

    !*FD Interpolate a global scalar field
    !*FD onto a projected grid. 
    !*FD 
    !*FD This uses a simple bilinear interpolation, which assumes
    !*FD that the global grid boxes are rectangular - i.e. it works
    !*FD in lat-lon space.
    !*FD
    !*FD Either localsp or localdp must be present (or both), depending
    !*FD which precision output is required.
  
    use glimmer_utils

    ! Argument declarations

    type(projection),        intent(in)           :: proj      !*FD Target map projection
    real(rk), dimension(:,:),intent(in)           :: global    !*FD Global field (input)
    real(rk), dimension(:),  intent(in)           :: lats      !*FD Latitudes of global gridpoints 
    real(rk), dimension(:),  intent(in)           :: lons      !*FD Longitudes of global gridpoints 
    type(downscale),         intent(in)           :: downs     !*FD Downscaling parameters
    real(sp),dimension(:,:),intent(out),optional :: localsp   !*FD Local field on projected grid (output) sp
    real(dp),dimension(:,:),intent(out),optional :: localdp   !*FD Local field on projected grid (output) dp

    ! Local variable declarations

    integer  :: nxg,nyg                      ! Dimensions of global fields
    integer  :: i,j                          ! Counter variables for main loop
    real(rk) :: dx,dy,x,y,ilon,jlat          ! Intermediate variables in interpolation formula
    real(rk),dimension(2,2) :: f             ! Temporary array holding the four points in the 
                                             ! interpolation domain.

    ! Retrieve the dimensions of the global domain

    nxg=size(global,1) ; nyg=size(global,2)

    ! Check to see that we have the right number of lats and lons

    if ((nxg/=size(lons)).or.(nyg/=size(lats))) then
      print*,'size mismatch in interp_to_local'
      stop
    endif

    ! check we have one output at least...

    if (.not.(present(localsp).or.present(localdp))) then
      print*, 'WARNING interp_to_local has no output'
    endif

    ! Main interpolation loop

    do i=1,proj%nx
      do j=1,proj%ny

        ! Find out where point i,j is in lat-lon space

        ilon=downs%llons(i,j)
        jlat=downs%llats(i,j)

        ! Compile the temporary array f from adjacent points 

        f(1,1)=global(downs%il(i,j),downs%jl(i,j))
        f(1,2)=global(downs%il(i,j),downs%jlp(i,j))
        f(2,1)=global(downs%ilp(i,j),downs%jl(i,j))
        f(2,2)=global(downs%ilp(i,j),downs%jlp(i,j))

        ! Calculate dx of interpolation domain and correct for boundary effects

        dx=lons(downs%ilp(i,j))-lons(downs%il(i,j))
        dx=loncorrect(dx)

        ! Calculate dy of interpolation domain

        dy=lats(downs%jl(i,j))-lats(downs%jlp(i,j))

        ! Calculate location of interpolation point within the rectangle

        x=ilon-lons(downs%il(i,j))
        x=loncorrect(x)
        y=lats(downs%jl(i,j))-jlat

        ! Apply the bilinear interpolation

        if (present(localsp)) localsp(i,j)=bilinear_interp(x,y,f,dx,dy)
        if (present(localdp)) localdp(i,j)=bilinear_interp(x,y,f,dx,dy)

      enddo
    enddo

  end subroutine interp_to_local

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine pointwise_to_global(proj,local,lons,lats,global)

    !*FD Upscale to global domain by
    !*FD pointwise sampling.
    !*FD
    !*FD Note that this is the mathematically inverse process of the 
    !*FD \texttt{interp\_to\_local} routine.

    ! Arguments

    type(projection),       intent(in)  :: proj      !*FD Projection to use
    real(rk),dimension(:,:),intent(in)  :: local     !*FD Local field (input)
    real(rk),dimension(:,:),intent(out) :: global    !*FD Global field (output)
    real(rk),dimension(:),  intent(in)  :: lats      !*FD Latitudes of grid-points (degrees)
    real(rk),dimension(:),  intent(in)  :: lons      !*FD Longitudes of grid-points (degrees)

    ! Internal variables

    real(rk),dimension(2,2) :: f
    integer :: nxg,nyg,nxl,nyl,i,j,xx,yy
    real(rk) :: x,y

    nxg=size(global,1) ; nyg=size(global,2)
    nxl=size(local,1)  ; nyl=size(local,2)

    do i=1,nxg
      do j=1,nyg
        call ll_to_xy(lons(i),lats(j),x,y,proj)
        xx=int(x) ; yy=int(y)
        if (nint(x)<=1.or.nint(x)>nxl-1.or.nint(y)<=1.or.nint(y)>nyl-1) then
          global(i,j)=0.0
        else
          f=local(xx:xx+1,yy:yy+1)
          global(i,j)=bilinear_interp(x-real(xx),y-real(yy),f,real(1.0,rk),real(1.0,rk))
        endif
      enddo
    enddo  

  end subroutine pointwise_to_global

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine mean_to_global(proj,local,gboxx,gboxy,gboxn,global)

    !*FD Upscale to global domain by
    !*FD areal averaging.
    !*FD
    !*FD Note that:
    !*FD \begin{itemize}
    !*FD \item \texttt{gboxx} and \texttt{gboxy} are the same size as \texttt{local}
    !*FD \item \texttt{gboxn} is the same size as \texttt{global}
    !*FD \item This method is \emph{not} the mathematical inverse of the
    !*FD \texttt{interp\_to\_local} routine.
    !*FD \end{itemize}

    ! Arguments

    type(projection),       intent(in)  :: proj   !*FD Projection of local grid.
    real(rk),dimension(:,:),intent(in)  :: local  !*FD Data on projected grid (input).
    integer, dimension(:,:),intent(in)  :: gboxx  !*FD $x$-indicies of global grid-box 
                                                  !*FD containing given local grid-box.
    integer, dimension(:,:),intent(in)  :: gboxy  !*FD $y$-indicies of global grid-box 
                                                  !*FD containing given local grid-box.
    integer, dimension(:,:),intent(in)  :: gboxn  !*FD Number of local grid-boxes contained in each global box
    real(rk),dimension(:,:),intent(out) :: global !*FD Data on global grid (output).

    ! Internal variables

    integer :: nxl,nyl,i,j

    ! Beginning of code

    nxl=size(local,1) ; nyl=size(local,2)

    global=0.0

    do i=1,nxl
      do j=1,nyl
        global(gboxx(i,j),gboxy(i,j))=global(gboxx(i,j),gboxy(i,j))+local(i,j)
       enddo
    enddo  

    where (gboxn.ne.0)
      global=global/gboxn
    elsewhere
      global=0.0
    endwhere  

  end subroutine mean_to_global

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(rk) function bilinear_interp(x,y,f,dx,dy)

    !*FD Performs bilinear interpolation 
    !*FD in a rectangular domain. Note that the bilinear interpolation formula is:
    !*FD  \[f_{\mathtt{x},\mathtt{y}}=(1-t)(1-u)f_{1,1}+t(1-u)f_{2,1}+tuf_{2,2}+(1-t)uf_{1,2}\]
    !*FD with
    !*FD \[t=\frac{\mathtt{x}}{\mathtt{dx}}\] \[u=\frac{\mathtt{y}}{\mathtt{dy}}.\]
    !*FDRV The value of \texttt{f} at \texttt{x,y}

    ! Argument declarations

    real(rk),dimension(2,2),intent(in) :: f     !*FD The interpolation domain;
                                                !*FD i.e. the four points surrounding the
                                                !*FD target.
    real(rk),               intent(in) :: dx    !*FD $x$-dimension of the rectangle contained in \texttt{f}.
    real(rk),               intent(in) :: dy    !*FD $y$-dimension of the rectangle contained in \texttt{f}.
    real(rk),               intent(in) :: x     !*FD The $x$-location of the target, offset from
                                                !*FD the location of \texttt{f(1,1)}. Note: \emph{not}
                                                !*FD the fractional displacement, but the \emph{actual}
                                                !*FD displacement. The fractional displacement is 
                                                !*FD in the function code.
    real(rk),               intent(in) :: y     !*FD The $y$-location of the target, offset from
                                                !*FD the location of \texttt{f(1,1)}. Note: \emph{not}
                                                !*FD the fractional displacement, but the \emph{actual}
                                                !*FD displacement. The fractional displacement is 
                                                !*FD in the function code.
    ! Local variables

    real(rk) :: t,u                         ! Fractional displacement of target point

    ! Calculate fractional displacement of target

    t=x/dx
    u=y/dy

    ! Apply bilinear interpolation formula

    bilinear_interp=(1-t)*(1-u)*f(1,1)+t*(1-u)*f(2,1)+t*u*f(2,2)+(1-t)*u*f(1,2)

  end function bilinear_interp

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine find_ll_index(il,jl,lon,lat,lons,lats)

    !*FD Find the global gridpoint at the first corner of the box surrounding
    !*FD a given location in lat-lon space.

    use glimmer_utils

    ! Arguments

    real(rk),             intent(in)  :: lon    !*FD Longitude of location to be indexed (input)
    real(rk),             intent(in)  :: lat    !*FD Latitude of location to be indexed (input)
    real(rk),dimension(:),intent(in)  :: lats   !*FD Latitudes of global grid points 
    real(rk),dimension(:),intent(in)  :: lons   !*FD Longitudes of global grid points 
    integer,              intent(out) :: il     !*FD $x$-gridpoint index (output)
    integer,              intent(out) :: jl     !*FD $y$-gridpoint index (output)
  
    ! Internal variables

    integer :: nx,ny

    nx=size(lons) ; ny=size(lats)

    if ((lon>lons(nx)).and.(lon<360.0)) then
      il=nx
    else
      il=1
      do while (lon>=array_bcs(lons,il))
        il=il+1
      enddo
      il=il-1
    endif

    if ((lat<lats(ny)).and.(lat>-90.0)) then
      jl=ny
    else
      jl=1
      do while (lat<=array_bcs_lats(lats,jl))
        jl=jl+1
      enddo
      jl=jl-1
    endif

  end subroutine find_ll_index

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine index_local_boxes(il,ilp,jl,jlp,lons,lats,proj)

    !*FD Indexes the corners of the
    !*FD global grid box in which each local grid box sits.

    use glimmer_utils

    ! Arguments

    integer, dimension(:,:),intent(out) :: il,ilp,jl,jlp !*FD Array of indicies (see \texttt{downscale} type)
    real(rk),dimension(:),  intent(in)  :: lons          !*FD Longitudinal positions of global grid points (degrees)
    real(rk),dimension(:),  intent(in)  :: lats          !*FD Latitudinal positions of global grid points (degrees)
    type(projection),       intent(in)  :: proj          !*FD Projection to be used

    ! Internal variables

    integer :: i,j,nxg,nyg
    real(rk) :: ilon,jlat

    nxg=size(lons) ; nyg=size(lats)

    do i=1,proj%nx
      do j=1,proj%ny

        ! Find out where point i,j is in lat-lon space

        call xy_to_ll(ilon,jlat,real(i,rk),real(j,rk),proj)

        ! Index that location onto the global grid

        call find_ll_index(il(i,j),jl(i,j),ilon,jlat,lons,lats)

        ! Compile the temporary array f from adjacent points 

        ilp(i,j)=il(i,j)+1
        jlp(i,j)=jl(i,j)+1

        call fix_bcs2d(il(i,j) ,jl(i,j) ,nxg,nyg)
        call fix_bcs2d(ilp(i,j),jlp(i,j),nxg,nyg)

      enddo
    enddo

  end subroutine index_local_boxes

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine downs_write_restart(downs,unit)

    !*FD Write the part of a
    !*FD restart file relating to the downscaling.
    !*FD Note that the logical file unit must already be open

    ! Arguments

    type(downscale),intent(in) :: downs  !*FD The downscaling parameters to be written
    integer,        intent(in) :: unit   !*FD The logical file unit to use

    ! temporary variable, this is need to fix an internal compiler error of the SUN WS f95 compiler.
    logical :: temp

    ! Beginning of code
    
    temp = associated(downs%il)
    write(unit) temp
    if (temp)  write(unit) downs%il

    temp = associated(downs%ilp)
    write(unit) temp
    if (temp)  write(unit) downs%il

    temp = associated(downs%jl)
    write(unit) temp
    if (temp)  write(unit) downs%il

    temp = associated(downs%jlp)
    write(unit) temp
    if (temp)  write(unit) downs%il

  end subroutine downs_write_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine downs_read_restart(downs,unit)

    !*FD Read the part of a
    !*FD restart file relating to the downscaling.
    !*FD Note that the logical file unit must already be open.

    type(downscale),intent(inout) :: downs !*FD The downscaling parameters to be read.
    integer,        intent(in)    :: unit  !*FD The logical file unit to use.

    logical :: tempflag

    read(unit) tempflag
    if(tempflag) read(unit) downs%il
    read(unit) tempflag
    if(tempflag) read(unit) downs%ilp
    read(unit) tempflag
    if(tempflag) read(unit) downs%jl
    read(unit) tempflag
    if(tempflag) read(unit) downs%jlp

  end subroutine downs_read_restart

end module glimmer_interp
