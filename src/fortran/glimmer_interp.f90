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

    real(rk),dimension(:,:),pointer :: llats => null() !*FD The latitude of each point in x-y space.
    real(rk),dimension(:,:),pointer :: llons => null() !*FD The longitude of each point in x-y space.

    integer, dimension(:,:,:),pointer :: xloc => null() !*FD The x-locations of the corner points of the
                                                        !*FD interpolation domain.
    integer, dimension(:,:,:),pointer :: yloc => null() !*FD The y-locations of the corner points of the
                                                        !*FD interpolation domain.
    real(rk),dimension(:,:),  pointer :: xfrac => null()
    real(rk),dimension(:,:),  pointer :: yfrac => null()

  end type downscale

  type upscale
  
    !*FD Derived type containing indexing information
    !*FD for upscaling by areal averaging.

    integer, dimension(:,:),pointer :: gboxx => null() !*FD $x$-indicies of global grid-box 
                                                       !*FD containing given local grid-box.
    integer, dimension(:,:),pointer :: gboxy => null() !*FD $y$-indicies of global grid-box 
                                                       !*FD containing given local grid-box.
    integer, dimension(:,:),pointer :: gboxn => null() !*FD Number of local grid-boxes 
                                                       !*FD contained in each global box.
    logical                         :: set = .false.   !*FD Set if the type has been initialised.
  end type upscale

  interface mean_to_global
    module procedure mean_to_global_sp,mean_to_global_dp
  end interface

contains

  subroutine new_downscale(downs,proj,grid)

    use glimmer_global_grid

    !*FD Initialises a downscale variable,
    !*FD according to given projected and global grids

    ! Arguments

    type(downscale),intent(out)      :: downs   !*FD Downscaling variable to be set
    type(projection),intent(in)      :: proj    !*FD Projection to use
    type(global_grid),intent(in)     :: grid    !*FD Global grid to use

    ! Internal variables

    real(rk) :: llat,llon
    integer :: i,j

    ! Allocate arrays

    allocate(downs%xloc (proj%nx,proj%ny,4))
    allocate(downs%yloc (proj%nx,proj%ny,4))
    allocate(downs%xfrac(proj%nx,proj%ny))
    allocate(downs%yfrac(proj%nx,proj%ny))
    allocate(downs%llons(proj%nx,proj%ny))
    allocate(downs%llats(proj%nx,proj%ny))

    ! index local boxes

    call index_local_boxes(downs%xloc,downs%yloc,downs%xfrac,downs%yfrac,grid,proj)

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

  subroutine interp_wind_to_local(proj,zonwind,merwind,downs,xwind,ywind)

    !*FD Interpolates a global wind field 
    !*FD (or any vector field) onto a given projected grid.

    use glimmer_utils

    ! Argument declarations

    type(projection),       intent(in)  :: proj             !*FD Target map projection
    real(rk),dimension(:,:),intent(in)  :: zonwind          !*FD Zonal component (input)
    real(rk),dimension(:,:),intent(in)  :: merwind          !*FD Meridional components (input)
    type(downscale),        intent(in)  :: downs            !*FD Downscaling parameters
    real(rk),dimension(:,:),intent(out) :: xwind,ywind      !*FD x and y components on the projected grid (output)

    ! Declare two temporary arrays to hold the interpolated zonal and meridional winds

    real(dp),dimension(size(xwind,1),size(xwind,2)) :: tempzw,tempmw

    ! Check input arrays are conformal to one another

    call check_conformal(zonwind,merwind,'interp_wind 1')
    call check_conformal(xwind,ywind,'interp_wind 2')

    ! Interpolate onto the projected grid

    call interp_to_local(proj,zonwind,downs,localdp=tempzw)
    call interp_to_local(proj,merwind,downs,localdp=tempmw)

    ! Apply rotation

    xwind=tempzw*proj%costheta-tempmw*proj%sintheta
    ywind=tempzw*proj%sintheta+tempmw*proj%costheta

  end subroutine interp_wind_to_local

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine interp_to_local(proj,global,downs,localsp,localdp,global_fn)

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
    type(downscale),         intent(in)           :: downs     !*FD Downscaling parameters
    real(sp),dimension(:,:), intent(out),optional :: localsp   !*FD Local field on projected grid (output) sp
    real(dp),dimension(:,:), intent(out),optional :: localdp   !*FD Local field on projected grid (output) dp
    real(sp),optional :: global_fn                             !*FD Function returning values in global field. This  
                                                               !*FD may be used as an alternative to passing the
                                                               !*FD whole array in \texttt{global} if, for instance the
                                                               !*FD data-set is in a large file, being accessed point by point.
                                                               !*FD In these circumstances, \texttt{global}
                                                               !*FD may be of any size, and its contents are irrelevant.

    ! Local variable declarations

    integer  :: i,j                          ! Counter variables for main loop
    real(rk),dimension(4) :: f               ! Temporary array holding the four points in the 
                                             ! interpolation domain.

    ! check we have one output at least...

    if (.not.(present(localsp).or.present(localdp))) then
      print*, 'WARNING interp_to_local has no output'
    endif

    ! Main interpolation loop

    do i=1,proj%nx
      do j=1,proj%ny

        ! Compile the temporary array f from adjacent points 

        if (present(global_fn)) then
          f(1)=global_fn(downs%xloc(i,j,1),downs%yloc(i,j,1))
          f(2)=global_fn(downs%xloc(i,j,2),downs%yloc(i,j,2))
          f(3)=global_fn(downs%xloc(i,j,3),downs%yloc(i,j,3))
          f(4)=global_fn(downs%xloc(i,j,4),downs%yloc(i,j,4))
        else
          f(1)=global(downs%xloc(i,j,1),downs%yloc(i,j,1))
          f(2)=global(downs%xloc(i,j,2),downs%yloc(i,j,2))
          f(3)=global(downs%xloc(i,j,3),downs%yloc(i,j,3))
          f(4)=global(downs%xloc(i,j,4),downs%yloc(i,j,4))
        end if

         ! Apply the bilinear interpolation

        if (present(localsp)) localsp(i,j)=bilinear_interp(downs%xfrac(i,j),downs%yfrac(i,j),f)
        if (present(localdp)) localdp(i,j)=bilinear_interp(downs%xfrac(i,j),downs%yfrac(i,j),f)

      enddo
    enddo

  end subroutine interp_to_local

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine mean_to_local(proj,global,lats,lons,localsp,localdp,global_fn)

    !*FD Average a high-resolution global field onto the projected grid
    !*FD This assumes that the global field is sufficiently high-resolution 
    !*FD compared with the local grid - it just averages the points contained 
    !*FD in each local grid-box.
 
    use glimmer_project
    use glimmer_utils

    ! Argument declarations

    type(projection),                intent(in)  :: proj      !*FD Target map projection
    real(rk),dimension(:,:),         intent(in)  :: global    !*FD Global field (input)
    real(rk),dimension(:),           intent(in)  :: lats      !*FD Latitudes of global gridpoints 
    real(rk),dimension(:),           intent(in)  :: lons      !*FD Longitudes of global gridpoints 
    real(sp),dimension(:,:),optional,intent(out) :: localsp   !*FD Local field on projected grid (output) sp
    real(dp),dimension(:,:),optional,intent(out) :: localdp   !*FD Local field on projected grid (output) dp
    real(sp),optional                            :: global_fn !*FD Function returning values in global field. This  
                                                              !*FD may be used as an alternative to passing the
                                                              !*FD whole array in \texttt{global} if, for instance the
                                                              !*FD data-set is in a large file, being accessed point by point.
                                                              !*FD In these circumstances, \texttt{global}
                                                              !*FD may be of any size, and its contents are irrelevant.

    integer :: nxg,nyg,i,j,xbox,ybox
    real(rk) :: lat,lon,x,y
    real(dp),dimension(proj%nx,proj%ny) :: temp_out
    integer,dimension(proj%nx,proj%ny) :: mean_count

    ! Find size of input field

    nxg=size(lons) ; nyg=size(lats)

    if (.not.present(global_fn)) then 
       if ((nxg/=size(lons)).or.(nyg/=size(lats))) then
          print*,'size mismatch in interp_to_local'
          stop
       end if
    end if

    ! check we have one output at least...

    if (.not.(present(localsp).or.present(localdp))) then
      print*, 'WARNING mean_to_local has no output'
    endif

    ! Zero some things

    mean_count=0
    temp_out=0.0

    ! Loop over all global points

    do i=1,nxg

       lon=lons(i)

       do j=1,nyg

          ! Find location in local coordinates

          lat=lats(j)  ! (Have already found lat above)
          call ll_to_xy(lon,lat,x,y,proj)
          xbox=nint(x)
          ybox=nint(y)
          
          ! Add to appropriate location and update count

          if (xbox.ge.1.and.xbox.le.proj%nx.and. &
              ybox.ge.1.and.ybox.le.proj%ny) then
             if (present(global_fn)) then
                temp_out(xbox,ybox)=temp_out(xbox,ybox)+global_fn(i,j)
             else
                temp_out(xbox,ybox)=temp_out(xbox,ybox)+global(i,j)
             end if
             mean_count(xbox,ybox)=mean_count(xbox,ybox)+1
          end if

       end do
    end do

    ! Divide by number of contributing points and copy to output

    if (present(localsp)) localsp=temp_out/real(mean_count,dp)
    if (present(localdp)) localdp=temp_out/real(mean_count,dp)
 
  end subroutine mean_to_local

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine pointwise_to_global(proj,local,lons,lats,global,mask)

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
    integer, dimension(:,:), intent(in),optional :: mask !*FD output mask for upscaling

    ! Internal variables

    real(rk),dimension(2,2) :: f
    integer :: nxg,nyg,nxl,nyl,i,j,xx,yy
    real(rk) :: x,y
    real(rk),dimension(size(local,1),size(local,2)) :: tempmask

    nxg=size(global,1) ; nyg=size(global,2)
    nxl=size(local,1)  ; nyl=size(local,2)

    do i=1,nxg
      do j=1,nyg
        call ll_to_xy(lons(i),lats(j),x,y,proj)
        xx=int(x) ; yy=int(y)
        if (nint(x)<=1.or.nint(x)>nxl-1.or.nint(y)<=1.or.nint(y)>nyl-1) then
          global(i,j)=0.0
        else
          f=local(xx:xx+1,yy:yy+1)*tempmask(xx:xx+1,yy:yy+1)
          global(i,j)=bilinear_interp((x-real(xx))/real(1.0,rk),(y-real(yy))/real(1.0,rk),f)
        endif
      enddo
    enddo  

  end subroutine pointwise_to_global

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine mean_to_global_sp(proj,ups,local,global,mask)

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
    type(upscale),          intent(in)  :: ups    !*FD Upscaling indexing data.
    real(sp),dimension(:,:),intent(in)  :: local  !*FD Data on projected grid (input).
    real(rk),dimension(:,:),intent(out) :: global !*FD Data on global grid (output).
    integer, dimension(:,:),intent(in),optional :: mask !*FD Mask for upscaling

    ! Internal variables

    integer :: nxl,nyl,i,j
    real(rk),dimension(size(local,1),size(local,2)) :: tempmask

    ! Beginning of code

    nxl=size(local,1) ; nyl=size(local,2)

    global=0.0

    if (present(mask)) then
      tempmask=mask
    else
      tempmask=1
    endif

    do i=1,nxl
      do j=1,nyl
        global(ups%gboxx(i,j),ups%gboxy(i,j))= &
               global(ups%gboxx(i,j),ups%gboxy(i,j))+local(i,j)*tempmask(i,j)
       enddo
    enddo  

    where (ups%gboxn.ne.0)
      global=global/ups%gboxn
    elsewhere
      global=0.0
    endwhere  

  end subroutine mean_to_global_sp

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine mean_to_global_dp(proj,ups,local,global,mask)

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
    type(upscale),          intent(in)  :: ups    !*FD Upscaling indexing data.
    real(dp),dimension(:,:),intent(in)  :: local  !*FD Data on projected grid (input).
    real(rk),dimension(:,:),intent(out) :: global !*FD Data on global grid (output).
    integer, dimension(:,:),intent(in),optional :: mask !*FD Mask for upscaling

    ! Internal variables

    integer :: nxl,nyl,i,j
    real(rk),dimension(size(local,1),size(local,2)) :: tempmask

    ! Beginning of code

    nxl=size(local,1) ; nyl=size(local,2)

    global=0.0

    if (present(mask)) then
      tempmask=mask
    else
      tempmask=1
    endif

    do i=1,nxl
      do j=1,nyl
        global(ups%gboxx(i,j),ups%gboxy(i,j))= &
               global(ups%gboxx(i,j),ups%gboxy(i,j))+local(i,j)*tempmask(i,j)
       enddo
    enddo  

    where (ups%gboxn.ne.0)
      global=global/ups%gboxn
    elsewhere
      global=0.0
    endwhere  

  end subroutine mean_to_global_dp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(rk) function bilinear_interp(xp,yp,f)

    !*FD Performs bilinear interpolation 
    !*FD in a rectangular domain. Note that the bilinear interpolation formula is:
    !*FD  \[f_{\mathtt{x},\mathtt{y}}=(1-X')(1-Y')f_{1}+X'(1-Y')f_{2}+X'Y'f_{3}+(1-X')Y'f_{4}\]
    !*FD where $X'$ and $Y'$ are the fractional displacements of the target point within the domain.
    !*RV The value of \texttt{f} at \texttt{x,y}

    ! Argument declarations

    real(rk),             intent(in) :: xp    !*FD The fractional $x$-displacement of the target.
    real(rk),             intent(in) :: yp    !*FD The fractional $y$-displacement of the target.
    real(rk),dimension(4),intent(in) :: f     !*FD The interpolation domain;
                                              !*FD i.e. the four points surrounding the
                                              !*FD target, presented anticlockwise from bottom-
                                              !*FD left 
    ! Apply bilinear interpolation formula

    bilinear_interp=(1-xp)*(1-yp)*f(1)+xp*(1-yp)*f(2)+xp*yp*f(3)+(1-xp)*yp*f(4)

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
    integer,dimension(1) :: loc

    nx=size(lons) ; ny=size(lats)

    if ((lon>maxval(lons)).and.(lon<360.0)) then
      loc=maxloc(lons)
      il=loc(1)
    else
      il=1
      do
        if (lon>=array_bcs(lons,il).and.lon<=array_bcs(lons,il+1)) exit
        il=il+1
      enddo
    endif

    if ((lat<lats(ny)).and.(lat>-90.0)) then
      jl=ny
      return
    endif

    if ((lat>lats(1)).and.(lat<90.0)) then
      jl=1
      return
    endif

    jl=1
    do 
      if (lat>lats(jl)) exit
      jl=jl+1
    enddo
    
  end subroutine find_ll_index

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine index_local_boxes(xloc,yloc,xfrac,yfrac,grid,proj)

    !*FD Indexes the corners of the
    !*FD global grid box in which each local grid box sits.

    use glimmer_utils
    use glimmer_global_grid

    ! Arguments

    integer, dimension(:,:,:),intent(out) :: xloc,yloc   !*FD Array of indicies (see \texttt{downscale} type)
    real(rk),dimension(:,:),  intent(out) :: xfrac,yfrac !*FD Fractional off-sets of grid points
    type(global_grid),        intent(in)  :: grid        !*FD Global grid to be used
    type(projection),         intent(in)  :: proj        !*FD Projection to be used

    ! Internal variables

    integer :: i,j,il,jl,temp
    real(rk) :: ilon,jlat,xa,ya,xb,yb,xc,yc,xd,yd

    do i=1,proj%nx
      do j=1,proj%ny

        ! Find out where point i,j is in lat-lon space

        call xy_to_ll(ilon,jlat,real(i,rk),real(j,rk),proj)

        ! Index that location onto the global grid

        call find_ll_index(il,jl,ilon,jlat,grid%lons,grid%lats)

        xloc(i,j,1)=il  ! This is the starting point - we now need to find
        yloc(i,j,1)=jl  ! three other points that enclose the interpolation target
        
 
        if (jlat>grid%lats(grid%ny)) then
          
          ! For all points except on the bottom row

          xloc(i,j,2)=il+1
          yloc(i,j,2)=jl

          xloc(i,j,3)=il+1
          yloc(i,j,3)=jl-1

          xloc(i,j,4)=il
          yloc(i,j,4)=jl-1

          call fix_bcs2d(xloc(i,j,2) ,yloc(i,j,2),grid%nx,grid%ny)
          call fix_bcs2d(xloc(i,j,3) ,yloc(i,j,3),grid%nx,grid%ny)
          call fix_bcs2d(xloc(i,j,4) ,yloc(i,j,4),grid%nx,grid%ny)

          if (jl==1) then
            temp=xloc(i,j,3)
            xloc(i,j,3)=xloc(i,j,4)
            xloc(i,j,4)=temp
          endif

        else

          ! The bottom row

          xloc(i,j,2)=il-1
          yloc(i,j,2)=jl

          xloc(i,j,3)=il-1
          yloc(i,j,3)=jl+1

          xloc(i,j,4)=il
          yloc(i,j,4)=jl+1
            
          call fix_bcs2d(xloc(i,j,2) ,yloc(i,j,2),grid%nx,grid%ny)
          call fix_bcs2d(xloc(i,j,3) ,yloc(i,j,3),grid%nx,grid%ny)
          call fix_bcs2d(xloc(i,j,4) ,yloc(i,j,4),grid%nx,grid%ny)

          temp=xloc(i,j,3)
          xloc(i,j,3)=xloc(i,j,4)
          xloc(i,j,4)=temp

        endif

        ! Now, find out where each of those points is on the projected
        ! grid, and calculate fractional displacements accordingly

        call ll_to_xy(grid%lons(xloc(i,j,1)),grid%lats(yloc(i,j,1)),xa,ya,proj)
        call ll_to_xy(grid%lons(xloc(i,j,2)),grid%lats(yloc(i,j,2)),xb,yb,proj)
        call ll_to_xy(grid%lons(xloc(i,j,3)),grid%lats(yloc(i,j,3)),xc,yc,proj)
        call ll_to_xy(grid%lons(xloc(i,j,4)),grid%lats(yloc(i,j,4)),xd,yd,proj)

        call calc_fractional(xfrac(i,j),yfrac(i,j),real(i,rk),real(j,rk), &
                             xa,ya,xb,yb,xc,yc,xd,yd)

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
    
!   temp = associated(downs%il)
!   write(unit) temp
!   if (temp)  write(unit) downs%il

!   temp = associated(downs%ilp)
!   write(unit) temp
!   if (temp)  write(unit) downs%il

!   temp = associated(downs%jl)
!   write(unit) temp
!   if (temp)  write(unit) downs%il

!   temp = associated(downs%jlp)
!   write(unit) temp
!   if (temp)  write(unit) downs%il

  end subroutine downs_write_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine downs_read_restart(downs,unit)

    !*FD Read the part of a
    !*FD restart file relating to the downscaling.
    !*FD Note that the logical file unit must already be open.

    type(downscale),intent(inout) :: downs !*FD The downscaling parameters to be read.
    integer,        intent(in)    :: unit  !*FD The logical file unit to use.

    logical :: tempflag

!    read(unit) tempflag
!    if(tempflag) read(unit) downs%il
!    read(unit) tempflag
!    if(tempflag) read(unit) downs%ilp
!    read(unit) tempflag
!    if(tempflag) read(unit) downs%jl
!    read(unit) tempflag
!    if(tempflag) read(unit) downs%jlp

  end subroutine downs_read_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine new_upscale(ups,grid,proj,mask)

    use glimmer_global_grid

    !*FD Compiles an index of which global grid box contains a given
    !*FD grid box on the projected grid, and sets derived type \texttt{ups}
    !*FD accordingly.

    ! Arguments

    type(upscale),         intent(out) :: ups        !*FD Upscaling type to be set
    type(global_grid),     intent(in)  :: grid       !*FD Global grid to be used
    type(projection),      intent(in)  :: proj       !*FD Projection being used
    integer,dimension(:,:),intent(in)  :: mask       !*FD Upscaling mask to be used
    
    ! Internal variables

    integer  :: i,j,ii,jj,nx,ny,gnx,gny
    real(rk) :: plon,plat
    
    ! Beginning of code

    if (associated(ups%gboxx)) deallocate(ups%gboxx)
    if (associated(ups%gboxy)) deallocate(ups%gboxy)
    if (associated(ups%gboxn)) deallocate(ups%gboxn)

    allocate(ups%gboxx(proj%nx,proj%ny))
    allocate(ups%gboxy(proj%nx,proj%ny))     
    allocate(ups%gboxn(grid%nx,grid%ny))

    gnx=grid%nx ; gny=grid%ny
    nx =proj%nx ; ny =proj%ny

    ups%gboxx=0 ; ups%gboxy=0

    do i=1,nx
      do j=1,ny
        call xy_to_ll(plon,plat,real(i,rk),real(j,rk),proj)
        ii=1 ; jj=1
        do
          ups%gboxx(i,j)=ii
          if (ii>gnx) then
            print*,'global index failure'
            stop
          endif  
          if (lon_between(grid%lon_bound(ii),grid%lon_bound(ii+1),plon)) exit
          ii=ii+1
        enddo

        jj=1

        do
          ups%gboxy(i,j)=jj
          if (jj>gny) then
            print*,'global index failure'
            stop
          endif  
          if ((grid%lat_bound(jj)>=plat).and.(plat>grid%lat_bound(jj+1))) exit
          jj=jj+1
        enddo

      enddo
    enddo

    ups%gboxn=0

    do i=1,nx
      do j=1,ny
        ups%gboxn(ups%gboxx(i,j),ups%gboxy(i,j))=ups%gboxn(ups%gboxx(i,j),ups%gboxy(i,j))+mask(i,j)
      enddo
    enddo

    ups%set=.true.

  end subroutine new_upscale

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine copy_upscale(in,out)

    type(upscale),intent(in)  :: in
    type(upscale),intent(out) :: out

    if (.not.in%set) then
      print*,'Attempt to copy un-initialised upscale type'
      stop
    endif

    if (associated(out%gboxx)) deallocate(out%gboxx)
    if (associated(out%gboxy)) deallocate(out%gboxy)
    if (associated(out%gboxn)) deallocate(out%gboxn)

    allocate(out%gboxx(size(in%gboxx,1),size(in%gboxx,2)))
    allocate(out%gboxy(size(in%gboxy,1),size(in%gboxy,2)))
    allocate(out%gboxn(size(in%gboxn,1),size(in%gboxn,2)))

    out%gboxx=in%gboxx
    out%gboxy=in%gboxy
    out%gboxn=in%gboxn

    out%set=.true.

  end subroutine copy_upscale

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  logical function lon_between(a,b,x)

    !*FD Checks to see whether a 
    !*FD longitudinal coordinate is between two bounds,
    !*FD taking into account the periodic boundary conditions.
    !*RV Returns \texttt{.true.} if $\mathtt{x}\geq \mathtt{a}$ and $\mathtt{x}<\mathtt{b}$.

    ! Arguments

    real(rk),intent(in) :: a  !*FD Lower bound on interval for checking
    real(rk),intent(in) :: b  !*FD Upper bound on interval for checking
    real(rk),intent(in) :: x  !*FD Test value (degrees)

    ! Internal variables

    real(rk) :: ta,tb

    ! Beginning of code

    if (a<b) then
      lon_between=((x>=a).and.(x<b))
    else
      if (x<a) then
        ta=a-360.0
        tb=b
      else 
        ta=a
        tb=b+360.0
      endif
      lon_between=((x>=ta).and.(x<tb))
    endif

  end function lon_between

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine calc_fractional(x,y,xp,yp,xa,ya,xb,yb,xc,yc,xd,yd)

    !*FD Performs a coordinate transformation to locate the point
    !*FD $(X',Y')$ fractionally within an arbitrary quadrilateral, 
    !*FD defined by the points $(x_A,y_A)$, $(x_B,y_B)$, 
    !*FD $(x_C,y_C)$ and $(x_D,y_D)$, which are ordered 
    !*FD anticlockwise.

    real(rk),intent(out) :: x !*FD The fractional $x$ location.
    real(rk),intent(out) :: y !*FD The fractional $y$ location.
    real(rk),intent(in)  :: xp,yp,xa,ya,xb,yb,xc,yc,xd,yd

    real(rk) :: a,b,c

    a=(yb-ya)*(xc-xd)-(yc-yd)*(xb-xa)

    b=xp*(yc-yd)-yp*(xc-xd) &
     +xd*(yb-ya)-yd*(xb-xa) &
     -xp*(yb-ya)+yp*(xb-xa) &
     -xa*(yc-yd)+ya*(xc-xd) 
         
    c=xp*(yd-ya)+yp*(xa-xd)+ya*xd-xa*yd

    if (a/=0.0) then
      x=(-b-sqrt(b**2-4*a*c))/(2*a)
    else
      x=-c/b
    endif

    y=(yp-ya-x*(yb-ya))/(yd+x*(yc-yd-yb+ya)-ya)

  end subroutine calc_fractional

end module glimmer_interp
