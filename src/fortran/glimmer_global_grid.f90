
module glimmer_global_grid

  use glimmer_global

  implicit none

  ! ------------------------------------------------------------
  ! GLOBAL_GRID derived type
  ! ------------------------------------------------------------

  type global_grid
    
    !*FD Contains parameters specifying the global grid configuration.

    ! Dimensions of grid ---------------------------------------

    integer :: nx = 0  !*FD Number of points in the $x$-direction.
    integer :: ny = 0  !*FD Number of points in the $y$-direction.

    ! Locations of grid-points ---------------------------------

    real(rk),pointer,dimension(:) :: lats      => null() !*FD Latitudinal locations of 
                                                         !*FD data-points in global fields (degrees)
    real(rk),pointer,dimension(:) :: lons      => null() !*FD Longitudinal locations of 
                                                         !*FD data-points in global fields (degrees)

    ! Locations of grid-box boundaries -------------------------

    real(rk),pointer,dimension(:) :: lat_bound => null() !*FD Latitudinal boundaries of 
                                                         !*FD data-points in global fields (degrees)
    real(rk),pointer,dimension(:) :: lon_bound => null() !*FD Longitudinal boundaries of 
                                                         !*FD data-points in global fields (degrees)

  end type global_grid

contains

  subroutine new_global_grid(grid,lons,lats,lonb,latb)

    !*FD Initialises a new global grid type

    type(global_grid),             intent(inout) :: grid !*FD The grid to be initialised
    real(rk),dimension(:),         intent(in)    :: lons !*FD Longitudinal positions of grid-points (degrees)
    real(rk),dimension(:),         intent(in)    :: lats !*FD Latitudinal positions of grid-points (degrees)
    real(rk),dimension(:),optional,intent(in)    :: lonb !*FD Longitudinal boundaries of grid-boxes (degrees)
    real(rk),dimension(:),optional,intent(in)    :: latb !*FD Latitudinal boundaries of grid-boxes (degrees)

    ! Check to see if things are allocated, and if so, deallocate them

    if (associated(grid%lats))      deallocate(grid%lats)
    if (associated(grid%lons))      deallocate(grid%lons)
    if (associated(grid%lat_bound)) deallocate(grid%lat_bound)
    if (associated(grid%lon_bound)) deallocate(grid%lon_bound)

    ! Find size of grid

    grid%nx=size(lons) ; grid%ny=size(lats)

    ! Allocate arrays

    allocate(grid%lons(grid%nx))
    allocate(grid%lats(grid%ny))
    allocate(grid%lon_bound(grid%nx+1))
    allocate(grid%lat_bound(grid%ny+1))

    ! Check dimensions of boundary arrays, if supplied

    if (present(lonb)) then
      if (.not.size(lonb)==grid%nx+1) then
        print*,'ERROR: new_global_grid. lonb mismatch'
        stop
      endif
    endif

    if (present(latb)) then
      if (.not.size(latb)==grid%ny+1) then
        print*,'ERROR: new_global_grid. latb mismatch'
        stop
      endif
    endif

    ! Copy lats and lons over

    grid%lats=lats
    grid%lons=lons

    ! Calculate boundaries if necessary

    if (present(lonb)) then
      grid%lon_bound=lonb
    else
      call calc_bounds_lon(lons,grid%lon_bound)
    endif

    if (present(latb)) then
      grid%lat_bound=latb
    else
      call calc_bounds_lat(lats,grid%lat_bound)
    endif

  end subroutine new_global_grid

!-----------------------------------------------------------------------------

  subroutine calc_bounds_lon(lons,lonb)
  
    !*FD Calculates the longitudinal boundaries between
    !*FD global grid-boxes. Note that we assume that the boundaries lie 
    !*FD half-way between the points, although 
    !*FD this isn't strictly true for a Gaussian grid.

    implicit none

    real(rk),dimension(:),intent(in)  :: lons !*FD locations of global grid-points (degrees)
    real(rk),dimension(:),intent(out) :: lonb !*FD boundaries of grid-boxes (degrees)

    integer :: nxg,i

    nxg=size(lons)

    ! Longitudes
 
    do i=1,nxg-1
      lonb(i+1)=mid_lon(lons(i),lons(i+1))
    enddo

    lonb(1)=mid_lon(lons(nxg),lons(1))
    lonb(nxg+1)=lonb(1)

  end subroutine calc_bounds_lon

!---------------------------------------------------------------------------------

  subroutine calc_bounds_lat(lat,latb)

  !*FD Calculates the boundaries between
  !*FD global grid-boxes. Note that we assume that the boundaries lie 
  !*FD half-way between the 
  !*FD points, both latitudinally and longitudinally, although 
  !*FD this isn't strictly true for a Gaussian grid.

    implicit none

    real(rk),dimension(:),intent(in)  :: lat !*FD locations of global grid-points (degrees)
    real(rk),dimension(:),intent(out) :: latb !*FD boundaries of grid-boxes (degrees)

    integer :: nyg,j

    nyg=size(lat)

    ! Latitudes first - we assume the boundaries of the first and 
    ! last boxes coincide with the poles. Not sure how to
    ! handle it if they don't...

    latb(1)=90.0
    latb(nyg+1)=-90.0

    do j=2,nyg
      latb(j)=lat(j-1)-(lat(j-1)-lat(j))/2.0
    enddo

  end subroutine calc_bounds_lat

!-------------------------------------------------------------

  real(rk) function mid_lon(a,b)

    !*FD Calculates the mid-point between two longitudes.
    !*FD \texttt{a} must be west of \texttt{b}.

    real(rk),intent(in) :: a,b

    real(rk) :: aa,bb,out

    aa=a ; bb=b

    if (aa>360.0.or.aa<0.0  .or. &
        bb>360.0.or.bb<0.0) then
        Print*,'Out of range in mid_lon'
        stop
    endif

    if (aa>bb) aa=aa-360.0

    out=aa+((bb-aa)/2.0)

    do
      if (out<=360.0) exit
      out=out-360.0
    end do

    do
      if (out>=0.0) exit
      out=out+360.0
    end do
    
    mid_lon=out

  end function mid_lon

end module glimmer_global_grid
