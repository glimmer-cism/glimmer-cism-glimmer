!
! $Id: smb_netcdf.f90,v 1.1.2.2 2005-05-12 10:09:24 icrutt Exp $
!

module smb_netcdf

  use netcdf
  use glimmer_global

  type basicAttInfo
     character(len=nf90_max_name) :: long_name
     character(len=nf90_max_name) :: units
     real                         :: missing_value
  end type basicAttInfo

  type dimInfo
     character(len=nf90_max_name) :: name ! dimension name
     integer                      :: len  ! length of dimension
     integer                      :: id   ! id assigned by lib call
  end type dimInfo

  ! used to store info about coordinates 
  type varInfo_1d
     character(len=nf90_max_name)         :: name   ! variable name
     integer                              :: len    ! length of variable
     integer                              :: xtype  ! type e.g. NF90_FLOAT
     integer                              :: ndims  ! num associated dims
     integer,dimension(nf90_max_var_dims) :: dimIDs ! IDs of assoc dims
     integer                              :: nAtts  ! num associated atts
     integer                              :: id     ! id assigned by lib call
     type(basicAttInfo)                   :: basicAtts
     real,pointer,dimension(:)            :: data
  end type varInfo_1d

  ! used to store 'actual' data
  type varInfo_2d
     character(len=nf90_max_name)         :: name   ! variable name
     integer                              :: len    ! length of variable
     integer                              :: xtype  ! type e.g. NF90_FLOAT
     integer                              :: ndims  ! num associated dims
     integer,dimension(nf90_max_var_dims) :: dimIDs ! IDs of assoc dims
     integer                              :: nAtts  ! num associated atts
     integer                              :: id     ! id assigned by lib call
     type(basicAttInfo)                   :: basicAtts
     real,pointer,dimension(:,:)          :: data
  end type varInfo_2d

  type varInfo_3d
     character(len=nf90_max_name)         :: name   ! variable name
     integer                              :: len    ! length of variable
     integer                              :: xtype  ! type e.g. NF90_FLOAT
     integer                              :: ndims  ! num associated dims
     integer,dimension(nf90_max_var_dims) :: dimIDs ! IDs of assoc dims
     integer                              :: nAtts  ! num associated atts
     integer                              :: id     ! id assigned by lib call
     type(basicAttInfo)                   :: basicAtts
     real,pointer,dimension(:,:,:)        :: data
  end type varInfo_3d

contains

  ! An example routine
  subroutine SMBExampleNetCDF(filename,nx,ny,dx,data)
    
    implicit none
    
    ! args
    character(len=100),intent(in)  :: filename  ! output filename
    integer,intent(in)             :: nx        ! num grid pts in x-direction
    integer,intent(in)             :: ny        ! num grid pts in y-direction
    integer,intent(in)             :: dx        ! grid resolution (meters)
    real(rk),intent(in),dimension(:,:) :: data  ! data to write out--one timeslice
    ! locals
    integer :: ii
    integer :: nDims  ! num dimensions (& coords)
    integer :: nVars  ! num vars
    integer :: ncid   ! dataset id
    type(dimInfo),pointer,dimension(:)    :: dims    ! dimension info
    type(varInfo_1d),pointer,dimension(:) :: coords  ! coord info
    type(varInfo_2d),pointer,dimension(:) :: vars    ! var info
    
    ! initialisation
    nDims = 3  ! x, y and time
    nVars = 1  ! just one data array
    
    allocate(dims(nDims))
    allocate(coords(nDims))
    allocate(vars(nVars))
    
    dims(1)%name = 'x-direction'
    dims(1)%len = nx
    dims(2)%name = 'y-direction'
    dims(2)%len = ny
    dims(3)%name = 'time'
    dims(3)%len = nf90_unlimited
    
    ! bit hard wired...
    allocate(coords(1)%data(nx))
    allocate(coords(2)%data(ny))
    ! NB Need to know num timesteps in advance
    ! ...recommend one per file
    allocate(coords(3)%data(1))
    allocate(vars(1)%data(nx,ny))
    
    coords(1)%name = 'x-direction'
    coords(1)%xtype = nf90_float
    coords(1)%basicAtts%units = 'meters'
    coords(2)%name = 'y-direction'
    coords(2)%xtype = nf90_float
    coords(2)%basicAtts%units = 'meters'
    coords(3)%name = 'time'
    coords(3)%xtype = nf90_float
    coords(3)%basicAtts%units = 'time units'
    
    vars(1)%name = 'model-data'
    vars(1)%xtype = nf90_float
    vars(1)%basicAtts%units = 'not sure'

    ! fill out the coordinate arrays--grid cell midpoints
    do ii=1,nx
       coords(1)%data(ii) = (ii*dx)-dx/2
    end do
    do ii=1,ny
       coords(2)%data(ii) = ii*dx-dx/2
    end do
    coords(3)%data(1) = 1
    vars(1)%data = data
    
    ! create a dataset
    call SMBCreateNetCDF(filename,ncid)
    call SMBInitNetCDF(ncid,dims,coords,vars)
    call SMBWriteCoordsNetCDF(ncid,coords)
    call SMBWrite2dVarsNetCDF(ncid,vars)
    ! A second call with appropriate timestep arg
    ! would augment dataset
    ! **NB time coord would need to be set apprpriately
    ! in advance
!  call SMBWrite2dVarsNetCDF(ncid,vars,2)
    call SMBCloseNetCDF(ncid)
    
    ! and read it
    !  call SMBReadNetCDF("foo.nc")
    
    ! clean up
    do ii=1,nDims
       if(associated(coords(ii)%data)) deallocate(coords(ii)%data)
    end do
    do ii=1,nVars
       if(associated(vars(ii)%data)) deallocate(vars(ii)%data)
    end do
    if(associated(dims))   deallocate(dims)
    if(associated(coords)) deallocate(coords)
    if(associated(vars))   deallocate(vars)
    
  end subroutine SMBExampleNetCDF

  ! Opens a new dataset & returns dataset id
  subroutine SMBCreateNetCDF(filename,ncid)

    implicit none

    ! args
    character(len=100),intent(in) :: filename
    integer,intent(out)           :: ncid
    ! locals
    integer                       :: status

    ! create the dataset
    status = nf90_create(filename,nf90_clobber,ncid)
    if (status /= nf90_noerr) call handle_err(status)

  end subroutine SMBCreateNetCDF

  ! Close an open dataset with given id
  subroutine SMBCloseNetCDF(ncid)

    implicit none

    ! args
    integer,intent(in)            :: ncid
    ! locals
    integer                       :: status

    ! close the dataset
    status = nf90_close(ncid)
    if (status /= nf90_noerr) call handle_err(status)

  end subroutine SMBCloseNetCDF

  ! Defines dims, coords and vars
  subroutine SMBInitNetCDF(ncid,dims,coords,vars)

    implicit none

    ! args
    integer,intent(in)                          :: ncid
    type(dimInfo),intent(inout),dimension(:)    :: dims
    type(varInfo_1d),intent(inout),dimension(:) :: coords
    type(varInfo_2d),intent(inout),dimension(:) :: vars
    ! locals
    integer :: ii,iid
    integer :: status
    integer,pointer,dimension(:) :: dimIDs

    ! check num dims = num coords
    if(size(dims) .ne. size(coords)) then
       print*, "ERROR: number of dimensions and coordinates differ: ", &
            size(dims), size(coords)
    end if

    ! allocate dimIDs array
    allocate(dimIDs(size(dims)))

    ! define the dimensions
    do ii=1,size(dims)
       status = nf90_def_dim(ncid,dims(ii)%name,dims(ii)%len,dims(ii)%id)
       if (status /= nf90_noerr) call handle_err(status)
    end do

    ! define the coordinates
    do ii=1,size(coords)
       iid=dims(ii)%id
       status = nf90_def_var(ncid,coords(ii)%name, &
            coords(ii)%xtype, (/iid/), coords(ii)%id)
       if (status /= nf90_noerr) call handle_err(status)
       dims(ii)%id=iid
       dimIDs(ii) = dims(ii)%id
    end do

    ! define the variables
    ! *NB* vars assumed to be defined over _all_ dims
    ! **A BIG assumption**
    do ii=1,size(vars)
       status = nf90_def_var(ncid,vars(ii)%name, &
            vars(ii)%xtype, dimIDs, vars(ii)%id)
       if (status /= nf90_noerr) call handle_err(status)
    end do

    ! set some attributes
    do ii=1,size(coords)
       status = nf90_put_att(ncid,coords(ii)%id,"units", &
            coords(ii)%basicAtts%units)
       if (status /= nf90_noerr) call handle_err(status)
    end do

    do ii=1,size(vars)
       status = nf90_put_att(ncid,vars(ii)%id,"units", &
            vars(ii)%basicAtts%units)
       if (status /= nf90_noerr) call handle_err(status)
    end do

    ! leave define mode--enter data mode
    status = nf90_enddef(ncid)
    if (status /= nf90_noerr) call handle_err(status)

    ! cleanup
    if(associated(dimIDs)) deallocate(dimIDs)

  end subroutine SMBInitNetCDF

  ! set values on coordinate axes
  subroutine SMBWriteCoordsNetCDF(ncid,coords)

    implicit none

    ! args
    integer,intent(in)                       :: ncid
    type(varInfo_1d),intent(in),dimension(:) :: coords
    ! locals
    integer :: ii
    integer :: status

    ! set the coordinate values
    do ii=1,size(coords)
       status = nf90_put_var(ncid,coords(ii)%id,coords(ii)%data)
       if (status /= nf90_noerr) call handle_err(status)
    end do
    
  end subroutine SMBWriteCoordsNetCDF

  ! set values for variables
  subroutine SMBWrite2dVarsNetCDF(ncid,vars,timestep)

    implicit none

    ! args
    integer,intent(in)                       :: ncid
    type(varInfo_2d),intent(in),dimension(:) :: vars
    integer,intent(in),optional              :: timestep
    ! locals
    integer :: ii
    integer :: status                        ! error status
    integer :: start

    if (present(timestep)) then
       start = timestep
    else
       start = 1
    end if

    ! set the actual data
    do ii=1,size(vars)
       status = nf90_put_var(ncid,vars(ii)%id,vars(ii)%data, &
            start=(/ 1, 1, start /))
       if (status /= nf90_noerr) call handle_err(status)
    end do

  end subroutine SMBWrite2dVarsNetCDF

  !! reading an unknown netcdf file
  subroutine SMBReadNetCDF(filename)

    implicit none

    ! args
    character(len=100),intent(in) :: filename

    ! locals
    integer :: ii,jj,kk
    integer :: ncid       ! dataset ID
    integer :: status     ! error status
    integer :: nDims      ! Num dimensions
    integer :: nVarDims   ! Num dims in a given var
    integer :: nVars      ! Num variables
    integer :: nGlobAtts  ! num gloabl attributes
    integer :: unlimDimID ! ID of the unlimited dimension
    type(dimInfo),pointer,dimension(:)    :: dims   ! dimension info
    type(varInfo_1d),pointer,dimension(:) :: coords ! info for 'coord' vars 
    type(varInfo_3d),pointer,dimension(:) :: vars   ! var info 

     ! open the dataset
    status = nf90_open(filename,nf90_nowrite,ncid)
    if (status /= nf90_noerr) call handle_err(status)

    ! preliminary inquiry
    status = nf90_inquire(ncid,nDims,nVars,nGlobAtts,unlimDimID)
    if (status /= nf90_noerr) call handle_err(status)

    ! allocate storage once numbers are known
    allocate(dims(nDims))
    allocate(coords(nDims))
    allocate(vars(nVars-nDims))  ! NB split vars: 'coords' & 'actual' vars

    ! dimension IDs are between 1 and nDims
    ! loop over IDs to glean information
    do ii=1,nDims
       status = nf90_inquire_dimension(ncid,ii, &
            dims(ii)%name,dims(ii)%len)
       if (status /= nf90_noerr) call handle_err(status)
    end do

    jj=1  ! coord counter
    kk=1  ! 'actual' var counter
    do ii=1,nVars
       ! is it a 'coord' or 'actual' var?
       status = nf90_inquire_variable(ncid,ii, &
            ndims=nVarDims)
       if (status /= nf90_noerr) call handle_err(status)
       ! load up a 'coordinate' variable
       if (nVarDims==1) then  
          status = nf90_inquire_variable(ncid,ii, &
            name=coords(jj)%name, xtype=coords(jj)%xtype, &
            ndims=coords(jj)%ndims, &
            dimids=coords(jj)%dimIDs, &
            nAtts=coords(jj)%nAtts)
          if (status /= nf90_noerr) call handle_err(status)
          allocate(coords(jj)%data(dims(jj)%len))
          status = nf90_get_var(ncid,ii,coords(jj)%data)
          if (status /= nf90_noerr) call handle_err(status)
!          print*,"coords: ", jj, coords(jj)%data(1)
          jj = jj+1
       else if (nVarDims==3) then
          status = nf90_inquire_variable(ncid,ii, &
               name=vars(kk)%name, xtype=vars(kk)%xtype, &
               ndims=vars(kk)%ndims, &
               dimids=vars(kk)%dimIDs, &
               nAtts=vars(kk)%nAtts)
          if (status /= nf90_noerr) call handle_err(status)
          allocate(vars(kk)%data(dims(vars(kk)%dimIDs(1))%len, &
               dims(vars(kk)%dimIDs(2))%len, &
               dims(vars(kk)%dimIDs(3))%len))
          status = nf90_get_var(ncid,ii,vars(kk)%data)
          if (status /= nf90_noerr) call handle_err(status)
 !         print*,"vars: ", kk, vars(kk)%data(1,1,1)
          kk = kk+1
       else
          print*, "num dims in var not as expected"
          stop "Stopped"
       end if
    end do

    ! close the dataset
    status = nf90_close(ncid)
    if (status /= nf90_noerr) call handle_err(status)

    ! deallocate
    do ii=1,nDims
       if(associated(coords(ii)%data)) deallocate(coords(ii)%data)
    end do
    do ii=1,nVars
       if(associated(vars(ii)%data)) deallocate(vars(ii)%data)
    end do
    if(associated(dims)) deallocate(dims)
    if(associated(coords)) deallocate(coords)
    if(associated(vars)) deallocate(vars)

  end subroutine SMBReadNetCDF

  !! simple error handler
  subroutine handle_err(status)
    
    implicit none
    
    integer,intent(in) :: status

    if(status /= nf90_noerr) then
       write (*,*) trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine handle_err

end module smb_netcdf
