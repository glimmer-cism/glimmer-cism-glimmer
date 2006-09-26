
module glimmer_restart_common

  use netcdf
  use glimmer_global, only: sp,dp

  implicit none

  integer,parameter :: restartflen = 80
  integer,parameter :: varnamelen = 50
  character(1),dimension(3) :: com_dims=(/'x','y','t'/)

  type restart_file
     character(restartflen) :: fname
     integer :: ncid
     logical :: def
     integer :: count
     logical :: write   !If false, we're open to read
  end type restart_file

  ! Interfaces ----------------------------------------
  ! Pointer scalar variables

!!$  interface write_pointscalvar
!!$
!!$  end interface
!!$
!!$  interface read_pointscalvar
!!$
!!$  end interface

  ! Allocatable arrays

!!$  interface write_allocarr
!!$     module procedure write_allocarr_int_1d, write_allocarr_realsp_1d, write_allocarr_realdp_1d
!!$     module procedure write_allocarr_int_2d, write_allocarr_realsp_2d, write_allocarr_realdp_2d
!!$     module procedure write_allocarr_int_3d, write_allocarr_realsp_3d, write_allocarr_realdp_3d
!!$  end interface
!!$
!!$  interface read_allocarr
!!$     module procedure read_allocarr_int_1d, read_allocarr_realsp_1d, read_allocarr_realdp_1d
!!$     module procedure read_allocarr_int_2d, read_allocarr_realsp_2d, read_allocarr_realdp_2d
!!$     module procedure read_allocarr_int_3d, read_allocarr_realsp_3d, read_allocarr_realdp_3d
!!$  end interface

  ! Pointer arrays

  interface write_pointarr
     module procedure write_pointarr_int_1d, write_pointarr_realsp_1d, write_pointarr_realdp_1d, write_pointarr_logical_1d
     module procedure write_pointarr_int_2d, write_pointarr_realsp_2d, write_pointarr_realdp_2d, write_pointarr_logical_2d
     module procedure write_pointarr_int_3d, write_pointarr_realsp_3d, write_pointarr_realdp_3d, write_pointarr_logical_3d
  end interface

  interface read_pointarr
     module procedure read_pointarr_int_1d, read_pointarr_realsp_1d, read_pointarr_realdp_1d!, read_pointarr_logical_1d
     module procedure read_pointarr_int_2d, read_pointarr_realsp_2d, read_pointarr_realdp_2d, read_pointarr_logical_2d
     module procedure read_pointarr_int_3d, read_pointarr_realsp_3d, read_pointarr_realdp_3d!, read_pointarr_logical_3d
  end interface

contains

  !-----------------------------------------------------------------
  ! Shared stuff
  !-----------------------------------------------------------------

  function new_restart_file(fname)

    type(restart_file) :: new_restart_file
    character(*),intent(in) :: fname
    integer :: status
    
    new_restart_file%fname=fname
    status = nf90_create(fname,NF90_CLOBBER,new_restart_file%ncid)
    if (status/=0) call ncdf_err(status,__LINE__,fname)
    new_restart_file%def = .true.
    new_restart_file%count = 0
    new_restart_file%write = .true.

  end function new_restart_file

  !-----------------------------------------------------------------

  function open_restart_file(fname)

    type(restart_file) :: open_restart_file
    character(*),intent(in) :: fname
    integer :: status
    
    open_restart_file%fname=fname
    status = nf90_open(fname,NF90_NOWRITE,open_restart_file%ncid)
    if (status/=0) call ncdf_err(status,__LINE__,fname)
    open_restart_file%def = .false.
    open_restart_file%count = 0
    open_restart_file%write = .false.

  end function open_restart_file

  !------------------------------------------------------------------

  subroutine close_restart_file(file)

    type(restart_file) :: file
    integer :: status

    status = nf90_close(file%ncid)
    if (status/=0) call ncdf_err(status,__LINE__,file%fname)

  end subroutine close_restart_file

  !------------------------------------------------------------------

  subroutine new_dimension(file,prefix,len,id)

    type(restart_file) :: file
    character(1),intent(in) :: prefix
    integer,intent(in) :: len
    integer,intent(out) :: id

    integer :: status
    character(10) :: dimname

    call set_define(file)

    write(dimname,'(I9)')len
    dimname=trim(prefix)//trim(adjustl(dimname))

    status=nf90_inq_dimid(file%ncid,dimname,id)
    if (status==NF90_NOERR) then
       return
    else if (status==NF90_EBADDIM) then
       status=nf90_def_dim(file%ncid,dimname,len,id)
       if (status==NF90_EBADDIM) call ncdf_err(status,__LINE__,file%fname,prefix)
    else
       call ncdf_err(status,__LINE__)
    end if

  end subroutine new_dimension

  !------------------------------------------------------------------

  subroutine set_define(file)

    type(restart_file),intent(inout) :: file
    integer :: status

    if (.not.file%def) then
       status = nf90_redef(file%ncid)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,file%fname)
       file%def=.true.
    end if

  end subroutine set_define

  !------------------------------------------------------------------

  subroutine end_define(file)

    type(restart_file),intent(inout) :: file
    integer :: status

    if (file%def) then
       status = nf90_enddef(file%ncid)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,file%fname)
       file%def=.false.
    end if

  end subroutine end_define

  !------------------------------------------------------------------

  subroutine check_read(file)
    
    use glimmer_log

    type(restart_file),intent(inout) :: file
    
    if (file%write) then
       call write_log('Attempted to read from restart file open for write',GM_FATAL)
    end if

  end subroutine check_read

  !------------------------------------------------------------------

  subroutine check_write(file)
    
    use glimmer_log

    type(restart_file),intent(inout) :: file
    
    if (.not.file%write) then
       call write_log('Attempted to write to restart file open for read',GM_FATAL)
    end if

  end subroutine check_write

  !------------------------------------------------------------------
  ! POINTER ARRAYS
  !------------------------------------------------------------------
  ! WRITE code - 1D
  !------------------------------------------------------------------

  subroutine write_pointarr_int_1d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,dimension(:),pointer     :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer,dimension(1) :: shp
    integer,dimension(1) :: lbounds

    if (associated(values)) then
       shp=shape(values)
       lbounds=lbound(values)
    end if

    call write_array_common(file,prefix,name,NF90_INT,associated(values),shp,lbounds,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    end if

  end subroutine write_pointarr_int_1d

  !------------------------------------------------------------------
    
  subroutine write_pointarr_realsp_1d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(sp),dimension(:),pointer    :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer,dimension(1) :: shp
    integer,dimension(1) :: lbounds

    if (associated(values)) then
       shp=shape(values)
       lbounds=lbound(values)
    end if

    call write_array_common(file,prefix,name,NF90_FLOAT,associated(values),shp,lbounds,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    end if

  end subroutine write_pointarr_realsp_1d

  !------------------------------------------------------------------
    
  subroutine write_pointarr_realdp_1d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(dp),dimension(:),pointer    :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer,dimension(1) :: shp
    integer,dimension(1) :: lbounds

    if (associated(values)) then
       shp=shape(values)
       lbounds=lbound(values)
    end if

    call write_array_common(file,prefix,name,NF90_DOUBLE,associated(values),shp,lbounds,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    end if

  end subroutine write_pointarr_realdp_1d

  !------------------------------------------------------------------
    
  subroutine write_pointarr_logical_1d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    logical,dimension(:),pointer    :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer :: nx
    integer,allocatable,dimension(:) :: tmp
    integer,dimension(1) :: shp
    integer,dimension(1) :: lbounds

    if (associated(values)) then
       shp=shape(values)
       lbounds=lbound(values)
    end if

    call write_array_common(file,prefix,name,NF90_INT,associated(values),shp,lbounds,varid)

    if (associated(values)) then
       allocate(tmp(size(values)))
       where (values)
          tmp=1
       elsewhere
          tmp=0
       end where
       status=nf90_put_var(file%ncid,varid,tmp)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    end if

  end subroutine write_pointarr_logical_1d

  !------------------------------------------------------------------
  ! POINTER ARRAYS
  !------------------------------------------------------------------
  ! WRITE code - 2D
  !------------------------------------------------------------------
    
  subroutine write_pointarr_int_2d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,dimension(:,:),pointer   :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer,dimension(2) :: shp
    integer,dimension(2) :: lbounds

    if (associated(values)) then
       shp=shape(values)
       lbounds=lbound(values)
    end if

    call write_array_common(file,prefix,name,NF90_INT,associated(values),shp,lbounds,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    end if

  end subroutine write_pointarr_int_2d

  !------------------------------------------------------------------
  
  subroutine write_pointarr_realsp_2d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(sp),dimension(:,:),pointer  :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer,dimension(2) :: shp
    integer,dimension(2) :: lbounds

    if (associated(values)) then
       shp=shape(values)
       lbounds=lbound(values)
    end if

    call write_array_common(file,prefix,name,NF90_FLOAT,associated(values),shp,lbounds,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    end if

  end subroutine write_pointarr_realsp_2d

  !------------------------------------------------------------------
 
  subroutine write_pointarr_realdp_2d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(dp),dimension(:,:),pointer  :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer,dimension(2) :: shp
    integer,dimension(2) :: lbounds

    if (associated(values)) then
       shp=shape(values)
       lbounds=lbound(values)
    end if

    call write_array_common(file,prefix,name,NF90_DOUBLE,associated(values),shp,lbounds,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    end if

  end subroutine write_pointarr_realdp_2d

  !------------------------------------------------------------------
    
  subroutine write_pointarr_logical_2d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    logical,dimension(:,:),pointer    :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer,allocatable,dimension(:,:) :: tmp
    integer,dimension(2) :: shp
    integer,dimension(2) :: lbounds

    if (associated(values)) then
       shp=shape(values)
       lbounds=lbound(values)
    end if

    call write_array_common(file,prefix,name,NF90_INT,associated(values),shp,lbounds,varid)

    if (associated(values)) then
       allocate(tmp(size(values,1),size(values,2)))
       where (values)
          tmp=1
       elsewhere
          tmp=0
       end where
       status=nf90_put_var(file%ncid,varid,tmp)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    end if

  end subroutine write_pointarr_logical_2d
    
  !------------------------------------------------------------------
  ! POINTER ARRAYS
  !------------------------------------------------------------------
  ! WRITE code - 3D
  !------------------------------------------------------------------

  subroutine write_pointarr_int_3d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,dimension(:,:,:),pointer :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer,dimension(3) :: shp
    integer,dimension(3) :: lbounds

    if (associated(values)) then
       shp=shape(values)
       lbounds=lbound(values)
    end if

    call write_array_common(file,prefix,name,NF90_INT,associated(values),shp,lbounds,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    end if

  end subroutine write_pointarr_int_3d

  !------------------------------------------------------------------
    
  subroutine write_pointarr_realsp_3d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(sp),dimension(:,:,:),pointer :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer,dimension(3) :: shp
    integer,dimension(3) :: lbounds

    if (associated(values)) then
       shp=shape(values)
       lbounds=lbound(values)
    end if

    call write_array_common(file,prefix,name,NF90_FLOAT,associated(values),shp,lbounds,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    end if

  end subroutine write_pointarr_realsp_3d

  !------------------------------------------------------------------
    
  subroutine write_pointarr_realdp_3d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(dp),dimension(:,:,:),pointer :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer,dimension(3) :: shp
    integer,dimension(3) :: lbounds

    if (associated(values)) then
       shp=shape(values)
       lbounds=lbound(values)
    end if

    call write_array_common(file,prefix,name,NF90_DOUBLE,associated(values),shp,lbounds,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    end if

  end subroutine write_pointarr_realdp_3d

  !------------------------------------------------------------------
    
  subroutine write_pointarr_logical_3d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    logical,dimension(:,:,:),pointer    :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer,allocatable,dimension(:,:,:) :: tmp
    integer,dimension(3) :: shp
    integer,dimension(3) :: lbounds

    if (associated(values)) then
       shp=shape(values)
       lbounds=lbound(values)
    end if

    call write_array_common(file,prefix,name,NF90_INT,associated(values),shp,lbounds,varid)

    if (associated(values)) then
       allocate(tmp(size(values,1),size(values,2),size(values,3)))
       where (values)
          tmp=1
       elsewhere
          tmp=0
       end where
       status=nf90_put_var(file%ncid,varid,tmp)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    end if

  end subroutine write_pointarr_logical_3d

  !------------------------------------------------------------------

  subroutine write_array_common(file,prefix,name,typecode,assoc,lens,lbs,varid)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,           intent(in)    :: typecode
    logical,           intent(in)    :: assoc
    integer,dimension(:),intent(in)  :: lens
    integer,dimension(:),intent(in)  :: lbs
    integer,           intent(out)   :: varid

    character(varnamelen) :: varname
    integer :: status
    integer,dimension(size(lens)) :: dimids
    integer :: ndims,i

    ndims=size(lens)

    call new_varname(varname,prefix,file%count)

    if (assoc) then
       do i=1,ndims
          call new_dimension(file,com_dims(i),lens(i),dimids(i))
       end do
       call set_define(file)
       status=nf90_def_var(file%ncid,varname,typecode,dimids,varid)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
       call write_lbounds(file%ncid,varid,lbs)
    else
       call set_define(file)
       status=nf90_def_var(file%ncid,varname,typecode,varid=varid)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
       status=nf90_put_att(file%ncid,varid,name='NULL',values='NULL')
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    endif

    status=nf90_put_att(file%ncid,varid,name='varname',values=name)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)

    call end_define(file)
    file%count=file%count+1

  end subroutine write_array_common

  !------------------------------------------------------------------
  ! POINTER ARRAYS
  !------------------------------------------------------------------
  ! READ code - 1D
  !------------------------------------------------------------------

  subroutine read_pointarr_int_1d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,dimension(:),pointer     :: values

    integer :: varid,status
    integer,dimension(1) :: dimlens
    integer,dimension(1) :: lbounds
    integer,dimension(1) :: ubounds
    logical :: nullt

    if(associated(values)) then
       deallocate(values)
       values => null()
    end if

    call read_array_common(file,prefix,name,varid,nullt,dimlens,lbounds)

    if (.not.nullt) then
       ubounds=lbounds+dimlens-1
       allocate(values(lbounds(1):ubounds(1)))
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name)
    end if

  end subroutine read_pointarr_int_1d

  !------------------------------------------------------------------

  subroutine read_pointarr_realsp_1d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(sp),dimension(:),pointer    :: values

    integer :: varid,status
    integer,dimension(1) :: dimlens
    integer,dimension(1) :: lbounds
    integer,dimension(1) :: ubounds
    logical :: nullt

    if(associated(values)) then
       deallocate(values)
       values => null()
    end if

    call read_array_common(file,prefix,name,varid,nullt,dimlens,lbounds)

    if (.not.nullt) then
       ubounds=lbounds+dimlens-1
       allocate(values(lbounds(1):ubounds(1)))
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name)
    end if

  end subroutine read_pointarr_realsp_1d

  !------------------------------------------------------------------

  subroutine read_pointarr_realdp_1d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(dp),dimension(:),pointer    :: values

    integer :: varid,status
    integer,dimension(1) :: dimlens
    integer,dimension(1) :: lbounds
    integer,dimension(1) :: ubounds
    logical :: nullt

    if(associated(values)) then
       deallocate(values)
       values => null()
    end if

    call read_array_common(file,prefix,name,varid,nullt,dimlens,lbounds)

    if (.not.nullt) then
       ubounds=lbounds+dimlens-1
       allocate(values(lbounds(1):ubounds(1)))
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name)
    end if

  end subroutine read_pointarr_realdp_1d

  !------------------------------------------------------------------
  ! POINTER ARRAYS
  !------------------------------------------------------------------
  ! READ code - 2D
  !------------------------------------------------------------------

  subroutine read_pointarr_int_2d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,dimension(:,:),pointer  :: values

    integer :: varid,status
    integer,dimension(2) :: dimlens
    integer,dimension(2) :: lbounds
    integer,dimension(2) :: ubounds
    logical :: nullt

    if(associated(values)) then
       deallocate(values)
       values => null()
    end if

    call read_array_common(file,prefix,name,varid,nullt,dimlens,lbounds)

    if (.not.nullt) then
       ubounds=lbounds+dimlens-1
       allocate(values(lbounds(1):ubounds(1),lbounds(2):ubounds(2)))
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name)
    end if

  end subroutine read_pointarr_int_2d

  !------------------------------------------------------------------

  subroutine read_pointarr_realsp_2d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(sp),dimension(:,:),pointer  :: values

    integer :: varid,status
    integer,dimension(2) :: dimlens
    integer,dimension(2) :: lbounds
    integer,dimension(2) :: ubounds
    logical :: nullt

    if(associated(values)) then
       deallocate(values)
       values => null()
    end if

    call read_array_common(file,prefix,name,varid,nullt,dimlens,lbounds)

    if (.not.nullt) then
       ubounds=lbounds+dimlens-1
       allocate(values(lbounds(1):ubounds(1),lbounds(2):ubounds(2)))
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name)
    end if

  end subroutine read_pointarr_realsp_2d

  !------------------------------------------------------------------

  subroutine read_pointarr_realdp_2d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(dp),dimension(:,:),pointer  :: values

    integer :: varid,status
    integer,dimension(2) :: dimlens
    integer,dimension(2) :: lbounds
    integer,dimension(2) :: ubounds
    logical :: nullt

    if(associated(values)) then
       deallocate(values)
       values => null()
    end if

    call read_array_common(file,prefix,name,varid,nullt,dimlens,lbounds)

    if (.not.nullt) then
       ubounds=lbounds+dimlens-1
       allocate(values(lbounds(1):ubounds(1),lbounds(2):ubounds(2)))
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name)
    end if

  end subroutine read_pointarr_realdp_2d
 
  !------------------------------------------------------------------

  subroutine read_pointarr_logical_2d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    logical,dimension(:,:),pointer  :: values

    integer :: varid,status
    integer,dimension(2) :: dimlens
    integer,dimension(2) :: lbounds
    integer,dimension(2) :: ubounds
    logical :: nullt
    integer,dimension(:,:),allocatable :: tmp

    if(associated(values)) then
       deallocate(values)
       values => null()
    end if

    call read_array_common(file,prefix,name,varid,nullt,dimlens,lbounds)

    if (.not.nullt) then
       ubounds=lbounds+dimlens-1
       allocate(values(lbounds(1):ubounds(1),lbounds(2):ubounds(2)))
       allocate(tmp(lbounds(1):ubounds(1),lbounds(2):ubounds(2)))
        status=nf90_get_var(file%ncid,varid,tmp)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name)
       where (tmp==1)
          values=.true.
       elsewhere
          values=.false.
       end where
    end if

  end subroutine read_pointarr_logical_2d
    
  !------------------------------------------------------------------
  ! POINTER ARRAYS
  !------------------------------------------------------------------
  ! READ code - 3D
  !------------------------------------------------------------------

  subroutine read_pointarr_int_3d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,dimension(:,:,:),pointer :: values

    integer :: varid,status
    integer,dimension(3) :: dimlens
    integer,dimension(3) :: lbounds
    integer,dimension(3) :: ubounds
    logical :: nullt

    if(associated(values)) then
       deallocate(values)
       values => null()
    end if

    call read_array_common(file,prefix,name,varid,nullt,dimlens,lbounds)

    if (.not.nullt) then
       ubounds=lbounds+dimlens-1
       allocate(values(lbounds(1):ubounds(1),lbounds(2):ubounds(2),lbounds(3):ubounds(3)))
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name)
    end if

  end subroutine read_pointarr_int_3d

  !------------------------------------------------------------------

  subroutine read_pointarr_realsp_3d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(sp),dimension(:,:,:),pointer :: values

    integer :: varid,status
    integer,dimension(3) :: dimlens
    integer,dimension(3) :: lbounds
    integer,dimension(3) :: ubounds
    logical :: nullt

    if(associated(values)) then
       deallocate(values)
       values => null()
    end if

    call read_array_common(file,prefix,name,varid,nullt,dimlens,lbounds)

    if (.not.nullt) then
       ubounds=lbounds+dimlens-1
       allocate(values(lbounds(1):ubounds(1),lbounds(2):ubounds(2),lbounds(3):ubounds(3)))
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name)
    end if

  end subroutine read_pointarr_realsp_3d

  !------------------------------------------------------------------

  subroutine read_pointarr_realdp_3d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(dp),dimension(:,:,:),pointer :: values

    integer :: varid,status
    integer,dimension(3) :: dimlens
    integer,dimension(3) :: lbounds
    integer,dimension(3) :: ubounds
    logical :: nullt

    if(associated(values)) then
       deallocate(values)
       values => null()
    end if

    call read_array_common(file,prefix,name,varid,nullt,dimlens,lbounds)

    if (.not.nullt) then
       ubounds=lbounds+dimlens-1
       allocate(values(lbounds(1):ubounds(1),lbounds(2):ubounds(2),lbounds(3):ubounds(3)))
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name)
    end if

  end subroutine read_pointarr_realdp_3d

  !------------------------------------------------------------------

  subroutine read_array_common(file,prefix,name,varid,nullt,dimlens,lbounds)

    use glimmer_log

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,           intent(out)   :: varid
    logical,           intent(out)   :: nullt
    integer,dimension(:),intent(out) :: dimlens
    integer,dimension(:),intent(out) :: lbounds

    character(varnamelen) :: varname,nametest,nulltest,dimname
    integer :: status,namelen,i
    integer,dimension(size(dimlens)) :: dimids

    call new_varname(varname,prefix,file%count)

    status=nf90_inq_varid(file%ncid,varname,varid)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)

    status=nf90_inquire_attribute(file%ncid,varid,'varname',len=namelen)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    status=nf90_get_att(file%ncid,varid,'varname',nametest)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    if (namelen>varnamelen) then
       call write_log('Variable name too long',GM_FATAL,__FILE__,__LINE__)
    end if
    nametest(namelen+1:)=repeat(' ',varnamelen-namelen)

    if (name/=nametest) then
       call write_log('Restart read mismatch: '//trim(varname)//', ' &
            //trim(name)//', '//trim(nametest),GM_FATAL)
    end if

    status=nf90_get_att(file%ncid,varid,'NULL',nulltest)
    select case(status)
    case(NF90_NOERR)
       nullt=.true.
    case(NF90_ENOTATT)
       nullt=.false.
       status=nf90_inquire_variable(file%ncid,varid,dimids=dimids)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
       do i = 1,size(dimlens)
          status=nf90_inquire_dimension(file%ncid,dimids(i),name=dimname,len=dimlens(i))
          if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
          if (dimname(1:1)/=com_dims(i)) then
             print*,'DIMENSION name error:',dimname,com_dims(i),i
          end if
       end do
       call read_lbounds(file%ncid,varid,lbounds)
    case default
       call ncdf_err(status,__LINE__,prefix,name,varname)
    end select

    file%count=file%count+1

  end subroutine read_array_common

  !------------------------------------------------------------------
  ! Utility code
  !------------------------------------------------------------------

  subroutine new_varname(varname,prefix,count)

    character(*),intent(out) :: varname
    character(*),intent(in)  :: prefix
    integer,     intent(in)  :: count

    write(varname,'(I9)')count
    varname=trim(prefix)//trim(adjustl(varname))

  end subroutine new_varname

  !------------------------------------------------------------------

  subroutine write_lbounds(ncid,varid,lb)

    integer,intent(in) :: ncid
    integer,intent(in) :: varid
    integer,dimension(:),intent(in) :: lb

    integer :: nbound,i,status
    character :: itxt

    nbound=size(lb)

    do i=1,nbound
       write(itxt,'(I1)')i
       status=nf90_put_att(ncid,varid,'lbound'//itxt,lb(i))
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,' in WRITE_LBOUNDS')
    end do

  end subroutine write_lbounds

  !------------------------------------------------------------------

  subroutine read_lbounds(ncid,varid,lb)

    use netcdf

    integer,intent(in)  :: ncid
    integer,intent(in)  :: varid
    integer,dimension(:),intent(out) :: lb

    integer :: nbound,i,status
    character :: itxt

    nbound=size(lb)

    do i=1,nbound
       write(itxt,'(I1)')i
       status=nf90_get_att(ncid,varid,'lbound'//itxt,lb(i))
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,' in READ_LBOUNDS')
    end do

  end subroutine read_lbounds

  !------------------------------------------------------------------

  subroutine write_allocatable(file,prefix,name,alloc)

    type(restart_file), intent(inout) :: file
    character(*) :: prefix
    character(*) :: name
    logical :: alloc

    character(varnamelen) :: varname
    integer :: varid,status

    call set_define(file)

    call new_varname(varname,prefix,file%count)
    status=nf90_def_var(file%ncid,varname,NF90_CHAR,varid=varid)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    status=nf90_put_att(file%ncid,varid,name='varname',values=name)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)

    if (.not.alloc) then
       status=nf90_put_att(file%ncid,varid,name='UNALLOC',values='UNALLOC')
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    end if

    call end_define(file)
    file%count=file%count+1

  end subroutine write_allocatable

  !------------------------------------------------------------------

  subroutine read_allocatable(file,prefix,name,alloc)

    use glimmer_log

    type(restart_file), intent(inout) :: file
    character(*) :: prefix
    character(*) :: name
    logical :: alloc

    character(varnamelen) :: varname,nametest,nulltest
    integer :: varid,status,namelen

    call new_varname(varname,prefix,file%count)

    status=nf90_inq_varid(file%ncid,varname,varid)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)

    status=nf90_inquire_attribute(file%ncid,varid,'varname',len=namelen)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    status=nf90_get_att(file%ncid,varid,'varname',nametest)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    if (namelen>varnamelen) then
       call write_log('Variable name too long',GM_FATAL,__FILE__,__LINE__)
    end if
    nametest(namelen+1:)=repeat(' ',varnamelen-namelen)

    if (name/=nametest) then
       call write_log('Restart read mismatch: '//trim(varname)//', ' &
            //trim(name)//', '//trim(nametest),GM_FATAL)
    end if

    status=nf90_get_att(file%ncid,varid,'UNALLOC',nulltest)
    select case(status)
    case(NF90_NOERR)
       alloc=.false.
    case(NF90_ENOTATT)
       alloc=.true.
    case default
       call ncdf_err(status,__LINE__,prefix,name,varname)
    end select

    file%count=file%count+1

  end subroutine read_allocatable
  
  !------------------------------------------------------------------

  subroutine write_pointer(file,prefix,name,assoc)

    type(restart_file), intent(inout) :: file
    character(*) :: prefix
    character(*) :: name
    logical :: assoc

    character(varnamelen) :: varname
    integer :: varid,status

    call set_define(file)

    call new_varname(varname,prefix,file%count)
    status=nf90_def_var(file%ncid,varname,NF90_CHAR,varid=varid)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    status=nf90_put_att(file%ncid,varid,name='varname',values=name)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)

    if (.not.assoc) then
       status=nf90_put_att(file%ncid,varid,name='NULL',values='NULL')
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    end if

    call end_define(file)
    file%count=file%count+1

  end subroutine write_pointer

  !------------------------------------------------------------------

  subroutine write_null_array_pointer(file,prefix,name)

    type(restart_file), intent(inout) :: file
    character(*) :: prefix
    character(*) :: name

    character(varnamelen) :: varname
    integer :: varid,status

    call set_define(file)

    call new_varname(varname,prefix,file%count)
    status=nf90_def_var(file%ncid,varname,NF90_CHAR,varid=varid)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    status=nf90_put_att(file%ncid,varid,name='varname',values=name)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    status=nf90_put_att(file%ncid,varid,name='NULL',values='NULL')
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)

    file%count=file%count+1

  end subroutine write_null_array_pointer

  !------------------------------------------------------------------

  subroutine write_array_pointer(file,prefix,name,sh)

    type(restart_file), intent(inout) :: file
    character(*) :: prefix
    character(*) :: name
    integer,dimension(:) :: sh

    character(varnamelen) :: varname
    integer :: varid,status,dimid

    call set_define(file)
    call new_dimension(file,com_dims(1),size(sh),dimid)

    call new_varname(varname,prefix,file%count)
    status=nf90_def_var(file%ncid,varname,NF90_INT,(/dimid/),varid)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    status=nf90_put_att(file%ncid,varid,name='varname',values=name)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)

    call end_define(file)
    status=nf90_put_var(file%ncid,varid,sh)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)

    file%count=file%count+1

  end subroutine write_array_pointer

  !------------------------------------------------------------------

  subroutine read_array_pointer(file,prefix,name,sh,assoc)

    use glimmer_log

    type(restart_file), intent(inout) :: file
    character(*) :: prefix
    character(*) :: name
    integer,dimension(:) :: sh
    logical,intent(out) :: assoc

    character(varnamelen) :: varname,nametest,nulltest
    integer :: varid,status,namelen

    call new_varname(varname,prefix,file%count)

    status=nf90_inq_varid(file%ncid,varname,varid)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)

    status=nf90_inquire_attribute(file%ncid,varid,'varname',len=namelen)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    status=nf90_get_att(file%ncid,varid,'varname',nametest)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    if (namelen>varnamelen) then
       call write_log('Variable name too long',GM_FATAL,__FILE__,__LINE__)
    end if
    nametest(namelen+1:)=repeat(' ',varnamelen-namelen)

    if (name/=nametest) then
       call write_log('Restart read mismatch: '//trim(varname)//', ' &
            //trim(name)//', '//trim(nametest),GM_FATAL)
    end if

    status=nf90_get_att(file%ncid,varid,'NULL',nulltest)
    select case(status)
    case(NF90_NOERR)
       assoc=.false.
    case(NF90_ENOTATT)
       assoc=.true.
       status=nf90_get_var(file%ncid,varid,sh)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    case default
       call ncdf_err(status,__LINE__,prefix,name,varname)
    end select

    file%count=file%count+1

  end subroutine read_array_pointer

  !------------------------------------------------------------------

  subroutine read_pointer(file,prefix,name,assoc)

    use glimmer_log

    type(restart_file), intent(inout) :: file
    character(*) :: prefix
    character(*) :: name
    logical :: assoc

    character(varnamelen) :: varname,nametest,nulltest
    integer :: varid,status,namelen

    call new_varname(varname,prefix,file%count)

    status=nf90_inq_varid(file%ncid,varname,varid)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)

    status=nf90_inquire_attribute(file%ncid,varid,'varname',len=namelen)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    status=nf90_get_att(file%ncid,varid,'varname',nametest)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__,prefix,name,varname)
    if (namelen>varnamelen) then
       call write_log('Variable name too long',GM_FATAL,__FILE__,__LINE__)
    end if
    nametest(namelen+1:)=repeat(' ',varnamelen-namelen)

    if (name/=nametest) then
       call write_log('Restart read mismatch: '//trim(varname)//', ' &
            //trim(name)//', '//trim(nametest),GM_FATAL)
    end if

    status=nf90_get_att(file%ncid,varid,'NULL',nulltest)
    select case(status)
    case(NF90_NOERR)
       assoc=.false.
    case(NF90_ENOTATT)
       assoc=.true.
    case default
       call ncdf_err(status,__LINE__,prefix,name,varname)
    end select

    file%count=file%count+1

  end subroutine read_pointer

  !------------------------------------------------------------------

  subroutine ncdf_err(status,line,char1,char2,char3)

    integer :: status
    integer :: line
    character(*),optional :: char1,char2,char3
    
    character(150) :: char

    char = ''
    if (present(char1)) char=trim(char)//trim(char1)
    if (present(char2)) char=trim(char)//' '//trim(char2)
    if (present(char3)) char=trim(char)//' '//trim(char3)

    print*,'NetCDF ERROR: ',trim(nf90_strerror(status)),' at ',line
    if (present(char1)) print*,'Diagnostics: ',char
    stop

  end subroutine ncdf_err

end module glimmer_restart_common
