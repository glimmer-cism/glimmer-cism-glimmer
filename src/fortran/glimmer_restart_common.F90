
module glimmer_restart_common

  use netcdf
  use glimmer_global, only: sp,dp

  implicit none

  integer,parameter :: restartflen = 80
  integer,parameter :: varnamelen = 50

  type restart_file
     character(restartflen) :: fname
     integer :: ncid
     logical :: def
     integer :: count
  end type restart_file

  interface write_var
     module procedure write_var_int_scalar, write_var_real_sp_scalar, write_var_real_dp_scalar
     module procedure write_var_int_1d, write_var_real_sp_1d, write_var_real_dp_1d
     module procedure write_var_int_2d, write_var_real_sp_2d, write_var_real_dp_2d
     module procedure write_var_int_3d, write_var_real_sp_3d, write_var_real_dp_3d
  end interface

  interface read_var
     module procedure read_var_int_scalar, read_var_real_sp_scalar, read_var_real_dp_scalar
     module procedure read_var_real_sp_2d
  end interface

contains

  function new_restart_file(fname)

    type(restart_file) :: new_restart_file
    character(*),intent(in) :: fname
    integer :: status
    
    new_restart_file%fname=fname
    status = nf90_create(fname,NF90_CLOBBER,new_restart_file%ncid)
    if (status/=0) call ncdf_err(status)
    new_restart_file%def = .true.
    new_restart_file%count = 0
    
  end function new_restart_file

  !-----------------------------------------------------------------

  function open_restart_file(fname)

    type(restart_file) :: open_restart_file
    character(*),intent(in) :: fname
    integer :: status
    
    open_restart_file%fname=fname
    status = nf90_open(fname,NF90_NOWRITE,open_restart_file%ncid)
    if (status/=0) call ncdf_err(status)
    open_restart_file%def = .false.
    open_restart_file%count = 0
    
  end function open_restart_file

  !------------------------------------------------------------------

  subroutine close_restart_file(file)

    type(restart_file) :: file
    integer :: status

    status = nf90_close(file%ncid)
    if (status/=0) call ncdf_err(status)

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
       if (status==NF90_EBADDIM) call ncdf_err(status)
    else
       call ncdf_err(status)
    end if

  end subroutine new_dimension

  !------------------------------------------------------------------

  subroutine set_define(file)

    type(restart_file),intent(inout) :: file
    integer :: status

    if (.not.file%def) then
       status = nf90_redef(file%ncid)
       if (status/=NF90_NOERR) call ncdf_err(status)
       file%def=.true.
    end if

  end subroutine set_define

  !------------------------------------------------------------------

  subroutine end_define(file)

    type(restart_file),intent(inout) :: file
    integer :: status

    if (file%def) then
       status = nf90_enddef(file%ncid)
       if (status/=NF90_NOERR) call ncdf_err(status)
       file%def=.false.
    end if

  end subroutine end_define
  
  !------------------------------------------------------------------
  ! Write_var specific subroutines
  !------------------------------------------------------------------
  ! Scalar variables
  !------------------------------------------------------------------

  subroutine write_var_int_scalar(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,           intent(in)    :: values

    integer :: status,varid

    call write_scalar_common(file,prefix,name,NF90_INT,varid)
    status=nf90_put_var(file%ncid,varid,values)
    if (status/=NF90_NOERR) call ncdf_err(status)

  end subroutine write_var_int_scalar

  !------------------------------------------------------------------

  subroutine write_var_real_sp_scalar(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(sp),          intent(in)    :: values

    integer :: status,varid

    call write_scalar_common(file,prefix,name,NF90_FLOAT,varid)
    status=nf90_put_var(file%ncid,varid,values)
    if (status/=NF90_NOERR) call ncdf_err(status)

  end subroutine write_var_real_sp_scalar

  !------------------------------------------------------------------

  subroutine write_var_real_dp_scalar(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(dp),          intent(in)    :: values

    integer :: status,varid

    call write_scalar_common(file,prefix,name,NF90_DOUBLE,varid)
    status=nf90_put_var(file%ncid,varid,values)
    if (status/=NF90_NOERR) call ncdf_err(status)

  end subroutine write_var_real_dp_scalar

  !------------------------------------------------------------------

  subroutine write_scalar_common(file,prefix,name,typecode,varid)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,           intent(in)    :: typecode
    integer,           intent(out)   :: varid

    character(varnamelen) :: varname
    integer :: status

    call new_varname(varname,prefix,file%count)

    ! Create new variable, and label it
    call set_define(file)
    status=nf90_def_var(file%ncid,varname,typecode,varid=varid)
    if (status/=NF90_NOERR) call ncdf_err(status)
    status=nf90_put_att(file%ncid,varid,name='varname',values=name)
    if (status/=NF90_NOERR) call ncdf_err(status)

    call end_define(file)
    file%count=file%count+1

  end subroutine write_scalar_common

  !------------------------------------------------------------------
  ! 1D pointer arrays
  !------------------------------------------------------------------
    
  subroutine write_var_int_1d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,dimension(:),pointer     :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer :: nx

    if (associated(values)) nx=size(values)

    call write_1d_common(file,prefix,name,NF90_INT,associated(values),nx,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status)
    end if

  end subroutine write_var_int_1d

  !------------------------------------------------------------------
    
  subroutine write_var_real_sp_1d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(sp),dimension(:),pointer    :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer :: nx

    if (associated(values)) nx=size(values)

    call write_1d_common(file,prefix,name,NF90_FLOAT,associated(values),nx,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status)
    end if

  end subroutine write_var_real_sp_1d

  !------------------------------------------------------------------
    
  subroutine write_var_real_dp_1d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(dp),dimension(:),pointer    :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer :: nx

    if (associated(values)) nx=size(values)

    call write_1d_common(file,prefix,name,NF90_DOUBLE,associated(values),nx,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status)
    end if

  end subroutine write_var_real_dp_1d

  !------------------------------------------------------------------

  subroutine write_1d_common(file,prefix,name,typecode,assoc,nx,varid)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,           intent(in)    :: typecode
    logical,           intent(in)    :: assoc
    integer,           intent(in)    :: nx
    integer,           intent(out)   :: varid

    character(varnamelen) :: varname
    integer :: status,dimid

    call new_varname(varname,prefix,file%count)

    if (assoc) then
       call new_dimension(file,'x',nx,dimid)
       call set_define(file)
       status=nf90_def_var(file%ncid,varname,typecode,(/dimid/),varid)
       if (status/=NF90_NOERR) call ncdf_err(status)
    else
       call set_define(file)
       status=nf90_def_var(file%ncid,varname,typecode,varid=varid)
       if (status/=NF90_NOERR) call ncdf_err(status)
       status=nf90_put_att(file%ncid,varid,name='NULL',values='NULL')
       if (status/=NF90_NOERR) call ncdf_err(status)
    endif

    status=nf90_put_att(file%ncid,varid,name='varname',values=name)
    if (status/=NF90_NOERR) call ncdf_err(status)

    call end_define(file)
    file%count=file%count+1

  end subroutine write_1d_common

  !------------------------------------------------------------------
  ! 2D pointer arrays
  !------------------------------------------------------------------
    
  subroutine write_var_int_2d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,dimension(:,:),pointer   :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer :: nx,ny

    if (associated(values)) then
       nx=size(values,1)
       ny=size(values,2)
    end if

    call write_2d_common(file,prefix,name,NF90_INT,associated(values),nx,ny,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status)
    end if

  end subroutine write_var_int_2d

  !------------------------------------------------------------------
  
  subroutine write_var_real_sp_2d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(sp),dimension(:,:),pointer  :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer :: nx,ny

    if (associated(values)) then
       nx=size(values,1)
       ny=size(values,2)
    end if

    call write_2d_common(file,prefix,name,NF90_FLOAT,associated(values),nx,ny,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status)
    end if

  end subroutine write_var_real_sp_2d

  !------------------------------------------------------------------
 
  subroutine write_var_real_dp_2d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(dp),dimension(:,:),pointer  :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer :: nx,ny

    if (associated(values)) then
       nx=size(values,1)
       ny=size(values,2)
    end if

    call write_2d_common(file,prefix,name,NF90_DOUBLE,associated(values),nx,ny,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status)
    end if

  end subroutine write_var_real_dp_2d

  !------------------------------------------------------------------

  subroutine write_2d_common(file,prefix,name,typecode,assoc,nx,ny,varid)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,           intent(in)    :: typecode
    logical,           intent(in)    :: assoc
    integer,           intent(in)    :: nx
    integer,           intent(in)    :: ny
    integer,           intent(out)   :: varid

    character(varnamelen) :: varname
    integer :: status,dimid1,dimid2

    call new_varname(varname,prefix,file%count)

    if (assoc) then
       call new_dimension(file,'x',nx,dimid1)
       call new_dimension(file,'y',ny,dimid2)
       call set_define(file)
       status=nf90_def_var(file%ncid,varname,typecode,(/dimid1,dimid2/),varid)
       if (status/=NF90_NOERR) call ncdf_err(status)
    else
       call set_define(file)
       status=nf90_def_var(file%ncid,varname,typecode,varid=varid)
       if (status/=NF90_NOERR) call ncdf_err(status)
       status=nf90_put_att(file%ncid,varid,name='NULL',values='NULL')
       if (status/=NF90_NOERR) call ncdf_err(status)
    endif

    status=nf90_put_att(file%ncid,varid,name='varname',values=name)
    if (status/=NF90_NOERR) call ncdf_err(status)

    call end_define(file)
    file%count=file%count+1

  end subroutine write_2d_common

  !------------------------------------------------------------------
  ! 3D pointer arrays
  !------------------------------------------------------------------
    
  subroutine write_var_int_3d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,dimension(:,:,:),pointer :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer :: nx,ny,nt

    if (associated(values)) then
       nx=size(values,1)
       ny=size(values,2)
       nt=size(values,3)
    end if

    call write_3d_common(file,prefix,name,NF90_INT,associated(values),nx,ny,nt,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status)
    end if

  end subroutine write_var_int_3d

  !------------------------------------------------------------------
    
  subroutine write_var_real_sp_3d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(sp),dimension(:,:,:),pointer :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer :: nx,ny,nt

    if (associated(values)) then
       nx=size(values,1)
       ny=size(values,2)
       nt=size(values,3)
    end if

    call write_3d_common(file,prefix,name,NF90_FLOAT,associated(values),nx,ny,nt,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status)
    end if

  end subroutine write_var_real_sp_3d

  !------------------------------------------------------------------
    
  subroutine write_var_real_dp_3d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(dp),dimension(:,:,:),pointer :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer :: nx,ny,nt

    if (associated(values)) then
       nx=size(values,1)
       ny=size(values,2)
       nt=size(values,3)
    end if

    call write_3d_common(file,prefix,name,NF90_DOUBLE,associated(values),nx,ny,nt,varid)

    if (associated(values)) then
       status=nf90_put_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status)
    end if

  end subroutine write_var_real_dp_3d

  !------------------------------------------------------------------

  subroutine write_3d_common(file,prefix,name,typecode,assoc,nx,ny,nt,varid)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,           intent(in)    :: typecode
    logical,           intent(in)    :: assoc
    integer,           intent(in)    :: nx
    integer,           intent(in)    :: ny
    integer,           intent(in)    :: nt
    integer,           intent(out)   :: varid

    character(varnamelen) :: varname
    integer :: status,dimid1,dimid2,dimid3

    call new_varname(varname,prefix,file%count)

    if (assoc) then
       call new_dimension(file,'x',nx,dimid1)
       call new_dimension(file,'y',ny,dimid2)
       call new_dimension(file,'t',nt,dimid3)
       call set_define(file)
       status=nf90_def_var(file%ncid,varname,typecode,(/dimid1,dimid2,dimid3/),varid)
       if (status/=NF90_NOERR) call ncdf_err(status)
    else
       call set_define(file)
       status=nf90_def_var(file%ncid,varname,typecode,varid=varid)
       if (status/=NF90_NOERR) call ncdf_err(status)
       status=nf90_put_att(file%ncid,varid,name='NULL',values='NULL')
       if (status/=NF90_NOERR) call ncdf_err(status)
    endif

    status=nf90_put_att(file%ncid,varid,name='varname',values=name)
    if (status/=NF90_NOERR) call ncdf_err(status)

    call end_define(file)
    file%count=file%count+1

  end subroutine write_3d_common

  !------------------------------------------------------------------
  ! read_var specific subroutines
  !------------------------------------------------------------------
  ! Scalar variables
  !------------------------------------------------------------------

  subroutine read_var_int_scalar(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,           intent(out)   :: values

    integer :: varid,status

    call read_scalar_common(file,prefix,name,varid)

    status=nf90_get_var(file%ncid,varid,values)
    if (status/=NF90_NOERR) call ncdf_err(status)

  end subroutine read_var_int_scalar

  !------------------------------------------------------------------

  subroutine read_var_real_sp_scalar(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(sp),          intent(out)   :: values

    integer :: varid,status

    call read_scalar_common(file,prefix,name,varid)

    status=nf90_get_var(file%ncid,varid,values)
    if (status/=NF90_NOERR) call ncdf_err(status)

  end subroutine read_var_real_sp_scalar

  !------------------------------------------------------------------

  subroutine read_var_real_dp_scalar(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(dp),          intent(out)   :: values

    integer :: varid,status

    call read_scalar_common(file,prefix,name,varid)

    status=nf90_get_var(file%ncid,varid,values)
    if (status/=NF90_NOERR) call ncdf_err(status)

  end subroutine read_var_real_dp_scalar

  !------------------------------------------------------------------

  subroutine read_scalar_common(file,prefix,name,varid)

    use glimmer_log

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,           intent(inout) :: varid

    character(varnamelen) :: varname,nametest
    integer :: status,namelen

    call new_varname(varname,prefix,file%count)

    status=nf90_inq_varid(file%ncid,varname,varid)
    if (status/=NF90_NOERR) call ncdf_err(status)

    status=nf90_inquire_attribute(file%ncid,varid,'varname',len=namelen)
    if (status/=NF90_NOERR) call ncdf_err(status)
    status=nf90_get_att(file%ncid,varid,'varname',nametest)
    if (status/=NF90_NOERR) call ncdf_err(status)
    if (namelen>varnamelen) then
       call write_log('Variable name too long',GM_FATAL,__FILE__,__LINE__)
    end if
    nametest(namelen+1:)=repeat(' ',varnamelen-namelen)

    if (name/=nametest) then
       call write_log('Restart read mismatch: '//trim(varname)//', ' &
            //trim(name)//', '//trim(nametest),GM_FATAL)
    end if

    file%count=file%count+1

  end subroutine read_scalar_common

  !------------------------------------------------------------------
  ! 2D pointer arrays
  !------------------------------------------------------------------

  subroutine read_var_int_2d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,dimension(:,:),pointer  :: values

    integer :: varid,status
    integer,dimension(2) :: dimlens
    logical :: nullt

    if(associated(values)) then
       deallocate(values)
       values => null()
    end if

    call read_array_common(file,prefix,name,varid,nullt,dimlens)

    if (.not.nullt) then
       allocate(values(dimlens(1),dimlens(2)))
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status)
    end if

  end subroutine read_var_int_2d

  !------------------------------------------------------------------

  subroutine read_var_real_sp_2d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(sp),dimension(:,:),pointer  :: values

    integer :: varid,status
    integer,dimension(2) :: dimlens
    logical :: nullt

    if(associated(values)) then
       deallocate(values)
       values => null()
    end if

    call read_array_common(file,prefix,name,varid,nullt,dimlens)

    if (.not.nullt) then
       allocate(values(dimlens(1),dimlens(2)))
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status)
    end if

  end subroutine read_var_real_sp_2d

  !------------------------------------------------------------------

  subroutine read_array_common(file,prefix,name,varid,nullt,dimlens)

    use glimmer_log

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,           intent(out)   :: varid
    logical,           intent(out)   :: nullt
    integer,dimension(:),intent(out) :: dimlens

    character(varnamelen) :: varname,nametest,nulltest
    integer :: status,namelen,i
    integer,dimension(size(dimlens)) :: dimids

    call new_varname(varname,prefix,file%count)

    status=nf90_inq_varid(file%ncid,varname,varid)
    if (status/=NF90_NOERR) call ncdf_err(status)

    status=nf90_inquire_attribute(file%ncid,varid,'varname',len=namelen)
    if (status/=NF90_NOERR) call ncdf_err(status)
    status=nf90_get_att(file%ncid,varid,'varname',nametest)
    if (status/=NF90_NOERR) call ncdf_err(status)
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
       if (status/=NF90_NOERR) call ncdf_err(status)
       do i = 1,size(dimlens)
          status=nf90_inquire_dimension(file%ncid,dimids(i),len=dimlens(i))
          if (status/=NF90_NOERR) call ncdf_err(status)
       end do
    case default
       call ncdf_err(status)
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

  subroutine ncdf_err(status)

    integer :: status

    print*,'NetCDF ERROR: ',trim(nf90_strerror(status))
    stop

  end subroutine ncdf_err

end module glimmer_restart_common
