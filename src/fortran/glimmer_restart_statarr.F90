
module glimmer_restart_statarr

  use glimmer_restart_common

  implicit none

  interface write_statarr
     module procedure write_statarr_int_1d, write_statarr_realsp_1d, write_statarr_realdp_1d, write_statarr_char_1d
     module procedure write_statarr_int_2d, write_statarr_realsp_2d, write_statarr_realdp_2d
     module procedure write_statarr_int_3d, write_statarr_realsp_3d, write_statarr_realdp_3d
  end interface

  interface read_statarr
     module procedure read_statarr_int_1d, read_statarr_realsp_1d, read_statarr_realdp_1d, read_statarr_char_1d
     module procedure read_statarr_int_2d, read_statarr_realsp_2d, read_statarr_realdp_2d
     module procedure read_statarr_int_3d, read_statarr_realsp_3d, read_statarr_realdp_3d
  end interface

  private
  public :: write_statarr, read_statarr

contains

  !------------------------------------------------------------------
  ! STATIC ARRAYS
  !------------------------------------------------------------------
  ! WRITE code - 1D
  !------------------------------------------------------------------

  subroutine write_statarr_int_1d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,dimension(:),intent(in)  :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer :: nx

    nx=size(values)

    call write_stat_common(file,prefix,name,NF90_INT,(/nx/),varid)

    status=nf90_put_var(file%ncid,varid,values)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
 
  end subroutine write_statarr_int_1d

  !------------------------------------------------------------------

  subroutine write_statarr_realsp_1d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(sp),dimension(:),intent(in) :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer :: nx

    nx=size(values)

    call write_stat_common(file,prefix,name,NF90_FLOAT,(/nx/),varid)

    status=nf90_put_var(file%ncid,varid,values)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
 
  end subroutine write_statarr_realsp_1d

  !------------------------------------------------------------------

  subroutine write_statarr_realdp_1d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(dp),dimension(:),intent(in) :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer :: nx

    nx=size(values)

    call write_stat_common(file,prefix,name,NF90_DOUBLE,(/nx/),varid)

    status=nf90_put_var(file%ncid,varid,values)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
 
  end subroutine write_statarr_realdp_1d

  !------------------------------------------------------------------

  subroutine write_statarr_char_1d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    character(*),dimension(:),intent(in) :: values

    integer :: status,varid
    character(varnamelen) :: varname
    integer :: nx

    print*,'WRITING CHARACTER ARRAY RESTART NOT IMPLEMENTED'

  end subroutine write_statarr_char_1d

  !------------------------------------------------------------------
  ! STATIC ARRAYS
  !------------------------------------------------------------------
  ! WRITE code - 2D
  !------------------------------------------------------------------

  subroutine write_statarr_int_2d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,dimension(:,:),intent(in)  :: values

    integer :: status,varid
    character(varnamelen) :: varname

    call write_stat_common(file,prefix,name,NF90_INT,shape(values),varid)

    status=nf90_put_var(file%ncid,varid,values)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
 
  end subroutine write_statarr_int_2d

  !------------------------------------------------------------------

  subroutine write_statarr_realsp_2d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(sp),dimension(:,:),intent(in)  :: values

    integer :: status,varid
    character(varnamelen) :: varname

    call write_stat_common(file,prefix,name,NF90_FLOAT,shape(values),varid)

    status=nf90_put_var(file%ncid,varid,values)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
 
  end subroutine write_statarr_realsp_2d

  !------------------------------------------------------------------

  subroutine write_statarr_realdp_2d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(dp),dimension(:,:),intent(in)  :: values

    integer :: status,varid
    character(varnamelen) :: varname

    call write_stat_common(file,prefix,name,NF90_DOUBLE,shape(values),varid)

    status=nf90_put_var(file%ncid,varid,values)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
 
  end subroutine write_statarr_realdp_2d

  !------------------------------------------------------------------
  ! STATIC ARRAYS
  !------------------------------------------------------------------
  ! WRITE code - 3D
  !------------------------------------------------------------------

  subroutine write_statarr_int_3d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,dimension(:,:,:),intent(in)  :: values

    integer :: status,varid
    character(varnamelen) :: varname

    call write_stat_common(file,prefix,name,NF90_INT,shape(values),varid)

    status=nf90_put_var(file%ncid,varid,values)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
 
  end subroutine write_statarr_int_3d

  !------------------------------------------------------------------

  subroutine write_statarr_realsp_3d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(sp),dimension(:,:,:),intent(in)  :: values

    integer :: status,varid
    character(varnamelen) :: varname

    call write_stat_common(file,prefix,name,NF90_FLOAT,shape(values),varid)

    status=nf90_put_var(file%ncid,varid,values)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
 
  end subroutine write_statarr_realsp_3d
  
  !------------------------------------------------------------------

  subroutine write_statarr_realdp_3d(file,prefix,name,values)

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    real(dp),dimension(:,:,:),intent(in)  :: values

    integer :: status,varid
    character(varnamelen) :: varname

    call write_stat_common(file,prefix,name,NF90_DOUBLE,shape(values),varid)

    status=nf90_put_var(file%ncid,varid,values)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
 
  end subroutine write_statarr_realdp_3d

  !------------------------------------------------------------------
  ! STATIC ARRAYS
  !------------------------------------------------------------------
  ! WRITE code - common
  !------------------------------------------------------------------

  subroutine write_stat_common(file,prefix,name,typecode,sizes,varid)

    type(restart_file),  intent(inout) :: file
    character(*),        intent(in)    :: prefix
    character(*),        intent(in)    :: name
    integer,             intent(in)    :: typecode
    integer,dimension(:),intent(in)    :: sizes
    integer,             intent(out)   :: varid

    character(varnamelen) :: varname
    integer :: status
    integer,dimension(size(sizes)) :: dimid
    integer :: rank,i

    rank = size(sizes)

    call new_varname(varname,prefix,file%count)

    do i=1,rank
       call new_dimension(file,com_dims(i),sizes(i),dimid(i))
    end do
    call set_define(file)
    status=nf90_def_var(file%ncid,varname,typecode,dimid,varid)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)

    status=nf90_put_att(file%ncid,varid,name='varname',values=name)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)

    call end_define(file)
    file%count=file%count+1

  end subroutine write_stat_common

  !------------------------------------------------------------------
  ! STATIC ARRAYS
  !------------------------------------------------------------------
  ! READ code - 1D
  !------------------------------------------------------------------

  subroutine read_statarr_int_1d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),  intent(inout) :: file
    character(*),        intent(in)    :: prefix
    character(*),        intent(in)    :: name
    integer,dimension(:),intent(out)   :: values

    integer :: varid,status
    integer,dimension(1) :: dimlens
 
    call read_statarr_common(file,prefix,name,varid,dimlens)

    if (all(dimlens==shape(values))) then
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
    end if

  end subroutine read_statarr_int_1d

  !------------------------------------------------------------------

  subroutine read_statarr_realsp_1d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),  intent(inout) :: file
    character(*),        intent(in)    :: prefix
    character(*),        intent(in)    :: name
    real(sp),dimension(:),intent(out)  :: values

    integer :: varid,status
    integer,dimension(1) :: dimlens
 
    call read_statarr_common(file,prefix,name,varid,dimlens)

    if (all(dimlens==shape(values))) then
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
    end if

  end subroutine read_statarr_realsp_1d

  !------------------------------------------------------------------

  subroutine read_statarr_realdp_1d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),  intent(inout) :: file
    character(*),        intent(in)    :: prefix
    character(*),        intent(in)    :: name
    real(dp),dimension(:),intent(out)  :: values

    integer :: varid,status
    integer,dimension(1) :: dimlens
 
    call read_statarr_common(file,prefix,name,varid,dimlens)

    if (all(dimlens==shape(values))) then
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
    end if

  end subroutine read_statarr_realdp_1d

  !------------------------------------------------------------------

  subroutine read_statarr_char_1d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),  intent(inout) :: file
    character(*),        intent(in)    :: prefix
    character(*),        intent(in)    :: name
    character(*),dimension(:),intent(out)  :: values

    integer :: varid,status
    integer,dimension(1) :: dimlens

    print*,'READING CHARACTER ARRAY RESTART NOT IMPLEMENTED'

  end subroutine read_statarr_char_1d

  !------------------------------------------------------------------
  ! STATIC ARRAYS
  !------------------------------------------------------------------
  ! READ code - 2D
  !------------------------------------------------------------------

  subroutine read_statarr_int_2d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),  intent(inout) :: file
    character(*),        intent(in)    :: prefix
    character(*),        intent(in)    :: name
    integer,dimension(:,:),intent(out) :: values

    integer :: varid,status
    integer,dimension(2) :: dimlens
 
    call read_statarr_common(file,prefix,name,varid,dimlens)

    if (all(dimlens==shape(values))) then
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
    end if

  end subroutine read_statarr_int_2d

  !------------------------------------------------------------------

  subroutine read_statarr_realsp_2d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),  intent(inout) :: file
    character(*),        intent(in)    :: prefix
    character(*),        intent(in)    :: name
    real(sp),dimension(:,:),intent(out) :: values

    integer :: varid,status
    integer,dimension(2) :: dimlens
 
    call read_statarr_common(file,prefix,name,varid,dimlens)

    if (all(dimlens==shape(values))) then
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
    end if

  end subroutine read_statarr_realsp_2d

  !------------------------------------------------------------------

  subroutine read_statarr_realdp_2d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),  intent(inout) :: file
    character(*),        intent(in)    :: prefix
    character(*),        intent(in)    :: name
    real(dp),dimension(:,:),intent(out)  :: values

    integer :: varid,status
    integer,dimension(2) :: dimlens
 
    call read_statarr_common(file,prefix,name,varid,dimlens)

    if (all(dimlens==shape(values))) then
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
    end if

  end subroutine read_statarr_realdp_2d

  !------------------------------------------------------------------
  ! STATIC ARRAYS
  !------------------------------------------------------------------
  ! READ code - 2D
  !------------------------------------------------------------------

  subroutine read_statarr_int_3d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),  intent(inout) :: file
    character(*),        intent(in)    :: prefix
    character(*),        intent(in)    :: name
    integer,dimension(:,:,:),intent(out) :: values

    integer :: varid,status
    integer,dimension(3) :: dimlens
 
    call read_statarr_common(file,prefix,name,varid,dimlens)

    if (all(dimlens==shape(values))) then
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
    end if

  end subroutine read_statarr_int_3d

  !------------------------------------------------------------------

  subroutine read_statarr_realsp_3d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),  intent(inout) :: file
    character(*),        intent(in)    :: prefix
    character(*),        intent(in)    :: name
    real(sp),dimension(:,:,:),intent(out) :: values

    integer :: varid,status
    integer,dimension(3) :: dimlens
 
    call read_statarr_common(file,prefix,name,varid,dimlens)

    if (all(dimlens==shape(values))) then
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
    end if

  end subroutine read_statarr_realsp_3d

  !------------------------------------------------------------------

  subroutine read_statarr_realdp_3d(file,prefix,name,values)

    use glimmer_log

    type(restart_file),  intent(inout) :: file
    character(*),        intent(in)    :: prefix
    character(*),        intent(in)    :: name
    real(dp),dimension(:,:,:),intent(out)  :: values

    integer :: varid,status
    integer,dimension(3) :: dimlens
 
    call read_statarr_common(file,prefix,name,varid,dimlens)

    if (all(dimlens==shape(values))) then
       status=nf90_get_var(file%ncid,varid,values)
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
    end if

  end subroutine read_statarr_realdp_3d

  !------------------------------------------------------------------
  ! STATIC ARRAYS
  !------------------------------------------------------------------
  ! READ code - common
  !------------------------------------------------------------------

  subroutine read_statarr_common(file,prefix,name,varid,dimlens)

    use glimmer_log

    type(restart_file),intent(inout) :: file
    character(*),      intent(in)    :: prefix
    character(*),      intent(in)    :: name
    integer,           intent(out)   :: varid
    integer,dimension(:),intent(out) :: dimlens

    character(varnamelen) :: varname,nametest,nulltest
    integer :: status,namelen,i
    integer,dimension(size(dimlens)) :: dimids

    call new_varname(varname,prefix,file%count)

    status=nf90_inq_varid(file%ncid,varname,varid)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)

    status=nf90_inquire_attribute(file%ncid,varid,'varname',len=namelen)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
    status=nf90_get_att(file%ncid,varid,'varname',nametest)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
    if (namelen>varnamelen) then
       call write_log('Variable name too long',GM_FATAL,__FILE__,__LINE__)
    end if
    nametest(namelen+1:)=repeat(' ',varnamelen-namelen)

    if (name/=nametest) then
       call write_log('Restart read mismatch: '//trim(varname)//', ' &
            //trim(name)//', '//trim(nametest),GM_FATAL)
    end if

    status=nf90_inquire_variable(file%ncid,varid,dimids=dimids)
    if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
    do i = 1,size(dimlens)
       status=nf90_inquire_dimension(file%ncid,dimids(i),len=dimlens(i))
       if (status/=NF90_NOERR) call ncdf_err(status,__LINE__)
    end do

    file%count=file%count+1

  end subroutine read_statarr_common


end module glimmer_restart_statarr
