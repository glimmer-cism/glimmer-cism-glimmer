! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_outp.f90 - part of the GLIMMER ice model         + 
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

#define NC outfile%nc

module bin2ncdf_mod
  !*FD module containing some utility routines used for converting old glimmer binary output
  !*FD to new glimmer netCDF output
  !*FD written by Magnus Hagdorn, May 2004

  use glimmer_global, only : sp
  
  integer numx0, numy0
  integer numx1, numy1
  integer numz
  integer numspots
  real(kind=sp), dimension(:), allocatable :: x0data, y0data
  real(kind=sp), dimension(:), allocatable :: x1data, y1data
  real(kind=sp), dimension(:), allocatable :: zdata
  real(kind=sp), dimension(:), allocatable :: x0_spots, x1_spots, y0_spots, y1_spots
  real(kind=sp), dimension(:,:), allocatable :: data0d
  real(kind=sp), dimension(:,:), allocatable :: data2d_0, data2d_1
  real(kind=sp), dimension(:,:,:), allocatable :: data3d_0, data3d_1
  real, dimension(:), allocatable :: times
  integer, parameter :: numvar0d=14
  integer, dimension(numvar0d) :: var0d_order

contains
  subroutine timeslice(time,outfile)
    !*FD move time counter
    use glimmer_ncdf
    implicit none
    real,intent(in) :: time
    type(glimmer_nc_output), pointer :: outfile

    ! local variables
    integer i, status, numt
    logical found

    status = nf90_inquire_dimension(NC%id,NC%timedim,len=numt)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (allocated(times)) then
       if (size(times) .lt. numt) then
          deallocate(times)
          allocate(times(numt))
          status = nf90_get_var(NC%id,NC%timevar,times)
          call nc_errorhandle(__FILE__,__LINE__,status)
       end if
       
       do i=1,numt
          found = times(i).eq.time
          if (found) then
             exit
          end if
       end do
       if (.not.found) then
          outfile%timecounter=numt+1
       end if
    else
       allocate(times(1))
       times(1) = time
       found = .false.
    end if
          
    if (.not.found) then
       status = nf90_put_var(NC%id, NC%timevar,time,(/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if   
  end subroutine timeslice

  subroutine open0d(fname,unit)
    !*FD open a 0D file
    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(inout) :: unit

    ! local vars
    integer i
    integer idum

#ifdef CVF
    open(unit,file=fname,form='unformatted',convert='BIG_ENDIAN')
    ! convert is a non-standard specifier and doesn't work with the
    ! Intel compiler.
    ! Substituting with
#else
    open(unit,file=fname,form='unformatted')
#endif

    ! get number of spots
    read(unit) numspots

    ! get coord definitions
    allocate(x0_spots(numspots))
    read(unit) (x0_spots(i),i=1,numspots)
    allocate(y0_spots(numspots))
    read(unit) (y0_spots(i),i=1,numspots)
    
    allocate(x1_spots(numspots))
    read(unit) (x1_spots(i),i=1,numspots)
    allocate(y1_spots(numspots))
    read(unit) (y1_spots(i),i=1,numspots)

    allocate(data0d(numspots,numvar0d))
  end subroutine open0d

  subroutine open2d(fname,unit)
    !*FD open a 2D file
    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(inout) :: unit

    ! local vars
    integer i
    integer idum

#ifdef CVF
    open(unit,file=fname,form='unformatted',convert='BIG_ENDIAN')
    ! convert is a non-standard specifier and doesn't work with the
    ! Intel compiler.
    ! Substituting with
#else
    open(unit,file=fname,form='unformatted')
#endif

    ! get coord definitions
    read(unit) idum,idum,numy0,numx0
    allocate(x0data(numx0))
    allocate(y0data(numy0))
    allocate(data2d_0(numx0,numy0))
    read(unit) (y0data(i),i=1,numy0)
    read(unit) (x0data(i),i=1,numx0)

    read(unit) idum,idum,numy1,numx1
    allocate(x1data(numx1))
    allocate(y1data(numy1))
    allocate(data2d_1(numx1,numy1))
    read(unit) (y1data(i),i=1,numy1)
    read(unit) (x1data(i),i=1,numx1)
  end subroutine open2d

  subroutine open3d(fname,unit)
    !*FD open a 3D file
    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(inout) :: unit

    ! local vars
    integer i
    integer idum

#ifdef CVF
    open(unit,file=fname,form='unformatted',convert='BIG_ENDIAN')
    ! convert is a non-standard specifier and doesn't work with the
    ! Intel compiler.
    ! Substituting with
#else
    open(unit,file=fname,form='unformatted')
#endif

    ! get coord definitions
    read(unit) idum,idum,numy0,numx0,numz
    allocate(x0data(numx0))
    allocate(y0data(numy0))
    allocate(data3d_0(numx0,numy0,numz))
    read(unit) (y0data(i),i=1,numy0)
    read(unit) (x0data(i),i=1,numx0)

    read(unit) idum,idum,numy1,numx1,numz
    allocate(x1data(numx1))
    allocate(y1data(numy1))
    allocate(data3d_1(numx1,numy1,numz))
    read(unit) (y1data(i),i=1,numy1)
    read(unit) (x1data(i),i=1,numx1)

    allocate(zdata(numz))
  end subroutine open3d

  subroutine close0d(unit)
    !*FD close a 2D file
    implicit none
    integer,intent(in) :: unit

    close(unit)
    deallocate(x0_spots, y0_spots)
    deallocate(x1_spots, y1_spots)
    deallocate(data0d)
  end subroutine close0d

  subroutine close2d(unit)
    !*FD close a 2D file
    implicit none
    integer,intent(in) :: unit

    close(unit)
    deallocate(x0data,y0data,data2d_0)
    deallocate(x1data,y1data,data2d_1)
  end subroutine close2d

  subroutine close3d(unit)
    !*FD close a 2D file
    implicit none
    integer,intent(in) :: unit

    close(unit)
    deallocate(x0data,y0data,data3d_0)
    deallocate(x1data,y1data,data3d_1)
    deallocate(zdata)
  end subroutine close3d

  function map2d(id,index)
    !*FD map id and index to netCDF variable id
    use glimmer_ncdf
    implicit none
    integer :: map2d
    integer, intent(in) :: id,index

    if (id==0) then
       select case(index)
       case(1)
          map2d = NC_B_UFLX
          return
       case(2)
          map2d = NC_B_VFLX
          return
       case(3)
          map2d = NC_B_DIFFU
          return
       case(4)
          map2d = NC_B_BTRC
          return
       case(5)
          map2d = NC_B_UBAS
          return
       case(6)
          map2d = NC_B_VBAS
          return
       end select
    else if (id==1) then
       select case(index)
       case(1)
          map2d = NC_B_THK
          return
       case(2)
          map2d = NC_B_USURF
          return
       case(3)
          map2d = NC_B_LSURF
          return
       case(4)
          map2d = NC_B_TOPG
          return
       case(5)
          map2d = NC_B_ACAB
          return
       case(6)
          map2d = NC_B_BMLT
          return
       case(7)
          map2d = NC_B_BWAT
          return
       case(8)
          map2d = NC_B_ARTM
          return
       case(9)
          map2d = NC_B_BTEMP
          return
       case(10)
          map2d = NC_B_ARNG
          return
       case(11)
          map2d = NC_B_PRCP
          return
       case(12)
          map2d = NC_B_ABLT
          return
       case(13)
          map2d = NC_B_DUSRFDTM
          return
       end select
    end if
    write(*,*) 'some error occured, we should not be here'
    stop
  end function map2d

  function map3d(id,index)
    !*FD map id and index to netCDF variable id
    use glimmer_ncdf
    implicit none
    integer :: map3d
    integer, intent(in) :: id,index

    if (id==0) then
       select case(index)
       case(2)
          map3d = NC_B_UVEL
          return
       case(3)
          map3d = NC_B_VVEL
          return
!!$       case(4)
!!$          map3d = NC_B_EFVS
!!$          return
!!$       case(5)
!!$          map3d = NC_B_TAU
!!$          return
!!$       case(6)
!!$          map3d = NC_B_TAUXZ
!!$          return
!!$       case(7)
!!$          map3d = NC_B_TAUYY
!!$          return
!!$       case(8)
!!$          map3d = NC_B_TAUXY
!!$          return
!!$       case(9)
!!$          map3d = NC_B_TAUXX
!!$          return
!!$       case(10)
!!$          map3d = NC_B_TAUYY
!!$          return
!!$       case(11)
!!$          map3d = NC_B_GDSX
!!$          return
!!$       case(12)
!!$          map3d = NC_B_GDSY
!!$          return
       end select
    else if (id==1) then
       select case(index)
       case(2)
          map3d = NC_B_WVEL
          return
       case(3)
          map3d = NC_B_WGRD
          return
       case(4)
          map3d = NC_B_FLWA
          return
       case(5)
          map3d = NC_B_TEMP
          return
       end select
    end if
    write(*,*) 'some error occured, we should not be here'
    stop
  end function map3d

  subroutine scan0d(vars)
    !*FD return which spots...
    use glimmer_ncdf
    implicit none
    logical, dimension(:) :: vars

    vars(NC_B_UFLX_SPOT) = .true.
    vars(NC_B_VFLX_SPOT) = .true.
    vars(NC_B_UBAS_SPOT) = .true.
    vars(NC_B_VBAS_SPOT) = .true.
    vars(NC_B_BTRC_SPOT) = .true.
    vars(NC_B_THK_SPOT) = .true.
    vars(NC_B_USURF_SPOT) = .true.
    vars(NC_B_LSURF_SPOT) = .true.
    vars(NC_B_TOPG_SPOT) = .true.
    vars(NC_B_ACAB_SPOT) = .true.
    vars(NC_B_BMLT_SPOT) = .true.
    vars(NC_B_BWAT_SPOT) = .true.
    vars(NC_B_ARTM_SPOT) = .true.
    vars(NC_B_BTEMP_SPOT) = .true.

    var0d_order(1) = NC_B_UFLX_SPOT
    var0d_order(2) = NC_B_VFLX_SPOT
    var0d_order(3) = NC_B_UBAS_SPOT
    var0d_order(4) = NC_B_VBAS_SPOT
    var0d_order(5) = NC_B_BTRC_SPOT
    var0d_order(6) = NC_B_THK_SPOT
    var0d_order(7) = NC_B_USURF_SPOT
    var0d_order(8) = NC_B_LSURF_SPOT
    var0d_order(9) = NC_B_TOPG_SPOT
    var0d_order(10) = NC_B_ACAB_SPOT
    var0d_order(11) = NC_B_BMLT_SPOT
    var0d_order(12) = NC_B_BWAT_SPOT
    var0d_order(13) = NC_B_ARTM_SPOT
    var0d_order(14) = NC_B_BTEMP_SPOT
  end subroutine scan0d

  subroutine scan2d(unit, vars)    
    !*FD scan a 2D file for variables
    use glimmer_ncdf
    implicit none
    integer, intent(in)   :: unit
    logical, dimension(:) :: vars

    ! local variables
    integer i,j,id,index
    integer :: ios=0
    real(kind=sp) :: rdummy

    do while (ios == 0)
       read(unit,iostat=ios) rdummy,index,id,rdummy
       if (id==0) then
          read(unit,iostat=ios) ((data2d_0(i,j),i=1,numx0),j=1,numy0)
       else
          read(unit,iostat=ios) ((data2d_1(i,j),i=1,numx1),j=1,numy1)
       end if
       vars(map2d(id,index)) = .true.
    end do
  end subroutine scan2d

  subroutine scan3d(unit, vars)    
    !*FD scan a 3D file for variables
    use glimmer_ncdf
    implicit none
    integer, intent(in)   :: unit
    logical, dimension(:) :: vars

    ! local variables
    integer i,j,k,z,id,index
    integer :: ios=0
    real(kind=sp) :: rdummy

    do while (ios == 0)
       read(unit,iostat=ios) rdummy,index,id,rdummy
       if (id==0) then
          read(unit,iostat=ios) (((data3d_0(i,j,k),k=1,numz),i=1,numx0),j=1,numy0)
       else
          read(unit,iostat=ios) (((data3d_1(i,j,k),k=1,numz),i=1,numx1),j=1,numy1)
       end if
       if (index.ne.1) then
          vars(map3d(id,index)) = .true.
       end if
    end do
  end subroutine scan3d

  subroutine copy0d(unit,outfile)
    !*FD copy 2D data to netCDF file
    use glimmer_ncdf
    implicit none
    integer, intent(in)   :: unit
    type(glimmer_nc_output), pointer :: outfile

    ! local variables
    integer :: ios=0
    integer i,v,status
    real(kind=sp) :: rdummy, time

    do while (ios == 0)
       write(*,*) 'nlub'
       read(unit,iostat=ios) time,rdummy,rdummy,rdummy, ((data0d(i,v),i=1,numspots),v=1,numvar0d)
       call timeslice(time,outfile)
       do v=1,numvar0d
          status = nf90_put_var(NC%id, NC%varids(var0d_order(v)), data0d(:,v), (/1,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end do
  end subroutine copy0d

  subroutine copy2d(unit,outfile)
    !*FD copy 2D data to netCDF file
    use glimmer_ncdf
    implicit none
    integer, intent(in)   :: unit
    type(glimmer_nc_output), pointer :: outfile

    ! local variables
    integer i,j,id,index,status
    integer :: ios=0
    real(kind=sp) :: rdummy, time
    
    do while (ios == 0)
       read(unit,iostat=ios) time,index,id,rdummy
       call timeslice(time,outfile)
       if (id==0) then
          read(unit,iostat=ios) ((data2d_0(i,j),i=1,numx0),j=1,numy0)
          status = nf90_put_var(NC%id, NC%varids(map2d(id,index)), data2d_0, (/1,1,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       else
          read(unit,iostat=ios) ((data2d_1(i,j),i=1,numx1),j=1,numy1)
          status = nf90_put_var(NC%id, NC%varids(map2d(id,index)), data2d_1, (/1,1,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end if
    end do
  end subroutine copy2d

  subroutine copy3d(unit,outfile)
    !*FD copy 3D data to netCDF file
    use glimmer_ncdf
    implicit none
    integer, intent(in)   :: unit
    type(glimmer_nc_output), pointer :: outfile

    ! local variables
    integer i,j,k,id,index,status
    integer :: ios=0
    real(kind=sp) :: rdummy, time
    
    do while (ios == 0)
       read(unit,iostat=ios) time,index,id,rdummy
       call timeslice(time,outfile)
       if (id==0) then
          read(unit,iostat=ios) (((data3d_0(i,j,k),k=1,numz),i=1,numx0),j=1,numy0)
          if (index.ne.1) then
             status = nf90_put_var(NC%id, NC%varids(map3d(id,index)), data3d_0, (/1,1,1,outfile%timecounter/))
             call nc_errorhandle(__FILE__,__LINE__,status)
          end if
       else
          read(unit,iostat=ios) (((data3d_1(i,j,k),k=1,numz),i=1,numx1),j=1,numy1)
          if (index.ne.1) then
             status = nf90_put_var(NC%id, NC%varids(map3d(id,index)), data3d_1, (/1,1,1,outfile%timecounter/))
             call nc_errorhandle(__FILE__,__LINE__,status)
          end if
       end if
    end do
  end subroutine copy3d
end module bin2ncdf_mod

program bin2ncdf
  use bin2ncdf_mod
  use glimmer_ncdf
  use glimmer_ncfile
  use glimmer_types
  use paramets, only : len0
  implicit none

  character(len=30) :: infile
  integer unit0d, unit2d, unit3d, status
  type(glimmer_nc_output), pointer :: outfile
  type(glimmer_global_type) :: model
  logical :: do0d,do2d,do3d

  allocate(outfile)

  write(*,*) 'Enter base name of input files'
  read(*,*) infile
  write(*,*) 'Enter name of output netCDF file'
  read(*,*) NC%filename

  ! first scan of input file to determine list of variables
  inquire (exist=do0d,file=trim(infile)//'.gl0')
  if (do0d) then
     write(*,*) 'Doing 0d...'
     call open0d(trim(infile)//'.gl0',unit0d)
     !call scan0d(NC%do_var)
     !NC%do_spot=.true.
  end if
  inquire (exist=do2d,file=trim(infile)//'.gl2')
  if (do0d) then
     write(*,*) 'Doing 2d...'
     call open2d(trim(infile)//'.gl2',unit2d)
     call scan2d(unit2d,NC%do_var)
     call close2d(unit2d)
  end if
  inquire (exist=do3d,file=trim(infile)//'.gl3')
  if (do3d) then
     write(*,*) 'Doing 3d...'
     call open3d(trim(infile)//'.gl3',unit3d)
     call scan3d(unit3d,NC%do_var)
     call close3d(unit3d)
  end if

  ! creating netCDF file
  model%general%upn = numz
  model%numerics%dew = (x1data(2)-x1data(1))/len0
  model%numerics%dns = (y1data(2)-y1data(1))/len0
  allocate(model%numerics%sigma(model%general%upn))
  model%numerics%sigma = 1
  model%general%ewn = numx1
  model%general%nsn = numy1
  call glimmer_nc_createfile(outfile, model)

  ! copy data
  if (do0d) then
     !call copy0d(unit0d,outfile)
     call close0d(unit0d)
  end if
  if (do2d) then
     call open2d(trim(infile)//'.gl2',unit2d)
     call copy2d(unit2d,outfile)
     call close2d(unit2d)
  end if
  if (do3d) then
     call open3d(trim(infile)//'.gl3',unit3d)
     call copy3d(unit3d,outfile)
     call close3d(unit3d) 
  end if

  status = nf90_close(NC%id)
end program bin2ncdf
  
  
