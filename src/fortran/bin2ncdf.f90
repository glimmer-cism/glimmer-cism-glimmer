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
  real(kind=sp), dimension(:), allocatable :: x0data, y0data
  real(kind=sp), dimension(:), allocatable :: x1data, y1data
  real(kind=sp), dimension(:,:), allocatable :: data2d_0, data2d_1
  real, dimension(:), allocatable :: times


contains
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

  subroutine close2d(unit)
    !*FD close a 2D file
    implicit none
    integer,intent(in) :: unit

    close(unit)
    deallocate(x0data,y0data,data2d_0)
    deallocate(x1data,y1data,data2d_1)
  end subroutine close2d

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
    write(*,*) 'some errer'
    stop
  end function map2d

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
       found = .false.
    end if
          
    if (.not.found) then
       status = nf90_put_var(NC%id, NC%timevar,time,(/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if   
  end subroutine timeslice

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

end module bin2ncdf_mod

program bin2ncdf
  use bin2ncdf_mod
  use glimmer_ncdf
  use glimmer_ncfile
  use glimmer_types
  implicit none

  character(len=30) :: file2d
  integer unit2d, status
  type(glimmer_nc_output), pointer :: outfile
  type(glimmer_global_type) :: model
  
  allocate(outfile)

  write(*,*) 'Enter name of 2D file'
  read(*,*) file2d
  write(*,*) 'Enter name of output netCDF file'
  read(*,*) NC%filename

  ! first scan of input file to determine list of variables
  call open2d(file2d,unit2d)
  call scan2d(unit2d,NC%do_var)
  call close2d(unit2d)

  ! creating netCDF file
  model%general%upn = 11
  model%numerics%dew = 1
  model%numerics%dns = 1
  allocate(model%numerics%sigma(model%general%upn))
  model%numerics%sigma = 1
  model%general%ewn = numx1
  model%general%nsn = numy1
  call glimmer_nc_createfile(outfile, model)

  ! copy data
  call open2d(file2d,unit2d)
  call copy2d(unit2d,outfile)
  call close2d(unit2d)

  status = nf90_close(NC%id)
end program bin2ncdf
  
  
