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

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#define NC outfile%nc

program input2ncdf
  !*FD convert glimmer input fields to netCDF
  !*FD written by Magnus Hagdorn, May 2004
  use glimmer_ncdf
  use glimmer_ncfile
  use glimmer_CFproj
  use glide_types

  use glimmer_paramets, only : len0
  implicit none

  ! file names
  character(len=30) :: latfile, topofile, rtopofile, sdfile
  character(len=30) :: surffile, precepfile, pressurf,mask
  integer unit, status, ptype
  logical :: polar
  real :: longitude_of_central_meridian
  real :: latitude_of_projection_origin
  real :: scale_factor_at_proj_origin
  real :: standard_parallel
  real :: standard_parallel_2
  real :: false_easting
  real :: false_northing  
  type(glimmer_nc_output), pointer :: outfile
  type(glide_global_type) :: model

  integer nx,ny
  real time,delta
  real(kind=sp), dimension(:,:), pointer :: data
  integer,       dimension(:,:), pointer :: maskdata

  ! start logging
  call open_log(unit=50)

  allocate(outfile)
  write(*,*) 'Enter name of latitude file.'
  read(*,*) latfile
  inquire(file=latfile,exist=NC%do_var(NC_B_LAT))
  write(*,*) 'Enter name of topography file.'
  read(*,*) topofile
  inquire(file=topofile,exist=NC%do_var(NC_B_TOPG))
  write(*,*) 'Enter name of relaxed topography file.'
  read(*,*) rtopofile
  inquire(file=rtopofile,exist=NC%do_var(NC_B_RELX))
  write(*,*) 'Enter name of standard deviation file.'
  read(*,*) sdfile
  inquire(file=sdfile,exist=NC%do_var(NC_B_STD_DEV))
  write(*,*) 'Enter name of surface file.'
  read(*,*) surffile
  inquire(file=surffile,exist=NC%do_var(NC_B_USURF))
  write(*,*) 'Enter name of present precipitation file.'
  read(*,*) precepfile
  inquire(file=precepfile,exist=NC%do_var(NC_B_PRESPRCP))
  write(*,*) 'Enter name of present ice surface file.'
  read(*,*) pressurf
  inquire(file=pressurf,exist=NC%do_var(NC_B_PRESUSRF))
  write(*,*) 'Enter name of interface mask file.'
  read(*,*) mask
  inquire(file=mask,exist=NC%do_var(NC_B_MASK))

  write(*,*) 'Enter name of output netCDF file.'
  read(*,*) NC%filename

  ! Define a projection 

  write(*,*) 'Enter the projection type:'
  write(*,*) CFP_LAEA,') Lambert Azimuthal Equal-Area'
  write(*,*) CFP_AEA ,') Albers Equal Area Conic'
  write(*,*) CFP_LCC ,') Lambert Conic Conformal'
  write(*,*) CFP_STERE,') Stereographic'
  read(*,*) ptype

  write(*,*) 'Enter longitude of central meridian:'
  read(*,*) longitude_of_central_meridian

  write(*,*) 'Enter latitude of projection origin:'
  read(*,*) latitude_of_projection_origin

  write(*,*) 'Enter false easting:'
  read(*,*) false_easting

  write(*,*) 'Enter false northing:'
  read(*,*) false_northing

  if (ptype.ne.CFP_LAEA) then
     write(*,*) 'Enter standard parallel'
     read(*,*) standard_parallel
  end if

  if (ptype.eq.CFP_AEA.or.ptype.eq.CFP_LCC) then
     write(*,*) 'Enter standard parallel 2'
     read(*,*) standard_parallel_2
  end if

  if (ptype.eq.CFP_STERE) then
     write(*,*) 'Polar? (true or false)'
     read(*,*) polar
     write(*,*) 'Enter scale factor at projection origin'
     read(*,*) scale_factor_at_proj_origin
  end if

  call CFproj_define(model%projection, &
       ptype, &
       polar, &
       longitude_of_central_meridian, &
       latitude_of_projection_origin, &
       scale_factor_at_proj_origin, &
       standard_parallel, &
       standard_parallel_2, &
       false_easting, &
       false_northing)

  if (.not.NC%do_var(NC_B_TOPG)) then
     write(*,*) 'No topo file, bailing out...'
     stop
  end if

  call readplan(topofile,data,time,nx,ny,delta)
  delta = delta/len0
  model%general%upn = 1
  model%numerics%dew = delta
  model%numerics%dns = delta
  allocate(model%numerics%sigma(model%general%upn))
  model%numerics%sigma = 1
  model%general%ewn = nx
  model%general%nsn = ny
  call glimmer_nc_createfile(outfile, model)

  status = nf90_put_var(NC%id, NC%varids(NC_B_TOPG), data, (/1,1,1/))
  call nc_errorhandle(__FILE__,__LINE__,status)

  status = nf90_put_var(NC%id, NC%timevar,time,(/1/))
  call nc_errorhandle(__FILE__,__LINE__,status)

  if (NC%do_var(NC_B_LAT)) then
     call readplan(latfile,data,time,nx,ny,delta)
     status = nf90_put_var(NC%id, NC%varids(NC_B_LAT), data, (/1,1,1/))
     call nc_errorhandle(__FILE__,__LINE__,status)
  end if

  if (NC%do_var(NC_B_RELX)) then
     call readplan(rtopofile,data,time,nx,ny,delta)
     status = nf90_put_var(NC%id, NC%varids(NC_B_RELX), data, (/1,1,1/))
     call nc_errorhandle(__FILE__,__LINE__,status)
  end if

  if (NC%do_var(NC_B_USURF)) then
     call readplan(surffile,data,time,nx,ny,delta)
     status = nf90_put_var(NC%id, NC%varids(NC_B_USURF), data, (/1,1,1/))
     call nc_errorhandle(__FILE__,__LINE__,status)
  end if

  if (NC%do_var(NC_B_PRESPRCP)) then
     call readplan(precepfile,data,time,nx,ny,delta)
     status = nf90_put_var(NC%id, NC%varids(NC_B_PRESPRCP), data, (/1,1,1/))
     call nc_errorhandle(__FILE__,__LINE__,status)
  end if

  if (NC%do_var(NC_B_PRESUSRF)) then
     call readplan(pressurf,data,time,nx,ny,delta)
     status = nf90_put_var(NC%id, NC%varids(NC_B_PRESUSRF), data, (/1,1,1/))
     call nc_errorhandle(__FILE__,__LINE__,status)
  end if

  if (NC%do_var(NC_B_STD_DEV)) then
     call readplan(sdfile,data,time,nx,ny,delta)
     status = nf90_put_var(NC%id, NC%varids(NC_B_STD_DEV), data, (/1,1,1/))
     call nc_errorhandle(__FILE__,__LINE__,status)
  end if

  if (NC%do_var(NC_B_MASK)) then
     call maskread(mask,maskdata,nx,ny)
     status = nf90_put_var(NC%id, NC%varids(NC_B_MASK), maskdata, (/1,1,1/))
     call nc_errorhandle(__FILE__,__LINE__,status)
  end if

  status = nf90_close(NC%id)
  call close_log
contains
  subroutine readplan(fname,data,time,ewn,nsn, delta)
    use glimmer_global, only : sp
    use glide_messages
    implicit none
    character(*),intent(in) :: fname
    real(kind=sp), dimension(:,:), pointer :: data
    integer, intent(out) :: ewn,nsn
    real(kind=sp), intent(out) :: time, delta

    ! local 
    integer :: unit=1
    character :: cdum*4 
    logical :: here
    integer i,j
    character(40) :: errtxt

    if (associated(data)) then
       deallocate(data)
    end if

    inquire(file=fname,exist=here)
    if ( .not. here ) then
      call write_log(GM_FATAL,'planform file '//trim(fname)//' not found',__FILE__,__LINE__)
    endif

#ifdef CVF
    open(unit,file=fname,form='unformatted',convert='BIG_ENDIAN')
    ! convert is a non-standard specifier and doesn't work with the
    ! Intel compiler.
    ! Substituting with

#else
    open(unit,file=fname,form='unformatted')
#endif

    read(unit) time, cdum, ewn,nsn, delta
    allocate(data(ewn,nsn))
    read(unit) ((data(i,j),i=1,ewn),j=1,nsn)

    close(unit)
  end subroutine readplan

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine maskread(filename,data,nx,ny)

    use glimmer_global, only : sp

    implicit none

    !*FD Reads a binary integer mask file

    character(*),intent(in) :: filename
    integer, dimension(:,:), pointer :: data
    integer, intent(in) :: nx,ny

	integer      :: unit=1    
    logical      :: there
    integer      :: i,j
    real(sp)     :: rdum
    character(4) :: cdum
    integer      :: idum

    if (associated(data)) then
       deallocate(data)
    end if
	allocate(data(nx,ny))

    inquire(file=filename,exist=there)
    if ( .not. there ) then
      call write_log(GM_FATAL,'mask file '//trim(filename)//' not found',__FILE__,__LINE__)
    endif

#ifdef CVF
    open(unit,file=filename,form='unformatted',convert='BIG_ENDIAN')
    ! convert is a non-standard specifier and doesn't work with the
    ! Intel compiler.
    ! Substituting with

#else

    open(unit,file=filename,form='unformatted')

#endif

    read(unit,end=10) rdum, cdum, idum, rdum 
    read(unit,end=10) ((data(i,j),i=1,nx),j=1,ny)

    close(unit)

    return

10    print*, 'read past end of mask file'
    Print*,'Reading file:',filename
    stop

  end subroutine maskread

end program input2ncdf
