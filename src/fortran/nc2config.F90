! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  nc2config.f90 - part of the GLIMMER ice model            + 
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

program nc2config

  !*FD Program to extract the config file data from a
  !*FD glimmer netcdf output file, and write it out, either
  !*FD to file or the screen. Command-line argument and 
  !*FD option support is provided, where available

  use netcdf
  use glimmer_global, only: endline

  implicit none

  character(100) :: infile, outfile, opt1, opt2, opt3
  logical :: stdout = .true.
  character(10000) :: config
  integer :: status,ncid,attlen,unit,i,ellen

#ifdef COMMAND_LINE

  if (nargs().eq.2) then
     call getarg(1,opt1)
     if (opt1=='-h') then
        call usage
        stop
     else
        infile=opt1
     end if
  else if (nargs().eq.4) then
     call getarg(1,opt1)
     call getarg(2,opt2)
     call getarg(3,opt3)
     if (opt1=='-o') then
        outfile=opt2
        infile=opt3
        stdout=.false.
     else if (opt2=='-o') then
        outfile=opt3
        infile=opt1
        stdout=.false.
     else
        write(0,*)'ERROR: Unknown argument combination'
        call usage
        stop
     end if
  else if (nargs().eq.1) then
#endif

     ! These are the default inputs
     Print*,'Enter name of input file:'
     read*,infile
     Print*,'Enter name of output file:'
     read*,outfile
     stdout=.false.

#ifdef COMMAND_LINE
  else
     write(0,*)'ERROR: Unknown argument combination'
     call usage
     stop
  endif
#endif

  ! Open file and look for appropriate attribute

  status=nf90_open(infile,0,ncid)
  if (status/=NF90_NOERR) call netcdf_error(status)
  status=nf90_inquire_attribute(ncid,NF90_GLOBAL,'configuration',len=attlen)
  if (status==NF90_ENOTATT) then
     write(0,*)'ERROR: No configuration data found in file ',trim(infile)
     stop
  else if (status/=NF90_NOERR) then
     call netcdf_error(status)
  end if
  status=nf90_get_att(ncid,NF90_GLOBAL,'configuration',config)
  if (status/=NF90_NOERR) call netcdf_error(status)
  config=config(:attlen)

  ! Format and print output

  ellen=len(endline)

  if (stdout) then
     unit=6
  else
     unit=20
     open(unit,file=outfile)
  end if

  i=1
  do 
     if (trim(config(i:))=='') exit
     if (config(i:i+ellen-1)==endline) then
        write(unit,'(A)')''
        i=i+ellen
     else if (config(i:i)=='[') then
        write(unit,'(A)')''
        write(unit,'(A)',advance='no')config(i:i)
        i=i+1
     else
        write(unit,'(A)',advance='no')config(i:i)
        i=i+1
     end if
  end do

contains

  subroutine usage

    print*,'usage: nc2config [-h | [-o outfile] [infile]]'
    print*,'    -h: display this message'
    print*,'    -o: specify output file. If absent, stdout is used'

  end subroutine usage

  subroutine netcdf_error(status)

    integer :: status

    print*,nf90_strerror(status)
    stop

  end subroutine netcdf_error

end program nc2config
