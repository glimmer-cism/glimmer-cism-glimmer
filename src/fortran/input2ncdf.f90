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

program input2ncdf
  !*FD convert glimmer input fields to netCDF
  !*FD written by Magnus Hagdorn, May 2004
  use glimmer_ncdf
  use glimmer_ncfile
  use glimmer_types
  implicit none

  ! file names
  character(len=30) :: latfile, topofile, rtopofile, surffile, precepfile, pressurf
  integer unit
  type(glimmer_nc_output), pointer :: outfile
  type(glimmer_global_type) :: model

  
  write(*,*) 'Enter name of latitude file.'
  read(*,*) latfile
  write(*,*) 'Enter name of topography file.'
  read(*,*) topofile
  write(*,*) 'Enter name of relaxed topography file.'
  read(*,*) rtopofile
  write(*,*) 'Enter name of surface file file.'
  read(*,*) surffile
  write(*,*) 'Enter name of present precipitation file.'
  read(*,*) precepfile
  write(*,*) 'Enter name of present ice surface file.'
  read(*,*) pressurf

  allocate(outfile)
  write(*,*) 'Enter name of output netCDF file.'
  read(*,*) NC%filename


end program input2ncdf
