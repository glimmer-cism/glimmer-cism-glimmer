
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_example.f90 - part of the GLIMMER ice model      + 
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

program glimmer_example

!*FD This program demonstrates the use of GLIMMER. It loads in
!*FD some example global fields and associated grid data,
!*FD Initialises the model, and then runs it for 1000 years.

  use glimmer_main

  implicit none

  type(glimmer_params) :: ice_sheet
  integer :: nx,ny,i,j
  character(fname_length) :: paramfile='example.glp'
  real(rk) :: time
  real(rk),dimension(:,:),allocatable :: temp,precip,zonwind,merwind,orog,coverage
  real(rk),dimension(:,:),allocatable :: albedo, orog_out, ice_frac, fw
  real(rk),dimension(:),allocatable :: lats,lons
  logical :: out
  
  nx=64 ; ny=32

  allocate(temp(nx,ny),precip(nx,ny),zonwind(nx,ny))
  allocate(merwind(nx,ny),lats(ny),lons(nx),orog(nx,ny))
  allocate(coverage(nx,ny),orog_out(nx,ny),albedo(nx,ny),ice_frac(nx,ny),fw(nx,ny))

  temp=0.0
  precip=0.0
  zonwind=0.0
  merwind=0.0
  lats=0.0
  lons=0.0
  orog=0.0
  albedo=0.0
  orog_out=0.0

! Load in latitudes

  open(20,file='latitudes.txt')
  do j=1,ny
    read(20,100)lats(j)
  enddo
  close(20)

! Load in longitudes

  open(20,file='longitudes.txt')
  do i=1,nx
    read(20,100)lons(i)
  enddo
  close(20)
 
! Load in sample temperature field

  open(20,file='tempfield.dat')
  do i=1,nx
    do j=1,ny
      read(20,101)temp(i,j)
    enddo
  enddo
  close(20)

! Load in sample precip field

  open(20,file='precipfield.dat')
  do i=1,nx
    do j=1,ny
      read(20,101)precip(i,j)
    enddo
  enddo
  close(20)

! Load in sample zonal wind field

  open(20,file='zonwindfield.dat')
  do i=1,nx
    do j=1,ny
      read(20,101)zonwind(i,j)
    enddo
  enddo
  close(20)

! Load in sample meridional wind field

  open(20,file='merwindfield.dat')
  do i=1,nx
    do j=1,ny
      read(20,101)merwind(i,j)
    enddo
  enddo
  close(20)

! Load large-scale orography

  open(20,file='largescale.orog')
  do i=1,nx
    do j=1,ny
      read(20,101)orog(i,j)
    enddo
  enddo
  close(20)

  precip=precip*24*365   ! Convert background precip rate from mm/h to mm/a

  time=24.0*360.0

  call initialise_glimmer(ice_sheet,lats,lons,paramfile)

!  call glimmer_read_restart(ice_sheet,25,'glimmer.restart')

! Set all inputs to zero...

!  temp=0.0
!  precip=0.0
!  zonwind=0.0
!  merwind=0.0

  if (glimmer_coverage_map(ice_sheet,coverage).ne.0) then
    Print*,'unable to get coverage map'
    stop
  endif

  do time=0.0*24.0*360.0,10000.0*24.0*360.0,24.0*360.0
    call glimmer(ice_sheet,time,temp,precip,zonwind,merwind,orog, &
                 orog_out=orog_out,albedo=albedo,output_flag=out, &
                 ice_frac=ice_frac,fw_flux=fw) 
    if(out) orog=coverage*orog_out+(1-coverage)*orog
  enddo
  call end_glimmer(ice_sheet)

! call glimmer_write_restart(ice_sheet,20,'glimmer.restart')

100 format(f9.5)
101 format(e12.5)

end program glimmer_example
