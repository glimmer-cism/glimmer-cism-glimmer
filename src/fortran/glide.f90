! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide.f90 - part of the GLIMMER ice model                + 
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

module glide
  !*FD the top-level GLIDE module

  use glide_types

  integer, private, parameter :: dummyunit=99

contains
  subroutine glide_initialise(model,fname)
    !*FD initialise GLIDE model instance
    use glide_setup
    use glimmer_ncparams
    use glimmer_ncfile
    use glimmer_ncinfile
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    character(len=*), intent(in) :: fname   !*FD name of paramter file

    ! read configuration file
    call glide_readconfig(model,fname)
    call glide_printconfig(6,model)
    ! scale parameters
    call glide_scale_params(model)
    ! allocate arrays
    call glide_allocarr(model)
    
    ! load sigma file
    call glide_load_sigma(model,dummyunit)

    ! netCDF I/O
    call ReadNCParams(model)
    ! open all input files
    call openall_in(model)
    ! and read first time slice
    call readall(model)

    ! open all output files
    call openall_out(model)
  end subroutine glide_initialise
  
  subroutine glide_finalise(model)
    !*FD finalise GLIDE model instance
    use glimmer_ncfile
    use glimmer_ncinfile
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    
    call closeall_in(model)
    call closeall_out(model)
    
    call glide_deallocarr(model)
  end subroutine glide_finalise
end module glide
