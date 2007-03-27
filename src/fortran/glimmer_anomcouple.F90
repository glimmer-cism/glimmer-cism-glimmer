! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_anomcouple.f90 - part of the GLIMMER ice model   + 
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
#include <config.inc>
#endif

module glimmer_anomcouple

  !*FD This module provides code to handle anomaly coupling. Although
  !*FD written for use with GLINT, it has general applicability

  use glimmer_global

  implicit none

  type anomaly_coupling
     character(fname_length) :: fname_reference !*FD File containing reference climate
     character(fname_length) :: fname_modelclim !*FD File containing mean model climate
     integer :: nslices                         !*FD Number of time-slices in climatologies
     real(dp),dimension(:,:,:),pointer :: temp_ref !*FD Reference climate (temperature)
     real(dp),dimension(:,:,:),pointer :: temp_mod !*FD Model climate (temperature)
     real(dp),dimension(:,:,:),pointer :: prcp_ref !*FD Reference climate (precip)
     real(dp),dimension(:,:,:),pointer :: prcp_mod !*FD Model climate (precip)
     real(dp),dimension(:)    ,pointer :: time     !*FD Time axis (fraction of year)
  end type anomaly_coupling

end module glimmer_anomcouple
