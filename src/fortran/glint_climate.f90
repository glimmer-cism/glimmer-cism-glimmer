! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_climate.f90 - part of the GLIMMER ice model        + 
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

module glint_climate
  !*FD handle glint climate

contains

  subroutine calcartm(usrf,artm,arng,g_orog,g_artm,g_arng,lapse_rate)

    !*FD Calculate the surface air temperature and mean annual range,
    !*FD according to various different models.

    use glimmer_global, only : dp, sp 
    use glint_mbal
    implicit none

    real(sp),dimension(:,:),         intent(in)    :: usrf   !*FD Surface elevation (km)
    real(sp),dimension(:,:),         intent(out)   :: artm   !*FD Surface annual mean air temperature ($^{\circ}$C)
    real(sp),dimension(:,:),         intent(out)   :: arng   !*FD Surface annual air tempurature range ($^{\circ}$C)
    real(dp),dimension(:,:),         intent(in)    :: g_orog !*FD Global orography on local grid (m)
    real(sp),dimension(:,:),         intent(in)    :: g_artm !*FD Supplied global air temperatures ($^{\circ}$C)
    real(sp),dimension(:,:),         intent(in)    :: g_arng !*FD Supplied global air temp range ($^{\circ}$C)
    real(sp),                        intent(in)    :: lapse_rate !*FD Lapse rate for temperature correction

    ! Copy the fields
    
    arng=g_arng
    artm=g_artm
    
    ! Reduce temperatures to sea-level
    
    call glint_lapserate(artm,real(g_orog,rk),real(-lapse_rate,rk))
    
    ! Raise them to high-res orography 
    
    call glint_lapserate(artm,real(usrf,rk),real(lapse_rate,rk))
                
  end subroutine calcartm

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine fix_acab(ablt,acab,prcp,thck,usrf)

    use glimmer_global, only : sp 

    real(sp),dimension(:,:),intent(inout) :: ablt
    real(sp),dimension(:,:),intent(inout) :: acab
    real(sp),dimension(:,:),intent(in)    :: prcp
    real(sp),dimension(:,:),intent(in)    :: thck
    real(sp),dimension(:,:),intent(in)    :: usrf

    integer  :: nx,ny

    nx=size(prcp,1) ; ny=size(prcp,2)

    ! Adjust ablation to be no greater than ice available for melting

    where (ablt>thck) 
      ablt=thck+prcp
      acab=prcp-ablt
    endwhere

    ! If the upper ice/land surface is at or below sea-level, set accumulation,
    ! ablation and mass-balance to zero. This is to prevent accumulation of ice below
    ! sea-level.

    where (usrf<=0.0)
      ablt=prcp
      acab=0.0
    end where

    ! Remove accumulation from domain edges

    ablt(:,1)=prcp(:,1)
    acab(:,1)=0.0
    ablt(:,ny)=prcp(:,ny)
    acab(:,ny)=0.0
    ablt(1,:)=prcp(1,:)
    acab(1,:)=0.0
    ablt(nx,:)=prcp(nx,:)
    acab(nx,:)=0.0

  end subroutine fix_acab

end module glint_climate
