! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  eismint3_forcing.f90 - part of the GLIMMER ice model     + 
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

module eismint3_forcing

  use glimmer_global, only: sp
  use eismint3_types

  implicit none

  !*FD Provides climate forcing for EISMINT3 Greeland scenario 1

contains

  subroutine eismint3_initialise(climate,config,model)

    !*FD initialise EISMINT3 climate

    use glide_types
    use glimmer_config
    use eismint3_io
    use glide_io
    use glide_setup

    implicit none

    type(eismint3_climate) :: climate      !*FD structure holding EISMINT3 climate
    type(ConfigSection), pointer :: config !*FD structure holding sections of configuration file   
    type(glide_global_type) :: model       !*FD model instance

    call coordsystem_allocate(model%general%ice_grid,climate%prcp)
    call coordsystem_allocate(model%general%ice_grid,climate%artm)
    call coordsystem_allocate(model%general%ice_grid,climate%arng)
    call coordsystem_allocate(model%general%ice_grid,climate%acab)
    call coordsystem_allocate(model%general%ice_grid,climate%usrf)
    call coordsystem_allocate(model%general%ice_grid,climate%ablt)
    call coordsystem_allocate(model%general%ice_grid,climate%presusurf)
    call coordsystem_allocate(model%general%ice_grid,climate%landsea)

    call eismint3_io_readall(climate,model)
    call glimmer_pdd_init(climate%pdd_scheme,config)
 
    ! Calculate initial thickness.

    where (climate%presusurf>0.0)
       climate%landsea=.true.
    elsewhere
       climate%landsea=.false.
    end where

    call eismint3_temp(climate%artm,climate%arng,climate%presusurf,model%climate%lati)
    call glimmer_pdd_mbal(climate%pdd_scheme,climate%artm,climate%arng,climate%prcp,climate%ablt,climate%acab,climate%landsea)

    ! Put it into glide

    call glide_set_thk(model,max(0.0,climate%acab*model%numerics%tinc))
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus,model%geometry%lsrf)
    model%geometry%usrf = model%geometry%thck + model%geometry%lsrf

  end subroutine eismint3_initialise

  subroutine eismint3_clim(climate,model)

    use glide_types
    use glide_io

    type(eismint3_climate) :: climate      !*FD structure holding EISMINT3 climate
    type(glide_global_type) :: model       !*FD model instance

    call glide_get_usurf(model,climate%usrf)
    where (climate%usrf>0.0)
       climate%landsea=.true.
    elsewhere
       climate%landsea=.false.
    end where

    call eismint3_temp(climate%artm,climate%arng,climate%usrf,model%climate%lati)
    call glimmer_pdd_mbal(climate%pdd_scheme,climate%artm,climate%arng,climate%prcp,climate%ablt,climate%acab,climate%landsea)

    where (.not.climate%landsea) climate%acab=0.0

    call glide_set_acab(model,climate%acab)
    call glide_set_artm(model,climate%artm)

  end subroutine eismint3_clim

  subroutine eismint3_temp(artm,arng,usrf,lati)

    real(sp),dimension(:,:),intent(out) :: artm,arng
    real(sp),dimension(:,:),intent(in)  :: usrf,lati

    artm=49.13-0.007992*max(usrf,20*(lati-65.0))-0.7576*lati
    arng=30.78-0.006277*usrf-0.3262*lati-artm

  end subroutine eismint3_temp

end module eismint3_forcing
