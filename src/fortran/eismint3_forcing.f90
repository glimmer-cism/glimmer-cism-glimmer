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
    call coordsystem_allocate(model%general%ice_grid,climate%landsea)

    call eismint3_io_readall(climate,model)
    call glimmer_pdd_init(climate%pdd_scheme,config)
 
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

    climate%artm=49.13-0.007992*max(climate%usrf,20*(model%climate%lati-65.0))-0.7576*model%climate%lati
    climate%arng=30.38-0.006277*climate%usrf-0.3262*model%climate%lati-climate%artm

    call glimmer_pdd_mbal(climate%pdd_scheme,climate%artm,climate%arng,climate%prcp,climate%ablt,climate%acab,climate%landsea)

    where (.not.climate%landsea) climate%acab=0.0

    call glide_set_acab(model,climate%acab)
    call glide_set_artm(model,climate%artm)

  end subroutine eismint3_clim

end module eismint3_forcing
