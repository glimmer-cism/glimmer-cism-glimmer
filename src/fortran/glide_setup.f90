! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_setup.f90 - part of the GLIMMER ice model        + 
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

module glide_setup

  !*FD Contains general routines for initialisation, etc, called
  !*FD from the top-level glimmer subroutines.

  private
  public :: readconfig, printconfig

contains

  subroutine readconfig(model,fname)
    !*FD read GLIDE configuration file
    use glide_types
    use glimmer_config
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    character(len=*), intent(in) :: fname   !*FD name of paramter file
    
    ! local variables
    type(ConfigSection), pointer :: config
    type(ConfigSection), pointer :: section
    type(glimmer_nc_output), pointer :: output => null()
    type(glimmer_nc_input), pointer :: input => null()

    model%funits%ncfile = fname

    ! read configuration
    call ConfigRead(fname,config)

    ! read grid size  parameters
    call GetSection(config,section,'grid')
    if (associated(section)) then
       call handle_grid(section, model)
    end if
    ! read time parameters
    call GetSection(config,section,'time')
    if (associated(section)) then
       call handle_time(section, model)
    end if
    ! read options parameters
    call GetSection(config,section,'options')
    if (associated(section)) then
       call handle_options(section, model)
    end if
    ! read parameters
    call GetSection(config,section,'parameters')
    if (associated(section)) then
       call handle_parameters(section, model)
    end if

  end subroutine readconfig


  subroutine printconfig(unit,model)
    !*FD print model configuration to unit
    use glide_types
    implicit none
    integer, intent(in) :: unit       !*FD unit to print to
    type(glide_global_type)  :: model !*FD model instance

    call print_grid(unit,model)
    call print_time(unit,model)
    call print_options(unit,model)
    call print_parameters(unit,model)
  end subroutine printconfig
    

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! private procedures
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! grid sizes
  subroutine handle_grid(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'ewn',model%general%ewn)
    call GetValue(section,'nsn',model%general%nsn)
    call GetValue(section,'upn',model%general%upn)
    call GetValue(section,'dew',model%numerics%dew)
    call GetValue(section,'dns',model%numerics%dns)
    call GetValue(section,'sigma_file',model%funits%sigfile)
  end subroutine handle_grid

  subroutine print_grid(unit,model)
    use glide_types
    implicit none
    integer, intent(in) :: unit
    type(glide_global_type)  :: model

    write(unit,*) 'Grid specification'
    write(unit,*) '------------------'
    write(unit,*) 'ewn             : ',model%general%ewn
    write(unit,*) 'nsn             : ',model%general%nsn
    write(unit,*) 'upn             : ',model%general%upn
    write(unit,*) 'EW grid spacing : ',model%numerics%dew
    write(unit,*) 'NS grid spacing : ',model%numerics%dns
    write(unit,*) 'sigma file      : ',trim(model%funits%sigfile)
    write(unit,*)
  end subroutine print_grid

  ! time
  subroutine handle_time(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'dt',model%numerics%tinc)
    call GetValue(section,'ntem',model%numerics%ntem)
    call GetValue(section,'nvel',model%numerics%nvel)
    call GetValue(section,'niso',model%numerics%niso)
  end subroutine handle_time
  
  subroutine print_time(unit,model)
    use glide_types
    implicit none
    integer, intent(in) :: unit
    type(glide_global_type)  :: model

    write(unit,*) 'Time steps'
    write(unit,*) '----------'
    write(unit,*) 'main time step    : ',model%numerics%tinc
    write(unit,*) 'thermal dt factor : ',model%numerics%ntem
    write(unit,*) 'velo dt factor    : ',model%numerics%nvel
    write(unit,*) 'isostasy dt factor: ',model%numerics%niso
    write(unit,*) ''
  end subroutine print_time

  ! options
  subroutine handle_options(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'ioparams',model%funits%ncfile)
    call GetValue(section,'temperature',model%options%whichtemp)
    call GetValue(section,'flow_law',model%options%whichflwa)
    call GetValue(section,'isostasy',model%options%whichisot)
    call GetValue(section,'sliding_law',model%options%whichslip)
    call GetValue(section,'basal_water',model%options%whichbwat)
    call GetValue(section,'marine_margin',model%options%whichmarn)
    call GetValue(section,'slip_coeff',model%options%whichbtrc)
    call GetValue(section,'stress_calc',model%options%whichstrs)
    call GetValue(section,'evolution',model%options%whichevol)
    call GetValue(section,'vertical_integration',model%options%whichwvel)
  end subroutine handle_options
  
  subroutine print_options(unit,model)
    use glide_types
    implicit none
    integer, intent(in) :: unit
    type(glide_global_type)  :: model

    ! local variables
    character(len=*), dimension(0:1), parameter :: temperature = (/ &
         'isothermal', &
         'full      '/)
    character(len=*), dimension(0:2), parameter :: flow_law = (/ &
         'Patterson and Budd               ', &
         'Patterson and Budd (temp=-10degC)', &
         'const 1e-16a^-1Pa^-n             ' /)
    character(len=*), dimension(0:2), parameter :: isostasy = (/ &
         'none   ', &
         'local  ', &
         'elastic' /)
    character(len=*), dimension(0:4), parameter :: sliding = (/ &
         'gravity', &
         'unknown', &
         'unknown', &
         'unknown', &
         'zero   ' /)
    character(len=*), dimension(0:2), parameter :: basal_water = (/ &
         'local water balance', &
         'local + const flux ', &
         'none               ' /)
    character(len=*), dimension(0:2), parameter :: marine_margin = (/ &
         'threshold   ', &
         'no ice shelf', &
         'none        ' /)
    character(len=*), dimension(0:1), parameter :: slip_coeff = (/ &
         '~basal water', &
         'zero        ' /)
    character(len=*), dimension(0:3), parameter :: stress = (/ &
         'zeroth-order                     ', &
         'first-order                      ', &
         'vertically-integrated first-order',&
         'none                             ' /)
    character(len=*), dimension(0:2), parameter :: evolution = (/ &
         'pseudo-diffusion', &
         'unknown         ', &
         'diffusion       ' /)
    character(len=*), dimension(0:1), parameter :: vertical_integration = (/ &
         'standard     ', &
         'obey upper BC' /)


    write(unit,*) 'GLIDE options'
    write(unit,*) '-------------'
    write(*,*) 'I/O parameter file      : ',trim(model%funits%ncfile)
    if (model%options%whichtemp.lt.0 .or. model%options%whichtemp.gt.size(temperature)) then
       write(*,*) 'Error, temperature out of range'
       stop
    end if
    write(unit,*) 'temperature calculation : ',model%options%whichtemp,temperature(model%options%whichtemp)
    if (model%options%whichflwa.lt.0 .or. model%options%whichflwa.gt.size(flow_law)) then
       write(*,*) 'Error, flow_law out of range'
       stop
    end if
    write(unit,*) 'flow law                : ',model%options%whichflwa,flow_law(model%options%whichflwa)
    if (model%options%whichisot.lt.0 .or. model%options%whichisot.gt.size(isostasy)) then
       write(*,*) 'Error, isostasy out of range'
       stop
    end if
    write(unit,*) 'isostasy                : ',model%options%whichisot,isostasy(model%options%whichisot)
    if (model%options%whichslip.lt.0 .or. model%options%whichslip.gt.size(sliding)) then
       write(*,*) 'Error, sliding_law out of range'
       stop
    end if
    write(unit,*) 'sliding_law             : ',model%options%whichslip, sliding(model%options%whichslip)
    if (model%options%whichbwat.lt.0 .or. model%options%whichbwat.gt.size(basal_water)) then
       write(*,*) 'Error, basal_water out of range'
       stop
    end if
    write(unit,*) 'basal_water             : ',model%options%whichbwat,basal_water(model%options%whichbwat)
    if (model%options%whichmarn.lt.0 .or. model%options%whichmarn.gt.size(marine_margin)) then
       write(*,*) 'Error, marine_margin out of range'
       stop
    end if
    write(unit,*) 'marine_margin           : ', model%options%whichmarn, marine_margin(model%options%whichmarn)
    if (model%options%whichbtrc.lt.0 .or. model%options%whichbtrc.gt.size(slip_coeff)) then
       write(*,*) 'Error, slip_coeff out of range'
       stop
    end if
    write(unit,*) 'slip_coeff              : ', model%options%whichbtrc, slip_coeff(model%options%whichbtrc)
    if (model%options%whichstrs.lt.0 .or. model%options%whichstrs.gt.size(stress)) then
       write(*,*) 'Error, stress_calc out of range'
       stop
    end if
    write(unit,*) 'stress_calc             : ', model%options%whichstrs, stress(model%options%whichstrs)
    if (model%options%whichevol.lt.0 .or. model%options%whichevol.gt.size(evolution)) then
       write(*,*) 'Error, evolution out of range'
       stop
    end if
    write(unit,*) 'evolution               : ', model%options%whichevol, evolution(model%options%whichevol)
    if (model%options%whichwvel.lt.0 .or. model%options%whichwvel.gt.size(vertical_integration)) then
       write(*,*) 'Error, vertical_integration out of range'
       stop
    end if
    write(unit,*) 'vertical_integration    : ',model%options%whichwvel,vertical_integration(model%options%whichwvel)
    write(unit,*) ''
  end subroutine print_options

  ! parameters
  subroutine handle_parameters(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model
    real, pointer, dimension(:) :: temp => NULL()

    call GetValue(section,'ice_limit',model%numerics%thklim)
    call GetValue(section,'marine_limit',model%numerics%mlimit)
    call GetValue(section,'geothermal',model%paramets%geot)
    call GetValue(section,'flow_factor',model%paramets%fiddle)
    call GetValue(section,'hydro_time',model%paramets%hydtim)
    call GetValue(section,'isos_time',model%paramets%isotim)
    call GetValue(section,'basal_tract',temp,5)
    if (associated(temp)) then
       model%paramets%bpar = temp
       deallocate(temp)
    end if
  end subroutine handle_parameters

  subroutine print_parameters(unit,model)
    use glide_types
    implicit none
    integer, intent(in) :: unit
    type(glide_global_type)  :: model

    write(unit,*) 'Parameters'
    write(unit,*) '----------'
    write(unit,*) 'ice limit             : ',model%numerics%thklim
    write(unit,*) 'marine depth limit    : ',model%numerics%mlimit
    write(unit,*) 'geothermal heat flux  : ',model%paramets%geot
    write(unit,*) 'flow enhancement      : ',model%paramets%fiddle
    write(unit,*) 'basal hydro time const: ',model%paramets%hydtim
    write(unit,*) 'isostasy time const   : ',model%paramets%isotim
    write(unit,*) 'basal tractio factors : ',model%paramets%bpar(1)
    write(unit,*) '                        ',model%paramets%bpar(2)
    write(unit,*) '                        ',model%paramets%bpar(3)
    write(unit,*) '                        ',model%paramets%bpar(4)
    write(unit,*) '                        ',model%paramets%bpar(5)
    write(unit,*)
  end subroutine print_parameters


end module glide_setup
