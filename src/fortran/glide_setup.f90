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
  public :: glide_readconfig, glide_printconfig, glide_scale_params, &
       glide_calclsrf, glide_marinlim, glide_load_sigma, glide_maskthck

contains

  subroutine glide_readconfig(model,config)
    !*FD read GLIDE configuration file
    use glide_types
    use glimmer_config
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file
    
    ! local variables
    type(ConfigSection), pointer :: section
    type(glimmer_nc_output), pointer :: output => null()
    type(glimmer_nc_input), pointer :: input => null()

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
  end subroutine glide_readconfig

  subroutine glide_printconfig(model)
    !*FD print model configuration to log
    use glimmer_log
    use glide_types
    implicit none
    type(glide_global_type)  :: model !*FD model instance

    call write_log_div
    call print_grid(model)
    call print_time(model)
    call print_options(model)
    call print_parameters(model)
  end subroutine glide_printconfig
    
  subroutine glide_scale_params(model)
    !*FD scale parameters
    use glide_types
    use physcon,  only: scyr, gn
    use paramets, only: thk0,tim0,len0, tau0, vel0, vis0
    implicit none
    type(glide_global_type)  :: model !*FD model instance

    tau0 = (vel0/(vis0*len0))**(1.0/gn)

    model%numerics%ntem = model%numerics%ntem * model%numerics%tinc
    model%numerics%nvel = model%numerics%nvel * model%numerics%tinc
    model%numerics%niso = model%numerics%niso * model%numerics%tinc

    model%numerics%dt     = model%numerics%tinc * scyr / tim0
    model%numerics%dttem  = model%numerics%ntem * scyr / tim0 
    model%numerics%thklim = model%numerics%thklim  / thk0

    model%numerics%dew = model%numerics%dew / len0
    model%numerics%dns = model%numerics%dns / len0

    model%paramets%isotim = model%paramets%isotim * scyr / tim0         
    model%numerics%mlimit = model%numerics%mlimit / thk0
  end subroutine glide_scale_params

!-------------------------------------------------------------------------

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glide_maskthck(whichthck,crita,critb,dom,pointno,totpts,empty)
    
    !*FD Calculates the contents of the mask array.

    use glimmer_global, only : dp, sp 

    implicit none

    !-------------------------------------------------------------------------
    ! Subroutine arguments
    !-------------------------------------------------------------------------

    integer,                intent(in)  :: whichthck  !*FD Option determining
                                                      !*FD which method to use.
    real(dp),dimension(:,:),intent(in)  :: crita      !*FD Ice thickness
    real(sp),dimension(:,:),intent(in)  :: critb      !*FD Mass balance
    integer, dimension(:),  intent(out) :: dom        
    integer, dimension(:,:),intent(out) :: pointno    !*FD Output mask
    integer,                intent(out) :: totpts     !*FD Total number of points
    logical,                intent(out) :: empty      !*FD Set if no mask points set.

    !-------------------------------------------------------------------------
    ! Internal variables
    !-------------------------------------------------------------------------

    integer,dimension(size(crita,2),2) :: band
    logical,dimension(size(crita,2))   :: full
    integer :: covtot 
    integer :: ew,ns,ewn,nsn

    !-------------------------------------------------------------------------

    ewn=size(crita,1) ; nsn=size(crita,2)

    pointno = 0
    covtot  = 0 

    !-------------------------------------------------------------------------

    do ns = 1,nsn

      full(ns) = .false.

      do ew = 1,ewn
        if ( thckcrit(crita(max(1,ew-1):min(ewn,ew+1),max(1,ns-1):min(nsn,ns+1)),critb(ew,ns)) ) then

          covtot = covtot + 1
          pointno(ew,ns) = covtot 

          if ( .not. full(ns) ) then
            band(ns,1) = ew
            full(ns)   = .true.
          else
            band(ns,2) = ew
          end if
               
        end if
      end do
    end do
  
    totpts = covtot
                                             
    dom(1:2) = (/ewn,1/); empty = .true.

    do ns = 1,nsn
           
      if (full(ns)) then

        if (empty) then
          empty  = .false.
          dom(3) = ns
        end if
        dom(4) = ns
        dom(1) = min0(dom(1),band(ns,1))
        dom(2) = max0(dom(2),band(ns,2))
      end if
    end do

  contains

    logical function thckcrit(ca,cb)

      implicit none

      real(dp),dimension(:,:),intent(in) :: ca 
      real(sp),               intent(in) :: cb

      select case (whichthck)
      case(5)

        ! whichthck=5 is not a 'known case'

        if ( ca(2,2) > 0.0d0 .or. cb > 0.0) then
          thckcrit = .true.
        else
          thckcrit = .false.
        end if

      case default

        ! If the thickness in the region under consideration
        ! or the mass balance is positive, thckcrit is .true.

        if ( any((ca(:,:) > 0.0d0)) .or. cb > 0.0 ) then
          thckcrit = .true.
        else
          thckcrit = .false.
        end if

      end select

    end function thckcrit

  end subroutine glide_maskthck

  subroutine glide_calclsrf(thck,topg,lsrf)

    !*FD Calculates the elevation of the lower surface of the ice, 
    !*FD by considering whether it is floating or not.

    use glimmer_global, only : dp
    use physcon, only : rhoi, rhoo

    implicit none

    real(dp), intent(in),  dimension(:,:) :: thck !*FD Ice thickness
    real(dp), intent(in),  dimension(:,:) :: topg !*FD Bedrock topography elevation
    real(dp), intent(out), dimension(:,:) :: lsrf !*FD Lower ice surface elevation

    real(dp), parameter :: con = - rhoi / rhoo

    where (topg < con * thck)
      lsrf = con * thck
    elsewhere
      lsrf = topg
    end where
  end subroutine glide_calclsrf

!-------------------------------------------------------------------------

  subroutine glide_marinlim(which,whicht,thck,usrf,relx,topg,lati,mlimit,eus)

    !*FD Removes non-grounded ice, according to one of two altenative
    !*FD criteria, and sets upper surface of non-ice-covered points 
    !*FD equal to the topographic height, or sea-level, whichever is higher.

    use glimmer_global, only : dp, sp
    use paramets, only : f  

    implicit none

    !---------------------------------------------------------------------
    ! Subroutine arguments
    !---------------------------------------------------------------------

    integer,                intent(in)    :: which   !*FD Option to choose ice-removal method
                                                     !*FD \begin{description}
                                                     !*FD \item[0] Set thickness to zero if 
                                                     !*FD relaxed bedrock is below a given level.
                                                     !*FD \item[1] Set thickness to zero if
                                                     !*FD ice is floating.
                                                     !*FD \end{description}
    integer,                intent(in)    :: whicht  !*FD Thickness calculation option. Only acted on
                                                     !*FD if equals six.
    real(dp),dimension(:,:),intent(inout) :: thck    !*FD Ice thickness (scaled)
    real(dp),dimension(:,:),intent(out)   :: usrf    !*FD Upper ice surface (scaled)
    real(dp),dimension(:,:),intent(in)    :: relx    !*FD Relaxed topography (scaled)
    real(dp),dimension(:,:),intent(in)    :: topg    !*FD Actual topography (scaled)
    real(sp),dimension(:,:),intent(in)    :: lati    !*FD Array of latitudes (only used if 
                                                     !*FD $\mathtt{whicht}=6$).
    real(dp)                              :: mlimit  !*FD Lower limit on topography elevation for
                                                     !*FD ice to be present (scaled). Used with 
                                                     !*FD $\mathtt{which}=0$.
    real, intent(in) :: eus                          !*FD eustatic sea level

    !---------------------------------------------------------------------

    select case (which)

    case(0) ! Set thickness to zero if relaxed bedrock is below a 
            ! given level

      where (relx < mlimit+eus)
        thck = 0.0d0
        usrf = max(0.0d0,topg)
      end where

    case(1) ! Set thickness to zero if ice is floating

      where (thck < f * topg)
        thck = 0.0d0
        usrf = max(0.0d0,topg)
      end where
 
    end select

    !---------------------------------------------------------------------
    ! Not sure what this option does - it's only used when whichthck=6
    !---------------------------------------------------------------------

    select case(whicht)
    case(6)
      where (lati == 0.0)
        thck = 0.0d0
        usrf = max(0.0d0,topg)
      end where
    end select
  end subroutine glide_marinlim

  subroutine glide_load_sigma(model,unit)

    !*FD Loads a file containing
    !*FD sigma vertical coordinates.
    use glide_types
    use glimmer_log
    implicit none

    ! Arguments
    type(glide_global_type),intent(inout) :: model !*FD Ice model to use
    integer,               intent(in)    :: unit  !*FD Logical file unit to use. 
                                                  !*FD The logical file unit specified 
                                                  !*FD must not already be in use

    ! Internal variables

    integer :: up,upn
    logical :: there

    ! Beginning of code

    upn=model%general%upn

    inquire (exist=there,file=model%funits%sigfile)
  
    if (there) then
       call write_log('Reading sigma file: '//model%funits%sigfile)
       open(unit,file=model%funits%sigfile)
       read(unit,'(f5.2)',err=10,end=10) (model%numerics%sigma(up), up=1,upn)
       close(unit)
       return
    else
       call write_log('Calculating sigma')
       do up=1,upn
          model%numerics%sigma(up) = f(real(up-1)/real(upn-1),2.)
       end do
       return
    end if
    
10  call error_log('something wrong with sigma coord file')
    
    stop
    
  contains
    function f(x,n)
      implicit none
      real :: f,x,n
      
      f = (1-(x+1)**(-n))/(1-2**(-n))
    end function f
  end subroutine glide_load_sigma


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

  subroutine print_grid(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message


    call write_log('Grid specification')
    call write_log('------------------')
    write(message,*) 'ewn             : ',model%general%ewn
    call write_log(trim(message))
    write(message,*) 'nsn             : ',model%general%nsn
    call write_log(trim(message))
    write(message,*) 'upn             : ',model%general%upn
    call write_log(trim(message))
    write(message,*) 'EW grid spacing : ',model%numerics%dew
    call write_log(trim(message))
    write(message,*) 'NS grid spacing : ',model%numerics%dns
    call write_log(trim(message))
    write(message,*) 'sigma file      : ',trim(model%funits%sigfile)
    call write_log(trim(message))
    call write_log('')
  end subroutine print_grid

  ! time
  subroutine handle_time(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'tstart',model%numerics%tstart)
    call GetValue(section,'tend',model%numerics%tend)
    call GetValue(section,'dt',model%numerics%tinc)
    call GetValue(section,'ntem',model%numerics%ntem)
    call GetValue(section,'nvel',model%numerics%nvel)
    call GetValue(section,'niso',model%numerics%niso)
  end subroutine handle_time
  
  subroutine print_time(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message

    call write_log('Time steps')
    call write_log('----------')
    write(message,*) 'start time        : ',model%numerics%tstart
    call write_log(message)
    write(message,*) 'end time          : ',model%numerics%tend
    call write_log(message)
    write(message,*) 'main time step    : ',model%numerics%tinc
    call write_log(message)
    write(message,*) 'thermal dt factor : ',model%numerics%ntem
    call write_log(message)
    write(message,*) 'velo dt factor    : ',model%numerics%nvel
    call write_log(message)
    write(message,*) 'isostasy dt factor: ',model%numerics%niso
    call write_log('')
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
    call GetValue(section,'topo_is_relaxed',model%options%whichrelaxed)
    call GetValue(section,'hotstart',model%options%hotstart)
  end subroutine handle_options
  
  subroutine print_options(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message

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
    character(len=*), dimension(0:3), parameter :: slip_coeff = (/ &
         'zero        ', &
         'const       ', &
         'const if T>0', &
         '~basal water' /)
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

    call write_log('GLIDE options')
    call write_log('-------------')
    write(message,*) 'I/O parameter file      : ',trim(model%funits%ncfile)
    call write_log(message)
    if (model%options%whichtemp.lt.0 .or. model%options%whichtemp.ge.size(temperature)) then
       call error_log('Error, temperature out of range')
       stop
    end if
    write(message,*) 'temperature calculation : ',model%options%whichtemp,temperature(model%options%whichtemp)
    call write_log(message)
    if (model%options%whichflwa.lt.0 .or. model%options%whichflwa.ge.size(flow_law)) then
       call error_log('Error, flow_law out of range')
       stop
    end if
    write(message,*) 'flow law                : ',model%options%whichflwa,flow_law(model%options%whichflwa)
    call write_log(message)
    if (model%options%whichisot.lt.0 .or. model%options%whichisot.ge.size(isostasy)) then
       call error_log('Error, isostasy out of range')
       stop
    end if
    write(message,*) 'isostasy                : ',model%options%whichisot,isostasy(model%options%whichisot)
    call write_log(message)
    if (model%options%whichslip.lt.0 .or. model%options%whichslip.ge.size(sliding)) then
       call error_log('Error, sliding_law out of range')
       stop
    end if
    write(message,*) 'sliding_law             : ',model%options%whichslip, sliding(model%options%whichslip)
    call write_log(message)
    if (model%options%whichbwat.lt.0 .or. model%options%whichbwat.ge.size(basal_water)) then
       call error_log('Error, basal_water out of range')
       stop
    end if
    write(message,*) 'basal_water             : ',model%options%whichbwat,basal_water(model%options%whichbwat)
    call write_log(message)
    if (model%options%whichmarn.lt.0 .or. model%options%whichmarn.ge.size(marine_margin)) then
       call error_log('Error, marine_margin out of range')
       stop
    end if
    write(message,*) 'marine_margin           : ', model%options%whichmarn, marine_margin(model%options%whichmarn)
    call write_log(message)
    if (model%options%whichbtrc.lt.0 .or. model%options%whichbtrc.ge.size(slip_coeff)) then
       write(*,*) 'Error, slip_coeff out of range'
       stop
    end if
    write(message,*) 'slip_coeff              : ', model%options%whichbtrc, slip_coeff(model%options%whichbtrc)
    call write_log(message)
    if (model%options%whichstrs.lt.0 .or. model%options%whichstrs.ge.size(stress)) then
       call error_log('Error, stress_calc out of range')
       stop
    end if
    write(message,*) 'stress_calc             : ', model%options%whichstrs, stress(model%options%whichstrs)
    call write_log(message)
    if (model%options%whichevol.lt.0 .or. model%options%whichevol.ge.size(evolution)) then
       call error_log('Error, evolution out of range')
       stop
    end if
    write(message,*) 'evolution               : ', model%options%whichevol, evolution(model%options%whichevol)
    call write_log(message)
    if (model%options%whichwvel.lt.0 .or. model%options%whichwvel.ge.size(vertical_integration)) then
       call error_log('Error, vertical_integration out of range')
       stop
    end if
    write(message,*) 'vertical_integration    : ',model%options%whichwvel,vertical_integration(model%options%whichwvel)
    call write_log(message)
    if (model%options%whichrelaxed.eq.1) then
       call write_log('First topo time slice is relaxed')
    end if
    if (model%options%hotstart.eq.1) then
       call write_log('Hotstarting model')
    end if
    call write_log('')
  end subroutine print_options

  ! parameters
  subroutine handle_parameters(section, model)
    use glimmer_config
    use glide_types
    use glimmer_log
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
       model%paramets%btrac_const=temp(1)
       deallocate(temp)
    end if
    call GetValue(section,'basal_tract_const',model%paramets%btrac_const)
  end subroutine handle_parameters

  subroutine print_parameters(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message

    call write_log('Parameters')
    call write_log('----------')
    write(message,*) 'ice limit             : ',model%numerics%thklim
    call write_log(message)
    write(message,*) 'marine depth limit    : ',model%numerics%mlimit
    call write_log(message)
    write(message,*) 'geothermal heat flux  : ',model%paramets%geot
    call write_log(message)
    write(message,*) 'flow enhancement      : ',model%paramets%fiddle
    call write_log(message)
    write(message,*) 'basal hydro time const: ',model%paramets%hydtim
    call write_log(message)
    write(message,*) 'isostasy time const   : ',model%paramets%isotim
    call write_log(message)
    write(message,*) 'basal traction param  : ',model%paramets%btrac_const
    call write_log(message)
    write(message,*) 'basal traction factors: ',model%paramets%bpar(1)
    call write_log(message)
    write(message,*) '                        ',model%paramets%bpar(2)
    call write_log(message)
    write(message,*) '                        ',model%paramets%bpar(3)
    call write_log(message)
    write(message,*) '                        ',model%paramets%bpar(4)
    call write_log(message)
    write(message,*) '                        ',model%paramets%bpar(5)
    call write_log(message)
    call write_log('')
  end subroutine print_parameters


end module glide_setup
