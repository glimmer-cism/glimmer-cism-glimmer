
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_setup.f90 - part of the GLIMMER ice model        + 
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

module glimmer_setup

  !*FD Contains general routines for initialisation, etc, called
  !*FD from the top-level glimmer subroutines.

contains

  subroutine read_config_file(unit,filename,model,proj)

    use physcon,  only: scyr
    use paramets, only: acc0,thk0,tim0,len0
    use glimmer_ncparams
    use glide_messages
    use glimmer_types
    use glimmer_project
    use glimmer_config
    implicit none

    !*FD Reads in configuration file for an individual ice model instance. 

    integer,                  intent(in)    :: unit     !*FD Logical file unit to use for reading.
    character(*),             intent(in)    :: filename !*FD Filename to read.
    type(glimmer_global_type),intent(inout) :: model    !*FD Model parameters to set.
    type(projection),         intent(inout) :: proj     !*FD Projection parameters to set.

    ! Internal variables/pointers

    type(ConfigSection), pointer :: config,section
    character(80) :: outtxt
    real(sp),dimension(:),pointer :: airt_temp => null()
    real(sp),dimension(:),pointer :: nmsb_temp => null()
    real(dp),dimension(:),pointer :: bpar_temp => null()

    ! Lists of allowed section/value names
    ! REMEMBER to UPDATE THESE if new variables are added,
    ! or their names in the config flie are changed

    character(20),dimension(11) :: allowed_sections = (/ &
         'output',            &
         'domain size',       &
         'projection',        &
         'sigma coordinates', &
         'options',           &
         'timesteps',         &
         'grid-lengths',      &
         'parameters',        &
         'forcing',           &
         'constants',         &
         'PDD scheme'/)

    character(20),dimension(1)  :: output_values = (/'configuration file'/)
    character(20),dimension(3)  :: domain_values = (/'east-west','north-south','vertical'/)
    character(20),dimension(6)  :: projection_values = (/ &
         'type',             &
         'false easting',    &
         'false northing',   &
         'parallel',         &
         'central meridian', &
         'standard parallel' /)
    character(20),dimension(1) :: sigma_values = (/'file'/)
    character(20),dimension(13) :: options_values = (/ &
         'whichtemp', &
         'whichartm', &
         'whichthck', &
         'whichflwa', &
         'whichisot', &
         'whichslip', &
         'whichbwat', &
         'whichmarn', &
         'whichbtrc', &
         'whichacab', &
         'whichevol', &
         'whichwvel', &
         'whichprecip' /)
    character(20),dimension(3) :: timesteps_values = (/ &
         'temperature','velocity','isostasy'/)
    character(20),dimension(2) :: grid_values = (/ &
         'east-west','north-south' /)
    character(20),dimension(9) :: params_values = (/ &
         'geot','fiddle','airt','nmsb','hydtim','isotim', &
         'bpar','thklim','mlimit'/)
    character(20),dimension(3) :: forcing_values = (/ &
         'file','trun','pfac'/)
    character(20),dimension(4) :: constants_values = (/ &
         'albedo','lapse rate','precip rate','air temp'/)
    character(20),dimension(3) :: pdd_values = (/ &
         'wmax','ice pdd factor','snow pdd factor'/)

    ! Open and read the configuration structure

    call ConfigRead(filename,config)

    ! Verify we have allowed section names

    if (ValidateSections(config,allowed_sections)/=0) then
       write(outtxt,*)'Unexpected sections in configuration file: ',trim(filename)
       call glide_msg(GM_FATAL,__FILE__,__LINE__,trim(outtxt))
    end if
     
    ! read Output section

    call GetSection(config,section,allowed_sections(1))
    if (associated(section)) then
       if (ValidateValueNames(section,output_values)/=0) then
          write(outtxt,*)'Unexpected value names in configuration file: ',trim(filename)
          call glide_msg(GM_FATAL,__FILE__,__LINE__,trim(outtxt))
       end if
       call GetValue(section,'configuration file',model%funits%ncfile)
    end if

    ! read domain size section

    call GetSection(config,section,allowed_sections(2))
    if (associated(section)) then
       if (ValidateValueNames(section,domain_values)/=0) then
          write(outtxt,*)'Unexpected value names in configuration file: ',trim(filename)
          call glide_msg(GM_FATAL,__FILE__,__LINE__,trim(outtxt))
       end if
       call GetValue(section,'east-west',  model%general%ewn)
       call GetValue(section,'north-south',model%general%nsn)
       call GetValue(section,'vertical',   model%general%upn)
    end if

    ! read projection section

    call GetSection(config,section,allowed_sections(3))
    if (associated(section)) then
       if (ValidateValueNames(section,projection_values)/=0) then
          write(outtxt,*)'Unexpected value names in configuration file: ',trim(filename)
          call glide_msg(GM_FATAL,__FILE__,__LINE__,trim(outtxt))
       end if
       call GetValue(section,'type',             proj%p_type)
       call GetValue(section,'false easting',    proj%cpx)
       call GetValue(section,'false northing',   proj%cpy)
       call GetValue(section,'parallel',         proj%latc)
       call GetValue(section,'central meridian', proj%lonc)
       call GetValue(section,'standard parallel',proj%std_par)
     end if

    ! read sigma coordinates section

    call GetSection(config,section,allowed_sections(4))
    if (associated(section)) then
       if (ValidateValueNames(section,sigma_values)/=0) then
          write(outtxt,*)'Unexpected value names in configuration file: ',trim(filename)
          call glide_msg(GM_FATAL,__FILE__,__LINE__,trim(outtxt))
       end if
       call GetValue(section,'file',model%funits%sigfile)
    end if

    ! read options section

    call GetSection(config,section,allowed_sections(5))
    if (associated(section)) then
       if (ValidateValueNames(section,options_values)/=0) then
          write(outtxt,*)'Unexpected value names in configuration file: ',trim(filename)
          call glide_msg(GM_FATAL,__FILE__,__LINE__,trim(outtxt))
       end if
       call GetValue(section,'whichtemp',model%options%whichtemp)
       call GetValue(section,'whichartm',model%options%whichartm)
       call GetValue(section,'whichthck',model%options%whichthck)
       call GetValue(section,'whichflwa',model%options%whichflwa)
       call GetValue(section,'whichisot',model%options%whichisot)
       call GetValue(section,'whichslip',model%options%whichslip)
       call GetValue(section,'whichbwat',model%options%whichbwat)
       call GetValue(section,'whichmarn',model%options%whichmarn)
       call GetValue(section,'whichbtrc',model%options%whichbtrc)
       call GetValue(section,'whichacab',model%options%whichacab)
       call GetValue(section,'whichevol',model%options%whichevol)
       call GetValue(section,'whichwvel',model%options%whichwvel)
       call GetValue(section,'whichprecip',model%options%whichprecip)
    end if

    ! read timesteps section

    call GetSection(config,section,allowed_sections(6))
    if (associated(section)) then
       if (ValidateValueNames(section,timesteps_values)/=0) then
          write(outtxt,*)'Unexpected value names in configuration file: ',trim(filename)
          call glide_msg(GM_FATAL,__FILE__,__LINE__,trim(outtxt))
       end if
       call GetValue(section,'temperature',model%numerics%ntem)
       call GetValue(section,'velocity',   model%numerics%nvel)
       call GetValue(section,'isostasy',   model%numerics%niso)
    end if

    ! read grid-lengths section

    call GetSection(config,section,allowed_sections(7))
    if (associated(section)) then
       if (ValidateValueNames(section,grid_values)/=0) then
          write(outtxt,*)'Unexpected value names in configuration file: ',trim(filename)
          call glide_msg(GM_FATAL,__FILE__,__LINE__,trim(outtxt))
       end if
       call GetValue(section,'east-west',  model%numerics%dew)
       call GetValue(section,'north-south',model%numerics%dns)
    end if

    ! read parameters section

    call GetSection(config,section,allowed_sections(8))
    if (associated(section)) then
       if (ValidateValueNames(section,params_values)/=0) then
          write(outtxt,*)'Unexpected value names in configuration file: ',trim(filename)
          call glide_msg(GM_FATAL,__FILE__,__LINE__,trim(outtxt))
       end if
       call GetValue(section,'geot',  model%paramets%geot)
       call GetValue(section,'fiddle',model%paramets%fiddle)
       call GetValue(section,'airt',  airt_temp,numval=size(model%paramets%airt))
       call GetValue(section,'nmsb',  nmsb_temp,numval=size(model%paramets%nmsb))
       call GetValue(section,'hydtim',model%paramets%hydtim)
       call GetValue(section,'isotim',model%paramets%isotim)
       call GetValue(section,'bpar',  bpar_temp,numval=size(model%paramets%bpar))
       call GetValue(section,'thklim',model%numerics%thklim)
       call GetValue(section,'mlimit',model%numerics%mlimit)

       if (associated(airt_temp)) then
          if (size(airt_temp).ne.size(model%paramets%airt)) then
             call glide_msg(GM_FATAL,__FILE__,__LINE__,'Wrong number of elements in airt')
          else
             model%paramets%airt=airt_temp
          end if
       endif

       if (associated(nmsb_temp)) then
          if (size(nmsb_temp).ne.size(model%paramets%nmsb)) then
             call glide_msg(GM_FATAL,__FILE__,__LINE__,'Wrong number of elements in nmsb')
          else
             model%paramets%nmsb=nmsb_temp
          end if
       end if

       if (associated(bpar_temp)) then
          if (size(bpar_temp).ne.size(model%paramets%bpar)) then
             call glide_msg(GM_FATAL,__FILE__,__LINE__,'Wrong number of elements in bpar')
          else
             model%paramets%bpar=bpar_temp
          end if
       end if

    end if

    ! read forcing section

    call GetSection(config,section,allowed_sections(9))
    if (associated(section)) then
       if (ValidateValueNames(section,forcing_values)/=0) then
          write(outtxt,*)'Unexpected value names in configuration file: ',trim(filename)
          call glide_msg(GM_FATAL,__FILE__,__LINE__,trim(outtxt))
       end if
       call GetValue(section,'file',model%funits%forcfile)
       call GetValue(section,'trun',model%forcdata%trun)
       call GetValue(section,'pfac',model%climate%pfac)
    end if

    ! read constants section

    call GetSection(config,section,allowed_sections(10))
    if (associated(section)) then
       if (ValidateValueNames(section,constants_values)/=0) then
          write(outtxt,*)'Unexpected value names in configuration file: ',trim(filename)
          call glide_msg(GM_FATAL,__FILE__,__LINE__,trim(outtxt))
       end if
       call GetValue(section,'albedo',     model%climate%ice_albedo)
       call GetValue(section,'lapse rate', model%climate%ulapse_rate)
       call GetValue(section,'precip rate',model%climate%uprecip_rate)
       call GetValue(section,'air temp',   model%climate%usurftemp)
    end if

    ! read PDD scheme section

    call GetSection(config,section,allowed_sections(11))
    if (associated(section)) then
       if (ValidateValueNames(section,pdd_values)/=0) then
          write(outtxt,*)'Unexpected value names in configuration file: ',trim(filename)
          call glide_msg(GM_FATAL,__FILE__,__LINE__,trim(outtxt))
       end if
       call GetValue(section,'wmax',model%pddcalc%wmax)
       call GetValue(section,'ice pdd factor',model%pddcalc%pddfac_ice)
       call GetValue(section,'snow pdd factor',model%pddcalc%pddfac_snow)
    end if

    ! Copy some parameters and derive others

    proj%nx=model%general%ewn
    proj%ny=model%general%nsn
    proj%dx=model%numerics%dew
    proj%dy=model%numerics%dns

    model%numerics%ntem = model%numerics%ntem * model%numerics%tinc
    model%numerics%nvel = model%numerics%nvel * model%numerics%tinc
    model%numerics%niso = model%numerics%niso * model%numerics%tinc

    model%numerics%dt     = model%numerics%tinc * scyr / tim0
    model%numerics%dttem  = model%numerics%ntem * scyr / tim0 
    model%numerics%thklim = model%numerics%thklim  / thk0

    model%numerics%dew = model%numerics%dew / len0
    model%numerics%dns = model%numerics%dns / len0

    model%paramets%nmsb(1) = model%paramets%nmsb(1) / (acc0 * scyr) 
    model%paramets%nmsb(2) = model%paramets%nmsb(2) / (acc0 * scyr)

    model%paramets%isotim = model%paramets%isotim * scyr / tim0         
    model%numerics%mlimit = model%numerics%mlimit / thk0

    ! Read output file configuration
    call ReadNCParams(model)

  end subroutine read_config_file

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine testinisthk(model,unit,first,newtemps,global_orog)

    !*FD Loads in relevant starting fields, and calculates the 
    !*FD present-day temperature field for the purposes of
    !*FD parameterised precip and temperature.

    ! Use statements for modules containing parameters 

    use paramets
    use physcon, only : scyr

    ! Use statements for modules containing subprograms

    use glimmer_outp
    use glimmer_temp
    use glimmer_thck
	use glimmer_velo
    use glide_messages
    use glimmer_ncinfile
    implicit none

    ! Subroutine arguments

    type(glimmer_global_type),intent(inout) :: model !*FD Model parameters being set
    integer,                  intent(in)    :: unit  !*FD File unit to use for read operations
    logical,                  intent(inout) :: first !*FD The `first time' flag for the instance
	logical,                  intent(inout) :: newtemps !*FD The new temperatures flag
	real(rk),dimension(:,:),  intent(in)    :: global_orog !*FD The global orography on the local grid

    ! Internal variables

    real(sp),dimension(:,:),allocatable :: arng
    type(glimmer_nc_input), pointer :: ic
    logical found_precip,found_presurf,found_usurf

	!-----------------------------------------------------
    ! open all netCDF input files
	!-----------------------------------------------------

    call openall_in(model)

    ic=>model%funits%in_first
    found_precip = .false.
    found_presurf = .false.
    found_usurf = .false.
    do while(associated(ic))
       ! read present-day precip file if required      
       if (model%options%whichprecip==3) then
          if (ic%nc%do_var(NC_B_PRESPRCP)) then
             found_precip = .true.
          end if
       else
          ic%nc%do_var(NC_B_PRESPRCP) = .false.
       end if

       if ((model%options%whichartm.eq.0).or. (model%options%whichartm.eq.1).or. &
            (model%options%whichartm.eq.2).or. (model%options%whichartm.eq.3).or. &
            (model%options%whichartm.eq.4)) then
          if (ic%nc%do_var(NC_B_PRESUSRF)) then
             found_presurf = .true.
          else
             ic%nc%do_var(NC_B_PRESUSRF) = .false.
          end if
       end if

       if (ic%nc%do_var(NC_B_USURF)) then
          found_usurf = .true.
       end if

       ic=>ic%next
    end do

    if (model%options%whichprecip==3 .and. .not. found_precip) then
       call glide_msg(GM_FATAL,__FILE__,__LINE__,"To run with whichprecip=3, you need to supply a forcing filename")
    endif

    call readall(model)

    ! -------------------------------------------------------------------
    ! Read in forcing file, if required, or set up simple forcing
    ! otherwise
    ! -------------------------------------------------------------------

    if (model%funits%forcfile .ne. 'none') then
      call redtsout(model%funits%forcfile,unit,model%forcdata)
    else
      allocate(model%forcdata%forcing(1,2))
      model%forcdata%forcing(1,:) = (/ 2 * model%forcdata%trun, 0.0 /)
      model%forcdata%flines = 1
    end if

    ! -------------------------------------------------------------------
    ! Read present-day surface elevation, if required, and calculate
    ! the present-day surface temperature from it.
    ! -------------------------------------------------------------------

    if (found_presurf) then

       ! Allocate arng array, passed to calcartm for air temperature range
       
       allocate(arng(size(model%climate%presartm,1),size(model%climate%presartm,2)))

       !----------------------------------------------------------------------
       ! Calculate the present-day mean air temperature and range, based on
       ! surface elevation and latitude
       !----------------------------------------------------------------------
       
       call calcartm(model, 3, &
            model%climate%presusrf, &
            model%climate%lati,     &
            model%climate%presartm, &  !** OUTPUT
            arng)                      !** OUTPUT
       
       !----------------------------------------------------------------------
       ! Calculate present-day mass-balance based on present-day elevation,
       ! temperature, temperature range, and PDD method.
       !----------------------------------------------------------------------
       
       call calcacab(model%numerics,         &
            model%paramets,         &
            model%pddcalc,          &
            1,                      &
            model%geometry%usrf,  &
            model%climate%presartm, &
            arng,                   &
            model%climate%presprcp, &
            model%climate%ablt,     &
            model%climate%lati,     &
            model%climate%acab,     &
            model%geometry%thck)     
       
       ! Set ice thickness to be mass-balance*time-step, where positive
       
       model%geometry%thck = max(0.0d0,model%climate%acab*model%numerics%dt)
       
       ! Calculate the elevation of the lower ice surface
       
       call calclsrf(model%geometry%thck,model%geometry%topg,model%geometry%lsrf)
       
       ! Calculate the elevation of the upper ice surface by adding thickness
       ! onto the lower surface elevation.
       
       model%geometry%usrf = model%geometry%thck + model%geometry%lsrf
       
       call calcartm(model, 3, &
            model%geometry%usrf, &
            model%climate%lati,     &
            model%climate%artm, &  !** OUTPUT
            arng)                      !** OUTPUT
       call timeevoltemp(model,0,model%climate%artm)     ! calculate initial temperature distribution
	   newtemps = .true.                          ! we have new temperatures


       first=.false.
       deallocate(arng) 
       
    else    
       
       ! -----------------------------------------------------------------
       ! Calculate the lower and upper surfaces of the ice-sheet 
       ! -----------------------------------------------------------------
       
       call calclsrf(model%geometry%thck,model%geometry%topg,model%geometry%lsrf)
       model%geometry%usrf = model%geometry%thck + model%geometry%lsrf
       
    endif

    ! -------------------------------------------------------------------
    ! Calculate the upper surface, if necessary, otherwise calculates
    ! the thickness
    ! -------------------------------------------------------------------
    
    if (found_usurf) then
       model%geometry%thck = model%geometry%usrf - model%geometry%lsrf
    endif
    
  end subroutine testinisthk

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine maskthck(whichthck,crita,critb,dom,pointno,totpts,empty)
    
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

  end subroutine maskthck

!-------------------------------------------------------------------------

  subroutine calclsrf(thck,topg,lsrf)

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

  end subroutine calclsrf

!-------------------------------------------------------------------------

  subroutine marinlim(which,whicht,thck,usrf,relx,topg,lati,mlimit,ablation_field)

    !*FD Removes non-grounded ice, according to one of two altenative
    !*FD criteria, and sets upper surface of non-ice-covered points 
    !*FD equal to the topographic height, or sea-level, whichever is higher.

    use glimmer_global, only : dp, sp,rk
    use physcon, only : f  

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
    real(rk),dimension(:,:),intent(inout) :: ablation_field 

    !---------------------------------------------------------------------

    select case (which)

    case(0) ! Set thickness to zero if relaxed bedrock is below a 
            ! given level

      where (relx < mlimit)
        ablation_field=ablation_field+thck
        thck = 0.0d0
        usrf = max(0.0d0,topg)
      end where

    case(1) ! Set thickness to zero if ice is floating

      where (thck < f * topg)
        ablation_field=ablation_field+thck
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

  end subroutine marinlim

!-------------------------------------------------------------------------

  subroutine allocarr(model)

    !*FD Allocates the model arrays, and initialises some of them to zero.
    !*FD These are the arrays allocated, and their dimensions:
    !*FD
    !*FD In \texttt{model\%temper}:
    !*FD \begin{itemize}
    !*FD \item \texttt{temp(upn,0:ewn+1,0:nsn+1))}
    !*FD \item \texttt{flwa(upn,ewn,nsn))}
    !*FD \item \texttt{bwat(ewn,nsn))}
    !*FD \item \texttt{bmlt(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%climate}:
    !*FD \begin{itemize}
    !*FD \item \texttt{acab(ewn,nsn))}
    !*FD \item \texttt{artm(ewn,nsn))}
    !*FD \item \texttt{arng(ewn,nsn))}
    !*FD \item \texttt{g\_artm(ewn,nsn))}
    !*FD \item \texttt{g\_arng(ewn,nsn))}
    !*FD \item \texttt{lati(ewn,nsn))}
    !*FD \item \texttt{ablt(ewn,nsn))}
    !*FD \item \texttt{prcp(ewn,nsn))}
    !*FD \item \texttt{presprcp(ewn,nsn))}
    !*FD \item \texttt{presartm(ewn,nsn))}
    !*FD \item \texttt{presusrf(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%velocity}:
    !*FD \begin{itemize}
    !*FD \item \texttt{uvel(upn,ewn-1,nsn-1))}
    !*FD \item \texttt{vvel(upn,ewn-1,nsn-1))}
    !*FD \item \texttt{wvel(upn,ewn,nsn))}
    !*FD \item \texttt{wgrd(upn,ewn,nsn))}
    !*FD \item \texttt{uflx(ewn-1,nsn-1))}
    !*FD \item \texttt{vflx(ewn-1,nsn-1))}
    !*FD \item \texttt{diffu(ewn,nsn))}
    !*FD \item \texttt{btrc(ewn,nsn))}
    !*FD \item \texttt{ubas(ewn,nsn))}
    !*FD \item \texttt{vbas(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%geomderv}:
    !*FD \begin{itemize}
    !*FD \item \texttt{dthckdew(ewn,nsn))}
    !*FD \item \texttt{dusrfdew(ewn,nsn))}
    !*FD \item \texttt{dthckdns(ewn,nsn))}
    !*FD \item \texttt{dusrfdns(ewn,nsn))}
    !*FD \item \texttt{dthckdtm(ewn,nsn))}
    !*FD \item \texttt{dusrfdtm(ewn,nsn))}
    !*FD \item \texttt{stagthck(ewn-1,nsn-1))}
    !*FD \end{itemize}
  
    !*FD In \texttt{model\%geometry}:
    !*FD \begin{itemize}
    !*FD \item \texttt{thck(ewn,nsn))}
    !*FD \item \texttt{usrf(ewn,nsn))}
    !*FD \item \texttt{lsrf(ewn,nsn))}
    !*FD \item \texttt{topg(ewn,nsn))}
    !*FD \item \texttt{relx(ewn,nsn))}
    !*FD \item \texttt{mask(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%thckwk}:
    !*FD \begin{itemize}
    !*FD \item \texttt{olds(ewn,nsn,thckwk\%nwhich))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%isotwk}:
    !*FD \begin{itemize}
    !*FD \item \texttt{load(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%numerics}:
    !*FD \begin{itemize}
    !*FD \item \texttt{sigma(upn))}
    !*FD \end{itemize}

    use glimmer_types

    implicit none

    type(glimmer_global_type),intent(inout) :: model

    integer :: ewn,nsn,upn

    ! for simplicity, copy these values...

    ewn=model%general%ewn
    nsn=model%general%nsn
    upn=model%general%upn

    ! Allocate appropriately

    allocate(model%temper%temp(upn,0:ewn+1,0:nsn+1)); model%temper%temp = 0.0
    allocate(model%temper%flwa(upn,ewn,nsn))   
    allocate(model%temper%bwat(ewn,nsn));             model%temper%bwat = 0.0
    allocate(model%temper%bmlt(ewn,nsn));             model%temper%bmlt = 0.0

    allocate(model%climate%acab(ewn,nsn));            model%climate%acab = 0.0
    allocate(model%climate%artm(ewn,nsn));            model%climate%artm = 0.0
    allocate(model%climate%arng(ewn,nsn));            model%climate%arng = 0.0
    allocate(model%climate%g_artm(ewn,nsn));          model%climate%g_artm = 0.0
    allocate(model%climate%g_arng(ewn,nsn));          model%climate%g_arng = 0.0
    allocate(model%climate%lati(ewn,nsn));            model%climate%lati = 0.0
    allocate(model%climate%loni(ewn,nsn));            model%climate%loni = 0.0
    allocate(model%climate%out_mask(ewn,nsn));        model%climate%out_mask = 1.0
    allocate(model%climate%ablt(ewn,nsn));            model%climate%ablt = 0.0
    allocate(model%climate%prcp(ewn,nsn));            model%climate%prcp = 0.0
    allocate(model%climate%presprcp(ewn,nsn));        model%climate%presprcp = 0.0
    allocate(model%climate%presartm(ewn,nsn));        model%climate%presartm = 0.0
    allocate(model%climate%presusrf(ewn,nsn));        model%climate%presusrf = 0.0

    allocate(model%velocity%uvel(upn,ewn-1,nsn-1));   model%velocity%uvel = 0.0d0
    allocate(model%velocity%vvel(upn,ewn-1,nsn-1));   model%velocity%vvel = 0.0d0
    allocate(model%velocity%wvel(upn,ewn,nsn));       model%velocity%wvel = 0.0d0
    allocate(model%velocity%wgrd(upn,ewn,nsn));       model%velocity%wgrd = 0.0d0
    allocate(model%velocity%uflx(ewn-1,nsn-1));       model%velocity%uflx = 0.0d0
    allocate(model%velocity%vflx(ewn-1,nsn-1));       model%velocity%vflx = 0.0d0
    allocate(model%velocity%diffu(ewn,nsn));          model%velocity%diffu = 0.0d0
    allocate(model%velocity%btrc(ewn,nsn));           model%velocity%btrc = 0.0d0
    allocate(model%velocity%ubas(ewn,nsn));           model%velocity%ubas = 0.0d0
    allocate(model%velocity%vbas(ewn,nsn));           model%velocity%vbas = 0.0d0

    allocate(model%geomderv%dthckdew(ewn,nsn));       model%geomderv%dthckdew = 0.0d0 
    allocate(model%geomderv%dusrfdew(ewn,nsn));       model%geomderv%dusrfdew = 0.0d0
    allocate(model%geomderv%dthckdns(ewn,nsn));       model%geomderv%dthckdns = 0.0d0
    allocate(model%geomderv%dusrfdns(ewn,nsn));       model%geomderv%dusrfdns = 0.0d0
    allocate(model%geomderv%dthckdtm(ewn,nsn));       model%geomderv%dthckdtm = 0.0d0
    allocate(model%geomderv%dusrfdtm(ewn,nsn));       model%geomderv%dusrfdtm = 0.0d0
    allocate(model%geomderv%stagthck(ewn-1,nsn-1));   model%geomderv%stagthck = 0.0d0
  
    allocate(model%geometry%thck(ewn,nsn));           model%geometry%thck = 0.0d0
    allocate(model%geometry%usrf(ewn,nsn));           model%geometry%usrf = 0.0d0
    allocate(model%geometry%lsrf(ewn,nsn));           model%geometry%lsrf = 0.0d0
    allocate(model%geometry%topg(ewn,nsn));           model%geometry%topg = 0.0d0
    allocate(model%geometry%relx(ewn,nsn));           model%geometry%relx = 0.0d0
    allocate(model%geometry%mask(ewn,nsn));           model%geometry%mask = 0
    allocate(model%geometry%std_dev(ewn,nsn));        model%geometry%std_dev = 0

    allocate(model%thckwk%olds(ewn,nsn,model%thckwk%nwhich))
                                                      model%thckwk%olds = 0.0d0
    allocate(model%isotwk%load(ewn,nsn));             model%isotwk%load = 0.0d0 
    allocate(model%numerics%sigma(upn))

  end subroutine allocarr

end module glimmer_setup
