
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

  subroutine initial(unit,nmlfile,model,proj)

    !*FD Reads in namelists for an individual ice model instance. This
    !*FD subroutine also holds the defaults for all the parameters
    !*FD that are read in from the namelist

    use physcon,  only: scyr
    use paramets, only: acc0,thk0,tim0,len0
    use glimmer_types
    use glimmer_project
    use glimmer_ncparams
    implicit none

    ! Subroutine arguments

    integer,                  intent(in)    :: unit     !*FD Logical file unit to use for reading.
    character(*),             intent(in)    :: nmlfile  !*FD Namelist filename to read.
    type(glimmer_global_type),intent(inout) :: model    !*FD Model parameters to set.
    type(projection),         intent(inout) :: proj     !*FD Projection parameters to set.

    ! -------------------------------------------------------------------
    ! Internal variables - these are needed, as namelists can't include
    ! components of derived types, or pointers, or anything else useful...
    ! -------------------------------------------------------------------

    logical :: there  ! Flag used when checking file is present.

    ! temporary variables...
    ! ...for general

    integer  :: ewn,nsn,upn

    ! ...for the projection

    integer  :: p_type
    real(sp) :: cpx,cpy,latc,lonc,std_par

    ! ...for the filenames

    character(fname_length) :: sigfile,ncconf

    ! ...for options

    integer :: whichtemp, whichartm, whichthck, whichflwa, &
               whichisot, whichslip, whichbwat, whichmarn, &
               whichbtrc, whichacab, whichstrs, & 
               whichevol, whichwvel, whichprecip

    ! ...for numerics

    real(sp) :: ntem, nvel, niso
    real(sp),dimension(3) :: nout
    real(sp) :: nstr
    real(dp) :: thklim, mlimit, dew, dns

    ! ...for funits

    character(fname_length) :: forcfile

    ! ...for other parameters

    real(dp) :: geot,fiddle,hydtim,isotim
    real(sp),dimension(2) :: airt
    real(sp),dimension(3) :: nmsb
    real(dp),dimension(5) :: bpar

    ! ...for constants

    real(sp) :: lapse_rate,precip_rate,air_temp,albedo

    ! ...for the forcing

    real(sp) :: trun ! The total time for the run.

    ! -------------------------------------------------------------------
    ! namelist definitions
    ! -------------------------------------------------------------------

    namelist / sizs / ewn, nsn, upn
    namelist / prj  / p_type, cpx, cpy, latc, lonc, std_par
    namelist / opts / whichtemp, whichartm, whichthck, whichflwa, &
                      whichisot, whichslip, whichbwat, whichmarn, &
                      whichbtrc, whichacab, whichstrs, & 
                      whichevol, whichwvel, whichprecip
    namelist / nums / ntem, nvel, niso, nout, nstr, thklim, mlimit, dew, dns 
    namelist / pars / geot, fiddle, airt, nmsb, hydtim, isotim, bpar
    namelist / dats / forcfile
    namelist / cons / lapse_rate,precip_rate,air_temp,albedo
    namelist / forc / trun

    ! -------------------------------------------------------------------
    ! Copy default values of namelist parameters here
    ! -------------------------------------------------------------------

    ! For general:

    ewn = model%general%ewn
    nsn = model%general%nsn
    upn = model%general%upn

    ! For the projection:

    p_type = proj%p_type
    cpx    = proj%cpx
    cpy    = proj%cpy
    latc   = proj%latc
    lonc   = proj%lonc
    std_par= proj%std_par

    ! For options:

    whichtemp   = model%options%whichtemp
    whichartm   = model%options%whichartm
    whichthck   = model%options%whichthck
    whichflwa   = model%options%whichflwa
    whichisot   = model%options%whichisot
    whichslip   = model%options%whichslip
    whichbwat   = model%options%whichbwat
    whichmarn   = model%options%whichmarn
    whichbtrc   = model%options%whichbtrc
    whichacab   = model%options%whichacab
    whichstrs   = model%options%whichstrs
    whichevol   = model%options%whichevol
    whichwvel   = model%options%whichwvel
    whichprecip = model%options%whichprecip

    ! For numerics:

    ntem   = model%numerics%ntem
    nvel   = model%numerics%nvel
    niso   = model%numerics%niso
    nout   = model%numerics%nout
    nstr   = model%numerics%nstr
    thklim = model%numerics%thklim
    mlimit = model%numerics%mlimit
    dew    = model%numerics%dew
    dns    = model%numerics%dns

    ! For funits

    forcfile     = model%funits%forcfile

    ! For other parameters

    geot   = model%paramets%geot
    fiddle = model%paramets%fiddle
    hydtim = model%paramets%hydtim
    isotim = model%paramets%isotim
    airt   = model%paramets%airt
    nmsb   = model%paramets%nmsb
    bpar   = model%paramets%bpar

    ! For constants

    lapse_rate  = model%climate%ulapse_rate
    precip_rate = model%climate%uprecip_rate
    air_temp    = model%climate%usurftemp
    albedo      = model%climate%ice_albedo

    ! For the forcing

    trun        = model%forcdata%trun

    ! -------------------------------------------------------------------
    ! Load the namelist, having first checked that it exists
    ! -------------------------------------------------------------------

    inquire (exist=there,file=nmlfile)         ! Check the namelist file exists
  
    if (there) then                            ! If it does exist, read it in
      open(unit,file=nmlfile)
      read(unit,100)ncconf
      read(unit,nml=sizs)
      read(unit,nml=prj)
      read(unit,100)sigfile
      read(unit,nml=opts)
      read(unit,nml=nums)
      read(unit,nml=pars)
      read(unit,nml=dats)
      read(unit,nml=cons)
      read(unit,nml=forc)
      close(unit)
    else
      print*,'Error opening namelist file ',trim(nmlfile),' - it doesn''t exist!'
      stop
    end if

    ! -------------------------------------------------------------------
    ! Copy data from temporary namelist variables to the model parameters
    ! -------------------------------------------------------------------

    ! copy projection data

    proj%p_type=p_type
    proj%nx=ewn
    proj%ny=nsn
    proj%cpx=cpx
    proj%cpy=cpy
    proj%dx=dew
    proj%dy=dns
    proj%latc=latc
    proj%lonc=lonc
    proj%std_par=std_par

    ! copy to options type

    model%options%whichtemp=whichtemp
    model%options%whichartm=whichartm
    model%options%whichthck=whichthck
    model%options%whichflwa=whichflwa
    model%options%whichisot=whichisot
    model%options%whichslip=whichslip
    model%options%whichbwat=whichbwat
    model%options%whichmarn=whichmarn
    model%options%whichbtrc=whichbtrc
    model%options%whichacab=whichacab
    model%options%whichstrs=whichstrs
    model%options%whichevol=whichevol
    model%options%whichwvel=whichwvel
    model%options%whichprecip=whichprecip

    ! copy general type

    model%general%nsn = nsn
    model%general%ewn = ewn
    model%general%upn = upn

    ! copy numerics type

    model%numerics%ntem   = ntem
    model%numerics%nvel   = nvel
    model%numerics%niso   = niso
    model%numerics%nout   = nout
    model%numerics%nstr   = nstr
    model%numerics%thklim = thklim
    model%numerics%mlimit = mlimit
    model%numerics%dew    = dew
    model%numerics%dns    = dns

    ! copy funits type

    model%funits%forcfile = forcfile

    ! filenames

    sigfile=adjustl(sigfile)

    model%funits%sigfile     = sigfile
    model%funits%ncfile = adjustl(ncconf)

    ! Constants

    model%climate%uprecip_rate = precip_rate
    model%climate%usurftemp    = air_temp
    model%climate%ice_albedo   = albedo
    model%climate%ulapse_rate  = lapse_rate

    ! Forcing

    model%forcdata%trun = trun

    ! back to original code

    model%numerics%ntem = model%numerics%ntem * model%numerics%tinc
    model%numerics%nvel = model%numerics%nvel * model%numerics%tinc
    model%numerics%niso = model%numerics%niso * model%numerics%tinc

    model%numerics%dt     = model%numerics%tinc * scyr / tim0
    model%numerics%dttem  = model%numerics%ntem * scyr / tim0 
    model%numerics%thklim = model%numerics%thklim  / thk0

    model%numerics%dew = model%numerics%dew / len0
    model%numerics%dns = model%numerics%dns / len0

    ! deal with parameters here (partly original code)

    model%paramets%geot   = geot
    model%paramets%fiddle = fiddle
    model%paramets%airt   = airt

    model%paramets%nmsb(1) = nmsb(1) / (acc0 * scyr) 
    model%paramets%nmsb(2) = nmsb(2) / (acc0 * scyr)
  
    model%paramets%hydtim = hydtim
    model%paramets%isotim = isotim * scyr / tim0         

    model%paramets%bpar   = bpar

    model%numerics%mlimit = model%numerics%mlimit / thk0

    ! Read output file configuration
    call ReadNCParams(model)
    
100 format(A)

  end subroutine initial

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
       print*,"To run with whichprecip=3, you need to supply a"
       print*,"forcing filename..."
       stop
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
            model%climate%acab)     
       
       ! Set ice thickness to be mass-balance*time-step, where positive
       
       model%geometry%thck = max(0.0d0,model%climate%acab*model%numerics%dt)
       
       ! Calculate the elevation of the lower ice surface
       
       call calclsrf(model%geometry%thck,model%geometry%topg,model%geometry%lsrf)
       
       ! Calculate the elevation of the upper ice surface by adding thickness
       ! onto the lower surface elevation.
       
       model%geometry%usrf = model%geometry%thck + model%geometry%lsrf
       
       call timeevoltemp(model,0,global_orog)     ! calculate initial temperature distribution
	   newtemps = .true.                          ! we have new temperatures

	   call calcflwa(model%numerics,        &              ! Calculate Glen's A
                     model%velowk,          &
                     model%paramets%fiddle, &
                     model%temper%flwa,     &
                     model%temper%temp,     &
                     model%geometry%thck,   &
                     model%options%whichflwa)

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

  subroutine marinlim(which,whicht,thck,usrf,relx,topg,lati,mlimit)

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

    !---------------------------------------------------------------------

    select case (which)

    case(0) ! Set thickness to zero if relaxed bedrock is below a 
            ! given level

      where (relx < mlimit)
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

    allocate(model%thckwk%olds(ewn,nsn,model%thckwk%nwhich))
                                                      model%thckwk%olds = 0.0d0
    allocate(model%isotwk%load(ewn,nsn));             model%isotwk%load = 0.0d0 
    allocate(model%numerics%sigma(upn))

  end subroutine allocarr

end module glimmer_setup
