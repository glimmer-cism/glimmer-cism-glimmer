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

  subroutine calcartm(instance,which,usrf,lati,artm,arng,g_orog,g_artm,g_arng)

    !*FD Calculate the surface air temperature and mean annual range,
    !*FD according to various different models.

    use glimmer_global, only : dp, sp 
    use paramets, only : len0, thk0
    use glint_type
    use glint_mbal
    use glimmer_log
    implicit none


    type(glint_instance),intent(inout)   :: instance         !*FD Instance whose elements are to be allocated.
    real(dp),dimension(:,:),         intent(in)    :: usrf   !*FD Surface elevation (km)
    real(sp),dimension(:,:),         intent(in)    :: lati   !*FD Array of latitudes
    real(sp),dimension(:,:),         intent(out)   :: artm   !*FD Surface annual mean air temperature ($^{\circ}$C)
    real(sp),dimension(:,:),         intent(out)   :: arng   !*FD Surface annual air tempurature range ($^{\circ}$C)
    real(dp),dimension(:,:),optional,intent(in)    :: g_orog !*FD Global orography on local grid (m)
    real(sp),dimension(:,:),optional,intent(in)    :: g_artm !*FD Supplied global air temperatures ($^{\circ}$C)
    real(sp),dimension(:,:),optional,intent(in)    :: g_arng !*FD Supplied global air temp range ($^{\circ}$C)

    integer, intent(in) :: which                    !*FD which method to use (see documentation for
                                                    !*FD allowed values of whichartm)

    real(sp) :: dist, ewct, nsct, inve, esurf   
    integer :: ns,ew

    !--------------------------------------------------------------------
    ! If this is the first call, set up some constants
    !--------------------------------------------------------------------

    if (instance%model%temper%first1) then
      instance%model%temper%grid = instance%model%numerics%dew * len0

      if (which == 0 .or. which == 4) then
          instance%climate%airt(2) = instance%climate%airt(2) * thk0
      end if

       instance%model%temper%first1 = .false.
    end if

    !--------------------------------------------------------------------
    ! Now calculate the temp and range according to different methods:
    !--------------------------------------------------------------------

    select case(which)

            ! ----------------------------------------------------------
    case(0) ! Linear decrease from sea-level according to lapse-rate
            ! ----------------------------------------------------------

      artm = instance%climate%airt(1) + usrf * instance%climate%airt(2)

            ! ----------------------------------------------------------
    case(1) ! 2d EISMINT test - Cubic function of distance from domain centre
            ! ----------------------------------------------------------

      ewct = real(instance%model%general%ewn+1) / 2.0
      nsct = real(instance%model%general%nsn+1) / 2.0

      do ns = 1,instance%model%general%nsn
        do ew = 1,instance%model%general%ewn
          dist = instance%model%temper%grid * sqrt((real(ew) - ewct)**2 + (real(ns) - nsct)**2) 
          artm(ew,ns) = instance%climate%airt(1) + dist**3 * instance%climate%airt(2)
        end do
      end do            

      ! version for 1d (ew) eismint test
      ! ** ewct = real(model%general%tewn+1) / 2.0
      ! ** do ns = 1,model%general%nsn; do ew = 1,model%general%ewn
      ! **   dist = grid * sqrt((real(ew) - ewct)**2) 
      ! **   artm(ew,ns) = model%paramets%airt(1) + dist**3 * model%paramets%airt(2)
      ! ** end do; end do            
            ! ----------------------------------------------------------
    case(2) ! Linear function of distance from domain centre
            ! ----------------------------------------------------------

      ewct = real(instance%model%general%ewn+1) / 2.0; nsct = real(instance%model%general%nsn+1) / 2.0

      do ns = 1,instance%model%general%nsn
        do ew = 1,instance%model%general%ewn
          dist = instance%model%temper%grid * sqrt((real(ew) - ewct)**2 + (real(ns) - nsct)**2) 
          artm(ew,ns) = instance%climate%airt(1) + dist * instance%climate%airt(2)
        end do
      end do            

            ! ----------------------------------------------------------
    case(3) ! Obtaining appropriate temperature forcing from file
            ! ----------------------------------------------------------

      do while ( instance%model%numerics%time .ge. instance%forcdata%forcing(instance%model%temper%tpt+1,1) .and. &
            instance%model%temper%tpt+1 .le. instance%forcdata%flines)
          instance%model%temper%tpt = instance%model%temper%tpt + 1
          instance%model%temper%perturb = instance%forcdata%forcing(instance%model%temper%tpt,2)
         ! ** print '(a26,f10.1,f10.3)', &
         ! **       '---> forcing read at      ', model%numerics%time, perturb 
      end do

      !   find air temps based on regr equation for greenland

      do ns = 1, instance%model%general%nsn
        do ew = 1, instance%model%general%ewn 
      
      !   find latitude dependent height of inversion layer  
 
          inve = max(0.0, 300 * (lati(ew,ns) - 65.0) / 15.0)
       
      !   make sure that surface is above sealevel

          esurf = max(0.0, sngl(usrf(ew,ns)*thk0))
     
      !   find mean annual temp (is surface above inversion layer?)

          if ( esurf <= inve ) then
            artm(ew,ns) = 49.13 - 0.007992 * inve - 0.7576 * lati(ew,ns) + instance%model%temper%perturb
          else
            artm(ew,ns) = 49.13 - 0.007992 * esurf - 0.7576 * lati(ew,ns) + instance%model%temper%perturb
          end if

      !   find July temp (express as annual half range for convenience)
    
          arng(ew,ns)= 30.38 - 0.006277 * esurf - 0.3262 * lati(ew,ns) + instance%model%temper%perturb - artm(ew,ns)

        end do
      end do

            ! ----------------------------------------------------------
    case(4) ! Air temperature is function of latitude and height
            ! ----------------------------------------------------------
    
      ! Note that we are using the array lati to hold sealevel air temperatures

      artm = lati + usrf * instance%climate%airt(2)

            ! ----------------------------------------------------------
    case(5) ! Uniform temperature, zero range
            ! ----------------------------------------------------------

      artm=instance%climate%usurftemp
      arng=0.0

            ! ----------------------------------------------------------
    case(6) ! Uniform temperature, lapse-rate corrected, zero range
            ! ----------------------------------------------------------

      artm=instance%climate%usurftemp
      arng=0.0
      call glint_lapserate(artm,real(usrf*thk0,rk),real(instance%climate%ulapse_rate,rk))

            ! ----------------------------------------------------------
    case(7) ! Supplied large-scale temperature and range
            ! ----------------------------------------------------------

      ! Check we have the necessary arguments first

      if (present(g_arng).and.present(g_artm).and.present(g_orog)) then

      ! Copy the fields

        arng=g_arng
        artm=g_artm

      ! Reduce temperatures to sea-level

        call glint_lapserate(artm,real(g_orog,rk),real(-instance%climate%ulapse_rate,rk))

      ! Raise them to high-res orography 

        call glint_lapserate(artm,real(usrf*thk0,rk),real(instance%climate%ulapse_rate,rk))

      else
          call write_log('ERROR: Error in arguments to CALCARTM - stopping', &
              GM_FATAL,__FILE__,__LINE__)

      endif

            ! ----------------------------------------------------------
    case(8) ! Leave everything alone
            ! ----------------------------------------------------------

                 ! -----------------------------------------------------
    case default ! Flag an error otherwise
                 ! -----------------------------------------------------

       call write_log('ERROR: Unsupported value of whichartm', &
           GM_FATAL,__FILE__,__LINE__)

    end select
                
  end subroutine calcartm

  subroutine calcacab(numerics,climate,pddcalc,which,usrf,artm,arng,prcp,ablt,lati,acab,thck)

    use glimmer_global, only : dp, sp 
    use paramets, only : len0
    use glint_pdd
    use glide_types
    use glint_type
    use glimmer_log

    implicit none

    !-----------------------------------------------------------------------------
    ! Subroutine arguments
    !-----------------------------------------------------------------------------

    type(glide_numerics), intent(in)    :: numerics
    type(glint_climate), intent(in)     :: climate
    type(glint_pdd_params),  intent(inout) :: pddcalc
    integer,                intent(in)    :: which
    real(dp),dimension(:,:),intent(in)    :: usrf
    real(sp),dimension(:,:),intent(in)    :: artm
    real(sp),dimension(:,:),intent(in)    :: arng
    real(sp),dimension(:,:),intent(inout) :: prcp
    real(sp),dimension(:,:),intent(out)   :: ablt
    real(sp),dimension(:,:),intent(in)    :: lati
    real(sp),dimension(:,:),intent(out)   :: acab
    real(dp),dimension(:,:),intent(in)    :: thck

    !-----------------------------------------------------------------------------
    ! Internal variables
    !-----------------------------------------------------------------------------

    real(sp) :: grid
    real(sp) :: dist, ewct, nsct
    integer  :: ns,ew,ewn,nsn

    !-----------------------------------------------------------------------------

    ewn=size(prcp,1) ; nsn=size(prcp,2)

    !-----------------------------------------------------------------------------

    grid = numerics%dew * len0

    select case(which)
    case(0)

       ewct = real(ewn+1) / 2.0
       nsct = real(nsn+1) / 2.0

       do ns = 1,nsn
          do ew = 1,ewn
             dist = grid * sqrt((real(ew) - ewct)**2 + (real(ns) - nsct)**2) 
             acab(ew,ns) = amin1(climate%nmsb(1), climate%nmsb(2) * (climate%nmsb(3) - dist))
          end do
       end do

    case(1)
       call glint_pdd_mbal(pddcalc,artm,arng,prcp,lati,ablt,acab) 
    case(2) 
       acab = prcp
    case(3)
!       call enmabal
       call glint_pdd_mbal(pddcalc,artm,arng,prcp,lati,ablt,acab)
    case default
       call write_log('Unrecognised value of whichacab',GM_FATAL,__FILE__,__LINE__)
    end select

    ! Adjust ablation to be no greater than ice available for melting

    where (ablt*numerics%dt>(thck+prcp*numerics%dt)) 
      ablt=(thck/numerics%dt)+prcp
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
    ablt(:,nsn)=prcp(:,nsn)
    acab(:,nsn)=0.0
    ablt(1,:)=prcp(1,:)
    acab(1,:)=0.0
    ablt(ewn,:)=prcp(ewn,:)
    acab(ewn,:)=0.0

  end subroutine calcacab

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine calcprcp(which,prcp,up_rate,f1,xwind,ywind,artm,orog,dx,dy,presprcp, &
                      presartm,pfac)

    use glint_mbal

    integer,                intent(in)    :: which    !*FD Precipitation option
    real(sp),dimension(:,:),intent(inout) :: prcp     !*FD Large-scale precip on input,
                                                      !*FD desired precip field on output
    real(sp),               intent(in)    :: up_rate  !*FD user-specified uniform precip rate
    real(rk),               intent(in)    :: f1       !*FD scaling parameter
    real(rk),dimension(:,:),intent(in)    :: xwind    !*FD x-wind
    real(rk),dimension(:,:),intent(in)    :: ywind    !*FD y-wind
    real(sp),dimension(:,:),intent(in)    :: artm     !*FD air temperature
    real(rk),dimension(:,:),intent(in)    :: orog     !*FD local orography
    real(rk),               intent(in)    :: dx       !*FD x-spacing
    real(rk),               intent(in)    :: dy       !*FD y-spacing
    real(sp),dimension(:,:),intent(in)    :: presprcp !*FD present-day precip
    real(sp),dimension(:,:),intent(in)    :: presartm !*FD present-day air temp
    real(sp),               intent(in)    :: pfac     !*FD precip parameterisation factor

   select case(which)

    case(0)   ! Uniform precipitation ----------------------------------------

      prcp=up_rate/f1

    case(1)   ! Large scale precip -------------------------------------------
              ! Note that we / by 1000 to convert to meters, and then by
              ! f1 for scaling in the model.

      prcp=prcp/(1000.0*f1)

    case(2)   ! Precip parameterization --------------------------------------

      call glint_precip(prcp,xwind,ywind,artm,orog,dx,dy,fixed_a=.true.)
      prcp=prcp/(1000.0*f1)

    case(3)   ! Prescribed small-scale precip from file, ---------------------
              ! adjusted for temperature

      prcp = presprcp * pfac ** (artm - presartm)

    end select

  end subroutine calcprcp

end module glint_climate
