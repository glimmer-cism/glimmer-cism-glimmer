
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_temp.f90 - part of the GLIMMER ice model         + 
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

module glimmer_temp

  use glimmer_types

contains

  subroutine timeevoltemp(model,which,artm)

    !*FD Calculates the ice temperature, according to one
    !*FD of several alternative methods.

    use glimmer_global, only : dp
    use paramets,       only : thk0
    use glimmer_velo
    use glimmer_thck
    use glide_messages

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    type(glimmer_global_type),intent(inout) :: model       !*FD Ice model parameters.
    integer,                  intent(in)    :: which       !*FD Flag to choose method.
    real(sp),dimension(:,:),  intent(in)    :: artm        !*FD Surface air temperature

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    real(dp),dimension(size(model%numerics%sigma)) :: subd, diag, supd, rhsd
    real(dp),dimension(size(model%numerics%sigma)) :: prevtemp, iteradvt, diagadvt
    real(dp) :: tempresid

    integer :: iter, again
    integer :: ew,ns

    real(dp),parameter :: tempthres = 0.001d0, floatlim = 10.0d0 / thk0
    integer, parameter :: mxit = 100
    integer, parameter :: ewbc = 1, nsbc = 1 

    logical, allocatable, dimension(:,:) :: floater
    character(80) :: outtxt

    !------------------------------------------------------------------------------------
    ! ewbc/nsbc set the type of boundary condition aplied at the end of
    ! the domain. a value of 0 implies zero gradient.
    !------------------------------------------------------------------------------------
    ! Calculate the ice thickness according to different methods
    !------------------------------------------------------------------------------------

    select case(which)

    case(0) ! Set column to surface air temperature -------------------------------------

       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             model%temper%temp(:,ew,ns) = dmin1(0.0d0,dble(artm(ew,ns)))
          end do
       end do

    case(1) ! Do full temperature solution ---------------------------------------------

       ! Work out where the ice is floating, and assign to mask array -------------------

       allocate(floater(model%general%ewn,model%general%nsn))
       floater = abs(model%geometry%topg-model%geometry%lsrf) > floatlim

       ! Calculate time-derivatives of thickness and upper surface elevation ------------

       call timeders(model%thckwk,            &
            model%geometry%thck,     &
            model%geomderv%dthckdtm, &
            model%geometry%mask,     &
            model%numerics%time,     &
            1)

       call timeders(model%thckwk,            &
            model%geometry%usrf,     &
            model%geomderv%dusrfdtm, &
            model%geometry%mask,     &
            model%numerics%time,     &
            2)

       ! Calculate the vertical velocity of the grid ------------------------------------

       call gridwvel(model%numerics%sigma,  &
            model%numerics%thklim, &
            model%velocity%uvel,   &
            model%velocity%vvel,   &
            model%geomderv,        &
            model%geometry%thck,   &
            model%velocity%wgrd)

       ! Calculate the actual vertical velocity; method depends on whichwvel ------------

       select case(model%options%whichwvel)
       case(0) 

          ! Usual vertical integration

          call wvelintg(model%velocity%uvel,                        &
               model%velocity%vvel,                        &
               model%geomderv,                             &
               model%numerics,                             &
               model%velowk,                               &
               model%velocity%wgrd(model%general%upn,:,:), &
               model%geometry%thck,                        &
               model%temper%bmlt,                          &
               model%velocity%wvel)

       case(1)

          ! Vertical integration constrained so kinematic upper BC obeyed.

          call wvelintg(model%velocity%uvel,                        &
               model%velocity%vvel,                        &
               model%geomderv,                             &
               model%numerics,                             &
               model%velowk,                               &
               model%velocity%wgrd(model%general%upn,:,:), &
               model%geometry%thck,                        &
               model%temper%  bmlt,                        &
               model%velocity%wvel)

          call chckwvel(model%numerics,                             &
               model%geomderv,                             &
               model%velocity%uvel(1,:,:),                 &
               model%velocity%vvel(1,:,:),                 &
               model%velocity%wvel,                        &
               model%geometry%thck,                        &
               model%climate% acab)

       case default

          call glide_msg(GM_FATAL,__FILE__,__LINE__,'Unrecognised value of whichwvel')

       end select

       ! Allocate some arrays and initialise ----------------------------------------------

       allocate(model%tempwk%initadvt(model%general%upn,model%general%ewn,model%general%nsn))
       allocate(model%tempwk%inittemp(model%general%upn,model%general%ewn,model%general%nsn))
       allocate(model%tempwk%dissip(model%general%upn,model%general%ewn,model%general%nsn))

       model%tempwk%inittemp = 0.0d0

       ! ----------------------------------------------------------------------------------

       call finddisp(model, &
            model%geometry%thck, &
            model%geomderv%stagthck, &
            model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns, &
            model%temper%flwa)

       call hadvall(model, &
            model%temper%temp, &
            model%velocity%uvel, &
            model%velocity%vvel, &
            model%geometry%thck)

       iter = 0
       again = 0

       do while (again .eq. 0)

          tempresid = 0.0d0

          do ns = 2,model%general%nsn-1
             do ew = 2,model%general%ewn-1
                if(model%geometry%thck(ew,ns)>model%numerics%thklim) then

                   call hadvpnt(model%tempwk,                           &
                        iteradvt,                               &
                        diagadvt,                               &
                        model%temper%temp(:,ew-2:ew+2,ns),      &
                        model%temper%temp(:,ew,ns-2:ns+2),      &
                        model%velocity%uvel(:,ew-1:ew,ns-1:ns), &
                        model%velocity%vvel(:,ew-1:ew,ns-1:ns))

                   call findvtri(model,iter,ew,ns,subd,diag,supd,rhsd,iteradvt,diagadvt, &
                        model%temper%temp(:,ew,ns), &
                        model%velocity%wgrd(:,ew,ns), &
                        model%velocity%wvel(:,ew,ns), &
                        model%geometry%thck(ew,ns), &
                        artm(ew,ns), &
                        floater(ew,ns))

                   prevtemp = model%temper%temp(:,ew,ns)

                   call tridag(subd(2:model%general%upn), &
                        diag(1:model%general%upn), &
                        supd(1:model%general%upn-1), &
                        rhsd(1:model%general%upn), &
                        model%temper%temp(1:model%general%upn,ew,ns))

                   call corrpmpt(model%temper%temp(:,ew,ns),model%geometry%thck(ew,ns),model%temper%bwat(ew,ns), &
                        model%numerics%sigma,model%general%upn)

                   tempresid = dmax1(tempresid,maxval(abs(model%temper%temp(:,ew,ns)-prevtemp(:))))

                endif
             end do
          end do

          if (0 == iter) then
             deallocate(model%tempwk%initadvt)
          end if

          iter = iter + 1

          if (tempresid > tempthres .and. iter < mxit) then
             again = 0
          else 
             again = 1
          end if

       end do

       model%temper%niter = model%temper%niter + iter 

       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             if (model%geometry%thck(ew,ns) <= model%numerics%thklim) then
                model%temper%temp(:,ew,ns) = dmin1(0.0d0,dble(artm(ew,ns)))
             end if
          end do
       end do

       deallocate(model%tempwk%inittemp)

       call swapbndt(ewbc, &
            model%temper%temp(:,1,:), &
            model%temper%temp(:,2,:), &
            model%temper%temp(:,model%general%ewn,:), &
            model%temper%temp(:,model%general%ewn-1,:), &
            0)
       call swapbndt(nsbc, &
            model%temper%temp(:,:,1), &
            model%temper%temp(:,:,2), &
            model%temper%temp(:,:,model%general%nsn), &
            model%temper%temp(:,:,model%general%nsn-1), &
            1)

       ! Calculate basal melt rate --------------------------------------------------

       call calcbmlt(model, &
            model%temper%temp, &
            model%geometry%thck, &
            model%temper%bwat, &
            model%geomderv%stagthck, &
            model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns, &
            model%velocity%ubas, &
            model%velocity%vbas, &
            model%temper%bmlt, &
            floater)

       ! Calculate basal water depth ------------------------------------------------

       call calcbwat(model, &
            model%options%whichbwat, &
            model%temper%bmlt, &
            model%temper%bwat, &
            model%geometry%thck, &
            model%geometry%topg, &
            model%temper%temp(model%general%upn,:,:), &
            floater) 

       ! Deallocate arrays ----------------------------------------------------------

       deallocate(model%tempwk%dissip,floater)

    case(2) ! Do something else, unspecified ---------------------------------------

       model%temper%bwat = 0.0d0

       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             model%temper%temp(:,ew,ns) = dmin1(0.0d0,dble(artm(ew,ns))) * (1.0d0 - model%numerics%sigma)
             call corrpmpt(model%temper%temp(:,ew,ns),model%geometry%thck(ew,ns),model%temper%bwat(ew,ns),&
                  model%numerics%sigma,model%general%upn)
          end do
       end do

    case default

       call glide_msg(GM_FATAL,__FILE__,__LINE__,'Unrecognised value of whichtemp')

    end select

    ! Output some information ----------------------------------------------------

    write(outtxt,*)"temp ", model%numerics%time, iter, model%temper%niter, &
         real(model%temper%temp(model%general%upn,model%general%ewn/2+1,model%general%nsn/2+1))
    call glide_msg(GM_TIMESTEP,__FILE__,__LINE__,trim(adjustl(outtxt)))

       ! Calculate Glenn's A --------------------------------------------------------

       call calcflwa(model%numerics,        &
            model%velowk,          &
            model%paramets%fiddle, &
            model%temper%flwa,     &
            model%temper%temp,     &
            model%geometry%thck,   &
            model%options%whichflwa) 

  end subroutine timeevoltemp

!-----------------------------------------------------------------------------------

  subroutine calcartm(model,which,usrf,lati,artm,arng,g_orog,g_artm,g_arng)

    !*FD Calculate the surface air temperature and mean annual range,
    !*FD according to various different models.

    use glimmer_global, only : dp, sp 
    use paramets, only : len0, thk0
    use glide_messages

    use glimmer_mbal

    implicit none

    type(glimmer_global_type),       intent(inout) :: model  !*FD Model parameters to be used
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

    if (model%temper%first1) then
      model%temper%grid = model%numerics%dew * len0

      if (which == 0 .or. which == 4) then
        model%paramets%airt(2) = model%paramets%airt(2) * thk0
      end if

      model%temper%first1 = .false.
    end if

    !--------------------------------------------------------------------
    ! Now calculate the temp and range according to different methods:
    !--------------------------------------------------------------------

    select case(which)

            ! ----------------------------------------------------------
    case(0) ! Linear decrease from sea-level according to lapse-rate
            ! ----------------------------------------------------------

      artm = model%paramets%airt(1) + usrf * model%paramets%airt(2)

            ! ----------------------------------------------------------
    case(1) ! 2d EISMINT test - Cubic function of distance from domain centre
            ! ----------------------------------------------------------

      ewct = real(model%general%ewn+1) / 2.0
      nsct = real(model%general%nsn+1) / 2.0

      do ns = 1,model%general%nsn
        do ew = 1,model%general%ewn
          dist = model%temper%grid * sqrt((real(ew) - ewct)**2 + (real(ns) - nsct)**2) 
          artm(ew,ns) = model%paramets%airt(1) + dist**3 * model%paramets%airt(2)
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

      ewct = real(model%general%ewn+1) / 2.0; nsct = real(model%general%nsn+1) / 2.0

      do ns = 1,model%general%nsn
        do ew = 1,model%general%ewn
          dist = model%temper%grid * sqrt((real(ew) - ewct)**2 + (real(ns) - nsct)**2) 
          artm(ew,ns) = model%paramets%airt(1) + dist * model%paramets%airt(2)
        end do
      end do            

            ! ----------------------------------------------------------
    case(3) ! Obtaining appropriate temperature forcing from file
            ! ----------------------------------------------------------

      do while ( model%numerics%time .ge. model%forcdata%forcing(model%temper%tpt+1,1) .and. &
                 model%temper%tpt+1 .le. model%forcdata%flines)
         model%temper%tpt = model%temper%tpt + 1
         model%temper%perturb = model%forcdata%forcing(model%temper%tpt,2)
         ! ** print '(a26,f10.1,f10.3)', &
         ! **       '---> forcing read at      ', model%numerics%time, perturb 
      end do

      !   find air temps based on regr equation for greenland

      do ns = 1, model%general%nsn
        do ew = 1, model%general%ewn 
      
      !   find latitude dependent height of inversion layer  
 
          inve = max(0.0, 300 * (lati(ew,ns) - 65.0) / 15.0)
       
      !   make sure that surface is above sealevel

          esurf = max(0.0, real(usrf(ew,ns)*thk0,sp))
     
      !   find mean annual temp (is surface above inversion layer?)

          if ( esurf <= inve ) then
            artm(ew,ns) = 49.13 - 0.007992 * inve - 0.7576 * lati(ew,ns) + model%temper%perturb
          else
            artm(ew,ns) = 49.13 - 0.007992 * esurf - 0.7576 * lati(ew,ns) + model%temper%perturb
          end if

      !   find July temp (express as annual half range for convenience)
    
          arng(ew,ns)= 30.38 - 0.006277 * esurf - 0.3262 * lati(ew,ns) + model%temper%perturb - artm(ew,ns)

        end do
      end do

            ! ----------------------------------------------------------
    case(4) ! Air temperature is function of latitude and height
            ! ----------------------------------------------------------
    
      ! Note that we are using the array lati to hold sealevel air temperatures

      artm = lati + usrf * model%paramets%airt(2)

            ! ----------------------------------------------------------
    case(5) ! Uniform temperature, zero range
            ! ----------------------------------------------------------

      artm=model%climate%usurftemp
      arng=0.0

            ! ----------------------------------------------------------
    case(6) ! Uniform temperature, lapse-rate corrected, zero range
            ! ----------------------------------------------------------

      artm=model%climate%usurftemp
      arng=0.0
      call glimmer_lapserate(artm,real(usrf*thk0,rk),real(model%climate%ulapse_rate,rk))

            ! ----------------------------------------------------------
    case(7) ! Supplied large-scale temperature and range
            ! ----------------------------------------------------------

      ! Check we have the necessary arguments first

      if (present(g_arng).and.present(g_artm).and.present(g_orog)) then

      ! Copy the fields

        arng=g_arng
        artm=g_artm

      ! Reduce temperatures to sea-level

        call glimmer_lapserate(artm,real(g_orog,rk),real(-model%climate%ulapse_rate,rk))

      ! Raise them to high-res orography 

        call glimmer_lapserate(artm,real(usrf*thk0,rk),real(model%climate%ulapse_rate,rk))

      else

        call glide_msg(GM_FATAL,__FILE__,__LINE__,'Error in arguments to CALCARTM - stopping')

      endif

            ! ----------------------------------------------------------
    case(8) ! Leave everything alone
            ! ----------------------------------------------------------

                 ! -----------------------------------------------------
    case default ! Flag an error otherwise
                 ! -----------------------------------------------------

      call glide_msg(GM_FATAL,__FILE__,__LINE__,'ERROR: Unsupported value of whichartm')

    end select
                
  end subroutine calcartm

!-------------------------------------------------------------------------

  subroutine hadvpnt(tempwk,iteradvt,diagadvt,tempx,tempy,uvel,vvel)

    use glimmer_global, only : dp
    use glimmer_utils, only: hsum

    implicit none

    type(glimmer_tempwk) :: tempwk
    real(dp), dimension(:,:,:), intent(in) :: uvel, vvel
    real(dp), dimension(:,:), intent(in) :: tempx, tempy
    real(dp), dimension(:), intent(out) :: iteradvt, diagadvt

    real(dp), dimension(size(iteradvt)) :: u, v

    iteradvt = 0.0d0
    diagadvt = 0.0d0

    u = tempwk%advconst(1) * hsum(uvel(:,:,:))
    v = tempwk%advconst(2) * hsum(vvel(:,:,:))
  
    if (u(1) > 0.0d0) then
      iteradvt = u * (- 4.0d0*tempx(:,2) + tempx(:,1))
      diagadvt = u * 3.0d0
    else if (u(1) < 0.0d0) then
      iteradvt = u * (4.0d0*tempx(:,4) - tempx(:,5))
      diagadvt = - u * 3.0d0
    end if

    if (v(1) > 0.0d0) then
      iteradvt = iteradvt + v * (- 4.0d0*tempy(:,2) + tempy(:,1))
      diagadvt = diagadvt + v * 3.0d0
    else if (v(1) < 0.0d0) then
      iteradvt = iteradvt + v * (4.0d0*tempy(:,4) - tempy(:,5))
      diagadvt = diagadvt - v * 3.0d0
    end if

  end subroutine hadvpnt

!-------------------------------------------------------------------------

  subroutine fohadvpnt(tempwk,iteradvt,diagadvt,tempx,tempy,uvel,vvel)

    use glimmer_global, only : dp
    use glimmer_utils, only: hsum

    implicit none

    type(glimmer_tempwk) :: tempwk
    real(dp), dimension(:,:,:), intent(in) :: uvel, vvel
    real(dp), dimension(:,:), intent(in) :: tempx, tempy
    real(dp), dimension(:), intent(out) :: iteradvt, diagadvt
  
    real(dp), dimension(size(iteradvt)) :: u, v

    iteradvt = 0.0d0
    diagadvt = 0.0d0

    u = tempwk%advconst(1) * hsum(uvel(:,:,:))
    v = tempwk%advconst(2) * hsum(vvel(:,:,:))
    
    if (u(1) > 0.0d0) then
      iteradvt = - u * 2.0d0 * tempx(:,1)
      diagadvt = 2.0d0 * u 
    else if (u(1) < 0.0d0) then
      iteradvt = u * 2.0d0 * tempx(:,3)
      diagadvt = - 2.0d0 * u 
    end if

    if (v(1) > 0.0d0) then
      iteradvt = iteradvt - v * 2.0d0 * tempy(:,1) 
      diagadvt = diagadvt + 2.0d0 * v 
    else if (v(1) < 0.0d0) then
      iteradvt = iteradvt + v * 2.0d0 * tempy(:,3)
      diagadvt = diagadvt - 2.0d0 * v 
    end if

  end subroutine fohadvpnt

!-------------------------------------------------------------------------

  subroutine hadvall(model,temp,uvel,vvel,thck)

    use glimmer_global, only : dp 

    implicit none

    type(glimmer_global_type) :: model
    real(dp), dimension(:,0:,0:), intent(in) :: temp
    real(dp), dimension(:,:,:), intent(in) :: uvel, vvel
    real(dp), dimension(:,:), intent(in) :: thck

    real(dp), dimension(size(temp,dim=1)) :: diagadvt

    integer :: ew,ns

    if (model%tempwk%first1) then
      model%tempwk%advconst(1) = model%numerics%dttem / (16.0d0 * model%numerics%dew)
      model%tempwk%advconst(2) = model%numerics%dttem / (16.0d0 * model%numerics%dns)
      model%tempwk%first1 = .false.
    end if

    model%tempwk%initadvt = 0.0d0

    do ns = 2,model%general%nsn-1
      do ew = 2,model%general%ewn-1
        if (thck(ew,ns) > model%numerics%thklim) then

          call hadvpnt(model%tempwk,                   &
                       model%tempwk%initadvt(:,ew,ns), &
                       diagadvt,                       &
                       temp(:,ew-2:ew+2,ns),           &
                       temp(:,ew,ns-2:ns+2),           &
                       uvel(:,ew-1:ew,ns-1:ns),        &
                       vvel(:,ew-1:ew,ns-1:ns))
        end if
      end do
    end do

  end subroutine hadvall

!-------------------------------------------------------------------------

  subroutine findvtri(model,iter,ew,ns,subd,diag,supd,rhsd,iteradvt,diagadvt,temp,wgrd,wvel,thck,artm,float)

    use glimmer_global, only : dp, sp 
    use physcon, only : rhoi, shci, coni, scyr 
    use paramets, only : tim0, thk0, acc0 !, geot

    implicit none

    type(glimmer_global_type) :: model
    integer, intent(in) :: ew, ns, iter
    real(dp), dimension(:), intent(in) :: temp, wgrd, wvel, iteradvt, diagadvt
    real(dp), intent(in) :: thck
    real(sp), intent(in) :: artm 
    real(dp), dimension(:), intent(out) :: subd, diag, supd, rhsd
    logical, intent(in) :: float
    integer :: up

    real(dp) :: fact(3), dupnp1
    real(dp), dimension(size(model%numerics%sigma)) :: weff

    if (model%tempwk%first2) then

      allocate(model%tempwk%dups(model%general%upn,3))
      
      model%tempwk%dups = 0.0d0

      do up = 2, model%general%upn-1
        model%tempwk%dups(up,1) = (model%numerics%sigma(up+1) - model%numerics%sigma(up-1)) * &
                                  (model%numerics%sigma(up)   - model%numerics%sigma(up-1))
        model%tempwk%dups(up,2) = (model%numerics%sigma(up+1) - model%numerics%sigma(up-1)) *  &
                                  (model%numerics%sigma(up+1) - model%numerics%sigma(up))
        model%tempwk%dups(up,3) = model%numerics%sigma(up+1)  - model%numerics%sigma(up-1)
      end do

      model%tempwk%zbed = 1.0d0 / thk0
      model%tempwk%dupn = model%numerics%sigma(model%general%upn) - model%numerics%sigma(model%general%upn-1)
      model%tempwk%wmax = 5.0d0 * tim0 / (scyr * thk0)

      model%tempwk%cons = (/ 2.0d0 * tim0 * model%numerics%dttem * coni / (2.0d0 * rhoi * shci * thk0**2), &
              model%numerics%dttem / 2.0d0, &
              2.0d0 * tim0 * model%numerics%dttem * model%paramets%geot / (thk0 * rhoi * shci), &
              tim0 * acc0 * model%numerics%dttem * model%paramets%geot / coni /)

      model%tempwk%first2 = .false.

    end if

    fact(1) = model%tempwk%cons(1) / thck**2
    fact(2) = model%tempwk%cons(2) / thck

    weff = wvel - wgrd

    if (maxval(abs(weff)) > model%tempwk%wmax) then
      weff = 0.0d0
    end if

    subd(2:model%general%upn-1) = fact(2) * weff(2:model%general%upn-1) / &
                                  model%tempwk%dups(2:model%general%upn-1,3)

    supd(2:model%general%upn-1) = - subd(2:model%general%upn-1) - fact(1) / &
                                  model%tempwk%dups(2:model%general%upn-1,2)

    subd(2:model%general%upn-1) = subd(2:model%general%upn-1) - fact(1) / &
                                  model%tempwk%dups(2:model%general%upn-1,1)

    diag(2:model%general%upn-1) = 1.0d0 - subd(2:model%general%upn-1) &
                                  - supd(2:model%general%upn-1) &
                                  + diagadvt(2:model%general%upn-1)

    rhsd(1) = artm
    supd(1) = 0.0d0
    subd(1) = 0.0d0
    diag(1) = 1.0d0

    dupnp1 = model%tempwk%zbed / thck  

    supd(model%general%upn) = 0.0d0 
    subd(model%general%upn) = - model%tempwk%cons(1) / (thck**2 * model%tempwk%dupn * dupnp1)
    diag(model%general%upn) = 1.0d0 - subd(model%general%upn) + diagadvt(model%general%upn)

    if (iter == 0) then

      model%tempwk%inittemp(2:model%general%upn-1,ew,ns) = temp(2:model%general%upn-1) * &
                                 (2.0d0 - diag(2:model%general%upn-1)) &
                              - temp(1:model%general%upn-2) * subd(2:model%general%upn-1) &
                              - temp(3:model%general%upn) * supd(2:model%general%upn-1) & 
                              - model%tempwk%initadvt(2:model%general%upn-1,ew,ns) &
                              + model%tempwk%dissip(2:model%general%upn-1,ew,ns)

      model%tempwk%inittemp(model%general%upn,ew,ns) = temp(model%general%upn) * &
                                 (2.0d0 - diag(model%general%upn)) &
                          - temp(model%general%upn-1) * subd(model%general%upn) &
                          - model%tempwk%cons(3) / (thck * dupnp1) &
                          + model%tempwk%cons(4) * weff(model%general%upn) & 
                          - model%tempwk%initadvt(model%general%upn,ew,ns)  &
                          + model%tempwk%dissip(model%general%upn,ew,ns)

    ! *tp* model%tempwk%inittemp(2:model%general%upn-1,ew,ns) = temp(2:model%general%upn-1) &
    ! *tp*                           - model%tempwk%initadvt(2:model%general%upn-1,ew,ns) + model%tempwk%dissip(2:model%general%upn-1,ew,ns)  

    ! *tp* model%tempwk%inittemp(model%general%upn,ew,ns) = temp(model%general%upn) &
    ! *tp*                       - cons(3) / (thck * dupnp1) &
    ! *tp*                       + cons(4) * (wvel(model%general%upn) - wgrd(model%general%upn)) &
    ! *tp*                       - model%tempwk%initadvt(model%general%upn,ew,ns) + model%tempwk%dissip(model%general%upn,ew,ns)  

    end if

! now do the basal boundary
! for grounded ice, a heat flux is applied
! for floating ice, temperature held constant

    if (float) then

      diag(model%general%upn) = 1.0d0  

      if (iter == 0) then

      model%tempwk%inittemp(model%general%upn,ew,ns) = temp(model%general%upn) 

      end if
    else 

      diag(model%general%upn) = 1.0d0 - subd(model%general%upn) + diagadvt(model%general%upn)

      if (iter == 0) then

        model%tempwk%inittemp(model%general%upn,ew,ns) = temp(model%general%upn) * (2.0d0 - diag(model%general%upn)) &
                            - temp(model%general%upn-1) * subd(model%general%upn) &
                            - model%tempwk%cons(3) / (thck * dupnp1) &
                            + model%tempwk%cons(4) * weff(model%general%upn) & 
                            - model%tempwk%initadvt(model%general%upn,ew,ns) + model%tempwk%dissip(model%general%upn,ew,ns)

      end if

    end if

    rhsd(2:model%general%upn) = model%tempwk%inittemp(2:model%general%upn,ew,ns) - iteradvt(2:model%general%upn)

  end subroutine findvtri

!-----------------------------------------------------------------------

  subroutine finddisp(model,thck,stagthck,dusrfdew,dusrfdns,flwa)

    use glimmer_global, only : dp
    use glimmer_utils, only : hsum, lsum
    use physcon, only : rhoi, grav, shci, gn
    use paramets, only : vis0, len0, thk0, tim0

    implicit none

    type(glimmer_global_type) :: model
    real(dp), dimension(:,:), intent(in) :: thck, stagthck, dusrfdew, dusrfdns
    real(dp), dimension(:,:,:), intent(in) :: flwa

    integer, parameter :: p1 = gn + 1  
    integer :: ew,ns

    real(dp) :: c2

    integer :: iew, ins

! two methods of doing this. 
! 1. find dissipation at u-pts and then average
! 2. find dissipation at H-pts by averaging quantities from u-pts
! 2. works best for eismint divide (symmetry) but 1 likely to be better for full expts

    if (.true.) then

      if (model%tempwk%first3) then
        allocate(model%tempwk%c1(model%general%upn))
        model%tempwk%c1 = (model%numerics%sigma * rhoi * grav * thk0**2 / len0)**p1 * &
                          2.0d0 * vis0 * model%numerics%dttem * tim0 / (16.0d0 * rhoi * shci)
        model%tempwk%first3 = .false.
      end if

     model%tempwk%dissip = 0.0d0

     do ns = 2, model%general%nsn-1
        do ew = 2, model%general%ewn-1
          if (thck(ew,ns) > model%numerics%thklim) then
  
            c2 = (sum(stagthck(ew-1:ew,ns-1:ns))/4.0d0 * dsqrt(sum((dusrfdew(ew-1:ew,ns-1:ns))/4.0d0)**2 &
                  + sum((dusrfdns(ew-1:ew,ns-1:ns))/4.0d0)**2))**p1                      

            model%tempwk%dissip(:,ew,ns) = c2 * model%tempwk%c1 * (hsum(flwa(:,ew-1:ew+1,ns-1:ns+1)) + flwa(:,ew,ns) &
                        + lsum(flwa(:,ew-1:ew+1,ns)) + lsum(flwa(:,ew,ns-1:ns+1)))

          end if
        end do
      end do

    else

! old method based on u pts

      if (model%tempwk%first3) then
        allocate(model%tempwk%c1(model%general%upn))
        model%tempwk%c1 = (model%numerics%sigma * rhoi * grav * thk0**2 / len0)**p1 * &
                          2.0d0 * vis0 * model%numerics%dttem * tim0 / (16.0d0 * rhoi * shci)
        model%tempwk%first3 = .false.
      end if

      model%tempwk%dissip = 0.0d0

      do ns = 2, model%general%nsn-1
        do ew = 2, model%general%ewn-1
          if (thck(ew,ns) > model%numerics%thklim) then

            do ins = ns-1,ns
              do iew = ew-1,ew
                c2 = (stagthck(iew,ins) * dsqrt(dusrfdew(iew,ins)**2 + dusrfdns(iew,ins)**2))**p1                        
                model%tempwk%dissip(:,ew,ns) = model%tempwk%dissip(:,ew,ns) + c2 * hsum(flwa(:,iew:iew+1,ins:ins+1))
              end do
            end do

            model%tempwk%dissip(:,ew,ns) = model%tempwk%c1 * model%tempwk%dissip(:,ew,ns)

          end if
        end do
      end do

    end if

  end subroutine finddisp

!-----------------------------------------------------------------------------------

  subroutine calcbmlt(model,temp,thck,bwat,stagthck,dusrfdew,dusrfdns,ubas,vbas,bmlt,floater)

    use glimmer_global, only : dp 
    use physcon, only : rhoi, grav, lhci, coni, shci
    use paramets, only : vel0, thk0, len0, tim0

    implicit none 

    type(glimmer_global_type) :: model
    real(dp), dimension(:,0:,0:), intent(in) :: temp
    real(dp), dimension(:,:), intent(in) :: thck, bwat, stagthck, dusrfdew, dusrfdns, ubas, vbas  
    real(dp), dimension(:,:), intent(out) :: bmlt
    logical, dimension(:,:), intent(in) :: floater

    real(dp), dimension(size(model%numerics%sigma)) :: pmptemp
    real(dp) :: slterm, newmlt

    integer :: ewp, nsp,up,ew,ns

    if (model%tempwk%first4) then

      allocate(model%tempwk%dupa(model%general%upn),model%tempwk%dupb(model%general%upn))
      allocate(model%tempwk%dupc(model%general%upn))

      model%tempwk%dupc = (/ (model%numerics%sigma(2) - model%numerics%sigma(1)) / 2.0d0, &
              ((model%numerics%sigma(up+1) - model%numerics%sigma(up-1)) / 2.0d0, &
              up=2,model%general%upn-1), (model%numerics%sigma(model%general%upn) - &
                           model%numerics%sigma(model%general%upn-1)) / 2.0d0  /)
      model%tempwk%dupa = (/ 0.0d0, 0.0d0, &
              ((model%numerics%sigma(up) - model%numerics%sigma(up-1)) / &
              ((model%numerics%sigma(up-2) - model%numerics%sigma(up-1)) * &
              (model%numerics%sigma(up-2) - model%numerics%sigma(up))), &
               up=3,model%general%upn)  /)
      model%tempwk%dupb = (/ 0.0d0, 0.0d0, &
              ((model%numerics%sigma(up) - model%numerics%sigma(up-2)) / &
              ((model%numerics%sigma(up-1) - model%numerics%sigma(up-2)) * &
              (model%numerics%sigma(up-1) - model%numerics%sigma(up))), &
              up=3,model%general%upn)  /)

      model%tempwk%f = (/ tim0 * coni / (thk0**2 * lhci * rhoi), &
                          tim0 * model%paramets%geot / &
                                 (thk0 * lhci * rhoi), &
                          tim0 * thk0 * rhoi * shci /  &
                                 (thk0 * tim0 * model%numerics%dttem * lhci * rhoi), &
                          tim0 * thk0**2 * vel0 * grav * rhoi / &
                                 (4.0d0 * thk0 * len0 * rhoi * lhci) /)

      model%tempwk%first4 = .false.

    end if

    do ns = 2, model%general%nsn-1
      do ew = 2, model%general%ewn-1
        if (thck(ew,ns) > model%numerics%thklim .and. .not. floater(ew,ns)) then

          call calcpmpt(pmptemp,thck(ew,ns),model%numerics%sigma,model%general%upn)

          if (abs(temp(model%general%upn,ew,ns)-pmptemp(model%general%upn)) .lt. 0.001) then

            slterm = 0.0d0

            do nsp = ns-1,ns
              do ewp = ew-1,ew
                slterm = slterm - stagthck(ewp,nsp) * &
                  (dusrfdew(ewp,nsp) * ubas(ewp,nsp) + dusrfdns(ewp,nsp) * vbas(ewp,nsp))
              end do
            end do

            bmlt(ew,ns) = 0.0d0
            newmlt = model%tempwk%f(4) * slterm - model%tempwk%f(2) + model%tempwk%f(3) * &
              model%tempwk%dupc(model%general%upn) * &
              thck(ew,ns) * model%tempwk%dissip(model%general%upn,ew,ns)

            up = model%general%upn - 1
 
            do while (abs(temp(up,ew,ns)-pmptemp(up)) .lt. 0.001 .and. up .ge. 3)
              bmlt(ew,ns) = bmlt(ew,ns) + newmlt
              newmlt = model%tempwk%f(3) * model%tempwk%dupc(up) * thck(ew,ns) * model%tempwk%dissip(up,ew,ns)
              up = up - 1
            end do

            up = up + 1

            if (up == model%general%upn) then
              bmlt(ew,ns) = newmlt - &
                model%tempwk%f(1) * ( (temp(up-2,ew,ns) - pmptemp(up-2)) * model%tempwk%dupa(up) &
                + (temp(up-1,ew,ns) - pmptemp(up-1)) * model%tempwk%dupb(up) ) / thck(ew,ns) 
            else
              bmlt(ew,ns) = bmlt(ew,ns) + max(0.0d0, newmlt - &
                model%tempwk%f(1) * ( (temp(up-2,ew,ns) - pmptemp(up-2)) * model%tempwk%dupa(up) &
                + (temp(up-1,ew,ns) - pmptemp(up-1)) * model%tempwk%dupb(up) ) / thck(ew,ns)) 
            end if

      ! first-order version
      ! newmlt + (f(1) * (temp(up-1,ew,ns) - pmptemp(up-1))) / ((model%numerics%sigma(up) - &
      ! model%numerics%sigma(up-1)) * thck(ew,ns))

          else

            bmlt(ew,ns) = 0.0d0

          end if

        else

          bmlt(ew,ns) = 0.0d0
                     
        end if
      end do
    end do                                       

  end subroutine calcbmlt

!-------------------------------------------------------------------

  subroutine calcbwat(model,which,bmlt,bwat,thck,topg,btem,floater)

    use glimmer_global, only : dp 
    use physcon, only : rhoi, grav, rhow, scyr
    use paramets, only : tim0, thk0, len0 
    use glide_messages

    implicit none

    type(glimmer_global_type) :: model
    integer, intent(in) :: which

    real(dp), dimension(:,:), intent(inout) :: bwat
    real(dp), dimension(:,:), intent(in) :: bmlt, thck, topg, btem
    logical, dimension(:,:), intent(in) :: floater

    real(dp), dimension(2), parameter :: &
         blim = (/ 0.00001 / thk0, 0.001 / thk0 /)

    real(dp), parameter :: smthf = 0.01d0 
    real(dp) :: estimate, dwphidew, dwphidns, dwphi, pmpt, bave
    real(dp), dimension(:,:), allocatable :: wphi, bwatu, bwatv, fluxew, fluxns, bint, smth

    integer :: t_wat,ns,ew

    character(80) :: outtxt

    select case (which)  ! This is the value of calcbwat
    case(0)

       if (model%tempwk%first5) then
          model%paramets%hydtim = tim0 / (model%paramets%hydtim * scyr)
          estimate = 0.2d0 / model%paramets%hydtim
          call find_dt_wat(model%numerics%dttem,estimate,model%tempwk%dt_wat,model%tempwk%nwat) 

          ! ** print *, model%numerics%dttem*tim0/scyr, model%tempwk%dt_wat*tim0/scyr, model%tempwk%nwat

          model%tempwk%c = (/ model%tempwk%dt_wat, 1.0d0 - 0.5d0 * model%tempwk%dt_wat * model%paramets%hydtim, &
               1.0d0 + 0.5d0 * model%tempwk%dt_wat * model%paramets%hydtim, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /) 

          model%tempwk%first5 = .false.
       end if

       do t_wat = 1, model%tempwk%nwat
          do ns = 2,model%general%nsn-1
             do ew = 2,model%general%ewn-1

                if (model%numerics%thklim < thck(ew,ns) .and. .not. floater(ew,ns)) then
                   bwat(ew,ns) = (model%tempwk%c(1) * bmlt(ew,ns) + model%tempwk%c(2) * bwat(ew,ns)) / &
                        model%tempwk%c(3)
                   if (blim(1) > bwat(ew,ns)) then
                      bwat(ew,ns) = 0.0d0
                   end if
                else
                   bwat(ew,ns) = 0.0d0
                end if

             end do
          end do
       end do

       allocate(smth(model%general%ewn,model%general%nsn))

       do ns = 2,model%general%nsn-1
          do ew = 2,model%general%ewn-1

             if (blim(2) < bwat(ew,ns)) then
                smth(ew,ns) = bwat(ew,ns) + smthf * &
                     (bwat(ew-1,ns) + bwat(ew+1,ns) + bwat(ew,ns-1) + bwat(ew,ns+1) - 4.0d0 * bwat(ew,ns))
             else 
                smth(ew,ns) = bwat(ew,ns)
             end if

          end do
       end do

       bwat(2:model%general%ewn-1,2:model%general%nsn-1) = smth(2:model%general%ewn-1,2:model%general%nsn-1)
       deallocate(smth)

    case(1)

       if (model%tempwk%first5) then

          model%tempwk%watvel = model%paramets%hydtim * tim0 / (scyr * len0)
          estimate = (0.2d0 * model%tempwk%watvel) / min(model%numerics%dew,model%numerics%dns)
          call find_dt_wat(model%numerics%dttem,estimate,model%tempwk%dt_wat,model%tempwk%nwat) 

          write(outtxt,*) model%numerics%dttem*tim0/scyr, model%tempwk%dt_wat*tim0/scyr, model%tempwk%nwat
          call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,trim(adjustl(outtxt)))

          model%tempwk%c = (/ rhow * grav, rhoi * grav, 2.0d0 * model%numerics%dew, 2.0d0 * model%numerics%dns, &
               0.25d0 * model%tempwk%dt_wat / model%numerics%dew, 0.25d0 * model%tempwk%dt_wat / model%numerics%dns, &
               0.5d0 * model%tempwk%dt_wat / model%numerics%dew, 0.5d0 * model%tempwk%dt_wat / model%numerics%dns /)

          model%tempwk%first5 = .false. 

       end if

       ! ** add any melt_water

       bwat = max(0.0d0,bwat + model%numerics%dttem * bmlt)

       allocate(wphi(model%general%ewn,model%general%nsn))
       allocate(bwatu(model%general%ewn,model%general%nsn))
       allocate(bwatv(model%general%ewn,model%general%nsn)) 
       allocate(fluxew(model%general%ewn,model%general%nsn))
       allocate(fluxns(model%general%ewn,model%general%nsn))
       allocate(bint(model%general%ewn-1,model%general%nsn-1)) 

       ! ** split time evolution into steps to avoid CFL problems

       do t_wat = 1,model%tempwk%nwat

          ! ** find potential surface using paterson p112, eq 4
          ! ** if no ice then set to sea level or land surface potential
          ! ** if frozen then set high 

          do ns = 1,model%general%nsn
             do ew = 1,model%general%ewn
                if (model%numerics%thklim < thck(ew,ns) .and. .not. floater(ew,ns)) then
                   call calcpmptb(pmpt,thck(ew,ns))
                   if (btem(ew,ns) == pmpt) then
                      wphi(ew,ns) = model%tempwk%c(1) * (topg(ew,ns) + bwat(ew,ns)) + model%tempwk%c(2) * thck(ew,ns)
                   else
                      wphi(ew,ns) = model%tempwk%c(1) * (topg(ew,ns) + thck(ew,ns))
                   end if
                else 
                   wphi(ew,ns) = max(model%tempwk%c(1) * topg(ew,ns),0.0d0)
                end if
             end do
          end do

          ! ** determine x,y components of water velocity assuming
          ! ** contstant velocity magnitude and using potential
          ! ** to determine direction

          do ns = 2,model%general%nsn-1
             do ew = 2,model%general%ewn-1
                if (thck(ew,ns) > model%numerics%thklim) then

                   dwphidew = (wphi(ew+1,ns) - wphi(ew-1,ns)) / model%tempwk%c(3)       
                   dwphidns = (wphi(ew,ns+1) - wphi(ew,ns-1)) / model%tempwk%c(4)  

                   dwphi = - model%tempwk%watvel / sqrt(dwphidew**2 + dwphidns**2)

                   bwatu(ew,ns) = dwphi * dwphidew  
                   bwatv(ew,ns) = dwphi * dwphidns  

                else
                   bwatu(ew,ns) = 0.0d0
                   bwatv(ew,ns) = 0.0d0
                end if
             end do
          end do

          ! ** use two-step law wendroff to solve dW/dt = -dF/dx - dF/dy

          ! ** 1. find fluxes F=uW

          fluxew = bwat * bwatu
          fluxns = bwat * bwatv

          ! ** 2. do 1st LW step on staggered grid for dt/2

          do ns = 1,model%general%nsn-1
             do ew = 1,model%general%ewn-1

                bave = 0.25 * sum(bwat(ew:ew+1,ns:ns+1))

                if (bave > 0.0d0) then

                   bint(ew,ns) = bave - &
                        model%tempwk%c(5) * (sum(fluxew(ew+1,ns:ns+1)) - sum(fluxew(ew,ns:ns+1))) - &
                        model%tempwk%c(6) * (sum(fluxns(ew:ew+1,ns+1)) - sum(fluxns(ew:ew+1,ns)))

                else
                   bint(ew,ns) = 0.0d0
                end if
             end do
          end do

          ! ** 3. find fluxes F=uW on staggered grid griven new Ws

          fluxew(1:model%general%ewn-1,1:model%general%nsn-1) = bint * 0.25 * &
               (bwatu(1:model%general%ewn-1,1:model%general%nsn-1) + &
               bwatu(2:model%general%ewn,1:model%general%nsn-1) + &
               bwatu(1:model%general%ewn-1,2:model%general%nsn) + &
               bwatu(2:model%general%ewn,2:model%general%nsn))
          fluxns(1:model%general%ewn-1,1:model%general%nsn-1) = bint * 0.25 * &
               (bwatv(1:model%general%ewn-1,1:model%general%nsn-1) + &
               bwatv(2:model%general%ewn,1:model%general%nsn-1) + &
               bwatv(1:model%general%ewn-1,2:model%general%nsn) + &
               bwatv(2:model%general%ewn,2:model%general%nsn))

          ! ** 4. finally do 2nd LW step to get back on to main grid

          do ns = 2,model%general%nsn-1
             do ew = 2,model%general%ewn-1
                if (bwat(ew,ns) > 0.0d0) then

                   bwat(ew,ns) = bwat(ew,ns) - &
                        model%tempwk%c(7) * (sum(fluxew(ew,ns-1:ns)) - sum(fluxew(ew-1,ns-1:ns))) - &
                        model%tempwk%c(8) * (sum(fluxns(ew-1:ew,ns)) - sum(fluxns(ew-1:ew,ns-1)))

                else
                   bwat(ew,ns) = 0.0d0
                end if
             end do
          end do
       end do

       where (blim(1) > bwat) 
          bwat = 0.0d0
       end where

       deallocate(bint,fluxew,fluxns,bwatu,bwatv,wphi)

    case(2)

       bwat = 0.0d0

    case default

       call glide_msg(GM_FATAL,__FILE__,__LINE__,'Unrecognised value of whichbwat')

    end select

    ! How to call the flow router.
    ! call advectflow(bwat,phi,bmlt,model%geometry%mask)

  contains

    subroutine find_dt_wat(dttem,estimate,dt_wat,nwat)

      implicit none

      real(dp), intent(out) :: dt_wat
      integer, intent(out) :: nwat
      real(dp), intent(in) :: dttem, estimate

      nwat = int(dttem/estimate) + 1
      dt_wat = dttem / nwat

    end subroutine find_dt_wat

  end subroutine calcbwat

!-------------------------------------------------------------------

  subroutine tridag(a,b,c,r,u)
 
    use glimmer_global, only : dp

    implicit none

    real(dp), dimension(:), intent(in) :: a,b,c,r
    real(dp), dimension(:), intent(out) :: u

    real(dp), dimension(size(b)) :: gam
    integer :: n, j
    real(dp) :: bet

    n=size(b)

    bet = b(1); u(1) = r(1)/bet

    do j = 2,n
      gam(j) = c(j-1) / bet
      bet = b(j) - a(j-1) * gam(j)
      u(j) = (r(j)- a(j-1) * u(j-1)) / bet
    end do

    do j = n-1,1,-1
      u(j) = u(j) - gam(j+1) * u(j+1)
    end do

  end subroutine tridag

!-------------------------------------------------------------------

  subroutine corrpmpt(temp,thck,bwat,sigma,upn)

    use glimmer_global, only : dp

    implicit none 

    real(dp), dimension(:), intent(inout) :: temp
    real(dp), intent(in) :: thck, bwat
    integer,intent(in) :: upn
    real(dp),dimension(:),intent(in) :: sigma

    real(dp), dimension(:) :: pmptemp(size(temp))

! corrects a temperature column for melting point effects
! 1. if temperature at any point in column is above pressure melting point then 
! set temperature to pressure melting point 
! 2. if bed is wet set basal temperature to pressure melting point 

    call calcpmpt(pmptemp,thck,sigma,upn)

    temp = dmin1(temp,pmptemp)

    if (bwat > 0.0d0) temp(upn) = pmptemp(upn)
                                
  end subroutine corrpmpt

!-------------------------------------------------------------------

  subroutine calcpmpt(pmptemp,thck,sigma,upn)

    use glimmer_global, only : dp !, upn
    use physcon, only : rhoi, grav, pmlt 
    use paramets, only : thk0

    implicit none 

    real(dp), dimension(:), intent(out) :: pmptemp
    real(dp), intent(in) :: thck
    real(dp),intent(in),dimension(:) :: sigma
    integer, intent(in) :: upn

    real(dp), parameter :: fact = - grav * rhoi * pmlt * thk0

    pmptemp = fact * thck * sigma

  end subroutine calcpmpt

!-------------------------------------------------------------------

  subroutine calcpmptb(pmptemp,thck)

    use glimmer_global, only : dp
    use physcon, only : rhoi, grav, pmlt 
    use paramets, only : thk0

    implicit none 

    real(dp), intent(out) :: pmptemp
    real(dp), intent(in) :: thck
  
    real(dp), parameter :: fact = - grav * rhoi * pmlt * thk0

    pmptemp = fact * thck 

  end subroutine calcpmptb

!-------------------------------------------------------------------

  subroutine swchpnt(a,b,c,d,e)

    implicit none 

    integer, intent(inout) :: a, b, c 
    integer, intent(in) :: d, e

    if (a == d) then
      a = e
      b = d
      c = -1
    else
      a = d
      b = e
      c = 1
    end if

  end subroutine swchpnt

!-------------------------------------------------------------------

  subroutine swapbndt(bc,a,b,c,d,flag)

    use glimmer_global, only : dp

    implicit none

    real(dp), intent(out), dimension(:,:) :: a, c
    real(dp), intent(in), dimension(:,:) :: b, d
    integer, intent(in) :: bc, flag

    if (bc == 0) then
      a = b
      c = d
    end if

  end subroutine swapbndt

!-------------------------------------------------------------------

  subroutine advectflow(model,bwat,phi,bmlt,mask,unit)

    use glimmer_global, only : dp 

    implicit none

    type(glimmer_global_type) :: model
    real(dp), intent(inout), dimension(:,:) :: phi, bwat, bmlt
    integer, intent(in), dimension(:,:) :: mask
    integer, intent(in) :: unit

    real(dp), allocatable, dimension(:,:) :: flats, r

    integer :: ewp, nsp, ct, pt, tct
    real(dp) :: dd(9), diffs(model%general%ewn,model%general%nsn,9)
    real(dp) :: sumd, dsq
    real(dp) :: da 
    integer :: ew,ns

    allocate(flats(model%general%ewn,model%general%nsn))

! get rid of pools and pits within the phi array c so that water
! is not collected
   
    call fillholes(model,phi,flats,mask)
 
! the main body of the calcs. finding the direction of the vector
! for each segment and then the angles and then the max vector

    dsq = sqrt(model%numerics%dew*model%numerics%dew + model%numerics%dns*model%numerics%dns); 
    dd = 1.0d0 / (/ dsq, model%numerics%dns, dsq, model%numerics%dew, 1.0d0, model%numerics%dew, &
                  dsq, model%numerics%dns, dsq /)
    dd(5) = 0.0d0
    dd((/1,3,7,9/)) = dd((/1,3,7,9/)) * sqrt(2.0d0) / 2.0d0

    da = model%numerics%dew * model%numerics%dns

    do ns = 2,model%general%nsn-1
      do ew = 2,model%general%ewn-1
        if (flats(ew,ns) == 0 .and. phi(ew,ns) /= model%tempwk%noflow) then

          ct = 1

          do nsp = ns-1,ns+1
            do ewp = ew-1,ew+1
              diffs(ew,ns,ct) = max(0.0d0,(phi(ew,ns) - phi(ewp,nsp)) * dd(ct)) ! fiddle factor ** 1.15
              ct = ct + 1
            end do
          end do

          sumd = sum(diffs(ew,ns,:))
          diffs(ew,ns,:) = diffs(ew,ns,:) / sumd
  
        end if
      end do
    end do

! Sort data into highest to lowest

    deallocate(flats)

    tct = count(phi(2:model%general%ewn-1,2:model%general%nsn-1) /= model%tempwk%noflow)
    allocate(r(tct,3))
    call phisorder(model,phi,r)
 
! calculate the water flow

    bwat = 0.0
   
    do pt = 1,tct

      ew = nint(r(pt,2))
      ns = nint(r(pt,3))

      bwat(ew,ns) = bwat(ew,ns) + da * bmlt(ew,ns)
   
      ct = 1
      do nsp = ns-1,ns+1
        do ewp = ew-1,ew+1
          bwat(ewp,nsp) = bwat(ewp,nsp) + bwat(ew,ns) * diffs(ew,ns,ct)
          ct = ct + 1
        end do
      end do

    end do

    open(unit,file='bwat')
    write(unit,'(g12.3)') bwat
    close(unit)

    deallocate(r)
      
  end subroutine advectflow

!-----------------------------------------------------------------------------------

  subroutine phisorder(model,var,r)

    use glimmer_global, only : dp
  
    implicit none

    type(glimmer_global_type) :: model
    real(dp), intent(in), dimension(:,:) :: var
    real(dp), intent(out), dimension(:,:) :: r
  
    integer :: ct,ns,ew

    ct = 1

    do ns = 2,model%general%nsn-1
      do ew = 2,model%general%ewn-1
        if (var(ew,ns) /= model%tempwk%noflow) then

          r(ct,1:3) = (/ var(ew,ns), dble(ew), dble(ns) /)
          ct = ct + 1

        end if
      end do
    end do

    call sorttp

  contains

    subroutine sorttp
  
      use glimmer_utils

      implicit none

      integer, dimension(:), allocatable :: index
      integer :: i

      allocate(index(size(r(:,1))))

      call indexx(-real(r(:,1),sp),index)

      do i = 2,3
        r(:,i) = r(index,i)
      end do
    
      deallocate(index)

    end subroutine sorttp
  
  end subroutine phisorder 

!-----------------------------------------------------------------------------------

  subroutine fillholes(model,phi,flats,mask)

    use glimmer_global, only : dp 

    implicit none

    type (glimmer_global_type) :: model
    real(dp), intent(inout), dimension(:,:) :: phi, flats
    integer, intent(in), dimension(:,:) :: mask

    real(dp), allocatable, dimension(:,:) :: old_phi
    integer, allocatable, dimension(:,:) :: pool

    real(dp) :: pvs(9), max_val
    real(dp), parameter :: null = 1e+20
    integer :: flag,ns,ew

    allocate(pool(model%general%ewn,model%general%nsn))
    allocate(old_phi(model%general%ewn,model%general%nsn))
    
    flag = 1

    do while (flag .eq. 1)

      flag = 0

      old_phi = phi

      do ew=2,model%general%ewn-1
        do ns=2,model%general%nsn-1; 

          flats(ew,ns) = 0

          if (mask(ew,ns) .eq. 1) then

            if (any(old_phi(ew-1:ew+1,ns-1:ns+1) < old_phi(ew,ns))) then
              pool(ew,ns) = 0
            else
              pool(ew,ns) = 1
            end if

            if (pool(ew,ns) .eq. 1) then
              flag = 1
              pvs = (/ old_phi(ew-1:ew+1,ns-1), old_phi(ew-1:ew+1,ns+1), old_phi(ew-1:ew+1,ns) /)
              where (pvs == old_phi(ew,ns))
                pvs = null
              end where

              max_val = minval(pvs)

              if (max_val .ne. null) then
                phi(ew,ns) = max_val
              else
                flag = 0
                flats(ew,ns) = 1
              end if
            end if
          end if
        end do
      end do
    end do

    deallocate(pool,old_phi)

  end subroutine fillholes

end module glimmer_temp

