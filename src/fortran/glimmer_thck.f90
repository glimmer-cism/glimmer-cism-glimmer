! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_thck.f90 - part of the GLIMMER ice model         + 
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

module glide_thck

  use glide_types

contains

  subroutine init_thck(model)
    !*FD initialise work data for ice thickness evolution
    implicit none
    type(glide_global_type) :: model

    allocate(model%thckwk%float(model%general%ewn,model%general%nsn))
    
    model%pcgdwk%fc2 = (/ model%numerics%alpha * model%numerics%dt / (2.0d0 * model%numerics%dew * model%numerics%dew), &
         model%numerics%dt, (1.0d0-model%numerics%alpha) / model%numerics%alpha, &
         1.0d0 / model%numerics%alpha, model%numerics%alpha * model%numerics%dt / &
         (2.0d0 * model%numerics%dns * model%numerics%dns), 0.0d0 /) 
  end subroutine init_thck

  subroutine thck_lin_evolve(model,newtemps,logunit)

    !*FD this subroutine solves the linearised ice thickness equation by computing the
    !*FD diffusivity from quantities of the previous time step

    use glide_velo
    implicit none
    ! subroutine arguments
    type(glide_global_type) :: model
    logical, intent(in) :: newtemps                     !*FD true when we should recalculate Glen's A
    integer,intent(in) :: logunit                       !*FD unit for logging

    if (model%geometry%empty) then

       model%geometry%thck = dmax1(0.0d0,model%geometry%thck + model%climate%acab * model%pcgdwk%fc2(2))
       print *, "* thck empty - net accumulation added", model%numerics%time

    else

       ! calculate basal velos
       if (newtemps) then
          call slipvelo(model%numerics,                &
               model%velowk,                  &
               model%geomderv,                &
               (/1,model%options%whichbtrc/), &
               model%temper%   bwat,          &
               model%velocity% btrc,          &
               model%geometry% relx,          &
               model%velocity% ubas,          &
               model%velocity% vbas)
       end if
       call slipvelo(model%numerics,                &
            model%velowk,                  &
            model%geomderv,                &
            (/2,model%options%whichbtrc/), &
            model%temper%   bwat,          &
            model%velocity% btrc,          &
            model%geometry% relx,          &
            model%velocity% ubas,          &
            model%velocity% vbas)

       ! calculate Glen's A if necessary
       if (newtemps) then
          call velo_integrate_flwa(model%velowk,model%geomderv%stagthck,model%temper%flwa)
       end if
       ! calculate diffusivity
       call velo_calc_diffu(model%velowk,model%geomderv%stagthck,model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns,model%velocity%diffu)

       ! get new thicknesses
       call thck_evolve(model,logunit,.true.,model%geometry%thck,model%geometry%thck)

       ! calculate horizontal velocity field
       call slipvelo(model%numerics,                &
            model%velowk,                  &
            model%geomderv,                &
            (/3,model%options%whichbtrc/), &
            model%temper%bwat,             &
            model%velocity%btrc,           &
            model%geometry%relx,           &
            model%velocity%ubas,           &
            model%velocity%vbas)
       call velo_calc_velo(model%velowk,model%geomderv%stagthck,model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns,model%temper%flwa,model%velocity%diffu,model%velocity%ubas, &
            model%velocity%vbas,model%velocity%uvel,model%velocity%vvel,model%velocity%uflx,model%velocity%vflx)
    end if
  end subroutine thck_lin_evolve

  subroutine thck_nonlin_evolve(model,newtemps,logunit)

    !*FD this subroutine solves the ice thickness equation by doing an outer, 
    !*FD non-linear iteration to update the diffusivities and in inner, linear
    !*FD iteration to calculate the new ice thickness distrib

    use glimmer_global, only : dp
    use glide_velo
    implicit none
    ! subroutine arguments
    type(glide_global_type) :: model
    logical, intent(in) :: newtemps                     !*FD true when we should recalculate Glen's A
    integer,intent(in) :: logunit                       !*FD unit for logging

    ! local variables
    integer, parameter :: pmax=50                       !*FD maximum Picard iterations
    real(kind=dp), parameter :: tol=1.0d-6
    real(kind=dp) :: residual
    integer p
    logical first_p


    if (model%geometry%empty) then

       model%geometry%thck = dmax1(0.0d0,model%geometry%thck + model%climate%acab * model%pcgdwk%fc2(2))
       print *, "* thck empty - net accumulation added", model%numerics%time

    else

       ! calculate basal velos
       if (newtemps) then
          call slipvelo(model%numerics,                &
               model%velowk,                  &
               model%geomderv,                &
               (/1,model%options%whichbtrc/), &
               model%temper%   bwat,          &
               model%velocity% btrc,          &
               model%geometry% relx,          &
               model%velocity% ubas,          &
               model%velocity% vbas)
          ! calculate Glen's A if necessary
          call velo_integrate_flwa(model%velowk,model%geomderv%stagthck,model%temper%flwa)
       end if

       first_p = .true.
       model%thckwk%oldthck = model%geometry%thck
       ! do Picard iteration
       do p=1,pmax
          model%thckwk%oldthck2 = model%thckwk%oldthck

          call stagvarb(model%geometry% thck, &
               model%geomderv% stagthck,&
               model%general%  ewn, &
               model%general%  nsn)

          call geomders(model%numerics, &
               model%geometry% usrf, &
               model%geomderv% stagthck,&
               model%geomderv% dusrfdew, &
               model%geomderv% dusrfdns)

          call geomders(model%numerics, &
               model%geometry% thck, &
               model%geomderv% stagthck,&
               model%geomderv% dthckdew, &
               model%geomderv% dthckdns)


          call slipvelo(model%numerics,                &
               model%velowk,                  &
               model%geomderv,                &
               (/2,model%options%whichbtrc/), &
               model%temper%   bwat,          &
               model%velocity% btrc,          &
               model%geometry% relx,          &
               model%velocity% ubas,          &
               model%velocity% vbas)

          ! calculate diffusivity
          call velo_calc_diffu(model%velowk,model%geomderv%stagthck,model%geomderv%dusrfdew, &
               model%geomderv%dusrfdns,model%velocity%diffu)

          ! get new thicknesses
          call thck_evolve(model,logunit,first_p,model%geometry%thck,model%thckwk%oldthck)

          first_p = .false.
          residual = maxval(abs(model%thckwk%oldthck-model%thckwk%oldthck2))
          #ifdef DEBUG
          !write(*,*) '*Picard iter ',p,residual
          #endif
          if (residual.le.tol) then
             exit
          end if

       end do

       model%geometry%thck = model%thckwk%oldthck

       ! calculate horizontal velocity field
       call slipvelo(model%numerics,                &
            model%velowk,                  &
            model%geomderv,                &
            (/3,model%options%whichbtrc/), &
            model%temper%bwat,             &
            model%velocity%btrc,           &
            model%geometry%relx,           &
            model%velocity%ubas,           &
            model%velocity%vbas)
       call velo_calc_velo(model%velowk,model%geomderv%stagthck,model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns,model%temper%flwa,model%velocity%diffu,model%velocity%ubas, &
            model%velocity%vbas,model%velocity%uvel,model%velocity%vvel,model%velocity%uflx,model%velocity%vflx)
    end if
  end subroutine thck_nonlin_evolve

  subroutine thck_evolve(model,logunit,calc_rhs,old_thck,new_thck)

    !*FD set up sparse matrix and solve matrix equation to find new ice thickness distribution
    !*FD this routine does not override the old thickness distribution

    use glimmer_global, only : dp
    use paramets, only : thk0,vel0
    implicit none
    ! subroutine arguments
    type(glide_global_type) :: model
    integer,intent(in) :: logunit                       !*FD unit for logging
    logical,intent(in) :: calc_rhs                      !*FD set to true when rhs should be calculated 
                                                        !*FD i.e. when doing lin solution or first picard iteration
    real(dp), intent(in), dimension(:,:) :: old_thck    !*FD contains ice thicknesses from previous time step
    real(dp), intent(inout), dimension(:,:) :: new_thck !*FD on entry contains first guess for new ice thicknesses
                                                        !*FD on exit contains ice thicknesses of new time step

    ! local variables
    real(dp), dimension(5) :: sumd 
    real(dp) :: err
    real(dp), parameter :: tolbnd = 1.0d-6
    integer :: thckflag = 0

    integer :: linit
    integer, parameter :: mxtbnd = 10, ewbc = 1, nsbc = 1

    ! ewbc/nsbc set the type of boundary condition aplied at the end of
    ! the domain. a value of 0 implies zero gradient.

    integer :: ew,ns

       !* the number of grid points and the number of nonzero 
       !* matrix elements (including bounary points)
       model%pcgdwk%pcgsize = (/ model%geometry%totpts, model%geometry%totpts * 5 /)

       model%pcgdwk%ct = 1

       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn 

             if (model%geometry%mask(ew,ns) /= 0) then

                if (ew == 1 .or. ew == model%general%ewn .or. &
                     ns == 1 .or. &
                     ns == model%general%nsn) then

                   call putpcgc(model%pcgdwk,&
                        1.0d0, &
                        model%geometry%mask(ew,ns), &
                        model%geometry%mask(ew,ns))           ! point (ew,ns)
                   if (calc_rhs) then
                      model%pcgdwk%rhsd(model%geometry%mask(ew,ns)) = old_thck(ew,ns) 
                   end if
                   model%pcgdwk%answ(model%geometry%mask(ew,ns)) = new_thck(ew,ns) 

                else if (any((model%geometry%thck(ew-1:ew+1,ns-1:ns+1) == 0.0d0)) .and. thckflag == 5) then

                   call putpcgc(model%pcgdwk, &
                        1.0d0, &
                        model%geometry%mask(ew,ns), &
                        model%geometry%mask(ew,ns))           ! point (ew,ns)

                   if (calc_rhs) then
                      model%pcgdwk%rhsd(model%geometry%mask(ew,ns)) = old_thck(ew,ns) + &
                           model%climate%acab(ew,ns) * model%pcgdwk%fc2(2) 
                   end if
                   model%pcgdwk%answ(model%geometry%mask(ew,ns)) = new_thck(ew,ns) 

                else

                   call findsums(ew,ns,model%pcgdwk%fc2)

                   call putpcgc(model%pcgdwk,sumd(1), &
                        model%geometry%mask(ew-1,ns), &
                        model%geometry%mask(ew,ns))       ! point (ew-1,ns)
                   call putpcgc(model%pcgdwk,sumd(2), &
                        model%geometry%mask(ew+1,ns), &
                        model%geometry%mask(ew,ns))       ! point (ew+1,ns)
                   call putpcgc(model%pcgdwk,sumd(3), &
                        model%geometry%mask(ew,ns-1), &
                        model%geometry%mask(ew,ns))       ! point (ew,ns-1)
                   call putpcgc(model%pcgdwk,sumd(4), &
                        model%geometry%mask(ew,ns+1), &
                        model%geometry%mask(ew,ns))       ! point (ew,ns+1)
                   call putpcgc(model%pcgdwk,1.0d0 + sumd(5), &
                        model%geometry%mask(ew,ns), &
                        model%geometry%mask(ew,ns)) ! point (ew,ns)

                   if (calc_rhs) then
                      model%pcgdwk%rhsd(model%geometry%mask(ew,ns)) = old_thck(ew,ns) * &
                           (1.0d0 - model%pcgdwk%fc2(3) * sumd(5)) - model%pcgdwk%fc2(3) * &
                           (old_thck(ew-1,ns) * sumd(1) + old_thck(ew+1,ns) * sumd(2) + &
                           old_thck(ew,ns-1) * sumd(3) + old_thck(ew,ns+1) * sumd(4)) - &
                           model%pcgdwk%fc2(4) * (model%geometry%lsrf(ew,ns) * sumd(5) + &
                           model%geometry%lsrf(ew-1,ns) * sumd(1) + model%geometry%lsrf(ew+1,ns) * sumd(2) + &
                           model%geometry%lsrf(ew,ns-1) * sumd(3) + model%geometry%lsrf(ew,ns+1) * sumd(4)) +  &
                           model%climate%acab(ew,ns) * model%pcgdwk%fc2(2)
                   end if
                   model%pcgdwk%answ(model%geometry%mask(ew,ns)) = new_thck(ew,ns) 

                end if
             end if
          end do
       end do

       model%pcgdwk%pcgsize(2) = model%pcgdwk%ct - 1 

       call slapsolv(model,.true.,linit,err,logunit)   

       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn 

             if (model%geometry%mask(ew,ns) /= 0) then
                new_thck(ew,ns) = model%pcgdwk%answ(model%geometry%mask(ew,ns))
             end if

          end do
       end do

       ! *tp+* implement bcs

       call swapbndh(ewbc, &
            new_thck(1,:), new_thck(2,:), &
            new_thck(model%general%ewn,:), new_thck(model%general%ewn-1,:))
       call swapbndh(nsbc, &
            new_thck(:,1), new_thck(:,2), &
            new_thck(:,model%general%nsn), new_thck(:,model%general%nsn-1))

       model%pcgdwk%tlinit = model%pcgdwk%tlinit + linit
       model%pcgdwk%mlinit = max(linit,model%pcgdwk%mlinit)

       new_thck = max(0.0d0, new_thck)
#ifdef DEBUG
       print *, "* thck ", model%numerics%time, linit, model%pcgdwk%mlinit, model%pcgdwk%tlinit, model%geometry%totpts, &
            real(thk0*new_thck(model%general%ewn/2+1,model%general%nsn/2+1)), &
            real(vel0*maxval(abs(model%velocity%ubas))), real(vel0*maxval(abs(model%velocity%vbas))) 
#endif

  contains
    
    subroutine findsums(ew,ns,fc)

      implicit none

      integer, intent(in) :: ew, ns
      real(dp),intent(in),dimension(6) :: fc

      sumd(1) = fc(1) * (sum(model%velocity%diffu(ew-1,ns-1:ns)) + &
           sum(model%velocity%ubas(ew-1,ns-1:ns)))
      sumd(2) = fc(1) * (sum(model%velocity%diffu(ew,ns-1:ns)) + &
           sum(model%velocity%ubas(ew,ns-1:ns)))
      sumd(3) = fc(5) * (sum(model%velocity%diffu(ew-1:ew,ns-1)) + &
           sum(model%velocity%ubas(ew-1:ew,ns-1)))
      sumd(4) = fc(5) * (sum(model%velocity%diffu(ew-1:ew,ns)) + &
           sum(model%velocity%ubas(ew-1:ew,ns)))
      sumd(5) = - sum(sumd(1:4))

    end subroutine findsums
  end subroutine thck_evolve

!---------------------------------------------------------------------------------

  subroutine slapsolv(model,first,iter,err,lunit)

    use glimmer_global, only : dp 
    use glide_stop
  !use pcgdwk 
  !use funits, only : ulog

    implicit none

    type(glide_global_type) :: model
    logical, intent(in) :: first
    integer, intent(out) :: iter
    real(dp), intent(out) :: err
    integer,intent(in) :: lunit

    real(dp), dimension(:), allocatable :: rwork
    integer, dimension(:), allocatable :: iwork

    real(dp), parameter :: tol = 1.0d-12

    integer, parameter :: isym = 0, itol = 2, itmax = 101
    integer :: ierr, mxnelt

    if (first) then
      call ds2y(model%pcgdwk%pcgsize(1),model%pcgdwk%pcgsize(2), &
                model%pcgdwk%pcgrow,model%pcgdwk%pcgcol, &
                model%pcgdwk%pcgval,isym)
    end if

    mxnelt = 20 * model%pcgdwk%pcgsize(1)

    allocate(rwork(mxnelt),iwork(mxnelt))

!**     solve the problem using the SLAP package routines     
!**     -------------------------------------------------
!**     n ... order of matrix a (in)
!**     b ... right hand side vector (in)                        
!**     x ... initial quess/final solution vector (in/out)                        
!**     nelt ... number of non-zeroes in A (in)
!**     ia, ja ... sparse matrix format of A (in)
!**     a ... matrix helt in SLAT column format (in)
!**     isym ... storage method (0 is complete) (in)
!**     itol ... convergence criteria (2 recommended) (in)                     
!**     tol ... criteria for convergence (in)
!**     itmax ... maximum number of iterations (in)
!**     iter ... returned number of iterations (out)
!**     err ... error estimate of solution (out)
!**     ierr ... returned error message (0 is ok) (out)
!**     iunit ... unit for error writes during iteration (0 no write) (in)
!**     rwork ... workspace for SLAP routines (in)
!**     mxnelt ... maximum array and vector sizes (in)
!**     iwork ... workspace for SLAP routines (in)

    call dslucs(model%pcgdwk%pcgsize(1),model%pcgdwk%rhsd,model%pcgdwk%answ,model%pcgdwk%pcgsize(2), &
                model%pcgdwk%pcgrow,model%pcgdwk%pcgcol,model%pcgdwk%pcgval, &
                isym,itol,tol,itmax,iter,err,ierr,0,rwork,mxnelt,iwork,mxnelt)

    if (ierr /= 0) then
      print *, 'pcg error ', ierr, itmax, iter
      call dcpplt(model%pcgdwk%pcgsize(1),model%pcgdwk%pcgsize(2),model%pcgdwk%pcgrow, &
                model%pcgdwk%pcgcol,model%pcgdwk%pcgval,isym,lunit)
      print *, model%pcgdwk%pcgval
      call glide_finalise(model,.true.)
      stop
    end if

    deallocate(rwork,iwork)

  end subroutine slapsolv 

!---------------------------------------------------------------------------------

  subroutine putpcgc(pcgdwk,value,col,row)

    use glimmer_global, only : dp
 ! use pcgdwk, only : pcgval, pcgcol, pcgrow, ct

    implicit none

    type(glide_pcgdwk) :: pcgdwk
    integer, intent(in) :: row, col
    real(dp), intent(in) :: value

    if (value /= 0.0d0) then
      pcgdwk%pcgval(pcgdwk%ct) = value
      pcgdwk%pcgcol(pcgdwk%ct) = col
      pcgdwk%pcgrow(pcgdwk%ct) = row
      pcgdwk%ct = pcgdwk%ct + 1
    end if

  end subroutine putpcgc
  
!---------------------------------------------------------------------------------

  subroutine geomders(numerics,ipvr,stagthck,opvrew,opvrns)

    use glimmer_global, only : dp

    implicit none 

    type(glide_numerics) :: numerics
    real(dp), intent(out), dimension(:,:) :: opvrew, opvrns
    real(dp), intent(in), dimension(:,:) :: ipvr, stagthck

    real(dp) :: dew2, dns2 
    integer :: ew,ns,ewn,nsn

    ! Obviously we don't need to do this every time,
    ! but will do so for the moment.
    dew2 = 1.d0/(2.0d0 * numerics%dew)
    dns2 = 1.d0/(2.0d0 * numerics%dns)
    ewn=size(ipvr,1)
    nsn=size(ipvr,2)

    do ns=1,nsn-1
       do ew = 1,ewn-1
          if (stagthck(ew,ns) /= 0.0d0) then
             opvrew(ew,ns) = (ipvr(ew+1,ns+1)+ipvr(ew+1,ns)-ipvr(ew,ns)-ipvr(ew,ns+1)) * dew2
             opvrns(ew,ns) = (ipvr(ew+1,ns+1)+ipvr(ew,ns+1)-ipvr(ew,ns)-ipvr(ew+1,ns)) * dns2
          else
             opvrew(ew,ns) = 0.
             opvrns(ew,ns) = 0.
          end if
       end do
    end do

    !opvrew = (cshift(cshift(ipvr,1,2),1,1) + cshift(ipvr,1,1) - ipvr - cshift(ipvr,1,2)) * dew2
    !opvrns = (cshift(cshift(ipvr,1,2),1,1) + cshift(ipvr,1,2) - ipvr - cshift(ipvr,1,1)) * dns2
    
  end subroutine geomders

!---------------------------------------------------------------------------------

  subroutine geom2ders(model,ipvr,stagthck,opvrew,opvrns)

    use glimmer_global, only : dp ! ew, ewn, ns, nsn
!    use numerics, only : dew, dns

    implicit none 

    type(glide_global_type) :: model
    real(dp), intent(out), dimension(:,:) :: opvrew, opvrns
    real(dp), intent(in), dimension(:,:) :: ipvr, stagthck

    real(dp) :: dewsq4, dnssq4
    integer :: ew,ns

    integer :: pt(2)

    dewsq4 = 4.0d0 * model%numerics%dew * model%numerics%dew
    dnssq4 = 4.0d0 * model%numerics%dns * model%numerics%dns
 
    do ns = 2, model%general%nsn-2
      do ew = 2, model%general%ewn-2
        if (stagthck(ew,ns) .gt. 0.0d0) then
          opvrew(ew,ns) = centerew(ew,ns)
          opvrns(ew,ns) = centerns(ew,ns)
        else
          opvrew(ew,ns) = 0.0d0
          opvrns(ew,ns) = 0.0d0
        end if
      end do
    end do

! *** 2nd order boundaries using upwinding

    do ew = 1, model%general%ewn-1, model%general%ewn-2

      pt = whichway(ew)

      do ns = 2, model%general%nsn-2 
        if (stagthck(ew,ns) .gt. 0.0d0) then
          opvrew(ew,ns) = boundyew(pt,ew,ns)
          opvrns(ew,ns) = centerns(ew,ns)
        else
          opvrew(ew,ns) = 0.0d0
          opvrns(ew,ns) = 0.0d0
        end if
      end do

    end do

    do ns = 1, model%general%nsn-1, model%general%nsn-2

      pt = whichway(ns)

      do ew = 2, model%general%ewn-2  
        if (stagthck(ew,ns) .gt. 0.0d0) then
          opvrew(ew,ns) = centerew(ew,ns)
          opvrns(ew,ns) = boundyns(pt,ew,ns)
        else
          opvrew(ew,ns) = 0.0d0
          opvrns(ew,ns) = 0.0d0
        end if
      end do

    end do

    do ns = 1, model%general%nsn-1, model%general%nsn-2
      do ew = 1, model%general%ewn-1, model%general%ewn-2
        if (stagthck(ew,ns) .gt. 0.0d0) then
          pt = whichway(ew)
          opvrew(ew,ns) = boundyew(pt,ew,ns)
          pt = whichway(ns)
          opvrns(ew,ns) = boundyns(pt,ew,ns)
        else
          opvrew(ew,ns) = 0.0d0
          opvrns(ew,ns) = 0.0d0
        end if
      end do
    end do

  contains

    function centerew(ew,ns)

      implicit none

      real(dp) :: centerew
      integer ns,ew

      centerew = (sum(ipvr(ew+2,ns:ns+1)) + sum(ipvr(ew-1,ns:ns+1)) - &
                  sum(ipvr(ew+1,ns:ns+1)) - sum(ipvr(ew,ns:ns+1))) / dewsq4
    
    end function centerew 

    function centerns(ew,ns)

      implicit none

      real(dp) :: centerns
      integer ns,ew

      centerns = (sum(ipvr(ew:ew+1,ns+2)) + sum(ipvr(ew:ew+1,ns-1)) - &
                  sum(ipvr(ew:ew+1,ns+1)) - sum(ipvr(ew:ew+1,ns))) / dnssq4
  
    end function centerns 

    function boundyew(pt,ew,ns)

      implicit none

      integer, intent(in) :: pt(2)
      real(dp) :: boundyew
      integer ns,ew

      boundyew = pt(1) * (3.0d0 * sum(ipvr(pt(2),ns:ns+1)) - 7.0d0 * sum(ipvr(pt(2)+pt(1),ns:ns+1)) + &
                 5.0d0 * sum(ipvr(pt(2)+2*pt(1),ns:ns+1)) - sum(ipvr(pt(2)+3*pt(1),ns:ns+1))) / dewsq4

    end function boundyew

    function boundyns(pt,ew,ns)

      implicit none

      integer, intent(in) :: pt(2)
      real(dp) :: boundyns
      integer ns,ew

      boundyns = pt(1) * (3.0d0 * sum(ipvr(ew:ew+1,pt(2))) - 7.0d0 * sum(ipvr(ew:ew+1,pt(2)+pt(1))) + &
                 5.0d0 * sum(ipvr(ew:ew+1,pt(2)+2*pt(1))) - sum(ipvr(ew:ew+1,pt(2)+3*pt(1)))) / dnssq4

    end function boundyns

    function whichway(i)

      implicit none

      integer, intent(in) :: i
      integer :: whichway(2) 

      if (i == 1) then 
        whichway = (/1,1/)
      else
        whichway = (/-1,i+1/)
      end if

    end function whichway

  end subroutine geom2ders

!---------------------------------------------------------------------------------

  subroutine stagvarb(ipvr,opvr,ewn,nsn)

    use glimmer_global, only : dp ! ewn, nsn
 
    implicit none 

    real(dp), intent(out), dimension(:,:) :: opvr 
    real(dp), intent(in), dimension(:,:) :: ipvr
    integer :: ewn,nsn

    opvr(1:ewn-1,1:nsn-1) = (ipvr(2:ewn,1:nsn-1) + ipvr(1:ewn-1,2:nsn) + &
                             ipvr(2:ewn,2:nsn)   + ipvr(1:ewn-1,1:nsn-1)) / 4.0d0

  end subroutine stagvarb

!---------------------------------------------------------------------------------

  subroutine timeders(thckwk,ipvr,opvr,mask,time,which)

    use glimmer_global, only : dp, sp
    use paramets, only : conv

    implicit none 

    type(glide_thckwk) :: thckwk
    real(dp), intent(out), dimension(:,:) :: opvr 
    real(dp), intent(in), dimension(:,:) :: ipvr
    real(sp), intent(in) :: time 
    integer, intent(in), dimension(:,:) :: mask
    integer, intent(in) :: which

    real(sp) :: factor

    factor = (time - thckwk%oldtime)
    if (factor .eq.0) then
       opvr = 0.0d0
    else
       factor = 1./factor
       where (mask /= 0)
          opvr = conv * (ipvr - thckwk%olds(:,:,which)) * factor
       elsewhere
          opvr = 0.0d0
       end where
    end if

    thckwk%olds(:,:,which) = ipvr

    if (which == thckwk%nwhich) then
      thckwk%oldtime = time
    end if

  end subroutine timeders

  subroutine filterthck(thck,ewn,nsn)

    use glimmer_global, only : dp ! ew, ewn, ns, nsn

    implicit none

    real(dp), dimension(:,:), intent(inout) :: thck
    real(dp), dimension(:,:), allocatable :: smth
    integer :: ewn,nsn

    real(dp), parameter :: f = 0.1d0 / 16.0d0
    integer :: count
    integer :: ns,ew

    allocate(smth(ewn,nsn))
    count = 1

    do ns = 3,nsn-2
      do ew = 3,ewn-2

        if (all((thck(ew-2:ew+2,ns) > 0.0d0)) .and. all((thck(ew,ns-2:ns+2) > 0.0d0))) then
          smth(ew,ns) =  thck(ew,ns) + f * &
                        (thck(ew-2,ns) - 4.0d0 * thck(ew-1,ns) + 12.0d0 * thck(ew,ns) - &
                         4.0d0 * thck(ew+1,ns) + thck(ew+2,ns) + &
                         thck(ew,ns-2) - 4.0d0 * thck(ew,ns-1) - &
                         4.0d0 * thck(ew,ns+1) + thck(ew,ns+2))
          count = count + 1
        else
          smth(ew,ns) = thck(ew,ns)
        end if

      end do
    end do

    thck(3:ewn-2,3:nsn-2) = smth(3:ewn-2,3:nsn-2)
    print *, count

    deallocate(smth)            

  end subroutine filterthck

!----------------------------------------------------------------------

  subroutine swapbndh(bc,a,b,c,d)

    use glimmer_global, only : dp

    implicit none

    real(dp), intent(out), dimension(:) :: a, c
    real(dp), intent(in), dimension(:) :: b, d
    integer, intent(in) :: bc

    if (bc == 0) then
      a = b
      c = d
    end if

  end subroutine swapbndh

!------------------------------------------------------------------------

  subroutine stagleapthck(model,uflx,vflx,thck)

    use glimmer_global, only : dp

    implicit none
 
    type(glide_global_type) :: model
    real(dp), intent(in), dimension(:,:) :: uflx, vflx
    real(dp), intent(inout), dimension(:,:) :: thck

    real(dp), dimension(:,:), allocatable :: newthck, bnduflx, bndvflx 
    integer :: ns,ew

    allocate(newthck(model%general%ewn,model%general%nsn))
    newthck = 0.0d0

    allocate(bnduflx(0:model%general%ewn,0:model%general%nsn))
    allocate(bndvflx(0:model%general%ewn,0:model%general%nsn))

    bnduflx(1:model%general%ewn-1,1:model%general%nsn-1) = uflx
    bndvflx(1:model%general%ewn-1,1:model%general%nsn-1) = vflx

    bnduflx(0:model%general%ewn:model%general%ewn,:) = 0.0d0
    bnduflx(:,0:model%general%nsn:model%general%nsn) = 0.0d0
    bndvflx(0:model%general%ewn:model%general%ewn,:) = 0.0d0
    bndvflx(:,0:model%general%nsn:model%general%nsn) = 0.0d0

    ! bnduflx(0,:) = (/0.0d0,uflx(1,:),0.0d0/)
    ! bndvflx(0,:) = (/0.0d0,vflx(1,:),0.0d0/) 
    ! bnduflx(model%general%ewn,:) = (/0.0d0,uflx(model%general%ewn-1,:),0.0d0/)
    ! bndvflx(model%general%ewn,:) = (/0.0d0,vflx(model%general%ewn-1,:),0.0d0/)
  
    bnduflx(:,0) = (/0.0d0,uflx(:,1),0.0d0/)
    bndvflx(:,0) = (/0.0d0,vflx(:,1),0.0d0/);
    ! bnduflx(:,model%general%nsn) = (/0.0d0,uflx(:,model%general%nsn-1),0.0d0/)
    ! bnduflx(:,model%general%nsn) = (/0.0d0,vflx(:,model%general%nsn-1),0.0d0/)

    do ns = 1,model%general%nsn
      do ew = 1,model%general%ewn

        newthck(ew,ns) = - &
        (sum(bnduflx(ew,ns-1:ns)) - sum(bnduflx(ew-1,ns-1:ns))) * model%thckwk%few - &
        (sum(bndvflx(ew-1:ew,ns)) - sum(bndvflx(ew-1:ew,ns-1))) * model%thckwk%fns 

      end do
    end do

    newthck = model%thckwk%oldthck  - model%thckwk%basestate  + newthck  

    newthck(1,:) = model%thckwk%oldthck(1,:)
    newthck(:,1) = model%thckwk%oldthck(:,1)
    newthck(model%general%ewn,:) = model%thckwk%oldthck(model%general%ewn,:)
    newthck(:,model%general%nsn) = model%thckwk%oldthck(:,model%general%nsn)

    model%thckwk%oldthck = thck

    if (model%thckwk%first1) then
      model%thckwk%basestate = newthck - model%thckwk%oldthck  
      newthck = newthck - model%thckwk%basestate  
      model%thckwk%first1 = .false.
    end if

    where (newthck .gt. 0.0d0) 
      thck = newthck
    elsewhere
      thck = 0.0d0
    end where

    deallocate(newthck,bnduflx,bndvflx)

  end subroutine stagleapthck

end module glide_thck

