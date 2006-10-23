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

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

module glide_thck

  use glide_types

  private

!lipscomb - added thck_remap_evolve
  public :: init_thck, thck_nonlin_evolve, thck_lin_evolve, thck_remap_evolve, &
            stagvarb, geomders, timeders, stagleapthck

#ifdef DEBUG_PICARD
  ! debugging Picard iteration
  integer, private, parameter :: picard_unit=101
  real, private, parameter    :: picard_interval=500.
  integer, private            :: picard_max=0
#endif

contains

  subroutine init_thck(model)
    !*FD initialise work data for ice thickness evolution
    use glimmer_log
    implicit none
    type(glide_global_type) :: model

    
    model%pcgdwk%fc2 = (/ model%numerics%alpha * model%numerics%dt / (2.0d0 * model%numerics%dew * model%numerics%dew), &
         model%numerics%dt, (1.0d0-model%numerics%alpha) / model%numerics%alpha, &
         1.0d0 / model%numerics%alpha, model%numerics%alpha * model%numerics%dt / &
         (2.0d0 * model%numerics%dns * model%numerics%dns), 0.0d0 /) 

#ifdef DEBUG_PICARD
    call write_log('Logging Picard iterations')
    open(picard_unit,name='picard_info.data',status='unknown')
    write(picard_unit,*) '#time    max_iter'
#endif

    ! allocate memory for ADI scheme
    if (model%options%whichevol.eq.1) then
       allocate(model%thckwk%alpha(max(model%general%ewn, model%general%nsn)))
       allocate(model%thckwk%beta (max(model%general%ewn, model%general%nsn)))
       allocate(model%thckwk%gamma(max(model%general%ewn, model%general%nsn)))
       allocate(model%thckwk%delta(max(model%general%ewn, model%general%nsn)))
    end if
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
#ifdef DEBUG       
       print *, "* thck empty - net accumulation added", model%numerics%time
#endif
    else

       ! calculate basal velos
       if (newtemps) then
          call slipvelo(model,                &
               1,                             &
               model%velocity% btrc,          &
               model%velocity% ubas,          &
               model%velocity% vbas)
          ! calculate Glen's A if necessary
          call velo_integrate_flwa(model%velowk,model%geomderv%stagthck,model%temper%flwa)
       end if
       call slipvelo(model,                &
            2,                             &
            model%velocity% btrc,          &
            model%velocity% ubas,          &
            model%velocity% vbas)

       ! calculate diffusivity
       call velo_calc_diffu(model%velowk,model%geomderv%stagthck,model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns,model%velocity%diffu)

       ! get new thicknesses
       call thck_evolve(model,logunit,.true.,model%geometry%thck,model%geometry%thck)

       ! calculate horizontal velocity field
       call slipvelo(model,                &
            3,                             &
            model%velocity%btrc,           &
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
    use glide_setup
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
#ifdef DEBUG
       print *, "* thck empty - net accumulation added", model%numerics%time
#endif
    else

       ! calculate basal velos
       if (newtemps) then
          call slipvelo(model,                &
               1,                             &
               model%velocity% btrc,          &
               model%velocity% ubas,          &
               model%velocity% vbas)
          ! calculate Glen's A if necessary
          call velo_integrate_flwa(model%velowk,model%geomderv%stagthck,model%temper%flwa)
       end if

       first_p = .true.
       model%thckwk%oldthck = model%geometry%thck
       ! do Picard iteration
       do p=1,pmax
          model%thckwk%oldthck2 = model%geometry%thck

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

          call slipvelo(model,                &
               2,                             &
               model%velocity% btrc,          &
               model%velocity% ubas,          &
               model%velocity% vbas)

          ! calculate diffusivity
          call velo_calc_diffu(model%velowk,model%geomderv%stagthck,model%geomderv%dusrfdew, &
               model%geomderv%dusrfdns,model%velocity%diffu)

          ! get new thicknesses
          call thck_evolve(model,logunit,first_p,model%thckwk%oldthck,model%geometry%thck)

          first_p = .false.
          residual = maxval(abs(model%geometry%thck-model%thckwk%oldthck2))
          if (residual.le.tol) then
             exit
          end if
          
       end do
#ifdef DEBUG_PICARD
       picard_max=max(picard_max,p)
       if (model%numerics%tinc > mod(model%numerics%time,picard_interval)) then
          write(picard_unit,*) model%numerics%time,p
          picard_max = 0
       end if
#endif

       ! calculate horizontal velocity field
       call slipvelo(model,                &
            3,                             &
            model%velocity%btrc,           &
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
    use glide_setup, only: glide_calclsrf
    use glimmer_global, only : dp
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

    integer :: linit
    integer, parameter :: mxtbnd = 10, ewbc = 1, nsbc = 1

    ! ewbc/nsbc set the type of boundary condition aplied at the end of
    ! the domain. a value of 0 implies zero gradient.

    integer :: ew,ns

       !* the number of grid points and the number of nonzero 
       !* matrix elements (including bounary points)
       model%pcgdwk%pcgsize = (/ model%geometry%totpts, model%geometry%totpts * 5 /)

       model%pcgdwk%ct = 1

       ! Boundary Conditions
       ! lower and upper BC
       do ew = 1,model%general%ewn
          ns=1
          if (model%geometry%mask(ew,ns) /= 0) then
             call putpcgc(model%pcgdwk,1.0d0, model%geometry%mask(ew,ns), model%geometry%mask(ew,ns))
             if (calc_rhs) then
                model%pcgdwk%rhsd(model%geometry%mask(ew,ns)) = old_thck(ew,ns) 
             end if
             model%pcgdwk%answ(model%geometry%mask(ew,ns)) = new_thck(ew,ns)
          end if
          ns=model%general%nsn
          if (model%geometry%mask(ew,ns) /= 0) then
             call putpcgc(model%pcgdwk,1.0d0, model%geometry%mask(ew,ns), model%geometry%mask(ew,ns))
             if (calc_rhs) then
                model%pcgdwk%rhsd(model%geometry%mask(ew,ns)) = old_thck(ew,ns) 
             end if
             model%pcgdwk%answ(model%geometry%mask(ew,ns)) = new_thck(ew,ns)
          end if
       end do          
       !left and right BC
       if (model%options%periodic_ew.eq.1) then
          do ns=2,model%general%nsn-1
             ew = 1
             if (model%geometry%mask(ew,ns) /= 0) then
                call findsums(model%general%ewn-2,model%general%ewn-1,ns-1,ns)
                call generate_row(model%general%ewn-2,ew,ew+1,ns-1,ns,ns+1)
             end if
             ew=model%general%ewn
             if (model%geometry%mask(ew,ns) /= 0) then
                call findsums(1,2,ns-1,ns)
                call generate_row(ew-1,ew,3,ns-1,ns,ns+1)
             end if
          end do
       else
          do ns=2,model%general%nsn-1
             ew=1
             if (model%geometry%mask(ew,ns) /= 0) then
                call putpcgc(model%pcgdwk,1.0d0, model%geometry%mask(ew,ns), model%geometry%mask(ew,ns))
                if (calc_rhs) then
                   model%pcgdwk%rhsd(model%geometry%mask(ew,ns)) = old_thck(ew,ns) 
                end if
                model%pcgdwk%answ(model%geometry%mask(ew,ns)) = new_thck(ew,ns)
             end if
             ew=model%general%ewn
             if (model%geometry%mask(ew,ns) /= 0) then
                call putpcgc(model%pcgdwk,1.0d0, model%geometry%mask(ew,ns), model%geometry%mask(ew,ns))
                if (calc_rhs) then
                   model%pcgdwk%rhsd(model%geometry%mask(ew,ns)) = old_thck(ew,ns) 
                end if
                model%pcgdwk%answ(model%geometry%mask(ew,ns)) = new_thck(ew,ns)
             end if
          end do
       end if
       ! ice body
       do ns = 2,model%general%nsn-1
          do ew = 2,model%general%ewn-1

             if (model%geometry%mask(ew,ns) /= 0) then
                
                call findsums(ew-1,ew,ns-1,ns)
                call generate_row(ew-1,ew,ew+1,ns-1,ns,ns+1)

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

       model%pcgdwk%tlinit = model%pcgdwk%tlinit + linit
       model%pcgdwk%mlinit = max(linit,model%pcgdwk%mlinit)

       new_thck = max(0.0d0, new_thck)
#ifdef DEBUG
       print *, "* thck ", model%numerics%time, linit, model%pcgdwk%mlinit, model%pcgdwk%tlinit, model%geometry%totpts, &
            real(thk0*new_thck(model%general%ewn/2+1,model%general%nsn/2+1)), &
            real(vel0*maxval(abs(model%velocity%ubas))), real(vel0*maxval(abs(model%velocity%vbas))) 
#endif
       
       !------------------------------------------------------------
       ! calculate upper and lower surface
       !------------------------------------------------------------
       call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
       model%geometry%usrf = max(0.d0,model%geometry%thck + model%geometry%lsrf)

  contains
    subroutine generate_row(ewm,ew,ewp,nsm,ns,nsp)
      ! calculate row of sparse matrix equation
      implicit none
      integer, intent(in) :: ewm,ew,ewp  ! ew index to left, central, right node
      integer, intent(in) :: nsm,ns,nsp  ! ns index to lower, central, upper node

      ! fill sparse matrix
      call putpcgc(model%pcgdwk,sumd(1), model%geometry%mask(ewm,ns), model%geometry%mask(ew,ns))       ! point (ew-1,ns)
      call putpcgc(model%pcgdwk,sumd(2), model%geometry%mask(ewp,ns), model%geometry%mask(ew,ns))       ! point (ew+1,ns)
      call putpcgc(model%pcgdwk,sumd(3), model%geometry%mask(ew,nsm), model%geometry%mask(ew,ns))       ! point (ew,ns-1)
      call putpcgc(model%pcgdwk,sumd(4), model%geometry%mask(ew,nsp), model%geometry%mask(ew,ns))       ! point (ew,ns+1)
      call putpcgc(model%pcgdwk,1.0d0 + sumd(5), model%geometry%mask(ew,ns), model%geometry%mask(ew,ns))! point (ew,ns)

      ! calculate RHS
      if (calc_rhs) then
         model%pcgdwk%rhsd(model%geometry%mask(ew,ns)) = old_thck(ew,ns) * &
              (1.0d0 - model%pcgdwk%fc2(3) * sumd(5)) - model%pcgdwk%fc2(3) * &
              (old_thck(ewm,ns) * sumd(1) + old_thck(ewp,ns) * sumd(2) + &
              old_thck(ew,nsm) * sumd(3) + old_thck(ew,nsp) * sumd(4)) - &
              model%pcgdwk%fc2(4) * (model%geometry%lsrf(ew,ns) * sumd(5) + &
              model%geometry%lsrf(ewm,ns) * sumd(1) + model%geometry%lsrf(ewp,ns) * sumd(2) + &
              model%geometry%lsrf(ew,nsm) * sumd(3) + model%geometry%lsrf(ew,nsp) * sumd(4)) +  &
              model%climate%acab(ew,ns) * model%pcgdwk%fc2(2)
      end if
      model%pcgdwk%answ(model%geometry%mask(ew,ns)) = new_thck(ew,ns)      
    end subroutine generate_row

    subroutine findsums(ewm,ew,nsm,ns)
      ! calculate diffusivities
      implicit none
      integer, intent(in) :: ewm,ew  ! ew index to left, right
      integer, intent(in) :: nsm,ns  ! ns index to lower, upper

      ! calculate sparse matrix elements
      sumd(1) = model%pcgdwk%fc2(1) * (&
           (model%velocity%diffu(ewm,nsm) + model%velocity%diffu(ewm,ns)) + &
           (model%velocity%ubas (ewm,nsm) + model%velocity%ubas (ewm,ns)))
      sumd(2) = model%pcgdwk%fc2(1) * (&
           (model%velocity%diffu(ew,nsm) + model%velocity%diffu(ew,ns)) + &
           (model%velocity%ubas (ew,nsm) + model%velocity%ubas (ew,ns)))
      sumd(3) = model%pcgdwk%fc2(5) * (&
           (model%velocity%diffu(ewm,nsm) + model%velocity%diffu(ew,nsm)) + &
           (model%velocity%ubas (ewm,nsm) + model%velocity%ubas (ew,nsm)))
      sumd(4) = model%pcgdwk%fc2(5) * (&
           (model%velocity%diffu(ewm,ns) + model%velocity%diffu(ew,ns)) + &
           (model%velocity%ubas (ewm,ns) + model%velocity%ubas (ew,ns)))
      sumd(5) = - (sumd(1) + sumd(2) + sumd(3) + sumd(4))
    end subroutine findsums
  end subroutine thck_evolve

!---------------------------------------------------------------------------------

  subroutine slapsolv(model,first,iter,err,lunit)

    use glimmer_global, only : dp 
    use glide_stop
    use glimmer_log
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
      write(*,*) model%numerics%time
      call glide_finalise(model,.true.)
      call close_log
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
          opvrew(ew,ns) = boundyew(pt,ns)
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
          opvrns(ew,ns) = boundyns(pt,ew)
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
          opvrew(ew,ns) = boundyew(pt,ns)
          pt = whichway(ns)
          opvrns(ew,ns) = boundyns(pt,ew)
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

    function boundyew(pt,ns)

      implicit none

      integer, intent(in) :: pt(2)
      real(dp) :: boundyew
      integer ns

      boundyew = pt(1) * (3.0d0 * sum(ipvr(pt(2),ns:ns+1)) - 7.0d0 * sum(ipvr(pt(2)+pt(1),ns:ns+1)) + &
                 5.0d0 * sum(ipvr(pt(2)+2*pt(1),ns:ns+1)) - sum(ipvr(pt(2)+3*pt(1),ns:ns+1))) / dewsq4

    end function boundyew

    function boundyns(pt,ew)

      implicit none

      integer, intent(in) :: pt(2)
      real(dp) :: boundyns
      integer ew

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

  !-----------------------------------------------------------------------------
  ! ADI routines
  !-----------------------------------------------------------------------------

  subroutine stagleapthck(model,newtemps,logunit)
    
    !*FD this subroutine solves the ice sheet thickness equation using the ADI scheme
    !*FD diffusivities are updated for each half time step

    use glide_setup, only: glide_calclsrf
    use glide_velo
    use glimmer_utils
    implicit none
    ! subroutine arguments
    type(glide_global_type) :: model
    logical, intent(in) :: newtemps                     !*FD true when we should recalculate Glen's A
    integer,intent(in) :: logunit                       !*FD unit for logging

    ! local variables
    integer ew,ns, n

    if (model%geometry%empty) then

       model%geometry%thck = dmax1(0.0d0,model%geometry%thck + model%climate%acab * model%pcgdwk%fc2(2))
#ifdef DEBUG       
       print *, "* thck empty - net accumulation added", model%numerics%time
#endif
    else

       ! calculate basal velos
       if (newtemps) then
          call slipvelo(model,                &
               1,                             &
               model%velocity% btrc,          &
               model%velocity% ubas,          &
               model%velocity% vbas)
          ! calculate Glen's A if necessary
          call velo_integrate_flwa(model%velowk,model%geomderv%stagthck,model%temper%flwa)
       end if
       call slipvelo(model,                &
            2,                             &
            model%velocity% btrc,          &
            model%velocity% ubas,          &
            model%velocity% vbas)

       ! calculate diffusivity
       call velo_calc_diffu(model%velowk,model%geomderv%stagthck,model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns,model%velocity%diffu)

       model%velocity%total_diffu(:,:) = model%velocity%diffu(:,:) + model%velocity%ubas(:,:)

       ! first ADI step, solve thickness equation along rows j
       n = model%general%ewn
       do ns=2,model%general%nsn-1
          call adi_tri ( model%thckwk%alpha, model%thckwk%beta, model%thckwk%gamma, model%thckwk%delta, &
               model%geometry%thck(:,ns), model%geometry%lsrf(:,ns), model%climate%acab(:,ns), &
               model%velocity%vflx(:,ns), model%velocity%vflx(:,ns-1), &
               model%velocity%total_diffu(:,ns),  model%velocity%total_diffu(:,ns-1), &
               model%numerics%dt, model%numerics%dew, model%numerics%dns )

          call tridag( model%thckwk%alpha(2:n), model%thckwk%beta(1:n), model%thckwk%gamma(1:n-1), model%thckwk%delta(1:n), &
               model%thckwk%oldthck(:,ns), n )
       end do

       model%thckwk%oldthck(:,:) = max(model%thckwk%oldthck(:,:), 0.d0)

       ! second ADI step, solve thickness equation along columns i
       n = model%general%nsn
       do ew=2,model%general%ewn-1
          call adi_tri ( model%thckwk%alpha, model%thckwk%beta, model%thckwk%gamma, model%thckwk%delta, &
               model%thckwk%oldthck(ew,:), model%geometry%lsrf(ew, :), model%climate%acab(ew, :), &
               model%velocity%uflx(ew,:), model%velocity%uflx(ew-1,:), &
               model%velocity%total_diffu(ew,:), model%velocity%total_diffu(ew-1,:), &
               model%numerics%dt, model%numerics%dns, model%numerics%dew )

          call tridag( model%thckwk%alpha(2:n), model%thckwk%beta(1:n), model%thckwk%gamma(1:n-1), model%thckwk%delta(1:n), &
               model%geometry%thck(ew, :), n)
       end do

       model%geometry%thck(:,:) = max(model%geometry%thck(:,:), 0.d0)

       ! calculate horizontal velocity field
       call slipvelo(model,                &
            3,                             &
            model%velocity%btrc,           &
            model%velocity%ubas,           &
            model%velocity%vbas)
       call velo_calc_velo(model%velowk,model%geomderv%stagthck,model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns,model%temper%flwa,model%velocity%diffu,model%velocity%ubas, &
            model%velocity%vbas,model%velocity%uvel,model%velocity%vvel,model%velocity%uflx,model%velocity%vflx)
    end if

    !------------------------------------------------------------
    ! calculate upper and lower surface
    !------------------------------------------------------------
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0,model%geometry%thck + model%geometry%lsrf)

  end subroutine stagleapthck


  subroutine adi_tri(a,b,c,d,thk,tpg,mb,flx_p,flx_m,dif_p,dif_m,dt,ds1, ds2)
    !*FD construct tri-diagonal matrix system for a column/row
    use glimmer_global, only : dp, sp
    implicit none
    
    real(dp), dimension(:), intent(out) :: a !*FD alpha (subdiagonal)
    real(dp), dimension(:), intent(out) :: b !*FD alpha (diagonal)
    real(dp), dimension(:), intent(out) :: c !*FD alpha (superdiagonal)
    real(dp), dimension(:), intent(out) :: d !*FD right-hand side
    
    real(dp), dimension(:), intent(in) :: thk   !*FD ice thickness
    real(dp), dimension(:), intent(in) :: tpg   !*FD lower surface of ice
    real(sp), dimension(:), intent(in) :: mb    !*FD mass balance
    real(dp), dimension(:), intent(in) :: flx_p !*FD flux +1/2
    real(dp), dimension(:), intent(in) :: flx_m !*FD flux -1/2
    real(dp), dimension(:), intent(in) :: dif_p !*FD diffusivity +1/2
    real(dp), dimension(:), intent(in) :: dif_m !*FD diffusivity -1/2
    
    real(dp), intent(in) :: dt !*FD time step
    real(dp), intent(in) :: ds1, ds2 !*FD spatial steps inline and transversal

    ! local variables
    real(dp) :: f1, f2, f3
    integer :: i,n
    
    n = size(thk)

    f1 = dt/(4*ds1*ds1)
    f2 = dt/(4*ds2)
    f3 = dt/2.

    a(:) = 0.
    b(:) = 0.
    c(:) = 0.
    d(:) = 0.

    a(1) = 0.
    do i=2,n
       a(i) = f1*(dif_m(i-1)+dif_p(i-1))
    end do
    do i=1,n-1
       c(i) = f1*(dif_m(i)+dif_p(i))
    end do
    c(n) = 0.
    b(:) = -(a(:)+c(:))

    ! calculate RHS
    do i=2,n-1
       d(i) = thk(i) - &
            f2 * (flx_p(i-1) + flx_p(i) - flx_m(i-1) - flx_m(i)) + &
            f3 * mb(i) - &
            a(i)*tpg(i-1) - b(i)*tpg(i) - c(i)*tpg(i+1)
    end do

    b(:) = 1.+b(:)

  end subroutine adi_tri

!---------------------------------------------------------------------
!lipscomb - new remapping subroutine, patterned after thck_lin_evolve
 
  subroutine thck_remap_evolve(model,newtemps,logunit,l_advtemp)
 
    ! Solve the continuity equation for ice thickness and temperature, as well
    !  as any tracers, using an incremental remapping transport scheme.
    ! This subroutine first computes the horizontal velocities, then calls the
    !  remapping driver.
 
    use glimmer_global, only : dp
    use physcon, only: scyr
    use paramets, only: thk0, tim0
    use glide_velo
    use glide_remap
 
    implicit none
 
    ! in-out arguments
 
    type(glide_global_type), intent(inout) :: model
    integer,intent(in) :: logunit      ! unit for logging
 
    logical, intent(in) ::   &
         newtemps,    &! true when recalculating Glen's A
         l_advtemp     ! if true, advect temperature and tracers
                       ! if false, advect thickness only
 
    ! local variables
 
    integer, parameter ::   &
          ntrace = 2         ! no. of tracer fields
                             ! 1 = temperature, 2 = ice age
 
    logical, parameter ::  &
        predict_correct = .false.      ! if true, use predicted thickness
                                       ! to correct the velocity 
    integer :: i, j, k, k1, k2, nt
 
    integer ::   &
         ewn,    &! model%general%ewn
         nsn,    &! model%general%nsn
         upn,    &! model%general%upn, number of vertical levels
         nlyr     ! number of ice layers, upn - 1
 
    real(dp), dimension(model%general%ewn, model%general%nsn) ::     &
         workh,      &! work array for thck
         worku,      &! work array for uvel
         workv,      &! work array for vvel
                      ! NOTE: These array bounds are different from uvel, vvel!
         e1, e2       ! internal energy before and after vertical remapping
 
    real(dp), dimension(model%general%ewn, model%general%nsn,   &
                        model%general%upn-1) ::     &
         workl        ! thickness of a single layer 
 
    real(dp), dimension(model%general%ewn, model%general%nsn,     &
                        ntrace,            0:model%general%upn) ::  &
         workt        ! tracer array (temperature, age)
                      ! Note: Array bounds are different from temper%temp
 
    real(dp), dimension(model%general%upn-1) :: &
         dsigma       ! model%velowk%dups, fractional thickness in layer k
      
    real(dp), dimension(model%general%ewn, model%general%nsn,   &
                        0:model%general%upn-1) ::     &
         zi1,        &! depths of layer interfaces after horizontal remapping
                      ! and accumulation/ablation
         zi2          ! depths of layer interfaces in sigma coordinates
 
    real (dp), dimension(model%general%ewn, model%general%nsn, ntrace) ::  &
         ht_init,    &! initial sum over ice sheet of h*tracer
         ht_final     ! final sum over ice sheet of h*tracer
 
    real(dp) ::   &
         dt,             &! time step (s)
         hlyr,           &! layer thickness
         hmlt,           &! thickness of melted ice
         dz,             &! difference of two depths
         grad             ! tracer gradient
 
    logical, parameter ::   &
         monotonicity_check = .true.       ,&! monotonicity check for remapping
         horiz_conservation_check = .true. ,&! bug check for horizontal remapping
         vert_conservation_check = .true.    ! bug check for vertical remapping
 
    !------------------------------------------------------------------
    ! Initialise
    !------------------------------------------------------------------
 
    ewn = model%general%ewn
    nsn = model%general%nsn
    upn = model%general%upn
    nlyr = upn - 1
    dsigma(:) = model%velowk%dups(:)
 
    worku(:,:) = c0
    workv(:,:) = c0
    workh(:,:) = c0
    workt(:,:,:,:) = c0

    if (model%geometry%empty) then
 
       ! surface ablation/accumulation
 
       model%geometry%thck(:,:) = dmax1(0.0d0,model%geometry%thck(:,:) &
                                 + model%climate%acab(:,:) * model%pcgdwk%fc2(2))
#ifdef DEBUG       
       print *, "* thck empty - net accumulation added", model%numerics%time
#endif
 
    else    ! not empty
 
    !------------------------------------------------------------------
    ! Calculate basal velocity
    !------------------------------------------------------------------
 
       !lipscomb - Set flag = 0 to compute basal velocity as linear function 
       !           of gravitational driving stress.  
       !lipscomb - Not sure why the flag is not simply set to
       !           model%glide_options%which_slip
 
       call slipvelo(model,                        &
                     0,                            &
                     model%velocity%btrc,          &
                     model%velocity%ubas,          &
                     model%velocity%vbas)
 
    !------------------------------------------------------------------
    ! Calculate horizontal velocity field
    !------------------------------------------------------------------
 
       !lipscomb - With flag = 0, this subroutine computes uvel and vvel,
       !           following Payne and Dongelmans eq. 8.  
 
       call zerovelo (model%velowk,             &
                      model%numerics%sigma,     &
                      0,                        &  ! Payne & Dongelmans
                      model%geomderv%stagthck,  &
                      model%geomderv%dusrfdew,  &
                      model%geomderv%dusrfdns,  &
                      model%temper%flwa,        &
                      model%velocity%ubas,      &
                      model%velocity%vbas,      &
                      model%velocity%uvel,      &  ! uvel, vvel are needed
                      model%velocity%vvel,      &
                      model%velocity%uflx,      &  ! uflx, vflx are not needed
                      model%velocity%vflx,      &  ! 
                      model%velocity%diffu)        ! diffu not needed
 
       if (predict_correct) then
 
    !------------------------------------------------------------------
    ! To reduce the time discretization error, there is the option
    ! of predicting the thickness field at time n + 1/2, then 
    ! correcting the velocity based on the mid-time thickness.
    !------------------------------------------------------------------
 
          workh(:,:) = model%geometry%thck(:,:)
          workt(:,:,1,1) = c1   ! dummy tracer field
 
          ! compute vertical average velocity
          ! Note: uvel and vvel have dimensions(ewn-1, nsn-1, upn)
          !       worku and workv have dimensions (ewn, nsn)
          !       worku and workv are initialized to zero above
 
          do k = 1, nlyr
             do j = 1, nsn-1
             do i = 1, ewn-1
                worku(i,j) = worku(i,j)   &
                           + model%velocity%uvel(k,i,j) * dsigma(k)
                workv(i,j) = workv(i,j)   &
                           + model%velocity%vvel(k,i,j) * dsigma(k)
             enddo
             enddo
          enddo
 
          ! predict thickness at time n + 1
 
          call remap_driver (ewn,                    &
                             nsn,                    &
                             1,                      &  ! no. of tracer fields
                             1,                      &  ! no. of ghost cells
                             model%gridwk,           &  ! grid quantities
                             model%numerics%dt,      &  ! time step (s)
                             worku,                  &  ! uvel
                             workv,                  &  ! vvel
                             workh,                  &  ! thck
                             workt(:,:,1,1),         &  ! dummy tracer
                             horiz_conservation_check, &
                             monotonicity_check)
 
          ! Given the predicted thickness at time n + 1, estimate the thickness
          ! at time n + 1/2
 
          workh(:,:) = p5 * (model%geometry%thck(:,:) + workh(:,:))
 
          ! Correct H and grad(S) at grid cell corners
 
          call stagvarb(workh,                         &
                        model%geomderv%stagthck,       &
                        ewn,                           &
                        nsn) 
 
          call geomders(model%numerics,                &
                        model%geometry%usrf,           &
                        model%geomderv%stagthck,       & 
                        model%geomderv%dusrfdew,       &
                        model%geomderv%dusrfdns)
 
          ! Correct the velocities
 
          call zerovelo (model%velowk,             &
                         model%numerics%sigma,     &
                         0,                        &
                         model%geomderv%stagthck,  &
                         model%geomderv%dusrfdew,  &
                         model%geomderv%dusrfdns,  &
                         model%temper%flwa,        &
                         model%velocity%ubas,      &
                         model%velocity%vbas,      &
                         model%velocity%uvel,      &
                         model%velocity%vvel,      &
                         model%velocity%uflx,      & ! uflx, vflx, and diffu
                         model%velocity%vflx,      & ! not needed for remapping
                         model%velocity%diffu)        
 
       endif   ! predict_correct
 
    !------------------------------------------------------------------
    ! Calculate flow factor A
    !------------------------------------------------------------------
 
       if (newtemps) then
          call velo_integrate_flwa (model%velowk,             &
                                    model%geomderv%stagthck,  &
                                    model%temper%flwa)
       endif
 
       if (l_advtemp) then
 
       !---------------------------------------------------------------
       ! Advance the thickness, temperature, and tracers in each ice layer.
       !---------------------------------------------------------------
 
          ! Increment ice age 
          ! Note: If l_advtemp is false, ice age is not updated or advected.
 
          model%geometry%age(:,:,:) = model%geometry%age(:,:,:)  &
                                    + model%numerics%dt
 
          do k = 1, nlyr
 
          ! Fill work arrays of the desired dimensions.
          ! Note: Work array dimensions are different from glimmer dimensions.
          !  model%velocity%uvel and vvel are (ewn-1, nsn-1, upn)
          !  model%temper%temp is (upn, 0:ewn+1, 0:nsn+1)
          !  worku and workv are (ewn, nsn),
          !  workt is (ewn, nwn, 0:upn)
 
             do j = 1, nsn
             do i = 1, ewn
                workl(i,j,k)   = model%geometry%thck(i,j) * dsigma(k)
                workt(i,j,1,k) = model%temper%temp(k,i,j)
                workt(i,j,2,k) = model%geometry%age(k,i,j)
             enddo
             enddo
 
             do j = 1, nsn-1
             do i = 1, ewn-1
                worku(i,j) = model%velocity%uvel(k,i,j)
                workv(i,j) = model%velocity%vvel(k,i,j)
             enddo
             enddo
 
          ! remap the thickness and tracer fields
 
!lipscomb - Pretend there is a layer of ghost cells so remapping loops 
!           are correct.
 
             call remap_driver(ewn,                    &
                               nsn,                    &
                               ntrace,                 &! no. of tracer fields
                               1,                      &! no. of ghost cells
                               model%gridwk,           &! grid quantities
                               model%numerics%dt,      &! time step (s)
                               worku(:,:),             &! uvel
                               workv(:,:),             &! vvel
                               workl(:,:,k),           &! layer thickness
                               workt(:,:,:,k),         &! tracers
                               horiz_conservation_check, &
                               monotonicity_check)
          enddo   ! k
 
          !------------------------------------------------------------
          ! Surface ablation/accumulation
          ! Compute thickness, temperature and age changes at surface.
          !------------------------------------------------------------
 
          workh(:,:) = model%climate%acab(:,:) * model%numerics%dt
 
          do j = 1, nsn
          do i = 1, ewn
             if (workh(i,j) > c0) then   ! surface accumulation
                k = 1
                hlyr = workl(i,j,k) + workh(i,j)
 
             ! layer 1 temperature 
                workt(i,j,1,k) = (workl(i,j,k)*workt(i,j,1,k)  &
                                + workh(i,j)*model%climate%artm(i,j)) / hlyr
 
             ! layer 1 age (sfc accumulation has mean age dt/2)
                workt(i,j,2,k) = (workl(i,j,k)*workt(i,j,2,k)  &
                                + workh(i,j)*p5*model%numerics%dt) / hlyr
 
                workl(i,j,k) = hlyr
 
             elseif (workh(i,j) < c0) then  ! surface ablation
                hmlt = -workh(i,j)
                do k = 1, nlyr
                   if (hmlt > workl(i,j,k)) then
                      hmlt = hmlt - workl(i,j,k)
                      workl(i,j,k) = c0  
                      workt(i,j,1,k) = c0  ! temperature
                      workt(i,j,2,k) = c0  ! age
                   elseif (hmlt > c0) then
                      workl(i,j,k) = workl(i,j,k) - hmlt
                      hmlt = c0
                      exit   ! k loop
                   endif
                enddo  ! k
             ! Note: If hmlt > 0 after k loop, energy is not conserved
             endif     ! workh > 0.
          enddo        ! i
          enddo        ! j
 
          !------------------------------------------------------------
          ! In each column of ice, remap the tracers to standard sigma layers.
          !------------------------------------------------------------

          ! Update ice thickness
          workh(:,:) = c0
          do k = 1, nlyr
             workh(:,:) = workh(:,:) + workl(:,:,k)
          enddo
       
          if (vert_conservation_check) then
             ! compute sum of h*tracer for each column 
             ht_init(:,:,:) = c0
             do k = 1, nlyr
                do nt = 1, ntrace
                   do j = 1, nsn
                   do i = 1, ewn
                      ht_init(i,j,nt) = ht_init(i,j,nt)   &
                                        + workl(i,j,k)*workt(i,j,nt,k)
                   enddo
                   enddo
                enddo
             enddo
          endif   ! vert_conservation_check
 
          !------------------------------------------------------------
          ! Compute depths at each level, before and after vert remapping.
          !------------------------------------------------------------
 
          zi1(:,:,0) = c0
          zi2(:,:,0) = c0
          
          do k = 1, nlyr-1
             zi1(:,:,k) = zi1(:,:,k-1) + workl(:,:,k)          ! before
             zi2(:,:,k) = zi2(:,:,k-1) + workh(:,:)*dsigma(k)  ! after
          enddo
          
          zi1(:,:,nlyr) = workh(:,:)
          zi2(:,:,nlyr) = workh(:,:)
 
          !------------------------------------------------------------
          ! Set tracer values at surface and base.
          !------------------------------------------------------------
 
          do j = 1, nsn
          do i = 1, ewn
 
             ! temperature
             workt(i,j,1,0)      = model%climate%artm(i,j)     ! sfc air temp
             workt(i,j,1,nlyr+1) = model%temper%temp(upn,i,j)  ! basal temp
             
             ! age
             workt(i,j,2,0) = c0    ! top sfc age = 0.
             dz = p5 * (zi1(i,j,nlyr) - zi1(i,j,nlyr-2))
             if (dz > eps) then     ! extrapolate gradient to base
                grad = (workt(i,j,2,nlyr) - workt(i,j,2,nlyr-1)) / dz
                workt(i,j,2,nlyr+1) = workt(i,j,2,nlyr)  &
                                    + p5*workl(i,j,nlyr)*grad
             else                   ! use bottom-layer value at base
                workt(i,j,2,nlyr+1) = workt(i,j,2,nlyr)
             endif
          enddo     ! i
          enddo     ! j
 
          !------------------------------------------------------------
          ! Conservative, 2nd-order accurate vertical remapping
          ! Reduces to 1st order where needed to preserve monotonicity
          !------------------------------------------------------------
 
          call vertical_remap(ewn,      nsn,       &
                              ntrace,   nlyr,      &
                              zi1,      zi2,       &
                              workt)
 
          if (vert_conservation_check) then
             ! compute sum of h*tracer for each column
             ht_final(:,:,:) = c0
             do k = 1, nlyr
                do nt = 1, ntrace
                   do j = 1, nsn
                   do i = 1, ewn
                      ht_final(i,j,nt) = ht_final(i,j,nt) &
                                       + workh(i,j)*dsigma(k)*workt(i,j,nt,k)
                   enddo
                   enddo
                enddo
             enddo
 
             ! compare to initial value
             do nt = 1, ntrace
                do j = 1, nsn
                do i = 1, ewn
                   if (abs(ht_init(i,j,nt)) >= eps) then 
                      if (abs( (ht_final(i,j,nt)-ht_init(i,j,nt))  &
                              / ht_init(i,j,nt) ) > eps) then
                         print*, 'Vertical remapping: tracer conservation error'
                         print*, 'i, j, nt, ht_init, ht_final =',  &
                                  i, j, nt, ht_init(i,j,nt), ht_final(i,j,nt)
                         stop
                      endif
                   endif
                enddo
                enddo
             enddo
 
          endif   ! vert_conservation_check
 
          !---------------------------------------------------------------
          ! update thickness and tracer arrays
          !---------------------------------------------------------------
 
          model%geometry%thck(:,:) = workh(:,:)
          
          do k = 1, nlyr
             do j = 1, nsn
                do i = 1, ewn
                   model%temper%temp(k,i,j) = workt(i,j,1,k)
                   model%geometry%age(k,i,j) = workt(i,j,2,k)
                enddo
             enddo
          enddo
          
          !------------------------------------------------------------
          ! Apply periodic ew BC to temperature field 
          !------------------------------------------------------------
          if (model%options%periodic_ew.eq.1) then 
             model%temper%temp(:,0,:)     = model%temper%temp(:,ewn-2,:)
             model%temper%temp(:,1,:)     = model%temper%temp(:,ewn-1,:)
             model%temper%temp(:,ewn,:)   = model%temper%temp(:,2,:)
             model%temper%temp(:,ewn+1,:) = model%temper%temp(:,3,:)
          end if
 
       else   ! advect thickness only
 
          workt(:,:,1,1) = c1   ! dummy tracer field
 
          ! compute vertical average velocity
          ! Note: uvel and vvel have dimensions(ewn-1, nsn-1, upn)
          !       worku and workv have dimensions (ewn, nsn)
 
          worku(:,:) = c0
          workv(:,:) = c0
          do k = 1, nlyr
             do j = 1, nsn-1
             do i = 1, ewn-1
                worku(i,j) = worku(i,j)   &
                           + model%velocity%uvel(k,i,j) * dsigma(k)
                workv(i,j) = workv(i,j)   &
                           + model%velocity%vvel(k,i,j) * dsigma(k)
             enddo
             enddo
          enddo
 
          !------------------------------------------------------------
          ! Advance the thickness.
          !------------------------------------------------------------
 
          call remap_driver (ewn,                    &
                             nsn,                    &
                             1,                      &  ! no. of tracer fields
                             1,                      &  ! no. of ghost cells
                             model%gridwk,           &  ! grid quantities
                             model%numerics%dt,      &  ! time step (s)
                             worku,                  &  ! vert average uvel
                             workv,                  &  ! vert average vvel
                             model%geometry%thck,    &  ! total thickness
                             workt(:,:,1,1),         &  ! dummy tracer
                             horiz_conservation_check, &!
                             monotonicity_check)
 
          !------------------------------------------------------------
          ! surface ablation/accumulation
          !------------------------------------------------------------
 
          workh(:,:) = model%climate%acab(:,:) * model%numerics%dt
 
          model%geometry%thck(:,:) = max(c0, model%geometry%thck(:,:) &
                                                         + workh(:,:))
 
       endif     ! l_advtemp

    end if   ! empty
 
  end subroutine thck_remap_evolve
 
!****************************************************************************
  subroutine vertical_remap(ewn,      nsn,       &
                            ntrace,   nlyr,      &
                            zi1,      zi2,       &
                            trcr)
 
    ! Conservative remapping of tracer fields from one set of vertical 
    ! coordinates to another.  Tracer fields are reconstructed linearly
    ! in each layer.  The remapping is second-order accurate except
    ! where tracer gradients are limited to preserve monotonicity.
 
    use glimmer_global, only : dp
    use paramets
    use glide_remap, only: c0, c1, p5   !lipscomb - move elsewhere?
 
    implicit none
 
    ! in-out arguments
 
    integer, intent(in) ::  &
         ewn, nsn,   &! number of cells in EW and NS directions
         ntrace,     &! number of tracer fields
         nlyr         ! number of vertical layers
 
    real(dp), dimension (ewn, nsn, 0:nlyr), intent(in) ::  &
         zi1,        &! layer interfaces in old coordinate system
                      ! zi1(0) = 0. = value at top surface
                      ! zi1(k) = depth of lower interface of layer k
         zi2          ! layer interfaces in new coordinate system
                      ! Note: zi1(0) = zi2(0)
                      !       zi1(nlyr) = zi2(nlyr)
 
    real(dp), dimension (ewn, nsn, ntrace, 0:nlyr+1), intent(inout) ::   &
         trcr         ! tracer field to be remapped
                      ! trcr(0) = value at top surface
                      ! tracer(k) = value at midpoint of layer k
                      ! trcr(nlyr+1) = value at bottom surface
 
    ! local variables
 
    integer :: i, j, k, k1, k2, nt
 
    real(dp), dimension(ewn,nsn,ntrace,nlyr) ::       &
         gradt,      &! tracer gradient within a layer
         htsum        ! sum of thickness*tracer in a layer
         
    real(dp), dimension(ewn, nsn, 0:nlyr+1) ::        &
         zmid,       &! depths of layer midpoints for zi1
         lmask        ! = 1. if layer thickness > 0, else = 0.
 
    real(dp), dimension(ewn, nsn, nlyr) ::        &
         hlyr         ! layer thickness
 
    real(dp) ::   &
         dz,             &! difference of two depths
         zlo,            &! lower z value of overlap region
         zhi,            &! higher z value of overlap region
         zav,            &! average of two depths
         tup, tdn,       &! tracer value in layers above and below
         tmax, tmin,     &! max/min tracer value in a layer and its neighbors
         tzmax, tzmin,   &! max/min value of reconstructed tracer in layer
         wk1, wk2         ! work variables
 
       ! initialise
       gradt(:,:,:,:) = c0
       htsum(:,:,:,:) = c0
 
       !-----------------------------------------------------------------
       ! Compute midpoint and thickness of each layer.
       ! Tracer value is assumed to reside at zmid.
       !-----------------------------------------------------------------
 
       zmid(:,:,0) = c0
       do k = 1, nlyr
          zmid(:,:,k) = p5 * (zi1(:,:,k-1) + zi1(:,:,k))
          hlyr(:,:,k) = zi1(:,:,k) - zi1(:,:,k-1)
       enddo
       zmid(:,:,nlyr+1) = zi1(:,:,nlyr)
 
       !-----------------------------------------------------------------
       ! Compute mask.
       ! Layers with zero thickness are masked out so that the 
       !  tracer values in these layers do not affect the gradient. 
       ! Set lmask = 1 at top and bottom boundaries.
       !-----------------------------------------------------------------
 
       do k = 1, nlyr
          do j = 1, nsn
             do i = 1, ewn
                if (hlyr(i,j,k) > c0) then
                   lmask(i,j,k) = c1
                else
                   lmask(i,j,k) = c0
                endif
             enddo
          enddo
       enddo
 
       lmask(:,:,0) = c1
       lmask(:,:,nlyr+1) = c1
 
       !-----------------------------------------------------------------
       ! Compute tracer gradients within each layer
       ! Note: Skipping this loop gives a 1st-order accurate remapping.
       !-----------------------------------------------------------------
 
       do k = 1, nlyr
          do nt = 1, ntrace
             do j = 1, nsn
                do i = 1, ewn
                   if (hlyr(i,j,k) > c0) then
 
                ! Load values in neighbor cells
                ! For cells of zero thickness, substitute the value in cell k
 
                      tup =     lmask(i,j,k-1)  * trcr(i,j,nt,k-1)  &
                          + (c1-lmask(i,j,k-1)) * trcr(i,j,nt,k) 
                      tdn =     lmask(i,j,k+1)  * trcr(i,j,nt,k+1)  &
                          + (c1-lmask(i,j,k+1)) * trcr(i,j,nt,k) 
 
                      dz = zmid(i,j,k+1) - zmid(i,j,k-1)
                      gradt(i,j,nt,k) = (tdn - tup) / dz
             
                ! max and min deviations of tracer values in layer k-1, k, k+1
                ! from value in layer k
                      tmax = max(tup,trcr(i,j,nt,k),tdn) - trcr(i,j,nt,k)
                      tmin = min(tup,trcr(i,j,nt,k),tdn) - trcr(i,j,nt,k)
 
                ! max and min deviation of reconstructed values from mean value
                      wk1 = gradt(i,j,nt,k)*(zi1(i,j,k-1) - zmid(i,j,k))
                      wk2 = gradt(i,j,nt,k)*(zi1(i,j,k) - zmid(i,j,k))
                      tzmax = max(wk1, wk2)
                      tzmin = min(wk1, wk2)
 
                ! Limit tracer gradients
                
                      if (abs(tzmin) > c0) then
                         wk1 = max(c0, tmin/tzmin)
                      else
                         wk1 = c1
                      endif
                      
                      if (abs(tzmax) > c0) then
                         wk2 = max(c0, tmax/tzmax)
                      else
                         wk2 = c1
                      endif
 
                      gradt(i,j,nt,k) = gradt(i,j,nt,k) * min(c1, wk1, wk2)
                      
                   endif  ! hlyr > c0
                enddo     ! i
             enddo        ! j
          enddo           ! nt
       enddo              ! k
       
       ! new layer thicknesses
       do k = 1, nlyr
          hlyr(:,:,k) = zi2(:,:,k) - zi2(:,:,k-1)
       enddo
 
       !-----------------------------------------------------------------
       ! Compute sum of h*T for each new layer (k2) by integrating
       ! over the regions of overlap with old layers (k1).
       ! The basic formula is as follows:
       !
       ! int_zlo^zhi [T(z) dz] = int_zlo^zhi [trcr + gradt*(z - zmid)] dz
       !                       = dz * [trcr + gradt*(zav - zmid)] 
       ! where dz = zhi - zlo, zav = (zhi+zlo)/2
       !
       ! Note: It might be worth trying a more efficient
       !       search algorithm if the number of layers is large.
       !       This algorithm scales as nlyr^2.
       !-----------------------------------------------------------------
 
       do k2 = 1, nlyr
          do k1 = 1, nlyr
             do nt = 1, ntrace
                do j = 1, nsn
                   do i = 1, ewn
                      zhi = min (zi1(i,j,k1),   zi2(i,j,k2)) 
                      zlo = max (zi1(i,j,k1-1), zi2(i,j,k2-1))
                      dz = max(zhi - zlo, c0)
                      zav = p5 * (zlo + zhi)
                      htsum(i,j,nt,k2) = htsum(i,j,nt,k2)    &
                                       + dz * (trcr(i,j,nt,k1) & 
                                       + gradt(i,j,nt,k1)*(zav - zmid(i,j,k1)))
                   enddo   ! i
                enddo      ! j
             enddo         ! nt
          enddo            ! k1
       enddo               ! k2
 
       !-----------------------------------------------------------------
       ! Compute tracers in new layers (zi2)
       !-----------------------------------------------------------------
 
       trcr(:,:,:,:) = c0
       do k = 1, nlyr
          do nt = 1, ntrace
             do j = 1, nsn
                do i = 1, ewn
                   if (hlyr(i,j,k) > c0) then
                      trcr(i,j,nt,k) = htsum(i,j,nt,k) / hlyr(i,j,k)
                   endif
                enddo   ! i
             enddo      ! j
          enddo         ! nt
       enddo            ! k
 
  end subroutine vertical_remap
        
!****************************************************************************

end module glide_thck

