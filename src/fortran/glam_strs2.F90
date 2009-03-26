module glam_strs2

!whl - Use statements modified for glimmer
!whl - NOTES:
!    - Do not need in, sprs2_dp, or sa
!    - rhoi = 917 in glam_physcon, 910 in glimmer_physcon
!    - Added evs0 and lambda0 to glimmer_paramets
!      Glimmer defines tau0 in terms of vis0; glam defines vis0 in terms of tau0 = rho*h*thk0

!whl - to do - Make sure that any hardwired constants in the code do not assume
!              a scaling different from the glimmer scaling.

!!use glam_general, only : dp, in, sprs2_dp
!!use glam_physcon, only : gn, rhoi, rhoo, grav, pi, scyr
!!use glam_paramets, only : thk0, len0, vel0, vis0, tim0, tau0, lambda0, evs0
!!use glam_funits, only : ulog, unin, betafile

use glimmer_paramets, only : dp
use glimmer_physcon,  only : gn, rhoi, rhoo, grav, pi, scyr
use glimmer_paramets, only : thk0, len0, vel0, vis0, tim0, tau0, lambda0, evs0
use glimmer_log,      only : write_log

implicit none

!whl - The following were moved up from glam_strs2 subroutines
  integer, save :: locplusup
  logical, save :: lateralboundry = .false.
  integer, dimension(6), save :: loc_latbc

  real (kind = dp), allocatable, dimension(:,:,:),     save :: flwafact
  real (kind = dp), allocatable, dimension(:),         save :: dups

!whl - The following are from module strswk
  real (kind = dp), allocatable, dimension(:,:,:,:,:), save :: corr
  real (kind = dp), allocatable, dimension(:,:,:,:),   save :: usav
  real (kind = dp), allocatable, dimension(:,:,:),     save :: tvel
  real (kind = dp), allocatable, dimension(:,:),       save :: valubbc
  real (kind = dp), allocatable, dimension(:),         save :: dup, dupm

  integer, dimension(:,:), allocatable :: typebbc
  integer, dimension(:,:), allocatable :: uindx, umask

  integer, parameter :: mnbdy = 1, noice = 0
  integer, parameter :: boundarys = 2
  integer, dimension(2), parameter :: strdomain = (/ mnbdy, boundarys /)

  real (kind = dp), parameter :: effstrminsq = (1.0e-20_dp * tim0)**2

! *sp* 'p' are later defined as variants on 'gn', or Glen's 'n'(=3),
!       e.g. p1=gn+1=4, etc ... 
  real (kind = dp) :: p1, p2, dew2, dns2, dew4, dns4

!whl - The following are from module strcalcs
!whl - Should these have the save attribute?

  real (kind = dp) :: cdxdy
  real (kind = dp), dimension(2) :: cdxdx
  real (kind = dp), dimension(:),   allocatable :: cdsds, cds
  real (kind = dp), dimension(:,:), allocatable :: cvert, cdsdx, fvert

  real (kind = dp), dimension(:), allocatable :: dsigmadew, dsigmadns
  real (kind = dp), dimension(:), allocatable :: d2sigmadew2, d2sigmadns2, d2sigmadewdns
  real (kind = dp) :: d2sigmadewdsigma, d2sigmadnsdsigma

  real (kind = dp), dimension(2), parameter ::   &
           oneorfour = (/ 1.0_dp, 4.0_dp /),     &
           fourorone = (/ 4.0_dp, 1.0_dp /),     &
           oneortwo  = (/ 1.0_dp, 2.0_dp /),     &
           twoorone  = (/ 2.0_dp, 1.0_dp /)

!*sp* coeff. for forward diff. template
  real (kind = dp), dimension(3), parameter ::   &
           onesideddiff = (/ -3.0_dp, 4.0_dp, -1.0_dp /)

!whl - The following are from module geomderv
  real (kind = dp), dimension(:,:), allocatable :: &
    d2thckdew2, d2usrfdew2, d2thckdns2, d2usrfdns2, d2thckdewdns, d2usrfdewdns

!whl - The following are from module pcgdwk.
!whl - save attribute?

  real (kind = dp), dimension(:), allocatable :: pcgval, rhsd, answ
  integer, dimension(:), allocatable :: pcgcol, pcgrow
  integer, dimension(2) :: pcgsize
  integer :: ct

!whl - The following is from glam_funits
  integer, parameter :: unin = 90

!***********************************************************************

contains

!***********************************************************************
!whl - added a subroutine to be called at initialization (in lieu of 'first')

subroutine glam_velo_fordsiapstr_init (ewn,   nsn,   upn,    &
                                       dew,   dns,           &
                                       sigma, stagsigma)

! Allocate arrays and initialize variables.

    implicit none

    integer, intent(in) :: ewn, nsn, upn
    real (kind = dp), intent(in) :: dew, dns

    real (kind = dp), dimension(:), intent(in)  :: sigma
    real (kind = dp), dimension(:), intent(out)  :: stagsigma

    integer :: up

!whl - to do - Many tasks currently done in glam.F90 should be done here or elsewhere.
! - Read in namelist values from glam.nml
! - Read in restart values?
! - Allocate arrays if not already done in glimmer
! - Define initial thickness
! - Mask the thickness
! - Initialize the horizontal remapping
! - Add ppm thickness routines

    allocate( dupm(upn) )
    allocate( cvert(upn,2) )
    allocate( cdsdx(upn,2) )
    allocate( cdsds(upn) )
    allocate( cds(upn) )
    allocate( fvert(upn,3) )

    ! *sp* Note that "dup" is defined as a vector and is read in from file like sigma coordinate 
    ! (NOT a scalar)
    dupm = - 0.25_dp / dup

    !eta = (/ (dup * real(up-1,dp), up = 1, upn) /)
!whl - eta is currently not used
!    eta = (/ (dup(up) * real(up-1,dp), up = 1, upn) /)

!whl - to do - Make sure sigma levels are evenly spaced

    stagsigma = (sigma(1:upn-1) + sigma(2:upn)) / 2.0_dp

    ! *sp*  p1 = -1/n   - used with rate factor in eff. visc. def.
    ! *sp*  p2 = (1-n)/2n   - used with eff. strain rate in eff. visc. def. 
    p1 = -1.0_dp / real(gn,dp)      
    p2 = (1.0_dp - real(gn,dp)) / (2.0_dp * real(gn,dp))

    dew2 = 2.0_dp * dew; dns2 = 2.0_dp * dns        ! *sp* 2x the standard grid spacing
    dew4 = 4.0_dp * dew; dns4 = 4.0_dp * dns        ! *sp* 4x the standard grid spacing

    allocate(dsigmadew(upn),  dsigmadns(upn))
    allocate(d2sigmadew2(upn),d2sigmadns2(upn),d2sigmadewdns(upn))

    allocate (d2thckdew2(ewn-1,nsn-1),d2thckdns2(ewn-1,nsn-1),d2thckdewdns(ewn-1,nsn-1), &
              d2usrfdew2(ewn-1,nsn-1),d2usrfdns2(ewn-1,nsn-1),d2usrfdewdns(ewn-1,nsn-1))

    allocate(valubbc(ewn-1,nsn-1),typebbc(ewn-1,nsn-1))
    allocate(umask(ewn-1,nsn-1)) ! this will be moved to main

!whl - moved from findefvsstr
    allocate(flwafact(1:upn-1,ewn,nsn))
    flwafact = 0.0_dp

! *sp* determine constants used in various FD calculations associated with 'findcoefst'   
! NOTE: there is some question about the definitions here vs. in write-up (see notes in subroutine)

!whl - moved from findcoefstr
     call calccoeffsinit(upn, dew, dns)

!whl - moved from vertintg
    allocate(dups(upn)) 
    dups = (/ (sigma(up+1) - sigma(up), up=1,upn-1), 0.0d0 /)

end subroutine glam_velo_fordsiapstr_init

!***********************************************************************

subroutine glam_velo_fordsiapstr(ewn,      nsn,    upn,  &
                                 dew,      dns,          &
                                 sigma,    stagsigma,    &
                                 thck,     usrf,         &
                                 lsrf,     topg,         &
                                 dthckdew, dthckdns,     &
                                 dusrfdew, dusrfdns,     & 
                                 dlsrfdew, dlsrfdns,     &
                                 stagthck, flwa,         & 
                                 mintauf,                & 
                                 whichbabc,              &
                                 whichefvs,              &
                                 whichresid,             &
                                 uvel,     vvel,         &
                                 uflx,     vflx,         &
                                 efvs,                   &
                                 gdsx,     gdsy)

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  real (kind = dp), intent(in) :: dew, dns

  real (kind = dp), dimension(:),     intent(in)  :: sigma, stagsigma
  real (kind = dp), dimension(:,:),   intent(in)  :: thck, usrf, lsrf, topg
  real (kind = dp), dimension(:,:),   intent(in)  :: dthckdew, dthckdns
  real (kind = dp), dimension(:,:),   intent(in)  :: dusrfdew, dusrfdns
  real (kind = dp), dimension(:,:),   intent(in)  :: dlsrfdew, dlsrfdns
  real (kind = dp), dimension(:,:),   intent(in)  :: stagthck
  real (kind = dp), dimension(:,:),   intent(in)  :: minTauf
  real (kind = dp), dimension(:,:,:), intent(in)  :: flwa

!whl - to do - Merge whichbabc with whichbtrc?
  integer, intent(in) :: whichbabc
  integer, intent(in) :: whichefvs
  integer, intent(in) :: whichresid

  real (kind = dp), dimension(:,:,:), intent(out) :: uvel, vvel
  real (kind = dp), dimension(:,:),   intent(out) :: uflx, vflx
  real (kind = dp), dimension(:,:,:), intent(out) :: gdsx, gdsy
  real (kind = dp), dimension(:,:,:), intent(out) :: efvs

  integer :: ew, ns, up

  real (kind = dp), parameter :: minres = 1.0d-5 
  real (kind = dp), save, dimension(2) :: resid  

  real (kind = dp), save :: etaorder = 1.0073e-6_dp  ! linear 1.0073e-6_dp 

  integer, parameter :: cmax = 100 
   
  integer :: count, linit 

  character(len=100) :: message

!whl - Moved initialization stuff to glam_velo_fordsiapstr_init

!whl - Took these out of initialization because these will change
!      for prognostic thickness

  ! *sp* geometric 1st deriv. for generic input variable 'ipvr',
  !      output as 'opvr' (includes 'upwinding' for boundary values)
  call geom2ders(ewn, nsn, dew, dns, usrf, stagthck, d2usrfdew2, d2usrfdns2)	
  call geom2ders(ewn, nsn, dew, dns, thck, stagthck, d2thckdew2, d2thckdns2)

  ! *sp* geometric (2nd) cross-deriv. for generic input variable 'ipvr', output as 'opvr'
  call geom2derscros(dew, dns, thck, stagthck, d2thckdewdns)
  call geom2derscros(dew, dns, usrf, stagthck, d2usrfdewdns)

  valubbc = 0.0_dp
  typebbc = 0.0_dp

  ! *sp* make a 2d array identifying if the associated point has zero thickness,
  !      has non-zero thickness and is interior, or has non-zero thickness
  !      and is along a boundary

  umask = maskvelostr(ewn, nsn, thck, stagthck)

  allocate(uindx(ewn-1,nsn-1))

  ! *sp* if a point from the 2d array 'mask' is associated with non-zero ice thickness,
  !      either a boundary or interior point, give it a unique number. If not, give it a zero			 
  uindx = indxvelostr(ewn, nsn, upn,  &
                      umask,pcgsize(1))

  allocate(tvel(upn,ewn-1,nsn-1)) 

  tvel = 0.0_dp 

  ! *sp* calculate scaled driving stresses in x,y directions	
  call calcgdststr(ewn,      nsn,        &
                   sigma,    stagsigma,  &
                   thck,                 &
                   dusrfdew, dusrfdns,   &
                   gdsx,     gdsy)

  ! *sp* allocate space for variables used by 'mindcrash' function
  allocate(corr(upn,ewn-1,nsn-1,2,2),usav(upn,ewn-1,nsn-1,2))

  ! *sp* an initial guess at the size of the sparse matrix
  pcgsize(2) = pcgsize(1) * 19

  ! *sp* allocate space matrix variables
  allocate (pcgrow(pcgsize(2)),pcgcol(pcgsize(2)),rhsd(pcgsize(1)), &
            answ(pcgsize(1)),pcgval(pcgsize(2)))

  !whl - Removed subroutine findbtrcstr; superseded by calcbetasquared
  
  resid = 1.0_dp
  count = 1
  linit = 0;

  ! *sp* main iteration on stress, vel, and eff. visc. solutions,

  do while ( maxval(resid) > minres .and. count < cmax)	
!  do while ( resid(1) > minres .and. count < cmax)  ! *sp* for 1d solutions (d*/dy=0) 

    ! *sp* effective viscosity calculation, based on previous estimate for vel. field
    call findefvsstr(ewn,  nsn,  upn,      &
                     stagsigma,  count,    &
                     whichefvs,  efvs,     &
                     uvel,       vvel,     &
                     flwa,       thck,     &
                     dusrfdew,   dthckdew, &
                     dusrfdns,   dthckdns, &
                     umask)

    ! *sp* calculation of coeff. for stress balance calc. 
    call findcoefstr(ewn,  nsn,   upn,            &
                     dew,  dns,   sigma,          &
                     2,           efvs,           &
                     vvel,        uvel,           &
                     thck,        dusrfdns,       &
                     dusrfdew,    dthckdew,       &
                     d2usrfdew2,  d2thckdew2,     &
                     dusrfdns,    dthckdns,       &
                     d2usrfdns2,  d2thckdns2,     &
                     d2usrfdewdns,d2thckdewdns,   &
                     dlsrfdew,    dlsrfdns,       &
                     stagthck,    whichbabc,      &
                     valubbc,     typebbc,        &
                     uindx,       umask,          &
                     lsrf,        topg,           &
                     minTauf,     flwa )

    ! *sp* solve 'Ax=b' for the across-flow velocity component
    tvel = slapsolvstr(ewn,  nsn,   upn,  &
                       vvel, uindx, linit)

!   !*sp* must be ACTIVE for test case (OR commented OUT for ice shelf test case) 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     ! periodic boundary conditions ...
!     tvel(:,1:ewn-1,nsn-1) = tvel(:,1:ewn-1,2)
!     tvel(:,1:ewn-2,1) = tvel(:,1:ewn-2,nsn-2)
!     tvel(:,ewn-1,1:nsn-1) = tvel(:,2,1:nsn-1)
!     tvel(:,1,1:nsn-1) = tvel(:,ewn-2,1:nsn-1)

!     tvel(:,1:ewn-1,nsn-2) = tvel(:,1:ewn-1,3)    ! if using rempping for dH/dt (domain def. is different) 
!     tvel(:,1:ewn-1,2) = tvel(:,1:ewn-1,nsn-3)
!     tvel(:,ewn-2,1:nsn-1) = tvel(:,3,1:nsn-1) 
!     tvel(:,2,1:nsn-1) = tvel(:,ewn-3,1:nsn-1)
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !  ... however, call it 'tvel' so that 'vvel', the solution from the previous iteration,
    !  gets passed to the solver for 'uvel' (rather than the new value of 'vvel' ... why?) 

    ! *sp* calculation of coeff. for stress balance calc. 

    ! - along-flow stress balance -
    call findcoefstr(ewn,  nsn,   upn,            &
                     dew,  dns,   sigma,          &
                     1,           efvs,           &
                     uvel,        vvel,           &
                     thck,        dusrfdew,       &
                     dusrfdew,    dthckdew,       &
                     d2usrfdew2,  d2thckdew2,     &
                     dusrfdns,    dthckdns,       &
                     d2usrfdns2,  d2thckdns2,     &
                     d2usrfdewdns,d2thckdewdns,   &
                     dlsrfdew,    dlsrfdns,       &
                     stagthck,    whichbabc,      &
                     valubbc,     typebbc,        &
                     uindx,       umask,          &
                     lsrf,        topg,           &
                     minTauf,     flwa )


    ! *sp* solve 'Ax=b' for along-flow velocity component
    uvel = slapsolvstr(ewn,  nsn,   upn,  &
                       uvel, uindx, linit)

    ! *sp* correct the velocity estimates using the "unstable manifold" correction scheme
    uvel = mindcrshstr(1,whichresid,uvel,count,resid(1))
    vvel = mindcrshstr(2,whichresid,tvel,count,resid(2))

!   !*sp* must be ACTIVE for test case (OR commented OUT for ice shelf test case) 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     ! periodic boundary conditions ...
!     uvel(:,1:ewn-1,nsn-1) = uvel(:,1:ewn-1,2)
!     uvel(:,1:ewn-1,1) = uvel(:,1:ewn-1,nsn-2)
!     uvel(:,1,1:nsn-1) = uvel(:,ewn-2,1:nsn-1)
!     uvel(:,ewn-1,1:nsn-1) = uvel(:,2,1:nsn-1)

!     uvel(:,1:ewn-1,nsn-2) = uvel(:,1:ewn-1,3)     ! if using rempping for dH/dt (domain def. is different) 
!     uvel(:,1:ewn-1,2) = uvel(:,1:ewn-1,nsn-3)
!     uvel(:,ewn-2,1:nsn-1) = uvel(:,3,1:nsn-1)
!     uvel(:,2,1:nsn-1) = uvel(:,ewn-3,1:nsn-1)
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    count = count + 1

    ! *sp* output status of iteration: iteration number, max residual, and location of max residual
    print '(" * strs ",i3,3g20.6)', count, resid(1), resid(2), minres

!whl - write this info to the log file
    write(message,'(" * strs ",i3,3g20.6)') count, resid(1), resid(2), minres
    call write_log (message)
  end do

!whl - to do - Compute the stresses outside this module?
!              The stresses are no longer passed as arguments.
  ! *sp* calc. stress fields from converged vel. fields
!  call calcstrsstr(ewn,  nsn,  upn,       &
!                   stagsigma,  thck,      &
!                   dusrfdew,   dusrfdns,  &
!                   dthckdew,   dthckdns,  &
!                   uvel,       vvel,      &
!                   efvs,                  &
!                   tauxx,      tauyy,     &
!                   tauxy,      tau,       &
!                   tauxz,      tauyz,     &
!                   umask)

  do ns = 1,nsn-1
  do ew = 1,ewn-1 
     ! *sp* calc. fluxes from converged vel. fields (for input to thickness evolution subroutine)
     if (umask(ew,ns) > 0) then	
        uflx(ew,ns) = vertintg(upn, sigma, uvel(:,ew,ns)) * stagthck(ew,ns)
        vflx(ew,ns) = vertintg(upn, sigma, vvel(:,ew,ns)) * stagthck(ew,ns)
     end if
  end do
  end do

  deallocate(tvel)
  deallocate(uindx,corr,usav)
  deallocate(pcgval,pcgrow,pcgcol,rhsd,answ)

end subroutine glam_velo_fordsiapstr

!***********************************************************************

function indxvelostr(ewn,  nsn,  upn,  &
                     mask, pointno)	

! *sp* if a point from the 2d array 'mask' is associated with non-zero ice thickness, 
!      either a boundary or interior point, give it a unique number. If not, give it a zero.

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  integer, intent(in), dimension(:,:) :: mask
  integer, intent(out) :: pointno

  integer :: ew, ns
  integer, dimension(size(mask,1),size(mask,2)) :: indxvelostr

  pointno = 1

  do ew = 1,ewn-1
  do ns = 1,nsn-1			
    if (any(strdomain == mask(ew,ns))) then
      indxvelostr(ew,ns) = pointno
      pointno = pointno + 1
    else
      indxvelostr(ew,ns) = 0
    end if
  end do
  end do

! add twop ghost points at upper and lower boundaries

  pointno = (pointno - 1) * (upn + 2)

  return

end function indxvelostr

!***********************************************************************

subroutine calcgdststr(ewn,      nsn,        &
                       sigma,    stagsigma,  &
                       thck,                 &
                       dusrfdew, dusrfdns,   &
                       gdsx,     gdsy)

! *sp* calculate scaled driving stresses in x,y directions

  implicit none

  integer, intent(in) :: ewn, nsn
  real (kind = dp), intent(in),  dimension(:)      :: sigma, stagsigma
  real (kind = dp), intent(in),  dimension(:,:)    :: thck, dusrfdew, dusrfdns
  real (kind = dp), intent(out), dimension(:,:,:)  :: gdsx, gdsy

  integer :: ew, ns
  real (kind = dp) :: fx, fy

! *sp* the factor of 1/4 is used to average sfc slopes over 4 grid points (below), 
! ... the minus sign is to give a positive value for the driving stress, 
! ... the terms in ( ) give the ratio of thickness to horiz length scales, H0/X0	

  real (kind = dp), parameter :: c =  -0.25_dp * (thk0 / len0)	

  do ns = 2, nsn-1
  do ew = 2, ewn-1

!whl - Change to thck > 0?
   if (thck(ew,ns) .ne. 0.0_dp) then 

    ! *sp* gives the scaled driving stress. e.g. rho*g*H*dh/dx * (1/Tau_0) * (H0/X0) = H'*dh'/dx' 
    fx = c * thck(ew,ns) * sum(dusrfdew(ew-1:ew,ns-1:ns))
    fy = c * thck(ew,ns) * sum(dusrfdns(ew-1:ew,ns-1:ns))

!whl - Reduce gdsx and gdsy to rank 2 arrays?

    ! *sp* make rank 3 arrays containing scaled driv. stresses in x, y, z directions 
    !      (constant in z direction)
    gdsx(:,ew,ns) = stagsigma * fx
    gdsy(:,ew,ns) = stagsigma * fy
  
  else

    gdsx(:,ew,ns) = 0.0_dp
    gdsy(:,ew,ns) = 0.0_dp

   end if  ! thck .ne. 0.0
  end do
  end do

  return

end subroutine calcgdststr

!***********************************************************************

subroutine findefvsstr(ewn,  nsn, upn,       &
                       stagsigma, count,     &
                       whichefvs, efvs,      &
                       uvel,      vvel,      &
                       flwa,      thck,      &
                       dusrfdew,  dthckdew,  &
                       dusrfdns,  dthckdns,  &
                       mask)

! *sp* calculate the eff. visc.	
! (NOTE: this version looks to AGREE with the version I came up with in the write-up, 
!        with the correct number of cross terms.)

  implicit none 

  integer, intent(in) :: ewn, nsn, upn 
  real (kind = dp), intent(in), dimension(:)     :: stagsigma
  real (kind = dp), intent(in), dimension(:,:,:) :: uvel, vvel, flwa
  real (kind = dp), intent(inout), dimension(:,:,:) :: efvs
  real (kind = dp), intent(in), dimension(:,:) :: thck, dthckdew, dusrfdew, & 
                                                  dusrfdns, dthckdns
  integer, intent(in), dimension(:,:) :: mask
  integer, intent(in) :: whichefvs, count
       
  integer :: ew, ns, up

  real (kind = dp), dimension(size(efvs,1)) :: effstr, ugradup, vgradup, &
                                               ugradew, ugradns, vgradew, vgradns

  integer, dimension(2) :: mew, mns

! *sp* this is the 1/4(X0/H0)^2 factor in front of the term ((dv/dz)^2+(du/dz)^2) 
  real (kind = dp), parameter :: f1 = 0.25_dp * (len0 / thk0)**2

  select case(whichefvs)

  case(0)			! *sp* calculate eff. visc. from eff. strain rate, etc

  if (1 == count) then
    do ns = 2,nsn-1; do ew = 2,ewn-1
    if (thck(ew,ns) > 0.0_dp) then
      ! *sp* term: 1/2*A^(-1/n)
      forall (up = 1:upn-1) flwafact(up,ew,ns) = 0.5_dp * (sum(flwa(up:up+1,ew,ns)) / 2.0_dp)**p1
    end if; end do; end do
  end if

  do ns = 2,nsn-1
  do ew = 2,ewn-1
   if (thck(ew,ns) > 0.0_dp) then

    ugradup = vertideriv(upn, hsum(uvel(:,ew-1:ew,ns-1:ns)), thck(ew,ns))
    vgradup = vertideriv(upn, hsum(vvel(:,ew-1:ew,ns-1:ns)), thck(ew,ns))

    ugradew = horizderiv(upn,  stagsigma,                &
                         sum(uvel(:,ew-1:ew,ns-1:ns),3), &
                         dew4, ugradup,                  &             
                         sum(dusrfdew(ew-1:ew,ns-1:ns)), &
                         sum(dthckdew(ew-1:ew,ns-1:ns)))

    vgradew = horizderiv(upn,  stagsigma,                &
                         sum(vvel(:,ew-1:ew,ns-1:ns),3), &
                         dew4, vgradup,                  &             
                         sum(dusrfdew(ew-1:ew,ns-1:ns)), &
                         sum(dthckdew(ew-1:ew,ns-1:ns)))

    ugradns = horizderiv(upn,  stagsigma,               &
                         sum(uvel(:,ew-1:ew,ns-1:ns),2), &
                         dns4, ugradup,                  &                              
                         sum(dusrfdns(ew-1:ew,ns-1:ns)), &
                         sum(dthckdns(ew-1:ew,ns-1:ns)))

    vgradns = horizderiv(upn,  stagsigma,               &
                         sum(vvel(:,ew-1:ew,ns-1:ns),2), &
                         dns4, vgradup,                  &                              
                         sum(dusrfdns(ew-1:ew,ns-1:ns)), &
                         sum(dthckdns(ew-1:ew,ns-1:ns)))

    ! *sp* eff. strain rate (squared)
    effstr = ugradew**2 + vgradns**2 + ugradew*vgradns + &
             0.25_dp * (vgradew + ugradns)**2 + &
             f1 * (ugradup**2 + vgradup**2)

    ! *sp* set eff. strain rate (squared) to some min value where
    !      it falls below some threshold value, 'effstrminsq'
    where (effstr < effstrminsq)
      effstr = effstrminsq
    end where

    ! *sp* p2 = (1-n)/2n, where the factor of 1/2 comes from taking 
    !      the sqr root of the squared eff. strain rate ...
    efvs(:,ew,ns) = flwafact(:,ew,ns) * effstr**p2	

    else  
      efvs(:,ew,ns) = effstrminsq
    end if
   end do
   end do

  case(1)			! *sp* set a constant eff. visc. 
	
    efvs = 1 / 100.0_dp    	

!whl - changed default to case 2
  case (2)

    efvs = effstrminsq
  
  end select
         
  return

end subroutine findefvsstr

!***********************************************************************

subroutine calcstrsstr(ewn,  nsn,  upn,       &
                       stagsigma,  thck,      &
                       dusrfdew,   dusrfdns,  &
                       dthckdew,   dthckdns,  &
                       uvel,       vvel,      &
                       efvs,                  &
                       tauxx,      tauyy,     &
                       tauxy,      tau,       &
                       tauxz,      tauyz,     &
                       mask)

! *sp* calc. stress based on vel. and eff. visc. solutions
!whl - no longer called from within glam_strs2

  implicit none
 
  integer, intent(in) :: ewn, nsn, upn 

  real (kind = dp), intent(in), dimension(:)     :: stagsigma
  real (kind = dp), intent(in), dimension(:,:,:) :: efvs, uvel, vvel 
  real (kind = dp), intent(in), dimension(:,:) :: thck, dusrfdew, &
                                                  dusrfdns, dthckdew, dthckdns 
  integer, intent(in), dimension(:,:) :: mask

  integer :: ew, ns, up
  real (kind = dp), intent(out), dimension(:,:,:) :: tauxx, tauyy, tauxy, &
                                                     tauxz, tauyz, tau

  real (kind = dp), parameter :: f1 = len0 / thk0

  integer, dimension(2) :: mew, mns
  
  do ns = 2,nsn-1
  do ew = 2,ewn-1; 
   if (thck(ew,ns) > 0.0_dp) then

    tauxz(:,ew,ns) = vertideriv(upn, hsum(uvel(:,ew-1:ew,ns-1:ns)), thck(ew,ns))

    tauyz(:,ew,ns) = vertideriv(upn, hsum(vvel(:,ew-1:ew,ns-1:ns)), thck(ew,ns))

    tauxx(:,ew,ns) = horizderiv(upn,  stagsigma,                &
                                sum(uvel(:,ew-1:ew,ns-1:ns),3), &
                                dew4, tauxz(:,ew,ns),           &             
                                sum(dusrfdew(ew-1:ew,ns-1:ns)), &
                                sum(dthckdew(ew-1:ew,ns-1:ns)))

    tauyy(:,ew,ns) = horizderiv(upn,  stagsigma,                &
                                sum(vvel(:,ew-1:ew,ns-1:ns),2), &
                                dns4, tauyz(:,ew,ns),           &             
                                sum(dusrfdns(ew-1:ew,ns-1:ns)), &
                                sum(dthckdns(ew-1:ew,ns-1:ns)))

    tauxy(:,ew,ns) = horizderiv(upn,  stagsigma,                &
                                sum(uvel(:,ew-1:ew,ns-1:ns),2), &
                                dns4, tauxz(:,ew,ns),           &             
                                sum(dusrfdns(ew-1:ew,ns-1:ns)), &
                                sum(dthckdns(ew-1:ew,ns-1:ns))) + &
                     horizderiv(upn,  stagsigma,                &
                                sum(vvel(:,ew-1:ew,ns-1:ns),3), &
                                dew4, tauyz(:,ew,ns),           &             
                                sum(dusrfdew(ew-1:ew,ns-1:ns)), &
                                sum(dthckdew(ew-1:ew,ns-1:ns)))
 
   else

    tauxz(:,ew,ns) = 0.0_dp
    tauyz(:,ew,ns) = 0.0_dp
    tauxx(:,ew,ns) = 0.0_dp
    tauyy(:,ew,ns) = 0.0_dp
    tauxy(:,ew,ns) = 0.0_dp

   end if
  end do
  end do

  tauxz = f1 * efvs * tauxz
  tauyz = f1 * efvs * tauyz
  tauxx = 2.0_dp * efvs * tauxx
  tauyy = 2.0_dp * efvs * tauyy
  tauxy = efvs * tauxy  

  tau = sqrt(tauxz**2 + tauyz**2 + tauxx**2 + tauyy**2 + tauxx*tauyy + tauxy**2)
         
  return

end subroutine calcstrsstr

!***********************************************************************

function vertideriv(upn, varb, thck)

  implicit none 

  integer, intent(in) :: upn
  real (kind = dp), intent(in), dimension(:) :: varb
  real (kind = dp), intent(in) :: thck    

  real (kind = dp), dimension(size(varb)-1) :: vertideriv

! *sp* 'dupm' defined as -1/(2*del_sigma), in which case it seems like 
!       there should be a '-' in front of this expression...

  vertideriv(1:upn-1) = dupm * (varb(2:upn) - varb(1:upn-1)) / thck


  return

end function vertideriv

!***********************************************************************

function horizderiv(upn,     stagsigma,   &
                    varb,    grid,        &
                    dvarbdz, dusrfdx, dthckdx)

  implicit none
  
  integer, intent(in) :: upn
  real (kind = dp), dimension(:), intent(in) :: stagsigma
  real (kind = dp), dimension(:,:), intent(in) :: varb
  real (kind = dp), dimension(:), intent(in) :: dvarbdz
  real (kind = dp), intent(in) :: dusrfdx, dthckdx, grid

  real (kind = dp) :: horizderiv(size(varb,1)-1)
  
  ! *sp* where does this factor of 1/4 come from ... averaging? 
  horizderiv = (varb(1:upn-1,2) + varb(2:upn,2) - varb(1:upn-1,1) - varb(2:upn,1)) / grid - &    
                dvarbdz * (dusrfdx - stagsigma * dthckdx) / 4.0_dp

  return

end function horizderiv

!***********************************************************************

function getlocrange(upn, indx)

  implicit none

  integer, intent(in) :: upn 
  integer, intent(in) :: indx
  integer, dimension(2) :: getlocrange

  getlocrange = (indx - 1) * (upn + 2) + 1 + (/ 1, upn /)

  return

end function getlocrange

!***********************************************************************

function getlocation(upn, indx)

!whl - This function is not currently used.

  implicit none

  integer, intent(in) :: upn
  integer, intent(in) :: indx
  integer :: getlocation

  getlocation = (indx - 1) * (upn + 2) + 1

  return

end function getlocation

!***********************************************************************

function getlocationarray(ewn, nsn, upn,  &
                          mask, mnbdy, boundarys)
    implicit none

    integer, intent(in) :: ewn, nsn, upn 
    integer, dimension(:,:), intent(in) :: mask
    integer, intent(in) :: mnbdy, boundarys
    integer, dimension(ewn-1,nsn-1) :: getlocationarray, temparray
    integer :: cumsum, ew, ns

    cumsum = 0

    do ew=1,ewn-1
        do ns=1,nsn-1
        !*sp* if in main body, there are upn+2 columns ...
        if ( mask(ew,ns) .eq. mnbdy .or. mask(ew,ns) .eq. boundarys ) then
            cumsum = cumsum + ( upn + 2 )
            getlocationarray(ew,ns) = cumsum
            temparray(ew,ns) = upn + 2
        else
            getlocationarray(ew,ns) = 0
            temparray(ew,ns) = 1
        end if
        end do
    end do

    getlocationarray = ( getlocationarray + 1 ) - temparray

    return

end function getlocationarray

!***********************************************************************

function slapsolvstr(ewn, nsn, upn, &
                     vel, uindx, its)

! *sp* routine to solve Ax=b sparse matrix problem 

  use glimmer_log, only: glimmer_get_logunit

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  real (kind = dp), dimension(:,:,:), intent(in) :: vel
  integer, dimension(:,:), intent(in) :: uindx
  real (kind = dp), dimension(size(vel,1),size(vel,2),size(vel,3)) :: slapsolvstr
  integer, intent(inout) :: its

  integer :: ew, ns 

  real (kind = dp), dimension(:), allocatable :: rwork
  integer, dimension(:), allocatable :: iwork

  real (kind = dp), parameter :: tol = 1.0e-12_dp
  real (kind = dp) :: err

  integer, parameter :: isym = 0, itol = 2, itmax = 100
  integer, dimension(2) :: loc
  integer :: iter, ierr, mxnelt

! ** move to values subr   

  pcgsize(2) = ct - 1
                    
  call ds2y(pcgsize(1),pcgsize(2),pcgrow,pcgcol,pcgval,isym)                   

!** plot the matrix to check that it has the correct form
!*sp* this outputs a .txt file of sparse matrixes

!call dcpplt(pcgsize(1),pcgsize(2),pcgrow,pcgcol,pcgval,isym,glimmer_get_logunit())

  mxnelt = 60 * pcgsize(1); allocate(rwork(mxnelt),iwork(mxnelt))

!**     solve the problem using the SLAP package routines     
!**     -------------------------------------------------
!**     n ... order of matrix a (in)
!**     b ... right hand side vector (in)                        
!**     x ... initial quess/final solution vector (in/out)                        
!**     nelt ... number of non-zeroes in A (in)
!**     ia, ja ... sparse matrix format of A (in)
!**     a ... matrix held in SLAT column format (in)
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

! *sp* initial estimate for vel. field?
  do ns = 1,nsn-1
  do ew = 1,ewn-1
   if (uindx(ew,ns) /= 0) then
    loc = getlocrange(upn, uindx(ew,ns))
    answ(loc(1):loc(2)) = vel(:,ew,ns)
    answ(loc(1)-1) = vel(1,ew,ns)
    answ(loc(2)+1) = vel(upn,ew,ns)
   end if
  end do
  end do

  call dslucs(pcgsize(1),rhsd,answ,pcgsize(2),pcgrow,pcgcol,pcgval, &
              isym,itol,tol,itmax,iter,err,ierr,0,rwork,mxnelt,iwork,mxnelt)

  if (ierr .ne. 0) then
    print *, 'pcg error ', ierr, itmax, iter, tol, err 
    ! stop
  end if

  deallocate(rwork,iwork)

  do ns = 1,nsn-1
  do ew = 1,ewn-1
     if (uindx(ew,ns) /= 0) then
       loc = getlocrange(upn, uindx(ew,ns))
       slapsolvstr(:,ew,ns) = answ(loc(1):loc(2))
     else 
       slapsolvstr(:,ew,ns) = 0.0d0
     end if
  end do
  end do

  its = its + iter

  return

end function slapsolvstr


!***********************************************************************

function mindcrshstr(pt,whichresid,vel,count,resid)

! *sp* function to perform 'unstable manifold correction' - 
!      corrects velocity fields u, v, from an initial guess at u, v, and eff. visc. 
!      iteratively to eventually reach a consistent solution in which 
!          A(u_old,v_old) u_new = b(u_old,v_old)
!      (coeff. matrix A is function of vel field through visc, 
!       rhs b is function of vel field through bcs) 

  implicit none

  real (kind = dp), intent(in), dimension(:,:,:) :: vel
  integer, intent(in) :: count, pt, whichresid 

  real (kind = dp), intent(out) :: resid

  real (kind = dp), dimension(size(vel,1),size(vel,2),size(vel,3)) :: mindcrshstr

  real (kind = dp), parameter :: ssthres = 5.0_dp * pi / 6.0_dp, &
                                 critlimit = 10.0_dp / (scyr * vel0), &
                                 small = 1.0e-16_dp

  real (kind = dp), intrinsic :: abs, acos
  
  integer, dimension(2), save :: new = 1, old = 2
  integer :: locat(3)

  integer :: nr
  integer, dimension(size(vel,1),size(vel,2),size(vel,3)) :: vel_ne_0

  if (count == 1) then
    usav(:,:,:,pt) = 0.0d0
  end if

  corr(:,:,:,new(pt),pt) = vel - usav(:,:,:,pt)           

  if (count > 1) then

    where (acos((corr(:,:,:,new(pt),pt) * corr(:,:,:,old(pt),pt)) / &
          (abs(corr(:,:,:,new(pt),pt)) * abs(corr(:,:,:,old(pt),pt)) + small)) > &
           ssthres .and. corr(:,:,:,new(pt),pt) - corr(:,:,:,old(pt),pt) /= 0.0_dp )

      mindcrshstr = usav(:,:,:,pt) + &
                    corr(:,:,:,new(pt),pt) * abs(corr(:,:,:,old(pt),pt)) / &
                    abs(corr(:,:,:,new(pt),pt) - corr(:,:,:,old(pt),pt)) 

    elsewhere

      mindcrshstr = vel;
 
    end where
    
  else 

    mindcrshstr = vel;
   
  end if

  if (new(pt) == 1) then; old(pt) = 1; new(pt) = 2; else; old(pt) = 1; new(pt) = 2; end if

  select case (whichresid)

! *sp* residual calculation method; 

!whl - to do - changed from default to case(0) and renumbered other cases
!      maxval (0); maxval ignoring basal vel (1); mean value (2)

   case(0)
    resid = maxval( abs((usav(:,:,:,pt) - vel ) / vel ), MASK = vel .ne. 0.0_dp)  
    locat = maxloc( abs((usav(:,:,:,pt) - vel ) / vel ), MASK = vel .ne. 0.0_dp)

   case(1)
    nr = size( vel, dim=1 ) ! number of grid points in vertical ...
    resid = maxval( abs((usav(1:nr-1,:,:,pt) - vel(1:nr-1,:,:) ) / vel(1:nr-1,:,:) ),  &
                        MASK = vel .ne. 0.0_dp)
    locat = maxloc( abs((usav(1:nr-1,:,:,pt) - vel(1:nr-1,:,:) ) / vel(1:nr-1,:,:) ),  &
            MASK = vel .ne. 0.0_dp)

   case(2)
    nr = size( vel, dim=1 )
    vel_ne_0 = 0
    where ( vel .ne. 0.0_dp ) vel_ne_0 = 1

    !*sp* include basal velocities in resid. calculation
    resid = sum( abs((usav(:,:,:,pt) - vel ) / vel ), &
            MASK = vel .ne. 0.0_dp) / sum( vel_ne_0 )

    !*sp* ignore basal velocities in resid. calculation
    !	resid = sum( abs((usav(1:nr-1,:,:,pt) - vel(1:nr-1,:,:) ) / vel(1:nr-1,:,:) ),   &
    !           MASK = vel .ne. 0.0_dp) / sum( vel_ne_0(1:nr-1,:,:) )

    ! *sp* note that the location of the max residual is somewhat irrelevent here,
    !      since we are using the mean resid for convergence testing
    locat = maxloc( abs((usav(:,:,:,pt) - vel ) / vel ), MASK = vel .ne. 0.0_dp)

  end select

    usav(:,:,:,pt) = vel

!  print '("* ",i3,g20.6,3i6,g20.6)', count, resid, locat, vel(locat(1),locat(2),locat(3))

  return

end function mindcrshstr

!***********************************************************************

subroutine findcoefstr(ewn,  nsn,   upn,            &
                       dew,  dns,   sigma,          &
                       pt,          efvs,           &
                       thisvel,     othervel,       &
                       thck,        thisdusrfdx,    &
                       dusrfdew,    dthckdew,       &
                       d2usrfdew2,  d2thckdew2,     &
                       dusrfdns,    dthckdns,       &
                       d2usrfdns2,  d2thckdns2,     &
                       d2usrfdewdns,d2thckdewdns,   &
                       dlsrfdew,    dlsrfdns,       &
                       stagthck,    whichbabc,      &
                       valubbc,     typebbc,        &
                       uindx,       mask,           &
                       lsrf,        topg,           &
                       minTauf,     flwa )

! *sp* find coeffecients in stress balance equation ...
! ... also appears to contain important bc information
! and puts that information into the sparse coeff. matrix

  implicit none

  integer, intent(in) :: ewn, nsn, upn 
  real (kind = dp), intent(in) :: dew, dns
  real (kind = dp), dimension(:), intent(in) :: sigma

  real (kind = dp), dimension(:,:,:), intent(in) :: efvs, thisvel, &
                                                    othervel
  real (kind = dp), dimension(:,:), intent(in) :: stagthck, valubbc, thisdusrfdx, &
                                                  dusrfdew,   dthckdew,      &
                                                  d2usrfdew2, d2thckdew2,    &
                                                  dusrfdns,   dthckdns,      &
                                                  d2usrfdns2, d2thckdns2,    &
                                                  d2usrfdewdns,d2thckdewdns, &
                                                  dlsrfdew,   dlsrfdns,      &
                                                  thck, lsrf, topg 

  real (kind = dp), dimension(:,:), intent(in) :: minTauf

  real (kind = dp), dimension(:,:,:), intent(in) :: flwa

  integer, dimension(:,:), intent(in) :: mask, uindx, typebbc
  integer, intent(in) :: pt, whichbabc

  real (kind = dp), dimension(ewn-1,nsn-1) :: betasquared
  real (kind = dp), dimension(2,2,2) :: localefvs   
  real (kind = dp), dimension(3,3,3) :: localothervel
  real (kind = dp), dimension(upn) :: boundaryvel
  real (kind = dp) :: flwabar

  integer, dimension(ewn-1,nsn-1) :: loc_array
  integer, dimension(6) :: loc
  integer, dimension(3) :: shift
  integer :: ew, ns, up

  ! *sp* 'localplusup' is used to find locations in the sparse matrix relative to
  ! the present location. In general, it is approx. the number of grid cells in 
  ! the vertical direction  

  ct = 1

! *sp* read in mask to ID where to apply floating ice bc 
!      open(unin,file=shlfmask,status='old')         
!      read(unin,*) ((shelfmask(ew,ns), ew=1,ewn-1), ns=1,nsn-1)
!      close(unin)

  ! *sp* calculate/specify value of 'betasquared', to be input to subroutine 'bodyset', below
  call calcbetasquared (whichbabc,              &
                        dew,        dns,        &
                        ewn,        nsn,        &
                        lsrf,       topg,       &
                        thck,                   &
                        thisvel(upn,:,:),       &   
                        othervel(upn,:,:),      &
                        minTauf,                &
                        betasquared )

  do ns = 1,nsn-1
    do ew = 1,ewn-1 

     ! *sp* depth-ave rate factor, needed for ice shelf b.c. (below)
!     flwabar = sum( ( flwa(:,ew,ns) + flwa(:,ew,ns+1) + flwa(:,ew+1,ns) + flwa(:,ew+1,ns+1) ) / 4.0_dp, 1 ) / real(upn)    
     flwabar = 3.171d-24 / vis0    ! isothermal, temperate value 
!     flwabar = 1.8075e-25 / vis0    ! EISMINT-ROSS test 3-4 value

!     ! ...or, calculate the depth-averaged value (complicated code so as not to include funny values at boundaries)
!     flwabar = ( sum( flwa(:,ew,ns), 1, flwa(1,ew,ns)*vis0 < 1.0d-10 )/real(upn) + &
!               sum( flwa(:,ew,ns+1), 1, flwa(1,ew,ns+1)*vis0 < 1.0d-10 )/real(upn)  + &
!               sum( flwa(:,ew+1,ns), 1, flwa(1,ew+1,ns)*vis0 < 1.0d-10 )/real(upn)  + &
!               sum( flwa(:,ew+1,ns+1), 1, flwa(1,ew+1,ns+1)*vis0 < 1.0d-10 )/real(upn) ) / &
!               ( sum( flwa(:,ew,ns)/flwa(:,ew,ns), 1, flwa(1,ew,ns)*vis0 < 1.0d-10 )/real(upn) + &
!               sum( flwa(:,ew,ns+1)/flwa(:,ew,ns+1), 1, flwa(1,ew,ns+1)*vis0 < 1.0d-10 )/real(upn) + &
!               sum( flwa(:,ew+1,ns)/flwa(:,ew+1,ns), 1, flwa(1,ew+1,ns)*vis0 < 1.0d-10 )/real(upn) + &
!               sum( flwa(:,ew+1,ns+1)/flwa(:,ew+1,ns+1), 1, flwa(1,ew+1,ns+1)*vis0 < 1.0d-10 )/real(upn) )

    if( ns == 1 .and. ew == 1 ) then
           loc_array = getlocationarray(ewn, nsn, upn,  &
                                        mask, mnbdy,boundarys)
    end if

    loc(1) = loc_array(ew,ns)

    if (mnbdy == mask(ew,ns)) then

        call calccoeffs( upn,               sigma,              &
                         stagthck(ew,ns),                       &
                         dusrfdew(ew,ns),   dusrfdns(ew,ns),    &
                         dthckdew(ew,ns),   dthckdns(ew,ns),    &
                         d2usrfdew2(ew,ns), d2usrfdns2(ew,ns),  &
                         d2usrfdewdns(ew,ns),                   &
                         d2thckdew2(ew,ns), d2thckdns2(ew,ns),  &
                         d2thckdewdns(ew,ns))

        ! get index of cardinal neighbours
        loc(2) = loc_array(ew+1,ns)
        loc(3) = loc_array(ew-1,ns)
        loc(4) = loc_array(ew,ns+1)
        loc(5) = loc_array(ew,ns-1)

        ! *sp* this loop fills coeff. for all vertical layers, including sfc. and bed bcs
        do up = 1,upn

            ! *sp* adjust indices at sfc and bed so that correct values of 'efvs' and 'othervel'
            !      are passed to function

            shift = indshift( 0, ew, ns, up, ewn, nsn, upn, loc_array, stagthck(ew-1:ew+1,ns-1:ns+1) )  

            !whl - added several arguments to bodyset call
            call bodyset(ew,  ns,  up,        &
                         ewn, nsn, upn,       &
                         dew,      dns,       &
                         pt,       loc_array, &
                         loc,      stagthck,  &
                         thisdusrfdx,         &
                         dusrfdew, dusrfdns,  &
                         dlsrfdew, dlsrfdns,  &
                         efvs(up-1+shift(1):up+shift(1),ew:ew+1,ns:ns+1),  &
                         othervel(up-1+shift(1):up+1+shift(1),  &
                         ew-1+shift(2):ew+1+shift(2),  &
                         ns-1+shift(3):ns+1+shift(3)), &
                         betasquared(ew,ns) ) 

        end do  ! upn

!   !*sp* must be ACTIVE for test cases 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    else if (boundarys == mask(ew,ns)) then
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!     else if (boundarys == mask(ew,ns) .and. shelfmask(ew,ns) == 1 )then

   !*sp* must be commented OUT for test case OR commented IN for ice shelf test case) 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        call calccoeffs( upn,               sigma,              &
                         stagthck(ew,ns),                       &
                         dusrfdew(ew,ns),   dusrfdns(ew,ns),    &
                         dthckdew(ew,ns),   dthckdns(ew,ns),    &
                         d2usrfdew2(ew,ns), d2usrfdns2(ew,ns),  &
                         d2usrfdewdns(ew,ns),                   &
                         d2thckdew2(ew,ns), d2thckdns2(ew,ns),  &
                         d2thckdewdns(ew,ns))

        do up = 1, upn

           lateralboundry = .true.

            shift = indshift(  1, ew, ns, up,                  &
                               ewn, nsn, upn,                  &
                               loc_array,                      & 
                               stagthck(ew-1:ew+1,ns-1:ns+1) )

            call bodyset(ew,  ns,  up,        &
                         ewn, nsn, upn,       &
                         dew,      dns,       &
                         pt,       loc_array, &
                         loc,      stagthck,  &
                         thisdusrfdx,         &
                         dusrfdew, dusrfdns,  &
                         dlsrfdew, dlsrfdns,  &
                         efvs(up-1+shift(1):up+shift(1),ew:ew+1,ns:ns+1),  &
                         othervel(up-1+shift(1):up+1+shift(1),  &
                         ew-1+shift(2):ew+1+shift(2),  &
                         ns-1+shift(3):ns+1+shift(3)), &
                         betasquared(ew,ns), abar=flwabar )        
        end do
        lateralboundry = .false.

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!     else 
!     else if (boundarys == mask(ew,ns) .and. shelfmask(ew,ns) == 0 )then 
!     else if (boundarys == mask(ew,ns) )then 

!    !*sp* must be ACTIVE for test case OR commented OUT for ice shelf test case 
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!        !*sp* puts zeros on rhs, coincident w/ location of the additional equation for the HO sfc and basal bcs
!        locplusup = loc(1)
!        call valueset(0.0_dp)
!        locplusup = loc(1) + upn + 1
!        call valueset(0.0_dp)
!        do up = 1,upn
!           locplusup = loc(1) + up
!           call valueset( thisvel(up,ew,ns) )     ! *sp* vel at margin set to 0
!        end do

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!     end if

!   !*sp* must be ACTIVE for test cases  
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    end if
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    end do;     ! ew 
  end do        ! ns

  return

end subroutine findcoefstr

!***********************************************************************

subroutine bodyset(ew,  ns,  up,           &
                   ewn, nsn, upn,          &
                   dew,      dns,          &
                   pt,       loc_array,    &
                   loc,      stagthck,     &
                   thisdusrfdx,            &
                   dusrfdew, dusrfdns,     &
                   dlsrfdew, dlsrfdns,     &
                   local_efvs,             &
                   local_othervel,         &
                   betasquared,            &
                   local_thisvel,          &
                   abar)

  implicit none

!whl - added several in/out arguments here

  integer, intent(in) :: ewn, nsn, upn
  integer, intent(in) :: ew, ns, up
  real (kind = dp), intent(in) :: dew, dns
  integer, intent(in) :: pt
  integer, dimension(ewn-1,nsn-1), intent(in) :: loc_array
  integer, dimension(6), intent(in) :: loc

  real (kind = dp), dimension(:,:), intent(in) :: stagthck
  real (kind = dp), dimension(:,:), intent(in) :: dusrfdew, dusrfdns
  real (kind = dp), dimension(:,:), intent(in) :: dlsrfdew, dlsrfdns
  real (kind = dp), dimension(:,:), intent(in) :: thisdusrfdx

  real (kind = dp), dimension(2,2,2), intent(in) :: local_efvs
  real (kind = dp), dimension(3,3,3), intent(in) :: local_othervel
! *sp* this vel compoenent is taken from the solution 
! at the end of the previous iteration, and is treated 
! as a known value (see below, where all rhs coef. are
! multiplied by this vel.)

  real (kind = dp), intent(in) :: betasquared
  real (kind = dp), intent(in), optional :: local_thisvel
  real (kind = dp), intent(in), optional :: abar

! *sp* needed for dirichlet basal bc (see note below)
  real (kind = dp), dimension(3,3,3) :: g

! *sp* source term on rhs of lateral boundary condition, 
! e.g. source = rho*g*H/(2*Neff) * ( 1 - rho_i / rho_w ) for ice shelf 
  real (kind = dp) :: source, slopex, slopey

! *sp* lateral boundary normal and vector to indicate use of forward 
! or bacward one-sided diff. for lateral bcs
  real (kind = dp), dimension(2) :: fwdorbwd, normal

  integer, dimension(2) :: bcflag  ! *sp* indicates choice of sfc and basal bcs ...

  locplusup = loc(1) + up       

  if( lateralboundry )then

  ! *********************************************************************************************
  ! lateral boundary conditions 
  
  ! if at sfc or bed, source due to seawater pressure is 0 and bc normal vector
  !  should contain components (ds/dx, ds/dy) or (db/dx, db,dy)
     source = 0.0_dp

     call getlatboundinfo( ew,  ns,  up,                            &
                           ewn, nsn, upn,                           & 
                           stagthck(ew-1:ew+1, ns-1:ns+1),          &
                           loc_array, fwdorbwd, normal, loc_latbc)
                
    if( up == 1 .or. up == upn )then
    
        if( up == 1 )then
           locplusup = loc(1) + up - 1  ! reverse the sparse matrix / rhs vector row index by 1 ...
           slopex = dusrfdew(ew,ns); slopey = dusrfdns(ew,ns)
        else
           locplusup = loc(1) + up + 1  ! advance the sparse matrix / rhs row vector index by 1 ...
           slopex = dlsrfdew(ew,ns); slopey = dlsrfdns(ew,ns)
        end if

        g = normhorizmainbc_lat(dew,           dns,             &
                                slopex,        slopey,          &
                                dsigmadew(up), dsigmadns(up),   &
                                pt,            2,               &
                                dup(up),                        &
                                oneorfour,     fourorone,       &
                                onesideddiff,                   &
                                normal,        fwdorbwd)

        ! add on coeff. associated w/ du/digma  
        g(:,3,3) = g(:,3,3) & 
                 + vertimainbc( stagthck(ew,ns),    bcflag,dup(up), &
                                local_efvs,         betasquared )

        ! put the coeff. for the b.c. equation in the same place as the prev. equation
        ! (w.r.t. cols), on a new row ...

        call fillsprsebndy( g, locplusup, loc_latbc, up, normal ) 

        rhsd(locplusup) = sum( croshorizmainbc_lat(dew,           dns,           &
                                                   slopex,        slopey,        &
                                                   dsigmadew(up), dsigmadns(up), &
                                                   pt,            2,             &
                                                   dup(up),       local_othervel,&
                                                   oneortwo,      twoorone,      &
                                                   onesideddiff,                 &
                                                   normal,fwdorbwd)              &         
                                                 * local_othervel )
                
    end if     ! up = 1 or up = upn (IF at lateral boundary and IF at surface or bed)

    ! if in main body and at ice/ocean boundary, calculate depth-averaged stress
    ! due to sea water, bc normal vector components should be boundary normal 
    locplusup = loc(1) + up
    slopex = normal(1)
    slopey = normal(2)

    ! Two options here, (1) use the 1d solution that involves the rate factor, 
    !                   (2) use the more general solution that involves the eff. visc.
    ! ... only one of these options should be active (comment the other lines out)

    ! (1) 
    ! source term (strain rate at shelf/ocean boundary) from Weertman's analytical solution 
    ! (eq. 2, Pattyn+, 2006, JGR v.111; eq. 8, Vieli&Payne, 2005, JGR v.110). Note that this 
    ! contains the 1d assumption that ice is not spreading lateraly !(assumes dv/dy = 0 for u along flow)
    ! Note that factor of 2 in front of 'stagthck' is NOT part of the standard bc. Here, it is used to 
    ! correct for the fact that the staggered thickness will be 1/2 of the normal thickness at a boundary 
!    source = abar*vis0 * ( 1.0_dp/4.0_dp * rhoi * grav * 2.0d0 * stagthck(ew,ns)*thk0 * ( 1.0_dp - rhoi/rhoo))**3.0_dp
!
!    ! multiply by 4 so that case where v=0, du/dy = 0, LHS gives: du/dx = du/dx|_shelf 
!    ! (i.e. LHS = 4*du/dx, requires 4*du/dx_shelf)
!    source = source * 4.0_dp
!
!    ! split source based on the boundary normal orientation and non-dimensinoalize
!    if( normal(1) .ne. 0.0d0 .and. normal(2) .ne. 0.0d0 )then
!    ! *sp* this necessary so that splitting of scalar source term on RHS is 1/2 to y-dir and 1/2 to x-dir 
!    ! (rather than 1/sqrt(2) to each dir)
!    ! Also note that this is not really appropriate to apply to 2d flow, since terms other than du/dx in 
!    ! eff. strain rate are ignored. For 2d flow, should use option (2) below. 
!       source = source * normal(pt) * 1 / sqrt( 2.0d0 ) * tim0
!    else
!       source = source * normal(pt) * tim0
!    end if


    ! (2)
    ! source term (strain rate at shelf/ocean boundary) from MacAyeal depth-ave solution. 
    ! As above, factor of 2 in front of 'stagthck' is not part of the formal solution but is used here to
    ! correct for the fact that the boundary thickness on the staggered grid will generally be ~1/2 of the 
    ! full thickness at the boundary (as a result of averaging to make 'stagthck' from 'thck'). 
    source = rhoi * grav * 2.0d0 * stagthck(ew,ns) * thk0 / 2.0_dp * ( 1.0_dp - rhoi / rhoo )

    ! terms after "/" below count number of non-zero efvs cells ... needed for averaging efvs at boundary 
    source = source / ( evs0 * sum(local_efvs, local_efvs .gt. 1.0d-10) / &
             sum( local_efvs/local_efvs,local_efvs .gt. 1.0d-10 ) )

    if( normal(1) .ne. 0.0d0 .and. normal(2) .ne. 0.0d0 )then
    ! *sp* this necessary so that splitting of scalar source term on RHS is 1/2 to y-dir and 1/2 to x-dir 
    ! (rather than 1/sqrt(2) to each dir)
       source = source * normal(pt) * 1 / sqrt( 2.0d0 ) * tim0
    else
       source = source * normal(pt) * tim0 ! non-dim
    end if
                        
    g = normhorizmainbc_lat(dew,           dns,        &
                            slopex,        slopey,     &
                            dsigmadew(up), dsigmadns(up),  &
                            pt,            1,          &
                            dup(up),                   &
                            oneorfour,     fourorone,  &
                            onesideddiff,              &
                            normal,        fwdorbwd)

    ! put the coeff. for the b.c. equation in the same place as the prev. equation
    ! (w.r.t. cols), on a new row ...
    call fillsprsebndy( g, locplusup, loc_latbc, up, normal )

    rhsd(locplusup) = sum( croshorizmainbc_lat(dew,           dns,            &
                                               slopex,        slopey,         &
                                               dsigmadew(up), dsigmadns(up),  &
                                               pt,            1,              &
                                               dup(up),       local_othervel, &
                                               oneortwo,      twoorone,       &
                                               onesideddiff,                  &
                                               normal,        fwdorbwd)       &
                                              * local_othervel ) + source
 
  else   ! NOT at a lateral boundary 

! *********************************************************************************************
! normal discretization for points inside of lateral boundary and inside main body of ice sheet
        
     g = normhorizmain(pt,up,local_efvs)
     g(:,2,2) = g(:,2,2) + vertimain(hsum(local_efvs),up)             
     call fillsprsemain(g,locplusup,loc,up) 
     rhsd(locplusup) = thisdusrfdx(ew,ns) - sum(croshorizmain(pt,up,local_efvs) * local_othervel)

  end if

! *********************************************************************************************
! higher-order sfc and bed boundary conditions in main body of ice sheet (NOT at lat. boundry)

  if(  ( up == upn  .or. up == 1 ) .and. .not. lateralboundry) then

     if( up == 1 )then                ! specify necessary variables for free sfc
        bcflag = (/1,0/)
        locplusup = loc(1) + up - 1   ! reverse the sparse matrix / rhs vector row index by 1 ...
        slopex = dusrfdew(ew,ns); slopey = dusrfdns(ew,ns)
     else                             ! specify necessary variables for basal bc
        !bcflag = (/0,0/)             ! flag for u=v=0 at bed; doesn't work well;
                                      ! better to specify very large value for betasquared below
        bcflag = (/1,1/)              ! flag for specififed stress at bed: Tau_zx = betasquared * u_bed,
                                      ! where betasquared is MacAyeal-type traction parameter
        locplusup = loc(1) + up + 1   ! advance the sparse matrix / rhs row vector index by 1 ...
        slopex = dlsrfdew(ew,ns); slopey = dlsrfdns(ew,ns)
     end if

     g = normhorizmainbc(dew,           dns,     &
                         slopex,        slopey,  &
                         dsigmadew(up), dsigmadns(up),  &
                         pt,            bcflag,  &
                         dup(up),                &
                         oneorfour,     fourorone)

     ! add on coeff. associated w/ du/digma
     g(:,2,2) = g(:,2,2)   &
              + vertimainbc(stagthck(ew,ns),bcflag,dup(up),local_efvs,betasquared)

     ! put the coeff. for the b.c. equation in the same place as the prev. equation
     ! (w.r.t. cols), on a new row ...
     call fillsprsemain(g,locplusup,loc,up)

     rhsd(locplusup) = sum( croshorizmainbc(dew,           dns,            &
                                            slopex,        slopey,         &
                                            dsigmadew(up), dsigmadns(up),  &
                                            pt,            bcflag,         &
                                            dup(up),       local_othervel, &
                                            oneortwo,      twoorone)       &
                                            * local_othervel )

  end if   ! (up = 1 or up = upn) and lateralboundry = F

! *********************************************************************************************

  return

end subroutine bodyset

!***********************************************************************

subroutine valueset(local_value)

! *sp* plugs values into the rhs vector for eventual soln. to Ax=rhs ...

  implicit none

  real (kind = dp), intent(in) :: local_value

  call putpcgc(1.0_dp,locplusup,locplusup)
  rhsd(locplusup) = local_value 

  return

end subroutine valueset

!***********************************************************************

subroutine calccoeffsinit (upn, dew, dns)

! *sp* determines constants used in various FD calculations associated with 'findcoefst'

  implicit none

  integer, intent(in) :: upn
  real (kind = dp), intent(in) :: dew, dns

! these are coefficients used in finite differences of vertical terms.

! *sp* the factor of 1/4 here may actually be a factor of 1/2 ... 
! Is it from the use of function 'hsum' to sum visc. in horiz. direction? 
! If so, this would then require dividing by 4 to get a mean (in map view, 
! each vel. point is surrounded by 4 visc. points)

  cvert(:,1) = (len0**2) / (4.0_dp * thk0**2 * dup**2)

! ... used for specified traction basal bc (OLD VERSION)
  cvert(:,2) = (len0) / (16.0_dp * thk0 * dup)


! these are coefficients used in finite differences of horizontal terms
! for d/dx(fdu/dx), d/dx(fdu/dy), d/dsigma(fdu/dx), d/dx(fdu/dsigma) and
! du/dsigma.  in some cases need separate coeffs for ew and ns dimensions.

! *sp* the following 4 terms are diff. then those in the write-up on p.44 
! by a factor of 1/2 (however, they are the same as those on p.34)

  cdxdx = (/ 0.25_dp / dew**2, 0.25_dp / dns**2 /)
  cdxdy = 0.0625_dp / (dew * dns)
  cdsdx(:,1) = 0.0625_dp / (dew * dup); cdsdx(:,2) = 0.0625_dp / (dns * dup);

! *sp* this term is diff. by a factor of 1/8 from that on p.44
  cdsds = 0.25_dp / (dup * dup)

! *sp* this term is diff. by a factor of 1/8 from that on p.44 (but is equal 
! to the def. given on p.34)	
  cds = 0.0625_dp / dup

  return
     
end subroutine calccoeffsinit

!***********************************************************************

subroutine calccoeffs(upn,        sigma,                    &
                      stagthck,                             &
                      dusrfdew,   dusrfdns,                 &
                      dthckdew,   dthckdns,                 &
                      d2usrfdew2, d2usrfdns2, d2usrfdewdns, &
                      d2thckdew2, d2thckdns2, d2thckdewdns)

! *sp* SUBROUTINE called from 'findcoefst' to find coefficients in 
!      stress balance equations
! Detemines coeficients needed for finite differencing.
! This is a column-based operation.  
! In general these coefficients refer to grid transformations and averaging 
!  of efvs to half grid points.

  implicit none

  integer, intent(in) :: upn 
  real (kind = dp), dimension(:), intent(in) :: sigma
  real (kind = dp), intent(in) :: stagthck, dusrfdew, dusrfdns, dthckdew, dthckdns, &
                                  d2usrfdew2, d2usrfdns2, d2usrfdewdns, &
                                  d2thckdew2, d2thckdns2, d2thckdewdns
     
! *sp* these next 3 values are not used in this subroutine ... 
! Are they saved to the module 'strscals' after being used here, 
!  or are they not used at all?

  fvert(:,1) = cvert(:,1) / stagthck**2     ! *sp* note that this is effected by the def. of 'cvert', which may be off by factor of 1/2

  fvert(:,2) = cvert(:,2) / stagthck        ! *sp* NOTE that these values may not be needed anymore ... only accessed by old bc functions?
  fvert(:,3) = 2.0_dp * fvert(:,1)          ! *sp* this also affected by def. of 'cvert' through def. of 'fvert(1)'

  dsigmadew = calcdsigmadx(upn, sigma, dusrfdew, dthckdew, stagthck)
  dsigmadns = calcdsigmadx(upn, sigma, dusrfdns, dthckdns, stagthck)

  d2sigmadew2 = calcd2sigmadxdy(upn,        sigma,      &
                                d2usrfdew2, d2thckdew2, &
                                dusrfdew,   dusrfdew,   &
                                dthckdew,   dthckdew,   &
                                stagthck)

  d2sigmadns2 = calcd2sigmadxdy(upn,        sigma,      &
                                d2usrfdns2, d2thckdns2, &
                                dusrfdns,   dusrfdns,   &
                                dthckdns,   dthckdns,   &
                                stagthck)

  d2sigmadewdns = calcd2sigmadxdy(upn,          sigma,         &
                                  d2usrfdewdns, d2thckdewdns,  &
                                  dusrfdew,     dusrfdns,      &
                                  dthckdew,     dthckdns,      &
                                  stagthck)

  d2sigmadewdsigma = calcd2sigmadxdsigma(dthckdew,stagthck)
  d2sigmadnsdsigma = calcd2sigmadxdsigma(dthckdns,stagthck)

  return

end subroutine calccoeffs

!***********************************************************************

function calcdsigmadx(upn,     sigma,    &
                      dusrfdx, dthckdx,  &
                      stagthck)

  implicit none

  integer, intent(in) :: upn  
  real (kind = dp), dimension(:), intent(in) :: sigma
  real (kind = dp), intent(in) :: stagthck, dusrfdx, dthckdx
  real (kind = dp), dimension(upn) :: calcdsigmadx

  calcdsigmadx = (dusrfdx - sigma * dthckdx) / stagthck

  return

end function calcdsigmadx

!***********************************************************************

function calcd2sigmadxdy(upn,        sigma,       &
                         d2usrfdxdy, d2thckdxdy,  &
                         dusrfdx,    dusrfdy,     &
                         dthckdx,    dthckdy,     &
                         stagthck)

  implicit none

  integer, intent(in) :: upn 
  real (kind = dp), dimension(:), intent(in) :: sigma
  real (kind = dp), intent(in) :: d2usrfdxdy, d2thckdxdy, dusrfdx, dusrfdy, &
                                  dthckdx, dthckdy, stagthck 
  real (kind = dp), dimension(upn) :: calcd2sigmadxdy

  calcd2sigmadxdy = (stagthck * d2usrfdxdy - &
                     dusrfdx * dthckdy - dusrfdy * dthckdx + &
                     sigma * (2.0_dp * dthckdx * dthckdy - &
                     stagthck * d2thckdxdy)) / stagthck**2

  return

end function calcd2sigmadxdy

!***********************************************************************

function calcd2sigmadxdsigma(dthckdx,stagthck)

  implicit none

  real (kind = dp), intent(in) :: dthckdx, stagthck 
  real (kind = dp) :: calcd2sigmadxdsigma

  calcd2sigmadxdsigma = - dthckdx / stagthck

  return

end function calcd2sigmadxdsigma

!***********************************************************************

function vertimain(efvs,up)

! *sp* function to come up with coeff. that correspond to the 'vertical' terms 
! on the LHS of the standard equation: (X/H)^2 * d/dz( du/dz) 
 
  implicit none

  real (kind = dp), dimension(2), intent(in) :: efvs

  real (kind = dp), dimension(3) :: vertimain

  integer, intent(in) :: up

  vertimain(3) = fvert(up,1) * efvs(2)  ! *sp* coeffs. for standard 2nd-order, centered diff.
  vertimain(1) = fvert(up,1) * efvs(1)
  vertimain(2) = - vertimain(3) - vertimain(1)

  return

end function vertimain

!***********************************************************************

function normhorizmain(which,up,efvs)

! *sp* FUNCTION called from 'findcoefst' to calculate normal-stress grad terms 
!      like: d/dx(f(du/dx)), d/dy(f(dv/dy)), etc.  
! ... calls FUNCTIONS: horiztermdxdx, horiztermdsdx, horiztermdxds,
!                      horiztermdsds, horiztermds 
! determines coefficients from d/dx(fdu/dx) and d/dy(fdu/dy)

  implicit none

  integer, intent(in) :: which, up
  real (kind = dp), dimension(:,:,:), intent(in) :: efvs

  real (kind = dp), dimension(3,3,3) :: normhorizmain
  real (kind = dp), dimension(3,3,3) :: g, h
  real (kind = dp), dimension(2) :: sumefvsup, sumefvsew, sumefvsns
  real (kind = dp) :: sumefvs

  g = 0.0_dp
  h = 0.0_dp

  sumefvsup = hsum(efvs)
  sumefvsew = sum(sum(efvs,3),1)
  sumefvsns = sum(sum(efvs,2),1)
  sumefvs = sum(efvs)

! *sp* here, coeff. for all the norm. horiz. terms are summed to come up with coeff. at 
! the center of a 3x3x3 block ... 

! for d(f.du/dx)/dx
   
  g(2,:,2) = horiztermdxdx(sumefvsew,cdxdx(1))
  g(:,1:3:2,2) = g(:,1:3:2,2) + horiztermdsdx(dsigmadew(up),sumefvsup,cdsdx(up,1)) 
  g(1:3:2,:,2) = g(1:3:2,:,2) + horiztermdxds(dsigmadew(up),sumefvsew,cdsdx(up,1)) 
  g(:,2,2) = g(:,2,2) + horiztermdsds(dsigmadew(up)**2,sumefvsup,cdsds(up)) 
  g(1:3:2,2,2) = g(1:3:2,2,2) + horiztermds(d2sigmadew2(up)+d2sigmadewdsigma*dsigmadew(up),sumefvs,cds(up))

! for d(f.du/dy)/dy 

  h(2,2,:) = horiztermdxdx(sumefvsns,cdxdx(2))
  h(:,2,1:3:2) = h(:,2,1:3:2) + horiztermdsdx(dsigmadns(up),sumefvsup,cdsdx(up,2)) 
  h(1:3:2,2,:) = h(1:3:2,2,:) + horiztermdxds(dsigmadns(up),sumefvsns,cdsdx(up,2)) 
  h(:,2,2) = h(:,2,2) + horiztermdsds(dsigmadns(up)**2,sumefvsup,cdsds(up))  
  h(1:3:2,2,2) = h(1:3:2,2,2) + horiztermds(d2sigmadns2(up)+d2sigmadnsdsigma*dsigmadns(up),sumefvs,cds(up))

  normhorizmain = g * fourorone(which) + h * oneorfour(which)

  return

end function normhorizmain

!***********************************************************************
   
function croshorizmain(which,up,efvs)
! *sp* FUNCTION called from 'findcoefst' to calculate cross-stress grad terms 
!      like: d/dx(f(du/dy)), d/dy(f(dv/dx)), etc.  
! ... calls FUNCTIONS: horiztermdxdy, horiztermdsdx, horiztermdxds, 
!                      horiztermdsds, horiztermds 
! determines coefficients from d/dx(fdu/dy) and d/dy(fdu/dx)

  implicit none

  integer, intent(in) :: which, up
  real (kind = dp), dimension(:,:,:), intent(in) :: efvs

  real (kind = dp), dimension(3,3,3) :: croshorizmain
  real (kind = dp), dimension(3,3,3) :: g = 0.0_dp, h = 0.0_dp
  real (kind = dp), dimension(2) :: sumefvsup, sumefvsew, sumefvsns
  real (kind = dp) :: sumefvs

  g = 0.0_dp
  h = 0.0_dp

  sumefvsup = hsum(efvs)
  sumefvsew = sum(sum(efvs,3),1)
  sumefvsns = sum(sum(efvs,2),1)
  sumefvs = sum(efvs)

! *sp* here, coeff. for all the cross horiz. terms are summed to come up with coeff. at 
! the center of a 3x3x3 block ... 

! for d(f.du/dy)/dx

  g(2,:,1:3:2) = horiztermdxdy(sumefvsew,cdxdy)
  g(:,2,1:3:2) = g(:,2,1:3:2) + horiztermdsdx(dsigmadew(up),sumefvsup,cdsdx(up,2))
  g(1:3:2,:,2) = g(1:3:2,:,2) + horiztermdxds(dsigmadns(up),sumefvsew,cdsdx(up,1))
  g(:,2,2) = g(:,2,2) + horiztermdsds(dsigmadew(up)*dsigmadns(up),sumefvsup,cdsds(up))   
  g(1:3:2,2,2) = g(1:3:2,2,2) + horiztermds(d2sigmadewdns(up)+d2sigmadnsdsigma*dsigmadew(up),sumefvs,cds(up))

! for d(f.du/dx)/dy 

  h(2,1:3:2,:) = transpose(horiztermdxdy(sumefvsns,cdxdy))
  h(:,1:3:2,2) = h(:,1:3:2,2) + horiztermdsdx(dsigmadns(up),sumefvsup,cdsdx(up,1))
  h(1:3:2,2,:) = h(1:3:2,2,:) + horiztermdxds(dsigmadew(up),sumefvsns,cdsdx(up,2))
  h(:,2,2) = h(:,2,2) + horiztermdsds(dsigmadew(up)*dsigmadns(up),sumefvsup,cdsds(up))
  h(1:3:2,2,2) = h(1:3:2,2,2) + horiztermds(d2sigmadewdns(up)+d2sigmadewdsigma*dsigmadns(up),sumefvs,cds(up))

  croshorizmain = g * twoorone(which) + h * oneortwo(which)

  return

end function croshorizmain

!***********************************************************************

! ***************************************************************************
! *sp* functions to deal with higher-order boundary conditions at sfc and bed
! ***************************************************************************

function vertimainbc(thck,bcflag,dup,efvs,betasquared)

! altered form of 'vertimain' that calculates coefficients for higher-order
! b.c. that go with the 'normhorizmain' term: -(X/H)^2 * dsigma/dzhat * du/dsigma 
   
    implicit none

	real (kind = dp), intent(in) :: dup, thck, betasquared
	real (kind = dp), intent(in), dimension(2,2,2) :: efvs
    integer, intent(in), dimension(2) :: bcflag

	real (kind = dp) :: c
	real (kind = dp), dimension(3) :: vertimainbc

    c = 0.0_dp

    ! for higher-order FREE SURFACE B.C. for x ('which'=1) or y ('which'=2) direction ...
    if( bcflag(1) == 1 )then

           c = - 1 / thck / (2*dup) * (len0**2 / thk0**2)   ! value of coefficient

           vertimainbc(:) = 0.0_dp
           vertimainbc(3) = -c 
           vertimainbc(1) = c
           vertimainbc(2) = vertimainbc(3) + vertimainbc(1) ! should = 0


    ! for higher-order BASAL B.C. w/ specified basal traction, add on the necessary source term ...
    if( bcflag(2) == 1 )then

            ! last set of terms is mean visc. of ice nearest to the bed
            vertimainbc(2) = vertimainbc(2)   &
                           + ( betasquared / ( sum( efvs(2,:,:) ) / 4.0_dp ) ) * (len0 / thk0)

    end if


    ! for higher-order BASAL B.C. U=V=0, in x ('which'=1) or y ('which'=2) direction ...
    else if( bcflag(1) == 0 )then

           ! if u,v set to 0, there are no coeff. assoc. with du/digma terms ...
           vertimainbc(:) = 0.0_dp

    end if

    return

end function vertimainbc

!***********************************************************************

function normhorizmainbc(dew,       dns,        &
                         dusrfdew,  dusrfdns,   &
                         dsigmadew, dsigmadns,  &
                         which,     bcflag,     &
                         dup,                   &
                         oneorfour, fourorone)

! determines higher-order surface and basal boundary conditions for LHS of equation ...
! ... that is, gives 3x3x3 coeff. array for either u or v component of velocity,
! depending on the value of the flag 'which'. Example of function call:
!
!  g = normhorizmainbc(dusrfew(ew,ns),dusrfnx(ew,ns),dsigmadew(ew,ns,up),dsigmadns(ew,ns,up),which,up,bcflag)	
!
! ... where g is a 3x3x3 array.
!
! 'dusrfdns' and 'dusrfdew' are n x m arrays, where n and m are the number of grid cells in the
! horiz x and y directions.
!
! 'dsigmadew' and 'dsigmadew' are vectors w/ dimensions of 'upn', the number of grid points in 
! the vertical.
!
! 'bcflag' is a 1 x 2 vector to indicate (1) which b.c. is being solved for (surface or bed) and 
! (2), if solving for the bed b.c., which type of b.c. to use (for now, only allows for u=v=0). 
! For example, bcflag = [ 0, 0 ] denotes free sfc bc; bcflag = [ 1, 0 ] denotes basal bc w/ u=v=0;
! bcflag = [ 1, 1 ] denotes basal bc w/ some other choice of basal b.c., etc...
!
! The vectors "fourorone" and "oneorfour" are given by: fourorone = [ 4 1 ]; oneorfour = [ 1 4 ].
! A single value is chosen from each vector and applied to the calculation of coefficients below.
! The "correct" value needed to satisfy the expression is chosen based on the "which" flag, which
! takes on a value of 1 for calculations in the x direction and a value of 2 for calculations in 
! the y direction. 

    implicit none

    real (kind = dp), intent(in) :: dew, dns
	real (kind = dp), intent(in) :: dusrfdew, dusrfdns, dsigmadew, dsigmadns, dup
	real (kind = dp), intent(in), dimension(2) :: oneorfour, fourorone
	real (kind = dp), dimension(3,3,3) :: normhorizmainbc
	real (kind = dp), dimension(3,3,3) :: g
	real (kind = dp) :: c

    integer, intent(in) :: which
    integer, intent(in), dimension(2) :: bcflag

    c = 0.0_dp
    g(:,:,:) = 0.0_dp

    ! for higher-order FREE SURFACE B.C. for x ('which'=1) or y ('which'=2) direction ...
    if( bcflag(1) == 1 )then

           ! first, coeff. that go with du/dsigma, and thus are associated 
           ! with u(1,2,2) and u(3,2,2) ...

           c = ( fourorone(which) * dusrfdew * dsigmadew   &
               + oneorfour(which) * dusrfdns * dsigmadns )/(2*dup)
           g(3,2,2) = -c 
           g(1,2,2) = c

           ! next, coeff. that go with du/dxhat and du/dyhat terms ...
           c = fourorone(which) * dusrfdew / (2*dew)
           g(2,3,2) = c
           g(2,1,2) = -c

           c = oneorfour(which) * dusrfdns / (2*dns)
           g(2,2,3) = c
           g(2,2,1) = -c

    ! for higher-order BASAL B.C. U=V=0, in x ('which'=1) or y ('which'=2) direction ...
    ! note that this requires that rhs(up) be set to 0 as well ...
    else if( bcflag(1) == 0 )then

           g(:,:,:) = 0.0_dp
           g(2,2,2) = 1.0_dp;

    end if

    normhorizmainbc = g

    return

end function normhorizmainbc

!***********************************************************************

function croshorizmainbc(dew,       dns,       &
                         dusrfdew,  dusrfdns,  &
                         dsigmadew, dsigmadns, &
                         which,     bcflag,    &
                         dup,       local_othervel,  &
                         oneortwo,  twoorone)

! The vectors "twoorone" and "oneortwo" are given by: twoorone = [ 2 1 ]; oneortwo = [ 1 2 ];
! The entire funciton is analagous to "normhorizmainbc" ... 
!  see additional comments there for more description.

    implicit none

    integer, intent(in) :: which
    integer, intent(in), dimension(2) :: bcflag

    real (kind = dp), intent(in) :: dew, dns
	real (kind = dp), intent(in), dimension(2) :: oneortwo, twoorone
	real (kind = dp), intent(in) :: dusrfdew, dusrfdns, dsigmadew, dsigmadns, dup
	real (kind = dp), intent(in), dimension(3,3,3) :: local_othervel
	real (kind = dp), dimension(3,3,3) :: g, croshorizmainbc
	real (kind = dp) :: c

    c = 0.0_dp
    g(:,:,:) = 0.0_dp

    ! for higher-order FREE SURFACE B.C. for x ('which'=1) or y ('which'=2) direction ...
    if( bcflag(1) == 1 )then

           ! first, coeff. that go with du/dsigma, and thus are associated
           ! with u(1,2,2) and u(3,2,2) ...

           c = ( - twoorone(which) * dusrfdew * dsigmadns   &
                 - oneortwo(which) * dusrfdns * dsigmadew )/(2*dup)
           g(3,2,2) = -c
           g(1,2,2) = c

           ! next, coeff. that go with du/dxhat and du/dyhat terms ...
           c = - oneortwo(which) * dusrfdns / (2*dew)
           g(2,3,2) = c
           g(2,1,2) = -c

           c = - twoorone(which) * dusrfdew / (2*dns)
           g(2,2,3) = c
           g(2,2,1) = -c

    ! for higher-order BASAL B.C. U=V=0, in x ('which'=1) or y ('which'=2) direction ...
    else if( bcflag(1) == 0 )then

       g(:,:,:) = 0.0_dp
       g(2,2,2) = 0.0_dp   ! specify the value of the boundary velocity in the rhs HERE

       ! this forces the multiplication by 'local_otherval' in the main program 
       ! to result in a value of 1, thus leaving the boundary vel. unchanged
       g = g / local_othervel  

    end if

    croshorizmainbc = g

    return

end function croshorizmainbc

!***********************************************************************

function normhorizmainbc_lat(dew,       dns,   &
                             dusrfdew,  dusrfdns,  &
                             dsigmadew, dsigmadns, &
                             which,     what,      &
                             dup,                  &
                             oneorfour, fourorone, &
                             onesideddiff,         &
                             normal,    fwdorbwd)

    implicit none

    real (kind = dp), intent(in) :: dew, dns
	real (kind = dp), intent(in) :: dusrfdew, dusrfdns, dsigmadew, dsigmadns, dup
	real (kind = dp), intent(in), dimension(2) :: oneorfour, fourorone, normal, fwdorbwd
	real (kind = dp), intent(in), dimension(3) :: onesideddiff

	real (kind = dp), dimension(3,3,3) :: normhorizmainbc_lat
	real (kind = dp), dimension(3,3,3) :: g
	real (kind = dp), dimension(2) :: whichbc
	real (kind = dp) :: c

    integer, intent(in) :: which, what

    c = 0.0_dp;	g(:,:,:) = 0.0_dp; whichbc = (/ 0.0_dp, 1.0_dp /)

    ! for higher-order FREE SURFACE B.C. for x ('which'=1) or y ('which'=2) direction ...

    ! first, coeff. that go with du/dsigma, and thus are associated with u(1,2,2) and u(3,2,2) ...
    ! ...note that these are stored in an empty column of 'g' (a corner column) so that we don't 
    ! overwrite these values in the case of fwd/bwd horiz. diffs., which require 3 spaces
    c = ( fourorone(which) * dusrfdew * dsigmadew    &
            + oneorfour(which) * dusrfdns * dsigmadns )/(2*dup)
    g(3,3,3) = -c * whichbc(what) 
    g(1,3,3) = c * whichbc(what)

    if( normal(1) .eq. 0.0_dp )then     ! centered in x ...

           c = fourorone(which) * dusrfdew / (2*dew)
           g(2,3,2) = c
           g(2,1,2) = -c

    elseif( normal(1) .ne. 0.0_dp )then     ! forward/backward in x ...

           c = fourorone(which) * fwdorbwd(1) * onesideddiff(1) * dusrfdew / (2*dew)
           g(2,2-int(fwdorbwd(1)),2) = c

           c = fourorone(which) * fwdorbwd(1) * onesideddiff(2) * dusrfdew / (2*dew)
           g(2,2,2) = c

           c = fourorone(which) * fwdorbwd(1) * onesideddiff(3) * dusrfdew / (2*dew)
           g(2,2+int(fwdorbwd(1)),2) = c

    end if


    if( normal(2) .eq. 0.0_dp ) then   ! centered in y ... 
                                       ! (NOTE that y coeff. are stored in g(1,:,:) )

           c = oneorfour(which) * dusrfdns / (2*dns)
           g(1,2,3) = c * whichbc(what)
           g(1,2,1) = -c * whichbc(what)

    elseif( normal(2) .ne. 0.0_dp) then ! forward/backward in y ...

           c = oneorfour(which) * fwdorbwd(2) * onesideddiff(1) * dusrfdns / (2*dns)
           g(1,2,2-int(fwdorbwd(2))) = c

           c = oneorfour(which) * fwdorbwd(2) * onesideddiff(2) * dusrfdns / (2*dns)
           g(1,2,2) = c

           c = oneorfour(which) * fwdorbwd(2) * onesideddiff(3) * dusrfdns / (2*dns)
           g(1,2,2+int(fwdorbwd(2))) = c

    end if

    normhorizmainbc_lat = g

    return

end function normhorizmainbc_lat

!***********************************************************************

function croshorizmainbc_lat (dew,       dns,       &
                              dusrfdew,  dusrfdns,  &
                              dsigmadew, dsigmadns, &
                              which,     what,      &
                              dup,       local_othervel,  &
                              oneortwo,  twoorone,  &
                              onesideddiff,         &
                              normal,    fwdorbwd)

! The vectors "twoorone" and "oneortwo" are given by: twoorone = [ 2 1 ]; oneortwo = [ 1 2 ];
! The entire funciton is analagous to "normhorizmainbc" ... 
! see additional comments there for more description.

    implicit none

    real (kind = dp), intent(in) :: dew, dns
	real (kind = dp), intent(in), dimension(2) :: oneortwo, twoorone, fwdorbwd, normal
	real (kind = dp), intent(in), dimension(3) :: onesideddiff
	real (kind = dp), intent(in) :: dusrfdew, dusrfdns, dsigmadew, dsigmadns, dup
	real (kind = dp), intent(in), dimension(3,3,3) :: local_othervel

    integer, intent(in) :: which, what

	real (kind = dp), dimension(3,3,3) :: g, croshorizmainbc_lat
	real (kind = dp), dimension(3) :: gvert
	real (kind = dp), dimension(2) :: whichbc
	real (kind = dp) :: c

    integer, dimension(2) :: inormal

    c = 0.0_dp
    g(:,:,:) = 0.0_dp
    gvert = 0.0_dp
    whichbc = (/ 0.0_dp, 1.0_dp /)
    croshorizmainbc_lat = 0.0_dp

    ! first, coeff. that go with du/dsigma, and thus are associated with u(1,2,2) and u(3,2,2) 
    ! ... note that these are stored in a separate vector (to avoid being overwritten if stored in normal 'g')	

    c = ( - twoorone(which) * dusrfdew * dsigmadns   &
              - oneortwo(which) * dusrfdns * dsigmadew )/(2*dup)
    gvert(3) = -c * whichbc(what)
    gvert(1) = c * whichbc(what)

    if( normal(1) .eq. 0.0_dp )then        ! centered in x ...

           c = -oneortwo(which) * dusrfdns / (2*dew)
           g(2,3,2) = c
           g(2,1,2) = -c

    elseif( normal(1) .ne. 0.0_dp )then    ! forward/backward in x ...
                                           ! (NOTE that x coeff. are stored in g(2,:,:) )

           c = -oneortwo(which) * fwdorbwd(1) * onesideddiff(1) * dusrfdns / (2*dew)
           g(2,2-int(fwdorbwd(1)),2) = c

           c = -oneortwo(which) * fwdorbwd(1) * onesideddiff(2) * dusrfdns / (2*dew)
           g(2,2,2) = c

           c = -oneortwo(which) * fwdorbwd(1) * onesideddiff(3) * dusrfdns / (2*dew)
           g(2,2+int(fwdorbwd(1)),2) = c

    end if

    if( normal(2) .eq. 0.0_dp )then    ! centered in y ...
                                       ! (NOTE that y coeff. are stored in g(1,:,:) )

           c = -twoorone(which) * dusrfdew / (2*dns)
           g(1,2,3) = c
           g(1,2,1) = -c

    elseif( normal(2) .ne. 0.0_dp )then ! forward/backward in y ...

           c = -twoorone(which) * fwdorbwd(2) * onesideddiff(1) * dusrfdew / (2*dns)
           g(1,2,2-int(fwdorbwd(2))) = c

           c = -twoorone(which) * fwdorbwd(2) * onesideddiff(2) * dusrfdew / (2*dns)
           g(1,2,2) = c

           c = -twoorone(which) * fwdorbwd(2) * onesideddiff(3) * dusrfdew / (2*dns)
           g(1,2,2+int(fwdorbwd(2))) = c

    end if

    ! Now rearrange position of coefficients in structure 'g' so that they are multiplied by 
    ! the correct velocity component in 'local_othervel' in 'bodyset' ...
    ! ... this can be done by using the boundary normal vector to shift the indices of the rows/columns
    ! in 'g', in the appropriate direction. First, convert the boundary normal to an integer index ...
    inormal(1) = int( normal(1)/abs(normal(1)) )
    inormal(2) = int( normal(2)/abs(normal(2)) )
    if( abs( inormal(1) ) .ne. 1 )then; inormal(1) = 0; end if
    if( abs( inormal(2) ) .ne. 1 )then; inormal(2) = 0; end if

    croshorizmainbc_lat(2,:,2+inormal(2)) = g(2,:,2)    ! move x-coeffs. appropriate amount
    croshorizmainbc_lat(1,2+inormal(1),:) = g(1,2,:)    ! move y-coeffs. appropriate amount

    ! sum coeffs. in same column and collapse so that all coeff. are on level (2,:,:)	
    croshorizmainbc_lat(2,:,:) = croshorizmainbc_lat(2,:,:) + croshorizmainbc_lat(1,:,:)    

    ! set remaining coeff. on this level to to 0 ...
    croshorizmainbc_lat(1,:,:) = 0.0_dp

    ! account for vertical terms stored seperately in 'gvert'
    croshorizmainbc_lat(1,2+inormal(1),2+inormal(2)) = gvert(1) * whichbc(what)
    croshorizmainbc_lat(3,2+inormal(1),2+inormal(2)) = gvert(3) * whichbc(what)

    return

end function croshorizmainbc_lat

!***********************************************************************

! ---> the following routines are for derivatives in the main body 

function horiztermdxdx(efvs,fact)

! this is the d/dx(f.du/dx) and d/dy(f.du/dy) terms

  implicit none

  real (kind = dp), dimension(2), intent(in) :: efvs
  real (kind = dp), intent(in) :: fact

  real (kind = dp), dimension(3) :: horiztermdxdx  

  horiztermdxdx(3) = efvs(2) * fact
  horiztermdxdx(1) = efvs(1) * fact
  horiztermdxdx(2) = - horiztermdxdx(3) - horiztermdxdx(1)

  return

end function horiztermdxdx

!***********************************************************************

function horiztermdxdy(efvs,fact)

! this is the d/dy(f.du/dx) and d/dx(f.du/dy) terms

  implicit none

  real (kind = dp), dimension(2), intent(in) :: efvs
  real (kind = dp), intent(in) :: fact

  real (kind = dp), dimension(3,2) :: horiztermdxdy

  horiztermdxdy(3,2) = efvs(2) * fact 
  horiztermdxdy(2,2) = horiztermdxdy(3,2)
  horiztermdxdy(3,1) = - horiztermdxdy(3,2)
  horiztermdxdy(2,1) = - horiztermdxdy(3,2)

  horiztermdxdy(1,2) = - efvs(1) * fact
  horiztermdxdy(2,2) = horiztermdxdy(2,2) + horiztermdxdy(1,2)
  horiztermdxdy(2,1) = horiztermdxdy(2,1) - horiztermdxdy(1,2)
  horiztermdxdy(1,1) = - horiztermdxdy(1,2)

  return

end function horiztermdxdy

!***********************************************************************

function horiztermdsdx(dsigmadxy,efvs,fact)

! this is the d/ds(f.du/dx) and d/ds(f.du/dy) terms

  implicit none

  real (kind = dp), dimension(2), intent(in) :: efvs
  real (kind = dp), intent(in) :: dsigmadxy, fact

  real (kind = dp), dimension(3,2) :: horiztermdsdx  

  horiztermdsdx(3,2) = dsigmadxy * efvs(2) * fact
  horiztermdsdx(2,2) = horiztermdsdx(3,2)
  horiztermdsdx(3,1) = - horiztermdsdx(3,2)
  horiztermdsdx(2,1) = - horiztermdsdx(3,2)

  horiztermdsdx(1,2) = - dsigmadxy * efvs(1) * fact
  horiztermdsdx(2,2) = horiztermdsdx(2,2) + horiztermdsdx(1,2)
  horiztermdsdx(2,1) = horiztermdsdx(2,1) - horiztermdsdx(1,2)
  horiztermdsdx(1,1) = - horiztermdsdx(1,2)

  return

end function horiztermdsdx

!***********************************************************************

function horiztermdxds(dsigmadxy,efvs,fact)

! this is the d/dx(f.du/ds) and d/dy(f.du/ds) terms

  implicit none

  real (kind = dp), dimension(2), intent(in) :: efvs
  real (kind = dp), intent(in) :: dsigmadxy, fact

  real (kind = dp), dimension(2,3) :: horiztermdxds

  horiztermdxds(2,3) = dsigmadxy * efvs(2) * fact
  horiztermdxds(2,2) = horiztermdxds(2,3)
  horiztermdxds(1,3) = - horiztermdxds(2,3)
  horiztermdxds(1,2) = - horiztermdxds(2,3)

  horiztermdxds(2,1) = - dsigmadxy * efvs(1) * fact
  horiztermdxds(2,2) = horiztermdxds(2,2) + horiztermdxds(2,1)
  horiztermdxds(1,2) = horiztermdxds(1,2) - horiztermdxds(2,1)
  horiztermdxds(1,1) = - horiztermdxds(2,1)

  return

end function horiztermdxds

!***********************************************************************

function horiztermdsds(dsigmadxysq,efvs,fact)

! this is the d/ds(f.du/ds) term

  implicit none

  real (kind = dp), dimension(2), intent(in) :: efvs
  real (kind = dp), intent(in) :: dsigmadxysq, fact

  real (kind = dp), dimension(3) :: horiztermdsds

  horiztermdsds(3) = dsigmadxysq * efvs(2) * fact
  horiztermdsds(1) = dsigmadxysq * efvs(1) * fact

  horiztermdsds(2) = - horiztermdsds(3) - horiztermdsds(1)

  return

end function horiztermdsds

!***********************************************************************

function horiztermds(d2sigmadxy2etc,efvs,fact)

! this is the f.du/ds term

  implicit none

  real (kind = dp), intent(in) :: efvs, d2sigmadxy2etc, fact

  real (kind = dp), dimension(2) :: horiztermds

  horiztermds(2) = d2sigmadxy2etc * efvs * fact
  horiztermds(1) = - horiztermds(2)

  return

end function horiztermds

! ---> end of routines for derivatives in the main body 

!***********************************************************************

subroutine fillsprsemain(inp,locplusup,ptindx,up)

  implicit none

  real (kind = dp), dimension(3,3,3), intent(in):: inp
  integer, intent(in) :: locplusup, up
  integer, dimension(6), intent(in) :: ptindx

! insert entries to these points that are on same level
  call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup)
  call putpcgc(inp(2,3,2),ptindx(2)+up,locplusup)
  call putpcgc(inp(2,1,2),ptindx(3)+up,locplusup)
  call putpcgc(inp(2,2,3),ptindx(4)+up,locplusup)  
  call putpcgc(inp(2,2,1),ptindx(5)+up,locplusup)

! add points for level above (that is, points in the local array with a LARGER first index,
! which correspond to grid points that are CLOSER TO THE BED than at current level)
  call putpcgc(inp(3,2,2),ptindx(1)+up+1,locplusup)
  call putpcgc(inp(3,3,2),ptindx(2)+up+1,locplusup)
  call putpcgc(inp(3,1,2),ptindx(3)+up+1,locplusup)
  call putpcgc(inp(3,2,3),ptindx(4)+up+1,locplusup)  
  call putpcgc(inp(3,2,1),ptindx(5)+up+1,locplusup)

! add points for level below (that is, points in the local array with a SMALLER first index,
! which correspond to grid points that are CLOSER TO THE SURFACE than at current level) 
  call putpcgc(inp(1,2,2),ptindx(1)+up-1,locplusup)
  call putpcgc(inp(1,3,2),ptindx(2)+up-1,locplusup)
  call putpcgc(inp(1,1,2),ptindx(3)+up-1,locplusup)
  call putpcgc(inp(1,2,3),ptindx(4)+up-1,locplusup)  
  call putpcgc(inp(1,2,1),ptindx(5)+up-1,locplusup)

  return

end subroutine fillsprsemain

!***********************************************************************

subroutine fillsprsebndy(inp,locplusup,ptindx,up,normal)

  implicit none

  integer, intent(in) :: locplusup, up
  integer, dimension(6), intent(in) :: ptindx
  real (kind = dp), dimension(3,3,3), intent(in) :: inp
  real (kind = dp), dimension(2), intent(in) :: normal

  ! at points where mixed centered and one-side diffs. would apply

  if( normal(1) == 0.0_dp )then         ! at boundary normal to y, centered diffs in x 
    if( normal(2) == -1.0_dp )then      ! at boundary w/ normal [0,-1]
           call putpcgc(inp(1,3,3),ptindx(5)+up-1,locplusup)
           call putpcgc( inp(2,3,3)+inp(1,2,1),ptindx(5)+up,locplusup)
           call putpcgc(inp(3,3,3),ptindx(5)+up+1,locplusup)
           call putpcgc(inp(1,2,3),ptindx(4)+up,locplusup)
    else                                ! at boundary w/ normal [0,1]
           call putpcgc(inp(1,3,3),ptindx(4)+up-1,locplusup)
           call putpcgc(inp(2,3,3)+inp(1,2,3),ptindx(4)+up,locplusup)
           call putpcgc(inp(3,3,3),ptindx(4)+up+1,locplusup)
           call putpcgc(inp(1,2,1),ptindx(5)+up,locplusup)
    end if
    call putpcgc(inp(1,2,2),ptindx(1)+up,locplusup)
  end if

  if( normal(2) == 0.0_dp )then            ! at boundary normal to x, centered diffs in y 
        if( normal(1) == -1.0_dp )then     ! at boundary w/ normal [-1,0]
           call putpcgc(inp(1,3,3),ptindx(3)+up-1,locplusup)
           call putpcgc( inp(2,3,3)+inp(2,1,2),ptindx(3)+up,locplusup)
           call putpcgc(inp(3,3,3),ptindx(3)+up+1,locplusup)
           call putpcgc(inp(2,3,2),ptindx(2)+up,locplusup)
        else                                 ! at boundary w/ normal [1,0]
           call putpcgc(inp(1,3,3),ptindx(2)+up-1,locplusup)
           call putpcgc( inp(2,3,3)+inp(2,3,2),ptindx(2)+up,locplusup)
           call putpcgc(inp(3,3,3),ptindx(2)+up+1,locplusup)
           call putpcgc(inp(2,1,2),ptindx(3)+up,locplusup)
    end if
    call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup)
  end if

  ! at corners where only one-side diffs. apply
  if( normal(1) .gt. 0.0_dp .and. normal(2) .ne. 0.0_dp )then
    if( normal(2) .gt. 0.0_dp )then      ! corner w/ normal [ 1/sqrt(2), 1/sqrt(2) ]
           call putpcgc(inp(1,3,3),ptindx(2)+up-1,locplusup)
           call putpcgc(inp(3,3,3),ptindx(2)+up+1,locplusup)
           call putpcgc(inp(2,3,3)+inp(2,3,2)+inp(1,2,3),ptindx(2)+up,locplusup)
           call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup)
           call putpcgc(inp(1,2,2),ptindx(6)+up,locplusup)
           call putpcgc(inp(1,2,1),ptindx(5)+up,locplusup)
           call putpcgc(inp(2,1,2),ptindx(3)+up,locplusup)
    else                                 ! corner w/ normal [ 1/sqrt(2), -1/sqrt(2) ]
           call putpcgc(inp(1,3,3),ptindx(2)+up-1,locplusup)
           call putpcgc(inp(3,3,3),ptindx(2)+up+1,locplusup)
           call putpcgc(inp(2,3,3)+inp(1,2,1)+inp(2,3,2),ptindx(2)+up,locplusup)
           call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup)
           call putpcgc(inp(2,1,2),ptindx(3)+up,locplusup)
           call putpcgc(inp(1,2,2),ptindx(6)+up,locplusup)
           call putpcgc(inp(1,2,3),ptindx(4)+up,locplusup)
    end if
  end if

  if( normal(1) .lt. 0.0_dp .and. normal(2) .ne. 0.0_dp )then
    if( normal(2) .gt. 0.0_dp )then       ! corner w/ normal [ -1/sqrt(2), 1/sqrt(2) ]
           call putpcgc(inp(1,3,3),ptindx(3)+up-1,locplusup)
           call putpcgc(inp(3,3,3),ptindx(3)+up+1,locplusup)
           call putpcgc(inp(2,3,3)+inp(1,2,3)+inp(2,1,2),ptindx(3)+up,locplusup)
           call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup)
           call putpcgc(inp(2,3,2),ptindx(2)+up,locplusup)
           call putpcgc(inp(1,2,2),ptindx(6)+up,locplusup)
           call putpcgc(inp(1,2,1),ptindx(5)+up,locplusup)
    else                                  ! corner w/ normal [ -1/sqrt(2), -1/sqrt(2) ]
           call putpcgc(inp(1,3,3),ptindx(3)+up-1,locplusup)
           call putpcgc(inp(3,3,3),ptindx(3)+up+1,locplusup)
           call putpcgc(inp(2,3,3)+inp(2,1,2)+inp(1,2,1),ptindx(3)+up,locplusup)
           call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup)
           call putpcgc(inp(1,2,2),ptindx(6)+up,locplusup)
           call putpcgc(inp(2,3,2),ptindx(2)+up,locplusup)
           call putpcgc(inp(1,2,3),ptindx(4)+up,locplusup)
    end if
  end if

  return

end subroutine fillsprsebndy

!***********************************************************************

subroutine getlatboundinfo( ew, ns, up, ewn, nsn, upn,  &  
                           thck, loc_array,             &
                           fwdorbwd, normal, loc_latbc)
  implicit none

  integer, intent(in) :: ew, ns, up
  integer, intent(in) :: ewn, nsn, upn
  integer, dimension(ewn-1,nsn-1), intent(in) :: loc_array
  real (kind = dp), dimension(3,3), intent(in) :: thck

  real (kind = dp), dimension(2), intent(out) :: fwdorbwd, normal
  integer, dimension(6), intent(out) :: loc_latbc

  real (kind = dp), dimension(3,3) :: mask, maskcorners
  real (kind = dp), dimension(3,3) :: thckmask
  real (kind = dp), dimension(3) :: testvect
  real (kind = dp) :: phi, deg2rad

  deg2rad = 3.141592654d0 / 180.0d0
  loc_latbc = 0; phi = 0
  mask(:,1) = (/ 0.0d0, 180.0d0, 0.0d0 /)
  mask(:,2) = (/ 270.0d0, 0.0d0, 90.0d0 /)
  mask(:,3) = (/ 0.0d0, 360.0d0, 0.0d0 /)
  maskcorners(:,1) = (/ 225.0d0, 0.0d0, 135.0d0 /)
  maskcorners(:,2) = (/ 0.0d0, 0.0d0, 0.0d0 /)
  maskcorners(:,3) = (/ 315.0d0, 0.0d0, 45.0d0 /)

  ! specify new value of 'loc' vector such that fwd/bwd diffs. are set up correctly in sparse matrix
  ! when function 'fillsprsebndy' is called. Also, specify appropriate values for the vectors 'normal'
  ! and 'fwdorbwd', which specify the orientation of the boundary normal and the direction of forward or
  ! backward differencing to be done in the lateral boundary condition functions 'normhorizmainbc_lat'
  ! and 'crosshorizmainbc_lat'

  ! following is algorithm for calculating boundary normal at 45 deg. increments, based on arbitray
  ! boundary shape

  where( thck .ne. 0.0d0 )
        thckmask = 0.0_dp
  elsewhere( thck .eq. 0.0d0 )
        thckmask = 1.0d0
  endwhere

  testvect = sum( thckmask * mask, 1 )

!if( up .eq. 3 )then ! temporary code for debugging
!  do i = 3,1,-1
!  print *, 'thck = ', thck(:,i)
!  end do
!  print *, ' '
!
!  do i = 3,1,-1
!      print *, 'thckmask = ', thckmask(:,i)
!  end do
!  print *, ' '
!
!  print *, 'testvect =  ', testvect
!  print *, ' '
!end if

    ! calculate the angle of the normal in cart. (x,y) system w/ 0 deg. at 12 O'clock, 90 deg. at 3 O'clock, etc.
    if( sum( sum( thckmask, 1 ) ) .eq. 1.0d0 )then
        phi = sum( sum( thckmask * maskcorners, 1 ) )
    else
        if( any( testvect .eq. 360.0d0 ) )then
            if( sum( testvect ) .eq. 450.0d0 )then
                phi = 45.0d0
            elseif( sum( testvect ) .eq. 630.0d0 )then
                phi = 315.0d0
            else
                phi = 0.0d0
            end if
        elseif( all( testvect .ne. 360 ) )then
            phi = sum( testvect ) / sum( testvect/testvect, testvect .ne. 0.0d0 )
        end if
    end if

!  print *, 'up,ew,ns = ', up,ew,ns 
!  print *, 'phi = ', phi
!  print *, ' ' 
!  pause 

    ! define normal vectors and change to loc_array based on this angle
    if( phi .eq. 0.0d0 )then
         loc_latbc(1) = loc_array(ew,ns-1); loc_latbc(4) = loc_array(ew,ns); loc_latbc(5) = loc_array(ew,ns-2)
         loc_latbc(2) = loc_array(ew+1,ns); loc_latbc(3) = loc_array(ew-1,ns)
         normal = (/ 0.0_dp, 1.0_dp /); fwdorbwd = (/ -1.0_dp, -1.0_dp /)
    elseif( phi .eq. 45.0d0 )then
         loc_latbc(1) = loc_array(ew-1,ns); loc_latbc(2) = loc_array(ew,ns); loc_latbc(3) = loc_array(ew-2,ns)
         loc_latbc(6) = loc_array(ew,ns-1); loc_latbc(4) = loc_array(ew,ns); loc_latbc(5) = loc_array(ew,ns-2)
         normal = (/ 1.0_dp/sqrt(2.0_dp), 1.0_dp/sqrt(2.0_dp) /); fwdorbwd = (/ -1.0_dp, -1.0_dp /)
    elseif( phi .eq. 90.0d0 )then
         loc_latbc(1) = loc_array(ew-1,ns); loc_latbc(2) = loc_array(ew,ns); loc_latbc(3) = loc_array(ew-2,ns)
         loc_latbc(4) = loc_array(ew,ns+1); loc_latbc(5) = loc_array(ew,ns-1)
         normal = (/ 1.0_dp, 0.0_dp /); fwdorbwd = (/ -1.0_dp, -1.0_dp /)
    elseif( phi .eq. 135.0d0 )then
         loc_latbc(1) = loc_array(ew-1,ns); loc_latbc(2) = loc_array(ew,ns); loc_latbc(3) = loc_array(ew-2,ns)
         loc_latbc(6) = loc_array(ew,ns+1); loc_latbc(4) = loc_array(ew,ns+2); loc_latbc(5) = loc_array(ew,ns)
         normal = (/ 1.0_dp/sqrt(2.0_dp), -1.0_dp/sqrt(2.0_dp) /); fwdorbwd = (/ -1.0_dp, 1.0_dp /)
    elseif( phi .eq. 180.0d0 )then
         loc_latbc(1) = loc_array(ew,ns+1); loc_latbc(4) = loc_array(ew,ns+2); loc_latbc(5) = loc_array(ew,ns)
         loc_latbc(2) = loc_array(ew+1,ns); loc_latbc(3) = loc_array(ew-1,ns)
         normal = (/ 0.0_dp, -1.0_dp /); fwdorbwd = (/ 1.0_dp, 1.0_dp /)
    elseif( phi .eq. 225.0d0 )then
         loc_latbc(1) = loc_array(ew+1,ns); loc_latbc(2) = loc_array(ew+2,ns); loc_latbc(3) = loc_array(ew,ns)
         loc_latbc(6) = loc_array(ew,ns+1); loc_latbc(4) = loc_array(ew,ns+2); loc_latbc(5) = loc_array(ew,ns);
         normal = (/ -1.0_dp/sqrt(2.0_dp), -1.0_dp/sqrt(2.0_dp) /); fwdorbwd = (/ 1.0_dp, 1.0_dp /)
    elseif( phi .eq. 270.0d0 )then
         loc_latbc(1) = loc_array(ew+1,ns); loc_latbc(2) = loc_array(ew+2,ns); loc_latbc(3) = loc_array(ew,ns)
         loc_latbc(4) = loc_array(ew,ns+1); loc_latbc(5) = loc_array(ew,ns-1)
         normal = (/ -1.0_dp, 0.0_dp /); fwdorbwd = (/ 1.0_dp, 1.0_dp /)
    else
         loc_latbc(1) = loc_array(ew+1,ns); loc_latbc(2) = loc_array(ew+2,ns); loc_latbc(3) = loc_array(ew,ns)
         loc_latbc(6) = loc_array(ew,ns-1); loc_latbc(4) = loc_array(ew,ns); loc_latbc(5) = loc_array(ew,ns-2)
         normal = (/ -1.0_dp/sqrt(2.0_dp), 1.0_dp/sqrt(2.0_dp) /); fwdorbwd = (/ 1.0_dp, -1.0_dp /)
    end if

  return

end subroutine getlatboundinfo

!***********************************************************************

function indshift( which, ew, ns, up, ewn, nsn, upn, loc_array, thck )  
! function output is a vector containing necessary index shifts for portions of 'othervel' 
! extracted near domain boundaries. NOTE that this contains duplication of some of the code in the 
! subroutine "getlatboundinfo", and the two could be combines at some point

  implicit none

  integer, intent(in) :: which
  integer, intent(in) :: ew, ns, up, ewn, nsn, upn 
  integer, dimension(ewn-1,nsn-1), intent(in) :: loc_array
  real (kind = dp), dimension(3,3), intent(in) :: thck

  integer, dimension(3) :: indshift
  integer :: upshift = 0, ewshift = 0, nsshift = 0

  real (kind = dp), dimension(3,3) :: mask, maskcorners
  real (kind = dp), dimension(3,3) :: thckmask
  real (kind = dp), dimension(3) :: testvect
  real (kind = dp) :: phi, deg2rad

  deg2rad = 3.141592654d0 / 180.0d0
  mask(:,1) = (/ 0.0d0, 180.0d0, 0.0d0 /)
  mask(:,2) = (/ 270.0d0, 0.0d0, 90.0d0 /)
  mask(:,3) = (/ 0.0d0, 360.0d0, 0.0d0 /)
  maskcorners(:,1) = (/ 225.0d0, 0.0d0, 135.0d0 /)
  maskcorners(:,2) = (/ 0.0d0, 0.0d0, 0.0d0 /)
  maskcorners(:,3) = (/ 315.0d0, 0.0d0, 45.0d0 /)

  if( up==1 )then   !! first treat bed/sfc, which aren't complicated
      upshift = 1
  elseif( up == upn )then
      upshift = -1
  else
      upshift = 0
  end if

  select case(which)

      case(0)   !! internal to lateral boundaries; no shift to ew,ns indices

          ewshift = 0; nsshift = 0;

      case(1)   !! at lateral boundaries; shift to ew,ns may be non-zero

          where( thck .ne. 0.0d0 ) 
            thckmask = 0.0_dp
          elsewhere( thck .eq. 0.0d0 )
            thckmask = 1.0d0
          endwhere

          testvect = sum( thckmask * mask, 1 )

        ! calculate the angle of the normal in cart. (x,y) system w/ 0 deg. at 12 O'clock, 90 deg. at 3 O'clock, etc.
        if( sum( sum( thckmask, 1 ) ) .eq. 1.0d0 )then
            phi = sum( sum( thckmask * maskcorners, 1 ) )
        else
            if( any( testvect .eq. 360.0d0 ) )then
                if( sum( testvect ) .eq. 450.0d0 )then
                    phi = 45.0d0
                elseif( sum( testvect ) .eq. 630.0d0 )then
                    phi = 315.0d0
                else
                    phi = 0.0d0
                end if
            elseif( all( testvect .ne. 360 ) )then
                phi = sum( testvect ) / sum( testvect/testvect, testvect .ne. 0.0d0 )
            end if
        end if

        ! define shift to indices based on this angle 
        if( phi .eq. 0.0d0 )then
            nsshift = -1; ewshift = 0
        elseif( phi .eq. 45.0d0 )then
            nsshift = -1; ewshift = -1
        elseif( phi .eq. 90.0d0 )then
            nsshift = 0; ewshift = -1
        elseif( phi .eq. 135.0d0 )then
            nsshift = 1; ewshift = -1
        elseif( phi .eq. 180.0d0 )then
            nsshift = 1; ewshift = 0
        elseif( phi .eq. 225.0d0 )then
            nsshift = 1; ewshift = 1
        elseif( phi .eq. 270.0d0 )then
            nsshift = 0; ewshift = 1
        elseif( phi .eq. 315.0d0 )then
            nsshift = -1; ewshift = 1
        end if

  end select

  indshift = (/ upshift, ewshift, nsshift /)

  return

end function indshift

!***********************************************************************

function maskvelostr(ewn,  nsn,   &
                     thck, stagthck)

! *sp* make 2d array containing a '0' (no ice present in associated cell), 
! '1' (cell contains ice and is in main body of ice sheet), 
! '2' (cell contains ice and is on a boundary)
! ... a 'cell' is an area bound by 4 grid points

  implicit none

  integer, intent(in) :: ewn, nsn 
  real (kind = dp), intent(in), dimension(:,:) :: thck, stagthck

  integer :: ew, ns
  integer, dimension(size(thck,1)-1,size(thck,2)-1) :: maskvelostr
  integer, dimension(3,3) :: template    ! *sp* added ...

  do ns = 1,nsn-1
  do ew = 1,ewn-1

! *sp* if all of the 4 points on the non-staggered grid contain non-zero ice thickness,
!       we are either in the main body of the ice or on a boundary

    if (all(thck(ew:ew+1,ns:ns+1) > 0.0_dp )) then

      ! *sp* if at the domain edges, define as a boundary
      if (ew == 1 .or. ew == ewn-1) then
        maskvelostr(ew,ns) = boundarys 
      else if (ns == 1 .or. ns == nsn-1) then
        maskvelostr(ew,ns) = boundarys 
      else      ! *sp* if not at domain edge, define as main body 
        maskvelostr(ew,ns) = mnbdy
      end if
    
    else if (any(thck(ew:ew+1,ns:ns+1) > 0.0_dp )) then

      maskvelostr(ew,ns) = boundarys

    else

      maskvelostr(ew,ns) = noice

    end if


!	if ( thck(ew,ns) > 20.0_dp / thk0 ) then
!		maskvelostr(ew,ns) = mnbdy
!	else
!		maskvelostr(ew,ns) = noice
!	end if
! end do; end do

!  do ns = 1,nsn-1; do ew = 1,ewn-1
!	if ( maskvelostr(ew-1,ns) .eq. noice .and. maskvelostr(ew,ns) .eq. mnbdy ) then
!		maskvelostr(ew,ns) = boundarys
!	else if ( maskvelostr(ew,ns) .eq. mnbdy .and. maskvelostr(ew+1,ns) .eq. noice ) then
!		maskvelostr(ew,ns) = boundarys
!	else if ( maskvelostr(ew,ns-1) .eq. noice .and. maskvelostr(ew,ns) .eq. mnbdy ) then
!		maskvelostr(ew,ns) = boundarys
!	else if ( maskvelostr(ew,ns+1) .eq. noice .and. maskvelostr(ew,ns) .eq. mnbdy ) then
!		maskvelostr(ew,ns) = boundarys
!	end if


  end do   ! ew
  end do   ! ns

  return

end function maskvelostr

!***********************************************************************

function viewsparse( dim, pcgrow, pcgcol, pcgval )

! function to construct sparse matrix from SLAPSOLV solver inputs 
! (for debugging test problems only!)

    implicit none

    integer, intent( in ), dimension( : ) :: pcgrow, pcgcol
    integer, intent( in ) :: dim
	real (kind = dp), intent( in ), dimension( : ) :: pcgval
	real (kind = dp), dimension(dim,dim) :: viewsparse 
    integer :: i, j, k, c

    c = 0; viewsparse = 0.0_dp

    do k = 1, size( pcgrow )
      if( pcgrow(k) .ne. 0 )then
        c = c + 1
      end if
    end do

    do k = 1, c
      i = pcgrow(k); j = pcgcol(k) 
      viewsparse( i, j ) = pcgval(k)
    end do
   
    return

end function viewsparse

!***********************************************************************

!whl - Currently this function is not called

function siabccalc(upn,      sigma,    &
                   which,    flwa,     &
                   thck,               & 
                   dusrfdew, dusrfdns, &
                   x,        y)

! calculate column of horiz. vel. using 0-order SIA, for lateral bc approx.

    implicit none

    integer, intent(in) :: upn 
    integer, intent(in) :: which, x, y
    real (kind = dp), dimension(:), intent(in) :: sigma
	real (kind = dp), intent(in) :: thck, dusrfdew, dusrfdns
    real (kind = dp), dimension(:,:,:), intent(in) :: flwa
	real (kind = dp), dimension(upn) :: shapefunc, siabccalc 
	real (kind = dp), dimension(2) :: slope
	real (kind = dp) :: u_surf, abar 
	real (kind = dp), dimension(1,2,2) :: abar0
    integer :: c, i, j

    c = 0; abar = 0.0_dp

    ! calculate dimensional sfc velocity from 0-order SIA
    slope = (/ dusrfdew, dusrfdns /) * (thk0/len0)

    abar0(1,:,:) = sum( flwa(:,x:(x+1),y:(y+1)), 1 ) / upn

    do i = 1,2; do j = 1,2
    if( abar0(1,i,j) .ne. 0.0_dp )then
      c = c + 1
      abar = abar + abar0(1,i,j)
    end if
    end do; end do

    abar = ( abar / c ) * vis0

    u_surf = ( 2.0_dp * abar ) / ( gn + 1.0_dp ) * ( -rhoi * grav * slope(which) )**gn & 
                * (thck * thk0)**(gn + 1.0_dp)

    shapefunc = 1.0_dp - ( sigma )**( gn + 1.0_dp )

    siabccalc = ( u_surf * shapefunc ) / vel0  !scale back to dimensionless value

    return

end function siabccalc

! *sp* function to specify map of betasquared
subroutine calcbetasquared (whichbabc,               & 
                            dew,         dns,        &
                            ewn,         nsn,        &
                            lsrf,        topg,       &
                            thck,                    &
                            thisvel,     othervel,   &
                            minTauf,                 &
                            betasquared, betafile) 

  integer, intent(in) :: whichbabc
  integer, intent(in) :: ewn, nsn

  real (kind = dp), intent(in) :: dew, dns
  real (kind = dp), intent(in), dimension(ewn,nsn) :: lsrf, topg, thck
  real (kind = dp), intent(in), dimension(ewn-1,nsn-1) :: thisvel, othervel, minTauf

  real (kind = dp), intent(out), dimension(ewn-1,nsn-1) :: betasquared

  character (len=30), intent(in), optional :: betafile
  real (kind = dp) :: smallnum = 1.0d-4
  real (kind = dp), dimension(ewn) :: grounded
  real (kind = dp) :: alpha, beta, dx, thck_gl, betalow, betahigh, roughness
  integer :: ew, ns

  select case(whichbabc)

    case(0)     ! constant value

      betasquared = 1.0d0

    case(1)     ! specify simple pattern, e.g. a strip of low B^2 to make an ice stream

      betasquared = 1.0d10

      do ew=1, ewn-1; do ns=10, nsn-10
        betasquared(ew,ns) = 1.0d0
      end do; end do

    case(2)     ! read in a map of betasquared from a file
                !whl - to do - Create a netCDF version of betafile

      open(unin,file=trim(betafile),status='old')
      read(unin,*) ((betasquared(ew,ns), ew=1,ewn-1), ns=1,nsn-1)
      close(unin)

      betasquared = betasquared / dsqrt( (thisvel*vel0*scyr)**2 + (othervel*vel0*scyr)**2 + (smallnum)**2 )
!      betasquared = betasquared * ( (thisvel*vel0*scyr)**2 + (othervel*vel0*scyr)**2 + (smallnum)**2 )**(-0.5d0)

    case(3)

      ! betasquared as proxy for till yield stress via Bueler iteration to achieve TauB ~ Tau0
      ! TauB = Tau0 * U_b * ( U_b^2 + U0 )^(-1/2)  
      ! ... where U0 is some small number to avoid div. by 0 (3rd term under radical)
      ! Here, betasquared is re-written to encompass the extra terms with the assumption that U_b in the radical
      ! can be taken from the value at the previous iteration.
      ! NOTE: previously assigned value for betasquared is now assumed equal to till yield stress, in Pa !!!

      betasquared = 1.0d10

      do ew=1, ewn-1; do ns=10, nsn-10
!        betasquared(ew,ns) = 1.05d0 * 7.0d3
        betasquared(ew,ns) = 5.0d3
      end do; end do

      betasquared = betasquared / dsqrt( (thisvel*vel0*scyr)**2 + (othervel*vel0*scyr)**2 + (smallnum)**2 )

    case(4)     ! same as case(3) but taking yield stress from basal processes model

      betasquared = minTauf
!      betasquared = minTauf / dsqrt( (thisvel*vel0*scyr)**2 + (othervel*vel0*scyr)**2 + (smallnum)**2 )
!      betasquared = betasquared * ( (thisvel*vel0*scyr)**2 + (othervel*vel0*scyr)**2 + (smallnum)**2 )**(-0.5d0)

    case(5)     ! simple 2d ice shelf

     betalow = 1.0d0
     betahigh = 1.0d5

     grounded = rhoo*( 0.0_dp - topg(:,3) ) / ( rhoi * thck(:,3) )

     do ew=1, ewn-1

      if( ( grounded(ew) .lt. 1.0_dp ) .and. ( grounded(ew+1) .lt. 1.0_dp ) )then

        betasquared(ew,:) = betahigh

      elseif( ( grounded(ew) .gt. 1.0_dp ) .and. ( grounded(ew+1) .gt. 1.0_dp ) )then

        betasquared(ew,:) = betalow

      elseif( grounded(ew) .lt. 1.0_dp .and. grounded(ew+1) .ge. 1.0_dp )then

        dx = ( 1.0_dp - grounded(ew) ) / ( ( grounded(ew+1) - grounded(ew) ) / dew )
        thck_gl = thck(ew,3) + ( thck(ew+1,3) - thck(ew,3) ) / dew * dx
        alpha = ( thck_gl - thck(ew,3) ) / ( thck(ew+1,3) - thck(ew,3) )
        betasquared(ew,:) = alpha*betahigh + ( 1.0_dp - alpha )*betalow

!        print *, 'x = ', (ew-1)*dew*len0
!        print *, 'dx = ', dx*len0
!        print *, 'x g.l. = ', ((ew-1)*dew+dx)*len0
!        print *, 'alpha = ', alpha
!        print *, 'betasquared at g.l. = ', betasquared(ew,3)

      end if

     end do

    case(6)     ! ISMIP-HOM experiment C; spatially periodic traction parameter

        do ns=1,nsn-1; do ew=1,ewn-1
            ! units as in ismip-hom document ( Pa * a * m^-1 )
            betasquared(ew,ns) = 1000 + 1000 * sin(2*pi/(lambda0*len0 )*(ns-1)*dns*len0) * &
                                 sin(2*pi/(lambda0*len0)*(ew-1)*dew*len0)
       end do; end do

    case(7)     ! circular ice shelf: set B^2 ~ 0 except for at center, where B^2 >> 0 to enforce u,v=0 there

      betasquared = 1.0d-5
      betasquared( (ewn-1)/2:(ewn-1)/2+1, (nsn-1)/2:(nsn-1)/2+1 ) = 1.0d10

!whl - change default to case 8
    case(8)    ! frozen (u=v=0) ice-bed interface

      betasquared = 1.0d10

  end select
  
  ! convert to dimensional model units ( Pa * s * m^-1 ) and then non-dimensionalize
  betasquared = ( betasquared * scyr ) / ( tau0 * tim0 / len0 )

end subroutine calcbetasquared


!***********************************************************************
!whl - copied this function from glam_velo.f90

function vertintg(upn, sigma, in) 
 
  implicit none 
 
  integer, intent(in) :: upn
  real (kind = dp), dimension(:), intent(in) :: sigma
  real (kind = dp), dimension(:), intent(in) :: in 
  real (kind = dp) :: vertintg 
 
  integer :: up 
 
  vertintg = 0.0d0 
 
  do up = upn-1, 1, -1 
    vertintg = vertintg + sum(in(up:up+1)) * dups(up) 
  end do 
 
  vertintg = vertintg / 2.0d0 
 
  return 
 
end function vertintg 

!***********************************************************************
!whl - copied this subroutine copied from glam_thck.f90

subroutine geom2derscros(dew,  dns,   &
                         ipvr, stagthck, opvrewns)        

! *sp* geometric (2nd) cross-deriv. for generic input variable 'ipvr', output as 'opvr'       
 
  implicit none 

  real (kind = dp), intent(in) :: dew, dns
  real (kind = dp), intent(out), dimension(:,:) :: opvrewns
  real (kind = dp), intent(in), dimension(:,:) :: ipvr, stagthck
 
!whl - Replace by a loop over ewn, nsn?  What is eoshift?
  where (stagthck /= 0.0d0)
    opvrewns = (eoshift(eoshift(ipvr,1,0.0_dp,2),1,0.0_dp,1) + ipvr   &
               - eoshift(ipvr,1,0.0_dp,1) - eoshift(ipvr,1,0.0_dp,2)) / (dew*dns)
  elsewhere
    opvrewns = 0.0d0
  end where
 
  return
 
end subroutine geom2derscros

!***********************************************************************
!whl - copied this subroutine from glam_thck.f90

subroutine geom2ders(ewn,    nsn,  &
                     dew,    dns,  &
                     ipvr,   stagthck,  &
                     opvrew, opvrns)       

! *sp* geometric 1st deriv. for generic input variable 'ipvr', 
!      output as 'opvr' (includes 'upwinding' for boundary values)

  implicit none 
 
  integer, intent(in) :: ewn, nsn 
  real (kind = dp), intent(in) :: dew, dns
  real (kind = dp), intent(out), dimension(:,:) :: opvrew, opvrns
  real (kind = dp), intent(in), dimension(:,:) :: ipvr, stagthck
 
  integer :: ew, ns 
  real (kind = dp) :: dewsq4, dnssq4
 
  integer :: pt(2)
 
  dewsq4 = 4.0d0 * dew * dew
  dnssq4 = 4.0d0 * dns * dns

!whl - Inline the functions below?
 
  do ns = 2, nsn-2 
  do ew = 2, ewn-2
    if (stagthck(ew,ns) .gt. 0.0d0) then
      opvrew(ew,ns) = centerew(ew,ns,ipvr,dewsq4)
      opvrns(ew,ns) = centerns(ew,ns,ipvr,dnssq4)
    else
      opvrew(ew,ns) = 0.0d0
      opvrns(ew,ns) = 0.0d0
    end if
  end do
  end do
 
! *** 2nd order boundaries using upwinding
 
  do ew = 1, ewn-1, ewn-2
 
    pt = whichway(ew)
 
    do ns = 2, nsn-2 
      if (stagthck(ew,ns) .gt. 0.0d0) then
        opvrew(ew,ns) = boundyew(ns,pt,ipvr,dewsq4)
        opvrns(ew,ns) = centerns(ew,ns,ipvr,dnssq4)
      else
        opvrew(ew,ns) = 0.0d0
        opvrns(ew,ns) = 0.0d0
      end if
    end do
 
  end do
 
  do ns = 1, nsn-1, nsn-2
 
    pt = whichway(ns)
 
    do ew = 2, ewn-2  
      if (stagthck(ew,ns) .gt. 0.0d0) then
        opvrew(ew,ns) = centerew(ew,ns,ipvr,dewsq4)
        opvrns(ew,ns) = boundyns(ew,pt,ipvr,dnssq4)
      else
        opvrew(ew,ns) = 0.0d0
        opvrns(ew,ns) = 0.0d0
      end if
    end do
 
  end do
 
  do ns = 1, nsn-1, nsn-2
    do ew = 1, ewn-1, ewn-2
      if (stagthck(ew,ns) .gt. 0.0d0) then
        pt = whichway(ew)
        opvrew(ew,ns) = boundyew(ns,pt,ipvr,dewsq4)
        pt = whichway(ns)
        opvrns(ew,ns) = boundyns(ew,pt,ipvr,dnssq4)
      else
        opvrew(ew,ns) = 0.0d0
        opvrns(ew,ns) = 0.0d0
      end if
    end do
  end do
 
end subroutine geom2ders
 
!***********************************************************************

  function centerew(ew, ns, ipvr, dewsq4)
 
    implicit none

    integer, intent(in) :: ew, ns 
    real(kind=dp), intent(in) :: ipvr(:,:)
    real(kind=dp), intent(in) :: dewsq4
    real (kind = dp) :: centerew
 
    centerew = (sum(ipvr(ew+2,ns:ns+1)) + sum(ipvr(ew-1,ns:ns+1)) - &
                sum(ipvr(ew+1,ns:ns+1)) - sum(ipvr(ew,ns:ns+1))) / dewsq4
 
    return
    
  end function centerew 
 
!***********************************************************************

  function centerns(ew, ns, ipvr, dnssq4)
 
    implicit none
 
    integer, intent(in) :: ew, ns 
    real(kind=dp), intent(in) :: ipvr(:,:)
    real(kind=dp), intent(in) :: dnssq4
    real (kind = dp) :: centerns
 
    centerns = (sum(ipvr(ew:ew+1,ns+2)) + sum(ipvr(ew:ew+1,ns-1)) - &
                sum(ipvr(ew:ew+1,ns+1)) - sum(ipvr(ew:ew+1,ns))) / dnssq4
 
    return
    
  end function centerns 
 
!***********************************************************************

  function boundyew(ns,pt,ipvr,dewsq4)
 
    implicit none

    integer, intent(in) :: ns
    integer, intent(in) :: pt(2)
    real(kind=dp), intent(in) :: ipvr(:,:)
    real(kind=dp), intent(in) :: dewsq4
    real (kind = dp) :: boundyew
 
    boundyew = pt(1) * (3.0d0 * sum(ipvr(pt(2),ns:ns+1)) - 7.0d0 * sum(ipvr(pt(2)+pt(1),ns:ns+1)) + &
               5.0d0 * sum(ipvr(pt(2)+2*pt(1),ns:ns+1)) - sum(ipvr(pt(2)+3*pt(1),ns:ns+1))) / dewsq4
 
    return
 
  end function boundyew
 
!***********************************************************************

  function boundyns(ew,pt,ipvr,dnssq4)
 
    implicit none
 
    integer, intent(in) :: ew
    integer, intent(in) :: pt(2)
    real(kind=dp), intent(in) :: ipvr(:,:)
    real(kind=dp), intent(in) :: dnssq4
    real (kind = dp) :: boundyns
 
    boundyns = pt(1) * (3.0d0 * sum(ipvr(ew:ew+1,pt(2))) - 7.0d0 * sum(ipvr(ew:ew+1,pt(2)+pt(1))) + &
               5.0d0 * sum(ipvr(ew:ew+1,pt(2)+2*pt(1))) - sum(ipvr(ew:ew+1,pt(2)+3*pt(1)))) / dnssq4
 
    return
 
  end function boundyns
 
!***********************************************************************

  function whichway(i)
 
    implicit none
 
    integer, intent(in) :: i
    integer :: whichway(2) 
 
    if (i == 1) then 
      whichway = (/1,1/)
    else
      whichway = (/-1,i+1/)
    end if
 
    return
 
  end function whichway
 
!***********************************************************************
!whl - copied this from module strswk
! CURRENTLY NOT BEING CALLED by 'findcoefstr' (commented out) 

!function vertibasl(whichbc,efvs,basalbc,rhs)            

! *sp* function called by 'findcoefstr' to stuff coeff. assoc. w/ deriv. 
!      "d/d_sigma ( _ d/d_sigma)" basal b.c. directly in sparse matrix 

!  implicit none 
!  integer, intent(in) :: whichbc
!  real (kind = dp), intent(in) :: efvs, basalbc
!  real (kind = dp), intent(inout) :: rhs
!  real (kind = dp), dimension(2) :: vertibasl 

!end function vertibasl 

!***********************************************************************
!whl - copied this from module general

    function hsum(inp) 
 
      implicit none
 
      real (kind = dp), dimension(:,:,:), intent(in) :: inp
      real (kind = dp), dimension(size(inp,dim=1)) :: hsum
 
      hsum = sum(sum(inp(:,:,:),dim=3),dim=2)
 
      return 
 
    end function hsum

!***********************************************************************

subroutine putpcgc(value,col,row) 
 
  implicit none
 
  integer, intent(in) :: row, col 
  real (kind = dp), intent(in) :: value 
 
  if (value /= 0.0d0) then
     pcgval(ct) = value
     pcgcol(ct) = col
     pcgrow(ct) = row
     ct = ct + 1
  end if
 
  return
 
end subroutine putpcgc 

!***********************************************************************

end module glam_strs2

!***********************************************************************
