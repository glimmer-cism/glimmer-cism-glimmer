! $Id: smb_wrappers.f90,v 1.1.2.1 2005-03-20 20:11:28 gethinw Exp $
!

!! Temporary wrapper for init driven from file
subroutine SMBInitForceWrapper(smb_p,nx,ny,dx,mindt)
  
  use smb_mecons
  use smb_force

  implicit none

  type(smb_params), intent(out) :: smb_p  ! instance of params object
  integer,intent(in)            :: nx     ! Number of grid point in x-direction
  integer,intent(in)            :: ny     ! Number of grid point in y-direction
  integer,intent(in)            :: dx     ! grid resolution (m)
  integer,intent(in)            :: mindt  ! timstep (minutes)

  ! Initialise some members of the derived type
  smb_p%nx    = nx
  smb_p%ny    = ny
  smb_p%dx    = dx
  smb_p%mindt = mindt

  call SMBInitForce(smb_p%nx,smb_p%ny,smb_p%mnl,smb_p%mnal,smb_p%dx,smb_p%dz0, &
       smb_p%zdeep,smb_p%dzdeep,smb_p%mindt,smb_p%outmap, &
       smb_p%albkind,smb_p%VanDus,smb_p%kindz0,smb_p%kstab, &
       smb_p%supice,smb_p%t0ini,smb_p%tempb,smb_p%rint,smb_p%denup, &
       smb_p%dendw,smb_p%dry,smb_p%alboldwet,smb_p%albolddry,smb_p%albfrs, &
       smb_p%albice,smb_p%albsup,smb_p%albwat,smb_p%decalb,smb_p%scswat, &
       smb_p%scsnow,smb_p%z0fs,smb_p%z0ms,smb_p%z0ice,smb_p%absice, &
       smb_p%abslow,smb_p%tsc0,smb_p%tsc1,smb_p%tsc10,smb_p%tscfac, &
       smb_p%zwind,smb_p%ztemp,smb_p%ist,smb_p%dt, &
       smb_p%x1,smb_p%x2,smb_p%y1,smb_p%y2, &
       smb_p%coefdz,smb_p%facstr,smb_p%bounl,smb_p%bounh, &
       smb_p%nal,smb_p%nl,smb_p%byear,smb_p%tatm,smb_p%wvel,&
       smb_p%prmass,smb_p%glr,smb_p%atmrad,smb_p%rlhf,smb_p%shf, &
       smb_p%terrad,smb_p%swrtop,smb_p%fluxrn,smb_p%enextr,smb_p%tempin,smb_p%t0, &
       smb_p%surfwt,smb_p%frsnow,smb_p%tsnwfall,smb_p%alb,smb_p%conden,smb_p%Psurf, &
       smb_p%eatm,smb_p%esat,smb_p%denair,smb_p%e0,smb_p%rhum,smb_p%frivel,smb_p%topsl, &
       smb_p%slope,smb_p%factin,smb_p%dfact,smb_p%zdum,smb_p%zbtm, &
       smb_p%zref,smb_p%hmin,smb_p%hlim,smb_p%hsurf,smb_p%zpack,smb_p%z0,smb_p%z0sn,smb_p%acc, &
       smb_p%rain,smb_p%condt,smb_p%runoff,smb_p%wetextent,smb_p%percozone, &
       smb_p%dryzone,smb_p%bareice,smb_p%den,smb_p%enst,smb_p%enfi,smb_p%totmfi, &
       smb_p%totmst,smb_p%dm,smb_p%absden,smb_p%absdm,smb_p%nk,smb_p%zbal,smb_p%ztal, &
       smb_p%jyear,smb_p%z,smb_p%dz,smb_p%rmass,smb_p%temp,smb_p%dens,smb_p%watcon, &
       smb_p%wetage,smb_p%rppart,smb_p%source)

end subroutine SMBInitForceWrapper

!! Wrapper for on-line init routine
subroutine SMBInitWrapper(smb_p,nx,ny,dx,mindt,paramFile)
  
  use smb_mecons
  use smb_debug

  implicit none

  type(smb_params),intent(out)  :: smb_p      ! params type instance
  integer,intent(in)            :: nx         ! Number of grid point in x-direction
  integer,intent(in)            :: ny         ! Number of grid point in y-direction
  integer,intent(in)            :: dx         ! grid resolution (m)
  integer,intent(in)            :: mindt      ! timstep (minutes)
  character(len=100),intent(in) :: paramFile  ! file holding parameters

  ! Initialise some members of the derived type
  smb_p%nx    = nx
  smb_p%ny    = ny
  smb_p%dx    = dx
  smb_p%mindt = mindt

  call SMBInit(smb_p%nx,smb_p%ny,smb_p%mnl,smb_p%mnal, &
       smb_p%dz0,smb_p%zdeep,smb_p%dzdeep, &
       smb_p%mindt,smb_p%outmap,paramFile,chcons, &
       smb_p%albkind,smb_p%VanDus,smb_p%kindz0,smb_p%kstab,smb_p%supice, &
       smb_p%alboldwet,smb_p%albolddry,smb_p%albfrs,smb_p%albice,smb_p%albsup,smb_p%albwat, &
       smb_p%decalb,smb_p%scswat,smb_p%scsnow, &
       smb_p%tsc0,smb_p%tsc1,smb_p%tsc10,smb_p%tscfac, &
       smb_p%z0ice,smb_p%zwind,smb_p%ztemp, &
       smb_p%t0ini,smb_p%tempb,smb_p%rint,smb_p%denup, &
       smb_p%dendw,smb_p%dry, &
       smb_p%z0fs,smb_p%z0ms,smb_p%absice,smb_p%abslow, &
       smb_p%ist,smb_p%dt, &
       smb_p%coefdz,smb_p%facstr,smb_p%bounl,smb_p%bounh, &
       smb_p%nal,smb_p%nl,smb_p%byear,smb_p%tatm,smb_p%wvel, &
       smb_p%prmass,smb_p%glr,smb_p%atmrad,smb_p%rlhf,smb_p%shf, &
       smb_p%terrad,smb_p%swrtop,smb_p%fluxrn,smb_p%enextr,smb_p%tempin,smb_p%Psurf,smb_p%t0, &
       smb_p%surfwt,smb_p%frsnow,smb_p%tsnwfall,smb_p%alb,smb_p%conden, &
       smb_p%eatm,smb_p%esat,smb_p%denair,smb_p%e0,smb_p%rhum,smb_p%frivel,smb_p%topsl, &
       smb_p%slope,smb_p%factin,smb_p%dfact,smb_p%zdum,smb_p%zbtm, &
       smb_p%zref,smb_p%hmin,smb_p%hlim,smb_p%hsurf,smb_p%zpack,smb_p%z0,smb_p%z0sn,smb_p%acc, &
       smb_p%rain,smb_p%condt,smb_p%runoff,smb_p%wetextent,smb_p%percozone, &
       smb_p%dryzone,smb_p%bareice,smb_p%den,smb_p%enst,smb_p%enfi,smb_p%totmfi, &
       smb_p%totmst,smb_p%dm,smb_p%absden,smb_p%absdm,smb_p%nk,smb_p%zbal,smb_p%ztal, &
       smb_p%jyear,smb_p%z,smb_p%dz,smb_p%rmass,smb_p%temp,smb_p%dens,smb_p%watcon, &
       smb_p%wetage,smb_p%rppart,smb_p%source)
  
end subroutine SMBInitWrapper

!! Wrapper for main loop routine which is
!! driven from file
subroutine SMBForceLoopWrapper(smb_p, &
     totacc,totrn,totcon,totoff,massbal)
  
  use smb_mecons
  use smb_force

  implicit none

  type(smb_params),intent(inout)        :: smb_p    !
  real(rk),dimension(smb_p%nx,smb_p%ny) :: totacc   ! total water equivalent snowfall
  real(rk),dimension(smb_p%nx,smb_p%ny) :: totrn    ! total rainfall
  real(rk),dimension(smb_p%nx,smb_p%ny) :: totcon   ! if -ve, evap, if +ve, condensation 
  real(rk),dimension(smb_p%nx,smb_p%ny) :: totoff   ! total runoff
  real(rk),dimension(smb_p%nx,smb_p%ny) :: massbal  ! totacc+totrn+totcon-totof

  call SMBForceLoop(smb_p%nx,smb_p%ny,smb_p%mnl,smb_p%mnal,smb_p%dz0, &
       smb_p%zdeep,smb_p%dzdeep,smb_p%mindt,smb_p%outmap, &
       smb_p%albkind, &
       smb_p%VanDus,smb_p%kindz0,smb_p%kstab, &
       smb_p%alboldwet,smb_p%albolddry,smb_p%albfrs,smb_p%albice, &
       smb_p%albsup,smb_p%albwat,smb_p%decalb,smb_p%scswat, &
       smb_p%scsnow,smb_p%z0fs,smb_p%z0ms,smb_p%z0ice,smb_p%abslow, &
       smb_p%zwind,smb_p%ztemp,smb_p%ist, &
       smb_p%dt,smb_p%iyear, &
       smb_p%x1,smb_p%x2,smb_p%y1,smb_p%y2, &
       smb_p%ratec,smb_p%cst,smb_p%dro,smb_p%coefdz,smb_p%facstr, &
       smb_p%bounl,smb_p%bounh,smb_p%nl, &
       smb_p%tatm,&
       smb_p%wvel,smb_p%prmass,smb_p%glr, &
       smb_p%atmrad,smb_p%rlhf,smb_p%shf,smb_p%terrad,smb_p%swrtop, &
       smb_p%fluxrn,smb_p%enextr,smb_p%tempin,smb_p%t0, &
       smb_p%surfwt,smb_p%frsnow,smb_p%tsnwfall,smb_p%alb,smb_p%conden, &
       smb_p%Psurf,smb_p%eatm,smb_p%esat,smb_p%denair,smb_p%e0, &
       smb_p%rhum,smb_p%frivel,smb_p%topsl, &
       smb_p%factin,smb_p%dfact,smb_p%zdum,smb_p%zbtm, &
       smb_p%zref,smb_p%hmin,smb_p%hlim,smb_p%hsurf, &
       smb_p%zpack,smb_p%z0,smb_p%z0sn, &
       smb_p%acc,smb_p%rain,smb_p%condt,smb_p%runoff, &
       smb_p%percozone,smb_p%dryzone,smb_p%bareice, &
       smb_p%den,smb_p%enst,smb_p%enfi,smb_p%totmfi, &
       smb_p%totmst,smb_p%dm,smb_p%absden,smb_p%absdm, &
       smb_p%nk,smb_p%jyear,smb_p%z,smb_p%dz,smb_p%rmass, &
       smb_p%temp,smb_p%dens,smb_p%watcon,smb_p%wetage, &
       smb_p%rppart,smb_p%source, &
       totacc,totrn,totcon,totoff,massbal, &
       smb_p%rkarma)

end subroutine SMBForceLoopWrapper

!! Wrapper for on-line main-loop routine
subroutine SMBStepWrapper(smb_p,   &
     totacc,totrn,totcon,totoff,massbal, &
     SLM,Ta,prec,U10m,V10m,hum,SWdown,LWdown,Psurf)   !MB wspeed replaced with U10m and V10m
  
  use smb_mecons
  use smb_const

  implicit none

  type(smb_params),intent(inout)                      :: smb_p   ! params for SMB scheme
  real(rk),intent(inout),dimension(smb_p%nx,smb_p%ny) :: totacc  ! total water equivalent snowfall
  real(rk),intent(inout),dimension(smb_p%nx,smb_p%ny) :: totrn   ! total rainfall
  real(rk),intent(inout),dimension(smb_p%nx,smb_p%ny) :: totcon  ! if -ve, evap, if +ve, condensation 
  real(rk),intent(inout),dimension(smb_p%nx,smb_p%ny) :: totoff  ! total runoff
  real(rk),intent(inout),dimension(smb_p%nx,smb_p%ny) :: massbal ! totacc+totrn+totcon-totof
  real(rk),intent(in),dimension(smb_p%nx,smb_p%ny)    :: SLM     ! **TODO**
  real(rk),intent(in),dimension(smb_p%nx,smb_p%ny)    :: Ta      ! air temp
  real(rk),intent(in),dimension(smb_p%nx,smb_p%ny)    :: prec    ! precip in right units? //MB No: see prmass
  real(rk),intent(in),dimension(smb_p%nx,smb_p%ny)    :: U10m    ! Dan has written convesrion
  real(rk),intent(in),dimension(smb_p%nx,smb_p%ny)    :: V10m    ! Dan has written convesrion 
  real(rk),intent(in),dimension(smb_p%nx,smb_p%ny)    :: hum     ! Relative humidity (%)
  real(rk),intent(in),dimension(smb_p%nx,smb_p%ny)    :: SWdown  ! Need to pick apart net// MB Dan should take care of this too.
  real(rk),intent(in),dimension(smb_p%nx,smb_p%ny)    :: LWdown  ! ditto
  real(rk),intent(in),dimension(smb_p%nx,smb_p%ny)    :: Psurf   ! surface pressure
  
  !local variables
  real(rk),dimension(smb_p%nx,smb_p%ny) :: prmass  ! precip in right units? //MB Yes  
  real(rk),dimension(smb_p%nx,smb_p%ny) :: wspeed  ! recalculated if U and V are provided instead of wspeed (call windcalc)


  ! Processing some input: get wind speed, precip in the right units, and humidity STILL TO FIX!
  call windcalc (U10m,V10m,wspeed)
  prmass=prec*(60*smb_p%mindt)        ! precip in mm w.e cumulated for 1 timestep


  call SMBStep(smb_p%mnl,smb_p%mnal,smb_p%dz0,smb_p%zdeep,smb_p%dzdeep, &
       smb_p%mindt,smb_p%albkind, &
       smb_p%VanDus,smb_p%kindz0,smb_p%kstab, &
       smb_p%alboldwet,smb_p%albolddry,smb_p%albfrs,smb_p%albice, &
       smb_p%albsup,smb_p%albwat,smb_p%decalb,smb_p%scswat, &
       smb_p%scsnow,smb_p%z0fs,smb_p%z0ms,smb_p%z0ice,smb_p%abslow, &
       smb_p%zwind,smb_p%ztemp,smb_p%ist,smb_p%dt, &
       smb_p%iyear, &
       smb_p%x1,smb_p%x2,smb_p%y1,smb_p%y2,smb_p%ratec,smb_p%cst, &
       smb_p%dro,smb_p%coefdz,smb_p%facstr, &
       smb_p%bounl,smb_p%bounh,smb_p%nl,SLM, &
       smb_p%rlhf,smb_p%shf,smb_p%terrad,smb_p%swrtop, &
       smb_p%fluxrn,smb_p%enextr,smb_p%tempin,smb_p%t0,smb_p%surfwt, &
       smb_p%frsnow,smb_p%tsnwfall,smb_p%alb,smb_p%conden, &
       smb_p%eatm,smb_p%esat,smb_p%denair,smb_p%e0, &
       smb_p%frivel,smb_p%topsl, &
       smb_p%factin,smb_p%dfact,smb_p%zdum,smb_p%zbtm,smb_p%zref, &
       smb_p%hmin,smb_p%hlim,smb_p%hsurf,smb_p%zpack,smb_p%z0,smb_p%z0sn, &
       smb_p%acc,smb_p%rain,smb_p%condt,smb_p%runoff,smb_p%percozone, &
       smb_p%dryzone,smb_p%bareice, &
       smb_p%den,smb_p%enst,smb_p%enfi,smb_p%totmfi,smb_p%totmst,smb_p%dm, &
       smb_p%absden,smb_p%absdm, &
       smb_p%nk,smb_p%jyear,smb_p%z,smb_p%dz,smb_p%rmass,smb_p%temp, &
       smb_p%dens,smb_p%watcon,smb_p%wetage, &
       smb_p%rppart,smb_p%source, &
       Ta,prmass,wspeed, &
       hum,SWdown,LWdown,Psurf,totacc, &
       totrn,totcon,totoff,massbal, &
       smb_p%rkarma)

end subroutine SMBStepWrapper

!! Wrapper for the cleanup routine 
subroutine SMBCleanupWrapper(smb_p)
  
  use smb_mecons

  implicit none

  type(smb_params) :: smb_p

  call SMBCleanup(smb_p%nal,smb_p%nl, &
       smb_p%tatm,smb_p%wvel, &
       smb_p%prmass,smb_p%glr,smb_p%atmrad,smb_p%rlhf,smb_p%shf, &
       smb_p%terrad,smb_p%swrtop,smb_p%fluxrn,smb_p%enextr,smb_p%tempin,smb_p%t0, &
       smb_p%surfwt,smb_p%frsnow,smb_p%tsnwfall,smb_p%alb,smb_p%conden,smb_p%Psurf, &
       smb_p%eatm,smb_p%esat,smb_p%denair,smb_p%e0,smb_p%rhum,smb_p%frivel,smb_p%topsl, &
       smb_p%slope,smb_p%factin,smb_p%dfact,smb_p%zdum,smb_p%zbtm, &
       smb_p%zref,smb_p%hmin,smb_p%hlim,smb_p%hsurf,smb_p%zpack,smb_p%z0,smb_p%z0sn,smb_p%acc, &
       smb_p%rain,smb_p%condt,smb_p%runoff,smb_p%wetextent,smb_p%percozone, &
       smb_p%dryzone,smb_p%bareice,smb_p%den,smb_p%enst,smb_p%enfi,smb_p%totmfi, &
       smb_p%totmst,smb_p%dm,smb_p%absden,smb_p%absdm,smb_p%nk,smb_p%zbal,smb_p%ztal, &
       smb_p%byear,smb_p%jyear,smb_p%z,smb_p%dz,smb_p%rmass,smb_p%temp,smb_p%dens,smb_p%watcon, &
       smb_p%wetage,smb_p%rppart,smb_p%source)

end subroutine SMBCleanupWrapper

