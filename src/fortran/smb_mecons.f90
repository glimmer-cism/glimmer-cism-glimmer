!
! $Id: smb_mecons.f90,v 1.1.2.1 2005-03-20 20:11:28 gethinw Exp $
!
!****************************************************************************
!
! MODULE: SMB
!
! SURFACE MASS BALANCE MODEL based largely on Greuell's code
!
!****************************************************************************

module smb_mecons

  use glimmer_global

  implicit none

  type smb_params

     ! Vars given default values--can be overwritten from file
     ! Called in file 'modsetup':
     integer  :: nx         = 75      ! Num grid pts in x-direction
     integer  :: ny         = 140     ! Num grid pts in y-direction
     integer  :: mnl        = 100     ! MAXIMUM number of vertical layer
     integer  :: mnal       = 3       ! MAXIMUM number of annual layer
     integer  :: dx         = 2000    ! grid resolution (m)
     real(rk) :: dz0        = 0.09    ! dz0 is the desired grid point size at the surface
     real(rk) :: zdeep      = 30.0    ! zdeep the depth of the model 
     real(rk) :: dzdeep     = 4.0     ! size of lowermost grid cell
     integer  :: mindt      = 30      ! timestep in min
     logical  :: outmap     = .false. !extensive output saving (includes saving map output, i.e. large files) 
     
     ! Called in file 'method':
     integer  :: albkind = 1         !choose between 3 different albedo schemes 
     logical  :: VanDus  = .false.   !Conductivity as in Van Dusen,1929 (true) or Sturm, 1997
     integer  :: kindz0  = 4         !roughness lengths (1:z0t=z0; 2:z0t=0.00001; 3:Brutsaert; 4:Andreas)
     integer  :: kstab   = 2         !stability functions(1=Businger; 2=Hogstrom; 3=Beljaars and Holtslag; 4=Dyer)
     logical  :: supice  = .true.    !run off according to Zuo and Oe (true) or all water on top of ice runs off (false)
     
     ! Called in file 'physpar':
     real(rk) :: t0ini     = -5.0    !Initial surface temperature (C)
     real(rk) :: tempb     = -10.0   !Temperature in the lower layer (C)
     real(rk) :: rint      = 1.00    !initial snow pack depth at elevation of measurements (m)
     real(rk) :: denup     = 300.0   !density of upper layer (kg/m3)
     real(rk) :: dendw     = 910.0   !density of lower layer (kg/m3)
     logical  :: dry       = .true.  !initial snow layer dry (true) or wet (false)
     real(rk) :: alboldwet = 0.50    !albedo of old wet snow
     real(rk) :: albolddry = 0.65    !albedo of old dry snow
     real(rk) :: albfrs    = 0.84    !albedo of fresh snow
     real(rk) :: albice    = 0.44    !albedo of ice
     real(rk) :: albsup    = 0.50    !albedo of superimposed ice
     real(rk) :: albwat    = 0.15    !albedo of water
     real(rk) :: decalb    = 15.0    !decay time for snow albedo
     real(rk) :: scswat    = 200.0   !scaling height (mm) for the effect of surficial water on the surface albedo 
     real(rk) :: scsnow    = 0.01    !scaling height for the effect of fresh snow
     real(rk) :: z0fs      = 0.00012 !momentum roughness length of frozen snow (m)
     real(rk) :: z0ms      = 0.0023  !momentum roughness length of melting snow (m)
     real(rk) :: z0ice     = 0.0044  !momentum roughness length of ice (m)
     real(rk) :: absice    = 2.8     !radiation penetration in ice (m-1)
     real(rk) :: abslow    = 20.0    !radiation penetration in snow (m-1)
     real(rk) :: tsc0      = 15.0    !time scale for runoff of surficial water on a horizontal surface (days) 
     real(rk) :: tsc1      = 0.0001  !time scale for runoff of surficial water on a surface with slope of 1 degree (days)
     real(rk) :: tsc10     = 0.0     !time scale for runoff of surficial water on a steep surface (days) 
     real(rk) :: tscfac    = 10.0    !ratio of time scales of runoff below and above surface
     real(rk) :: zwind     = 10.0    !height of wind velocity input (m) for bulk method or highest level profile method
     real(rk) :: ztemp     = 2.0     !height of temperature input (m) for bulk method or highest level profile method

     integer  :: iyear     = 1       ! current year
     
     !Model variables
     integer  :: ist                 ! loop stuff
     real(rk) :: dt                  ! timestep in sec
     integer  :: x1                  ! coordinates target 
     integer  :: x2                  ! coordinates target 
     integer  :: y1                  ! coordinates target 
     integer  :: y2                  ! coordinates target 
     real(rk) :: ratec               ! use for dry snow densification
     real(rk) :: cst                 ! use for dry snow densification
     real(rk) :: dro                 ! use for dry snow densification
     real(rk) :: rkarma              ! von karman constant
     
     ! to calculate and maintain the vertical grid
     real(rk) :: coefdz              !
     real(rk) :: facstr              !
     real(rk) :: bounl               !
     real(rk) :: bounh               !
     
     !----------ALLOCATABLE ARRAYS---------------------------
     ! *NB* initialisation using '=> null()' if fortran 95 only

     integer,pointer,dimension(:,:)    :: nal => null()             !nb of annual layer
     integer,pointer,dimension(:,:)    :: nl => null()              !nb of layer
     integer,pointer,dimension(:,:)    :: byear => null()           !year at the beginning of the simulation
     real(rk),pointer,dimension(:,:)   :: tatm => null()            !atmospheric input
     real(rk),pointer,dimension(:,:)   :: wvel => null()            !atmospheric input
     real(rk),pointer,dimension(:,:)   :: prmass => null()          !atmospheric input
     real(rk),pointer,dimension(:,:)   :: glr => null()             !atmospheric input
     real(rk),pointer,dimension(:,:)   :: atmrad => null()          !atmospheric input
     real(rk),pointer,dimension(:,:)   :: rlhf => null()            !energy related
     real(rk),pointer,dimension(:,:)   :: shf => null()             !energy related
     real(rk),pointer,dimension(:,:)   :: terrad => null()          !energy related
     real(rk),pointer,dimension(:,:)   :: swrtop => null()          !energy related
     real(rk),pointer,dimension(:,:)   :: fluxrn => null()          !energy related
     real(rk),pointer,dimension(:,:)   :: enextr => null()          !energy related
     real(rk),pointer,dimension(:,:)   :: tempin => null()          !temperature
     real(rk),pointer,dimension(:,:)   :: t0 => null()              !temperature
     real(rk),pointer,dimension(:,:)   :: surfwt => null()          !albedo related
     real(rk),pointer,dimension(:,:)   :: frsnow => null()          !albedo related
     real(rk),pointer,dimension(:,:)   :: tsnwfall => null()        !albedo related
     real(rk),pointer,dimension(:,:)   :: alb => null()             !albedo related
     real(rk),pointer,dimension(:,:)   :: conden => null()          !Turb HF related
     real(rk),pointer,dimension(:,:)   :: Psurf => null()           !Turb HF related
     real(rk),pointer,dimension(:,:)   :: eatm => null()            !Turb HF related
     real(rk),pointer,dimension(:,:)   :: esat => null()            !Turb HF related
     real(rk),pointer,dimension(:,:)   :: denair => null()          !Turb HF related
     real(rk),pointer,dimension(:,:)   :: e0 => null()              !Turb HF related
     real(rk),pointer,dimension(:,:)   :: rhum => null()            !Turb HF related
     real(rk),pointer,dimension(:,:)   :: frivel => null()          !Turb HF related
     real(rk),pointer,dimension(:,:)   :: topsl => null()           !runoff related
     real(rk),pointer,dimension(:,:)   :: slope => null()           !runoff related
     real(rk),pointer,dimension(:,:)   :: factin => null()          !runoff related
     real(rk),pointer,dimension(:,:)   :: dfact => null()           !runoff related
     real(rk),pointer,dimension(:,:)   :: zdum => null()            !grid adjustment
     real(rk),pointer,dimension(:,:)   :: zbtm => null()            !grid adjustment
     real(rk),pointer,dimension(:,:)   :: zref => null()            !geometry
     real(rk),pointer,dimension(:,:)   :: hmin => null()            !geometry
     real(rk),pointer,dimension(:,:)   :: hlim => null()            !geometry
     real(rk),pointer,dimension(:,:)   :: hsurf => null()           !geometry
     real(rk),pointer,dimension(:,:)   :: zpack => null()           !geometry
     real(rk),pointer,dimension(:,:)   :: z0 => null()              !geometry
     real(rk),pointer,dimension(:,:)   :: z0sn => null()            !geometry
     real(rk),pointer,dimension(:,:)   :: acc => null()             !mass balance
     real(rk),pointer,dimension(:,:)   :: rain => null()            !mass balance
     real(rk),pointer,dimension(:,:)   :: condt => null()           !mass balance
     real(rk),pointer,dimension(:,:)   :: runoff => null()          !mass balance
     real(rk),pointer,dimension(:,:)   :: wetextent => null()       !defined zones
     real(rk),pointer,dimension(:,:)   :: percozone => null()       !defined zones
     real(rk),pointer,dimension(:,:)   :: dryzone => null()         !defined zones
     real(rk),pointer,dimension(:,:)   :: bareice => null()         !defined zones
     real(rk),pointer,dimension(:,:)   :: den => null()             !mass-energy cons.
     real(rk),pointer,dimension(:,:)   :: enst => null()            !mass-energy cons.
     real(rk),pointer,dimension(:,:)   :: enfi => null()            !mass-energy cons.
     real(rk),pointer,dimension(:,:)   :: totmfi => null()          !mass-energy cons.
     real(rk),pointer,dimension(:,:)   :: totmst => null()          !mass-energy cons.
     real(rk),pointer,dimension(:,:)   :: dm => null()              !mass-energy cons.
     real(rk),pointer,dimension(:,:)   :: absden => null()          !mass-energy cons.
     real(rk),pointer,dimension(:,:)   :: absdm => null()           !mass-energy cons.
     integer,pointer,dimension(:,:,:)  :: nk => null()              !Nb of nl per annual layer
     integer,pointer,dimension(:,:,:)  :: jyear => null()           !year for deposition of layer nk
     real(rk),pointer,dimension(:,:,:) :: zbal => null()            !depth of btm of annual lyr 
     real(rk),pointer,dimension(:,:,:) :: ztal => null()            !depth of top of annual lyr 
     real(rk),pointer,dimension(:,:,:) :: z => null()               !grid variables
     real(rk),pointer,dimension(:,:,:) :: dz => null()              !grid variables
     real(rk),pointer,dimension(:,:,:) :: rmass => null()           !englacial var.
     real(rk),pointer,dimension(:,:,:) :: temp => null()            !englacial var.
     real(rk),pointer,dimension(:,:,:) :: dens => null()            !englacial var.
     real(rk),pointer,dimension(:,:,:) :: watcon => null()          !englacial var.
     real(rk),pointer,dimension(:,:,:) :: wetage => null()          !englacial var.
     real(rk),pointer,dimension(:,:,:) :: rppart => null()          !energy penetration
     real(rk),pointer,dimension(:,:,:) :: source => null()          !energy penetration

  end type smb_params
  
contains

  subroutine SMBInit(nx,ny,mnl,mnal, &
       dz0,zdeep,dzdeep, &
       mindt,outmap,paramFile,chcons, &
       albkind,VanDus,kindz0,kstab,supice, &
       alboldwet,albolddry,albfrs,albice,albsup,albwat, &
       decalb,scswat,scsnow, &
       tsc0,tsc1,tsc10,tscfac, &
       z0ice,zwind,ztemp, &
       t0ini,tempb,rint,denup, &
       dendw,dry, &
       z0fs,z0ms,absice,abslow, &
       ist,dt, &
       coefdz,facstr,bounl,bounh, &
       nal,nl,byear,tatm,wvel, &
       prmass,glr,atmrad,rlhf,shf, &
       terrad,swrtop,fluxrn,enextr,tempin,Psurf,t0, &
       surfwt,frsnow,tsnwfall,alb,conden, &
       eatm,esat,denair,e0,rhum,frivel,topsl, &
       slope,factin,dfact,zdum,zbtm, &
       zref,hmin,hlim,hsurf,zpack,z0,z0sn,acc, &
       rain,condt,runoff,wetextent,percozone, &
       dryzone,bareice,den,enst,enfi,totmfi, &
       totmst,dm,absden,absdm,nk,zbal,ztal, &
       jyear,z,dz,rmass,temp,dens,watcon, &
       wetage,rppart,source)

    use smb_const

    implicit none

    ! Parmaters

    character(len=100) :: paramFile  ! file holding parameters

    logical  :: dry      !initial snow layer dry (true) or wet (false)
    logical  :: outmap   !extensive output saving (includes saving map output, i.e. large files) 
    logical  :: chcons   !check for mass convervation (true) or not
    logical  :: VanDus   !Conductivity as in Van Dusen,1929 (true) or Sturm, 1997
    logical  :: supice   !run off according to Zuo and Oe (true) or all water on top of ice runs off (false)

    integer  :: nx       !Number of grid point in x-direction
    integer  :: ny       !Number of grid point in y-direction
    integer  :: mnl      !MAXIMUM number of vertical layer
    integer  :: mnal     !MAXIMUM number of annual layer
    integer  :: mindt    !timestep in min
    integer  :: ist      ! loop stuff
    integer  :: albkind  !choose between 3 different albedo schemes 
    integer  :: kindz0   !roughness lengths (1:z0t=z0; 2:z0t=0.00001; 3:Brutsaert; 4:Andreas)
    integer  :: kstab    !stability functions(1=Businger; 2=Hogstrom; 3=Beljaars and Holtslag; 4=Dyer)

    real(rk)   :: alboldwet   !albedo of old wet snow
    real(rk)   :: albolddry   !albedo of old dry snow
    real(rk)   :: albfrs      !albedo of fresh snow
    real(rk)   :: albice      !albedo of ice
    real(rk)   :: albsup      !albedo of superimposed ice
    real(rk)   :: albwat      !albedo of water
    real(rk)   :: decalb      !decay time for snow albedo
    real(rk)   :: scswat      !scaling height (mm) for the effect of surficial water on the surface albedo 
    real(rk)   :: scsnow      !scaling height for the effect of fresh snow
    real(rk)   :: tsc0        !time scale for runoff of surficial water on a horizontal surface (days) 
    real(rk)   :: tsc1        !time scale for runoff of surficial water on a surface with slope of 1 degree (days)
    real(rk)   :: tsc10       !time scale for runoff of surficial water on a steep surface (days) 
    real(rk)   :: tscfac      !ratio of time scales of runoff below and above surface
    real(rk)   :: z0fs        !momentum roughness length of frozen snow (m)
    real(rk)   :: z0ms        !momentum roughness length of melting snow (m)
    real(rk)   :: z0ice       !momentum roughness length of ice (m)
    real(rk)   :: zwind       !height of wind velocity input (m) for bulk method or highest level profile method
    real(rk)   :: ztemp       !height of temperature input (m) for bulk method or highest level profile method

    real(rk) :: dz0      !dz0 is the desired grid point size at the surface
    real(rk) :: zdeep    !zdeep the depth of the model 
    real(rk) :: dzdeep   !size of lowermost grid cell
    real(rk) :: t0ini    !Initial surface temperature (C)
    real(rk) :: tempb    !Temperature in the lower layer (C)
    real(rk) :: rint     !initial snow pack depth at elevation of measurements (m)
    real(rk) :: denup    !density of upper layer (kg/m3)
    real(rk) :: dendw    !density of lower layer (kg/m3)
    real(rk) :: absice   !radiation penetration in ice (m-1)
    real(rk) :: abslow   !radiation penetration in snow (m-1)
    real(rk) :: dt       ! timestep in sec
    real(rk) :: coefdz   !
    real(rk) :: facstr   !
    real(rk) :: bounl    !
    real(rk) :: bounh    !

    integer,pointer,dimension(:,:)    :: nal         !nb of annual layer
    integer,pointer,dimension(:,:)    :: nl          !nb of layer
    integer,pointer,dimension(:,:)    :: byear       !year at the beginning of the experiment
    real(rk),pointer,dimension(:,:)   :: tatm        !atmospheric input
    real(rk),pointer,dimension(:,:)   :: wvel        !atmospheric input
    real(rk),pointer,dimension(:,:)   :: prmass      !atmospheric input
    real(rk),pointer,dimension(:,:)   :: glr         !atmospheric input
    real(rk),pointer,dimension(:,:)   :: atmrad      !atmospheric input
    real(rk),pointer,dimension(:,:)   :: rlhf        !energy related
    real(rk),pointer,dimension(:,:)   :: shf         !energy related
    real(rk),pointer,dimension(:,:)   :: terrad      !energy related
    real(rk),pointer,dimension(:,:)   :: swrtop      !energy related
    real(rk),pointer,dimension(:,:)   :: fluxrn      !energy related
    real(rk),pointer,dimension(:,:)   :: enextr      !energy related
    real(rk),pointer,dimension(:,:)   :: tempin      !temperature
    real(rk),pointer,dimension(:,:)   :: t0          !temperature
    real(rk),pointer,dimension(:,:)   :: surfwt      !albedo related
    real(rk),pointer,dimension(:,:)   :: frsnow      !albedo related
    real(rk),pointer,dimension(:,:)   :: tsnwfall    !albedo related
    real(rk),pointer,dimension(:,:)   :: alb         !albedo related
    real(rk),pointer,dimension(:,:)   :: conden      !Turb HF related
    real(rk),pointer,dimension(:,:)   :: Psurf       !Turb HF related
    real(rk),pointer,dimension(:,:)   :: eatm        !Turb HF related
    real(rk),pointer,dimension(:,:)   :: esat        !Turb HF related
    real(rk),pointer,dimension(:,:)   :: denair      !Turb HF related
    real(rk),pointer,dimension(:,:)   :: e0          !Turb HF related
    real(rk),pointer,dimension(:,:)   :: rhum        !Turb HF related
    real(rk),pointer,dimension(:,:)   :: frivel      !Turb HF related
    real(rk),pointer,dimension(:,:)   :: topsl       !runoff related
    real(rk),pointer,dimension(:,:)   :: slope       !runoff related
    real(rk),pointer,dimension(:,:)   :: factin      !runoff related
    real(rk),pointer,dimension(:,:)   :: dfact       !runoff related
    real(rk),pointer,dimension(:,:)   :: zdum        !grid adjustment
    real(rk),pointer,dimension(:,:)   :: zbtm        !grid adjustment
    real(rk),pointer,dimension(:,:)   :: zref        !geometry
    real(rk),pointer,dimension(:,:)   :: hmin        !geometry
    real(rk),pointer,dimension(:,:)   :: hlim        !geometry
    real(rk),pointer,dimension(:,:)   :: hsurf       !geometry
    real(rk),pointer,dimension(:,:)   :: zpack       !geometry
    real(rk),pointer,dimension(:,:)   :: z0          !geometry
    real(rk),pointer,dimension(:,:)   :: z0sn        !geometry
    real(rk),pointer,dimension(:,:)   :: acc         !mass balance
    real(rk),pointer,dimension(:,:)   :: rain        !mass balance
    real(rk),pointer,dimension(:,:)   :: condt       !mass balance
    real(rk),pointer,dimension(:,:)   :: runoff      !mass balance
    real(rk),pointer,dimension(:,:)   :: wetextent   !defined zones
    real(rk),pointer,dimension(:,:)   :: percozone   !defined zones
    real(rk),pointer,dimension(:,:)   :: dryzone     !defined zones
    real(rk),pointer,dimension(:,:)   :: bareice     !defined zones
    real(rk),pointer,dimension(:,:)   :: den         !mass-energy cons.
    real(rk),pointer,dimension(:,:)   :: enst        !mass-energy cons.
    real(rk),pointer,dimension(:,:)   :: enfi        !mass-energy cons.
    real(rk),pointer,dimension(:,:)   :: totmfi      !mass-energy cons.
    real(rk),pointer,dimension(:,:)   :: totmst      !mass-energy cons.
    real(rk),pointer,dimension(:,:)   :: dm          !mass-energy cons.
    real(rk),pointer,dimension(:,:)   :: absden      !mass-energy cons.
    real(rk),pointer,dimension(:,:)   :: absdm       !mass-energy cons.
    integer,pointer,dimension(:,:,:)  :: nk          !Nb of nl per annual layer
    integer,pointer,dimension(:,:,:)  :: jyear       !year for deposition of layer nk
    real(rk),pointer,dimension(:,:,:) :: zbal        !depth of btm of annual lyr 
    real(rk),pointer,dimension(:,:,:) :: ztal        !depth of top of annual lyr 
    real(rk),pointer,dimension(:,:,:) :: z           !grid variables
    real(rk),pointer,dimension(:,:,:) :: dz          !grid variables
    real(rk),pointer,dimension(:,:,:) :: rmass       !englacial var.
    real(rk),pointer,dimension(:,:,:) :: temp        !englacial var.
    real(rk),pointer,dimension(:,:,:) :: dens        !englacial var.
    real(rk),pointer,dimension(:,:,:) :: watcon      !englacial var.
    real(rk),pointer,dimension(:,:,:) :: wetage      !englacial var.
    real(rk),pointer,dimension(:,:,:) :: rppart      !energy penetration
    real(rk),pointer,dimension(:,:,:) :: source      !energy penetration

    ! Local Variables
    integer :: i
    integer :: x
    integer :: y

    ! Read params from file
    call SMBReadParams(paramFile, &
         mnl,mnal,dz0,zdeep,dzdeep, &
         outmap,chcons,t0ini,tempb,rint,denup,dendw,dry, &
         alboldwet,albolddry,albfrs,albice,albsup,albwat, &
         decalb,scswat,scsnow,z0fs,z0ms,z0ice,absice, &
         abslow,tsc0,tsc1,tsc10,tscfac,zwind,ztemp, &
         albkind,VanDus,kindz0,kstab,supice)
        
    !computed from values read from file
    dt=mindt*60.0                          !timestep in sec

    ! some coefficients needed to calculate and maintain the vertical grid are computed 
    coefdz=(dzdeep-dz0)/zdeep
    bounl=blc*dz0
    bounh=bhc*dz0
    facstr=(absice-abslow)/(denice-denlow) !all terms defined in myin/inppar1 (GK94)

    ! Size of following arrays depend upon number or grid points in
    ! x- and y-directions and upon the maximum number of vertical
    ! and anual layers
    ! Initialise too for safety
    allocate(nal(nx,ny)); nal=0
    allocate(nl(nx,ny)); nl=0
    allocate(tatm(nx,ny)); tatm=0
    allocate(wvel(nx,ny)); wvel=0
    allocate(prmass(nx,ny)); prmass=0
    allocate(glr(nx,ny)); glr=0
    allocate(atmrad(nx,ny)); atmrad=0
    allocate(rlhf(nx,ny)); rlhf=0
    allocate(shf(nx,ny)); shf=0
    allocate(terrad(nx,ny)); terrad=0
    allocate(swrtop(nx,ny)); swrtop=0
    allocate(fluxrn(nx,ny)); fluxrn=0
    allocate(enextr(nx,ny)); enextr=0
    allocate(tempin(nx,ny)); tempin=0
    allocate(t0(nx,ny)); t0=0
    allocate(surfwt(nx,ny)); surfwt=0
    allocate(frsnow(nx,ny)); frsnow=0
    allocate(tsnwfall(nx,ny)); tsnwfall=0
    allocate(alb(nx,ny)); alb=0
    allocate(conden(nx,ny)); conden=0
    allocate(Psurf(nx,ny)); Psurf=0
    allocate(eatm(nx,ny)); eatm=0
    allocate(esat(nx,ny)); esat=0
    allocate(denair(nx,ny)); denair=0
    allocate(e0(nx,ny)); e0=0
    allocate(rhum(nx,ny)); rhum=0
    allocate(frivel(nx,ny)); frivel=0
    allocate(topsl(nx,ny)); topsl=0
    allocate(slope(nx,ny)); slope=0
    allocate(factin(nx,ny)); factin=0
    allocate(dfact(nx,ny)); dfact=0
    allocate(zdum(nx,ny)); zdum=0
    allocate(zbtm(nx,ny)); zbtm=0
    allocate(zref(nx,ny)); zref=0
    allocate(hmin(nx,ny)); hmin=0
    allocate(hlim(nx,ny)); hlim=0
    allocate(hsurf(nx,ny)); hsurf=0
    allocate(zpack(nx,ny)); zpack=0
    allocate(z0(nx,ny)); z0=0
    allocate(z0sn(nx,ny)); z0sn=0
    allocate(acc(nx,ny)); acc=0
    allocate(rain(nx,ny)); rain=0
    allocate(condt(nx,ny)); condt=0
    allocate(runoff(nx,ny)); runoff=0
    allocate(wetextent(nx,ny)); wetextent=0
    allocate(percozone(nx,ny)); percozone=0
    allocate(dryzone(nx,ny)); dryzone=0
    allocate(bareice(nx,ny)); bareice=0
    allocate(den(nx,ny)); den=0
    allocate(enst(nx,ny)); enst=0
    allocate(enfi(nx,ny)); enfi=0
    allocate(totmfi(nx,ny)); totmfi=0
    allocate(totmst(nx,ny)); totmst=0
    allocate(dm(nx,ny)); dm=0
    allocate(absden(nx,ny)); absden=0
    allocate(absdm(nx,ny)); absdm=0
    ! depends upon max number of annual layers
    allocate(nk(nx,ny,mnal)); nk=0
    allocate(zbal(nx,ny,mnal)); zbal=0
    allocate(ztal(nx,ny,mnal)); ztal=0 
    ! depends upon max number of vertical layers
    allocate(byear(nx,ny)); byear=0
    allocate(jyear(nx,ny,mnl)); jyear=0
    allocate(z(nx,ny,mnl)); z=0
    allocate(dz(nx,ny,mnl)); dz=0
    allocate(rmass(nx,ny,mnl)); rmass=0
    allocate(temp(nx,ny,mnl)); temp=0
    allocate(dens(nx,ny,mnl)); dens=0
    allocate(watcon(nx,ny,mnl)); watcon=0
    allocate(wetage(nx,ny,mnl)); wetage=0
    allocate(rppart(nx,ny,mnl)); rppart=0
    allocate(source(nx,ny,mnl)); source=0

    ! initialise some loop variables
    nl=0
    ist=0

    ! compute the position of the top and the bottom of each annual layer (initial situation)
    ! nal is the number of annual layers; the lowest layer is number one
    ! zbal is the depth of the bottom of each annual layer 
    ! ztal is the depth of the top of each annual layer 
    ! rint is the snow pack depth
  
    nal=2
    zbal(:,:,1)=zdeep
    zbal(:,:,2)=rint

      do x=1,nx
         do y=1,ny
            if (nal(x,y).ne.1) then
               do i=1,nal(x,y)-1
                  ztal(x,y,i)=zbal(x,y,i+1)
               enddo
            endif
            ztal(x,y,nal(x,y))=0.0
         enddo
      enddo

      ! make the initial grid
      ! in: dz0,coefdz,zbal,ztal,nal,mnal,mnl
      ! out: nk,nl,z,dz
      do x=1,nx
         do y=1,ny
            call grid (mnal,mnl,dz(x,y,:),z(x,y,:),nl(x,y),dz0,coefdz,                   &
                 zbal(x,y,:),ztal(x,y,:),nk(x,y,:),nal(x,y))
            ! out:temp, dens, watcon, rmass, wetage,surfwt,tempin
            call initiate (mnal,mnl,tempb,nk(x,y,:),temp(x,y,:),t0ini,z(x,y,:),rint,     &
                 dry,dens(x,y,:),denup,denice,jyear(x,y,:),byear(x,y),            &
                 dendw,epsden,wetage(x,y,:),nl(x,y),rmass(x,y,:),dz(x,y,:),              &
                 watcon(x,y,:),tempin(x,y),surfwt(x,y))
         enddo
      enddo

      ! the top of the ice at the beginning of the calculations is the reference height
      zref=0.0
      do x=1,nx
         do y=1,ny
            do i=nk(x,y,nal(x,y))+1,nl(x,y)
               zref(x,y)=zref(x,y)+dz(x,y,i)
            enddo
         enddo
      enddo

      ! the momentum roughness length of snow depends on whether the snow is wet or dry      
      do x=1,nx
         do y=1,ny
            if (watcon(x,y,1).eq.0.0) then
               z0sn(x,y)=z0fs
            else
               z0sn(x,y)=z0ms
            endif
            hmin(x,y)=height(mnl,dz(x,y,:),zref(x,y),nl(x,y))
         enddo
      enddo

    end subroutine SMBInit


    ! Main time loop subroutine
    subroutine SMBStep(mnl,mnal,dz0,zdeep,dzdeep, &
         mindt,albkind, &
         VanDus,kindz0,kstab, &
         alboldwet,albolddry,albfrs,albice,albsup,albwat,decalb,scswat, &
         scsnow,z0fs,z0ms,z0ice,abslow, &
         zwind,ztemp,ist,dt, &
         iyear, &
         x1,x2,y1,y2,ratec,cst,dro,coefdz,facstr, &
         bounl,bounh,nl,SLM,&
         rlhf,shf,terrad,swrtop, &
         fluxrn,enextr,tempin,t0,surfwt,frsnow,tsnwfall,alb,conden, &
         eatm,esat,denair,e0,frivel,topsl, &
         factin,dfact,zdum,zbtm,zref,hmin,hlim,hsurf,zpack,z0,z0sn, &
         acc,rain,condt,runoff,percozone,dryzone,bareice, &
         den,enst,enfi,totmfi,totmst,dm,absden,absdm, &
         nk,jyear,z,dz,rmass,temp,dens,watcon,wetage, &
         rppart,source, &
         tatm,prmass,wvel, &
         rhum,glr,atmrad,Psurf,totacc, &
         totrn,totcon,totoff,massbal, &
         rkarma)
      
      use smb_const
      use smb_debug
      
      implicit none
      
      ! Arguments
      logical,intent(in)    :: VanDus  !Conductivity as in Van Dusen,1929 (true) or Sturm, 1997

      integer,intent(in)    :: mnl     !MAXIMUM number of vertical layer
      integer,intent(in)    :: mnal    !MAXIMUM number of annual layer
      integer,intent(in)    :: mindt   !timestep in min
      integer,intent(in)    :: albkind !choose between 3 different albedo schemes 
      integer,intent(in)    :: kindz0  !roughness lengths (1:z0t=z0; 2:z0t=0.00001; 3:Brutsaert; 4:Andreas)
      integer,intent(in)    :: kstab   !stability functions(1=Businger; 2=Hogstrom; 3=Beljaars and Holtslag; 4=Dyer)
      integer,intent(inout) :: ist     ! loop stuff
      integer,intent(in)    :: iyear   ! current year
      integer,intent(in)    :: x1      ! coordinates target 
      integer,intent(in)    :: x2      ! coordinates target 
      integer,intent(in)    :: y1      ! coordinates target 
      integer,intent(in)    :: y2      ! coordinates target 

      real(rk),intent(in)   :: dz0         !dz0 is the desired grid point size at the surface
      real(rk),intent(in)   :: zdeep       !zdeep the depth of the model 
      real(rk),intent(in)   :: dzdeep      !size of lowermost grid cell
      real(rk),intent(in)   :: alboldwet   !albedo of old wet snow
      real(rk),intent(in)   :: albolddry   !albedo of old dry snow
      real(rk),intent(in)   :: albfrs      !albedo of fresh snow
      real(rk),intent(in)   :: albice      !albedo of ice
      real(rk),intent(in)   :: albsup      !albedo of superimposed ice
      real(rk),intent(in)   :: albwat      !albedo of water
      real(rk),intent(in)   :: decalb      !decay time for snow albedo
      real(rk),intent(in)   :: scswat      !scaling height (mm) for the effect of surficial water on the surface albedo 
      real(rk),intent(in)   :: scsnow      !scaling height for the effect of fresh snow
      real(rk),intent(in)   :: z0fs        !momentum roughness length of frozen snow (m)
      real(rk),intent(in)   :: z0ms        !momentum roughness length of melting snow (m)
      real(rk),intent(in)   :: z0ice       !momentum roughness length of ice (m)
      real(rk),intent(in)   :: abslow      !radiation penetration in snow (m-1)
      real(rk),intent(in)   :: zwind       !height of wind velocity input (m) for bulk method or highest level profile method
      real(rk),intent(in)   :: ztemp       !height of temperature input (m) for bulk method or highest level profile method
      real(rk),intent(in)   :: dt          ! timestep in sec
      real(rk),intent(in)   :: coefdz      !
      real(rk),intent(in)   :: facstr      !
      real(rk),intent(in)   :: bounl       !
      real(rk),intent(in)   :: bounh       !
      real(rk),intent(out)  :: ratec       ! use for dry snow densification
      real(rk),intent(out)  :: cst         ! use for dry snow densification
      real(rk),intent(out)  :: dro         ! use for dry snow densification
      real(rk),intent(out)  :: rkarma      ! von karman constant

      integer,intent(inout),dimension(:,:)  :: nl              !nb of layer
      real(rk),intent(in),dimension(:,:) :: SLM             !sea-land mask
      real(rk),dimension(:,:) :: tatm            !atmospheric input
      real(rk),dimension(:,:) :: wvel            !atmospheric input
      real(rk),dimension(:,:) :: prmass          !atmospheric input
      real(rk),dimension(:,:) :: glr             !atmospheric input
      real(rk),dimension(:,:) :: atmrad          !atmospheric input
      real(rk),dimension(:,:) :: rlhf            !energy related
      real(rk),dimension(:,:) :: shf             !energy related
      real(rk),dimension(:,:) :: terrad          !energy related
      real(rk),dimension(:,:) :: swrtop          !energy related
      real(rk),dimension(:,:) :: fluxrn          !energy related
      real(rk),dimension(:,:) :: enextr          !energy related
      real(rk),dimension(:,:) :: tempin          !temperature
      real(rk),dimension(:,:) :: t0              !temperature
      real(rk),dimension(:,:) :: surfwt          !albedo related
      real(rk),dimension(:,:) :: frsnow          !albedo related
      real(rk),dimension(:,:) :: tsnwfall        !albedo related
      real(rk),dimension(:,:) :: alb             !albedo related
      real(rk),dimension(:,:) :: conden          !Turb HF related
      real(rk),dimension(:,:) :: Psurf           !Turb HF related
      real(rk),dimension(:,:) :: eatm            !Turb HF related
      real(rk),dimension(:,:) :: esat            !Turb HF related
      real(rk),dimension(:,:) :: denair          !Turb HF related
      real(rk),dimension(:,:) :: e0              !Turb HF related
      real(rk),dimension(:,:) :: rhum            !Turb HF related
      real(rk),dimension(:,:) :: frivel          !Turb HF related
      real(rk),dimension(:,:) :: topsl           !runoff related
      real(rk),dimension(:,:) :: factin          !runoff related
      real(rk),dimension(:,:) :: dfact           !runoff related
      real(rk),dimension(:,:) :: zdum            !grid adjustment
      real(rk),dimension(:,:) :: zbtm            !grid adjustment
      real(rk),dimension(:,:) :: zref            !geometry
      real(rk),dimension(:,:) :: hmin            !geometry
      real(rk),dimension(:,:) :: hlim            !geometry
      real(rk),dimension(:,:) :: hsurf           !geometry
      real(rk),dimension(:,:) :: zpack           !geometry
      real(rk),dimension(:,:) :: z0              !geometry
      real(rk),dimension(:,:) :: z0sn            !geometry
      real(rk),dimension(:,:) :: acc             !mass balance
      real(rk),dimension(:,:) :: rain            !mass balance
      real(rk),dimension(:,:) :: condt           !mass balance
      real(rk),dimension(:,:) :: runoff          !mass balance
      real(rk),dimension(:,:) :: percozone       !defined zones
      real(rk),dimension(:,:) :: dryzone         !defined zones
      real(rk),dimension(:,:) :: bareice         !defined zones
      real(rk),dimension(:,:) :: den             !mass-energy cons.
      real(rk),dimension(:,:) :: enst            !mass-energy cons.
      real(rk),dimension(:,:) :: enfi            !mass-energy cons.
      real(rk),dimension(:,:) :: totmfi          !mass-energy cons.
      real(rk),dimension(:,:) :: totmst          !mass-energy cons.
      real(rk),dimension(:,:) :: dm              !mass-energy cons.
      real(rk),dimension(:,:) :: absden          !mass-energy cons.
      real(rk),dimension(:,:) :: absdm           !mass-energy cons.
      real(rk),dimension(:,:) :: totacc          !mass balance
      real(rk),dimension(:,:) :: totrn           !mass balance
      real(rk),dimension(:,:) :: totcon          !mass balance
      real(rk),dimension(:,:) :: totoff          !mass balance
      real(rk),dimension(:,:) :: massbal         !mass balance
      integer,dimension(:,:,:) :: nk          !Nb of nl per annual layer
      integer,dimension(:,:,:) :: jyear       !year for deposition of layer nk
      real(rk),dimension(:,:,:) :: z           !grid variables
      real(rk),dimension(:,:,:) :: dz          !grid variables
      real(rk),dimension(:,:,:) :: rmass       !englacial var.
      real(rk),dimension(:,:,:) :: temp        !englacial var.
      real(rk),dimension(:,:,:) :: dens        !englacial var.
      real(rk),dimension(:,:,:) :: watcon      !englacial var.
      real(rk),dimension(:,:,:) :: wetage      !englacial var.
      real(rk),dimension(:,:,:) :: rppart      !energy penetration
      real(rk),dimension(:,:,:) :: source      !energy penetration


      ! Local Variables
      integer :: i  ! loop stuff
      integer :: x  ! loop stuff
      integer :: y  ! loop stuff

      ! Reset the accumulators every time-step
      totoff=0.0
      totacc=0.0
      totrn=0.0
      totcon=0.0
     
      ist=ist+1
            
      do x=x1,x2
         do y=y1,y2

            if (SLM(x,y).ge.25.0) then

               ! if the grid becomes shorter than zdeep, another layer is added
               zbtm(x,y)=z(x,y,nl(x,y))+0.5*dz(x,y,nl(x,y))
               if (zbtm(x,y).lt.zdeep) then
                  nl(x,y)=nl(x,y)+1
                  dz(x,y,nl(x,y))=dzdeep
                  zref(x,y)=zref(x,y)+dzdeep
                  z(x,y,nl(x,y))=zbtm(x,y)+0.5*dz(x,y,nl(x,y))
                  temp(x,y,nl(x,y))=tempin(x,y)
                  watcon(x,y,nl(x,y))=0.0
                  dens(x,y,nl(x,y))=dens(x,y,nl(x,y)-1)
                  rmass(x,y,nl(x,y))=dens(x,y,nl(x,y))*dz(x,y,nl(x,y))
                  if (dens(x,y,nl(x,y)).lt.denice-epsden) then
                     wetage(x,y,nl(x,y))=0.0
                  else
                     wetage(x,y,nl(x,y))=-1.0   
                  endif
                  nk(x,y,1)=nk(x,y,1)+1
                  jyear(x,y,nl(x,y))=jyear(x,y,nl(x,y)-1)
                  enst(x,y)=enst(x,y)+cp*rmass(x,y,nl(x,y))*temp(x,y,nl(x,y))+                 &
                       watcon(x,y,nl(x,y))*fus
                  totmst(x,y)=totmst(x,y)+rmass(x,y,nl(x,y))+watcon(x,y,nl(x,y))
               endif
               
               ! if the lowest grid point is below zdeep, its temperature is put equal
               !    to tempin
               do i=nl(x,y),1,-1
                  if (z(x,y,i).gt.zdeep) then
                     enst(x,y)=enst(x,y)+cp*rmass(x,y,i)*(tempin(x,y)-temp(x,y,i))
                     temp(x,y,i)=tempin(x,y)
                  else
                     exit
                  endif
               enddo
               
               ! if the uppermost grid boxes gets too thin or too thick, the grid must be adapted
               !OUT:z,dz,nl,jyear,nk,temp,rmass,dens,watcon,wetage
               if (dz(x,y,1).lt.bounl.or.dz(x,y,1).gt.bounh) then
                  call chgrid (mnal,mnl,frsnow(x,y),cp,denice,fus,densn,dz0,coefdz,blc,         &
                       bhc,iyear,nl(x,y),nk(x,y,1:mnal),jyear(x,y,1:mnl),dz(x,y,1:mnl),         &
                       z(x,y,1:mnl),wetage(x,y,1:mnl),temp(x,y,1:mnl),                          &
                       dens(x,y,1:mnl),watcon(x,y,1:mnl),rmass(x,y,1:mnl))
               endif
               
               ! the grid is also adapted when one of the other grid boxes becomes too small      
               do i=2,nl(x,y)
                  if (dz(x,y,i).lt.bounl) then 
                     call chgrid (mnal,mnl,frsnow(x,y),cp,denice,fus,densn,dz0,coefdz,blc,       &
                          bhc,iyear,nl(x,y),nk(x,y,1:mnal),jyear(x,y,1:mnl),dz(x,y,1:mnl),       &
                          z(x,y,1:mnl),wetage(x,y,1:mnl),temp(x,y,1:mnl),                        &
                          dens(x,y,1:mnl),watcon(x,y,1:mnl),rmass(x,y,1:mnl))
                     exit
                  endif
               enddo
               
               ! THIS IS FOR ENERGY CONSERVATION CHECK
               enst(x,y)=toten(mnl,nl(x,y),cp,rmass(x,y,:),temp(x,y,:),                          &
                    watcon(x,y,:),fus,surfwt(x,y))
               ! THIS IS FOR MASS CONSERVATION CHECK
               totmst(x,y)=totmss(mnl,nl(x,y),rmass(x,y,:),watcon(x,y,:),surfwt(x,y))      
               
               ! PRECIPITATION
               ! whether precipitation falls as snow or rain depends on the temperature of
               ! the free atmosphere
               
               if (tatm(x,y).gt.2) then 
                  rain(x,y)=prmass(x,y)
                  acc(x,y)=0.0
               else
                  rain(x,y)=0.0
                  acc(x,y)=prmass(x,y)
               endif
               
               if (acc(x,y).gt.0.0) then
                  tsnwfall(x,y)=0.0                       ! reset time since last snowfall
               else
                  tsnwfall(x,y)=tsnwfall(x,y)+mindt
               endif

               totacc(x,y)=totacc(x,y)+acc(x,y)
               totrn(x,y)=totrn(x,y)+rain(x,y)
               fluxrn(x,y)=rain(x,y)/dt*cpwat*tatm(x,y)
               
               if (prmass(x,y).ne.0.0) then
                  ! changes in water content, mass, density and fresh-snow amount due to snow and rain fall     
                  watcon(x,y,1)=watcon(x,y,1)+rain(x,y)
                  enst(x,y)=enst(x,y)+rain(x,y)*fus
                  rmass(x,y,1)=rmass(x,y,1)+acc(x,y)
                  dz(x,y,1)=dz(x,y,1)+acc(x,y)/densn
                  dens(x,y,1)=rmass(x,y,1)/dz(x,y,1)
                  enst(x,y)=enst(x,y)+acc(x,y)*temp(x,y,1)*cp
                  frsnow(x,y)=frsnow(x,y)+acc(x,y)
               endif
               
               !  RADIATION PENETRATION and SURFACE ENERGY BALANCE
               
               ! calculate penetration of solar radiation
               ! OUT: rppart
               call strpen (nl(x,y),rppart(x,y,1:mnl),dens(x,y,1:mnl),                            &
                    denlow,abslow,facstr,dz(x,y,1:mnl),mnl)
               
               ! calculation of the momentum roughness length of snow   
               ! hsurf is the height of the surface relative to the height of the ice surface 
               ! at the beginning of the simulation
               hsurf(x,y)=height(mnl,dz(x,y,1:mnl),zref(x,y),nl(x,y))
               hlim(x,y)=hmin(x,y)+0.01
               ! if the surface is higher by 1 cm than before, fresh snow gives z0 its fresh-snow value
               if (hsurf(x,y).gt.hlim(x,y)) then
                  z0sn(x,y)=z0fs
                  hmin(x,y)=hsurf(x,y)
               elseif (hsurf(x,y).lt.hmin(x,y)) then
                  hmin(x,y)=hsurf(x,y)
               endif
               ! but this is replaced by the wet-snow value if there is water in the uppermost grid box
               if (watcon(x,y,1).gt.0.0) z0sn(x,y)=z0ms
               
               ! Surface albedo
               
               if (albkind.eq.1) then
                  !OUT: alb
                  call albedo1 (wetage(x,y,1),frsnow(x,y),densn,alb(x,y),alboldwet,                 &
                       albolddry,albfrs,albice,decalb,dens(x,y,1),denice,epsden,albsup,             &
                       albwat,surfwt(x,y),scswat,scsnow,tatm(x,y),tsnwfall(x,y))  
                  
               elseif (albkind.eq.2) then
                  !OUT:alb
                  call albedo2(tatm(x,y),tsnwfall(x,y),alb(x,y))
                  
               elseif (albkind.eq.3) then
                  !OUT:alb
                  call albedo3(tatm(x,y),tsnwfall(x,y),alb(x,y),surfwt(x,y),                        &
                       scswat,frsnow(x,y),densn,wetage(x,y,1),scsnow,dens(x,y,1),                   &
                       denice,epsden,albsup,alboldwet,albice,albfrs,albwat,decalb)
               endif
               
               ! Shortwave radiation at the surface
               swrtop(x,y)=glr(x,y)*(1.0-alb(x,y))
               
               ! source= energy input from the atmosphere at grid point i
               do i=1,nl(x,y)
                  source(x,y,i)=rppart(x,y,i)*swrtop(x,y)
               enddo
               
               ! Outgoing longwave
               ! surf. temp. (t0) and surf. water vapour pressure (e0)
               ! t0 defined when calling stemp. If t0>0 deg, then t0 is set back to 0 deg.
               t0(x,y)=stemp(temp(x,y,1),temp(x,y,2),dz(x,y,1),dz(x,y,2))
               e0(x,y)=satwvp(t0(x,y),1)
               
               terrad(x,y)=sb*(t0(x,y)+tkel)**4
               
               ! Calculate air density and roughness length for wind (used in turbhf)
               ! Roughness length depends now on snow density, as described in WG94
    
               denair(x,y)=Psurf(x,y)/(gascda*(tatm(x,y)+tkel))
               
               if (dens(x,y,1).lt.denz0) then
                  z0(x,y)=z0sn(x,y)
               else
                  z0(x,y)=z0sn(x,y)+(dens(x,y,1)-denz0)/(denice-denz0)*(z0ice-z0sn(x,y))
               endif
            
               !saturation water vapour pressure:
               esat(x,y)=satwvp(tatm(x,y),2)
               eatm(x,y)=rhum(x,y)*esat(x,y)
              
               ! Turbulent heat fluxes, using the BULK METHOD
               ! OUT: frivel,conden,shf,rlhf
               call turbhf (ist,eatm(x,y),tatm(x,y),e0(x,y),wvel(x,y),t0(x,y),zwind,               & 
                    ztemp,z0(x,y),denair(x,y),shf(x,y),conden(x,y),rlhf(x,y),                      &
                    kstab,kindz0,frivel(x,y),rkarma)
               
               ! Calculate extra energy due to condensation
               condt(x,y)=conden(x,y)*dt                    !condt is in kg m-2
               totcon(x,y)=totcon(x,y)+condt(x,y)
               if (t0(x,y).eq.0.0) then                     !==>MELTING surface: 
                  !                                         !condt (>0 or <0) is added to watcon
                  !                                         !if resulting watcon is still >0, then it represents
                  !                                         !the new water content of uppermost grid box
                  !                                         !if it is <0, it is removed from snow mass, and 
                  !                                         !reset to 0.
                  watcon(x,y,1)=watcon(x,y,1)+condt(x,y)
                  enextr(x,y)=condt(x,y)*fus                 ! in J m-2, 
                  if (watcon(x,y,1).lt.0.0) then
                     rmass(x,y,1)=rmass(x,y,1)+watcon(x,y,1)
                     dz(x,y,1)=rmass(x,y,1)/dens(x,y,1)
                     dens(x,y,1)=rmass(x,y,1)/dz(x,y,1)
                     enextr(x,y)=enextr(x,y)-watcon(x,y,1)*(fus-cp*temp(x,y,1))
                     rlhf(x,y)=rlhf(x,y)+watcon(x,y,1)/dt*fus
                     watcon(x,y,1)=0.0
                  endif
               else         
                  !                                          !==>FROZEN surface: 
                  !                                          !condt is either deposition or sublimation
                  !                                          !and is simply added/removed from snow mass 
                  rmass(x,y,1)=rmass(x,y,1)+condt(x,y)
                  dz(x,y,1)=rmass(x,y,1)/dens(x,y,1)
                  dens(x,y,1)=rmass(x,y,1)/dz(x,y,1)
                  enextr(x,y)=condt(x,y)*temp(x,y,1)*cp
               endif
               enst(x,y)=enst(x,y)+enextr(x,y)
               
               ! source(x,y,1)= energy at the upper grid point 
               source(x,y,1)=source(x,y,1)+atmrad(x,y)-terrad(x,y)+shf(x,y)+                        &
                    rlhf(x,y)+fluxrn(x,y)
               
               !For energy conservation check
               do i=1,nl(x,y)
                  enst(x,y)=enst(x,y)+source(x,y,i)*dt
               enddo
               
               !Subglacial processes:
               !OUT:temp,dens,watcon,rmass,frsnow,surfwt,runoff,topsl
               !    bareice,percozone,dryzone
               call entemp (dt,dz(x,y,1:mnl),z(x,y,1:mnl),nl(x,y),temp(x,y,1:mnl),                 &
                    dens(x,y,1:mnl),watcon(x,y,1:mnl),rmass(x,y,1:mnl),frsnow(x,y),cp,             &
                    denice,fus,denwat,runoff(x,y),source(x,y,1:mnl),bareice(x,y),topsl(x,y),       &
                    surfwt(x,y),factin(x,y),dfact(x,y),percozone(x,y),dryzone(x,y),mnl,VanDus)
               
               totoff(x,y)=totoff(x,y)+runoff(x,y)
               
               ! calculate the depth of snow pack (in m.we)
               zpack(x,y)=0.0
               do i=1,nl(x,y)
                  if (dens(x,y,i).lt.830.0) zpack(x,y)=zpack(x,y)+(dz(x,y,i)*dens(x,y,i))*1000.0
               enddo
               
               ! the grid is updated
               zdum(x,y)=0.0
               do i=1,nl(x,y)
                  z(x,y,i)=zdum(x,y)+0.5*dz(x,y,i)
                  zdum(x,y)=zdum(x,y)+dz(x,y,i)
               enddo
               
               ! wet snow age
               do i=1,nl(x,y)
                  if (watcon(x,y,i).gt.0.0.and.wetage(x,y,i).ne.-1.0)                              &
                       wetage(x,y,i)=wetage(x,y,i)+1.0/(24.0*(60.0/mindt))
               enddo
               
               ! densification of dry snow, as in GK94
               do i=1,nl(x,y)
                  if (dens(x,y,i).lt.550.0) then
                     ratec=11.0*exp(-1222.0/(temp(x,y,i)+tkel))
                     cst=1.0
                  else
                     ratec=575.0*exp(-2574.0/(temp(x,y,i)+tkel))
                     cst=0.5
                  endif
                  dro=ratec*(accden/1000.0)**cst*(denice-dens(x,y,i))/(ndyear*(60.0/mindt))
                  dens(x,y,i)=dens(x,y,i)+dro
               enddo
               
               ! ==BEGIN=DEBUG==
               ! ENERGY + MASS  CONSERVATION CHECK
               ! it is checked whether no energy or mass got lost
               if (chcons) then
                  enfi(x,y)=toten(mnl,nl(x,y),cp,rmass(x,y,1:mnl),temp(x,y,1:mnl),              &
                       watcon(x,y,1:mnl),fus,surfwt(x,y))+runoff(x,y)*fus
                  den(x,y)=enfi(x,y)-enst(x,y)
                  absden(x,y)=abs(den(x,y))
                  if (abs(den(x,y)).gt.1.0e-2) then
                     print 111,den(x,y)
111                  format (' you have no conservation of energy ',                            &
                          '(den = ',f12.4,' Joule; one of many possible ',                      &
                          'reasons is that your grid boxes are too small')
                     print*,'x,y',x,y
                     stop ! can change to return 
                  endif
                  
                  totmfi=totmss(mnl,nl(x,y),rmass(x,y,1:mnl),watcon(x,y,1:mnl),surfwt(x,y))
                  dm(x,y)=totmfi(x,y)+runoff(x,y)-(totmst(x,y)+condt(x,y)+prmass(x,y))
                  absdm(x,y)=abs(dm(x,y))
                  if (abs(dm(x,y)).gt.1.0e-4) then
                     print 112,dm(x,y)
112                  format (' you have no mass conservation (dm = ',f12.4,')')
                     print*,'x,y',x,y
                     stop ! can change to return 
                  endif
               endif
               ! ==END=DEBUG=====
               
            endif
         enddo   ! loop on y
      enddo   ! loop on x

      ! this time-step's mass balance
      massbal=totacc+totrn+totcon-totoff
      
    end subroutine SMBStep
    
    ! Cleanup routine
    subroutine SMBCleanup(nal,nl,tatm,wvel, &
         prmass,glr,atmrad,rlhf,shf, &
         terrad,swrtop,fluxrn,enextr,tempin,t0, &
         surfwt,frsnow,tsnwfall,alb,conden,Psurf, &
         eatm,esat,denair,e0,rhum,frivel,topsl, &
         slope,factin,dfact,zdum,zbtm, &
         zref,hmin,hlim,hsurf,zpack,z0,z0sn,acc, &
         rain,condt,runoff,wetextent,percozone, &
         dryzone,bareice,den,enst,enfi,totmfi, &
         totmst,dm,absden,absdm,nk,zbal,ztal, &
         byear,jyear,z,dz,rmass,temp,dens,watcon, &
         wetage,rppart,source)
      
      implicit none
      
      integer,pointer,dimension(:,:)    :: nal         !nb of annual layer
      integer,pointer,dimension(:,:)    :: nl          !nb of layer
      integer,pointer,dimension(:,:)    :: byear       !year at the beginning of the simulation

      real(rk),pointer,dimension(:,:)   :: tatm        !atmospheric input
      real(rk),pointer,dimension(:,:)   :: wvel        !atmospheric input
      real(rk),pointer,dimension(:,:)   :: prmass      !atmospheric input
      real(rk),pointer,dimension(:,:)   :: glr         !atmospheric input
      real(rk),pointer,dimension(:,:)   :: atmrad      !atmospheric input
      real(rk),pointer,dimension(:,:)   :: rlhf        !energy related
      real(rk),pointer,dimension(:,:)   :: shf         !energy related
      real(rk),pointer,dimension(:,:)   :: terrad      !energy related
      real(rk),pointer,dimension(:,:)   :: swrtop      !energy related
      real(rk),pointer,dimension(:,:)   :: fluxrn      !energy related
      real(rk),pointer,dimension(:,:)   :: enextr      !energy related
      real(rk),pointer,dimension(:,:)   :: tempin      !temperature
      real(rk),pointer,dimension(:,:)   :: t0          !temperature
      real(rk),pointer,dimension(:,:)   :: surfwt      !albedo related
      real(rk),pointer,dimension(:,:)   :: frsnow      !albedo related
      real(rk),pointer,dimension(:,:)   :: tsnwfall    !albedo related
      real(rk),pointer,dimension(:,:)   :: alb         !albedo related
      real(rk),pointer,dimension(:,:)   :: conden      !Turb HF related
      real(rk),pointer,dimension(:,:)   :: Psurf       !Turb HF related
      real(rk),pointer,dimension(:,:)   :: eatm        !Turb HF related
      real(rk),pointer,dimension(:,:)   :: esat        !Turb HF related
      real(rk),pointer,dimension(:,:)   :: denair      !Turb HF related
      real(rk),pointer,dimension(:,:)   :: e0          !Turb HF related
      real(rk),pointer,dimension(:,:)   :: rhum        !Turb HF related
      real(rk),pointer,dimension(:,:)   :: frivel      !Turb HF related
      real(rk),pointer,dimension(:,:)   :: topsl       !runoff related
      real(rk),pointer,dimension(:,:)   :: slope       !runoff related
      real(rk),pointer,dimension(:,:)   :: factin      !runoff related
      real(rk),pointer,dimension(:,:)   :: dfact       !runoff related
      real(rk),pointer,dimension(:,:)   :: zdum        !grid adjustment
      real(rk),pointer,dimension(:,:)   :: zbtm        !grid adjustment
      real(rk),pointer,dimension(:,:)   :: zref        !geometry
      real(rk),pointer,dimension(:,:)   :: hmin        !geometry
      real(rk),pointer,dimension(:,:)   :: hlim        !geometry
      real(rk),pointer,dimension(:,:)   :: hsurf       !geometry
      real(rk),pointer,dimension(:,:)   :: zpack       !geometry
      real(rk),pointer,dimension(:,:)   :: z0          !geometry
      real(rk),pointer,dimension(:,:)   :: z0sn        !geometry
      real(rk),pointer,dimension(:,:)   :: acc         !mass balance
      real(rk),pointer,dimension(:,:)   :: rain        !mass balance
      real(rk),pointer,dimension(:,:)   :: condt       !mass balance
      real(rk),pointer,dimension(:,:)   :: runoff      !mass balance
      real(rk),pointer,dimension(:,:)   :: wetextent   !defined zones
      real(rk),pointer,dimension(:,:)   :: percozone   !defined zones
      real(rk),pointer,dimension(:,:)   :: dryzone     !defined zones
      real(rk),pointer,dimension(:,:)   :: bareice     !defined zones
      real(rk),pointer,dimension(:,:)   :: den         !mass-energy cons.
      real(rk),pointer,dimension(:,:)   :: enst        !mass-energy cons.
      real(rk),pointer,dimension(:,:)   :: enfi        !mass-energy cons.
      real(rk),pointer,dimension(:,:)   :: totmfi      !mass-energy cons.
      real(rk),pointer,dimension(:,:)   :: totmst      !mass-energy cons.
      real(rk),pointer,dimension(:,:)   :: dm          !mass-energy cons.
      real(rk),pointer,dimension(:,:)   :: absden      !mass-energy cons.
      real(rk),pointer,dimension(:,:)   :: absdm       !mass-energy cons.

      integer,pointer,dimension(:,:,:)  :: nk          !Nb of nl per annual layer
      integer,pointer,dimension(:,:,:)  :: jyear       !year for deposition of layer nk

      real(rk),pointer,dimension(:,:,:) :: zbal        !depth of btm of annual lyr 
      real(rk),pointer,dimension(:,:,:) :: ztal        !depth of top of annual lyr 
      real(rk),pointer,dimension(:,:,:) :: z           !grid variables
      real(rk),pointer,dimension(:,:,:) :: dz          !grid variables
      real(rk),pointer,dimension(:,:,:) :: rmass       !englacial var.
      real(rk),pointer,dimension(:,:,:) :: temp        !englacial var.
      real(rk),pointer,dimension(:,:,:) :: dens        !englacial var.
      real(rk),pointer,dimension(:,:,:) :: watcon      !englacial var.
      real(rk),pointer,dimension(:,:,:) :: wetage      !englacial var.
      real(rk),pointer,dimension(:,:,:) :: rppart      !energy penetration
      real(rk),pointer,dimension(:,:,:) :: source      !energy penetration

      if(associated(nal)) deallocate(nal)
      if(associated(nl)) deallocate(nl)
      if(associated(tatm)) deallocate(tatm)
      if(associated(wvel)) deallocate(wvel)
      if(associated(prmass)) deallocate(prmass)
      if(associated(glr)) deallocate(glr)
      if(associated(atmrad)) deallocate(atmrad)
      if(associated(rlhf)) deallocate(rlhf)
      if(associated(shf)) deallocate(shf)
      if(associated(terrad)) deallocate(terrad)
      if(associated(swrtop)) deallocate(swrtop)
      if(associated(fluxrn)) deallocate(fluxrn)
      if(associated(enextr)) deallocate(enextr)
      if(associated(tempin)) deallocate(tempin)
      if(associated(t0)) deallocate(t0)
      if(associated(surfwt)) deallocate(surfwt)
      if(associated(frsnow)) deallocate(frsnow)
      if(associated(tsnwfall)) deallocate(tsnwfall)
      if(associated(alb)) deallocate(alb)
      if(associated(conden)) deallocate(conden)
      if(associated(Psurf)) deallocate(Psurf)
      if(associated(eatm)) deallocate(eatm)
      if(associated(esat)) deallocate(esat)
      if(associated(denair)) deallocate(denair)
      if(associated(e0)) deallocate(e0)
      if(associated(rhum)) deallocate(rhum)
      if(associated(frivel)) deallocate(frivel)
      if(associated(topsl)) deallocate(topsl)
      if(associated(slope)) deallocate(slope)
      if(associated(factin)) deallocate(factin)
      if(associated(dfact)) deallocate(dfact)
      if(associated(zdum)) deallocate(zdum)
      if(associated(zbtm)) deallocate(zbtm)
      if(associated(zref)) deallocate(zref)
      if(associated(hmin)) deallocate(hmin)
      if(associated(hlim)) deallocate(hlim)
      if(associated(hsurf)) deallocate(hsurf)
      if(associated(zpack)) deallocate(zpack)
      if(associated(z0)) deallocate(z0)
      if(associated(z0sn)) deallocate(z0sn)
      if(associated(acc)) deallocate(acc)
      if(associated(rain)) deallocate(rain)
      if(associated(condt)) deallocate(condt)
      if(associated(runoff)) deallocate(runoff)
      if(associated(wetextent)) deallocate(wetextent)
      if(associated(percozone)) deallocate(percozone)
      if(associated(dryzone)) deallocate(dryzone)
      if(associated(bareice)) deallocate(bareice)
      if(associated(den)) deallocate(den)
      if(associated(enst)) deallocate(enst)
      if(associated(enfi)) deallocate(enfi)
      if(associated(totmfi)) deallocate(totmfi)
      if(associated(totmst)) deallocate(totmst)
      if(associated(dm)) deallocate(dm)
      if(associated(absden)) deallocate(absden)
      if(associated(absdm)) deallocate(absdm)
      if(associated(nk)) deallocate(nk)
      if(associated(zbal)) deallocate(zbal)
      if(associated(ztal)) deallocate(ztal)
      if(associated(byear)) deallocate(byear)
      if(associated(jyear)) deallocate(jyear)
      if(associated(z)) deallocate(z)
      if(associated(dz)) deallocate(dz)
      if(associated(rmass)) deallocate(rmass)
      if(associated(temp)) deallocate(temp)
      if(associated(dens)) deallocate(dens)
      if(associated(watcon)) deallocate(watcon)
      if(associated(wetage)) deallocate(wetage)
      if(associated(rppart)) deallocate(rppart)
      if(associated(source)) deallocate(source)
      
    end subroutine SMBCleanup
    
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCC            SUBROUTINES                  CCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    !!
    subroutine SMBWriteReal3dToFile(outFile,outData)

      ! write a 3d array of real-valued data to file 

      implicit none

      character(len=100)         :: outFile ! file to write output to
      real(rk), dimension(:,:,:) :: outData ! 3d, real-valued data array

      ! locals
      integer :: j,iNumRecords,nx,ny
      integer :: ios = 1
      logical :: lConnected = .true.
      integer :: iUnit = 500

      ! determine extent of data array
      nx = size(outData,1)
      ny = size(outData,2)
      iNumRecords = size(outData,3)

      ! Check that unit number is available
      inquire(UNIT=iUnit,OPENED=lConnected)
      if (lConnected==.true.) then
         print *, "ERROR: unit already taken:", iUnit
         stop
      end if

      ! Open the file
      open (FILE=trim(outFile),UNIT=iUnit,IOSTAT=ios,STATUS='unknown', &
           ACCESS='direct',FORM='unformatted',RECL=nx*ny*4)
      if (ios/=0) then
         print *, "ERROR: opening file:", outFile
         stop
      end if

      ! Write the data
      do j=1,iNumRecords
         write(iUnit,rec=j) outData(:,:,j)
      enddo

      ! Close the file
      close(iUnit)

    end subroutine SMBWriteReal3dToFile

    !!
    subroutine SMBWriteReal2dToFile(outFile,outData)

      ! write a 2d array of real-valued data to file 

      implicit none

      character(len=100)       :: outFile ! file to write output to
      real(rk), dimension(:,:) :: outData ! 2d, real-valued data array

      ! locals
      integer :: nx,ny
      integer :: ios = 1
      logical :: lConnected = .true.
      integer :: iUnit = 500

      ! determine extent of data array
      nx = size(outData,1)
      ny = size(outData,2)

      ! Check that unit number is available
      inquire(UNIT=iUnit,OPENED=lConnected)
      if (lConnected==.true.) then
         print *, "ERROR: unit already taken:", iUnit
         stop
      end if

      ! Open the file
      open (FILE=trim(outFile),UNIT=iUnit,IOSTAT=ios,STATUS='unknown', &
           ACCESS='direct',FORM='unformatted',RECL=nx*ny*4)
      if (ios/=0) then
         print *, "ERROR: opening file:", outFile
         stop
      end if

      ! Write the data
      write(iUnit,rec=1) outData(:,:)

      ! Close the file
      close(iUnit)

    end subroutine SMBWriteReal2dToFile

    !!
    subroutine SMBReadParams(paramFile, &
         mnl,mnal,dz0,zdeep,dzdeep, &
         outmap,chcons,t0ini,tempb,rint,denup,dendw,dry, &
         alboldwet,albolddry,albfrs,albice,albsup,albwat, &
         decalb,scswat,scsnow,z0fs,z0ms,z0ice,absice, &
         abslow,tsc0,tsc1,tsc10,tscfac,zwind,ztemp, &
         albkind,VanDus,kindz0,kstab,supice)

      ! Read model parameters from file

      implicit none

      character(len=100) :: paramFile  ! file holding parameters

!      integer  :: nx         !Num grid pts in x-direction 
!      integer  :: ny         !Num grid pts in y-direction 
      integer  :: mnl        !MAXIMUM number of vertical layer
      integer  :: mnal       !MAXIMUM number of annual layer
!      integer  :: dx         !grid resolution (m)
      real(rk) :: dz0        !dz0 is the desired grid point size at the surface
      real(rk) :: zdeep      !zdeep the depth of the model 
      real(rk) :: dzdeep     !size of lowermost grid cell
!      integer  :: mindt      !timestep in min
      logical  :: outmap     !extensive output saving (includes saving map output, i.e. large files) 
      logical  :: chcons     !check for mass convervation (true) or not
      real(rk) :: t0ini      !Initial surface temperature (C)
      real(rk) :: tempb      !Temperature in the lower layer (C)
      real(rk) :: rint       !initial snow pack depth at elevation of measurements (m)
      real(rk) :: denup      !density of upper layer (kg/m3)
      real(rk) :: dendw      !density of lower layer (kg/m3)
      logical  :: dry        !initial snow layer dry (true) or wet (false)
      real(rk) :: alboldwet  !albedo of old wet snow
      real(rk) :: albolddry  !albedo of old dry snow
      real(rk) :: albfrs     !albedo of fresh snow
      real(rk) :: albice     !albedo of ice
      real(rk) :: albsup     !albedo of superimposed ice
      real(rk) :: albwat     !albedo of water
      real(rk) :: decalb     !decay time for snow albedo
      real(rk) :: scswat     !scaling height (mm) for the effect of surficial water on the surface albedo 
      real(rk) :: scsnow     !scaling height for the effect of fresh snow
      real(rk) :: z0fs       !momentum roughness length of frozen snow (m)
      real(rk) :: z0ms       !momentum roughness length of melting snow (m)
      real(rk) :: z0ice      !momentum roughness length of ice (m)
      real(rk) :: absice     !radiation penetration in ice (m-1)
      real(rk) :: abslow     !radiation penetration in snow (m-1)
      real(rk) :: tsc0       !time scale for runoff of surficial water on a horizontal surface (days) 
      real(rk) :: tsc1       !time scale for runoff of surficial water on a surface with slope of 1 degree (days)
      real(rk) :: tsc10      !time scale for runoff of surficial water on a steep surface (days) 
      real(rk) :: tscfac     !ratio of time scales of runoff below and above surface
      real(rk) :: zwind      !height of wind velocity input (m) for bulk method or highest level profile method
      real(rk) :: ztemp      !height of temperature input (m) for bulk method or highest level profile method
      integer  :: albkind    !choose between 3 different albedo schemes 
      logical  :: VanDus     !Conductivity as in Van Dusen,1929 (true) or Sturm, 1997
      integer  :: kindz0     !roughness lengths (1:z0t=z0; 2:z0t=0.00001; 3:Brutsaert; 4:Andreas)
      integer  :: kstab      !stability functions(1=Businger; 2=Hogstrom; 3=Beljaars and Holtslag; 4=Dyer)
      logical  :: supice     !run off according to Zuo and Oe (true) or all water on top of ice runs off (false)

      ! local variables
      logical :: lConnected = .true.
      integer :: ios = 1
      integer :: iUnit = 20

      ! Check that unit number is available
      inquire(UNIT=iUnit,OPENED=lConnected)
      if (lConnected==.true.) then
         print *, "ERROR: unit already taken"
         stop
      end if

      ! Open the file
      open (FILE=trim(paramFile),UNIT=iUnit,IOSTAT=ios,STATUS='old')
      if (ios/=0) then
         print *, "Error: opening file"
         stop
      end if

      ! Read the contents
      read (iUnit,*) mnl          
      read (iUnit,*) mnal     
      read (iUnit,*) dz0     
      read (iUnit,*) zdeep  
      read (iUnit,*) dzdeep
      read (iUnit,*) outmap    
      read (iUnit,*) chcons    
      read (iUnit,*) t0ini
      read (iUnit,*) tempb
      read (iUnit,*) rint
      read (iUnit,*) denup
      read (iUnit,*) dendw
      read (iUnit,*) dry
      read (iUnit,*) alboldwet
      read (iUnit,*) albolddry
      read (iUnit,*) albfrs
      read (iUnit,*) albice
      read (iUnit,*) albsup
      read (iUnit,*) albwat
      read (iUnit,*) decalb
      read (iUnit,*) scswat
      read (iUnit,*) scsnow
      read (iUnit,*) z0fs
      read (iUnit,*) z0ms
      read (iUnit,*) z0ice
      read (iUnit,*) absice
      read (iUnit,*) abslow
      read (iUnit,*) tsc0
      read (iUnit,*) tsc1
      read (iUnit,*) tsc10
      read (iUnit,*) tscfac
      read (iUnit,*) zwind
      read (iUnit,*) ztemp
      read (iUnit,*) albkind
      read (iUnit,*) VanDus
      read (iUnit,*) kindz0
      read (iUnit,*) kstab
      read (iUnit,*) supice
      close(iUnit)

    end subroutine SMBReadParams

    !!
    subroutine windcalc (U10m,V10m,wspeed)

      !calculate the wind speed from the 2 components U and V

      implicit none

      real(rk), intent(in), dimension(:,:)  :: U10m
      real(rk), intent(in), dimension(:,:)  :: V10m
      real(rk), intent(out), dimension(:,:) ::wspeed

      wspeed=sqrt(U10m**2+V10m**2)

    end subroutine windcalc

    !!
    subroutine SMBUpdateFromOrog(nx,ny,dx,orog,slope, &
         supice,tsc10,tsc0,tsc1,mindt,tscfac,dfact,factin)

      ! Update slope and runoff related stuff given a new
      ! orography

      implicit none

      ! Arguments
      integer, intent(in) :: nx           ! Number of grid point in x-direction
      integer, intent(in) :: ny           ! Number of grid point in y-direction
      integer, intent(in) :: dx           ! grid resolution (m)
      real(rk), intent(in), dimension(:,:) :: orog      !orography
      real(rk), intent(inout), dimension(:,:) :: slope  !runoff related
      logical, intent(in) :: supice       !run off according to Zuo and Oe (true) or all water on top of ice runs off (false)
      real(rk), intent(in) :: tsc0        !time scale for runoff of surficial water on a horizontal surface (days) 
      real(rk), intent(in) :: tsc1        !time scale for runoff of surficial water on a surface with slope of 1 degree (days)
      real(rk), intent(in) :: tsc10       !time scale for runoff of surficial water on a steep surface (days) 
      integer, intent(in) :: mindt        !timestep in min
      real(rk), intent(in) :: tscfac      !ratio of time scales of runoff below and above surface
      real(rk), intent(out), dimension(:,:) :: dfact    !runoff related
      real(rk),intent(out), dimension(:,:) :: factin    !runoff related
      ! Local variables
      integer :: x,y
      
      ! CALLED IN THE INITIALIZATION SO FAR.
      ! SHOULD BE CALLED IN THE TIMELOOP FOR ON-LINE RUNS
      
      !Compute the slope. 
      !OUT: slope
      call calc_slope(nx,ny,dx,orog,slope)
      !Compute the parameters necessary for the timescale of runoff 
      !OUT:dfact,factin
      do x=1,nx
         do y=1,ny
            call runfac (supice,tsc10,tsc0,tsc1,slope(x,y),mindt,                   &
                 tscfac,dfact(x,y),factin(x,y))
         enddo
      enddo
      
    end subroutine SMBUpdateFromOrog    

  !!
  ! this routines produces the initial grid
  subroutine grid(mnal,mnl,dz,z,nl,dz0,coefdz,zbal,ztal,nk,nal)
    
    implicit none

    ! arguments
    integer, intent(in) :: mnal,mnl
    real(rk), intent(out), dimension(mnl) :: dz,z
    integer, intent(out) :: nl
    real(rk), intent(in) :: dz0,coefdz
    real(rk), intent(in), dimension(mnal) :: zbal,ztal
    integer, intent(out), dimension(mnal) :: nk
    integer, intent(in) :: nal
    ! local variables
    integer :: i,j,k,ip
    real(rk) :: zdum,boun,dzk,factor
    real(rk), dimension(mnl) :: dzd

    j=0
    ip=0

    do i=nal,1,-1
       k=0
       zdum=ztal(i-ip)
       
       ! dzd is the desired thickness of each layer
       do
          k=k+1
          dzd(k)=dz0+coefdz*zdum
          zdum=zdum+dzd(k)
          ! we stop after traversing the bottom of each layer
          if (zdum.ge.zbal(i)) exit
       enddo

       ! set the number of grid points in the annual layer: nk
       boun=zdum-0.5*dzd(k)
       if (zbal(i).gt.boun) then
          nk(i)=k
          dzk=zdum-ztal(i-ip)
       else
          nk(i)=k-1
          dzk=zdum-dzd(k)-ztal(i-ip)
       endif

       ! ip is the number of annual layers that is too thin to form a single grid point
       if (nk(i).eq.0) then
          ip=ip+1
       else
          ! give the ultimate thickness to each grid point
          ! j is the total number of grid points until now
          factor=(zbal(i)-ztal(i-ip))/dzk
          
          zdum=ztal(i-ip)
          do k=1,nk(i)
             j=j+1
             dz(j)=dzd(k)*factor
             z(j)=zdum+0.5*dz(j)
             zdum=zdum+dz(j)
          enddo
          
          ip=0
       endif
    enddo

    nl=j
  
    return
  end subroutine grid


  ! this routine is called when the sub-surface grid is changed
  subroutine chgrid (mnal,mnl,frsnow,cp,denice,fus,densn,dz0,coefdz,blc,     &
       bhc,iyear,nl,nk,jyear,dz,z,wetage,temp,dens,watcon,rmass)

    implicit none
    
    ! Constants
    real(rk), parameter :: epsden=1.0
    ! Arguments
    integer, intent(in) :: mnal,mnl
    real(rk), intent(inout) :: frsnow  ! not sure...I think
    real(rk), intent(in) :: cp,denice,fus,densn
    real(rk), intent(in) :: dz0,coefdz,blc,bhc
    integer, intent(in) :: iyear
    integer, intent(inout) :: nl
    integer, intent(inout), dimension(mnal) :: nk
    integer, intent(inout), dimension(mnl) :: jyear
    real(rk), intent(inout), dimension(mnl) :: dz
    real(rk), intent(out), dimension(mnl) :: z
    real(rk), intent(inout), dimension(mnl) :: wetage,temp,dens,watcon,rmass
    ! Variables
    real(rk) :: dzdes,zdum,bounl,bounh
    real(rk), dimension(mnl) :: dznew,rmnew,tnew,dnew,wcnew,agenew
    integer, dimension(mnl) :: jynew
    real(rk) :: absbl,enptem,enpwat,enmin,entp,dzfrsn,dzold,rmold,denold
    integer :: i,j,k,ik,ial
    logical :: skiptag


    absbl=blc*dz0  !ABSOLUTE MIN THICKNESS

    zdum=0.0
    j=0

    do i=1,nl

       skiptag=.false.

       dzdes=dz0+coefdz*zdum
       z(i)=zdum+0.5*dz(i)
       zdum=zdum+dz(i)
       bounl=blc*dzdes    !blc=0.5 
       bounh=bhc*dzdes    !bhc=2
       
       ! fusion of 2 layers
       
       if (dz(i).lt.bounl.and.i.ne.nl) then
          
          ! under certain conditions no fusion should take place, that is 
          !    (provided the grid point is thicker than the absolute minimum)
          !    1) a grid point with a density below that of ice is positioned above
          !       a grid point having the density of ice
          !    2) between 2 grid points deposited during different years
          if (dz(i).ge.absbl.and.                                                     &
               ((dens(i).lt.denice-epsden.and.dens(i+1).ge.denice-epsden)             &
               .or.(jyear(i).ne.jyear(i+1)))) then          
             j=j+1
             dznew(j)=dz(i)
             rmnew(j)=rmass(i)
             tnew(j)=temp(i)
             dnew(j)=dens(i)
             wcnew(j)=watcon(i)
             jynew(j)=jyear(i)
             agenew(j)=wetage(i)
             skiptag=.true.
          endif
          
          if (.not.skiptag) then
             temp(i+1)=(rmass(i+1)*temp(i+1)+rmass(i)*temp(i))/(rmass(i+1)+rmass(i))
          
             if (wetage(i).ne.-1.0.and.wetage(i+1).ne.-1.0) then
                wetage(i+1)=(rmass(i+1)*wetage(i+1)+rmass(i)*wetage(i))/(rmass(i+1)+rmass(i))
             else
                wetage(i+1)=-1.0
             endif
             watcon(i+1)=watcon(i+1)+watcon(i)
             rmass(i+1)=rmass(i+1)+rmass(i)
          
             ! this for the case a thinning grid point at the surface is fused with an
             ! underlying grid point with the density of ice
             if (i.eq.1.and.dens(2).eq.denice) then
                dens(i+1)=denice
                dz(i+1)=rmass(i+1)/dens(i+1)
             else
                dz(i+1)=dz(i+1)+dz(i)
                dens(i+1)=rmass(i+1)/dz(i+1)
             endif
             
             if (watcon(i+1).gt.0.0.and.temp(i+1).lt.0.0) then
                enptem=abs(temp(i+1))*cp*rmass(i+1)
                enpwat=watcon(i+1)*fus
                enmin=min(enpwat,enptem)
                entp=temp(i+1)*cp*rmass(i+1)
                watcon(i+1)=watcon(i+1)-enmin/fus
                rmass(i+1)=rmass(i+1)+enmin/fus
                temp(i+1)=(entp+enmin)/(rmass(i+1)*cp)
             endif
          endif

             ! no change
             
       elseif (dz(i).lt.bounh) then
          j=j+1         
          dznew(j)=dz(i)
          rmnew(j)=rmass(i)
          tnew(j)=temp(i)
          dnew(j)=dens(i)
          wcnew(j)=watcon(i)
          jynew(j)=jyear(i)
          agenew(j)=wetage(i)
          
          ! layer is split into 2 parts
          
       else
          if (i.eq.1) then  
             dzfrsn=frsnow/densn
             dzold=dz(1)-dzfrsn
             rmold=rmass(1)-frsnow
             denold=rmold/dzold
             dznew(1)=dzfrsn
             
             if (dznew(1).lt.bounl) dznew(1)=bounl
             
             dznew(2)=dz(1)-dznew(1)
             
             if (dznew(2).ge.bounl) then         
                dnew(2)=denold
                rmnew(2)=dnew(2)*dznew(2)
                rmnew(1)=rmass(1)-rmnew(2)
                dnew(1)=rmnew(1)/dznew(1)
             else
                dznew(2)=bounl
                dznew(1)=dz(1)-dznew(2)
                rmnew(1)=dznew(1)*densn
                rmnew(2)=rmass(1)-rmnew(1)
                dnew(1)=densn
                dnew(2)=rmnew(2)/dznew(2)
             endif
             
             if (dnew(2).gt.denice) then
                dnew(2)=denice
                dznew(2)=rmnew(2)/dnew(2)
             endif
             
             wcnew(1)=(dznew(1)/dz(1))*watcon(1)
             wcnew(2)=(dznew(2)/dz(1))*watcon(1)
             
             tnew(1)=temp(1)
             tnew(2)=temp(1)
             
             if (frsnow.gt.0.0) then
                agenew(1)=0.0
             else
                agenew(1)=wetage(1)
             endif
             agenew(2)=wetage(1)
             
             frsnow=0.0
             j=2
             
          else
             k=j+1   
             do ik=k,k+1
                j=ik
                dznew(j)=0.5*dz(i)
                rmnew(j)=0.5*rmass(i)
                tnew(j)=temp(i)
                dnew(j)=dens(i)
                wcnew(j)=0.5*watcon(i)
                agenew(j)=wetage(i)
             enddo
          endif
          
          if (i.eq.1) then
             jynew(1)=iyear
             jynew(2)=jyear(1)
          else
             jynew(j-1)=jyear(i)
             jynew(j)=jyear(i)
          endif
       endif
       
    enddo
       
       nl=j
       
       zdum=0.0
       do i=1,nl
          dz(i)=dznew(i)
          z(i)=zdum+0.5*dz(i)
          zdum=zdum+dz(i)
          rmass(i)=rmnew(i)
          temp(i)=tnew(i)
          dens(i)=dnew(i)
          watcon(i)=wcnew(i)
          jyear(i)=jynew(i)
          wetage(i)=agenew(i)
       enddo
       
       do j=nl+1,mnl
          dz(j)=0.0
          z(j)=0.0
          rmass(j)=0.0
          temp(j)=0.0
          dens(j)=0.0
          watcon(j)=0.0
          wetage(j)=0.0
       enddo
       
       do ial=1,mnal
          nk(ial)=0
       enddo
       
       do i=nl,1,-1
          if (i.eq.nl) then
             ial=1
          else
             if (jyear(i).ne.jyear(i+1)) ial=ial+1
          endif
          nk(ial)=nk(ial)+1
       enddo
       
       return
     end subroutine chgrid


     subroutine initiate (mnal,mnl,tempb,nk,temp,t0ini,z,rint,dry,dens,denup,       &
          denice,jyear,byear,dendw,epsden,wetage,nl,rmass,dz,watcon,tempin,surfwt)
       
       !initiate the englacial variables: temp,dens,watcon,rmass,tempin
       ! plus variable related to presence of water: wetage and surfwt
       
       implicit none
       
       ! Arguments 
       integer, intent(in) :: mnal,mnl
       real(rk), intent(in) :: tempb
       integer, intent(in), dimension(mnal) :: nk
       real(rk), intent(out), dimension(mnl) :: temp
       real(rk), intent(in) :: t0ini
       real(rk), intent(in), dimension(mnl) :: z
       real(rk), intent(in) :: rint
       logical, intent(in) :: dry
       real(rk), intent(out), dimension(mnl) :: dens
       real(rk), intent(in) :: denup,denice
       integer, intent(out), dimension(mnl) :: jyear
       integer, intent(in) :: byear
       real(rk), intent(in) :: dendw,epsden
       real(rk), intent(out), dimension(mnl) :: wetage
       integer, intent(in) :: nl
       real(rk), intent(out), dimension(mnl) :: rmass,watcon
       real(rk), intent(in), dimension(mnl) :: dz
       real(rk), intent(out) :: tempin,surfwt
       
       ! Variables
       real(rk) :: dencap
       integer :: i
 
       tempin=tempb
       if (tempin.gt.0.0) tempin=0.0

       ! initialisation of the top layer	 
       do i=1,nk(2)
          temp(i)=t0ini+z(i)/rint*(tempin-t0ini)
          if (dry) then
             dens(i)=denup
             dencap=0.0
          else
             call inidensandwat (denup,denice,dencap,dens(i))
          endif
          
          jyear(i)=byear
          
          watcon(i)=dencap*dz(i)
          rmass(i)=dens(i)*dz(i)
       enddo
       
       ! initialisation of the 2nd layer, I think (nk(2) is the number of layer in upper annual layer)
       do i=nk(2)+1,nl
          temp(i)=tempin
          if (dry) then
             dens(i)=dendw
             dencap=0.0
          else
             call inidensandwat (dendw,denice,dencap,dens(i))
          endif
          jyear(i)=0
          watcon(i)=dencap*dz(i)
          rmass(i)=dens(i)*dz(i)
       enddo
       
       ! the wetage of ice is -1.0
       do i=1,nl
          if (dens(i).lt.denice-epsden) then
             wetage(i)=0.0
          else
             wetage(i)=-1.0   
          endif
       enddo
       
       surfwt=0.0
       
       return
     end subroutine initiate


    subroutine inidensandwat (denmea,denice,dencap,densol)

        ! computes density of capillary water and of solid material (snow/ice)
        ! from measured density (which includes water + snow/ice)
        ! the density of the capillary water is computed according to Coleou
        ! routine has been checked on 04-04-03      

      implicit none

      ! Constants
      real(rk), parameter :: c1=0.057
      real(rk), parameter :: c2=0.017
      ! Arguments
      real(rk), intent(in) :: denmea,denice
      real(rk), intent(out) :: dencap,densol
      ! Varaiables
      real(rk) :: a,b,c,determ

      if (denmea.lt.900.0) then
         a=1.0      
         b=(c1-c2-1.0)*denmea
         c=c1*denice*denmea+denmea**2*(c2-c1)
         determ=b**2-4.0*a*c
         if (determ.gt.0.0) then
            dencap=(-b-determ**0.5)/(2.0*a)
         else
            dencap=-b/(2.0*a)
         endif
      else
         dencap=0.0
      endif
      densol=denmea-dencap
      
      return
    end subroutine inidensandwat
            

    subroutine strpen (nl,rppart,dens,denlow,abslow,facstr,dz,mnl)

      ! routine for calculating penetration of short-wave radiation
      ! Follows GK94

      implicit none

      ! Constants
      real(rk), parameter :: eps=1.0e-5
      real(rk), parameter :: pabs=0.36
      ! Arguments
      integer, intent(in) :: nl,mnl
      real(rk), intent(out), dimension(mnl) :: rppart  
      real(rk), intent(in), dimension(mnl) :: dens,dz
      real(rk), intent(in) :: denlow,abslow,facstr
      ! Variables
      real(rk) :: abscf,swrbt,swrtp
      integer :: i

      do i=1,nl
         rppart(i)=0.0
      enddo

      swrtp=1.0-pabs                             !short wave radiation TOP layer
      do i=1,nl
         if (dens(i).lt.denlow) then
            rppart(i)=swrtp
            exit
         else
            abscf=abslow+(dens(i)-denlow)*facstr    !exctinction coeff
            swrbt=swrtp*exp(-abscf*dz(i))           !short wave radiation BOTTOM lyr
            rppart(i)=swrtp-swrbt                   !energy absorbed at layer i
            swrtp=swrbt                             !radiation top becomes rad bottom
            if (swrtp.lt.eps) then
               rppart(1)=rppart(1)+swrtp            !for energy conservation?
               exit
            endif
         endif
      enddo
      
      rppart(1)=rppart(1)+pabs
      
      return
    end subroutine strpen

    !!
    subroutine albedo1 (wetsur,frsnow,densn,alb,alboldwet,albolddry,albfrs,  &
         albice,decalb,densur,denice,epsden,albsup,albwat,      &
         surfwt,scswat,scsnow,tatm,tsnwfall)
      
      ! Albedo param is based on Greuell's original param
      
      ! the albedo is a function of the time that the uppermost grid point
      ! has been wet, given by wetsur (days)
      ! wetsur is negative if the uppermost grid point is ice
      ! for snow the albedo decreases exponentially with wetsur from albfrs
      ! to alboldwet with a time constant decalb (days)
      ! fresh snow (frsnow in mm w.e.), which has a depth dzfrsn (m), on top
      ! of the rest causes an extra increase in albedo (depth scale: scsnow (m))   

      implicit none

      ! Arguments
      real(rk), intent(in) :: wetsur,frsnow,densn
      real(rk), intent(out) :: alb
      real(rk), intent(in) :: alboldwet,albolddry,albfrs,albice,decalb
      real(rk), intent(in) :: densur,denice,epsden,albsup
      real(rk), intent(in) :: albwat,scswat,scsnow
      real(rk), intent(in) :: surfwt, tatm, tsnwfall
      ! Variables
      real(rk) :: dzfrsn

      if (wetsur.gt.0.0) then
         ! the following case must be superimposed ice
         if (densur.ge.denice-epsden) then
            alb=albsup
         else
            alb=alboldwet+(albfrs-alboldwet)*exp(-wetsur/decalb)
         endif
      elseif (wetsur.eq.0.0.and.tatm.gt.(-10)) then
         !decay for dry snow
         alb=albolddry+(albfrs-albolddry)*exp(-((tsnwfall)/(12.0*60.*24.0))*(tatm+10.0)/10.0)
      elseif (wetsur.eq.0.0.and.tatm.le.(-10)) then
         alb=albfrs   
      else   !(wetsur=-1 corresponds to ice)
         alb=albice     
      endif
      
      dzfrsn=frsnow/densn
      if (dzfrsn.gt.0.0)  alb=albfrs-(albfrs-alb)*exp(-dzfrsn/scsnow)
      if (surfwt.gt.0.0)  alb=albwat-(albwat-alb)*exp(-surfwt/scswat)         
      
      return
    end subroutine albedo1


    subroutine albedo2(tatm,tsnwfall,alb)

      implicit none
     
      ! Parameters
      real(rk), intent(in) :: tatm,tsnwfall
      real(rk), intent(out) :: alb

      if (tatm.le.-10.0) then
         alb=0.88
      elseif (tatm.lt.0.0.and.tatm.gt.-10.0) then
         alb=0.82-0.006*tatm
      elseif (tatm.ge.0.0.and.tatm.lt.8.0) then
         alb=0.82-0.03*tatm-1.74e-3*(tatm**2)-1.14e-4*(tatm**3)
      else
         alb=0.44
      endif
      
      if (tatm.lt.0.0) then
         alb=alb-0.0061*tsnwfall
      else
         alb=alb-0.0015*tsnwfall
      endif
      
      return
    end subroutine albedo2

    subroutine albedo3(tatm,tsnwfall,alb,surfwt,scswat,frsnow,densn,wetsur,scsnow,     &
         densur,denice,epsden,albsup,alboldwet,albice,albfrs,albwat,decalb)

      implicit none
      
      ! Arguments
      real(rk), intent(in) :: tatm,tsnwfall
      real(rk), intent(out) :: alb
      real(rk), intent(in) :: surfwt,scswat,frsnow,densn,wetsur,scsnow
      real(rk), intent(in) :: densur,denice,epsden
      real(rk), intent(in) :: albsup,alboldwet,albice,albfrs,albwat,decalb
      ! Variables
      integer:: reset
      real(rk) :: dzfrsn,critfrsn
      logical :: skiptag
      
      skiptag=.false.
      dzfrsn=frsnow/densn
      critfrsn=(1e-3*surfwt)       

      if (dzfrsn.gt.critfrsn.and.tatm.lt.2.0) then
         reset=1
      else
         reset=0
      endif

      if (tatm.le.-10.0) then
         alb=0.88
      elseif (tatm.lt.0.0.and.tatm.gt.-10.0) then
         alb=0.82-0.006*tatm
      elseif (tatm.ge.0.0.and.tatm.lt.8.0.and.surfwt.le.0.0) then
         alb=0.82-0.03*tatm-1.74e-3*(tatm**2)-1.14e-4*(tatm**3)
      elseif (tatm.ge.0.0.and.tatm.lt.8.0.and.surfwt.gt.0.0) then
           
         if (wetsur.ge.0.0) then
            ! the following case must be superimposed ice
            if (densur.ge.denice-epsden) then
               alb=albsup
            else
               alb=alboldwet+(albfrs-alboldwet)*exp(-wetsur/decalb)
            endif
         else
            alb=albice
         endif
          
         alb=albfrs-(albfrs-alb)*exp(-dzfrsn/scsnow)
         
         if (surfwt.gt.0.0.and.reset.eq.0) then
            alb=albwat-(albwat-alb)*exp(-surfwt/scswat)
            skiptag=.true.
         endif
         
      else
         alb=0.44
      endif

      if (.not.skiptag) then
         if (tatm.lt.0.0) then
            alb=alb-0.0061*tsnwfall
         else
            alb=alb-0.0015*tsnwfall
         endif
      endif

 
      return
    end subroutine albedo3


    subroutine turbhf (ist,eatm,tatm,e0,wvel,t0,zwind,ztemp,z0,denair,shf,           &
         conden,rlhf,kstab,kindz0,frivel,rkarma)
      
      implicit none
      
      ! Constants
      real(rk), parameter :: cpair=1005.0
      real(rk), parameter :: gamma=0.0098
      real(rk), parameter :: tkel=273.15
      real(rk), parameter :: grav=9.81
      real(rk), parameter :: vap=2.5e6
      real(rk), parameter :: crit=0.01
      real(rk), parameter :: fus=0.335e6
      real(rk), parameter :: visc=1.461e-5
      real(rk), parameter :: gascwv=461.89
      real(rk), parameter :: prand=0.71
      real(rk), parameter :: schmid=0.6
      integer, parameter :: miter=30
      ! Arguments
      integer, intent(in) :: ist
      real(rk), intent(in) :: eatm,tatm,e0,wvel,t0
      real(rk), intent(in) :: zwind,ztemp,z0,denair
      integer, intent(in) :: kstab,kindz0
      real(rk), intent(out) :: conden,rlhf,frivel,shf,rkarma
      ! Variables
      logical :: stop1
      real(rk) :: rmom1, rmom2
      integer :: iter
      real(rk) :: obukhl,tp,cap,obuold,denom,shfn,rlhfn
      real(rk) :: wind1,temp1,wvp1,zwind1,z0t,z0w,arg
      real(rk) :: reyst,arg1,arg2,ztemp1
      real(rk) :: sen1,sen2,rlat1,rlat2,scatem,scavpr
      real(rk) :: dshf,drlhf,zwvp1

      ! INPUT
      ! choose stability functions
      ! kstab = 1: Businger
      ! kstab = 2: Hogstrom          
      ! kstab = 3: Beljaars and Holtslag
      ! kstab = 4: Dyer
      
      ! kindz0 = 1: z0t=z0w=z0
      ! kindz0 = 2: z0t=z0w=0.00001
      ! kindz0 = 3: z0t and z0w computed according to Brutsaert
      ! kindz0 = 4: z0t and z0w computed according to Andreas

      ! ist: should be 1 the first time the routine is called
      !      At later call should be > 1
      ! denair: density of the air (kg/m3)

      ! eatm: water-vapour pressure (Pa) at height ztemp                
      ! tatm: temperature (ûC) at height ztemp                
      ! wvel: wind speed (m/s) at height zwind 
      ! zwind: height above surface of wvel (m)   GIVEN in in/inppar1            
      ! ztemp: height above surface of eatm and tatm (m)  GIVEN in in/inppar1

      ! Bulk method:               
      ! e0: water-vapour pressure (Pa) at height z0w                
      ! t0: temperature (ûC)  at height z0t
      ! z0: momentum roughness length (m)

      ! OUTPUT
      ! shf: sensible heat flux (W/m2)
      ! rlhf: latent heat flux (W/m2)
      ! conden: condensation (kg/m2/s)
      ! frivel: friction velocity (m/s)
      ! rkarma: 
      ! z0t: roughness length for temperature (m)   --only needed in this subroutine             
      ! z0w: roughness length for water vapour (m)  --only needed in this subroutine
       

      if (kstab.eq.1) then
         rkarma=0.35
      elseif (kstab.eq.2) then
         rkarma=0.40
      elseif (kstab.eq.3) then
         rkarma=0.40
      elseif (kstab.eq.4) then
         rkarma=0.41
      else
         rkarma=0.4
      endif

      ! count number of iteration steps
      iter=0
      stop1=.false.
      obukhl=0.0

      ! temperature is converted to potential temperature
      tp=tatm+gamma*ztemp
      
      ! the latent heat flux depends on whether water or snow (ice) evaporates
      !     (or condensates)
      if (t0.eq.0.0) then
         cap=vap
      else
         cap=vap+fus
      endif
      
      ! the first time the subroutine is called, the following variables
      !    are set equal to 0; at later calls the variables take their
      !    values calculated during the previous call
      if (ist.eq.1) then   
         frivel=0.0
         shf=0.0
         conden=0.0
      endif
      
      if (abs(shf).lt.crit.and.abs(rlhf).lt.crit) then
         frivel=rkarma*wvel/log(zwind/z0)
      endif
      
      do
         obuold=obukhl
         if (stop1) then
            obukhl=1.0e10
         else
            denom=rkarma*grav*(shf/((tatm+tkel)*cpair)-0.61*conden)
            if (denom.ne.0.0) then
               obukhl=frivel**3*denair/denom                          
            else
               obukhl=1.0e10
            endif
         endif
         ! the iteration does not converge if the Obukhov length switches from
         !    positive to negative and vice versa at each consequetive step      
         if (iter.gt.20) then
            if ((obukhl.gt.0.0.and.obuold.le.0.0).or.                      &
                 (obukhl.le.0.0.and.obuold.gt.0.0)) stop1=.true. 
         endif
         
         if (obukhl.eq.0.0) then
            shfn=0.0 
            conden=0.0
            rlhfn=0.0
         else
            ! heights and meteorological values of the lowest level are set
            wind1=0.0    
            temp1=t0      
            wvp1=e0     
            zwind1=z0       
            if (frivel.eq.0.0) then
               z0t=z0
               z0w=z0
            else
               ! there are 4 different ways to compute the roughness length for temperature
               !    and water vapour
               if (kindz0.eq.1) then
                  z0t=z0
                  z0w=z0
               elseif (kindz0.eq.2) then
                  z0t=1.0e-5
                  z0w=1.0e-5
               elseif (kindz0.eq.3) then
                  arg=7.3*(frivel*z0/visc)**0.25
                  z0t=z0*exp(-rkarma*(arg*prand**0.5-5.0))
                  z0w=z0*exp(-rkarma*(arg*schmid**0.5-5.0))
               elseif (kindz0.eq.4) then
                  reyst=frivel*z0/visc
                  arg1=0.317-0.565*log(reyst)-0.183*log(reyst)**2
                  z0t=0.001*exp(arg1)
                  arg2=0.396-0.512*log(reyst)-0.180*log(reyst)**2
                  z0w=0.001*exp(arg2)
               endif
            endif
            ztemp1=z0t    
            zwvp1=z0w     
            
            ! calculation of the primitive of the stability function divided by z
            !    at the upper and the lower level (1)
            if (obukhl.lt.0.0) then
               rmom1=insmom(zwind1,obukhl,kstab)
               rmom2=insmom(zwind,obukhl,kstab)
               sen1=inssen(ztemp1,obukhl,kstab)
               sen2=inssen(ztemp,obukhl,kstab)
               rlat1=inssen(zwvp1,obukhl,kstab)
               rlat2=inssen(ztemp,obukhl,kstab)
            else
               rmom1=stamom(zwind1,obukhl,kstab)
               rmom2=stamom(zwind,obukhl,kstab)
               sen1=stasen(ztemp1,obukhl,kstab)
               sen2=stasen(ztemp,obukhl,kstab)
               rlat1=stasen(zwvp1,obukhl,kstab)
               rlat2=stasen(ztemp,obukhl,kstab)
            endif
            
            ! frivel, scatem, scavp are velo (u*), temp(teta*), and humidity (q*) scale
            frivel=rkarma*(wvel-wind1)/(rmom2-rmom1)
            scatem=rkarma*(tp-temp1)/(sen2-sen1)
            scavpr=rkarma*(eatm-wvp1)/(rlat2-rlat1)
            ! sensible heat flux, condensation, latent heat flux:
            shfn=denair*cpair*frivel*scatem
            conden=frivel*scavpr/(gascwv*(tatm+tkel))
            rlhfn=conden*cap
         endif
         
         iter=iter+1

         
         dshf=shfn-shf      
         drlhf=rlhfn-rlhf 
         if (iter.le.miter) then
            shf=shfn
            rlhf=rlhfn
         else
            shf=shf+0.5*dshf
            rlhf=rlhf+0.5*drlhf
         endif
         
         ! the iteration continues when
         ! 1) the sensible or the latent heat flux does not converge within crit
         ! 2) less than miter (MAX number of iteration)  steps have been made
         if ((abs(dshf).gt.crit.or.abs(drlhf).gt.crit).and.                &
              iter.le.miter) then
            continue  ! i.e. go round the do loop again
         else
            exit  ! exit the do loop
         endif
      enddo
      
      return
    end subroutine turbhf

    !!
    !!Computes timescale for runoff factors: dfact and factin
    subroutine runfac (supice,tsc10,tsc0,tsc1,slope,mindt,          &
         tscfac,dfact,factin)
      
      implicit none
      ! Arguments
      logical, intent(in) :: supice
      real(rk), intent(in) :: tsc10,tsc0,tsc1,slope,tscfac
      integer, intent(in) :: mindt
      real(rk), intent(out) :: dfact,factin
      ! Variables
      real(rk) :: c1wat,c2wat,c3wat,tscex,factex
      
      if (supice) then
         c1wat=tsc10
         c2wat=tsc0-tsc10
         c3wat=-log((tsc1-c1wat)/c2wat)/tan(3.1415/180.0)
         tscex=c1wat+c2wat*exp(-c3wat*slope)
         factex=mindt/24.0/60.0/tscex
         if (factex.gt.1.0) factex=1.0
         factin=factex/tscfac
      else
         factex=1.0
         factin=1.0
      endif
      dfact=factex-factin
      
      return
    end subroutine runfac

    ! this routine computes changes in the subsurface temperature, water 
    !   content and density
    subroutine entemp (dt,dz,z,nl,temp,dens,watcon,rmass,frsnow,cp,denice,             &
         fus,denwat,runoff,source,bareice,topsl,surfwt,factin,dfact,                   &
         percozone,dryzone,mnl,VanDus)

      implicit none

      ! Arguments
      integer, intent(in) :: mnl,nl
      real(rk), intent(in) :: dt,cp,denice,fus,denwat
      real(rk), intent(in) :: factin,dfact
      real(rk), intent(out) :: bareice,runoff
      real(rk), intent(out) :: topsl,percozone,dryzone
      real(rk), intent(inout) :: surfwt,frsnow
      real(rk), intent(in), dimension(mnl) :: z,source
      real(rk), intent(inout), dimension(mnl) :: dz,temp,watcon,rmass,dens
      logical, intent(in) :: VanDus
      ! Variables
      real(rk) :: wcmax,roffs,roffin,remfw,dwat,wccap,entrem
      real(rk) :: dmass,enmin,enpden,enptem,enpwat,dzmelt
      real(rk) :: totcap,totfrt,totfri,totfw,totmel
      real(rk) :: dzsl,frsl
      real(rk), dimension(mnl) :: encon,cond
      logical :: skiptag1
      logical :: skiptag2
      integer :: i
      
      ! cond    conductivity of the grid points
      ! cp      specific heat capacity of ice
      ! ddens   difference between the density of ice and the actual density
      ! denice  density of ice
      ! dens    dry density of the grid points
      ! dfact   difference between factex and factin
      ! dmass   mass of the frozen water
      ! dryzone  =1 initially, assumes that there is no melting at all.
      !         If not true, will be set to 0.
      ! dwat    the difference in mass between the water in a grid point and
      !         the maximum possible amount of water in the grid point
      ! dzmelt  mass of melt water formed in a grid point
      ! dzsl    depth of slush layer in uppermost grid point with slush
      ! dt      time step (in seconds)
      ! encon   energy transferred from grid point i to grid point i+1
      !         by conduction
      ! enpden  the energy needed to fill freeze the water filling all the
      !         pore space
      ! enptem  the energy available in cold snow for freezing water
      ! enpwat  the energy that is necessary to freeze all the water
      ! entrem  energy remaining in negative temperature after freezing
      ! frsnow  mass of fresh snow
      ! fus     latent heat of fusion
      ! the factors factex and factin determine which part of the free water
      !         runs off in one time step. factex applies to the water on top of the
      !         surface. factin applies to the englacial water. Both factors are a
      !         function of the local slope (slope)
      ! nl      number of layers
      ! percozone  =1 initially, assumes that all water that melts will immediately
      !         refreeze (belongs to perco zone). If not true, will be set to 0.
      ! poros   porosity
      ! remfw   remaining amoun t of free water during the procedure of
      !         filling up all of the pore space
      ! rmass   mass of the snow in a grid point (water excluded)
      ! roffin  runoff within the snow pack
      ! roffs   runoff over the surface
      ! runoff  runoff
      ! slope   local slope (dimensionless; tangent of angle)
      ! source  energy input from the atmosphere at grid point i
      ! surfwt  mass of water on top of the ice (mm w.e.)
      ! temp    temperature of grid point i in degrees Celsius
      ! topsl   depth of the top of the slush layer
      ! totcap  total amount of capillary water over entire column
      ! totfri  total amount of internal freezing over the entire column
      ! totfrt  total amount of freezing on top of the surface    
      ! totfw   total mass of free water in the column (this excludes
      !         capillary water)
      ! totmel  total amount of melt over the entire column
      ! tscin   time scale for runoff of free water within the snow pack (in days)
      ! watcon  mass of water in grid point i
      ! wcmax   maximum possible mass of water in a grid point
      ! wcmaxp  wcmax given by the available pore space
      ! wcmaxc  wcmax according to equation of coleou which describes capillary water
      ! wmi     irreducible water content in % of mass according to coleou (mass of
      !         water devided by sum of masses of water and snow)
      ! z       depth below the surface of the centre of grid point i (m)   
      
      ! compute temperature change due to conduction and energy exchange
      !    with the atmosphere
      percozone=1
      dryzone=1
      bareice=0

      
      do i=1,nl
         cond(i)=conduc(dens(i),VanDus)
      enddo

      do i=1,nl-1
         encon(i)=0.5*(cond(i)+cond(i+1))*(temp(i)-temp(i+1))/             &
              (z(i+1)-z(i))
      enddo
      encon(nl)=0.0
      
      temp(1)=temp(1)+(source(1)-encon(1))*dt/(cp*rmass(1))
      do i=2,nl
         temp(i)=temp(i)+(source(i)+encon(i-1)-encon(i))*dt/(cp*rmass(i))
      enddo
      
      ! melt              
      totmel=0.0
      do i=1,nl
         if (temp(i).gt.0.0) then
            dzmelt=temp(i)*rmass(i)*cp/fus
            watcon(i)=watcon(i)+dzmelt
            rmass(i)=rmass(i)-dzmelt
            dz(i)=rmass(i)/dens(i)
            temp(i)=0.0
            totmel=totmel+dzmelt
            dryzone=0
            ! compute change to the mass of fresh snow
            if (i.eq.1) then
               frsnow=frsnow-dzmelt
               if (frsnow.lt.0.0) frsnow=0.0
            endif
            
         endif
      enddo
      
      ! freezing of water, percolation and retention by capillary forces  (do-loop 11)
      ! water on top of the ice is added to the water content of the uppermost grid point
      watcon(1)=watcon(1)+surfwt
      totfw=0.0
      
      totfri=0.0
      totfrt=0.0
      totcap=0.0
      
      do i=1,nl
         if (watcon(i).gt.0.0) then
            ! compute how much water freezes
            ! the amount of refreezing is limited by:
            ! either temp cannot be raised above melting point (ENPTEM)
            ! or available amount of water (ENPWAT)
            ! or available pore space (ENPDEN)
            if (temp(i).eq.0.0.or.dens(i).eq.denice) then
               continue
            else
               enpwat=watcon(i)*fus
               enptem=abs(temp(i))*rmass(i)*cp
               enpden=(denice-dens(i))*dz(i)*fus
               enmin=min(enpwat,enptem,enpden)
               dmass=enmin/fus        !mass of frozen water
               
               ! compute changes in variables due to freezing      
               watcon(i)=watcon(i)-dmass
               rmass(i)=rmass(i)+dmass
               dens(i)=rmass(i)/dz(i)
               entrem=enptem-enmin
               temp(i)=-1.0*entrem/(rmass(i)*cp)
               totfri=totfri+dmass
            endif

            ! surficial water may freeze to the top layer
            if (i.eq.1.and.temp(1).lt.0.0.and.watcon(1).gt.0.0) then
               enpwat=watcon(i)*fus
               enptem=abs(temp(i))*rmass(i)*cp
               enmin=min(enpwat,enptem)
               dmass=enmin/fus
               
               watcon(i)=watcon(i)-dmass
               rmass(i)=rmass(i)+dmass
	       
               dz(i)=dz(i)+dmass/denice
               dens(i)=rmass(i)/dz(i)
               entrem=enptem-enmin
               temp(i)=-1.0*entrem/(rmass(i)*cp)
               totfrt=totfrt+dmass
            endif
            
            ! percolation and retention of water
            ! compute the maximum amount of water that can be retained in the grid
            !    point
            wccap=dencapfun(dens(i),denice,denwat)*dz(i)
            
            ! remaining water percolates downwards and, if it is the lowest grid
            !    point and cannot be retained by capillary forces, it is free water
            if (watcon(i).gt.wccap) then
               dwat=watcon(i)-wccap
               if (i.ne.nl) then
                  watcon(i+1)=watcon(i+1)+dwat
               else
                  totfw=dwat
               endif
               watcon(i)=watcon(i)-dwat
            endif
            totcap=totcap+watcon(i)
         endif
      enddo

      ! If all the water has frozen again, call it part of the perco zone
      ! as defined by van der veen, p.353
      do i=1,nl   
         if (watcon(i).gt.0.0) then
            percozone=0
            exit
         endif
      enddo

      ! Check if the surface is bare-ice or not (use pore-closure density=830 kg/m3)
      if (dens(1).gt.830.0) bareice=1

      skiptag1=.false.
      skiptag2=.false.
      
      if (totfw.eq.0.0) then
         runoff=0.0
         topsl=-999.0
         surfwt=0.0
         skiptag1=.true.
      endif

      
      if (.not.skiptag1) then
         ! part of the free water runs off
         roffin=totfw*factin
         totfw=totfw-roffin
         remfw=totfw
      
         if (totfw.eq.0.0) then
            runoff=roffin
            topsl=-999.0
            surfwt=0.0
            skiptag2=.true.
         endif

      ! the remainder fills all of the pores
         if (.not.skiptag2) then
            do i=nl,1,-1
         
               wcmax=(denice-dens(i))/denice*dz(i)*denwat
               watcon(i)=watcon(i)+remfw
            
               if (watcon(i).gt.wcmax) then
                  remfw=watcon(i)-wcmax
                  watcon(i)=wcmax
               
                  ! when all the pores are filled, the rest is put on top of the surface
                  ! there it partly runs off (more easily than within the snow pack)
                  if (i.eq.1) then
                     surfwt=remfw
                     totfw=totfw-surfwt
                     topsl=0.0
                     roffs=surfwt*dfact
                     surfwt=surfwt-roffs
                  endif
               
               else
                  ! compute the position of the top of the slush layer
                  surfwt=0.0
                  wccap=dencapfun(dens(i),denice,denwat)*dz(i)
                  if (wccap.lt.wcmax) then
                     frsl=(watcon(i)-wccap)/(wcmax-wccap)
                  else
                     frsl=watcon(i)/wcmax
                  endif
                  dzsl=frsl*dz(i)
                  topsl=z(i)+0.5*dz(i)-dzsl
                  roffs=0.0
                  exit
               endif
            enddo
                  
            runoff=roffin+roffs
         endif
      endif
      
       
     
      return
    end subroutine entemp
    
    !!
    !SLOPE CALCULATIONS (Budd and Warner, 1996)
    subroutine calc_slope(nx,ny,dx,H,slope)

      implicit none

      ! Arguments
      integer, intent(in) :: nx,ny,dx
      real(rk), intent(in), dimension(nx,ny) :: H
      real(rk), intent(out), dimension(nx,ny) :: slope
      ! Variables
      integer :: x,y,nx1,ny1
      real(rk), dimension(nx,ny) :: sx,sy

      nx1=nx-1
      ny1=ny-1
      
      sx=0.0
      sy=0.0
      
      do x=2,nx1
         do y=2,ny1
            sx(x,y)=(H(x-1,y)-H(x+1,y))/(2*dx)
            sy(x,y)=(H(x,y-1)-H(x,y+1))/(2*dx)
         enddo
      enddo

      slope=sqrt(sx**2+sy**2)
      
      return
    end subroutine calc_slope

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccc!              FUNCTIONS                    ccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    real(rk) function height(mnl,dz,zref,nl)
      
      implicit none

      integer, intent(in) :: mnl,nl
      real(rk), intent(in), dimension(mnl) :: dz
      real(rk), intent(in) :: zref

      ! local variables
      integer :: i
      
      height=0.0
      do i=1,nl
         height=height+dz(i)
      end do
      height=height-zref

      return
    end function height

    real(rk) function dencapfun(densit,denice,denwat)

      implicit none

      ! calculation of the density of the amount of capillary water (kg/m3)
      !     according to Coleou
      ! has been tested on 03-04-03 and was ok for all input densities
      
      ! dencapfun  = computed density of capillary water (kg/m3)
      ! dencol  = maximum amount of density of capillary water according to Coleou (kg/m3)
      ! denice  = density of pure ice (kg/m3)
      ! denpor  = density of water when pores would be filled entirely (kg/m3)
      ! densit  = density of the dry snow in the grid box (so water is excluded) (kg/m3)
      ! denwat  = density of water (kg/m3)
      ! poros   = porosity

      ! Arguments
      real(rk), intent(in) :: densit,denice,denwat
      ! Variables
      real(rk) :: ddens,poros,wmi,dencol,denpor
        
      ddens=denice-densit
      if (abs(ddens).lt.1.0) then
         dencapfun=0.0
      else
         poros=ddens/denice
         wmi=0.057*poros/(1.0-poros)+0.017
         dencol=wmi/(1.0-wmi)*densit
         denpor=poros*denwat
         dencapfun=min(dencol,denpor)
      endif
      
      return
    end function dencapfun
    
    !!
    real(rk) function satwvp (t,kind) 

      implicit none

      real(rk), intent(in) :: t
      integer, intent(in) :: kind
      
      ! kind=1: saturation vapour pressure for ice
      ! kind=2: saturation vapour pressure for water
      
      if (kind.eq.1) then
         satwvp=(6.1091+t*(5.0347E-01+t*(1.8860E-02                        &
              +t*(4.1762E-04+t*(5.8247E-06+t*(4.8388E-08+                  &
           1.8388E-10*t))))))*100.0
      else
         satwvp=(6.1078+t*(4.4365E-01+t*(1.4289E-02                        &
              +t*(2.6506E-04+t*(3.0312E-06+t*(2.0341E-08+                  &
              6.1368E-11*t))))))*100.0
      endif
      
      return
    end function satwvp

    !!
    real(rk) function stemp (temp1,temp2,dz1,dz2)

      implicit none

      ! Arguments
      real(rk), intent(in) :: temp1,temp2,dz1,dz2
      ! Variables
      real(rk) :: tgrad
      
      tgrad=2.0*(temp2-temp1)/(dz1+dz2)
      stemp=temp1-tgrad*0.5*dz1
      
      if (stemp.gt.0.0) stemp=0.0
      
      return
    end function stemp


    real(rk) function conduc(dens,VanDus)

      implicit none
      !Argument"
      logical, intent(in)::VanDus
      !Variable
      real(rk), intent(in) :: dens
  
      
      if (VanDus) then
    ! the following function computes conductivity as a function of
    !   the density of the snow or ice: VON DUSEN'S APPROACH (1929)
         conduc=0.21e-01+0.42e-03*dens+0.22e-08*dens**3
      else
    ! Other possibility: equation of on Sturm et al. 1997,based on 488 measurments
    ! Also used in Bassford thesis
         conduc=0.138-1.01e-03*dens+3.233e-06*dens**2
      endif

      return
    end function conduc


    ! NEXT 4 FUNCTIONS: calculate the primitive of stability function,
    ! insmom and inssen: for unstable stratification (L<0)
    ! stamom and stasen: for stable stratificatin (L>0)
    !!
    real(rk) function insmom(z,obl,kstab)

      implicit none

      ! Arguments
      real(rk), intent(in) :: z, obl
      integer, intent(in) :: kstab
      ! Variables
      real(rk) :: fiw,coeff
        
      if (kstab.eq.1) then
         coeff=15.0
      elseif (kstab.eq.2) then
	 coeff=19.0
      elseif (kstab.eq.3) then
	 coeff=19.0
      elseif (kstab.eq.4) then
	 coeff=16.0
      endif
      
      fiw=(1.0-coeff*z/obl)**0.25
      insmom=log(z)-2.*log((1.+fiw)/2.)-log((1.+fiw**2)/2.)             &
           +2.0*atan(fiw)-3.1415/2.0
      
      return
    end function insmom

    !!
    real(rk) function inssen(z,obl,kstab)
        
      implicit none
      
      ! Arguments
      real(rk), intent(in) :: z,obl
      integer, intent(in) :: kstab
      ! Variables
      real(rk) :: fit,coeff

      if (kstab.eq.1) then
         coeff=9.0
      elseif (kstab.eq.2) then
         coeff=11.6
      elseif (kstab.eq.3) then
         coeff=11.6
      elseif (kstab.eq.4) then
         coeff=16.0
      endif
      
      fit=(1.0-coeff*z/obl)**0.5
      inssen=log(z)-2.*log((1.+fit)/2.)
      
      return
    end function inssen
    
    !!
    real(rk) function stamom(z,obl,kstab)
      
      implicit none
     
      ! Constants
      real(rk), parameter :: a=1.0,b=0.667,c=5.0,d=0.35
      ! Arguments
      real(rk), intent(in) :: z,obl
      integer, intent(in) :: kstab
      
      if (kstab.eq.1) then
         stamom=log(z)+4.7*z/obl
      elseif (kstab.eq.2) then
         stamom=log(z)+6.0*z/obl
      elseif (kstab.eq.3) then
         stamom=log(z)+a*z/obl-b/d*(1.0+c)*exp(-d*z/obl)-               &
               b*d/obl**2*(-z*obl/d-(obl/d)**2)*exp(-d*z/obl)
      elseif (kstab.eq.4) then
         stamom=log(z)+5.0*z/obl
      endif

      return
    end function stamom
      
    !!
    real(rk) function stasen(z,obl,kstab)
     
      implicit none

      ! Constants
      real(rk), parameter :: a=1.0,b=0.667,c=5.0,d=0.35
      ! Arguments
      real(rk), intent(in) :: z,obl
      integer, intent(in) :: kstab

      if (kstab.eq.1) then
         stasen=0.74*log(z)+4.7*z/obl
      elseif (kstab.eq.2) then
         stasen=0.95*log(z)+7.8*z/obl
      elseif (kstab.eq.3) then
         stasen=log(z)+(1.0+0.6667*a*z/obl)**1.5                        &
             -b/d*(1.0+c)*exp(-d*z/obl)-                                &
               b*d/obl**2*(-z*obl/d-(obl/d)**2)*exp(-d*z/obl)           
      elseif (kstab.eq.4) then
         stasen=log(z)+5.0*z/obl
      endif

      return
    end function stasen

    !!
    real(rk) function toten(mnl,nl,cp,rmass,temp,watcon,fus,surfwt)

      implicit none

      integer, intent(in) :: mnl,nl
      real(rk), intent(in) :: cp,fus,surfwt
      real(rk), intent(in), dimension(mnl) :: temp,watcon,rmass

      ! local variables
      integer :: i

      toten=surfwt*fus
      do i=1,nl
         toten=toten+cp*rmass(i)*temp(i)+watcon(i)*fus
      end do
      
      return
    end function toten

    !! 
    real(rk) function totmss(mnl,nl,rmass,watcon,surfwt)
      
      implicit none

      integer, intent(in) :: mnl, nl
      real(rk), intent(in), dimension(mnl) :: rmass,watcon
      real(rk), intent(in) :: surfwt

      ! local variables
      integer :: i
 
      totmss=surfwt
      do i=1,nl
         totmss=totmss+rmass(i)+watcon(i)
      end do

      return
    end function totmss
 
  end module smb_mecons
