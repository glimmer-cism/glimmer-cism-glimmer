module smb_force

  use glimmer_global

  implicit none

  ! Constants relating to forcing
  integer, parameter :: nt1     = 1460    !Number of record in 6-hourly input file (normal year)
  integer, parameter :: nt2     = 1464    !Number of record in 6-hourly input file (leap year)
  ! mout=number of hours of the experiment (8754 for 1 yr minus 6h
  ! -the last record is in fact 31 dec. 1991, at 18h00, in ERA40)
  integer, parameter :: mout1   = nt1*6
  integer, parameter :: mout2   = nt2*6
  integer, parameter :: xtarg   = 16          !coordinate target for single cell calculation
  integer, parameter :: ytarg   = 33         !coordinate target
  ! Variables not in the derived type
  character(len=36)  :: commpath          !find file locatio
  character(len=100) :: file1             !input filename of current yr
  character(len=100) :: file2             !input filename of current yr
  character(len=100) :: file3             !input filename of current yr
  character(len=100) :: file4             !input filename of current yr
  character(len=100) :: file5             !input filename of current yr
  character(len=100) :: file6             !input filename of current yr
  character(len=100) :: file7             !input filename of current yr
  character(len=100) :: file11            !input filename of follow. yr
  character(len=100) :: file22            !input filename of follow. yr
  character(len=100) :: file33            !input filename of follow. yr
  character(len=100) :: file44            !input filename of follow. yr
  character(len=100) :: file55            !input filename of follow. yr
  character(len=100) :: file66            !input filename of follow. yr
  character(len=100) :: file77            !input filename of follow. yr
  character(len=100) :: tmpFile
  !
  real(rk),pointer,dimension(:,:) :: orog  => null()           !orography
  integer,pointer,dimension(:,:)  :: drain => null()           !map of drainzones
  real(rk),pointer,dimension(:,:) :: SLM   => null()           !sea-land mask 
  integer  :: nin                 ! counters for interpolation
  integer  :: rinter              ! counters for interpolation
  real(rk),pointer,dimension(:,:,:) :: Ta => null()            !input 
  real(rk),pointer,dimension(:,:,:) :: prec => null()          !input 
  real(rk),pointer,dimension(:,:,:) :: wspeed => null()        !input 
  real(rk),pointer,dimension(:,:,:) :: Tdew => null()          !input 
  real(rk),pointer,dimension(:,:,:) :: SWdown => null()        !input
  real(rk),pointer,dimension(:,:,:) :: LWdown => null()        !input
  real(rk),pointer,dimension(:,:,:) :: press => null()       !input
  !
  real(rk),pointer,dimension(:)     :: quickocean => null()      !runoff location+volume
  real(rk),pointer,dimension(:,:)   :: quickoceantime => null()  !history of runoff loc. + vol.        
  real(rk),pointer,dimension(:,:,:) :: runoffmap => null()       !daily map output
  real(rk),pointer,dimension(:,:,:) :: albmap => null()          !daily map output
  real(rk),pointer,dimension(:,:,:) :: massbalmap => null()      !daily map output
  real(rk),pointer,dimension(:,:,:) :: drymap => null()          !daily map output
  real(rk),pointer,dimension(:,:)   :: dayroff => null()         !runoff related
  ! for output
  character(len=4) :: chyear              !current yr
  character(len=4) :: chyear2             !follow. yr
  ! Called in file 'modsetup':
  integer  :: whichtest = 1        !Test value (1: run for one point, 2: run for the whole grid)
  integer  :: count               ! other time counters
  integer  :: daycount            ! daily counter
  integer  :: daycount1           !
  integer  :: day                 ! other time counters
  integer  :: timecall            ! other time counters
  integer  :: ntotal              ! related to lenght of current year
  integer  :: yr                  ! current year
  integer  :: qoceansize = 6      ! see 'trackrunoff'

contains

  subroutine SMBInitForce(nx,ny,mnl,mnal,dx,dz0, &
       zdeep,dzdeep,mindt,outmap, &
       albkind,VanDus,kindz0,kstab,supice, &
       t0ini,tempb,rint,denup, &
       dendw,dry,alboldwet,albolddry,albfrs,albice,albsup,albwat,decalb,scswat, &
       scsnow,z0fs,z0ms,z0ice,absice,abslow,tsc0,tsc1,tsc10,tscfac, &
       zwind,ztemp,ist,dt, &
       x1,x2,y1,y2,coefdz,facstr,bounl,bounh, &
       nal,nl,byear,tatm,wvel,&
       prmass,glr,atmrad,rlhf,shf, &
       terrad,swrtop,fluxrn,enextr,tempin,t0, &
       surfwt,frsnow,tsnwfall,alb,conden,Psurf, &
       eatm,esat,denair,e0,rhum,frivel,topsl, &
       slope,factin,dfact,zdum,zbtm, &
       zref,hmin,hlim,hsurf,zpack,z0,z0sn,acc, &
       rain,condt,runoff,wetextent,percozone, &
       dryzone,bareice,den,enst,enfi,totmfi, &
       totmst,dm,absden,absdm,nk,zbal,ztal, &
       jyear,z,dz,rmass,temp,dens,watcon, &
       wetage,rppart,source)

    use smb_debug
    use smb_mecons, only: SMBInit,SMBUpdateFromOrog
    
    implicit none

    ! Parmaters
    integer :: nx       !Number of grid point in x-direction
    integer :: ny       !Number of grid point in y-direction
    integer :: mnl      !MAXIMUM number of vertical layer
    integer :: mnal     !MAXIMUM number of annual layer
    integer :: dx       !grid resolution (m)
    real(rk) :: dz0         !dz0 is the desired grid point size at the surface
    real(rk) :: zdeep       !zdeep the depth of the model 
    real(rk) :: dzdeep       !size of lowermost grid cell
    integer :: mindt    !timestep in min
    logical  :: outmap  !extensive output saving (includes saving map output, i.e. large files) 
    integer :: albkind  !choose between 3 different albedo schemes 
    logical :: VanDus   !Conductivity as in Van Dusen,1929 (true) or Sturm, 1997
    integer :: kindz0   !roughness lengths (1:z0t=z0; 2:z0t=0.00001; 3:Brutsaert; 4:Andreas)
    integer :: kstab    !stability functions(1=Businger; 2=Hogstrom; 3=Beljaars and Holtslag; 4=Dyer)
    logical :: supice   !run off according to Zuo and Oe (true) or all water on top of ice runs off (false)
    real(rk) :: t0ini       !Initial surface temperature (C)
    real(rk) :: tempb       !Temperature in the lower layer (C)
    real(rk) :: rint        !initial snow pack depth at elevation of measurements (m)
    real(rk) :: denup       !density of upper layer (kg/m3)
    real(rk) :: dendw       !density of lower layer (kg/m3)
    logical :: dry      !initial snow layer dry (true) or wet (false)
    real(rk) :: alboldwet   !albedo of old wet snow
    real(rk) :: albolddry   !albedo of old dry snow
    real(rk) :: albfrs      !albedo of fresh snow
    real(rk) :: albice      !albedo of ice
    real(rk) :: albsup      !albedo of superimposed ice
    real(rk) :: albwat      !albedo of water
    real(rk) :: decalb      !decay time for snow albedo
    real(rk) :: scswat      !scaling height (mm) for the effect of surficial water on the surface albedo 
    real(rk) :: scsnow      !scaling height for the effect of fresh snow
    real(rk) :: z0fs        !momentum roughness length of frozen snow (m)
    real(rk) :: z0ms        !momentum roughness length of melting snow (m)
    real(rk) :: z0ice       !momentum roughness length of ice (m)
    real(rk) :: absice      !radiation penetration in ice (m-1)
    real(rk) :: abslow      !radiation penetration in snow (m-1)
    real(rk) :: tsc0        !time scale for runoff of surficial water on a horizontal surface (days) 
    real(rk) :: tsc1        !time scale for runoff of surficial water on a surface with slope of 1 degree (days)
    real(rk) :: tsc10       !time scale for runoff of surficial water on a steep surface (days) 
    real(rk) :: tscfac      !ratio of time scales of runoff below and above surface
    real(rk) :: zwind       !height of wind velocity input (m) for bulk method or highest level profile method
    real(rk) :: ztemp       !height of temperature input (m) for bulk method or highest level profile method
    integer :: ist              ! loop stuff
    real(rk) :: dt            ! timestep in sec
    integer :: x1            ! coordinates target 
    integer :: x2            ! coordinates target 
    integer :: y1            ! coordinates target 
    integer :: y2            ! coordinates target 
    real(rk) :: coefdz 
    real(rk) :: facstr
    real(rk) :: bounl
    real(rk) :: bounh

    integer,pointer,dimension(:,:):: nal             !nb of annual layer
    integer,pointer,dimension(:,:):: nl              !nb of layer
    real(rk),pointer,dimension(:,:):: tatm            !atmospheric input
    real(rk),pointer,dimension(:,:):: wvel            !atmospheric input
    real(rk),pointer,dimension(:,:):: prmass          !atmospheric input
    real(rk),pointer,dimension(:,:):: glr             !atmospheric input
    real(rk),pointer,dimension(:,:):: atmrad          !atmospheric input
    real(rk),pointer,dimension(:,:):: rlhf            !energy related
    real(rk),pointer,dimension(:,:):: shf             !energy related
    real(rk),pointer,dimension(:,:):: terrad          !energy related
    real(rk),pointer,dimension(:,:):: swrtop          !energy related
    real(rk),pointer,dimension(:,:):: fluxrn          !energy related
    real(rk),pointer,dimension(:,:):: enextr          !energy related
    real(rk),pointer,dimension(:,:):: tempin          !temperature
    real(rk),pointer,dimension(:,:):: t0              !temperature
    real(rk),pointer,dimension(:,:):: surfwt          !albedo related
    real(rk),pointer,dimension(:,:):: frsnow          !albedo related
    real(rk),pointer,dimension(:,:):: tsnwfall        !albedo related
    real(rk),pointer,dimension(:,:):: alb             !albedo related
    real(rk),pointer,dimension(:,:):: conden          !Turb HF related
    real(rk),pointer,dimension(:,:):: Psurf           !Turb HF related
    real(rk),pointer,dimension(:,:):: eatm            !Turb HF related
    real(rk),pointer,dimension(:,:):: esat            !Turb HF related
    real(rk),pointer,dimension(:,:):: denair          !Turb HF related
    real(rk),pointer,dimension(:,:):: e0              !Turb HF related
    real(rk),pointer,dimension(:,:):: rhum            !Turb HF related
    real(rk),pointer,dimension(:,:):: frivel          !Turb HF related
    real(rk),pointer,dimension(:,:):: topsl           !runoff related
    real(rk),pointer,dimension(:,:):: slope           !runoff related
    real(rk),pointer,dimension(:,:):: factin          !runoff related
    real(rk),pointer,dimension(:,:):: dfact           !runoff related
    real(rk),pointer,dimension(:,:):: zdum            !grid adjustment
    real(rk),pointer,dimension(:,:):: zbtm            !grid adjustment
    real(rk),pointer,dimension(:,:):: zref            !geometry
    real(rk),pointer,dimension(:,:):: hmin            !geometry
    real(rk),pointer,dimension(:,:):: hlim            !geometry
    real(rk),pointer,dimension(:,:):: hsurf           !geometry
    real(rk),pointer,dimension(:,:):: zpack           !geometry
    real(rk),pointer,dimension(:,:):: z0              !geometry
    real(rk),pointer,dimension(:,:):: z0sn            !geometry
    real(rk),pointer,dimension(:,:):: acc             !mass balance
    real(rk),pointer,dimension(:,:):: rain            !mass balance
    real(rk),pointer,dimension(:,:):: condt           !mass balance
    real(rk),pointer,dimension(:,:):: runoff          !mass balance
    real(rk),pointer,dimension(:,:):: wetextent       !defined zones
    real(rk),pointer,dimension(:,:):: percozone       !defined zones
    real(rk),pointer,dimension(:,:):: dryzone         !defined zones
    real(rk),pointer,dimension(:,:):: bareice         !defined zones
    real(rk),pointer,dimension(:,:):: den             !mass-energy cons.
    real(rk),pointer,dimension(:,:):: enst            !mass-energy cons.
    real(rk),pointer,dimension(:,:):: enfi            !mass-energy cons.
    real(rk),pointer,dimension(:,:):: totmfi          !mass-energy cons.
    real(rk),pointer,dimension(:,:):: totmst          !mass-energy cons.
    real(rk),pointer,dimension(:,:):: dm              !mass-energy cons.
    real(rk),pointer,dimension(:,:):: absden          !mass-energy cons.
    real(rk),pointer,dimension(:,:):: absdm           !mass-energy cons.
    integer,pointer,dimension(:,:,:)::nk          !Nb of nl per annual layer
    real(rk),pointer,dimension(:,:,:)::zbal        !depth of btm of annual lyr 
    real(rk),pointer,dimension(:,:,:)::ztal        !depth of top of annual lyr 
    integer,pointer,dimension(:,:):: byear       !
    integer,pointer,dimension(:,:,:):: jyear       !year for deposition of layer nk
    real(rk),pointer,dimension(:,:,:):: z           !grid variables
    real(rk),pointer,dimension(:,:,:):: dz          !grid variables
    real(rk),pointer,dimension(:,:,:):: rmass       !englacial var.
    real(rk),pointer,dimension(:,:,:):: temp        !englacial var.
    real(rk),pointer,dimension(:,:,:):: dens        !englacial var.
    real(rk),pointer,dimension(:,:,:):: watcon      !englacial var.
    real(rk),pointer,dimension(:,:,:):: wetage      !englacial var.
    real(rk),pointer,dimension(:,:,:):: rppart      !energy penetration
    real(rk),pointer,dimension(:,:,:):: source      !energy penetration



    !The next three files contains the tunable parameters:
    open (file='smb_config/offline',unit=10,status='old')
    read (10,*) qoceansize
    read (10,*) whichtest  
    close(10)

    if (whichtest.eq.1) then
       x1=xtarg
       y1=ytarg
       x2=xtarg
       y2=ytarg
       print*,'Test option 1, short run'
    elseif (whichtest.eq.2) then
       x1=1
       y1=1
       x2=nx
       y2=ny
       print*,'Test option 2, long run!!'
    else
       print*,'choose option 1 or 2 for test...exiting'
       stop  ! NB we can pass back an error cond using alt form of return
    endif

    rinter=360/mindt                       !CHANGE THIS VALUE IF NOT 6-HOURLY INPUT 
    nin=1+rinter                           !number of interpolated input required
    daycount=(60/mindt)*24                                   
    daycount1=(60/mindt)*24+1 

    ! In force module
    allocate(quickocean(qoceansize)); quickocean=0.0
    allocate(dayroff(nx,ny)); dayroff=0
    allocate(drain(nx,ny)); drain=0
    allocate(SLM(nx,ny));   SLM=0
    allocate(orog(nx,ny));  orog=0
    allocate(Ta(nx,ny,nin)); Ta=0
    allocate(prec(nx,ny,nin)); prec=0
    allocate(wspeed(nx,ny,nin)); wspeed=0
    allocate(Tdew(nx,ny,nin)); Tdew=0
    allocate(SWdown(nx,ny,nin)); SWdown=0
    allocate(LWdown(nx,ny,nin)); LWdown=0
    allocate(press(nx,ny,nin)); press=0

    yr=1                                   !CAN BE REMOVED WITH ON-LINE RUNS

    ! Load the Sea-Land mask, orography, and mask used for runoff location
    open(10,file='smb_data/grmask20.txt')
    read(10,*) SLM
    close(10)
    SLM = SLM * 8.36  ! just scale it up so I can make the test .ge. 3.0 -> 25.0

    open(15,file='smb_data/grelev20.txt')
    read(15,*) orog
    close(15)

    open(10,file='smb_data/drainzone.txt')
    read(10,*) drain
    close(10)

    call SMBInit(nx,ny,mnl,mnal,dz0,zdeep,dzdeep, &
         mindt,outmap,'smb_config/online',chcons, &
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

    ! Available to be called by outside world, everytime
    ! orog changes
    ! GW:  Need to think of a way to ensure this?
    ! MB when using ERA 40, no need to worry about updating this,
    ! as long as the runs remains  about 10 years long
    call SMBUpdateFromOrog(nx,ny,dx,orog,slope, &
         supice,tsc10,tsc0,tsc1,mindt,tscfac,dfact,factin)

  end subroutine SMBInitForce

    ! NB TEMPORARY WRAPPER SUBROUTINE 
    subroutine SMBForceLoop(nx,ny,mnl,mnal,dz0,zdeep,dzdeep, &
         mindt,outmap,albkind, &
         VanDus,kindz0,kstab, &
         alboldwet,albolddry,albfrs,albice,albsup,albwat,decalb,scswat, &
         scsnow,z0fs,z0ms,z0ice,abslow, &
         zwind,ztemp,ist,dt, &
         iyear, &
         x1,x2,y1,y2,ratec,cst,dro,coefdz,facstr, &
         bounl,bounh,nl,tatm, &
         wvel,prmass,glr,atmrad,rlhf,shf,terrad,swrtop, &
         fluxrn,enextr,tempin,t0,surfwt,frsnow,tsnwfall,alb,conden, &
         Psurf,eatm,esat,denair,e0,rhum,frivel,topsl, &
         factin,dfact,zdum,zbtm,zref,hmin,hlim,hsurf,zpack,z0,z0sn, &
         acc,rain,condt,runoff,percozone,dryzone,bareice, &
         den,enst,enfi,totmfi,totmst,dm,absden,absdm, &
         nk,jyear,z,dz,rmass,temp,dens,watcon,wetage, &
         rppart,source, &
         totacc,totrn,totcon,totoff,massbal,rkarma)

      use smb_const
      use smb_mecons, only: SMBStep,SMBWriteReal3dToFile
      use smb_debug
      
      implicit none
      
      ! Parameters
      integer :: nx       !Number of grid point in x-direction
      integer :: ny       !Number of grid point in y-direction
      integer :: mnl      !MAXIMUM number of vertical layer
      integer :: mnal     !MAXIMUM number of annual layer
      real(rk) :: dz0         !dz0 is the desired grid point size at the surface
      real(rk) :: zdeep       !zdeep the depth of the model 
      real(rk) :: dzdeep       !size of lowermost grid cell
      integer :: mindt    !timestep in min
      logical  :: outmap  !extensive output saving (includes saving map output, i.e. large files) 
      integer :: albkind  !choose between 3 different albedo schemes 
      logical :: VanDus   !Conductivity as in Van Dusen,1929 (true) or Sturm, 1997
      integer :: kindz0   !roughness lengths (1:z0t=z0; 2:z0t=0.00001; 3:Brutsaert; 4:Andreas)
      integer :: kstab    !stability functions(1=Businger; 2=Hogstrom; 3=Beljaars and Holtslag; 4=Dyer)
      real(rk) :: alboldwet   !albedo of old wet snow
      real(rk) :: albolddry   !albedo of old dry snow
      real(rk) :: albfrs      !albedo of fresh snow
      real(rk) :: albice      !albedo of ice
      real(rk) :: albsup      !albedo of superimposed ice
      real(rk) :: albwat      !albedo of water
      real(rk) :: decalb      !decay time for snow albedo
      real(rk) :: scswat      !scaling height (mm) for the effect of surficial water on the surface albedo 
      real(rk) :: scsnow      !scaling height for the effect of fresh snow
      real(rk) :: z0fs        !momentum roughness length of frozen snow (m)
      real(rk) :: z0ms        !momentum roughness length of melting snow (m)
      real(rk) :: z0ice       !momentum roughness length of ice (m)
      real(rk) :: abslow      !radiation penetration in snow (m-1)
      real(rk) :: zwind       !height of wind velocity input (m) for bulk method or highest level profile method
      real(rk) :: ztemp       !height of temperature input (m) for bulk method or highest level profile method
      integer :: ist              ! loop stuff
      real(rk) :: dt            ! timestep in sec
      integer :: iyear         ! current year
      integer :: x1            ! coordinates target 
      integer :: x2            ! coordinates target 
      integer :: y1            ! coordinates target 
      integer :: y2            ! coordinates target 
      real(rk) :: ratec            ! use for dry snow densification
      real(rk) :: cst              ! use for dry snow densification
      real(rk) :: dro              ! use for dry snow densification
      real(rk) :: rkarma           ! von karman constant
      real(rk) :: coefdz 
      real(rk) :: facstr
      real(rk) :: bounl
      real(rk) :: bounh
      integer,dimension(:,:) :: nl              !nb of layer
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
      real(rk),dimension(:,:) :: totacc          !mass balance
      real(rk),dimension(:,:) :: totrn           !mass balance
      real(rk),dimension(:,:) :: totcon          !mass balance
      real(rk),dimension(:,:) :: totoff          !mass balance
      real(rk),dimension(:,:) :: massbal         !mass balance

      ! Local Variables
      integer :: j                ! loop stuff
      integer :: t                ! loop stuff
      integer :: m                ! loop stuff
      integer :: nt               ! loop stuff

      real(rk),dimension(nx,ny) :: esat_dew        !Turb HF related
      real(rk),dimension(nx,ny) :: esat_air        !Turb HF related
      real(rk),dimension(nx,ny) :: Tdewp           !atmospheric input


      ! Adjust some parameters, for long simulations (10 years)
      ! THIS BIT HAS TO DO WITH THE FILES NEEDED FOR OFF-LINE RUNS.

1233  if (yr.eq.1) then
         chyear='1991'
         chyear2='1992'
         iyear=1991
         goto 1234
      elseif (yr.eq.1461) then
         chyear='1992'
         chyear2='1993'
         iyear=1992
         goto 1235
      elseif (yr.eq.2925) then
         chyear='1993'
         chyear2='1994'
         iyear=1993
         goto 1234
      elseif (yr.eq.4385) then
         chyear='1994'
         chyear2='1995'
         iyear=1994
         goto 1234
      elseif (yr.eq.5845) then
         chyear='1995'
         chyear2='1996' 
         iyear=1995
         goto 1234
      elseif (yr.eq.7305) then
         chyear='1996'
         chyear2='1997' 
         iyear=1996
         goto 1235
      elseif (yr.eq.8769) then
         chyear='1997'
         chyear2='1998'
         iyear=1997
         goto 1234
      elseif (yr.eq.10229) then
         chyear='1998'
         chyear2='1999'  
         iyear=1998
         goto 1234
      elseif (yr.eq.11689) then
         chyear='1999'
         chyear2='2000'
         iyear=1999
         goto 1234
      elseif (yr.eq.13149) then
         chyear='2000'
         chyear2='2001'
         iyear=2000
         goto 1235
      elseif (yr.eq.14613) then
         chyear='2001'
         chyear2='2001'
         iyear=2001
         goto 1234
      endif

1234  ntotal=mout1*(60.0/mindt)   !total number of iterations
      nt=nt1                      !number of record in the bin.file 
      goto 1236

1235  ntotal=mout2*(60.0/mindt)
      nt=nt2
 
      ! Allocation of arrays for the output
1236  allocate (quickoceantime(qoceansize,ndyear),massbalmap(nx,ny,ndyear))
      if (outmap) then
         allocate (runoffmap(nx,ny,ndyear),albmap(nx,ny,ndyear),     &
              drymap(nx,ny,ndyear))  
      endif
      
      commpath=trim('/export/dryas/array-02/ggxmb/newint/')
      file1=commpath//'LWdown'//chyear//'.bin'
      file2=commpath//'SWdown'//chyear//'.bin'
      file3=commpath//'Prec'//chyear//'.bin'
      file4=commpath//'SLpress'//chyear//'.bin'
      file5=commpath//'T2mcorr'//chyear//'-5k.bin'
      file6=commpath//'Wspeed'//chyear//'.bin'
      file7=commpath//'Tdewcorr'//chyear//'-5k.bin'
      
      file11=commpath//'LWdown'//chyear2//'.bin' 
      file22=commpath//'SWdown'//chyear2//'.bin'
      file33=commpath//'Prec'//chyear2//'.bin'
      file44=commpath//'SLpress'//chyear2//'.bin'
      file55=commpath//'T2mcorr'//chyear2//'-5k.bin'
      file66=commpath//'Wspeed'//chyear2//'.bin'
      file77=commpath//'Tdewcorr'//chyear2//'-5k.bin'
         
      print*,'file of current year =',file1
      print*,'file of next year =',file11

      day=1
      count=1
      timecall=rinter

      if (outmap) then
         albmap=0.0
         runoffmap=0.0
         drymap=0.0
      endif
      massbalmap=0.0

      print*,'INIT done for year ',chyear,', starting time loop'
      print*,'parameters are: outmap=',outmap
      print*,'mass and energy conservation',chcons

      ! ==MAIN TIME LOOP==
      !    do t=1,ntotal
      do t=7248,10176   !run from 1st june to july 31st


   ! Call the input files
         !OUT: LWdown/SWdown/prec/SLpress/Ta/wspeed/Tdew
         if (timecall.eq.rinter.and.(1+t/timecall).le.(nt-1)) then
            timecall=1+t/timecall
            call inputdata(nx,ny,timecall,trim(file1),nin,LWdown)  
            call inputdata(nx,ny,timecall,trim(file2),nin,SWdown)
            call inputdata(nx,ny,timecall,trim(file3),nin,prec)
            call inputdata(nx,ny,timecall,trim(file4),nin,press)
            call inputdata(nx,ny,timecall,trim(file5),nin,Ta)
            call inputdata(nx,ny,timecall,trim(file6),nin,wspeed)
            call inputdata(nx,ny,timecall,trim(file7),nin,Tdew)
            m=1
            timecall=1
         elseif (timecall.eq.rinter.and.(1+t/timecall).eq.nt.and.yr.lt.14613) then
            timecall=1+t/timecall
            call inputdata2(nx,ny,timecall,trim(file1),trim(file11),nin,LWdown)  
            call inputdata2(nx,ny,timecall,trim(file2),trim(file22),nin,SWdown)
            call inputdata2(nx,ny,timecall,trim(file3),trim(file33),nin,prec)
            call inputdata2(nx,ny,timecall,trim(file4),trim(file44),nin,press)
            call inputdata2(nx,ny,timecall,trim(file5),trim(file55),nin,Ta)
            call inputdata2(nx,ny,timecall,trim(file6),trim(file66),nin,wspeed)
            call inputdata2(nx,ny,timecall,trim(file7),trim(file77),nin,Tdew)
            m=1
            timecall=1
         else
            timecall=timecall+1
            m=m+1
         endif

         tatm=Ta(:,:,m)-tkel                  !Temperature in deg C
         prmass=prec(:,:,m)*1000/(360.0/mindt)! precip in mm for 30 min
         wvel=wspeed(:,:,m)                   ! wind speed, in m s-1
         Tdewp=Tdew(:,:,m)-tkel                ! Dew point Temp, in deg C provided by ERA40
         glr=SWdown(:,:,m)/21600              ! in W m-2 (21600 sec in 6 hours)
         atmrad=LWdown(:,:,m)/21600           ! in W m-2    
         Psurf=press(:,:,m)*exp(-grav*orog/(gascda*(tav+tkel))) !Surface pressure from sea level P.
         call satwvp2D(Tdewp,2,esat_dew)          ! Saturated Water vapour pressure at tdew, 2m
         call satwvp2D(tatm,2,esat_air)           ! satur. water vapour pressure at tatm, 2m
         rhum=esat_dew/esat_air                    ! relative humidity at 2m


         call SMBStep(mnl,mnal,dz0,zdeep,dzdeep, &
              mindt,albkind, &
              VanDus,kindz0,kstab, &
              alboldwet,albolddry,albfrs,albice,albsup,albwat,decalb,scswat, &
              scsnow,z0fs,z0ms,z0ice,abslow, &
              zwind,ztemp,ist,dt, &
              iyear, &
              x1,x2,y1,y2,ratec,cst,dro,coefdz,facstr, &
              bounl,bounh,nl,SLM, &
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

         !SUM UP RUNOFF, USE THIS ONE IN THE TRACKRUNOFF subroutine
         dayroff=dayroff+runoff
      
         if (sum(dayroff).gt.0.0.and.count.eq.daycount) then
            ! TRACK RUNOFF -IF ANY- - ONCE PER DAY - need to check if 'bareice' and 'dens'
            ! varies a lot or not.
            call trackrunoff(nx,ny,qoceansize,dayroff,SLM,quickocean,drain)
            quickoceantime(:,day)=quickocean  ! GW--need to excise this...
            dayroff=0.0
         endif
         
         if (outmap) then
            runoffmap(:,:,day)=runoffmap(:,:,day)+runoff
            albmap(:,:,day)=albmap(:,:,day)+alb
            drymap(:,:,day)=drymap(:,:,day)+dryzone
         endif

         ! accumulate the mass balance at every time step
         massbalmap(:,:,day)=massbalmap(:,:,day)+massbal 

         ! Once a day...
         !(60/mindt)*24+1=97 if dt=15min, 73 if dt=20min, 49 for dt=30min
         if (count.eq.daycount1) then  
            day=day+1
            count=1
         endif
         count=count+1

         !     if (t.eq.ntotal) then
         if (t.eq.10176) then
1111        print*,'WRITING some output'
            print*,'mass balance at test point', sum(massbalmap(xtarg,ytarg,:))
            print*,'mass balance of the whole ice sheet',sum(massbalmap)*20*20/1e6

            if (outmap) then

               do j=1,ndyear
                  drymap(:,:,j) = drymap(:,:,j)/daycount
                  albmap(:,:,j) = albmap(:,:,j)/daycount
               enddo
               
               tmpFile = 'RES/soc/drymap.'//chyear
               call SMBWriteReal3dToFile(tmpFile,drymap)
               tmpFile = 'RES/soc/massbalmap.'//chyear
               call SMBWriteReal3dToFile(tmpFile,massbalmap)
               tmpFile = 'RES/soc/runoffmap.'//chyear
               call SMBWriteReal3dToFile(tmpFile,runoffmap)
               tmpFile = 'RES/soc/albmap.'//chyear
               call SMBWriteReal3dToFile(tmpFile,albmap)

            endif
         endif
      end do
      ! =================

      if (outmap) then
         if(associated(runoffmap)) deallocate(runoffmap)
         if(associated(albmap))    deallocate(albmap)
         if(associated(drymap))    deallocate(drymap)
      endif
      if(associated(quickoceantime)) deallocate(quickoceantime)
      if(associated(massbalmap))     deallocate(massbalmap)

      if (yr.eq.1) then
         return ! REMOVE THAT TO RUN FOR MORE THAN 1 YR
         yr=1461
         goto 1233
      elseif (yr.eq.1461) then
         yr=2925
         goto 1233
      elseif (yr.eq.2925) then
         yr=4385
         goto 1233
      elseif (yr.eq.4385) then
         yr=5845
         goto 1233
      elseif (yr.eq.5845) then
         yr=7305
         goto 1233
      elseif (yr.eq.7305) then
         yr=8769
         goto 1233
      elseif (yr.eq.8769) then
         yr=10229
         goto 1233
      elseif (yr.eq.10229) then
         yr=11689
         goto 1233
      elseif (yr.eq.11689) then
         yr=13149
         goto 1233
      elseif (yr.eq.13149) then
         yr=14613
         goto 1233
      endif

    end subroutine SMBForceLoop

    subroutine SMBCleanupForce

      implicit none

      if(associated(drain))   deallocate(drain)
      if(associated(SLM))     deallocate(SLM)
      if(associated(orog))    deallocate(orog)
      if(associated(Ta))      deallocate(Ta)
      if(associated(Tdew))    deallocate(Tdew)
      if(associated(wspeed))  deallocate(wspeed)
      if(associated(press))   deallocate(press)
      if(associated(SWdown))  deallocate(SWdown)
      if(associated(LWdown))  deallocate(LWdown)
      if(associated(prec))    deallocate(prec)
      if(associated(quickocean)) deallocate(quickocean)
      if(associated(dayroff)) deallocate(dayroff)

    end subroutine SMBCleanupForce

  !!
  ! This subroutines open the input data
  subroutine inputdata(nx,ny,t,datain,nin,interdat)

    implicit none

    integer, parameter :: blocksize=4
    
    ! Arguments
    integer, intent(in) :: nx,ny,t,nin
    character(len=*),intent(in) :: datain 
    real(rk), intent(out),dimension(nx,ny,nin) ::interdat 
    ! Variables
    integer :: i,ios
    real(rk), dimension(nx,ny) :: dataout,dataout1,add
    real(blocksize), dimension(nx,ny) :: readin,readin1

    open(unit=20,file=datain,status='unknown',                        &
         access='direct',form='unformatted',                          &
         iostat=ios,recl=nx*ny)

    read(unit=20,rec=t,iostat=ios) readin
    read(unit=20,rec=t+1) readin1

    dataout(:,:)=readin(:,:)
    dataout1(:,:)=readin1(:,:)

    close(unit=20)
   
    add(:,:)=(dataout1(:,:)-dataout(:,:))/(nin-1)
    do i=1,nin
       interdat(:,:,i)=dataout(:,:)+(i-1)*add(:,:)
    enddo

    return
  end subroutine inputdata

  ! This subroutines open the input data at the limit of 2 different years
  subroutine inputdata2(nx,ny,t,datain,datain2,nin,interdat)

    implicit none

    integer, parameter :: blocksize=4

    ! Arguments
    integer, intent(in) :: nx,ny,t,nin
    real(rk), intent(out),dimension(nx,ny,nin) :: interdat  
    character(len=*),intent (in) :: datain,datain2
    ! Variables
    integer :: i
    real(rk),dimension(nx,ny)::dataout,dataout1,add
    real(blocksize), dimension(nx,ny)::readin,readin1

    open(unit=20,file=datain,status='unknown',                        &
         access='direct',form='unformatted',recl=nx*ny)
    open(unit=21,file=datain2,status='unknown',                        &
         access='direct',form='unformatted',recl=nx*ny)

    read(unit=20,rec=t) readin
    read(unit=21,rec=1) readin1

    dataout(:,:)=readin(:,:)
    dataout1(:,:)=readin1(:,:)
    
    close(unit=20)
    close(unit=21)
    
    add(:,:)=(dataout1(:,:)-dataout(:,:))/(nin-1)
    do i=1,nin
       interdat(:,:,i)=dataout(:,:)+(i-1)*add(:,:)
    enddo

    return
  end subroutine inputdata2

  !!
  !   SUBROUTINE TO TRACK RUNOFF
  subroutine trackrunoff(nx,ny,qoceansize,roff,SLM,quickocean,drain)
    
    implicit none
    
    ! Arguments
    integer, intent(in) :: nx,ny,qoceansize
    real(rk), intent(in), dimension(nx,ny) :: roff
    real(rk), intent(in), dimension(nx,ny) :: SLM
    integer, intent(in), dimension(nx,ny) :: drain
    real(rk), intent(out), dimension(qoceansize) :: quickocean
    ! Variables
    integer :: x,y
    
    quickocean=0.0
    
    do x=1,nx
       do y=1,ny
          if (SLM(x,y).lt.25.0.or.(roff(x,y).eq.0.0)) then
             continue  ! effectively empty
          else
             quickocean(drain(x,y))=roff(x,y)+quickocean(drain(x,y))
          endif
       enddo   ! end main y loop
    enddo   ! end main x loop
    
    return
  end subroutine trackrunoff



    subroutine satwvp2D(t,kind,eatm) 

      implicit none

      real(rk), intent(in), dimension(:,:) :: t
      real(rk), intent(out), dimension(:,:):: eatm
      integer, intent(in) :: kind
      
      ! kind=1: saturation vapour pressure for ice
      ! kind=2: saturation vapour pressure for water
      
      if (kind.eq.1) then
         eatm=(6.1091+t*(5.0347E-01+t*(1.8860E-02                        &
              +t*(4.1762E-04+t*(5.8247E-06+t*(4.8388E-08+                  &
           1.8388E-10*t))))))*100.0
      else
         eatm=(6.1078+t*(4.4365E-01+t*(1.4289E-02                        &
              +t*(2.6506E-04+t*(3.0312E-06+t*(2.0341E-08+                  &
              6.1368E-11*t))))))*100.0
      endif
      
      return
    end subroutine satwvp2D

end module smb_force
