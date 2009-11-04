module plume_namelist

contains

subroutine namelist(nl_filename,varoutrat,seedtype,dt1,dt2,dtswtim,phi,labtim, &
                         tottim,outtim,fouttim,louttim,snottim,lnottim)

! set program parameters
! where appropriate, read from namelist file
! derive some parameters at end of routine, after read from file

  use plume_global

  implicit none

! local variables
  character(len = *),intent(in) :: nl_filename
  integer varoutrat,seedtype,l

  real(kind=kdp) dt1,dt2,dtswtim,phi,labtim
  real(kind=kdp) tottim,outtim,fouttim,louttim,snottim,lnottim
  
  real(kind=kdp) infloain,infloein,cseedfix,cinffix
      

  namelist /plume_nml/ &
       mixlayer &
       ,  restart &
       ,  frazil &
       ,  nonlin &
       ,  horturb &
       ,  entrain &
       ,  entype &
       ,  basmelt &
       ,  rholinear &
       ,  thermobar &
       ,  intrace &
       ,  vardrag &
       ,  topedit &
       ,  tangle &
       ,  negfrz &
       ,  tottim &
       ,  outtim &
       ,  labtim &
       ,  snottim &
       ,  lnottim &
       ,  dt1 &
       ,  m & 
       ,  n &
       ,  hx &
       ,  hy &
       ,  gldep &
       ,  ifdep &
       ,  wcdep &
       ,  context &
       ,  bathtype &
       ,  bsmoothit &
       ,  depinit &
       ,  depinffix &
       ,  meltinf &
       ,  salttop &
       ,  infloain &
       ,  infloein &
       ,  saltbot &
       ,  temptop &
       ,  tempbot &
       ,  phi &
       ,  ah &
       ,  kh &
       ,  cdb &
       ,  cl & 
       ,  ef &
       ,  ti &
       ,  nus &
       ,  nbar &
       ,  nice &
       ,  seedtype &
       ,  cseedfix &
       ,  cinffix

! ++++++++++++++++++
! set default values
! ++++++++++++++++++

! set ascii output (screen and file) properties 

      varoutrat = 1       ! output resolution varies:
                          ! 1 - in both directions
                          ! 2 - according to extent of plume in y-direction
                          ! 3 - according to extent of plume in x-direction

! set switches 
! ------------
      mixlayer    = .false. ! model a mixed-layer, rather than a plume
                            ! (i.e. whole domain is given initial thickness)
      restart     = .false. ! restart from previous model dump
      frazil      = .false. ! include frazil ice
      nonlin      = .true.  ! advection terms included
      horturb     = .false. ! horizontal diffusion included
      entrain     = .true.  ! entrainment included
      entype      = 1       ! entrainment parameterisation:
                            ! 1 - full kochergin formulation
                            ! 2 - fractional kochergin formulation (using ef)
                            ! 3 - halved pedersen formulation
                            ! 4 - halved and modified pedersen formulation
      basmelt     = .true.  ! include direct basal melting and freezing
      rholinear   = .true.  ! use linear equation of state instead of unesco
      thermobar   = .false. ! thermobaricity included (code needs checking)
      intrace     = .false. ! include tracers tracking the fate of each inflow
      vardrag     = .false. ! include spatially-varying drag coefficient
      topedit     = .false. ! apply hand-editing to read-in topography
      tangle      = .false. ! use turning angle when applying drag
                            ! (probably wise to check or reprogram)
      negfrz      = .false. ! artificially keep frazil concentrations positive

! set simulation time parameters
! ------------------------------

      tottim  = 0050.0d0  ! total simulation time in days
      outtim  = 050.0d0   ! file output frequency in days
      labtim  = 001.0d0   ! units in which to name files in days
      snottim = 005.0d0   ! short note output frequency in days
      lnottim = 050.0d0   ! long note output frequency in days
 
      dt1     = 0360.0d0  ! long first timestep in seconds

! test values
!c      dt1     = 150.0d0            
!c      tottim  = dt1           
!c      outtim  = dt1          
!c      labtim  = dt1         
!c      snottim = dt1        
!c      lnottim = dt1 

! set grid and domain properties
! ------------------------------
      m = 0600           ! number of east-west cells
      n = 0600           ! number of north-south cells
      hx = 5000.d0       ! east-west cell dimension (m) for uniform grid
      hy = hx            ! north-south cell dimension (m) for uniform grid
      gldep = 1400.d0    ! thickness of ice shelf at grounding line
      ifdep = 285.d0     ! thickness of ice shelf at ice front (or plateau)
      wcdep = 1600.d0    ! depth of water column beneath grounding line 

! set topography properties
! -------------------------
      context = "isomip" ! topography choices to use if bathtype = 0
                         ! choices are currently fris,larsen, or isomip
                         ! note isomip topography incomplete

      bathtype = 01     ! ice shelf bathymetry scenario:
                        ! 0 => read bathymetry from file
                        ! 1 => no walls, rises northwards
                        ! 2 => right-angle west wall, rises northwards
                        ! 3 => 45 degree west wall, rises northwards
                        ! 4 => right-angle west wall, rises || to wall
                        ! 5 => 45 degree west wall, rises || to wall
                        ! 6 => carlson inlet
                        ! 7 => rutford ice stream 
                        ! 8 => carlson inlet test 1
                        ! 9 => carlson inlet test 2
                        ! 10=> evans ice stream
                        ! 11=> whole filchner-ronne simplified
                        ! 12=> analytical ice shelf

      kcorn = int(135.d0*1000.d0/hx) ! dist. of corner from inflow (km)
      rad = int(35.d0*1000.d0/hx)    ! radius of rounding on corner (km)

      cweight = 4.d0         ! smoothing weight of central point 
      nweight = 1.d0         ! smoothing weight of neighbour point
      bsmoothit = 00         ! its of smoothing to apply throughout domain

! smooth region flags and iterations (including bsmoothit) if variable smoothing
! 1) fris
      smflag(01)   = .false. ! rutford ice stream
      smoothit(01) = 20
      smflag(02)   = .false. ! support force glacier
      smoothit(02) = 10
      smflag(03)   = .false. ! spare slot
      smoothit(03) = 00
! 2) larsen
      smflag(01)   = .false. ! choyce point
      smoothit(01) = 15
      smflag(02)   = .false. ! thuronyi bluff
      smoothit(02) = 15
      smflag(03)   = .false. ! spare slot
      smoothit(03) = 00

! set inflow properties and initial thickness
! -------------------------------------------
      ninfmin = 01             ! lowest-numbered inflow turned on
      ninfmax = 01             ! highest-numbered inflow turned on 

      depinit = 01.0d0    ! initial plume thickness if mixed-layer model
      depinffix = 01.d0   ! plume thickness assigned to inflow regions
      meltinf = 9.8d-1    ! how far along gade line inflow properties are:
                          ! 1 => properties at shelf, 0 => properties of ambient

! inflow extent if using single inflow on simple bathymetry
      infloain = 290.d0   ! western boundary (km)
      infloein = 310.d0   ! eastern boundary (km)
      knfloa = 1          ! southern boundary (cells)
      knfloe = 3          ! northern boundary (cells)

! set ambient fluid properties
! ----------------------------
      namb = 302           ! increments in ambient water column (minimum +2)
      dzincr = 10.d0       ! depth of each increment 
                           ! (so total (namb - 2)*dzincr = gldep + wcdep)
      salttop = 34.500d0   ! sea surface salinity 
      saltbot = 34.950d0   ! sea bed salinity (gldep + wcdep)
      temptop = -1.900d0   ! sea surface temperature 
      tempbot = -2.500d0   ! sea bed temperature (gldep + wcdep)

! set physical and geographical parameters
! ----------------------------------------
      g = 9.81d0       ! gravitational constant
      phi = +00.d0     ! latitude (in degrees, negative => s. hemisphere)
                       ! (set to zero to remove rotation)
      ah = 0000.d0     ! horizontal eddy viscosity
      kh = ah          ! horizontal eddy diffusivity
      cdb = 2.5d-3     ! bottom drag coefficient (background value if variable)
      cdbvar = 5.0d-3  ! bottom drag coefficient in varied areas (if variable)
      cl = 1.775d-2    ! kochergin entrainment coefficient (should be 2.75d-2)
      ef = 5.0d-1      ! kochergin entrainment factor (1=> match jung & back)
      rho0 = 1030.d0   ! reference seawater density
      rhoi = 920.d0    ! density of ice

! drag region flags used if using variable drag
      drflag(01) = .false.  !east of henry ice rise
      drflag(02) = .false.  !east of berkner island
      drflag(03) = .false.  !spare slot

! set ice-related physical parameters
! -----------------------------------
      lat = 335000.d0  ! latent heat of ice fusion
      c0 = 3974.d0     ! specific heat capacity of weddell seawater
      ci = 2009.d0     ! specific heat capacity of ice shelf
      nu0 = 1.95d-6    ! kinematic viscosity of seawater
      pr = 13.8d0      ! molecular prandtl number of seawater
      sc = 2432.d0     ! molecular schmidt number of seawater
      fta = -5.73d-2   ! slope of liquidus for seawater (tf decrease with s)
      ftb = 8.32d-2    ! offset of liquidus for seawater (tf at s=0, z=surface)
      ftc = -7.61d-4   ! freezing temp change with depth (tf decrease with zeta)
      ti = -25.d0      ! temperature of shelf (heat conduction during melting)
      si = 0.d0        ! salinity of shelf (salt trapped in ice during freezing)
      nus = +1.d0      ! nusselt number: 
                       ! >=0 - that constant value
                       ! -1  - correct full variable (not working)
                       ! -2  - correct variable (no turbulent part for large)
                       ! -3  - incorrect full hammar&shen variable (not working)
                       ! -4  - incorrect h&s variable (no turbulent for large)
      kt = 1.4d-7      ! molecular thermal diffusivity of seawater
      ks = 8.0d-10     ! molecular haline diffusivity of seawater
      ar = 2.0d-2      ! aspect ratio of frazil discs
      eps = 7.4d-6     ! turbulent dissipation rate
      nbar = 1.0d3     ! cap on number of crystals for secondary nucleation

! set frazil model parameters
! ---------------------------
      nice = 10        ! number of frazil size classes

      r(01)  = 0.01d-3 ! radius of frazil discs
      r(02)  = 0.05d-3    
      r(03)  = 0.15d-3    
      r(04)  = 0.3d-3    
      r(05)  = 0.4d-3    
      r(06)  = 0.5d-3    
      r(07)  = 0.6d-3    
      r(08)  = 0.8d-3    
      r(09)  = 1.0d-3    
      r(10)  = 2.0d-3    

      seedtype = 2      ! seeding strategy:
                        ! 1 => seed at southernmost supercooled cells 
                        !      (smedsrud & jenkins)
                        ! 2 => seed newly-supercooled cells if c(i) < cseed(i) 
                        !      (holland & feltham)

      cseedfix = 1.0d-7 ! seed concentration in class i
      cinffix  = 0.0d0  ! inflow concentration in class i

! set modelling parameters
! ------------------------
      edepth = 1.0d-3  ! critical min depth for entrainment
      mdepth = edepth  ! critical min depth for basal melting
      fdepth = edepth  ! critical min depth for frazil dynamics
      dcr    = edepth  ! critical min plume depth
      septol = 1.0d-2  ! smallest plume-ambient density anomaly tolerated before
                       ! plume separation is flagged
      small = 1.0d-15  ! smallest number

! +++++++++++++++++++++
! read values from file
! +++++++++++++++++++++

      open(21,file=nl_filename,status='old')

        read(21,plume_nml)
        
      close(21)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
! some parameters dependent upon namelist-read variables
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++

! convert timesteps to seconds

      tottim  = tottim*24.d0*3600.d0  ! total simulation time in seconds
      outtim  = outtim*24.d0*3600.d0  ! file output frequency in seconds
      labtim  = labtim*24.d0*3600.d0  ! units in which to name files
      snottim = snottim*24.d0*3600.d0 ! short note output frequency in seconds
      lnottim = lnottim*24.d0*3600.d0 ! long note output frequency in seconds
      fouttim = outtim                ! first file output in seconds
      louttim = tottim                ! last file output in seconds

! uncomment to do short test runs of a few timesteps
!      tottim  = dt1*60*24*10
!      tottim  = dt1
!      outtim  = tottim
!      labtim  = tottim
!      snottim = tottim
!      lnottim = tottim
!      fouttim = tottim
!      louttim = tottim

      dt2     = dt1                   ! short second timestep in seconds
      dtswtim = tottim                ! timestep switch time in seconds

! derive ideal inflow extents

      infloa = int(infloain*1000.d0/hx)! western boundary (cells)
      infloe = int(infloein*1000.d0/hx)! eastern boundary (cells)

! set case-specific inflow flags if using realistic bathymetry

      if (context.eq."fris") then

        inflag(01) = .true.   ! evans ice stream
        inflag(02) = .false.  ! carlson inlet
        inflag(03) = .false.  ! rutford ice stream
        inflag(04) = .false.  ! institute ice stream
        inflag(05) = .false.  ! mollereisstrom
        inflag(06) = .false.  ! foundation ice stream
        inflag(07) = .false.  ! support force glacier
        inflag(08) = .false.  ! recovery glacier
        inflag(09) = .false.  ! slessor glacier

      end if

      if (context.eq."larsen") then

        inflag(01) = .false.   ! attlee glacier
        inflag(02) = .true.    ! whole larsen b
        inflag(03) = .true.    ! whole larsen c

      end if

      if (context.eq."isomip") then

        inflag(01) = .true.    ! whole southern boundary

      end if

! determine ambient fluid gradients

      sgrad = (saltbot - salttop)/(gldep + wcdep) ! rate of s change with depth
      tgrad = (tempbot - temptop)/(gldep + wcdep) ! rate of t change with depth

! set frazil seeds and inflows

      do l = 1,nice     ! frazil seed population 
        cseed(l) = cseedfix
      end do
      cinftot = 0.d0    ! inflow frazil
      do l = 1,nice
        cinf(l) = cinffix
        cinftot = cinftot + cinf(l)
      end do

      return
end subroutine


end module


