
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_modu.f90 - part of the GLIMMER ice model         + 
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

module glimmer_types

  !*FD Holds type definitions for the derived types used by each 
  !*FD instance of the ice model. Originally, each of these types
  !*FD was a module containing variables, which were used as containers
  !*FD for global variables. However, the need to allow for multiple
  !*FD ice model instances meant that the nested derived types were instituted
  !*FD instead. However, there is probably one too many levels in this scheme. 
  !*FD It would be better if the different types here were contained in the 
  !*FD higher-level instance type (\texttt{glimmer\_params}), rather than 
  !*FD the intermediate model type (\texttt{glimmer\_global\_type}).

  use glimmer_global
  use glimmer_ncdf

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_general

    !*FD Holds fundamental parameters of the ice model geometry.

    integer :: ewn   !*FD The number of grid-points in the E-W direction.
    integer :: nsn   !*FD The number of grid-points in the N-S direction.
    integer :: upn   !*FD The number of vertical levels in the model.

  end type glimmer_general

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_options

    !*FD Holds user options controlling the methods used in the ice-model
    !*FD integration.

    integer :: whichtemp

    !*FD Method of ice temperature calculation:
    !*FD \begin{description} 
    !*FD \item[0] Set column to surface air temperature
    !*FD \item[1] Do full temperature solution (also find vertical velocity
    !*FD and apparent vertical velocity)
    !*FD \end{description}

    integer :: whichartm

    !*FD Method of calculation of surface air temperature:
    !*FD \begin{description} 
    !*FD \item[0] Linear function of surface elevation 
    !*FD \item[1] Cubic function of distance from domain centre
    !*FD \item[2] Linear function of distance from domain centre
    !*FD \item[3] Greenland conditions (function of surface elevation 
    !*FD and latitude) including forcing
    !*FD \item[4] Antarctic conditions (sea-level air temperature 
    !*FD -- function of position)
    !*FD \item[5] Uniform temperature, zero range (temperature set in
    !*FD \texttt{cons} namelist) 
    !*FD \item[6] Uniform temperature, corrected for height, zero range.
    !*FD \item[7] Use large-scale temperature and range. 
    !*FD \end{description}

    integer :: whichthck

    !*FD Source of initial conditions: 
    !*FD \begin{description} 
    !*FD \item[1] Read from file 
    !*FD \item[2] Set equal to one time-step of net accumulation
    !*FD (where positive)
    !*FD \item[3] Stepped, linear function of distance from domain centre 
    !*FD \item[4] Read from file 
    !*FD \item[5--7] Unknown 
    !*FD \end{description}

    integer :: whichflwa

    !*FD Method for calculating flow factor $A$:
    !*FD \begin{description} 
    !*FD \item[0] \emph{Patterson and Budd} relationship 
    !*FD \item[1] \emph{Patterson and Budd} relationship, 
    !*FD with temperature set to $-10^{\circ}\mathrm{C}$ 
    !*FD \item[2] Set equal to $1\times 10^{-16}\,\mathrm{yr}^{-1}
    !*FD \,\mathrm{Pa}^{-n}$
    !*FD \end{description}

    integer :: whichisot

    !*FD Bedrock elevation: 
    !*FD \begin{description} 
    !*FD \item[0] Fixed at input values
    !*FD \item[1] Local function of ice loading history (ODE)
    !*FD \item[2] Local function of ice loading history (ODE) with flexure
    !*FD \end{description}

    integer :: whichslip

    !*FD Horizontal bed velocity: 
    !*FD \begin{description} 
    !*FD \item[0] Linear function of gravitational driving stress 
    !*FD \item[1--3] Unknown 
    !*FD \item[4] Set to zero everywhere 
    !*FD \end{description}

    integer :: whichbwat

    !*FD Basal water depth: 
    !*FD \begin{description} 
    !*FD \item[0] Calculated from local basal water balance 
    !*FD \item[1] as {\bf 0}, including constant horizontal flow 
    !*FD \item[2] Set to zero everywhere 
    !*FD \end{description}

    integer :: whichmarn

    !*FD Ice thickness: 
    !*FD \begin{description} 
    !*FD \item[0] Set thickness to zero if relaxed bedrock is more 
    !*FD than certain water depth (??) 
    !*FD \item[1] Set thickness to zero if floating 
    !*FD \item[2] No action 
    !*FD \end{description}

    integer :: whichbtrc

    !*FD Basal slip coefficient: 
    !*FD \begin{description}
    !*FD \item[0] \texttt{tanh} function of basal water depth 
    !*FD \item[1] Set equal to zero everywhere 
    !*FD \end{description}

    integer :: whichacab

    !*FD Net accumulation: 
    !*FD \begin{description}
    !*FD \item[0] EISMINT moving margin 
    !*FD \item[1] PDD mass-balance model [recommended] 
    !*FD \item[2] Accumulation only 
    !*FD \end{description}

    integer :: whichstrs

    !*FD Stress solution: 
    !*FD \begin{description}
    !*FD \item[0] Zeroth-order 
    !*FD \item[1] First-order 
    !*FD \item[2] Vertically-integrated first-order 
    !*FD \item[3] No action (use when velocity found elsewhere) 
    !*FD \end{description}

    integer :: whichevol

    !*FD Thickness evolution method:
    !*FD \begin{description}
    !*FD \item[0] Pseudo-diffusion approach 
    !*FD \item[2] Diffusion approach (also calculates velocities) 
    !*FD \end{description}

    integer :: whichwvel

    !*FD Vertical velocities: 
    !*FD \begin{description}
    !*FD \item[0] Usual vertical integration 
    !*FD \item[1] Vertical integration constrained so that 
    !*FD upper kinematic B.C. obeyed 
    !*FD \end{description}

    integer :: whichprecip

    !*FD Source of precipitation:
    !*FD \begin{description}
    !*FD \item[0] Uniform precipitation rate (set internally 
    !*FD at present) 
    !*FD \item[1] Use large-scale precipitation rate 
    !*FD \item[2] Use parameterization of \emph{Roe and Lindzen} 
    !*FD \end{description}

  end type glimmer_options

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_geometry

    !*FD Holds fields and other information relating to the
    !*FD geometry of the ice sheet and bedrock.

    real(dp),dimension(:,:),pointer :: thck 
    
    !*FD The thickness of the ice, divided by \texttt{thck0}.

    real(dp),dimension(:,:),pointer :: usrf
    
    !*FD The elevation of the upper ice surface, divided by \texttt{thck0}.

    real(dp),dimension(:,:),pointer :: lsrf 
    
    !*FD The elevation of the lower ice surface, divided by \texttt{thck0}.

    real(dp),dimension(:,:),pointer :: topg 
    
    !*FD The elevation of the topography, divided by \texttt{thck0}.

    real(dp),dimension(:,:),pointer :: relx 
    
    !*FD The elevation of the relaxed topography, by \texttt{thck0}.

    integer, dimension(:,:),pointer :: mask 
    
    !*FD Set to zero for all points where $\mathtt{thck}=0$, otherwise non-zero.
    !*FD the non-zero points are numbered in sequence from the bottom left to the 
    !*FD top right, going along the rows.

    integer :: totpts 
    
    !*FD The total number of points with non-zero thickness

    integer, dimension(4) :: dom            !*FD I have no idea what this is for.
    logical :: empty                        !*FD I have no idea what this is for.

  end type glimmer_geometry

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_geomderv

    !*FD Holds the horizontal and temporal derivatives of the thickness and
    !*FD upper surface elevation, as well as the thickness on the staggered grid.

    real(dp),dimension(:,:),pointer :: dthckdew !*FD E-W derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdew !*FD E-W derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: dthckdns !*FD N-S derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdns !*FD N-S derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: dthckdtm !*FD Temporal derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdtm !*FD Temporal derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: stagthck !*FD Thickness averaged onto the staggered grid.

  end type glimmer_geomderv

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_velocity
    real(dp),dimension(:,:,:),pointer :: uvel, vvel, wvel, wgrd
    real(dp),dimension(:,:)  ,pointer :: uflx, vflx, diffu
    real(dp),dimension(:,:)  ,pointer :: ubas, vbas, btrc
  end type glimmer_velocity

  type glimmer_climate
    real(sp),dimension(:,:),pointer :: acab, artm
    real(sp),dimension(:,:),pointer :: arng, lati, ablt
    real(sp),dimension(:,:),pointer :: g_arng,g_artm  ! Global fields as input
    real(sp),dimension(:,:),pointer :: prcp
    real(sp),dimension(:,:),pointer :: presprcp  !*FD Present-day annual precip (mm)
    real(sp),dimension(:,:),pointer :: presartm  !*FD Present-day mean air-temp ($^{\circ}$C)
    real(dp),dimension(:,:),pointer :: presusrf  !*FD Present-day upper surface (km)
    real(sp) :: uprecip_rate   !*FD Uniform precipitaion rate in m/a
    real(sp) :: usurftemp      !*FD Uniform surface temperature in deg C
    real(sp) :: ice_albedo     !*FD Ice albedo. (fraction)
    real(sp) :: ulapse_rate    !*FD Uniform lapse rate in deg C/km
  end type glimmer_climate

  type glimmer_temper
    real(dp),dimension(:,:,:),pointer :: temp
    real(dp),dimension(:,:,:),pointer :: flwa  
    real(dp),dimension(:,:),  pointer :: bwat, bmlt
    integer :: niter
    real(sp) :: perturb 
    real(sp) :: grid
    integer :: tpt         ! pointer to time series data
    logical :: first1
  end type glimmer_temper

  type glimmer_forcdata
    real(sp),pointer,dimension(:,:) :: forcing
    integer :: flines
    real(sp) :: trun
  end type glimmer_forcdata

  type glimmer_stress
    real(dp),dimension(:,:,:),pointer :: efvs
    real(dp),dimension(:,:,:),pointer :: tauxz, tauyz, tauxy, tauxx, tauyy, tau, gdsx, gdsy
  end type glimmer_stress

  type glimmer_funits
    integer, dimension(8) :: indices0dx
    integer, dimension(8) :: indices0dy
    integer, dimension(n2d) :: which2d
    integer, dimension(n3d) :: which3d 
    character(fname_length) :: usrffile, topgfile, relxfile    ! input filenames
    character(fname_length) :: sigfile                         !*FD sigma coordinates file
    character(fname_length) :: prcpfile                        !*FD Precipitation file
    character(fname_length) :: presusrffile                    !*FD Present-day surface file
    character(fname_length) :: forcfile                        !*FD Temperature forcing file
    character(fname_length) :: output_stem                     !*FD output filename stem
    character(fname_length) :: ncfile                          !*FD configuration file for netCDF I/O
    type(glimmer_nc_output),pointer :: out_first=>NULL()       !*FD first element of linked list defining netCDF outputs
    type(glimmer_nc_input), pointer :: in_first=>NULL()        !*FD first element of linked list defining netCDF inputs
    character(fname_length) :: latifile                        !*FD File containing latitudes of all points
  end type glimmer_funits

  type glimmer_numerics
    real(sp) :: time         !*FD main loop counter in years
    real(sp) :: tinc         !*FD time step of main loop in years 
    real(sp) :: ntem         !*FD temperature time step in years
    real(sp) :: nvel         !*FD velocity time step in years
    real(sp) :: niso         !*FD flexure time step in years
    real(sp),dimension(3) :: nout !*FD output time step in years 
                                  !*FD for time series, horiz. and full data                        
    real(sp) :: nstr         !*FD output start time in years
    real(dp) :: alpha        !*FD richard suggests 1.5 - was a parameter in original
    real(dp) :: alphas       !*FD was a parameter in the original
    real(dp) :: thklim    
    real(dp) :: mlimit  
    real(dp) :: dew 
    real(dp) :: dns 
    real(dp) :: dt       
    real(dp) :: dttem 
    real(dp),dimension(:),pointer :: sigma 
    real(sp) :: nshlf 
  end type glimmer_numerics

  type glimmer_pddcalc

    !*FD Holds parameters for positive-degree-day mass-balance
    !*FD calculation. The table has two axes - the $x$ axis is the
    !*FD difference between mean annual and July temps, while the
    !*FD $y$- axis is the mean annual temp

    integer :: dx !*FD Spacing of values in x-direction ($^{\circ}$C) --- default=1
    integer :: dy !*FD Spacing of values in y-direction ($^{\circ}$C) --- default=1
    integer :: ix !*FD Lower bound of $x$-axis ($^{\circ}$C) --- default=0
    integer :: iy !*FD Lower bound of $y$-axis ($^{\circ}$C) --- default=-50
    integer :: nx !*FD Number of values in x-direction
    integer :: ny !*FD Number of values in y-direction

    real(sp),dimension(:,:),pointer :: pddtab !*FD PDD table - must be 
                                              !*FD allocated with dimensions nx,ny.

    logical  :: pt_alloc !*FD set \texttt{.true.} if \texttt{pddtab} has been allocated.
    real(sp) :: dailytemp, tma, tmj, dtmj
    real(sp) :: dd_sigma
    logical  :: first
  end type glimmer_pddcalc

  type glimmer_isotwk
    real(dp),pointer,dimension(:,:) :: load
    real(dp),dimension(4) :: fact
    logical :: first1
    integer :: nflx
    real(sp),dimension(:,:),pointer :: dflct
    logical :: first2
  end type glimmer_isotwk

  type glimmer_velowk
    real(dp),dimension(:),pointer :: depth, dupsw, depthw
    real(dp),dimension(:),pointer :: suvel, svvel
    real(dp),dimension(:,:),pointer :: fslip
    logical :: first1
    real(dp),dimension(:,:),pointer :: dintflwa
    logical :: first2
    real(dp),dimension(:),pointer :: dups
    logical :: first3
    logical :: first4
    real(dp),dimension(4) :: fact
    logical :: first5
    real(dp) :: watwd, watct, trc0
    real(dp) :: trcmin, marine, trcmax
    real(dp),dimension(4) :: c 
    logical :: first6
  end type glimmer_velowk

  type glimmer_pcgdwk
    real(dp), pointer, dimension(:) :: pcgval, rhsd, answ
    integer, pointer, dimension(:) :: pcgcol, pcgrow
    integer, dimension(2) :: pcgsize
    integer :: ct
    real(dp), dimension(4) :: fc
    real(dp), dimension(6) :: fc2
    integer :: mlinit
    integer :: tlinit
    logical :: first1
  end type glimmer_pcgdwk

  type glimmer_thckwk
    integer :: nwhich
    real(dp), pointer, dimension(:,:,:) :: olds
    real(sp) :: oldtime
    real(dp), dimension(:,:), pointer :: oldthck, basestate 
    real(dp) :: few, fns
    logical :: first1
  end type glimmer_thckwk

  type glimmer_tempwk
    real(dp), pointer, dimension(:,:,:) :: inittemp, dissip, initadvt
    real(dp) :: advconst(2)
    real(dp), dimension(:), pointer :: dupa, dupb, dupc, c1
    real(dp), dimension(:,:),pointer :: dups
    real(dp) :: noflow
    logical :: first1
    real(dp),dimension(4) :: cons 
    real(dp) :: zbed, dupn, wmax 
    logical :: first2
    logical :: first3
    real(dp), dimension(4) :: f
    logical :: first4
    real(dp) :: dt_wat
    integer  :: nwat
    real(dp), dimension(8) :: c 
    logical  :: first5
  end type glimmer_tempwk

  type glimmer_paramets
    real(dp) :: geot     ! W m^{-2}
    real(dp) :: fiddle   ! -
    real(sp),dimension(2) :: airt       ! K, K km^{-3}
    real(sp),dimension(3) :: nmsb  ! m yr^{-1}, yr^{-1}, m                
    real(dp) :: hydtim  ! yr^{-1} converted to s^{-1} and scaled, 
                                    ! 0 if no drainage = 0.0d0 * tim0 / scyr
    real(dp) :: isotim
    real(dp),dimension(5) :: bpar 
  end type glimmer_paramets

  type glimmer_global_type
    type(glimmer_general)  :: general
    type(glimmer_options)  :: options
    type(glimmer_geometry) :: geometry
    type(glimmer_geomderv) :: geomderv
    type(glimmer_velocity) :: velocity
    type(glimmer_climate)  :: climate
    type(glimmer_temper)   :: temper
    type(glimmer_forcdata) :: forcdata
    type(glimmer_stress)   :: stress
    type(glimmer_funits)   :: funits
    type(glimmer_numerics) :: numerics
    type(glimmer_pddcalc)  :: pddcalc
    type(glimmer_isotwk)   :: isotwk
    type(glimmer_velowk)   :: velowk
    type(glimmer_pcgdwk)   :: pcgdwk
    type(glimmer_thckwk)   :: thckwk
    type(glimmer_tempwk)   :: tempwk
    type(glimmer_paramets) :: paramets
  end type glimmer_global_type

contains

  subroutine glimmer_global_type_initialise(params)

    !*FD Initialise the model parameters. Note that this is
    !*FD \emph{not} the place to put the default values of
    !*FD parameters read in from the namelist. Those defaults
    !*FD are set in the \texttt{initial} subroutine.

    type(glimmer_global_type),intent(inout) :: params

    params%general%ewn  = 0
    params%general%nsn  = 0
    params%general%upn  = 0

    params%options%whichtemp     = 1
    params%options%whichartm     = 3
    params%options%whichthck     = 4
    params%options%whichflwa     = 0
    params%options%whichisot     = 1
    params%options%whichslip     = 4
    params%options%whichbwat     = 2
    params%options%whichmarn     = 0
    params%options%whichbtrc     = 1
    params%options%whichacab     = 2
    params%options%whichstrs     = 2
    params%options%whichevol     = 0
    params%options%whichwvel     = 0
    params%options%whichprecip   = 0

    params%geometry%dom          = 0
    params%geometry%totpts       = 0
    params%geometry%empty        = .true.

    params%climate%uprecip_rate  = 0.0
    params%climate%usurftemp     = 0.0
    params%climate%ice_albedo    = 0.0
    params%climate%ulapse_rate   = 0.0

    params%temper%niter          = 0
    params%temper%perturb        = 0.0
    params%temper%grid           = 0.0 
    params%temper%tpt            = 0
    params%temper%first1         = .true.

    params%forcdata%flines       = 0

    params%funits%indices0dx     = 0
    params%funits%indices0dy     = 0
    params%funits%which2d        = 1 
    params%funits%which3d        = 1
    params%funits%usrffile       = 'none'
    params%funits%topgfile       = 'none'
    params%funits%relxfile       = 'none'
    params%funits%sigfile        = ''
    params%funits%output_stem    = 'untitled'
    params%funits%ncfile         = ''

    params%numerics%time         = 0.0
!    params%numerics%tinc         = 20.0  ! We don't want to touch this...
    params%numerics%ntem         = 5.0  
    params%numerics%nvel         = 10.0 
    params%numerics%niso         = 5.0 
    params%numerics%nout         = (/ 1.0, 1.0e1, 1.0e1 /) 
    params%numerics%nstr         = 0.0
    params%numerics%alpha        = 0.5d0
    params%numerics%alphas       = 0.5d0
    params%numerics%thklim       = 100.0d0
    params%numerics%mlimit       = -200.0d0
    params%numerics%dew          = 20.0d3
    params%numerics%dns          = 20.0d3
    params%numerics%dt           = 0.0
    params%numerics%dttem        = 0.0
    params%numerics%nshlf        = 0.0

    params%pddcalc%dx = 1
    params%pddcalc%dy = 1
    params%pddcalc%ix = 0
    params%pddcalc%iy = -50
    params%pddcalc%nx = 31
    params%pddcalc%ny = 71
    params%pddcalc%pt_alloc=.false.
    params%pddcalc%dailytemp = 0.0
    params%pddcalc%tma =0.0 
    params%pddcalc%tmj =0.0
    params%pddcalc%dtmj =0.0
    params%pddcalc%dd_sigma = 5.0
    params%pddcalc%first=.true.

    params%isotwk%first1 = .true.
    params%isotwk%fact   = 0.0
    params%isotwk%nflx   = 0
    params%isotwk%first2 = .true.

    params%velowk%first1 = .true.
    params%velowk%first2 = .true.
    params%velowk%first3 = .true.
    params%velowk%first4 = .true.
    params%velowk%first5 = .true.
    params%velowk%watwd = 3.0d0
    params%velowk%watct = 10.0d0
    params%velowk%trc0  = 0.0
    params%velowk%trcmin = 0.0d0
    params%velowk%marine = 1.0d0
    params%velowk%trcmax = 10.0d0
    params%velowk%c      = 0.0
    params%velowk%first6 = .true.

    params%pcgdwk%ct  = 0
    params%pcgdwk%fc  = 0.0
    params%pcgdwk%fc2 = 0.0
    params%pcgdwk%mlinit = 0
    params%pcgdwk%tlinit = 0
    params%pcgdwk%first1 = .true.

    params%thckwk%nwhich = 2
    params%thckwk%oldtime = 0.0
    params%thckwk%few = 0.0
    params%thckwk%fns = 0.0
    params%thckwk%first1=.true.

    params%tempwk%advconst = 0.0
    params%tempwk%cons = 0.0
    params%tempwk%zbed = 0.0
    params%tempwk%dupn = 0.0
    params%tempwk%wmax = 0.0
    params%tempwk%f = 0.0
    params%tempwk%dt_wat = 0.0
    params%tempwk%nwat = 0
    params%tempwk%c = 0.0
    params%tempwk%noflow = -1
    params%tempwk%first1 = .true.
    params%tempwk%first2 = .true.
    params%tempwk%first3 = .true.
    params%tempwk%first4 = .true.
    params%tempwk%first5 = .true.

    params%paramets%geot = -5.0d-2
    params%paramets%fiddle = 3.0d0
    params%paramets%airt = (/ -3.150, -1.0e-2 /)
    params%paramets%nmsb = (/ 0.5, 1.05e-5, 450.0e3 /)
    params%paramets%hydtim = 1000.0d0
    params%paramets%isotim = 3000.0d0
    params%paramets%bpar = (/ 2.0d0, 10.0d0, 10.0d0, 0.0d0, 1.0d0 /)

  end subroutine glimmer_global_type_initialise

end module glimmer_types

