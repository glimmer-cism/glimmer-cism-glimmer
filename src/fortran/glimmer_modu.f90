
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
  !*FD 
  !*FD Note that this \emph{is} now where the defaults are defined for these
  !*FD variables.
 
  use glimmer_global
  use glimmer_ncdf

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_general

    !*FD Holds fundamental parameters of the ice model geometry.

    integer :: ewn = 0  !*FD The number of grid-points in the E-W direction.
    integer :: nsn = 0  !*FD The number of grid-points in the N-S direction.
    integer :: upn = 1  !*FD The number of vertical levels in the model.

  end type glimmer_general

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_options

    !*FD Holds user options controlling the methods used in the ice-model
    !*FD integration.

    integer :: whichtemp = 1

    !*FD Method of ice temperature calculation:
    !*FD \begin{description} 
    !*FD \item[0] Set column to surface air temperature
    !*FD \item[1] Do full temperature solution (also find vertical velocity
    !*FD and apparent vertical velocity)
    !*FD \end{description}

    integer :: whichartm = 3

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

    integer :: whichthck = 4

    !*FD Source of initial conditions: 
    !*FD \begin{description} 
    !*FD \item[1] Read from file 
    !*FD \item[2] Set equal to one time-step of net accumulation
    !*FD (where positive)
    !*FD \item[3] Stepped, linear function of distance from domain centre 
    !*FD \item[4] Read from file 
    !*FD \item[5--7] Unknown 
    !*FD \end{description}

    integer :: whichflwa = 0

    !*FD Method for calculating flow factor $A$:
    !*FD \begin{description} 
    !*FD \item[0] \emph{Patterson and Budd} relationship 
    !*FD \item[1] \emph{Patterson and Budd} relationship, 
    !*FD with temperature set to $-10^{\circ}\mathrm{C}$ 
    !*FD \item[2] Set equal to $1\times 10^{-16}\,\mathrm{yr}^{-1}
    !*FD \,\mathrm{Pa}^{-n}$
    !*FD \end{description}

    integer :: whichisot = 1

    !*FD Bedrock elevation: 
    !*FD \begin{description} 
    !*FD \item[0] Fixed at input values
    !*FD \item[1] Local function of ice loading history (ODE)
    !*FD \item[2] Local function of ice loading history (ODE) with flexure
    !*FD \end{description}

    integer :: whichslip = 4

    !*FD Horizontal bed velocity: 
    !*FD \begin{description} 
    !*FD \item[0] Linear function of gravitational driving stress 
    !*FD \item[1--3] Unknown 
    !*FD \item[4] Set to zero everywhere 
    !*FD \end{description}

    integer :: whichbwat = 2

    !*FD Basal water depth: 
    !*FD \begin{description} 
    !*FD \item[0] Calculated from local basal water balance 
    !*FD \item[1] as {\bf 0}, including constant horizontal flow 
    !*FD \item[2] Set to zero everywhere 
    !*FD \end{description}

    integer :: whichmarn = 0

    !*FD Ice thickness: 
    !*FD \begin{description} 
    !*FD \item[0] Set thickness to zero if relaxed bedrock is more 
    !*FD than certain water depth (??) 
    !*FD \item[1] Set thickness to zero if floating 
    !*FD \item[2] No action 
    !*FD \end{description}

    integer :: whichbtrc = 1

    !*FD Basal slip coefficient: 
    !*FD \begin{description}
    !*FD \item[0] \texttt{tanh} function of basal water depth 
    !*FD \item[1] Set equal to zero everywhere 
    !*FD \end{description}

    integer :: whichacab = 2

    !*FD Net accumulation: 
    !*FD \begin{description}
    !*FD \item[0] EISMINT moving margin 
    !*FD \item[1] PDD mass-balance model [recommended] 
    !*FD \item[2] Accumulation only 
    !*FD \end{description}

    integer :: whichstrs = 2

    !*FD Stress solution: 
    !*FD \begin{description}
    !*FD \item[0] Zeroth-order 
    !*FD \item[1] First-order 
    !*FD \item[2] Vertically-integrated first-order 
    !*FD \item[3] No action (use when velocity found elsewhere) 
    !*FD \end{description}

    integer :: whichevol = 0

    !*FD Thickness evolution method:
    !*FD \begin{description}
    !*FD \item[0] Pseudo-diffusion approach 
    !*FD \item[2] Diffusion approach (also calculates velocities) 
    !*FD \end{description}

    integer :: whichwvel = 0

    !*FD Vertical velocities: 
    !*FD \begin{description}
    !*FD \item[0] Usual vertical integration 
    !*FD \item[1] Vertical integration constrained so that 
    !*FD upper kinematic B.C. obeyed 
    !*FD \end{description}

    integer :: whichprecip = 0

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

    real(dp),dimension(:,:),pointer :: thck => null()
    
    !*FD The thickness of the ice, divided by \texttt{thck0}.

    real(dp),dimension(:,:),pointer :: usrf => null()
    
    !*FD The elevation of the upper ice surface, divided by \texttt{thck0}.

    real(dp),dimension(:,:),pointer :: lsrf => null() 
    
    !*FD The elevation of the lower ice surface, divided by \texttt{thck0}.

    real(dp),dimension(:,:),pointer :: topg => null() 
    
    !*FD The elevation of the topography, divided by \texttt{thck0}.

    real(dp),dimension(:,:),pointer :: relx => null()
     
    !*FD The elevation of the relaxed topography, by \texttt{thck0}.

    integer, dimension(:,:),pointer :: mask => null()
    
    !*FD Set to zero for all points where $\mathtt{thck}=0$, otherwise non-zero.
    !*FD the non-zero points are numbered in sequence from the bottom left to the 
    !*FD top right, going along the rows.

    integer :: totpts = 0
    
    !*FD The total number of points with non-zero thickness

    integer, dimension(4) :: dom   = 0      !*FD I have no idea what this is for.
    logical               :: empty = .true. !*FD I have no idea what this is for.

  end type glimmer_geometry

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_geomderv

    !*FD Holds the horizontal and temporal derivatives of the thickness and
    !*FD upper surface elevation, as well as the thickness on the staggered grid.

    real(dp),dimension(:,:),pointer :: dthckdew => null() !*FD E-W derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdew => null() !*FD E-W derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: dthckdns => null() !*FD N-S derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdns => null() !*FD N-S derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: dthckdtm => null() !*FD Temporal derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdtm => null() !*FD Temporal derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: stagthck => null() !*FD Thickness averaged onto the staggered grid.

  end type glimmer_geomderv

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_velocity

    !*FD Holds the velocity fields in 2D and 3D. At least some of these fields
    !*FD are stored on the displaced grid.

    real(dp),dimension(:,:,:),pointer :: uvel  => null() !*FD 3D $x$-velocity.
    real(dp),dimension(:,:,:),pointer :: vvel  => null() !*FD 3D $y$-velocity.
    real(dp),dimension(:,:,:),pointer :: wvel  => null() !*FD 3D $z$-velocity.
    real(dp),dimension(:,:,:),pointer :: wgrd  => null() !*FD 3D grid vertical velocity.
    real(dp),dimension(:,:)  ,pointer :: uflx  => null() !*FD 
    real(dp),dimension(:,:)  ,pointer :: vflx  => null() !*FD 
    real(dp),dimension(:,:)  ,pointer :: diffu => null() !*FD 
    real(dp),dimension(:,:)  ,pointer :: ubas  => null() !*FD 
    real(dp),dimension(:,:)  ,pointer :: vbas  => null() !*FD 
    real(dp),dimension(:,:)  ,pointer :: btrc  => null() !*FD 
  end type glimmer_velocity

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_climate

    !*FD Holds fields relating to the atmospheric climate of the model.

    real(sp),dimension(:,:),pointer :: acab     => null() !*FD Annual mass balance.
    real(sp),dimension(:,:),pointer :: artm     => null() !*FD Annual mean air temperature
    real(sp),dimension(:,:),pointer :: arng     => null() !*FD Annual air temperature range
    real(sp),dimension(:,:),pointer :: lati     => null() !*FD Latitudes of model grid points
    real(sp),dimension(:,:),pointer :: ablt     => null() !*FD Annual ablation.
    real(sp),dimension(:,:),pointer :: g_arng   => null() !*FD Global annual air temperature range
    real(sp),dimension(:,:),pointer :: g_artm   => null() !*FD Global annual mean air temperature
    real(sp),dimension(:,:),pointer :: prcp     => null() !*FD Annual precipitation.
    real(sp),dimension(:,:),pointer :: presprcp => null() !*FD Present-day annual precip (mm)
    real(sp),dimension(:,:),pointer :: presartm => null() !*FD Present-day mean air-temp ($^{\circ}$C)
    real(dp),dimension(:,:),pointer :: presusrf => null() !*FD Present-day upper surface (km)
    integer, dimension(:,:),pointer :: out_mask => null() !*FD Array indicating whether a point 
                                                          !*FD should be considered or ignored 
                                                          !*FD when upscaling data for output. 
                                                          !*FD 1 means use, 0 means ignore.
    real(sp) :: uprecip_rate =   0.5 !*FD Uniform precipitaion rate in m/a
    real(sp) :: usurftemp    = -20.0 !*FD Uniform surface temperature in $^{\circ}$C.
    real(sp) :: ice_albedo   =   0.4 !*FD Ice albedo. (fraction)
    real(sp) :: ulapse_rate  =  -8.0 !*FD Uniform lapse rate in deg C/km 
                                     !*FD (N.B. This should be \emph{negative}!)
  end type glimmer_climate

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_temper

    !*FD Holds fields relating to temperature.

    real(dp),dimension(:,:,:),pointer :: temp => null() !*FD Three-dimensional temperature field.
    real(dp),dimension(:,:,:),pointer :: flwa => null() !*FD Glenn's $A$.
    real(dp),dimension(:,:),  pointer :: bwat => null() !*FD Basal water depth(?)
    real(dp),dimension(:,:),  pointer :: bmlt => null() !*FD Basal melt-rate(?)
    integer  :: niter   = 0      !*FD
    real(sp) :: perturb = 0.0    !*FD
    real(sp) :: grid    = 0.0    !*FD
    integer  :: tpt     = 0      !*FD Pointer to time series data
    logical  :: first1  = .true. !*FD
  end type glimmer_temper

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_forcdata

    !*FD Holds information relating to time-series forcing.

    real(sp),dimension(:,:),pointer :: forcing => null() !*FD Forcing data(?)
    integer  :: flines = 0
    real(sp) :: trun   = 10000
  end type glimmer_forcdata

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_funits
    character(fname_length) :: sigfile                         !*FD sigma coordinates file
    character(fname_length) :: forcfile                        !*FD Temperature forcing file
    character(fname_length) :: ncfile                          !*FD configuration file for netCDF I/O
    type(glimmer_nc_output),pointer :: out_first=>NULL()       !*FD first element of linked list defining netCDF outputs
    type(glimmer_nc_input), pointer :: in_first=>NULL()        !*FD first element of linked list defining netCDF inputs
  end type glimmer_funits

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_numerics

    !*FD Parameters relating to the model numerics.

    real(sp) :: time   =    0.0   !*FD main loop counter in years
    real(sp) :: tinc   =   20.0   !*FD time step of main loop in years 
    real(sp) :: ntem   =    1.0   !*FD temperature time step in years
    real(sp) :: nvel   =    1.0   !*FD velocity time step in years
    real(sp) :: niso   =    1.0   !*FD flexure time step in years
    real(sp) :: nstr   =    0.0   !*FD output start time in years
    real(dp) :: alpha  =    0.5d0 !*FD richard suggests 1.5 - was a parameter in original
    real(dp) :: alphas =    0.5d0 !*FD was a parameter in the original
    real(dp) :: thklim =  100.0   
    real(dp) :: mlimit = -200.0d0
    real(dp) :: dew    =   20.0d3
    real(dp) :: dns    =   20.0d3
    real(dp) :: dt     =    0.0
    real(dp) :: dttem  =    0.0
    real(sp) :: nshlf  =    0.0

    ! Data output frequency -------------------------------------------------

    real(sp),dimension(3)         :: nout = (/1.0,10.0,10.0/)
    
    !*FD output time step in years for time series, horiz. and full data 
                                          
    ! Vertical coordinate ---------------------------------------------------
                                                               
    real(dp),dimension(:),pointer :: sigma => null() !*FD Sigma values for 
                                                     !*FD vertical spacing of 
                                                     !*FD model levels

  end type glimmer_numerics

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_pddcalc

    !*FD Holds parameters for positive-degree-day mass-balance
    !*FD calculation. The table has two axes - the $x$ axis is the
    !*FD difference between mean annual and July temps, while the
    !*FD $y$- axis is the mean annual temp

    integer  :: dx        = 1   !*FD Spacing of values in x-direction ($^{\circ}$C)
    integer  :: dy        = 1   !*FD Spacing of values in y-direction ($^{\circ}$C)
    integer  :: ix        = 0   !*FD Lower bound of $x$-axis ($^{\circ}$C)
    integer  :: iy        = -50 !*FD Lower bound of $y$-axis ($^{\circ}$C)
    integer  :: nx        = 31  !*FD Number of values in x-direction
    integer  :: ny        = 71  !*FD Number of values in y-direction
    real(sp) :: dailytemp = 0.0 
    real(sp) :: tma       = 0.0
    real(sp) :: tmj       = 0.0
    real(sp) :: dtmj      = 0.0
    real(sp) :: dd_sigma  = 5.0
    logical  :: first     = .true.
    logical  :: pt_alloc  = .false. !*FD set \texttt{.true.} if \texttt{pddtab}
                                    !*FD has been allocated.
 
    ! The actual PDD table ---------------------------------------------

    real(sp),dimension(:,:),pointer :: pddtab  => null() 
    
    !*FD PDD table - must be allocated with dimensions nx,ny.

 end type glimmer_pddcalc

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_isotwk
    real(dp),dimension(:,:),pointer :: load  => null()
    real(sp),dimension(:,:),pointer :: dflct => null()
    real(dp),dimension(4)           :: fact  =  0.0
    integer :: nflx   = 0
    logical :: first1 = .true.
    logical :: first2 = .true.
  end type glimmer_isotwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_velowk
    real(dp),dimension(:),  pointer :: depth    => null()
    real(dp),dimension(:),  pointer :: dupsw    => null()
    real(dp),dimension(:),  pointer :: depthw   => null()
    real(dp),dimension(:),  pointer :: suvel    => null()
    real(dp),dimension(:),  pointer :: svvel    => null()
    real(dp),dimension(:,:),pointer :: fslip    => null()
    real(dp),dimension(:,:),pointer :: dintflwa => null()
    real(dp),dimension(:),  pointer :: dups     => null()
    real(dp),dimension(4) :: fact
    real(dp),dimension(4) :: c    = 0.0
    real(dp) :: watwd  = 3.0d0
    real(dp) :: watct  = 10.0d0
    real(dp) :: trc0   = 0.0
    real(dp) :: trcmin = 0.0d0
    real(dp) :: marine = 1.0d0
    real(dp) :: trcmax = 10.0d0

    ! Initialisation flags ---------------------------------------

    logical :: first1 = .true.
    logical :: first2 = .true.
    logical :: first3 = .true.
    logical :: first4 = .true.
    logical :: first5 = .true.
    logical :: first6 = .true.
  end type glimmer_velowk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_pcgdwk
    real(dp),dimension(:),pointer :: pcgval  => null()
    real(dp),dimension(:),pointer :: rhsd    => null()
    real(dp),dimension(:),pointer :: answ    => null()
    integer, dimension(:),pointer :: pcgcol  => null()
    integer, dimension(:),pointer :: pcgrow  => null()
    integer, dimension(2)         :: pcgsize = 0
    real(dp),dimension(4)         :: fc      = 0.0
    real(dp),dimension(6)         :: fc2     = 0.0
    integer :: ct     = 0
    integer :: mlinit = 0
    integer :: tlinit = 0
    logical :: first1 = .true.
  end type glimmer_pcgdwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_thckwk
    real(dp),dimension(:,:),  pointer :: oldthck   => null()
    real(dp),dimension(:,:),  pointer :: basestate => null()
    real(dp),dimension(:,:,:),pointer :: olds      => null()
    integer  :: nwhich  = 2
    real(sp) :: oldtime = 0.0
    real(dp) :: few     = 0.0
    real(dp) :: fns     = 0.0
    logical  :: first1  =.true.
  end type glimmer_thckwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_tempwk
    real(dp),dimension(:,:,:),pointer :: inittemp => null()
    real(dp),dimension(:,:,:),pointer :: dissip   => null()
    real(dp),dimension(:,:,:),pointer :: initadvt => null()
    real(dp),dimension(:),    pointer :: dupa     => null()
    real(dp),dimension(:),    pointer :: dupb     => null()
    real(dp),dimension(:),    pointer :: dupc     => null()
    real(dp),dimension(:),    pointer :: c1       => null()
    real(dp),dimension(:,:),  pointer :: dups     => null()
    real(dp),dimension(4)             :: cons     = 0.0
    real(dp),dimension(4)             :: f        = 0.0
    real(dp),dimension(8)             :: c        = 0.0
    real(dp) :: noflow      = -1
    real(dp) :: advconst(2) = 0.0
    real(dp) :: zbed        = 0.0
    real(dp) :: dupn        = 0.0
    real(dp) :: wmax        = 0.0
    real(dp) :: dt_wat      = 0.0
    integer  :: nwat        = 0
    logical  :: first1      = .true.
    logical  :: first2      = .true.
    logical  :: first3      = .true.
    logical  :: first4      = .true.
    logical  :: first5      = .true.
  end type glimmer_tempwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_paramets
    real(sp),dimension(2) :: airt = (/ -3.150, -1.0e-2 /)       ! K, K km^{-3}
    real(sp),dimension(3) :: nmsb = (/ 0.5, 1.05e-5, 450.0e3 /) ! m yr^{-1}, yr^{-1}, m                
    real(dp),dimension(5) :: bpar = (/ 2.0d0, 10.0d0, 10.0d0, 0.0d0, 1.0d0 /)
    real(dp) :: geot   = -5.0d-2  ! W m^{-2}
    real(dp) :: fiddle = 3.0d0    ! -
    real(dp) :: hydtim = 1000.0d0 ! yr^{-1} converted to s^{-1} and scaled, 
                                  ! 0 if no drainage = 0.0d0 * tim0 / scyr
    real(dp) :: isotim = 3000.0d0
  end type glimmer_paramets

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glimmer_global_type
    type(glimmer_general)  :: general
    type(glimmer_options)  :: options
    type(glimmer_geometry) :: geometry
    type(glimmer_geomderv) :: geomderv
    type(glimmer_velocity) :: velocity
    type(glimmer_climate)  :: climate
    type(glimmer_temper)   :: temper
    type(glimmer_forcdata) :: forcdata
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

end module glimmer_types

