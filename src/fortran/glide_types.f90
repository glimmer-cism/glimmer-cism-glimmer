! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_types.f90 - part of the GLIMMER ice model         + 
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

module glide_types

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
  use glimmer_cfproj, only : CFproj_projection

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_general

    !*FD Holds fundamental parameters of the ice model geometry.

    integer :: ewn = 0  !*FD The number of grid-points in the E-W direction.
    integer :: nsn = 0  !*FD The number of grid-points in the N-S direction.
    integer :: upn = 1  !*FD The number of vertical levels in the model.

  end type glide_general

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_options

    !*FD Holds user options controlling the methods used in the ice-model
    !*FD integration.

    integer :: whichtemp = 1

    !*FD Method of ice temperature calculation:
    !*FD \begin{description} 
    !*FD \item[0] Set column to surface air temperature
    !*FD \item[1] Do full temperature solution (also find vertical velocity
    !*FD and apparent vertical velocity)
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

    integer :: whichbtrc = 0

    !*FD Basal slip coefficient: 
    !*FD \begin{description}
    !*FD \item[0] Set equal to zero everywhere
    !*FD \item[1] Set (non--zero) constant
    !*FD \item[2] Set to (non--zero) constant where where temperature is at pressure melting point of ice, otherwise to zero
    !*FD \item[3] \texttt{tanh} function of basal water depth 
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

    integer :: whichrelaxed = 0
    !*FD relaxed topography:
    !*FD \begin{description}
    !*FD \item[0] get relaxed topo from separate variable
    !*FD \item[1] first time slice of input topo is relaxed
    !*FD \end{description}

    integer :: hotstart = 0
    !*FD hotstart the model
    !*FD \begin{description}
    !*FD \item[0] normal start-up
    !*FD \item[1] hotstart model from previous run
    !*FD \end{description}

  end type glide_options

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_geometry

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

    integer, dimension(:,:),pointer :: thkmask => null()
    !*FD see glide_mask.f90 for possible values

    integer :: totpts = 0
    !*FD The total number of points with non-zero thickness

    integer, dimension(4) :: dom   = 0      !*FD I have no idea what this is for.
    logical               :: empty = .true. !*FD I have no idea what this is for.

  end type glide_geometry

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_geomderv

    !*FD Holds the horizontal and temporal derivatives of the thickness and
    !*FD upper surface elevation, as well as the thickness on the staggered grid.

    real(dp),dimension(:,:),pointer :: dthckdew => null() !*FD E-W derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdew => null() !*FD E-W derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: dthckdns => null() !*FD N-S derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdns => null() !*FD N-S derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: dthckdtm => null() !*FD Temporal derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdtm => null() !*FD Temporal derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: stagthck => null() !*FD Thickness averaged onto the staggered grid.

  end type glide_geomderv

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_velocity

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
  end type glide_velocity

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_climate
     !*FD Holds fields used to drive the model
     real(sp),dimension(:,:),pointer :: acab     => null() !*FD Annual mass balance.
     real(sp),dimension(:,:),pointer :: artm     => null() !*FD Annual mean air temperature
     real(sp),dimension(:,:),pointer :: lati     => null() !*FD Latitudes of model grid points
     real(sp),dimension(:,:),pointer :: loni     => null() !*FD Longitudes of model grid points
     real(sp) :: eus = 0.                                  !*FD eustatic sea level
  end type glide_climate

  type glide_temper

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
  end type glide_temper

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_funits
    character(fname_length) :: sigfile=''                      !*FD sigma coordinates file
    character(fname_length) :: ncfile=''                       !*FD configuration file for netCDF I/O
    type(glimmer_nc_output),pointer :: out_first=>NULL()       !*FD first element of linked list defining netCDF outputs
    type(glimmer_nc_input), pointer :: in_first=>NULL()        !*FD first element of linked list defining netCDF inputs
  end type glide_funits

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_numerics

    !*FD Parameters relating to the model numerics.
    real(sp) :: tstart = 0.0      !*FD starting time
    real(sp) :: tend   = 20000.0  !*FD end time
    real(sp) :: time   =    0.0   !*FD main loop counter in years
    real(sp) :: tinc   =   20.0   !*FD time step of main loop in years 
    real(sp) :: ntem   =    1.0   !*FD temperature time step (multiplier of main time step)
    real(sp) :: nvel   =    1.0   !*FD velocity time step (multiplier of main time step)
    real(sp) :: niso   =    1.0   !*FD flexure time step (multiplier of main time step)
    real(dp) :: alpha  =    0.5d0 !*FD richard suggests 1.5 - was a parameter in original
    real(dp) :: alphas =    0.5d0 !*FD was a parameter in the original
    real(dp) :: thklim =  100.0   
    real(dp) :: mlimit = -200.0d0
    real(dp) :: dew    =   20.0d3
    real(dp) :: dns    =   20.0d3
    real(dp) :: dt     =    0.0
    real(dp) :: dttem  =    0.0
    real(sp) :: nshlf  =    0.0
    
    ! Vertical coordinate ---------------------------------------------------
                                                               
    real(dp),dimension(:),pointer :: sigma => null() !*FD Sigma values for 
                                                     !*FD vertical spacing of 
                                                     !*FD model levels

  end type glide_numerics

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_isotwk
    real(dp),dimension(:,:),pointer :: load  => null()
    real(sp),dimension(:,:),pointer :: dflct => null()
    real(dp),dimension(4)           :: fact  =  0.0
    integer :: nflx   = 0
  end type glide_isotwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_velowk
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
    real(dp) :: btrac_const = 0.0d0
  end type glide_velowk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_pcgdwk
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
  end type glide_pcgdwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_thckwk
    real(dp),dimension(:,:),  pointer :: oldthck   => null()
    real(dp),dimension(:,:),  pointer :: basestate => null()
    real(dp),dimension(:,:),pointer :: float => null()
    real(dp),dimension(:,:,:),pointer :: olds      => null()
    integer  :: nwhich  = 2
    real(sp) :: oldtime = 0.0
    real(dp) :: few     = 0.0
    real(dp) :: fns     = 0.0
    logical  :: first1  =.true.
  end type glide_thckwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_tempwk
    logical, dimension(:,:), pointer  :: floater  => null() 
    real(dp),dimension(:,:,:),pointer :: inittemp => null()
    real(dp),dimension(:,:,:),pointer :: dissip   => null()
    real(dp),dimension(:,:,:),pointer :: initadvt => null()
    real(dp),dimension(:),    pointer :: dupa     => null()
    real(dp),dimension(:),    pointer :: dupb     => null()
    real(dp),dimension(:),    pointer :: dupc     => null()
    real(dp),dimension(:),    pointer :: c1       => null()
    real(dp),dimension(:,:),  pointer :: dups     => null()
    real(dp),dimension(:,:),  pointer :: wphi     => null()
    real(dp),dimension(:,:),  pointer :: bwatu    => null()
    real(dp),dimension(:,:),  pointer :: bwatv    => null()
    real(dp),dimension(:,:),  pointer :: fluxew   => null()
    real(dp),dimension(:,:),  pointer :: fluxns   => null()
    real(dp),dimension(:,:),  pointer :: bint     => null()
    real(dp),dimension(:,:),  pointer :: smth     => null()
    real(dp),dimension(4)             :: cons     = 0.0
    real(dp),dimension(4)             :: f        = 0.0
    real(dp),dimension(8)             :: c        = 0.0
    real(dp) :: noflow      = -1
    real(dp) :: advconst(2) = 0.0
    real(dp) :: zbed        = 0.0
    real(dp) :: dupn        = 0.0
    real(dp) :: wmax        = 0.0
    real(dp) :: dt_wat      = 0.0
    real(dp) :: watvel      = 0.0
    integer  :: nwat        = 0
  end type glide_tempwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_paramets
    real(dp),dimension(5) :: bpar = (/ 2.0d0, 10.0d0, 10.0d0, 0.0d0, 1.0d0 /)
    real(dp) :: btrac_const = 0.d0 ! m yr^{-1} Pa^{-1} (gets scaled during init)
    real(dp) :: geot   = -5.0d-2  ! W m^{-2}
    real(dp) :: fiddle = 3.0d0    ! -
    real(dp) :: hydtim = 1000.0d0 ! yr^{-1} converted to s^{-1} and scaled, 
                                  ! 0 if no drainage = 0.0d0 * tim0 / scyr
    real(dp) :: isotim = 3000.0d0
  end type glide_paramets

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_global_type
    type(glide_general)  :: general
    type(glide_options)  :: options
    type(glide_geometry) :: geometry
    type(glide_geomderv) :: geomderv
    type(glide_velocity) :: velocity
    type(glide_climate)  :: climate
    type(glide_temper)   :: temper
    type(glide_funits)   :: funits
    type(glide_numerics) :: numerics
    type(glide_isotwk)   :: isotwk
    type(glide_velowk)   :: velowk
    type(glide_pcgdwk)   :: pcgdwk
    type(glide_thckwk)   :: thckwk
    type(glide_tempwk)   :: tempwk
    type(glide_paramets) :: paramets
    type(CFproj_projection) :: projection
  end type glide_global_type

contains
  
  subroutine glide_allocarr(model)
    
    !*FD Allocates the model arrays, and initialises some of them to zero.
    !*FD These are the arrays allocated, and their dimensions:
    !*FD
    !*FD In \texttt{model\%temper}:
    !*FD \begin{itemize}
    !*FD \item \texttt{temp(upn,0:ewn+1,0:nsn+1))}
    !*FD \item \texttt{flwa(upn,ewn,nsn))}
    !*FD \item \texttt{bwat(ewn,nsn))}
    !*FD \item \texttt{bmlt(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%velocity}:
    !*FD \begin{itemize}
    !*FD \item \texttt{uvel(upn,ewn-1,nsn-1))}
    !*FD \item \texttt{vvel(upn,ewn-1,nsn-1))}
    !*FD \item \texttt{wvel(upn,ewn,nsn))}
    !*FD \item \texttt{wgrd(upn,ewn,nsn))}
    !*FD \item \texttt{uflx(ewn-1,nsn-1))}
    !*FD \item \texttt{vflx(ewn-1,nsn-1))}
    !*FD \item \texttt{diffu(ewn,nsn))}
    !*FD \item \texttt{btrc(ewn,nsn))}
    !*FD \item \texttt{ubas(ewn,nsn))}
    !*FD \item \texttt{vbas(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%climate}:
    !*FD \begin{itemize}
    !*FD \item \texttt{acab(ewn,nsn))}
    !*FD \item \texttt{artm(ewn,nsn))}
    !*FD \item \texttt{lati(ewn,nsn))}
    !*FD \item \texttt{loni(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%geomderv}:
    !*FD \begin{itemize}
    !*FD \item \texttt{dthckdew(ewn,nsn))}
    !*FD \item \texttt{dusrfdew(ewn,nsn))}
    !*FD \item \texttt{dthckdns(ewn,nsn))}
    !*FD \item \texttt{dusrfdns(ewn,nsn))}
    !*FD \item \texttt{dthckdtm(ewn,nsn))}
    !*FD \item \texttt{dusrfdtm(ewn,nsn))}
    !*FD \item \texttt{stagthck(ewn-1,nsn-1))}
    !*FD \end{itemize}
  
    !*FD In \texttt{model\%geometry}:
    !*FD \begin{itemize}
    !*FD \item \texttt{thck(ewn,nsn))}
    !*FD \item \texttt{usrf(ewn,nsn))}
    !*FD \item \texttt{lsrf(ewn,nsn))}
    !*FD \item \texttt{topg(ewn,nsn))}
    !*FD \item \texttt{relx(ewn,nsn))}
    !*FD \item \texttt{mask(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%thckwk}:
    !*FD \begin{itemize}
    !*FD \item \texttt{olds(ewn,nsn,thckwk\%nwhich))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%isotwk}:
    !*FD \begin{itemize}
    !*FD \item \texttt{load(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%numerics}:
    !*FD \begin{itemize}
    !*FD \item \texttt{sigma(upn))}
    !*FD \end{itemize}

    implicit none

    type(glide_global_type),intent(inout) :: model

    integer :: ewn,nsn,upn

    ! for simplicity, copy these values...

    ewn=model%general%ewn
    nsn=model%general%nsn
    upn=model%general%upn

    ! Allocate appropriately

    allocate(model%temper%temp(upn,0:ewn+1,0:nsn+1)); model%temper%temp = 0.0
    allocate(model%temper%flwa(upn,ewn,nsn))   
    allocate(model%temper%bwat(ewn,nsn));             model%temper%bwat = 0.0
    allocate(model%temper%bmlt(ewn,nsn));             model%temper%bmlt = 0.0

    allocate(model%velocity%uvel(upn,ewn-1,nsn-1));   model%velocity%uvel = 0.0d0
    allocate(model%velocity%vvel(upn,ewn-1,nsn-1));   model%velocity%vvel = 0.0d0
    allocate(model%velocity%wvel(upn,ewn,nsn));       model%velocity%wvel = 0.0d0
    allocate(model%velocity%wgrd(upn,ewn,nsn));       model%velocity%wgrd = 0.0d0
    allocate(model%velocity%uflx(ewn-1,nsn-1));       model%velocity%uflx = 0.0d0
    allocate(model%velocity%vflx(ewn-1,nsn-1));       model%velocity%vflx = 0.0d0
    allocate(model%velocity%diffu(ewn,nsn));          model%velocity%diffu = 0.0d0
    allocate(model%velocity%btrc(ewn,nsn));           model%velocity%btrc = 0.0d0
    allocate(model%velocity%ubas(ewn,nsn));           model%velocity%ubas = 0.0d0
    allocate(model%velocity%vbas(ewn,nsn));           model%velocity%vbas = 0.0d0

    allocate(model%climate%acab(ewn,nsn));            model%climate%acab = 0.0
    allocate(model%climate%artm(ewn,nsn));            model%climate%artm = 0.0
    allocate(model%climate%lati(ewn,nsn));            model%climate%lati = 0.0
    allocate(model%climate%loni(ewn,nsn));            model%climate%loni = 0.0

    allocate(model%geomderv%dthckdew(ewn,nsn));       model%geomderv%dthckdew = 0.0d0 
    allocate(model%geomderv%dusrfdew(ewn,nsn));       model%geomderv%dusrfdew = 0.0d0
    allocate(model%geomderv%dthckdns(ewn,nsn));       model%geomderv%dthckdns = 0.0d0
    allocate(model%geomderv%dusrfdns(ewn,nsn));       model%geomderv%dusrfdns = 0.0d0
    allocate(model%geomderv%dthckdtm(ewn,nsn));       model%geomderv%dthckdtm = 0.0d0
    allocate(model%geomderv%dusrfdtm(ewn,nsn));       model%geomderv%dusrfdtm = 0.0d0
    allocate(model%geomderv%stagthck(ewn-1,nsn-1));   model%geomderv%stagthck = 0.0d0
  
    allocate(model%geometry%thck(ewn,nsn));           model%geometry%thck = 0.0d0
    allocate(model%geometry%usrf(ewn,nsn));           model%geometry%usrf = 0.0d0
    allocate(model%geometry%lsrf(ewn,nsn));           model%geometry%lsrf = 0.0d0
    allocate(model%geometry%topg(ewn,nsn));           model%geometry%topg = 0.0d0
    allocate(model%geometry%relx(ewn,nsn));           model%geometry%relx = 0.0d0
    allocate(model%geometry%mask(ewn,nsn));           model%geometry%mask = 0
    allocate(model%geometry%thkmask(ewn,nsn));        model%geometry%thkmask = 0

    allocate(model%thckwk%olds(ewn,nsn,model%thckwk%nwhich))
                                                      model%thckwk%olds = 0.0d0
    allocate(model%thckwk%oldthck(ewn,nsn));          model%thckwk%oldthck = 0.0d0
    allocate(model%thckwk%basestate(ewn,nsn));        model%thckwk%basestate = 0.0d0
    allocate(model%isotwk%load(ewn,nsn));             model%isotwk%load = 0.0d0 
    allocate(model%numerics%sigma(upn))

  end subroutine glide_allocarr

  subroutine glide_deallocarr(model)
    !*FD deallocate model arrays
    implicit none
    type(glide_global_type),intent(inout) :: model

    deallocate(model%temper%temp)
    deallocate(model%temper%flwa)
    deallocate(model%temper%bwat)
    deallocate(model%temper%bmlt)

    deallocate(model%velocity%uvel)
    deallocate(model%velocity%vvel)
    deallocate(model%velocity%wvel)
    deallocate(model%velocity%wgrd)
    deallocate(model%velocity%uflx)
    deallocate(model%velocity%vflx)
    deallocate(model%velocity%diffu)
    deallocate(model%velocity%btrc)
    deallocate(model%velocity%ubas)
    deallocate(model%velocity%vbas)

    deallocate(model%climate%acab)
    deallocate(model%climate%artm)
    deallocate(model%climate%lati)
    deallocate(model%climate%loni)

    deallocate(model%geomderv%dthckdew)
    deallocate(model%geomderv%dusrfdew)
    deallocate(model%geomderv%dthckdns)
    deallocate(model%geomderv%dusrfdns)
    deallocate(model%geomderv%dthckdtm)
    deallocate(model%geomderv%dusrfdtm)
    deallocate(model%geomderv%stagthck)
  
    deallocate(model%geometry%thck)
    deallocate(model%geometry%usrf)
    deallocate(model%geometry%lsrf)
    deallocate(model%geometry%topg)
    deallocate(model%geometry%relx)
    deallocate(model%geometry%mask)
    deallocate(model%geometry%thkmask)

    deallocate(model%thckwk%olds)
    deallocate(model%thckwk%oldthck)
    deallocate(model%thckwk%basestate)
    deallocate(model%isotwk%load)
    deallocate(model%numerics%sigma)
  end subroutine glide_deallocarr
end module glide_types

