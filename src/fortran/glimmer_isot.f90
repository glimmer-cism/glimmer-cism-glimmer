! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_isot.f90 - part of the GLIMMER ice model         + 
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

module glide_isot

  !*FD This module contains subroutines to calculate the isostatic
  !*FD depression/rebound due to the weight of ice and/or water on the
  !*FD Earth's crust.

  use glide_types

  private flextopg

contains

!----------------------------------------------------------------------
! PUBLIC subroutines
!----------------------------------------------------------------------

  subroutine isosevol(numerics,paramets,isotwk,which,thck,evol,rtop)

    !*FD Calculates the elevation of the bedrock topography, by considering
    !*FD the isostatic depression caused by the over-lying ice and water.

    use glimmer_global, only : dp
    use physcon, only : rhoi, rhom

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
 
    type(glide_numerics),   intent(in)    :: numerics !*FD Numerical parameters
    type(glide_paramets),   intent(in)    :: paramets !*FD Other model parameters
    type(glide_isotwk),     intent(inout) :: isotwk   !*FD Isostasy fields and parameters
    real(dp), dimension(:,:), intent(inout) :: thck     !*FD Ice thickness field (scaled by $1/\mathtt{thck0}$).
    real(dp), dimension(:,:), intent(inout) :: evol     !*FD The current elevation of the 
                                                        !*FD topography/bathymetry (scaled by $1/\mathtt{thck0}$).
    real(dp), dimension(:,:), intent(inout) :: rtop     !*FD The elevation of the topography
                                                        !*FD in a relaxed state (scaled by $1/\mathtt{thck0}$).
    integer,                  intent(in)    :: which    !*FD Option determining which method
                                                        !*FD to use in the calculation.

    !*FD
    !*FD The possible range of values of \texttt{which} is as follows:
    !*FD \begin{itemize}
    !*FD \item \texttt{which=0} The topography is fixed. Somewhat bizarrely, this is accomplished by
    !*FD setting the relaxed topography to be equal to the present state.
    !*FD \item \texttt{which=1} Use a local function of ice thickness.
    !*FD \item \texttt{which=2} Use a (different) local function of ice thickness, incorporating
    !*FD flexure.
    !*FD \end{itemize}

    ! -----------------------------------------------------------------------------
    ! Initialisation, if this is the first call
    ! -----------------------------------------------------------------------------

    if (isotwk%first1) then
      isotwk%fact = (/ numerics%dt / paramets%isotim, &
                       1.0d0 - 0.5d0 * numerics%dt / paramets%isotim, &
                       1.0d0 + 0.5d0 * numerics%dt / paramets%isotim, &
                       rhoi/rhom /)
      isotwk%first1 = .false.
    end if

    ! -----------------------------------------------------------------------------
    ! Choose one of three possible methods
    ! -----------------------------------------------------------------------------

    select case(which)
    case(0)
      rtop = evol    
    case(1)
      evol = (isotwk%fact(1) * (rtop - isotwk%fact(4)*thck) + &
              isotwk%fact(2) * evol) / isotwk%fact(3)
    case(2) 
      if ( 0 == mod(numerics%time,numerics%niso) ) then
        call flextopg(isotwk,numerics,isotwk%load,thck)
      end if
      evol = (isotwk%fact(1) * (rtop + isotwk%load) + &
              isotwk%fact(2) * evol) / isotwk%fact(3)
    end select 

  end subroutine isosevol

!---------------------------------------------------------------
! PRIVATE subroutines
!---------------------------------------------------------------

  subroutine flextopg(isotwk,numerics,flex,thck)

    !*FD This is the subroutine called for the case \texttt{whichisot}=2.
    !*FD It calculates the load due to ice (and possibly water)

    use glimmer_global, only : dp, sp
    use paramets, only : len0
    use physcon, only : rhoi, rhom, grav, pi

    implicit none

    !-----------------------------------------------------------
    ! Subroutine arguments
    !-----------------------------------------------------------

    type(glide_isotwk),     intent(inout) :: isotwk    !*FD Isostasy work arrays
    type(glide_numerics),   intent(in)    :: numerics  !*FD Model numerics parameters
    real(dp), dimension(:,:), intent(out)   :: flex      !*FD The load due to ice and
                                                         !*FD water (?)
    real(dp), dimension(:,:), intent(in)    :: thck      !*FD Current ice thickness (scaled)
    
    !-----------------------------------------------------------
    ! Internal variables
    !-----------------------------------------------------------

    integer :: ns,ew,nsn,ewn
    integer :: ewpt, nspt, ewflx, nsflx, ikelv
    integer, parameter :: nkelv = 82 
      
    ! ----------------------------------------------------------
    ! ** constants used in calc
    ! **
    ! ** Youngs modulus (Nm^-2)

    real(sp),parameter :: youngs = 8.35e10

    ! ** thickness of lithosphere (m)

    real(sp),parameter :: thklith = 110.0e3

    ! ** radius of lithosphere (m)

    real(sp),parameter :: radlith = 6.244e6

    ! ** Poissons ratio

    real(sp),parameter :: poiss = 0.25

    ! ----------------------------------------------------------
    ! ** Zero order kelvin function (for every dkelv from 0.0 to 8.0)

    real(sp),parameter ::  dkelv = 0.1 

    ! ** Look-up table of values:

    real(sp),parameter, dimension(nkelv) :: &
    kelvin0 = (/ -0.785, -0.777, -0.758, -0.733, -0.704, &
                 -0.672, -0.637, -0.602, -0.566, -0.531, &
                 -0.495, -0.460, -0.426, -0.393, -0.362, &
                 -0.331, -0.303, -0.275, -0.249, -0.225, &
                 -0.202, -0.181, -0.161, -0.143, -0.126, &
                 -0.111, -0.096, -0.083, -0.072, -0.061, &
                 -0.051, -0.042, -0.035, -0.028, -0.021, &
                 -0.016, -0.011, -0.007, -0.003,  0.000, & 
                  0.002,  0.004,  0.006,  0.008,  0.009, & 
                  0.010,  0.010,  0.011,  0.011,  0.011, & 
                  0.011,  0.011,  0.011,  0.011,  0.010, & 
                  0.010,  0.009,  0.009,  0.008,  0.008, & 
                  0.007,  0.007,  0.006,  0.006,  0.005, & 
                  0.005,  0.004,  0.004,  0.003,  0.003, & 
                  0.003,  0.002,  0.002,  0.002,  0.002, & 
                  0.001,  0.001,  0.001,  0.001,  0.001, & 
                  0.000,  0.000 /)

    ! ----------------------------------------------------------
    ! ** quantities calculated
    ! **
    ! ** flexural rigidity
    ! ** radius of stiffness
    ! ** multiplier for loads
     
    real(sp),save :: rigid, alpha, multi, dist, load

    ! ----------------------------------------------------------

    ewn=size(thck,1) ; nsn=size(thck,2)

    if (isotwk%first2) then                                                  

      rigid = (youngs * thklith**3) / (12.0 * (1.0 - poiss**2))
      alpha = (rigid / ((youngs * thklith / radlith**2) + rhom * grav))**0.25
      multi = grav * (numerics%dew * len0)**2 * alpha**2 / (2.0 * pi * rigid)
      isotwk%nflx = 7 * int(alpha / (numerics%dew * len0)) + 1

      allocate(isotwk%dflct(isotwk%nflx,isotwk%nflx))

      do nsflx = 1,isotwk%nflx
        do ewflx = 1,isotwk%nflx
          dist = len0 * numerics%dew * sqrt(real(ewflx-1)**2 + real(nsflx-1)**2) / alpha
          ikelv = min(nkelv-1,int(dist/dkelv) + 1)
          isotwk%dflct(ewflx,nsflx) = multi * &
                              (kelvin0(ikelv) + &
                              (dist - dkelv * (ikelv-1)) * &
                              (kelvin0(ikelv+1)- kelvin0(ikelv)) / dkelv)
        end do
      end do

      isotwk%first2 = .false.
         
    end if

    flex = 0.0

    ! ** now loop through all of the points on the model grid

    do ns = 1,nsn
      do ew = 1,ewn
        if (thck(ew,ns) > numerics%thklim) then

          ! ** for each point find the load it is imposing which
          ! ** depends on whether there it is ice covered or ocean
          ! ** covered

          load = rhoi * thck(ew,ns)

          ! ** now apply the calculated deflection field using the 
          ! ** ice/ocean load at that point and the function of
          ! ** how it will affect its neighbours

          ! ** the effect is linear so that we can sum the deflection
          ! ** from all imposed loads

          ! ** only do this if there is a load to impose and 
          ! ** be careful not to extend past grid domain
                 
          do nsflx = max(1.0,real(ns)-isotwk%nflx+1), min(nsn,nint(real(ns)+isotwk%nflx-1.0))
            do ewflx = max(1.0,real(ew)-isotwk%nflx+1), min(ewn,nint(real(ew)+isotwk%nflx-1.0))

              ! ** find the correct function value to use
              ! ** the array dflct is one quadrant of the a square
              ! ** centered at the point imposing the load

              nspt = abs(ns - nsflx) + 1; ewpt = abs(ew - ewflx) + 1
              flex(ewflx,nsflx) = load * isotwk%dflct(ewpt,nspt) + flex(ewflx,nsflx)
            end do
          end do

        end if

      end do
    end do
         
  end subroutine flextopg

end module glide_isot
