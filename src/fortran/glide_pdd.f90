
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_pdd.f90 - part of the GLIMMER ice model            + 
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

module glint_pdd

  !*FD The GLINT annual positive degree day mass-balance scheme.
  !*FD Based on the original pdd mass-balance code from Tony's model.
  !*FD 
  !*FD {\bf N.B.} the module variables in this module are used for back-door 
  !*FD message passing, to make the integration of the PDD table look more 
  !*FD comprehensible, and avoid the need to have two customised copies of 
  !*FD the integration code.

  use glimmer_global

  implicit none

  type glint_pdd_params

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
    real(sp) :: dd_sigma  = 5.0 !*FD Standard deviation of daily temperature (K)
    logical  :: first     = .true.
    logical  :: pt_alloc  = .false. !*FD set \texttt{.true.} if \texttt{pddtab}
                                    !*FD has been allocated.
 
    ! The actual PDD table ---------------------------------------------

    real(sp),dimension(:,:),pointer :: pddtab  => null() 
    
    !*FD PDD table - must be allocated with dimensions nx,ny.

    ! Parameters for the PDD calculation

    real(sp) :: pddfs               !*FD Later set to \texttt{(rhow / rhoi) * pddfac\_snow / (acc0 * scyr)}
    real(sp) :: pddfi               !*FD Later set to \texttt{(rhow / rhoi) * pddfac\_ice  / (acc0 * scyr)}
    real(sp) :: wmax        = 0.6   !*FD Fraction of melted snow that refreezes
    real(dp) :: pddfac_ice  = 0.008 !*FD PDD factor for ice (m day$^{-1}$ $^{\circ}C$^{-1}$)
    real(dp) :: pddfac_snow = 0.003 !*FD PDD factor for snow (m day$^{-1}$ $^{\circ}C$^{-1}$)

    ! andreas values for pdd coeffs
    ! pddfs 0.005 and pddfi 0.016
    ! eismint pddfs 0.003 and pddfi 0.008
    ! ritz et al 1997 pddfs 0.005 and pddfi 0.016

  end type glint_pdd_params
 	 
  real(sp) :: dd_sigma            !*FD The value of $\sigma$ in the PDD integral
  real(sp) :: t_a_prime           !*FD The value of $T'_{a}$ in the PDD integral
  real(sp) :: mean_annual_temp    !*FD Mean annual temperature
  real(sp) :: mean_july_temp      !*FD Mean july temperature

  private
  public :: glint_pdd_params, glint_pdd_init, glint_pdd_mbal

contains

!-------------------------------------------------------------------------------
! PUBLIC subroutines
!-------------------------------------------------------------------------------

  subroutine glint_pdd_init(params)

    type(glint_pdd_params),intent(inout) :: params !*FD The positive-degree-day parameters

    ! Deallocate arrays

    if (associated(params%pddtab)) deallocate(params%pddtab)
    
    ! Allocate pdd table

    allocate(params%pddtab(params%nx,params%ny))
    call pddtabgrn(params)  

  end subroutine glint_pdd_init

!-------------------------------------------------------------------------------

  subroutine glint_pdd_mbal(pddcalc,artm,arng,prcp,lati,ablt,acab)

    !*FD Calculates mass-balance over the ice model domain, by the
    !*FD positive-degree-day method.

    implicit none 
 
    type(glint_pdd_params),   intent(inout) :: pddcalc !*FD The positive-degree-day parameters
    real(sp), dimension(:,:), intent(in)    :: artm    !*FD Annual mean air-temperature 
                                                       !*FD ($^{\circ}$C)
    real(sp), dimension(:,:), intent(in)    :: arng    !*FD Annual temerature half-range ($^{\circ}$C)
    real(sp), dimension(:,:), intent(in)    :: prcp    !*FD Annual accumulated precipitation 
                                                       !*FD (mm water equivalent)
    real(sp), dimension(:,:), intent(in)    :: lati    !*FD Latitudes of each point in the 
                                                       !*FD domain ($^{\circ}$N)
    real(sp), dimension(:,:), intent(out)   :: ablt    !*FD Annual ablation (mm water equivalent)
    real(sp), dimension(:,:), intent(out)   :: acab    !*FD Annual mass-balance (mm water equivalent)

    ! Internal variables

    real(sp) :: wfrac, pablt, tx, ty, pdd
    integer  :: ns,ew,nsn,ewn,kx,ky,jx,jy

    ! Get size of arrays. All arrays should be the same size as this.

    ewn=size(artm,1) ; nsn=size(artm,2)

    !-----------------------------------------------------------------------
    ! Main loop
    !-----------------------------------------------------------------------

    do ns = 1, nsn
      do ew = 1, ewn 

        if (lati(ew,ns) .gt. 0.0) then

          ! Find the no. of pdd from the mean annual temp and its range

          ky = int((artm(ew,ns)-pddcalc%iy)/pddcalc%dy)
          kx = int((arng(ew,ns)-pddcalc%ix)/pddcalc%dx) 
    
          ! Check to see if indicies are in range

          if ( kx < 0 ) then 
            tx = 0
            jx = 2
            kx = 1
          else if ( kx > pddcalc%nx-2 ) then
            tx = 1.0
            jx = pddcalc%nx
            kx = pddcalc%nx-1
          else
            tx = arng(ew,ns) - kx * pddcalc%dx - pddcalc%ix
            jx = kx + 2
            kx = kx + 1
          end if

          if ( ky < 0 ) then 
            ty = 0.0
            jy = 2
            ky = 1
          else if ( ky > pddcalc%ny-2 ) then
            ty = 1.0
            jy = pddcalc%ny
            ky = pddcalc%ny-1
          else
            ty = artm(ew,ns) - ky * pddcalc%dy - pddcalc%iy;
            jy = ky + 2
            ky = ky + 1
          end if
            
          ! this is done using a look-up table constructed earlier

          pdd = pddcalc%pddtab(kx,ky)*(1.0-tx)*(1.0-ty) + &
                pddcalc%pddtab(jx,ky) * tx * (1.0 - ty) + &
                pddcalc%pddtab(jx,jy) * tx * ty +         &
                pddcalc%pddtab(kx,jy) * (1.0 - tx) * ty

          ! now start to find the actual net annual accumulation
          ! correct prcpitation for changes in air temperature
          ! REMOVED as we are taking precip as an input

          ! prcp(ew,ns) = climate%presprcp(ew,ns) * &
          !              pfac ** (artm(ew,ns) - climate%presartm(ew,ns))
 
          ! this is the depth of superimposed ice that would need to be
          ! melted before runoff can occur (prcp is already scaled)

          wfrac = pddcalc%wmax * prcp(ew,ns)

          ! this is the total potential ablation of SNOW
          ! note we convert to scaled ablation
    
          pablt = pdd * pddcalc%pddfs

          ! if the total snow ablation is less than the depth of 
          ! superimposed ice - no runoff occurs

          ! else if the total snow ablation is more than the depth
          ! of superimposed ice BUT less than the total amount of
          ! prcpitation - runoff occurs (at a rate equal to the
          ! total potential snowmelt minus that which forms superimposed ice)

          ! else if the total snow ablation is more than the amount
          ! of prcpitation - all snow that is not superimposed ice is lost 
          ! and the potential ablation not used on snow is used on ice
          ! (including the superimposed ice)

          ! there is a change in the pddfi term, replaced wfrac with prcp
          ! error spotted by jonathan 18-04-00

          if ( pablt <= wfrac ) then
            ablt(ew,ns) = 0.0
          else if(pablt > wfrac .and.pablt <= prcp(ew,ns)) then   
            ablt(ew,ns) = pablt - wfrac 
          else
            ablt(ew,ns) = prcp(ew,ns) - wfrac + pddcalc%pddfi*(pdd-prcp(ew,ns)/pddcalc%pddfs) 
          end if

          ! Finally, mass-balance is difference between accumulation and
          ! ablation.

          acab(ew,ns) = prcp(ew,ns) - ablt(ew,ns)
        else
          acab(ew,ns) = 0.0
          ablt(ew,ns) = 0.0
        end if
      end do
    end do

  end subroutine glint_pdd_mbal

!-------------------------------------------------------------------------------
! PRIVATE subroutines and functions
!-------------------------------------------------------------------------------

  subroutine pddtabgrn(pddcalc)

    !*FD Initialises the positive-degree-day-table.

    use glimmer_global, only: sp
    use paramets, only: scyr, acc0
    use physcon, only: rhoi,rhow
    use glimmer_log

    implicit none

    type(glint_pdd_params),intent(inout) :: pddcalc !*FD PDD parameters

    ! Internal variables

    real(sp)           :: tma,dtmj
    real(sp),parameter :: twopi = 3.1416 * 2.0 
    real(sp)           :: fac
    integer            :: mfin
    integer  :: kx,ky

    ! Initialise a couple of constants

    pddcalc%pddfs = (rhow / rhoi) * pddcalc%pddfac_snow / (acc0 * scyr)
    pddcalc%pddfi = (rhow / rhoi) * pddcalc%pddfac_ice  / (acc0 * scyr)

    !--------------------------------------------------------------------
    ! Main loops:
    !  tma  -- the mean annual temperature (y-axis of pdd table)
    !  dtmj -- difference from mean july temperature (x-axis of table)
    !  tmj -- the actual july temperature
    !--------------------------------------------------------------------

    call write_log('Calculating PDD table...',GM_DIAGNOSTIC)

    do tma = pddcalc%iy, pddcalc%iy+(pddcalc%ny-1)*pddcalc%dy, pddcalc%dy

       ky = findgrid(tma,real(pddcalc%iy),real(pddcalc%dy))

       do dtmj = pddcalc%ix, pddcalc%ix+(pddcalc%nx-1)*pddcalc%dx, pddcalc%dx

          mean_july_temp = tma + dtmj   
          kx  = findgrid(dtmj,real(pddcalc%ix),real(pddcalc%dx)) 

          ! need these lines to take account of the backdoor message passing used here

          mean_annual_temp=tma
          dd_sigma=pddcalc%dd_sigma

          pddcalc%pddtab(kx,ky)=(1.0/(dd_sigma*sqrt(twopi)))*romberg_int(inner_integral,0.0,twopi)

          ! convert to days     

          pddcalc%pddtab(kx,ky) = 365.0 * pddcalc%pddtab(kx,ky) / twopi

       end do
    end do

    call write_log('   ...done.',GM_DIAGNOSTIC)

  end subroutine pddtabgrn

!-------------------------------------------------------------------------------
        
  real(sp) function inner_integral(day)

    !*FD Calculates the value of the inner integral, i.e.
    !*FD \begin{equation}
    !*FD \int^{T_{a}'+2.5\sigma}_{0}T_{a}\times
    !*FD \exp\left(\frac{-(T_a-T_{a}')^2}{2\sigma^2}\right)\,dT
    !*FD \end{equation}

    real(sp) :: day !*FD The `day', in radians, so that a year is $2\pi$ long.

    real(sp) :: upper_limit,fac
    integer  :: mfin

    t_a_prime=mean_annual_temp+(mean_july_temp-mean_annual_temp)*cos(day)

    upper_limit=t_a_prime+2.5*dd_sigma

    if (upper_limit<=0.0) then
      inner_integral=0.0
    else
      inner_integral=romberg_int(pdd_integrand,0.0,upper_limit)
    endif

  end function inner_integral

!-------------------------------------------------------------------------------
        
  real(sp) function pdd_integrand(artm)

    !*FD The expression to be integrated in the calculation of the PDD table. The whole
    !*FD integral is:
    !*FD \begin{equation}
    !*FD D=\frac{1}{\sigma\sqrt{2\pi}}\int^{A}_{0}\int^{T_{a}'+2.5\sigma}_{0}T_{a}\times
    !*FD \exp\left(\frac{-(T_a-T_{a}')^2}{2\sigma^2}\right)\,dTdt
    !*FD \end{equation}

    implicit none

    real(sp), intent(in) :: artm      !*FD The annual mean air temperature (degC)

     pdd_integrand = artm *  exp(- (artm - t_a_prime)**2 / (2.0 * dd_sigma**2))

  end function pdd_integrand

!--------------------------------------------------------------------------------------

  recursive real(sp) function romberg_int(fct,lgr,rgr)

    !*FD Function to perform Romberg Integration on function \texttt{fct}, between
    !*FD limits \texttt{lgr} and \texttt{rgr}. The precision of the routine is 
    !*FD determined by the value of \texttt{ord}, an internal variable. 
    !*FD
    !*FD This routine is an implementation of ACM algorithm 60, by F. L. Bauer.
    !*FD (Comm. ACM, vol. 4, issue 6, June 1961).

    implicit none

    real(sp)            :: fct    !*FD Function to be integrated
    real(sp),intent(in) :: lgr    !*FD Lower bound
    real(sp),intent(in) :: rgr    !*FD Upper bound
    integer,parameter :: ord = 6

    real(sp),dimension(ord+1) :: t
    real(sp) :: l,u,m
    integer :: f,h,j,n

    external fct

    l=rgr-lgr
    t(1)=(fct(lgr)+fct(rgr))/2.0
    n=1

    do h=1,ord
       u=0
       m=l/(2*n)

       do j=1,2*n-1,2
          u=u+fct(lgr+j*m)
       end do

       t(h+1)=((u/n)+t(h))/2.0
       f=1
       
       do j=h,1,-1
          f=f*4
          t(j)=t(j+1)+(t(j+1)-t(j))/(f-1)
       end do

       n=2*n

    end do

    romberg_int=t(1)*l

  end function romberg_int

!-------------------------------------------------------------------------------

  integer function findgrid(rin,init,step)

    !*FD Calculates which row or column of the pdd table corresponds
    !*FD to a given value on the appropriate axis, so that:
    !*FD \[
    !*FD \mathtt{findgrid}=\frac{\mathtt{rin}-\mathtt{init}}{\mathtt{step}+1}
    !*FD \] 
    !*RV The relevant array index.

    use glimmer_global, only : sp
    
    implicit none
    
    real(sp), intent(in) :: rin  !*FD Value of axis variable at current point.
    real(sp), intent(in) :: init !*FD Value of axis variable at first point.
    real(sp), intent(in) :: step !*FD Grid spacing.
    
    findgrid = (rin - init) / step + 1

  end function findgrid
 
end module glint_pdd
