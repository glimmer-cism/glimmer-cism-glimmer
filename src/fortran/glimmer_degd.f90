
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_degd.f90 - part of the GLIMMER ice model         + 
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

module glimmer_degd

  !*FD Contains subroutines required to calculate the mass-balance
  !*FD of the ice-sheet. Formerly, this module included code specific to
  !*FD Antarctica, but this has been removed for the time being. The current
  !*FD code was originally specific to Greenland, but is in fact applied
  !*FD regardless of which region the ice model is covering.

  private pddtabgrn,qromb,qromb2,polint,trapzd
  private trapzd2,pdd1stint,pdd2ndint,findgrid

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

contains

!-------------------------------------------------------------------------------
! PUBLIC subroutine
!-------------------------------------------------------------------------------

  subroutine masbgrn(pddcalc,artm,arng,prcp,lati,ablt,acab)

    !*FD Calculates mass-balance over the ice model domain, by the
    !*FD positive-degree-day method.

    use glimmer_global, only : sp
    use paramets, only : pddfs, pddfi, wmax

    implicit none 
 
    type(glimmer_pddcalc),    intent(inout) :: pddcalc !*FD The positive-degree-day parameters
    real(sp), dimension(:,:), intent(in)    :: artm    !*FD Annual mean air-temperature 
                                                       !*FD ($^{\circ}$C)
    real(sp), dimension(:,:), intent(in)    :: arng    !*FD Annual temerature range ($^{\circ}$C)
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

    ! Check to see if pdd table is allocated. If not, then allocate it.

    if (.not.pddcalc%pt_alloc) then
      allocate(pddcalc%pddtab(pddcalc%nx,pddcalc%ny))
      pddcalc%pt_alloc=.true.
    endif

    ! If this is the first call to the subroutine, perform initialisataion
    ! of pdd table

    if (pddcalc%first) then 
      call pddtabgrn(pddcalc)  
      pddcalc%first = .false.
    end if

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

          wfrac = wmax * prcp(ew,ns)

          ! this is the total potential ablation of SNOW
          ! note we convert to scaled ablation
    
          pablt = pdd * pddfs

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
            ablt(ew,ns) = prcp(ew,ns) - wfrac + pddfi*(pdd-prcp(ew,ns)/pddfs) 
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

  end subroutine masbgrn                                  

!-------------------------------------------------------------------------------
! PRIVATE subroutines and functions
!-------------------------------------------------------------------------------

  subroutine pddtabgrn(pddcalc)

    !*FD Initialises the positive-degree-day-table.

    use glimmer_global, only: sp   
 
    implicit none

    type(glimmer_pddcalc),intent(inout) :: pddcalc !*FD PDD parameters

    ! Internal variables

    real(sp)           :: tma,dtmj
    real(sp),parameter :: twopi = 3.1416 * 2.0 
    integer  :: kx,ky
    
    !--------------------------------------------------------------------
    ! Main loops:
    !  tma  -- the mean annual temperature (y-axis of pdd table)
    !  dtmj -- difference from mean july temperature (x-axis of table)
    !  tmj -- the actual july temperature
    !--------------------------------------------------------------------

    do tma = pddcalc%iy, pddcalc%iy+(pddcalc%ny-1)*pddcalc%dy, pddcalc%dy

      ky = findgrid(tma,real(pddcalc%iy),real(pddcalc%dy))

      do dtmj = pddcalc%ix, pddcalc%ix+(pddcalc%nx-1)*pddcalc%dx, pddcalc%dx

        pddcalc%tmj = tma + dtmj   
        kx  = findgrid(dtmj,real(pddcalc%ix),real(pddcalc%dx)) 

        ! need these lines to take account of the backdoor message passing used here

        pddcalc%tma=tma
        pddcalc%dtmj=dtmj
     
        call qromb(pddcalc,0.0,twopi,pddcalc%pddtab(kx,ky))

        pddcalc%pddtab(kx,ky) = pddcalc%pddtab(kx,ky) / (pddcalc%dd_sigma * sqrt(twopi))
     
        ! convert to days     

        pddcalc%pddtab(kx,ky) = 365.0 * pddcalc%pddtab(kx,ky) / twopi

      end do
    end do

  end subroutine pddtabgrn

!-------------------------------------------------------------------------------

  subroutine qromb(pddcalc,a,b,ss)

    use glimmer_global, only: sp
      
    implicit none
    
    type(glimmer_pddcalc),intent(inout) :: pddcalc      
    real(sp),          intent(in)    :: a,b
    real(sp),          intent(inout) :: ss

    real(sp), parameter :: eps = 1.0e-6
    integer,  parameter :: jmax=20, jmaxp=jmax+1, k=5, km=k-1

    real(sp) :: dss
    integer j
    real(sp),dimension(jmaxp):: h,s
    
    h(1)=1.0
    
    do j = 1,jmax
    
      call trapzd(pddcalc,a,b,s(j),j)
     
      if (j>=k) then
        call polint(h(j-km),s(j-km),k,0.,ss,dss)
        if (abs(dss)<= eps*abs(ss)) return
      endif
     
      s(j+1)=s(j)
      h(j+1)=0.25*h(j)
     
    end do

    stop 'too many steps in qromb'

  end subroutine qromb

!-------------------------------------------------------------------------------
    
  subroutine qromb2(pddcalc,a,b,ss)
    
    use glimmer_global, only: sp
      
    implicit none
          
    type(glimmer_pddcalc),intent(inout) :: pddcalc      
    real(sp),          intent(in)    :: a,b
    real(sp),          intent(inout) :: ss

    real(sp), parameter :: eps=1.e-6
    integer,  parameter :: jmax=20, jmaxp=jmax+1, k=5, km=k-1

    real(sp) :: dss
    integer j
    real(sp), dimension (jmaxp):: h,s
    
    h(1)=1.0
    
    do j = 1,jmax
    
      call trapzd2(pddcalc,a,b,s(j),j)
     
      if (j>=k) then
        call polint(h(j-km),s(j-km),k,0.,ss,dss)
        if (abs(dss)<= eps*abs(ss)) return
      endif
     
      s(j+1)=s(j)
      h(j+1)=0.25*h(j)
     
    end do

    stop 'too many steps in qromb'

  end subroutine qromb2
    
!-------------------------------------------------------------------------------

  subroutine polint(xa,ya,n,x,y,dy)

     use glimmer_global, only: sp
      
    implicit none

    integer, intent(in) :: n
    integer:: i, m, ns
    integer, parameter :: nmax = 10
 
    real(sp), intent(in), dimension(n)::xa, ya     
    real(sp), intent(in):: x
    real(sp), intent(out):: dy,y
    real(sp) :: den, dif, dift, ho, hp, w
    real(sp), dimension(nmax):: c, d
    
    ns=1
    dif=abs(x-xa(1))
    
    do i=1, n
    
      dift=abs(x-xa(i))
     
      if (dift < dif) then
        ns=i
        dif=dift
      end if
     
      c(i)=ya(i)
      d(i)=ya(i)
     
    end do

    y=ya(ns)
    ns=ns-1
    
    do m=1,n-1
      do i=1,n-m

        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
    
        if (den == 0.) then
           write(*,*) 'failure in polint'
           stop
        end if

        den=w/den
        d(i)=hp*den
        c(i)=ho*den
      
      enddo
        
      if (2*ns < n-m)then
        dy=c(ns+1)
      else
        dy=d(ns)
        ns=ns-1
      endif
     
      y=y+dy

    enddo

  end subroutine polint

!-------------------------------------------------------------------------------
    
  subroutine trapzd(pddcalc,a,b,s,n)

    use glimmer_global, only: sp
      
    implicit none
          
    type(glimmer_pddcalc) :: pddcalc      
    integer, intent(in) :: n
    real(sp), intent(in) :: a, b
    real(sp), intent(inout) :: s

    integer :: it,j
    real(sp) :: del,sum,tnm,x
    
    if (n == 1) then
    
      s=0.5*(b-a)*(pdd1stint(pddcalc,a)+pdd1stint(pddcalc,b))
     
    else 
    
      it=2**(n-2)
      tnm=it
      del=(b-a)/tnm
      x=a+0.5*del
      sum=0.
    
      do j=1,it
        sum=sum+pdd1stint(pddcalc,x)
        x=x+del
      enddo

      s=0.5*(s+(b-a)*sum/tnm)
    
    endif
      
  end subroutine trapzd

!-------------------------------------------------------------------------------

  subroutine trapzd2(pddcalc,a,b,s,n)

    use glimmer_global, only : sp
      
    implicit none
          
    type(glimmer_pddcalc) :: pddcalc      
    integer, intent(in) :: n
    real(sp), intent(in) :: a, b
    real(sp), intent(inout) :: s
    integer :: it, j
    real  (sp) :: del, sum, tnm, x
    
    if (n == 1) then
    
      s=0.5*(b-a)*(pdd2ndint(pddcalc,a)+pdd2ndint(pddcalc,b))
     
    else 
    
      it=2**(n-2)
      tnm=it
      del=(b-a)/tnm
      x=a+0.5*del
      sum=0.
    
      do j=1,it
        sum=sum+pdd2ndint(pddcalc,x)
        x=x+del
      enddo

      s=0.5*(s+(b-a)*sum/tnm)
    
    endif

  end subroutine trapzd2

!-------------------------------------------------------------------------------

  real function pdd1stint(pddcalc,day)
    
    use glimmer_global, only : sp
      
    implicit none
    
    type(glimmer_pddcalc) :: pddcalc      
    real(sp), intent(in):: day 

    real(sp) :: upplim

    pddcalc%dailytemp = pddcalc%tma + (pddcalc%tmj - pddcalc%tma) * cos(day)  
    
    upplim = pddcalc%dailytemp + 2.5 * pddcalc%dd_sigma
    
    if ( upplim <= 0.0 ) then
      pdd1stint = 0.0
    else
      call qromb2(pddcalc,0.0,upplim,pdd1stint)
    end if

  end function pdd1stint

!-------------------------------------------------------------------------------
        
  real function pdd2ndint(pddcalc,artm)

    use glimmer_global, only : sp
         
    implicit none
    
    type(glimmer_pddcalc) :: pddcalc      
    real(sp), intent(in):: artm
    
    pdd2ndint = artm *  exp(- (artm - pddcalc%dailytemp)**2 / (2.0 * pddcalc%dd_sigma**2))

  end function pdd2ndint

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
  
end module glimmer_degd

