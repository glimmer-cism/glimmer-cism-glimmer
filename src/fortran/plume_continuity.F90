module plume_continuity

  use plume_global
  implicit none

contains

subroutine continuity(icalcan,kcalcan,icalcen,kcalcen)      

! evaluate entrainment and continuity equation to find interface position

     implicit none

! local variables

      integer i,k,l,icalcan,kcalcan,icalcen,kcalcen,mflag
      real(kind=kdp) pdepc(lm,ln),bspeed(lm,ln),speed(lm,ln)
      real(kind=kdp) rhopac,rhoa,delrho,rhoq,redg,tt,vmid,umid
      real(kind=kdp) rich,sm,arg,delta,iold,ites1,ites2
      real(kind=kdp) dragrt,prden,scden,gambt,gambs,tfb,tfi,c1,c2,c3
      real(kind=kdp) deltam(lm,ln),fppntot,ucrit,ucl


! 0. preliminaries
! ----------------
      prden = 12.5d0*pr**(2.d0/3.d0) - 9.d0
      scden = 12.5d0*sc**(2.d0/3.d0) - 9.d0

      do i = icalcan,icalcen
        do k = kcalcan,kcalcen
          deltam(i,k) = 0.d0
          pdepc(i,k) = bpos(i,k) - ipos(i,k)
          ! find flow speed on scalar grid
          tt = 5.0d-1*dy(k)*rdyv(k)
          vmid = tt*sv(i,k-1) + (1.-tt)*sv(i,k) 
          tt = 5.0d-1*dx(i-1)*rdxu(i)
          umid = tt*su(i,k) + (1.-tt)*su(i-1,k)
          speed(i,k) = umid**2 + vmid**2 + small
          bspeed(i,k) = dsqrt(speed(i,k))
        end do
      end do

! 1. entrainment
! --------------

      if (entrain) then
        do i = icalcan,icalcen
          do k = kcalcan,kcalcen
            if(pdepc(i,k).gt.edepth)then          

              rhopac = rhop(i,k)
              rhoa = rhoamb(i,k)
              delrho = rhoa - rhopac
              rhoq = 5.0d-1*(rhopac+rhoa)
              redg = delrho/rhoq
              rich = g*redg*pdepc(i,k)/speed(i,k)
              rich = dmax1(5.0d-2,rich)    


! full kochergin entrainment
              if (entype.eq.1) then
                sm = rich/(7.25d-1*(rich + 1.86d-1 - dsqrt(rich*rich- &
     &             3.16d-1*rich + 3.46d-2)))
                arg = dmax1(small,speed(i,k) + g*redg*pdepc(i,k)/sm)
                entr(i,k) = cl**2*dsqrt(arg)/sm                
              end if

! reduced kochergin entrainment matched to pedersen
              if (entype.eq.2) then
                sm = rich/(7.25d-1*(rich + 1.86d-1 - dsqrt(rich*rich- &
     &             3.16d-1*rich + 3.46d-2)))
                arg = dmax1(small,speed(i,k) + g*redg*pdepc(i,k)/sm)
                entr(i,k) = ef*cl**2*(drag(i,k)/3.0d-3)*dsqrt(arg)/sm                
              end if

! half pedersen entrainment
              if (entype.eq.3) then
                entr(i,k) = 3.6d-2*bspeed(i,k)*drag(i,k)/rich
              end if

! half modified pedersen entrainment
              if (entype.eq.4) then
                entr(i,k) = 3.6d-2*bspeed(i,k)*dsin(1.0d-3)
              end if

              deltam(i,k) = deltam(i,k) + entr(i,k)

            end if
          end do
        end do
      end if

! 2. basal melting
! ----------------

      if (basmelt) then
        do i = icalcan,icalcen
          do k = kcalcan,kcalcen
            if (pdepc(i,k).gt.mdepth) then         
              dragrt = dsqrt(drag(i,k))

! find turbulent exchange coefficients
              gambt = (dragrt*bspeed(i,k))/ &
     &       (2.12d0*dlog((dragrt*bspeed(i,k)*pdepc(i,k))/nu0) + prden)
              gambs = (dragrt*bspeed(i,k))/ &
     &       (2.12d0*dlog((dragrt*bspeed(i,k)*pdepc(i,k))/nu0) + scden)

! calculate freezing point of plume at shelf base, decide if melt (mflag = 1) 
! or freeze and calculate freezing point of ice at shelf base
              tfb = fta*salt(i,k) + ftb + ftc*(gldep + wcdep -bpos(i,k))
              mflag = (1 + int(sign(1.d0,temp(i,k) - tfb)))/2
              tfi = (1 - mflag)*fta*si  &
     &                           + ftb + ftc*(gldep + wcdep - bpos(i,k))

! calculate coefficients in quadratic to be solved
              c1 = lat/c0 + mflag*(ci/c0)*(tfi-ti)
              c2 = gambs*(lat/c0 + mflag*(ci/c0)*(tfb-ti))  &
     &                                           + gambt*(tfi-temp(i,k))
              c3 = gambs*gambt*(tfb - temp(i,k))

! calculate melt rate
              bmelt(i,k) = -(c2 - dsqrt(c2*c2 - 4.d0*c1*c3))/(2.d0*c1)

!! multiply by 10 if freezing
!              if (mflag.eq.0) then
!                bmelt(i,k) = bmelt(i,k)*10.d0
!              end if

              deltam(i,k) = deltam(i,k) + bmelt(i,k)

            end if
          end do
        end do
      end if

! 3.frazil precipitation
! ----------------------

      if (frazil) then
        do i = icalcan,icalcen
          do k = kcalcan,kcalcen
            if (pdepc(i,k).gt.fdepth) then         

              fppntot = 0.d0
              do l=1,nice
                if (c(i,k,l).gt.0.d0) then

! calculate precipitation of frazil
                  ucl = &
     &              dsqrt((1.0d-1*g*re(l)*(rho0-rhoi))/(rho0*drag(i,k)))
                  ucrit = (1.d0 - bspeed(i,k)*bspeed(i,k)/(ucl*ucl))
                  fppn(i,k,l) = &
     &                     dmin1(-(rhoi/rho0)*c(i,k,l)*wi(l)*ucrit,0.d0)

                  fppntot = fppntot + fppn(i,k,l)
                end if
              end do

              deltam(i,k) = deltam(i,k) + fppntot

            end if
          end do
        end do
      end if

! final divergence plus entrainment and melting
! ---------------------------------------------
!    
      do i = icalcan,icalcen
        do k = kcalcan,kcalcen
          delta = dt*(rdxu(i)*(u(i-1,k)-u(i,k)) &
                   + rdyv(k)*(v(i,k-1)-v(i,k))) &
                   + deltam(i,k)*dt

! update interface position
          iold = ipos(i,k)
          ipos(i,k) = ipos(i,k) - delta

! check for negative predicted depth of ambient fluid 
          jcd_negdep(i,k) = 0
          ites1 = iold - 2.d0*delta
          ites1 = bpos(i,k) - ites1
          if (ites1 .lt. 0.d0) jcd_negdep(i,k) = 2
          ites2 = iold - 3.d0*delta
          ites2 = bpos(i,k) - ites2 
          if (ites2 .lt. 0.d0) jcd_negdep(i,k) = 1

        end do
      end do

      return

end subroutine continuity

end module
