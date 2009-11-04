module plume_scalar

contains

subroutine scalar(icalcan,kcalcan,icalcen,kcalcen,frzcut,seedtype)             

! calculates transport and entrainment of all scalars, evolution of 
! frazil ice, basal melting and freezing, and inflow tracking

  use plume_global
  use plume_functions

  implicit none

! local variables

      integer i,k,l,icalcan,kcalcan,icalcen,kcalcen,seedtype,izo,izu
      integer idel,kdel,idx,kdy,ihilf,khilf,mflag,seedindex

      real(kind=kdp) deltat(lm,ln),deltas(lm,ln)
      real(kind=kdp) deltac(lm,ln,lice),deltatr(lm,ln,linf)
      real(kind=kdp) pdepc(lm,ln),pdepcp(lm,ln),vmid(lm,ln)
      real(kind=kdp) umid(lm,ln),speed(lm,ln),bspeed(lm,ln)
      real(kind=kdp) depth(lm,ln)
      real(kind=kdp) one,slon,sloe,slos,slow,sumslo
      real(kind=kdp) pressure,ttt,difu,difo,rhoa,redg
      real(kind=kdp) tt,rich,sm,arg
      real(kind=kdp) dxx,dyy,sx,sy,sxy,r1,r2
      real(kind=kdp) dife,difw,difn,difs
      real(kind=kdp) dragrt,prden,scden,gambt,gambs,tfb,tfi,c1,c2,c3
      real(kind=kdp) gamct,gamcs,ucrit,ucl
      real(kind=kdp) fmelttot,wturb,mfac1,mfac2,gi,gim1,mi,mip1
      real(kind=kdp) frzcut(lice)

      one = 1.d0
      prden = 12.5d0*pr**(2.d0/3.d0) - 9.d0
      scden = 12.5d0*sc**(2.d0/3.d0) - 9.d0

! calculate values to be used for newly-wet cells using weights 
! of interface gradient (thus tempa etc is temp to be used for whole plume)

      do i = icalcan,icalcen
        do k = kcalcan,kcalcen
          if (jcd_fl(i,k) .eq. 1) then
            slon = dmax1(0.d0,ipos(i,k)-ipos(i,k+1))*jcd(i,k+1)
            sloe = dmax1(0.d0,ipos(i,k)-ipos(i+1,k))*jcd(i+1,k)
            slos = dmax1(0.d0,ipos(i,k)-ipos(i,k-1))*jcd(i,k-1)
            slow = dmax1(0.d0,ipos(i,k)-ipos(i-1,k))*jcd(i-1,k)
            sumslo = slon + sloe + slos + slow        

            if (sumslo .gt. dcr) then
              tempa(i,k) = (temp(i,k+1)*slon + temp(i+1,k)*sloe &
                   + temp(i,k-1)*slos + temp(i-1,k)*slow)/sumslo
              salta(i,k) = (salt(i,k+1)*slon + salt(i+1,k)*sloe &
                   + salt(i,k-1)*slos + salt(i-1,k)*slow)/sumslo

              if (rholinear) then
                rhop(i,k) = rho_func_linear(tempa(i,k),salta(i,k))
              else
                if (thermobar) then
                  pressure = 1.0d-1*(wcdep + gldep - ipos(i,k))
                  ttt = tinsitu_func(tempa(i,k),salta(i,k),pressure)
                  rhop(i,k) =rho_func_nonlinear(ttt,salta(i,k),pressure)
                else
                  rhop(i,k) =  &
                         rho_func_nonlinear(tempa(i,k),salta(i,k),0.d0)
                end if
              end if

              if (frazil) then
                ctota(i,k) = 0.d0
                do l = 1,nice
                  ca(i,k,l) = (c(i,k+1,l)*slon + c(i+1,k,l)*sloe &
                   + c(i,k-1,l)*slos + c(i-1,k,l)*slow)/sumslo
                  ctota(i,k) = ctota(i,k) + ca(i,k,l)
                end do
                rhop(i,k) = (1.d0-ctota(i,k))*rhop(i,k)+ctota(i,k)*rhoi
              end if

              if (intrace) then
                do l = ninfmin,ninfmax
                  intra(i,k,l) =  &
                    (intr(i,k+1,l)*slon + intr(i+1,k,l)*sloe  &
                   + intr(i,k-1,l)*slos + intr(i-1,k,l)*slow)/sumslo
                end do
              end if

            end if
          end if

        end do
      end do

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 0)preliminaries
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      do i = icalcan,icalcen
        do k = kcalcan,kcalcen
          deltat(i,k) = 0.d0
          deltas(i,k) = 0.d0

          pdepc(i,k) = bpos(i,k) - ipos(i,k)
          pdepcp(i,k) = pdepc(i,k) + dcr
          tt = 5.0d-1*dx(i)*rdxu(i)
          umid(i,k) = tt*su(i-1,k) + (1.d0-tt)*su(i,k)
          tt = 5.0d-1*dy(k)*rdyv(k)
          vmid(i,k) = tt*sv(i,k-1) + (1.d0-tt)*sv(i,k) 
          speed(i,k) = umid(i,k)**2 + vmid(i,k)**2 + small
          bspeed(i,k) = dsqrt(speed(i,k))

! calculate freezing temperature at plume mid-depth (even if no frazil)
          depth(i,k) = gldep + wcdep - 5.0d-1*(ipos(i,k) + bpos(i,k))
          tf(i,k) = freezing_temp_func(salta(i,k),depth(i,k))
        end do
      end do

      if (frazil) then
        do i = icalcan,icalcen
          do k = kcalcan,kcalcen
            do l = 1,nice
              deltac(i,k,l) = 0.d0
            end do
          end do
        end do
      end if

      if (intrace) then
        do i = icalcan,icalcen
          do k = kcalcan,kcalcen
            do l = ninfmin,ninfmax
              deltatr(i,k,l) = 0.d0
            end do
          end do
        end do
      end if

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 1)entrainment
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if (entrain) then
        do i = icalcan,icalcen
          do k = kcalcan,kcalcen
            if (pdepc(i,k).gt.edepth) then         

! interpolate ambient values to interface and calculate
! quantities for entrainment
              izo = int((wcdep + gldep - ipos(i,k))/dzincr) + 1
              izu = izo + 1
              difu = dble(izo)*dzincr - (wcdep + gldep - ipos(i,k))
              difo = dzincr - difu
              atemp(i,k) = (difu*tamb(izo) + difo*tamb(izu))/dzincr
              asalt(i,k) = (difu*samb(izo) + difo*samb(izu))/dzincr
              rhoa = (difu*rhovf(izo) + difo*rhovf(izu))/dzincr
              redg = (rhoa - rhop(i,k))/(5.0d-1*(rhoa + rhop(i,k)))
              rich = g*redg*pdepc(i,k)/speed(i,k)
              rich = dmax1(5.0d-02,rich)    

! full kochergin entrainment
              if (entype.eq.1) then
                sm = rich/(7.25d-1*(rich + 1.86d-1 - dsqrt(rich*rich- &
                3.16d-1*rich + 3.46d-2)))
                arg = dmax1(small,speed(i,k) + g*redg*pdepc(i,k)/sm)
                entr(i,k) = cl**2*dsqrt(arg)/sm                
              end if

! reduced kochergin entrainment matched to pedersen
              if (entype.eq.2) then
                sm = rich/(7.25d-1*(rich + 1.86d-1 - dsqrt(rich*rich- &
                3.16d-1*rich + 3.46d-2)))
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

              deltat(i,k) = deltat(i,k)  &
                   + dt*entr(i,k)*(atemp(i,k) - tempa(i,k))/pdepcp(i,k)
              deltas(i,k) = deltas(i,k)  &
                   + dt*entr(i,k)*(asalt(i,k) - salta(i,k))/pdepcp(i,k)

            endif
          end do
        end do

        if (frazil) then
          do i = icalcan,icalcen
            do k = kcalcan,kcalcen
              if (pdepc(i,k).gt.edepth) then         
                do l=1,nice
                  deltac(i,k,l) = deltac(i,k,l)  &
                                + dt*entr(i,k)*(-ca(i,k,l))/pdepcp(i,k)
                end do
              end if
            end do
          end do
        end if

      end if

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 2)advection (selective vector upstream)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if (nonlin) then
        do i = icalcan,icalcen
          do k = kcalcan,kcalcen

!  choose upstream quarter
            idel = -int(sign(one,umid(i,k)))
            kdel = -int(sign(one,vmid(i,k)))
            ihilf = i + idel
            khilf = k + kdel
!  choose the relevant upstream cell dimensions 
            idx = i + (idel-1)/2
            kdy = k + (kdel-1)/2
            dxx = dx(idx)
            dyy = dy(kdy)
!  weighting coefficients
            sx = abs(umid(i,k))*dt/dxx
            sy = abs(vmid(i,k))*dt/dyy
            sxy = sx - sy
            r1 = (dsign(one,sxy) + 1.)*5.0d-1
            r2 = one - r1

            if ((sx.gt.one).or.(sy.gt.one)) then
              write(*,*) 'error: tsc courant exceeded at i=',i,' k=',k
              write(11,*) 'error: tsc courant exceeded at i=',i,' k=',k
            end if

            deltat(i,k) = deltat(i,k)  &
                   + (r1*sy+r2*sx)*(tempa(ihilf,khilf)*jcd(ihilf,khilf) &
                   + (1 - jcd(ihilf,khilf))*tempa(i,k))  &
                   + r1*sxy*(tempa(ihilf,k)*jcd(ihilf,k) &
                   + (1-jcd(ihilf,k))*tempa(i,k)) &
                   - r2*sxy*(tempa(i,khilf)*jcd(i,khilf) &
                   + (1-jcd(i,khilf))*tempa(i,k)) &
                   - (r2*sy+r1*sx)*tempa(i,k)                         

            deltas(i,k) = deltas(i,k)  &
                   + (r1*sy+r2*sx)*(salta(ihilf,khilf)*jcd(ihilf,khilf) &
                   + (1 - jcd(ihilf,khilf))*salta(i,k))  &
                   + r1*sxy*(salta(ihilf,k)*jcd(ihilf,k) &
                   + (1-jcd(ihilf,k))*salta(i,k)) &
                   - r2*sxy*(salta(i,khilf)*jcd(i,khilf) &
                   + (1-jcd(i,khilf))*salta(i,k)) &
                   - (r2*sy+r1*sx)*salta(i,k)                         

            if (frazil) then
              do l = 1,nice
                deltac(i,k,l) = deltac(i,k,l) &
                    + (r1*sy+r2*sx)*(ca(ihilf,khilf,l)*jcd(ihilf,khilf) &
                    + (1 - jcd(ihilf,khilf))*ca(i,k,l))  &
                    + r1*sxy*(ca(ihilf,k,l)*jcd(ihilf,k)  &
                    + (1-jcd(ihilf,k))*ca(i,k,l)) &
                    - r2*sxy*(ca(i,khilf,l)*jcd(i,khilf) &
                    + (1-jcd(i,khilf))*ca(i,k,l)) &
                    - (r2*sy+r1*sx)*ca(i,k,l)      
              end do
            end if

            if (intrace) then
              do l = ninfmin,ninfmax
                deltatr(i,k,l) = deltatr(i,k,l) &
                 + (r1*sy+r2*sx)*(intra(ihilf,khilf,l)*jcd(ihilf,khilf) &
                    + (1 - jcd(ihilf,khilf))*intra(i,k,l)) &
                    + r1*sxy*(intra(ihilf,k,l)*jcd(ihilf,k)  &
                    + (1-jcd(ihilf,k))*intra(i,k,l)) &
                    - r2*sxy*(intra(i,khilf,l)*jcd(i,khilf) &
                    + (1-jcd(i,khilf))*intra(i,k,l)) &
                    - (r2*sy+r1*sx)*intra(i,k,l)      
              end do
            end if

          end do
        end do

      end if

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 3)horizontal diffusion
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if (horturb) then

        do i = icalcan,icalcen
          do k = kcalcan,kcalcen

            dife = (tempa(i+1,k) - tempa(i,k))*jcd(i+1,k)*kh*rdx(i)
            difw = (tempa(i,k) - tempa(i-1,k))*jcd(i-1,k)*kh*rdx(i-1)
            difn = (tempa(i,k+1) - tempa(i,k))*jcd(i,k+1)*kh*rdy(k)
            difs = (tempa(i,k) - tempa(i,k-1))*jcd(i,k-1)*kh*rdy(k-1)
            deltat(i,k) = deltat(i,k) &
                     + ((dife - difw)*rdxu(i) + (difn-difs)*rdyv(k))*dt 

            dife = (salta(i+1,k) - salta(i,k))*jcd(i+1,k)*kh*rdx(i)
            difw = (salta(i,k) - salta(i-1,k))*jcd(i-1,k)*kh*rdx(i-1)
            difn = (salta(i,k+1) - salta(i,k))*jcd(i,k+1)*kh*rdy(k)
            difs = (salta(i,k) - salta(i,k-1))*jcd(i,k-1)*kh*rdy(k-1)
            deltas(i,k) = deltas(i,k)  &
                     + ((dife - difw)*rdxu(i) + (difn-difs)*rdyv(k))*dt

          end do 
        end do

        if (frazil) then
          do i = icalcan,icalcen
            do k = kcalcan,kcalcen
              do l = 1,nice
                dife = (ca(i+1,k,l) - ca(i,k,l))*jcd(i+1,k)*kh*rdx(i)
                difw = (ca(i,k,l) - ca(i-1,k,l))*jcd(i-1,k)*kh*rdx(i-1)
                difn = (ca(i,k+1,l) - ca(i,k,l))*jcd(i,k+1)*kh*rdy(k)
                difs = (ca(i,k,l) - ca(i,k-1,l))*jcd(i,k-1)*kh*rdy(k-1)
                deltac(i,k,l) = deltac(i,k,l)  &
                     + ((dife - difw)*rdxu(i) + (difn-difs)*rdyv(k))*dt
              end do
            end do
          end do
        end if

        if (intrace) then
          do i = icalcan,icalcen
            do k = kcalcan,kcalcen
              do l = ninfmin,ninfmax
                dife =  &
                   (intra(i+1,k,l) - intra(i,k,l))*jcd(i+1,k)*kh*rdx(i)
                difw =  &
                 (intra(i,k,l) - intra(i-1,k,l))*jcd(i-1,k)*kh*rdx(i-1)
                difn =  &
                  (intra(i,k+1,l) - intra(i,k,l))*jcd(i,k+1)*kh*rdy(k)
                difs =  &
                 (intra(i,k,l) - intra(i,k-1,l))*jcd(i,k-1)*kh*rdy(k-1)
                deltatr(i,k,l) = deltatr(i,k,l)  &
                     + ((dife - difw)*rdxu(i) + (difn-difs)*rdyv(k))*dt
              end do
            end do
          end do
        end if

      end if

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 4)basal melting and freezing
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if (basmelt) then
        do i = icalcan,icalcen
          do k = kcalcan,kcalcen
            if (pdepc(i,k).gt.mdepth) then         
              dragrt = dsqrt(drag(i,k))

! find turbulent exchange coefficients
              gambt = (dragrt*bspeed(i,k))/ &
            (2.12d0*dlog((dragrt*bspeed(i,k)*pdepc(i,k))/nu0) + prden)
              gambs = (dragrt*bspeed(i,k))/ &
            (2.12d0*dlog((dragrt*bspeed(i,k)*pdepc(i,k))/nu0) + scden)

! calculate freezing point of plume at shelf base, decide if melt (mflag = 1) 
! or freeze and calculate freezing point of ice at shelf base.  note that heat 
! conduction into the ice shelf is taken into account here if melting (as removal 
! of warmest ice steepens gradient) but is not represented in heat source terms
              tfb = fta*salta(i,k) + ftb + ftc*(gldep+wcdep-bpos(i,k))
              mflag = (1 + int(sign(1.d0,tempa(i,k) - tfb)))/2
              tfi = (1 - mflag)*fta*si  &
                                + ftb + ftc*(gldep + wcdep - bpos(i,k))

! calculate coefficients in quadratic to be solved
              c1 = lat/c0 + mflag*(ci/c0)*(tfi-ti)
              c2 = gambs*(lat/c0 + mflag*(ci/c0)*(tfb-ti))  &
                                               + gambt*(tfi-tempa(i,k))
              c3 = gambs*gambt*(tfb - tempa(i,k))

! calculate melt rate
              bmelt(i,k) = -(c2 - dsqrt(c2*c2 - 4.d0*c1*c3))/(2.d0*c1)

!! multiply by 10 if freezing
!              if (mflag.eq.0) then
!                bmelt(i,k) = bmelt(i,k)*10.d0
!              end if

! calculate basal temperature and salinity
              btemp(i,k) = (gambt*tempa(i,k)+mflag*(ci/c0)*bmelt(i,k)*ti &
                        - (lat/c0)*bmelt(i,k) )/  &
                        (gambt + mflag*(ci/c0)*bmelt(i,k))
              bsalt(i,k) = (btemp(i,k) - ftb  &
                                 - ftc*(gldep + wcdep - bpos(i,k)))/fta

! calculate change of heat, salt, and frazil due to meltwater/freezewater flux
! (ice shelf has zero salinity and no frazil)
              deltat(i,k) = deltat(i,k)  &
                  + dt*bmelt(i,k)*(btemp(i,k) - tempa(i,k))/pdepcp(i,k)
              deltas(i,k) = deltas(i,k)  &
                                 - dt*bmelt(i,k)*salta(i,k)/pdepcp(i,k)

! add diffusive flux of latent heat release/uptake
              deltat(i,k) = deltat(i,k)  &
                                   - dt*(lat/c0)*bmelt(i,k)/pdepcp(i,k)

            end if
          end do
        end do

        if (frazil) then
          do i = icalcan,icalcen
            do k = kcalcan,kcalcen
              if (pdepc(i,k).gt.mdepth) then         
                do l=1,nice
                  deltac(i,k,l) = deltac(i,k,l)  &
                                  - dt*bmelt(i,k)*ca(i,k,l)/pdepcp(i,k)
                end do
              end if
            end do
          end do
        end if

      end if

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 5)frazil nucleation, growth, melting, and precipitation
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if (frazil) then
        do i = icalcan,icalcen
          do k = kcalcan,kcalcen
            if (pdepc(i,k).gt.fdepth) then         
!              
              fmelttot = 0.d0

              if (nice.eq.1) then
! ----------------------------
! jenkins and bombosch version
! ----------------------------
                if (ca(i,k,1).ge.cmin(1)) then

! a) calculate growth/melting of frazil 
                  gamct = (nuss(1)*kt)/(ar*r(1))
                  gamcs = (nuss(1)*ks)/(ar*r(1))

                  c1 = 5.0d-1*lat*r(1)/(c0*ca(i,k,1)*pdepc(i,k)*gamct)
                  c2 = ftb + ftc*depth(i,k) - tempa(i,k)  &
                                                 + lat*gamcs/(c0*gamct)
                  c3 = (2.0d0*ctota(i,k)*pdepc(i,k)*gamcs/r(1)) &
                   *(fta*salta(i,k)+ ftb + ftc*depth(i,k) - tempa(i,k))

                  fmelt(i,k,1) = &
                            -(c2 - dsqrt(c2*c2 - 4.d0*c1*c3))/(2.d0*c1)
                  ctempd(i,k) = tempa(i,k) - fmelt(i,k,1)* &
                     (5.0d-1*lat*r(1)/(c0*ctota(i,k)*pdepc(i,k)*gamct))

! b) calculate precipitation of frazil
                  ucl =  &
                   dsqrt((1.0d-1*g*re(1)*(rho0-rhoi))/(rho0*drag(i,k)))
                  ucrit = (1.d0 - bspeed(i,k)*bspeed(i,k)/(ucl*ucl))
                  fppn(i,k,1) =  &
                         dmin1(-(rhoi/rho0)*ca(i,k,1)*wi(1)*ucrit,0.d0)

! c) calculate overall frazil source for this size class and 
!    total melt for ts source
                  deltac(i,k,1) = deltac(i,k,1)  &
                                 - dt*fppn(i,k,1)*ca(i,k,1)/pdepcp(i,k) &
              + (rho0/rhoi)*dt*(fppn(i,k,1) - fmelt(i,k,1))/pdepcp(i,k)
                  fmelttot = fmelt(i,k,1)
                end if

              else
! ----------------------------
! smedsrud and jenkins version
! ----------------------------
                if (ctota(i,k).ge.cmin(1)) then

! a) calculate secondary nucleation of frazil

                  fnuc(i,k,1) = 0.d0
                  do l=2,nice
                    wturb = dsqrt(wi(l)*wi(l)  &
                                + ((4.d0*eps)/(15.d0*nu0))*re(l)*re(l))
                    fnuc(i,k,l) =  &
                           - (rhoi/rho0)*pdepc(i,k)*pi*dble(nbar)*wturb &
                                            *(re(1)**3/re(l))*ca(i,k,l)
                    fnuc(i,k,1) = fnuc(i,k,1) - fnuc(i,k,l)
                  end do

! b) calculate growth/melting of frazil 
                  mfac1 = (2.d0*c0*kt/lat)*(tf(i,k) - tempa(i,k))
                  mfac2 = (rhoi/rho0)*pdepc(i,k)

                  if (tempa(i,k).lt.tf(i,k)) then
                    gi = mfac1*nuss(1)*ca(i,k,1)/(r(1)*r(1))
                    fmelt(i,k,1) = mfac2*(vol(1)/(vol(2) - vol(1)))*gi

                    do l=2,nice-1
                      gi = mfac1*nuss(l)*ca(i,k,l)/(r(l)*r(l))
                      gim1 = mfac1*nuss(l-1)*ca(i,k,l-1)/(r(l-1)*r(l-1))
                      fmelt(i,k,l) = mfac2*( &
                           (vol(l)/(vol(l+1) - vol(l)))*gi  &           
                                   - (vol(l)/(vol(l) - vol(l-1)))*gim1)             
                    end do

                    gim1 = mfac1*nuss(nice-1)*ca(i,k,nice-1) &
                                                 /(r(nice-1)*r(nice-1))
                    fmelt(i,k,nice) = - mfac2* &
                             (vol(nice)/(vol(nice) - vol(nice-1)))*gim1
                  else
                    mi = mfac1*nuss(1)*ca(i,k,1) &
                                        *(1.d0/r(1) + 1.d0/thi(1))/r(1)
                    mip1 = mfac1*nuss(2)*ca(i,k,2) &
                                        *(1.d0/r(2) + 1.d0/thi(2))/r(2)
                    fmelt(i,k,1) = mfac2*  &
                                 ((vol(1)/(vol(2) - vol(1)))*mip1 - mi)

                    do l=2,nice-1
                      mi = mfac1*nuss(l)*ca(i,k,l) &
                                        *(1.d0/r(l) + 1.d0/thi(l))/r(l)
                      mip1 = mfac1*nuss(l+1)*ca(i,k,l+1)* &
                                   (1.d0/r(l+1) + 1.d0/thi(l+1))/r(l+1)
                      fmelt(i,k,l) = mfac2*( &
                           (vol(l)/(vol(l+1) - vol(l)))*mip1      &       
                                     - (vol(l)/(vol(l) - vol(l-1)))*mi)  
                    end do

                    mi = mfac1*nuss(nice)*ca(i,k,nice)* &
                                (1.d0/r(nice) + 1.d0/thi(nice))/r(nice)
                    fmelt(i,k,nice) = - mfac2* &
                               (vol(nice)/(vol(nice) - vol(nice-1)))*mi
                  end if

! c) calculate precipitation of frazil
                  do l=1,nice
                    ucl =  &
                   dsqrt((1.0d-1*g*re(l)*(rho0-rhoi))/(rho0*drag(i,k)))
                    ucrit = (1.d0-bspeed(i,k)*bspeed(i,k)/(ucl*ucl))
                    fppn(i,k,l) =  &
                         dmin1(-(rhoi/rho0)*ca(i,k,l)*wi(l)*ucrit,0.d0)
                  end do

! d) calculate overall frazil source for this size class and 
!    total melt for ts source
                  do l=1,nice
                    deltac(i,k,l) = deltac(i,k,l) &
                                 - dt*fppn(i,k,l)*ca(i,k,l)/pdepcp(i,k) &
                          + (rho0/rhoi)*dt*(fppn(i,k,l) + fnuc(i,k,l)  &
                                            - fmelt(i,k,l))/pdepcp(i,k)
                    fmelttot = fmelttot + fmelt(i,k,l)
                  end do
                end if

              end if

! -------------
! both versions
! -------------
! calculate change of heat and salt due to frazil meltwater/freezewater flux
              deltat(i,k) = deltat(i,k)  &
                       + dt*fmelttot*(tf(i,k) - tempa(i,k))/pdepcp(i,k)
              deltas(i,k) = deltas(i,k)  &
                                   - dt*fmelttot*salta(i,k)/pdepcp(i,k)

! add diffusive flux of latent heat release/uptake
              deltat(i,k) = deltat(i,k)-dt*(lat/c0)*fmelttot/pdepcp(i,k)

            end if
          end do
        end do

      end if

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! final
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      do i = icalcan,icalcen
        do k = kcalcan,kcalcen
          temp(i,k) = tempa(i,k) + deltat(i,k)
          salt(i,k) = salta(i,k) + deltas(i,k)
        end do
      end do

      if (frazil) then

! if desired, check pre-emptively for negative frazil concentrations 
! and set to zero if found.  keep track of total negative concentration
        if (negfrz) then

          do i = icalcan,icalcen
            do k = kcalcan,kcalcen
              do l = 1,nice
                if (ca(i,k,l) + deltac(i,k,l).ge.0.d0) then
                  c(i,k,l) = ca(i,k,l) + deltac(i,k,l)
                else
                  c(i,k,l) = 0.d0
                  frzcut(l) = frzcut(l) - ca(i,k,l) - deltac(i,k,l)
                end if
              end do
            end do
          end do

        else

          do i = icalcan,icalcen
            do k = kcalcan,kcalcen
              do l = 1,nice
                c(i,k,l) = ca(i,k,l) + deltac(i,k,l)
              end do
            end do
          end do

        end if

! total up frazil concentrations
        do i = icalcan,icalcen
          do k = kcalcan,kcalcen
            ctot(i,k) = 0.d0
            do l = 1,nice
              ctot(i,k) = ctot(i,k) + c(i,k,l) 
            end do
          end do
        end do

      end if

      if (intrace) then
        do i = icalcan,icalcen
          do k = kcalcan,kcalcen
            do l = ninfmin,ninfmax
              intr(i,k,l) = intra(i,k,l) + deltatr(i,k,l)
            end do
          end do
        end do
      end if

! seed frazil ice
! ---------------
      if (frazil) then 
!           
! smedsrud and jenkins seding strategy:
! (find and then seed southernmost supercooled cells          
!           
        if (seedtype.eq.1) then
!     
          seedindex = n
          do i = icalcan,icalcen
            do k = kcalcan,kcalcen
              if (pdepc(i,k).gt.fdepth) then
                if (tempa(i,k).lt.tf(i,k)) then
                  seedindex = min(k,seedindex)
                end if
              end if
            end do
          end do

          if (seedindex.lt.n) then
            do i = icalcan,icalcen
              if (pdepc(i,seedindex).gt.fdepth) then
                ctot(i,seedindex) = 0.d0
                do l=1,nice
                  c(i,seedindex,l) = cseed(l)
                  ctot(i,seedindex) = ctota(i,seedindex) + cseed(l)
                end do
              end if
            end do
          end if

        end if

! holland and feltham seding strategy:
! (if cell is newly supercooled (i.e. was superheated on last timestep), seed
! with frazil if advection and diffusion have not exceeded seed population)
!           
        if (seedtype.eq.2) then

          do i = icalcan,icalcen
            do k = kcalcan,kcalcen
              if (pdepc(i,k).gt.fdepth) then         
                if (tempa(i,k).gt.tf(i,k)) then
                  jcd_fseed(i,k) = 0
                else
                  if (jcd_fseed(i,k).eq.0) then
                    ctot(i,k) = 0.d0
                    do l = 1,nice
                      if (c(i,k,l).lt.cseed(l)) c(i,k,l) = cseed(l)
                      ctot(i,k) = ctot(i,k) + c(i,k,l)
                    end do
                    jcd_fseed(i,k) = 1     
                  end if
                end if
              end if
            end do
          end do
!           
        end if
!           
      end if
!           
! tidy up variables
! -----------------
      do i = icalcan,icalcen
        do k = kcalcan,kcalcen
          tempa(i,k) = temp(i,k)
          salta(i,k) = salt(i,k)
        end do
      end do
!     
      if (frazil) then
        do i = icalcan,icalcen
          do k = kcalcan,kcalcen
            ctota(i,k) = ctot(i,k)
            do l = 1,nice
              ca(i,k,l) = c(i,k,l)
            end do
          end do
        end do
      end if
!     
      if (intrace) then
        do i = icalcan,icalcen
          do k = kcalcan,kcalcen
            do l = ninfmin,ninfmax
              intra(i,k,l) = intr(i,k,l)
            end do
          end do
        end do
      end if

      return
end subroutine

end module
