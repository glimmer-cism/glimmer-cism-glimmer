module plume_momentum

contains

subroutine momentum(icalcan,kcalcan,icalcen,kcalcen, &
                                                   cosfdt,sinfdt,f)

! calculates velocity components within plume

! all lower-case code is turning angle theory by alex wilchinsky
! (wilchinsky, feltham and holland, submitted to j. phys. oceanogr.)
! it is probably wise to check or reprogram this before use

  use plume_global
  use plume_functions	
  implicit none

! local variables

      integer i,k,icalcan,kcalcan,icalcen,kcalcen,jcvfac
      integer izo,izu,idel,kdel,iidx,kkdy,ihilf,khilf

      real(kind=kdp) cosfdt,sinfdt,one,termnl,termnl2,corx
      real(kind=kdp) redgu,slorho,islope
      real(kind=kdp) pdepu,zu,zum,arfac,tt,uu,umid,vmid
      real(kind=kdp) speed,tbotm,difu,difo,rhoa,tu,salu
      real(kind=kdp) rhoc,rhoe,rhou,rhoq,dxx,dyy,sx,sy,sxy,r1,r2
      real(kind=kdp) tlate,tlatw,tlats,tlatn,hordif,cory,redgv
      real(kind=kdp) pdepv,zv,zvm,tv,salv,rhon,rhov
      real(kind=kdp) delta,itest,ctotu,ctotv,dragu,dragv

      logical olddrag,norotation,variableekman,draginmelt 
      logical newudrag(lm,ln)
      real(kind=kdp) av,ekthick,kq,costang,sintang,thickratio,f
      real(kind=kdp) ugriddrag(lm,ln)

      sintang = 0.0
      kq = 0.0

! set switches for turning angle model

      olddrag       = .false. ! use the old drag magnitude
      norotation    = .false. ! switch off drag rotation
      variableekman = .false. ! use variable vertical eddy visc (hence ek thick)
      draginmelt    = .true.  ! write drag array so that drag used in melting
                              ! (only works with spatially constant drag coeff)

      one = 1.d0

! calculate ekman layer thickness for constant vertical eddy diffusivity

      if (.not.variableekman) then
        av = 5.d-4
        ekthick = dsqrt(-2.d0*av/f) 
      end if

! start main loop

      do i = icalcan,icalcen
      do k = kcalcan,kcalcen

! skip cell if dry land
      if (jc(i,k).le.0) goto 100

!*******************************************************************
!*** u-component  (eastern - cross-slope) **************************
!*******************************************************************

      termnl = 0.d0
      termnl2 = 0.d0
      hordif = 0.d0
      corx = 0.d0
      slorho = 0.d0
      islope = 0.d0
      redgu = 0.d0  

! plume thickness on the u-grid
      pdepu = 5.0d-1*(pdep(i,k) + pdep(i+1,k))

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! skip u-component if eastern cell is dry
      if (jc(i+1,k).le.0) goto 50

      islope = ipos(i+1,k) - ipos(i,k)

! final wet/dry logic
      if ((jcd(i,k).le.0).and.(jcd(i+1,k).le.0)) then
        goto 50
      else if ((jcd(i,k).le.0).and.(islope.gt.0.d0)) then
        goto 50
      else if ((jcd(i+1,k).le.0).and.(islope.le.0.d0)) then
        goto 50
      endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! control of negative depth effects by enhanced friction
      jcvfac = 0
      arfac = 1.d0
      jcvfac = max0(jcd_negdep(i,k),jcd_negdep(i+1,k))
      if (jcvfac .ge. 1) arfac = 75.d0*jcvfac  

! velocity components on the u-grid
      tt = 5.0d-1 
      uu = 5.0d-1*dy(k)*rdyv(k)
      vmid = tt*uu*(sv(i,k-1) + sv(i+1,k-1)) +  &
            tt*(1.d0-uu)*(sv(i,k) + sv(i+1,k))
      umid = su(i,k)

      speed = umid**2+vmid**2 + small

! drag coefficient on the u-grid
      dragu = 5.0d-1*(drag(i,k) + drag(i+1,k))

! 1)calculation of friction parameter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      tbotm = 1.d0/(1.d0+arfac*dragu*dsqrt(speed)*dt/(pdepu+dcr))

! turning angle theory

      if (tangle) then 

! set up stuff for rewriting drag if necessary

        if (draginmelt) then
          dragu = cdb
          ugriddrag(i,k) = 5.0d-1*(drag(i,k) + drag(i+1,k))
          newudrag(i,k) = .false.
        endif

! find ekman layer thickness if necessary and thickness ratio

        if (variableekman) then
          av = abs(0.16d0*0.4d0*dragu*speed / (f*2.71828d0))
          ekthick = dsqrt(-2.d0*av/f) 
        end if
        thickratio = 2.d0*pdepu/ekthick

! calculate velocities at base of ekman layer (eq. 32)

        kq = (ekthick/av)*dragu   &
              * dsqrt(speed*((u0a(i,k) + one)**2 + v0a(i,k)**2))
        u0(i,k) = - (one + kq)*kq / ((one + kq)**2 + one) 
        v0(i,k) = - kq / ((one + kq)**2 + one)

! calculate geostrophic velocity magnitude (eq. 36)

        kq = thickratio / dsqrt( &
                (thickratio + u0(i,k) - v0(i,k))**2 &
                           + (u0(i,k) + v0(i,k))**2  )  

! calculate final drag magnitude (eq. 42)

        kq = (kq*av/ekthick) * &
                 dsqrt( 2.d0*(u0(i,k)**2 + v0(i,k)**2) )
!         
        if (olddrag) kq = dragu*dsqrt(speed) 
!        
        if (draginmelt) then
          if (speed > small) then
            ugriddrag(i,k) = kq/dsqrt(speed)
          else
            ugriddrag(i,k) = cdb
          endif
          newudrag(i,k) = .true.
        endif
!           
! calculate turning angle          
!           
        tang(i,k) = atan2( u0(i,k)-v0(i,k) , -u0(i,k)-v0(i,k) ) &
       - atan2( u0(i,k)+v0(i,k) , thickratio+u0(i,k)-v0(i,k) )
!      
        if (norotation) tang(i,k) = 0.0 
!         
        costang = dcos(tang(i,k))
        sintang = dsin(tang(i,k))
        tbotm = one / (one + arfac*kq*costang*dt / (pdepu+dcr))

      end if    

! 2)baroclinic pressure gradient
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      zu = wcdep + gldep - 5.0d-1*(ipos(i,k) + ipos(i+1,k))
      zum = zu - 5.0d-1*pdepu      
      izo = int(zu/dzincr) + 1
      izu = izo + 1
      difu = dble(izo)*dzincr - zu
      difo = dzincr - difu
      rhoa =(difu*rhovf(izo)+difo*rhovf(izu))/dzincr

! density of water fraction
      salu = 5.0d-1*(salt(i,k)+salt(i+1,k))
      if (rholinear) then
        tu = 5.0d-1*(temp(i,k)+temp(i+1,k))
        rhoc = rho_func_linear(temp(i,k),salt(i,k))
        rhoe = rho_func_linear(temp(i+1,k),salt(i+1,k))
        rhou = rho_func_linear(tu,salu)
      else
        if (thermobar) then
          tu = 5.0d-1*(tins(i,k)+tins(i+1,k))
          rhoc = rho_func_nonlinear(temp(i,k),salt(i,k),zum*1.0d-1)
          rhoe = rho_func_nonlinear(temp(i+1,k),salt(i+1,k),zum*1.d-1)
          rhou = rho_func_nonlinear(tu,salu,zu*1.0d-1)
        else
          tu = 5.0d-1*(temp(i,k)+temp(i+1,k))
          rhoc = rho_func_nonlinear(temp(i,k),salt(i,k),0.d0)
          rhoe = rho_func_nonlinear(temp(i+1,k),salt(i+1,k),0.d0)
          rhou = rho_func_nonlinear(tu,salu,0.d0)
        end if
      end if 

! effects of frazil
      if (frazil) then
        ctotu = 5.0d-1*(ctot(i,k) + ctot(i+1,k))
        rhoc = (1.d0 - ctot(i,k))*rhoc + ctot(i,k)*rhoi
        rhoe = (1.d0 - ctot(i+1,k))*rhoe + ctot(i+1,k)*rhoi
        rhou = (1.d0 - ctotu)*rhou + ctotu*rhoi  
      end if

      rhoq = 5.0d-1*(rhou + rhoa)
      redgu = (rhoa - rhou)/rhoq
      slorho =5.0d-1*(rhoe-rhoc)*gdt*rdx(i)*pdepu*pdepu/rhoq

! 3)barotropic pressure gradient
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      islope = islope*rdx(i)*gdt*redgu*pdepu

! 4)coriolis
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      corx = cosfdt*ua(i,k) + pdepu*vmid*sinfdt

! if turning angle add other part of drag here because uses v not u

      if (tangle) then
        corx = corx + arfac*kq*sintang*dt*vmid
      end if

! 5&6)nonlinear terms (selective vector upstream)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (nonlin) then

!  choose upstream quarter
        idel = - int(sign(one,umid))
        kdel = - int(sign(one,vmid))
        ihilf = i + idel
        khilf = k + kdel
!  choose the relevant upstream cell dimensions 
        iidx = i + (idel+1)/2
        kkdy = k + (kdel-1)/2
        dxx = dxu(iidx)  
        dyy = dy(kkdy)
!  weighting coefficients
        sx = abs(umid)*dt/dxx
        sy = abs(vmid)*dt/dyy

        if ((sx.gt.one).or.(sy.gt.one)) then
          write(*,*) 'error: u courant exceeded at i=',i,' k=',k
          write(11,*) 'error: u courant exceeded at i=',i,' k=',k
        end if

        sxy = sx -sy
        r1 = (sign(one,sxy) + 1.d0)*5.0d-1
        r2 = one - r1
        termnl=(r1*sy + r2*sx)*(ua(ihilf,khilf)*dble(jcd_u(ihilf,khilf)) &
                + dble(1 - jcd_u(ihilf,khilf))*ua(i,k)) &
                + r1*sxy*(ua(ihilf,k)*dble(jcd_u(ihilf,k)) &
                + dble(1 - jcd_u(ihilf,k))*ua(i,k)) &
                - r2*sxy*(ua(i,khilf)*dble(jcd_u(i,khilf)) &
                + dble(1 - jcd_u(i,khilf))*ua(i,k)) &
                - (r2*sy + r1*sx)*ua(i,k)

        termnl2 = - ua(i,k)*dt*((su(ihilf,k) - su(i,k))*idel/dxx  &
                         + (sv(iidx,k) - sv(iidx,k-1))*rdyv(k))
      end if
!          
! 7)lateral shear stress terms (east component)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (horturb) then

        tlate = ahdxu(i+1)*(su(i+1,k) - su(i,k))
        tlatw = ahdxu(i)*(su(i,k) - su(i-1,k))
        tlats = ahdy(k-1)*(su(i,k) - su(i,k-1))
        tlatn = ahdy(k)*(su(i,k+1) - su(i,k))

        hordif = pdepu*((tlate-tlatw)*rdx(i) + (tlatn-tlats)*rdyv(k))

      end if

! final
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      u(i,k) = (corx+islope+slorho+termnl+termnl2+hordif)*tbotm

50    continue
!*******************************************************************
!*** v-component  (northern - up-slope) ****************************
!*******************************************************************

      termnl = 0.d0
      termnl2 = 0.d0
      cory  = 0.d0
      slorho = 0.d0
      islope = 0.d0
      redgv = 0.d0

! plume thickness on the v-grid
      pdepv = 5.0d-1*(pdep(i,k) + pdep(i,k+1))

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! skip u-component if northern cell is dry
      if (jc(i,k+1).le.0) goto 100

      islope = ipos(i,k+1) - ipos(i,k)

! final wet/dry logic
      if ((jcd(i,k).le.0).and.(jcd(i,k+1).le.0)) then
        goto 100
      else if ((jcd(i,k).le.0).and.(islope.gt.0.d0)) then
        goto 100
      else if ((jcd(i,k+1).le.0).and.(islope.le.0.d0)) then
        goto 100
      endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! control of negative depth effects by enhanced friction
      jcvfac = 0
      arfac = 1.d0
      jcvfac = max0(jcd_negdep(i,k),jcd_negdep(i,k+1))
      if (jcvfac .ge. 1) arfac = 75.d0*jcvfac  

! velocity components on the v-grid 
      uu = 5.0d-1*dx(i)*rdxu(i)
      tt = 5.0d-1
      umid = tt*uu*(su(i-1,k) + su(i-1,k+1)) + &
            tt*(1.d0-uu)*(su(i,k) + su(i,k+1))
      vmid = sv(i,k)

      speed = umid**2+vmid**2 + small

! drag coefficient on the v-grid
      dragv = 5.0d-1*(drag(i,k) + drag(i,k))

! 1)calculation of friction parameter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      tbotm = 1.d0/(1.d0+arfac*dragv*dsqrt(speed)*dt/(pdepv+dcr))

! turning angle theory

      if (tangle) then 

        if (draginmelt) dragv = cdb

! find ekman layer thickness if necessary and thickness ratio

        if (variableekman) then
          av = abs(0.16d0*0.4d0*dragv*speed / (f*2.71828d0))
          ekthick = dsqrt(-2.d0*av/f)
        end if
        thickratio = 2.d0*pdepv/ekthick

! calculate velocities at base of ekman layer (eq. 32)

        kq = (ekthick/av)*dragv &
              * dsqrt(speed*((u0a(i,k) + one)**2 + v0a(i,k)**2))
        u0(i,k) = - (one + kq)*kq / ((one + kq)**2 + one) 
        v0(i,k) = - kq / ((one + kq)**2 + one)

! calculate geostrophic velocity magnitude (eq. 36)

        kq = thickratio / dsqrt( &
                (thickratio + u0(i,k) - v0(i,k))**2    &
                           + (u0(i,k) + v0(i,k))**2)  

! calculate final drag magnitude (eq. 42)

        kq = (kq*av/ekthick) * &
                 dsqrt( 2.d0*(u0(i,k)**2 + v0(i,k)**2) )

        if (olddrag) kq = dragv*dsqrt(speed) 

! calculate turning angle          
!           
        tang(i,k) = atan2( u0(i,k)-v0(i,k) , -u0(i,k)-v0(i,k) ) &
       - atan2( u0(i,k)+v0(i,k) , thickratio+u0(i,k)-v0(i,k) )

        if (norotation) tang(i,k) = 0.0 

        costang = dcos(tang(i,k))
        sintang = dsin(tang(i,k))
        tbotm = one / (one + arfac*kq*costang*dt / (pdepv+dcr))

      end if      

! 2)baroclinic pressure gradient
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      zv = wcdep + gldep - 5.0d-1*(ipos(i,k) + ipos(i,k+1))
      zvm = zv - 5.0d-1*pdepv     
      izo = int(zv/dzincr) + 1
      izu = izo + 1
      difu = dble(izo)*dzincr - zv
      difo = dzincr - difu
      rhoa =(difu*rhovf(izo)+difo*rhovf(izu))/dzincr

! density of water fraction
      salv = 5.0d-1*(salt(i,k)+salt(i,k+1))
      if (rholinear) then
        tv = 5.0d-1*(temp(i,k)+temp(i,k+1))
        rhoc = rho_func_linear(temp(i,k),salt(i,k))
        rhon = rho_func_linear(temp(i,k+1),salt(i,k+1))
        rhov = rho_func_linear(tv,salv)
      else
        if (thermobar) then
          tv = 5.0d-1*(tins(i,k)+tins(i,k+1))
          rhoc = rho_func_nonlinear(temp(i,k),salt(i,k),zvm*1.0d-1)
          rhon = rho_func_nonlinear(temp(i,k+1),salt(i,k+1),zvm*1.d-1)
          rhov = rho_func_nonlinear(tv,salv,zv*1.0d-1)
        else
          tv = 5.0d-1*(temp(i,k)+temp(i,k+1))
          rhoc = rho_func_nonlinear(temp(i,k),salt(i,k),0.d0)
          rhon = rho_func_nonlinear(temp(i,k+1),salt(i,k+1),0.d0)
          rhov = rho_func_nonlinear(tv,salv,0.d0)
        end if
      end if 

! effects of frazil
      if (frazil) then
        ctotv = 5.0d-1*(ctot(i,k) + ctot(i,k+1))
        rhoc = (1.d0 - ctot(i,k))*rhoc + ctot(i,k)*rhoi
        rhon = (1.d0 - ctot(i,k+1))*rhon + ctot(i,k+1)*rhoi
        rhov = (1.d0 - ctotv)*rhov + ctotv*rhoi  
      end if

      rhoq = 5.0d-1*(rhov + rhoa)
      redgv = (rhoa - rhov)/rhoq
      slorho = 5.0d-1*(rhon-rhoc)*gdt*rdy(k)*pdepv*pdepv/rhoq

! 3)barotropic pressure gradient
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      islope = islope*rdy(k)*gdt*redgv*pdepv

! 4)coriolis
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      cory = cosfdt*va(i,k) - pdepv*umid*sinfdt

! if turning angle add other part of drag here because uses u not v

      if (tangle) then
        cory = cory - arfac*kq*sintang*dt*umid
      end if      

! 5&6)nonlinear terms (selective vector upstream)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (nonlin) then

!  choose upstream quarter
        idel = - int(sign(one,umid))
        kdel = - int(sign(one,vmid))
        ihilf = i + idel
        khilf = k + kdel
!  choose the relevant upstream cell dimensions 
        iidx = i + (idel-1)/2
        kkdy = k + (kdel+1)/2
        ihilf = idel + i
        khilf = kdel + k
        dxx = dx(iidx)
        dyy = dyv(kkdy)
! weighting coefficients
        sx = abs(umid)*dt/dxx
        sy = abs(vmid)*dt/dyy
! test courant number
        if ((sx.gt.one).or.(sy.gt.one)) then
          write(*,*) 'error: v courant exceeded at i=',i,' k=',k
          write(11,*) 'error: v courant exceeded at i=',i,' k=',k
        end if
! calculate nonlinear terms
        sxy = sx -sy
        r1 = (dsign(one,sxy) + 1.d0)*5.0d-1
        r2 = one - r1
        termnl=(r1*sy + r2*sx)*(va(ihilf,khilf)*dble(jcd_v(ihilf,khilf)) &
                + dble(1 - jcd_v(ihilf,khilf))*va(i,k)) &
                + r1*sxy*(va(ihilf,k)*dble(jcd_v(ihilf,k)) &
                + dble(1 - jcd_v(ihilf,k))*va(i,k)) &
                - r2*sxy*(va(i,khilf)*dble(jcd_v(i,khilf)) &
                + dble(1 - jcd_v(i,khilf))*va(i,k)) &
                - (r2*sy + r1*sx)*va(i,k)

         termnl2 = - va(i,k)*dt*((sv(i,khilf) - sv(i,k))*kdel/dyy  &
                          + (su(i,kkdy) - su(i-1,kkdy))*rdxu(i))
      end if
!          
! 7)lateral shear stress terms (north component)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (horturb) then

        tlate = ahdx(i)*(sv(i+1,k) - sv(i,k))
        tlatw = ahdx(i-1)*(sv(i,k) - sv(i-1,k))
        tlats = ahdyv(k)*(sv(i,k) - sv(i,k-1))
        tlatn = ahdyv(k+1)*(sv(i,k+1) - sv(i,k))

        hordif = pdepv*((tlate-tlatw)*rdxu(i) + (tlatn-tlats)*rdy(k))

      end if
!          
! final
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      v(i,k) = (cory+islope+slorho+termnl+termnl2+hordif)*tbotm

100   continue
      end do
      end do
!     
! write new drag on scalar grid if necessary 
!     
      if (tangle.and.draginmelt) then
!     
        do i = icalcan+1,icalcen
          do k = kcalcan,kcalcen
            if (newudrag(i-1,k).or.newudrag(i,k)) then
              drag(i,k) = 0.5d0*(ugriddrag(i-1,k) + ugriddrag(i,k))
            end if
          end do
        end do
!     
      endif
!     
!*******************************************************************
!*** final velocities **********************************************
!*******************************************************************

      do i = icalcan,icalcen
        do k = kcalcan,kcalcen
          pdepu = 5.0d-1*(pdep(i,k) + pdep(i+1,k)) + small
          pdepv = 5.0d-1*(pdep(i,k) + pdep(i,k+1)) + small
          su(i,k) = u(i,k)/pdepu
          sv(i,k) = v(i,k)/pdepv
! test for negative depths and set flow to zero if so
          delta = dt*(rdxu(i)*(u(i-1,k)-u(i,k)) &
                   + rdyv(k)*(v(i,k-1)-v(i,k)))
          itest = bpos(i,k) - ipos(i,k) + delta 
          if (itest .lt. 0.d0) then
            u(i,k) = 0.d0
            su(i,k) = 0.d0
            u(i-1,k) = 0.d0
            su(i-1,k) = 0.d0
            v(i,k) = 0.d0
            sv(i,k) = 0.d0
            v(i,k-1) = 0.d0
            sv(i,k-1) = 0.d0
          endif                 
          ua(i,k) = u(i,k)
          va(i,k) = v(i,k)
          u0a(i,k) = u0(i,k)
          v0a(i,k) = v0(i,k)
        end do
      end do            

      return
end subroutine

end module
