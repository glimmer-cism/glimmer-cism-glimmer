module plume_update

  use plume_global

contains

subroutine update(iwetmin,iwetmax,kwetmin,kwetmax, &
                                icalcan,icalcen,kcalcan,kcalcen,negdep)

! update depths, depth-averaged velocities (not depth-integrated),
! and wetted area

  implicit none

! local variables

      integer i,k,iwetmin,iwetmax,kwetmin,kwetmax
      integer icalcan,icalcen,kcalcan,kcalcen,izo,izu

      real(kind=kdp) iconti,error,pdepold,pdepc,pdepe
      real(kind=kdp) dunew,pdepn,dvnew,zd,difu,difo,negdep

! update interface position


      do i = icalcan,icalcen
        do k = kcalcan,kcalcen

! interface position corrected for negative depths at new timestep
          iconti = ipos(i,k)
          ipos(i,k) = dmin1(bpos(i,k),ipos(i,k))

! report any errors
          error = abs(iconti - ipos(i,k))
          negdep = negdep+error               
          if (error.gt.5.0d-1) then
            write(*,*) 'error: negative depth ',error,' at i=', &
                                                        i,' k=',k
            write(11,*) 'error: negative depth ',error,' at i=', &
                                                        i,' k=',k
          end if

          sv(i,k) = 0.d0
          su(i,k) = 0.d0

! plume thickness and wetted area (main update of pdep)
          jcd_fl(i,k) = 0      
          pdepold = pdep(i,k)
          jcd(i,k) = 0
          pdep(i,k) = bpos(i,k) - ipos(i,k)

! find newly wet area
          if (pdepold .eq. 0.d0 .and. pdep(i,k) .gt. 0.d0) then
            jcd_fl(i,k) = 1 
          end if

! find wetted area boundaries (main update of jcd)
          if (pdep(i,k).ge.dcr) then
            jcd(i,k) = 1
            if (i.lt.iwetmin) then
              iwetmin = i
            else 
              if (i.gt.iwetmax) iwetmax = i
            end if
            if (k.lt.kwetmin) then
              kwetmin = k
            else 
              if (k.gt.kwetmax) kwetmax = k
            end if
          end if

        end do
      end do

! update velocities
      do i = icalcan,icalcen
        do k = kcalcan,kcalcen
          jcd_u(i,k) = 0
          jcd_v(i,k) = 0
          pdepc = bpos(i,k) - ipos(i,k)
! u-component
          pdepe = bpos(i+1,k) - ipos(i+1,k)
          dunew = 5.0d-1*(pdepe + pdepc)
          if (dunew.gt.small) then
            su(i,k) = u(i,k)/dunew
            jcd_u(i,k) = 1
          endif
! v-component
          pdepn = bpos(i,k+1) - ipos(i,k+1)
          dvnew = 5.0d-1*(pdepn + pdepc)
          if (dvnew.gt.small) then
            sv(i,k) = v(i,k)/dvnew
            jcd_v(i,k) = 1
          endif
        end do
      end do

! interpolate ambient density field experienced
      do i = icalcan,icalcen
        do k = kcalcan,kcalcen
          zd = gldep + wcdep - ipos(i,k)
          izo = int(zd/dzincr) + 1
          izu = izo + 1
          difu = dble(izo)*dzincr - zd
          difo = dzincr - difu
          rhoamb(i,k) =(difu*rhovf(izo)+difo*rhovf(izu))/dzincr
        end do
      end do

      return
end subroutine

end module
