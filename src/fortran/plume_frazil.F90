module plume_frazil

  use plume_global

contains

subroutine frazil_calc

  implicit none

! set up various parameter arrays used in the frazil calculations

! local variables

      integer l
      real(kind=kdp) cdc,winew,rey,mstar

      do l=1,nice

! calculate effective radii

        re(l) = r(l)*(1.5d0*ar)**(1.d0/3.d0)

! calculate crystal thicknesses

        thi(l) = 2.d0*r(l)*ar

! calculate crystal volumes

        vol(l) = pi*r(l)*r(l)*thi(l)

! iteratively calculate frazil rising velocity

        cdc = 10.d0
        wi(l) = 0.d0
        winew = dsqrt((4.d0*g*ar*r(l)*(rho0-rhoi))/(rho0*cdc))
   40   if (abs(wi(l)-winew).gt.small) then
          wi(l) = winew
          rey = 2.d0*winew*r(l)/nu0
          cdc = 10**(1.386d0 - 0.892d0*dlog10(rey) &
     &                                       + 0.111d0*(dlog10(rey))**2)
          winew = dsqrt((4.d0*g*ar*r(l)*(rho0-rhoi))/(rho0*cdc))
          goto 40
        end if
        wi(l) = winew

! set the minimum of each size class to be one crystal

        cmin(l) = vol(l)

! set the nusselt number for each size class according to option chosen

! if nusselt number set to positive, use that constant value
        nuss(l) = nus

! correct full variable formulation
        if (nus.eq.-1.d0) then
          write(*,*) 'error: nusselt number option not coded yet'
        end if

! correct variable formulation with no turbulent part for large crystals
        if (nus.eq.-2.d0) then
          mstar = r(l)/(nu0**3.d0/eps)**2.5d-1
          if (mstar.lt.(1.d0/dsqrt(pr))) then
            nuss(l) = 1.d0 + mstar*1.7d-1*dsqrt(pr)
          else
            nuss(l) = 1.d0 + mstar*5.5d-1*(pr/mstar)**(1.d0/3.d0)
          end if 
        end if

! incorrect full hammar & shen formulation
        if (nus.eq.-3.d0) then
          write(*,*) 'error: nusselt number option not coded yet'
        end if

! incorrect hammar & shen formulation with no turbulent part for large crystals
        if (nus.eq.-4.d0) then
          mstar = r(l)/(nu0**3.d0/eps)**2.5d-1
          if (mstar.lt.(1.d0/dsqrt(pr))) then
            nuss(l) = 1.d0/mstar + 1.7d-1*dsqrt(pr)
          else
            nuss(l) = 1.d0/mstar + 5.5d-1*(pr/mstar)**(1.d0/3.d0)
          end if 
        end if

      end do

      return
end subroutine

end module
