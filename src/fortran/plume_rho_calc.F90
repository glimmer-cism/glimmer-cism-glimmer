module plume_rho_calc

  use plume_global
  use plume_functions

contains

subroutine rho_calc(icalcan,icalcen,kcalcan,kcalcen,sepflag)

! calculate density and warn if plume has reached separation point

  implicit none

! local variables

      logical sepflag

      integer i,k,icalcan,icalcen,kcalcan,kcalcen,izo,izu

      real(kind=kdp)p,rhoatmp,depth,difu,difo


! get density of water fraction. if thermobaricity included, compute
! in situ temp from potential temp. for future reference
      if (rholinear) then
        do i = icalcan,icalcen
          do k = kcalcan,kcalcen
            rhop(i,k) = rho_func_linear(temp(i,k),salta(i,k)) 
          end do
        end do
      else
        if (thermobar) then
          do i = icalcan,icalcen
            do k = kcalcan,kcalcen
              p = (wcdep + gldep - ipos(i,k) - 5.0d-1*pdep(i,k))*1.0d-1
              tins(i,k) = tinsitu_func(temp(i,k),salta(i,k),p)  
              rhop(i,k) = rho_func_nonlinear(tins(i,k),salta(i,k),p) 
            end do
          end do
        else
          do i = icalcan,icalcen
            do k = kcalcan,kcalcen
              rhop(i,k) = rho_func_nonlinear(temp(i,k),salta(i,k),0.d0) 
            end do
          end do
        end if
      end if

! add density of ice fraction if appropriate
      if (frazil) then
        do i = icalcan,icalcen
          do k = kcalcan,kcalcen
            rhop(i,k) = (1.d0 - ctot(i,k))*rhop(i,k) + ctot(i,k)*rhoi
          end do
        end do
      end if

! calculate density of ambient fluid and test for separation
      do i = icalcan,icalcen
        do k = kcalcan,kcalcen
          if (pdep(i,k).gt.dcr) then
            depth = wcdep + gldep - ipos(i,k) - 5.0d-1*pdep(i,k)
            izo = int(depth/dzincr) + 1
            izu = izo + 1
            difu = dble(izo)*dzincr - depth
            difo = dzincr - difu
            rhoatmp = (difu*rhovf(izo)+difo*rhovf(izu))/dzincr

            if (((rhoatmp - rhop(i,k)).le.septol).and. &
                                               (.not.sepflag)) then
              write(*,*) 'warning: plume separated at i=',i,' k=',k
              write(11,*) 'warning: plume separated at i=',i,' k=',k
              sepflag = .true.
            end if
!           
          end if
        end do
      end do

      return
end subroutine

end module

