module plume_inflow

use plume_global
use plume_functions

contains

subroutine inflow_calc(icalcan,kcalcan,icalcen,kcalcen)


  implicit none

! local variables

      integer i,k,l,icalcan,kcalcan,icalcen,kcalcen

      real(kind=kdp) pressure,ttt

      do i = icalcan,icalcen
        do k = kcalcan,kcalcen
          if (depinf(i,k).gt.0.d0) then

! reset inflow depth

            ipos(i,k) = bpos(i,k) - depinf(i,k)         
            pdep(i,k) = depinf(i,k)

            if (pdep(i,k).ge.dcr) then           

! reset inflow properties 

              jcd(i,k) = 1
              salt(i,k) = saltinf(i,k)
              salta(i,k) = saltinf(i,k)
              temp(i,k) = tempinf(i,k)
              tempa(i,k) = tempinf(i,k)

! calculate density of inflow water fraction

              if (rholinear) then
                rhop(i,k) = rho_func_linear(tempinf(i,k),saltinf(i,k))
              else
                if (thermobar) then
                  pressure = 1.0d-1*(gldep + wcdep - ipos(i,k))
                  ttt = tinsitu_func(tempinf(i,k),saltinf(i,k),pressure)
                  rhop(i,k) =  &
     &                   rho_func_nonlinear(ttt,saltinf(i,k),pressure)
                else
                  rhop(i,k) = &
     &              rho_func_nonlinear(tempinf(i,k),saltinf(i,k),0.d0)
                end if
              end if

! set frazil and add density of ice fraction

              if (frazil) then
                do l = 1,nice
                  c(i,k,l) = cinf(l)
                  ca(i,k,l) = cinf(l)
                end do
                ctot(i,k) = cinftot
                ctota(i,k) = cinftot
                rhop(i,k) = (1.d0 - cinftot)*rhop(i,k) + cinftot*rhoi
              end if

! set inflow tracers

              if (intrace) then
                do l = ninfmin,ninfmax
                  if (intrin(i,k).eq.l) then
                    intr(i,k,l) = 1.d0
                    intra(i,k,l) = 1.d0
                  else
                    intr(i,k,l) = 0.d0
                    intra(i,k,l) = 0.d0
                  end if 
                end do
              end if

            end if

          end if
        end do
      end do

      return
end subroutine

subroutine inflow_set_isomip()

  implicit none

! set isomip inflow depths by hand if using real topography

! local variables

      integer i,k

! -----------------------------
! set up inflow regions by hand
! -----------------------------
! unless otherwise indicated the inflows are set up by defining a square
! region and then setting all ocean domain within that region to inflow if it
! is deeper than a given depth.

! the inflows are:
!  01 - whole southern boundary

! whole southern boundary
! -----------------------

      if (inflag(01)) then

        do i = m/2-5,m/2+5
          do k = n/2-5,n/2+5
            intrin(i,k) = 1
            depinf(i,k) = depinffix
          end do
        end do

      end if

! ----------------------------------------------
! mask regions with grounding line and ice front
! ----------------------------------------------

      do i = 1,m
        do k = 1,n
          if ((bpos(i,k).le.0.d0).or.(bpos(i,k).ge.(gldep+wcdep))) then
            intrin(i,k) = 0
            depinf(i,k) = 0.d0
          end if
        end do
      end do

      return
end subroutine

subroutine inflow_set_fris()

  implicit none

! set fris inflow depths by hand if using real topography

! local variables

      integer i,k,imin,imax,kmin,kmax

      real(kind=kdp) cutdepth


! -----------------------------
! set up inflow regions by hand
! -----------------------------
! unless otherwise indicated the inflows are set up by defining a square
! region and then setting all ocean domain within that region to inflow if it
! is deeper than a given depth.

! the inflows are:
!  01 - evans ice stream
!  02 - carlson inlet
!  03 - rutford ice stream
!  04 - institute ice stream
!  05 - mollereisstrom
!  06 - foundation ice stream
!  07 - support force glacier
!  08 - recovery glacier
!  09 - slessor glacier

! evans ice stream
! ----------------

      if (inflag(01)) then

        cutdepth = 0900.d0
        imin = int(0050.d0*1000.d0/hx)
        imax = int(0150.d0*1000.d0/hx) 
        kmin = int(0200.d0*1000.d0/hy)
        kmax = int(0300.d0*1000.d0/hy)

        do i = imin,imax
          do k = kmin,kmax
            if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
              intrin(i,k) = 1
              depinf(i,k) = depinffix
            end if
          end do
        end do

      end if

! carlson inlet
! -------------

      if (inflag(02)) then

        cutdepth = 1400.d0
        imin = int(0150.d0*1000.d0/hx)
        imax = int(0210.d0*1000.d0/hx)
        kmin = int(0050.d0*1000.d0/hy)
        kmax = int(0150.d0*1000.d0/hy)

        do i = imin,imax
          do k = kmin,kmax
            if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
              intrin(i,k) = 2
              depinf(i,k) = depinffix
            end if
          end do
        end do

      end if

! rutford ice stream
! ------------------

      if (inflag(03)) then

        cutdepth = 1500.d0
        imin = int(0200.d0*1000.d0/hx)
        imax = int(0300.d0*1000.d0/hx)
        kmin = int(0000.d0*1000.d0/hy)
        kmax = int(0070.d0*1000.d0/hy)

        do i = imin,imax
          do k = kmin,kmax
            if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
              intrin(i,k) = 3
              depinf(i,k) = depinffix
            end if
          end do
        end do

      end if

! institute ice stream
! --------------------

      if (inflag(04)) then

        cutdepth = 1000.d0
        imin = int(0450.d0*1000.d0/hx)
        imax = int(0600.d0*1000.d0/hx)
        kmin = int(0100.d0*1000.d0/hy)
        kmax = int(0200.d0*1000.d0/hy)

        do i = imin,imax
          do k = kmin,kmax
            if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
              intrin(i,k) = 4
              depinf(i,k) = depinffix
            end if
          end do
        end do

      end if

! mollereisstrom
! --------------

      if (inflag(05)) then

        cutdepth = 1050.d0
        imin = int(0650.d0*1000.d0/hx)
        imax = int(0750.d0*1000.d0/hx)
        kmin = int(0200.d0*1000.d0/hy)
        kmax = int(0300.d0*1000.d0/hy)

        do i = imin,imax
          do k = kmin,kmax
            if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
              intrin(i,k) = 5
              depinf(i,k) = depinffix
            end if
          end do
        end do

      end if

! foundation ice stream
! ---------------------

      if (inflag(06)) then

        cutdepth = 1500.d0
        imin = int(0800.d0*1000.d0/hx)
        imax = int(0900.d0*1000.d0/hx)
        kmin = int(0200.d0*1000.d0/hy)
        kmax = int(0300.d0*1000.d0/hy)

        do i = imin,imax
          do k = kmin,kmax
            if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
              intrin(i,k) = 6
              depinf(i,k) = depinffix
            end if
          end do
        end do

      end if

! support force glacier
! ---------------------

      if (inflag(07)) then

        cutdepth = 1300.d0
        imin = int(0850.d0*1000.d0/hx)
        imax = int(0950.d0*1000.d0/hx)
        kmin = int(0450.d0*1000.d0/hy)
        kmax = int(0550.d0*1000.d0/hy)

        do i = imin,imax
          do k = kmin,kmax
            if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
              intrin(i,k) = 7
              depinf(i,k) = depinffix
            end if
          end do
        end do

      end if

! recovery glacier
! ----------------

      if (inflag(08)) then

        cutdepth = 1000.d0
        imin = int(0850.d0*1000.d0/hx)
        imax = int(0950.d0*1000.d0/hx)
        kmin = int(0630.d0*1000.d0/hy)
        kmax = int(0750.d0*1000.d0/hy)

        do i = imin,imax
          do k = kmin,kmax
            if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
              intrin(i,k) = 8
              depinf(i,k) = depinffix
            end if
          end do
        end do

      end if

! slessor glacier
! ---------------

      if (inflag(09)) then

        cutdepth = 1200.d0
        imin = int(0900.d0*1000.d0/hx)
        imax = int(1000.d0*1000.d0/hx)
        kmin = int(0800.d0*1000.d0/hy)
        kmax = int(0900.d0*1000.d0/hy)

        do i = imin,imax
          do k = kmin,kmax
            if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
              intrin(i,k) = 9
              depinf(i,k) = depinffix
            end if
          end do
        end do

      end if

! ----------------------------------------------
! mask regions with grounding line and ice front
! ----------------------------------------------

      do i = 1,m
        do k = 1,n
          if ((bpos(i,k).le.0.d0).or.(bpos(i,k).ge.(gldep+wcdep))) then
            intrin(i,k) = 0
            depinf(i,k) = 0.d0
          end if
        end do
      end do

      return
 end subroutine

subroutine inflow_set_larsen()

  implicit none

! set larsen inflow depths by hand if using real topography

! local variables

      integer i,k,imin,imax,kmin,kmax

      real(kind=kdp) cutdepth

! -----------------------------
! set up inflow regions by hand
! -----------------------------
! unless otherwise indicated the inflows are set up by defining a square
! region and then setting all ocean domain within that region to inflow if it
! is deeper than a given depth.

! the inflows are:
!  01 - attlee glacier
!  02 - whole larsen b
!  03 - whole larsen c

! attlee glacier
! --------------

      if (inflag(01)) then

        cutdepth = 1000.d0
        imin = int(0065.d0*1000.d0/hx)
        imax = int(0080.d0*1000.d0/hx) 
        kmin = int(0290.d0*1000.d0/hy)
        kmax = int(0310.d0*1000.d0/hy)

        do i = imin,imax
          do k = kmin,kmax
            if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
              intrin(i,k) = 1
              depinf(i,k) = depinffix
            end if
          end do
        end do

      end if

! whole larsen B
! --------------

      if (inflag(02)) then

        cutdepth = 0400.d0

        ! south
        imin = int(0130.d0*1000.d0/hx)
        imax = int(0160.d0*1000.d0/hx) 
        kmin = int(0320.d0*1000.d0/hy)
        kmax = int(0425.d0*1000.d0/hy)
        do i = imin,imax
          do k = kmin,kmax
            if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
              intrin(i,k) = 2
              depinf(i,k) = depinffix
            end if
          end do
        end do

        ! north
        imin = int(0130.d0*1000.d0/hx)
        imax = int(0200.d0*1000.d0/hx) 
        kmin = int(0425.d0*1000.d0/hy)
        kmax = int(0440.d0*1000.d0/hy)
        do i = imin,imax
          do k = kmin,kmax
            if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
              intrin(i,k) = 2
              depinf(i,k) = depinffix
            end if
          end do
        end do

      end if

! whole larsen C (set all>cutdepth to inflow, then remove promontories)
! ---------------------------------------------------------------------

      if (inflag(03)) then

        cutdepth = 1000.d0

        ! set all to inflow west of jason peninsula
        imin = int(0000.d0*1000.d0/hx)
        imax = int(0150.d0*1000.d0/hx) 
        kmin = int(0000.d0*1000.d0/hy)
        kmax = int(0310.d0*1000.d0/hy)
        do i = imin,imax
          do k = kmin,kmax
            if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
              intrin(i,k) = 3
              depinf(i,k) = depinffix
            end if
          end do
        end do

        ! unset cape alexander
        imin = int(0120.d0*1000.d0/hx)
        imax = int(0140.d0*1000.d0/hx) 
        kmin = int(0240.d0*1000.d0/hy)
        kmax = int(0270.d0*1000.d0/hy)
        do i = imin,imax
          do k = kmin,kmax
            intrin(i,k) = 0
            depinf(i,k) = 0.d0
          end do
        end do

        ! unset cole peninsula
        imin = int(0075.d0*1000.d0/hx)
        imax = int(0090.d0*1000.d0/hx) 
        kmin = int(0210.d0*1000.d0/hy)
        kmax = int(0240.d0*1000.d0/hy)
        do i = imin,imax
          do k = kmin,kmax
            intrin(i,k) = 0
            depinf(i,k) = 0.d0
          end do
        end do

        ! unset kenyon peninsula
        imin = int(0070.d0*1000.d0/hx)
        imax = int(0140.d0*1000.d0/hx) 
        kmin = int(0020.d0*1000.d0/hy)
        kmax = int(0060.d0*1000.d0/hy)
        do i = imin,imax
          do k = kmin,kmax
            intrin(i,k) = 0
            depinf(i,k) = 0.d0
          end do
        end do

      end if

! ----------------------------------------------
! mask regions with grounding line and ice front
! ----------------------------------------------

      do i = 1,m
        do k = 1,n
          if ((bpos(i,k).le.0.d0).or.(bpos(i,k).ge.(gldep+wcdep))) then
            intrin(i,k) = 0
            depinf(i,k) = 0.d0
          end if
        end do
      end do

      return
end subroutine

end module
