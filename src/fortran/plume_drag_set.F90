module plume_drag_set

  use plume_global

contains

subroutine drag_set()

  implicit none

! set variable drag-coefficient areas by hand 

! local variables

      integer i,k,imin,imax,kmin,kmax,icen,kcen,ii,kk

      real(kind=kdp) rotfac

! ---------------------------
! set up drag regions by hand
! ---------------------------
! unless otherwise indicated the drag regions are defined such that all cells 
! within a rotated box (defined by limits and then rotated about centre 
! according to rotfac) have a different drag to the background value 

! east of henry ice rise
! ----------------------

      if (drflag(01)) then

!c        rotfac = +055.0d0*pi/180.d0
!c        imin = int(0560.d0*1000.d0/hx)
!c        imax = int(0620.d0*1000.d0/hx)
!c        kmin = int(0370.d0*1000.d0/hy)
!c        kmax = int(0500.d0*1000.d0/hy)
        rotfac = +000.0d0*pi/180.d0
        imin = int(0225.d0*1000.d0/hx)
        imax = int(0400.d0*1000.d0/hx)
        kmin = int(0000.d0*1000.d0/hy)
        kmax = int(0110.d0*1000.d0/hy)

        icen = imin + nint(dble(imax-imin)/2.d0)
        kcen = kmin + nint(dble(kmax-kmin)/2.d0)

        do i = 1,m
          do k = 1,n

            ii = nint(dble(i-icen)*cos(rotfac)) &
     &              + nint(dble(k-kcen)*sin(rotfac)) + icen
            kk = nint(dble(k-kcen)*cos(rotfac)) &
     &              - nint(dble(i-icen)*sin(rotfac)) + kcen

            if (((ii.ge.imin).and.(ii.le.imax)) &
     &                   .and.((kk.ge.kmin).and.(kk.le.kmax))) then
              drag(i,k) = cdbvar
            end if

          end do
        end do

      end if

! east of berkner island
! ----------------------

      if (drflag(02)) then

        rotfac = +045.0d0*pi/180.d0
        imin = int(0600.d0*1000.d0/hx)
        imax = int(0730.d0*1000.d0/hx)
        kmin = int(0700.d0*1000.d0/hy)
        kmax = int(0830.d0*1000.d0/hy)

        icen = imin + nint(dble(imax-imin)/2.d0)
        kcen = kmin + nint(dble(kmax-kmin)/2.d0)

        do i = 1,m
          do k = 1,n

            ii = nint(dble(i-icen)*cos(rotfac)) &
     &              + nint(dble(k-kcen)*sin(rotfac)) + icen
            kk = nint(dble(k-kcen)*cos(rotfac)) &
     &              - nint(dble(i-icen)*sin(rotfac)) + kcen

            if (((ii.ge.imin).and.(ii.le.imax)) &
     &                   .and.((kk.ge.kmin).and.(kk.le.kmax))) then
              drag(i,k) = cdbvar
            end if

          end do
        end do

      end if

! spare slot
! ----------

      if (drflag(03)) then

        write(*,*) 'error: drag area not set'

      end if

! mask regions with grounding line and ice front
! ----------------------------------------------

      do i = 1,m
        do k = 1,n
          if ((bpos(i,k).le.0.d0).or.(bpos(i,k).ge.(gldep+wcdep))) then
            drag(i,k) = 0.d0
          end if
        end do
      end do

      return
end subroutine

end module
