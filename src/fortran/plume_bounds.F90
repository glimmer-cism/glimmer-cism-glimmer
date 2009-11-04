module plume_bounds

  use plume_global
  implicit none

contains

  subroutine bound_u_v(icalcan,icalcen,kcalcan,kcalcen)

    implicit none

    ! calculates boundary values for velocity and transports
    ! from neumann conditions. 

    ! local variables

    integer i,k,icalcan,icalcen,kcalcan,kcalcen

    ! southern boundary
    if (kcalcan.le.2) then
       do i = 1,m-1
          su(i,1) = su(i,2)
          u(i,1) = su(i,2)*5.0d-1*(pdep(i,2) + pdep(i+1,2))
          sv(i,1) = sv(i,2)
          v(i,1) = sv(i,2)*pdep(i,2)
       end do
    end if

    ! northern boundary
    if (kcalcen.ge.n-1) then
       do i = 1,m-1
          su(i,n) = su(i,n-1)
          u(i,n) = su(i,n-1)*5.0d-1*(pdep(i,n-1) + pdep(i+1,n-1))
          sv(i,n) = sv(i,n-1)
          v(i,n) = sv(i,n-1)*pdep(i,n-1)
       end do
    end if

    ! western boundary
    if (icalcan.le.2) then
       do k = 1,n-1
          su(1,k) = su(2,k)
          u(1,k) = su(2,k)*pdep(2,k)
          sv(1,k) = sv(2,k)
          v(1,k) = sv(2,k)*5.0d-1*(pdep(2,k)+pdep(2,k+1))
       end do
    end if

    ! eastern boundary
    if (icalcen.ge.m-1) then
       do k = 1,n-1
          su(m,k) = su(m-1,k)
          u(m,k) = su(m-1,k)*pdep(m-1,k)
          sv(m,k) = sv(m-1,k)
          v(m,k) = sv(m-1,k)*5.0d-1*(pdep(m-1,k)+pdep(m-1,k+1))
       end do
    end if

    return
  end subroutine bound_u_v

  !*********************************************************************
  !*********************************************************************
  !*********************************************************************

  subroutine bound_interface(icalcan,icalcen,kcalcan,kcalcen)

    ! calculates boundary values for interface depth
    ! from neumann conditions. 

    implicit none

    ! local variables

    integer :: i,k,icalcan,icalcen,kcalcan,kcalcen
    real(kind=kdp):: pdepold

    ! southern boundary
    if(kcalcan.le.2) then
       do i = 1,m
          jcd_fl(i,1) = 0       
          jcd(i,1) = 0
          ipos(i,1) = min(bpos(i,1),ipos(i,1))
          pdepold = pdep(i,1)
          pdep(i,1) = bpos(i,1) - ipos(i,1)
          if (pdep(i,1).ge.dcr) jcd(i,1) = 1
          if (pdepold.lt.small .and. pdep(i,1) .ge. small) &
               jcd_fl(i,1) = 1      
       end do
    end if

    ! northern boundary
    if(kcalcen.ge.n-1) then
       do i = 1,m
          jcd_fl(i,n) = 0       
          jcd(i,n) = 0
          ipos(i,n) = min(bpos(i,n),ipos(i,n))
          pdepold = pdep(i,n)
          pdep(i,n) = bpos(i,n) - ipos(i,n)
          if (pdep(i,n).ge.dcr) jcd(i,n) = 1
          if (pdepold.lt.small .and. pdep(i,n) .ge. small) &
               jcd_fl(i,n) = 1      
       end do
    end if

    ! western boundary
    if (icalcan.le.2) then
       do k = 1,n
          jcd(1,k) = jcd(2,k)
          jcd_fl(1,k) = jcd_fl(2,k)
          ipos(1,k) = ipos(2,k)
          pdep(1,k) = pdep(2,k)
       end do
    end if

    ! eastern boundary
    if (icalcen.ge.m-1) then
       do k = 1,n
          jcd(m,k) = jcd(m-1,k)
          jcd_fl(m,k) = jcd_fl(m-1,k)
          ipos(m,k) = ipos(m-1,k)
          pdep(m,k) = pdep(m-1,k)
       end do
    end if

    return
  end subroutine bound_interface

  !*********************************************************************
  !*********************************************************************
  !*********************************************************************
  subroutine bound_tsdc(icalcan,icalcen,kcalcan,kcalcen)

    ! calculates boundary values for temperature, salinity, frazil, density, and
    ! inflow tracers from neumann conditions. 

    implicit none

    ! local variables

    integer i,k,l,icalcan,icalcen,kcalcan,kcalcen


    ! southern boundary

    if(kcalcan.le.2) then

       do i = 1,m
          temp(i,1) = temp(i,2)                          
          tempa(i,1) = tempa(i,2)                          
          salt(i,1) = salt(i,2)                          
          salta(i,1) = salta(i,2)                          
          rhop(i,1) = rhop(i,2)                          
          tf(i,1) = tf(i,2)                          
       end do

       if (frazil) then
          do i = 1,m
             ctot(i,1) = ctot(i,2)                          
             do l = 1,nice
                c(i,1,l) = c(i,2,l)
                ca(i,1,l) = ca(i,2,l)
             end do
          end do
       end if

       if (intrace) then
          do i = 1,m
             do l = ninfmin,ninfmax
                intr(i,1,l) = intr(i,2,l)
                intra(i,1,l) = intra(i,2,l)
             end do
          end do
       end if

    end if

    ! northern boundary

    if(kcalcen.ge.n-1) then

       do i = 1,m
          temp(i,n) = temp(i,n-1)                          
          tempa(i,n) = tempa(i,n-1)                          
          salt(i,n) = salt(i,n-1)                          
          salta(i,n) = salta(i,n-1)                          
          rhop(i,n) = rhop(i,n-1)                          
          tf(i,n) = tf(i,n-1)                          
       end do

       if (frazil) then
          do i = 1,m
             ctot(i,n) = ctot(i,n-1)                          
             do l = 1,nice
                c(i,n,l) = c(i,n-1,l)
                ca(i,n,l) = ca(i,n-1,l)
             end do
          end do
       end if

       if (intrace) then
          do i = 1,m
             do l = ninfmin,ninfmax
                intr(i,n,l) = intr(i,n-1,l)
                intra(i,n,l) = intra(i,n-1,l)
             end do
          end do
       end if

    end if

    ! western boundary

    if (icalcan.le.2) then

       do k = 1,n
          temp(1,k) = temp(2,k)                          
          tempa(1,k) = tempa(2,k)                          
          salt(1,k) = salt(2,k)                          
          salta(1,k) = salta(2,k)                          
          rhop(1,k) = rhop(2,k)                          
          tf(1,k) = tf(2,k)                          
       end do

       if (frazil) then
          do k = 1,n
             ctot(1,k) = ctot(2,k)                          
             do l = 1,nice
                c(1,k,l) = c(2,k,l)
                ca(1,k,l) = ca(2,k,l)
             end do
          end do
       end if

       if (intrace) then
          do k = 1,n
             do l = ninfmin,ninfmax
                intr(1,k,l) = intr(2,k,l)
                intra(1,k,l) = intra(2,k,l)
             end do
          end do
       end if

    end if

    ! eastern boundary

    if (icalcen.ge.m-1) then

       do k = 1,n
          temp(m,k) = temp(m-1,k)                          
          tempa(m,k) = tempa(m-1,k)                          
          salt(m,k) = salt(m-1,k)                          
          salta(m,k) = salta(m-1,k)                          
          rhop(m,k) = rhop(m-1,k)                          
          tf(m,k) = tf(m-1,k)                          
       end do

       if (frazil) then
          do k = 1,n
             ctot(m,k) = ctot(m-1,k)                          
             do l = 1,nice
                c(m,k,l) = c(m-1,k,l)
                ca(m,k,l) = ca(m-1,k,l)
             end do
          end do
       end if

       if (intrace) then
          do k = 1,n
             do l = ninfmin,ninfmax
                intr(m,k,l) = intr(m-1,k,l)
                intra(m,k,l) = intra(m-1,k,l)
             end do
          end do
       end if

    end if

    return
  end subroutine bound_tsdc

end module plume_bounds
