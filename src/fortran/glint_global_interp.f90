
module glint_global_interp

  implicit none

contains

  subroutine global_interp (idl,il,alon,jl,alat,a,mask,idlo,ilo,alono,jlo,alato,ao,masko,ier)
    !
    !          this subroutine does an area weighted average from one grid,
    !          on a spherical earth, to another.  logical masks may be assigned
    !          for each grid, and only those grid boxes which are masked true
    !          on both grids will be used.  a value of amm   will be assigned
    !          to all nodes of the new grid which are initially false or have
    !          no data from the old grid.  the new mask will also be changed to
    !          false where no data is available.
    !
    !          restrictions:  longitude must be the first dimension and it
    !                         be monotonically increasing (west to east).
    !
    !                         latitude must be the second dimension and it
    !                         must be monotonic.
    !
    !                         values for longitude and latitude must be in
    !                         degrees.
    !
    !                         arrays that wrap around must repeat longitudes
    !                         with a 360 degree increment.  it will be assumed
    !                         that values in the wrapped input and mask arrays
    !                         will also repeat (wrapped values in these arrays
    !                         will not be used).
    !
    !        input
    !
    !   integer   idl    first dimension of input a and mask.
    !   integer   il     number of grid boxes in longitude for a and mask.
    !   real      alon   longitude (deg) limits of grid boxes for a and mask.
    !   integer   jl     number of grid boxes in latitude for a and mask.
    !   real      alat   latitude (deg) limits of grid boxes for a and mask.
    !   real      a      array of input data.
    !   logical   mask   mask for input data (.false. to mask out data).
    !
    !        output
    !
    !   integer   idlo   first dimension of output ao and masko.
    !   integer   ilo    number of grid boxes in longitude for ao and masko.
    !   real      alono  longitude (deg) limits of grid boxes for ao and masko.
    !   integer   jlo    number of grid boxes in latitude for ao and masko.
    !   real      alato  latitude (deg) limits of grid boxes for ao and masko.
    !   real      ao     array of output data.
    !   logical   masko  mask for output data (.false. to mask out data).
    !   integer   ier    error indication:
    !                    (values may be summed for multiple errors)
    !                    0  no errors
    !                    1  input longitude dimension and/or length <=0.
    !                    2  output dimension and/or length <=0.
    !                    4  input latititude dimension <=0.
    !                    8  output latitude dimension <=0.
    !                   16  wrap-around on input longitude grid doesn't
    !                       repeat (+360).
    !                   32  wrap-around on output longitude grid doesn't
    !                       repeat (+360).
    !                   64  longitude of input is not monotonic increasing.
    !                  128  longitude of output is not monotonic increasing.
    !                  256  latitude of input is not monotonic.
    !                  512  latitude of output is not monotonic.
    !                 1024  input longitude wraps but doesn't repeat identically.
    !                 2048  output longitude wraps but doesn't repeat identically.
    !                   -1  output mask is changed.
    !                   -2  output mask contains all false values.

    ! --------------------------------------------------------
    ! Subroutine arguments
    ! --------------------------------------------------------

    integer                     :: idl,il,jl,idlo,ilo,jlo,ier
    real,   dimension(il+1)     :: alon
    real,   dimension(jl+1)     :: alat
    real,   dimension(ilo+1)    :: alono
    real,   dimension(jlo+1)    :: alato
    real,   dimension(idl,jl)   :: a
    real,   dimension(idlo,jlo) :: ao
    logical,dimension(idl,jl)   :: mask
    logical,dimension(idlo,jlo) :: masko

    ! --------------------------------------------------------
    ! Internal variables
    ! --------------------------------------------------------

    real :: api=3.1415926536
    real :: amm,almx,almn,sgn,al,dln,almxo,almno,dlno,amnlto
    real :: amxlto,amnlt,amxlt,amnlno,amxlno,amnln,amxln,wt,avg
    real :: slatmx,wlat,slatmn,slon,slonp,slonmx,slonmn,delon
    integer :: i,j,iil,iilo,j1,j2,jj,i1,i2,k,ii,iii,iip

    ! --------------------------------------------------------

    amm=50.
    ier=0

    ! Check that the sizes of the arrays given are sensible --

    if (idl.lt.il.or.il.le.0) ier=1
    if (idlo.lt.ilo.or.ilo.le.0) ier=ier+2
    if (jl.le.0)  ier=ier+4
    if (jlo.le.0) ier=ier+8
    if (ier.gt.0) return

    ! Check monotonic increasing input longitudes ------------

    do i=2,il
       if (alon(i).le.alon(i-1)) then
          ier=ier+64
          exit
       endif
    end do

    ! Check monotonic increasing output longitudes -----------

    do i=2,ilo 
       if (alono(i).le.alono(i-1)) then
          ier=ier+128
          exit
       endif
    end do

    ! Check monotonicity of input latitudes ------------------

    sgn=(alat(2)-alat(1))
    do j=2,jl
       if (sgn.lt.0.0) then
          if (alat(j)-alat(j-1).ge.0) then
             ier=ier+256
             exit
          endif
       else if (sgn.gt.0.0) then
          if (alat(j)-alat(j-1).le.0.0) then
             ier=ier+256
             exit
          endif
       else
          ier=ier+256
          exit
       endif
    end do

    ! Check monotonicity of output latitudes ------------------

    sgn=(alato(2)-alato(1))
    do j=2,jlo 
       if (sgn.lt.0.0) then
          if (alato(j)-alato(j-1).ge.0.0) then
             ier=ier+512
             exit
          endif
       else if (sgn.gt.0.0) then
          if (alato(j)-alato(j-1).le.0.0) then
             ier=ier+512
             exit
          endif
       else
          ier=ier+512
          exit
       endif
    end do

    ! Find wrap around of input grid, if it exists ------------

    iil=il
    almx=alon(1)
    almn=alon(1)

    do i=2,il+1
       almx=max(almx,alon(i))
       almn=min(almn,alon(i))
       al=abs(alon(i)-alon(1))-360.0
       if (abs(al).le.1.e-4) then
          iil=i-1
          exit
       else if (al.gt.0.0) then
          ier=ier+1024
          go to 12
       endif
    end do

    dln=0.0
    if (almn.lt.0.0) then
       dln=int(-almn/360.0+.001)*360.0
    else if (almn.gt.360.0) then
       dln=-int(almn/360.0+.001)*360.0
    endif
12  continue

    ! Find wrap around of output grid, if it exists -----------

    iilo=ilo
    almxo=alono(1)
    almno=alono(1)
    do i=2,ilo+1
       almxo=max(almxo,alono(i))
       almno=min(almno,alono(i))
       al=abs(alono(i)-alono(1))-360.0
       if (abs(al).le.1.e-4) then
          iilo=i-1
          exit
       else if (al.gt.0.0) then
          ier=ier+2048
          go to 15
       endif
    end do

    dlno=0.0
    if (almno.lt.0.0) then
       dlno=int(-almno/360.0+.001)*360.0
    else if (almno.gt.360.0) then
       dlno=-int(almno/360.0+.001)*360.0
    endif
15  continue

    ! Test for errors.  return if any --------------------------

    if (ier.ne.0) return

    ! The output grid needs to begin with or after the input grid.

    if (almno+dlno.lt.almn+dln) dlno=dlno+360.0

    do j=1,jlo ! loop 200 - over output latitudes
       ! find index limits in latitude to cover the new grid.
       j1=jl+1
       j2=0
       amnlto=min(alato(j),alato(j+1))
       amxlto=max(alato(j),alato(j+1))

       ! search for index limits in j.

       do jj=1,jl
          amnlt=min(alat(jj),alat(jj+1))
          amxlt=max(alat(jj),alat(jj+1))
          ! find jj limits
          if (amxlt.gt.amnlto.and.amnlt.lt.amxlto) then
             j1=min(jj,j1)
             j2=max(jj,j2)
          endif
       end do

       ! if input grid doesn't at least partially cover the
       ! output grid box, no values will be assigned.  mask out
       ! all values for the latitude.

       if (j2.lt.j1) then
          do i=1,iilo 
             ao(i,j)=amm
             if (masko(i,j)) ier=-1
             masko(i,j)=.false.
          end do
!          go to 200
          cycle
       endif

       do i=1,iilo ! loop 100

          ! no need to compute if it is masked out.
          if (.not.masko(i,j)) cycle !go to 100

          ! find index limits in longitude to cover the new grid.
          i1=3*il+1
          i2=0
          amnlno=min(alono(i),alono(i+1))+dlno
          amxlno=max(alono(i),alono(i+1))+dlno

          ! search for index limits in i.
          ! because of wrap around it is necessary to
          ! look through the data twice.
          ! the output grid longitudes have been adjusted
          ! (using dlno) such that the first longitude in
          ! the output grid is greater than the first
          ! longitude on the input grid.

          do k=0,1
             do ii=1,iil 
                amnln=min(alon(ii),alon(ii+1))+dln+k*360.0
                amxln=max(alon(ii),alon(ii+1))+dln+k*360.0
                ! find ii limits
                if (amxln.gt.amnlno.and.amnln.lt.amxlno) then
                   i1=min(ii+k*il,i1)
                   i2=max(ii+k*il,i2)
                endif
             end do
          end do

          ! if input grid doesn't partially cover the output
          ! grid box, no values will be assigned.  mask out
          ! the grid box.

          if (i2.lt.i1) then
             ao(i,j)=amm
             if (masko(i,j)) ier=-1
             masko(i,j)=.false.
!             go to 100
             cycle
          endif

          wt=0.0
          avg=0.0

          do jj=j1,j2
             slatmx=max(alat(jj),alat(jj+1))
             slatmn=min(alat(jj),alat(jj+1))
             wlat=max(sin(min(amxlto,slatmx)*api/180.)-sin(max(amnlto,slatmn)*api/180.),0.0)
             if (wlat.ne.0.0) then
                do iii=i1,i2
                   slon=dln
                   slonp=dln
                   if (iii.gt.iil) then
                      slon=slon+360.
                      slonp=slonp+360.
                   endif
                   ii=mod(iii-1,iil)+1
                   iip=ii+1
                   if (mask(ii,jj)) then
                      slon=slon+alon(ii)
                      slonp=slonp+alon(iip)
                      slonmx=max(slon,slonp)
                      slonmn=min(slon,slonp)
                      delon=max(min(amxlno,slonmx)-max(amnlno,slonmn),0.0)
                      wt=wt+wlat*delon
                      avg=avg+a(ii,jj)*wlat*delon
                   endif
                end do
             endif
          end do

          if (wt.gt.0.0) then
             ao(i,j)=avg/wt
          else
             ao(i,j)=amm
             if (masko(i,j)) ier=-1
             masko(i,j)=.false.
          endif
100       continue
       end do
200    continue
    end do

    ! Finish filling the output array from wrap-around.
    
    if (iilo.lt.ilo) then
       do j=1,jlo 
          do i=iilo+1,ilo
             ao(i,j)=ao(i-iilo,j)
             masko(i,j)=masko(i-iilo,j)
          end do
       end do
    endif

    ! Check if output masko is all false.

    if (all(masko==.false.)) ier=-2

  end subroutine global_interp

end module glint_global_interp
