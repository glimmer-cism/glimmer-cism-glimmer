module plume_plotascii

contains

subroutine plotascii(icalcan,icalcen,kcalcan,kcalcen,varoutrat)

! plot ascii output of a field

  use plume_global

  implicit none

! local variables

      character*1 a(m,n),deca(ldec),dry,land,wet,sep,inf

      integer i,k,jci,icalcan,icalcen,kcalcan,kcalcen,varoutrat
      integer minc,ninc,msize,nsize,izo,izu

      real(kind=kdp) tmp,xsig,depth,difu,difo,rhoatmp


! data definitions

      data dry,land,wet,sep,inf/'.','l','w','*','>'/
      data deca/'w','a','b','c','d','e','f','g','h','i','j','k', &
       'l','m','n','o','p','q','r','s','t','u','v','x','y','z'/

      xsig = 5.0d-1

! if plume is too large for screen then use coarser resolution
! (varible grid step in each direction chosen to fit 50/70 box

      minc = 1
      ninc = 1
      msize = icalcen - icalcan + 1
      nsize = kcalcen - kcalcan + 1
      if (msize.gt.50) minc = msize/50 + 1 
      if (nsize.gt.74) ninc = nsize/74 + 1 

! choose output resolution fixed solely on extent of plume in one direction

      if (varoutrat.eq.2) then
        minc = ninc
      end if

      if (varoutrat.eq.3) then
        ninc = minc
      end if

! calculate array of plume/land/ambient values

      do i = 1,m
        do k = 1,n
          a(i,k) = land
          if (jc(i,k).le.0) goto 10
          a(i,k) = dry
          if (jcd(i,k).le.0) goto 10
          a(i,k) = wet
          tmp = pdep(i,k)
          if (abs(tmp).lt. dcr) goto 10
          tmp = tmp + sign(xsig,tmp)
          jci = int(abs(tmp)/10.d0) + 1
          jci = max0(jci,1)
          jci = min0(jci,ldec)
          a(i,k) = deca(jci)
10        continue
        end do
      end do

! over-write plume thickness in inflow region

      do i = 1,m
        do k = 1,n
          if (depinf(i,k).gt.0.d0) a(i,k) = inf
        end do
      end do

! over-write plume thickness where it has separated

      do i = 1,m
        do k = 1,n
          depth = wcdep + gldep - ipos(i,k) - 5.0d-1*pdep(i,k)
          izo = int(depth/dzincr) + 1
          izu = izo + 1
          difu = dble(izo)*dzincr - depth
          difo = dzincr - difu
          rhoatmp = (difu*rhovf(izo)+difo*rhovf(izu))/dzincr
          if ((jcd(i,k).eq.1).and.((rhoatmp - rhop(i,k)).le.septol)) &
                                                       a(i,k) = sep
        end do
      end do

! draw horizontal scale

      write(*,*)' '
      write(11,*)' '

      write(*,50) (k/100,k = kcalcan,kcalcen,2*ninc)
      write(*,50) ((k - (k/100*100))/10,k = kcalcan,kcalcen,2*ninc)
      write(*,50) (k - (k/100*100) - (k - (k/100*100))/10*10 &
                                      ,k = kcalcan,kcalcen,2*ninc)
      write(11,50) (k/100,k = kcalcan,kcalcen,2*ninc)
      write(11,50) ((k - (k/100*100))/10,k = kcalcan,kcalcen,2*ninc)
      write(11,50) (k - (k/100*100) - (k - (k/100*100))/10*10 &
                                       ,k = kcalcan,kcalcen,2*ninc)

! draw plume and vertical scale

      write(*,*)' '
      write(11,*)' '

      do i = icalcan,icalcen,minc
        write(*,100) i,(a(i,k),k = kcalcan,kcalcen,ninc)
        write(11,100) i,(a(i,k),k = kcalcan,kcalcen,ninc)
      end do

      write(*,*)' '
      write(11,*)' '

  50  format(5x,37(i1,' '))
  100 format(1x,i3,1x,74a1)

      return
end subroutine

end module

