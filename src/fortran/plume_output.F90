module plume_output

contains

subroutine output(runtim,labtim, &
                    jobid,output_dir)                    
       
  ! produces surface output (single file version)
  
  use plume_global
  use plume_functions

  implicit none

! local variables

      character(len=4) outftim
      character(len=3) jobid
      character(len=*) output_dir

      integer i,k,l,izo,izu

      real(kind=kdp) runtim,labtim
      real(kind=kdp) odep,wtoi
      real(kind=kdp) rhopp(lm,ln),rhoap(lm,ln)
      real(kind=kdp) zc,difu,difo,tambz,sambz


! calculate origin depth

      odep=wcdep+gldep

! calculate water m/s to ice m/y conversion factor for output

      wtoi=31557600.d0*(rho0/rhoi)

! calculate potential density for ambient and plume (to give differences)

      do i = 1,m
        do k = 1,n
          zc =  wcdep + gldep - ipos(i,k) 
          izo = int(zc/dzincr) + 1
          izu = izo + 1
          difu = dble(izo)*dzincr - zc
          difo = dzincr - difu
          tambz = (difu*tamb(izo) + difo*tamb(izu))/dzincr
          sambz = (difu*samb(izo) + difo*samb(izu))/dzincr

          if (rholinear) then
            rhoap(i,k) = rho_func_linear(tambz,sambz)
            rhopp(i,k) = rho_func_linear(temp(i,k),salt(i,k))
          else
            rhoap(i,k) = rho_func_nonlinear(tambz,sambz,0.d0)
            rhopp(i,k) = rho_func_nonlinear(temp(i,k),salt(i,k),0.d0)
          end if

          if (frazil) then
            rhopp(i,k) = (1.d0 - ctot(i,k))*rhopp(i,k) + ctot(i,k)*rhoi
          end if

        end do
      end do

! ------------------------
! output in separate files
! ------------------------

! make file name in units of labtim

      write(outftim,1000) int(runtim/labtim)      

! write header file
! -----------------

      open(12,file=trim(output_dir)//jobid//'_'//outftim//'_header')

! write switches

      write(12,2000) frazil

! write dimensions

      write(12,2010) m,n,nice

! write grid

      do i = 1,m
        write(12,2020) dx(i)
      end do
      do k = 1,n
        write(12,2020) dy(k)
      end do

      close(12)

! write plume state file
! ----------------------

! note output formats differ from model variables:
!    draft and plume base position are output relative to sea level
!    all rates are output in m/y and all melt rates are of ice

      open(13,file=trim(output_dir)//jobid//'_'//outftim//'_states')

      do i = 1,m
        do k = 1,n
          write(13,3000) real(jcd(i,k)),su(i,k),sv(i,k),temp(i,k),salt(i,k), &
                        pdep(i,k),bpos(i,k)-odep,ipos(i,k)-odep,  &
                        rhoamb(i,k),rhop(i,k),rhoap(i,k),rhopp(i,k), &
                        entr(i,k)*31557600.d0,atemp(i,k),asalt(i,k), &
                        bmelt(i,k)*wtoi,btemp(i,k),bsalt(i,k), &
                        tf(i,k),depinf(i,k)
        end do
      end do

      close(13)

! write frazil file (if required)
! -------------------------------

! note output formats differ from model variables:
!    all ice interaction rates are output in m/y of ice

      if (frazil) then

        open(14,file=trim(output_dir)//jobid//'_'//outftim//'_frazil')

! write total fields

        do i = 1,m
          do k = 1,n
            write(14,4000) ctot(i,k),ctempd(i,k)
          end do
        end do

! write individual size-class fields

        do i = 1,m
          do k = 1,n
            do l = 1,nice
              write(14,4010) c(i,k,l),fmelt(i,k,l)*wtoi, &
                            fppn(i,k,l)*wtoi,fnuc(i,k,l)*wtoi
            end do
          end do
        end do

        close(14)

      end if

 1000 format('',i4.4)
 2000 format(l1)
 2010 format(3i4)
 2020 format(d13.7)
 3000 format(20(' ',d13.7))
 4000 format(2(' ',d13.7))
 4010 format(4(' ',d13.7))

      return
end  subroutine

end module


!
! the following code is deprecated but should not be deleted.
! it contains output to one large file in a specific format.
! this is readable by old graphics routines so needs to be kept
! for backward-compatibility, but also contains some extra output
! fields that could be added to the newer files in future.
!
!
!! ------------------------
!! output in one large file
!! ------------------------
!
!! open file (file name is in units of labtim)
!
!      write(outftim,1000) int(runtim/labtim)      
!      open(12,file=jobid//'_d_'//outftim)
!
!! write preamble
!
!      write(12,*) runtim,m,n,vardrag,frazil,intrace,tangle,negfrz,nice,
!     &            ninfmin,ninfmax,gldep+wcdep,cdb,iwetmin,iwetmax,
!     &                                      kwetmin,kwetmax,negdep
!
!! write grid
!
!      do i = 1,m
!        write(12,*) dx(i)
!      end do
!      do k = 1,n
!        write(12,*) dy(k)
!      end do
!
!! write general data fields
!
!      do i = 1,m
!        do k = 1,n
!          write(12,*) jcd(i,k),su(i,k),sv(i,k),temp(i,k),salt(i,k),
!     &                pdep(i,k),bpos(i,k),ipos(i,k), 
!     &                rhoamb(i,k),rhop(i,k),rhoap(i,k),rhopp(i,k),
!     &                entr(i,k),atemp(i,k),asalt(i,k),
!     &                bmelt(i,k),btemp(i,k),bsalt(i,k),
!     &                tf(i,k)
!        end do
!      end do
!
!! write variable drag field if included
!
!      if (vardrag) then
!        do i = 1,m
!          do k = 1,n
!            write(12,*) drag(i,k)
!          end do
!        end do
!      end if
!
!! write frazil fields if included
!
!      if (frazil) then
!        do i = 1,m
!          do k = 1,n
!            write(12,*) ctot(i,k),ctempd(i,k)
!            do l = 1,nice
!              write(12,*) c(i,k,l),fmelt(i,k,l),fppn(i,k,l),fnuc(i,k,l)
!            end do
!          end do
!        end do
!      end if
!
!! write inflow tracer fields if included
!
!      if (intrace) then
!        do i = 1,m
!          do k = 1,n
!            write(12,*) intrin(i,k)
!            do l = ninfmin,ninfmax
!              write(12,*) intr(i,k,l)
!            end do
!          end do
!        end do
!      end if
!
!! write turning angle fields if included
!
!      if (tangle) then
!        do i = 1,m
!          do k = 1,n
!            write(12,*) u0(i,k),v0(i,k),u0a(i,k),v0a(i,k),tang(i,k)
!          end do
!        end do
!      end if
!
!! close file
!
!      close(12)


