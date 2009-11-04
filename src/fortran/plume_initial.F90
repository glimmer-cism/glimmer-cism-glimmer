module plume_initial

  use plume_global
  use plume_topog
  use plume_inflow
  use plume_ambient
  use plume_drag_set
  use plume_grid_set
  use plume_functions

contains

subroutine initial(iwetmin,iwetmax,kwetmin,kwetmax)
	  
  implicit none

! read data and set all initial fields

! local variables

      integer i,k,l,izo,izu,iwetmin,iwetmax,kwetmin,kwetmax

      real(kind=kdp) zd,difu,difo,tambz,sambz
      real(kind=kdp) depth,sambindep,tambindep,c1,c2,c3,rhoa


! initialise all arrays
! ---------------------

! m and n are defined in limits module

      do k = 1,n
        do i = 1,m
          u(i,k) = 0.d0
          v(i,k) = 0.d0
          ua(i,k) = 0.d0
          va(i,k) = 0.d0
          su(i,k) = 0.d0
          sv(i,k) = 0.d0
          u0(i,k) = 0.d0
          v0(i,k) = 0.d0
          u0a(i,k) = 0.d0
          v0a(i,k) = 0.d0
          tang(i,k) = 0.d0
          pdep(i,k) = 0.d0
          ipos(i,k) = 0.d0
          bpos(i,k) = 0.d0
          jc(i,k) = 0
          jcd(i,k) = 0
          jcd_u(i,k) = 0
          jcd_v(i,k) = 0
          jcd_fl(i,k) = 0       
          jcd_negdep(i,k) = 0
          jcd_fseed(i,k) = 0
          ctot(i,k) = 0.d0
          tf(i,k) = 0.d0
          entr(i,k) = 0.d0
          atemp(i,k) = 0.d0
          asalt(i,k) = 0.d0
          drag(i,k) = 0.d0
          bmelt(i,k) = 0.d0
          btemp(i,k) = 0.d0
          bsalt(i,k) = 0.d0
          ctempd(i,k) = 0.d0
! initialise ice to zero even when frazil is off so that densities are correct
          do l = 1,nice
            c(i,k,l) = 0.d0
            ca(i,k,l) = 0.d0
            fmelt(i,k,l) = 0.d0
            fppn(i,k,l) = 0.d0
            fnuc(i,k,l) = 0.d0
          end do
          tempinf(i,k) = 0.d0
          saltinf(i,k) = 0.d0
          depinf(i,k) = 0.d0
        end do
      end do

      if (intrace) then
        do k = 1,n
          do i = 1,m
            intrin(i,k) = 0
            do l = ninfmin,ninfmax
              intr(i,k,l) = 0.d0
            end do
          end do
        end do
      end if

! get cell dimensions
! -------------------

      call grid_set

! get topography and inflow regions, either from data or settings
! ---------------------------------------------------------------

      if (bathtype.gt.0) then
        call topog_inflow_set
      else
        call topog_read_edit
        call topog_smooth
        if (context.eq."fris")   call inflow_set_fris 
        if (context.eq."isomip") call inflow_set_isomip 
        if (context.eq."larsen") call inflow_set_larsen 
      end if

! set ambient properties
! ----------------------

      call ambient

! set drag coefficient
! --------------------

      do k = 1,n
        do i = 1,m
          drag(i,k) = cdb
        end do
      end do

      if (vardrag) call drag_set

! set inflow properties according to gade (jpo 1979) lines  
! --------------------------------------------------------
! (depth of plume mid-point so that t=tf initially if meltinf = 1)

      do k = 1,n
        do i = 1,m
          if (depinf(i,k).gt.0.d0) then

            depth = wcdep + gldep - bpos(i,k) + depinf(i,k)/2.d0

! interpolate ambient values 
            izo = int(depth/dzincr) + 1
            izu = izo + 1
            difu = dble(izo)*dzincr - depth 
            difo = dzincr - difu
            tambindep = (difu*tamb(izo) + difo*tamb(izu))/dzincr
            sambindep = (difu*samb(izo) + difo*samb(izu))/dzincr

! calculate gade values
            c1 = fta*(1.d0 - ci/c0)
            c2 = (ftb + ftc*depth)*(1.d0 - ci/c0) &
     &               - tambindep - lat/c0 + (ci/c0)*(fta*sambindep + ti)
            c3 = sambindep*(lat/c0 + (ci/c0)*(ftb + ftc*depth - ti))
            saltinf(i,k) = - (c2 + dsqrt(c2*c2 - 4.d0*c1*c3))/(2.d0*c1)
            tempinf(i,k) = freezing_temp_func(saltinf(i,k),depth) 

! calculate plume properties from proportions of ambient and melt waters

            saltinf(i,k) = sambindep &
     &                              + meltinf*(saltinf(i,k) - sambindep)
            tempinf(i,k) = tambindep &
     &                              + meltinf*(tempinf(i,k) - tambindep)

          end if
        end do
      end do

! initialise whole domain to wet and set plume initial thickness if mixed-layer
! model, but otherwise initialise wet area to inflow points only
! -----------------------------------------------------------------------------

      if (mixlayer) then

        iwetmin = 1
        iwetmax = m
        kwetmin = 1
        kwetmax = n

        do k = 1,n
          do i = 1,m
            if ((bpos(i,k).gt.0.d0).and.(bpos(i,k).lt.(wcdep+gldep))) &
     &                                            pdep(i,k) = depinit
          end do
        end do

      else

        iwetmin = m
        iwetmax = 1
        kwetmin = n
        kwetmax = 1

        do k = 1,n
          do i = 1,m
            if (depinf(i,k).gt.0.d0) then
              iwetmin = min(iwetmin,i)
              iwetmax = max(iwetmax,i)
              kwetmin = min(kwetmin,k)
              kwetmax = max(kwetmax,k)
            end if
          end do
        end do

      end if

! initialise flags, depths, and scalars
! -------------------------------------

      do k = 1,n
        do i = 1,m
          jc(i,k) = 1
          ipos(i,k) = bpos(i,k) - pdep(i,k)

! set index for solid ice points (no water column)
          if (bpos(i,k).le.dcr) then
            ipos(i,k) = bpos(i,k)
            jc(i,k) = 0
          endif
          pdep(i,k) = bpos(i,k) - ipos(i,k)

! set field for wet/dry points (jcd)
          if (pdep(i,k).ge.dcr) jcd(i,k) = 1

! set fields to ambient-fluid properties for depth
! (used when considering newly-wet cells)

          zd = dmax1(0.d0,wcdep + gldep - bpos(i,k))
          izo = int(zd/dzincr) + 1
          izu = izo + 1
          difu = dble(izo)*dzincr - zd
          difo = dzincr - difu
          tambz = (difu*tamb(izo) + difo*tamb(izu))/dzincr
          sambz = (difu*samb(izo) + difo*samb(izu))/dzincr
          rhoa =(difu*rhovf(izo)+difo*rhovf(izu))/dzincr

          tempa(i,k) = tambz                 
          temp(i,k) = tambz                 
          salta(i,k) = sambz                 
          salt(i,k) = sambz                 
          rhoamb(i,k) = rhoa
          rhop(i,k) = rhoa

        end do
      end do

      return
end subroutine

end module
