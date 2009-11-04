
module plume

  ! 2d depth-averaged plume model from model of
  ! johann jungclaus & jan backhaus, 1994
  ! (journal of geophysical research, 99, 12375) 

  ! converted to model of an ice shelf water plume using formulation of
  ! adrian jenkins and andreas bombosch, 1995
  ! (journal of geophysical research, 100, 6967)
  ! and
  ! lars smedsrud and adrian jenkins, 2004
  ! (journal of geophysical research, 109, c03025)

  ! by paul holland, 2004 onwards.  basic reference is
  ! paul holland and daniel feltham, 2006
  ! (journal of physical oceanography, 36, 2312)

  ! model axes in code are 
  !   x/i/m - east (parallel to isobaths of ice shelf base, perp. to inflow)
  !   y/k/n - north (perp. to isobaths of ice shelf base, parallel to inflow)
  !   z     - upwards from the sea bed (defined as gldep+wcdep)

  ! set array sizes to maximum (data storage actually used set later)
  ! lamb = maximum number of vertical levels for ambient fluid t,s,rho 
  ! lice = maximum no of frazil size classes, linf = maximum no of inflows
  ! ldr = maximum number of regions in which drag coefficient is varied
  ! lsm = maximum number of regions in which smoothing is varied

  use plume_global

  use plume_inflow
  use plume_bounds
  use plume_continuity
  use plume_momentum
  use plume_frazil
  use plume_plotascii
  use plume_initial
  use plume_namelist
  use plume_output
  use plume_rho_calc
  use plume_scalar
  use plume_update

  implicit none


  ! local variables


  logical sepflag,ltmp
  character(len=3) jobid
  character(len=128) nl_filename
  character(len=128) output_dir

  integer,dimension(8) :: idt_values,dt_values
  integer :: count,count_rate,count_max
  real(kind=ksp)::sys_tim_a,sys_tim_b

  integer i,k,l,nsteps,istep
  integer iwetmin,iwetmax,kwetmin,kwetmax,ikcplus
  integer icalcan,kcalcan,icalcen,kcalcen,itmp
  integer ndays,nhours,nmin,nsec,varoutrat,seedtype

  real(kind=kdp) dt1,dt2,dtswtim
  real(kind=kdp) phi,tottim,outtim,fouttim,louttim,snottim,lnottim
  real(kind=kdp) radian,f,cosfdt,sinfdt,pressure,ttt
  real(kind=kdp) dtmp,runtim,labtim,totvol,negdep
  real(kind=kdp) frzcut(lice),systim



contains

  subroutine plume_initialize()



    ! *******************************************************************
    ! *** set parameters ************************************************
    ! *******************************************************************

    call namelist(nl_filename,varoutrat,seedtype,dt1,dt2,dtswtim,phi,labtim, &
         &              tottim,outtim,fouttim,louttim,snottim,lnottim)


    ! *******************************************************************
    ! *** preliminary operations ****************************************
    ! *******************************************************************

    ! open output file 
    open(11,file=trim(output_dir)//jobid//'_output')   


    ! calculate general derived quantities

    pi = 4.d0*datan(1.d0)
    nsteps = int(dtswtim/dt1 + (tottim - dtswtim)/dt2)
    radian = pi/180.d0
    f = 4.d0*pi*dsin(phi*radian)/86164.d0 ! coriolis parameter 
    dt = dt1
    gdt = g*dt1
    fdt = f*dt1
    cosfdt = 1.d0 ! used to be dcos(fdt)
    sinfdt = fdt  ! used to be dsin(fdt)

    ! set up frazil parameters

    if (frazil) call frazil_calc()

    ! initialise fields, including topography and grids

    call initial(iwetmin,iwetmax,kwetmin,kwetmax)

    ! initialise wet area to surround all inflow points

    ikcplus = 2  ! set additional area in which to check wetness
    icalcan = max0(2,iwetmin-ikcplus)
    kcalcan = max0(2,kwetmin-ikcplus)
    icalcen = min0(m-1,iwetmax+ikcplus)
    kcalcen = min0(n-1,kwetmax+ikcplus)

    ! initialise time, total of negative depths, and separation and negative
    ! frazil warning counters

    runtim = 0.d0
    negdep = 0.d0
    sepflag = .false.
    do l=1,nice
       frzcut(l) = 0.d0
    end do

    ! initialise elapsed time and set system time and write to output
    ! (these are gnu fortran routines and may not work on all compilers)

    systim = 0.d0

    call system_clock(count,count_rate,count_max)
    sys_tim_a = real(count) / count_rate

    call date_and_time(VALUES=idt_values)
    

    write(*,*) ' '
    write(11,*) ' '
    write(*,6000) idt_values(3),idt_values(2),idt_values(1),idt_values(5), &
                  idt_values(6),idt_values(7)
    write(11,6000) idt_values(3),idt_values(2),idt_values(1),idt_values(5),&
                   idt_values(6),idt_values(7)
    write(*,*) ' '
    write(11,*) ' '

    ! if restarting from an old run
    ! -----------------------------
    if (restart) then

       ! read old fields

       open(31,file = trim(output_dir)//jobid//'/start.dat')

       read(31,*) runtim,itmp,itmp,ltmp,ltmp,ltmp,ltmp,itmp, &
            &     itmp,itmp,dtmp,dtmp,iwetmin,iwetmax,kwetmin,kwetmax,negdep

       do i = 1,m
          read(31,*) dtmp
       end do
       do k = 1,n
          read(31,*) dtmp
       end do

       do i = 1,m
          do k = 1,n
             read(31,*) jcd(i,k),su(i,k),sv(i,k),temp(i,k),salt(i,k), &
                  pdep(i,k),bpos(i,k),ipos(i,k), &
                  rhoamb(i,k),rhop(i,k),dtmp,dtmp, &
                  entr(i,k),atemp(i,k),asalt(i,k), &
                  bmelt(i,k),btemp(i,k),bsalt(i,k), &
                  tf(i,k)
          end do
       end do

       if (vardrag) then
          do i = 1,m
             do k = 1,n
                read(31,*) drag(i,k)
             end do
          end do
       end if

       if (frazil) then
          do i = 1,m
             do k = 1,n
                read(31,*) ctot(i,k),dtmp
                do l = 1,nice
                   read(31,*) c(i,k,l),fmelt(i,k,l),fppn(i,k,l),fnuc(i,k,l)
                   ca(i,k,l) = c(i,k,l)
                end do
             end do
          end do
       end if

       if (intrace) then
          do i = 1,m
             do k = 1,n
                read(31,*) intrin(i,k)
                do l = ninfmin,ninfmax
                   read(31,*) intr(i,k,l)
                   intra(i,k,l) = intr(i,k,l)
                end do
             end do
          end do
       end if

       if (tangle) then
          do i = 1,m
             do k = 1,n
                read(31,*) u0(i,k),v0(i,k),u0a(i,k),v0a(i,k),tang(i,k)
             end do
          end do
       end if

       close(31)

       ! populate other arrays

       do i = 1,m
          do k = 1,n
             ua(i,k) = su(i,k)*5.0d-1*(pdep(i,k) + pdep(i+1,k))
             va(i,k) = sv(i,k)*5.0d-1*(pdep(i,k) + pdep(i,k+1))
             tempa(i,k) = temp(i,k)
             salta(i,k) = salt(i,k)
             ctota(i,k) = ctot(i,k)
          end do
       end do

       ! sort out timestep

       if (runtim.lt.dtswtim) then
          nsteps = int((dtswtim - runtim)/dt1 + (tottim - dtswtim)/dt2)
          dt = dt1
       else
          nsteps = int((tottim - runtim)/dt2)
          dt = dt2
       end if

    end if


6000 format(' start system date ',i2.2,'/',i2.2,'/',i4.4, &
         &                                 '; time ',i2.2,':',i2.2,':',i2.2)

  end subroutine plume_initialize


  subroutine plume_runstep()

    ! change timestep and recalculate basic quantities if necessary

    if (runtim.eq.dtswtim) then
       dt = dt2
       gdt = g*dt2
       fdt = f*dt2
       cosfdt = 1.d0 ! used to be dcos(fdt)
       sinfdt = fdt  ! used to be dsin(fdt)
    end if

    ! update time in seconds

    runtim = runtim + dt

    ! find indices of current wetted area

    icalcan = max0(2,iwetmin-ikcplus)
    kcalcan = max0(2,kwetmin-ikcplus)
    icalcen = min0(m-1,iwetmax+ikcplus)
    kcalcen = min0(n-1,kwetmax+ikcplus)

    ! --------------------
    ! calculate velocities
    ! --------------------

    call momentum(icalcan,kcalcan,icalcen,kcalcen,cosfdt,sinfdt,f)

    ! set velocity components on the open boundaries if the plume 
    ! has reached that far

    call bound_u_v(icalcan,icalcen,kcalcan,kcalcen)

    ! ---------------------------------------------------
    ! evaluate continuity equation for interface position
    ! ---------------------------------------------------

    call continuity(icalcan,kcalcan,icalcen,kcalcen)   

    ! interface and scalars passed forward on open boundary

    call inflow_calc(icalcan,kcalcan,icalcen,kcalcen)

    ! ----------------------------------------------------------
    ! update plume thickness, wet/dry boundaries, and velocities
    ! ----------------------------------------------------------

    call update(iwetmin,iwetmax,kwetmin,kwetmax,    &
         icalcan,icalcen,kcalcan,kcalcen,negdep)

    ! update interface position and flow on open boundaries

    call bound_interface(icalcan,icalcen,kcalcan,kcalcen)
    call bound_u_v(icalcan,icalcen,kcalcan,kcalcen)

    ! ---------------------------------------------------------------
    ! solve transport equations for temperature, salinity, and frazil
    ! ---------------------------------------------------------------

    call scalar(icalcan,kcalcan,icalcen,kcalcen,frzcut,seedtype)             

    ! -----------------------------------
    ! calculate density and update bounds
    ! -----------------------------------

    call rho_calc(icalcan,icalcen,kcalcan,kcalcen,sepflag)

    call bound_tsdc(icalcan,icalcen,kcalcan,kcalcen)

    ! write short step output and reset separation warning flag

    if (mod(runtim,snottim).eq.0) then
       sepflag = .false.
       ndays=int(runtim/86400.d0)
       nhours=int((runtim-ndays*86400.d0)/3600.d0)
       nmin=int((runtim-ndays*86400.d0-nhours*3600.d0)/60.d0)
       write(11,6100) ndays,nhours,nmin,int(runtim/dt)
       write(*,6100) ndays,nhours,nmin,int(runtim/dt)
    end if

    ! write long step output 

    if (mod(runtim,lnottim).eq.0) then

       ! write time and elapsed system time

       call date_and_time(VALUES=dt_values)
       write(*,*) ' '
       write(11,*) ' '
       write(*,6200) dt_values(3),dt_values(2),dt_values(1),dt_values(5),&
                     dt_values(6),dt_values(7)
       write(11,6200) dt_values(3),dt_values(2),dt_values(1),dt_values(5),&
                     dt_values(6),dt_values(7) 

       call system_clock(count,count_rate,count_max)
       sys_tim_b = real(count) / count_rate
       systim = systim + (sys_tim_b-sys_tim_a)
       sys_tim_a = sys_tim_b

       ndays=int(systim/86400.d0)
       nhours=int((systim-ndays*86400.d0)/3600.d0)
       nmin=int((systim-ndays*86400.d0-nhours*3600.d0)/60.d0)
       nsec=int(systim-ndays*86400.d0-nhours*3600.d0-nmin*60.d0)
       write(11,6300) ndays,nhours,nmin,nsec
       write(*,6300) ndays,nhours,nmin,nsec

       ! write totals of plume volume, negative depths, and negative frazil

       totvol = 0.d0
       do i = 1,m
          do k = 1,n
             totvol = totvol + pdep(i,k)*dx(i)*dy(k)
          end do
       end do
       write(11,*) 'total plume volume = ',totvol,' m^3'
       write(*,*) 'total plume volume = ',totvol,' m^3'

       if (negdep.ne.0.d0) then
          write(11,*) 'warning: negdep = ',negdep,' (cumulative)'  
          write(*,*) 'warning: negdep = ',negdep,' (cumulative)'
       end if

       do l=1,nice
          if (frzcut(l).gt.0.d0) then
             write(11,6350) ' warning: frzcut(',l,') = ',frzcut(l), &
                  ' (cumulative)'
             write(*,6350) ' warning: frzcut(',l,') = ',frzcut(l), &
                  ' (cumulative)'
          end if
       end do

       ! write ascii rendering of plume thickness
       !    
       call plotascii(icalcan,icalcen,kcalcan,kcalcen,varoutrat)

    end if

    ! write surface output to file

    if ((mod(runtim,outtim).eq.0).and. &
         ((runtim.ge.fouttim).and.(runtim.le.louttim))) then

       ndays=int(runtim/86400.d0)
       nhours=int((runtim-ndays*86400.d0)/3600.d0)
       nmin=int((runtim-ndays*86400.d0-nhours*3600.d0)/60.d0)

       call output(runtim,labtim,jobid,output_dir)

       write(11,6400) ndays,nhours,nmin,int(runtim/dt)
       write(*,6400) ndays,nhours,nmin,int(runtim/dt)
       write(*,*) ' '
       write(11,*) ' '

    end if

6100 format(' calculated ',i4,' days, ',i2,' hours and ',i2, &
         &                                         ' minutes (step ',i7,')')
6200 format(' current system date ',i2.2,'/',i2.2,'/',i4.4, &
         &                                 '; time ',i2.2,':',i2.2,':',i2.2)
6300 format(' elapsed system time ',i4,' days ',i2,' hours ',i2, &
         &                                        ' minutes ',i2,' seconds')
6350 format(a,i3,a,d10.3,a)
6400 format(' file output at ',i4,' days, ',i2,' hours and ',i2, &
         &                                         ' minutes (step ',i7,')')

  end subroutine plume_runstep

  subroutine plume_finalize()

    ! if program executed whole loop, report positively

    call system_clock(count,count_rate,count_max)
    sys_tim_b = real(count) / count_rate
    systim = systim + (sys_tim_b -  sys_tim_a)
    sys_tim_a = sys_tim_b

    ndays=int(systim/86400.d0)
    nhours=int((systim-ndays*86400.d0)/3600.d0)
    nmin=int((systim-ndays*86400.d0-nhours*3600.d0)/60.d0)
    nsec=int(systim-ndays*86400.d0-nhours*3600.d0-nmin*60.d0)
    write(11,6500) ndays,nhours,nmin
    write(*,6500) ndays,nhours,nmin
    write(*,*) ' '
    write(11,*) ' '

    ! *******************************************************************
    ! *** end of main loop **********************************************
    ! *******************************************************************

    ! close output file

    close(11)

6500 format(' :-) total system time ',i4,' days ',i2,' hours ',i2, &
         &                                        ' minutes ',i2,' seconds')

  end subroutine plume_finalize

  !subroutine plume_writeout()

  !end subroutine


end module plume

