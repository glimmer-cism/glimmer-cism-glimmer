! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_timestep.f90 - part of the GLIMMER ice model       + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004 GLIMMER contributors - see COPYRIGHT file 
! for list of contributors.
!
! This program is free software; you can redistribute it and/or 
! modify it under the terms of the GNU General Public License as 
! published by the Free Software Foundation; either version 2 of 
! the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License 
! along with this program; if not, write to the Free Software 
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 
! 02111-1307 USA
!
! GLIMMER is maintained by:
!
! Ian Rutt
! School of Geographical Sciences
! University of Bristol
! University Road
! Bristol
! BS8 1SS
! UK
!
! email: <i.c.rutt@bristol.ac.uk> or <ian.rutt@physics.org>
!
! GLIMMER is hosted on NeSCForge:
!
! http://forge.nesc.ac.uk/projects/glimmer/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module glint_timestep
  !*FD timestep of a GLINT instance

  use glint_type
  private
  public glint_i_tstep

contains

  subroutine glint_i_tstep(time,instance,lats,lons,g_temp,g_temp_range, &
                          g_precip,g_zonwind,g_merwind,g_orog,g_orog_out,g_albedo,g_ice_frac,&
                          g_water_in,g_water_out,t_win,t_wout,ice_vol,out_f)

    !*FD Performs time-step of an ice model instance. Note that this 
    !*FD code will need to be altered to take account of the 
    !*FD energy-balance mass-balance model when it is completed.
    !*FD
    !*FD Note also that input quantities here are accumulated totals since the
    !*FD last call.
    use glide
    use glide_setup
    use glint_mbal
    use paramets
    use glint_io
    use glint_climate
    use glint_routing
    implicit none

    ! ------------------------------------------------------------------------  
    ! Arguments
    ! ------------------------------------------------------------------------  

    real(rk),               intent(in)   :: time         !*FD Current time in years
    type(glint_instance), intent(inout)  :: instance     !*FD Model instance
    real(rk),dimension(:),  intent(in)   :: lats         !*FD Latitudes of global grid points (degrees north)
    real(rk),dimension(:),  intent(in)   :: lons         !*FD Longitudes of global grid points (degrees east)
    real(rk),dimension(:,:),intent(in)   :: g_temp       !*FD Global mean surface temperature field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_temp_range !*FD Global surface temperature half-range field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_precip     !*FD Global precip field total (mm)
    real(rk),dimension(:,:),intent(in)   :: g_zonwind    !*FD Global mean surface zonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_merwind    !*FD Global mean surface meridonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_orog       !*FD Input global orography (m)
    real(rk),dimension(:,:),intent(out)  :: g_orog_out   !*FD Output orography (m)
    real(rk),dimension(:,:),intent(out)  :: g_albedo     !*FD Output surface albedo 
    real(rk),dimension(:,:),intent(out)  :: g_ice_frac   !*FD Output ice fraction
    real(rk),dimension(:,:),intent(out)  :: g_water_in   !*FD Input water flux (mm)
    real(rk),dimension(:,:),intent(out)  :: g_water_out  !*FD Output water flux (mm)
    real(rk),               intent(out)  :: t_win        !*FD Total water input (kg)
    real(rk),               intent(out)  :: t_wout       !*FD Total water output (kg)
    real(rk),               intent(out)  :: ice_vol      !*FD Output ice volume (m$^3$)
    type(output_flags),     intent(in)   :: out_f        !*FD Flags to tell us whether to do output       

    ! ------------------------------------------------------------------------  
    ! Internal variables
    ! ------------------------------------------------------------------------  

    real(rk),dimension(:,:),allocatable :: upscale_temp  ! temporary array for upscaling
    real(rk),dimension(:,:),allocatable :: routing_temp  ! temporary array for flow routing
    real(rk),dimension(:,:),allocatable :: accum_temp    ! temporary array for accumulation
    real(rk),dimension(:,:),allocatable :: ablat_temp    ! temporary array for ablation
    integer, dimension(:,:),allocatable :: fudge_mask    ! temporary array for fudging
    real(rk) :: f1 ! Scaling factor for converting precip and run-off amounts.
    character(40) :: timetxt
    real(rk) :: start_volume,end_volume,flux_fudge

    ! Set scaling factor for water flux calculations -------------------------

    f1 = scyr * thk0 / tim0

    ! ------------------------------------------------------------------------  
    ! Update internal clock
    ! ------------------------------------------------------------------------  

    instance%model%numerics%time=time  

    ! Set mass-balance accumulation start if necessary -----------------------

    if (instance%first_accum) then
      instance%accum_start=nint(time-instance%climate%tinc_mbal)
      instance%first_accum=.false.
    endif

    ! ------------------------------------------------------------------------  
    ! Downscale input fields, but only if required by options selected
    ! ------------------------------------------------------------------------  

    ! Temperature downscaling ------------------------------------------------

    select case(instance%climate%whichartm)
    case(7)
      call interp_to_local(instance%proj,  &
                           g_temp,         &
                           instance%downs, &
                           localsp=instance%climate%g_artm)

      call interp_to_local(instance%proj,  &
                           g_temp_range,   &
                           instance%downs, &
                           localsp=instance%climate%g_arng)
    end select

    ! Precip downscaling -----------------------------------------------------

    select case(instance%climate%whichprecip)
    case(1,2)
      call interp_to_local(instance%proj,  &
                           g_precip,       &
                           instance%downs, &
                           localsp=instance%climate%prcp)
    end select

    ! Orography downscaling --------------------------------------------------

    call interp_to_local(instance%proj,    &
                         g_orog,           &
                         instance%downs,   &
                         localdp=instance%global_orog)

    ! Wind downscaling -------------------------------------------------------

    select case(instance%climate%whichprecip)
    case(2)
      call interp_wind_to_local(instance%proj,    &
                                g_zonwind,        &
                                g_merwind,        &
                                instance%downs,   &
                                instance%xwind,instance%ywind)
    end select

    ! ------------------------------------------------------------------------  
    ! Sort out some local orography - scale orography by correct amount
    ! ------------------------------------------------------------------------  

    instance%local_orog=instance%model%geometry%usrf*thk0

    ! Remove bathymetry. This relies on the point 1,1 being underwater.
    ! However, it's a better method than just setting all points < 0.0 to zero

    call glint_remove_bath(instance%local_orog,1,1)

    ! ------------------------------------------------------------------------  
    ! Calculate or process the surface temperature and half-range
    ! ------------------------------------------------------------------------  

    call calcartm(instance,                       &
                  instance%climate%whichartm,     &
                  instance%model%geometry%usrf,         &
                  instance%model%climate%lati,          &
                  instance%model%climate%artm,          &
                  instance%climate%arng,          &
                  g_orog=instance%global_orog,          &
                  g_artm=instance%climate%g_artm, &
                  g_arng=instance%climate%g_arng)

    ! ------------------------------------------------------------------------  
    ! Calculate or process the precipitation field 
    ! ------------------------------------------------------------------------  

    call calcprcp(instance%climate%whichprecip, &
                  instance%climate%prcp, &
                  instance%climate%uprecip_rate, &
                  f1, &
                  instance%xwind, &
                  instance%ywind, &
                  instance%model%climate%artm, &
                  instance%local_orog, &
                  instance%proj%dx, &
                  instance%proj%dy, &
                  instance%climate%presprcp, &
                  instance%climate%presartm, &
                  instance%climate%pfac)

    ! ------------------------------------------------------------------------ 
    ! Calculate ablation, and thus mass-balance
    ! ------------------------------------------------------------------------ 

    call calcacab(instance%model%numerics, &
         instance%climate, &
         instance%pddcalc,  &
         instance%climate%whichacab, &
         instance%model%geometry% usrf,      &
         instance%model%climate%  artm,      &
         instance%climate%  arng,      &
         instance%climate%  prcp,      &
         instance%climate%  ablt,      &
         instance%model%climate%  lati,      &
         instance%model%climate%  acab,      &
         instance%model%geometry%  thck)

    ! Accumulate mass-balance if necessary

    instance%climate%prcp_save = instance%climate%prcp_save + instance%climate%prcp
    instance%climate%ablt_save = instance%climate%ablt_save + instance%climate%ablt
    instance%climate%acab_save = instance%climate%acab_save + instance%model%climate%acab

    ! If it's not time for a dynamics/temp/velocity timestep, return

    if (time-instance%accum_start.lt.instance%model%numerics%tinc) return

    ! ------------------------------------------------------------------------  
    ! ICE TIMESTEP begins HERE ***********************************************
    ! ------------------------------------------------------------------------  
 
    ! Copy and zero accumulated totals ---------------------------------------

    instance%climate%prcp = instance%climate%prcp_save
    instance%climate%ablt = instance%climate%ablt_save
    instance%model%climate%acab = instance%climate%acab_save

    instance%climate%prcp_save = 0.0
    instance%climate%ablt_save = 0.0
    instance%climate%acab_save = 0.0

    instance%accum_start=time

    ! Calculate the initial ice volume (scaled) ------------------------------

    start_volume=sum(instance%model%geometry%thck)

    ! Allocate temporary upscaling array -------------------------------------

    allocate(upscale_temp(instance%model%general%ewn,instance%model%general%nsn))
    allocate(routing_temp(instance%model%general%ewn,instance%model%general%nsn))
    allocate(accum_temp(instance%model%general%ewn,instance%model%general%nsn))
    allocate(ablat_temp(instance%model%general%ewn,instance%model%general%nsn))
    allocate(fudge_mask(instance%model%general%ewn,instance%model%general%nsn))

    accum_temp=0.0 ; ablat_temp=0.0

    ! Do water budget accounting ---------------------------------------------

    if (out_f%water_out.or.out_f%total_wout.or.out_f%water_in .or.out_f%total_win) then
      accum_temp=instance%climate%prcp*instance%model%numerics%dt
      ablat_temp=instance%climate%ablt*instance%model%numerics%dt
    endif


    ! ------------------------------------------------------------------------  
    ! If first step, use as seed...
    ! ------------------------------------------------------------------------  

    if (instance%first) then
      where (instance%local_orog>0.0)
        instance%model%geometry%thck = max(0.0d0, &
                                         instance%model%climate%acab* &
                                         instance%model%numerics%dt)
      endwhere

      call glide_calclsrf(instance%model%geometry%thck, &
                    instance%model%geometry%topg, &
                    instance%model%climate%eus, &
                    instance%model%geometry%lsrf)
      instance%model%geometry%usrf = instance%model%geometry%thck + &
                                     instance%model%geometry%lsrf

      ! Reset mass-balance variables to zero, as we've already used these

      instance%climate%prcp = 0.0
      instance%climate%ablt = 0.0
      instance%model%climate%acab = 0.0

    end if

    ! ------------------------------------------------------------------------
    ! do the first part of the glide time step
    ! ------------------------------------------------------------------------

    call glide_tstep_p1(instance%model,time)
    
    ! write glint data to file
    call glint_io_writeall(instance,instance%model)

    ! ------------------------------------------------------------------------
    ! do the second part of the glide time step
    ! ------------------------------------------------------------------------

    call glide_tstep_p2(instance%model)

    ! ------------------------------------------------------------------------
    ! do the final part of the glide time step (isostasy)
    ! ------------------------------------------------------------------------

    call glide_tstep_p3(instance%model)

    ! Calculate flux fudge factor --------------------------------------------

    if (out_f%water_out.or.out_f%total_wout.or.out_f%water_in .or.out_f%total_win) then

      end_volume=sum(instance%model%geometry%thck)

      where (instance%model%geometry%thck>0.0)
        fudge_mask=1
      elsewhere
        fudge_mask=0
      endwhere

      flux_fudge=(start_volume+sum(accum_temp)-sum(ablat_temp)-end_volume)/sum(fudge_mask)

      ! Apply fudge_factor

      where(instance%model%geometry%thck>0.0)
        ablat_temp=ablat_temp+flux_fudge
      endwhere

    endif

    ! ------------------------------------------------------------------------ 
    ! Upscaling of output
    ! ------------------------------------------------------------------------ 

    ! Upscale the output orography field, and re-dimensionalise --------------

    if (out_f%orog) call get_i_upscaled_fields(instance,orog=g_orog_out)

    ! Use thickness to calculate albedo and ice fraction ---------------------

    if (out_f%albedo.or.out_f%ice_frac.or.out_f%water_out) then
      call get_i_upscaled_fields(instance,albedo=g_albedo,ice_frac=g_ice_frac)
    endif

    ! Calculate ice volume ---------------------------------------------------

    if (out_f%ice_vol) then
      ice_vol=f1*1000.0*sum(instance%model%geometry%thck)*instance%proj%dx*instance%proj%dy
    endif

    ! Upscale water flux fields ----------------------------------------------
    ! First water input (i.e. mass balance + ablation)

    if (out_f%water_in) then

      where (instance%model%geometry%thck>0.0)
        upscale_temp=f1*1000.0*accum_temp
      elsewhere
        upscale_temp=0.0
      endwhere

      call mean_to_global(instance%proj,  &
                          instance%ups,   &
                          upscale_temp,   &
                          g_water_in,     &
                          instance%climate%out_mask)

    endif

    ! Now water output (i.e. ablation) - and do routing

    if (out_f%water_out) then

      where (instance%model%geometry%thck>0.0)
        upscale_temp=f1*1000.0*ablat_temp
      elsewhere
        upscale_temp=0.0
      endwhere

      call flow_router(instance%model%geometry%usrf, &
                       upscale_temp, &
                       routing_temp, &
                       instance%climate%out_mask, &
                       instance%proj%dx, &
                       instance%proj%dy)

      call mean_to_global(instance%proj,  &
                          instance%ups,   &
                          routing_temp,   &
                          g_water_out,    &
                          instance%climate%out_mask)

    endif

    ! Sum water fluxes and convert if necessary ------------------------------

    if (out_f%total_win) then
      t_win  = sum(accum_temp)*f1*1000.0*instance%proj%dx*instance%proj%dy
    endif

    if (out_f%total_wout) then
      t_wout = sum(ablat_temp)*f1*1000.0*instance%proj%dx*instance%proj%dy
    endif

    ! Tidy up ----------------------------------------------------------------

    deallocate(upscale_temp,accum_temp,ablat_temp,fudge_mask,routing_temp)
    instance%first=.false.
 
  end subroutine glint_i_tstep

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_i_upscaled_fields(instance,orog,albedo,ice_frac)

    !*FD Upscales and returns certain fields, according to the 
    !*FD arguments supplied

    use paramets

    type(glint_instance),            intent(in)  :: instance
    real(rk),dimension(:,:),optional,intent(out) :: orog
    real(rk),dimension(:,:),optional,intent(out) :: albedo
    real(rk),dimension(:,:),optional,intent(out) :: ice_frac

    real(rk),dimension(:,:),allocatable :: if_temp,upscale_temp

	  ! Calculate orography

    if (present(orog)) then
      call mean_to_global(instance%proj, &
                          instance%ups_orog, &
                          instance%model%geometry%usrf, &
                          orog,    &
                          instance%climate%out_mask)
      orog=thk0*orog
    endif

    if (present(albedo).or.present(ice_frac)) then

      if (present(albedo)) then
        allocate(if_temp(size(albedo,1),size(albedo,2)))
      else
        allocate(if_temp(size(ice_frac,1),size(ice_frac,2)))
      endif
      allocate(upscale_temp(instance%model%general%ewn,instance%model%general%nsn))

      ! First, ice coverage on local grid 
  
      where (instance%model%geometry%thck>0.0)
        upscale_temp=1.0
      elsewhere
        upscale_temp=0.0
      endwhere

      ! Upscale it...

      call mean_to_global(instance%proj, &
                          instance%ups, &
                          upscale_temp, &
                          if_temp,    &
                          instance%climate%out_mask)

      if (present(ice_frac)) ice_frac=if_temp

    endif

    ! Calculate albedo -------------------------------------------------------

    if (present(albedo)) then 
      where (if_temp>0.0)
        albedo=instance%climate%ice_albedo
      elsewhere
        albedo=0.0
      endwhere
    endif

    if (allocated(if_temp)) deallocate(if_temp)
    if (allocated(upscale_temp)) deallocate(upscale_temp)

  end subroutine get_i_upscaled_fields

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_remove_bath(orog,x,y)

    !*FD Sets ocean areas to zero height, working recursively from
    !*FD a known ocean point.

    real(rk),dimension(:,:),intent(inout) :: orog !*FD Orography --- used for input and output
    integer,                intent(in)    :: x,y  !*FD Location of starting point (index)

    integer :: nx,ny

    nx=size(orog,1) ; ny=size(orog,2)

    if (orog(x,y).lt.0.0) orog(x,y)=0.0
    call glint_find_bath(orog,x,y,nx,ny)

  end subroutine glint_remove_bath

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  recursive subroutine glint_find_bath(orog,x,y,nx,ny)

    !*FD Recursive subroutine called by {\tt glimmer\_remove\_bath}.

    real(rk),dimension(:,:),intent(inout) :: orog  !*FD Orography --- used for input and output
    integer,                intent(in)    :: x,y   !*FD Starting point
    integer,                intent(in)    :: nx,ny !*FD Size of array {\tt orography}

    integer,dimension(4) :: xi=(/ -1,1,0,0 /)
    integer,dimension(4) :: yi=(/ 0,0,-1,1 /)
    integer :: ns=4,i

    do i=1,ns
      if (x+xi(i).le.nx.and.x+xi(i).gt.0.and. &
          y+yi(i).le.ny.and.y+yi(i).gt.0) then
        if (orog(x+xi(i),y+yi(i)).lt.0.0) then
          orog(x+xi(i),y+yi(i))=0.0
          call glint_find_bath(orog,x+xi(i),y+yi(i),nx,ny)
        endif
      endif
    enddo

  end subroutine glint_find_bath

end module glint_timestep
  
  
