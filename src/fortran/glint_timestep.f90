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
    implicit none

    ! ------------------------------------------------------------------------  
    ! Arguments
    ! ------------------------------------------------------------------------  

    real(rk),               intent(in)   :: time         !*FD Current time in years
    type(glint_instance), intent(inout)  :: instance     !*FD Model instance
    real(rk),dimension(:),  intent(in)   :: lats         !*FD Latitudes of global grid points (degrees north)
    real(rk),dimension(:),  intent(in)   :: lons         !*FD Longitudes of global grid points (degrees east)
    real(rk),dimension(:,:),intent(in)   :: g_temp       !*FD Global mean surface temperature field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_temp_range !*FD Global surface temperature range field ($^{\circ}$C)
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
    real(rk),dimension(:,:),allocatable :: accum_temp    ! temporary array for accumulation
    real(rk),dimension(:,:),allocatable :: ablat_temp    ! temporary array for ablation
    real(rk) :: f1 ! Scaling factor for converting precip and run-off amounts.

    f1 = scyr * thk0 / tim0

    ! Allocate temporary upscaling array -------------------------------------

    allocate(upscale_temp(instance%model%general%ewn,instance%model%general%nsn))
    allocate(accum_temp(instance%model%general%ewn,instance%model%general%nsn))
    allocate(ablat_temp(instance%model%general%ewn,instance%model%general%nsn))

    ! ------------------------------------------------------------------------  
    ! Update internal clock
    ! ------------------------------------------------------------------------  

    instance%model%numerics%time=time  

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
    ! Calculate the precipitation field 
    ! ------------------------------------------------------------------------  

    select case(instance%climate%whichprecip)

    case(0)   ! Uniform precipitation ----------------------------------------

       instance%climate%prcp=instance%climate%uprecip_rate/f1

    case(1)   ! Large scale precip -------------------------------------------
       ! Note that we / by 1000 to convert to meters, and then by
       ! f1 for scaling in the model.

       instance%climate%prcp=instance%climate%prcp/(1000.0*f1)

    case(2)   ! Precip parameterization --------------------------------------

       call glint_precip(instance%climate%prcp, &
            instance%xwind,&
            instance%ywind,&
            instance%model%climate%artm,&
            instance%local_orog,&
            instance%proj%dx,&
            instance%proj%dy,&
            fixed_a=.true.)
       instance%climate%prcp=instance%climate%prcp/(1000.0*f1)

    case(3)   ! Prescribed small-scale precip from file, ---------------------
       ! adjusted for temperature

       instance%climate%prcp = instance%climate%presprcp * &
            pfac ** (instance%model%climate%artm - &
            instance%climate%presartm)

    end select

    ! ------------------------------------------------------------------------  
    ! If first step, use as seed...
    ! ------------------------------------------------------------------------  

    if (instance%first) then
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
            instance%model%climate%  acab)
       instance%model%geometry%thck = max(0.0d0, instance%model%climate%acab* instance%model%numerics%dt)

       call glide_calclsrf(instance%model%geometry%thck, &
            instance%model%geometry%topg, &
            instance%model%geometry%lsrf)
       instance%model%geometry%usrf = instance%model%geometry%thck + instance%model%geometry%lsrf
       !instance%first=.false.

       ! If this is the first time-setp, we need to make sure the
       ! accumulated ice is accounted for

       accum_temp=instance%model%geometry%thck

       ! Reset mass-balance variables to zero, as we've already used these

       instance%climate%prcp = 0.0
       instance%climate%ablt = 0.0
       instance%model%climate%acab = 0.0
    end if

    ! write glint data to file
    call glint_writeall(instance)

    ! ------------------------------------------------------------------------
    ! At this point we need to do the water flux output calculations, as the 
    ! limits of the ice-sheet may change before we get to the end of the
    ! subroutine, where the other outputs are calculated
    !
    ! Both these fields are also scaled to be in mm
    ! ------------------------------------------------------------------------
    ! Begin by calculating the actual mass balance (as opposed to the possible
    ! ablation/accumulation
    ! ------------------------------------------------------------------------

    if (out_f%water_out.or.out_f%total_wout.or. &
         out_f%water_in .or.out_f%total_win) then

       if (.not.instance%first) then
          accum_temp=instance%climate%prcp*instance%model%numerics%dt
          ablat_temp=min(instance%model%geometry%thck+accum_temp, &
               instance%climate%ablt*instance%model%numerics%dt)
       else
          ablat_temp=0.0
       endif

    endif

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

    ! Now water output (i.e. ablation)

    if (out_f%water_out) then

       where (instance%model%geometry%thck>0.0)
          upscale_temp=f1*1000.0*ablat_temp
       elsewhere
          upscale_temp=0.0
       endwhere

       call mean_to_global(instance%proj,  &
            instance%ups,   &
            upscale_temp,   &
            g_water_out,    &
            instance%climate%out_mask)

    endif

    ! ------------------------------------------------------------------------
    ! Sum water fluxes and convert if necessary
    ! ------------------------------------------------------------------------

    if (out_f%total_win) then
       t_win  = sum(accum_temp)*f1*1000.0*instance%proj%dx*instance%proj%dy
    endif

    if (out_f%total_wout) then
       t_wout = sum(ablat_temp)*f1*1000.0*instance%proj%dx*instance%proj%dy
    endif

    ! do the glide time step
    call glide_tstep(instance%model,time)

    ! ------------------------------------------------------------------------ 
    ! Upscaling of output
    ! ------------------------------------------------------------------------ 

    ! Upscale the output orography field, and re-dimensionalise --------------

    if (out_f%orog) then
       call mean_to_global(instance%proj, &
            instance%ups_orog, &
            instance%model%geometry%usrf, &
            g_orog_out,    &
            instance%climate%out_mask)
       g_orog_out=thk0*g_orog_out
    endif

    ! Use thickness to calculate albedo and ice fraction ---------------------

    if (out_f%albedo.or.out_f%ice_frac.or.out_f%water_out) then

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
            g_albedo,    &
            instance%climate%out_mask)

       ! Copy to ice fraction

       g_ice_frac=g_albedo

    endif

    ! Calculate albedo -------------------------------------------------------

    if (out_f%albedo) then 
       where (g_ice_frac>0.0)
          g_albedo=instance%climate%ice_albedo
       elsewhere
          g_albedo=0.0
       endwhere
    endif

    ! Calculate ice volume ---------------------------------------------------

    if (out_f%ice_vol) then
       ice_vol=f1*sum(instance%model%geometry%thck)*instance%proj%dx*instance%proj%dy
    endif

    ! Tidy up ----------------------------------------------------------------

    deallocate(upscale_temp,accum_temp,ablat_temp)
    instance%first=.false.

  end subroutine glint_i_tstep

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
  
  
