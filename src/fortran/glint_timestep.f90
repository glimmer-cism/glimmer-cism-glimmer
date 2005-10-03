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
  use glint_constants
  private
  public glint_i_tstep

contains

  subroutine glint_i_tstep(time,instance,g_temp,g_temp_range, &
       g_precip,g_zonwind,g_merwind,g_humid,g_lwdown,g_swdown,g_airpress, &
       g_orog,g_orog_out,g_albedo,g_ice_frac,g_veg_frac,g_snowice_frac,g_snowveg_frac,&
       g_snow_depth,g_water_in,g_water_out,t_win,&
       t_wout,ice_vol,out_f,orogflag,mbal_skip,ice_tstep)

    !*FD Performs time-step of an ice model instance. Note that this 
    !*FD code will need to be altered to take account of the 
    !*FD energy-balance mass-balance model when it is completed.
    !*FD
    !*FD Note also that input quantities here are accumulated/average totals since the
    !*FD last call.
    use glide
    use glide_setup
    use glide_io
    use paramets
    use glint_io
    use glint_mbal_io
    use glint_climate
    use glint_routing
    use glimmer_log
    implicit none

    ! ------------------------------------------------------------------------  
    ! Arguments
    ! ------------------------------------------------------------------------  

    integer,                intent(in)   :: time         !*FD Current time in hours
    type(glint_instance), intent(inout)  :: instance     !*FD Model instance
    real(rk),dimension(:,:),intent(in)   :: g_temp       !*FD Global mean surface temperature field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_temp_range !*FD Global surface temperature half-range field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_precip     !*FD Global precip field total (mm)
    real(rk),dimension(:,:),intent(in)   :: g_zonwind    !*FD Global mean surface zonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_merwind    !*FD Global mean surface meridonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_humid      !*FD Global surface humidity (%)
    real(rk),dimension(:,:),intent(in)   :: g_lwdown     !*FD Global downwelling longwave (W/m^2)
    real(rk),dimension(:,:),intent(in)   :: g_swdown     !*FD Global downwelling shortwave (W/m^2)
    real(rk),dimension(:,:),intent(in)   :: g_airpress   !*FD Global surface air pressure (Pa)
    real(rk),dimension(:,:),intent(in)   :: g_orog       !*FD Input global orography (m)
    real(rk),dimension(:,:),intent(out)  :: g_orog_out   !*FD Output orography (m)
    real(rk),dimension(:,:),intent(out)  :: g_albedo     !*FD Output surface albedo 
    real(rk),dimension(:,:),intent(out)  :: g_ice_frac   !*FD Output ice fraction
    real(rk),dimension(:,:),intent(out)  :: g_veg_frac   !*FD Output veg fraction
    real(rk),dimension(:,:),intent(out)  :: g_snowice_frac !*FD Output snow-ice fraction
    real(rk),dimension(:,:),intent(out)  :: g_snowveg_frac !*FD Output snow-veg fraction
    real(rk),dimension(:,:),intent(out)  :: g_snow_depth !*FD Output snow depth (m)
    real(rk),dimension(:,:),intent(out)  :: g_water_in   !*FD Input water flux (mm)
    real(rk),dimension(:,:),intent(out)  :: g_water_out  !*FD Output water flux (mm)
    real(rk),               intent(out)  :: t_win        !*FD Total water input (kg)
    real(rk),               intent(out)  :: t_wout       !*FD Total water output (kg)
    real(rk),               intent(out)  :: ice_vol      !*FD Output ice volume (m$^3$)
    type(output_flags),     intent(in)   :: out_f        !*FD Flags to tell us whether to do output   
    logical,                intent(in)   :: orogflag     !*FD Set if we have new global orog
    logical,                intent(in)   :: mbal_skip    !*FD set if we are to skip mass-balance accumulation
    logical,                intent(out)  :: ice_tstep    !*FD Set if we have done an ice time step

    ! ------------------------------------------------------------------------  
    ! Internal variables
    ! ------------------------------------------------------------------------  

    real(rk),dimension(:,:),allocatable :: upscale_temp  ! temporary array for upscaling
    real(rk),dimension(:,:),allocatable :: routing_temp  ! temporary array for flow routing
    real(rk),dimension(:,:),allocatable :: accum_temp    ! temporary array for accumulation
    real(rk),dimension(:,:),allocatable :: ablat_temp    ! temporary array for ablation
    integer, dimension(:,:),allocatable :: fudge_mask    ! temporary array for fudging
    real(sp),dimension(:,:),allocatable :: thck_temp     ! temporary array for volume calcs
    real(rk) :: start_volume,end_volume,flux_fudge

    ! Assume we always need this, as it's too complicated to work out when we do and don't

    allocate(thck_temp(instance%proj%nx,instance%proj%ny))
    ice_tstep=.false.

    if (.not.mbal_skip) then

       ! Downscale input fields -------------------------------------------------

       call glint_downscaling(instance,g_temp,g_temp_range,g_precip,g_orog,g_zonwind,g_merwind, &
            g_humid,g_lwdown,g_swdown,g_airpress,orogflag)

       ! ------------------------------------------------------------------------  
       ! Sort out some local orography and remove bathymetry. This relies on the 
       ! point 1,1 being underwater. However, it's a better method than just 
       ! setting all points < 0.0 to zero
       ! ------------------------------------------------------------------------  

       call glide_get_usurf(instance%model,instance%local_orog)
       call glint_remove_bath(instance%local_orog,1,1)

       ! ------------------------------------------------------------------------  
       ! Adjust the surface temperatures using the lapse-rate, by reducing to
       ! sea-level and then back up to high-res orography
       ! ------------------------------------------------------------------------  

       call glint_lapserate(instance%artm,real(instance%global_orog,rk),real(-instance%lapse_rate,rk))
       call glint_lapserate(instance%artm,real(instance%local_orog,rk), real(instance%lapse_rate,rk))

       ! Process the precipitation field if necessary ---------------------------
       ! and convert from mm to m

       call glint_calc_precip(instance)

       ! Get ice thickness, if necessary ----------------------------------------

       if (instance%whichacab==3) then
          call glide_get_thk(instance%model,thck_temp)
       endif

       ! Do accumulation --------------------------------------------------------

       call glint_accumulate(instance%mbal_accum,instance%artm,instance%arng,instance%prcp, &
            instance%snowd,instance%siced,instance%xwind,instance%ywind, &
            instance%local_orog,real(thck_temp,rk),instance%humid,instance%swdown,instance%lwdown, &
            instance%airpress)

    end if

    ! Write output if necessary ----------------------------------------------

    call set_time(instance%model,real(time)*real(hours2years))

    ! Initialise water budget quantities to zero. These will be over-ridden if
    ! there's an ice-model time-step

    t_win=0.0       ; t_wout=0.0
    g_water_out=0.0 ; g_water_in=0.0

    ! ------------------------------------------------------------------------  
    ! ICE TIMESTEP begins HERE ***********************************************
    ! ------------------------------------------------------------------------  

    if (time-instance%last_timestep.ge.instance%ice_tstep) then

       ice_tstep=.true.
       instance%last_timestep=time

       ! Get the mass-balance ------------------------------------------------

       call glint_get_mbal(instance%mbal_accum,instance%artm,instance%prcp,instance%ablt, &
            instance%acab,instance%snowd,instance%siced)

       ! Calculate the initial ice volume (scaled) ------------------------------

       call glide_get_thk(instance%model,thck_temp)
       start_volume=sum(thck_temp)

       ! Constrain accumulation according to topography and domain edges -----

       call fix_acab(instance%ablt,instance%acab,instance%prcp,thck_temp,instance%local_orog)

       ! Put climate inputs in the appropriate places, with conversion ----------

       call glide_set_acab(instance%model,instance%acab/real(instance%ice_tstep*hours2years,sp))
       call glide_set_artm(instance%model,instance%artm)

       ! Adjust acab for output. This is done here otherwise roundoff error leaves a 
       ! residual thickness when it's taken away from the current ice

       where (instance%acab<-thck_temp)
          instance%acab=-thck_temp
       end where

       ! Do water budget accounting ---------------------------------------------

       if (out_f%water_out.or.out_f%total_wout.or.out_f%water_in .or.out_f%total_win) then
          allocate(accum_temp(instance%proj%nx,instance%proj%ny))
          allocate(ablat_temp(instance%proj%nx,instance%proj%ny))
          accum_temp=instance%prcp
          ablat_temp=instance%ablt
       endif

       ! ---------------------------------------------------------------------
       ! do the different parts of the glint timestep
       ! ---------------------------------------------------------------------

       call glide_tstep_p1(instance%model,real(time,rk)*hours2years)
       call glide_tstep_p2(instance%model)
       call glide_tstep_p3(instance%model)

       ! Calculate flux fudge factor --------------------------------------------

       if (out_f%water_out.or.out_f%total_wout.or.out_f%water_in .or.out_f%total_win) then

          allocate(fudge_mask(instance%proj%nx,instance%proj%ny))

          call glide_get_thk(instance%model,thck_temp)
          end_volume=sum(thck_temp)

          where (thck_temp>0.0)
             fudge_mask=1
          elsewhere
             fudge_mask=0
          endwhere

          flux_fudge=(start_volume+sum(accum_temp)-sum(ablat_temp)-end_volume)/sum(fudge_mask)

          ! Apply fudge_factor

          where(thck_temp>0.0)
             ablat_temp=ablat_temp+flux_fudge
          endwhere

       endif

       ! Upscale water flux fields ----------------------------------------------
       ! First water input (i.e. mass balance + ablation)

       if (out_f%water_in) then
          allocate(upscale_temp(instance%proj%nx,instance%proj%ny))

          where (thck_temp>0.0)
             upscale_temp=accum_temp
          elsewhere
             upscale_temp=0.0
          endwhere

          call mean_to_global(instance%ups,   &
               upscale_temp,   &
               g_water_in,     &
               instance%out_mask)
          deallocate(upscale_temp)
       endif

       ! Now water output (i.e. ablation) - and do routing

       if (out_f%water_out) then
          allocate(upscale_temp(instance%proj%nx,instance%proj%ny))
          allocate(routing_temp(instance%proj%nx,instance%proj%ny))

          where (thck_temp>0.0)
             upscale_temp=ablat_temp
          elsewhere
             upscale_temp=0.0
          endwhere

          call glide_get_usurf(instance%model,instance%local_orog)
          call flow_router(instance%local_orog, &
               upscale_temp, &
               routing_temp, &
               instance%out_mask, &
               instance%proj%dx, &
               instance%proj%dy)

          call mean_to_global(instance%ups,   &
               routing_temp,   &
               g_water_out,    &
               instance%out_mask)
          deallocate(upscale_temp,routing_temp)

       endif

       ! Sum water fluxes and convert if necessary ------------------------------

       if (out_f%total_win) then
          t_win  = sum(accum_temp)*instance%proj%dx*instance%proj%dy
       endif

       if (out_f%total_wout) then
          t_wout = sum(ablat_temp)*instance%proj%dx*instance%proj%dy
       endif
    else
       call glide_io_writeall(instance%model,instance%model)
    end if

    ! ------------------------------------------------------------------------ 
    ! Upscaling of output
    ! ------------------------------------------------------------------------ 

    ! We now upscale all fields at once...

    call get_i_upscaled_fields(instance,g_orog_out,g_albedo,g_ice_frac,g_veg_frac, &
         g_snowice_frac,g_snowveg_frac,g_snow_depth)

    ! Calculate ice volume ---------------------------------------------------

    if (out_f%ice_vol) then
       call glide_get_thk(instance%model,thck_temp)
       ice_vol=sum(thck_temp)*instance%proj%dx*instance%proj%dy
    endif

    ! Write glint data to file -----------------------------------------------

    call glint_io_writeall(instance,instance%model)
    call glint_mbal_io_writeall(instance%mbal_accum,instance%model)

    ! Tidy up ----------------------------------------------------------------

    if (allocated(accum_temp)) deallocate(accum_temp)
    if (allocated(ablat_temp)) deallocate(ablat_temp)
    deallocate(thck_temp)

  end subroutine glint_i_tstep

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_remove_bath(orog,x,y)

    !*FD Sets ocean areas to zero height, working recursively from
    !*FD a known ocean point.

    real(sp),dimension(:,:),intent(inout) :: orog !*FD Orography --- used for input and output
    integer,                intent(in)    :: x,y  !*FD Location of starting point (index)

    integer :: nx,ny

    nx=size(orog,1) ; ny=size(orog,2)

    if (orog(x,y).lt.0.0) orog(x,y)=0.0
    call glint_find_bath(orog,x,y,nx,ny)

  end subroutine glint_remove_bath

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  recursive subroutine glint_find_bath(orog,x,y,nx,ny)

    !*FD Recursive subroutine called by {\tt glimmer\_remove\_bath}.

    real(sp),dimension(:,:),intent(inout) :: orog  !*FD Orography --- used for input and output
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


