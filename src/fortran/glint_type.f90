! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_type.f90 - part of the GLIMMER ice model           + 
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

module glint_type

  !*FD contains type definitions for GLINT

  use glimmer_global
  use glint_proj
  use glint_interp
  use glint_mbal
  use glide_types

  implicit none

  type glint_instance

    !*FD Derived type holding information about ice model instance. 

    type(projection)                 :: proj               !*FD The projection definition of the instance.
    type(downscale)                  :: downs              !*FD Downscaling parameters.
    type(upscale)                    :: ups                !*FD Upscaling parameters
    type(upscale)                    :: ups_orog           !*FD Upscaling parameters for orography (to cope
                                                           !*FD with need to convert to spectral form).
    type(glide_global_type)          :: model              !*FD The instance and all its arrays.
    type(glint_mbal_params)          :: mbal_params        !*FD mass balance scheme parameters
    character(fname_length)          :: paramfile          !*FD The name of the configuration file.
    integer                          :: ice_tstep          !*FD Ice timestep in hours
    integer                          :: mbal_tstep         !*FD Mass-balance timestep in hours

    ! Climate inputs from global model --------------------------

    real(sp),dimension(:,:),pointer :: artm        => null() !*FD Annual mean air temperature
    real(sp),dimension(:,:),pointer :: arng        => null() !*FD Annual air temperature half-range
    real(sp),dimension(:,:),pointer :: prcp        => null() !*FD Precipitation (mm or m)
    real(sp),dimension(:,:),pointer :: snowd       => null() !*FD Snow depth (m)
    real(sp),dimension(:,:),pointer :: siced       => null() !*FD Superimposed ice depth (m)
    real(rk),dimension(:,:),pointer :: xwind       => null() !*FD $x$-component of surface winds (m/s)
    real(rk),dimension(:,:),pointer :: ywind       => null() !*FD $y$-component of surface winds (m/s)
    real(dp),dimension(:,:),pointer :: global_orog => null() !*FD Global orography (m)
    real(sp),dimension(:,:),pointer :: local_orog  => null() !*FD Local orography (m)
 
    ! Locally calculated climate/mass-balance fields ------------

    real(sp),dimension(:,:),pointer :: ablt => null() !*FD Annual ablation.
    real(sp),dimension(:,:),pointer :: acab => null() !*FD Annual mass-balance.

    ! Arrays to accumulate mass-balance quantities --------------

    real(sp),dimension(:,:),pointer :: prcp_save => null() !*FD used to accumulate precip
    real(sp),dimension(:,:),pointer :: ablt_save => null() !*FD used to accumulate ablation
    real(sp),dimension(:,:),pointer :: acab_save => null() !*FD used to accumulate mass-balance
    real(sp),dimension(:,:),pointer :: artm_save => null() !*FD used to average air-temperature

    ! Fractional coverage information ---------------------------
    
    real(rk) ,dimension(:,:),pointer :: frac_coverage => null() !*FD Fractional coverage of each 
                                                                !*FD global gridbox by the projected grid.
    real(rk) ,dimension(:,:),pointer :: frac_cov_orog => null() !*FD Fractional coverage of each 
                                                                !*FD global gridbox by the projected grid 
                                                                !*FD (orography).
    ! Output masking --------------------------------------------

    integer, dimension(:,:),pointer :: out_mask => null() 
    
    !*FD Array indicating whether a point should be considered or ignored 
    !*FD when upscaling data for output. 1 means use, 0 means ignore.

    ! Accumulation information ----------------------------------

    real(rk) :: accum_start = 0.0    !*FD Time when mass-balance accumulation started
    logical  :: first_accum = .true. !*FD First time accumulation flag

    ! Climate options -------------------------------------------

    integer :: whichacab = 1

    !*FD Which mass-balance scheme: 
    !*FD \begin{description}
    !*FD \item[1] PDD mass-balance model
    !*FD \item[2] Accumulation only 
    !*FD \item[3] RAPID energy balance model
    !*FD \end{description}

    integer :: whichprecip = 1

    !*FD Source of precipitation:
    !*FD \begin{description}
    !*FD \item[1] Use large-scale precip as is.
    !*FD \item[2] Use parameterization of \emph{Roe and Lindzen} 
    !*FD \end{description}

    ! Climate parameters ----------------------------------------------------------

    real(sp) :: ice_albedo   =   0.4 !*FD Ice albedo. (fraction)
    real(sp) :: lapse_rate   =   8.0 !*FD Uniform lapse rate in deg C/km 
                                     !*FD (N.B. This should be \emph{positive} for 
                                     !*FD temperature falling with height!)

    ! Counter for averaging temperature input --------------------------------------

    integer  :: av_count = 0 !*FD Counter for averaging temperature input

  end type glint_instance

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type output_flags
    !*FD A derived type used internally to communicate the outputs which need
    !*FD to be upscaled, thus avoiding unnecessary calculation

    logical :: orog         !*FD Set if we need to upscale the orography
    logical :: albedo       !*FD Set if we need to upscale the albedo
    logical :: ice_frac     !*FD Set if we need to upscale the ice fraction
    logical :: water_in     !*FD Set if we need to upscale the input water flux
    logical :: water_out    !*FD Set if we need to upscale the output water flux
    logical :: total_win    !*FD Set if we need to sum the total water taken up by ice sheet
    logical :: total_wout   !*FD Set if we need to sum the total ablation by the ice sheet
    logical :: ice_vol      !*FD Set if we need to calculate the total ice volume
  end type output_flags

contains

  subroutine glint_i_allocate(instance,nxg,nyg,nxgo,nygo)

    !*FD Allocate top-level arrays in
    !*FD the model instance, and ice model arrays.

    implicit none

    type(glint_instance),intent(inout) :: instance  !*FD Instance whose elements are to be allocated.
    integer,             intent(in)    :: nxg       !*FD Longitudinal size of global grid (grid-points).
    integer,             intent(in)    :: nyg       !*FD Latitudinal size of global grid (grid-points).
    integer,             intent(in)    :: nxgo       !*FD Longitudinal size of global orog grid (grid-points).
    integer,             intent(in)    :: nygo       !*FD Latitudinal size of global orog grid (grid-points).

    integer ewn,nsn

    ewn=get_ewn(instance%model)
    nsn=get_nsn(instance%model)

    ! First deallocate if necessary
    ! Downscaled global arrays

    if (associated(instance%artm))          deallocate(instance%artm)
    if (associated(instance%arng))          deallocate(instance%arng)
    if (associated(instance%prcp))          deallocate(instance%prcp)
    if (associated(instance%snowd))         deallocate(instance%snowd)
    if (associated(instance%siced))         deallocate(instance%siced)
    if (associated(instance%xwind))         deallocate(instance%xwind)
    if (associated(instance%ywind))         deallocate(instance%ywind)
    if (associated(instance%global_orog))   deallocate(instance%global_orog) 
    if (associated(instance%local_orog))    deallocate(instance%local_orog)

    ! Local climate arrays

    if (associated(instance%ablt))          deallocate(instance%ablt)
    if (associated(instance%acab))          deallocate(instance%acab)

    ! Accumulation arrays

    if (associated(instance%prcp_save))     deallocate(instance%prcp_save)
    if (associated(instance%ablt_save))     deallocate(instance%ablt_save)
    if (associated(instance%acab_save))     deallocate(instance%acab_save)
    if (associated(instance%artm_save))     deallocate(instance%artm_save)

    ! Fractional coverage

    if (associated(instance%frac_coverage)) deallocate(instance%frac_coverage)
    if (associated(instance%frac_cov_orog)) deallocate(instance%frac_cov_orog)

    ! Output mask

    if (associated(instance%out_mask))      deallocate(instance%out_mask)

    ! Then reallocate and zero...
    ! Global input fields

    allocate(instance%artm(ewn,nsn));        instance%artm = 0.0
    allocate(instance%arng(ewn,nsn));        instance%arng = 0.0
    allocate(instance%prcp(ewn,nsn));        instance%prcp        = 0.0
    allocate(instance%snowd(ewn,nsn));       instance%snowd       = 0.0
    allocate(instance%siced(ewn,nsn));       instance%siced       = 0.0
    allocate(instance%xwind(ewn,nsn));       instance%xwind       = 0.0
    allocate(instance%ywind(ewn,nsn));       instance%ywind       = 0.0
    allocate(instance%global_orog(ewn,nsn)); instance%global_orog = 0.0
    allocate(instance%local_orog(ewn,nsn));  instance%local_orog  = 0.0

    ! Local fields

    allocate(instance%ablt(ewn,nsn)); instance%ablt = 0.0
    allocate(instance%acab(ewn,nsn)); instance%acab = 0.0

    ! Accumulation arrays

    allocate(instance%prcp_save(ewn,nsn)); instance%prcp_save = 0.0
    allocate(instance%ablt_save(ewn,nsn)); instance%ablt_save = 0.0
    allocate(instance%acab_save(ewn,nsn)); instance%acab_save = 0.0
    allocate(instance%artm_save(ewn,nsn)); instance%artm_save = 0.0

    ! Fractional coverage map

    allocate(instance%frac_coverage(nxg,nyg)); instance%frac_coverage = 0.0
    allocate(instance%frac_cov_orog(nxgo,nygo)); instance%frac_cov_orog = 0.0

    ! Output mask

    allocate(instance%out_mask(ewn,nsn)); instance%out_mask = 1.0

  end subroutine glint_i_allocate

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_i_readconfig(instance,config)

    !*FD read glint configuration

    use glimmer_config

    implicit none

    ! Arguments

    type(ConfigSection), pointer       :: config      !*FD structure holding sections of configuration file   
    type(glint_instance),intent(inout) :: instance    !*FD The instance being initialised.

    ! Internals

    type(ConfigSection), pointer :: section

    ! GLINT projection parameters
    call proj_readconfig(instance%proj,config)         ! read glint projection configuration
    
    call GetSection(config,section,'GLINT climate')
    if (associated(section)) then
       call GetValue(section,'precip_mode',instance%whichprecip)
       call GetValue(section,'acab_mode',instance%whichacab)
       call GetValue(section,'ice_albedo',instance%ice_albedo)
       call GetValue(section,'lapse_rate',instance%lapse_rate)
    end if

  end subroutine glint_i_readconfig

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_i_printconfig(instance)

    use glimmer_log

    implicit none

    ! Argument

    type(glint_instance),intent(inout) :: instance    !*FD The instance to be printed

    ! Internal

    character(len=100) :: message

    call proj_printconfig(instance%proj)  ! Print projection info

    call write_log('GLINT climate')
    call write_log('-------------')
    write(message,*) 'precip_mode ',instance%whichprecip
    call write_log(message)
    write(message,*) 'acab_mode   ',instance%whichacab
    call write_log(message)
    write(message,*) 'ice_albedo  ',instance%ice_albedo
    call write_log(message)
    write(message,*) 'lapse_rate  ',instance%lapse_rate
    call write_log(message)
    call write_log('')

  end subroutine glint_i_printconfig

end module glint_type
