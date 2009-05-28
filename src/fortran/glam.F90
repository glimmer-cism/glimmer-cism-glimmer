
!***********************************************************************
module glam         
!***********************************************************************


    ! *sfp** 1st-order ice sheet dynamics from Payne/Price solver and dH/dt from LANL incremental remapping
    ! ... analagous in form to "glissade" ?

    ! *sfp** contains only stubs for now

    use glide_types
    use glimmer_paramets, only : vis0, vis0_glam 
    use glimmer_physcon, only :
    use glide_mask

    use glam_strs2, only: glam_velo_fordsiapstr, umask
    use remap_advection
    use remap_glamutils

    implicit none
    private

    public :: glam_driver

    ! *sfp** note that initializtion routines for "glam_velo_fordsiapstr" and "remap_advection
    ! have been moved to initialization portion of "glide.F90"

    ! *sfp** driver subroutine for Payne/Price HO dynamics and LANL inc. remapping for dH/dt
    ! ... called from 'glide'

    contains

    subroutine glam_driver( model )

        type(glide_global_type), intent(inout) :: model

        ! arrays passed in and out of remapping routine
        real (kind = dp), allocatable, dimension(:,:,:) ::   &
          thck_ir,            &
          dew_ir,   dns_ir,   &
          dewt_ir,  dnst_ir,  &
          dewu_ir,  dnsu_ir,  &
          hm_ir,    tarea_ir, &
          ubar_ir,  vbar_ir

        real (kind = dp), allocatable, dimension(:,:,:,:) :: trace_ir
        real (kind = dp) :: dt_ir
        integer :: ewn, nsn

        integer ::         & 
          ntrace_ir       ,&! number of tracers to be remapped
          nghost_ir         ! number of ghost cells

        ! *sfp** specify subroutine arguments here that are not already in the 
        ! model derived type. These are just dummy values for now
        ! to get things compiling ... 

        real (kind=dp), dimension(model%general%ewn-1,model%general%nsn-1) :: minTauf
        minTauf = 0.0d0

        ewn = model%general%ewn
        nsn = model%general%nsn

        ! *sfp* calculate mask for staggered thickness grid, using CISM subroutines. This will eventually
        ! take the place of the mask subroutines that are internal to 'glam_strs2'
        call glide_set_mask(model%numerics, model%geomderv%stagthck, model%geomderv%stagtopg, &
                             model%general%ewn, model%general%nsn, model%climate%eus, &
                             umask ) 



        ! Compute the higher-order velocities using the method of Payne and Price

        !whl - to do - Make sure that the sigma field passed to glam is consistent with glam numerics.
        ! Note that the argument 'eta' was removed from the call, as it is not used. 

        ! *sfp** note that the variables 'dlsrfdew', 'dlsrfdns' appear to be derived from their related
        ! usrf and thck derivs. rather than specified directly as in 'glam'. However, shouldn't they be
        ! e.g. 'dlsrfdew  = dusrfdew - dthckdew' rather than 'dlsrfdew = dthckdew - dusrfdew' ???   

        call glam_velo_fordsiapstr( model%general%ewn,       model%general%nsn,                 &
                                    model%general%upn,                                          &
                                    model%numerics%dew,      model%numerics%dns,                &
                                    model%numerics%sigma,    model%numerics%stagsigma,          &
                                    model%geometry%thck,     model%geometry%usrf,               &
                                    model%geometry%lsrf,     model%geometry%topg,               &
                                    model%geomderv%dthckdew, model%geomderv%dthckdns,           &
                                    model%geomderv%dusrfdew, model%geomderv%dusrfdns,           &
                                    model%geomderv%dusrfdew-model%geomderv%dthckdew,            &
                                    model%geomderv%dusrfdns-model%geomderv%dthckdns,            & 
                                    model%geomderv%stagthck, model%temper%flwa*vis0/vis0_glam,  &
                                    minTauf, umask,                                             &
                                    model%options%which_ho_babc,                                &
                                    model%options%which_ho_efvs,                                &
                                    model%options%which_ho_resid,                               &
                                    model%options%periodic_ew,                                  &
                                    model%options%periodic_ns,                                  &
                                    model%velocity_hom%beta,                                    & 
                                    model%velocity_hom%uvel, model%velocity_hom%vvel,           &
                                    model%velocity_hom%uflx, model%velocity_hom%vflx,           &
                                    model%velocity_hom%efvs )
                                    !model%velocity_hom%efvs,                                    & 
                                    !model%velocity_hom%kinematic_bc_u,                          &
                                    !model%velocity_hom%kinematic_bc_v )

        ! *sfp** put necessary variables in format for inc. remapping

!         call horizontal_remap_in(model%numerics%dt,       model%geometry%thck(1:ewn-1,1:nsn-1),  &
!                                  ntrace_ir,               nghost_ir,                             &
!                                  model%numerics%dew,      model%numerics%dns,                    &
!                                  model%velocity_hom%uflx, model%velocity_hom%vflx,               &
!                                  model%geomderv%stagthck, thck_ir,                      &
!                                  dew_ir,                  dns_ir,                       &
!                                  dewt_ir,                 dnst_ir,                      &
!                                  dewu_ir,                 dnsu_ir,                      &
!                                  hm_ir,                   tarea_ir,                     &
!                                  ubar_ir,                 vbar_ir,                      &
!                                  trace_ir,                dt_ir )

        ! *sfp** call remapping code

!         call horizontal_remap  ( dt_ir,                                  & 
!                                  ewn-1,               nsn-1,             &
!                                  ntrace_ir,           nghost_ir,         &
!                                  ubar_ir,             vbar_ir,           &
!                                  thck_ir,             trace_ir,          &
!                                  dew_ir,              dns_ir,            &
!                                  dewt_ir,             dnst_ir,           &
!                                  dewu_ir,             dnsu_ir,           &
!                                  hm_ir,               tarea_ir )


        ! *sfp** put variables back into format to be used by glam

!         call horizontal_remap_out (thck_ir,            model%geometry%thck,    &
!                                    model%climate%acab, model%numerics%dt )

!       These to be moved elsewhere ... somewhere in "glide_stop.F90"?
!
!        ! *sfp** finalization routine, to be written
!        call glam_velo_fordsiapstr_final( )
!
!        ! *sfp** finalization routine for remapping, exists in 'remap_utils' 
!        call horizontal_remap_final( )


    end subroutine glam_driver 


!***********************************************************************
end module glam
!***********************************************************************

