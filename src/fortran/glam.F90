
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

    use remap_advection
    use remap_glamutils

    use glide_velo_higher
    use glide_thck

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

        ! Compute the higher-order velocities using the method of Payne and Price

        !whl - to do - Make sure that the sigma field passed to glam is consistent with glam numerics.
        ! Note that the argument 'eta' was removed from the call, as it is not used. 

        ! *sfp** note that the variables 'dlsrfdew', 'dlsrfdns' appear to be derived from their related
        ! usrf and thck derivs. rather than specified directly as in 'glam'. However, shouldn't they be
        ! e.g. 'dlsrfdew  = dusrfdew - dthckdew' rather than 'dlsrfdew = dthckdew - dusrfdew' ???   


        ! *tjb** Moved the PP call to glide_velo_higher.  This needs to prompt some more rethinking
        !        regarding what's "glide", what's "glam", and what's "glissade" so that we have
        !        more logical separation between high-level modules
        call run_ho_diagnostic(model)

        ! *sfp** put necessary variables in format for inc. remapping

        call horizontal_remap_in(model%numerics%dt,       model%geometry%thck(1:ewn-1,1:nsn-1),  &
                                  ntrace_ir,               nghost_ir,                             &
                                  model%numerics%dew,      model%numerics%dns,                    &
                                  model%velocity_hom%uflx, model%velocity_hom%vflx,               &
                                  model%geomderv%stagthck, thck_ir,                      &
                                  dew_ir,                  dns_ir,                       &
                                  dewt_ir,                 dnst_ir,                      &
                                  dewu_ir,                 dnsu_ir,                      &
                                  hm_ir,                   tarea_ir,                     &
                                  ubar_ir,                 vbar_ir,                      &
                                  trace_ir,                dt_ir )

        ! *sfp** call remapping code

         call horizontal_remap  ( dt_ir,                                  & 
                                  ewn-1,               nsn-1,             &
                                  ntrace_ir,           nghost_ir,         &
                                  ubar_ir,             vbar_ir,           &
                                  thck_ir,             trace_ir,          &
                                  dew_ir,              dns_ir,            &
                                  dewt_ir,             dnst_ir,           &
                                  dewu_ir,             dnsu_ir,           &
                                  hm_ir,               tarea_ir )


        ! *sfp** put variables back into format to be used by glam

         call horizontal_remap_out (thck_ir,            mask_ir, model%geometry%thck,    &
                                    model%climate%acab, model%numerics%dt )

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

