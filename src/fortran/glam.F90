
!***********************************************************************
module glam         
!***********************************************************************


    ! *sp* 1st-order ice sheet dynamics from Payne/Price solver and dH/dt from LANL incremental remapping
    ! ... analagous in form to "glissade" ?

    ! *sp* contains only stubs for now

    use glide_types

    use glam_strs2
    use remap_advection
    use remap_glamutils

    implicit none
    private

    public :: glam_driver

    contains

    ! *sp* note that initializtion routines for "glam_velo_fordsiapstr" and "remap_advection
    ! have been moved to initialization portion of "glide.F90"

    ! *sp* driver subroutine for Payne/Price HO dynamics and LANL inc. remapping for dH/dt
    ! ... called from 'glide'
    subroutine glam_driver( model )

        type(glide_global_type), intent(inout) :: model


        ! *sp* specify a few subroutine arguments here that are not already in the 
        ! model derived type (e.g. 'which*" flags). These are just dummy values for now
        ! to get things compiling ... 
        integer :: whichbabc, whichefvs, whichresid
        real (kind=dp), dimension(model%general%upn) :: sigma
        real (kind=dp), dimension(model%general%upn-1) :: stagsigma
        real (kind=dp), dimension(model%general%ewn-1,model%general%nsn-1) :: minTauf
        whichbabc = 9; whichefvs = 9; whichresid = 9
        sigma = 0.0d0; stagsigma = 0.0d0
        minTauf = 0.0d0

        ! *sp* also not consistent amoung the two models at present is "sigma" and "stagsigam". For now
        ! these will not be passed in. 'minTauf' avoided for now as well. 

        ! Note that the argument 'eta' was removed from the call, as it is not used. 

        ! *sp* note that the variables 'dlsrfdew', 'dlsrfdns' appear to be derived from their related
        ! usrf and thck derivs. rather than specified directly as in 'glam'. However, shouldn't they be
        ! e.g. 'dlsrfdew  = dusrfdew - dthckdew' rather than 'dlsrfdew = dthckdew - dusrfdew' ???   

        ! *sp* stub to 1st order solution 
        call glam_velo_fordsiapstr( model%general%ewn, model%general%nsn, model%general%upn,    &
                                    model%numerics%dew, model%numerics%dns,                     &
                                    sigma, stagsigma,                                           &                                    
                                    model%geometry%thck, model%geometry%usrf,                   &                                    
                                    model%geometry%lsrf, model%geometry%topg,                   &
                                    model%geomderv%dthckdew, model%geomderv%dthckdns,           &
                                    model%geomderv%dusrfdew, model%geomderv%dusrfdns,           &
                                    model%geomderv%dthckdew-model%geomderv%dusrfdew,            &
                                    model%geomderv%dthckdns-model%geomderv%dusrfdns,            & 
                                    model%geomderv%stagthck,                                    &
                                    minTauf,                                                    &
                                    model%temper%flwa,                                          &
                                    whichbabc, whichefvs, whichresid,                           &
                                    model%velocity_hom%uvel, model%velocity_hom%vvel,           &                                    
                                    model%velocity_hom%uflx, model%velocity_hom%vflx,           &                                                        
                                    model%velocity_hom%efvs, model%velocity_hom%tau,            &                                    
                                    model%velocity_hom%gdsx, model%velocity_hom%gdsy )

!   *sp* I think these are now redunant
!                                    ubas, vbas,                                                 &
!                                    tauxz, tauyz,                                               &
!                                    tauxx, tauyy,                                               &
!                                    tauxy,                                                      &


    ! *sp* This is the original subroutine signature from 'glam' code 
    !
    !    call glam_velo_fordsiapstr (ewn,      nsn,    upn,  &
    !                                dew,      dns,    &
    !                               sigma,    stagsigma,    &
    !                               thck,     usrf,         &
    !                               lsrf,     topg,         &
    !                               dthckdew, dthckdns,     &
    !                               dusrfdew, dusrfdns,     &
    !                               dlsrfdew, dlsrfdns,     &
    !                               stagthck, minTauf,      &
    !                               flwa,                   &
    !                               whichbabc,              &
    !                               whichefvs,              &
    !                               whichresid,             &
    !                               uvel,     vvel,         &
    !                               uflx,     vflx,         &
    !                               ubas,     vbas,         &
    !                               efvs,     tau,          &
    !                               tauxz,    tauyz,        &
    !                               tauxx,    tauyy,        &
    !                               tauxy,                  &
    !                               gdsx,     gdsy)


        ! *sp* stubs to remapping code (eventually replace w/ full calls below)

        ! *sp* put necessary variables in format for inc. remapping
        call horizontal_remap_in( ) 

        ! *sp* call remapping code
        call horizontal_remap( ) 

        ! *sp* put variables back into format to be used by glam
        call horizontal_remap_out( )


    !    call horizontal_remap_in( dt, thck(1:ewn-1,1:nsn-1), dew, dns, uflx, vflx, stagthck, &
    !              thck_ir, dew_ir, dns_ir, dewt_ir, dnst_ir, dewu_ir, dnsu_ir, hm_ir, tarea_ir, &
    !              ubar_ir, vbar_ir, trace_ir, dt_ir )
    !    call horizontal_remap( dt_ir, 2, ewn-1, nsn-1, ubar_ir, vbar_ir, thck_ir, trace_ir, dew_ir, &
    !              dns_ir, dewt_ir, dnst_ir, dewu_ir, dnsu_ir, hm_ir, tarea_ir )
    !    call horizontal_remap_out( thck_ir, thck, acab, dt )


!       These to be moved elsewhere ... somewhere in "glide_stop.F90"?
!
!        ! *sp* finalization routine, to be written
!        call glam_velo_fordsiapstr_final( )
!
!        ! *sp* finalization routine for remapping, exists in 'remap_utils' 
!        call horizontal_remap_final( )


    end subroutine glam_driver 


!***********************************************************************
end module glam
!***********************************************************************

