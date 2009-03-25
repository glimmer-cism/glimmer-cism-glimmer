!***********************************************************************
module glam_strs2
!***********************************************************************

    ! *sp* Payne/Price 1st order dynamics

    use glide_types 
    use glimmer_global, only : dp

    implicit none
    private

    public :: glam_velo_fordsiapstr_init, glam_velo_fordsiapstr, glam_velo_fordsiapstr_final

    contains

        subroutine glam_velo_fordsiapstr_init( )

            implicit none

        end subroutine glam_velo_fordsiapstr_init
   
 

        subroutine glam_velo_fordsiapstr(ewn,      nsn,    upn,  &
                                         dew,      dns,          &
                                         sigma,    stagsigma,    &
                                         thck,     usrf,         &
                                         lsrf,     topg,         &
                                         dthckdew, dthckdns,     &
                                         dusrfdew, dusrfdns,     & 
                                         dlsrfdew,               & 
                                         dlsrfdns,               &
                                         stagthck,               &
                                         minTauf,                & 
                                         flwa,                   & 
                                         whichbabc,              &
                                         whichefvs,              &
                                         whichresid,             &
                                         uvel,     vvel,         &
                                         uflx,     vflx,         &
                                         efvs,     tau,          &
                                         gdsx,     gdsy)

!   *sp* I think these are now redunant
!                                    ubas, vbas,                                                 &
!                                    tauxz, tauyz,                                               &
!                                    tauxx, tauyy,                                               &
!                                    tauxy,                                                      &

        implicit none

        integer, intent(in) :: ewn, nsn, upn
        real (kind = dp), intent(in) :: dew, dns

        real (kind = dp), dimension(:),     intent(in)  :: sigma
        !whl - change to intent(in) when initialized elsewhere
        real (kind = dp), dimension(:),     intent(inout) :: stagsigma

        real (kind = dp), dimension(:,:),   intent(in)  :: thck, usrf, lsrf, topg
        real (kind = dp), dimension(:,:),   intent(in)  :: dthckdew, dthckdns
        real (kind = dp), dimension(:,:),   intent(in)  :: dusrfdew, dusrfdns
        real (kind = dp), dimension(:,:),   intent(in)  :: dlsrfdew, dlsrfdns
        real (kind = dp), dimension(:,:),   intent(in)  :: stagthck
        real (kind = dp), dimension(:,:),   intent(in)  :: minTauf
        real (kind = dp), dimension(:,:,:), intent(in)  :: flwa

        ! *sp* various flags
        integer, intent(in) :: whichbabc
        integer, intent(in) :: whichefvs
        integer, intent(in) :: whichresid

        real (kind = dp), dimension(:,:,:), intent(out) :: uvel, vvel
        real (kind = dp), dimension(:,:),   intent(out) :: uflx, vflx 
        real (kind = dp), dimension(:,:,:), intent(out) :: efvs

        ! *sp* Note that we need these as 3d, not 2d arrays (that is, they need to
        ! include the driving stress at various levels, not just at the sfc)
        real (kind = dp), dimension(:,:,:), intent(out) :: gdsx, gdsy

        ! *sp* this will need to be fixed throughout, as we assume that 'tau' is a 3d array 
!        real (kind = dp), dimension(:,:,:), intent(out) :: tau
        type(glide_tensor), intent(out) :: tau

         !*sp* I think these are all redundant and/or handeled eleswhere in cism now
!        real (kind = dp), dimension(:,:),   intent(out) :: ubas, vbas
!        real (kind = dp), dimension(:,:,:), intent(out) :: tauxz, tauyz
!        real (kind = dp), dimension(:,:,:), intent(out) :: tauxx, tauyy
!        real (kind = dp), dimension(:,:,:), intent(out) :: tauxy

        end subroutine glam_velo_fordsiapstr


        subroutine glam_velo_fordsiapstr_final( )

        end subroutine glam_velo_fordsiapstr_final


!***********************************************************************
end module glam_strs2
!***********************************************************************
