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
                                         thck,     usrf,         &
                                         lsrf,     topg,         &
                                         dthckdew, dthckdns,     &
                                         dusrfdew, dusrfdns,     & 
                                         dlsrfdew,               & 
                                         dlsrfdns,               &
                                         stagthck,               &
                                         flwa,                   & 
                                         whichbabc,              &
                                         whichefvs,              &
                                         whichresid,             &
!                                         uvel,      vvel )
                                         uvel,     vvel,         &
                                         uflx,     vflx,         &
                                         efvs,     tau,          &
                                         gdsx,     gdsy)

!   *sp* these need to be dealt w/ still
!                                         sigma,    stagsigma,    &
!                                         minTauf,      & 
!                                         ubas,     vbas,         &
!
!   *sp* I think these are now redundant
!                                         tauxz,    tauyz,        & 
!                                         tauxx,    tauyy,        & 
!                                         tauxy,                  & 

        implicit none

        integer, intent(in) :: ewn, nsn, upn
        real (kind = dp), intent(in) :: dew, dns

!        real (kind = dp), dimension(:),     intent(in)  :: sigma
        !whl - change to intent(in) when initialized elsewhere
!        real (kind = dp), dimension(:),     intent(inout) :: stagsigma

        real (kind = dp), dimension(:,:),   intent(in)  :: thck, usrf, lsrf, topg
        real (kind = dp), dimension(:,:),   intent(in)  :: dthckdew, dthckdns
        real (kind = dp), dimension(:,:),   intent(in)  :: dusrfdew, dusrfdns
        real (kind = dp), dimension(:,:),   intent(in)  :: dlsrfdew, dlsrfdns
        real (kind = dp), dimension(:,:),   intent(in)  :: stagthck
!        real (kind = dp), dimension(:,:),   intent(in)  :: minTauf
        real (kind = dp), dimension(:,:,:), intent(in)  :: flwa

        ! *sp* various flags
        integer, intent(in) :: whichbabc
        integer, intent(in) :: whichefvs
        integer, intent(in) :: whichresid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! *sp*
!
! NOTE that all of the following vars should have the intent(out) 
! attribute, but that this is not currently working due to them also
! having the pointer attribute. In lieu of a proper fix to this, I'm
! leaving the intent statements off for now, just so we can get the 
! stubs compiling properly and committed. Will have to come back to 
! this and figure out a fix later on.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        real (kind = dp), dimension(:,:,:), intent(out) :: uvel, vvel
        real (kind = dp), dimension(:,:,:) :: uvel, vvel

        ! *sp* NOTE that these are currently delcared elsewhere as 3d arrays, but they 
        ! SHOULD be 2d arrays (Tim will fix at some point)
!        real (kind = dp), dimension(:,:),   intent(out) :: uflx, vflx 
!        real (kind = dp), dimension(:,:,:),   intent(out) :: uflx, vflx 
        real (kind = dp), dimension(:,:,:) :: uflx, vflx 

        ! *sp* I think these relate to old code and don't need passing out anymore 
!        real (kind = dp), dimension(:,:),   intent(out) :: ubas, vbas

        ! *sp* Note that we need these as 3d, not 2d arrays (that is, they need to
        ! include the driving stress at various levels, not just at the sfc)
!        real (kind = dp), dimension(:,:,:), intent(out) :: gdsx, gdsy
!        real (kind = dp), dimension(:,:), intent(out) :: gdsx, gdsy
        real (kind = dp), dimension(:,:) :: gdsx, gdsy

!        real (kind = dp), dimension(:,:,:), intent(out) :: efvs
        real (kind = dp), dimension(:,:,:) :: efvs

        ! *sp* this will need to be fixed throughout, as we assume that 'tau' is a 3d array 
!        real (kind = dp), dimension(:,:,:), intent(out) :: tau
!        type(glide_tensor), intent(out) :: tau
        type(glide_tensor) :: tau

        ! *sp* These are calculated and passed out of the current version of the 1st order solver,
        ! but I'm guessing there is something analagous in existing CISM code, making these redundant
!        real (kind = dp), dimension(:,:,:), intent(out) :: tauxz, tauyz
!        real (kind = dp), dimension(:,:,:), intent(out) :: tauxx, tauyy
!        real (kind = dp), dimension(:,:,:), intent(out) :: tauxy

        end subroutine glam_velo_fordsiapstr


        subroutine glam_velo_fordsiapstr_final( )

        end subroutine glam_velo_fordsiapstr_final


!***********************************************************************
end module glam_strs2
!***********************************************************************
