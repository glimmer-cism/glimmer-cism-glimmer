
!***********************************************************************
module remap_glamutils      
!***********************************************************************

    ! *sp* contains various subroutines needed when using LANL incremental remapping code
    ! for thickness evolution in glam/glimmer codes

!    use glam_general, only : sp, dp
!    use glam_paramets, only : len0, thk0, tim0, vel0

    implicit none
    private

    public :: horizontal_remap_init, horizontal_remap_final, &
              horizontal_remap_in, horizontal_remap_out


!    ! *sp* arrays needed to pass GLAM variables to/from inc. remapping solver
!    real (kind = dp), allocatable, dimension(:,:,:) :: thck_ir, dew_ir, dns_ir, dewt_ir, dnst_ir, &
!                                                dewu_ir, dnsu_ir, hm_ir, tarea_ir, ubar_ir, vbar_ir
!    real (kind = dp), allocatable, dimension(:,:,:,:) :: trace_ir
!    real (kind = dp) :: dt_ir

!    ! *sp* mask to apply for domains where initial ice thickness limits are equivalent to the 
!    ! domain edge (e.g. the locations where bcs are applied)
!    real (kind = dp), allocatable, dimension(:,:) :: mask_ir


    contains


    ! *sp* initialize variables for use in inc. remapping code   
    subroutine horizontal_remap_init( ewn, nsn )

      implicit none

      integer, intent(in) :: ewn, nsn

!      dt_ir = 0.0_dp

!      ! allocate arrays/vars 
!      allocate( thck_ir(1:ewn-1,1:nsn-1,1) ); thck_ir = 0.0_dp
!      allocate( dew_ir(1:ewn-1,1:nsn-1,1) ); dew_ir = 0.0_dp
!      allocate( dns_ir(1:ewn-1,1:nsn-1,1) ); dns_ir = 0.0_dp
!      allocate( dewt_ir(1:ewn-1,1:nsn-1,1) ); dewt_ir = 0.0_dp
!      allocate( dnst_ir(1:ewn-1,1:nsn-1,1) ); dnst_ir = 0.0_dp
!      allocate( dewu_ir(1:ewn-1,1:nsn-1,1) ); dewu_ir = 0.0_dp
!      allocate( dnsu_ir(1:ewn-1,1:nsn-1,1) ); dnsu_ir = 0.0_dp
!      allocate( hm_ir(1:ewn-1,1:nsn-1,1) ); hm_ir = 0.0_dp
!      allocate( tarea_ir(1:ewn-1,1:nsn-1,1) ); tarea_ir = 0.0_dp
!      allocate( ubar_ir(1:ewn-1,1:nsn-1,1) ); ubar_ir = 0.0_dp
!      allocate( vbar_ir(1:ewn-1,1:nsn-1,1) ); vbar_ir = 0.0_dp
!      allocate( trace_ir(1:ewn-1,1:nsn-1,1,1) ); trace_ir = 0.0_dp
!      allocate( mask_ir(1:ewn,1:nsn) ); mask_ir = 0.0_dp

    end subroutine horizontal_remap_init


    ! *sp* get GLAM variables in order for use in inc. remapping code   
    subroutine horizontal_remap_in( ) 
!    subroutine horizontal_remap_in( dt, thck, dew, dns, uflx, vflx, stagthck, thck_ir, dew_ir, dns_ir, &
!                                dewt_ir, dnst_ir, dewu_ir, dnsu_ir, hm_ir, tarea_ir, ubar_ir, vbar_ir, &
!                                trace_ir, dt_ir )

        implicit none

!        real (kind = dp), dimension(:,:), intent(in) :: thck, uflx, vflx, stagthck
!        real (kind = dp), intent(in) :: dew, dns, dt
!        real (kind = dp), dimension(:,:,:), intent(out) :: thck_ir, dew_ir, dns_ir, dewt_ir, dnst_ir, &
!                                       dewu_ir, dnsu_ir, hm_ir, tarea_ir, ubar_ir, vbar_ir
!        real (kind = dp), dimension(:,:,:,:), intent(out) :: trace_ir
!        real (kind = dp), intent(out) :: dt_ir

!        thck_ir(:,:,1) = thck(:,:)*thk0
!        dew_ir(:,:,1) = dew*len0; dns_ir(:,:,1) = dns*len0
!        dewt_ir(:,:,1) = dew*len0; dnst_ir(:,:,1) = dns*len0
!        dewu_ir(:,:,1) = dew*len0; dnsu_ir(:,:,1) = dns*len0
!        hm_ir(:,:,1) = 1.0_dp
!        tarea_ir = 1.0_dp / ( dew_ir * dns_ir )

!        where( stagthck > 0.0_dp )
!            ubar_ir(:,:,1) = uflx/stagthck*vel0;
!            vbar_ir(:,:,1) = vflx/stagthck*vel0;
!        elsewhere
!            ubar_ir(:,:,1) = 0.0_dp
!            vbar_ir(:,:,1) = 0.0_dp
!        endwhere

!        trace_ir(:,:,1,1) = 1.0_dp;
!        dt_ir = dt * tim0

!        where( thck > 0.0_dp )
!            mask_ir = 1.0_dp
!        end where

    end subroutine horizontal_remap_in


    ! *sp* take output from inc. remapping and put back in GLAM format
    subroutine horizontal_remap_out( )
!    subroutine horizontal_remap_out( thck_ir, thck, acab, dt )

        implicit none

!        real (kind = dp), intent(in) :: dt
!        real (kind = sp), intent(in), dimension(:,:) :: acab
!        real (kind = dp), dimension(:,:,:), intent(in) :: thck_ir
!        real (kind = dp), dimension(:,:), intent(inout) :: thck

!        thck(:,:) = ( thck_ir(:,:,1) / thk0 + acab(:,:)*dt ) * mask_ir

    end subroutine horizontal_remap_out



    ! *sp* deallocate variables for use in inc. remapping code   
    subroutine horizontal_remap_final( )

      implicit none

!      ! deallocate arrays 
!      deallocate( thck_ir )
!      deallocate( dew_ir ); deallocate( dns_ir )
!      deallocate( dewt_ir ); deallocate( dnst_ir )
!      deallocate( dewu_ir ); deallocate( dnsu_ir )
!      deallocate( hm_ir )
!      deallocate( tarea_ir )
!      deallocate( ubar_ir ); deallocate( vbar_ir )
!      deallocate( trace_ir )
!      deallocate( mask_ir )

    end subroutine horizontal_remap_final



!***********************************************************************
end module remap_glamutils
!***********************************************************************

