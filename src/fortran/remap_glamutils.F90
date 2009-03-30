
module remap_glamutils      

  ! *sfp** contains various subroutines needed when using LANL incremental remapping code
  ! for thickness evolution in glam/glimmer codes

    use glimmer_paramets, only: sp, dp, len0, thk0, tim0, vel0

    ! *sfp** arrays needed to pass GLAM variables to/from inc. remapping solver
    real (kind = dp), allocatable, dimension(:,:,:) ::   &
          thck_ir,            &
          dew_ir,   dns_ir,   &
          dewt_ir,  dnst_ir,  &
          dewu_ir,  dnsu_ir,  &
          hm_ir,    tarea_ir, &
          ubar_ir,  vbar_ir
    real (kind = dp), allocatable, dimension(:,:,:,:) :: trace_ir
    real (kind = dp) :: dt_ir

    ! *sfp** mask to apply for domains where initial ice thickness limits are equivalent
    ! to the domain edge (e.g. the locations where bcs are applied)

    real (kind = dp), allocatable, dimension(:,:) :: mask_ir
    
    contains

!----------------------------------------------------------------------

    subroutine horizontal_remap_init (ewn, nsn)

    ! *sfp** initialize variables for use in inc. remapping code   

      implicit none

      integer, intent(in) :: ewn, nsn   ! horizontal dimensions

!whl - to do - Set ntrace to the actual number of tracers we want to transport
!              (e.g., ice temperature, ice age)
!              Hardwired ntrace = 1 for now

      integer, parameter  :: ntrace = 1    ! number of tracers

      dt_ir = 0.0_dp     ! time step

      ! allocate arrays/vars 
      allocate( thck_ir(1:ewn-1,1:nsn-1,1) ); thck_ir = 0.0_dp
      allocate( dew_ir(1:ewn-1,1:nsn-1,1) ); dew_ir = 0.0_dp
      allocate( dns_ir(1:ewn-1,1:nsn-1,1) ); dns_ir = 0.0_dp
      allocate( dewt_ir(1:ewn-1,1:nsn-1,1) ); dewt_ir = 0.0_dp
      allocate( dnst_ir(1:ewn-1,1:nsn-1,1) ); dnst_ir = 0.0_dp
      allocate( dewu_ir(1:ewn-1,1:nsn-1,1) ); dewu_ir = 0.0_dp
      allocate( dnsu_ir(1:ewn-1,1:nsn-1,1) ); dnsu_ir = 0.0_dp
      allocate( hm_ir(1:ewn-1,1:nsn-1,1) ); hm_ir = 0.0_dp
      allocate( tarea_ir(1:ewn-1,1:nsn-1,1) ); tarea_ir = 0.0_dp
      allocate( ubar_ir(1:ewn-1,1:nsn-1,1) ); ubar_ir = 0.0_dp
      allocate( vbar_ir(1:ewn-1,1:nsn-1,1) ); vbar_ir = 0.0_dp
      allocate( trace_ir(1:ewn-1,1:nsn-1,ntrace,1) ); trace_ir = 0.0_dp
      allocate( mask_ir(1:ewn,1:nsn) ); mask_ir = 0.0_dp

    end subroutine horizontal_remap_init

!----------------------------------------------------------------------

    subroutine horizontal_remap_final( )

    ! *sfp** deallocate variables for use in inc. remapping code   

      implicit none

      ! deallocate arrays 
      deallocate( thck_ir )
      deallocate( dew_ir ); deallocate( dns_ir )
      deallocate( dewt_ir ); deallocate( dnst_ir )
      deallocate( dewu_ir ); deallocate( dnsu_ir )
      deallocate( hm_ir )
      deallocate( tarea_ir )
      deallocate( ubar_ir ); deallocate( vbar_ir )
      deallocate( trace_ir )
      deallocate( mask_ir )

    end subroutine horizontal_remap_final

!----------------------------------------------------------------------

    subroutine horizontal_remap_in( dt,       thck,     &
                                    ntrace,   nghost,   &
                                    dew,      dns,      &
                                    uflx,     vflx,     &
                                    stagthck, thck_ir,  &
                                    dew_ir,   dns_ir,   &
                                    dewt_ir,  dnst_ir,  &
                                    dewu_ir,  dnsu_ir,  &
                                    hm_ir,    tarea_ir, &
                                    ubar_ir,  vbar_ir, &
                                    trace_ir, dt_ir )

    ! *sfp** get GLAM variables in order for use in inc. remapping code   

    implicit none

!whl - to do - add comments describing the arguments

    integer, intent(out) ::   &
         ntrace           ,&! number of tracer arrays to be remapped
         nghost             ! number of ghost cells

    real (kind = dp), dimension(:,:), intent(in) :: thck, uflx, vflx, stagthck
    real (kind = dp), intent(in) :: dew, dns, dt
    real (kind = dp), dimension(:,:,:), intent(out) ::   & 
         thck_ir,               &
         dew_ir,    dns_ir,     &
         dewt_ir,   dnst_ir,    &
         dewu_ir,   dnsu_ir,    &
         hm_ir,     tarea_ir,   &
         ubar_ir,   vbar_ir

    real (kind = dp), dimension(:,:,:,:), intent(out) :: trace_ir
    real (kind = dp), intent(out) :: dt_ir

!whl - Hardwire ntrace and nghost for now
!      Initially, no tracers are actually remapped, but the remapping routine
!       requires ntrace >= 1

    ntrace = 1

!whl - The number of ghost cells is applicable only when we have a parallel code.
!      For remapping, nghost = 2 is ideal because then no boundary updates are needed.
!      
    nghost = 2

    thck_ir(:,:,1) = thck(:,:)*thk0
    dew_ir(:,:,1)  = dew*len0; dns_ir(:,:,1) = dns*len0
    dewt_ir(:,:,1) = dew*len0; dnst_ir(:,:,1) = dns*len0
    dewu_ir(:,:,1) = dew*len0; dnsu_ir(:,:,1) = dns*len0
    hm_ir(:,:,1) = 1.0_dp
    tarea_ir = 1.0_dp / ( dew_ir * dns_ir )

    where( stagthck > 0.0_dp )
        ubar_ir(:,:,1) = uflx/stagthck*vel0;
        vbar_ir(:,:,1) = vflx/stagthck*vel0;
    elsewhere
        ubar_ir(:,:,1) = 0.0_dp
        vbar_ir(:,:,1) = 0.0_dp
    endwhere

!whl - to do - Fill the tracer array with ice temperature and other tracers
    trace_ir(:,:,:,1) = 1.0_dp;
    dt_ir = dt * tim0

    where( thck > 0.0_dp )
        mask_ir = 1.0_dp
    end where

    end subroutine horizontal_remap_in

!----------------------------------------------------------------------

    subroutine horizontal_remap_out( thck_ir, thck, acab, dt )

    ! *sfp** take output from inc. remapping and put back in GLAM format

    implicit none

    real (kind = dp), intent(in) :: dt
    real (kind = sp), intent(in), dimension(:,:) :: acab
    real (kind = dp), dimension(:,:,:), intent(in) :: thck_ir
    real (kind = dp), dimension(:,:), intent(inout) :: thck

    thck(:,:) = ( thck_ir(:,:,1) / thk0 + acab(:,:)*dt ) * mask_ir

    end subroutine horizontal_remap_out

!----------------------------------------------------------------------

end module remap_glamutils

