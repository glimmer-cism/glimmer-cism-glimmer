
module remap_glamutils      

  ! *sfp** contains various subroutines needed when using LANL incremental remapping code
  ! for thickness evolution in glam/glimmer codes

    use glimmer_paramets, only: sp, dp, len0, thk0, tim0, vel0
    use glide_grids,      only: periodic_boundaries, periodic_boundaries_3d
    use xls

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

    integer :: ewn_ir, nsn_ir
    
    contains

!----------------------------------------------------------------------

    subroutine horizontal_remap_init (ewn, nsn, periodic_ew, periodic_ns)

    ! *sfp** initialize variables for use in inc. remapping code   

      implicit none

      integer, intent(in) :: ewn, nsn   ! horizontal dimensions
      logical, intent(in) :: periodic_ew, periodic_ns

!whl - to do - Set ntrace to the actual number of tracers we want to transport
!              (e.g., ice temperature, ice age)
!              Hardwired ntrace = 1 for now

      integer, parameter  :: ntrace = 1    ! number of tracers

      dt_ir = 0.0_dp     ! time step

      if (periodic_ew) then
        ewn_ir = ewn + 2
     else
        ewn_ir = ewn + 4
      end if

      if (periodic_ns) then
        nsn_ir = nsn + 2
      else
        nsn_ir = nsn + 4
      end if

      ! allocate arrays/vars 
      allocate( thck_ir(1:ewn_ir,1:nsn_ir,1) ); thck_ir = 0.0_dp
      allocate( dew_ir(1:ewn_ir,1:nsn_ir,1) ); dew_ir = 0.0_dp
      allocate( dns_ir(1:ewn_ir,1:nsn_ir,1) ); dns_ir = 0.0_dp
      allocate( dewt_ir(1:ewn_ir,1:nsn_ir,1) ); dewt_ir = 0.0_dp
      allocate( dnst_ir(1:ewn_ir,1:nsn_ir,1) ); dnst_ir = 0.0_dp
      allocate( dewu_ir(1:ewn_ir,1:nsn_ir,1) ); dewu_ir = 0.0_dp
      allocate( dnsu_ir(1:ewn_ir,1:nsn_ir,1) ); dnsu_ir = 0.0_dp
      allocate( hm_ir(1:ewn_ir,1:nsn_ir,1) ); hm_ir = 0.0_dp
      allocate( tarea_ir(1:ewn_ir,1:nsn_ir,1) ); tarea_ir = 0.0_dp
      allocate( ubar_ir(1:ewn_ir,1:nsn_ir,1) ); ubar_ir = 0.0_dp
      allocate( vbar_ir(1:ewn_ir,1:nsn_ir,1) ); vbar_ir = 0.0_dp
      allocate( trace_ir(1:ewn_ir,1:nsn_ir,ntrace,1) ); trace_ir = 0.0_dp
      allocate( mask_ir(1:ewn_ir,1:nsn_ir) ); mask_ir = 0.0_dp

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
                                    trace_ir, dt_ir, &
                                    periodic_ew, periodic_ns)

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

    logical, intent(in) :: periodic_ew, periodic_ns

    integer  :: ewn,  nsn
    integer  :: ngew, ngns !Number of *extra* ghost cells (depending on whether periodic bcs are enabled) in ew and ns


!whl - Hardwire ntrace and nghost for now
!      Initially, no tracers are actually remapped, but the remapping routine
!       requires ntrace >= 1

    ntrace = 1

!whl - The number of ghost cells is applicable only when we have a parallel code.
!      For remapping, nghost = 2 is ideal because then no boundary updates are needed.
!      
    nghost = 2

    ewn = size(thck, 1)
    nsn = size(thck, 2)

    ngew = (ewn_ir - ewn)/2
    ngns = (nsn_ir - nsn)/2
    
    thck_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = thck(:,:)*thk0
    dew_ir(:,:,1)  = dew*len0; dns_ir(:,:,1) = dns*len0
    dewt_ir(:,:,1) = dew*len0; dnst_ir(:,:,1) = dns*len0
    dewu_ir(:,:,1) = dew*len0; dnsu_ir(:,:,1) = dns*len0
    hm_ir(:,:,1) = 1.0_dp
    tarea_ir = 1.0_dp / ( dew_ir * dns_ir )

    call write_xls("uflx.txt", uflx)
    call write_xls("vflx.txt", vflx)

    !Copy fluxes over
    !*tjb** - Some rejiggering is needed here, because, while IR and Glide both place
    !velocities on a B-grid, the IR B-grid is the same size as the A-grid
    !(that is, there is an extra row of cells on the extreme right and bottom
    !of the domain).
    !If periodic BCs are used, then the fluxes already have a row of ghost cells
    !on the left and right.  This means that, like thickness, a second row
    !is needed on the left.  However, unlike thickness, *two* extra rows are needed
    !on the right, to account for the extra B-grid row.  Same goes for top and bottom.
    where( stagthck > 0.0_dp )
        ubar_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = uflx/stagthck*vel0;
        vbar_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = vflx/stagthck*vel0;
    elsewhere
        ubar_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = 0.0_dp
        vbar_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = 0.0_dp
    endwhere

    call write_xls("ubar_ir.txt", ubar_ir(:,:,1))
    call write_xls("vbar_ir.txt", vbar_ir(:,:,1))

    call periodic_boundaries(thck_ir(:,:,1), periodic_ew, periodic_ns, 2)
    call periodic_boundaries(ubar_ir(:ewn_ir-1,:nsn_ir-1,1), periodic_ew, periodic_ns, 2)
    call periodic_boundaries(vbar_ir(:ewn_ir-1,:nsn_ir-1,1), periodic_ew, periodic_ns, 2)

    !Copy the extra set of ghost cells over
    !Hard coded 5 as the source for these ghost cells,
    !because we go in two rows for the low-end ghost cells,
    !then go in two more for the source of the first two
    !high-end ghost cells.
    ubar_ir(ewn_ir,:,:) = ubar_ir(5,:,:)
    ubar_ir(:,nsn_ir,:) = ubar_ir(:,5,:)
    vbar_ir(ewn_ir,:,:) = ubar_ir(5,:,:)
    vbar_ir(:,nsn_ir,:) = ubar_ir(:,5,:)
     

    call write_xls("ubar_ir_pbc.txt", ubar_ir(:,:,1))
    call write_xls("vbar_ir_pbc.txt", vbar_ir(:,:,1))    

!whl - to do - Fill the tracer array with ice temperature and other tracers
    trace_ir(:,:,:,1) = 1.0_dp;
    dt_ir = dt * tim0

    where( thck_ir(:,:,1) > 0.0_dp )
        mask_ir  = 1.0_dp
    end where

    end subroutine horizontal_remap_in

!----------------------------------------------------------------------

    subroutine horizontal_remap_out( thck_ir, mask_ir, thck, acab, dt, periodic_ew, periodic_ns )
    use xls
    ! *sfp** take output from inc. remapping and put back in GLAM format

    implicit none

    real (kind = dp), intent(in) :: dt
    real (kind = sp), intent(in), dimension(:,:) :: acab
    real (kind = dp), dimension(:,:,:), intent(inout) :: thck_ir
    real (kind = dp), dimension(:,:), intent(in) :: mask_ir
    real (kind = dp), dimension(:,:), intent(inout) :: thck
    logical, intent(in) :: periodic_ew, periodic_ns

    integer :: ewn, nsn, ngew, ngns

    ewn = size(thck, 1)
    nsn = size(thck, 2)

    ngew = (ewn_ir - ewn)/2
    ngns = (nsn_ir - nsn)/2

    call write_xls("thck_ir.txt", thck_ir(:,:,1))
    call periodic_boundaries(thck_ir(:,:,1), periodic_ew, periodic_ns, 2)
    call write_xls("thck_ir_pbc.txt", thck_ir(:,:,1))


    !Map from IR thickness field back to Glide thickness field
    thck = thck_ir(1+ngew:ngew+ewn, 1+ngns:ngns+nsn,1) / thk0
    
    !Apply accumulation
    thck = thck + acab
    
    !Remove thickness from previously masked out locations
    thck = thck * mask_ir(1+ngew:ngew+ewn, 1+ngns:ngns+nsn)
    

    call write_xls("thck.txt", thck)
    end subroutine horizontal_remap_out

!----------------------------------------------------------------------

end module remap_glamutils

