module fo_upwind_advect

!----------------------------------------------------------------------

    ! init, finalize, and driver subroutines for mass advection based on 1st order upwinding

    use glimmer_paramets, only: sp, dp, len0, thk0, tim0, vel0, tim0, acc0, scyr
    use glide_types
    use glide_velo_higher  

    private
    public :: fo_upwind_advect_init, fo_upwind_advect_driver, fo_upwind_advect_final

    ! allocatable work arrays
    real (kind = dp), allocatable, dimension(:,:) ::    &
                       ubar, vbar,                      & 
                       ubar_grid, vbar_grid,            &
                       flux_net, thck_grid,             & 
                       mask, thck_old 

    contains

!----------------------------------------------------------------------

    subroutine fo_upwind_advect_init( ewn, nsn )
    
    ! initialization for 1st-order upwinding mass advection

    implicit none
   
    integer, intent(in) :: ewn, nsn   ! horizontal grid dimensions
 
    ! allocate work arrays
    allocate( ubar(ewn-1,nsn-1) ); ubar = 0.0_dp
    allocate( vbar(ewn-1,nsn-1) ); vbar = 0.0_dp
    allocate( ubar_grid(ewn+1,nsn+1) ); ubar_grid = 0.0_dp
    allocate( vbar_grid(ewn+1,nsn+1) ); vbar_grid = 0.0_dp
    allocate( thck_grid(ewn+2,nsn+2) ); thck_grid = 0.0_dp
    allocate( flux_net(ewn,nsn) ); flux_net = 0.0_dp 
    allocate( mask(ewn,nsn) ); mask = 0.0_dp
    allocate( thck_old(ewn,nsn) ); thck_old = 0.0_dp

    end subroutine fo_upwind_advect_init

!----------------------------------------------------------------------

    subroutine fo_upwind_advect_final( )

    ! finalization for 1st-order upwinding mass advection

    implicit none   

    ! deallocate work arrays
    deallocate( ubar )
    deallocate( vbar )
    deallocate( ubar_grid )
    deallocate( vbar_grid )
    deallocate( thck_grid )
    deallocate( flux_net )
    deallocate( mask )
    deallocate( thck_old )

    end subroutine fo_upwind_advect_final

!----------------------------------------------------------------------

    subroutine fo_upwind_advect_driver( model )

    ! driver routine for the 1st order, upwind mass transport scheme

        type(glide_global_type), intent(inout) :: model

        call run_ho_diagnostic(model)   ! get velocities and fluxes from HO dynamic subroutines

        call fo_upwind_advect_main( model%geometry%thck,    model%geomderv%stagthck,    &
                                    model%climate%acab,     model%numerics%dt,          &
                                    model%velocity_hom%uflx,model%velocity_hom%vflx,    &
                                    model%general%ewn,      model%general%nsn,          &
                                    model%numerics%dew,     model%numerics%dns )

    end subroutine fo_upwind_advect_driver


!----------------------------------------------------------------------

    subroutine fo_upwind_advect_main( thck, stagthck, acab, dt, uflx, vflx, ewn, nsn, dew, dns )

    ! 1st-order upwinding mass advection that uses a finite-volume like scheme for 
    ! mass conservation. Velocities from the staggered grid (B-grid) are averaged onto the 
    ! faces of the non-staggered grid (i.e. faces of the grid where scalers like thickness live).  
    ! Thus, the averaged velocities exist on a C-grid, allowing mass transport to be treated  
    ! in a finite-volume manner; depth averaged velocities give the fluxes out of each cell 
    ! centered on a thickness point and the thickness advected is chosen according to upwinding.
    ! 
    ! Note that this works at the calving front because a non-zero staggered thickness there 
    ! defines the velocities there. These velocites can be used to define the velocity at
    ! the face of the last non-zero thickness cell (on the normal grid) which corresponds to
    ! the location of the calving front. 
    !
    ! Note also that this code has NOT been tested extensively for anything other than a 
    ! simple shelf configuration (i.e. the "confined-shelf" in the tests/ directory.

    implicit none

    real (kind = dp), intent(in) :: dt
    real (kind = dp), dimension(:,:), intent(inout) :: thck
    real (kind = dp), dimension(:,:), intent(in) :: stagthck
    real (kind = sp), dimension(:,:), intent(in) :: acab
    real (kind = dp), dimension(:,:), intent(in) :: uflx, vflx
    real (kind = dp), intent(in) :: dew, dns  
    integer, intent(in)  :: ewn, nsn

    real (kind = dp) :: He, Hw, Hn, Hs, ue, uw, vn, vs  ! upwinding variables and interface velocities

    integer :: ew, ns 

    where( thck > 0.0_dp )      ! mask for eventually removing flux outside of the original domain
        mask = 1.0_dp           ! (i.e. stuff that moves past the calving front goes away)
    else where
        mask = 0.0_dp
    end where

    where( stagthck > 0.0_dp )  ! calculate the depth-ave velocities
        ubar = uflx / stagthck
        vbar = vflx / stagthck  
    end where

    thck_old = thck             ! save the old thickness for debugging purposes

    ! fill in the interior values on the extended velocity grid (extended B-grid)
    ubar_grid(2:ewn,2:nsn) = ubar    
    vbar_grid(2:ewn,2:nsn) = vbar    

    ! fill in the interior values on the extended thickness grid
    thck_grid(2:ewn+1,2:nsn+1) = thck(:,:)

    ! calculate the interface velocities from the extended B-grid, then use upwinding
    ! criterion to advect thickness in or out of cells (NOTE that parts of this could
    ! probably be vectorized at some point)
    do ns = 1, nsn
        do ew = 1, ewn

            ! interface depth-ave velocities
            ue = ( ubar_grid(ew+1,ns+1) + ubar_grid(ew+1,ns) ) / 2.0d0
            uw = ( ubar_grid(ew,ns+1) + ubar_grid(ew,ns) ) / 2.0d0
            vn = ( vbar_grid(ew,ns+1) + vbar_grid(ew+1,ns+1) ) / 2.0d0            
            vs = ( vbar_grid(ew,ns) + vbar_grid(ew+1,ns) ) / 2.0d0            

            ! choose thickness to advect based on upwinding
            if( ue > 0.0d0 )then
                He = - thck_grid(ew+1,ns+1) ! negative signs necessary so that flux to the east
            else                            ! results in mass loss in this volume (and vice versa)
                He = - thck_grid(ew+2,ns+1)
            end if
            if( uw > 0.0d0 )then
                Hw = thck_grid(ew,ns+1)
            else
                Hw = thck_grid(ew+1,ns+1)
            end if
            if( vn > 0.0d0 )then
                Hn = - thck_grid(ew+1,ns+1) ! negative signs here as above for ue, and He
            else
                Hn = - thck_grid(ew+1,ns+2)
            end if
            if( vs > 0.0d0 )then
                Hs = thck_grid(ew+1,ns)
            else
                Hs = thck_grid(ew+1,ns+1)
            end if

            ! net flux into/out of each cell
            flux_net(ew,ns) = ( ue*He*dns + uw*Hw*dns + vn*Hn*dew + vs*Hs*dew ) 

        end do
    end do

    thck = thck_old + ( 1 / (dns * dew) * flux_net ) * dt + (acab * dt)

    ! debugging
    print *, ' '
    print *, 'net volume change = ', sum( (thck-thck_old)*mask )*thk0 *dew*dns*len0**2 
    print *, 'net calving flux = ', sum( thck * (1.0d0-mask) )*thk0*dew*dns*len0**2
    print *, '(for the confined shelf experiment, the above two should sum to ~0)'
    print *, 'mean accum/ablat rate = ', sum( acab * (1.0d0-mask) ) / sum(mask) / (dt*tim0) * scyr
    print *, 'mean dH/dt = ', sum( (thck-thck_old)*mask )*thk0 / sum(mask) / (dt*tim0) * scyr
    print *, 'sum of flux change (should be ~0) = ', sum( flux_net*vel0*thk0*len0 ) 
    print *, ' '
!    pause

    thck = thck * mask               ! remove any mass advected outside of initial domain

    where( thck < 0.0_dp )           ! gaurd against thickness going negative
        thck = 0.0_dp
    end where

    end subroutine fo_upwind_advect_main

!----------------------------------------------------------------------

end module fo_upwind_advect
