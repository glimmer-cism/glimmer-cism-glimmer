#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

#include "glide_nan.inc"

#define shapedbg(x) write(*,*) "x", shape(x)

module glide_velo_higher
    use ice3d_lib
    use glimmer_global, only : dp
    use glide_types
    use glide_vertint
    use glide_grids, only: stagvarb
    use glimmer_physcon, only: gn
    use glide_mask
    implicit none
    
    !TODO: Parameterize the following globals
    real(dp), parameter :: VEL2ERR  = 1e-4
    real(dp), parameter :: TOLER    = 1e-6
    integer,  parameter :: CONVT    = 4
    real(dp), parameter :: SHTUNE   = 1.D-16
    integer,  parameter :: UPSTREAM = 0
    integer,  parameter :: MANIFOLD = 1
contains
        
    subroutine init_velo_hom_pattyn(model)
        type(glide_global_type) :: model

        !Init the beta field based on the selected option
        !If we selected to use 1/btrc, this needs to be computed later
        !If we selected to read the beta field, it's been done already and no
        !action is needed
        select case (model%options%which_ho_beta_in)
            case(HO_BETA_ALL_NAN)
                model%velocity_hom%beta = NAN
            case(HO_BETA_USE_SOFT)
                where (model%velocity%bed_softness /= 0)
                    model%velocity_hom%beta = 1 / model%velocity%bed_softness
                elsewhere
                    model%velocity_hom%beta = NAN
                end where
        end select

    end subroutine

    subroutine calc_slip_ratio(flwa, thick, slip_ratio_coef, beta)
        real(dp), dimension(:,:), intent(in) :: flwa !*FD Glen's A at the base
        real(dp), dimension(:,:), intent(in) :: thick !*FD Ice thickness
        real(dp),                 intent(in) :: slip_ratio_coef
        real(dp), dimension(:,:), intent(out):: beta

        beta = 1.0D0/(slip_ratio_coef*flwa*thick)
    end subroutine  

    !This is a temporary wrapper function to get all the HO setup and calling
    !code out of glide_thck so it keeps the clutter down.  I'm just passing the
    !model right now because the analysis on what needs to be passed has yet to
    !be done.
    subroutine run_ho_diagnostic(model)
        use glide_thckmask
        use glide_mask
        
        type(glide_global_type),intent(inout) :: model
        !For HO masking
        logical :: empty
        integer :: totpts
        real(sp), dimension(model%general%ewn-1, model%general%nsn-1) :: stagmassb

        !TEMPORARY arrays, these should at some point be placed in Model
        !probably
        integer, dimension(model%general%ewn-1, model%general%nsn-1)  :: geom_mask_stag
        real(dp), dimension(model%general%ewn-1, model%general%nsn-1) :: latbc_norms_stag

        !Beta field computations that change in time
        if (model%options%which_ho_beta_in == HO_BETA_USE_BTRC) then
           where (model%velocity%btrc /= 0)
               model%velocity_hom%beta = 1/model%velocity%btrc
            elsewhere
                model%velocity_hom%beta = NAN
            end where
        else if (model%options%which_ho_beta_in == HO_BETA_SLIP_RATIO) then
            call calc_slip_ratio(model%temper%flwa(model%general%upn,:,:), model%geomderv%stagthck, &
                                 model%paramets%slip_ratio, model%velocity_hom%beta)
        end if

            
        if (model%options%which_ho_diagnostic == HO_DIAG_PATTYN_STAGGERED) then
            write(*,*)"Running Pattyn staggered"
         
            !Compute the "point mask" (mask that assigns a unique value to each
            !grid point on which ice dynamics are enabled) for the
            !staggered grid.  
            stagmassb = 0

            
            call glide_maskthck(model%geomderv%stagthck, stagmassb, .true., model%geometry%dom, &
                                model%velocity_hom%velmask, totpts, empty)
                 
            !Compute the "geometry mask" (type of square) for the staggered grid

            call glide_set_mask(model%numerics, model%geomderv%stagthck, model%geomderv%stagtopg, &
                                model%general%ewn-1, model%general%nsn-1, model%climate%eus, &
                                geom_mask_stag) 

            !Compute the normal vectors to the marine margin for the staggered grid
            call glide_marine_margin_normal(model%geomderv%stagthck, geom_mask_stag, latbc_norms_stag)
                
            call velo_hom_pattyn(model%general%ewn, model%general%nsn, model%general%upn, &
                                 model%numerics%dew, model%numerics%dns, model%numerics%sigma, &
                                 model%geomderv%stagthck, model%geomderv%stagusrf, model%geomderv%staglsrf, &
                                 model%geomderv%dthckdew, model%geomderv%dthckdns, &
                                 model%geomderv%dusrfdew, model%geomderv%dusrfdns, &
                                 model%geomderv%dlsrfdew, model%geomderv%dlsrfdns, & 
                                 model%geomderv%d2usrfdew2, model%geomderv%d2usrfdns2, &
                                 model%geomderv%d2thckdew2, model%geomderv%d2thckdns2, &
                                 model%velocity_hom%velmask, totpts, &
                                 geom_mask_stag, &
                                 model%temper%flwa, real(gn, dp), model%velocity_hom%beta, &
                                 model%options%which_ho_bstress,&
                                 model%options%which_ho_efvs, &
                                 model%options%periodic_ew .eq. 1, &
                                 model%options%periodic_ns .eq. 1,&
                                 model%velocity_hom%kinematic_bc_u, model%velocity_hom%kinematic_bc_v, &
                                 latbc_norms_stag,&
                                 model%velocity_hom%uvel, model%velocity_hom%vvel, &
                                 model%velocity_hom%is_velocity_valid, &
                                 model%velocity_hom%uflx, model%velocity_hom%vflx, &
                                 model%velocity_hom%efvs, model%velocity_hom%tau, &
                                 model%velocity_hom%gdsx, model%velocity_hom%gdsy)
        else if (model%options%which_ho_diagnostic == HO_DIAG_PATTYN_UNSTAGGERED) then
            write(*,*)"Running Pattyn unstaggered"
            call velo_hom_pattyn_nonstag(model%general%ewn, model%general%nsn, model%general%upn, &
                                         model%numerics%dew, model%numerics%dns, model%numerics%sigma, &
                                         model%geometry%thck, model%geometry%usrf, model%geometry%lsrf, &
                                         model%geomderv%dusrfdew_unstag, model%geomderv%dusrfdns_unstag, &
                                         model%geomderv%dthckdew_unstag, model%geomderv%dthckdns_unstag, &
                                         model%geomderv%dlsrfdew_unstag, model%geomderv%dlsrfdns_unstag, &
                                         model%geomderv%d2usrfdew2_unstag, model%geomderv%d2usrfdns2_unstag, &
                                         model%geomderv%d2thckdew2_unstag, model%geomderv%d2thckdns2_unstag, &
                                         model%geometry%mask, model%geometry%totpts, &
                                         model%geometry%thkmask, &
                                         model%temper%flwa, real(gn, dp), &
                                         model%velocity_hom%beta, model%options%which_ho_bstress, &
                                         model%options%which_ho_efvs, &
                                         model%options%periodic_ew .eq. 1, &
                                         model%options%periodic_ns .eq. 1, &
                                         model%velocity_hom%kinematic_bc_u, model%velocity_hom%kinematic_bc_v, &
                                         model%geometry%marine_bc_normal, &
                                         model%velocity_hom%uvel, model%velocity_hom%vvel, &
                                         model%velocity_hom%is_velocity_valid)
        end if

        model%velocity_hom%is_velocity_valid = .true.
        
    end subroutine

    subroutine velo_hom_pattyn(ewn, nsn, upn, dew, dns, sigma, &
                               thck, usrf, lsrf, dthckdew, dthckdns, dusrfdew, dusrfdns, &
                               dlsrfdew, dlsrfdns, &
                               d2zdx2, d2zdy2, d2hdx2, d2hdy2, &
                               point_mask, totpts, geometry_mask, flwa, flwn, btrc, &
                               which_sliding_law, which_efvs, &
                               periodic_ew, periodic_ns, kinematic_bc_u, kinematic_bc_v, &
                               marine_bc_normal, &
                               uvel, vvel, valid_initial_guess, uflx, vflx, efvs, tau, gdsx, gdsy)
                           
        integer, intent(in)  :: ewn !*FD Number of cells X
        integer, intent(in)  :: nsn !*FD Number of cells Y
        integer, intent(in)  :: upn !*FD Number of cells Z
        real(dp), intent(in) :: dew
        real(dp), intent(in) :: dns
        real(dp), dimension(:), intent(in) :: sigma !*FD Sigma coord for rescaled Z dimension
        real(dp), dimension(:,:), intent(in) :: thck !*FD Thickness, on staggered grid
	real(dp), dimension(:,:), intent(in) :: usrf !*FD Upper surface profile
        real(dp), dimension(:,:), intent(in) :: lsrf !*FD Lower surface profile
	real(dp), dimension(:,:), intent(in) :: dthckdew !*FD X thickness gradient
	real(dp), dimension(:,:), intent(in) :: dthckdns !*FD Y thickness gradient
	real(dp), dimension(:,:), intent(in) :: dusrfdew !*FD X surface gradient
	real(dp), dimension(:,:), intent(in) :: dusrfdns !*FD Y surface gradient
	real(dp), dimension(:,:), intent(in) :: dlsrfdew !*FD X bed gradient
	real(dp), dimension(:,:), intent(in) :: dlsrfdns !*FD Y bed gradient
        real(dp), dimension(:,:), intent(in) :: d2zdx2, d2zdy2, d2hdx2, d2hdy2
        integer,  dimension(:,:), intent(in) :: point_mask     !*FD Numbers points in the staggered grid that are included in computation
        integer, intent(in) :: totpts
        integer, dimension(:,:), intent(in) :: geometry_mask
        real(dp), dimension(:,:,:), intent(in) :: flwa !*FD Glen's A (rate factor) - Used for thermomechanical coupling
        real(dp), dimension(:,:), intent(in)   :: btrc !*FD Basal Traction, either betasquared or tau0
        real(dp), dimension(:,:,:), intent(in)   :: kinematic_bc_u, kinematic_bc_v
        real(dp), dimension(:,:), intent(in)   :: marine_bc_normal
        real(dp), intent(in) :: flwn !*FD Exponent in Glenn power law
        integer, intent(in) :: which_sliding_law
        integer, intent(in) :: which_efvs
        logical, intent(in) :: periodic_ew !*Whether to use periodic boundary conditions
        logical, intent(in) :: periodic_ns
        real(dp), dimension(:,:,:), intent(inout) :: uvel 
        real(dp), dimension(:,:,:), intent(inout) :: vvel
        logical, intent(in) :: valid_initial_guess !*Whether or not the given uvel or vvel are appropriate initial guesses.  If not we'll have to roll our own.
        real(dp), dimension(:,:), intent(out) :: uflx
        real(dp), dimension(:,:), intent(out) :: vflx
        real(dp), dimension(:,:,:), intent(out) :: efvs !*FD Effective viscosity
        type(glide_tensor), intent(inout)         :: tau
        real(dp), dimension(:,:,:), intent(out) :: gdsx !*FD X driving stress
        real(dp), dimension(:,:,:), intent(out) :: gdsy !*FD Y driving stress

        !Second derivative of surface
                !Arrays for rescaled coordinate parameters
        real(dp), dimension(nsn-1, ewn-1, upn) :: ax, ay, bx, by, cxy
        
        real(dp), dimension(nsn-1, ewn-1, upn) :: kinematic_bc_u_t, kinematic_bc_v_t

        real(dp), dimension(nsn, ewn)::u,v

        !Arrays to hold transposed data
        real(dp), dimension(nsn-1, ewn-1) :: dlsrfdew_t, dlsrfdns_t, dusrfdew_t, & 
        dusrfdns_t, dthckdew_t, dthckdns_t, btrc_t, staglsrf_t, stagusrf_t, stagthck_t, &
        d2zdx2_t, d2zdy2_t, d2hdx2_t, d2hdy2_t
       
        real(dp), dimension(nsn, ewn, upn) :: flwa_t
        real(dp), dimension(nsn-1, ewn-1, upn) :: uvel_t, vvel_t, mu_t, flwa_t_stag, &
                                                  tau_xz_t, tau_yz_t, tau_xx_t, tau_yy_t, tau_xy_t

        integer,  dimension(nsn-1, ewn-1) :: point_mask_t
        integer,  dimension(nsn-1, ewn-1) :: geometry_mask_t
        real(dp), dimension(nsn-1, ewn-1) :: marine_bc_normal_t       
        real(dp), dimension(nsn-1, ewn-1) :: direction_x, direction_y
        
        direction_x = 0
        direction_y = 0

        !Move everything into transposed coordinates
        !This means that we need to transpose every data field with x and y that
        !we have thus far.
        !It also means that we need to swap x and y derivatives
        stagthck_t = transpose(thck)
        staglsrf_t = transpose(lsrf)
        stagusrf_t = transpose(usrf)
        btrc_t = transpose(btrc)

        !TODO: Recompute these derivatives?
        dthckdew_t = transpose(dthckdew)
        dthckdns_t = transpose(dthckdns)
        dlsrfdew_t = transpose(dlsrfdew)
        dlsrfdns_t = transpose(dlsrfdns)
        dusrfdew_t = transpose(dusrfdew)
        dusrfdns_t = transpose(dusrfdns)

        d2zdx2_t = transpose(d2zdx2)
        d2zdy2_t = transpose(d2zdy2)
        d2hdx2_t = transpose(d2hdx2)
        d2hdy2_t = transpose(d2hdy2)

        point_mask_t = transpose(point_mask)
        geometry_mask_t = transpose(geometry_mask)
        marine_bc_normal_t = transpose(marine_bc_normal)

        call glimToIce3d_3d(uvel,uvel_t,ewn-1,nsn-1,upn)
        call glimToIce3d_3d(vvel,vvel_t,ewn-1,nsn-1,upn)

        call glimToIce3d_3d(kinematic_bc_u, kinematic_bc_u_t, ewn-1, nsn-1, upn)
        call glimToIce3d_3d(kinematic_bc_v, kinematic_bc_v_t, ewn-1, nsn-1, upn)

        !The flow law field needs to be transposed AND staggered, which is a lot
        !of fun!!!
        call glimToIce3d_3d(flwa,flwa_t,ewn,nsn,upn)
        call staggered_field_3d(flwa_t, flwa_t_stag)

        !Compute a new geometry mask based on staggered geometry
       
        !Compute rescaled coordinate parameters (needed because Pattyn uses an
        !irregular Z grid and scales so that 0 is the surface, 1 is the bed)
        call init_rescaled_coordinates(dthckdew_t,dlsrfdew_t,dthckdns_t,dlsrfdns_t,stagusrf_t,stagthck_t,staglsrf_t,&
                                               dusrfdew_t,dusrfdns_t,d2zdx2_t,d2zdy2_t,d2hdx2_t,d2hdy2_t,&
                                               sigma,ax,ay,bx,by,cxy,dew,dns,direction_x,direction_y)
        !"Spin up" estimate with Pattyn's SIA model runs if we don't already
        !have a good initial guess
        if (.not. valid_initial_guess) then
            call veloc1(dusrfdew_t, dusrfdns_t, stagthck_t, flwa_t_stag, sigma, uvel_t, vvel_t, u, v, nsn-1, ewn-1, upn, &
                        FLWN, periodic_ew, periodic_ns)
            !If we are performing the plastic bed iteration, the SIA is not
            !enough and we need to spin up a better estimate by shoehorning the
            !tau0 values into a linear bed estimate
            if (which_sliding_law == 1) then
                call veloc2(mu_t, uvel_t, vvel_t, flwa_t_stag, dusrfdew_t, dusrfdns_t, stagthck_t, ax, ay, &
                        sigma, bx, by, cxy, btrc_t/100, dlsrfdew_t, dlsrfdns_t, FLWN, ZIP, VEL2ERR, MANIFOLD,&
                        TOLER, periodic_ew,periodic_ns, 0, which_efvs, .true., dew, dns,point_mask_t,totpts,geometry_mask_t, &
                        kinematic_bc_u_t, kinematic_bc_v_t, marine_bc_normal_t)
            end if
        end if

        !Higher order velocity estimation
        !I am assuming that efvs (effective viscosity) is the same as mu
        !A NOTE ON COORDINATE TRANSPOSITION:
        !Because of the transposition, ewn=maxy and nsn=maxx.  However, veloc2
        !passes maxy in *first*, so these really get passed in the same order
        !that they normally would.
        call veloc2(mu_t, uvel_t, vvel_t, flwa_t, dusrfdew_t, dusrfdns_t, stagthck_t, ax, ay, &
                    sigma, bx, by, cxy, btrc_t, dlsrfdew_t, dlsrfdns_t, FLWN, ZIP, VEL2ERR, MANIFOLD,&
                    TOLER, periodic_ew,periodic_ns, which_sliding_law, which_efvs, .true., dew,&
                    dns,point_mask_t,totpts, geometry_mask_t, kinematic_bc_u_t, kinematic_bc_v_t, marine_bc_normal_t)
       
        !Final computation of stress field for output
        call stressf(mu_t, uvel_t, vvel_t, flwa_t, stagthck_t, ax, ay, dew, dns, sigma, & 
                     tau_xz_t, tau_yz_t, tau_xx_t, tau_yy_t, tau_xy_t, flwn, zip, periodic_ew, periodic_ns) 
    
        !Transpose from the ice3d coordinate system (y,x,z) to the glimmer
        !coordinate system (z,x,y).  We need to do this for all the 3D outputs
        !that Glimmer expects
        call ice3dToGlim_3d(uvel_t, uvel, ewn-1, nsn-1, upn)
        call ice3dToGlim_3d(vvel_t, vvel, ewn-1, nsn-1, upn)
        call ice3dToGlim_3d(mu_t, efvs, ewn-1, nsn-1, upn)
        call ice3dToGlim_3d(tau_xz_t, tau%xz, ewn-1, nsn-1, upn)
        call ice3dToGlim_3d(tau_yz_t, tau%yz, ewn-1, nsn-1, upn)
        call ice3dToGlim_3d(tau_xx_t, tau%xx, ewn-1, nsn-1, upn)
        call ice3dToGlim_3d(tau_yy_t, tau%yy, ewn-1, nsn-1, upn)
        call ice3dToGlim_3d(tau_xy_t, tau%xy, ewn-1, nsn-1, upn)
    end subroutine velo_hom_pattyn

    subroutine hom_diffusion_pattyn(uvel_hom, vvel_hom, stagthck, dusrfdew, dusrfdns, sigma, diffu_x, diffu_y)
        !*FD Estimate of higher-order diffusivity vector given Pattyn's approach
        !*FD Implements (54) in Pattyn 2003, including vertical integration
        real(dp), dimension(:,:,:), intent(in) :: uvel_hom
        real(dp), dimension(:,:,:), intent(in) :: vvel_hom
        real(dp), dimension(:,:),   intent(in) :: stagthck
        real(dp), dimension(:,:),   intent(in) :: dusrfdew
        real(dp), dimension(:,:),   intent(in) :: dusrfdns
        real(dp), dimension(:),     intent(in) :: sigma
        real(dp), dimension(:,:),   intent(out):: diffu_x
        real(dp), dimension(:,:),   intent(out):: diffu_y

        !Local Variables
        real(dp), dimension(size(uvel_hom, 2), size(uvel_hom, 3)) :: u_int !*FD Vertically integrated u velocity
        real(dp), dimension(size(uvel_hom, 2), size(uvel_hom, 3)) :: v_int !*FD Vertically integrated v velocity

        call vertint_output2d(uvel_hom, u_int, sigma)
        call vertint_output2d(vvel_hom, v_int, sigma)

        diffu_x = u_int*stagthck*dusrfdew
        diffu_y = u_int*stagthck*dusrfdns

    end subroutine hom_diffusion_pattyn


    !*FD Version of higher order velo call that computes velocities on ice grid
    !*FD initially and then staggers them.  This might introduce less error than
    !computing velocities on a staggered grid; it is thought that the averaging
    !there causes problems.
    subroutine velo_hom_pattyn_nonstag(ewn, nsn, upn, dew, dns, sigma, thck, usrf, lsrf, &
                              dusrfdew, dusrfdns, dthckdew, dthckdns, dlsrfdew, dlsrfdns, &
                              d2zdx2, d2zdy2, d2hdx2, d2hdy2, &
                              point_mask, totpts, geometry_mask, flwa, flwn, btrc, which_sliding_law, &
                              which_efvs, periodic_ew, periodic_ns, kinematic_bc_u, kinematic_bc_v, marine_bc_normal,&
                              uvel, vvel, valid_initial_guess)
                            
        integer, intent(in)  :: ewn !*FD Number of cells X
        integer, intent(in)  :: nsn !*FD Number of cells Y
        integer, intent(in)  :: upn !*FD Number of cells Z
        real(dp), intent(in) :: dew !*FD Grid spacing X
	real(dp), intent(in) :: dns !*FD Grid spacing Y
	real(dp), dimension(:), intent(in) :: sigma !*FD Sigma coord for rescaled Z dimension
	real(dp), dimension(:,:), intent(in) :: thck !*FD Thickness, on non-staggered grid
	real(dp), dimension(:,:), intent(in) :: usrf !*FD Upper surface profile
        real(dp), dimension(:,:), intent(in) :: lsrf !*FD Lower surface profile
	real(dp), dimension(:,:), intent(in) :: dthckdew !*FD X thickness gradient
	real(dp), dimension(:,:), intent(in) :: dthckdns !*FD Y thickness gradient
	real(dp), dimension(:,:), intent(in) :: dusrfdew !*FD X surface gradient
	real(dp), dimension(:,:), intent(in) :: dusrfdns !*FD Y surface gradient
	real(dp), dimension(:,:), intent(in) :: dlsrfdew !*FD X bed gradient
	real(dp), dimension(:,:), intent(in) :: dlsrfdns !*FD Y bed gradient
        real(dp), dimension(:,:), intent(in) :: d2zdx2, d2zdy2, d2hdx2, d2hdy2
        integer,  dimension(:,:), intent(in) :: point_mask     !*FD Numbers points in the staggered grid that are included in computation
        integer, intent(in) :: totpts
        integer, dimension(:,:), intent(in) :: geometry_mask
        real(dp), dimension(:,:,:), intent(in) :: flwa !*FD Glen's A (rate factor) - Used for thermomechanical coupling
        real(dp), dimension(:,:), intent(in)   :: btrc !*FD Basal Traction, either betasquared or tau0
        real(dp), dimension(:,:,:), intent(in)   :: kinematic_bc_u, kinematic_bc_v
        real(dp), dimension(:,:), intent(in)   :: marine_bc_normal
        real(dp), intent(in) :: flwn !*FD Exponent in Glenn power law
        integer, intent(in) :: which_sliding_law
        integer, intent(in) :: which_efvs
        logical, intent(in) :: periodic_ew !*Whether to use periodic boundary conditions
        logical :: periodic_ns
        real(dp), dimension(:,:,:), intent(out) :: uvel 
        real(dp), dimension(:,:,:), intent(out) :: vvel
        logical, intent(in) :: valid_initial_guess !*Whether or not the given uvel or vvel are appropriate initial guesses.  If not we'll have to roll our own.
        
        !Arrays for rescaled coordinate parameters
        real(dp), dimension(nsn, ewn, upn) :: ax, ay, bx, by, cxy
        
        real(dp), dimension(nsn, ewn)::u,v

        !whether to upwind or downwind derivatives.
        real(dp), dimension(ewn, nsn) :: direction_x, direction_y
        
        !Arrays to hold transposed data
        real(dp), dimension(nsn, ewn) :: dlsrfdew_t, dlsrfdns_t, dusrfdew_t, & 
        dusrfdns_t, dthckdew_t, dthckdns_t, btrc_t_unstag,&
        d2zdx2_t, d2zdy2_t, d2hdx2_t, d2hdy2_t,usrf_t, lsrf_t, thck_t

        real(dp), dimension(nsn-1, ewn-1) :: btrc_t 
        real(dp), dimension(nsn, ewn, upn) :: flwa_t
        real(dp), dimension(nsn-1, ewn-1, upn) :: uvel_t, vvel_t
        real(dp), dimension(nsn, ewn, upn) :: mu_t
        real(dp), dimension(nsn, ewn, upn) :: uvel_t_unstag, vvel_t_unstag
        integer, dimension(nsn, ewn) :: point_mask_t
        integer, dimension(nsn, ewn) :: geometry_mask_t

        real(dp), dimension(nsn-1, ewn-1, upn) :: kinematic_bc_u_t
        real(dp), dimension(nsn-1, ewn-1, upn) :: kinematic_bc_v_t

        real(dp), dimension(nsn, ewn, upn) :: kinematic_bc_u_t_unstag
        real(dp), dimension(nsn, ewn, upn) :: kinematic_bc_v_t_unstag

        real(dp), dimension(nsn,ewn) :: marine_bc_normal_t

        !Determine whether to upwind or downwind derivatives at points on the
        !interior of the model domain (this is mainly important for the marine
        !margin
        !Since this is a call to Pattyn's model, derivatives need to be
        !transposed.
        direction_y = 0
        direction_x = 0
        call upwind_from_mask(geometry_mask, direction_y, direction_x)

        !Move everything into transposed coordinates
        !This means that we need to transpose every data field with x and y that
        !we have thus far.
        !It also means that we need to swap x and y derivatives
        thck_t = transpose(thck)
        lsrf_t = transpose(lsrf)
        usrf_t = transpose(usrf)
        btrc_t = transpose(btrc)
        
        dthckdew_t = transpose(dthckdew)
        dthckdns_t = transpose(dthckdns)
        dlsrfdew_t = transpose(dlsrfdew)
        dlsrfdns_t = transpose(dlsrfdns)
        dusrfdew_t = transpose(dusrfdew)
        dusrfdns_t = transpose(dusrfdns)

        d2zdx2_t = transpose(d2zdx2)
        d2zdy2_t = transpose(d2zdy2)
        d2hdx2_t = transpose(d2hdx2)
        d2hdy2_t = transpose(d2hdy2)

        point_mask_t = transpose(point_mask)
        geometry_mask_t = transpose(geometry_mask)

        marine_bc_normal_t = transpose(marine_bc_normal)
       
        call glimToIce3d_3d(uvel,uvel_t,ewn-1,nsn-1,upn)
        call glimToIce3d_3d(vvel,vvel_t,ewn-1,nsn-1,upn)
        
        call glimToIce3d_3d(kinematic_bc_u,kinematic_bc_u_t,ewn-1,nsn-1,upn)
        call glimToIce3d_3d(kinematic_bc_v,kinematic_bc_v_t,ewn-1,nsn-1,upn)

        call glimToIce3d_3d(flwa,flwa_t,ewn,nsn,upn)
      
        call unstagger_field_3d(vvel_t, uvel_t_unstag, periodic_ew, periodic_ns)
        call unstagger_field_3d(vvel_t, uvel_t_unstag, periodic_ew, periodic_ns)
       
        call unstagger_field_3d(kinematic_bc_u_t,kinematic_bc_u_t_unstag, periodic_ew, periodic_ns)
        call unstagger_field_3d(kinematic_bc_v_t,kinematic_bc_v_t_unstag, periodic_ew, periodic_ns)
        
        !In unstaggering the boundary condition fields, we need to remove points
        !that aren't on the boundary
        !do i = 1,nsn
        !    do j = 1, ewn
        !        if (i > 1 .and. i < nsn .and. j > 1 .and. j < ewn) then 
        !            if (thck(i+1, j) /= 0 .and. thck(i-1, j) /= 0 .and. thck(i, j+1) /= 0 .and. thck(i, j-1) /= 0) then
        !                kinematic_bc_u_t_unstag(i,j,:) = NaN
        !                kinematic_bc_v_t_unstag(i,j,:) = NaN
        !            end if
        !        end if
        !    end do
        !end do

        call unstagger_field_2d(btrc_t, btrc_t_unstag, periodic_ew, periodic_ns)

        !Compute rescaled coordinate parameters (needed because Pattyn uses an
        !irregular Z grid and scales so that 0 is the surface, 1 is the bed)
        call init_rescaled_coordinates(dthckdew_t,dlsrfdew_t,dthckdns_t,dlsrfdns_t,usrf_t,thck_t,lsrf_t,&
                                               dusrfdew_t,dusrfdns_t,d2zdx2_t,d2zdy2_t,d2hdx2_t,d2hdy2_t,&
                                               sigma,ax,ay,bx,by,cxy,dew,dns,direction_x,direction_y)
       
        !"Spin up" estimate with Pattyn's SIA model runs if we don't already
        !have a good initial guess
        if (.not. valid_initial_guess) then
            call veloc1(dusrfdew_t, dusrfdns_t, thck_t, flwa_t, sigma, uvel_t_unstag, vvel_t_unstag, u, v, nsn, ewn, upn, &
                        FLWN, periodic_ew, periodic_ns)
            !If we are performing the plastic bed iteration, the SIA is not
            !enough and we need to spin up a better estimate by shoehorning the
            !tau0 values into a linear bed estimate
            if (which_sliding_law == 1) then
                call veloc2(mu_t, uvel_t, vvel_t, flwa_t, dusrfdew_t, dusrfdns_t, thck_t, ax, ay, &
                        sigma, bx, by, cxy, btrc_t_unstag, dlsrfdew_t, dlsrfdns_t, FLWN, ZIP, VEL2ERR, MANIFOLD,&
                        TOLER, periodic_ew,periodic_ns, 0, which_efvs, .true., dew, dns,point_mask_t,totpts,geometry_mask_t,&
                        kinematic_bc_u_t_unstag, kinematic_bc_v_t_unstag, marine_bc_normal_t)
            end if
        end if
        
        !Higher order velocity estimation
        !I am assuming that efvs (effective viscosity) is the same as mu
        !A NOTE ON COORDINATE TRANSPOSITION:
        !Because of the transposition, ewn=maxy and nsn=maxx.  However, veloc2
        !passes maxy in *first*, so these really get passed in the same order
        !that they normally would.
        call veloc2(mu_t, uvel_t_unstag, vvel_t_unstag, flwa_t, dusrfdew_t, dusrfdns_t, thck_t, ax, ay, &
                    sigma, bx, by, cxy, btrc_t_unstag, dlsrfdew_t, dlsrfdns_t, FLWN, ZIP, VEL2ERR, MANIFOLD,&
                    TOLER, periodic_ew,periodic_ns, which_sliding_law, which_efvs, .true., dew,&
                    dns,point_mask_t,totpts,geometry_mask_t,kinematic_bc_u_t_unstag, kinematic_bc_v_t_unstag, marine_bc_normal_t)
       
        !Final computation of stress field for output
        !call stressf(mu_t, uvel_t, vvel_t, flwa_t, stagthck_t, ax, ay, dew, dns, sigma, & 
        !             tau_xz_t, tau_yz_t, tau_xx_t, tau_yy_t, tau_xy_t, flwn, zip, periodic_ew, periodic_ns) 
        
        
        call staggered_field_3d(uvel_t_unstag, uvel_t)
        call staggered_field_3d(vvel_t_unstag, vvel_t)
        
        !Transpose from the ice3d coordinate system (y,x,z) to the glimmer
        !coordinate system (z,x,y).  We need to do this for all the 3D outputs
        !that Glimmer expects
        call ice3dToGlim_3d(uvel_t, uvel, ewn-1, nsn-1, upn)
        call ice3dToGlim_3d(vvel_t, vvel, ewn-1, nsn-1, upn)
        !call ice3dToGlim_3d(mu_t, efvs, ewn-1, nsn-1, upn)
        
        !call ice3dToGlim_3d(tau_xz_t, tau%xz, ewn-1, nsn-1, upn)
        !call ice3dToGlim_3d(tau_yz_t, tau%yz, ewn-1, nsn-1, upn)
        !call ice3dToGlim_3d(tau_xx_t, tau%xx, ewn-1, nsn-1, upn)
        !call ice3dToGlim_3d(tau_yy_t, tau%yy, ewn-1, nsn-1, upn)
        !call ice3dToGlim_3d(tau_xy_t, tau%xy, ewn-1, nsn-1, upn)
    end subroutine velo_hom_pattyn_nonstag


    !Transposes from Glimmer 3D fields, stored (z,x,y), to Ice3D fields, stored
    !(y,x,z)
    subroutine glimToIce3d_3d(field, dest, nx, ny, nz)
        real(dp), dimension(:, :, :) :: field
        real(dp), dimension(:, :, :) :: dest
        integer :: nx, ny, nz, i, j, k
 
        do i=1,nx
            do j=1,ny
                do k=1,nz
                    dest(j,i,k) = field(k,i,j)
                end do
            end do
        end do
    end subroutine glimToIce3d_3d

    subroutine ice3dToGlim_3d(field, dest, nx, ny, nz)
        real(dp), dimension(:, :, :) :: field
        real(dp), dimension(:, :, :) :: dest
        integer :: nx, ny, nz, i, j, k
        do i=1,nx
            do j=1,ny
                do k=1,nz
                    dest(k,i,j) = field(j,i,k)
                end do
            end do
        end do
    end subroutine ice3dToGlim_3d

end module glide_velo_higher
