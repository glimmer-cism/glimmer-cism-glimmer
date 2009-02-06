#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

#include "glide_nan.inc"

module glide_velo_higher
    use ice3d_lib
    use glimmer_global, only : dp
    use glide_types
    use glide_vertint
    implicit none
    
    !TODO: Parameterize the following globals
    real(dp), parameter :: VEL2ERR  = 4e-3
    real(dp), parameter :: TOLER    = 1e-6
    integer,  parameter :: CONVT    = 4
    real(dp), parameter :: SHTUNE   = 1.D-16
    integer,  parameter :: UPSTREAM = 0
    integer,  parameter :: MANIFOLD = 1
contains
        subroutine setup_ismip_hom_a(surf, bed, thickness, delta_x, delta_y)
                double precision, dimension(:,:) :: surf, bed, thickness
                double precision, intent(in) :: delta_x, delta_y
                double precision :: m
                
                integer :: i, j, maxx, maxy
                
                maxx = size(surf, 2)
                maxy = size(surf, 1)
                
                !set m=.5 degree downward slope
                m = -tan(0.5*pi/180.0)
                write(*,*) delta_x, delta_y
                write(*,*)"m=",m 
                do i = 1, maxy
                        do j = 1, maxx
                    surf(i,j)=2600.+m*(j-1)*delta_x
                    thickness(i,j)=1000-500.*sin(2.*PI*(j-1)/(MAXX-2))*sin(2.*PI*(i-1)/(MAXY-2))
                    bed(i,j)=surf(i,j)-thickness(i,j)
                        end do
                end do
                
        end subroutine setup_ismip_hom_a
        
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

    subroutine velo_hom_pattyn(ewn, nsn, upn, dew, dns, sigma, &
                               thck, usrf, dthckdew, dthckdns, dusrfdew, dusrfdns, &
                               dlsrfdew, dlsrfdns, stagthck, mask, totpts, flwa, flwn, btrc, which_sliding_law, &
                               periodic_ew, periodic_ns, &
                               uvel, vvel, valid_initial_guess, uflx, vflx, efvs, tau, gdsx, gdsy)
                            
        integer  :: ewn !*FD Number of cells X
        integer  :: nsn !*FD Number of cells Y
        integer  :: upn !*FD Number of cells Z
        real(dp) :: dew !*FD Grid spacing X
	real(dp) :: dns !*FD Grid spacing Y
	real(dp), dimension(upn) :: sigma !*FD Sigma coord for rescaled Z dimension
	real(dp), dimension(ewn,nsn) :: thck !*FD Thickness, on non-staggered grid
	real(dp), dimension(ewn,nsn) :: usrf !*FD Upper surface profile
	real(dp), dimension(ewn-1,nsn-1) :: dthckdew !*FD X thickness gradient
	real(dp), dimension(ewn-1,nsn-1) :: dthckdns !*FD Y thickness gradient
	real(dp), dimension(ewn-1,nsn-1) :: dusrfdew !*FD X surface gradient
	real(dp), dimension(ewn-1,nsn-1) :: dusrfdns !*FD Y surface gradient
	real(dp), dimension(ewn-1,nsn-1) :: dlsrfdew !*FD X bed gradient
	real(dp), dimension(ewn-1,nsn-1) :: dlsrfdns !*FD Y bed gradient
	real(dp), dimension(ewn-1,nsn-1) :: stagthck !*FD Staggered thickness
        integer,  dimension(ewn-1,nsn-1) :: mask     !*FD Numbers points in the staggered grid that are included in computation
        integer :: totpts
	real(dp), dimension(:,:,:) :: flwa !*FD Glen's A (rate factor) - Used for thermomechanical coupling
        real(dp), dimension(ewn-1,nsn-1)   :: btrc !*FD Basal Traction, either betasquared or tau0
        real(dp) :: flwn !*FD Exponent in Glenn power law
        integer :: which_sliding_law
        logical :: periodic_ew !*Whether to use periodic boundary conditions
        logical :: periodic_ns
        real(dp), dimension(upn,ewn-1,nsn-1) :: uvel 
        real(dp), dimension(upn,ewn-1,nsn-1) :: vvel
        logical :: valid_initial_guess !*Whether or not the given uvel or vvel are appropriate initial guesses.  If not we'll have to roll our own.
        real(dp), dimension(:,:,:) :: uflx
        real(dp), dimension(:,:,:) :: vflx
        real(dp), dimension(:,:,:) :: efvs !*FD Effective viscosity
        type(glide_tensor)         :: tau
        real(dp), dimension(:,:) :: gdsx !*FD X driving stress
        real(dp), dimension(:,:) :: gdsy !*FD Y driving stress
        integer :: i, k

        !Second derivative of surface
        real(dp), dimension(ewn-1,nsn-1) :: d2zdx2, d2zdy2, d2hdx2, d2hdy2

        !Arrays for rescaled coordinate parameters
        real(dp), dimension(nsn-1, ewn-1, upn) :: ax, ay, bx, by, cxy
        !Arrays for staggered versions of quantities that were passed unstaggered
        real(dp), dimension(ewn-1, nsn-1) :: stagusrf, staglsrf
        
        real(dp), dimension(nsn, ewn)::u,v

        !Total number of grid points in 3 dimensions
        integer :: ijktot

        !Arrays to hold transposed data
        real(dp), dimension(nsn-1, ewn-1) :: dlsrfdew_t, dlsrfdns_t, dusrfdew_t, & 
        dusrfdns_t, dthckdew_t, dthckdns_t, btrc_t, staglsrf_t, stagusrf_t, stagthck_t, &
        d2zdx2_t, d2zdy2_t, d2hdx2_t, d2hdy2_t
       
        real(dp), dimension(nsn, ewn, upn) :: flwa_t
        real(dp), dimension(nsn-1, ewn-1, upn) :: uvel_t, vvel_t, mu_t, flwa_t_stag, &
                                                  tau_xz_t, tau_yz_t, tau_xx_t, tau_yy_t, tau_xy_t

        integer, dimension(nsn-1, ewn-1) :: mask_t
        
        ijktot = (nsn-1)*(ewn-1)*upn

        !Put the surface, bed, and thickness onto staggered grids
        !We do this because Ice3d is by nature unstaggered.
        !Note that we are already passed the staggered thickness
        call staggered_field_2d(usrf, stagusrf)
        staglsrf = stagusrf - stagthck

        !Compute second derivatives of thickness and surface, these are needed
        !for the rescaled coordinate parameters
        call d2f_field_stag(usrf, dew, dns, d2zdx2, d2zdy2, .false., .false.)
        call d2f_field_stag(thck, dew, dns, d2hdx2, d2hdy2, .false., .false.)
        
        !Move everything into transposed coordinates
        !This means that we need to transpose every data field with x and y that
        !we have thus far.
        !It also means that we need to swap x and y derivatives
        stagthck_t = transpose(stagthck)
        staglsrf_t = transpose(staglsrf)
        stagusrf_t = transpose(stagusrf)
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

        mask_t = transpose(mask)

        call glimToIce3d_3d(uvel,uvel_t,ewn-1,nsn-1,upn)
        call glimToIce3d_3d(vvel,vvel_t,ewn-1,nsn-1,upn)

        !The flow law field needs to be transposed AND staggered, which is a lot
        !of fun!!!
        call glimToIce3d_3d(flwa,flwa_t,ewn,nsn,upn)
        call staggered_field_3d(flwa_t, flwa_t_stag)
      
        !DEBUG!!!!
        call setup_ismip_hom_a(stagusrf_t, staglsrf_t, stagthck_t,dns,dew)
        call df_field_2d(stagusrf_t, dns, dew, dusrfdns_t, dusrfdew_t, .false., .false.)
        call df_field_2d(staglsrf_t, dns, dew, dlsrfdns_t, dlsrfdew_t, .false., .false.)
        call df_field_2d(stagthck_t, dns, dew, dthckdns_t, dthckdew_t, .false., .false.)

        !Compute rescaled coordinate parameters (needed because Pattyn uses an
        !irregular Z grid and scales so that 0 is the surface, 1 is the bed)
        call init_rescaled_coordinates(dthckdew_t,dlsrfdew_t,dthckdns_t,dlsrfdns_t,stagusrf_t,stagthck_t,staglsrf_t,&
                                               dusrfdew_t,dusrfdns_t,d2zdx2_t,d2zdy2_t,d2hdx2_t,d2hdy2_t,&
                                               sigma,ax,ay,bx,by,cxy,dew,dns)
       
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
                        TOLER, periodic_ew,periodic_ns, 0, dew, dns,mask_t,totpts)
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
                    TOLER, periodic_ew,periodic_ns, which_sliding_law, dew,&
                    dns,mask_t,totpts)
       
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
    subroutine velo_hom_pattyn_nonstag(ewn, nsn, upn, dew, dns, sigma, &
                               thck, usrf, lsrf, mask, totpts, flwa, flwn, btrc, which_sliding_law, &
                               periodic_ew, periodic_ns, &
                               uvel, vvel, valid_initial_guess)
                            
        integer  :: ewn !*FD Number of cells X
        integer  :: nsn !*FD Number of cells Y
        integer  :: upn !*FD Number of cells Z
        real(dp) :: dew !*FD Grid spacing X
	real(dp) :: dns !*FD Grid spacing Y
	real(dp), dimension(upn) :: sigma !*FD Sigma coord for rescaled Z dimension
	real(dp), dimension(ewn,nsn) :: thck !*FD Thickness, on non-staggered grid
	real(dp), dimension(ewn,nsn) :: usrf !*FD Upper surface profile
        real(dp), dimension(ewn,nsn) :: lsrf !*FD Lower surface profile
        integer,  dimension(ewn,nsn) :: mask     !*FD Numbers points in the staggered grid that are included in computation
        integer :: totpts
	real(dp), dimension(:,:,:) :: flwa !*FD Glen's A (rate factor) - Used for thermomechanical coupling
        real(dp), dimension(:,:)   :: btrc !*FD Basal Traction, either betasquared or tau0
        real(dp) :: flwn !*FD Exponent in Glenn power law
        integer :: which_sliding_law
        logical :: periodic_ew !*Whether to use periodic boundary conditions
        logical :: periodic_ns
        real(dp), dimension(upn,ewn-1,nsn-1) :: uvel 
        real(dp), dimension(upn,ewn-1,nsn-1) :: vvel
        logical :: valid_initial_guess !*Whether or not the given uvel or vvel are appropriate initial guesses.  If not we'll have to roll our own.
        integer :: i, k

        !Second derivative of surface
        real(dp), dimension(ewn,nsn) :: d2zdx2, d2zdy2, d2hdx2, d2hdy2

        !Arrays for rescaled coordinate parameters
        real(dp), dimension(nsn, ewn, upn) :: ax, ay, bx, by, cxy
        
        real(dp), dimension(nsn, ewn)::u,v

        !Total number of grid points in 3 dimensions
        integer :: ijktot
        
        real(dp), dimension(ewn, nsn) :: dusrfdew, dusrfdns, dlsrfdew, dlsrfdns, dthckdew, dthckdns


        !Arrays to hold transposed data
        real(dp), dimension(nsn, ewn) :: dlsrfdew_t, dlsrfdns_t, dusrfdew_t, & 
        dusrfdns_t, dthckdew_t, dthckdns_t, btrc_t_unstag,&
        d2zdx2_t, d2zdy2_t, d2hdx2_t, d2hdy2_t,usrf_t, lsrf_t, thck_t

        real(dp), dimension(nsn-1, ewn-1) :: btrc_t 
        real(dp), dimension(nsn, ewn, upn) :: flwa_t
        real(dp), dimension(nsn-1, ewn-1, upn) :: uvel_t, vvel_t, flwa_t_stag
        real(dp), dimension(nsn, ewn, upn) :: mu_t
        real(dp), dimension(nsn, ewn, upn) :: uvel_t_unstag, vvel_t_unstag
        integer, dimension(nsn, ewn) :: mask_t
        
        ijktot = (nsn)*(ewn)*upn

        !call fudge_mask(mask, totpts)

        !Compute first derivatives of geometry
        call df_field_2d(usrf, dew, dns, dusrfdew, dusrfdns, .false., .false.)
        call df_field_2d(lsrf, dew, dns, dlsrfdew, dlsrfdns, .false., .false.)
        call df_field_2d(thck, dew, dns, dthckdew, dthckdns, .false., .false.)

        !Compute second derivatives of thickness and surface, these are needed
        !for the rescaled coordinate parameters
        call d2f_field(usrf, dew, dns, d2zdx2, d2zdy2, .false., .false.)
        call d2f_field(thck, dew, dns, d2hdx2, d2hdy2, .false., .false.)
        
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

        mask_t = transpose(mask)

        call write_xls_int("mask.txt",mask_t)
        
        call glimToIce3d_3d(uvel,uvel_t,ewn-1,nsn-1,upn)
        call glimToIce3d_3d(vvel,vvel_t,ewn-1,nsn-1,upn)
        call glimToIce3d_3d(flwa,flwa_t,ewn,nsn,upn)
      
        call unstagger_field_3d_periodic(vvel_t, uvel_t_unstag)
        call unstagger_field_3d_periodic(vvel_t, uvel_t_unstag)
        call unstagger_field_2d_periodic(btrc_t, btrc_t_unstag)
        write(*,*) ewn, dew, nsn, dns, upn
        write(*,*) shape(btrc)
        !Compute rescaled coordinate parameters (needed because Pattyn uses an
        !irregular Z grid and scales so that 0 is the surface, 1 is the bed)
        call init_rescaled_coordinates(dthckdew_t,dlsrfdew_t,dthckdns_t,dlsrfdns_t,usrf_t,thck_t,lsrf_t,&
                                               dusrfdew_t,dusrfdns_t,d2zdx2_t,d2zdy2_t,d2hdx2_t,d2hdy2_t,&
                                               sigma,ax,ay,bx,by,cxy,dew,dns)
       
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
                        TOLER, periodic_ew,periodic_ns, 0, dew, dns,mask_t,totpts)
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
                    TOLER, periodic_ew,periodic_ns, which_sliding_law, dew,&
                    dns,mask_t,totpts)
       
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
        integer :: nx, ny, nz, ni, nj, nk, i, j, k
        logical :: to_glimmer
 
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
        integer :: nx, ny, nz, ni, nj, nk, i, j, k
        logical :: to_glimmer

        do i=1,nx
            do j=1,ny
                do k=1,nz
                    dest(k,i,j) = field(j,i,k)
                end do
            end do
        end do
    end subroutine ice3dToGlim_3d

    subroutine fudge_mask(mask, totpts)
        integer, dimension(:,:) :: mask
        integer, intent(out) :: totpts

        integer :: i, j, ni, nj,k

        ni = size(mask,1)
        nj = size(mask,2)
        k = 1
        totpts = ni*nj

        do i=1,ni
            do j=1,nj
                mask(i,j) = k
                k = k+1
            end do
        end do
    end subroutine

end module glide_velo_higher
