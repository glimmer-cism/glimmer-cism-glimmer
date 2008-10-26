#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

#define print_size(var) write(*,*)"size of var", shape(var)
#define assert(expr) call assert_impl("expr",(expr),__FILE__, __LINE__)

module glide_velo_higher
    use ice3d_lib
    use glimmer_global, only : dp
    use glide_types
    implicit none
    
    !TODO: Parameterize the following globals
    real(dp), parameter :: FLOWN    = 3. !TODO: This is passed in, but as an array
    real(dp), parameter :: VEL2ERR  = 4e-3
    real(dp), parameter :: TOLER    = 1e-6
    integer,  parameter :: CONVT    = 4
    real(dp), parameter :: SHTUNE   = 1.D-16
    integer,  parameter :: UPSTREAM = 0
    integer,  parameter :: MANIFOLD = 1
contains
    subroutine assert_impl(cond, res, filename, line)
        character(256) :: cond
        character(64) :: filename
        integer :: line
        logical :: res

        if (.not.res) then
            write(*,*) "Assertion Failed:", cond, "at", filename, ":", line
            stop
        end if
    end subroutine
    
    subroutine velo_hom_pattyn(ewn, nsn, upn, dew, dns, sigma, &
                               thck, usrf, dthckdew, dthckdns, dusrfdew, dusrfdns, &
                               dlsrfdew, dlsrfdns, stagthck, flwa, btrc, which_sliding_law, &
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
	real(dp), dimension(:,:,:) :: flwa !*FD Glen's A (rate factor)
        real(dp), dimension(ewn-1,nsn-1)   :: btrc !*FD Basal Traction, either betasquared or tau0
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
        real(dp), dimension(nsn-1, ewn-1, upn) :: ax, ay, bx, by, cxy,arrh
        !Arrays for staggered versions of quantities that were passed unstaggered
        real(dp), dimension(ewn-1, nsn-1) :: stagusrf, staglsrf
        
        real(dp), dimension(nsn, ewn)::u,v

        !Total number of grid points in 3 dimensions
        integer :: ijktot

        !Arrays to hold transposed data
        real(dp), dimension(nsn-1, ewn-1) :: dlsrfdew_t, dlsrfdns_t, dusrfdew_t, & 
        dusrfdns_t, dthckdew_t, dthckdns_t, btrc_t, staglsrf_t, stagusrf_t, stagthck_t, &
        d2zdx2_t, d2zdy2_t, d2hdx2_t, d2hdy2_t
       
        real(dp), dimension(nsn-1, ewn-1, upn) :: uvel_t, vvel_t, mu_t

        write(*,*)"BEGIN velo_hom_pattyn"
        write(*,*)"Sliding law = ", which_sliding_law
        ijktot = (nsn-1)*(ewn-1)*upn
        arrh = SHTUNE

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
        call glimToIce3d_3d(uvel,uvel_t,ewn-1,nsn-1,upn-1)
        call glimToIce3d_3d(vvel,vvel_t,ewn-1,nsn-1,upn-1)
       
        !Compute rescaled coordinate parameters (needed because Pattyn uses an
        !irregular Z grid and scales so that 0 is the surface, 1 is the bed)
        call init_rescaled_coordinates(dthckdew_t,dlsrfdew_t,dthckdns_t,dlsrfdns_t,stagusrf_t,stagthck_t,staglsrf_t,&
                                               dusrfdew_t,dusrfdns_t,d2zdx2_t,d2zdy2_t,d2hdx2_t,d2hdy2_t,&
                                               sigma,ax,ay,bx,by,cxy,dew,dns)
        
        !"Spin up" estimate with Pattyn's SIA model runs if we don't already
        !have a good initial guess
        if (.not. valid_initial_guess) then
            write(*,*)"SIA spinup"
            call veloc1(dusrfdew_t, dusrfdns_t, stagthck_t, arrh, sigma, uvel_t, vvel_t, u, v, nsn-1, ewn-1, upn, &
                        FLOWN, periodic_ew, periodic_ns)
            !If we are performing the plastic bed iteration, the SIA is not
            !enough and we need to spin up a better estimate by shoehorning the
            !tau0 values into a linear bed estimate
            if (which_sliding_law == 1) then
                write(*,*)"Linear bed spinup"
                call veloc2(mu_t, uvel_t, vvel_t, arrh, dusrfdew_t, dusrfdns_t, stagthck_t, ax, ay, &
                        sigma, bx, by, cxy, btrc_t/100, dlsrfdew_t, dlsrfdns_t, FLOWN, ZIP, VEL2ERR, MANIFOLD,&
                        TOLER, periodic_ew,periodic_ns, 0, dew, dns)
            end if
            write(*,*)"Spinup complete"
        end if
        !Higher order velocity estimation
        !I am assuming that efvs (effective viscosity) is the same as mu
        !A NOTE ON COORDINATE TRANSPOSITION:
        !Because of the transposition, ewn=maxy and nsn=maxx.  However, veloc2
        !passes maxy in *first*, so these really get passed in the same order
        !that they normally would.
        call veloc2(mu_t, uvel_t, vvel_t, arrh, dusrfdew_t, dusrfdns_t, stagthck_t, ax, ay, &
                    sigma, bx, by, cxy, btrc_t, dlsrfdew_t, dlsrfdns_t, FLOWN, ZIP, VEL2ERR, MANIFOLD,&
                    TOLER, periodic_ew,periodic_ns, which_sliding_law, dew, dns)
        call write_xls_3d("uvel_result.txt",uvel_t)
        call write_xls_3d("vvel_result.txt",vvel_t)
        !Transpose from the ice3d coordinate system (y,x,z) to the glimmer
        !coordinate system (z,x,y).  We need to do this for all the 3D outputs
        !that Glimmer expects.
        call ice3dToGlim_3d(uvel_t, uvel, ewn-1, nsn-1, upn-1)
        call ice3dToGlim_3d(vvel_t, vvel, ewn-1, nsn-1, upn-1)
        call ice3dToGlim_3d(mu_t, efvs, ewn-1, nsn-1, upn-1)

        !TODO: Stress field computation and output

        write(*,*)"END velo_hom_pattyn"
    end subroutine velo_hom_pattyn

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

end module glide_velo_higher
