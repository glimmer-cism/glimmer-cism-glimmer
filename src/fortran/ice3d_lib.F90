!--------------------------------------------------
!ice3d_lib.F90
!This file is adapted from Frank Pattyn's Ice3D model
!Code updated to Fortran 90 and integrated with Glimmer
!by Tim Bocek.
!See F. Pattyn, "A new three-dimensional higher-order thermomechanical
!        ice sheet model: basic sensitivity, ice stream development,
!        and ice flow across subglacial lakes", Journal of Geophysical
!        Research, Volume 108, no. B8, 2003.
!----------------------------------------------------

#include "glide_nan.inc"
#include "glide_mask.inc"

!The following defines contain debug options that aid in examining the sparse
!matrix.  They should be disabled unless the higher-order code needs debugging.

!Define to write the rows of the sparse matrix for each vel. component then stop.
!#define OUTPUT_SPARSE_MATRIX

!Define to output a NetCDF file of the partial iterations
#define OUTPUT_PARTIAL_ITERATIONS

!If defined, a vertically averaged pressure will be used across the ice front.
!Otherwise, a vertically explicit pressure will be used.
!#define AVERAGED_PRESSURE

!The maximum number of iterations to use in the unstable manifold loop
#define NUMBER_OF_ITERATIONS 1000

!#define VERY_VERBOSE

module ice3d_lib
    use glimmer_global
    use glimmer_physcon, only: pi, grav, rhoi, rhoo, scyr
    use glide_deriv
    use glide_types
    use glimmer_sparse_type
    use glimmer_sparse
    use glimmer_log
    use xls
    use glide_nonlin
    implicit none
    real(dp) :: small, zip, toler_adjust_factor
    PARAMETER(SMALL=1.D-10,ZIP=1.D-30)
    PARAMETER(toler_adjust_factor=1.1)

    real(dp) :: plastic_bed_regularization = 1e-2

    integer, dimension(:,:), allocatable :: pointtype

#ifdef VERY_VERBOSE
    logical, parameter :: sparverbose = .true.
#else
    logical, parameter :: sparverbose = .false.
#endif
    !Mu is recomputed only when the error tolerance reaches this level
    !This allows nonlinearities resulting from equation separation and 
    !the beta^2 computation to "quiet down" before attempting another
    !viscosity iteration.  Set to a very high number to turn this feature off
    !entirely
    real(dp) :: recompute_mu_toler = 1e10
    real(dp) :: recompute_mu_adjust = 1

    integer, parameter :: ITER_UNIT = 42
!------------------------------------------------------
!   lookup coordinates in sparse matrix
!
!   Reference:
!     1   ->   i - 1, j - 1, k
!x    2   ->   i - 1, j, k - 1
!x    3   ->   i - 1, j, k
!x    4   ->   i - 1, j, k + 1
!     5   ->   i - 1, j + 1, k
!x    6   ->   i, j - 1, k - 1
!x    7   ->   i, j - 1, k
!x    8   ->   i, j - 1, k + 1
!x    9   ->   i, j, k - 2
!x   10   ->   i, j, k - 1
!x   11   ->   i, j, k
!x   12   ->   i, j, k + 1
!x   13   ->   i, j, k + 2
!x   14   ->   i, j + 1, k - 1
!x   15   ->   i, j + 1, k
!x   16   ->   i, j + 1, k + 1
!    17   ->   i + 1, j - 1, k
!x   18   ->   i + 1, j, k - 1
!x   19   ->   i + 1, j, k
!x   20   ->   i + 1, j, k + 1
!    21   ->   i + 1, j + 1, k
!------------------------------------------------------
    integer, parameter :: IM1_JM1_K = 1
    integer, parameter :: IM1_J_KM1 = 2
    integer, parameter :: IM1_J_K = 3
    integer, parameter :: IM1_J_KP1 = 4
    integer, parameter :: IM1_JP1_K = 5
    integer, parameter :: I_JM1_KM1 = 6
    integer, parameter :: I_JM1_K = 7
    integer, parameter :: I_JM1_KP1 = 8
    integer, parameter :: I_J_KM2 = 9
    integer, parameter :: I_J_KM1 = 10
    integer, parameter :: I_J_K = 11
    integer, parameter :: I_J_KP1 = 12
    integer, parameter :: I_J_KP2 = 13
    integer, parameter :: I_JP1_KM1 = 14
    integer, parameter :: I_JP1_K = 15
    integer, parameter :: I_JP1_KP1 = 16
    integer, parameter :: IP1_JM1_K = 17
    integer, parameter :: IP1_J_KM1 = 18
    integer, parameter :: IP1_J_K = 19
    integer, parameter :: IP1_J_KP1 = 20
    integer, parameter :: IP1_JP1_K = 21

    !Extra stencil locations introduced to allow for 2nd order upwinding and downwinding
    !lateral derivs.
    integer, parameter :: IM2_J_K = 22
    integer, parameter :: IP2_J_K = 23
    integer, parameter :: I_JM2_K = 24
    integer, parameter :: I_JP2_K = 25

    integer, parameter :: STENCIL_SIZE = 25
contains
    !Initialize rescaled coordinate coefficients
    subroutine init_rescaled_coordinates(dhdx,dhbdx,dhdy,dhbdy,surf,h,hb,&
               dzdx,dzdy,d2zdx2, d2zdy2, d2hdx2, d2hdy2, zeta,ax,ay,bx,by,cxy,dx,dy, &
               direction_x, direction_y)
!
        INTEGER MAXX,MAXY,NZETA
        real(dp), dimension(:,:), intent(in) :: dhdx
        real(dp), dimension(:,:), intent(in) :: dhbdx
        real(dp), dimension(:,:), intent(in) :: dhdy
        real(dp), dimension(:,:), intent(in) :: dhbdy
        real(dp), dimension(:,:), intent(in) :: surf
        real(dp), dimension(:,:), intent(in) :: h
        real(dp), dimension(:,:), intent(in) :: hb
        real(dp), dimension(:,:), intent(in) :: dzdx
        real(dp), dimension(:,:), intent(in) :: dzdy
        real(dp), dimension(:),   intent(in) :: zeta
        real(dp), dimension(:,:,:), intent(inout) :: ax
        real(dp), dimension(:,:,:), intent(inout) :: ay
        real(dp), dimension(:,:,:), intent(inout) :: bx
        real(dp), dimension(:,:,:), intent(inout) :: by
        real(dp), dimension(:,:,:), intent(inout) :: cxy
        
        real(dp), intent(in) :: dx
        real(dp), intent(in) :: dy
        
        real(dp), dimension(:,:), intent(in) :: d2zdx2,d2zdy2,d2hdx2,d2hdy2
        
        real(dp), dimension(:,:), intent(in) :: direction_x, direction_y

        INTEGER :: i,j,k

#ifndef NO_RESCALE
        write(*,*) "Pattyn-Bocek model does not work with scaling yet.  Please re-compile with the flag -DNO_RESCALE."
        stop
#endif

#if 0
        call write_xls("h.txt",h)
        call write_xls("hb.txt",hb)
        call write_xls("surf.txt",surf)
        call write_xls("dhdx.txt",dhdx)
        call write_xls("dhdy.txt",dhdy)
        call write_xls("dhbdx.txt",dhbdx)
        call write_xls("dhbdy.txt",dhbdy)
        call write_xls("dzdx.txt",dzdx)
        call write_xls("dzdy.txt",dzdy)
        call write_xls("d2zdx2.txt",d2zdx2)
        call write_xls("d2zdy2.txt",d2zdy2)
        call write_xls("d2hdx2.txt",d2hdx2)
        call write_xls("d2hdy2.txt",d2hdy2)
#endif
        MAXY = size(h, 1)
        MAXX = size(h, 2)
        NZETA = size(zeta, 1)

        !Get a field of 2nd derivatives of the surface onto a staggered
        !grid
        !Determination of scaling factors ax, ay, bx, by, cxy
        !These are given by (40) through (42) in Pattyn 2002.
        !At boundaries or when ice thickness (h) is 0, these
        !factors are set to 0.
        !WARNING: In this block of code, dfdx really means dfdy
        !and vice versa.  This is because it is operating in
        !transposed coordinates instead of native Fortran
        !coordinates.
        do i = 1, MAXY
            do j = 1, MAXX
                if ( (i > 1) .and. (i < MAXY) .and. (j > 1) .and. (j < MAXX)&
                     .and. (h(i,j) > 0.)) then
                    do k = 1, NZETA
                        ax(i,j,k) = (dzdx(i, j) - zeta(k) * dhdx(i, j)) / h(i, j)
                        ay(i,j,k) = (dzdy(i, j) - zeta(k) * dhdy(i, j)) / h(i, j)
              
                        bx(i,j,k) = ( d2zdx2(i, j) - &
                                      zeta(k) * d2hdx2(i, j) - & 
                                      2. * ax(i, j, k) * dhdx(i, j) ) / h(i,j)
                         
                        by(i,j,k) = ( d2zdy2(i, j) - &
                                      zeta(k) * d2hdy2(i, j) - & 
                                      2. * ay(i, j, k) * dhdy(i, j) ) / h(i,j)
                         
                        cxy(i,j,k) = ( dfdy_2d(dzdy, i, j, dx) - & 
                                      zeta(k) * dfdy_2d(dhdy, i, j, dx) - &
                                      ax(i, j, k)*dhdy(i, j) - &
                                      ay(i, j, k)*dhdx(i, j) ) / h(i,j)
                
                    end do
                else
                    ax(i,j,:)=0.
                    ay(i,j,:)=0.
                    bx(i,j,:)=0.
                    by(i,j,:)=0.
                    cxy(i,j,:)=0.
                endif
            end do
        end do
#if 0
        call write_xls_3d("ax.txt",ax)
        call write_xls_3d("ay.txt",ay)
        call write_xls_3d("bx.txt",bx)
        call write_xls_3d("by.txt",by)
        call write_xls_3d("cxy.txt",cxy)
#endif
    end subroutine
!
!
!------------------------------------------------------
!   Velocity estimation according to 0th order model (SIA)
!------------------------------------------------------
!
      SUBROUTINE veloc1(dzdx,dzdy,h,arrh,zeta,uvel,vvel,u,v,&
        MAXY,MAXX,NZETA,FLOWN,PERIODIC_X,PERIODIC_Y)
!
        INTEGER MAXX,MAXY,NZETA
        real(dp) :: dzdx(MAXY,MAXX)
        real(dp) :: dzdy(MAXY,MAXX),h(MAXY,MAXX),zeta(NZETA)
        real(dp) :: arrh(MAXY,MAXX,NZETA),uvel(MAXY,MAXX,NZETA)
        real(dp) :: vvel(MAXY,MAXX,NZETA),u(MAXY,MAXX)
        real(dp) :: v(MAXY,MAXX),FLOWN
        logical :: PERIODIC_X, PERIODIC_Y
!
      INTEGER i,j,k
      real(dp) :: diff1(NZETA),d,grad,z,diffus,us,vs
!
      diff1 = 0
     
      do 10 i=2,MAXY-1
        do 20 j=2,MAXX-1
          grad=(dzdx(i,j)**2)+(dzdy(i,j)**2)
          us=0.
          vs=0.
          d = (RHOI*GRAV*h(i,j))**FLOWN
          !Accumulated vertical integration using trapezoid rule??
          do 30 k=(NZETA-1),1,-1
            z = (0.5*(zeta(k+1)+zeta(k)))**FLOWN
            diff1(k)=d*grad*h(i,j)*(arrh(i,j,k+1)+arrh(i,j,k))*z*&
              (zeta(k)-zeta(k+1))+diff1(k+1)
   30     CONTINUE
          diffus=0.
          !Vertical averaging of diffusivity (this works because the
          !differences between successive vertical layers must sum to 1!)
          do 40 k=2,NZETA
            diffus=diffus+0.5*(diff1(k)+diff1(k-1))*&
              (zeta(k)-zeta(k-1))
   40     CONTINUE
          !Estimation of velocity from SIA diffusivity
          do 50 k=1,NZETA
            uvel(i,j,k)=diff1(k)*dzdx(i,j)+us
            vvel(i,j,k)=diff1(k)*dzdy(i,j)+vs
   50     CONTINUE
        !Estimation of vertical average velocity from vertical average diff.
        u(i,j)=diffus*dzdx(i,j)
        v(i,j)=diffus*dzdy(i,j)
   20   CONTINUE
   10 CONTINUE

      !Periodic boundary conditions
      if (PERIODIC_X) then
        do 100 i=1,MAXY
          u(i,MAXX)=u(i,2)
          u(i,1)=u(i,MAXX-1)
          v(i,MAXX)=v(i,2)
          v(i,1)=v(i,MAXX-1)
          do 110 k=1,NZETA
            uvel(i,MAXX,k)=uvel(i,2,k)
            uvel(i,1,k)=uvel(i,MAXX-1,k)
            vvel(i,MAXX,k)=vvel(i,2,k)
            vvel(i,1,k)=vvel(i,MAXX-1,k)
  110     CONTINUE
  100   CONTINUE
      end if
      if (PERIODIC_Y) then
        do 120 j=1,MAXX
          u(MAXY,j)=u(2,j)
          u(1,j)=u(MAXY-1,j)
          v(MAXY,j)=v(2,j)
          v(1,j)=v(MAXY-1,j)
          do 130 k=1,NZETA
            uvel(MAXY,j,k)=uvel(2,j,k)
            uvel(1,j,k)=uvel(MAXY-1,j,k)
            vvel(MAXY,j,k)=vvel(2,j,k)
            vvel(1,j,k)=vvel(MAXY-1,j,k)
  130     CONTINUE
  120   CONTINUE
      endif
      return
      END subroutine
     !

    subroutine periodic_boundaries(m, apply_to_x, apply_to_y)
        real(dp), dimension(:,:), intent(inout) :: m
        integer :: maxx, maxy
        logical :: apply_to_x, apply_to_y
        maxx = size(m, 2)
        maxy = size(m, 1)

        if(apply_to_x) then
            m(:,maxx) = m(:,2)
            m(:,1) = m(:,maxx-1)
        end if

        if(apply_to_y) then
            m(maxy,:) = m(2,:)
            m(1,:) = m(maxy-1,:)
        end if
    end subroutine periodic_boundaries
    
    subroutine periodic_boundaries_stag(m, apply_to_x, apply_to_y)
        real(dp), dimension(:,:), intent(inout) :: m
        integer :: maxx, maxy
        logical, intent(in) :: apply_to_x, apply_to_y

        maxx = size(m, 2)
        maxy = size(m, 1)
        
        !UNCOMMENT for "standard" periodic boundary conditions (this might
        !(not be right)
    
        !See the discussion of periodic boundary conditions on a staggered
        !grid at the top of this code file
        !Note that we aren't using the last row and column of the staggered grid
        !right now!!!
        if(apply_to_x) then
            m(:,maxx) = m(:,2)
            m(:,1) = m(:,maxx-1)
        end if

        if(apply_to_y) then
            m(maxy,:) = m(2,:)
            m(1,:) = m(maxy-1,:)
        end if

    end subroutine periodic_boundaries_stag
    
    subroutine periodic_boundaries_3d(m, apply_to_x, apply_to_y)
        real(dp), dimension(:,:,:) :: m
        integer :: maxx, maxy, maxz , k
        logical, intent(in) :: apply_to_x, apply_to_y
                
        maxx = size(m,2)
        maxy = size(m,1)
        maxz = size(m,3)
                
    
        do k = 1, maxz
            call periodic_boundaries(m(:,:,k),apply_to_x, apply_to_y)
        end do    
    end subroutine periodic_boundaries_3d
    
    subroutine periodic_boundaries_3d_stag(m,apply_to_x,apply_to_y)
        real(dp), dimension(:,:,:) :: m
        integer :: k
        logical, intent(in) :: apply_to_x, apply_to_y
    
        do k=1, size(m, 3)
            call periodic_boundaries_stag(m(:,:,k),apply_to_x, apply_to_y)
        end do
    end subroutine periodic_boundaries_3d_stag
!
!------------------------------------------------------
!   Velocity estimation according to higher order model
!------------------------------------------------------
!
    SUBROUTINE veloc2(mu,uvel,vvel,arrh,dzdx,dzdy,h,ax,ay,&
                zeta,bx,by,cxy,beta,&
                dhbdx,dhbdy,FLOWN,ZIP,VEL2ERR,&
                TOLER, options, STAGGERED, delta_x, delta_y, &
                point_mask, active_points, geometry_mask, kinematic_bc_u, kinematic_bc_v, &
                marine_bc_normal)
                

        INTEGER MAXY,MAXX,NZETA,MANIFOLD

        real(dp), dimension(:,:,:), intent(inout) :: mu
        real(dp), dimension(:,:,:), intent(inout) :: uvel
        real(dp), dimension(:,:,:), intent(inout) :: vvel
        real(dp), dimension(:,:,:), intent(in) :: arrh
        real(dp), dimension(:,:), intent(in) :: dzdx
        real(dp), dimension(:,:), intent(in) :: dzdy
        real(dp), dimension(:,:), intent(in) :: h
        real(dp), dimension(:,:,:), intent(in) :: ax
        real(dp), dimension(:,:,:), intent(in) :: ay
        real(dp), dimension(:,:,:), intent(in) :: bx
        real(dp), dimension(:,:,:), intent(in) :: by
        real(dp), dimension(:,:,:), intent(in) :: cxy
        real(dp), dimension(:), intent(in) :: zeta
        real(dp), dimension(:,:), intent(in) :: dhbdx
        real(dp), dimension(:,:), intent(in) :: dhbdy
        real(dp), dimension(:,:), intent(in) :: beta
        type(glide_options), intent(in) :: options 
        real(dp), intent(in) :: FLOWN,ZIP,VEL2ERR,TOLER,delta_x, delta_y
        logical, intent(in) :: STAGGERED !Whether the model is run on a staggered grid or a colocated grid.
                                         !This is passed so that corrections arising from averaging can be made.
        
        integer, dimension(:,:), intent(in) :: point_mask 
        !*FD Contains a unique nonzero number on each point of ice that should be computed or that is
        !*FD ajdacent to a point that should be computed. Other numbers contain
        !*FD zeros
        integer, intent(in) :: active_points !*FD The number of points that should be computed.  This is equivalent to the number of nonzero entries in point_mask.

        integer, intent(in), dimension(:,:) :: geometry_mask !*FD point_mask field as described in glide_mask.inc

        !Contains NaN everywhere except where a kinematic boundary is to be
        !applied, in which case contains the value at the boundary
        real(dp), dimension(:,:,:), intent(inout) :: kinematic_bc_u, kinematic_bc_v
        
        !Contains NaN everywhere except on the marine ice edge of an ice shelf,
        !where this should contain the angle of the normal to the marine edge,
        !in radians, with 0=12 o'clock, pi/2=3 o'clock, etc.
        real(dp), dimension(:,:), intent(in)   :: marine_bc_normal

        INTEGER l,lacc,m,maxiter,iter,DU1,DU2,DV1,DV2,ijktot,k
        real(dp) :: error,tot,teta

        PARAMETER (DU1=1,DU2=2,DV1=3,DV2=4)

        real(dp), dimension(4,size(mu,1)*size(mu,2)*size(mu,3)) :: em

        real(dp), dimension(2*size(mu,1)*size(mu,2)*size(mu,3)) :: correction_vec
        !Velocity estimates computed for the *current* iteration.  uvel and
        !vvel, comparitively, hold the velocity estimates for the *last*
        !iteration.
        real(dp), dimension(size(mu,1),size(mu,2),size(mu,3)) :: ustar, vstar
        real(dp), dimension(size(mu,1),size(mu,2)) :: tau !Basal traction, to be computed from the provided beta parameter
       
        real(dp), dimension(size(mu,1),size(mu,2),size(mu,3)) :: dudx, dudy, dudz
        real(dp), dimension(size(mu,1),size(mu,2),size(mu,3)) :: dvdx, dvdy, dvdz
        real(dp), dimension(size(mu,1),size(mu,2)) :: direction_x, direction_y

        logical :: cont

        type(sparse_matrix_type) :: matrix
        type(sparse_solver_workspace) :: matrix_workspace
        type(sparse_solver_options) :: matrix_options

#ifdef OUTPUT_PARTIAL_ITERATIONS 
        integer :: ncid_debug
#endif
        !For timing the algorithm
        real(dp) :: solve_start_time, solve_end_time, iter_start_time, iter_end_time

        call cpu_time(solve_start_time)

        allocate(pointtype(size(mu, 1), size(mu, 2)))

        !Get the sizes from the mu field.  Calling code must make sure that
        !these all agree.
        maxy = size(mu,1)
        maxx = size(mu,2)
        nzeta = size(mu,3)
        ijktot = active_points*nzeta

        uvel = 0
        vvel = 0

        maxiter=NUMBER_OF_ITERATIONS
        error=VEL2ERR
        m=1
        lacc=0
        em = 0
        correction_vec = 0
        tot = 0
        !Set up sparse matrix options
        call sparse_solver_default_options(options%which_ho_sparse, matrix_options)
        matrix_options%base%tolerance=TOLER
        matrix_options%base%maxiters  = 100000
        
        if (options%which_ho_efvs == 1) then
            write(*,*) "Using linear rheology"
        else
            write(*,*) "Using full stress calculation"
        end if

        !Create the sparse matrix
        call new_sparse_matrix(ijktot, ijktot*STENCIL_SIZE, matrix)
        call sparse_allocate_workspace(matrix, matrix_options, matrix_workspace, ijktot*STENCIL_SIZE)

#if 1
        call write_xls_3d("arrh.txt",arrh)
        call write_xls("dzdx.txt",dzdx)
        call write_xls("dzdy.txt",dzdy)
        call write_xls_3d("ax.txt",ax)
        call write_xls_3d("ay.txt",ay)
        call write_xls_3d("bx.txt",bx)
        call write_xls_3d("by.txt",by)
        call write_xls_3d("cxy.txt",cxy)
        call write_xls("h.txt",h)
        call write_xls_3d("uvel_sia.txt",uvel)
        call write_xls_3d("vvel_sia.txt",vvel)
        call write_xls("beta.txt",beta)
        call write_xls("dhbdx.txt",dhbdx)
        call write_xls("dhbdy.txt",dhbdy)
        call write_xls("latbc.txt",marine_bc_normal)
        call write_xls_int("geometry_mask.txt",geometry_mask)
        call write_xls_3d("kinematic_bc_u.txt",kinematic_bc_u)
        call write_xls_3d("kinematic_bc_v.txt",kinematic_bc_v)
        call write_xls("marine_bc_norms.txt",marine_bc_normal)
        call write_xls_direction_guide("direction_guide.txt",3,3)
        call write_xls("normal_x.txt", sin(marine_bc_normal))
        call write_xls("normal_y.txt", -cos(marine_bc_normal))
        write(*,*) "ZETA=",zeta
#endif

        !Copy the velocity estimate from the previous iteration as the current velocity
        ustar=uvel
        vstar=vvel

        !Deal with periodic boundary conditions
        call periodic_boundaries_3d_stag(uvel,options%periodic_ew,options%periodic_ns)
        call periodic_boundaries_3d_stag(vvel,options%periodic_ew,options%periodic_ns)
        !Compute basal traction

        !Force a radially symmetric initial guess
        !call fill_radial_symmetry(uvel, 1000.0_dp, 1)
        !call fill_radial_symmetry(vvel, 1000.0_dp, 2)
        
        call comp_type(kinematic_bc_u(:,:,1), geometry_mask, point_mask, "comp_type_before.txt")
        !call mask_out_edges(geometry_mask, point_mask, kinematic_bc_u, kinematic_bc_v)
        !do k=1,nzeta
        !    where (h > 0 .and. point_mask == 0)
        !        kinematic_bc_u(:,:,k) = 0
        !        kinematic_bc_v(:,:,k) = 0
        !    endwhere
        !end do
        call comp_type(kinematic_bc_u(:,:,1), geometry_mask, point_mask, "comp_type_after.txt")

#ifdef OUTPUT_PARTIAL_ITERATIONS        
        ncid_debug = begin_iteration_debug(maxx, maxy, nzeta)
        call iteration_debug_step(ncid_debug, 0, mu, uvel, vvel, geometry_mask)
#endif
        nonlinear_iteration: do l=1,maxiter !Loop until we have reached the number of iterations allowed
            call cpu_time(iter_start_time)
            
            lacc=lacc+1
           
            !Every 10 iterations, we relax the error requirement somewhat (1/2
            !an order of magnitude) under the assumption that we won't converge 
            !given the current error tolerance.
            !(TODO: This leads to exponential increase in the tolerance - can we
            !get by without it?  Is this too liberal of a relaxation?)
            if (l == m*10) then
                error = error*toler_adjust_factor
#if DEBUG
                write(*,*) "Error tolerance is now", error
#endif
                m = m+1
            endif
             
            !Compute the basal traction.  This is done by either passing through
            !the beta parameter, or by using Shoof's plastic bed model (Schoof
            !2004).
            if (options%which_ho_bstress == HO_BSTRESS_LINEAR) then
                tau = beta
            else
                call plastic_bed(tau, beta, uvel(:,:,nzeta), vvel(:,:,nzeta))
            end if

            !Compute velocity derivatives
            call velderivs(uvel, vvel, dzdx, dzdy, geometry_mask, &
                           delta_x, delta_y, zeta, .false., &
                           direction_x, direction_y, &
                           dudx, dudy, dudz, dvdx, dvdy, dvdz)
#if 0
            call write_xls_3d("dudx.txt",dudx)
            call write_xls_3d("dudy.txt",dudy)
            call write_xls_3d("dudz.txt",dudz)
            call write_xls_3d("dvdx.txt",dvdx)
            call write_xls_3d("dvdy.txt",dvdy)
            call write_xls_3d("dvdz.txt",dvdz)
#endif
            !Compute viscosity
            if (options%which_ho_efvs == HO_EFVS_FULL) then
                if (tot < recompute_mu_toler) then
                    call muterm(mu,arrh,h,ax,ay,FLOWN,ZIP,dudx, dudy, dudz, dvdx, dvdy, dvdz)
                    recompute_mu_toler = recompute_mu_toler / recompute_mu_adjust
                    recompute_mu_toler = max(recompute_mu_toler, error*2)
                end if
            else if (options%which_ho_efvs == HO_EFVS_CONSTANT) then
                mu = 1d6
            end if

            call write_xls_3d("mu.txt",mu)
            !Apply periodic boundary conditions to the viscosity
            call periodic_boundaries_3d_stag(mu,options%periodic_ew,options%periodic_ns)
            !call write_xls_3d("mu.txt",mu)
            
            !Sparse matrix routine for determining velocities.  The new
            !velocities will get spit into ustar and vstar, while uvel and vvel
            !will still hold the old velocities.

#ifdef OUTPUT_SPARSE_MATRIX
            call open_sparse_output(l, ITER_UNIT)
#endif
            iter=sparuv(mu,dzdx,dzdy,ax,ay,bx,by,cxy,h,&
                uvel,vvel,dudx,dudy,dudz,dvdx,dvdy,dvdz,&
                ustar,vstar,tau,dhbdx,dhbdy,ijktot,MAXY,&
                MAXX,NZETA,TOLER, options%which_ho_source, delta_x, delta_y, zeta, point_mask, &
                geometry_mask,matrix, matrix_workspace, matrix_options, &
                kinematic_bc_u, kinematic_bc_v, &
                marine_bc_normal,direction_x,direction_y, STAGGERED, options%ho_include_thinice)

#ifdef OUTPUT_SPARSE_MATRIX
           close(ITER_UNIT)
#endif

            !Apply periodic boundary conditions to the computed velocity
            call periodic_boundaries_3d_stag(ustar,options%periodic_ew,options%periodic_ns)
            call periodic_boundaries_3d_stag(vstar,options%periodic_ew,options%periodic_ns)

#ifdef OUTPUT_PARTIAL_ITERATIONS
            call iteration_debug_step(ncid_debug, l, mu, ustar, vstar, geometry_mask)
            call iterdebug_vel_derivs(ncid_debug, l, dudx, dudy, dvdx, dvdy)
#endif

            !Apply unstable manifold correction.  This function returns
            !true if we need to keep iterating, false if we reached convergence
            cont = umc_correct_vels(ustar, uvel, &
                                                vstar, vvel, correction_vec, &
                                                maxy, maxx, nzeta, error, &
                                                tot, teta)    
            call cpu_time(iter_end_time)
#if DEBUG 
            write(*,*) l, iter, tot, teta, &
                       max( maxval(abs(ustar)), maxval(abs(vstar)) ),&
                       iter_end_time - iter_start_time
#endif
            !Check whether we have reached convergance
            if (.not. cont) exit nonlinear_iteration
     end do nonlinear_iteration

#ifdef OUTPUT_PARTIAL_ITERATIONS
      call end_debug_iteration(ncid_debug)
#endif
#if 0
      call write_xls_3d("uvel.txt",uvel)
      call write_xls_3d("vvel.txt",vvel)
      call write_xls("uvel_surf.txt",uvel(:,:,1))
      call write_xls("vvel_surf.txt",vvel(:,:,1))
#endif
      call sparse_destroy_workspace(matrix, matrix_options, matrix_workspace)
      call del_sparse_matrix(matrix) 

      call cpu_time(solve_end_time)
      deallocate(pointtype)
      write(*,*) "Pattyn higher-order solve took",solve_end_time - solve_start_time
      return
      END subroutine

      subroutine comp_type(kinematic_bc, geometry_mask, point_mask, outfilename)
            real(dp), dimension(:,:) :: kinematic_bc
            integer,  dimension(:,:) :: geometry_mask, point_mask
            character(*) outfilename

            integer, dimension(size(kinematic_bc,1), size(kinematic_bc,2)) :: types

            types = 1

            where(.not. GLIDE_HAS_ICE(geometry_mask))
                types = 0
            endwhere

            where(.not. IS_NAN(kinematic_bc))
                types = 2
            endwhere

            where(GLIDE_IS_CALVING(geometry_mask))
                types = 4
            endwhere
            
            where(GLIDE_IS_CALVING(geometry_mask) &
            .and. .not. IS_NAN(kinematic_bc))
                types = 6
            endwhere

            call write_xls_int(outfilename, types)
      end subroutine
  
!----------------------------------------------------------------------------------------
! BOP
!
! !IROUTINE: unstable_manifold_correction
!
! !INTERFACE:
    function umc_correct_vels(u_new, u_old, &
                                          v_new, v_old, vec_correction, &
                                          maxy, maxx, nzeta, toler, &
                                          tot_out, theta_out)
        ! !RETURN VALUE:
        logical :: umc_correct_vels !whether another iteration step is needed
 
! !PARAMETERS:
        real(dp), dimension(:,:,:), intent(in) :: u_new  !Computed u component from this iteration
        real(dp), dimension(:,:,:), intent(inout) :: u_old !Computed u component from last iteration
        real(dp), dimension(:,:,:), intent(in) :: v_new !Computed v component from this iteration
        real(dp), dimension(:,:,:), intent(inout) :: v_old !Computed v component from last iteration
        real(dp), dimension(:), intent(inout) :: vec_correction !Old correction vector for v
        integer, intent(in) :: maxy, maxx, nzeta !Grid size
        real(dp), intent(in) :: toler !Error tolerance for the iteration
        real(dp), intent(out):: tot_out !Optional output of error
        real(dp), intent(out):: theta_out !Optional output of angle

        real(dp), dimension(maxy*maxx*nzeta*2) :: vec_new, vec_old
        integer :: linearize_idx
        
        linearize_idx = 1
        call linearize_3d(vec_new, linearize_idx, u_new)
        call linearize_3d(vec_new, linearize_idx, v_new)
        
        linearize_idx = 1
        call linearize_3d(vec_old, linearize_idx, u_old)
        call linearize_3d(vec_old, linearize_idx, v_old)


        umc_correct_vels = unstable_manifold_correction(vec_new, &
                vec_old, vec_correction, maxy*maxx*nzeta*2, toler, tot_out, theta_out)

        linearize_idx = 1
        call delinearize_3d(vec_old, linearize_idx, u_old)
        call delinearize_3d(vec_old, linearize_idx, v_old)
        
    end function umc_correct_vels


    !Reduces a 3d higher order velocity estimate to a 2d velocity field.
    subroutine vel_2d_from_3d(uvel, vvel, u, v, zeta)
        real(dp), dimension(:,:,:), intent(in) :: uvel
        real(dp), dimension(:,:,:), intent(in) :: vvel
        real(dp), dimension(:,:),  intent(out) :: u
        real(dp), dimension(:,:),  intent(out) :: v
        real(dp), dimension(:), intent(in) :: zeta
    
        integer :: maxx, maxy, nzeta, i, j, k

        maxx = size(uvel,2)
        maxy = size(uvel,1)
        nzeta = size(uvel,3)

        u = 0
        v = 0
        do i=1,MAXY
            do j=1,MAXX
                do k=2,NZETA
                    !Average velocity for layer (k-1...k) * thickness of layer (k-1..k)
                    !Does this lead to a vertical velocity average???
                    u(i,j) = u(i,j) + (uvel(i,j,k)+uvel(i,j,k-1))*(zeta(k)-zeta(k-1))*0.5
                    v(i,j) = v(i,j) + (vvel(i,j,k)+vvel(i,j,k-1))*(zeta(k)-zeta(k-1))*0.5
                end do
            end do
        end do
    end subroutine vel_2d_from_3d


!-----------------------------------------------------
!   plastic bed computation
!-----------------------------------------------------
    subroutine plastic_bed(tau, tau0, ubas, vbas)
        real(dp), dimension(:,:), intent(out) :: tau
        real(dp), dimension(:,:), intent(in)  :: tau0
        real(dp), dimension(:,:), intent(in)  :: ubas
        real(dp), dimension(:,:), intent(in)  :: vbas
        
        integer :: maxy, maxx, i, j

        maxy = size(ubas,1)
        maxx = size(ubas,2)
        !TODO: Vectorize
        do i = 1,maxy
            do j = 1,maxx
                tau(i,j) = tau0(i,j) / (sqrt(ubas(i,j)**2 + vbas(i,j)**2 + plastic_bed_regularization ** 2))
            end do
        end do
    end subroutine plastic_bed


    !Compute x,y,z derivative fields of u and v, upwinding if necessary at the
    !ice shelf front
    subroutine velderivs(uvel, vvel, dzdx, dzdy, geometry_mask, & 
                         dx, dy, levels, UPSTREAM, &
                         direction_x, direction_y, &
                         dudx, dudy, dudz, dvdx, dvdy, dvdz)
        use glide_mask, only: upwind_from_mask

        real(dp), dimension(:,:,:), intent(in) :: uvel
        real(dp), dimension(:,:,:), intent(in) :: vvel
        real(dp), dimension(:,:), intent(in) :: dzdx
        real(dp), dimension(:,:), intent(in) :: dzdy
        integer, dimension(:,:), intent(in) :: geometry_mask
        real(dp), intent(in) :: dx
        real(dp), intent(in) :: dy
        real(dp), dimension(:), intent(in) :: levels
        logical, intent(in) :: UPSTREAM
        real(dp), dimension(:,:), intent(out) :: direction_x
        real(dp), dimension(:,:), intent(out) :: direction_y

        real(dp), dimension(:,:,:), intent(out) :: dudx
        real(dp), dimension(:,:,:), intent(out) :: dudy
        real(dp), dimension(:,:,:), intent(out) :: dudz
        real(dp), dimension(:,:,:), intent(out) :: dvdx        
        real(dp), dimension(:,:,:), intent(out) :: dvdy        
        real(dp), dimension(:,:,:), intent(out) :: dvdz
 
        real(dp), dimension(size(uvel,1),size(uvel,2),size(uvel,3)) :: uvel_test, vvel_test

        integer :: k

        !TEST of upwinding shelf boundary derivatives:
        !Everywhere where there's no ice, we set the velocity
        !to NaN.  If a location with no ice is ever used, this should
        !propegate.
        uvel_test = uvel
        vvel_test = vvel

#if 0
        do k=1,size(uvel,3)
            where(.not. GLIDE_HAS_ICE(geometry_mask))
                uvel_test(:,:,k) = NaN
                vvel_test(:,:,k) = NaN
            endwhere
        end do
#endif

        !TODO: Implement option for upstream differencing
        if (UPSTREAM) then
            !Upstream the number 
            call df_field_3d(uvel, dy, dx, levels, dudy, dudx, dudz, dzdy, dzdx)
            call df_field_3d(vvel, dy, dx, levels, dvdy, dvdx, dvdz, dzdy, dzdx)
        else
            direction_x = 0
            direction_y = 0
            call upwind_from_mask(geometry_mask, direction_y, direction_x)
            call write_xls("direction_x.txt",direction_x)
            call write_xls("direction_y.txt",direction_y)
            
            call df_field_3d(uvel_test, dy, dx, levels, dudy, dudx, dudz, direction_y, direction_x)
            call df_field_3d(vvel_test, dy, dx, levels, dvdy, dvdx, dvdz, direction_y, direction_x)
         end if

    end subroutine

!
!
!------------------------------------------------------
!   nonlinear viscosity term
!------------------------------------------------------
!
      subroutine muterm(mu, arrh, h, ax, ay, FLOWN, ZIP, &
                        dudx, dudy, dudz, dvdx, dvdy, dvdz)
!
        real(dp),                   intent(in)  :: FLOWN,ZIP
        real(dp), dimension(:,:,:), intent(out) :: mu
        real(dp), dimension(:,:,:), intent(in)  :: arrh
        real(dp), dimension(:,:),   intent(in)  :: h
        real(dp), dimension(:,:,:), intent(in)  :: ax
        real(dp), dimension(:,:,:), intent(in)  :: ay
        
        real(dp), dimension(:,:,:), intent(in) :: dudx, dudy, dudz
        real(dp), dimension(:,:,:), intent(in) :: dvdx, dvdy, dvdz
       !
        INTEGER :: i,j,k, MAXX, MAXY, NZETA
        real(dp) :: macht
        real(dp) :: exx,eyy,exy,exz,eyz,eeff
!
        MAXY = size(mu, 1)
        MAXX = size(mu, 2)
        NZETA = size(mu, 3)

        macht=(1.-FLOWN)/(2.*FLOWN)
     
        !Compute mu term from the acceleration fields
        do i=1,MAXY
            do j=1,MAXX
                do k=1,NZETA
          
                    if (h(i,j) > 0.) then
                        exx=dudx(i,j,k)+ax(i,j,k)*dudz(i,j,k)
                        eyy=dvdy(i,j,k)+ay(i,j,k)*dvdz(i,j,k)
                        exy=0.5*((dudy(i,j,k)+ay(i,j,k)*dudz(i,j,k))+(dvdx(i,j,k)+ax(i,j,k)*dvdz(i,j,k)))
                        exz=-0.5*dudz(i,j,k)/h(i,j)
                        eyz=-0.5*dvdz(i,j,k)/h(i,j)
                    else
                        exx=dudx(i,j,k)
                        eyy=dvdy(i,j,k)
                        exy=0.5*(dudy(i,j,k)+dvdx(i,j,k))
                        exz=0.
                        eyz=0.
                    endif
            
                    eeff=(exx*exx)+(eyy*eyy)+(exx*eyy)+(exy*exy)+(exz*exz)+(eyz*eyz)
                    mu(i,j,k)=0.5 * (arrh(i,j,k)**(-1./FLOWN)) * ((eeff+ZIP)**macht)
                end do
            end do
        end do
        return
    end subroutine

    function sparuv(mu,dzdx,dzdy,ax,ay,bx,by,cxy,h,uvel,vvel,dudx,dudy,dudz,dvdx,dvdy,dvdz,&
                    ustar,vstar,beta,dhbdx,dhbdy,&
                    IJKTOT,MAXY,MAXX,NZETA,TOLER,WHICH_SOURCE,GRIDX,GRIDY,zeta, point_mask, geometry_mask,&
                    matrix, matrix_workspace, matrix_options, kinematic_bc_u, kinematic_bc_v,latbc_normal, &
                    direction_x, direction_y, STAGGERED, INCLUDE_THIN)
        INTEGER IJKTOT,MAXY,MAXX,NZETA
        real(dp), dimension(:,:,:), intent(in) :: mu
        real(dp), dimension(:,:), intent(in) :: dzdx
        real(dp), dimension(:,:), intent(in) :: dzdy
        real(dp), dimension(:,:,:), intent(in) :: ax
        real(dp), dimension(:,:,:), intent(in) :: ay
        real(dp), dimension(:,:,:), intent(in) :: bx
        real(dp), dimension(:,:,:), intent(in) :: by
        real(dp), dimension(:,:,:), intent(in) :: cxy
        real(dp), dimension(:,:), intent(in) :: h
        real(dp), dimension(:,:,:), intent(in), target :: uvel
        real(dp), dimension(:,:,:), intent(in), target :: vvel
        
        real(dp), dimension(:,:,:), intent(in) :: dudx,dudy,dudz,dvdx,dvdy,dvdz
        
        real(dp), dimension(:,:,:), intent(out), target :: ustar
        real(dp), dimension(:,:,:), intent(out), target :: vstar
        real(dp), dimension(:,:), intent(in) :: dhbdx
        real(dp), dimension(:,:), intent(in) :: dhbdy
        real(dp), dimension(:,:), intent(in) :: beta
        real(dp), dimension(:), intent(in) :: zeta
        integer, dimension(:,:), intent(in) :: geometry_mask
        real(dp), dimension(:,:,:), intent(in), target :: kinematic_bc_u
        real(dp), dimension(:,:,:), intent(in), target :: kinematic_bc_v
        real(dp), dimension(:,:), intent(in)   :: latbc_normal !On the marine ice front, this is the angle of the normal to the ice front
        real(dp), intent(in) :: toler
        real(dp), intent(in) :: gridx
        real(dp), intent(in) :: gridy
        real(dp), dimension(:,:), intent(in) :: direction_x,direction_y

        !Sparse matrix variables.  These are passed in so that allocation can be
        !done once per velocity solve instead of once per iteration
        type(sparse_matrix_type), intent(inout) :: matrix
        type(sparse_solver_workspace), intent(inout) :: matrix_workspace
        type(sparse_solver_options), intent(inout) :: matrix_options

        
        logical, intent(in)::STAGGERED
        integer, intent(in)::WHICH_SOURCE
        logical, intent(in)::INCLUDE_THIN
        real(dp) :: rhs
        
        INTEGER i,j,k,m,sparuv,iter, ierr
        real(dp) :: d(IJKTOT),x(IJKTOT),coef(STENCIL_SIZE),err
      
        integer, dimension(:,:) :: point_mask
        
        integer :: stencil_center_idx
        integer :: si, sj, sk
        real(dp), dimension(:,:,:), pointer :: velpara, velperp, velpara_star, kinematic_bc_para
        integer :: whichcomponent
        character(1) :: componentstr

        pointtype = 0
!
!-------  velocity u
!
#ifdef OUTPUT_SPARSE_MATRIX
        write(ITER_UNIT,*) " Component Y X Z h mu point_type ", &
                   "im1jm1 im1km1 im1 im1kp1 im1jp1 jm1km1 jm1 jm1kp1 ", &
                   "km2 km1 center kp1 kp2 jp1km1 jp1 jp1kp1 ", &
                   "ip1jm1 ip1km1 ip1 ip1kp1 ip1jp1 im2 ip2 jm2 jp2 stencil_total rhs "
#endif

        sparuv = 0
        do whichcomponent = 1,2
            if (whichcomponent == 1) then !Set up to compute u component
                velpara => uvel
                velperp => vvel
                velpara_star => ustar
                componentstr = "u"
                kinematic_bc_para => kinematic_bc_u
            else !Compute v component
                velpara => vvel
                velperp => uvel
                velpara_star => vstar
                componentstr = "v"
                kinematic_bc_para => kinematic_bc_v
            end if
            !Initialize sparse matrix & vectors
            d=0
            x=0
            call sparse_clear(matrix)
#ifdef VERY_VERBOSE
            write(*,*)"Begin Matrix Assembly"
#endif
            do i=1,MAXY
                do j=1,MAXX
                    if (point_mask(i,j) /= 0) then
                        do k=1,NZETA
                            coef = 0
                            stencil_center_idx = csp_masked(I_J_K,i,j,k,point_mask,NZETA) 
                            if (.not. GLIDE_HAS_ICE( geometry_mask(i,j) ) .or. &
                                      GLIDE_IS_THIN( geometry_mask(i,j) ) &
                                         .and. .not. INCLUDE_THIN) then
                                !No ice - "pass through"
                                coef(I_J_K)=1.
                                !Normally, we use uvel(i,j,k) as our initial guess.
                                !However, in this case, we know the velocity
                                !should be 0, so we'll replace our uvel with that
                                velpara(i,j,k) = 0
                                rhs = 0
                                pointtype(i,j) = 1
                            else if (.not. IS_NAN(kinematic_bc_para(i,j,k))) then
                                !If a kinematic boundary condition was specified
                                !for this location, hold the location at the
                                !specified value.
                                coef(I_J_K)=1.
                                rhs = kinematic_bc_para(i,j,k)
                                pointtype(i,j) = 2
                            else if ((i.eq.1).or.(i.eq.MAXY).or.(j.eq.1).or.(j.eq.MAXX)) then

                                !Boundary condition at the edges of the domain.
                                !If we don't have a kinematic boundary already
                                !specified, we just "pass through" the initial
                                !guess.  If a kinematic b.c. is specified, we
                                !handle it in sparse_setup rather than here.
                                coef(I_J_K)=1.
                                rhs=velpara(i,j,k)
                                pointtype(i,j) = 3
                            else
                                pointtype(i,j) = 4
                                call sparse_setup(componentstr, i, j, k, mu, dzdx, dzdy, ax, ay, bx, by, cxy, &
                                     h, gridx, gridy, zeta, uvel, vvel, dudx, dudy, dudx, dvdx, dvdy, dvdz, &
                                     dhbdx, dhbdy, beta, geometry_mask, &
                                     latbc_normal, maxx, maxy, Nzeta, coef, rhs, direction_x, direction_y, STAGGERED, &
                                     WHICH_SOURCE)
                            endif
                            d(stencil_center_idx)=rhs
                            !Preliminary benchmarks indicate the we actually reach
                            !convergance faster if we use a 0 initial guess rather than
                            !the previous velocity.  More testing is needed.
                            !To use the current velocity as the guess, uncomment this
                            !line here and and in the Y velocity section of this
                            !function.
                            !x=csp(I_J_K,i,j,k,MAXX,NZETA))=uvel(i,j,k) !Use current velocity as guess of solution
                            if (abs(coef(I_J_K)) < SMALL) then
                                write(*,*) "WARNING: 0 on diagonal at position",i,j,k
                                write(*,*) "component:",componentstr
                                write(*,*) "Thickness at this point:", h(i,j)
                                write(*,*) "Mask at this point:", geometry_mask(i,j)
                                stop
                            end if

                            do m=1,STENCIL_SIZE
                                if (abs(coef(m)) > SMALL) then
                                    si = stencil_i(m,i)
                                    sj = stencil_j(m,j)
                                    sk = stencil_k(m,k)

                                    if (si > 0 .and. si <= maxy .and. &
                                        sj > 0 .and. sj <= maxx .and. &
                                        sk > 0 .and. sk <= nzeta) then
                                            if (point_mask(si,sj) == 0) then
                                                write(*,*) "ERROR: point is off mask."
                                                write(*,*) "component:",componentstr
                                                write(*,*) "location:",i,j,k
                                                write(*,*) "stencil:",si,sj,sk
                                                write(*,*) "position and coefficient in stencil:", m, coef(m)
                                                write(*,*) "h for stencil center and bad point:",h(i,j),h(si,sj)
                                                write(*,*) "Mask for stencil center and bad point:",&
                                                           geometry_mask(i,j), geometry_mask(si,sj)
                                                write(*,*) "Coefficient:", coef(m)
                                            end if
                                            call sparse_insert_val(matrix, &
                                                    stencil_center_idx, &
                                                    csp_stenciled(si,sj,sk,point_mask,NZETA), &
                                                    coef(m)) 
                                    end if
                                end if
                            end do !End loop through stencil
                        end do !END k loop
                    end if !End mask check
                end do !End j loop
            end do !End i loop
#ifdef VERY_VERBOSE
            write(*,*)"End Matrix Assembly"
            write(*,*)"Begin Matrix Solve"
#endif
            call write_xls_int("point_type.txt", pointtype)
            call sparse_solver_preprocess(matrix, matrix_options, matrix_workspace)        
            ierr = sparse_solve(matrix, d, x, matrix_options, matrix_workspace,  err, iter, verbose=sparverbose)
            sparuv = sparuv + iter
        
            call handle_sparse_error(matrix, matrix_options, ierr, __FILE__, __LINE__)      
            call sparse_solver_postprocess(matrix, matrix_options, matrix_workspace)
#ifdef VERY_VERBOSE
            write(*,*)"End Matrix Solve"
#endif
            !Delinearize the solution
            do i=1,MAXY
                do j=1,MAXX
                    do k=1,NZETA
                        if (point_mask(i,j) /= 0) then
                            !write(*,*)csp_masked(11,i,j,k,point_mask,nzeta)
                            velpara_star(i,j,k)=x(csp_masked(11,i,j,k,point_mask,nzeta))
                        else
                            velpara_star(i,j,k)=0
                        end if
                    end do
                end do
            end do
        end do !END whichcomponent loop
    end function sparuv
!
!
!
!------------------------------------------------------
!   lookup coordinates in sparse matrix
!
!   Reference:
!     1   ->   i - 1, j - 1, k
!     2   ->   i - 1, j, k - 1
!     3   ->   i - 1, j, k
!     4   ->   i - 1, j, k + 1
!     5   ->   i - 1, j + 1, k
!     6   ->   i, j - 1, k - 1
!     7   ->   i, j - 1, k
!     8   ->   i, j - 1, k + 1
!     9   ->   i, j, k - 2
!    10   ->   i, j, k - 1
!    11   ->   i, j, k
!    12   ->   i, j, k + 1
!    13   ->   i, j, k + 2
!    14   ->   i, j + 1, k - 1
!    15   ->   i, j + 1, k
!    16   ->   i, j + 1, k + 1
!    17   ->   i + 1, j - 1, k
!    18   ->   i + 1, j, k - 1
!    19   ->   i + 1, j, k
!    20   ->   i + 1, j, k + 1
!    21   ->   i + 1, j + 1, k
!------------------------------------------------------
!
      !Given a point and a stencil location,
      !the following three functions return the equivalent location
      function stencil_i(pos, i)
        integer, intent(in) :: pos, i
        integer :: stencil_i
        if (pos <= 5) then
            stencil_i = i - 1
        else if (pos <= 16 .or. pos == 24 .or. pos == 25) then
            stencil_i = i
        else if (pos == 22) then
            stencil_i = i - 2
        else if (pos == 23) then
            stencil_i = i + 2
        else
            stencil_i = i + 1
        end if
      end function

      function stencil_j(pos, j)
        integer, intent(in) :: pos, j
        integer :: stencil_j
        select case(pos)
            case(1, 6, 7, 8, 17)
                stencil_j = j - 1
            case(5, 14, 15, 16, 21)
                stencil_j = j + 1
            case(24)
                stencil_j = j - 2
            case(25)
                stencil_j = j + 2
            case default
                stencil_j = j
        end select
      end function

      function stencil_k(pos, k)
        integer, intent(in) :: pos, k
        integer :: stencil_k
        select case(pos)
            case(2, 6, 10, 14, 18)
                stencil_k = k - 1
            case(4, 8, 12, 16, 20)
                stencil_k = k + 1
            case(9)
                stencil_k = k - 2
            case(13)
                stencil_k = k + 2
            case default
                stencil_k = k
        end select
      end function
 !
      function csp_masked(pos, i, j, k, point_mask, nzeta)
         integer, intent(in) :: pos
         integer, intent(in) :: i !Y coordinate
         integer, intent(in) :: j !X coordinate
         integer, intent(in) :: k !Z coordinate
         integer, dimension(:,:), intent(in) :: point_mask
         integer, intent(in) :: nzeta

         integer :: csp_masked

         csp_masked = (point_mask(stencil_i(pos,i), stencil_j(pos,j)) - 1) * nzeta &
                      + stencil_k(pos, k)
      end function

      !Returns the coordinate without performing stencil lookups (this assumes
      !those have already been done)
      function csp_stenciled(si, sj, sk, point_mask, nzeta)
            integer, intent(in) :: si, sj, sk
            integer, dimension(:,:), intent(in) :: point_mask
            integer, intent(in) :: nzeta

            integer :: csp_stenciled

            csp_stenciled = (point_mask(si, sj) - 1)*nzeta + sk
      end function

!
!------------------------------------------------------
!   central difference calculation for uvel
!------------------------------------------------------
!

    !This differencing scheme function uses the following nomenclature.
    !This reduces the code duplication that would have otherwise happened!!!
    !u - Refers to velocity in the x direction
    !v - Refers to velocity in the y direction
    !para - Refers to a vector that is parallel to the component of the
    !       velocity being computed
    !perp - Refers to a vector that is perpendicular to the component of
    !       the velocity being computed
    !Under this nomenclature, d_para_d_para means either "dudx" or "dvdy".
    !d_perp_d_perp means "dvdy" or "dudx" (in each case, for computing u and v)
    !Note: x and y might be referred to directly when they don't change
    !      depending on the differencing scheme being used
    !
    !The "component" argument should be a single character indicating
    !whether the "u" or "v" component is being requested.
    !DON'T PANIC, this works!  It's been tested!  Really!  The bug is somewhere
    !else!
    subroutine sparse_setup(component, i,j,k,mu,dzdx,dzdy,ax,ay,bx,by,cxy,&
        h, dx, dy, dz, uvel, vvel, dudx_field, dudy_field, dudz_field, &
        dvdx_field, dvdy_field, dvdz_field, &
        dhbdx,dhbdy,beta,geometry_mask,latbc_normal,MAXX,MAXY,Ndz,&
        coef, rhs, direction_x, direction_y, STAGGERED, WHICH_SOURCE)
!
        integer, intent(in) :: i,j,k, MAXY, MAXX, Ndz

        !Output array
        real(dp), dimension(STENCIL_SIZE), intent(out) :: coef

        !Output RHS value
        real(dp), intent(out) :: rhs

        !Viscosity
        real(dp), dimension(:,:,:), intent(in) :: mu

        real(dp), dimension(:,:,:), intent(in) :: dudx_field, dudy_field, dudz_field 
        real(dp), dimension(:,:,:), intent(in) :: dvdx_field, dvdy_field, dvdz_field

        !Surface Gradients
        real(dp), dimension(:,:), target, intent(in) :: dzdx, dzdy

        !Rescaling Factors
        real(dp), dimension(:,:,:), target, intent(in) :: ax, ay, bx, by, cxy
        
        !Ice Thickness
        real(dp), dimension(:,:), intent(in) :: h

        !Grid Spacing (Z is an irregular grid)
        real(dp), intent(in) :: dx, dy
        real(dp), dimension(:), intent(in) :: dz

        !Current velocities
        real(dp), dimension(:,:,:), target, intent(in) :: uvel, vvel

        !Bedrock Gradients
        real(dp), dimension(:,:), target, intent(in) :: dhbdx, dhbdy

        !Basal Traction
        real(dp), dimension(:,:), intent(in) :: beta

        integer, dimension(:,:), intent(in) :: geometry_mask

        !Contains the angle of the normal to the marine margine, NaN
        !everywhere not on the margin
        real(dp), dimension(:,:), intent(in) :: latbc_normal

        logical, intent(in) :: STAGGERED

        integer, intent(in) :: WHICH_SOURCE

        real(dp), dimension(:,:),intent(in) :: direction_x, direction_y

        character(1), intent(in) :: component
!
        !Temporary calculations done in the ABSOLUTE (using u/x, v/y) coordinate
        !system
        real(dp), target :: dmudx, dmudy, dmudz, dmudx2, dmudy2
        real(dp), target :: dudx, dudy, dudz, dudx2, dudy2, dudz2, dudxz, dudyz, dudxy
        real(dp), target :: dvdx, dvdy, dvdz, dvdx2, dvdy2, dvdz2, dvdxz, dvdyz, dvdxy

        !Temporary calculations done in the RELATIVE (using para, perp)
        !coordinate system
        real(dp), pointer :: dpara_dpara, dpara_dperp, dpara_dx, dpara_dy, dpara_dz
        real(dp), pointer :: dpara_dpara2, dpara_dparaz, dpara_dz2, dpara_dperp2, dpara_dperpz
        real(dp), pointer :: dpara_dyz, dpara_dxz, dpara_dx2, dpara_dy2
        real(dp), pointer :: dperp_dpara, dperp_dperp, dperp_dx, dperp_dy, dperp_dz
        real(dp), pointer :: dperp_dxy, dperp_dxz, dperp_dyz, dperp_dz2
        real(dp), pointer :: dmu_dpara2, dmu_dperp2

        !Fields done in the RELATIVE coordinate system
        real(dp), dimension(:,:,:), pointer :: a_para, a_perp, b_para, b_perp
        real(dp), dimension(:,:,:), pointer :: vel_perp, vel_para
        real(dp), dimension(:,:), pointer :: dz_dpara, dz_dperp, dhb_dpara, dhb_dperp

        !Cached difference calculations for the nonstaggered Z grid
        real(dp) dz_down1, dz_down2, dz_down3, dz_up1, dz_up2, dz_up3
        real(dp) dz_cen1, dz_cen2, dz_cen3, dz_sec1, dz_sec2, dz_sec3
               !call write_xls_3d("uvel.txt",uvel)
        !call write_xls_3d("vvel.txt",vvel)
        character(20) :: point_type

        !Set up a system of pointers to encapsulate the nomenclature described
        !above
        !TODO: Go through the code and figure out which of these bindings
        !      can just become the variable that gets assigned to
        if (component == "u") then
            dpara_dpara => dudx
            dpara_dperp => dudy
            dpara_dx    => dudx
            dpara_dy    => dudy
            dpara_dz    => dudz
            
            dpara_dpara2 => dudx2
            dpara_dparaz => dudxz
            dpara_dz2    => dudz2
            dpara_dx2    => dudx2
            dpara_dy2    => dudy2
            dpara_dperp2 => dudy2
            dpara_dperpz => dudyz
            dpara_dyz    => dudyz
            dpara_dxz    => dudxz

            dperp_dpara => dvdx
            dperp_dperp => dvdy
            dperp_dx    => dvdx
            dperp_dy    => dvdy
            dperp_dz    => dvdz

            dperp_dxy => dvdxy
            dperp_dxz => dvdxz
            dperp_dyz => dvdyz
            dperp_dz2 => dvdz2

            dmu_dpara2 => dmudx2
            dmu_dperp2 => dmudy2

            a_para => ax
            a_perp => ay
            b_para => bx
            b_perp => by

            dz_dpara  => dzdx
            dz_dperp  => dzdy
            dhb_dpara => dhbdx
            dhb_dperp => dhbdy

            vel_perp => vvel
            vel_para => uvel
        else if (component == "v") then
            dpara_dpara => dvdy
            dpara_dperp => dvdx
            dpara_dx    => dvdx
            dpara_dy    => dvdy
            dpara_dz    => dvdz
            
            dpara_dpara2 => dvdy2
            dpara_dparaz => dvdyz
            dpara_dx2    => dvdx2
            dpara_dy2    => dvdy2
            dpara_dz2    => dvdz2
            dpara_dperp2 => dvdx2
            dpara_dperpz => dvdxz
            dpara_dyz    => dvdyz
            dpara_dxz    => dvdxz

            dperp_dpara => dudy
            dperp_dperp => dudx
            dperp_dx    => dudx
            dperp_dy    => dudy
            dperp_dz    => dudz

            dperp_dxy => dudxy
            dperp_dxz => dudxz
            dperp_dyz => dudyz
            dperp_dz2 => dudz2

            dmu_dpara2 => dmudy2
            dmu_dperp2 => dmudx2

            a_para => ay
            a_perp => ax
            b_para => by
            b_perp => bx

            dz_dpara  => dzdy
            dz_dperp  => dzdx
            dhb_dpara => dhbdy
            dhb_dperp => dhbdx

            vel_perp => uvel
            vel_para => vvel
        else
            write(*,*)"FATAL ERROR: sparse_setup called with invalid component"
        end if
        
        if (GLIDE_IS_CALVING(geometry_mask(i,j)) .and. &
            WHICH_SOURCE /= HO_SOURCE_DISABLED) then !Marine margin dynamic (Neumann) boundary condition
            call sparse_marine_margin(component,i,j,k,h,latbc_normal, vel_perp, mu, dx, dy, ax, ay, dz,coef, rhs, &
                                      dudx_field,dudy_field,dudz_field,dvdx_field,dvdy_field,dvdz_field,&
                                      direction_x, direction_y, geometry_mask, STAGGERED, WHICH_SOURCE)
            point_type = "lateral"
        else if (k.eq.1) then !Upper boundary condition (stress-free surface)
            point_type = "surface"
            !Finite difference coefficients for an irregular Z grid, downwinded
            dz_down1=(2.*dz(k)-dz(k+1)-dz(k+2))/(dz(k+1)-dz(k))/(dz(k+2)-dz(k))
            dz_down2=(dz(k+2)-dz(k))/(dz(k+2)-dz(k+1))/(dz(k+1)-dz(k))
            dz_down3=(dz(k)-dz(k+1))/(dz(k+2)-dz(k+1))/(dz(k+2)-dz(k))
 
            dpara_dpara = 4*dz_dpara(i,j)
            dpara_dz = 4*a_para(i,j,k)*dz_dpara(i,j) + a_perp(i,j,k)*dz_dperp(i,j) + 1./h(i,j)
            dpara_dperp = dz_dperp(i,j)

            dperp_dpara = dz_dperp(i,j)
            dperp_dperp = 2*dz_dpara(i,j)
            dperp_dz = 2.*a_perp(i,j,k)*dz_dpara(i,j)+a_para(i,j,k)*dz_dperp(i,j)

            coef(IM1_J_K) = -.5*dpara_dy/dy
            coef(I_JM1_K) = -.5*dpara_dx/dx
            coef(I_J_K) = dpara_dz*dz_down1
            coef(I_J_KP1) = dpara_dz*dz_down2
            coef(I_J_KP2) = dpara_dz*dz_down3
            coef(I_JP1_K) = dpara_dx*.5/dx
            coef(IP1_J_K) = dpara_dy*.5/dy
            !Note the transposition below (dfdx => dfdy, vice versa)
            rhs = -dperp_dx * dfdy_3d(vel_perp,i,j,k,dx) &
                -  dperp_dy * dfdx_3d(vel_perp,i,j,k,dy) &
                -  dperp_dz * dfdz_3d_downwind_irregular(vel_perp,i,j,k,dz)
                
            !"Blend" with lateral boundary condition if we are on the
            !boundary.  This functionally is just summing the equations.
            !Note that sparse_marine_margin is written in such a way as
            !to always add to coef, not overwrite it.
            !if (.not. IS_NAN(latbc_normal(i,j))) then !Marine margin dynamic (Neumann) boundary condition
            !    rhs2 = 0
            !    call sparse_marine_margin(component,i,j,k,h,latbc_normal, vel_perp, mu, dx, dy, ax, ay, dz,coef, rhs2)
            !    rhs = rhs + rhs2
            !end if
        else if (k.eq.Ndz) then !Lower boundary condition (Basal sliding)
            !Finite difference coefficients for an irregular Z grid, upwinded
            dz_up1=(dz(k)-dz(k-1))/(dz(k-1)-dz(k-2))/(dz(k)-dz(k-2))
            dz_up2=(dz(k-2)-dz(k))/(dz(k)-dz(k-1))/(dz(k-1)-dz(k-2))
            dz_up3=(2.*dz(k)-dz(k-1)-dz(k-2))/(dz(k)-dz(k-1))/(dz(k)-dz(k-2))

            if (IS_NAN(beta(i,j))) then !No sliding at the base
                coef(I_J_K)=1.
                rhs=0.
                point_type = "base-no-slip"
            else !Compute sliding 
                point_type = "base-stress-free"
                dpara_dpara = 4.*dhb_dpara(i,j)
                dpara_dz = 4.*a_para(i,j,k)*dhb_dpara(i,j) + a_perp(i,j,k)*dhb_dperp(i,j) + 1./h(i,j)
                dpara_dperp = dhb_dperp(i,j)

                dperp_dpara = dhb_dperp(i,j)
                dperp_dperp = 2.*dhb_dpara(i,j)
                dperp_dz = 2.*a_perp(i,j,k)*dhb_dpara(i,j) + a_para(i,j,k)*dhb_dperp(i,j)

                !If we have a non-zero friction, then we need to scale the
                !strain rates by the viscosity.  If there is no friction, then
                !we have a stress-free base and these coefficients disappear.
                if (beta(i,j) >= SMALL) then
                   point_type = "base-beta-squared"
                   dpara_dpara = dpara_dpara*mu(i,j,k)
                   dpara_dperp = dpara_dperp*mu(i,j,k)
                   dpara_dz = dpara_dz*mu(i,j,k)
                   dperp_dpara = dperp_dpara*mu(i,j,k)
                   dperp_dperp = dperp_dperp*mu(i,j,k)
                   dperp_dz = dperp_dz*mu(i,j,k)
                end if

                coef(IM1_J_K) = -.5*dpara_dy/dy
                coef(I_JM1_K) = -.5*dpara_dx/dx
                coef(I_J_KM2) = dpara_dz*dz_up1
                coef(I_J_KM1) = dpara_dz*dz_up2
                coef(I_J_K)   = dpara_dz*dz_up3
                coef(I_JP1_K) = .5*dpara_dx/dx
                coef(IP1_J_K) = .5*dpara_dy/dy
                !Transposition of derivatives below
                rhs = -dperp_dx * dfdy_3d(vel_perp, i, j, k, dx) &
                    -  dperp_dy * dfdx_3d(vel_perp, i, j, k, dy) &
                    -  dperp_dz * dfdz_3d_upwind_irregular(vel_perp, i, j, k, dz)

                !Adjust center of stencil for the basal friction parameter (this
                !occurs in all cases b/c the ice shelf case necessarily has beta = 0
                !CORRECTION: Adding in normalization term that Pattyn has
                !neglected.
                coef(11) = coef(11) + beta(i,j)! * sqrt(1 + dhbdx(i,j)**2 + dhbdy(i,j)**2)
                
                !"Blend" with lateral boundary condition if we are on the
                !boundary.  This functionally is just summing the equations.
                !Note that sparse_marine_margin is written in such a way as
                !to always add to coef, not overwrite it.
                !if (.not. IS_NAN(latbc_normal(i,j))) then !Marine margin dynamic (Neumann) boundary condition
                !    rhs2 = 0
                !    call sparse_marine_margin(component,i,j,k,h,latbc_normal, vel_perp, mu, dx, dy, ax, ay, dz,coef, rhs2)
                !    rhs = rhs + rhs2
                !end if
            endif
        else !Interior of the ice (e.g. not a boundary condition)
            !Finite difference coefficients for an irregular Z grid
            dz_cen1 = (dz(k)-dz(k+1))/(dz(k)-dz(k-1))/(dz(k+1)-dz(k-1))
            dz_cen2 = (dz(k+1)-2.*dz(k)+dz(k-1))/(dz(k)-dz(k-1))/(dz(k+1)-dz(k))
            dz_cen3 = (dz(k)-dz(k-1))/(dz(k+1)-dz(k))/(dz(k+1)-dz(k-1))
            dz_sec1 = 2./(dz(k)-dz(k-1))/(dz(k+1)-dz(k-1))
            dz_sec2 = 2./(dz(k)-dz(k+1))/(dz(k)-dz(k-1))
            dz_sec3 = 2./(dz(k+1)-dz(k))/(dz(k+1)-dz(k-1))

            !Derivative transposition below
            dmudx=dfdy_3d(mu,i,j,k,dx)
            dmudy=dfdx_3d(mu,i,j,k,dy)
            dmudz=dfdz_3d_irregular(mu,i,j,k,dz)
            dmudx2=dmudx+ax(i,j,k)*dmudz
            dmudy2=dmudy+ay(i,j,k)*dmudz

            dpara_dpara = 4.*dmu_dpara2
            dpara_dz = 4.*a_para(i,j,k)*dmu_dpara2 + 4.*b_para(i,j,k)*mu(i,j,k) &
                     +    a_perp(i,j,k)*dmu_dperp2 +    b_perp(i,j,k)*mu(i,j,k) &
                     +    dmudz/(h(i,j)*h(i,j))
            dpara_dpara2 = 4.*mu(i,j,k)
            dpara_dz2 = mu(i,j,k)*(4.*a_para(i,j,k)**2 + a_perp(i,j,k)**2 + 1/h(i,j)**2)
            dpara_dparaz = 8.*a_para(i,j,k)*mu(i,j,k)
            dpara_dperp=dmu_dperp2
            dpara_dperp2=mu(i,j,k)
            dpara_dperpz=2.*mu(i,j,k)*a_perp(i,j,k)

            dperp_dpara = dmu_dperp2
            dperp_dperp = 2 * dmu_dpara2
            dperp_dz  = 2 * a_perp(i,j,k) * dmu_dpara2  + a_para(i,j,k) * dmu_dperp2 + 3 * mu(i,j,k) * cxy(i,j,k) 
            dperp_dxy = 3 * mu(i,j,k)
            dperp_dxz = 3 * mu(i,j,k) * ay(i,j,k)
            dperp_dyz = 3 * mu(i,j,k) * ax(i,j,k)
            dperp_dz2 = 3 * mu(i,j,k) * ax(i,j,k)*ay(i,j,k)

            coef(IM1_J_KM1)  = dpara_dyz*dz_cen1*(-.5/dy)
            coef(IM1_J_K)  = (dpara_dyz*dz_cen2 + dpara_dy)*(-.5/dy) + dpara_dy2/(dy**2)
            coef(IM1_J_KP1)  = dpara_dyz*dz_cen3*(-.5/dy)
          
            coef(I_JM1_KM1)  = dpara_dxz*dz_cen1*(-.5/dx)
            coef(I_JM1_K)  = (dpara_dx + dpara_dxz*dz_cen2)*(-.5/dx) + dpara_dx2 / dx**2
            coef(I_JM1_KP1)  = dpara_dxz*dz_cen3*(-.5/dx)
          
            coef(I_J_KM1) = dpara_dz*dz_cen1 + dpara_dz2*dz_sec1
            coef(I_J_K) = dpara_dz*dz_cen2 + dpara_dz2*dz_sec2 + dpara_dx2*(-2/dx**2) + dpara_dy2*(-2/dy**2)
            coef(I_J_KP1) = dpara_dz*dz_cen3 + dpara_dz2*dz_sec3
          
            coef(I_JP1_KM1) = dpara_dxz * dz_cen1 * (.5/dx)
            coef(I_JP1_K) = (dpara_dx + dpara_dxz*dz_cen2) * (.5/dx) + dpara_dx2 / dx**2
            coef(I_JP1_KP1) = dpara_dxz * dz_cen3 * (.5/dx)
          
            coef(IP1_J_KM1) = dpara_dyz * dz_cen1 * (.5/dy)
            coef(IP1_J_K) = (dpara_dyz*dz_cen2 + dpara_dy) * (.5/dy) + dpara_dy2 / dy**2
            coef(IP1_J_KP1) = dpara_dyz * dz_cen3 * (.5/dy)
          
            !Transposition of derivatives below
            rhs =RHOI * GRAV * dz_dpara(i,j) &
                - dperp_dx  * dfdy_3d(vel_perp,i,j,k,dx)&
                - dperp_dy  * dfdx_3d(vel_perp,i,j,k,dy)&
                - dperp_dz  * dfdz_3d_irregular(vel_perp,i,j,k,dz)&
                - dperp_dxy * d2fdxy_3d(vel_perp,i,j,k,dx,dy)&
                - dperp_dxz * d2fdyz_3d(vel_perp,i,j,k,dx,dz)&
                - dperp_dyz * d2fdxz_3d(vel_perp,i,j,k,dy,dz)&
                - dperp_dz2 * d2fdz2_3d_irregular(vel_perp,i,j,k,dz)
        
            point_type = "interior"
        end if

#ifdef OUTPUT_SPARSE_MATRIX
        write(ITER_UNIT,*)component,i,j,k,h(i,j),mu(i,j,k), &
                  point_type, coef,sum(coef),rhs
#endif

    end subroutine sparse_setup

    !Computes finite differences for the marine margin
    subroutine sparse_marine_margin(component,i,j,k,h,normals, vel_perp, mu, dx, dy, ax, ay, zeta,coef, rhs, &
                                    dudx,dudy,dudz,dvdx,dvdy,dvdz,direction_x,direction_y, geometry_mask, STAGGERED, &
                                    WHICH_SOURCE)
        character(*), intent(in) :: component !*FD Either "u" or "v"
        integer, intent(in) :: i,j,k !*FD Point that the boundary condition is computed for
        real(dp), dimension(:,:), intent(in) :: h !*FD Ice thickness field
        real(dp), dimension(:,:), intent(in) :: normals !*FD Pre-computed angles that are normal to the marine margin
        real(dp), dimension(:,:,:), intent(in) :: vel_perp !*FD Velocity component that is perpendicular to the one being computed
        real(dp), dimension(:,:,:), intent(in) :: mu !*FD effective viscosity

        real(dp), intent(in) :: dx, dy !*FD grid spacing in x and y directions
        
        real(dp), intent(in), dimension(:,:,:), target :: ax, ay !*FD coefficients introduced during vertical rescaling
        
        real(dp), dimension(:), intent(in) :: zeta !*FD Irregular grid levels in vertical coordinate
        
        real(dp), dimension(STENCIL_SIZE), intent(inout) :: coef !*FD Output array of stencil coefficients
        real(dp),                intent(out)   :: rhs !*FD Element of the right-hand side vector corresponding to this point
        
        real(dp), dimension(:,:,:), intent(in) :: dudx,dudy,dudz,dvdx,dvdy,dvdz
        
        integer, dimension(:,:), intent(in) :: geometry_mask

        real(dp),dimension(:,:) :: direction_x, direction_y

        logical, intent(in) :: STAGGERED
        integer, intent(in) :: WHICH_SOURCE

        !Whether or not we need to upwind or downwind lateral derivatives
        !Upwind is for when there's no ice on positive side,
        !downwind is for when there's no ice on negative side
        logical :: upwind_x, upwind_y, downwind_x, downwind_y

        !x and y components of the normal vector with unit length
        real(dp), target :: n_x, n_y

        !Hydrostatic pressure at this level
        real(dp) :: source
        
        !Derivative of the perpendicular component (v if computing u and vice
        !versa) with respect to the parallel direction (x if computing u, y if
        !computing v)
        real(dp), target :: dperp_dx
        real(dp), target :: dperp_dy
        real(dp), target :: dperp_dz

        !The derivatives of the velocity component being computed with respect
        !to each dimension.  These end up as coefficients in the stencil
        !definition
        real(dp), target :: dx_coeff, dy_coeff, dz_coeff

        !Cached difference calculations for the nonstaggered Z grid
        real(dp) dz_down1, dz_down2, dz_down3, dz_up1, dz_up2, dz_up3
        real(dp) dz_cen1, dz_cen2, dz_cen3 
        real(dp), dimension(:,:,:), pointer :: a_para, a_perp

        real(dp), pointer :: n_para, n_perp, dperp_dpara, dperp_dperp
        real(dp), pointer :: para_coeff, perp_coeff

        real(dp) :: ntot = 1d0

        !Get the normal unit vector
        !The angles are defined in radians, with 0 = 12 o' clock, pi/2 = 3 o' clock,
        !etc.
        !We will need to convert this into an angle system with 0 = 3 o' clock,
        !pi/2 = 12 o' clock, etc.
        !At first I thought:
        !theta = normals(i,j) + 2*(PI/4 - normals(i,j))
        !But really this reduces to swapping sine and cosine (easy proof to work
        !through).
        !n_x =  (sin(normals(i,j)))
        !n_y = -(cos(normals(i,j)))

        n_x = sin(normals(i,j))
        n_y = -cos(normals(i,j))

       !Clamp n_x and n_y so that small values register as zero
        !(The need to do this arises from the fact that, sin(k*pi)
        !cannot be exactly 0)
        if (abs(n_x) < 1d-10) then
            n_x = 0
        end if

        if (abs(n_y) < 1d-10) then
            n_y = 0
        end if
   
        !Facilitates proper splitting of the source term so the it adds to 1
        !In the case of a 45 deg. boundary, this will result in changing
        !1/sqrt(2) to 1/2.  It should do nothing in the case of a boundary
        !aligned with the axis.
        !ntot = abs(n_x) + abs(n_y)


        !Divide the source term by viscosity, less error-prone to do it here
        !rather than multiply everything else by it (though it may give
        !iterative solvers more fits maybe?)
        source = source_term(i, j, k, H, zeta, WHICH_SOURCE, STAGGERED)/mu(i,j,k)

        !Determine whether to use upwinded or downwinded derivatives for the
        !horizontal directions.
        upwind_x = .false.
        upwind_y = .false.
        downwind_x = .false.
        downwind_y = .false.

        pointtype(i,j) = 5

        if (component=="u") then
            n_para => n_x
            n_perp => n_y

            dperp_dpara => dperp_dx
            dperp_dperp => dperp_dy
            
            a_para => ax
            a_perp => ay

            para_coeff => dx_coeff
            perp_coeff => dy_coeff
        
            dperp_dx = dvdx(i,j,k)
            dperp_dy = dvdy(i,j,k)
            dperp_dz = dvdz(i,j,k)
        else if (component=="v") then
            n_para => n_y
            n_perp => n_x
            
            dperp_dpara => dperp_dy
            dperp_dperp => dperp_dx
            
            a_para => ay
            a_perp => ax

            para_coeff => dy_coeff
            perp_coeff => dx_coeff

            dperp_dx = dudx(i,j,k)
            dperp_dy = dudy(i,j,k)
            dperp_dz = dudz(i,j,k)
        end if

        if (direction_y(i,j) > 0) then
           downwind_y = .true.
        else if (direction_y(i,j) < 0) then
            upwind_y = .true.
        end if

        if (direction_x(i,j) > 0) then
            downwind_x = .true.
        else if (direction_x(i,j) < 0) then
            upwind_x = .true.
        end if
       
        para_coeff = 4*n_para
        perp_coeff =   n_perp
        
        dz_coeff   =   (4*a_para(i,j,k)*n_para + a_perp(i,j,k)*n_perp)

        rhs = source/ntot * n_para &
                - 2*n_para*dperp_dperp & 
                -   n_perp*dperp_dpara &
                -   (2*a_para(i,j,k)*n_perp + a_perp(i,j,k)*n_para)*dperp_dz

        !Enter the z component into the finite difference scheme.
        !If we are on the top of bottom of the ice shelf, we will need
        !to upwind or downwind respectively
        !These terms arise from the vertical rescaling of the system;
        !any vertical explicitness of the system its self has been
        !thrown out based on the assumption of plug flow.
        if (k==1) then !Top of the ice shelf
            !Finite difference coefficients for an irregular Z grid, downwinded
            dz_down1=(2.*zeta(k)-zeta(k+1)-zeta(k+2))/(zeta(k+1)-zeta(k))/(zeta(k+2)-zeta(k))
            dz_down2=(zeta(k+2)-zeta(k))/(zeta(k+2)-zeta(k+1))/(zeta(k+1)-zeta(k))
            dz_down3=(zeta(k)-zeta(k+1))/(zeta(k+2)-zeta(k+1))/(zeta(k+2)-zeta(k))
        
            coef(I_J_KP2) = coef(I_J_KP2) + dz_coeff*dz_down3
            coef(I_J_KP1) = coef(I_J_KP1) + dz_coeff*dz_down2
            coef(I_J_K)   = coef(I_J_K)   + dz_coeff*dz_down1
        else if (k==size(zeta)) then !Bottom of the ice shelf
            !Finite difference coefficients for an irregular Z grid, upwinded
            dz_up1=(zeta(k)-zeta(k-1))/(zeta(k-1)-zeta(k-2))/(zeta(k)-zeta(k-2))
            dz_up2=(zeta(k-2)-zeta(k))/(zeta(k)-zeta(k-1))/(zeta(k-1)-zeta(k-2))
            dz_up3=(2.*zeta(k)-zeta(k-1)-zeta(k-2))/(zeta(k)-zeta(k-1))/(zeta(k)-zeta(k-2))

            coef(I_J_KM2) = coef(I_J_KM2) + dz_coeff*dz_up1
            coef(I_J_KM1) = coef(I_J_KM1) + dz_coeff*dz_up2
            coef(I_J_K)   = coef(I_J_K)   + dz_coeff*dz_up3
        else !Middle of the ice shelf, use centered difference
            dz_cen1 = (zeta(k)-zeta(k+1))/(zeta(k)-zeta(k-1))/(zeta(k+1)-zeta(k-1))
            dz_cen2 = (zeta(k+1)-2.*zeta(k)+zeta(k-1))/(zeta(k)-zeta(k-1))/(zeta(k+1)-zeta(k))
            dz_cen3 = (zeta(k)-zeta(k-1))/(zeta(k+1)-zeta(k))/(zeta(k+1)-zeta(k-1))

            coef(I_J_KM1) = coef(I_J_KM1) + dz_coeff*dz_cen1
            coef(I_J_K)   = coef(I_J_K)   + dz_coeff*dz_cen2
            coef(I_J_KP1) = coef(I_J_KP1) + dz_coeff*dz_cen3
        end if
        
        !Apply finite difference approximations of the x and y derivatives.  The
        !coefficients have already been computed above, so we just need to make
        !sure to use the correct difference form (upwind or downwind, etc.)
        if (downwind_y) then
            coef(IP2_J_K) = coef(IP2_J_K) - 0.5 * dy_coeff/dy
            coef(IP1_J_K) = coef(IP1_J_K) + 2.0 * dy_coeff/dy
            coef(I_J_K)   = coef(I_J_K)   - 1.5 * dy_coeff/dy
        else if (upwind_y) then
            coef(I_J_K)   = coef(I_J_K)   + 1.5 * dy_coeff/dy
            coef(IM1_J_K) = coef(IM1_J_K) - 2.0 * dy_coeff/dy
            coef(IM2_J_K) = coef(IM2_J_K) + 0.5 * dy_coeff/dy
        else
            coef(IP1_J_K) = coef(IP1_J_K) + .5*dy_coeff/dy
            coef(IM1_J_K) = coef(IM1_J_K) - .5*dy_coeff/dy  
        end if

        if (downwind_x) then
            coef(I_JP2_K) = coef(I_JP2_K) - 0.5 * dx_coeff/dx
            coef(I_JP1_K) = coef(I_JP1_K) + 2.0 * dx_coeff/dx
            coef(I_J_K)   = coef(I_J_K)   - 1.5 * dx_coeff/dx 
        else if (upwind_x) then
            coef(I_J_K)   = coef(I_J_K)   + 1.5 * dx_coeff/dx
            coef(I_JM1_K) = coef(I_JM1_K) - 2.0 * dx_coeff/dx
            coef(I_JM2_K) = coef(I_JM2_K) + 0.5 * dx_coeff/dx
        else
            coef(I_JP1_K) = coef(I_JP1_K) + 0.5 * dx_coeff/dx
            coef(I_JM1_K) = coef(I_JM1_K) - 0.5 * dx_coeff/dx  
        end if

#if 0
        write(*,*) &
                   "im1jm1 im1km1 im1 im1kp1 im1jp1 jm1km1 jm1 jm1kp1 ", &
                   "km2 km1 center kp1 kp2 jp1km1 jp1 jp1kp1 ", &
                   "ip1jm1 ip1km1 ip1 ip1kp1 ip1jp1 im2 ip2 jm2 jp2 rhs "
        write(*,*)coef,rhs 
#endif
    end subroutine

    !Computes the source term for an ice shelf from hydrostatic and cryostatic pressure
    function source_term(i, j, k, H, zeta, which_source, staggered)
        integer, intent(in) :: i,j,k
        real(dp),dimension(:,:), intent(in) :: H
        real(dp), dimension(:), intent(in) :: zeta
        
        integer, intent(in) :: which_source
        logical, intent(in) :: staggered

        real(dp) :: source_term
        
        real(dp) :: integrate_from

        !Determine the hydrostatic pressure at the current location
        !If we are at the part of the ice shelf above the water, 
        !then we are at 1ATM pressure which we just assume to be 0 
        !throughout the model anyway


        if (which_source == HO_SOURCE_AVERAGED) then
            !Vertically averaged model
            source_term = .5 * rhoi * grav * h(i,j) * (1 - rhoi/rhoo)
        else if (which_source == HO_SOURCE_EXPLICIT) then
#if 1
            !In order to get the correct source term, I integrate from the level above this
            !one to the current level so that all of the pressure is accounted for.
            if (k == 1) then
                source_term = 0
            else
                source_term = .5 * rhoi * grav * H(i,j) * (zeta(k)**2 - zeta(k-1)**2)

                !When we consider the contribution of hydrostatic pressure to the 
                !source term, we only want to integrate that part of the current
                !layer that is below sea level.
                integrate_from = max(1-rhoi/rhoo, zeta(k-1))

                if (zeta(k) > (1 - rhoi/rhoo)) then
                    source_term = source_term + rhoo * grav * H(i,j) * &
                        ( ( 1-rhoi/rhoo) * (zeta(k) - integrate_from) - &
                        .5*(zeta(k)**2 - integrate_from**2) )
                end if
                source_term = source_term / (zeta(k) - zeta(k-1))
            end if
#else
            !Cryostatic component
            source_term = rhoi * grav * h(i,j) * zeta(k)
            !Hydrostatic component, only if we're below the surface
            if (zeta(k) > (1 - rhoi/rhoo)) then
                source_term = source_term + rhoo * grav * h(i,j) * (1 - rhoi/rhoo - zeta(k))
            end if
#endif
        end if  

        !If we are running on a staggered grid, then we need to correct for the averaging that has
        !taken place on the boundary.  Since the thickness is, on average, halved, we double
        !the source term
        !(technique from Stephen Price)

    end function

!
!
!------------------------------------------------------
!   Calculation of stress field based on velocities
!------------------------------------------------------
!
      SUBROUTINE stressf(mu,uvel,vvel,arrh,h,ax,ay,dx,dy,&
        dz,txz,tyz,txx,tyy,txy,FLOWN,ZIP,PERIODIC_X, PERIODIC_Y)
!
        INTEGER MAXY,MAXX,NZETA
        real(dp) :: FLOWN,ZIP
        real(dp), dimension(:,:,:), intent(inout) :: mu
        real(dp), dimension(:,:,:), intent(in) :: uvel
        real(dp), dimension(:,:,:), intent(in) :: vvel
        real(dp), dimension(:,:,:), intent(in) :: arrh
        real(dp), dimension(:,:),   intent(in) :: h
        real(dp), dimension(:,:,:), intent(in) :: ax
        real(dp), dimension(:,:,:), intent(in) :: ay
        real(dp), dimension(:,:,:), intent(out) :: txz 
        real(dp), dimension(:,:,:), intent(out) :: tyz
        real(dp), dimension(:,:,:), intent(out) :: txx
        real(dp), dimension(:,:,:), intent(out) :: tyy
        real(dp), dimension(:,:,:), intent(out) :: txy
        real(dp),                   intent(in) :: dx
        real(dp),                   intent(in) :: dy
        real(dp), dimension(:),     intent(in) :: dz
        
        logical :: PERIODIC_X, PERIODIC_Y

      INTEGER i,j,k
      real(dp), dimension(:,:,:), allocatable :: dudz, dvdz, dudx, dvdx, dudy, dvdy

      real(dp) :: exx,eyy,exy,exz,eyz,eeff
      real(dp) :: macht
!
      
      macht=(1.-FLOWN)/(2.*FLOWN)
      
      MAXY = size(uvel, 1)
      MAXX = size(uvel, 2)
      NZETA = size(dz)
      
            !Allocate temporary derivative data
      allocate(dudx(MAXX, MAXY, NZETA))
      allocate(dudy(MAXX, MAXY, NZETA))
      allocate(dudz(MAXX, MAXY, NZETA))
      allocate(dvdx(MAXX, MAXY, NZETA))
      allocate(dvdy(MAXX, MAXY, NZETA))
      allocate(dvdz(MAXX, MAXY, NZETA))
      
      
      !TODO: Figure out what the UPSTREAM flag does!!
      !Because the values being derived (uvel and vvel) already
      !exist on the staggered grid, we use the unstaggered version
      !of the field derivative function
      call df_field_3d(uvel, dx, dy, dz, dudy, dudx, dudz)
      call df_field_3d(vvel, dx, dy, dz, dvdy, dvdx, dvdz)
      !TODO: Can we somehow avoid the code duplication with muterm?
      do  i=1,MAXY
        do  j=1,MAXX
          do  k=1,NZETA
            if (h(i,j).gt.0.) then
              exx=dudx(j,i,k)+ax(i,j,k)*dudz(j,i,k)
              eyy=dvdy(j,i,k)+ay(i,j,k)*dvdz(j,i,k)
              exy=0.5*((dudy(j,i,k)+ay(i,j,k)*dudz(j,i,k))+(dvdx(j,i,k)+ax(i,j,k)*dvdz(j,i,k)))
              exz=-0.5*dudz(j,i,k)/h(i,j)
              eyz=-0.5*dvdz(j,i,k)/h(i,j)
            else
              exx=dudx(j,i,k)
              eyy=dvdy(j,i,k)
              exy=0.5*(dudy(j,i,k)+dvdx(j,i,k))
              exz=0.
              eyz=0.
            endif
            eeff=(exx*exx)+(eyy*eyy)+(exx*eyy)+(exy*exy)+(exz*exz)+&
              (eyz*eyz)
            if (eeff < ZIP) eeff=ZIP
            mu(i,j,k)=0.5*(arrh(i,j,k)**(-1./FLOWN))*(eeff**macht)
            txz(i,j,k)=2.0*mu(i,j,k)*exz
            tyz(i,j,k)=2.0*mu(i,j,k)*eyz
            txx(i,j,k)=2.0*mu(i,j,k)*exx
            tyy(i,j,k)=2.0*mu(i,j,k)*eyy
            txy(i,j,k)=2.0*mu(i,j,k)*exy
        end do
      end do
    end do


      !Adjust boundaries if using periodic boundary conditions.
      !REFACTORED CODE BELOW
      call periodic_boundaries_3d_stag(txz,PERIODIC_X,PERIODIC_Y)
      call periodic_boundaries_3d_stag(txy,PERIODIC_X,PERIODIC_Y)
      call periodic_boundaries_3d_stag(txx,PERIODIC_X,PERIODIC_Y)
      call periodic_boundaries_3d_stag(tyy,PERIODIC_X,PERIODIC_Y)
      call periodic_boundaries_3d_stag(txy,PERIODIC_X,PERIODIC_Y)
      
            !Free temporary derivative data
      deallocate(dudx)
      deallocate(dudy)
      deallocate(dudz)
      deallocate(dvdx)
      deallocate(dvdy)
      deallocate(dvdz)
      
      return
      END subroutine
      
    function staggered_point_2d(f, i, j)
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i, j
        real(dp) :: staggered_point_2d
        
        staggered_point_2d = (f(i, j) + f(i, j+1) + f(i+1, j) + f(i+1, j+1))/4
        
    end function staggered_point_2d
    
    function staggered_point_3d(f,i,j,k)
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i, j, k
        real(dp) :: staggered_point_3d
        
        staggered_point_3d = ( f(i, j, k) + f(i, j+1, k) + f(i+1, j, k) + f(i+1, j+1, k) )/4
    
    !*FD Moves a field to a staggered grid version of that field
    !*FD This is done simply by averaging each point
    !*FD f_stag must have 1 fewer element in each dimension
    end function staggered_point_3d
    
    subroutine staggered_field_2d(f, f_stag)
        !See the note on staggered grids and periodic boundary conditions at the top
        !of this code file.
        !We will assume that f is the same size as f_stag.  If periodic boundary
        !conditions are not used, then the last row and column of f_stag will be
        !set to zero.  Otherwise, the staggered grid values will be set for them
        !such that the last row and column represents the staggered values
        !that occur between periodic repetitions of the nonstaggered grid.
        real(dp), dimension(:,:), intent(in)  :: f
        real(dp), dimension(:,:), intent(out) :: f_stag
        integer :: nx, ny, i, j
        
        nx = size(f, 1)
        ny = size(f, 2)
        
        !f_stag = (f(1:nx-1, 1:ny-1) + f(2:nx, 1:ny-1) + f(1:nx-1, 2:ny) + f(2:nx, 2:ny))/4
        do i = 1,nx-1
            do j = 1,ny-1
                f_stag(i,j) = staggered_point_2d(f, i, j)
            end do
        end do
        
!       if (periodic_x) then
!               do j = 1, ny-1
!                       f_stag(nx, j) = (f(nx, j) + f(nx, j+1) + f(1, j) +f(1,j+1))/4
!               end do
!       else
!               !DANGER: Did I transpose X and Y?
!               f_stag(nx,:) = 0
!       end if
!       
!       if (periodic_y) then
!               do i = 1, nx-1
!                       f_stag(i, ny) = (f(i, ny) + f(i+1, ny) + f(i, 1) + f(i+1, 1))/4
!               end do
!       else
!               f_stag(:, ny) = 0
!       end if
!       
!       !Do the corner only if both periodic boundaries are used
!       if (periodic_x .and. periodic_y) then
!               f_stag(nx, ny) = (f(nx, ny) + f(1, ny) + f(nx, 1) + f(1,1))/4
!       else
!               f_stag(nx, ny) = 0
!       end if
!       
        
    end subroutine staggered_field_2d
    
    !*FD Moves a field to a staggered grid version of that field
    !*FD This is done simply by averaging each point
    !*FD f_stag must have 1 fewer element in the x and y dimensions, and the same
    !*FD number of elements in the z dimension
    subroutine staggered_field_3d(f, f_stag)
    
        real(dp), dimension(:,:,:), intent(in)  :: f
        real(dp), dimension(:,:,:), intent(out) :: f_stag
        integer :: nz, k

        nz = size(f, 3)
    
        do k = 1, nz
            call staggered_field_2d(f(:,:,k), f_stag(:,:,k))
        end do
    end subroutine staggered_field_3d
    
    !*FD Copies a staggered grid onto a nonstaggered grid.  This verion
    !*FD assumes periodic boundary conditions.
    subroutine unstagger_field_2d(f_stag, f, periodic_x, periodic_y)
        real(dp), dimension(:,:), intent(in) :: f_stag
        real(dp), dimension(:,:), intent(out) :: f
        logical, intent(in) :: periodic_x, periodic_y

        real(dp), dimension(4) :: pts

        real(dp) :: s,n

        integer :: i,j, k,i1, i2, j1, j2, ni, nj
        
        ni = size(f, 1)
        nj = size(f, 2)

        do i = 1, size(f, 1)
            do j = 1, size(f, 2)
                s = 0
                n = 0
                
                i1 = i-1
                i2 = i

                if (i1 == 0) then
                    if (periodic_y) then
                        i1 = ni - 1
                    else
                        i1 = 1
                    end if
                end if
    
                if (i2 == ni) then
                    if (periodic_y) then
                        i2 = 1
                    else
                        i2 = ni - 1
                    end if
                end if
    
                j1 = j-1
                j2 = j
    
                if (j1 == 0) then
                    if (periodic_x) then
                        j1 = nj - 1
                    else
                        j1 = 1
                    end if
                end if
    
                if (j2 == nj) then
                    if (periodic_x) then
                        j2 = 1
                    else
                        j2 = nj - 1
                    end if
                end if
                
                !Place the points into an array, loop over them, and average
                !all the points that AREN'T NaN.
                pts = (/f_stag(i1, j1), f_stag(i2, j1), f_stag(i1, j2), f_stag(i2, j2)/)
            
                do k=1,4
                    if (.not. (IS_NAN(pts(k)))) then
                        s = s + pts(k)
                        n = n + 1
                    end if
                end do
                if (n /= 0) then
                    f(i,j) = s/n
                else
                    f(i,j) = NaN
                end if
            end do
        end do
    
    end subroutine unstagger_field_2d

    subroutine unstagger_field_3d(f, f_stag, periodic_x, periodic_y)
        real(dp), dimension(:,:,:) :: f, f_stag
        logical, intent(in) :: periodic_x, periodic_y

        integer :: i

        do i = 1,size(f,3)
            call unstagger_field_2d(f(:,:,i), f_stag(:,:,i), periodic_x, periodic_y)
        end do
        
    end subroutine unstagger_field_3d

#ifdef OUTPUT_PARTIAL_ITERATIONS
!------------------------------------------------------------------------------
!The following are methods to write iterations to a NetCDF file.  
!The NetCDF file used is not managed through Glimmer's NetCDF interface
!I know this is code duplication, it is for DEBUG PURPOSES ONLY
!------------------------------------------------------------------------------
function begin_iteration_debug(nx, ny, nz)
    use netcdf
    use glimmer_ncdf, only: nc_errorhandle
    integer :: nx, ny, nz
    integer :: err
    integer :: ncid

    integer :: xdim, ydim, zdim, iterdim, dims(4), dims2d(3), varid
    integer :: begin_iteration_debug

    ncid = 0
    err = nf90_create("iterdebug.nc", NF90_CLOBBER, ncid)
    call nc_errorhandle(__FILE__, __LINE__, err)
    err = nf90_def_dim(ncid, "x", nx, xdim)
    call nc_errorhandle(__FILE__, __LINE__, err)
    err = nf90_def_dim(ncid, "y", ny, ydim)
    call nc_errorhandle(__FILE__, __LINE__, err)
    err = nf90_def_dim(ncid, "z", nz, zdim)
    call nc_errorhandle(__FILE__, __LINE__, err)
    err = nf90_def_dim(ncid, "iter", NF90_UNLIMITED, iterdim)
    call nc_errorhandle(__FILE__, __LINE__, err)

    dims = (/xdim, ydim, zdim, iterdim/)
    dims2d = (/xdim, ydim, iterdim/)

    err = nf90_def_var(ncid, "mu", NF90_DOUBLE, dims, varid)
    call nc_errorhandle(__FILE__, __LINE__, err)
    err = nf90_def_var(ncid, "uvel", NF90_DOUBLE, dims, varid)
    call nc_errorhandle(__FILE__, __LINE__, err)
    err = nf90_def_var(ncid, "vvel", NF90_DOUBLE, dims, varid)
    call nc_errorhandle(__FILE__, __LINE__, err) 
    err = nf90_def_var(ncid, "velnorm", NF90_DOUBLE, dims, varid)
    call nc_errorhandle(__FILE__, __LINE__, err) 
    err = nf90_def_var(ncid, "mask", NF90_INT, dims2d, varid)
    call nc_errorhandle(__FILE__, __LINE__, err) 

    err = nf90_def_var(ncid, "dudx", NF90_DOUBLE, dims, varid)
    call nc_errorhandle(__FILE__, __LINE__, err) 
    err = nf90_def_var(ncid, "dudy", NF90_DOUBLE, dims, varid)
    call nc_errorhandle(__FILE__, __LINE__, err) 
    err = nf90_def_var(ncid, "dvdx", NF90_DOUBLE, dims, varid)
    call nc_errorhandle(__FILE__, __LINE__, err) 
    err = nf90_def_var(ncid, "dvdy", NF90_DOUBLE, dims, varid)
    call nc_errorhandle(__FILE__, __LINE__, err) 
  
    err = nf90_enddef(ncid)
    call nc_errorhandle(__FILE__, __LINE__, err)

    begin_iteration_debug=ncid


end function

subroutine iteration_debug_step(ncid, iter, mu, uvel, vvel, geometry_mask)
    use netcdf
    use glimmer_ncdf, only: nc_errorhandle

    integer :: ncid, iter
    real(dp), dimension(:,:,:) :: mu, uvel, vvel
    integer, dimension(:,:) :: geometry_mask

    integer :: varid, err

    integer :: nx, ny, nz
    integer :: i,j,k, start(4), count(4)
    
    nx = size(mu, 2)
    ny = size(mu, 1)
    nz = size(mu, 3)

    count = (/1,1,1,1/)
    do i = 1,ny
        do j = 1,nx
            do k = 1,nz
                start=(/j,i,k,iter+1/)
       err = nf90_inq_varid(ncid, "mu", varid)
       call nc_errorhandle(__FILE__, __LINE__, err)
       
       err = nf90_put_var(ncid, varid, mu(i,j,k), start)
       call nc_errorhandle(__FILE__, __LINE__, err)
       
       err = nf90_inq_varid(ncid, "uvel", varid)
       call nc_errorhandle(__FILE__, __LINE__, err) 
       
       err = nf90_put_var(ncid, varid, uvel(i,j,k), start)
       
       err = nf90_inq_varid(ncid, "vvel", varid)
       call nc_errorhandle(__FILE__, __LINE__, err)
       
       err = nf90_put_var(ncid, varid, vvel(i,j,k), start)
       call nc_errorhandle(__FILE__, __LINE__, err)
 
       err = nf90_inq_varid(ncid, "velnorm", varid)
       call nc_errorhandle(__FILE__, __LINE__, err)
       
       err = nf90_put_var(ncid, varid, sqrt(uvel(i,j,k)**2 + vvel(i,j,k)**2), start)
       call nc_errorhandle(__FILE__, __LINE__, err)

 
       err = nf90_inq_varid(ncid, "mask", varid)
       call nc_errorhandle(__FILE__, __LINE__, err)
       
       err = nf90_put_var(ncid, varid, geometry_mask(i,j), (/j,i,iter+1/))
       call nc_errorhandle(__FILE__, __LINE__, err)
            end do
        end do
    end do
   
    !Close and open the dataset so that changes get written to it
    call end_debug_iteration(ncid)

    err = nf90_open("iterdebug.nc", ior(NF90_WRITE,NF90_SHARE), ncid)
    call nc_errorhandle(__FILE__, __LINE__, err)

end subroutine

subroutine iterdebug_vel_derivs(ncid, iter, dudx, dudy, dvdx, dvdy)
    use netcdf
    use glimmer_ncdf, only: nc_errorhandle

    integer :: ncid, iter
    real(dp), dimension(:,:,:) :: dudx, dudy, dvdx, dvdy

    integer :: varid, err

    integer :: nx, ny, nz
    integer :: i,j,k, start(4), count(4)
    
    nx = size(dudx, 2)
    ny = size(dudx, 1)
    nz = size(dudx, 3)

    count = (/1,1,1,1/)
    do i = 1,ny
        do j = 1,nx
            do k = 1,nz
                start=(/j,i,k,iter+1/)
                err = nf90_inq_varid(ncid, "dudx", varid)
                call nc_errorhandle(__FILE__, __LINE__, err)
       
                err = nf90_put_var(ncid, varid, dudx(i,j,k), start)
                call nc_errorhandle(__FILE__, __LINE__, err)

                err = nf90_inq_varid(ncid, "dudy", varid)
                call nc_errorhandle(__FILE__, __LINE__, err)
       
                err = nf90_put_var(ncid, varid, dudy(i,j,k), start)
                call nc_errorhandle(__FILE__, __LINE__, err)
       
                err = nf90_inq_varid(ncid, "dvdx", varid)
                call nc_errorhandle(__FILE__, __LINE__, err)
       
                err = nf90_put_var(ncid, varid, dvdx(i,j,k), start)
                call nc_errorhandle(__FILE__, __LINE__, err)
       
                err = nf90_inq_varid(ncid, "dvdy", varid)
                call nc_errorhandle(__FILE__, __LINE__, err)
       
                err = nf90_put_var(ncid, varid, dvdy(i,j,k), start)
                call nc_errorhandle(__FILE__, __LINE__, err)
       
 
            end do
        end do
    end do
    


end subroutine

subroutine end_debug_iteration(ncid)
    use netcdf
    use glimmer_ncdf, only: nc_errorhandle

    integer :: ncid

    integer :: err

    err = nf90_close(ncid)
    call nc_errorhandle(__FILE__, __LINE__, err)
end subroutine
#endif

subroutine write_xls_direction_guide(filename, nx, ny)
    integer :: nx, ny
    character(*) :: filename
    
    real(dp), dimension(nx, ny) :: field

    field = 0
    field(1,1) = 1
    field(2,1) = 1

    call write_xls(filename, field)

end subroutine

subroutine open_sparse_output(iteration, iunit)
    integer :: iteration
    integer :: iunit

    character(32) filename
    character(2)  istring

    write(istring, '(i2)') iteration
    write(filename,*) "sparse"//istring//".txt"
    open(unit=iunit, file=filename)
end subroutine

!Detects locations that are surrounded entirely by 
subroutine mask_out_edges(geometry_mask, point_mask, uvelbc, vvelbc)
    integer, dimension(:,:), intent(in) :: geometry_mask
    integer, dimension(:,:), intent(in) :: point_mask
    real(dp), dimension(:,:,:), intent(inout) :: uvelbc, vvelbc
    
    logical, dimension(size(uvelbc, 1), size(uvelbc, 2)) :: is_potential_point

    integer :: i, j

    !PASS 1: Look for points that potentially need to be masked out (i.e. surrounded by no ice
    !or by dirichlet boundary)
    do i = 2, size(uvelbc, 1)-1
        do j = 2, size(uvelbc, 1)-1
            is_potential_point(i,j) = needs_mask(i,j,1)
        end do
    end do

    !PASS 2: If a point is surrounded entirely by points that might need to be masked out,
    !points with no ice, or dirichlet boundary conditions, then set a 0-value dirichlet
    !condition there as it should not contribute to the solution.
    do i = 2, size(uvelbc, 1)-1
        do j = 2, size(uvelbc, 1)-1
            if (needs_mask(i,j,2)) then
                uvelbc(i,j,:) = 0
                vvelbc(i,j,:) = 0
            end if
        end do
    end do
contains
    function needs_mask(i,j,pass)
        integer, intent(in) :: i, j
        integer, intent(in) :: pass

        integer :: i2, j2
        
        integer, dimension(-1:1, -1:1) :: local

        logical :: needs_mask

        local = 0
        do i2 = -1, 1
            do j2 = -1, 1
                if (i2 /= 0 .and. j2 /= 0) then
                    if (.not. IS_NAN(uvelbc(i+i2, j+j2,1))) then
                        local(i2, j2) = 1
                    end if

                    if (point_mask(i+i2, j+j2) == 0) then
                        local(i2, j2) = 1
                    end if
                end if

                if (pass == 2 .and. is_potential_point(i+i2, j+j2)) then
                    local(i2, j2) = 1
                end if
            end do
        end do

        needs_mask = .false.
        if (pass == 1 .and. sum(local) > 0) then
            needs_mask = .true.
        end if

        if (pass == 2 .and. sum(local) == 9) then
            needs_mask = .true.
        end if
    end function
end subroutine

!
!---------------------------------------------------
!
!                 END OF SUBROUTINES
!
!---------------------------------------------------
end module ice3d_lib

