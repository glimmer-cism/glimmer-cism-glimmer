module ice3d_lib
    use glimmer_global
    use glimmer_physcon, only: pi, grav,rhoi
    use glide_deriv
    use glimmer_sparse
    use glimmer_sparse_solver
    use xls
    implicit none
    double precision :: small, zip, notdef
    PARAMETER(SMALL=1.D-10,ZIP=1.D-30)
    PARAMETER(NOTDEF=-999999999.)
    integer,parameter :: dslucs_dbg_unit = 6

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

contains
!
!
!     *****************************************
!     *                                       *
!     *         Start of subroutines          *
!     *                                       *
!     *****************************************
!
!
!------------------------------------------------------
!   Copying matrices. This subroutine copies the newly
!   calculated ice thickness on the old ones
!------------------------------------------------------
!
      SUBROUTINE mapmatrix(h,hnew,MAXY,MAXX)
!
        INTEGER MAXX,MAXY
        double precision :: h(MAXY,MAXX),hnew(MAXY,MAXX)
!
      INTEGER i,j
!
      do 10 i=1,MAXY
        do 20 j=1,MAXX
          h(i,j)=hnew(i,j)
   20   CONTINUE
   10 CONTINUE
      return
      END subroutine
!
!
!------------------------------------------------------
!   Surface mass balance parameterization. This routine
!   should be altered by the user
!------------------------------------------------------
!
      SUBROUTINE environment(massb,MAXY,MAXX)
!
        INTEGER MAXX,MAXY
        double precision :: massb(MAXY,MAXX)
!
      INTEGER i,j
!
      do 10 i=1,MAXY
        do 20 j=1,MAXX
          massb(i,j)=0.3
   20   CONTINUE
   10 CONTINUE
      return
      END subroutine
!
!
!------------------------------------------------------
!   Thermomechanical coupling: determination of 'arrh'.
!   In this example this factor is kept constant as
!   temperature coupling is not considered
!------------------------------------------------------
!
      SUBROUTINE tcoupling(arrh,MAXY,MAXX,NZETA,SHTUNE)
!
        INTEGER MAXX,MAXY,NZETA
        double precision :: arrh(MAXY,MAXX,NZETA),SHTUNE
!
      INTEGER i,j,k
!
      do 10 i=1,MAXY
        do 20 j=1,MAXX
          do 30 k=1,NZETA
            arrh(i,j,k)=SHTUNE
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
      return
      END subroutine
      
subroutine init_zeta(zeta)
    double precision, dimension(:), intent(out) :: zeta
    integer :: nzeta, k

    nzeta = size(zeta,1)

    do k=1,NZETA
        !kprime = 41D0*k / nzeta
        !zeta(k)=(-2.5641025641D-4)*(kprime-1.)**2+3.5256410256D-2*(kprime-1.)-&
        !  8.0047080075D-13
        zeta(K) = (k-1D0)/(nzeta-1D0)
    end do
    zeta(1)=0.
    zeta(NZETA)=1
end subroutine init_zeta
 
!
!
!------------------------------------------------------
!   Initialization of geometric vectors and gradients.
!   h, hb, and surf are given on unstaggered grids, while
!   the derivatives are given on staggered grids
!------------------------------------------------------
!
    SUBROUTINE init_geometry_vectors(dhdx,dhbdx,dhdy,dhbdy,surf,h,hb,dzdx,dzdy,dx,dy)
!
        INTEGER MAXX,MAXY
        double precision, dimension(:,:), intent(out) :: dhdx
        double precision, dimension(:,:), intent(out) :: dhbdx
        double precision, dimension(:,:), intent(out) :: dhdy
        double precision, dimension(:,:), intent(out) :: dhbdy
        double precision, dimension(:,:), intent(out) :: surf
        double precision, dimension(:,:), intent(in)  :: h
        double precision, dimension(:,:), intent(in)  :: hb
        double precision, dimension(:,:), intent(out) :: dzdx
        double precision, dimension(:,:), intent(out) :: dzdy
        
        double precision, intent(in) :: dx
        double precision, intent(in) :: dy
        
        MAXX = size(h, 1)
        MAXY = size(h, 2)
        write(*,*) "BEGIN init_geometry_vectors"
        !Compute derivative fields for bed and thickness onto
        !a staggered grid.
        !Note that, because of the different indexing schemes, we switch the
        !locations where we place dhdy and dhdx.  We also need to switch the
        !x and y grid sizes given
        call df_field_2d_staggered(h,  dy, dx, dhdy,  dhdx, .true., .true.)
        call df_field_2d_staggered(hb, dy, dx, dhbdy, dhbdx, .true., .true.)
      

        !Compute surface, derivatives of surface
        !Because of transpose issues, dzdy has to be the derivatives
        !w.r.t x added, and vice versa.
        !I won't pretend to understand why!
        surf = h+hb
        dzdx = dhdx+dhbdx
        dzdy = dhdy+dhbdy !TODO: Still getting ~25% error here!!  Looks like floating point roundoff 
        write(*,*) "END init_geometry_vectors" 
    END subroutine
      
      
    !Initialize rescaled coordinate coefficients
    subroutine init_rescaled_coordinates(dhdx,dhbdx,dhdy,dhbdy,surf,h,hb,&
               dzdx,dzdy,d2zdx2, d2zdy2, d2hdx2, d2hdy2, zeta,ax,ay,bx,by,cxy,dx,dy)
!
        INTEGER MAXX,MAXY,NZETA
        double precision, dimension(:,:), intent(in) :: dhdx
        double precision, dimension(:,:), intent(in) :: dhbdx
        double precision, dimension(:,:), intent(in) :: dhdy
        double precision, dimension(:,:), intent(in) :: dhbdy
        double precision, dimension(:,:), intent(in) :: surf
        double precision, dimension(:,:), intent(in) :: h
        double precision, dimension(:,:), intent(in) :: hb
        double precision, dimension(:,:), intent(in) :: dzdx
        double precision, dimension(:,:), intent(in) :: dzdy
        double precision, dimension(:),   intent(in) :: zeta
        double precision, dimension(:,:,:) :: ax
        double precision, dimension(:,:,:) :: ay
        double precision, dimension(:,:,:) :: bx
        double precision, dimension(:,:,:) :: by
        double precision, dimension(:,:,:) :: cxy
        
        double precision, intent(in) :: dx
        double precision, intent(in) :: dy
        
        double precision, dimension(:,:), intent(in) :: d2zdx2,d2zdy2,d2hdx2,d2hdy2
        
        INTEGER :: i,j,k
#if 1
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
                                      zeta(k) * d2hdy2(i, j) - & 
                                      2. * ax(i, j, k) * dhdx(i, j) ) / h(i,j)
                         
                        by(i,j,k) = ( d2zdy2(i, j) - &
                                      zeta(k) * d2hdx2(i, j) - & 
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
        stop
#endif
    end subroutine
!
!
!------------------------------------------------------
!   ice-sheet equation
!------------------------------------------------------
!
    SUBROUTINE spice2d(h,hnew,hb,diffus,u,v,massb,MAXY,MAXX,&
        IJTOT,GRIDX,GRIDY,TSTEP,SMALL,PERIODIC,DIFTYPE)
!
        INTEGER MAXX,MAXY,IJTOT,PERIODIC,DIFTYPE
        double precision :: h(MAXY,MAXX),hnew(MAXY,MAXX),hb(MAXY,MAXX)
        double precision :: diffus(MAXY,MAXX),massb(MAXY,MAXX)
        double precision :: u(MAXY,MAXX),v(MAXY,MAXX)
        double precision :: GRIDX,GRIDY,TSTEP,SMALL
!
        INTEGER i,j,k,m,iter,cspice
        double precision :: d(IJTOT),coef(6),x(IJTOT),hmean
        double precision :: err,dipx,dipy,dimx,dimy,dtdx,dtdy,tol
        double precision :: alfa,beta,gamma,delta,admx,adpx,admy,adpy
        double precision :: dtdx2,dtdy2
      
      
        double precision, dimension(:), allocatable :: sa
        integer, dimension(:), allocatable :: ija

        !Allocate sparse matrix storage
        allocate(sa(ijtot*22))
        allocate(ija(ijtot*22))
!
   
        ija = 0
        sa  = 0
        d   = 0
        x   = 0
   
        tol=1.e-5
        k=IJTOT+2;
        dtdx=TSTEP/(2.*GRIDX*GRIDX)
        dtdy=TSTEP/(2.*GRIDY*GRIDY)
        dtdx2=TSTEP/(4.*GRIDX)
        dtdy2=TSTEP/(4.*GRIDY)
        if (DIFTYPE.eq.1) then
            alfa=1.
            beta=1.
            gamma=1.
            delta=1.
        else
            alfa=0.
            beta=0.
            gamma=0.
            delta=0.
        endif
        
        do 10 i=1,MAXY
            do 20 j=1,MAXX
                coef = 0
          
                if (i.eq.1.or.i.eq.MAXY.or.j.eq.1.or.j.eq.MAXX) then
                    coef(3)=1.
                    coef(6)=h(i,j)
                else
                    dipx=dtdx*(diffus(i,j+1)+diffus(i,j))
                    dimx=dtdx*(diffus(i,j)+diffus(i,j-1))
                    dipy=dtdy*(diffus(i+1,j)+diffus(i,j))
                    dimy=dtdy*(diffus(i,j)+diffus(i-1,j))
                    adpx=dtdx2*(u(i,j)+u(i,j+1))
                    admx=dtdx2*(u(i,j)+u(i,j-1))
                    adpy=dtdy2*(v(i,j)+v(i+1,j))
                    admy=dtdy2*(v(i,j)+v(i-1,j))
                    coef(1)=(-delta*dimy)+(delta-1.)*admy
                    coef(2)=(-beta*dimx)+(beta-1.)*admx
                    coef(3)=1.+alfa*dipx+beta*dimx+gamma*dipy+delta*dimy-&
                        (alfa-1.)*adpx+(beta-1.)*admx-(gamma-1.)*adpy+&
                        (delta-1.)*admy
                    coef(IM1_J_KP1)=(-alfa*dipx)-(alfa-1.)*adpx
                    coef(5)=(-gamma*dipy)-(gamma-1.)*adpy
                    coef(6)=h(i,j)+alfa*dipx*(hb(i,j+1)-hb(i,j))-beta*&
                        dimx*(hb(i,j)-hb(i,j-1))+gamma*dipy*(hb(i+1,j)-&
                        hb(i,j))-delta*dimy*(hb(i,j)-hb(i-1,j))+&
                        massb(i,j)*TSTEP
                endif
                sa(cspice(3,i,j,MAXX))=coef(3)
                d(cspice(3,i,j,MAXX))=coef(6)
                x(cspice(3,i,j,MAXX))=h(i,j)
                ija(cspice(3,i,j,MAXX))=k
                do 40 m=1,5
                    if (m.ne.3.and.abs(coef(m)).gt.SMALL) then
                        sa(k)=coef(m)
                        ija(k)=cspice(m,i,j,MAXX)
                        k=k+1
                    endif
   40           continue
   20       continue
   10   continue
        ija(MAXX*MAXY+1)=k
        !CALL linbcg(IJTOT,d,x,1,tol,200,iter,err,sa,ija)
        k=0
        hmean=0.
        do 50 i=1,MAXY
            do 60 j=1,MAXX
                k=k+1
                hnew(i,j)=x(k)
                if (hnew(i,j).lt.0.)hnew(i,j)=0.
                hmean=hmean+hnew(i,j)
   60       continue
   50   continue
        if (PERIODIC.eq.1) then
            do 100 i=1,MAXY
                hnew(i,1)=hnew(i,MAXX-1)
                hnew(i,MAXX)=hnew(i,2)
  100       CONTINUE
            do 110 j=1,MAXX
                hnew(1,j)=hnew(MAXY-1,j)
                hnew(MAXY,j)=hnew(2,j)
  110       CONTINUE
        endif
        write(*,*)iter,hmean/IJTOT,massb(10,10)
      
        deallocate(sa)
        deallocate(ija)
      
        return
    end subroutine
!
!
!------------------------------------------------------
!   Diffusivities according to 0th order model (SIA)
!------------------------------------------------------
!
    SUBROUTINE diffus1(dzdx,dzdy,h,diff1,zeta,arrh,diffus,&
        MAXY,MAXX,NZETA,FLOWN,PERIODIC_X, PERIODIC_Y)

        INTEGER MAXX,MAXY,NZETA
        double precision :: dzdx(MAXY,MAXX),dzdy(MAXY,MAXX)
        double precision :: h(MAXY,MAXX),diff1(MAXY,MAXX,NZETA)
        double precision :: zeta(NZETA),arrh(MAXY,MAXX,NZETA)
        double precision :: diffus(MAXY,MAXX),FLOWN

        logical :: periodic_x, periodic_y
        INTEGER i,j,k
        double precision :: grad,d,z
   
        diff1 = 0

        do 10 i=2,MAXY-1
            do 20 j=2,MAXX-1
                grad=dzdx(i,j)*dzdx(i,j)+dzdy(i,j)*dzdy(i,j)
                if (h(i,j).gt.0.) then
                    d=(RHOI*GRAV*h(i,j))**FLOWN
                else
                    d=0.
                endif
                do 30 k=(NZETA-1),1,-1
                    z=(0.5*(zeta(k+1)+zeta(k)))**FLOWN
                    diff1(i,j,k)=d*h(i,j)*grad*(arrh(i,j,k+1)+arrh(i,j,k))*&
                        z*(zeta(k+1)-zeta(k))+diff1(i,j,k+1)
   30           CONTINUE
                diffus(i,j)=0.
                do 40 k=2,NZETA
                    diffus(i,j)=diffus(i,j)+0.5*(diff1(i,j,k)+&
                        diff1(i,j,k-1))*(zeta(k)-zeta(k-1))
   40           CONTINUE
                diffus(i,j)=diffus(i,j)*h(i,j)
   20       CONTINUE
   10   CONTINUE

        if (PERIODIC_X) then
            do 100 i=1,MAXY
                diffus(i,MAXX)=diffus(i,2)
                diffus(i,1)=diffus(i,MAXX-1)
                do 110 k=1,NZETA
                    diff1(i,MAXX,k)=diff1(i,2,k)
                    diff1(i,1,k)=diff1(i,MAXX-1,k)
  110           CONTINUE
  100       CONTINUE
        else
            do 140 i=1,MAXY
                diffus(i,MAXX)=0.
                diffus(i,1)=0.
                do 150 k=1,NZETA
                    diff1(i,MAXX,k)=0.
                    diff1(i,1,k)=0.
  150           CONTINUE
  140       CONTINUE
        endif

        if (PERIODIC_Y) then
            do 120 j=1,MAXX
                diffus(MAXY,j)=diffus(2,j)
                diffus(1,j)=diffus(MAXY-1,j)
                do 130 k=1,NZETA
                    diff1(MAXY,j,k)=diff1(2,j,k)
                    diff1(1,j,k)=diff1(MAXY-1,j,k)
  130           CONTINUE
  120       CONTINUE
        else
            do 160 j=1,MAXX
                diffus(MAXY,j)=0.
                diffus(1,j)=0.
                do 170 k=1,NZETA
                    diff1(MAXY,j,k)=0.
                    diff1(1,j,k)=0.
  170           CONTINUE
  160       CONTINUE
        endif
 
        return
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
        double precision :: dzdx(MAXY,MAXX)
        double precision :: dzdy(MAXY,MAXX),h(MAXY,MAXX),zeta(NZETA)
        double precision :: arrh(MAXY,MAXX,NZETA),uvel(MAXY,MAXX,NZETA)
        double precision :: vvel(MAXY,MAXX,NZETA),u(MAXY,MAXX)
        double precision :: v(MAXY,MAXX),FLOWN
        logical :: PERIODIC_X, PERIODIC_Y
!
      INTEGER i,j,k
      double precision :: diff1(NZETA),d,grad,z,diffus,us,vs
!
      diff1 = 0
     
      do 10 i=2,MAXY-1
        do 20 j=2,MAXX-1
          grad=(dzdx(i,j)**2)+(dzdy(i,j)**2)
          us=0.
          vs=0.
          d = (RHOI*GRAV*h(i,j))**FLOWN
          do 30 k=(NZETA-1),1,-1
            z = (0.5*(zeta(k+1)+zeta(k)))**FLOWN
            diff1(k)=d*grad*h(i,j)*(arrh(i,j,k+1)+arrh(i,j,k))*z*&
              (zeta(k)-zeta(k+1))+diff1(k+1)
   30     CONTINUE
          diffus=0.
          do 40 k=2,NZETA
            diffus=diffus+0.5*(diff1(k)+diff1(k-1))*&
              (zeta(k)-zeta(k-1))
   40     CONTINUE
          do 50 k=1,NZETA
            uvel(i,j,k)=diff1(k)*dzdx(i,j)+us
            vvel(i,j,k)=diff1(k)*dzdy(i,j)+vs
   50     CONTINUE
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
        double precision, dimension(:,:), intent(inout) :: m
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
        double precision, dimension(:,:), intent(inout) :: m
        integer :: maxx, maxy
        logical :: apply_to_x, apply_to_y

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

        !UNCOMMENT for periodic boundary conditions set so that they reflect
        !"standard" treatment of the periodic conditions on the unstaggered grid
    
!               do i = 1, maxy
!                       
!                       !Using the identities: 
!                       !x_2 = (x_1.5 + x_2.5)/2 
!                       !x_{n-1} = (x_{n-.5} + x_{n-1.5})/2
!                       !x'_1.5 = x'_{n-.5} = (x_2 + x_{n-1})/2
!                       !                   = (x_1.5 + x_2.5 + x_{n-1.5} + x_{n-.5})/4
!                       new_boundary = (m(i,1) + m(i,2) + m(i,maxx-1) + m(i,maxx))/4
!                       
!                       !Because the periodic boundary conditions on the nonstaggered
!                       !grid influence the first and last two elements on the staggered
!                       !grid, we need to adjust the second boundary condition
!                       !Starting with the identity x_2.5 = (x_2 + x_3)/4
!                       != (x_1.5 + 2x_2.5 + x_3.5)/4
!                       !And x'_2.5 = (x'_1.5 + 2x'_2.5 + x_3.5)/4
!                       !We can solve x'_2.5 = x_2.5 + (x'_1.5 - x_1.5)/4
!                       m(i,2) = m(i,2) + (new_boundary - m(i,1))/4
!                       m(i,maxx-1) = m(i,maxx-1) + (new_boundary - m(i,maxx))/4
!                       
!                       !Now, place the new boundary conditions in the first and last
!                       !elements of this row
!                       m(i,1)     = new_boundary
!                       m(i, maxx) = new_boundary
!                       
!                        
!               end do
!               
!               !Do the same as the above for the y dimension
!               do i = 1, maxx
!                       new_boundary = (m(1,i) + m(2,i) + m(maxy-1,i) + m(maxy,i))/4
!                       
!                       m(2,i) = m(2,i) + (new_boundary - m(1,i))/4
!                       m(maxy-1, i) = m(maxy-1,i) + (new_boundary - m(maxy,i))/4
!                       
!                       m(1,i)    = new_boundary
!                       m(maxy,i) = new_boundary
!               end do

    !m(:,1) = (m(:,1) + m(:,2) + m(:,maxx-2) + m(:,maxx-1))/4
    !m(:,maxx-1) = m(:,1)

    !m(1,:) = (m(1,:) + m(2,:) + m(maxy-2,:) + m(maxy-1,:))/4
    !m(maxy-1,:) = m(1,:)
    
    
!               
    end subroutine periodic_boundaries_stag


!       subroutine constant_boundaries_stag(m, minx, maxx, miny, maxy)
!               double precision, dimension(:,:) :: m
!               double precision :: minx, maxx, miny, maxy
!               
!               !In order to apply a boundary condition to the constant boundaries,
!               !we need to pretend that we are applying the boundary condition on
!               !the unstaggered grid.  However, to avoid over-averaging the
!               !boundary conditions, we will control for the current "hidden" value
!               !of the boundary on the unstaggered grid.
!               !If we take x'_1.5 = (x'_1 + x_2)/2 and x_1.5 = (x_1 + x_2)/2,
!               !then x'_1.5 = x_1.5 + (x'_1 - x_1) / 2program HOM3D
!               !What we need is an expression for x_1, which can be found as
!               !x_1 = 2x_1.5 - x_2 (rearranging the identity for x_1.5 above)
!               !Substituting this, we get that
!               !x'_1.5 = x_1.5 + (x'_1 - 2x_1.5 + x2)/2
!               !However, we also don't know x2!  Substituting an interpolated
!               !value for x2, we get:
!               !x'_1.5 = 1/4*x_1.5 + 1/2*x'_1 + 1/4*x_2.5
!               m(1,:) = .5*minx + .25*m(1,:) + .25*m(2,:)
!               m(size(m,1),:) = .5*maxx + .25*m(size(m,1),:) + .25*m(size(m,1)-1,:)
!               
!               m(:,1) = .5*miny + .25*m(:,1) + .25*m(:,2)
!               m(size(m,1),:) = .5*maxy + .25*m(:,size(m,2)) + .25*m(:,size(m,2)-1)
!       
!       end subroutine constant_boundaries_stag
    
    subroutine periodic_boundaries_3d(m, apply_to_x, apply_to_y)
        double precision, dimension(:,:,:) :: m
        integer :: maxx, maxy, maxz , k
        logical :: apply_to_x, apply_to_y
                
        maxx = size(m,2)
        maxy = size(m,1)
        maxz = size(m,3)
                
    
        do k = 1, maxz
            call periodic_boundaries(m(:,:,k),apply_to_x, apply_to_y)
        end do    
    end subroutine periodic_boundaries_3d
    
    subroutine periodic_boundaries_3d_stag(m,apply_to_x,apply_to_y)
        double precision, dimension(:,:,:) :: m
        integer :: k
        logical :: apply_to_x, apply_to_y
    
        do k=1, size(m, 3)
            call periodic_boundaries_stag(m(:,:,k),apply_to_x, apply_to_y)
        end do
    end subroutine periodic_boundaries_3d_stag

!       subroutine constant_boundaries_3d_stag(m, minx,miny,maxx,maxy)
!               double precision, dimension(:,:,:) :: m
!               double precision :: minx,miny,maxx,maxy
!               integer :: k
!               
!               
!               do k = 1, size(m,3)
!                       call constant_boundaries_stag(m(:,:,k),minx,miny,maxx,maxy)
!               end do
!       
!       end subroutine

!
!
!------------------------------------------------------
!   Velocity estimation according to higher order model
!------------------------------------------------------
!
    SUBROUTINE veloc2(mu,uvel,vvel,arrh,dzdx,dzdy,h,ax,ay,&
                zeta,bx,by,cxy,beta,&
                dhbdx,dhbdy,FLOWN,ZIP,VEL2ERR,&
                MANIFOLD,TOLER,PERIODIC_X, PERIODIC_Y, PLASTIC, delta_x, delta_y)

        INTEGER MAXY,MAXX,NZETA,MANIFOLD,PLASTIC

        double precision, dimension(:,:,:) :: mu
        double precision, dimension(:,:,:) :: uvel
        double precision, dimension(:,:,:) :: vvel
        double precision, dimension(:,:,:) :: arrh
        double precision, dimension(:,:) :: dzdx
        double precision, dimension(:,:) :: dzdy
        double precision, dimension(:,:) :: h
        double precision, dimension(:,:,:) :: ax
        double precision, dimension(:,:,:) :: ay
        double precision, dimension(:,:,:) :: bx
        double precision, dimension(:,:,:) :: by
        double precision, dimension(:,:,:) :: cxy
        double precision, dimension(:) :: zeta
        double precision, dimension(:,:) :: dhbdx
        double precision, dimension(:,:) :: dhbdy
        double precision, dimension(:,:) :: beta

        double precision :: FLOWN,ZIP,VEL2ERR,TOLER,delta_x, delta_y

        INTEGER i,j,k,l,lacc,m,n,maxiter,iter,DU1,DU2,DV1,DV2,&
            sparuv
        double precision :: error,tot,alfa,norm1,norm2,norm3,norm4,&
            norm5,teta

        logical :: periodic_x, periodic_y
        PARAMETER (DU1=1,DU2=2,DV1=3,DV2=4)

        double precision, dimension(:,:), allocatable :: em

        !Velocity estimates computed for the *current* iteration.  uvel and
        !vvel, comparitively, hold the velocity estimates for the *last*
        !iteration.
        double precision, dimension(:,:,:), allocatable :: ustar, vstar
        double precision, dimension(:,:), allocatable :: tau !Basal traction, to be computed from the provided beta parameter
       
        !Get the sizes from the mu field.  Calling code must make sure that
        !these all agree.
        maxy = size(mu,1)
        maxx = size(mu,2)
        nzeta = size(mu,3)

        allocate(em(4, maxx*maxy*nzeta))
        allocate(ustar(maxy, maxx, nzeta))
        allocate(vstar(maxy, maxx, nzeta))
        allocate(tau(maxy, maxx))
        maxiter=100
        error=VEL2ERR
        m=1
        lacc=0
        em = 0

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
        write(*,*)zeta
#endif


        !Copy the velocity estimate from the previous iteration as the current velocity
        ustar=uvel
        vstar=vvel
      
        !Deal with periodic boundary conditions
        call periodic_boundaries_3d_stag(uvel,periodic_x,periodic_y)
        call periodic_boundaries_3d_stag(vvel,periodic_x,periodic_y)
            !Compute basal traction

      
        nonlinear_iteration: do l=1,maxiter !Loop until we have reached the number of iterations allowed
            lacc=lacc+1
            
            !Every 10 iterations, we relax the error requirement somewhat (1/2
            !an order of magnitude) under the assumption that we won't converge 
            !given the current error tolerance.
            !(TODO: This leads to exponential increase in the tolerance - can we
            !get by without it?)
            if (l == m*10) then
                error = error*5.
                m = m+1
            endif
            
            !Compute the basal traction.  This is done by either passing through
            !the beta parameter, or by using Shoof's plastic bed model (Schoof
            !2004).
            if (PLASTIC == 0) then
                tau = beta
            else
                call plastic_bed(tau, beta, uvel, vvel)
            end if

       
            !Compute viscosity
            call muterm(mu,uvel,vvel,arrh,h,ax,ay,delta_x,delta_y,zeta,FLOWN,ZIP)

            !Sparse matrix routine for determining velocities.  The new
            !velocities will get spit into ustar and vstar, while uvel and vvel
            !will still hold the old velocities.
            iter=sparuv(mu,dzdx,dzdy,ax,ay,bx,by,cxy,h,&
                uvel,vvel,ustar,vstar,tau,dhbdx, dhbdy,MAXX*MAXY*NZETA,MAXY,&
                MAXX,NZETA,TOLER, delta_x, delta_y, zeta)
          
            norm1=0.
            norm2=0.
            norm3=0.
            norm4=0.
            norm5=0.
            n=0
        
            !em is a vector of distances between successive elements of uvel and vvel
            !This is the correction vector referred to in section 5.1.2
            !As we compute em, we also compute various L2 (Euclidean) norms for use
            !in later computations
            do i = 1, MAXY
                do j = 1, MAXX
                    do k = 1, NZETA
                        n = n + 1
              
                        em(DU2,n)=ustar(i,j,k)-uvel(i,j,k)
                        em(DV2,n)=vstar(i,j,k)-vvel(i,j,k)
              
                        !\| c^{l-1}\|
                        norm1 = norm1 + ((em(DU2,n)-em(DU1,n))**2)+((em(DV2,n)-em(DV1,n))**2)
                        !\|c^l - c^{l-1}\|
                        norm2 = norm2 + (em(DU1,n)**2)+(em(DV1,n)**2)
                        !\|c^l\|
                        norm3 = norm3 + (em(DU2,n)**2)+(em(DV2,n)**2)
                        !\|c^l \cdot c^{l-1}\|
                        norm4 = norm4 + em(DU1,n)*em(DU2,n)+em(DV1,n)*em(DV2,n)
              
                        norm5 = norm5 + (ustar(i,j,k)**2)+(vstar(i,j,k)**2)
                end do
            end do
        end do
   
        !Compute the angle between successive correction vectors
        if ((abs(norm2) < SMALL) .or. (abs(norm3) < SMALL)) then
            teta=PI/2.
        else
            teta=acos(norm4/sqrt(norm2*norm3))
        endif
        
        if ( (teta <= (5.*PI/6.) ) .and. (MANIFOLD == 1) ) then
            !We've requested unstable manifold correction, and the angle is
            !small (less than 5pi/6, a value identified by Hindmarsh and Payne
            !to work well).   If this is the case, we compute and apply
            !a correction vector.
            
            !Compute the error between the last two *correction vectors* (not
            !the last two iteration values!)  (See (51) in Pattyn's paper)
            if (abs(norm2) > 0.) then !We're just avoiding a divide by 0 here
                alfa=sqrt(norm1/norm2)
            else
                alfa=1.
            endif
            
            if (alfa < 1.e-6) then
                lacc=2*maxiter !This will end the iteration below - it occurs if
                               !correction vector didn't change much.
            else
                !Update the previous guess of the velocity with the correction
                !vector.  This throws out the current iteration's computed
                !velocity, and instead uses the computed correction vector.
                n=0
                do i=1,MAXY
                    do j=1,MAXX
                        do k=1,NZETA
                            n=n+1
                            uvel(i,j,k)=uvel(i,j,k)+em(DU2,n)/alfa
                            em(DU1,n)=em(DU2,n)
                            vvel(i,j,k)=vvel(i,j,k)+em(DV2,n)/alfa
                            em(DV1,n)=em(DV2,n)
                        end do
                    end do
                end do
            endif
          
        else
            n=0
          
            !Copy this iteration's new values to the old values
            !for the next iteration - because the angle between correction
            !vectors is large we do not want to apply a correction.
            uvel=ustar
            vvel=vstar
            em(DU1,:) = em(DU2,:)
            em(DV1,:) = em(DV2,:)
          
        endif
        
        tot=sqrt(norm3/norm5)
        write(*,*) l,iter,tot,(teta*180./PI)
        
        if (tot.lt.error) lacc=2*maxiter
       
        if (lacc.gt.maxiter) exit nonlinear_iteration
        
      end do nonlinear_iteration

      call write_xls_3d("uvel.txt",uvel)
      call write_xls_3d("vvel.txt",vvel)
      call write_xls("uvel_surf.txt",uvel(:,:,1))
      call write_xls("vvel_surf.txt",vvel(:,:,1))

      deallocate(em)
      deallocate(ustar)
      deallocate(vstar)
      return
      END subroutine
      
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
    subroutine plastic_bed(tau, tau0, uvel, vvel)
        double precision, dimension(:,:), intent(out)  :: tau
        double precision, dimension(:,:), intent(in)   :: tau0
        double precision, dimension(:,:,:), intent(in) :: uvel
        double precision, dimension(:,:,:), intent(in) :: vvel
        
        integer :: maxy, maxx, nzeta, i, j, k

        maxy = size(uvel,1)
        maxx = size(uvel,2)
        nzeta = size(uvel,3)
        do i = 1,maxy
            do j = 1,maxx
                tau(i,j) = tau0(i,j) / (sqrt(uvel(i,j,nzeta)**2 + vvel(i,j,nzeta)**2) + ZIP)
            end do
        end do
    end subroutine plastic_bed

!
!
!------------------------------------------------------
!   nonlinear viscosity term
!------------------------------------------------------
!
      subroutine muterm(mu,uvel,vvel,arrh,h,ax,ay,dx,dy,dz,FLOWN,ZIP)
!
        double precision,                   intent(in)  :: FLOWN,ZIP
        double precision, dimension(:,:,:), intent(out) :: mu
        double precision, dimension(:,:,:), intent(in)  :: uvel
        double precision, dimension(:,:,:), intent(in)  :: vvel
        double precision, dimension(:,:,:), intent(in)  :: arrh
        double precision, dimension(:,:),   intent(in)  :: h
        double precision, dimension(:,:,:), intent(in)  :: ax
        double precision, dimension(:,:,:), intent(in)  :: ay
        double precision,                   intent(in)  :: dx
        double precision,                   intent(in)  :: dy
        double precision, dimension(:),     intent(in)  :: dz
        
        !Derivatives of velocity.  These will be allocated and torn down during this function.
        !TODO: Make this more efficient, so that we are only allocating these once per run!!
        double precision, dimension(:,:,:), allocatable :: dudx
        double precision, dimension(:,:,:), allocatable :: dudy
        double precision, dimension(:,:,:), allocatable :: dudz
        double precision, dimension(:,:,:), allocatable :: dvdx
        double precision, dimension(:,:,:), allocatable :: dvdy
        double precision, dimension(:,:,:), allocatable :: dvdz
!
        INTEGER :: i,j,k, MAXX, MAXY, NZETA
        double precision :: macht
        double precision :: exx,eyy,exy,exz,eyz,eeff
!
        MAXY = size(uvel, 1)
        MAXX = size(uvel, 2)
        NZETA = size(dz)

        macht=(1.-FLOWN)/(2.*FLOWN)
      
        !Allocate temporary derivative data
        allocate(dudx(MAXY, MAXX, NZETA))
        allocate(dudy(MAXY, MAXX, NZETA))
        allocate(dudz(MAXY, MAXX, NZETA))
        allocate(dvdx(MAXY, MAXX, NZETA))
        allocate(dvdy(MAXY, MAXX, NZETA))
        allocate(dvdz(MAXY, MAXX, NZETA))

        !TODO: Implement option for upstream differencing
        call df_field_3d(uvel, dy, dx, dz, dudy, dudx, dudz)
        call df_field_3d(vvel, dy, dx, dz, dvdy, dvdx, dvdz)
      
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

        !Free temporary derivative data
        deallocate(dudx)
        deallocate(dudy)
        deallocate(dudz)
        deallocate(dvdx)
        deallocate(dvdy)
        deallocate(dvdz)
        return
    end subroutine

!
!------------------------------------------------------
!   determination of velocity u and v
!   with higher order model using sparse
!   matrices
!------------------------------------------------------
!
    function sparuv(mu,dzdx,dzdy,ax,ay,bx,by,cxy,h,uvel,vvel,ustar,vstar,beta,dhbdx,dhbdy,&
                    IJKTOT,MAXY,MAXX,NZETA,TOLER,GRIDX,GRIDY,zeta)
        INTEGER IJKTOT,MAXY,MAXX,NZETA
        double precision, dimension(:,:,:) :: mu
        double precision, dimension(:,:) :: dzdx
        double precision, dimension(:,:) :: dzdy
        double precision, dimension(:,:,:) :: ax
        double precision, dimension(:,:,:) :: ay
        double precision, dimension(:,:,:) :: bx
        double precision, dimension(:,:,:) :: by
        double precision, dimension(:,:,:) :: cxy
        double precision, dimension(:,:) :: h
        double precision, dimension(:,:,:) :: uvel
        double precision, dimension(:,:,:) :: vvel
        double precision, dimension(:,:,:) :: ustar
        double precision, dimension(:,:,:) :: vstar
        double precision, dimension(:,:) :: dhbdx
        double precision, dimension(:,:) :: dhbdy
        double precision, dimension(:,:) :: beta
        double precision, dimension(:) :: zeta
        double precision :: toler
        double precision :: gridx
        double precision :: gridy
        double precision :: rhs
          
        INTEGER i,j,k,l,m,sparuv,iter, ierr
        double precision :: d(IJKTOT),x(IJKTOT),coef(22),err
      
        type(sparse_matrix_type) :: matrix
        type(sparse_solver_workspace) :: workspace
        type(sparse_solver_options) :: options

       
        !Set up sparse matrix options
        call sparse_solver_default_options(options)
        options%tolerance=TOLER

        !Create the sparse matrix
        call new_sparse_matrix(ijktot, ijktot*22, matrix)
        call sparse_allocate_workspace(matrix, options, workspace)

!
!-------  velocity u
!
        !Initialize sparse matrix & vectors
        d=0
        x=0

        do i=1,MAXY
            do j=1,MAXX
                do k=1,NZETA 
                    coef = 0
                    if (h(i,j).lt.SMALL) then
                        !No ice - "pass through"
                        coef(I_J_K)=1.
                        !Normally, we use uvel(i,j,k) as our initial guess.
                        !However, in this case, we know the velocity
                        !should be 0, so we'll replace our uvel with that
                        uvel(i,j,k) = 0
                    else if ((i.eq.1).or.(i.eq.MAXY).or.(j.eq.1).or.(j.eq.MAXX)) then
                        !Boundary condition
                        coef(I_J_K)=1.
                        rhs=uvel(i,j,k)
                    else

                        call sparse_setup("u", i, j, k, mu, dzdx, dzdy, ax, ay, bx, by, cxy, &
                                  h, gridx, gridy, zeta, uvel, vvel, dhbdx, dhbdy, beta, &
                                  maxx, maxy, Nzeta, coef, rhs)
                    endif
                    
                    d(csp(I_J_K,i,j,k,MAXX,NZETA))=rhs
                    !Preliminary benchmarks indicate the we actually reach
                    !convergance faster if we use a 0 initial guess rather than
                    !the previous velocity.  More testing is needed.
                    !To use the current velocity as the guess, uncomment this
                    !line here and and in the Y velocity section of this
                    !function.
                    !x(csp(I_J_K,i,j,k,MAXX,NZETA))=uvel(i,j,k) !Use current velocity as guess of solution
                        
                    do m=1,21
                      call sparse_insert_val(matrix,csp(11,i,j,k,MAXX,NZETA),csp(m,i,j,k,MAXX,NZETA),coef(m)) 
                      !write(*,*)i,j,k,m,MAXX,NZETA,csp(m,i,j,k,MAXX,NZETA)
                    end do

                end do
             end do
        end do

        call sparse_solver_preprocess(matrix, options, workspace)
        
        ierr = sparse_solve(matrix, d, x, options, workspace,  err, iter, verbose=.false.)
        if (ierr /= 0) then
            write(*,*)"Error in SLAP dslucs:",ierr
            call sparse_interpret_error(ierr)
            write(*,*)"This is considered fatal."
            stop
        end if
        
        l=0
        do i=1,MAXY
            do j=1,MAXX
                do k=1,NZETA
                    l=l+1
                    ustar(i,j,k)=x(l)
                end do
            end do
        end do
!
!-------  velocity v
!
        !Initialize sparse matrix & vectors
        d=0
        x=0
        call sparse_clear(matrix)
     
        do i=1,MAXY
            do j=1,MAXX
                do k=1,NZETA
                    coef = 0
                    if (h(i,j).lt.SMALL) then
                        coef(I_J_K)=1.
                        uvel(i,j,k) = 0
                    else if ((i.eq.1).or.(i.eq.MAXY).or.(j.eq.1).or.(j.eq.MAXX)) then
                        coef(I_J_K)=1.
                        rhs=vvel(i,j,k)
                    else
                        call sparse_setup("v", i, j, k, mu, dzdx, dzdy, ax, ay, bx, by, cxy, &
                                  h, gridx, gridy, zeta, uvel, vvel, dhbdx, dhbdy, beta, &
                                  maxx, maxy, Nzeta, coef, rhs)
                     endif
                    coef(22) = rhs
                    !11 corresponds to position (i,j,k)
                    d(csp(11,i,j,k,MAXX,NZETA))=rhs
                    !x(csp(11,i,j,k,MAXX,NZETA))=vvel(i,j,k) !Use current velocity as guess of solution
            
                    do m=1,21
                      call sparse_insert_val(matrix,csp(11,i,j,k,MAXX,NZETA),csp(m,i,j,k,MAXX,NZETA),coef(m))
                    end do
                end do
            end do
        end do

#if 1
        call write_sparse_system("uvel_triad.txt",matrix%row, matrix%col, matrix%val, matrix%n)
        open(37,file="uvel_rhs.txt")
        write(37,*) d
        close(37)
        !stop
#endif

        ierr = sparse_solve(matrix, d, x, options, workspace,  err, iter, verbose=.false.)
        if (ierr /= 0) then
            write(*,*)"Error in SLAP dslucs:",ierr
            call sparse_interpret_error(ierr)
            write(*,*)"This is considered fatal."
            stop
        end if

        l=0
        do i=1,MAXY
            do j=1,MAXX
                do k=1,NZETA
                    l=l+1
                    vstar(i,j,k)=x(l)
                end do
            end do
        end do
        sparuv=iter
        
        call sparse_destroy_workspace(matrix, options, workspace)
        call del_sparse_matrix(matrix) 
    end function sparuv
!
!
!------------------------------------------------------
!   lookup coordinates in sparse matrix for 2D matrix
!
!   Reference:
!     1   ->   i - 1, j
!     2   ->   i, j - 1
!     3   ->   i, j
!     4   ->   i, j + 1
!     5   ->   i + 1, j
!------------------------------------------------------
!
      FUNCTION cspice(pos,i,j,MAXX)
!
        INTEGER pos,i,j,cspice,MAXX
!
      if (pos.eq.1) cspice=(i-2)*MAXX+j
      if (pos.eq.2) cspice=(i-1)*MAXX+j-1
      if (pos.eq.3) cspice=(i-1)*MAXX+j
      if (pos.eq.4) cspice=(i-1)*MAXX+j+1
      if (pos.eq.5) cspice=i*MAXX+j
      return
      END function
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
      FUNCTION csp(pos,i,j,k,MAXX,NZETA)
!
        INTEGER pos,i,j,k,csp,MAXX,NZETA
!
      if (pos.eq.1) csp=(i-2)*MAXX*NZETA+(j-2)*NZETA+k
!
      if (pos.eq.2) csp=(i-2)*MAXX*NZETA+(j-1)*NZETA+k-1
      if (pos.eq.3) csp=(i-2)*MAXX*NZETA+(j-1)*NZETA+k
      if (pos.eq.4) csp=(i-2)*MAXX*NZETA+(j-1)*NZETA+k+1
!
      if (pos.eq.5) csp=(i-2)*MAXX*NZETA+j*NZETA+k
!
      if (pos.eq.6) csp=(i-1)*MAXX*NZETA+(j-2)*NZETA+k-1
      if (pos.eq.7) csp=(i-1)*MAXX*NZETA+(j-2)*NZETA+k
      if (pos.eq.8) csp=(i-1)*MAXX*NZETA+(j-2)*NZETA+k+1
!
      if (pos.eq.9) csp=(i-1)*MAXX*NZETA+(j-1)*NZETA+k-2
      if (pos.eq.10) csp=(i-1)*MAXX*NZETA+(j-1)*NZETA+k-1
      if (pos.eq.11) csp=(i-1)*MAXX*NZETA+(j-1)*NZETA+k
      if (pos.eq.12) csp=(i-1)*MAXX*NZETA+(j-1)*NZETA+k+1
      if (pos.eq.13) csp=(i-1)*MAXX*NZETA+(j-1)*NZETA+k+2
!
      if (pos.eq.14) csp=(i-1)*MAXX*NZETA+j*NZETA+k-1
      if (pos.eq.15) csp=(i-1)*MAXX*NZETA+j*NZETA+k
      if (pos.eq.16) csp=(i-1)*MAXX*NZETA+j*NZETA+k+1
!
      if (pos.eq.17) csp=i*MAXX*NZETA+(j-2)*NZETA+k
!
      if (pos.eq.18) csp=i*MAXX*NZETA+(j-1)*NZETA+k-1
      if (pos.eq.19) csp=i*MAXX*NZETA+(j-1)*NZETA+k
      if (pos.eq.20) csp=i*MAXX*NZETA+(j-1)*NZETA+k+1
!
      if (pos.eq.21) csp=i*MAXX*NZETA+j*NZETA+k
      return
      END function
      
      
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
    subroutine sparse_setup(component, i,j,k,mu,dzdx,dzdy,ax,ay,bx,by,cxy,&
        h,dx,dy,dz,uvel,vvel,dhbdx,dhbdy,beta,MAXX,MAXY,Ndz,&
        coef, rhs)
!
        integer :: i,j,k, MAXY, MAXX, Ndz

        !Output array
        double precision, dimension(22), intent(out) :: coef

        !Output RHS value
        double precision, intent(out) :: rhs

        !Viscosity
        double precision, dimension(:,:,:) :: mu

        !Surface Gradients
        double precision, dimension(:,:), target :: dzdx, dzdy

        !Rescaling Factors
        double precision, dimension(:,:,:), target :: ax, ay, bx, by, cxy
        
        !Ice Thickness
        double precision, dimension(:,:) :: h

        !Grid Spacing (Z is an irregular grid)
        double precision :: dx, dy
        double precision, dimension(:) :: dz

        !Current velocities
        double precision, dimension(:,:,:), target :: uvel, vvel

        !Bedrock Gradients
        double precision, dimension(:,:), target :: dhbdx, dhbdy

        !Basal Traction
        double precision, dimension(:,:) :: beta

        character(1) :: component
!
        !Temporary calculations done in the ABSOLUTE (using u/x, v/y) coordinate
        !system
        double precision, target :: dmudx, dmudy, dmudz, dmudx2, dmudy2
        double precision, target :: dudx, dudy, dudz, dudx2, dudy2, dudz2, dudxz, dudyz, dudxy
        double precision, target :: dvdx, dvdy, dvdz, dvdx2, dvdy2, dvdz2, dvdxz, dvdyz, dvdxy

        !Temporary calculations done in the RELATIVE (using para, perp)
        !coordinate system
        double precision, pointer :: dpara_dpara, dpara_dperp, dpara_dx, dpara_dy, dpara_dz
        double precision, pointer :: dpara_dpara2, dpara_dparaz, dpara_dz2, dpara_dperp2, dpara_dperpz
        double precision, pointer :: dpara_dyz, dpara_dxz, dpara_dx2, dpara_dy2
        double precision, pointer :: dperp_dpara, dperp_dperp, dperp_dx, dperp_dy, dperp_dz
        double precision, pointer :: dperp_dxy, dperp_dxz, dperp_dyz, dperp_dz2
        double precision, pointer :: dmu_dpara2, dmu_dperp2

        !Fields done in the RELATIVE coordinate system
        double precision, dimension(:,:,:), pointer :: a_para, a_perp, b_para, b_perp
        double precision, dimension(:,:,:), pointer :: vel_perp
        double precision, dimension(:,:), pointer :: dz_dpara, dz_dperp, dhb_dpara, dhb_dperp

        !Cached difference calculations for the nonstaggered Z grid
        double precision dz_down1, dz_down2, dz_down3, dz_up1, dz_up2, dz_up3
        double precision dz_cen1, dz_cen2, dz_cen3, dz_sec1, dz_sec2, dz_sec3
        
        !call write_xls_3d("uvel.txt",uvel)
        !call write_xls_3d("vvel.txt",vvel)

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
        else
            write(*,*)"FATAL ERROR: sparse_setup called with invalid component"
        end if

         
     
      
        !Set up a system of pointers

        if (k.eq.1) then !Upper boundary condition (stress-free surface)
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

            coef(3) = -.5*dpara_dy/dy
            coef(7) = -.5*dpara_dx/dx
            coef(11) = dpara_dz*dz_down1
            coef(12) = dpara_dz*dz_down2
            coef(13) = dpara_dz*dz_down3
            coef(15) = dpara_dx*.5/dx
            coef(19) = dpara_dy*.5/dy
            !Note the transposition below (dfdx => dfdy, vice versa)
            rhs = -dperp_dx * dfdy_3d(vel_perp,i,j,k,dx) &
                -  dperp_dy * dfdx_3d(vel_perp,i,j,k,dy) &
                -  dperp_dz * dfdz_3d_downwind_irregular(vel_perp,i,j,k,dz)

        else if (k.eq.Ndz) then !Lower boundary condition (Basal sliding)
            !Finite difference coefficients for an irregular Z grid, upwinded
            dz_up1=(dz(k)-dz(k-1))/(dz(k-1)-dz(k-2))/(dz(k)-dz(k-2))
            dz_up2=(dz(k-2)-dz(k))/(dz(k)-dz(k-1))/(dz(k-1)-dz(k-2))
            dz_up3=(2.*dz(k)-dz(k-1)-dz(k-2))/(dz(k)-dz(k-1))/(dz(k)-dz(k-2))

            if (beta(i,j).le.NOTDEF) then !No sliding at the base
                coef(I_J_K)=1.
                rhs=0.
            else !Compute sliding 
                dpara_dpara = 4.*dhb_dpara(i,j)
                dpara_dz = 4.*a_para(i,j,k)*dhb_dpara(i,j) + a_perp(i,j,k)*dhb_dperp(i,j) + 1./h(i,j)
                dpara_dperp = dhb_dperp(i,j)

                dperp_dpara = dhb_dperp(i,j)
                dperp_dperp = 2.*dhb_dpara(i,j)
                dperp_dz = 2.*a_perp(i,j,k)*dhb_dpara(i,j) + a_para(i,j,k)*dhb_dperp(i,j)

                !If we have a non-zero friction, then apply that friction
                !parameter
                !Otherwise, we have a spot of zero traction or an ice shelf\
                !and we leave everything unscaled
                if (beta(i,j) >= SMALL) then
                   dpara_dpara = dpara_dpara*mu(i,j,k)
                   dpara_dperp = dpara_dperp*mu(i,j,k)
                   dpara_dz = dpara_dz*mu(i,j,k)
                   dperp_dpara = dperp_dpara*mu(i,j,k)
                   dperp_dperp = dperp_dperp*mu(i,j,k)
                   dperp_dz = dperp_dz*mu(i,j,k)
                end if

                coef(3) = -.5*dpara_dy/dy
                coef(7) = -.5*dpara_dx/dx
                coef(9) = dpara_dz*dz_up1
                coef(10)= dpara_dz*dz_up2
                coef(11)= dpara_dz*dz_up3
                coef(15)= .5*dpara_dx/dx
                coef(19)= .5*dpara_dy/dy
                !Transposition of derivatives below
                rhs = -dperp_dx * dfdy_3d(vel_perp, i, j, k, dx) &
                    -  dperp_dy * dfdx_3d(vel_perp, i, j, k, dy) &
                    -  dperp_dz * dfdz_3d_upwind_irregular(vel_perp, i, j, k, dz)

                !Adjust center of stencil for the basal traction if it is
                !relevant
                if (beta(i,j) >= SMALL) then
                    coef(11) = coef(11) + beta(i,j)
                end if
            
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
            dperp_dz  = 2 * a_perp(i,j,k) * dmu_dpara2 + 3*mu(i,j,k) * cxy(i,j,k) + a_para(i,j,k) * dmu_dperp2
            dperp_dxy = 3 * mu(i,j,k)
            dperp_dxz = 3 * mu(i,j,k) * ay(i,j,k)
            dperp_dyz = 3 * mu(i,j,k) * ax(i,j,k)
            dperp_dz2 = 3 * mu(i,j,k) * ax(i,j,k)*ay(i,j,k)

            coef(2)  = dpara_dyz*dz_cen1*(-.5/dy)
            coef(3)  = (dpara_dyz*dz_cen2 + dpara_dy)*(-.5/dy) + dpara_dy2/(dy**2)
            coef(4)  = dpara_dyz*dz_cen3*(-.5/dy)
          
            coef(6)  = dpara_dxz*dz_cen1*(-.5/dx)
            coef(7)  = (dpara_dx + dpara_dxz*dz_cen2)*(-.5/dx) + dpara_dx2 / dx**2
            coef(8)  = dpara_dxz*dz_cen3*(-.5/dx)
          
            coef(10) = dpara_dz*dz_cen1 + dpara_dz2*dz_sec1
            coef(11) = dpara_dz*dz_cen2 + dpara_dz2*dz_sec2 + dpara_dx2*(-2/dx**2) + dpara_dy2*(-2/dy**2)
            coef(12) = dpara_dz*dz_cen3 + dpara_dz2*dz_sec3
          
            coef(14) = dpara_dxz * dz_cen1 * (.5/dx)
            coef(15) = (dpara_dx + dpara_dxz*dz_cen2) * (.5/dx) + dpara_dx2 / dx**2
            coef(16) = dpara_dxz * dz_cen3 * (.5/dx)
          
            coef(18) = dpara_dyz * dz_cen1 * (.5/dy)
            coef(19) = (dpara_dyz*dz_cen2 + dpara_dy) * (.5/dy) + dpara_dy2 / dy**2
            coef(20) = dpara_dyz * dz_cen3 * (.5/dy)
          
            !Transposition of derivatives below
            rhs = RHOI * GRAV * dz_dpara(i,j) &
                - dperp_dx  * dfdy_3d(vel_perp,i,j,k,dx)&
                - dperp_dy  * dfdx_3d(vel_perp,i,j,k,dy)&
                - dperp_dz  * dfdz_3d_irregular(vel_perp,i,j,k,dz)&
                - dperp_dxy * d2fdxy_3d(vel_perp,i,j,k,dx,dy)&
                - dperp_dxz * d2fdyz_3d(vel_perp,i,j,k,dx,dz)&
                - dperp_dyz * d2fdxz_3d(vel_perp,i,j,k,dy,dz)&
                - dperp_dz2 * d2fdz2_3d_irregular(vel_perp,i,j,k,dz)
            !write(*,*)i,j,k,dz_dpara(i,j),& !4
            !      dperp_dx,dfdy_3d(vel_perp,i,j,k,dx),& !6
            !      dperp_dy,dfdx_3d(vel_perp,i,j,k,dy),& !8
            !      dperp_dz,dfdz_3d_irregular(vel_perp,i,j,k,dz),& !10
            !      dperp_dxy,d2fdxy_3d(vel_perp,i,j,k,dx,dy),&
            !      dperp_dxz,d2fdyz_3d(vel_perp,i,j,k,dx,dz),&
            !      dperp_dyz,d2fdxz_3d(vel_perp,i,j,k,dy,dz),&
            !      dperp_dz2,d2fdz2_3d_irregular(vel_perp,i,j,k,dz)
        end if
    end subroutine sparse_setup

!
!
!------------------------------------------------------
!   Diffusive equation for higher-order model
!   (EISMINT TYPE II)
!------------------------------------------------------
!
      SUBROUTINE diffus2(h,u,v,uvel,vvel,dzdx,dzdy,diff1,&
        diffus,arrh,zeta,MAXY,MAXX,NZETA,SMALL,FLOWN,&
        PERIODIC_X, PERIODIC_Y)
!
        double precision :: h(MAXY,MAXX),u(MAXY,MAXX),v(MAXY,MAXX),&
          uvel(MAXY,MAXX,NZETA),&
          vvel(MAXY,MAXX,NZETA),dzdx(MAXY,MAXX),dzdy(MAXY,MAXX),&
          diff1(MAXY,MAXX,NZETA),diffus(MAXY,MAXX),zeta(NZETA),&
          arrh(MAXY,MAXX,NZETA),SMALL,FLOWN
        INTEGER MAXY,MAXX,NZETA
        logical :: PERIODIC_X,PERIODIC_Y
      INTEGER i,j,k
      double precision :: grad,vel,d,z
!
      do 10 i=1,MAXY
        do 20 j=1,MAXX
          if (h(i,j).le.0.) then
            u(i,j)=0.
            v(i,j)=0.
            do 30 k=1,NZETA
              uvel(i,j,k)=0.
              vvel(i,j,k)=0.
   30       CONTINUE
          endif
   20   CONTINUE
   10 CONTINUE
!
      do 40 i=2,MAXY-1
        do 50 j=2,MAXX-1
          grad=dzdx(i,j)*dzdx(i,j)+dzdy(i,j)*dzdy(i,j)
          if (sqrt(grad).gt.SMALL) then
            do 60 k=1,NZETA
              vel=sqrt(uvel(i,j,k)*uvel(i,j,k)+vvel(i,j,k)*&
                vvel(i,j,k))
              diff1(i,j,k)=vel/sqrt(grad)
   60       CONTINUE
          else
            d=(RHOI*GRAV*h(i,j))**FLOWN
            diff1(i,j,k)=0.
            do 61 k=(NZETA-1),1,-1
              z=(0.5*(zeta(k+1)+zeta(k)))**FLOWN
              diff1(i,j,k)=d*grad*h(i,j)*(arrh(i,j,k+1)+arrh(i,j,k))*&
                z*(zeta(k+1)-zeta(k))+diff1(i,j,k+1)
   61       CONTINUE
          endif
          diffus(i,j)=0.
          do 62 k=2,NZETA
            diffus(i,j)=diffus(i,j)+0.5*(diff1(i,j,k)+&
              diff1(i,j,k-1))*(zeta(k)-zeta(k-1))
   62     CONTINUE
          diffus(i,j)=diffus(i,j)*h(i,j)
   50   CONTINUE
   40 CONTINUE
   
        call periodic_boundaries_stag(diffus,periodic_x,periodic_y)
        call periodic_boundaries_3d_stag(diff1,periodic_x,periodic_y)
        if (.not. PERIODIC_X) then
            diffus(:,1) = 0
            diffus(:,maxx) = 0
            diff1(:,1,:) = 0
            diff1(:,maxx,:) = 0
        end if

        if (.not. PERIODIC_Y) then
            diffus(1,:) = 0
            diffus(maxy,:) = 0
            diff1(1,:,:) = 0
            diff1(maxy,:,:) = 0
        end if
 
      return
      END subroutine
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
        double precision :: FLOWN,ZIP
        double precision, dimension(:,:,:), intent(inout) :: mu
        double precision, dimension(:,:,:), intent(in) :: uvel
        double precision, dimension(:,:,:), intent(in) :: vvel
        double precision, dimension(:,:,:), intent(in) :: arrh
        double precision, dimension(:,:),   intent(in) :: h
        double precision, dimension(:,:,:), intent(in) :: ax
        double precision, dimension(:,:,:), intent(in) :: ay
        double precision, dimension(:,:,:), intent(out) :: txz 
        double precision, dimension(:,:,:), intent(out) :: tyz
        double precision, dimension(:,:,:), intent(out) :: txx
        double precision, dimension(:,:,:), intent(out) :: tyy
        double precision, dimension(:,:,:), intent(out) :: txy
        double precision,                   intent(in) :: dx
        double precision,                   intent(in) :: dy
        double precision, dimension(:),     intent(in) :: dz
        
        logical :: PERIODIC_X, PERIODIC_Y

      INTEGER i,j,k
      double precision, dimension(:,:,:), allocatable :: dudz, dvdz, dudx, dvdx, dudy, dvdy

      double precision :: exx,eyy,exy,exz,eyz,eeff
      double precision :: macht
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
!
!
!------------------------------------------------------
!   U-velocity estimation from diffusivities (TYPE II)
!------------------------------------------------------
!
      SUBROUTINE ufield(h,hb,uflux,diff1,uvel,u,zeta,MAXY,MAXX,&
        NZETA,SMALL,GRIDX,PERIODIC_X, PERIODIC_Y,DIFTYPE)
!
        INTEGER MAXY,MAXX,NZETA,DIFTYPE
        double precision :: h(MAXY,MAXX),uflux(MAXY,MAXX,NZETA)
        double precision :: diff1(MAXY,MAXX,NZETA),u(MAXY,MAXX)
        double precision :: uvel(MAXY,MAXX,NZETA), zeta(NZETA)
        double precision :: hb(MAXY,MAXX),SMALL,GRIDX
!
        logical :: periodic_x, periodic_y
      INTEGER i,j,k
      double precision :: grad
!
      do 10 i=1,MAXY-1
        do 20 j=1,MAXX-1
          if ((h(i,j).gt.SMALL).or.(h(i,j+1).gt.SMALL)) then
            grad=(h(i,j+1)+hb(i,j+1)-h(i,j)-hb(i,j))/&
              (GRIDX*(h(i,j)+h(i,j+1)))
          else
            grad=0.
          endif
          do 30 k=1,NZETA
            if (DIFTYPE.eq.1) then
              uflux(i,j,k)=-(h(i,j)*diff1(i,j,k)+h(i,j+1)*&
                diff1(i,j+1,k))*grad
            else
              uflux(i,j,k)=(uvel(i,j,k)+uvel(i,j+1,k))/2.
            endif
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
      do 11 i=1,MAXY-1
        do 21 j=2,MAXX-1
          do 31 k=1,NZETA
            uvel(i,j,k)=0.5*(uflux(i,j,k)+uflux(i,j-1,k))
   31     CONTINUE
          u(i,j)=0.
          do 40 k=2,NZETA
            u(i,j)=u(i,j)+(uvel(i,j,k)+uvel(i,j,k-1))*&
              (zeta(k)-zeta(k-1))*0.5
   40     CONTINUE
   21   CONTINUE
   11 CONTINUE
   
        call periodic_boundaries_stag(u,PERIODIC_X,PERIODIC_Y)
        call periodic_boundaries_3d_stag(uvel,PERIODIC_X,PERIODIC_Y)
        if (.not. PERIODIC_X) then
            u(:,1) = 0
            u(:,maxx) = 0
            uvel(:,1,:) = 0
            uvel(:,maxx,:) = 0
        end if

        if (.not. PERIODIC_Y) then
            u(1,:) = 0
            u(maxy,:) = 0
            uvel(1,:,:) = 0
            uvel(maxy,:,:) = 0
        end if
      
      return
      END subroutine
!
!
!------------------------------------------------------
!   V-velocity estimation from diffusivities (TYPE II)
!------------------------------------------------------
!
      SUBROUTINE vfield(h,hb,vflux,diff1,vvel,v,zeta,MAXY,MAXX,&
        NZETA,SMALL,GRIDY,PERIODIC_X,PERIODIC_Y,DIFTYPE)
!
        INTEGER MAXY,MAXX,NZETA,DIFTYPE
        double precision :: h(MAXY,MAXX),vflux(MAXY,MAXX,NZETA)
        double precision :: diff1(MAXY,MAXX,NZETA),v(MAXY,MAXX)
        double precision :: vvel(MAXY,MAXX,NZETA), zeta(NZETA)
        double precision :: hb(MAXY,MAXX),SMALL,GRIDY
!
        logical :: PERIODIC_X, PERIODIC_Y
      INTEGER i,j,k
      double precision :: grad
!
      do 10 i=1,MAXY-1
        do 20 j=1,MAXX-1
          if ((h(i,j).gt.SMALL).or.(h(i+1,j).gt.SMALL)) then
            grad=(h(i+1,j)+hb(i+1,j)-h(i,j)-hb(i,j))/&
              (GRIDY*(h(i,j)+h(i+1,j)))
          else
            grad=0.
          endif
          do 30 k=1,NZETA
            if (DIFTYPE.eq.1) then
              vflux(i,j,k)=-(h(i,j)*diff1(i,j,k)+h(i+1,j)*&
                diff1(i+1,j,k))*grad
            else
              vflux(i,j,k)=(vvel(i,j,k)+vvel(i+1,j,k))/2.
            endif
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
      do 11 i=2,MAXY-1
        do 21 j=1,MAXX-1
          do 31 k=1,NZETA
            vvel(i,j,k)=0.5*(vflux(i,j,k)+vflux(i-1,j,k))
   31     CONTINUE
          v(i,j)=0.
          do 40 k=2,NZETA
            v(i,j)=v(i,j)+(vvel(i,j,k)+vvel(i,j,k-1))*&
              (zeta(k)-zeta(k-1))*0.5
   40     CONTINUE
   21   CONTINUE
   11 CONTINUE
   
        call periodic_boundaries_stag(v,PERIODIC_X,PERIODIC_Y)
        call periodic_boundaries_3d_stag(vvel,PERIODIC_X,PERIODIC_Y)
        if (.not. PERIODIC_X) then
            v(:,1) = 0
            v(:,maxx) = 0
            vvel(:,1,:) = 0
            vvel(:,maxx,:) = 0
        end if

        if (.not. PERIODIC_Y) then
            v(1,:) = 0
            v(maxy,:) = 0
            vvel(1,:,:) = 0
            vvel(maxy,:,:) = 0
        end if
    
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
    subroutine unstagger_field_2d_periodic(f_stag, f)
        real(dp), dimension(:,:), intent(in) :: f_stag
        real(dp), dimension(:,:), intent(out) :: f

        integer :: i,j, i1, i2, j1, j2, ni, nj
        
        ni = size(f, 1)
        nj = size(f, 2)

        do i = 1, size(f, 1)
            do j = 1, size(f, 2)
                i1 = i-1
                i2 = i
    
                if (i1 == 0) then
                    i1 = ni - 1
                end if
    
                if (i2 == ni) then
                    i2 = 1
                end if
    
                j1 = j-1
                j2 = j
    
                if (j1 == 0) then
                    j1 = nj - 1
                end if
    
                if (j2 == nj) then
                    j2 = 1
                end if
    
                f(i,j) = (f_stag(i1, j1) + f_stag(i2, j1) + f_stag(i1, j2) + f_stag(i2, j2))/2
            end do
        end do
    
    end subroutine unstagger_field_2d_periodic

    function isinf(n)
        double precision::n
        logical :: isinf
        
        isinf = isnan(0*n)
    end function isinf
    
    function isnan(n)
        double precision ::n
        logical :: isnan
        
        isnan = (n /= n)
    end function isnan
    
    function checkrow(row)
        double precision, dimension(:) :: row
        integer :: i
        logical :: checkrow
        checkrow=.false.
        do i = 1,size(row,1)
            if (row(i) /= 0) then
                checkrow=.true.
            end if
        end do
    end function checkrow

!
!
!---------------------------------------------------
!
!                 END OF SUBROUTINES
!
!---------------------------------------------------
end module ice3d_lib

