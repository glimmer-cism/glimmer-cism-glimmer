module plume_grid_set

  use plume_global

contains

subroutine grid_set()

  implicit none

! read grid cell data

! local variables

      integer i,k,imid,kmid

      real(kind=kdp) rmdx,adx,adxu
      real(kind=kdp) rmdxu,rmdy,rmdyv,ady,adyv

      imid=(m+1)/2
      kmid=(n+1)/2

! ----------------
! grid for scalars
! ----------------

! set east-west cell dimensions for uniform scalar grid
      do i = 1,m
        dx(i) = hx
      end do

! set north-south cell dimensions for uniform scalar grid
      do k = 1,n
        dy(k) = hy
      end do

! -----------------------------------------------
! grid for velocities and turbulence coefficients
! -----------------------------------------------

! calculate cell dimensions for velocity grid
      do k=2,n
        dyv(k) = 5.0d-1*(dy(k) + dy(k-1))
        rdy(k)=1.d0/dy(k)
        rdyv(k)=1.d0/dyv(k)
      end do
      dyv(1)=dyv(2)
      rdy(1)=rdy(2)
      rdyv(1)=rdyv(2)

      do i=2,m
        dxu(i) = 5.0d-1*(dx(i) + dx(i-1))
        rdx(i)= 1.d0/dx(i)
        rdxu(i)= 1.d0/dxu(i)
      end do
      dxu(1)=dxu(2)
      rdx(1)=rdx(2)
      rdxu(1)=rdxu(2)

! calculate coefficients for horizontal eddy viscosity
      rmdx = rdx(imid)
      rmdxu = rdxu(imid)
      adx = dt*ah*rmdx**2
      adxu = dt*ah*rmdxu**2
      do i = 1,m
        ahdx(i) = adx/rdx(i)
        ahdxu(i) = adxu/rdxu(i)
      end do

      rmdy = rdy(kmid)
      rmdyv = rdyv(kmid)
      ady = dt*ah*rmdy**2
      adyv = dt*ah*rmdyv**2
      do k = 1,n
        ahdy(k) = ady/rdy(k)
        ahdyv(k) = adyv/rdyv(k)
      end do

      return
end subroutine

end  module
