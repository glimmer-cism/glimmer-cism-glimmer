!
! $Id: smb_caller.f90,v 1.1.2.1 2005-03-20 20:11:28 gethinw Exp $
!

program SMBProg

  use smb_mecons
  use smb_force
  use smb_netcdf

  implicit none

  !! Variables passed in from calling context (e.g. GCM)
  integer :: nx                              ! Number of grid point in x-direction
  integer :: ny                              ! Number of grid point in y-direction
  integer :: dx                              ! grid resolution (m)
  integer :: mindt                           ! timestep in min
  real(rk),pointer,dimension(:,:) :: totacc  ! total water equivalent snowfall
  real(rk),pointer,dimension(:,:) :: totrn   ! total rainfall
  real(rk),pointer,dimension(:,:) :: totcon  ! if -ve, evap, if +ve, condensation 
  real(rk),pointer,dimension(:,:) :: totoff  ! total runoff
  real(rk),pointer,dimension(:,:) :: massbal ! totacc+totrn+totcon-totoff

  !! derived type holding vars for SMB routine
  type(smb_params) :: smb_p

  !! Intialisation
  nx    = 75
  ny    = 140
  dx    = 20000
  mindt = 30

  !! Allocations & init to zero
  allocate(totacc(nx,ny));  totacc =0.0
  allocate(totrn(nx,ny));   totrn  =0.0
  allocate(totcon(nx,ny));  totcon =0.0
  allocate(totoff(nx,ny));  totoff =0.0
  allocate(massbal(nx,ny)); massbal=0.0

  !! these two driven from file
  call SMBInitForceWrapper(smb_p,nx,ny,dx,mindt)
  call SMBForceLoopWrapper(smb_p, &
       totacc,totrn,totcon,totoff,massbal)

  !! write a NetCDF file with final massbal value
  !! NB massbal contains '6 hour' data in this case
  call SMBExampleNetCDF("foo.nc",nx,ny,dx,massbal)

  !! temp
  call SMBCleanupForce

  !! mainly free arrays allocated by init
  call SMBCleanupWrapper(smb_p)

  !! De-allocate
  if(associated(totacc))  deallocate(totacc)
  if(associated(totrn))   deallocate(totrn)
  if(associated(totcon))  deallocate(totcon)
  if(associated(totoff))  deallocate(totoff)
  if(associated(massbal)) deallocate(massbal)

end program SMBProg
