
 module params_mod     ! *sp* used by LANL incremental remapping code, 'remap_advection'

      use kinds_mod

      implicit none
      save
      public

     integer (i4), parameter :: &
      stdin     =  5,         &! reserved unit for standard input
      stdout    =  6,         &! reserved unit for standard output
      stderr    =  6           ! reserved unit for standard error

     integer (i4), parameter :: &
      POP_maxBlocksClinic = 1,  &! number of blocks
      POP_nt = 1                 ! number of tracers

     ! *sp* commented out 3 lines below. These now passed as args. to driver of 'remap_advection'.

!     integer (i4), parameter :: &
!      nghost = 2,               &! number of ghost cells 
!      nx_block = 33,            &! number of cells in x direction 
!      ny_block = 33              ! number of cells in y direction 

     ! *sp* commented out line below. Now defined in 'remap_advection' after nx_block, 
     ! ny_block defined (passed in to driver) 

!     integer (i4), parameter :: &
!      nx = nx_block, ny = ny_block,      

     integer (i4), parameter :: &
      nt = POP_nt, nblock = POP_maxBlocksClinic

     integer (i4), parameter :: &
      my_task = 0,              &!
      master_task = 0

 end module params_mod
