 module exit_mod     ! *sp* used by LANL incremental remapping code, 'remap_advection'

   use kinds_mod
   use constants

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: exit_POP

! !DEFINED PARAMETERS:

   integer (int_kind), parameter, public :: &
      sigExit  =  0,    &! signal for normal exit
      sigAbort = -1      ! signal for aborting (exit due to error)

 contains

!***********************************************************************
!BOP
! !IROUTINE: exit_POP
! !INTERFACE:

 subroutine exit_POP(exit_mode, exit_message)

! !DESCRIPTION:
!  This routine prints a message, exits any message environment
!  and cleans up before stopping

! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
     exit_mode    ! method for exiting (normal exit or abort)

   character (*), intent(in) :: &
     exit_message ! message to print before stopping

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: ierr  ! error flag

!-----------------------------------------------------------------------
!
!  print message - must use unit 6 in place of stdout to
!  prevent circular dependence with io module
!
!-----------------------------------------------------------------------

      write (6,delim_fmt)
      write (6,blank_fmt)
      write (6,*) exit_message
      write (6,blank_fmt)
      write (6,delim_fmt)

!-----------------------------------------------------------------------
!
!  now we can stop
!
!-----------------------------------------------------------------------

   stop

!-----------------------------------------------------------------------
!EOC

 end subroutine exit_POP

!***********************************************************************





 end module exit_mod
