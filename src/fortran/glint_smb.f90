
module smb_mecons

  ! This module provides a dummy, hopefully warning-free interface
  ! in place of the Energy-balance mass-balance scheme. If either
  ! subroutine is called, a fatal error is flagged.

  use glimmer_global

  implicit none

  type smb_params
     integer       :: dummyint
     real(rk)      :: dummyreal
     character(40) :: dummypath
  end type smb_params

contains

  subroutine SMBInitWrapper(params,nx,ny,dxr,tstep,path)

    use glimmer_log

    type(smb_params) :: params
    integer :: nx,ny,dxr,tstep
    character(*) :: path

    ! Fatal error

    call write_log('Glimmer not compiled with EBMB mass-balance scheme',GM_FATAL,__FILE__,__LINE__)

    ! Need these lines to avoid warnings, though they are never executed

    params%dummyint=nx
    params%dummyint=ny
    params%dummyint=dxr
    params%dummyint=tstep
    params%dummypath=path

  end subroutine SMBInitWrapper

  !---------------------------------------------------------------------------------------------

  subroutine SMBStepWrapper(params,temp,thck,artm,prcp,U10m,V10m,humidity,SWdown,LWdown,Psurf)

    use glimmer_log

    type(smb_params)        :: params
    real(rk),dimension(:,:) :: temp,thck,artm,prcp,U10m,V10m,humidity,SWdown,LWdown,Psurf

    ! Fatal error

    call write_log('Glimmer not compiled with EBMB mass-balance scheme',GM_FATAL,__FILE__,__LINE__)

    ! Need this line to avoid warnings, though it is never executed

    params%dummyreal=sum(temp+thck+artm+prcp+U10m+V10m+humidity+SWdown+LWdown+Psurf)

  end subroutine SMBStepWrapper

end module smb_mecons
