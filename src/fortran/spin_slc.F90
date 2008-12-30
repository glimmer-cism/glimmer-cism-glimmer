! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  spin_slc.f90 - part of the GLIMMER ice model              + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

module spin_slc
  !*FD this module handles the sea level component

  use glimmer_ts
  use glimmer_global, only : fname_length

  type spin_slc_type
     !*FD parameters for the EIS sea level forcing
     character(len=fname_length) :: fname=''     !*FD name of file containing ts
     type(glimmer_tseries) :: slc_ts             !*FD time series 
  end type spin_slc_type
  
contains  
  subroutine spin_slc_config(config,slc)
    !*FD get SLC configuration from config file
    use glimmer_config
    implicit none
    type(spin_slc_type)           :: slc     !*FD slc data
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   
    ! local variables
    type(ConfigSection), pointer :: section

    call GetSection(config,section,'SPIN SLC')
    if (associated(section)) then
       call GetValue(section,'slc_file',slc%fname)
    end if
  end subroutine spin_slc_config

  subroutine spin_slc_printconfig(slc)
    !*FD print configuration to log
    use glimmer_log
    implicit none
    type(spin_slc_type)           :: slc     !*FD slc data
    ! local variables
    character(len=100) :: message

    call write_log('SPIN SLC')
    call write_log('-------')
    write(message,*) 'SLC file: ',trim(slc%fname)
    call write_log(message)
    call write_log('')
  end subroutine spin_slc_printconfig

  subroutine spin_init_slc(slc)
    !*FD initialise SLC forcing
    use glimmer_paramets, only: thk0
    implicit none
    type(spin_slc_type)           :: slc     !*FD slc data
    
    call glimmer_read_ts(slc%slc_ts,slc%fname)
    ! scale parameters
    slc%slc_ts%values = slc%slc_ts%values/thk0
  end subroutine spin_init_slc

  subroutine spin_eus(slc,model,time)
    !*FD calculate mass balance
    use glide_types
    use glimmer_global, only : rk
    implicit none
    type(spin_slc_type)        :: slc   !*FD slc data
    type(glide_global_type)   :: model !*FD model instance
    real(kind=rk), intent(in) :: time  !*FD current time

    call glimmer_ts_linear(slc%slc_ts,real(time),model%climate%eus)
  end subroutine spin_eus
end module spin_slc
