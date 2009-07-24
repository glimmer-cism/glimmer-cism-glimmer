! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  spin_slc.f90 - part of the GLIMMER ice model              + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! The Spin drivers are intended for model spin-up for Antarctica and Greenland.  The
! model_type variable differentiates between the two models, (0) for Antarctica,
! and (1) for Greenland.  These can be set inside the configuration file.  Also,
! included in the drivers are a earlier versions of the mass balance and
! temperature and are meant to be the same as the ones used in the EISMINT model
! intercomparison experiment.  More information can be found at
! http://homepages.vub.ac.be/~phuybrec/eismint.html.  These can be specified by
! setting the variable use_simple = 1 in the configuration file under temperature
! and massbalance separately.  These drivers also vary for Antarctica and
! Greenland.  The other drivers available in these drivers are based on the
! current papers of Philippe Huybrechts:  

! P. Huybrechts and J. de Wolde. The dynamic response of the greenland and
! antarctic ice sheets to multiple-century climatic warming. Journal of Climate,
! 12(8):2169-2188, 1999.

! I. Janssens and P. Huybrects. The treatment of meltwater retention in
! massbalance parameterisations of the greenland ice sheet. Annals of Glaciology,
! 31:133-140, 2000.

! P. Huybrechts. Sea-level changes at the LGM from ice-dynamic reconstructions
! of the Greenland and Antarctic ice sheets during the glacial cycles. Quater-
! nary Science Reviews, 21:203-231, 2002.

! P. Huybrechts, O. Rybak, F. Pattyn, U. Ruth, and D. Steinhage. Ice thinning,
! upstream advection, and non-climatic biases for the upper ice sheet. Climate
! of the Past, 3:577-589, 2007.

! The Spin series of drivers were developed by Brian Hand, at the University of
! Montana, 2008.
#ifdef HAVE_CONFIG_H
#include "config.inc"
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
