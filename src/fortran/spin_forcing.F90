! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! + spin_forcing.f90 - part of the GLIMMER ice model          + 
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

module spin_forcing
  
  use spin_types
  
contains
  subroutine spin_initialise(climate,config,model,pdd_scheme)
    !*FD initialise climate
    use glide_types
    use glimmer_config
    use spin_io
    use glimmer_pdd
    implicit none
    type(spin_climate_type) :: climate     !*FD structure holding climate
    type(ConfigSection), pointer :: config !*FD structure holding sections of configuration file   
    type(glide_global_type) :: model       !*FD model instance
    type(glimmer_pdd_params) :: pdd_scheme !structure for positive degree day
    !create the fields needed to run the spin climate model 
    ! read config
    call spin_readconfig(climate,config)
    ! print config
    call spin_printconfig(climate)
    ! create spin variables
    call spin_io_createall(model)
    ! initialise subsystems
    call spin_init_mb(climate%mb,model)
    call spin_init_temp(climate%temp, model)
    call spin_init_slc(climate%slc)

    call glimmer_pdd_init(pdd_scheme, config)
    ! and read first time slice
    call spin_io_readall(climate,model)

  end subroutine spin_initialise

  subroutine spin_readconfig(climate,config)
    !*FD read spin configuration
    use glimmer_config
    implicit none
    type(spin_climate_type) :: climate    !*FD structure holding climate
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   
 
    call spin_mb_config(config,climate%mb)
    call spin_temp_config(config,climate%temp)
    call spin_slc_config(config,climate%slc)
  end subroutine spin_readconfig

  subroutine spin_printconfig(climate)
    !*FD print spin configuration
    use glimmer_log
    implicit none
    type(spin_climate_type) :: climate  !*FD structure holding climate
    
    call write_log_div
    call write_log('Spin Model')
    call spin_mb_printconfig(climate%mb)
    call spin_temp_printconfig(climate%temp)
    call spin_slc_printconfig(climate%slc)
  end subroutine spin_printconfig

  subroutine spin_climate(climate,model,time,pdd_scheme)
    !*FD do the spin climate forcing
    use glide_types
    use glimmer_global, only : rk    
    use glimmer_pdd
    implicit none
    type(spin_climate_type) :: climate  !*FD structure holding climate
    type(glide_global_type)   :: model !*FD model instance
    real(kind=rk), intent(in) :: time  !*FD current time
    type(glimmer_pdd_params) :: pdd_scheme
    
    call spin_eus(climate%slc,model,time)
    call spin_temp_routine(climate%temp,model,time)
    call spin_mb_routine(climate%mb,climate%temp,model,time,pdd_scheme)
  end subroutine spin_climate
    
end module spin_forcing
