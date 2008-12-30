! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! + spin_forcing.f90 - part of the GLIMMER ice model          + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include <config.inc>
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
    call write_log('Ant Spin Model')
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
