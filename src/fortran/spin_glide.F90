! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  spin_glide.f90 - part of the GLIMMER ice model            + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

program spin_glide
  use glimmer_global, only:rk
  use glide
  use spin_forcing
  use spin_io
  use glimmer_log
  use glimmer_config
  use spin_types
  use glimmer_pdd
  implicit none

  type(glimmer_pdd_params) :: pdd_scheme
  type(glide_global_type) :: model        ! model instance
  type(spin_climate_type) :: climate       ! climate
  type(ConfigSection), pointer :: config  ! configuration stuff
  character(len=fname_length) :: fname   ! name of paramter file
  real(kind=rk) time
  


  write(*,*) 'Enter name of GLIDE configuration file to be read'
  read(*,*) fname
  
  ! start logging
  call open_log(unit=50, fname=logname(fname))

  ! read configuration
  call ConfigRead(fname,config)

  ! initialise GLIDE
  call glide_config(model,config)
  call glide_initialise(model)
  call spin_initialise(climate,config,model,pdd_scheme)
  !call spin_initialise(climate,config,model)
  call CheckSections(config)
  ! fill dimension variables
  call glide_nc_fillall(model)

  time = model%numerics%tstart
 ! call spin_climate(climate,model,time)
  call spin_climate(climate,model,time, pdd_scheme)
  call spinup_lithot(model)

  do while(time.le.model%numerics%tend)    
     call glide_tstep_p1(model,time)
     call spin_io_writeall(climate,model)
     call glide_tstep_p2(model)
     call glide_tstep_p3(model)
     ! override masking stuff for now
     time = time + model%numerics%tinc
     call spin_climate(climate,model,time, pdd_scheme)
     !call spin_climate(climate,model,time)
  end do

  ! finalise GLIDE
  call glide_finalise(model)
  call close_log

end program spin_glide
