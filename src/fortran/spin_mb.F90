! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  spin_mb.f90 - part of the GLIMMER ice model              + 
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
#include <config.inc>
#endif
!The mass balance model
module spin_mb
   
  use glimmer_ts
  use glimmer_global, only : dp, sp, fname_length
  implicit none 
  type spin_mb_type
     integer :: model_type = 0
     integer :: ny, nx
        
     character(len=fname_length) :: fname=''     !*FD name of file containing
                                                 ! oxygen isotope
     integer :: use_simple = 0  !Use simple Eisment forcing (1) or more advanced
                                !Huybrecht forcing techniques(0)
      
     real(sp),dimension(:,:),pointer :: ablt !*FD Ablation
     real(sp),dimension(:,:),pointer :: presprcp  !*FD Present-day precip (water-equivalent)
     logical,dimension(:,:),pointer :: landsea !*FD Land-sea mask
     real(sp),dimension(:,:),pointer :: prcp !*FD Annual accumulation
     real :: scale2d_f1
     real :: oisotope = 0.0
     type(glimmer_tseries) :: mb_ts             !*FD ELA time series 
  end type spin_mb_type


contains
  subroutine spin_mb_config(config,mb)
    !*FD get configuration from config file
    use glimmer_config
    implicit none
    type(spin_mb_type)           :: mb     !*FD mb data
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   
    ! local variables
    type(ConfigSection), pointer :: section

    call GetSection(config,section,'SPIN MB')
    if (associated(section)) then
       call GetValue(section,'oxygen_isotope_file',mb%fname)
       call GetValue(section,'ny', mb%ny)
       call GetValue(section,'nx', mb%nx)
       call GetValue(section,'model_type', mb%model_type)
       call GetVAlue(section,'use_simple', mb%use_simple)
    end if

  end subroutine spin_mb_config
    
  subroutine spin_mb_printconfig(mb)
    !*FD print configuration to log
    use glimmer_log
    implicit none
    type(spin_mb_type)      :: mb   !*FD mb data
    ! local variables
    character(len=100) :: message
    call write_log('Spin MB')
    call write_log('-------')
    
    call write_log(message)
    if (mb%model_type == 0) then
      write(message,*) 'model type: ', mb%model_type, ' Antarctica Model'
      else if (mb%model_type == 1) then
        write(message,*) 'model type: ', mb%model_type, ' Greenland Model'
        write(message,*) 'oxygen isotope file : ',trim(mb%fname)
      else
        write(message,*) 'Model type is default: ', mb%model_type, &
        ' Antarctica Model'
    end if
    call write_log(message)
    if(mb%use_simple == 0) then
      write(message,*) 'use simple forcing: ', mb%use_simple, ' off'
    else
      write(message,*) 'use simple forcing: ', mb%use_simple, ' on'
    end if 
    call write_log(message)
    call write_log('')
  end subroutine spin_mb_printconfig

  subroutine spin_init_mb(mb, model)
    use glide_types
    use glimmer_paramets, only: thk0, scyr, tim0 
    implicit none 
    type(spin_mb_type)        :: mb   !*FD mb data
    type(glide_global_type)   :: model !*FD model instance
    !scaling
    mb%scale2d_f1 = scyr * thk0/tim0 
    
    call coordsystem_allocate(model%general%ice_grid, mb%prcp)
    call coordsystem_allocate(model%general%ice_grid, mb%ablt)
    call coordsystem_allocate(model%general%ice_grid, mb%landsea) 
    call coordsystem_allocate(model%general%ice_grid, mb%presprcp)

    if (mb%model_type == 1) then
      call glimmer_read_ts(mb%mb_ts,mb%fname)
    end if
 end subroutine spin_init_mb
    
  
  subroutine spin_mb_routine(mb,temp, model, time, pdd_scheme)
    !picks either the simple or advanced routine based on the use_simple variable
    use glimmer_pdd
    use glimmer_config
    use glide_types
    use spin_temp
    !use glimmer_paramets only: thk0, scyr, tim0 
    implicit none
    type(spin_mb_type)   :: mb !the mass balance data
    type(spin_temp_type) :: temp !the temp data
    type(glide_global_type)   :: model !*FD model instance
    type(glimmer_pdd_params)  :: pdd_scheme
    real(kind=rk), intent(in) :: time  !*FD current time
    !write(*,*), model%general%x0 
    if (mb%model_type == 1) then
      !get the value for the oxygen isotope in the case of Greenland
      call glimmer_ts_linear(mb%mb_ts,real(time),mb%oisotope)
    end if 
    if(mb%use_simple == 1) then
      call spin_simple_massbalance(pdd_scheme, &
                                   mb%model_type, &
                                   model%geometry%usrf, &
                                   mb%landsea, &
                                   temp%presartm, &
                                   model%climate%artm, &
                                   mb%prcp, &
                                   mb%presprcp, &
                                   model%climate%acab, &
                                   mb%ablt, &
                                   temp%arng)
                                   
                                  
    
    else
      call spin_massbalance(pdd_scheme, &
                            mb%model_type, &
                            model%geometry%usrf, &
                            mb%landsea, &
                            temp%tinvp, &
                            temp%presartm, &
                            model%climate%artm, &
                            mb%prcp, &
                            mb%presprcp, &
                            model%climate%acab, &
                            mb%ablt, &
                            temp%arng, &
                            model%general%ewn, &
                            model%general%nsn, &
                            mb%oisotope)
                            
                           

    end if  
  
    !glimmer requires rescaling for acab
    model%climate%acab = model%climate%acab/mb%scale2d_f1

  end subroutine spin_mb_routine 


  subroutine spin_simple_massbalance(pdd_scheme, model_type, usrf, landsea, &
                                     presartm, artm, prcp, presprcp, acab, &
                                     ablt, arng)
    !routine to calculate massbalance for Antarctica and Greenland using the
    !simple eismint drivers
    use glimmer_pdd 
    implicit none
    type(glimmer_pdd_params) :: pdd_scheme
    real(dp), dimension(:,:), intent(inout) :: usrf !ice elevation
    real(sp), dimension(:,:), intent(out) :: artm !mean annual temp var   
    real(sp), dimension(:,:), intent(out) :: arng !temp half range var
    real(sp), dimension(:,:), intent(out) :: ablt !ablation
    real(sp), dimension(:,:), intent(out) :: acab !accumlation/ablation
    real(sp), dimension(:,:), intent(inout) :: prcp !precipitation field var
    logical, dimension(:,:), intent(inout)  :: landsea !land/sea mask
    real(sp), dimension(:,:), intent(in) :: presartm !present mean annual temp
    real(sp), dimension(:,:), intent(in) :: presprcp !present precipitation
    integer, intent(in) :: model_type !(0) is Antarcitca, (1) is Greenland
    !real, intent(in) :: eus !eustatic sea level
    !real(dp), intent(in) :: topg !topgraphy field
    real(sp) :: pfac=1.0533 !*FD Precip enhancement factor (default is supposed EISMINT value)
    !calculate the landsea matrix which determines where the pdd method will be
    !used
    where(usrf > 0.0) 
      landsea =.True.
    elsewhere
      landsea = .False.
    end where 
   
    select case(model_type)
    case(0)
      !precipitation calculation for Antarctica
      prcp = 1.5 * 2**(artm/10)
    
    case(1)
      !precipitation calculataion for Greenland
      prcp = presprcp * pfac**(artm - presartm)
    
    end select 
    !send into the pdd calculation
    call glimmer_pdd_mbal(pdd_scheme, &
                          artm, &
                          arng, &
                          prcp, &
                          ablt, & 
                          acab, &
                          landsea) 
  end subroutine spin_simple_massbalance
  
  subroutine spin_massbalance(pdd_scheme, model_type, usrf, landsea,          &
                            tinvp, presartm, artm, prcp, presprcp, acab, ablt,&
                            arng,ewn, nsn, oisotope)

    use glimmer_pdd
    type(glimmer_pdd_params)              :: pdd_scheme !pdd type, holds
                                                        !parameters for pos deg day
    integer, intent(in)                   :: model_type !(0) is Antarctica
                                                        !(1) is Greenland
    logical, dimension(:,:),intent(inout)  :: landsea !land/sea mask
    real(sp),dimension(:,:),intent(out)   :: tinvp !present inversion temp
    real(sp),dimension(:,:),intent(out)   :: acab !accumlation/ablation
    real(sp),dimension(:,:),intent(out)   :: ablt !ablation
    real(sp),dimension(:,:),intent(in)      :: presartm !present mean annual temp
    real(sp),dimension(:,:),intent(in)      :: artm !mean annual temp field var
    real(sp),dimension(:,:),intent(in)      :: presprcp !present precipitation
    real(sp),dimension(:,:),intent(in)      :: arng !temp half range
    real(sp),dimension(:,:),intent(inout)   :: prcp !precipitation field var
    integer, intent(in)                   :: ewn, nsn !grid points
    real(sp), dimension(ewn,nsn)          :: tinv !inversion temp field var
    real(sp), dimension(ewn,nsn)          :: z
    real(dp), dimension(:,:),intent(inout)     :: usrf !ice elevation
    real                                  :: tzero = 273.16 !Kelvin 
    real, intent(in)                      :: oisotope !oxygen-isotope
    !real, intent(in) :: eus !eustatic sea level
    !real(dp), intent(in) :: topg !topgraphy field
    
    !set up the landsea matrix as a true/false map of land/sea
    where(usrf > 0.0) 
      landsea =.True.
    elsewhere
      landsea = .False.
    end where 
    select case(model_type)
    case(0) !Antarctica Model, model_type = 0

    
      !calculate the present inversion temp, add
      !tzero to get the right temp in Kelvin
      tinvp = 0.67 * presartm + 88.9 + tzero
    
      !calculate the new inversion temp through time, add tzero to get 
      !the right temp in Kelvin
      tinv = 0.67*artm + 88.9 + tzero
    
      !calculate the precip field
      prcp = presprcp * exp(22.47*(tzero/tinvp - &
      tzero/tinv))*((tinvp/tinv)**2)*(1 + 0.046*(tinv - tinvp)) 
    
    !end of Antarctica model type calculations

    case(1) !start of Greenland model calculations
    
    !calculate the perturbed precipitation using a file containing
    !oxygen isotopes instead of temperatures
      prcp = presprcp * exp(0.169*(oisotope + 34.83)) 
    !end of Green model calculations
    end select     
       
    !call the glimmer pdd scheme and send in the correct fields
    call glimmer_pdd_mbal(pdd_scheme, &
                          artm, &
                          arng, &
                          prcp, &
                          ablt, & 
                          acab, &
                          landsea) 
    
  end subroutine spin_massbalance
  
  




  
end module spin_mb
