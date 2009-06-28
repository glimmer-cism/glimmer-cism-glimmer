! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  spin_temp.f90 - part of the GLIMMER ice model             + 
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

module spin_temp
  
  !*FD temperature forcing

  use glimmer_ts
  use glimmer_global, only : sp, dp, fname_length

  type spin_temp_type
     integer :: model_type = 0  !determines the model type Antarctica(0) or Greenland(1)
     integer :: use_simple = 0  !Use simple Eisment forcing (1) or more advanced
                                !Huybrecht forcing techniques(0)
     character(len=fname_length) :: fname=''     !*FD name of file containing temperature ts
     character(len=fname_length) :: ele_fname=''
     type(glimmer_tseries) :: temp_ts            !*FD temperature time series 
     type(glimmer_tseries) :: ele_ts
     real, dimension(:),pointer :: tperturb        !*FD temperature value
     !real :: tvalue
     integer:: torder
     real(sp),dimension(:,:),pointer :: arng !*FD Surface temp half-range
     real(sp),dimension(:,:),pointer :: presartm  !*FD Present-day surface temperature
     real(sp),dimension(:,:),pointer :: tinvp ! Inversion point temperature 
     real(dp),dimension(:,:),pointer :: local_usrf !local ice surface field for
                                                   !use in scaling
     integer,dimension(2) :: core_x_y !x and y of the core location to
                                           !calc the elevation change
     real(sp) :: temp_ele = 0.0 !temperature correction due to elevation changes at
                          !the site of the ice core
     real :: ele_correction = 0.0 !elevation value for correction 
  end type spin_temp_type

contains

subroutine spin_temp_config(config,temp)
    !*FD get temperature configuration from config file
    use glimmer_config
    implicit none
    type(spin_temp_type)           :: temp     !*FD mb data
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   
    ! local variables
    type(ConfigSection), pointer :: section
    integer, dimension(:), pointer :: dummy

    call GetSection(config,section,'SPIN Temperature')
    if (associated(section)) then
       call GetValue(section,'temp_file',temp%fname)
       call GetValue(section,'model_type', temp%model_type)
       call GetValue(section,'use_simple', temp%use_simple)
    end if
    
    call GetSection(config,section,'SPIN Elevation Correction')
    if (associated(section)) then
       call GetValue(section,'ele_file',temp%ele_fname) 
    end if
    
  
  end subroutine spin_temp_config

  subroutine spin_temp_printconfig(temp)
    !*FD print configuration to log
    use glimmer_log
    implicit none
    type(spin_temp_type)      :: temp   !*FD temperature data
    ! local variables
    character(len=100) :: message
    call write_log('SPIN Temperature')
    call write_log('---------------')
    write(message,*) 'temperature file  : ',trim(temp%fname)
    call write_log(message)
    if (temp%model_type == 0) then
      write(message,*) 'model type: ', temp%model_type, ' Antarctica Model'
    else if (temp%model_type == 1) then
      write(message,*) 'model type: ', temp%model_type, ' Greenland Model &
      --Huybrechts Temperature Equations'
    else if (temp%model_type == 2) then
      write(message,*) 'model type: ', temp%model_type, ' Greenland Model & 
      --Fausto Temperature Equations'
    else
      write(message,*) 'Model type is default: ', temp%model_type, ' Antarctica &
                       Model'
    end if
    call write_log(message)
    if(temp%use_simple == 0) then
      write(message,*) 'use simple: ', temp%use_simple, ' off'
    else
      write(message,*) 'use simple: ', temp%use_simple, ' on'
    end if 
    call write_log(message)
    call write_log('')
  end subroutine spin_temp_printconfig

  subroutine spin_init_temp(temp, model)
    !*FD initialise temperature forcing
    use glide_types
    implicit none
    type(spin_temp_type)     :: temp  !*FD mb data
    type(glide_global_type)  :: model !model data
    
    call glimmer_read_ts(temp%temp_ts,temp%fname,1)
    allocate(temp%tperturb(1))
    call coordsystem_allocate(model%general%ice_grid, temp%local_usrf) 
    call coordsystem_allocate(model%general%ice_grid, temp%presartm)
    call coordsystem_allocate(model%general%ice_grid, temp%arng)
    call coordsystem_allocate(model%general%ice_grid, temp%tinvp)
    
    if (len(trim(temp%ele_fname)) .ne. 0) then
      call glimmer_read_ts(temp%ele_ts,temp%ele_fname,1)
    end if
  end subroutine spin_init_temp

  subroutine spin_temp_routine(temp, model, time)
    ! choose the routine to calculate surface temperature
    use glide_types
    use glimmer_global, only : sp,dp,rk
    use glimmer_paramets, only: thk0
    implicit none
    type(spin_temp_type)      :: temp  !*FD temperature data
    type(glide_global_type)   :: model !*FD model instance
    real(kind=rk), intent(in) :: time  !*FD current time
   
    
    temp%local_usrf = model%geometry%usrf*thk0
    call glimmer_ts_linear(temp%temp_ts,real(time),temp%tperturb)
    
    if (len(trim(temp%ele_fname)) .ne. 0) then
      call glimmer_ts_linear(temp%ele_ts,real(time),temp%ele_correction)
    end if

    select case(temp%use_simple)

    case(0) !use the more advanced temperature forcing as given by Huybrechts
      call spin_surftemp(temp%model_type, &
                         temp%arng, &
                         model%climate%artm, &
                         temp%local_usrf, &
                         model%climate%lati, &
                         model%climate%eus, &
                         temp%tperturb(1), &
                         model%general%ewn, &
                         model%general%nsn, &
                         model%climate%loni, &
                         temp%core_x_y, &
                         temp%temp_ele, &
                         temp%ele_correction)


    case(1)  !the simple Eismint type temperature forcing use_simple = 0
      call spin_simple_surftemp(temp%model_type, &
                         temp%arng, &
                         model%climate%artm, &
                         temp%local_usrf, &
                         model%climate%lati, &
                         model%climate%eus, &
                         temp%tperturb(1), &
                         model%general%ewn, &
                         model%general%nsn)
    end select
  
  model%climate%tempanmly = temp%tperturb(1)
  end subroutine spin_temp_routine
  
  
  subroutine spin_simple_surftemp(model_type,arng, artm, usrf, lati, eus,tperturb, &
                                  ewn, nsn)
    integer, intent(in) :: model_type !*FDreferences which model to run (0)
                                      !*FD Antarctica, (1) Greenland
    real(sp), dimension(:,:), intent(inout) :: arng !*FD temperature half-range
    real(sp), dimension(:,:), intent(inout) :: artm !*FD mean annual temperature
    real(dp), dimension(:,:), intent(in) :: usrf !*FD surface elevation
    real(sp), dimension(:,:), intent(in) :: lati !*FD latitude
    real :: eus !eustatic sea level
    real, intent(in) :: tperturb !temperature forcing over time
    integer, intent(in) :: ewn, nsn !# of grid points 
    real, dimension(ewn, nsn) :: glandhinv
    select case(model_type)
    
    case(0) !Antarctica Temperature Model, model_type = 0
    
      !calculate the artm for Antarctica following the eismint  method
      artm = 34.46  - 0.00914*usrf + 0.68775*lati + tperturb
    
      !calculate the temperature half range by subtracting the artm from the 
      !summer temperature
      arng = 16.81 - 0.00692*usrf +  0.27937*lati + 2.0*tperturb -  artm 
    
    case(1) !Greenland Temperature Model, model_type = 1
      
      !create and calculate the inversion temperature elevation

      glandhinv = 20*(lati - 65)


      !calculate the mean annual temperature
      where (usrf >= glandhinv)
    
        artm = 49.13 - 0.007992 * usrf - 0.7576 * lati + tperturb

      elsewhere
        artm = 49.13 - 0.007992 * glandhinv - 0.7576 * lati + tperturb
      end where 
    
      !calculate the temperature half-range
      arng = 30.38 - 0.006277 * glandhinv - 0.3262 * lati - artm + 2.0*tperturb 
    end select 
  
  
  end subroutine spin_simple_surftemp

  subroutine spin_surftemp(model_type, arng, artm, usrf, lati, eus, tperturb, &
                           ewn, nsn, loni, core_x_y, temp_ele, ele_correction)
    
    integer, intent(in) :: model_type !references which model to run (0)
                                      !Antarctica, (1) Greenland
    real(sp), dimension(:,:), intent(out) :: arng !temperature half-range
    real(sp), dimension(:,:), intent(out) :: artm !mean annual temperature
    real(dp), dimension(:,:), intent(in) :: usrf !surface elevation
    real(sp), dimension(:,:), intent(in) :: lati !latitude
    real :: eus !eustatic sea level
    real, intent(inout) :: tperturb !temperature forcing over time
    integer, intent(in) :: ewn, nsn !# of grid points 
    real(sp), dimension(:,:), intent(in) :: loni
    real, dimension(ewn, nsn) :: glandhinv
    integer ns, ew   
    integer, dimension(2), intent(in)  :: core_x_y !the x,y grid values of the
                                                    !elevation at the site of the ice core
    real(sp), intent(inout) :: temp_ele !The correction to temperature due to changes in elevation
                                        !at the site of the ice core
    real, intent(in) :: ele_correction
    !make adjustments for elevation changes at the core site here
    if (ele_correction .ne. 0.0) then
      call calc_elevation_change(ele_correction, temp_ele,model_type)
      tperturb = tperturb + temp_ele
    end if

    select case(model_type)
    case(0) !Antarctica Temperature Model, model_type = 0
      !calculate the artm for Antarctica following Huybrechts method
      where (usrf <= 1500.)
        artm = 34.46  - 0.005102*usrf + 0.68775*lati + tperturb
    
      elsewhere    
        artm = 34.46  - 0.014285*usrf + 0.68775*lati + tperturb

      end where
    
      !calculate the temperature half range by subtracting the artm from the 
      !summer temperature
      arng = 16.81 - 0.00692*usrf + 0.27937*lati - artm + 2.0*tperturb 
    case(1) !Greenland Temp Model model_type 1 and 2
      !Huybrecht's temperature method for Greenland
      !create and calculate the inversion temperature elevation
      glandhinv = 20*(lati - 65)


      !calculate the mean annual temperature
      where (usrf >= glandhinv)
        artm = 49.13 - 0.007992 * usrf - 0.7576 * lati + tperturb
      elsewhere
        artm = 49.13 - 0.007992 * glandhinv - 0.7576 * lati + tperturb
      end where 
      
      !calculate the temperature half-range
      arng = 30.78  - 0.006277 * usrf - 0.3262 * lati - artm + 2.0*tperturb 
    
    case(2) !Fausto Temperature 
    
    
      !The paramaterization based on Fausto et al 2009
      !Note: the parameterization accounts for longitudinal changes in
      !temperature, but uses degrees W, whereas the current data set uses degrees E
    
     
      !without land stations
      !artm = 40.96 - (6.988*usrf/1000.0) - 0.6901*lati + 0.0742*(-loni) + tperturb
      !arng = 13.46 - (5.172*usrf/1000.0) - 0.1809*lati + 0.1049*(-loni) - artm + 2.0* tperturb

      !with land stations
      artm = 41.83 - (6.309*usrf/1000.0) - 0.7189*lati + 0.0672*(-loni) + tperturb
      arng = 14.70 - (5.426*usrf/1000.0) - 0.1585*lati + 0.0518*(-loni) - artm + 2.0* tperturb
    end select
  end subroutine spin_surftemp

  subroutine calc_elevation_change(core_elevation, temp_ele,model_type)
    implicit none
    real, intent(out) :: temp_ele  !The correction to temperature due to changes in elevation
                                  !at the site of the ice core
    real, intent(in) :: core_elevation !surface elevation at the core
    real  :: core_elevation_last = 0.0 !core elevation at the last time step
    real :: elevation_diff !difference between the present and last time
                               !step elevation
    integer, intent(in) :: model_type !references which model to run (0)
                                      !Antarctica, (1) Greenland
    !calculate the difference in surface elevation at each time step to
    !correct the temperature
    if (core_elevation_last .eq. 0.0) then 
      elevation_diff = 0.0
    else
      elevation_diff = elevation_diff +(core_elevation - core_elevation_last)
    end if 
    if (model_type .eq. 0) then
    
      temp_ele = temp_ele + elevation_diff * 0.014285 !our assumed lapse rate
    else  
      temp_ele = temp_ele + elevation_diff * 0.00792 !our assumed lapse rate
    end if
    core_elevation_last = core_elevation
  end subroutine calc_elevation_change
end module spin_temp
