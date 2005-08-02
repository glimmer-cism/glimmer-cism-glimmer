
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_daily_pdd.f90 - part of the GLIMMER ice model    + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004 GLIMMER contributors - see COPYRIGHT file 
! for list of contributors.
!
! This program is free software; you can redistribute it and/or 
! modify it under the terms of the GNU General Public License as 
! published by the Free Software Foundation; either version 2 of 
! the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License 
! along with this program; if not, write to the Free Software 
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 
! 02111-1307 USA
!
! GLIMMER is maintained by:
!
! Ian Rutt
! School of Geographical Sciences
! University of Bristol
! University Road
! Bristol
! BS8 1SS
! UK
!
! email: <i.c.rutt@bristol.ac.uk> or <ian.rutt@physics.org>
!
! GLIMMER is hosted on NeSCForge:
!
! http://forge.nesc.ac.uk/projects/glimmer/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module glimmer_daily_pdd

  use glimmer_global
  use physcon, only : pi

  implicit none

  type glimmer_daily_pdd_params

     !*FD Holds parameters for daily positive-degree-day mass-balance
     !*FD calculation. 

     real(sp) :: pddfs                  !*FD Later set to \texttt{(rhow / rhoi) * pddfac\_snow}
     real(sp) :: pddfi                  !*FD Later set to \texttt{(rhow / rhoi) * pddfac\_ice}
     real(sp) :: wmax        = 0.6_sp   !*FD Fraction of melted snow that refreezes
     real(sp) :: pddfac_ice  = 0.008_sp !*FD PDD factor for ice (m day$^{-1}$ $^{\circ}C$^{-1}$)
     real(sp) :: pddfac_snow = 0.003_sp !*FD PDD factor for snow (m day$^{-1}$ $^{\circ}C$^{-1}$)
     real(sp) :: rain_threshold = 1.0_sp !*FD Threshold for precip melting (degC)
     integer  :: whichrain = 1  !*FD method for determining whether precip is rain or snow.

  end type glimmer_daily_pdd_params

  real(sp),parameter :: one_over_pi=1.0_sp/pi

  private
  public :: glimmer_daily_pdd_params, glimmer_daily_pdd_init, glimmer_daily_pdd_mbal

contains

  subroutine glimmer_daily_pdd_init(params,config)

    use physcon, only : rhow, rhoi
    use glimmer_config
  
    type(glimmer_daily_pdd_params) :: params
    type(ConfigSection), pointer         :: config !*FD structure holding sections of configuration file

    ! Read the config file and output to log

    call daily_pdd_readconfig(params,config)
    call daily_pdd_printconfig(params)

    params%pddfs = (rhow / rhoi) * params%pddfac_snow
    params%pddfi = (rhow / rhoi) * params%pddfac_ice

  end subroutine glimmer_daily_pdd_init

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine daily_pdd_readconfig(params,config)

    !*FD Reads in configuration data for the annual PDD scheme.

    use glimmer_config

    type(glimmer_daily_pdd_params),intent(inout) :: params !*FD The positive-degree-day parameters
    type(ConfigSection), pointer         :: config !*FD structure holding sections of configuration file   

    ! local variables
    type(ConfigSection), pointer :: section
    
    call GetSection(config,section,'GLIMMER daily pdd')
    if (associated(section)) then
       call GetValue(section,'wmax',params%wmax)
       call GetValue(section,'pddfac_ice',params%pddfac_ice)
       call GetValue(section,'pddfac_snow',params%pddfac_snow)
       call GetValue(section,'rain_threshold',params%rain_threshold)
       call GetValue(section,'whichrain',params%whichrain)
    end if

  end subroutine daily_pdd_readconfig

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine daily_pdd_printconfig(params)

    use glimmer_log

    type(glimmer_daily_pdd_params),intent(inout) :: params !*FD The positive-degree-day parameters
    character(len=100) :: message

    call write_log_div

    call write_log('GLIMMER daily PDD Scheme parameters:')
    call write_log('-----------------------------------')
    write(message,*) 'Snow refreezing fraction',params%wmax
    call write_log(message)
    write(message,*) 'PDD factor for ice',params%pddfac_ice
    call write_log(message)
    write(message,*) 'PDD factor for snow',params%pddfac_snow
    call write_log(message)
    write(message,*) 'Rain threshold temperature',params%rain_threshold,' degC'
    call write_log(message)
    write(message,*) 'Rain/snow partition method',params%whichrain
    call write_log('')

  end subroutine daily_pdd_printconfig

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glimmer_daily_pdd_mbal(params,artm,arng,prcp,snowd,siced,ablt,acab,landsea)

    type(glimmer_daily_pdd_params)        :: params !*FD Daily PDD scheme parameters
    real(sp),dimension(:,:),intent(in)    :: artm   !*FD Daily mean air-temperature ($^{\circ}$C)
    real(sp),dimension(:,:),intent(in)    :: arng   !*FD Daily temperature half-range ($^{\circ}$C)
    real(sp),dimension(:,:),intent(in)    :: prcp   !*FD Daily precipitation (m)
    real(sp),dimension(:,:),intent(inout) :: snowd  !*FD Snow depth (m)
    real(sp),dimension(:,:),intent(inout) :: siced  !*FD Superimposed ice depth (m)
    real(sp),dimension(:,:),intent(out)   :: ablt   !*FD Daily ablation (m)
    real(sp),dimension(:,:),intent(out)   :: acab   !*FD Daily mass-balance (m)
    logical, dimension(:,:),intent(in)    :: landsea !*FD Land-sea mask (land is TRUE)

    real(sp),dimension(size(prcp,1),size(prcp,2)) :: rain ! Daily rain
    real(sp),dimension(size(prcp,1),size(prcp,2)) :: degd ! Degree-day
    real(sp),dimension(size(prcp,1),size(prcp,2)) :: giced ! Temporary array for glacial ice
    real(sp),dimension(size(prcp,1),size(prcp,2)) :: old_snow,old_sice

    integer :: nx,ny,i,j

    nx=size(prcp,1) ; ny=size(prcp,2)

    old_snow=snowd
    old_sice=siced

    rain=rainorsnw(params%whichrain,artm,arng,prcp,params%rain_threshold)
    degd=finddegdays(artm,arng)
    giced=0.0_sp

    do i=1,nx
       do j=1,ny
          if (landsea(i,j)) then
             call degdaymodel(params,snowd(i,j),siced(i,j),giced(i,j),degd(i,j),rain(i,j),prcp(i,j)) 
             acab(i,j)=snowd(i,j)+siced(i,j)+giced(i,j)-old_snow(i,j)-old_sice(i,j)
             ablt(i,j)=prcp(i,j)-acab(i,j)
          else
             ablt(i,j)=prcp(i,j)
             acab(i,j)=0.0
             snowd(i,j)=0.0
             siced(i,j)=0.0
          end if
       end do
    end do

  end subroutine glimmer_daily_pdd_mbal

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  elemental real(sp) function finddegdays(localartm,localarng)

    !*FD Finds the degree-day as a function of mean daily temperature and
    !*FD half-range. The calculation is made on the assumption that
    !*FD daily temperatures vary as $T = T_{0} - \Delta T * \cos(\theta)$.

    real(sp), intent(in) :: localartm !*FD Mean daily temperature (degC)
    real(sp), intent(in) :: localarng !*FD Daily temperture half-range (degC)

    real(sp) :: time

    if (localartm - localarng > 0.0_sp) then
       finddegdays = localartm
    else if (localartm + localarng < 0.0_sp) then  
       finddegdays = 0.0_sp
    else
       time = acos(localartm / localarng) 
       finddegdays = (localartm*(pi-time)+localarng*sin(time))*one_over_pi
    end if

  end function finddegdays

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  elemental real(sp) function rainorsnw(which,localartm,localarng,localprcp,threshold)

    !*FD Determines a value for snow precipitation, dependent on air temperature and half-range.
    !*FD Takes precipitation as input and returns amount of rain.

    integer, intent(in) :: which     !*FD Selects method of calculation
    real(sp),intent(in) :: localartm !*FD Air temperature (degC)
    real(sp),intent(in) :: localarng !*FD Temperature half-range (degC)
    real(sp),intent(in) :: localprcp !*FD Precipitation (arbitrary units)
    real(sp),intent(in) :: threshold !*FD Snow/rain threshold (degC)

    select case(which)

    case(1)

       ! Assume temperature evolution is sinusoidal

       if (localartm - localarng > threshold) then
          rainorsnw = localprcp
       else if (localartm + localarng < threshold) then  
          rainorsnw = 0.0_sp
       else
          rainorsnw = localprcp * (1.0_sp-one_over_pi * acos((localartm - threshold) / localarng))
       end if

    case default

       ! Just use mean temperature to determine if precip is snow or rain

       if (localartm > threshold) then
          rainorsnw = localprcp
       else
          rainorsnw = 0.0_sp
       end if

    end select

  end function rainorsnw

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine degdaymodel(params,snowdepth,sicedepth,gicedepth,degd,rain,prcp)

    !*FD Applies the positive degree-day model.
    !*FD The output of this subroutine is in the variables \texttt{snowdepth},
    !*FD \texttt{sicedepth}, and \texttt{gicedepth}, which give the new depths
    !*FD of snow, superimposed ice and glacial ice, respectively.

    type(glimmer_daily_pdd_params) :: params !*FD PDD parameters
    real(sp), intent(inout) :: snowdepth     !*FD Snow depth (m)
    real(sp), intent(inout) :: sicedepth     !*FD Superimposed ice depth (m)
    real(sp), intent(inout) :: gicedepth     !*FD Glacial ice depth (m)
    real(sp), intent(in)    :: degd          !*FD The degree-day (degC day)
    real(sp), intent(in)    :: prcp          !*FD Total precip (m)
    real(sp), intent(in)    :: rain          !*FD Rain (m)

    real(sp) :: potablt, wfrac

    !-------------------------------------------------------------------------

    ! add snow to snow depth

    snowdepth = snowdepth + prcp - rain

    ! assume that rain has gone into superimposed ice

    sicedepth = sicedepth + rain

    ! this is the depth of superimposed ice that would need to be
    ! melted before runoff can occur 

    wfrac = params%wmax * (snowdepth + sicedepth)

    ! this is the total potential ablation of SNOW

    potablt = degd * params%pddfs

    ! if not enough snow to cause sice/firn > wfrac
    ! no net mass loss but melted snow becomes sice

    if ( potablt + sicedepth <= wfrac ) then

       snowdepth = snowdepth - potablt
       sicedepth = sicedepth + potablt

    else 

       ! have exceeded ablation threshold. now deduct
       ! potential ablation needed to get sice up to wfrac
       ! and set sice to this threshold.

       potablt = potablt - (wfrac - sicedepth)

       snowdepth = snowdepth - (wfrac - sicedepth)
       sicedepth = wfrac

       ! now determine what to do with remaining potential
       ! ablation.  if this is not sufficient to get rid of
       ! all remaining snow then 

       if (potablt < snowdepth) then

          snowdepth = snowdepth - potablt 

       else

          ! now melting ice so change reminaing potential to ice
          ! start to melt superimposed ice and continue to melt glacier
          ! ice if necessary


          potablt = params%pddfi * (potablt - snowdepth) / params%pddfs
          snowdepth = 0.0_sp 

	  sicedepth = sicedepth - potablt

          if (sicedepth < 0.0_sp) then

             gicedepth = gicedepth + sicedepth
             sicedepth = 0.0_sp
          end if

       end if
    end if

  end subroutine degdaymodel

end module glimmer_daily_pdd
