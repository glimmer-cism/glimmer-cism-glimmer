
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

     real(dp) :: pddfs                  !*FD Later set to \texttt{(rhow / rhoi) * pddfac\_snow}
     real(dp) :: pddfi                  !*FD Later set to \texttt{(rhow / rhoi) * pddfac\_ice}
     real(dp) :: wmax        = 0.6_dp   !*FD Fraction of melted snow that refreezes
     real(dp) :: pddfac_ice  = 0.008_dp !*FD PDD factor for ice (m day$^{-1}$ $^{\circ}C$^{-1}$)
     real(dp) :: pddfac_snow = 0.003_dp !*FD PDD factor for snow (m day$^{-1}$ $^{\circ}C$^{-1}$)
     real(dp) :: rain_threshold = 1.0_dp !*FD Threshold for precip melting (degC)

  end type glimmer_daily_pdd_params

  real(dp),parameter :: one_over_pi=1.0_dp/pi

contains

  subroutine glimmer_daily_pdd_init(params)

    use physcon, only : rhow, rhoi
  
    type(glimmer_daily_pdd_params) :: params

    params%pddfs = (rhow / rhoi) * 0.003_dp 
    params%pddfi = (rhow / rhoi) * 0.008_dp 

  end subroutine glimmer_daily_pdd_init

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glimmer_daily_pdd_mbal(params,artm,arng,prcp,snowd,siced,ablt,acab)

    type(glimmer_daily_pdd_params)        :: params !*FD Daily PDD scheme parameters
    real(dp),dimension(:,:),intent(in)    :: artm   !*FD Daily mean air-temperature ($^{\circ}$C)
    real(dp),dimension(:,:),intent(in)    :: arng   !*FD Daily temperature half-range ($^{\circ}$C)
    real(dp),dimension(:,:),intent(in)    :: prcp   !*FD Daily precipitation (m)
    real(dp),dimension(:,:),intent(inout) :: snowd  !*FD Snow depth (m)
    real(dp),dimension(:,:),intent(inout) :: siced  !*FD Superimposed ice depth (m)
    real(dp),dimension(:,:),intent(out)   :: ablt   !*FD Daily ablation (m)
    real(dp),dimension(:,:),intent(out)   :: acab   !*FD Daily mass-balance (m)

    real(dp),dimension(size(prcp,1),size(prcp,2)) :: rain ! Daily rain
    real(dp),dimension(size(prcp,1),size(prcp,2)) :: degd ! Degree-day
    real(dp),dimension(size(prcp,1),size(prcp,2)) :: giced ! Temporary array for glacial ice
    real(dp),dimension(size(prcp,1),size(prcp,2)) :: old_snow,old_sice

    integer :: nx,ny,i,j

    nx=size(prcp,1) ; ny=size(prcp,2)

    old_snow=snowd
    old_sice=siced

    rain=rainorsnw(1,artm,arng,prcp,params%rain_threshold)
    degd=finddegdays(artm,arng)
    giced=0.0_dp

    do i=1,nx
       do j=1,ny
          call degdaymodel(params,snowd(i,j),siced(i,j),giced(i,j),degd(i,j),rain(i,j),prcp(i,j)) 
       end do
    end do

    acab=snowd+siced+giced-old_snow-old_sice
    ablt=prcp-acab

  end subroutine glimmer_daily_pdd_mbal

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  elemental real(dp) function finddegdays(localartm,localarng)

    !*FD Finds the degree-day as a function of mean daily temperature and
    !*FD half-range. The calculation is made on the assumption that
    !*FD daily temperatures vary as $T = T_{0} - \Delta T * \cos(\theta)$.

    real(dp), intent(in) :: localartm !*FD Mean daily temperature (degC)
    real(dp), intent(in) :: localarng !*FD Daily temperture half-range (degC)

    real(dp) :: time

    if (localartm - localarng > 0.0_dp) then
       finddegdays = localartm
    else if (localartm + localarng < 0.0_dp) then  
       finddegdays = 0.0_dp
    else
       time = acos(localartm / localarng) 
       finddegdays = (localartm*(pi-time)+localarng*sin(time))*one_over_pi
    end if

  end function finddegdays

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  elemental real(dp) function rainorsnw(which,localartm,localarng,localprcp,threshold)

    !*FD Determines a value for snow precipitation, dependent on air temperature and half-range.
    !*FD Takes precipitation as input and returns amount of rain.

    integer, intent(in) :: which     !*FD Selects method of calculation
    real(dp),intent(in) :: localartm !*FD Air temperature (degC)
    real(dp),intent(in) :: localarng !*FD Temperature half-range (degC)
    real(dp),intent(in) :: localprcp !*FD Precipitation (arbitrary units)
    real(dp),intent(in) :: threshold !*FD Snow/rain threshold (degC)

    select case(which)

    case(1)

       ! Assume temperature evolution is sinusoidal

       if (localartm - localarng > threshold) then
          rainorsnw = localprcp
       else if (localartm + localarng < threshold) then  
          rainorsnw = 0.0_dp
       else
          rainorsnw = localprcp * (1.0_dp-one_over_pi * acos((localartm - threshold) / localarng))
       end if

    case default

       ! Just use mean temperature to determine if precip is snow or rain

       if (localartm > threshold) then
          rainorsnw = localprcp
       else
          rainorsnw = 0.0_dp
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
    real(dp), intent(inout) :: snowdepth     !*FD Snow depth (m)
    real(dp), intent(inout) :: sicedepth     !*FD Superimposed ice depth (m)
    real(dp), intent(inout) :: gicedepth     !*FD Glacial ice depth (m)
    real(dp), intent(in)    :: degd          !*FD The degree-day (degC day)
    real(dp), intent(in)    :: prcp          !*FD Total precip (m)
    real(dp), intent(in)    :: rain          !*FD Rain (m)

    real(dp) :: potablt, wfrac

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
          snowdepth = 0.0_dp 

	  sicedepth = sicedepth - potablt

          if (sicedepth < 0.0_dp) then

             gicedepth = gicedepth + sicedepth
             sicedepth = 0.0_dp
          end if

       end if
    end if

  end subroutine degdaymodel

end module glimmer_daily_pdd
