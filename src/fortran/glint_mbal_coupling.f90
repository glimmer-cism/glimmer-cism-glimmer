
module glint_mbal_coupling

  use glint_mbal
  use glint_proj
  use glimmer_config

  !*FD Module to handle the accumulation of inputs and calculation of mass-balance

  type glint_mbc
    real(sp),dimension(:,:),pointer :: prcp_save => null() !*FD used to accumulate precip
    real(sp),dimension(:,:),pointer :: ablt_save => null() !*FD used to accumulate ablation
    real(sp),dimension(:,:),pointer :: acab_save => null() !*FD used to accumulate mass-balance
    real(sp),dimension(:,:),pointer :: artm_save => null() !*FD used to average air-temperature
    real(sp),dimension(:,:),pointer :: snowd     => null() !*FD Keeps track of snow depth
    real(sp),dimension(:,:),pointer :: siced     => null() !*FD Keeps track of superimposed ice depth 
    real(sp),dimension(:,:),pointer :: snowd_save => null() !*FD Saves snow depth
    real(sp),dimension(:,:),pointer :: siced_save => null() !*FD Saves superimposed ice depth 
    integer :: whichprcp=1 !*FD Option for precip calculation
    integer :: av_count =0 !*FD Counter for averaging temperature input
    logical :: new_accum=.true.
    type(glint_mbal_params) :: mbal
  end type glint_mbc

contains

  subroutine glint_mbc_init(params,proj,config,whichacab,whichprcp)

    type(glint_mbc)  :: params
    type(projection) :: proj
    type(ConfigSection), pointer :: config !*FD structure holding sections of configuration file
    integer          :: whichacab
    integer          :: whichprcp

    ! Deallocate if necessary

    if (associated(params%prcp_save))  deallocate(params%prcp_save)
    if (associated(params%ablt_save))  deallocate(params%ablt_save)
    if (associated(params%acab_save))  deallocate(params%acab_save)
    if (associated(params%artm_save))  deallocate(params%artm_save)
    if (associated(params%snowd))      deallocate(params%snowd)
    if (associated(params%siced))      deallocate(params%siced)
    if (associated(params%snowd_save)) deallocate(params%snowd_save)
    if (associated(params%siced_save)) deallocate(params%siced_save)

    ! Allocate arrays and zero

    allocate(params%prcp_save(proj%nx,proj%ny));  params%prcp_save = 0.0
    allocate(params%ablt_save(proj%nx,proj%ny));  params%ablt_save = 0.0
    allocate(params%acab_save(proj%nx,proj%ny));  params%acab_save = 0.0
    allocate(params%artm_save(proj%nx,proj%ny));  params%artm_save = 0.0
    allocate(params%snowd(proj%nx,proj%ny));      params%snowd = 0.0
    allocate(params%siced(proj%nx,proj%ny));      params%siced = 0.0
    allocate(params%snowd_save(proj%nx,proj%ny)); params%snowd_save = 0.0
    allocate(params%siced_save(proj%nx,proj%ny)); params%siced_save = 0.0

    ! Initialise the mass-balance scheme and other components

    call glint_mbal_init(params%mbal,config,whichacab)
    params%whichprcp=whichprcp

  end subroutine glint_mbc_init

  ! +++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_accumulate(params,artm,arng,prcp,snowd,siced,xwind,ywind,global_orog,local_orog)

    type(glint_mbc)  :: params
    real(sp),dimension(:,:),intent(inout) :: artm      !*FD Mean air temperature (degC)
    real(sp),dimension(:,:),intent(in) :: arng         !*FD Air temperature half-range (degC)
    real(sp),dimension(:,:),intent(inout) :: prcp      !*FD Precipitation (m)
    real(sp),dimension(:,:),intent(in) :: snowd        !*FD Snow depth (m)
    real(sp),dimension(:,:),intent(in) :: siced        !*FD Superimposed ice depth (m)
    real(rk),dimension(:,:),intent(in) :: xwind        !*FD $x$-component of surface winds (m/s)
    real(rk),dimension(:,:),intent(in) :: ywind        !*FD $y$-component of surface winds (m/s)
    real(dp),dimension(:,:),intent(in) :: global_orog  !*FD Global orography (m)
    real(sp),dimension(:,:),intent(in) :: local_orog   !*FD Local orography (m)

    real(sp),dimension(size(artm,1),size(artm,2)) :: ablt,acab
    
    ! Things to do the first time

    if (params%new_accum) then

       params%new_accum=.false.
       params%av_count =0

       ! Initialise 

       params%snowd=snowd
       params%siced=siced
       params%snowd_save=snowd
       params%siced_save=siced

       params%prcp_save=0.0
       params%ablt_save=0.0
       params%acab_save=0.0
       params%artm_save=0.0

    end if

    params%av_count=params%av_count+1

    ! Call mass-balance

    call glint_mbal_calc(params%mbal,artm,arng,prcp,params%snowd,params%siced,ablt,acab) 

    ! Accumulate

    params%prcp_save = params%prcp_save + prcp
    params%ablt_save = params%ablt_save + ablt
    params%acab_save = params%acab_save + acab
    params%artm_save = params%artm_save + artm

  end subroutine glint_accumulate

  ! +++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_get_mbal(params,artm,prcp,ablt,acab,snowd,siced)

    type(glint_mbc)  :: params
    real(sp),dimension(:,:),intent(out)   :: artm   !*FD Mean air temperature (degC)
    real(sp),dimension(:,:),intent(out)   :: prcp   !*FD Precipitation (m)
    real(sp),dimension(:,:),intent(out)   :: ablt   !*FD Ablation
    real(sp),dimension(:,:),intent(out)   :: acab   !*FD Mass-balance
    real(sp),dimension(:,:),intent(inout) :: snowd  !*FD Snow depth (m)
    real(sp),dimension(:,:),intent(inout) :: siced  !*FD Superimposed ice depth (m)

    if (.not.params%new_accum) then
       params%snowd=params%snowd-params%snowd_save
       params%siced=params%siced-params%siced_save
       params%artm_save=params%artm_save/real(params%av_count)
    end if

    params%new_accum=.true.

    artm=params%artm_save
    prcp=params%prcp_save
    ablt=params%ablt_save
    acab=params%acab_save
    snowd=snowd+params%snowd
    siced=siced+params%siced

    where (snowd<0.0) snowd=0.0
    where (siced<0.0) siced=0.0

  end subroutine glint_get_mbal

end module glint_mbal_coupling
