! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_initialise.f90 - part of the GLIMMER ice model     + 
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

module glint_initialise

  !*FD Initialise GLINT model instance

  use glint_type

  private
  public glint_i_initialise, glint_i_end, calc_coverage, get_i_upscaled_fields

contains

  subroutine glint_i_initialise(config,instance,radea,grid,time_step,start_time)

    !*FD Initialise a GLINT ice model instance

    use glimmer_config
    use glint_global_grid
    use glint_io
    use glide
    implicit none

    ! Arguments
    type(ConfigSection), pointer         :: config      !*FD structure holding sections of configuration file   
    type(glint_instance),  intent(inout) :: instance    !*FD The instance being initialised.
    real(rk),              intent(in)    :: radea       !*FD Radius of the earth (m).
    type(global_grid),     intent(in)    :: grid        !*FD Global grid to use
    real(rk),              intent(in)    :: time_step   !*FD Model time-step (years).
    real(rk),optional,     intent(in)    :: start_time  !*FD Start time of model (years).

    ! initialise model

    call glide_initialise(instance%model,config)

    ! create glint variables

    call glint_io_createall(instance%model)

    instance%model%numerics%tinc=time_step             ! Initialise the model time step

    ! read glint configuration

    call glint_i_readconfig(instance,config)    

    ! initialise the pdd scheme (does this regardless of need at the moment)

    call glint_pdd_init(instance%pddcalc,config)

    ! setting nx,ny of proj

    instance%proj%nx = instance%model%general%ewn
    instance%proj%ny = instance%model%general%nsn
    instance%proj%dx = instance%model%numerics%dew
    instance%proj%dy = instance%model%numerics%dns

    call new_proj(instance%proj,radea)                        ! Initialise the projection
    call new_downscale(instance%downs,instance%proj,grid)     ! Initialise the downscaling
    call glint_i_allocate(instance,grid%nx,grid%ny)           ! Allocate arrays appropriately

    call glint_i_readdata(instance)

    call new_upscale(instance%ups,grid,instance%proj, &
                     instance%climate%out_mask)           ! Initialise upscaling parameters

    call calc_coverage(instance%proj, &                         ! Calculate coverage map
                       instance%ups,  &             
                       grid,          &
                       radea,         &
                       instance%climate%out_mask, &
                       instance%frac_coverage)

    call copy_upscale(instance%ups,instance%ups_orog)    ! Set upscaling for orog output to same as for 
                                                         ! other fields.
    instance%frac_cov_orog=instance%frac_coverage        ! Set fractional coverage for orog to be same as
                                                         ! for other fields.

    if (present(start_time)) then
      instance%model%numerics%time = start_time       ! Initialise the counter.
    else                                              ! Despite being in the GLIMMER framework,
      instance%model%numerics%time = 0.0              ! each instance has a copy of the counter
    endif                                             ! for simplicity.

  end subroutine glint_i_initialise

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_i_end(instance)

    !*FD Performs tidying-up for an ice model. 

    use glide
    implicit none
    type(glint_instance),  intent(inout) :: instance    !*FD The instance being initialised.

    call glide_finalise(instance%model)

  end subroutine glint_i_end

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_i_readconfig(instance,config)
    !*FD read glint configuration
    use glimmer_config
    use glimmer_log
    implicit none
    type(ConfigSection), pointer         :: config      !*FD structure holding sections of configuration file   
    type(glint_instance),  intent(inout) :: instance    !*FD The instance being initialised.

    type(ConfigSection), pointer :: section
    character(len=100) :: message

    ! GLINT projection parameters
    call proj_readconfig(instance%proj,config)         ! read glint projection configuration
    call proj_printconfig(instance%proj)               ! and print it
    
    call GetSection(config,section,'GLINT climate')
    if (associated(section)) then
       call GetValue(section,'temperature_file',instance%forcdata%forcfile)
       call GetValue(section,'artm_mode',instance%climate%whichartm)
       call GetValue(section,'precip_mode',instance%climate%whichprecip)
       call GetValue(section,'acab_mode',instance%climate%whichacab)
    end if

    ! Print some configuration
    call write_log('GLINT climate')
    call write_log('-------------')
    write(message,*) 'external temperature file ', trim(instance%forcdata%forcfile)
    call write_log(message)
    write(message,*) 'artm_mode   ',instance%climate%whichartm
    call write_log(message)
    write(message,*) 'precip_mode ',instance%climate%whichprecip
    call write_log(message)
    write(message,*) 'acab_mode   ',instance%climate%whichacab
    call write_log(message)
    call write_log('')
  end subroutine glint_i_readconfig

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_i_readdata(instance)
    !*FD read data from netCDF file and initialise climate
    use glimmer_log
    use glimmer_ncdf
    use glint_climate
    use glint_io
    use glide_setup
    use glide_temp
    implicit none

    type(glint_instance),intent(inout)   :: instance    !*FD Instance whose elements are to be allocated.
    
    real(sp),dimension(:,:),allocatable :: arng
    type(glimmer_nc_input), pointer :: ic
    logical found_precip,found_presurf,found_usurf
    
    ! read data
    call glint_io_readall(instance,instance%model)

    if (instance%forcdata%forcfile .ne. '') then
      call redtsout(instance%forcdata%forcfile,50,instance%forcdata)
    else
      allocate(instance%forcdata%forcing(1,2))
      instance%forcdata%forcing(1,:) = (/ 2 * instance%forcdata%trun, 0.0 /)
      instance%forcdata%flines = 1
    end if

    ! -------------------------------------------------------------------
    ! Read present-day surface elevation, if required, and calculate
    ! the present-day surface temperature from it.
    ! -------------------------------------------------------------------

    if (found_presurf) then

       ! Allocate arng array, passed to calcartm for air temperature range
       
       allocate(arng(size(instance%climate%presartm,1),size(instance%climate%presartm,2)))

       !----------------------------------------------------------------------
       ! Calculate the present-day mean air temperature and range, based on
       ! surface elevation and latitude
       !----------------------------------------------------------------------
       
       call calcartm(instance, 3, &
            instance%climate%presusrf, &
            instance%model%climate%lati,     &
            instance%climate%presartm, &  !** OUTPUT
            arng)                      !** OUTPUT
       
       !----------------------------------------------------------------------
       ! Calculate present-day mass-balance based on present-day elevation,
       ! temperature, temperature range, and PDD method.
       !----------------------------------------------------------------------
       
       call calcacab(instance%model%numerics,         &
            instance%climate,         &
            instance%pddcalc,          &
            1,                      &
            instance%model%geometry%usrf,  &
            instance%climate%presartm, &
            arng,                   &
            instance%climate%presprcp, &
            instance%climate%ablt,     &
            instance%model%climate%lati,     &
            instance%model%climate%acab,     &
            instance%model%geometry%thck)     
       
       ! Set ice thickness to be mass-balance*time-step, where positive
       
       instance%model%geometry%thck = max(0.0d0,instance%model%climate%acab*instance%model%numerics%dt)
       
       ! Calculate the elevation of the lower ice surface
       
       call glide_calclsrf(instance%model%geometry%thck,instance%model%geometry%topg, &
            instance%model%climate%eus,instance%model%geometry%lsrf)
       
       ! Calculate the elevation of the upper ice surface by adding thickness
       ! onto the lower surface elevation.
       
       instance%model%geometry%usrf = instance%model%geometry%thck + instance%model%geometry%lsrf
       
       call timeevoltemp(instance%model,0)     ! calculate initial temperature distribution

       deallocate(arng) 
       
    else    
       
       ! -----------------------------------------------------------------
       ! Calculate the lower and upper surfaces of the ice-sheet 
       ! -----------------------------------------------------------------
       
       call glide_calclsrf(instance%model%geometry%thck,instance%model%geometry%topg, &
            instance%model%climate%eus,instance%model%geometry%lsrf)
       instance%model%geometry%usrf = instance%model%geometry%thck + instance%model%geometry%lsrf
       
    endif

  end subroutine glint_i_readdata

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine redtsout(fname,unit,forcdata)

    !*FD Sets up forcing data from time-series file (e.g. GRIP data)
    implicit none

    character(*), intent(in)          :: fname     !*FD filename to use
    integer,      intent(in)          :: unit      !*FD File unit to use
    type(glint_forcdata),intent(inout) :: forcdata  !*FD Parameters to be set

    ! Internal variables

    integer :: count, ios = 0
    logical :: there

    ! ---------------------------------------------------------------
    ! Check to see whether file exists
    ! ---------------------------------------------------------------

    inquire(file=fname,exist=there)    

    if ( .not. there ) then
       print*,'ERROR: Time series file not found'
       stop
    endif

    ! ---------------------------------------------------------------
    ! Read in the whole file so we know how many lines there are
    ! ---------------------------------------------------------------

    open(unit,file=fname,form='formatted')

    forcdata%flines = 0

    do while (ios == 0)
       forcdata%flines = forcdata%flines + 1
       read(unit,*,iostat=ios)  
    end do

    forcdata%flines = forcdata%flines - 1

    ! ---------------------------------------------------------------
    ! Allocate array appropriately, then read in data
    ! ---------------------------------------------------------------

    allocate(forcdata%forcing(forcdata%flines,2))

    rewind(unit)

    do count = 1, forcdata%flines
       read(unit,*) forcdata%forcing(count,:)
    end do

    close(unit)

  end subroutine redtsout

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine calc_coverage(proj,ups,grid,radea,mask,frac_coverage)

    !*FD Calculates the fractional
    !*FD coverage of the global grid-boxes by the ice model
    !*FD domain.

    use glint_proj
    use glint_global_grid
    use gmt, only: D2R

    ! Arguments

    type(projection),       intent(in)  :: proj          !*FD Projection to be used
    type(upscale),          intent(in)  :: ups           !*FD Upscaling used
    type(global_grid),      intent(in)  :: grid          !*FD Global grid used
    real(rk),               intent(in)  :: radea         !*FD Radius of the earth (m)
    integer, dimension(:,:),intent(in)  :: mask          !*FD Mask of points for upscaling
    real(rk),dimension(:,:),intent(out) :: frac_coverage !*FD Map of fractional 
                                                         !*FD coverage of global by local grid-boxes.
    ! Internal variables

    integer,dimension(grid%nx,grid%ny) :: tempcount
    integer :: i,j

    ! Beginning of code

    tempcount=0

    do i=1,proj%nx
      do j=1,proj%ny
        tempcount(ups%gboxx(i,j),ups%gboxy(i,j))=tempcount(ups%gboxx(i,j),ups%gboxy(i,j))+mask(i,j)
      enddo
    enddo

    do i=1,grid%nx
      do j=1,grid%ny
        if (tempcount(i,j)==0) then
          frac_coverage(i,j)=0.0
        else
          frac_coverage(i,j)=(tempcount(i,j)*proj%dx*proj%dy)/ &
                 (lon_diff(grid%lon_bound(i+1),grid%lon_bound(i))*D2R*radea**2*    &
                 (sin(grid%lat_bound(j)*D2R)-sin(grid%lat_bound(j+1)*D2R)))
        endif
      enddo
    enddo

    ! Fix points that should be 1.0 by checking their surroundings

    ! Interior points first

    do i=2,grid%nx-1
      do j=2,grid%ny-1
        if ((frac_coverage(i,j).ne.0).and. &
            (frac_coverage(i+1,j).ne.0).and. &
            (frac_coverage(i,j+1).ne.0).and. &
            (frac_coverage(i-1,j).ne.0).and. &
            (frac_coverage(i,j-1).ne.0)) &
                        frac_coverage(i,j)=1.0
      enddo
    enddo

    ! top and bottom edges

    do i=2,grid%nx/2
      if ((frac_coverage(i,1).ne.0).and. &
          (frac_coverage(i+1,1).ne.0).and. &
          (frac_coverage(i,2).ne.0).and. &
          (frac_coverage(i-1,1).ne.0).and. &
          (frac_coverage(i+grid%nx/2,1).ne.0)) &
                      frac_coverage(i,1)=1.0
    enddo

    do i=grid%nx/2+1,grid%nx-1
      if ((frac_coverage(i,1).ne.0).and. &
          (frac_coverage(i+1,1).ne.0).and. &
          (frac_coverage(i,2).ne.0).and. &
          (frac_coverage(i-1,1).ne.0).and. &
          (frac_coverage(i-grid%nx/2,1).ne.0)) &
                      frac_coverage(i,1)=1.0
    enddo

    do i=2,grid%nx/2
      if ((frac_coverage(i,grid%ny).ne.0).and. &
          (frac_coverage(i+1,grid%ny).ne.0).and. &
          (frac_coverage(i+grid%nx/2,grid%ny).ne.0).and. &
          (frac_coverage(i-1,grid%ny).ne.0).and. &
          (frac_coverage(i,grid%ny-1).ne.0)) &
                      frac_coverage(i,grid%ny)=1.0
    enddo

    do i=grid%nx/2+1,grid%nx-1
      if ((frac_coverage(i,grid%ny).ne.0).and. &
          (frac_coverage(i+1,grid%ny).ne.0).and. &
          (frac_coverage(i-grid%nx/2,grid%ny).ne.0).and. &
          (frac_coverage(i-1,grid%ny).ne.0).and. &
          (frac_coverage(i,grid%ny-1).ne.0)) &
                      frac_coverage(i,grid%ny)=1.0
    enddo
 
    ! left and right edges

    do j=2,grid%ny-1
      if ((frac_coverage(1,j).ne.0).and. &
          (frac_coverage(2,j).ne.0).and. &
          (frac_coverage(1,j+1).ne.0).and. &
          (frac_coverage(grid%nx,j).ne.0).and. &
          (frac_coverage(1,j-1).ne.0)) &
                      frac_coverage(1,j)=1.0
      if ((frac_coverage(grid%nx,j).ne.0).and. &
          (frac_coverage(1,j).ne.0).and. &
          (frac_coverage(grid%nx,j+1).ne.0).and. &
          (frac_coverage(grid%nx-1,j).ne.0).and. &
          (frac_coverage(grid%nx,j-1).ne.0)) &
                      frac_coverage(grid%nx,j)=1.0
    enddo

    ! corners

    if ((frac_coverage(1,1).ne.0).and. &
        (frac_coverage(2,1).ne.0).and. &
        (frac_coverage(1,2).ne.0).and. &
        (frac_coverage(grid%nx,1).ne.0).and. &
        (frac_coverage(grid%nx/2+1,1).ne.0)) &
                    frac_coverage(1,1)=1.0

    if ((frac_coverage(1,grid%ny).ne.0).and. &
        (frac_coverage(2,grid%ny).ne.0).and. &
        (frac_coverage(grid%nx/2+1,grid%ny).ne.0).and. &
        (frac_coverage(grid%nx,grid%ny).ne.0).and. &
        (frac_coverage(1,grid%ny-1).ne.0)) &
                    frac_coverage(1,grid%ny)=1.0

    if ((frac_coverage(grid%nx,1).ne.0).and. &
        (frac_coverage(1,1).ne.0).and. &
        (frac_coverage(grid%nx,2).ne.0).and. &
        (frac_coverage(grid%nx-1,1).ne.0).and. &
        (frac_coverage(grid%nx/2,1).ne.0)) &
                   frac_coverage(grid%nx,1)=1.0

    if ((frac_coverage(grid%nx,grid%ny).ne.0).and. &
        (frac_coverage(1,grid%ny).ne.0).and. &
        (frac_coverage(grid%nx/2,grid%ny).ne.0).and. &
        (frac_coverage(grid%nx-1,grid%ny).ne.0).and. &
        (frac_coverage(grid%nx,grid%ny-1).ne.0)) &
                   frac_coverage(grid%nx,grid%ny)=1.0

    ! Finally fix any rogue points > 1.0

    where (frac_coverage>1.0) frac_coverage=1.0

  end subroutine calc_coverage

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(rk) function lon_diff(a,b)

    implicit none

    real(rk),intent(in) :: a,b
    real(rk) :: aa,bb

    aa=a ; bb=b

      do
        if (aa>bb) exit
        aa=aa+360.0
      enddo

    lon_diff=aa-bb

  end function lon_diff

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_i_upscaled_fields(instance,orog,albedo,ice_frac)

    !*FD Upscales and returns certain fields, according to the 
    !*FD arguments supplied

    use paramets

    type(glint_instance),            intent(in)  :: instance
    real(rk),dimension(:,:),optional,intent(out) :: orog
    real(rk),dimension(:,:),optional,intent(out) :: albedo
    real(rk),dimension(:,:),optional,intent(out) :: ice_frac

    real(rk),dimension(:,:),allocatable :: if_temp,upscale_temp

	  ! Calculate orography

    if (present(orog)) then
      call mean_to_global(instance%proj, &
                          instance%ups_orog, &
                          instance%model%geometry%usrf, &
                          orog,    &
                          instance%climate%out_mask)
      orog=thk0*orog
    endif

    if (present(albedo).or.present(ice_frac)) then

      if (present(albedo)) then
        allocate(if_temp(size(albedo,1),size(albedo,2)))
      else
        allocate(if_temp(size(ice_frac,1),size(ice_frac,2)))
      endif
      allocate(upscale_temp(instance%model%general%ewn,instance%model%general%nsn))

      ! First, ice coverage on local grid 
  
      where (instance%model%geometry%thck>0.0)
        upscale_temp=1.0
      elsewhere
        upscale_temp=0.0
      endwhere

      ! Upscale it...

      call mean_to_global(instance%proj, &
                          instance%ups, &
                          upscale_temp, &
                          if_temp,    &
                          instance%climate%out_mask)

      if (present(ice_frac)) ice_frac=if_temp

    endif

    ! Calculate albedo -------------------------------------------------------

    if (present(albedo)) then 
      where (if_temp>0.0)
        albedo=instance%climate%ice_albedo
      elsewhere
        albedo=0.0
      endwhere
    endif

    if (allocated(if_temp)) deallocate(if_temp)
    if (allocated(upscale_temp)) deallocate(upscale_temp)

  end subroutine get_i_upscaled_fields

end module glint_initialise
