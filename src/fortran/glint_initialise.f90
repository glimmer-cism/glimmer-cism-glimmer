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

  subroutine glint_i_initialise(config,instance,grid,grid_orog,mbts,need_winds)

    !*FD Initialise a GLINT ice model instance

    use glimmer_config
    use glint_global_grid
    use glint_io
    use glint_mbal_io
    use glide
    use glint_constants
    implicit none

    ! Arguments
    type(ConfigSection), pointer         :: config      !*FD structure holding sections of configuration file   
    type(glint_instance),  intent(inout) :: instance    !*FD The instance being initialised.
    type(global_grid),     intent(in)    :: grid        !*FD Global grid to use
    type(global_grid),     intent(in)    :: grid_orog   !*FD Global grid to use for orography
    integer,               intent(out)   :: mbts        !*FD mass-balance time-step (hours)
    logical,               intent(inout) :: need_winds  !*FD Set if this instance needs wind input

    ! initialise model

    call glide_initialise(instance%model,config)
    instance%ice_tstep=get_tinc(instance%model)*years2hours

    ! create glint variables
    call glint_io_createall(instance%model)
    call glint_mbal_io_createall(instance%model)

    ! fill dimension variables
    call glide_nc_fillall(instance%model)

    ! read glint configuration

    call glint_i_readconfig(instance,config)    
    call glint_i_printconfig(instance)    

    ! Check we've used all the config sections

    call CheckSections(config)

    ! New projection and downscaling

    call new_proj(instance%proj, &                       ! Initialise the projection
                  nx=get_ewn(instance%model), &
                  ny=get_nsn(instance%model), &
                  dx=get_dew(instance%model), &
                  dy=get_dns(instance%model))
    call new_downscale(instance%downs,instance%proj,grid)     ! Initialise the downscaling

    call glint_i_allocate(instance,grid%nx,grid%ny,grid_orog%nx,grid_orog%ny)           ! Allocate arrays appropriately
    call glint_i_readdata(instance)
    call new_upscale(instance%ups,grid,instance%proj,instance%out_mask) ! Initialise upscaling parameters
    call new_upscale(instance%ups_orog,grid_orog,instance%proj,instance%out_mask) ! Initialise upscaling parameters

    call calc_coverage(instance%proj, &                         ! Calculate coverage map
                       instance%ups,  &             
                       grid,          &
                       instance%out_mask, &
                       instance%frac_coverage)

    call calc_coverage(instance%proj, &                         ! Calculate coverage map for orog
                       instance%ups_orog,  &             
                       grid_orog,     &
                       instance%out_mask, &
                       instance%frac_cov_orog)

    ! initialise the mass-balance accumulation

    call glint_mbc_init(instance%mbal_accum,instance%proj,config,instance%whichacab)
    instance%mbal_tstep=instance%mbal_accum%mbal%tstep
    mbts=instance%mbal_tstep

    call glint_mbal_io_writeall(instance%mbal_accum,instance%model)
    call glide_io_writeall(instance%model,instance%model)
    call glint_io_writeall(instance,instance%model)

    if (instance%whichprecip==2) need_winds=.true.

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
     
    ! read data
    call glint_io_readall(instance,instance%model)
      
    call glide_calclsrf(instance%model%geometry%thck,instance%model%geometry%topg, &
         instance%model%climate%eus,instance%model%geometry%lsrf)
    instance%model%geometry%usrf = instance%model%geometry%thck + instance%model%geometry%lsrf
 
  end subroutine glint_i_readdata

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine calc_coverage(proj,ups,grid,mask,frac_coverage)

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
                 (lon_diff(grid%lon_bound(i+1),grid%lon_bound(i))*D2R*proj%radea**2*    &
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
      call mean_to_global(instance%ups_orog, &
                          instance%model%geometry%usrf, &
                          orog,    &
                          instance%out_mask)
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

      call mean_to_global(instance%ups, &
                          upscale_temp, &
                          if_temp,    &
                          instance%out_mask)

      if (present(ice_frac)) ice_frac=if_temp

    endif

    ! Calculate albedo -------------------------------------------------------

    if (present(albedo)) then 
      where (if_temp>0.0)
        albedo=instance%ice_albedo
      elsewhere
        albedo=0.0
      endwhere
    endif

    if (allocated(if_temp)) deallocate(if_temp)
    if (allocated(upscale_temp)) deallocate(upscale_temp)

  end subroutine get_i_upscaled_fields

end module glint_initialise
