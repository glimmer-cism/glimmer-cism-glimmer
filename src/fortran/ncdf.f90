!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WARNING: this file was automatically generated on
! Wed, 12 Jan 2005 11:34:19 +0000
! from ncdf.f90.in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  ncdf.f90 - part of the GLIMMER ice model                 + 
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

module glimmer_ncdf
  !*FD Data structures and utility functions for netCDF I/O
  !*FD written by Magnus Hagdorn, 2004

  use glimmer_global, only: fname_length
  use netcdf

  integer, private, parameter :: num_vars = 68
  integer, private, parameter :: meta_len = 100

  integer, parameter :: NC_B_ABLT = 1 ! ablation
  integer, parameter :: NC_B_ABLT_SPOT = 2 ! ablation
  integer, parameter :: NC_B_ACAB = 3 ! accumulation, ablation rate
  integer, parameter :: NC_B_ACAB_SPOT = 4 ! accumulation, ablation rate
  integer, parameter :: NC_B_ARNG = 5 ! annual temperature range
  integer, parameter :: NC_B_ARNG_SPOT = 6 ! annual temperature range
  integer, parameter :: NC_B_ARTM = 7 ! annual mean air temperature
  integer, parameter :: NC_B_ARTM_SPOT = 8 ! annual mean air temperature
  integer, parameter :: NC_B_BMLT = 9 ! basal melt rate
  integer, parameter :: NC_B_BMLT_SPOT = 10 ! basal melt rate
  integer, parameter :: NC_B_BTEMP = 11 ! basal ice temperature
  integer, parameter :: NC_B_BTEMP_SPOT = 12 ! basal ice temperature
  integer, parameter :: NC_B_BTRC = 13 ! basal slip coefficient
  integer, parameter :: NC_B_BTRC_SPOT = 14 ! basal slip coefficient
  integer, parameter :: NC_B_BWAT = 15 ! basal water depth
  integer, parameter :: NC_B_BWAT_SPOT = 16 ! basal water depth
  integer, parameter :: NC_B_DIFFU = 17 ! apparent diffusivity
  integer, parameter :: NC_B_DIFFU_SPOT = 18 ! apparent diffusivity
  integer, parameter :: NC_B_DUSRFDTM = 19 ! rate of upper ice surface elevation change
  integer, parameter :: NC_B_DUSRFDTM_SPOT = 20 ! rate of upper ice surface elevation change
  integer, parameter :: NC_B_FLWA = 21 ! ??
  integer, parameter :: NC_B_FLWA_SPOT = 22 ! ??
  integer, parameter :: NC_B_LAT = 23 ! Latitude
  integer, parameter :: NC_B_LAT_SPOT = 24 ! Latitude
  integer, parameter :: NC_B_LON = 25 ! Longitude
  integer, parameter :: NC_B_LON_SPOT = 26 ! Longitude
  integer, parameter :: NC_B_LSURF = 27 ! ice lower surface elevation
  integer, parameter :: NC_B_LSURF_SPOT = 28 ! ice lower surface elevation
  integer, parameter :: NC_B_MASK = 29 ! upscaling and downscaling mask
  integer, parameter :: NC_B_MASK_SPOT = 30 ! upscaling and downscaling mask
  integer, parameter :: NC_B_PRCP = 31 ! precipitation
  integer, parameter :: NC_B_PRCP_SPOT = 32 ! precipitation
  integer, parameter :: NC_B_PRESPRCP = 33 ! present day precipitation
  integer, parameter :: NC_B_PRESPRCP_SPOT = 34 ! present day precipitation
  integer, parameter :: NC_B_PRESUSRF = 35 ! present day surface of the ice-sheet
  integer, parameter :: NC_B_PRESUSRF_SPOT = 36 ! present day surface of the ice-sheet
  integer, parameter :: NC_B_RELX = 37 ! relaxed bedrock topography
  integer, parameter :: NC_B_RELX_SPOT = 38 ! relaxed bedrock topography
  integer, parameter :: NC_B_STD_DEV = 39 ! standard deviation of sub-grid topography
  integer, parameter :: NC_B_STD_DEV_SPOT = 40 ! standard deviation of sub-grid topography
  integer, parameter :: NC_B_TEMP = 41 ! ice temperature
  integer, parameter :: NC_B_TEMP_SPOT = 42 ! ice temperature
  integer, parameter :: NC_B_THK = 43 ! ice thickness
  integer, parameter :: NC_B_THK_SPOT = 44 ! ice thickness
  integer, parameter :: NC_B_TOPG = 45 ! bedrock topography
  integer, parameter :: NC_B_TOPG_SPOT = 46 ! bedrock topography
  integer, parameter :: NC_B_UBAS = 47 ! basal slip velocity in x direction
  integer, parameter :: NC_B_UBAS_SPOT = 48 ! basal slip velocity in x direction
  integer, parameter :: NC_B_UFLX = 49 ! flux in x direction
  integer, parameter :: NC_B_UFLX_SPOT = 50 ! flux in x direction
  integer, parameter :: NC_B_USURF = 51 ! ice upper surface elevation
  integer, parameter :: NC_B_USURF_SPOT = 52 ! ice upper surface elevation
  integer, parameter :: NC_B_UVEL = 53 ! ice velocity in x direction
  integer, parameter :: NC_B_UVEL_SPOT = 54 ! ice velocity in x direction
  integer, parameter :: NC_B_VBAS = 55 ! basal slip velocity in y direction
  integer, parameter :: NC_B_VBAS_SPOT = 56 ! basal slip velocity in y direction
  integer, parameter :: NC_B_VFLX = 57 ! flux in x direction
  integer, parameter :: NC_B_VFLX_SPOT = 58 ! flux in x direction
  integer, parameter :: NC_B_VVEL = 59 ! ice velocity in y direction
  integer, parameter :: NC_B_VVEL_SPOT = 60 ! ice velocity in y direction
  integer, parameter :: NC_B_WGRD = 61 ! ?? some velo ??
  integer, parameter :: NC_B_WGRD_SPOT = 62 ! ?? some velo ??
  integer, parameter :: NC_B_WVEL = 63 ! vertical ice velocity
  integer, parameter :: NC_B_WVEL_SPOT = 64 ! vertical ice velocity

  type glimmer_nc_stat
     !*FD Data structure holding netCDF file description

     character(len=fname_length) :: filename = " "
     !*FD name of netCDF file
     logical, dimension(num_vars) :: do_var
     !*FD array specifying which variables should be written to netCDF file
     logical :: do_spot = .false.
     !*FD write spot data
     integer id
     !*FD id of netCDF file

     integer x0dim
     !*FD id of x0 dimension
     integer y0dim
     !*FD id of y0 dimension
     integer x1dim
     !*FD id of x1 dimension
     integer y1dim
     !*FD id of y1 dimension
     integer leveldim
     !*FD id of sigma level dimension
     integer timedim
     !*FD id of time dimension
     integer spotdim
     !*FD id of spot index dimensions
     integer x0var, x0_spotvar
     !*FD id of x0 variable
     integer y0var, y0_spotvar
     !*FD id of y0 variable
     integer x1var, x1_spotvar
     !*FD id of x1 variable
     integer y1var, y1_spotvar
     !*FD id of y1 variable
     integer levelvar
     !*FD id of sigma level variable
     integer timevar
     !*FD id of time variable 
     integer, dimension(num_vars) :: varids
     !*FD array holding variable ids
  end type glimmer_nc_stat

  type glimmer_nc_meta
     !*FD Data structure holding netCDF meta data, see CF user guide
     
     character(len=meta_len) :: title = ''
     !*FD title of netCDF file
     character(len=meta_len) :: institution = ''
     !*FD where the data was produced
     character(len=meta_len) :: references = ''
     !*FD list of references
     character(len=meta_len) :: source = ''
     !*FD this string will hold the GLIMMER version
     character(len=meta_len) :: history = ''
     !*FD netCDF file history string
     character(len=meta_len) :: comment = ''
     !*FD some comments
  end type glimmer_nc_meta

  type glimmer_nc_output
     !*FD element of linked list describing netCDF output file

     type(glimmer_nc_stat) :: nc
     !*FD structure containg file info
     integer :: freq=1000
     !*FD frequency at which data is written to file
     integer :: next_write=0
     !*FD next time step at which data is dumped
     integer :: timecounter=1
     !*FD time counter
     integer, pointer, dimension(:) :: spotx=>NULL()
     !*FD array containg spot x-index
     integer, pointer, dimension(:) :: spoty=>NULL()
     !*FD array containg spot y-index
     
     type(glimmer_nc_meta) :: metadata
     !*FD structure holding metadata

     type(glimmer_nc_output), pointer :: next=>NULL()
     !*FD next element in list
     type(glimmer_nc_output), pointer :: previous=>NULL()
     !*FD previous element in list
  end type glimmer_nc_output

  type glimmer_nc_input
     !*FD element of linked list describing netCDF input file

     type(glimmer_nc_stat) :: nc
     !*FD structure containg file info
     integer, pointer, dimension(:) :: times => NULL()     
     !*FD pointer to array holding times
     integer                        :: nt, current_time=1
     !*FDnumber of elements in times and current time index
     integer                        :: get_time_slice = 1     
     !*FD -1 if all times should be loaded, > 0 to load particular slice and then close file

     type(glimmer_nc_input), pointer :: next=>NULL()
     !*FD next element in list
     type(glimmer_nc_input), pointer :: previous=>NULL()
     !*FD previous element in list
  end type glimmer_nc_input

  interface delete
     module procedure delete_output, delete_input
  end interface

  interface add
     module procedure add_output, add_input
  end interface

contains
  function delete_output(oc, cf)
    !*FD remove element from linked list
	use glide_messages
    implicit none
    type(glimmer_nc_output), pointer :: delete_output
    type(glimmer_nc_output), pointer :: oc
    logical, intent(in), optional :: cf
    ! local variables
    logical closefile
    integer status

    if (present(cf)) then
       closefile = cf
    else
       closefile = .true.
    end if

    if (associated(oc)) then
       if (associated(oc%previous)) then
          oc%previous%next => oc%next
       end if
       if (associated(oc%next)) then
          oc%next%previous => oc%previous
          delete_output => oc%next
       else
          delete_output => NULL()
       end if
       if (closefile) then
          status = nf90_close(oc%nc%id)
          call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Closing output file '//trim(oc%nc%filename))
       end if
       deallocate(oc)
    end if
  end function delete_output
  
  function delete_input(ic,cf)
    !*FD remove element from linked list
	use glide_messages
    implicit none
    type(glimmer_nc_input), pointer :: delete_input
    type(glimmer_nc_input), pointer :: ic
    logical, intent(in), optional :: cf

    ! local variables
    logical closefile
    integer status

    if (present(cf)) then
       closefile = cf
    else
       closefile = .true.
    end if

    if (associated(ic)) then
       if (associated(ic%previous)) then
          ic%previous%next => ic%next
       end if
       if (associated(ic%next)) then
          ic%next%previous => ic%previous
          delete_input => ic%next
       else
          delete_input => NULL()
       end if
       if (closefile) then
          status = nf90_close(ic%nc%id)
          call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Closing input file '//trim(ic%nc%filename))
       end if
       deallocate(ic%times)
       deallocate(ic)
    end if
  end function delete_input

  function add_output(oc)
    !*FD add new element to linked list
    implicit none
    type(glimmer_nc_output), pointer :: add_output
    type(glimmer_nc_output), pointer :: oc

    allocate(add_output)
    add_output%nc%do_var(:) = .false.

    if (associated(oc)) then
       add_output%previous => oc
       if (associated(oc%next)) then
          add_output%next => oc%next
          oc%next%previous => add_output
       end if
       oc%next => add_output
    end if
  end function add_output

  function add_input(ic)
    !*FD add new element to linked list
    implicit none
    type(glimmer_nc_input), pointer :: add_input
    type(glimmer_nc_input), pointer :: ic

    allocate(add_input)
    add_input%nc%do_var(:) = .false.

    if (associated(ic)) then
       add_input%previous => ic
       if (associated(ic%next)) then
          add_input%next => ic%next
          ic%next%previous => add_input
       end if
       ic%next => add_input
    end if
  end function add_input

  subroutine check_vars(nc,unit)
    !*FD print netCDF variables handled by nc to unit
    implicit none
    type(glimmer_nc_stat) :: nc
    !*FD netCDF file descriptor
    integer, intent(in) :: unit
    !*FD file unit to be written to

    if (nc%do_var(NC_B_ABLT)) then
       write(unit,*) 'ablt'
    end if
    if (nc%do_var(NC_B_ABLT_SPOT)) then
       write(unit,*) 'ablt_spot'
    end if
    if (nc%do_var(NC_B_ACAB)) then
       write(unit,*) 'acab'
    end if
    if (nc%do_var(NC_B_ACAB_SPOT)) then
       write(unit,*) 'acab_spot'
    end if
    if (nc%do_var(NC_B_ARNG)) then
       write(unit,*) 'arng'
    end if
    if (nc%do_var(NC_B_ARNG_SPOT)) then
       write(unit,*) 'arng_spot'
    end if
    if (nc%do_var(NC_B_ARTM)) then
       write(unit,*) 'artm'
    end if
    if (nc%do_var(NC_B_ARTM_SPOT)) then
       write(unit,*) 'artm_spot'
    end if
    if (nc%do_var(NC_B_BMLT)) then
       write(unit,*) 'bmlt'
    end if
    if (nc%do_var(NC_B_BMLT_SPOT)) then
       write(unit,*) 'bmlt_spot'
    end if
    if (nc%do_var(NC_B_BTEMP)) then
       write(unit,*) 'btemp'
    end if
    if (nc%do_var(NC_B_BTEMP_SPOT)) then
       write(unit,*) 'btemp_spot'
    end if
    if (nc%do_var(NC_B_BTRC)) then
       write(unit,*) 'btrc'
    end if
    if (nc%do_var(NC_B_BTRC_SPOT)) then
       write(unit,*) 'btrc_spot'
    end if
    if (nc%do_var(NC_B_BWAT)) then
       write(unit,*) 'bwat'
    end if
    if (nc%do_var(NC_B_BWAT_SPOT)) then
       write(unit,*) 'bwat_spot'
    end if
    if (nc%do_var(NC_B_DIFFU)) then
       write(unit,*) 'diffu'
    end if
    if (nc%do_var(NC_B_DIFFU_SPOT)) then
       write(unit,*) 'diffu_spot'
    end if
    if (nc%do_var(NC_B_DUSRFDTM)) then
       write(unit,*) 'dusrfdtm'
    end if
    if (nc%do_var(NC_B_DUSRFDTM_SPOT)) then
       write(unit,*) 'dusrfdtm_spot'
    end if
    if (nc%do_var(NC_B_FLWA)) then
       write(unit,*) 'flwa'
    end if
    if (nc%do_var(NC_B_FLWA_SPOT)) then
       write(unit,*) 'flwa_spot'
    end if
    if (nc%do_var(NC_B_LAT)) then
       write(unit,*) 'lat'
    end if
    if (nc%do_var(NC_B_LAT_SPOT)) then
       write(unit,*) 'lat_spot'
    end if
    if (nc%do_var(NC_B_LON)) then
       write(unit,*) 'lon'
    end if
    if (nc%do_var(NC_B_LON_SPOT)) then
       write(unit,*) 'lon_spot'
    end if
    if (nc%do_var(NC_B_LSURF)) then
       write(unit,*) 'lsurf'
    end if
    if (nc%do_var(NC_B_LSURF_SPOT)) then
       write(unit,*) 'lsurf_spot'
    end if
    if (nc%do_var(NC_B_MASK)) then
       write(unit,*) 'mask'
    end if
    if (nc%do_var(NC_B_MASK_SPOT)) then
       write(unit,*) 'mask_spot'
    end if
    if (nc%do_var(NC_B_PRCP)) then
       write(unit,*) 'prcp'
    end if
    if (nc%do_var(NC_B_PRCP_SPOT)) then
       write(unit,*) 'prcp_spot'
    end if
    if (nc%do_var(NC_B_PRESPRCP)) then
       write(unit,*) 'presprcp'
    end if
    if (nc%do_var(NC_B_PRESPRCP_SPOT)) then
       write(unit,*) 'presprcp_spot'
    end if
    if (nc%do_var(NC_B_PRESUSRF)) then
       write(unit,*) 'presusrf'
    end if
    if (nc%do_var(NC_B_PRESUSRF_SPOT)) then
       write(unit,*) 'presusrf_spot'
    end if
    if (nc%do_var(NC_B_RELX)) then
       write(unit,*) 'relx'
    end if
    if (nc%do_var(NC_B_RELX_SPOT)) then
       write(unit,*) 'relx_spot'
    end if
    if (nc%do_var(NC_B_STD_DEV)) then
       write(unit,*) 'std_dev'
    end if
    if (nc%do_var(NC_B_STD_DEV_SPOT)) then
       write(unit,*) 'std_dev_spot'
    end if
    if (nc%do_var(NC_B_TEMP)) then
       write(unit,*) 'temp'
    end if
    if (nc%do_var(NC_B_TEMP_SPOT)) then
       write(unit,*) 'temp_spot'
    end if
    if (nc%do_var(NC_B_THK)) then
       write(unit,*) 'thk'
    end if
    if (nc%do_var(NC_B_THK_SPOT)) then
       write(unit,*) 'thk_spot'
    end if
    if (nc%do_var(NC_B_TOPG)) then
       write(unit,*) 'topg'
    end if
    if (nc%do_var(NC_B_TOPG_SPOT)) then
       write(unit,*) 'topg_spot'
    end if
    if (nc%do_var(NC_B_UBAS)) then
       write(unit,*) 'ubas'
    end if
    if (nc%do_var(NC_B_UBAS_SPOT)) then
       write(unit,*) 'ubas_spot'
    end if
    if (nc%do_var(NC_B_UFLX)) then
       write(unit,*) 'uflx'
    end if
    if (nc%do_var(NC_B_UFLX_SPOT)) then
       write(unit,*) 'uflx_spot'
    end if
    if (nc%do_var(NC_B_USURF)) then
       write(unit,*) 'usurf'
    end if
    if (nc%do_var(NC_B_USURF_SPOT)) then
       write(unit,*) 'usurf_spot'
    end if
    if (nc%do_var(NC_B_UVEL)) then
       write(unit,*) 'uvel'
    end if
    if (nc%do_var(NC_B_UVEL_SPOT)) then
       write(unit,*) 'uvel_spot'
    end if
    if (nc%do_var(NC_B_VBAS)) then
       write(unit,*) 'vbas'
    end if
    if (nc%do_var(NC_B_VBAS_SPOT)) then
       write(unit,*) 'vbas_spot'
    end if
    if (nc%do_var(NC_B_VFLX)) then
       write(unit,*) 'vflx'
    end if
    if (nc%do_var(NC_B_VFLX_SPOT)) then
       write(unit,*) 'vflx_spot'
    end if
    if (nc%do_var(NC_B_VVEL)) then
       write(unit,*) 'vvel'
    end if
    if (nc%do_var(NC_B_VVEL_SPOT)) then
       write(unit,*) 'vvel_spot'
    end if
    if (nc%do_var(NC_B_WGRD)) then
       write(unit,*) 'wgrd'
    end if
    if (nc%do_var(NC_B_WGRD_SPOT)) then
       write(unit,*) 'wgrd_spot'
    end if
    if (nc%do_var(NC_B_WVEL)) then
       write(unit,*) 'wvel'
    end if
    if (nc%do_var(NC_B_WVEL_SPOT)) then
       write(unit,*) 'wvel_spot'
    end if
    
  end subroutine check_vars

end module glimmer_ncdf

module glimmer_scales
  !*FD this module holds scales for various fields

  use glimmer_global, only : dp

  real(dp) :: scale2d_f1, scale2d_f2, scale2d_f3, scale2d_f4, scale2d_f5, scale2d_f6, scale2d_f7, scale2d_f8
  real(dp) :: scale3d_f1, scale3d_f2, scale3d_f3, scale3d_f4, scale3d_f5, scale3d_f6, scale3d_f7, scale3d_f8

contains
  subroutine glimmer_init_scales
    !*FD calculate scale factors (can't have non-integer powers)
    use physcon, only : scyr, gn
    use paramets, only : thk0, tim0, vel0, vis0, len0, tau0
    implicit none

    scale2d_f1 = scyr * thk0 / tim0
    scale2d_f2 = scyr * vel0 * thk0
    scale2d_f3 = vel0 / (vis0 * len0)
    scale2d_f4 = vel0 * scyr * len0
    scale2d_f5 = scyr * vel0
    scale2d_f6 = scyr * vel0 * len0 / (thk0**2)
    scale2d_f7 = tau0
    scale2d_f8 = tau0 * len0 / (scyr * vel0)

    scale3d_f1 = scyr * vel0
    scale3d_f2 = vis0 * (vel0/len0)**(gn - 1)
    scale3d_f3 = scyr * thk0
    scale3d_f4 = vel0/(vis0*len0)
    scale3d_f5 = 1.0d0/scale3d_f2**(1.0/gn)
    scale3d_f6 = scale3d_f4**(1.0/gn)
    scale3d_f7 = scyr * thk0/tim0
    scale3d_f8 = vis0*scyr
  end subroutine glimmer_init_scales
end module glimmer_scales

subroutine nc_errorhandle(file,line,status)
  !*FD handle netCDF error
  use netcdf
  implicit none
  character(len=*), intent(in) :: file
  !*FD name of f90 file error occured in
  integer, intent(in) :: line
  !*FD line number error occured at
  integer, intent(in) :: status
  !*FD netCDF return value
  
  if (status.ne.NF90_NOERR) then
     write(*,*) 'NETCDF Error (',file,line,'): ', nf90_strerror(status)
     stop
  end if
end subroutine nc_errorhandle
