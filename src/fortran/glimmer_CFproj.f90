! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_CFproj.f90 - part of the GLIMMER ice model       + 
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

module glimmer_CFproj

  !*FD Holds derived types and subroutines
  !*FD necessary for handling map projections. Most of the component
  !*FD names of the various derived types are self-explanatory.
  !*FD Note that this doesn't currently interface with the proj4
  !*FD library in anyway, it simply handles NetCDF data and projection
  !*FD parameters in an appropriate format.

  private
  public  CFproj_projection,CFproj_proj4,CFproj_GetProj,CFproj_PutProj,CFproj_define
  public  CFP_LAEA,CFP_AEA,CFP_LCC,CFP_STERE

  integer, parameter :: proj4len=100
  integer, parameter :: CFP_LAEA=1
  integer, parameter :: CFP_AEA=2
  integer, parameter :: CFP_LCC=3
  integer, parameter :: CFP_STERE=4

  ! Type definitions for various projections -------------------

  type CFproj_projection

     !*FD CF projection type. This type contains a pointer to an
     !*FD instance of each type of projection, all initialised to
     !*FD \texttt{NULL}.

     private 
     type(CFproj_laea),  pointer :: laea  => NULL() !*FD Pointer to Lambert azimuthal equal area type
     type(CFproj_aea),   pointer :: aea   => NULL() !*FD Pointer to Albers equal area conic type
     type(CFproj_lcc),   pointer :: lcc   => NULL() !*FD Pointer to Lambert conic conformal type
     type(CFproj_stere), pointer :: stere => NULL() !*FD Pointer to Stereographic type

  end type CFproj_projection

  !-------------------------------------------------------------

  type CFproj_stere
     !*FD Stereographic projection derived type 
     logical :: polar                               !*FD Polar projection?
     real :: longitude_of_central_meridian          
     real :: latitude_of_projection_origin          
     real :: scale_factor_at_projection_origin = 0. 
     real :: standard_parallel = 0.                 
     real :: false_easting          
     real :: false_northing
  end type CFproj_stere

  !-------------------------------------------------------------

  type CFproj_laea
     !*FD Lambert Azimuthal Equal Area
     real :: longitude_of_central_meridian
     real :: latitude_of_projection_origin
     real :: false_easting
     real :: false_northing
  end type CFproj_laea

  !-------------------------------------------------------------

  type CFproj_aea
     !*FD Albers Equal-Area Conic
     real,dimension(2) :: standard_parallel
     real :: longitude_of_central_meridian
     real :: latitude_of_projection_origin
     real :: false_easting
     real :: false_northing
  end type CFproj_aea

  !-------------------------------------------------------------

  type CFproj_lcc
     !*FD Lambert Conic Conformal
     real,dimension(2) :: standard_parallel
     real :: longitude_of_central_meridian
     real :: latitude_of_projection_origin
     real :: false_easting
     real :: false_northing
  end type CFproj_lcc

contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! public functions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function CFproj_GetProj(ncid)

    !*FD Read projection from a given netCDF file, returning
    !*FD an instance of type \texttt{CFproj\_projection}.
    !*FDRV Derived type instance containing projection parameters

    use netcdf
    use glide_messages

    implicit none

    type(CFproj_projection) :: CFproj_GetProj
    integer, intent(in) :: ncid                !*FD Handle of the file to be read.
    
    !local variables
    integer status
    integer nvars, varid
    integer natts, attid
    logical found_map
    character(len=50) :: attname,mapname

    ! getting variables
    status = nf90_inquire(ncid,nvariables=nvars)
    call nc_errorhandle(__FILE__,__LINE__,status)    
    
    ! looping over variables
    found_map=.false.
    do varid=1,nvars
       status = nf90_inquire_variable(ncid,varid,natts=natts)
       ! and loop over attributes
       do attid=1,natts
          status = nf90_inq_attname(ncid,varid,attid,attname)
          if (trim(attname) == 'grid_mapping_name') then
             found_map = .true.
             status = nf90_get_att(ncid,varid,attname,mapname)
             mapname = adjustl(mapname)
             call nc_errorhandle(__FILE__,__LINE__,status)
             exit
          end if
       end do
       if (found_map) exit
    end do

    if (found_map) then
       if (index(mapname,'lambert_azimuthal_equal_area').ne.0) then
          CFproj_GetProj%laea => CFproj_get_laea(ncid,varid)
          return
       else if (index(mapname,'albers_conical_equal_area').ne.0) then
          CFproj_GetProj%aea => CFproj_get_aea(ncid,varid)
          return
       else if (index(mapname,'lambert_conformal_conic').ne.0) then
          CFproj_GetProj%lcc => CFproj_get_lcc(ncid,varid)
       else if (index(mapname,'polar_stereographic').ne.0) then
          CFproj_GetProj%stere => CFproj_get_stere_polar(ncid,varid)
       else if (index(mapname,'stereographic').ne.0) then
          CFproj_GetProj%stere => CFproj_get_stere(ncid,varid)
       else
          call glide_msg(GM_WARNING,__FILE__,__LINE__,'Do not know about this projection: '//trim(mapname))
       end if
    else
       call glide_msg(GM_WARNING,__FILE__,__LINE__,'No map projection found')
    end if
  end function CFproj_GetProj

  !-------------------------------------------------------------------------

  function CFproj_proj4(projection)
    
    !*FD Returns a proj4 parameter string for a given set of projection parameters
    !*FDRV Pointer to array of projection parameter strings

    use glide_messages
    implicit none
    character(len=proj4len), dimension(:), pointer :: CFproj_proj4
    type(CFproj_projection) :: projection !*FD Projection of interest

    if (associated(projection%laea)) then
       CFproj_proj4 => CFproj_proj4_laea(projection%laea)
       return
    else if (associated(projection%aea)) then
       CFproj_proj4 => CFproj_proj4_aea(projection%aea)
       return
    else if (associated(projection%lcc)) then
       CFproj_proj4 => CFproj_proj4_lcc(projection%lcc)
       return
    else if (associated(projection%stere)) then
       CFproj_proj4 => CFproj_proj4_stere(projection%stere)
       return
    else
       call glide_msg(GM_WARNING,__FILE__,__LINE__,'No projection found!')
    end if
  end function CFproj_proj4

  !-------------------------------------------------------------------------

  subroutine CFproj_PutProj(ncid,mapid,projection)

    !*FD write projection to a netCDF file.

    use netcdf
    use glide_messages

    implicit none

    type(CFproj_projection) :: projection   !*FD Projection to be written.
    integer, intent(in) :: ncid             !*FD Handle of netCDF file.
    integer, intent(in) :: mapid            !*FD Handle of map projection in netCDF file.

    if (associated(projection%laea)) then
       call CFproj_put_laea(ncid,mapid,projection%laea)
       return
    else if (associated(projection%aea)) then
       call CFproj_put_aea(ncid,mapid,projection%aea)
       return
    else if (associated(projection%lcc)) then
       call CFproj_put_lcc(ncid,mapid,projection%lcc)
       return
    else if (associated(projection%stere)) then
       call CFproj_put_stere(ncid,mapid,projection%stere)
       return
    else
       call glide_msg(GM_FATAL,__FILE__,__LINE__,'Could not find any projection')
    end if
  end subroutine CFproj_PutProj

  !-------------------------------------------------------------------------

  subroutine CFproj_define(cfp,ptype,polar, &
       longitude_of_central_meridian, &
       latitude_of_projection_origin, &
       scale_factor_at_projection_origin, &
       standard_parallel, &
       standard_parallel_2, &
       false_easting, &
       false_northing)

    use glide_messages

    type(CFproj_projection),intent(inout) :: cfp 
    integer,intent(in) :: ptype
    logical,optional,intent(in) :: polar
    real,optional,intent(in) :: longitude_of_central_meridian
    real,optional,intent(in) :: latitude_of_projection_origin
    real,optional,intent(in) :: scale_factor_at_projection_origin
    real,optional,intent(in) :: standard_parallel
    real,optional,intent(in) :: standard_parallel_2
    real,optional,intent(in) :: false_easting
    real,optional,intent(in) :: false_northing

    if (associated(cfp%laea))  deallocate(cfp%laea)
    if (associated(cfp%aea))   deallocate(cfp%aea)
    if (associated(cfp%lcc))   deallocate(cfp%lcc)
    if (associated(cfp%stere)) deallocate(cfp%stere)

    select case(ptype)
    case(CFP_LAEA)
       allocate(cfp%laea)
       if(present(longitude_of_central_meridian)) &
            cfp%laea%longitude_of_central_meridian = longitude_of_central_meridian
       if(present(latitude_of_projection_origin)) &
            cfp%laea%latitude_of_projection_origin = latitude_of_projection_origin
       if(present(false_easting)) &
            cfp%laea%false_easting = false_easting
       if(present(false_northing)) &
            cfp%laea%false_northing = false_northing
    case(CFP_AEA)
       allocate(cfp%aea)
       if(present(longitude_of_central_meridian)) & 
            cfp%aea%longitude_of_central_meridian = longitude_of_central_meridian
       if(present(latitude_of_projection_origin)) &
            cfp%aea%latitude_of_projection_origin = latitude_of_projection_origin
       if(present(false_easting)) &
            cfp%aea%false_easting = false_easting
       if(present(false_northing))&
            cfp%aea%false_northing = false_northing
       if(present(standard_parallel).and.present(standard_parallel_2)) &
            cfp%aea%standard_parallel = (/ standard_parallel,standard_parallel_2 /)
    case(CFP_LCC)
       allocate(cfp%lcc)
       if (present(longitude_of_central_meridian)) &
            cfp%lcc%longitude_of_central_meridian = longitude_of_central_meridian
       if (present(latitude_of_projection_origin)) &
            cfp%lcc%latitude_of_projection_origin = latitude_of_projection_origin
       if (present(false_easting)) &
            cfp%lcc%false_easting = false_easting
       if (present(false_northing)) &
            cfp%lcc%false_northing = false_northing
       if (present(standard_parallel).and.present(standard_parallel_2)) &
            cfp%lcc%standard_parallel = (/ standard_parallel,standard_parallel_2 /)
    case(CFP_STERE)
       allocate(cfp%stere)
       if(present(polar)) cfp%stere%polar = polar
       if(present(longitude_of_central_meridian)) &
            cfp%stere%longitude_of_central_meridian = longitude_of_central_meridian
       if(present(latitude_of_projection_origin)) &
            cfp%stere%latitude_of_projection_origin = latitude_of_projection_origin
       if(present(scale_factor_at_projection_origin)) &
            cfp%stere%scale_factor_at_projection_origin = scale_factor_at_projection_origin
       if(present(standard_parallel)) &
            cfp%stere%standard_parallel = standard_parallel
       if(present(false_easting)) &
            cfp%stere%false_easting = false_easting
       if(present(false_northing)) &
            cfp%stere%false_northing = false_northing
    case default
       call glide_msg(GM_FATAL,__FILE__,__LINE__,'Unrecognised projection type')
    end select

  end subroutine CFproj_define

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! private readers
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function CFproj_get_stere(ncid,mapid)
    use netcdf
    implicit none
    type(CFproj_stere), pointer :: CFproj_get_stere
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid
    
    integer status

    allocate(CFproj_get_stere)
    CFproj_get_stere%polar = .false.
    status = nf90_get_att(ncid,mapid,'false_easting',CFproj_get_stere%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'false_northing',CFproj_get_stere%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'longitude_of_central_meridian',CFproj_get_stere%longitude_of_central_meridian)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'latitude_of_projection_origin',CFproj_get_stere%latitude_of_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'scale_factor_at_projection_origin',CFproj_get_stere%scale_factor_at_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)

  end function CFproj_get_stere

  function CFproj_get_stere_polar(ncid,mapid)
    use netcdf
    use glide_messages
    implicit none
    type(CFproj_stere), pointer :: CFproj_get_stere_polar
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid
    
    integer status
    real dummy

    allocate(CFproj_get_stere_polar)
    CFproj_get_stere_polar%polar = .true.
    status = nf90_get_att(ncid,mapid,'false_easting',CFproj_get_stere_polar%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'false_northing',CFproj_get_stere_polar%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'straight_vertical_longitude_from_pole',CFproj_get_stere_polar%longitude_of_central_meridian)
    call nc_errorhandle(__FILE__,__LINE__,status)
    CFproj_get_stere_polar%latitude_of_projection_origin = 90.
    status = nf90_get_att(ncid,mapid,'scale_factor_at_projection_origin',dummy)
    if (status.eq.NF90_NOERR) then
       CFproj_get_stere_polar%scale_factor_at_projection_origin = dummy
    end if
    status = nf90_get_att(ncid,mapid,'standard_parallel',dummy)
    if (status.eq.NF90_NOERR) then
       CFproj_get_stere_polar%standard_parallel = dummy
    end if
    if (CFproj_get_stere_polar%standard_parallel.ne.0 .and. CFproj_get_stere_polar%scale_factor_at_projection_origin.ne.0.) then
       call glide_msg(GM_FATAL,__FILE__,__LINE__,'(stereographic projection), can only handle either'// &
          ' standard_parallel or scale_factor_at_projection_origin')
    end if
  end function CFproj_get_stere_polar

  function CFproj_get_laea(ncid,mapid)
    use netcdf
    implicit none
    type(CFproj_laea), pointer :: CFproj_get_laea
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid
    
    integer status
    allocate(CFproj_get_laea)
    status = nf90_get_att(ncid,mapid,'false_easting',CFproj_get_laea%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'false_northing',CFproj_get_laea%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'longitude_of_central_meridian',CFproj_get_laea%longitude_of_central_meridian)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'latitude_of_projection_origin',CFproj_get_laea%latitude_of_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)
  end function CFproj_get_laea

  function CFproj_get_aea(ncid,mapid)
    use netcdf
    implicit none
    type(CFproj_aea), pointer :: CFproj_get_aea
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid
    
    integer status
    allocate(CFproj_get_aea)
    status = nf90_get_att(ncid,mapid,'false_easting',CFproj_get_aea%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'false_northing',CFproj_get_aea%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'longitude_of_central_meridian',CFproj_get_aea%longitude_of_central_meridian)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'latitude_of_projection_origin',CFproj_get_aea%latitude_of_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'standard_parallel',CFproj_get_aea%standard_parallel)
    call nc_errorhandle(__FILE__,__LINE__,status)
  end function CFproj_get_aea

  function CFproj_get_lcc(ncid,mapid)
    use netcdf
    implicit none
    type(CFproj_lcc), pointer :: CFproj_get_lcc
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid
    
    integer status
    allocate(CFproj_get_lcc)
    status = nf90_get_att(ncid,mapid,'false_easting',CFproj_get_lcc%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'false_northing',CFproj_get_lcc%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'longitude_of_central_meridian',CFproj_get_lcc%longitude_of_central_meridian)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'latitude_of_projection_origin',CFproj_get_lcc%latitude_of_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'standard_parallel',CFproj_get_lcc%standard_parallel)
    call nc_errorhandle(__FILE__,__LINE__,status)
  end function CFproj_get_lcc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! private converters to proj4 strings
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function CFproj_proj4_stere(stere)
    implicit none
    character(len=proj4len), dimension(:), pointer :: CFproj_proj4_stere
    type(CFproj_stere) :: stere

    allocate(CFproj_proj4_stere(6))
    write(CFproj_proj4_stere(1),*) 'proj=stere'
    write(CFproj_proj4_stere(2),*) 'lon_0=',stere%longitude_of_central_meridian
    write(CFproj_proj4_stere(3),*) 'lat_0=',stere%latitude_of_projection_origin
    if (stere%polar) then
       if (stere%standard_parallel.ne.0) then
          write(CFproj_proj4_stere(4),*) 'lat_ts=',stere%standard_parallel
       else
          write(CFproj_proj4_stere(4),*) 'k_0=',stere%scale_factor_at_projection_origin
       end if
    else
       write(CFproj_proj4_stere(4),*) 'k_0=',stere%scale_factor_at_projection_origin
    end if
    write(CFproj_proj4_stere(5),*) 'x_0=',stere%false_easting
    write(CFproj_proj4_stere(6),*) 'y_0=',stere%false_northing
  end function CFproj_proj4_stere

  function CFproj_proj4_laea(laea)
    implicit none
    character(len=proj4len), dimension(:), pointer :: CFproj_proj4_laea
    type(CFproj_laea) :: laea

    allocate(CFproj_proj4_laea(5))
    write(CFproj_proj4_laea(1),*) 'proj=laea'
    write(CFproj_proj4_laea(2),*) 'lon_0=',laea%longitude_of_central_meridian
    write(CFproj_proj4_laea(3),*) 'lat_0=',laea%latitude_of_projection_origin
    write(CFproj_proj4_laea(4),*) 'x_0=',laea%false_easting
    write(CFproj_proj4_laea(5),*) 'y_0=',laea%false_northing
  end function CFproj_proj4_laea

  function CFproj_proj4_aea(aea)
    implicit none
    character(len=proj4len), dimension(:), pointer :: CFproj_proj4_aea
    type(CFproj_aea) :: aea

    allocate(CFproj_proj4_aea(7))
    write(CFproj_proj4_aea(1),*) 'proj=aea'
    write(CFproj_proj4_aea(2),*) 'lon_0=',aea%longitude_of_central_meridian
    write(CFproj_proj4_aea(3),*) 'lat_0=',aea%latitude_of_projection_origin
    write(CFproj_proj4_aea(4),*) 'lat_1=',aea%standard_parallel(1)
    write(CFproj_proj4_aea(5),*) 'lat_2=',aea%standard_parallel(2)
    write(CFproj_proj4_aea(6),*) 'x_0=',aea%false_easting
    write(CFproj_proj4_aea(7),*) 'y_0=',aea%false_northing
  end function CFproj_proj4_aea

  function CFproj_proj4_lcc(lcc)
    implicit none
    character(len=proj4len), dimension(:), pointer :: CFproj_proj4_lcc
    type(CFproj_lcc) :: lcc

    allocate(CFproj_proj4_lcc(7))
    write(CFproj_proj4_lcc(1),*) 'proj=lcc'
    write(CFproj_proj4_lcc(2),*) 'lon_0=',lcc%longitude_of_central_meridian
    write(CFproj_proj4_lcc(3),*) 'lat_0=',lcc%latitude_of_projection_origin
    write(CFproj_proj4_lcc(4),*) 'lat_1=',lcc%standard_parallel(1)
    write(CFproj_proj4_lcc(5),*) 'lat_2=',lcc%standard_parallel(2)
    write(CFproj_proj4_lcc(6),*) 'x_0=',lcc%false_easting
    write(CFproj_proj4_lcc(7),*) 'y_0=',lcc%false_northing
  end function CFproj_proj4_lcc  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! private subroutines to write projection info
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CFproj_put_stere(ncid,mapid,stere)
    use netcdf
    implicit none
    type(CFproj_stere), pointer :: stere
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid

    integer status

    if (stere%polar) then
       status = nf90_put_att(ncid,mapid,'grid_mapping_name','polar_stereographic')
    else
       status = nf90_put_att(ncid,mapid,'grid_mapping_name','stereographic')
    end if
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'false_easting',stere%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'false_northing',stere%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'longitude_of_central_meridian',stere%longitude_of_central_meridian)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'latitude_of_projection_origin',stere%latitude_of_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (stere%polar) then
       if (stere%standard_parallel.ne.0) then
          status = nf90_put_att(ncid,mapid,'standard_parallel',stere%standard_parallel)
       else
          status = nf90_put_att(ncid,mapid,'scale_factor_at_projection_origin',stere%scale_factor_at_projection_origin)
       end if
    else
       status = nf90_put_att(ncid,mapid,'scale_factor_at_projection_origin',stere%scale_factor_at_projection_origin)
    end if
    call nc_errorhandle(__FILE__,__LINE__,status)
  end subroutine CFproj_put_stere

  subroutine CFproj_put_laea(ncid,mapid,laea)
    use netcdf
    implicit none
    type(CFproj_laea), pointer :: laea
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid

    integer status

    status = nf90_put_att(ncid,mapid,'grid_mapping_name','lambert_azimuthal_equal_area')
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'false_easting',laea%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'false_northing',laea%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'longitude_of_central_meridian',laea%longitude_of_central_meridian)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'latitude_of_projection_origin',laea%latitude_of_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)
  end subroutine CFproj_put_laea

  subroutine CFproj_put_aea(ncid,mapid,aea)
    use netcdf
    implicit none
    type(CFproj_aea), pointer :: aea
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid

    integer status

    status = nf90_put_att(ncid,mapid,'grid_mapping_name','albers_conical_equal_area')
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'false_easting',aea%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'false_northing',aea%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'longitude_of_central_meridian',aea%longitude_of_central_meridian)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'latitude_of_projection_origin',aea%latitude_of_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'standard_parallel',aea%standard_parallel)
    call nc_errorhandle(__FILE__,__LINE__,status)
  end subroutine CFproj_put_aea

  subroutine CFproj_put_lcc(ncid,mapid,lcc)
    use netcdf
    implicit none
    type(CFproj_lcc), pointer :: lcc
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid

    integer status

    status = nf90_put_att(ncid,mapid,'grid_mapping_name','lambert_conformal_conic')
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'false_easting',lcc%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'false_northing',lcc%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'longitude_of_central_meridian',lcc%longitude_of_central_meridian)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'latitude_of_projection_origin',lcc%latitude_of_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'standard_parallel',lcc%standard_parallel)
    call nc_errorhandle(__FILE__,__LINE__,status)
  end subroutine CFproj_put_lcc

end module glimmer_CFproj
