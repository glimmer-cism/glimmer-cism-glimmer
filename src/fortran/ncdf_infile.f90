
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  ncdf_file.f90 - part of the GLIMMER ice model            + 
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

! N.B. This hard-wired file replaces the automatically-generated
! NetCDF handling routines previously used, in order to remove the
! dependency on Python 2.3

#define NC infile%nc

module glimmer_ncinfile
  !*FD routines for GLIMMER netCDF file input
  !*FD written by Magnus Hagdorn, 2004

contains
  subroutine openall_in(model)
    !*FD open all netCDF files for input
    use glimmer_types
    use glimmer_ncdf
    implicit none
    type(glimmer_global_type) :: model
    
    ! local variables
    type(glimmer_nc_input), pointer :: ic

    ic=>model%funits%in_first
    do while(associated(ic))
       call glimmer_nc_openfile(ic,model)
       ic=>ic%next
    end do
  end subroutine openall_in

  subroutine readall(model)
    !*FD read from netCDF file
    use glimmer_types
    use glimmer_ncdf
    implicit none
    type(glimmer_global_type) :: model

    ! local variables
    type(glimmer_nc_input), pointer :: ic    

    ic=>model%funits%in_first
    do while(associated(ic))
       if (ic%current_time.le.ic%nt) then
          call glimmer_nc_read(ic,model)
       end if
       ic=>ic%next
    end do
  end subroutine readall

  subroutine closeall_in(model)
    !*FD close all netCDF files for input
    use glimmer_types
    use glimmer_ncdf
    implicit none
    type(glimmer_global_type) :: model
    
    ! local variables
    type(glimmer_nc_input), pointer :: ic

    ic=>model%funits%in_first
    do while(associated(ic))
       ic=>delete(ic)
    end do
    model%funits%in_first=>NULL()
  end subroutine closeall_in

  subroutine glimmer_nc_openfile(infile,model)
    !*FD open an existing netCDF file
    use glimmer_ncdf
    use glimmer_types
    use glimmer_cfproj
	use glide_messages
    implicit none
    type(glimmer_nc_input), pointer :: infile
    !*FD structure containg input netCDF descriptor
    type(glimmer_global_type) :: model
    !*FD the model instance

    ! local variables
    character(len=50) varname
    integer nvars
    integer i, dimsize
    integer status
	character(80) :: errtxt    
    
    ! open netCDF file
    status = nf90_open(NC%filename,NF90_NOWRITE,NC%id)
    call nc_errorhandle(__FILE__,__LINE__,status)
	call glide_stars
    call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'opening file '//trim(NC%filename)//' for input')

    ! getting projections
    model%projection = CFproj_GetProj(NC%id)

    ! reading dimensions
    ! getting ids
    status = nf90_inq_dimid(NC%id, 'x0', NC%x0dim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_inq_dimid(NC%id, 'y0', NC%y0dim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_inq_dimid(NC%id, 'x1', NC%x1dim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_inq_dimid(NC%id, 'y1', NC%y1dim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_inq_dimid(NC%id, 'level', NC%leveldim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_inq_dimid(NC%id, 'time', NC%timedim)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! getting lengths of x, y and level dimensions
    ! bail out if they do not match
    status = nf90_inquire_dimension(NC%id,NC%x0dim,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (dimsize.ne.model%general%ewn-1) then
	   write(errtxt,*)dimsize,model%general%ewn-1
       call glide_msg(GM_FATAL,__FILE__,__LINE__,'reading file '//trim(NC%filename)//&
	     ' size x0dim does not match: '//trim(errtxt))
    end if
    status = nf90_inquire_dimension(NC%id,NC%y0dim,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (dimsize.ne.model%general%nsn-1) then
       write(errtxt,*)dimsize,model%general%nsn-1
       call glide_msg(GM_FATAL,__FILE__,__LINE__,'reading file '//trim(NC%filename)//&
         ' size y0dim does not match: '//trim(errtxt))
    end if
    status = nf90_inquire_dimension(NC%id,NC%x1dim,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (dimsize.ne.model%general%ewn) then
       write(errtxt,*)dimsize,model%general%ewn
       call glide_msg(GM_FATAL,__FILE__,__LINE__,'reading file '//trim(NC%filename)//&
         ' size x1dim does not match: '//trim(errtxt))
    end if
    status = nf90_inquire_dimension(NC%id,NC%y1dim,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (dimsize.ne.model%general%nsn) then
       write(errtxt,*)dimsize,model%general%nsn
       call glide_msg(GM_FATAL,__FILE__,__LINE__,'reading file '//trim(NC%filename)//&
         ' size y1dim does not match: '//trim(errtxt))
    end if  

    ! getting variables
    status = nf90_inquire(NC%id,nvariables=nvars)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! get id of time variable
    status = nf90_inq_varid(NC%id,'time',NC%timevar)
    call nc_errorhandle(__FILE__,__LINE__,status)

    do i=1,nvars
       status = nf90_inquire_variable(NC%id,i,name=varname)
       call nc_errorhandle(__FILE__,__LINE__,status)
       select case(varname)
       case('lat')
          NC%do_var(NC_B_LAT) = .true.
          NC%varids(NC_B_LAT) = i
       case('lon')
          NC%do_var(NC_B_LON) = .true.
          NC%varids(NC_B_LON) = i
       case('mask')
          NC%do_var(NC_B_MASK) = .true.
          NC%varids(NC_B_MASK) = i
       case('presprcp')
          NC%do_var(NC_B_PRESPRCP) = .true.
          NC%varids(NC_B_PRESPRCP) = i
       case('presusrf')
          NC%do_var(NC_B_PRESUSRF) = .true.
          NC%varids(NC_B_PRESUSRF) = i
       case('relx')
          NC%do_var(NC_B_RELX) = .true.
          NC%varids(NC_B_RELX) = i
       case('topg')
          NC%do_var(NC_B_TOPG) = .true.
          NC%varids(NC_B_TOPG) = i
       case('usurf')
          NC%do_var(NC_B_USURF) = .true.
          NC%varids(NC_B_USURF) = i
       end select
    end do

    ! getting length of time dimension and allocating memory for array containing times
    status = nf90_inquire_dimension(NC%id,NC%timedim,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    allocate(infile%times(dimsize))
    infile%nt=dimsize
    status = nf90_get_var(NC%id,NC%timevar,infile%times)
  end subroutine glimmer_nc_openfile

  subroutine glimmer_nc_read(infile,model,scale_vars)
    !*FD read variables from a netCDF file
    use glimmer_ncdf
    use glimmer_types
    use glide_messages
    use glimmer_global, only : dp
    use physcon, only : scyr
    use paramets, only : thk0, acc0
    use glimmer_scales
    implicit none
    type(glimmer_nc_input), pointer :: infile
    !*FD structure containg output netCDF descriptor
    type(glimmer_global_type) :: model
    !*FD the model instance
    logical,optional :: scale_vars
    !*FD Specifies whether fields should be scaled by factors when read in.

    ! local variables
    integer status
    integer up
    integer :: ewnv, nsnv
    logical :: scale=.true.
    character(40) :: outtxt1,outtxt2,outtxt3

    ! Deal with optional argument

    if (present(scale_vars)) scale=scale_vars

    ! Set up various constants. These were originally only
    ! done once, but are done each time now, for safety.

    ewnv = model%general%ewn - 1
    nsnv = model%general%nsn - 1

	call glide_stars
    write(outtxt1,*)infile%current_time
    write(outtxt2,*)infile%times(infile%current_time)
    write(outtxt3,*)model%numerics%time
    call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Reading time slice '//trim(adjustl(outtxt1))//&
        '('//trim(adjustl(outtxt2))//') from file '//trim(NC%filename)//' at time '//trim(adjustl(outtxt3)))
    
    ! read variables
    if (NC%do_var(NC_B_LAT)) then
      call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Loading lat')
       status = nf90_get_var(NC%id, NC%varids(NC_B_LAT), &
            model%climate%lati, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_LON)) then
      call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Loading lon')
       status = nf90_get_var(NC%id, NC%varids(NC_B_LON), &
            model%climate%loni, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_MASK)) then
      call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Loading mask')
       status = nf90_get_var(NC%id, NC%varids(NC_B_MASK), &
            model%climate%out_mask, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_PRESPRCP)) then
      call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Loading presprcp')
       status = nf90_get_var(NC%id, NC%varids(NC_B_PRESPRCP), &
            model%climate%presprcp, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       if (scale) model%climate%presprcp = model%climate%presprcp/(scyr * acc0)
    end if

    if (NC%do_var(NC_B_PRESUSRF)) then
      call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Loading presusrf')
       status = nf90_get_var(NC%id, NC%varids(NC_B_PRESUSRF), &
            model%climate%presusrf, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       if (scale) model%climate%presusrf = model%climate%presusrf/(thk0)
    end if

    if (NC%do_var(NC_B_RELX)) then
      call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Loading relx')
       status = nf90_get_var(NC%id, NC%varids(NC_B_RELX), &
            model%geometry%relx, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       if (scale) model%geometry%relx = model%geometry%relx/(thk0)
    end if

    if (NC%do_var(NC_B_TOPG)) then
      call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Loading topg')
       status = nf90_get_var(NC%id, NC%varids(NC_B_TOPG), &
            model%geometry%topg, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       if (scale) model%geometry%topg = model%geometry%topg/(thk0)
    end if

    if (NC%do_var(NC_B_USURF)) then
      call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Loading usurf')
       status = nf90_get_var(NC%id, NC%varids(NC_B_USURF), &
            model%geometry%usrf, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       if (scale) model%geometry%usrf = model%geometry%usrf/(thk0)
    end if

    infile%current_time = infile%current_time + 1
  end subroutine glimmer_nc_read

end module glimmer_ncinfile
