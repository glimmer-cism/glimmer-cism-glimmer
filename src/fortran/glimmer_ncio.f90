! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_ncio.f90 - part of the GLIMMER ice model         + 
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

#define NCO outfile%nc
#define NCI infile%nc

module glimmer_ncio
  !*FD module for common netCDF I/O
  !*FD written by Magnus Hagdorn, 2004

  use glimmer_ncdf
  integer,parameter,private :: msglen=150
  
contains
  !*****************************************************************************
  ! netCDF output
  !*****************************************************************************  
  subroutine openall_out(model)
    !*FD open all netCDF files for output
    use glide_types
    use glimmer_ncdf
    implicit none
    type(glide_global_type) :: model
    
    ! local variables
    type(glimmer_nc_output), pointer :: oc

    oc=>model%funits%out_first
    do while(associated(oc))
       call glimmer_nc_createfile(oc,model)
       oc=>oc%next
    end do
  end subroutine openall_out

  subroutine closeall_out(model)
    !*FD close all netCDF files for output
    use glide_types
    use glimmer_ncdf
    implicit none
    type(glide_global_type) :: model
    
    ! local variables
    type(glimmer_nc_output), pointer :: oc

    oc=>model%funits%out_first
    do while(associated(oc))
       oc=>delete(oc)
    end do
    model%funits%out_first=>NULL()
  end subroutine closeall_out

  subroutine glimmer_nc_createfile(outfile,model)
    !*FD create a new netCDF file
    use glimmer_log
    use glide_types
    use glimmer_CFproj
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    !*FD structure containg output netCDF descriptor
    type(glide_global_type) :: model
    !*FD the model instance

    ! local variables
    integer status
    integer i,mapid
    character(len=msglen) message

    ! create new netCDF file
    status = nf90_create(NCO%filename,NF90_CLOBBER,NCO%id)
    call nc_errorhandle(__FILE__,__LINE__,status)
    call write_log_div
    write(message,*) 'Opening file ',trim(NCO%filename),' for output; '
    write(message,*) '  Starting output at ',outfile%next_write,' and write every ',outfile%freq,' years'
    call write_log(trim(message))

    ! writing meta data
    status = nf90_put_att(NCO%id, NF90_GLOBAL, 'Conventions', "CF-1.0")
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NCO%id, NF90_GLOBAL,'title',trim(outfile%metadata%title))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NCO%id, NF90_GLOBAL,'institution',trim(outfile%metadata%institution))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NCO%id, NF90_GLOBAL,'source',trim(outfile%metadata%source))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NCO%id, NF90_GLOBAL,'history',trim(outfile%metadata%history))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NCO%id, NF90_GLOBAL,'references',trim(outfile%metadata%references))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NCO%id, NF90_GLOBAL,'comment',trim(outfile%metadata%comment))
    call nc_errorhandle(__FILE__,__LINE__,status)
  
    ! defining time dimension and variable
    status = nf90_def_dim(NCO%id,'time',NF90_UNLIMITED,NCO%timedim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    !     time -- Model time
    call write_log('Creating variable time')
    status = nf90_def_var(NCO%id,'time',NF90_FLOAT,(/NCO%timedim/),NCO%timevar)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NCO%id, NCO%timevar, 'long_name', 'Model time')
    status = nf90_put_att(NCO%id, NCO%timevar, 'standard_name', 'time')
    status = nf90_put_att(NCO%id, NCO%timevar, 'units', 'year since 1-1-1 0:0:0')
    status = nf90_put_att(NCO%id, NCO%timevar, 'calendar', 'none')

    ! adding projection info
    if (CFproj_allocated(model%projection)) then
       status = nf90_def_var(NCO%id,glimmer_nc_mapvarname,NF90_CHAR,mapid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       call CFproj_PutProj(NCO%id,mapid,model%projection)
    end if

    ! setting the size of the level dimension
    NCO%nlevel = model%general%upn
  end subroutine glimmer_nc_createfile

  subroutine glimmer_nc_checkwrite(outfile,model,forcewrite)
    !*FD check if we should write to file
    use glimmer_log
    use glide_types
    implicit none
    type(glimmer_nc_output), pointer :: outfile    
    type(glide_global_type) :: model
    logical forcewrite

    character(len=msglen) :: message
    integer status

    ! check if we are still in define mode and if so leave it
    if (NCO%define_mode) then
       status = nf90_enddef(NCO%id)
       call nc_errorhandle(__FILE__,__LINE__,status)
       NCO%define_mode = .FALSE.
    end if

    if (model%numerics%time.gt.NCO%processsed_time) then
       if (NCO%just_processed) then
          ! finished writing during last time step, need to increase counter...
          
          outfile%timecounter = outfile%timecounter + 1
          status = nf90_sync(NCO%id)
          NCO%just_processed = .FALSE.
       end if
    end if
    if (model%numerics%time.ge.outfile%next_write .or. (forcewrite.and.model%numerics%time.gt.outfile%next_write-outfile%freq)) then
       if (.not.NCO%just_processed) then
          call write_log_div
          write(message,*) 'Writing to file ', trim(NCO%filename), ' at time ', model%numerics%time
          call write_log(trim(message))
          ! increase next_write
          outfile%next_write=outfile%next_write+outfile%freq
          NCO%processsed_time = model%numerics%time
          ! write time
          status = nf90_put_var(NCO%id,NCO%timevar,model%numerics%time,(/outfile%timecounter/))
          NCO%just_processed = .TRUE.         
       end if
    end if
  end subroutine glimmer_nc_checkwrite

  !*****************************************************************************
  ! netCDF input
  !*****************************************************************************  
  subroutine openall_in(model)
    !*FD open all netCDF files for input
    use glide_types
    use glimmer_ncdf
    implicit none
    type(glide_global_type) :: model
    
    ! local variables
    type(glimmer_nc_input), pointer :: ic

    ic=>model%funits%in_first
    do while(associated(ic))
       call glimmer_nc_openfile(ic,model)
       ic=>ic%next
    end do
  end subroutine openall_in

  subroutine closeall_in(model)
    !*FD close all netCDF files for input
    use glide_types
    use glimmer_ncdf
    implicit none
    type(glide_global_type) :: model
    
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
    use glide_types
    use glimmer_cfproj
    use glimmer_log
    implicit none
    type(glimmer_nc_input), pointer :: infile
    !*FD structure containg input netCDF descriptor
    type(glide_global_type) :: model
    !*FD the model instance

    ! local variables
    character(len=50) varname
    integer nvars
    integer i, dimsize
    integer status    
    character(len=msglen) message
    
    ! open netCDF file
    status = nf90_open(NCI%filename,NF90_NOWRITE,NCI%id)
    call nc_errorhandle(__FILE__,__LINE__,status)
    call write_log_div
    call write_log('opening file '//trim(NCI%filename)//' for input')

    ! getting projections
    model%projection = CFproj_GetProj(NCI%id)

    ! getting time dimension
    status = nf90_inq_dimid(NCI%id, 'time', NCI%timedim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    ! get id of time variable
    status = nf90_inq_varid(NCI%id,'time',NCI%timevar)
    call nc_errorhandle(__FILE__,__LINE__,status)
    
    ! getting length of time dimension and allocating memory for array containing times
    status = nf90_inquire_dimension(NCI%id,NCI%timedim,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    allocate(infile%times(dimsize))
    infile%nt=dimsize
    status = nf90_get_var(NCI%id,NCI%timevar,infile%times)

    ! setting the size of the level dimension
    NCI%nlevel = model%general%upn
  end subroutine glimmer_nc_openfile

  subroutine glimmer_nc_checkread(infile,model)
    !*FD check if we should read from file
    use glimmer_log
    use glide_types
    implicit none
    type(glimmer_nc_input), pointer :: infile
    !*FD structure containg output netCDF descriptor
    type(glide_global_type) :: model
    !*FD the model instance

    character(len=msglen) :: message

    if (infile%current_time.le.infile%nt) then
       if (.not.NCI%just_processed) then
          call write_log_div
          write(message,*) 'Reading time slice ',infile%current_time,'(',infile%times(infile%current_time),') from file ', &
               trim(NCI%filename), ' at time ', model%numerics%time
          call write_log(message)
          NCI%just_processed = .TRUE.
          NCI%processsed_time = model%numerics%time
       end if
    end if
    if (model%numerics%time.gt.NCI%processsed_time) then
       if (NCI%just_processed) then
          ! finished reading during last time step, need to increase counter...
          infile%current_time = infile%current_time + 1
          NCI%just_processed = .FALSE.
       end if
    end if
  end subroutine glimmer_nc_checkread

end module glimmer_ncio

