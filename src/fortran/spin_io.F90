!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WARNING: this file was automatically generated on
! Mon, 29 Dec 2008 08:35:15 +0000
! from ncdf_template.F90.in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  ncdf_template.f90 - part of the GLIMMER ice model        + 
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


module spin_io
  !*FD template for creating subsystem specific I/O routines
  !*FD written by Magnus Hagdorn, 2004

  character(len=*),private,parameter :: hotvars = ''

contains

  !*****************************************************************************
  ! netCDF output
  !*****************************************************************************
  subroutine spin_io_createall(model,data,outfiles)
    !*FD open all netCDF files for output
    use spin_types
    use glide_types
    use glimmer_ncdf
    use glimmer_ncio
    implicit none
    type(glide_global_type) :: model
    type(spin_climate_type), optional :: data
    type(glimmer_nc_output),optional,pointer :: outfiles
    
    ! local variables
    type(glimmer_nc_output), pointer :: oc

    if (present(outfiles)) then
       oc => outfiles
    else
       oc=>model%funits%out_first
    end if

    do while(associated(oc))
       if (present(data)) then
          call spin_io_create(oc,model,data)
       else
          call spin_io_create(oc,model)
       end if
       oc=>oc%next
    end do
  end subroutine spin_io_createall

  subroutine spin_io_writeall(data,model,atend,outfiles,time)
    !*FD if necessary write to netCDF files
    use spin_types
    use glide_types
    use glimmer_ncdf
    use glimmer_ncio
    implicit none
    type(spin_climate_type) :: data
    type(glide_global_type) :: model
    logical, optional :: atend
    type(glimmer_nc_output),optional,pointer :: outfiles
    real(sp),optional :: time

    ! local variables
    type(glimmer_nc_output), pointer :: oc
    logical :: forcewrite=.false.

    if (present(outfiles)) then
       oc => outfiles
    else
       oc=>model%funits%out_first
    end if

    if (present(atend)) then
       forcewrite = atend
    end if

    do while(associated(oc))
#ifdef HAVE_AVG
       if (oc%do_averages) then
          call spin_avg_accumulate(oc,data,model)
       end if
#endif
       call glimmer_nc_checkwrite(oc,model,forcewrite,time)
       if (oc%nc%just_processed) then
          ! write standard variables
          call spin_io_write(oc,data)
#ifdef HAVE_AVG
          if (oc%do_averages) then
             call spin_avg_reset(oc,data)
          end if
#endif
       end if
       oc=>oc%next
    end do
  end subroutine spin_io_writeall
  
  subroutine spin_io_create(outfile,model,data)
    use glide_types
    use spin_types
    use glimmer_ncdf
    use glimmer_map_types
    use glimmer_log
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    type(glide_global_type) :: model
    type(spin_climate_type), optional :: data

    integer status,varid,pos

    integer :: time_dimid
    integer :: x1_dimid
    integer :: y1_dimid

    ! defining dimensions
    status = nf90_inq_dimid(NCO%id,'time',time_dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_inq_dimid(NCO%id,'x1',x1_dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_inq_dimid(NCO%id,'y1',y1_dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)

    NCO%vars = ' '//trim(NCO%vars)//' '
    ! expanding hotstart variables
    pos = index(NCO%vars,' hot ') 
    if (pos.ne.0) then
       NCO%vars = NCO%vars(:pos)//NCO%vars(pos+4:)
       NCO%hotstart = .true.
    end if
    if (NCO%hotstart) then
       NCO%vars = trim(NCO%vars)//hotvars
    end if
    ! checking if we need to handle time averages
    pos = index(NCO%vars,"_tavg")
    if (pos.ne.0) then
       outfile%do_averages = .True.
    end if    

    !     ablt -- ablation
    pos = index(NCO%vars,' ablt ')
    status = nf90_inq_varid(NCO%id,'ablt',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
       call write_log('Creating variable ablt')
       status = nf90_def_var(NCO%id,'ablt',NF90_FLOAT,(/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       status = nf90_put_att(NCO%id, varid, 'long_name', 'ablation')
       status = nf90_put_att(NCO%id, varid, 'standard_name', 'ablation')
       status = nf90_put_att(NCO%id, varid, 'units', 'meters/year')
       if (glimmap_allocated(model%projection)) then
          status = nf90_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     end if

    !     arng -- half range temp
    pos = index(NCO%vars,' arng ')
    status = nf90_inq_varid(NCO%id,'arng',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
       call write_log('Creating variable arng')
       status = nf90_def_var(NCO%id,'arng',NF90_FLOAT,(/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       status = nf90_put_att(NCO%id, varid, 'long_name', 'half range temp')
       status = nf90_put_att(NCO%id, varid, 'standard_name', 'degree half range')
       status = nf90_put_att(NCO%id, varid, 'units', 'degree_Celsius')
       if (glimmap_allocated(model%projection)) then
          status = nf90_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     end if

    !     prcp -- precipitation
    pos = index(NCO%vars,' prcp ')
    status = nf90_inq_varid(NCO%id,'prcp',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
       call write_log('Creating variable prcp')
       status = nf90_def_var(NCO%id,'prcp',NF90_FLOAT,(/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       status = nf90_put_att(NCO%id, varid, 'long_name', 'precipitation')
       status = nf90_put_att(NCO%id, varid, 'standard_name', 'precipitation')
       status = nf90_put_att(NCO%id, varid, 'units', 'meters/year')
       if (glimmap_allocated(model%projection)) then
          status = nf90_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     end if

    !     presartm -- annual mean air temperature
    pos = index(NCO%vars,' presartm ')
    status = nf90_inq_varid(NCO%id,'presartm',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
       call write_log('Creating variable presartm')
       status = nf90_def_var(NCO%id,'presartm',NF90_FLOAT,(/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       status = nf90_put_att(NCO%id, varid, 'long_name', 'annual mean air temperature')
       status = nf90_put_att(NCO%id, varid, 'standard_name', 'surface_temperature')
       status = nf90_put_att(NCO%id, varid, 'cell_methods', 'time: mean')
       status = nf90_put_att(NCO%id, varid, 'units', 'degree_Celsius')
       if (glimmap_allocated(model%projection)) then
          status = nf90_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     end if

    !     presprcp -- present precipitation
    pos = index(NCO%vars,' presprcp ')
    status = nf90_inq_varid(NCO%id,'presprcp',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
       call write_log('Creating variable presprcp')
       status = nf90_def_var(NCO%id,'presprcp',NF90_FLOAT,(/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NCO%id, varid, 'coordinates', 'lon lat')
       status = nf90_put_att(NCO%id, varid, 'long_name', 'present precipitation')
       status = nf90_put_att(NCO%id, varid, 'standard_name', 'present precipitation')
       status = nf90_put_att(NCO%id, varid, 'units', 'meters/year')
       if (glimmap_allocated(model%projection)) then
          status = nf90_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     end if

  end subroutine spin_io_create

  subroutine spin_io_write(outfile,data)
    use spin_types
    use glimmer_ncdf
    use glimmer_paramets
    use glimmer_scales
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    !*FD structure containg output netCDF descriptor
    type(spin_climate_type) :: data
    !*FD the model instance

    ! local variables
    real tavgf
    integer status, varid
    integer up
     
    tavgf = outfile%total_time
    if (tavgf.ne.0.) then
       tavgf = 1./tavgf
    end if

    ! write variables
    status = nf90_inq_varid(NCO%id,'ablt',varid)
    if (status .eq. nf90_noerr) then
       status = nf90_put_var(NCO%id, varid, &
            data%mb%ablt, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = nf90_inq_varid(NCO%id,'arng',varid)
    if (status .eq. nf90_noerr) then
       status = nf90_put_var(NCO%id, varid, &
            data%temp%arng, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = nf90_inq_varid(NCO%id,'prcp',varid)
    if (status .eq. nf90_noerr) then
       status = nf90_put_var(NCO%id, varid, &
            data%mb%prcp, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = nf90_inq_varid(NCO%id,'presartm',varid)
    if (status .eq. nf90_noerr) then
       status = nf90_put_var(NCO%id, varid, &
            data%temp%presartm, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = nf90_inq_varid(NCO%id,'presprcp',varid)
    if (status .eq. nf90_noerr) then
       status = nf90_put_var(NCO%id, varid, &
            data%mb%presprcp, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

  end subroutine spin_io_write

  !*****************************************************************************
  ! netCDF input
  !*****************************************************************************  
  subroutine spin_io_readall(data,model)
    !*FD read from netCDF file
    use spin_types
    use glide_types
    use glimmer_ncio
    use glimmer_ncdf
    implicit none
    type(spin_climate_type) :: data
    type(glide_global_type) :: model

    ! local variables
    type(glimmer_nc_input), pointer :: ic    

    ic=>model%funits%in_first
    do while(associated(ic))
       call glimmer_nc_checkread(ic,model)
       if (ic%nc%just_processed) then
          call spin_io_read(ic,data)
       end if
       ic=>ic%next
    end do
  end subroutine spin_io_readall

  subroutine spin_io_read(infile,data,scale_vars)
    !*FD read variables from a netCDF file
    use glimmer_log
    use glimmer_ncdf
    use spin_types
    use glimmer_paramets
    use glimmer_scales
    implicit none
    type(glimmer_nc_input), pointer :: infile
    !*FD structure containg output netCDF descriptor
    type(spin_climate_type) :: data
    !*FD the model instance
    logical,optional :: scale_vars
    !*FD Specifies whether fields should be scaled by factors when read in.

    ! local variables
    integer status,varid
    integer up
    logical :: scale=.true.

    ! Deal with optional argument
    if (present(scale_vars)) scale=scale_vars
   
    ! read variables
    status = nf90_inq_varid(NCI%id,'presartm',varid)
    if (status .eq. nf90_noerr) then
       call write_log('  Loading presartm')
       status = nf90_get_var(NCI%id, varid, &
            data%temp%presartm, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = nf90_inq_varid(NCI%id,'presprcp',varid)
    if (status .eq. nf90_noerr) then
       call write_log('  Loading presprcp')
       status = nf90_get_var(NCI%id, varid, &
            data%mb%presprcp, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

  end subroutine spin_io_read

  subroutine spin_io_checkdim(infile,model,data)
    !*FD check if dimension sizes in file match dims of model
    use glimmer_log
    use glimmer_ncdf
    use glide_types
    use spin_types
    implicit none
    type(glimmer_nc_input), pointer :: infile
    !*FD structure containg output netCDF descriptor
    type(glide_global_type) :: model
    type(spin_climate_type), optional :: data

    integer status,dimid,dimsize
    character(len=150) message

    ! check dimensions
  end subroutine spin_io_checkdim

  !*****************************************************************************
  ! calculating time averages
  !*****************************************************************************  
#ifdef HAVE_AVG
  subroutine spin_avg_accumulate(outfile,data,model)
    use glide_types
    use spin_types
    use glimmer_ncdf
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    !*FD structure containg output netCDF descriptor
    type(glide_global_type) :: model
    type(spin_climate_type) :: data

    ! local variables
    real :: factor
    integer status, varid

    ! increase total time
    outfile%total_time = outfile%total_time + model%numerics%tinc
    factor = model%numerics%tinc

  end subroutine spin_avg_accumulate

  subroutine spin_avg_reset(outfile,data)
    use spin_types
    use glimmer_ncdf
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    !*FD structure containg output netCDF descriptor
    type(spin_climate_type) :: data

    ! local variables
    integer status, varid

    ! reset total time
    outfile%total_time = 0.

  end subroutine spin_avg_reset
#endif

  !*********************************************************************
  ! lots of accessor subroutines follow
  !*********************************************************************
  subroutine spin_get_ablt(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use spin_types
    implicit none
    type(spin_climate_type) :: data
    real, dimension(:,:), intent(out) :: outarray

    outarray = data%mb%ablt
  end subroutine spin_get_ablt

  subroutine spin_set_ablt(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use spin_types
    implicit none
    type(spin_climate_type) :: data
    real, dimension(:,:), intent(in) :: inarray

    data%mb%ablt = inarray
  end subroutine spin_set_ablt

  subroutine spin_get_arng(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use spin_types
    implicit none
    type(spin_climate_type) :: data
    real, dimension(:,:), intent(out) :: outarray

    outarray = data%temp%arng
  end subroutine spin_get_arng

  subroutine spin_set_arng(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use spin_types
    implicit none
    type(spin_climate_type) :: data
    real, dimension(:,:), intent(in) :: inarray

    data%temp%arng = inarray
  end subroutine spin_set_arng

  subroutine spin_get_prcp(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use spin_types
    implicit none
    type(spin_climate_type) :: data
    real, dimension(:,:), intent(out) :: outarray

    outarray = data%mb%prcp
  end subroutine spin_get_prcp

  subroutine spin_set_prcp(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use spin_types
    implicit none
    type(spin_climate_type) :: data
    real, dimension(:,:), intent(in) :: inarray

    data%mb%prcp = inarray
  end subroutine spin_set_prcp

  subroutine spin_get_presartm(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use spin_types
    implicit none
    type(spin_climate_type) :: data
    real, dimension(:,:), intent(out) :: outarray

    outarray = data%temp%presartm
  end subroutine spin_get_presartm

  subroutine spin_set_presartm(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use spin_types
    implicit none
    type(spin_climate_type) :: data
    real, dimension(:,:), intent(in) :: inarray

    data%temp%presartm = inarray
  end subroutine spin_set_presartm

  subroutine spin_get_presprcp(data,outarray)
    use glimmer_scales
    use glimmer_paramets
    use spin_types
    implicit none
    type(spin_climate_type) :: data
    real, dimension(:,:), intent(out) :: outarray

    outarray = data%mb%presprcp
  end subroutine spin_get_presprcp

  subroutine spin_set_presprcp(data,inarray)
    use glimmer_scales
    use glimmer_paramets
    use spin_types
    implicit none
    type(spin_climate_type) :: data
    real, dimension(:,:), intent(in) :: inarray

    data%mb%presprcp = inarray
  end subroutine spin_set_presprcp


end module spin_io
