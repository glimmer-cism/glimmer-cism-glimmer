!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WARNING: this file was automatically generated on
! Wed, 12 Jan 2005 11:34:19 +0000
! from ncdf_file.f90.in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

#define NC outfile%nc

module glimmer_ncfile
  !*FD routines for GLIMMER netCDF file I/O
  !*FD written by Magnus Hagdorn, 2004

  character(len=*), private, parameter :: mapvarname = 'mapping'

contains
  subroutine openall_out(model)
    !*FD open all netCDF files for output
    use glimmer_types
    use glimmer_ncdf
    implicit none
    type(glimmer_global_type) :: model
    
    ! local variables
    type(glimmer_nc_output), pointer :: oc

    oc=>model%funits%out_first
    do while(associated(oc))
       call glimmer_nc_createfile(oc,model)
       oc%next_write=model%numerics%time
       oc=>oc%next
    end do
  end subroutine openall_out

  subroutine writeall(model,atend)
    !*FD if necessary write to netCDF files
    use glimmer_types
    use glimmer_ncdf
    implicit none
    type(glimmer_global_type) :: model
    logical, optional :: atend

    ! local variables
    type(glimmer_nc_output), pointer :: oc
    logical :: forcewrite=.false.

    if (present(atend)) then
       forcewrite = atend
    end if

    oc=>model%funits%out_first
    do while(associated(oc))
       if (model%numerics%time.ge.oc%next_write .or. (forcewrite.and.model%numerics%time.gt.oc%next_write-oc%freq)) then
          call glimmer_nc_write(oc,model)
       end if
       oc=>oc%next
    end do
  end subroutine writeall

  subroutine closeall_out(model)
    !*FD close all netCDF files for output
    use glimmer_types
    use glimmer_ncdf
    implicit none
    type(glimmer_global_type) :: model
    
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
    use glimmer_ncdf
    use glimmer_types
	use glide_messages
    use glimmer_CFproj
    use paramets, only : len0
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    !*FD structure containg output netCDF descriptor
    type(glimmer_global_type) :: model
    !*FD the model instance

    ! local variables
    integer status
    integer i,mapid
	character(40) :: freqtxt
    
    ! create new netCDF file
    status = nf90_create(NC%filename,NF90_CLOBBER,NC%id)
    call nc_errorhandle(__FILE__,__LINE__,status)
	call glide_stars
	write(freqtxt,*)outfile%freq
    call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,&
		'Opening file '//trim(NC%filename)//' for output; write period every '//trim(adjustl(freqtxt))//' years')

    ! writing meta data
    status = nf90_put_att(NC%id, NF90_GLOBAL, 'Conventions', "CF-1.0")
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NC%id, NF90_GLOBAL,'title',trim(outfile%metadata%title))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NC%id, NF90_GLOBAL,'institution',trim(outfile%metadata%institution))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NC%id, NF90_GLOBAL,'source',trim(outfile%metadata%source))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NC%id, NF90_GLOBAL,'history',trim(outfile%metadata%history))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NC%id, NF90_GLOBAL,'references',trim(outfile%metadata%references))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NC%id, NF90_GLOBAL,'comment',trim(outfile%metadata%comment))
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! defining dimensions
    status = nf90_def_dim(NC%id,'x0',model%general%ewn-1,NC%x0dim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_def_dim(NC%id,'y0',model%general%nsn-1,NC%y0dim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_def_dim(NC%id,'x1',model%general%ewn,NC%x1dim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_def_dim(NC%id,'y1',model%general%nsn,NC%y1dim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_def_dim(NC%id,'level',model%general%upn,NC%leveldim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (NC%do_spot) then	 
       status = nf90_def_dim(NC%id,'spot',size(outfile%spotx),NC%spotdim)	 
       call nc_errorhandle(__FILE__,__LINE__,status)	 
    end if
    status = nf90_def_dim(NC%id,'time',NF90_UNLIMITED,NC%timedim)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! defining variables
    !     level -- sigma layers
    call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable level')
    status = nf90_def_var(NC%id,'level',NF90_FLOAT,(/NC%leveldim/),NC%levelvar)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NC%id, NC%levelvar, 'units', &
           '1')
    status = nf90_put_att(NC%id, NC%levelvar, 'long_name', &
           'sigma layers')

    !     time -- Model time
    call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable time')
    status = nf90_def_var(NC%id,'time',NF90_FLOAT,(/NC%timedim/),NC%timevar)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NC%id, NC%timevar, 'units', &
           'year')
    status = nf90_put_att(NC%id, NC%timevar, 'long_name', &
           'Model time')

    !     x0 -- Cartisian x-coordinate, midpoint
    call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable x0')
    status = nf90_def_var(NC%id,'x0',NF90_FLOAT,(/NC%x0dim/),NC%x0var)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NC%id, NC%x0var, 'units', &
           'meter')
    status = nf90_put_att(NC%id, NC%x0var, 'long_name', &
           'Cartisian x-coordinate, midpoint')

    !     x0_spot -- Cartisian x-coordinate, midpoint
    if (NC%do_spot) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable x0_spot')
       status = nf90_def_var(NC%id,'x0_spot',NF90_FLOAT,(/NC%spotdim/),NC%x0_spotvar)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%x0_spotvar, 'units', &
              'meter')
       status = nf90_put_att(NC%id, NC%x0_spotvar, 'long_name', &
              'Cartisian x-coordinate, midpoint')
    end if

    !     x1 -- Cartisian x-coordinate
    call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable x1')
    status = nf90_def_var(NC%id,'x1',NF90_FLOAT,(/NC%x1dim/),NC%x1var)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NC%id, NC%x1var, 'units', &
           'meter')
    status = nf90_put_att(NC%id, NC%x1var, 'long_name', &
           'Cartisian x-coordinate')

    !     x1_spot -- Cartisian x-coordinate
    if (NC%do_spot) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable x1_spot')
       status = nf90_def_var(NC%id,'x1_spot',NF90_FLOAT,(/NC%spotdim/),NC%x1_spotvar)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%x1_spotvar, 'units', &
              'meter')
       status = nf90_put_att(NC%id, NC%x1_spotvar, 'long_name', &
              'Cartisian x-coordinate')
    end if

    !     y0 -- Cartisian y-coordinate, midpoint
    call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable y0')
    status = nf90_def_var(NC%id,'y0',NF90_FLOAT,(/NC%y0dim/),NC%y0var)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NC%id, NC%y0var, 'units', &
           'meter')
    status = nf90_put_att(NC%id, NC%y0var, 'long_name', &
           'Cartisian y-coordinate, midpoint')

    !     y0_spot -- Cartisian y-coordinate, midpoint
    if (NC%do_spot) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable y0_spot')
       status = nf90_def_var(NC%id,'y0_spot',NF90_FLOAT,(/NC%spotdim/),NC%y0_spotvar)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%y0_spotvar, 'units', &
              'meter')
       status = nf90_put_att(NC%id, NC%y0_spotvar, 'long_name', &
              'Cartisian y-coordinate, midpoint')
    end if

    !     y1 -- Cartisian y-coordinate
    call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable y1')
    status = nf90_def_var(NC%id,'y1',NF90_FLOAT,(/NC%y1dim/),NC%y1var)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NC%id, NC%y1var, 'units', &
           'meter')
    status = nf90_put_att(NC%id, NC%y1var, 'long_name', &
           'Cartisian y-coordinate')

    !     y1_spot -- Cartisian y-coordinate
    if (NC%do_spot) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable y1_spot')
       status = nf90_def_var(NC%id,'y1_spot',NF90_FLOAT,(/NC%spotdim/),NC%y1_spotvar)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%y1_spotvar, 'units', &
              'meter')
       status = nf90_put_att(NC%id, NC%y1_spotvar, 'long_name', &
              'Cartisian y-coordinate')
    end if

    !     ablt -- ablation
    if (NC%do_var(NC_B_ABLT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable ablt')
       status = nf90_def_var(NC%id,'ablt',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_ABLT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_ABLT), 'long_name', &
              'ablation')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ABLT), 'units', &
              'meter/year')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ABLT), 'grid_mapping',mapvarname)
    end if

    !     ablt_spot -- ablation
    if (NC%do_var(NC_B_ABLT_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable ablt_spot')
       status = nf90_def_var(NC%id,'ablt_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_ABLT_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_ABLT_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ABLT_SPOT), 'long_name', &
              'ablation')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ABLT_SPOT), 'units', &
              'meter/year')
    end if

    !     acab -- accumulation, ablation rate
    if (NC%do_var(NC_B_ACAB)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable acab')
       status = nf90_def_var(NC%id,'acab',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_ACAB))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_ACAB), 'long_name', &
              'accumulation, ablation rate')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ACAB), 'standard_name', &
              'land_ice_surface_mass_balance')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ACAB), 'units', &
              'meter/year')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ACAB), 'grid_mapping',mapvarname)
    end if

    !     acab_spot -- accumulation, ablation rate
    if (NC%do_var(NC_B_ACAB_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable acab_spot')
       status = nf90_def_var(NC%id,'acab_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_ACAB_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_ACAB_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ACAB_SPOT), 'long_name', &
              'accumulation, ablation rate')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ACAB_SPOT), 'standard_name', &
              'land_ice_surface_mass_balance')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ACAB_SPOT), 'units', &
              'meter/year')
    end if

    !     arng -- annual temperature range
    if (NC%do_var(NC_B_ARNG)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable arng')
       status = nf90_def_var(NC%id,'arng',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_ARNG))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_ARNG), 'units', &
              'degree_Celcius')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ARNG), 'long_name', &
              'annual temperature range')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ARNG), 'grid_mapping',mapvarname)
    end if

    !     arng_spot -- annual temperature range
    if (NC%do_var(NC_B_ARNG_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable arng_spot')
       status = nf90_def_var(NC%id,'arng_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_ARNG_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_ARNG_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ARNG_SPOT), 'long_name', &
              'annual temperature range')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ARNG_SPOT), 'units', &
              'degree_Celcius')
    end if

    !     artm -- annual mean air temperature
    if (NC%do_var(NC_B_ARTM)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable artm')
       status = nf90_def_var(NC%id,'artm',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_ARTM))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_ARTM), 'long_name', &
              'annual mean air temperature')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ARTM), 'standard_name', &
              'surface_temperature')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ARTM), 'cell_methods', &
              'time: mean')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ARTM), 'units', &
              'degree_Celcius')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ARTM), 'grid_mapping',mapvarname)
    end if

    !     artm_spot -- annual mean air temperature
    if (NC%do_var(NC_B_ARTM_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable artm_spot')
       status = nf90_def_var(NC%id,'artm_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_ARTM_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_ARTM_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ARTM_SPOT), 'long_name', &
              'annual mean air temperature')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ARTM_SPOT), 'standard_name', &
              'surface_temperature')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ARTM_SPOT), 'cell_methods', &
              'time: mean')
       status = nf90_put_att(NC%id, NC%varids(NC_B_ARTM_SPOT), 'units', &
              'degree_Celcius')
    end if

    !     bmlt -- basal melt rate
    if (NC%do_var(NC_B_BMLT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable bmlt')
       status = nf90_def_var(NC%id,'bmlt',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_BMLT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_BMLT), 'long_name', &
              'basal melt rate')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BMLT), 'standard_name', &
              'land_ice_basal_melt_rate')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BMLT), 'units', &
              'meter/year')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BMLT), 'grid_mapping',mapvarname)
    end if

    !     bmlt_spot -- basal melt rate
    if (NC%do_var(NC_B_BMLT_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable bmlt_spot')
       status = nf90_def_var(NC%id,'bmlt_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_BMLT_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_BMLT_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BMLT_SPOT), 'long_name', &
              'basal melt rate')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BMLT_SPOT), 'standard_name', &
              'land_ice_basal_melt_rate')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BMLT_SPOT), 'units', &
              'meter/year')
    end if

    !     btemp -- basal ice temperature
    if (NC%do_var(NC_B_BTEMP)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable btemp')
       status = nf90_def_var(NC%id,'btemp',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_BTEMP))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_BTEMP), 'units', &
              'degree_Celcius')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BTEMP), 'long_name', &
              'basal ice temperature')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BTEMP), 'grid_mapping',mapvarname)
    end if

    !     btemp_spot -- basal ice temperature
    if (NC%do_var(NC_B_BTEMP_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable btemp_spot')
       status = nf90_def_var(NC%id,'btemp_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_BTEMP_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_BTEMP_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BTEMP_SPOT), 'long_name', &
              'basal ice temperature')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BTEMP_SPOT), 'units', &
              'degree_Celcius')
    end if

    !     btrc -- basal slip coefficient
    if (NC%do_var(NC_B_BTRC)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable btrc')
       status = nf90_def_var(NC%id,'btrc',NF90_FLOAT,(/NC%x0dim, NC%y0dim, NC%timedim/),NC%varids(NC_B_BTRC))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_BTRC), 'long_name', &
              'basal slip coefficient')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BTRC), 'units', &
              'meter/pascal/year')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BTRC), 'grid_mapping',mapvarname)
    end if

    !     btrc_spot -- basal slip coefficient
    if (NC%do_var(NC_B_BTRC_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable btrc_spot')
       status = nf90_def_var(NC%id,'btrc_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_BTRC_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_BTRC_SPOT), 'coordinates', &
              'y0_spot x0_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BTRC_SPOT), 'long_name', &
              'basal slip coefficient')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BTRC_SPOT), 'units', &
              'meter/pascal/year')
    end if

    !     bwat -- basal water depth
    if (NC%do_var(NC_B_BWAT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable bwat')
       status = nf90_def_var(NC%id,'bwat',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_BWAT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_BWAT), 'long_name', &
              'basal water depth')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BWAT), 'units', &
              'meter')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BWAT), 'grid_mapping',mapvarname)
    end if

    !     bwat_spot -- basal water depth
    if (NC%do_var(NC_B_BWAT_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable bwat_spot')
       status = nf90_def_var(NC%id,'bwat_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_BWAT_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_BWAT_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BWAT_SPOT), 'long_name', &
              'basal water depth')
       status = nf90_put_att(NC%id, NC%varids(NC_B_BWAT_SPOT), 'units', &
              'meter')
    end if

    !     diffu -- apparent diffusivity
    if (NC%do_var(NC_B_DIFFU)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable diffu')
       status = nf90_def_var(NC%id,'diffu',NF90_FLOAT,(/NC%x0dim, NC%y0dim, NC%timedim/),NC%varids(NC_B_DIFFU))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_DIFFU), 'long_name', &
              'apparent diffusivity')
       status = nf90_put_att(NC%id, NC%varids(NC_B_DIFFU), 'units', &
              'meter2/year')
       status = nf90_put_att(NC%id, NC%varids(NC_B_DIFFU), 'grid_mapping',mapvarname)
    end if

    !     diffu_spot -- apparent diffusivity
    if (NC%do_var(NC_B_DIFFU_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable diffu_spot')
       status = nf90_def_var(NC%id,'diffu_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_DIFFU_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_DIFFU_SPOT), 'coordinates', &
              'y0_spot x0_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_DIFFU_SPOT), 'long_name', &
              'apparent diffusivity')
       status = nf90_put_att(NC%id, NC%varids(NC_B_DIFFU_SPOT), 'units', &
              'meter2/year')
    end if

    !     dusrfdtm -- rate of upper ice surface elevation change
    if (NC%do_var(NC_B_DUSRFDTM)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable dusrfdtm')
       status = nf90_def_var(NC%id,'dusrfdtm',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_DUSRFDTM))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_DUSRFDTM), 'long_name', &
              'rate of upper ice surface elevation change')
       status = nf90_put_att(NC%id, NC%varids(NC_B_DUSRFDTM), 'units', &
              'meter/year')
       status = nf90_put_att(NC%id, NC%varids(NC_B_DUSRFDTM), 'grid_mapping',mapvarname)
    end if

    !     dusrfdtm_spot -- rate of upper ice surface elevation change
    if (NC%do_var(NC_B_DUSRFDTM_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable dusrfdtm_spot')
       status = nf90_def_var(NC%id,'dusrfdtm_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_DUSRFDTM_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_DUSRFDTM_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_DUSRFDTM_SPOT), 'long_name', &
              'rate of upper ice surface elevation change')
       status = nf90_put_att(NC%id, NC%varids(NC_B_DUSRFDTM_SPOT), 'units', &
              'meter/year')
    end if

    !     flwa -- ??
    if (NC%do_var(NC_B_FLWA)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable flwa')
       status = nf90_def_var(NC%id,'flwa',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%leveldim, NC%timedim/),NC%varids(NC_B_FLWA))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_FLWA), 'long_name', &
              '??')
       status = nf90_put_att(NC%id, NC%varids(NC_B_FLWA), 'units', &
              '??')
       status = nf90_put_att(NC%id, NC%varids(NC_B_FLWA), 'grid_mapping',mapvarname)
    end if

    !     flwa_spot -- ??
    if (NC%do_var(NC_B_FLWA_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable flwa_spot')
       status = nf90_def_var(NC%id,'flwa_spot',NF90_FLOAT,(/NC%spotdim, NC%leveldim, NC%timedim/),NC%varids(NC_B_FLWA_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_FLWA_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_FLWA_SPOT), 'long_name', &
              '??')
       status = nf90_put_att(NC%id, NC%varids(NC_B_FLWA_SPOT), 'units', &
              '??')
    end if

    !     lat -- Latitude
    if (NC%do_var(NC_B_LAT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable lat')
       status = nf90_def_var(NC%id,'lat',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_LAT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_LAT), 'long_name', &
              'Latitude')
       status = nf90_put_att(NC%id, NC%varids(NC_B_LAT), 'standard_name', &
              'latitude')
       status = nf90_put_att(NC%id, NC%varids(NC_B_LAT), 'units', &
              'degreeN')
       status = nf90_put_att(NC%id, NC%varids(NC_B_LAT), 'grid_mapping',mapvarname)
    end if

    !     lat_spot -- Latitude
    if (NC%do_var(NC_B_LAT_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable lat_spot')
       status = nf90_def_var(NC%id,'lat_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_LAT_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_LAT_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_LAT_SPOT), 'long_name', &
              'Latitude')
       status = nf90_put_att(NC%id, NC%varids(NC_B_LAT_SPOT), 'standard_name', &
              'latitude')
       status = nf90_put_att(NC%id, NC%varids(NC_B_LAT_SPOT), 'units', &
              'degreeN')
    end if

    !     lon -- Longitude
    if (NC%do_var(NC_B_LON)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable lon')
       status = nf90_def_var(NC%id,'lon',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_LON))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_LON), 'long_name', &
              'Longitude')
       status = nf90_put_att(NC%id, NC%varids(NC_B_LON), 'standard_name', &
              'longitude')
       status = nf90_put_att(NC%id, NC%varids(NC_B_LON), 'units', &
              'degreeE')
       status = nf90_put_att(NC%id, NC%varids(NC_B_LON), 'grid_mapping',mapvarname)
    end if

    !     lon_spot -- Longitude
    if (NC%do_var(NC_B_LON_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable lon_spot')
       status = nf90_def_var(NC%id,'lon_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_LON_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_LON_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_LON_SPOT), 'long_name', &
              'Longitude')
       status = nf90_put_att(NC%id, NC%varids(NC_B_LON_SPOT), 'standard_name', &
              'longitude')
       status = nf90_put_att(NC%id, NC%varids(NC_B_LON_SPOT), 'units', &
              'degreeE')
    end if

    !     lsurf -- ice lower surface elevation
    if (NC%do_var(NC_B_LSURF)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable lsurf')
       status = nf90_def_var(NC%id,'lsurf',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_LSURF))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_LSURF), 'long_name', &
              'ice lower surface elevation')
       status = nf90_put_att(NC%id, NC%varids(NC_B_LSURF), 'units', &
              'meter')
       status = nf90_put_att(NC%id, NC%varids(NC_B_LSURF), 'grid_mapping',mapvarname)
    end if

    !     lsurf_spot -- ice lower surface elevation
    if (NC%do_var(NC_B_LSURF_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable lsurf_spot')
       status = nf90_def_var(NC%id,'lsurf_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_LSURF_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_LSURF_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_LSURF_SPOT), 'long_name', &
              'ice lower surface elevation')
       status = nf90_put_att(NC%id, NC%varids(NC_B_LSURF_SPOT), 'units', &
              'meter')
    end if

    !     mask -- upscaling and downscaling mask
    if (NC%do_var(NC_B_MASK)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable mask')
       status = nf90_def_var(NC%id,'mask',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_MASK))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_MASK), 'long_name', &
              'upscaling and downscaling mask')
       status = nf90_put_att(NC%id, NC%varids(NC_B_MASK), 'units', &
              '1')
       status = nf90_put_att(NC%id, NC%varids(NC_B_MASK), 'grid_mapping',mapvarname)
    end if

    !     mask_spot -- upscaling and downscaling mask
    if (NC%do_var(NC_B_MASK_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable mask_spot')
       status = nf90_def_var(NC%id,'mask_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_MASK_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_MASK_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_MASK_SPOT), 'long_name', &
              'upscaling and downscaling mask')
       status = nf90_put_att(NC%id, NC%varids(NC_B_MASK_SPOT), 'units', &
              '1')
    end if

    !     prcp -- precipitation
    if (NC%do_var(NC_B_PRCP)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable prcp')
       status = nf90_def_var(NC%id,'prcp',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_PRCP))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRCP), 'long_name', &
              'precipitation')
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRCP), 'standard_name', &
              'lwe_precipitation_rate')
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRCP), 'units', &
              'meter/year')
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRCP), 'grid_mapping',mapvarname)
    end if

    !     prcp_spot -- precipitation
    if (NC%do_var(NC_B_PRCP_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable prcp_spot')
       status = nf90_def_var(NC%id,'prcp_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_PRCP_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRCP_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRCP_SPOT), 'long_name', &
              'precipitation')
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRCP_SPOT), 'standard_name', &
              'lwe_precipitation_rate')
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRCP_SPOT), 'units', &
              'meter/year')
    end if

    !     presprcp -- present day precipitation
    if (NC%do_var(NC_B_PRESPRCP)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable presprcp')
       status = nf90_def_var(NC%id,'presprcp',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_PRESPRCP))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRESPRCP), 'long_name', &
              'present day precipitation')
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRESPRCP), 'units', &
              'meter/year')
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRESPRCP), 'grid_mapping',mapvarname)
    end if

    !     presprcp_spot -- present day precipitation
    if (NC%do_var(NC_B_PRESPRCP_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable presprcp_spot')
       status = nf90_def_var(NC%id,'presprcp_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_PRESPRCP_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRESPRCP_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRESPRCP_SPOT), 'long_name', &
              'present day precipitation')
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRESPRCP_SPOT), 'units', &
              'meter/year')
    end if

    !     presusrf -- present day surface of the ice-sheet
    if (NC%do_var(NC_B_PRESUSRF)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable presusrf')
       status = nf90_def_var(NC%id,'presusrf',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_PRESUSRF))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRESUSRF), 'long_name', &
              'present day surface of the ice-sheet')
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRESUSRF), 'units', &
              'meter')
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRESUSRF), 'grid_mapping',mapvarname)
    end if

    !     presusrf_spot -- present day surface of the ice-sheet
    if (NC%do_var(NC_B_PRESUSRF_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable presusrf_spot')
       status = nf90_def_var(NC%id,'presusrf_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_PRESUSRF_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRESUSRF_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRESUSRF_SPOT), 'long_name', &
              'present day surface of the ice-sheet')
       status = nf90_put_att(NC%id, NC%varids(NC_B_PRESUSRF_SPOT), 'units', &
              'meter')
    end if

    !     relx -- relaxed bedrock topography
    if (NC%do_var(NC_B_RELX)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable relx')
       status = nf90_def_var(NC%id,'relx',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_RELX))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_RELX), 'long_name', &
              'relaxed bedrock topography')
       status = nf90_put_att(NC%id, NC%varids(NC_B_RELX), 'units', &
              'meter')
       status = nf90_put_att(NC%id, NC%varids(NC_B_RELX), 'grid_mapping',mapvarname)
    end if

    !     relx_spot -- relaxed bedrock topography
    if (NC%do_var(NC_B_RELX_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable relx_spot')
       status = nf90_def_var(NC%id,'relx_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_RELX_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_RELX_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_RELX_SPOT), 'long_name', &
              'relaxed bedrock topography')
       status = nf90_put_att(NC%id, NC%varids(NC_B_RELX_SPOT), 'units', &
              'meter')
    end if

    !     std_dev -- standard deviation of sub-grid topography
    if (NC%do_var(NC_B_STD_DEV)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable std_dev')
       status = nf90_def_var(NC%id,'std_dev',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_STD_DEV))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_STD_DEV), 'long_name', &
              'standard deviation of sub-grid topography')
       status = nf90_put_att(NC%id, NC%varids(NC_B_STD_DEV), 'units', &
              'meter')
       status = nf90_put_att(NC%id, NC%varids(NC_B_STD_DEV), 'grid_mapping',mapvarname)
    end if

    !     std_dev_spot -- standard deviation of sub-grid topography
    if (NC%do_var(NC_B_STD_DEV_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable std_dev_spot')
       status = nf90_def_var(NC%id,'std_dev_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_STD_DEV_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_STD_DEV_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_STD_DEV_SPOT), 'long_name', &
              'standard deviation of sub-grid topography')
       status = nf90_put_att(NC%id, NC%varids(NC_B_STD_DEV_SPOT), 'units', &
              'meter')
    end if

    !     temp -- ice temperature
    if (NC%do_var(NC_B_TEMP)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable temp')
       status = nf90_def_var(NC%id,'temp',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%leveldim, NC%timedim/),NC%varids(NC_B_TEMP))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_TEMP), 'long_name', &
              'ice temperature')
       status = nf90_put_att(NC%id, NC%varids(NC_B_TEMP), 'standard_name', &
              'land_ice_temperature')
       status = nf90_put_att(NC%id, NC%varids(NC_B_TEMP), 'units', &
              'degree_Celcius')
       status = nf90_put_att(NC%id, NC%varids(NC_B_TEMP), 'grid_mapping',mapvarname)
    end if

    !     temp_spot -- ice temperature
    if (NC%do_var(NC_B_TEMP_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable temp_spot')
       status = nf90_def_var(NC%id,'temp_spot',NF90_FLOAT,(/NC%spotdim, NC%leveldim, NC%timedim/),NC%varids(NC_B_TEMP_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_TEMP_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_TEMP_SPOT), 'long_name', &
              'ice temperature')
       status = nf90_put_att(NC%id, NC%varids(NC_B_TEMP_SPOT), 'standard_name', &
              'land_ice_temperature')
       status = nf90_put_att(NC%id, NC%varids(NC_B_TEMP_SPOT), 'units', &
              'degree_Celcius')
    end if

    !     thk -- ice thickness
    if (NC%do_var(NC_B_THK)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable thk')
       status = nf90_def_var(NC%id,'thk',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_THK))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_THK), 'long_name', &
              'ice thickness')
       status = nf90_put_att(NC%id, NC%varids(NC_B_THK), 'standard_name', &
              'land_ice_thickness')
       status = nf90_put_att(NC%id, NC%varids(NC_B_THK), 'units', &
              'meter')
       status = nf90_put_att(NC%id, NC%varids(NC_B_THK), 'grid_mapping',mapvarname)
    end if

    !     thk_spot -- ice thickness
    if (NC%do_var(NC_B_THK_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable thk_spot')
       status = nf90_def_var(NC%id,'thk_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_THK_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_THK_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_THK_SPOT), 'long_name', &
              'ice thickness')
       status = nf90_put_att(NC%id, NC%varids(NC_B_THK_SPOT), 'standard_name', &
              'land_ice_thickness')
       status = nf90_put_att(NC%id, NC%varids(NC_B_THK_SPOT), 'units', &
              'meter')
    end if

    !     topg -- bedrock topography
    if (NC%do_var(NC_B_TOPG)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable topg')
       status = nf90_def_var(NC%id,'topg',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_TOPG))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_TOPG), 'long_name', &
              'bedrock topography')
       status = nf90_put_att(NC%id, NC%varids(NC_B_TOPG), 'standard_name', &
              'bedrock_altitude')
       status = nf90_put_att(NC%id, NC%varids(NC_B_TOPG), 'units', &
              'meter')
       status = nf90_put_att(NC%id, NC%varids(NC_B_TOPG), 'grid_mapping',mapvarname)
    end if

    !     topg_spot -- bedrock topography
    if (NC%do_var(NC_B_TOPG_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable topg_spot')
       status = nf90_def_var(NC%id,'topg_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_TOPG_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_TOPG_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_TOPG_SPOT), 'long_name', &
              'bedrock topography')
       status = nf90_put_att(NC%id, NC%varids(NC_B_TOPG_SPOT), 'standard_name', &
              'bedrock_altitude')
       status = nf90_put_att(NC%id, NC%varids(NC_B_TOPG_SPOT), 'units', &
              'meter')
    end if

    !     ubas -- basal slip velocity in x direction
    if (NC%do_var(NC_B_UBAS)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable ubas')
       status = nf90_def_var(NC%id,'ubas',NF90_FLOAT,(/NC%x0dim, NC%y0dim, NC%timedim/),NC%varids(NC_B_UBAS))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_UBAS), 'long_name', &
              'basal slip velocity in x direction')
       status = nf90_put_att(NC%id, NC%varids(NC_B_UBAS), 'standard_name', &
              'land_ice_basal_x_velocity')
       status = nf90_put_att(NC%id, NC%varids(NC_B_UBAS), 'units', &
              'meter/year')
       status = nf90_put_att(NC%id, NC%varids(NC_B_UBAS), 'grid_mapping',mapvarname)
    end if

    !     ubas_spot -- basal slip velocity in x direction
    if (NC%do_var(NC_B_UBAS_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable ubas_spot')
       status = nf90_def_var(NC%id,'ubas_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_UBAS_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_UBAS_SPOT), 'coordinates', &
              'y0_spot x0_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_UBAS_SPOT), 'long_name', &
              'basal slip velocity in x direction')
       status = nf90_put_att(NC%id, NC%varids(NC_B_UBAS_SPOT), 'standard_name', &
              'land_ice_basal_x_velocity')
       status = nf90_put_att(NC%id, NC%varids(NC_B_UBAS_SPOT), 'units', &
              'meter/year')
    end if

    !     uflx -- flux in x direction
    if (NC%do_var(NC_B_UFLX)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable uflx')
       status = nf90_def_var(NC%id,'uflx',NF90_FLOAT,(/NC%x0dim, NC%y0dim, NC%timedim/),NC%varids(NC_B_UFLX))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_UFLX), 'long_name', &
              'flux in x direction')
       status = nf90_put_att(NC%id, NC%varids(NC_B_UFLX), 'units', &
              'meter2/year')
       status = nf90_put_att(NC%id, NC%varids(NC_B_UFLX), 'grid_mapping',mapvarname)
    end if

    !     uflx_spot -- flux in x direction
    if (NC%do_var(NC_B_UFLX_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable uflx_spot')
       status = nf90_def_var(NC%id,'uflx_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_UFLX_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_UFLX_SPOT), 'coordinates', &
              'y0_spot x0_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_UFLX_SPOT), 'long_name', &
              'flux in x direction')
       status = nf90_put_att(NC%id, NC%varids(NC_B_UFLX_SPOT), 'units', &
              'meter2/year')
    end if

    !     usurf -- ice upper surface elevation
    if (NC%do_var(NC_B_USURF)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable usurf')
       status = nf90_def_var(NC%id,'usurf',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%timedim/),NC%varids(NC_B_USURF))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_USURF), 'long_name', &
              'ice upper surface elevation')
       status = nf90_put_att(NC%id, NC%varids(NC_B_USURF), 'units', &
              'meter')
       status = nf90_put_att(NC%id, NC%varids(NC_B_USURF), 'grid_mapping',mapvarname)
    end if

    !     usurf_spot -- ice upper surface elevation
    if (NC%do_var(NC_B_USURF_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable usurf_spot')
       status = nf90_def_var(NC%id,'usurf_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_USURF_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_USURF_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_USURF_SPOT), 'long_name', &
              'ice upper surface elevation')
       status = nf90_put_att(NC%id, NC%varids(NC_B_USURF_SPOT), 'units', &
              'meter')
    end if

    !     uvel -- ice velocity in x direction
    if (NC%do_var(NC_B_UVEL)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable uvel')
       status = nf90_def_var(NC%id,'uvel',NF90_FLOAT,(/NC%x0dim, NC%y0dim, NC%leveldim, NC%timedim/),NC%varids(NC_B_UVEL))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_UVEL), 'long_name', &
              'ice velocity in x direction')
       status = nf90_put_att(NC%id, NC%varids(NC_B_UVEL), 'standard_name', &
              'land_ice_x_velocity')
       status = nf90_put_att(NC%id, NC%varids(NC_B_UVEL), 'units', &
              'meter/year')
       status = nf90_put_att(NC%id, NC%varids(NC_B_UVEL), 'grid_mapping',mapvarname)
    end if

    !     uvel_spot -- ice velocity in x direction
    if (NC%do_var(NC_B_UVEL_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable uvel_spot')
       status = nf90_def_var(NC%id,'uvel_spot',NF90_FLOAT,(/NC%spotdim, NC%leveldim, NC%timedim/),NC%varids(NC_B_UVEL_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_UVEL_SPOT), 'coordinates', &
              'y0_spot x0_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_UVEL_SPOT), 'long_name', &
              'ice velocity in x direction')
       status = nf90_put_att(NC%id, NC%varids(NC_B_UVEL_SPOT), 'standard_name', &
              'land_ice_x_velocity')
       status = nf90_put_att(NC%id, NC%varids(NC_B_UVEL_SPOT), 'units', &
              'meter/year')
    end if

    !     vbas -- basal slip velocity in y direction
    if (NC%do_var(NC_B_VBAS)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable vbas')
       status = nf90_def_var(NC%id,'vbas',NF90_FLOAT,(/NC%x0dim, NC%y0dim, NC%timedim/),NC%varids(NC_B_VBAS))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_VBAS), 'long_name', &
              'basal slip velocity in y direction')
       status = nf90_put_att(NC%id, NC%varids(NC_B_VBAS), 'standard_name', &
              'land_ice_basal_y_velocity')
       status = nf90_put_att(NC%id, NC%varids(NC_B_VBAS), 'units', &
              'meter/year')
       status = nf90_put_att(NC%id, NC%varids(NC_B_VBAS), 'grid_mapping',mapvarname)
    end if

    !     vbas_spot -- basal slip velocity in y direction
    if (NC%do_var(NC_B_VBAS_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable vbas_spot')
       status = nf90_def_var(NC%id,'vbas_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_VBAS_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_VBAS_SPOT), 'coordinates', &
              'y0_spot x0_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_VBAS_SPOT), 'long_name', &
              'basal slip velocity in y direction')
       status = nf90_put_att(NC%id, NC%varids(NC_B_VBAS_SPOT), 'standard_name', &
              'land_ice_basal_y_velocity')
       status = nf90_put_att(NC%id, NC%varids(NC_B_VBAS_SPOT), 'units', &
              'meter/year')
    end if

    !     vflx -- flux in x direction
    if (NC%do_var(NC_B_VFLX)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable vflx')
       status = nf90_def_var(NC%id,'vflx',NF90_FLOAT,(/NC%x0dim, NC%y0dim, NC%timedim/),NC%varids(NC_B_VFLX))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_VFLX), 'long_name', &
              'flux in x direction')
       status = nf90_put_att(NC%id, NC%varids(NC_B_VFLX), 'units', &
              'meter2/year')
       status = nf90_put_att(NC%id, NC%varids(NC_B_VFLX), 'grid_mapping',mapvarname)
    end if

    !     vflx_spot -- flux in x direction
    if (NC%do_var(NC_B_VFLX_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable vflx_spot')
       status = nf90_def_var(NC%id,'vflx_spot',NF90_FLOAT,(/NC%spotdim, NC%timedim/),NC%varids(NC_B_VFLX_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_VFLX_SPOT), 'coordinates', &
              'y0_spot x0_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_VFLX_SPOT), 'long_name', &
              'flux in x direction')
       status = nf90_put_att(NC%id, NC%varids(NC_B_VFLX_SPOT), 'units', &
              'meter2/year')
    end if

    !     vvel -- ice velocity in y direction
    if (NC%do_var(NC_B_VVEL)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable vvel')
       status = nf90_def_var(NC%id,'vvel',NF90_FLOAT,(/NC%x0dim, NC%y0dim, NC%leveldim, NC%timedim/),NC%varids(NC_B_VVEL))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_VVEL), 'long_name', &
              'ice velocity in y direction')
       status = nf90_put_att(NC%id, NC%varids(NC_B_VVEL), 'standard_name', &
              'land_ice_y_velocity')
       status = nf90_put_att(NC%id, NC%varids(NC_B_VVEL), 'units', &
              'meter/year')
       status = nf90_put_att(NC%id, NC%varids(NC_B_VVEL), 'grid_mapping',mapvarname)
    end if

    !     vvel_spot -- ice velocity in y direction
    if (NC%do_var(NC_B_VVEL_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable vvel_spot')
       status = nf90_def_var(NC%id,'vvel_spot',NF90_FLOAT,(/NC%spotdim, NC%leveldim, NC%timedim/),NC%varids(NC_B_VVEL_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_VVEL_SPOT), 'coordinates', &
              'y0_spot x0_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_VVEL_SPOT), 'long_name', &
              'ice velocity in y direction')
       status = nf90_put_att(NC%id, NC%varids(NC_B_VVEL_SPOT), 'standard_name', &
              'land_ice_y_velocity')
       status = nf90_put_att(NC%id, NC%varids(NC_B_VVEL_SPOT), 'units', &
              'meter/year')
    end if

    !     wgrd -- ?? some velo ??
    if (NC%do_var(NC_B_WGRD)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable wgrd')
       status = nf90_def_var(NC%id,'wgrd',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%leveldim, NC%timedim/),NC%varids(NC_B_WGRD))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_WGRD), 'long_name', &
              '?? some velo ??')
       status = nf90_put_att(NC%id, NC%varids(NC_B_WGRD), 'units', &
              'meter/year')
       status = nf90_put_att(NC%id, NC%varids(NC_B_WGRD), 'grid_mapping',mapvarname)
    end if

    !     wgrd_spot -- ?? some velo ??
    if (NC%do_var(NC_B_WGRD_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable wgrd_spot')
       status = nf90_def_var(NC%id,'wgrd_spot',NF90_FLOAT,(/NC%spotdim, NC%leveldim, NC%timedim/),NC%varids(NC_B_WGRD_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_WGRD_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_WGRD_SPOT), 'long_name', &
              '?? some velo ??')
       status = nf90_put_att(NC%id, NC%varids(NC_B_WGRD_SPOT), 'units', &
              'meter/year')
    end if

    !     wvel -- vertical ice velocity
    if (NC%do_var(NC_B_WVEL)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable wvel')
       status = nf90_def_var(NC%id,'wvel',NF90_FLOAT,(/NC%x1dim, NC%y1dim, NC%leveldim, NC%timedim/),NC%varids(NC_B_WVEL))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_WVEL), 'long_name', &
              'vertical ice velocity')
       status = nf90_put_att(NC%id, NC%varids(NC_B_WVEL), 'standard_name', &
              'land_ice_z_velocity')
       status = nf90_put_att(NC%id, NC%varids(NC_B_WVEL), 'units', &
              'meter/year')
       status = nf90_put_att(NC%id, NC%varids(NC_B_WVEL), 'grid_mapping',mapvarname)
    end if

    !     wvel_spot -- vertical ice velocity
    if (NC%do_var(NC_B_WVEL_SPOT)) then
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'Creating variable wvel_spot')
       status = nf90_def_var(NC%id,'wvel_spot',NF90_FLOAT,(/NC%spotdim, NC%leveldim, NC%timedim/),NC%varids(NC_B_WVEL_SPOT))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = nf90_put_att(NC%id, NC%varids(NC_B_WVEL_SPOT), 'coordinates', &
              'y1_spot x1_spot')
       status = nf90_put_att(NC%id, NC%varids(NC_B_WVEL_SPOT), 'long_name', &
              'vertical ice velocity')
       status = nf90_put_att(NC%id, NC%varids(NC_B_WVEL_SPOT), 'standard_name', &
              'land_ice_z_velocity')
       status = nf90_put_att(NC%id, NC%varids(NC_B_WVEL_SPOT), 'units', &
              'meter/year')
    end if


    ! adding projection info
    status = nf90_def_var(NC%id,mapvarname,NF90_CHAR,mapid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    call CFproj_PutProj(NC%id,mapid,model%projection)

    ! leaving define mode
    status = nf90_enddef(NC%id)
    call nc_errorhandle(__FILE__,__LINE__,status)
    
    ! filling coordinate variables
    do i=1, model%general%ewn-1
       status=nf90_put_var(NC%id,NC%x0var,((i-0.5)*model%numerics%dew*len0),(/i/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end do
    do i=1, model%general%nsn-1
       status=nf90_put_var(NC%id,NC%y0var,(i-0.5)*model%numerics%dns*len0,(/i/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end do
    do i=1, model%general%ewn
       status=nf90_put_var(NC%id,NC%x1var,(i-1.)*model%numerics%dew*len0,(/i/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end do
    do i=1, model%general%nsn
       status=nf90_put_var(NC%id,NC%y1var,(i-1.)*model%numerics%dns*len0,(/i/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end do    
    status=nf90_put_var(NC%id,NC%levelvar,model%numerics%sigma)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (NC%do_spot) then
       do i=1,size(outfile%spotx)
          status=nf90_put_var(NC%id,NC%x0_spotvar, &
               (real(outfile%spotx(i))-0.5)*real(model%numerics%dew*len0),(/i/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status=nf90_put_var(NC%id,NC%y0_spotvar, &
               (real(outfile%spoty(i))-0.5)*real(model%numerics%dns*len0),(/i/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          
          status=nf90_put_var(NC%id,NC%x1_spotvar, &
               (real(outfile%spotx(i))-1.0)*real(model%numerics%dew*len0),(/i/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status=nf90_put_var(NC%id,NC%y1_spotvar, &
               (real(outfile%spoty(i))-1.0)*real(model%numerics%dns*len0),(/i/))
          call nc_errorhandle(__FILE__,__LINE__,status)    
       end do
    end if
  end subroutine glimmer_nc_createfile

  subroutine glimmer_nc_write(outfile,model)
    !*FD write variables to a netCDF file
    use glimmer_ncdf
    use glimmer_types
    use glimmer_global, only : dp
	use glide_messages
    use physcon, only : scyr
    use paramets, only : thk0, vis0, acc0
    use glimmer_scales
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    !*FD structure containg output netCDF descriptor
    type(glimmer_global_type) :: model
    !*FD the model instance

    ! local variables
    integer status
    integer up, spot
    integer :: ewnv, nsnv
	character(40) :: timetxt

    ! Set up various constants. These were originally only
    ! done once, but are done each time now, for safety.

    ewnv = model%general%ewn - 1
    nsnv = model%general%nsn - 1

	call glide_stars
	write(timetxt,*)model%numerics%time
    call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__, &
	    'Writing to file '//trim(NC%filename)//' at time '//trim(adjustl(timetxt)))

    ! increase next_write 
    outfile%next_write=outfile%next_write+outfile%freq
    ! write time
    status = nf90_put_var(NC%id,NC%timevar,model%numerics%time,(/outfile%timecounter/))

    ! write variables
    if (NC%do_var(NC_B_ABLT)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_ABLT), &
            (scale2d_f1)*(model%climate%ablt), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_ABLT_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_ABLT_SPOT), &
               (scale2d_f1)*(model%climate%ablt(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_ACAB)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_ACAB), &
            (scale2d_f1)*(model%climate%acab), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_ACAB_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_ACAB_SPOT), &
               (scale2d_f1)*(model%climate%acab(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_ARNG)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_ARNG), &
            model%climate%arng, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_ARNG_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_ARNG_SPOT), &
               model%climate%arng(outfile%spotx(spot),outfile%spoty(spot)), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_ARTM)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_ARTM), &
            model%climate%artm, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_ARTM_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_ARTM_SPOT), &
               model%climate%artm(outfile%spotx(spot),outfile%spoty(spot)), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_BMLT)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_BMLT), &
            (scale2d_f1)*(model%temper%bmlt), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_BMLT_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_BMLT_SPOT), &
               (scale2d_f1)*(model%temper%bmlt(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_BTEMP)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_BTEMP), &
            model%temper%temp(model%general%upn,1:model%general%ewn,1:model%general%nsn), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_BTEMP_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_BTEMP_SPOT), &
               model%temper%temp(model%general%upn,outfile%spotx(spot),outfile%spoty(spot)), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_BTRC)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_BTRC), &
            (scale2d_f6)*(model%velocity%btrc(1:ewnv,1:nsnv)), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_BTRC_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_BTRC_SPOT), &
               (scale2d_f6)*(model%velocity%btrc(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_BWAT)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_BWAT), &
            (thk0)*(model%temper%bwat), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_BWAT_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_BWAT_SPOT), &
               (thk0)*(model%temper%bwat(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_DIFFU)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_DIFFU), &
            (scale2d_f4)*(model%velocity%diffu(1:ewnv,1:nsnv)), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_DIFFU_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_DIFFU_SPOT), &
               (scale2d_f4)*(model%velocity%diffu(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_DUSRFDTM)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_DUSRFDTM), &
            (scale2d_f1)*(model%geomderv%dusrfdtm), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_DUSRFDTM_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_DUSRFDTM_SPOT), &
               (scale2d_f1)*(model%geomderv%dusrfdtm(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_FLWA)) then
       do up=1,model%general%upn
          status = nf90_put_var(NC%id, NC%varids(NC_B_FLWA), &
               (vis0*scyr)*(model%temper%flwa(up,:,:)), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_FLWA_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_FLWA_SPOT), &
               (vis0*scyr)*(model%temper%flwa(:,outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,1,outfile%timecounter/),(/1,model%general%upn,1/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_LAT)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_LAT), &
            model%climate%lati, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_LAT_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_LAT_SPOT), &
               model%climate%lati(outfile%spotx(spot),outfile%spoty(spot)), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_LON)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_LON), &
            model%climate%loni, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_LON_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_LON_SPOT), &
               model%climate%loni(outfile%spotx(spot),outfile%spoty(spot)), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_LSURF)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_LSURF), &
            (thk0)*(model%geometry%lsrf), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_LSURF_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_LSURF_SPOT), &
               (thk0)*(model%geometry%lsrf(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_MASK)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_MASK), &
            model%climate%out_mask, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_MASK_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_MASK_SPOT), &
               model%climate%out_mask(outfile%spotx(spot),outfile%spoty(spot)), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_PRCP)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_PRCP), &
            (scale2d_f1)*(model%climate%prcp), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_PRCP_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_PRCP_SPOT), &
               (scale2d_f1)*(model%climate%prcp(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_PRESPRCP)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_PRESPRCP), &
            (scyr * acc0)*(model%climate%presprcp), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_PRESPRCP_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_PRESPRCP_SPOT), &
               (scyr * acc0)*(model%climate%presprcp(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_PRESUSRF)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_PRESUSRF), &
            (thk0)*(model%climate%presusrf), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_PRESUSRF_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_PRESUSRF_SPOT), &
               (thk0)*(model%climate%presusrf(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_RELX)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_RELX), &
            (thk0)*(model%geometry%relx), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_RELX_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_RELX_SPOT), &
               (thk0)*(model%geometry%relx(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_STD_DEV)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_STD_DEV), &
            (thk0)*(model%geometry%std_dev), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_STD_DEV_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_STD_DEV_SPOT), &
               (thk0)*(model%geometry%std_dev(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_TEMP)) then
       do up=1,model%general%upn
          status = nf90_put_var(NC%id, NC%varids(NC_B_TEMP), &
               model%temper%temp(up,1:model%general%ewn,1:model%general%nsn), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_TEMP_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_TEMP_SPOT), &
               model%temper%temp(:,outfile%spotx(spot),outfile%spoty(spot)), &
               (/spot,1,outfile%timecounter/),(/1,model%general%upn,1/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_THK)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_THK), &
            (thk0)*(model%geometry%thck), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_THK_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_THK_SPOT), &
               (thk0)*(model%geometry%thck(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_TOPG)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_TOPG), &
            (thk0)*(model%geometry%topg), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_TOPG_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_TOPG_SPOT), &
               (thk0)*(model%geometry%topg(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_UBAS)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_UBAS), &
            (scale2d_f5)*(model%velocity%ubas(1:ewnv,1:nsnv)), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_UBAS_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_UBAS_SPOT), &
               (scale2d_f5)*(model%velocity%ubas(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_UFLX)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_UFLX), &
            (scale2d_f2)*(model%velocity%uflx(1:ewnv,1:nsnv)), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_UFLX_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_UFLX_SPOT), &
               (scale2d_f2)*(model%velocity%uflx(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_USURF)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_USURF), &
            (thk0)*(model%geometry%usrf), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_USURF_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_USURF_SPOT), &
               (thk0)*(model%geometry%usrf(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_UVEL)) then
       do up=1,model%general%upn
          status = nf90_put_var(NC%id, NC%varids(NC_B_UVEL), &
               (scale3d_f1)*(model%velocity%uvel(up,1:ewnv,1:nsnv)), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_UVEL_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_UVEL_SPOT), &
               (scale3d_f1)*(model%velocity%uvel(:,outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,1,outfile%timecounter/),(/1,model%general%upn,1/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_VBAS)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_VBAS), &
            (scale2d_f5)*(model%velocity%vbas(1:ewnv,1:nsnv)), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_VBAS_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_VBAS_SPOT), &
               (scale2d_f5)*(model%velocity%vbas(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_VFLX)) then
       status = nf90_put_var(NC%id, NC%varids(NC_B_VFLX), &
            (scale2d_f2)*(model%velocity%vflx(1:ewnv,1:nsnv)), (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    if (NC%do_var(NC_B_VFLX_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_VFLX_SPOT), &
               (scale2d_f2)*(model%velocity%vflx(outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_VVEL)) then
       do up=1,model%general%upn
          status = nf90_put_var(NC%id, NC%varids(NC_B_VVEL), &
               (scale3d_f1)*(model%velocity%vvel(up,1:ewnv,1:nsnv)), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_VVEL_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_VVEL_SPOT), &
               (scale3d_f1)*(model%velocity%vvel(:,outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,1,outfile%timecounter/),(/1,model%general%upn,1/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_WGRD)) then
       do up=1,model%general%upn
          status = nf90_put_var(NC%id, NC%varids(NC_B_WGRD), &
               (scale3d_f7)*(model%velocity%wgrd(up,:,:)), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_WGRD_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_WGRD_SPOT), &
               (scale3d_f7)*(model%velocity%wgrd(:,outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,1,outfile%timecounter/),(/1,model%general%upn,1/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_WVEL)) then
       do up=1,model%general%upn
          status = nf90_put_var(NC%id, NC%varids(NC_B_WVEL), &
               (scale3d_f7)*(model%velocity%wvel(up,:,:)), (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    if (NC%do_var(NC_B_WVEL_SPOT)) then
       do spot=1,size(outfile%spotx)
          status = nf90_put_var(NC%id, NC%varids(NC_B_WVEL_SPOT), &
               (scale3d_f7)*(model%velocity%wvel(:,outfile%spotx(spot),outfile%spoty(spot))), &
               (/spot,1,outfile%timecounter/),(/1,model%general%upn,1/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if


    ! increase time counter
    outfile%timecounter = outfile%timecounter + 1

    status = nf90_sync(NC%id)

  end subroutine glimmer_nc_write
end module glimmer_ncfile
