
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_outp.f90 - part of the GLIMMER ice model         + 
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

module glimmer_outp

  use glimmer_types

contains

  subroutine initout(model,unit)

    use glimmer_global, only: n2d, n3d
    use paramets, only : len0 

    implicit none

    type(glimmer_global_type) :: model
    integer,intent(in) :: unit
    logical :: there
    integer :: ew,ns
    integer,dimension(n2d) :: which2d
    integer,dimension(n3d) :: which3d

    namelist / vars / which2d, which3d

    model%funits%b0dfile  = trim(model%funits%output_stem) // '.gl0'
    model%funits%b2dfile  = trim(model%funits%output_stem) // '.gl2'
    model%funits%b3dfile  = trim(model%funits%output_stem) // '.gl3'
    model%funits%outfile  = trim(model%funits%output_stem) // '.glw'

    ! Open 0d file and write out header information.

    open(unit,file=model%funits%b0dfile ,form='unformatted')

    write(unit)   size(model%funits%indices0dx)  
    write(unit) ((real(model%funits%indices0dx(ew))-0.5)*real(model%numerics%dew*len0), &
                  ew = 1,size(model%funits%indices0dx))
    write(unit) ((real(model%funits%indices0dy(ew))-0.5)*real(model%numerics%dns*len0), &
                  ew = 1,size(model%funits%indices0dy))
    write(unit) ((real(model%funits%indices0dx(ew))-1.0)*real(model%numerics%dew*len0), &
                  ew = 1,size(model%funits%indices0dx))
    write(unit) ((real(model%funits%indices0dy(ew))-1.0)*real(model%numerics%dns*len0), &
                  ew = 1,size(model%funits%indices0dy))
    close(unit)

    ! Open 3d file and write out header information.

    open(unit,file=model%funits%b3dfile,form='unformatted')

    write(unit) 3, 0, model%general%nsn-1, model%general%ewn-1, model%general%upn
    write(unit) (real((ns-0.5)*model%numerics%dns*len0), ns = 1, model%general%nsn-1)
    write(unit) (real((ew-0.5)*model%numerics%dew*len0), ew = 1, model%general%ewn-1)

    write(unit) 3, 1, model%general%nsn, model%general%ewn, model%general%upn
    write(unit) (real((ns-1.0)*model%numerics%dns*len0), ns = 1, model%general%nsn)
    write(unit) (real((ew-1.0)*model%numerics%dew*len0), ew = 1, model%general%ewn)

    close(unit)

    ! Open 2d file and write out header information.

    open(unit,file=model%funits%b2dfile,form='unformatted')

    write(unit) 2, 0, model%general%nsn-1, model%general%ewn-1 
    write(unit) (real((ns-0.5)*model%numerics%dns*len0), ns = 1, model%general%nsn-1)
    write(unit) (real((ew-0.5)*model%numerics%dew*len0), ew = 1, model%general%ewn-1)
   
    write(unit) 2, 1, model%general%nsn, model%general%ewn 
    write(unit) (real((ns-1.0)*model%numerics%dns*len0), ns = 1, model%general%nsn)
    write(unit) (real((ew-1.0)*model%numerics%dew*len0), ew = 1, model%general%ewn)

    close(unit)

    ! Check to see if an output configuration file exists.
    ! If so, read it in.

    inquire (exist=there,file=model%funits%outfile)
  
    if (there) then
      open(unit,file=model%funits%outfile)
      read(unit,nml=vars)
      close(unit)
      model%funits%which2d=which2d
      model%funits%which3d=which3d
    end if

  end subroutine initout

!-------------------------------------------------------------------------

  subroutine writ2dvr(model,unit)

    ! *** output for 2d fields units:

    ! *** positions      m
    ! *** thickness etc    m
    ! *** horizontal fluxes   m^2 yr^-1
    ! *** accumulation etc     m yr^-1
    ! *** basal melt    m yr^-1
    ! *** basal water    m

    use glimmer_global, only : dp
    use physcon, only : scyr
    use paramets, only : thk0, tim0, vel0, vis0, len0, tau0
   
    implicit none

    type(glimmer_global_type),intent(in) :: model
    integer,intent(in) :: unit

    real(dp):: f1, f2, f3, f4, f5, f6, f7, f8
    integer :: ewnv, nsnv

    ! Set up various constants. These were originally only
    ! done once, but are done each time now, for safety.

    ewnv = model%general%ewn - 1
    nsnv = model%general%nsn - 1
    f1 = scyr * thk0 / tim0
    f2 = scyr * vel0 * thk0
    f3 = vel0 / (vis0 * len0)
    f4 = vel0 * scyr * len0
    f5 = scyr * vel0
    f6 = scyr * vel0 * len0 / (thk0**2)
    f7 = tau0
    f8 = tau0 * len0 / (scyr * vel0) 

    print*,'******* DOING 2D OUTPUT at ',model%numerics%time,' *******'

    open(unit,file=model%funits%b2dfile,position='append',form='unformatted')

    if (model%funits%which2d(1) == 1) then
      write(unit) model%numerics%time, 1, 0, real(f2)
      write(unit) real(f2*model%velocity%uflx(1:ewnv,1:nsnv))
    end if

    if (model%funits%which2d(2) == 1) then
      write(unit) model%numerics%time, 2, 0, real(f2)
      write(unit) real(f2*model%velocity%vflx(1:ewnv,1:nsnv))
    end if

    if (model%funits%which2d(3) == 1) then
      write(unit) model%numerics%time, 3, 0, real(f4)
      write(unit) real(f4*model%velocity%diffu(1:ewnv,1:nsnv))
    end if

    if (model%funits%which2d(4) == 1) then
      write(unit) model%numerics%time, 4, 0, real(f6)
      write(unit) real(f6*model%velocity%btrc(1:ewnv,1:nsnv))
    end if       

    if (model%funits%which2d(5) == 1) then
      write(unit) model%numerics%time, 5, 0, real(f5)
      write(unit) real(f5*model%velocity%ubas(1:ewnv,1:nsnv))
    end if

    if (model%funits%which2d(6) == 1) then
      write(unit) model%numerics%time, 6, 0, real(f5)
      write(unit) real(f5*model%velocity%vbas(1:ewnv,1:nsnv))
    end if
  
    if (model%funits%which2d(7) == 1) then
      write(unit) model%numerics%time, 1, 1, real(thk0)
      write(unit) real(thk0*model%geometry%thck)
    end if
    
    if (model%funits%which2d(8) == 1) then
      write(unit) model%numerics%time, 2, 1, real(thk0)
      write(unit) real(thk0*model%geometry%usrf)
    end if

    if (model%funits%which2d(9) == 1) then
      write(unit) model%numerics%time, 3, 1, real(thk0)
      write(unit) real(thk0*model%geometry%lsrf)
    end if

    if (model%funits%which2d(10) == 1) then
      write(unit) model%numerics%time, 4, 1, real(thk0)
      write(unit) real(thk0*model%geometry%topg)
    end if

    if (model%funits%which2d(11) == 1) then
      write(unit) model%numerics%time, 5, 1, real(f1)
      write(unit) real(f1*model%climate%acab)
    end if

    if (model%funits%which2d(12) == 1) then
      write(unit) model%numerics%time, 6, 1, real(f1)
      write(unit) real(f1*model%temper%bmlt)
    end if

    if (model%funits%which2d(13) == 1) then
      write(unit) model%numerics%time, 7, 1, real(thk0)
      write(unit) real(thk0*model%temper%bwat)
    end if

    if (model%funits%which2d(14) == 1) then
      write(unit) model%numerics%time, 8, 1, 1.0
      write(unit) real(model%climate%artm)
    end if

    if (model%funits%which2d(15) == 1) then
      write(unit) model%numerics%time, 9, 1, 1.0
      write(unit) real(model%temper%temp(model%general%upn,1:model%general%ewn,1:model%general%nsn))
    end if

    if (model%funits%which2d(16) == 1) then
      write(unit) model%numerics%time, 10, 1, 1.0
      write(unit) real(model%climate%arng)
    end if

    if (model%funits%which2d(17) == 1) then
      write(unit) model%numerics%time, 11, 1, real(f1)
      write(unit) real(f1*model%climate%prcp)
    end if

    if (model%funits%which2d(18) == 1) then
      write(unit) model%numerics%time, 12, 1, real(f1)
      write(unit) real(f1*model%climate%ablt)
    end if

    if (model%funits%which2d(19) == 1) then
      write(unit) model%numerics%time, 13, 1, real(f1)
      write(unit) real(f1*model%geomderv%dusrfdtm)
    end if

    close(unit)

  end subroutine writ2dvr

!-------------------------------------------------------------------------

  subroutine writ3dvr(model,unit)

    ! *** output for 3d fields units:

    ! *** positions      m
    ! *** horizontal velocities   m yr^-1
    ! *** effective viscocity  Pa s
    ! *** vertical velocities   m yr^-1
    ! *** temperatures    deg. C
    ! *** flow constant    Pa^-3 yr^-1

    use glimmer_global, only : dp
    use physcon, only : scyr, gn
    use paramets, only : thk0, tim0, vel0, vis0, len0 

    implicit none

    type(glimmer_global_type) :: model
    integer,intent(in) :: unit

    real(sp), allocatable, dimension(:,:,:) :: z 
    real(dp) :: f1, f2, f3, f4, f5, f6
    integer  :: ewnv, nsnv
    integer :: ns,ew,up

    ! These calculations don't really need to be done
    ! every time, but are done so for safety.

    ewnv = model%general%ewn - 1
    nsnv = model%general%nsn - 1
    f1 = scyr * vel0
    f2 = vis0 * (vel0/len0)**(gn - 1)
    f3 = scyr * thk0
    f4 = vel0/(vis0*len0)
    f5 = 1.0d0/f2**(1.0/gn)
    f6 = f4**(1.0/gn); 

    allocate(z(model%general%upn,model%general%ewn,model%general%nsn))
         
    do ns = 1,nsnv
      do ew = 1,ewnv
        do up = 1,model%general%upn
          z(up,ew,ns) = thk0 * (sum(model%geometry%usrf(ew:ew+1,ns:ns+1))/4.0 - model%numerics%sigma(up) * &
                        sum(model%geometry%thck(ew:ew+1,ns:ns+1))/4.0) 
        end do
      end do
    end do

    open(unit,file=model%funits%b3dfile,position='append',form='unformatted')

    write(unit) model%numerics%time, 1, 0, 1.0
    write(unit) z(:,1:ewnv,1:nsnv)

    if (model%funits%which3d(1) == 1) then
      write(unit) model%numerics%time, 2, 0, real(f1)
      write(unit) real(f1*model%velocity%uvel(:,1:ewnv,1:nsnv))
    end if

    if (model%funits%which3d(2) == 1) then
      write(unit) model%numerics%time, 3, 0, real(f1)
      write(unit) real(f1*model%velocity%vvel(:,1:ewnv,1:nsnv))
    end if

    do ns = 1,model%general%nsn 
      do ew = 1, model%general%ewn
        do up = 1, model%general%upn
          z(up,ew,ns) = thk0 * (model%geometry%usrf(ew,ns) - model%numerics%sigma(up) * &
                        model%geometry%thck(ew,ns)) 
        end do
      end do
    end do

    write(unit) model%numerics%time, 1, 1, 1.0
    write(unit) z 

    if (model%funits%which3d(3) == 1) then
      write(unit) model%numerics%time, 2, 1, real(f3/tim0)
      write(unit) real(f3*model%velocity%wvel/tim0)
    end if

    if (model%funits%which3d(4) == 1) then
      write(unit) model%numerics%time, 3, 1, real(f3/tim0)
      write(unit) real(f3*model%velocity%wgrd/tim0)
    end if

    if (model%funits%which3d(5) == 1) then
      write(unit) model%numerics%time, 4, 1, real(vis0*scyr)
      write(unit) real(model%temper%flwa*vis0*scyr)
    end if

    if (model%funits%which3d(6) == 1) then
      write(unit) model%numerics%time, 5, 1, 1.0
      write(unit) real(model%temper%temp(:,1:model%general%ewn,1:model%general%nsn))
    end if

    deallocate(z)

    close(unit)

  end subroutine writ3dvr

!-------------------------------------------------------------------------

  subroutine writ0dvr(model,unit) 

    use glimmer_global, only : sp
    use physcon, only : scyr 
    use paramets, only : thk0, vel0, len0, tim0 

    implicit none

    type(glimmer_global_type) :: model
    integer,intent(in) :: unit
    real(sp) :: f1, f2, f3, f4, f5, f6
    integer :: i

    ! Again, don't need to do this every time, but safer
    ! this way...

    f1 = thk0 * len0**2 * model%numerics%dew * model%numerics%dns
    f2 = len0**2 * model%numerics%dew * model%numerics%dns 
    f3 = scyr * thk0 / tim0
    f4 = scyr * vel0 * thk0
    f5 = scyr * vel0
    f6 = scyr * vel0 * len0 / (thk0**2)

    open(unit,file=model%funits%b0dfile,position='append',form='unformatted')
      
    write(unit) &
      model%numerics%time, &
      real(f1*sum(model%geometry%thck)), &
      real(f2*count(model%geometry%thck>0.0d0)), &
      real(f2*count(model%temper%bwat>0.0d0)), &
     (real(f4*  model%velocity%uflx(model%funits%indices0dx(i),model%funits%indices0dy(i))), &
                i=1,size(model%funits%indices0dx)), & 
     (real(f4*  model%velocity%vflx(model%funits%indices0dx(i),model%funits%indices0dy(i))), &
                i=1,size(model%funits%indices0dx)), &
     (real(f5*  model%velocity%ubas(model%funits%indices0dx(i),model%funits%indices0dy(i))), &
                i=1,size(model%funits%indices0dx)), & 
     (real(f5*  model%velocity%vbas(model%funits%indices0dx(i),model%funits%indices0dy(i))), &
                i=1,size(model%funits%indices0dx)), &
     (real(f6*  model%velocity%btrc(model%funits%indices0dx(i),model%funits%indices0dy(i))), &
                i=1,size(model%funits%indices0dx)), &
     (real(thk0*model%geometry%thck(model%funits%indices0dx(i),model%funits%indices0dy(i))), &
                i=1,size(model%funits%indices0dx)), &
     (real(thk0*model%geometry%usrf(model%funits%indices0dx(i),model%funits%indices0dy(i))), &
                i=1,size(model%funits%indices0dx)), &
     (real(thk0*model%geometry%lsrf(model%funits%indices0dx(i),model%funits%indices0dy(i))), &
                i=1,size(model%funits%indices0dx)), &
     (real(thk0*model%geometry%topg(model%funits%indices0dx(i),model%funits%indices0dy(i))), &
                i=1,size(model%funits%indices0dx)), &
     (real(f3*  model%climate% acab(model%funits%indices0dx(i),model%funits%indices0dy(i))), &
                i=1,size(model%funits%indices0dx)), &
     (real(f3*  model%temper%  bmlt(model%funits%indices0dx(i),model%funits%indices0dy(i))), &
                i=1,size(model%funits%indices0dx)), &
     (real(thk0*model%temper%  bwat(model%funits%indices0dx(i),model%funits%indices0dy(i))), &
                i=1,size(model%funits%indices0dx)), &
     (real(     model%climate% artm(model%funits%indices0dx(i),model%funits%indices0dy(i))), &
                i=1,size(model%funits%indices0dx)), &
     (real(     model%temper%  temp(model%general%upn,model%funits%indices0dx(i),model%funits%indices0dy(i))),  &
                i=1,size(model%funits%indices0dx))

    close(unit)

  end subroutine writ0dvr

!-------------------------------------------------------------------------

  function rplanout(fname,unit,ewn,nsn)

    ! opens and reads planform file

    use glimmer_global, only : sp
      
    implicit none

    character(*),intent(in) :: fname
    integer,intent(in) :: unit,ewn,nsn
    real (kind = sp), dimension(ewn,nsn) :: rplanout 

    ! locals

    integer,dimension(2) :: idum   
    real (kind = sp) :: rdum     
    character :: cdum*4 
    integer :: ns,ew

    logical::there

    inquire(file=fname,exist=there)
    if ( .not. there ) then
      print*,'planform file ',fname,' not found'
      stop
    endif

#ifdef CVF
    open(unit,file=fname,form='unformatted',convert='BIG_ENDIAN')
    ! convert is a non-standard specifier and doesn't work with the
    ! Intel compiler.
    ! Substituting with

#else

    open(unit,file=fname,form='unformatted')

#endif

    read(unit,end=10) rdum, cdum, idum, rdum 
    read(unit,end=10) ((rplanout(ew,ns),ew=1,ewn),ns=1,nsn)

    close(unit)

    return

10    print*, 'read past end of planform file'
    Print*,'Reading file:',fname
    stop

  end function rplanout

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine redtsout(fname,unit,forcdata)

    !*FD Sets up forcing data from time-series file (e.g. GRIP data)

    implicit none

    character(*), intent(in)          :: fname     !*FD filename to use
    integer,      intent(in)          :: unit      !*FD File unit to use
    type(glimmer_forcdata),intent(inout) :: forcdata  !*FD Parameters to be set
  
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

end module glimmer_outp
