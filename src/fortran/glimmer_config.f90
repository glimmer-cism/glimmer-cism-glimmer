! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This code is taken from the Generic Mapping Tools and was 
! converted into Fortran 90 by Ian Rutt.
!
! Original code (in C) Copyright (c) 1991-2003 by P. Wessel and W. H. F. Smith
!
! Partial translation into Fortran 90 (c) 2004 Ian C. Rutt
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
!
! The Generic Mapping Tools are maintained by Paul Wessel and 
! Walter H. F. Smith. The GMT homepage is:
!
! http://gmt.soest.hawaii.edu/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module glimmer_config
  !*FD configuration file parser
  !*FD written by Magnus Hagdorn, May 2004
  !*FD everything is a singly linked list

  use glimmer_global, only : dp
  private :: handle_section, handle_value, InsertSection, InsertValue, dp

  integer, parameter :: namelen=50
  integer, parameter :: valuelen=200
  integer, parameter :: linelen=250

  type ConfigValue
     character(len=namelen) :: name
     character(len=valuelen) :: value
     type(ConfigValue), pointer :: next=>NULL()
  end type ConfigValue

  type ConfigSection
     character(len=namelen) :: name
     type(ConfigValue), pointer :: values=>NULL()
     type(ConfigSection), pointer :: next=>NULL()
  end type ConfigSection

  interface GetValue
     module procedure GetValueDouble, GetValueReal, GetValueInt, GetValueChar, &
          GetValueDoubleArray, GetValueRealArray, GetValueIntArray, GetValueCharArray
  end interface

contains
  subroutine ConfigRead(fname,config)
    !*FD read configuration file
    use glimmer_log
    implicit none
    character(len=*), intent(in) :: fname
    !*FD name of configuration file
    type(ConfigSection), pointer :: config

    ! local variables
    type(ConfigSection), pointer :: this_section
    type(ConfigValue), pointer ::  this_value
    logical there
    integer unit,ios,linenr
    character(len=linelen) :: line
    character(len=100) :: message

    inquire (exist=there,file=fname)
    if (.not.there) then
       call write_log('Cannot open configuration file '//trim(fname),GM_FATAL)
    end if
    
    unit=99
    open(unit,file=trim(fname),status='old')
    ios=0
    linenr=0
    config=>NULL()
    this_section=>NULL()
    do while(ios.eq.0)
       read(unit,fmt='(a250)',iostat=ios) line
       line = adjustl(line)
       if (ios.ne.0) then
          exit
       end if
       if (.not.(line(1:1).eq.'!' .or. line(1:1).eq.'#' .or. line(1:1).eq.';' .or. line(1:1).eq.' ')) then
          ! handle comments
          if (line(1:1).eq.'[') then
             ! new section
             call handle_section(linenr,line,this_section)
             this_value=>NULL()
             if (.not.associated(config)) then
                ! this is the first section in config file
                config=>this_section
             end if
          else
             ! handle value
             if (.not.associated(this_section)) then
                call write_log('No section defined yet',GM_ERROR)
                write(message,*) trim(adjustl(fname)), linenr
                call write_log(message,GM_FATAL)
             end if
             call handle_value(linenr,line,this_value)
             if (.not.associated(this_section%values)) then
                this_section%values => this_value
             end if
          end if
       end if
       linenr = linenr + 1
    end do
    close(unit)
    return
  end subroutine ConfigRead

  subroutine PrintConfig(config)
    implicit none
    type(ConfigSection), pointer :: config

    type(ConfigSection), pointer :: sec
    type(ConfigValue), pointer ::  val
    
    sec=>config
    do while(associated(sec))
       write(*,*) sec%name
       val=>sec%values
       do while(associated(val))
          write(*,*) '  ',trim(val%name),' == ', trim(val%value)
          val=>val%next
       end do
       write(*,*)
       sec=>sec%next
    end do
  end subroutine PrintConfig

  subroutine GetSection(config,found,name)
    !*FD Find and return section with name
    implicit none
    type(ConfigSection), pointer :: config
    type(ConfigSection), pointer :: found
    character(len=*),intent(in) :: name

    found=>config
    do while(associated(found))
       if (name.eq.trim(found%name)) then
          return
       end if
       found=>found%next
    end do
  end subroutine GetSection

  subroutine GetValueDoubleArray(section,name,val,numval)
    !*FD get real array value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    real(kind=dp), pointer, dimension(:) :: val
    integer,intent(in), optional :: numval

    ! local variables
    character(len=valuelen) :: value
    real(kind=dp), dimension(:),allocatable :: tempval
    integer i,numv

    if (present(numval)) then
       numv=numval
    else
       numv=100
    end if
    allocate(tempval(numv))
    value=''
    call GetValueChar(section,name,value)
    if (value.eq.'') return
    read(value,*,end=10) (tempval(i),i=1,numv)
10  i=i-1
    if (i.ge.1) then
       if (associated(val)) then
          deallocate(val)
       end if
       allocate(val(i))
       val = tempval(1:i)
    end if
  end subroutine GetValueDoubleArray

  subroutine GetValueRealArray(section,name,val,numval)
    !*FD get real array value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    real, pointer, dimension(:) :: val
    integer,intent(in), optional :: numval

    ! local variables
    character(len=valuelen) :: value
    real, dimension(:),allocatable :: tempval
    integer i,numv

    if (present(numval)) then
       numv=numval
    else
       numv=100
    end if
    allocate(tempval(numv))
    value=''
    call GetValueChar(section,name,value)
    if (value.eq.'') return
    read(value,*,end=10) (tempval(i),i=1,numv)
10  i=i-1
    if (i.ge.1) then
       if (associated(val)) then
          deallocate(val)
       end if
       allocate(val(i))
       val = tempval(1:i)
    end if
  end subroutine GetValueRealArray

  subroutine GetValueIntArray(section,name,val,numval)
    !*FD get integer array value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    integer, pointer, dimension(:) :: val
    integer,intent(in), optional :: numval

    ! local variables
    character(len=valuelen) :: value
    integer, dimension(:),allocatable :: tempval
    integer i,numv

    if (present(numval)) then
       numv=numval
    else
       numv=100
    end if
    allocate(tempval(numv))
    value=''
    call GetValueChar(section,name,value)
    if (value.eq.'') return
    read(value,*,end=20) (tempval(i),i=1,numv)
20  i=i-1

    if (i.ge.1) then
       if (associated(val)) then
          deallocate(val)
       end if
       allocate(val(i))
       val = tempval(1:i)
    end if
  end subroutine GetValueIntArray

  subroutine GetValueCharArray(section,name,val,numval)
    !*FD get character array value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    character(*), pointer, dimension(:) :: val
    integer,intent(in), optional :: numval

    ! local variables
    character(len=valuelen) :: value
    character(80), dimension(:),allocatable :: tempval
    integer i,numv

    if (present(numval)) then
       numv=numval
    else
       numv=100
    end if
    allocate(tempval(numv))
    value=''
    call GetValueChar(section,name,value)
    if (value.eq.'') return
    read(value,*,end=20) (tempval(i),i=1,numv)
20  i=i-1

    if (i.ge.1) then
       if (associated(val)) then
          deallocate(val)
       end if
       allocate(val(i))
       val = tempval(1:i)
    end if
  end subroutine GetValueCharArray

  subroutine GetValueReal(section,name,val)
    !*FD get real value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    real :: val

    ! local variables
    character(len=valuelen) :: value
    real temp
    integer ios

    value=''
    call GetValueChar(section,name,value)

    read(value,*,iostat=ios) temp
    if (ios==0) then
       val = temp
    end if
  end subroutine GetValueReal

  subroutine GetValueDouble(section,name,val)
    !*FD get double value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    real(kind=dp) :: val

    ! local variables
    character(len=valuelen) :: value
    real temp
    integer ios

    value=''
    call GetValueChar(section,name,value)

    read(value,*,iostat=ios) temp
    if (ios==0) then
       val = temp
    end if
  end subroutine GetValueDouble

  subroutine GetValueInt(section,name,val)
    !*FD get integer value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    integer :: val

    ! local variables
    character(len=valuelen) :: value
    integer temp
    integer ios

    value=''
    call GetValueChar(section,name,value)

    read(value,*,iostat=ios) temp
    if (ios==0) then
       val = temp
    end if
  end subroutine GetValueInt

  subroutine GetValueChar(section,name,val)
    !*FD get character value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    character(len=*) :: val
    
    type(ConfigValue), pointer :: value

    value=>section%values
    do while(associated(value))
       if (name.eq.trim(value%name)) then
          val = value%value
          return
       end if
       value=>value%next
    end do
  end subroutine GetValueChar

  !==================================================================================
  ! private procedures
  !==================================================================================

  subroutine handle_section(linenr,line,section)
    use glimmer_log
    implicit none
    integer, intent(in) :: linenr
    character(len=*) :: line
    type(ConfigSection), pointer :: section

    ! local variables
    integer i
    character(len=100) :: message

    do i=1,linelen
       if (line(i:i).eq.']') then
          exit
       end if
    end do
    if (line(i:i).ne.']') then
       write(message,*) 'Cannot find end of section ',linenr
       call write_log(message,GM_FATAL)
    end if

    call InsertSection(trim(adjustl(line(2:i-1))),section)
  end subroutine handle_section
  
  subroutine handle_value(linenr,line,value)
    use glimmer_log
    implicit none
    integer, intent(in) :: linenr
    character(len=*) :: line
    type(ConfigValue), pointer :: value

    ! local variables
    integer i
    character(len=100) :: message
    do i=1,linelen
       if (line(i:i).eq.'=' .or. line(i:i).eq.':') then
          exit
       end if
    end do
    if (.not.(line(i:i).eq.'=' .or. line(i:i).eq.':')) then
       write(message,*) 'Cannot find = or : ',linenr
       call write_log(message,GM_FATAL)
    end if

    call InsertValue(trim(adjustl(line(:i-1))), trim(adjustl(line(i+1:))),value)
  end subroutine handle_value

  subroutine InsertSection(name,section)
    !*FD add a new section
    implicit none
    character(len=*), intent(in) :: name
    type(ConfigSection), pointer :: section
    type(ConfigSection), pointer :: new_sec

    allocate(new_sec)
    new_sec%name = name

    if (associated(section)) then
       if (associated(section%next)) then
          new_sec%next => section%next
       end if
       section%next=>new_sec
    end if
    section=>new_sec
  end subroutine InsertSection

  subroutine InsertValue(name,val,value)
    !*FD add a new value
    implicit none
    character(len=*), intent(in) :: name, val
    type(ConfigValue), pointer :: value
    type(ConfigValue), pointer :: new_value

    allocate(new_value)
    
    new_value%name = name
    new_value%value = val

    if(associated(value)) then
       if (associated(value%next)) then
          new_value%next => value%next
       end if
       value%next => new_value
    end if
    value=>new_value
  end subroutine InsertValue
end module glimmer_config
