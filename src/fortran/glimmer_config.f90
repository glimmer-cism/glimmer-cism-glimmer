
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_config.f90 - part of the GLIMMER ice model       + 
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

module glimmer_config

  use glimmer_global

  !*FD Configuration file parser,
  !*FD written by Magnus Hagdorn, May 2004.
  !*FD Everything is a singly linked list

  private :: handle_section, handle_value, InsertSection, InsertValue

  integer, parameter :: namelen=20
  integer, parameter :: valuelen=200
  integer, parameter :: linelen=250

  type ConfigValue
     !*FD `Value' element of configuration linked list
     character(len=namelen) :: name              !*FD Name of value element
     character(len=valuelen) :: value            !*FD Value of this element
     type(ConfigValue), pointer :: next=>NULL()  !*FD Pointer to next element
  end type ConfigValue

  type ConfigSection
     !*FD `Section' element of configuration linked list
     character(len=namelen) :: name               !*FD Name of this section
     type(ConfigValue), pointer :: values=>NULL() !*FD Pointer to list of values
     type(ConfigSection), pointer :: next=>NULL() !*FD Pointer to next section
  end type ConfigSection

  interface GetValue
     module procedure GetValueReal, GetValueInt, GetValueChar, GetValueRealArray, GetValueIntArray, &
          GetValueDouble, GetValueDoubleArray, GetValueCharArray
  end interface

contains

  subroutine ConfigRead(fname,config)

    !*FD Read a configuration file.

    use glide_messages
    implicit none

    character(len=*), intent(in) :: fname    !*FD name of configuration file
    type(ConfigSection), pointer :: config   !*FD pointer to first section

    ! local variables
    type(ConfigSection), pointer :: this_section
    type(ConfigValue), pointer ::  this_value
    character(40) :: errtxt
    logical there
    integer unit,ios,linenr
    character(len=linelen) :: line

    inquire (exist=there,file=fname)
    if (.not.there) then
       call glide_msg(GM_FATAL,__FILE__,__LINE__,'Cannot open configuration file '//trim(fname))
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
                write(errtxt,*)linenr
                call glide_msg(GM_FATAL,__FILE__,__LINE__,'Error, no section defined yet. '// &
                  trim(adjustl(fname))//errtxt)
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

  !--------------------------------------------------------------------------

  subroutine PrintConfig(config)

    !*FD Prints the contents of a configuration linked-list to the screen.

    use glide_messages
    implicit none

    type(ConfigSection), pointer :: config !*FD Pointer to list to be printed.

    type(ConfigSection), pointer :: sec
    type(ConfigValue),   pointer :: val
    
    sec=>config
    do while(associated(sec))
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,trim(sec%name))
       val=>sec%values
       do while(associated(val))
          call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'  '//trim(val%name)//' == '//trim(val%value))
          val=>val%next
       end do
       call glide_msg(GM_DIAGNOSTIC,__FILE__,__LINE__,'')
       sec=>sec%next
    end do
  end subroutine PrintConfig

  !--------------------------------------------------------------------------

  integer function ValidateSections(config,section_list)

    !*FD Checks to see that all the sections in a config file are also
    !*FD present in a supplied list. Prints a message returns 1 if not.
    !*FDRV 0: No errors, 1: unexpected section name found.

    use glide_messages
    implicit none

    type(ConfigSection), pointer :: config    !*FD Pointer to list to be checked
    character(*),dimension(:) :: section_list !*FD List of allowed section names

    type(ConfigSection), pointer :: sec ! Keeps track of where we are
    character(80) :: outtxt

    ValidateSections=0
    sec=>config
    do 
       if (.not.any(sec%name==section_list)) then
          write(outtxt,*)'Unrecognised section name: ',trim(sec%name)
          call glide_msg(GM_ERROR,__FILE__,__LINE__,trim(outtxt))
          ValidateSections=1
       end if
       sec=>sec%next
       if (.not.associated(sec)) exit
    end do

  end function ValidateSections

  !--------------------------------------------------------------------------

  integer function ValidateValueNames(config,section_list)

    !*FD Checks to see that all the sections in a config file are also
    !*FD present in a supplied list. Prints a message returns 1 if not.
    !*FDRV 0: No errors, 1: unexpected section name found.

    use glide_messages
    implicit none

    type(ConfigSection), pointer :: config    !*FD Pointer to list to be checked
    character(*),dimension(:) :: section_list !*FD List of allowed section names

    type(ConfigValue), pointer :: val ! Keeps track of where we are
    character(80) :: outtxt

    ValidateValueNames=0
    val=>config%values
    do 
       if (.not.any(val%name==section_list)) then
          write(outtxt,*)'Unrecognised value name: ',trim(val%name)
          call glide_msg(GM_ERROR,__FILE__,__LINE__,trim(outtxt))
          ValidateValueNames=1
       end if
       val=>val%next
       if (.not.associated(val)) exit
    end do

  end function ValidateValueNames

  !--------------------------------------------------------------------------

  subroutine GetSection(config,found,name)

    !*FD Find and return section with name. The found section is
    !*FD returned in \texttt{found}.

    implicit none

    type(ConfigSection), pointer :: config !*FD The start of the linked list
    type(ConfigSection), pointer :: found  !*FD The found section
    character(len=*),intent(in) :: name    !*FD The name of the required section

    found=>config
    do while(associated(found))
       if (name.eq.trim(found%name)) then
          return
       end if
       found=>found%next
    end do
  end subroutine GetSection

  !--------------------------------------------------------------------------

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
    integer ios,i,numv

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

  !--------------------------------------------------------------------------

  subroutine GetValueDoubleArray(section,name,val,numval)
    !*FD get real(dp) array value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    real(dp), pointer, dimension(:) :: val
    integer,intent(in), optional :: numval

    ! local variables
    character(len=valuelen) :: value
    real, dimension(:),allocatable :: tempval
    integer ios,i,numv

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

  !--------------------------------------------------------------------------

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
    integer ios,i,numv

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

  !--------------------------------------------------------------------------

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
    integer ios,i,numv

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

  !--------------------------------------------------------------------------

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

  !--------------------------------------------------------------------------

  subroutine GetValueDouble(section,name,val)
    !*FD get real(dp) value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    real(dp) :: val

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

  !--------------------------------------------------------------------------

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

  !--------------------------------------------------------------------------

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
    use glide_messages
    implicit none
    integer, intent(in) :: linenr
    character(len=*) :: line
    type(ConfigSection), pointer :: section

    ! local variables
    integer i
    character(40) :: errtxt

    do i=1,linelen
       if (line(i:i).eq.']') then
          exit
       end if
    end do
    if (line(i:i).ne.']') then
       write(errtxt,*)linenr
       call glide_msg(GM_FATAL,__FILE__,__LINE__,'Cannot find end of section'//trim(errtxt))
    end if

    call InsertSection(trim(adjustl(line(2:i-1))),section)
  end subroutine handle_section
  
  subroutine handle_value(linenr,line,value)
    use glide_messages
    implicit none
    integer, intent(in) :: linenr
    character(len=*) :: line
    type(ConfigValue), pointer :: value

    ! local variables
    integer i
    character(40) :: errtxt

    do i=1,linelen
       if (line(i:i).eq.'=' .or. line(i:i).eq.':') then
          exit
       end if
    end do
    if (.not.(line(i:i).eq.'=' .or. line(i:i).eq.':')) then
       write(errtxt,*)linenr
       call glide_msg(GM_FATAL,__FILE__,__LINE__,'Cannot find = or : '//trim(errtxt))
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
