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

  private :: handle_section, handle_value, InsertSection, InsertValue

  integer, parameter :: namelen=20
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
     module procedure GetValueReal, GetValueInt, GetValueChar
  end interface

contains
  function ConfigRead(fname)
    !*FD read configuration file
    implicit none
    character(len=*), intent(in) :: fname
    !*FD name of configuration file
    type(ConfigSection), pointer :: ConfigRead

    ! local variables
    type(ConfigSection), pointer :: this_section
    type(ConfigValue), pointer ::  this_value
    logical there
    integer unit,ios,linenr
    character(len=linelen) :: line

    inquire (exist=there,file=fname)
    if (.not.there) then
       write(*,*) 'Cannot open configuration file ',trim(fname)
       stop
    end if

    open(unit,file=fname)
    ios=0
    linenr=0
    ConfigRead=>NULL()
    this_section=>NULL()
    do while(ios.eq.0)
       read(unit,fmt='(a250)',iostat=ios) line
       line = adjustl(line)
       if (.not.(line(1:1).eq.'!' .or. line(1:1).eq.'#' .or. line(1:1).eq.';' .or. line(1:1).eq.' ')) then
          ! handle comments
          if (line(1:1).eq.'[') then
             ! new section
             call handle_section(linenr,line,this_section)
             this_value=>NULL()
             if (.not.associated(ConfigRead)) then
                ! this is the first section in config file
                ConfigRead=>this_section
             end if
          else
             ! handle value
             if (.not.associated(ConfigRead)) then
                write(*,*) 'Error, no section defined yet'
                write(*,*) trim(adjustl(fname)), linenr
                stop
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
  end function ConfigRead

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

  function GetSection(config,name)
    !*FD Find and return section with name
    implicit none
    type(ConfigSection), pointer :: config
    type(ConfigSection), pointer :: GetSection
    character(len=*),intent(in) :: name

    GetSection=>config
    
    do while(associated(GetSection))
       if (name.eq.trim(GetSection%name)) then
          return
       end if
       GetSection=>GetSection%next
    end do
  end function GetSection

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

    read(value,*) temp
    if (ios==0) then
       val = temp
    end if
  end subroutine GetValueReal

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
    implicit none
    integer, intent(in) :: linenr
    character(len=*) :: line
    type(ConfigSection), pointer :: section

    ! local variables
    integer i

    do i=1,linelen
       if (line(i:i).eq.']') then
          exit
       end if
    end do
    if (line(i:i).ne.']') then
       write(*,*) 'Cannot find end of section'
       write(*,*) linenr
       stop
    end if

    section => InsertSection(trim(adjustl(line(2:i-1))),section)
  end subroutine handle_section
  
  subroutine handle_value(linenr,line,value)
    implicit none
    integer, intent(in) :: linenr
    character(len=*) :: line
    type(ConfigValue), pointer :: value

    ! local variables
    integer i

    do i=1,linelen
       if (line(i:i).eq.'=' .or. line(i:i).eq.':') then
          exit
       end if
    end do
    if (.not.(line(i:i).eq.'=' .or. line(i:i).eq.':')) then
       write(*,*) 'Cannot find = or : '
       write(*,*) linenr
       stop
    end if

    value => InsertValue(trim(adjustl(line(:i-1))), trim(adjustl(line(i+1:))),value)
  end subroutine handle_value

  function InsertSection(name,here)
    !*FD add a new section
    implicit none
    character(len=*), intent(in) :: name
    type(ConfigSection), pointer :: here, InsertSection

    allocate(InsertSection)
    InsertSection%name = name

    if (associated(here)) then
       if (associated(here%next)) then
          InsertSection%next => here%next
       end if
       here%next => InsertSection
    end if
  end function InsertSection

  function InsertValue(name,value,here)
    !*FD add a new value
    implicit none
    character(len=*), intent(in) :: name, value
    type(ConfigValue), pointer :: here, InsertValue

    allocate(InsertValue)
    
    InsertValue%name = name
    InsertValue%value = value

    if(associated(here)) then
       if (associated(here%next)) then
          InsertValue%next => here%next
       end if
       here%next => InsertValue
    end if
  end function InsertValue
end module glimmer_config
