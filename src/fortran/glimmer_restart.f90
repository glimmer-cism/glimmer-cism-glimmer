
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_restart.f90 - part of the GLIMMER ice model      + 
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

module glimmer_restart

  !*FD Reads and writes parts of restart
  !*FD files pertaining to individual ice model instances

contains

  subroutine glimmer_i_write_restart(instance,unit)

    !*FD Writes a restart for an ice model instance. Note that
    !*FD the logical file unit to be used must be open
    !*FD on entry, and remains open on exit.

    use glimmer_object

    implicit none

    type(glimmer_instance),intent(in) :: instance !*FD The ice model instance to be written
    integer,            intent(in) :: unit     !*FD The logical file unit to use

    ! temporary variable, this is need to fix an internal compiler error of the SUN WS f95 compiler.
    logical :: temp

    ! top level variables

    call proj_write_restart(instance%proj,unit)
    call downs_write_restart(instance%downs,unit)

    ! Domain size

    write(unit) instance%model%general% ewn
    write(unit) instance%model%general% nsn
    write(unit) instance%model%general% upn
    
    ! Remaining top level variables

    write(unit) instance% paramfile
    write(unit) instance% newtemps
    write(unit) instance% xwind
    write(unit) instance% ywind
    write(unit) instance% global_orog
    write(unit) instance% local_orog
    write(unit) instance% gboxx
    write(unit) instance% gboxy
    write(unit) instance% gboxn
    write(unit) instance% frac_coverage
    write(unit) instance% first

    ! Now the model components:

    ! options

    write(unit) instance%model%options% whichtemp
    write(unit) instance%model%options% whichartm
    write(unit) instance%model%options% whichthck
    write(unit) instance%model%options% whichflwa
    write(unit) instance%model%options% whichisot
    write(unit) instance%model%options% whichslip
    write(unit) instance%model%options% whichbwat
    write(unit) instance%model%options% whichmarn
    write(unit) instance%model%options% whichbtrc
    write(unit) instance%model%options% whichacab
    write(unit) instance%model%options% whichstrs
    write(unit) instance%model%options% whichevol
    write(unit) instance%model%options% whichwvel
    write(unit) instance%model%options% whichprecip

    ! geometry

    write(unit) instance%model%geometry% thck
    write(unit) instance%model%geometry% usrf
    write(unit) instance%model%geometry% lsrf
    write(unit) instance%model%geometry% topg
    write(unit) instance%model%geometry% relx 
    write(unit) instance%model%geometry% mask 
    write(unit) instance%model%geometry% dom
    write(unit) instance%model%geometry% totpts
    write(unit) instance%model%geometry% empty

    ! geomderv

    write(unit) instance%model%geomderv% dthckdew
    write(unit) instance%model%geomderv% dusrfdew
    write(unit) instance%model%geomderv% dthckdns
    write(unit) instance%model%geomderv% dusrfdns 
    write(unit) instance%model%geomderv% dthckdtm
    write(unit) instance%model%geomderv% dusrfdtm
    write(unit) instance%model%geomderv% stagthck

    ! velocity

    write(unit) instance%model%velocity% uvel
    write(unit) instance%model%velocity% vvel
    write(unit) instance%model%velocity% wvel
    write(unit) instance%model%velocity% wgrd
    write(unit) instance%model%velocity% uflx
    write(unit) instance%model%velocity% vflx
    write(unit) instance%model%velocity% diffu
    write(unit) instance%model%velocity% ubas
    write(unit) instance%model%velocity% vbas
    write(unit) instance%model%velocity% btrc

    ! climate

    write(unit) instance%model%climate% acab
    write(unit) instance%model%climate% artm
    write(unit) instance%model%climate% arng
    write(unit) instance%model%climate% lati
    write(unit) instance%model%climate% ablt
    write(unit) instance%model%climate% prcp
    write(unit) instance%model%climate% uprecip_rate
    write(unit) instance%model%climate% usurftemp
    write(unit) instance%model%climate% ice_albedo
    write(unit) instance%model%climate% ulapse_rate

    ! temper

    write(unit) instance%model%temper% temp
    write(unit) instance%model%temper% flwa  
    write(unit) instance%model%temper% bwat
    write(unit) instance%model%temper% bmlt
    write(unit) instance%model%temper% niter
    write(unit) instance%model%temper% perturb
    write(unit) instance%model%temper% grid
    write(unit) instance%model%temper% tpt
    write(unit) instance%model%temper% first1

    ! forcdata

    temp = associated(instance%model%forcdata%forcing)
    write(unit) temp
    write(unit) instance%model%forcdata% flines
    if (temp) write(unit) instance%model%forcdata% forcing

    ! funits

    write(unit) instance%model%funits% usrffile
    write(unit) instance%model%funits% topgfile
    write(unit) instance%model%funits% relxfile
    write(unit) instance%model%funits% sigfile
    write(unit) instance%model%funits% output_stem
    write(unit) instance%model%funits% ncfile

    ! numerics

    write(unit) instance%model%numerics% time
    write(unit) instance%model%numerics% tinc
    write(unit) instance%model%numerics% ntem
    write(unit) instance%model%numerics% nvel
    write(unit) instance%model%numerics% niso
    write(unit) instance%model%numerics% nout
    write(unit) instance%model%numerics% nstr
    write(unit) instance%model%numerics% alpha
    write(unit) instance%model%numerics% alphas
    write(unit) instance%model%numerics% thklim
    write(unit) instance%model%numerics% mlimit
    write(unit) instance%model%numerics% dew
    write(unit) instance%model%numerics% dns
    write(unit) instance%model%numerics% dt       
    write(unit) instance%model%numerics% dttem 
    write(unit) instance%model%numerics% sigma 
    write(unit) instance%model%numerics% nshlf 

    ! pddcalc

    write(unit) instance%model%pddcalc% dx
    write(unit) instance%model%pddcalc% dy
    write(unit) instance%model%pddcalc% ix
    write(unit) instance%model%pddcalc% iy
    write(unit) instance%model%pddcalc% nx
    write(unit) instance%model%pddcalc% ny
    write(unit) instance%model%pddcalc% pddtab
    write(unit) instance%model%pddcalc% pt_alloc
    write(unit) instance%model%pddcalc% dailytemp
    write(unit) instance%model%pddcalc% tma
    write(unit) instance%model%pddcalc% tmj
    write(unit) instance%model%pddcalc% dtmj
    write(unit) instance%model%pddcalc% dd_sigma
    write(unit) instance%model%pddcalc% first

    ! isotwk

    temp = associated(instance%model%isotwk%dflct)
    write(unit) temp
    write(unit) instance%model%isotwk% load
    write(unit) instance%model%isotwk% fact
    write(unit) instance%model%isotwk% first1
    write(unit) instance%model%isotwk% nflx
    if (temp) write(unit) instance%model%isotwk%dflct
    write(unit) instance%model%isotwk% first2

    ! velowk

    temp = associated(instance%model%velowk% dupsw)
    write(unit) instance%model%velowk% depth
    write(unit) temp
    if (temp) write(unit) instance%model%velowk% dupsw

    temp = associated(instance%model%velowk% depthw)
    write(unit) temp
    if (temp) write(unit) instance%model%velowk% depthw

    temp = associated(instance%model%velowk% suvel)
    write(unit) temp
    if (temp) write(unit) instance%model%velowk% suvel

    temp = associated(instance%model%velowk% svvel)
    write(unit) temp
    if (temp) write(unit) instance%model%velowk% svvel

    temp = associated(instance%model%velowk% fslip)
    write(unit) temp
    if (temp) write(unit) instance%model%velowk% fslip

    temp = associated(instance%model%velowk% dintflwa)
    write(unit) temp
    if (temp) write(unit) instance%model%velowk% dintflwa

    temp = associated(instance%model%velowk% dups)
    write(unit) temp
    if (temp) write(unit) instance%model%velowk% dups

    write(unit) instance%model%velowk% first1
    write(unit) instance%model%velowk% first2
    write(unit) instance%model%velowk% first3
    write(unit) instance%model%velowk% first4
    write(unit) instance%model%velowk% fact
    write(unit) instance%model%velowk% first5
    write(unit) instance%model%velowk% watwd 
    write(unit) instance%model%velowk% watct
    write(unit) instance%model%velowk% trc0
    write(unit) instance%model%velowk% trcmin
    write(unit) instance%model%velowk% marine
    write(unit) instance%model%velowk% trcmax
    write(unit) instance%model%velowk% c 
    write(unit) instance%model%velowk% first6

   ! pcgdwk

    write(unit) instance%model%pcgdwk% pcgsize
    write(unit) instance%model%pcgdwk% ct
    write(unit) instance%model%pcgdwk% fc
    write(unit) instance%model%pcgdwk% fc2
    write(unit) instance%model%pcgdwk% mlinit
    write(unit) instance%model%pcgdwk% tlinit
    write(unit) instance%model%pcgdwk% first1

    ! thckwk

    write(unit) instance%model%thckwk% nwhich
    write(unit) instance%model%thckwk% olds
    write(unit) instance%model%thckwk% oldtime

    temp = associated(instance%model%thckwk%oldthck)
    write(unit) temp
    if (temp) write(unit) instance%model%thckwk%oldthck
    temp = associated(instance%model%thckwk% basestate)
    write(unit) temp
    if (temp) write(unit) instance%model%thckwk% basestate 

    write(unit) instance%model%thckwk% few
    write(unit) instance%model%thckwk% fns
    write(unit) instance%model%thckwk% first1

   ! tempwk
    temp = associated(instance%model%tempwk% inittemp)
    write(unit) temp
    if (temp) write(unit) instance%model%tempwk% inittemp
    temp = associated(instance%model%tempwk% dupa)
    write(unit) temp
    if (temp) write(unit) instance%model%tempwk% dupa
    temp = associated(instance%model%tempwk% dupb)
    write(unit) temp
    if (temp) write(unit) instance%model%tempwk% dupb
    temp = associated(instance%model%tempwk% dupc)
    write(unit) temp
    if (temp) write(unit) instance%model%tempwk% dupc
    temp = associated(instance%model%tempwk% c1)
    write(unit) temp
    if (temp) write(unit) instance%model%tempwk% c1
    temp = associated(instance%model%tempwk% dups)
    write(unit) temp
    if (temp) write(unit) instance%model%tempwk% dups

    write(unit) instance%model%tempwk% advconst
    write(unit) instance%model%tempwk% noflow
    write(unit) instance%model%tempwk% first1
    write(unit) instance%model%tempwk% cons 
    write(unit) instance%model%tempwk% zbed
    write(unit) instance%model%tempwk% dupn
    write(unit) instance%model%tempwk% wmax 
    write(unit) instance%model%tempwk% first2
    write(unit) instance%model%tempwk% first3
    write(unit) instance%model%tempwk% f
    write(unit) instance%model%tempwk% first4
    write(unit) instance%model%tempwk% dt_wat
    write(unit) instance%model%tempwk% nwat
    write(unit) instance%model%tempwk% c
    write(unit) instance%model%tempwk% first5

   ! paramets

    write(unit) instance%model%paramets% geot
    write(unit) instance%model%paramets% fiddle
    write(unit) instance%model%paramets% airt
    write(unit) instance%model%paramets% nmsb           
    write(unit) instance%model%paramets% hydtim
    write(unit) instance%model%paramets% isotim
    write(unit) instance%model%paramets% bpar

  end subroutine glimmer_i_write_restart

  subroutine glimmer_i_read_restart(instance,unit,nxg,nyg)

    !*FD Reads a restart for an ice model instance.
    !*FD Note that the logical file unit to be used must be open
    !*FD on entry, and remains open on exit.

    use glimmer_object
    use glimmer_setup

    implicit none

    type(glimmer_instance),intent(inout) :: instance !*FD The ice model instance to be written
    integer,            intent(in)    :: unit     !*FD The logical file unit to use

    integer :: nxg,nyg  !*FD The size of the global fields
    logical :: force_flag,dflct_flag,flag

    ! First deallocate everything in sight...
    ! Removed as null initialisation isn't supported by f90

!    if (associated(instance%model%forcdata%forcing)) deallocate(instance%model%forcdata%forcing)
!    if (associated(instance%model%pddcalc%pddtab))   deallocate(instance%model%pddcalc%pddtab)
!    if (associated(instance%model%isotwk%dflct))     deallocate(instance%model%isotwk%dflct)
!    if (associated(instance%model%velowk%depth))     deallocate(instance%model%velowk%depth)
!    if (associated(instance%model%velowk%dupsw))     deallocate(instance%model%velowk%dupsw)
!    if (associated(instance%model%velowk%depthw))    deallocate(instance%model%velowk%depthw)
!    if (associated(instance%model%velowk%suvel))     deallocate(instance%model%velowk%suvel)
!    if (associated(instance%model%velowk%svvel))     deallocate(instance%model%velowk%svvel)
!    if (associated(instance%model%velowk%fslip))     deallocate(instance%model%velowk%fslip)
!    if (associated(instance%model%velowk%dintflwa))  deallocate(instance%model%velowk%dintflwa)
!    if (associated(instance%model%velowk%dups))      deallocate(instance%model%velowk%dups)
!    if (associated(instance%model%thckwk%oldthck))   deallocate(instance%model%thckwk%oldthck)
!    if (associated(instance%model%thckwk%basestate)) deallocate(instance%model%thckwk%basestate)
!    if (associated(instance%model%tempwk% inittemp)) deallocate(instance%model%tempwk% inittemp)
!    if (associated(instance%model%tempwk% dupa))     deallocate(instance%model%tempwk% dupa)
!    if (associated(instance%model%tempwk%dupb))      deallocate(instance%model%tempwk%dupb)
!    if (associated(instance%model%tempwk%dupc))      deallocate(instance%model%tempwk%dupc)
!    if (associated(instance%model%tempwk%c1))        deallocate(instance%model%tempwk%c1)
!    if (associated(instance%model%tempwk%dups))      deallocate(instance%model%tempwk%dups)

    call proj_read_restart(instance%proj,unit)
    call downs_read_restart(instance%downs,unit)

    ! Domain size

    read(unit) instance%model%general% ewn
    read(unit) instance%model%general% nsn
    read(unit) instance%model%general% upn

    ! Allocate arrays

    call glimmer_i_allocate(instance,nxg,nyg)

    ! top level variables

    read(unit) instance% paramfile
    read(unit) instance% newtemps
    read(unit) instance% xwind
    read(unit) instance% ywind
    read(unit) instance% global_orog
    read(unit) instance% local_orog
    read(unit) instance% gboxx
    read(unit) instance% gboxy
    read(unit) instance% gboxn
    read(unit) instance% frac_coverage
    read(unit) instance% first

    ! Now the model components:

    ! options

    read(unit) instance%model%options% whichtemp
    read(unit) instance%model%options% whichartm
    read(unit) instance%model%options% whichthck
    read(unit) instance%model%options% whichflwa
    read(unit) instance%model%options% whichisot
    read(unit) instance%model%options% whichslip
    read(unit) instance%model%options% whichbwat
    read(unit) instance%model%options% whichmarn
    read(unit) instance%model%options% whichbtrc
    read(unit) instance%model%options% whichacab
    read(unit) instance%model%options% whichstrs
    read(unit) instance%model%options% whichevol
    read(unit) instance%model%options% whichwvel
    read(unit) instance%model%options% whichprecip

    ! geometry

    read(unit) instance%model%geometry% thck
    read(unit) instance%model%geometry% usrf
    read(unit) instance%model%geometry% lsrf
    read(unit) instance%model%geometry% topg
    read(unit) instance%model%geometry% relx 
    read(unit) instance%model%geometry% mask 
    read(unit) instance%model%geometry% dom
    read(unit) instance%model%geometry% totpts
    read(unit) instance%model%geometry% empty

    ! geomderv

    read(unit) instance%model%geomderv% dthckdew
    read(unit) instance%model%geomderv% dusrfdew
    read(unit) instance%model%geomderv% dthckdns
    read(unit) instance%model%geomderv% dusrfdns 
    read(unit) instance%model%geomderv% dthckdtm
    read(unit) instance%model%geomderv% dusrfdtm
    read(unit) instance%model%geomderv% stagthck

    ! velocity

    read(unit) instance%model%velocity% uvel
    read(unit) instance%model%velocity% vvel
    read(unit) instance%model%velocity% wvel
    read(unit) instance%model%velocity% wgrd
    read(unit) instance%model%velocity% uflx
    read(unit) instance%model%velocity% vflx
    read(unit) instance%model%velocity% diffu
    read(unit) instance%model%velocity% ubas
    read(unit) instance%model%velocity% vbas
    read(unit) instance%model%velocity% btrc

    ! climate

    read(unit) instance%model%climate% acab
    read(unit) instance%model%climate% artm
    read(unit) instance%model%climate% arng
    read(unit) instance%model%climate% lati
    read(unit) instance%model%climate% ablt
    read(unit) instance%model%climate% prcp
    read(unit) instance%model%climate% uprecip_rate
    read(unit) instance%model%climate% usurftemp
    read(unit) instance%model%climate% ice_albedo
    read(unit) instance%model%climate% ulapse_rate

    ! temper

    read(unit) instance%model%temper% temp
    read(unit) instance%model%temper% flwa  
    read(unit) instance%model%temper% bwat
    read(unit) instance%model%temper% bmlt
    read(unit) instance%model%temper% niter
    read(unit) instance%model%temper% perturb
    read(unit) instance%model%temper% grid
    read(unit) instance%model%temper% tpt
    read(unit) instance%model%temper% first1
 
    ! forcdata

    read(unit) force_flag
    read(unit) instance%model%forcdata% flines

    if (force_flag) then
      allocate(instance%model%forcdata%forcing(instance%model%forcdata%flines,2)) 
      read(unit) instance%model%forcdata% forcing
    endif

    ! funits

    read(unit) instance%model%funits% usrffile
    read(unit) instance%model%funits% topgfile
    read(unit) instance%model%funits% relxfile
    read(unit) instance%model%funits% sigfile
    read(unit) instance%model%funits% output_stem
    read(unit) instance%model%funits% ncfile

    ! numerics

    read(unit) instance%model%numerics% time
    read(unit) instance%model%numerics% tinc
    read(unit) instance%model%numerics% ntem
    read(unit) instance%model%numerics% nvel
    read(unit) instance%model%numerics% niso
    read(unit) instance%model%numerics% nout
    read(unit) instance%model%numerics% nstr
    read(unit) instance%model%numerics% alpha
    read(unit) instance%model%numerics% alphas
    read(unit) instance%model%numerics% thklim
    read(unit) instance%model%numerics% mlimit
    read(unit) instance%model%numerics% dew
    read(unit) instance%model%numerics% dns
    read(unit) instance%model%numerics% dt       
    read(unit) instance%model%numerics% dttem 
    read(unit) instance%model%numerics% sigma 
    read(unit) instance%model%numerics% nshlf 

    ! pddcalc

    read(unit) instance%model%pddcalc% dx
    read(unit) instance%model%pddcalc% dy
    read(unit) instance%model%pddcalc% ix
    read(unit) instance%model%pddcalc% iy
    read(unit) instance%model%pddcalc% nx
    read(unit) instance%model%pddcalc% ny

    allocate(instance%model%pddcalc%pddtab(instance%model%pddcalc%nx, &
                                           instance%model%pddcalc%ny))

    read(unit) instance%model%pddcalc% pddtab
    read(unit) instance%model%pddcalc% pt_alloc
    read(unit) instance%model%pddcalc% dailytemp
    read(unit) instance%model%pddcalc% tma
    read(unit) instance%model%pddcalc% tmj
    read(unit) instance%model%pddcalc% dtmj
    read(unit) instance%model%pddcalc% dd_sigma
    read(unit) instance%model%pddcalc% first

    ! isotwk

    read(unit) dflct_flag
    read(unit) instance%model%isotwk% load
    read(unit) instance%model%isotwk% fact
    read(unit) instance%model%isotwk% first1
    read(unit) instance%model%isotwk% nflx

    if (dflct_flag) then 
      allocate(instance%model%isotwk%dflct &
                (instance%model%isotwk%nflx, &
                 instance%model%isotwk%nflx))
      read(unit) instance%model%isotwk% dflct
    endif

    read(unit) instance%model%isotwk% first2

    ! velowk

    allocate(instance%model%velowk%depth(instance%model%general%upn))
    read(unit) instance%model%velowk% depth

    read(unit) flag
    if (flag) then
      allocate(instance%model%velowk%dupsw(instance%model%general%upn))
      read(unit) instance%model%velowk% dupsw
    endif

    read(unit) flag
    if (flag) then
      allocate(instance%model%velowk%depthw(instance%model%general%upn))    
      read(unit) instance%model%velowk% depthw
    endif

    read(unit) flag
    if (flag) then
      allocate(instance%model%velowk%suvel(instance%model%general%upn))
      read(unit) instance%model%velowk% suvel
    endif

    read(unit) flag
    if (flag) then
      allocate(instance%model%velowk%svvel(instance%model%general%upn))
      read(unit) instance%model%velowk% svvel
    endif

    read(unit) flag
    if (flag) then
      allocate(instance%model%velowk%fslip &
                (instance%model%general%ewn, &
                 instance%model%general%nsn))
      read(unit) instance%model%velowk% fslip
    endif

    read(unit) flag
    if (flag) then
      allocate(instance%model%velowk%dintflwa &
                 (instance%model%general%ewn, &
                  instance%model%general%nsn))    
      read(unit) instance%model%velowk% dintflwa
    endif

    read(unit) flag
    if (flag) then
      allocate(instance%model%velowk%dups(instance%model%general%upn)) 
      read(unit) instance%model%velowk% dups     
    endif   

    read(unit) instance%model%velowk% first1
    read(unit) instance%model%velowk% first2
    read(unit) instance%model%velowk% first3
    read(unit) instance%model%velowk% first4
    read(unit) instance%model%velowk% fact
    read(unit) instance%model%velowk% first5
    read(unit) instance%model%velowk% watwd
    read(unit) instance%model%velowk% watct
    read(unit) instance%model%velowk% trc0
    read(unit) instance%model%velowk% trcmin
    read(unit) instance%model%velowk% marine
    read(unit) instance%model%velowk% trcmax
    read(unit) instance%model%velowk% c 
    read(unit) instance%model%velowk% first6

   ! pcgdwk

    read(unit) instance%model%pcgdwk% pcgsize
    read(unit) instance%model%pcgdwk% ct
    read(unit) instance%model%pcgdwk% fc
    read(unit) instance%model%pcgdwk% fc2
    read(unit) instance%model%pcgdwk% mlinit
    read(unit) instance%model%pcgdwk% tlinit
    read(unit) instance%model%pcgdwk% first1

    ! thckwk

    read(unit) instance%model%thckwk% nwhich
    read(unit) instance%model%thckwk% olds
    read(unit) instance%model%thckwk% oldtime


    read(unit) flag
    if(flag) then
      allocate(instance%model%thckwk%oldthck &
                (instance%model%general%ewn, &
                 instance%model%general%nsn))
      read(unit) instance%model%thckwk% oldthck
    endif

    read(unit) flag
    if(flag) then
      allocate(instance%model%thckwk%basestate &
                (instance%model%general%ewn, &
                 instance%model%general%nsn))
      read(unit) instance%model%thckwk% basestate 
    endif

    read(unit) instance%model%thckwk% few
    read(unit) instance%model%thckwk% fns
    read(unit) instance%model%thckwk% first1

   ! tempwk

    read(unit) flag
    if (flag) then
      allocate(instance%model%tempwk% inittemp(instance%model%general%upn,&
                                               instance%model%general%ewn,&
                                               instance%model%general%nsn))
      read(unit) instance%model%tempwk% inittemp
    endif

    read(unit) flag
    if (flag) then
      allocate(instance%model%tempwk% dupa(instance%model%general%upn))
      read(unit) instance%model%tempwk% dupa
    endif

    read(unit) flag
    if (flag) then
      allocate(instance%model%tempwk%dupb(instance%model%general%upn))
      read(unit)instance%model%tempwk%dupb
    endif

    read(unit) flag
    if (flag) then
      allocate(instance%model%tempwk%dupc(instance%model%general%upn))
      read(unit) instance%model%tempwk%dupc
    endif

    read(unit) flag
    if (flag) then
      allocate(instance%model%tempwk%c1(instance%model%general%upn))
      read(unit) instance%model%tempwk%c1
    endif

    read(unit) flag
    if (flag) then
      allocate(instance%model%tempwk%dups(instance%model%general%upn,3))
      read(unit) instance%model%tempwk%dups
    endif

    read(unit) instance%model%tempwk% advconst
    read(unit) instance%model%tempwk% noflow
    read(unit) instance%model%tempwk% first1
    read(unit) instance%model%tempwk% cons 
    read(unit) instance%model%tempwk% zbed
    read(unit) instance%model%tempwk% dupn
    read(unit) instance%model%tempwk% wmax 
    read(unit) instance%model%tempwk% first2
    read(unit) instance%model%tempwk% first3
    read(unit) instance%model%tempwk% f
    read(unit) instance%model%tempwk% first4
    read(unit) instance%model%tempwk% dt_wat
    read(unit) instance%model%tempwk% nwat
    read(unit) instance%model%tempwk% c
    read(unit) instance%model%tempwk% first5

   ! paramets

    read(unit) instance%model%paramets% geot
    read(unit) instance%model%paramets% fiddle
    read(unit) instance%model%paramets% airt
    read(unit) instance%model%paramets% nmsb           
    read(unit) instance%model%paramets% hydtim
    read(unit) instance%model%paramets% isotim
    read(unit) instance%model%paramets% bpar

  end subroutine glimmer_i_read_restart

end module glimmer_restart
