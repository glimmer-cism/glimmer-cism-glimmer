
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_messages.f90 - part of the GLIMMER ice model       + 
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

module glide_messages

  !*FD Module containing error/message handling code for GLIMMER.
  !*FD Six levels of message/error are defined:
  !*FD \begin{itemize}
  !*FD \item Diagnostic messages
  !*FD \item Timestep enumeration and related information
  !*FD \item Information messages
  !*FD \item Warning messages
  !*FD \item Error messages
  !*FD \item Fatal error messages
  !*FD \end{itemize}
  !*FD These are numbered 1--6, with increasing severity, and the level of
  !*FD message output may be set to output all messages, only those above a particular 
  !*FD severity, or none at all. It should be noted that even if all messages are
  !*FD turned off, the model will still halt if it encounters a fatal
  !*FD error!
  !*FD 
  !*FD The other point to note is that when calling the messaging routines,
  !*FD the numerical identifier of a message level should be replaced by the
  !*FD appropriate parameter:
  !*FD \begin{itemize}
  !*FD \item \texttt{GM\_DIAGNOSTIC}
  !*FD \item \texttt{GM\_TIMESTEP}
  !*FD \item \texttt{GM\_INFO}
  !*FD \item \texttt{GM\_WARNING}
  !*FD \item \texttt{GM\_ERROR}
  !*FD \item \texttt{GM\_FATAL}
  !*FD \end{itemize}

  integer,parameter :: GM_DIAGNOSTIC = 1 !*FD Numerical identifier for diagnostic messages.
  integer,parameter :: GM_TIMESTEP   = 2 !*FD Numerical identifier for timestep messages.
  integer,parameter :: GM_INFO       = 3 !*FD Numerical identifier for information messages.
  integer,parameter :: GM_WARNING    = 4 !*FD Numerical identifier for warning messages.
  integer,parameter :: GM_ERROR      = 5 !*FD Numerical identifier for (non-fatal) error messages.
  integer,parameter :: GM_FATAL      = 6 !*FD Numerical identifier for fatal error messages.

  integer,parameter            :: GM_levels = 6
  logical,dimension(GM_levels) :: gm_show = .true.
 
contains

  subroutine glide_msg(type,file,line,text)

    !*FD Generate a message, and output according to the value of \texttt{GM\_levels}.
    !*FD N.B. it is expected that the \texttt{\_\_FILE\_\_} and \texttt{\_\_LINE\_\_}
    !*FD preprocessor macros will be used when calling this subroutine.

    integer,intent(in)      :: type !*FD Type of error to be generated (see list above).
    character(*),intent(in) :: file !*FD Name of file containing code from which \texttt{glide\_msg}
                                    !*FD was called.
    integer,intent(in)      :: line !*FD Line number of call to \texttt{glide\_msg}.
    character(*),intent(in) :: text !*FD The line of text to be output.

    character(6) :: linetxt

    write(linetxt,'(I6)')line

    select case(type)
    case(GM_DIAGNOSTIC)
      if (gm_show(GM_DIAGNOSTIC)) write(*,*)'* ',text
    case(GM_TIMESTEP)
      if (gm_show(GM_TIMESTEP))   write(*,*)'* ',text
    case(GM_INFO)
      if (gm_show(GM_INFO))       write(*,*)'* MESSAGE: ',text
    case(GM_WARNING)
      if (gm_show(GM_WARNING))    write(*,*)'* WARNING: ',text
    case(GM_ERROR)
      if (gm_show(GM_ERROR))      write(*,*)'* ERROR: ',text
    case(GM_FATAL)
      if (gm_show(GM_FATAL)) then
        write(*,*)'* FATAL ERROR (',trim(file),':',trim(adjustl(linetxt)),') ',text
      end if
      stop
    case default
      write(*,*)'* Error in call to GLIDE_MSG, in ',trim(file),', at line ',trim(linetxt),':'
      write(*,*)'* Type ',type,' is not recognised'
    end select

  end subroutine glide_msg

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glide_stars

    !*FD Outputs a row of stars, only if all messages are enabled.

      if (gm_show(GM_DIAGNOSTIC)) write(*,*)&
        '**********************************************************************'

  end subroutine glide_stars

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glide_set_msg_level(level)

    !*FD Sets the output message level.

    integer, intent(in) :: level !*FD The message level (6 is all messages; 0 is no messages). 
    integer :: i

    do i=1,GM_levels
      if (i>(6-level)) then
        gm_show(i)=.true.
      else
        gm_show(i)=.false.
      endif
    enddo

  end subroutine glide_set_msg_level

end module glide_messages
