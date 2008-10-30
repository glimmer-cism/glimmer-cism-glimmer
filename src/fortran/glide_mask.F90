! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_mask.f90 - part of the GLIMMER ice model           + 
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

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

#include <glide_mask.inc>

#define MASK model%geometry%thkmask

module glide_mask
  !*FD masking ice thicknesses

contains
  subroutine glide_set_mask(model)
    use glide_types
    use glimmer_global, only : dp
    use glimmer_physcon, only : rhoi, rhoo
    implicit none
    type(glide_global_type) :: model        !*FD model instance

    ! local variables
    integer ew,ns
    real(dp), parameter :: con = - rhoi / rhoo

    MASK = 0
    model%geometry%iarea = 0.
    model%geometry%ivol = 0.
    do ns=1,model%general%nsn
       do ew = 1,model%general%ewn
          
          if (model%geometry%thck(ew,ns) .eq. 0.) then                               ! no ice
             if (model%geometry%topg(ew,ns) .lt. model%climate%eus) then             ! below SL
                MASK(ew,ns) = GLIDE_MASK_OCEAN
             else                                                                    ! above SL
                MASK(ew,ns) = GLIDE_MASK_LAND
             end if
          else
             model%geometry%iarea = model%geometry%iarea + 1.
             model%geometry%ivol = model%geometry%ivol + model%geometry%thck(ew,ns)
             if (model%geometry%topg(ew,ns) - model%climate%eus &                    ! ice
                  < con * model%geometry%thck(ew,ns)) then                           ! floating ice
                MASK(ew,ns) = GLIDE_MASK_SHELF
             else                                                                    ! grounded ice
                MASK(ew,ns) = GLIDE_MASK_INTERIOR
             end if
             if (model%geometry%thck(ew,ns) .le. model%numerics%thklim) then         ! ice below dynamic limit
                MASK(ew,ns) = ior(MASK(ew,ns),GLIDE_MASK_THIN_ICE)
             end if
          end if

       end do
    end do
    model%geometry%iarea = model%geometry%iarea * model%numerics%dew * model%numerics%dns
    model%geometry%ivol = model%geometry%ivol * model%numerics%dew * model%numerics%dns

    ! finding boundaries
    do ns=2,model%general%nsn-1
       do ew = 2,model%general%ewn-1
          if (GLIDE_IS_FLOAT(MASK(ew,ns))) then
             ! shelf front
             if (GLIDE_IS_OCEAN(MASK(ew-1,ns)) .or. GLIDE_IS_OCEAN(MASK(ew+1,ns)) .or. &
                  GLIDE_IS_OCEAN(MASK(ew,ns-1)) .or. GLIDE_IS_OCEAN(MASK(ew,ns+1))) then
                MASK(ew,ns) = ior(MASK(ew,ns),GLIDE_MASK_SHELF_FRONT)
             end if
          else if (GLIDE_IS_GROUND(MASK(ew,ns))) then
             ! land margin
             if (GLIDE_IS_LAND(MASK(ew-1,ns)) .or. GLIDE_IS_LAND(MASK(ew+1,ns)) .or. &
                  GLIDE_IS_LAND(MASK(ew,ns-1)) .or. GLIDE_IS_LAND(MASK(ew,ns+1))) then
                MASK(ew,ns) = ior(MASK(ew,ns),GLIDE_MASK_LAND_MARGIN)
             end if
             ! grounding line
             if (GLIDE_IS_FLOAT(MASK(ew-1,ns)) .or. &
                  GLIDE_IS_FLOAT(MASK(ew+1,ns)) .or. &
                  GLIDE_IS_FLOAT(MASK(ew,ns-1)) .or. & 
                  GLIDE_IS_FLOAT(MASK(ew,ns+1))) then
                MASK(ew,ns) = ior(MASK(ew,ns),GLIDE_MASK_GROUNDING_LINE)
             end if
          end if
          ! Edge of marine ice, whether floating or not
          if ((model%geometry%topg(ew,ns) .lt. model%climate%eus.and.&
               model%geometry%thck(ew,ns)>0.0).and. &
               (GLIDE_IS_OCEAN(MASK(ew-1,ns)) .or. GLIDE_IS_OCEAN(MASK(ew+1,ns)) .or. &
               GLIDE_IS_OCEAN(MASK(ew,ns-1)) .or. GLIDE_IS_OCEAN(MASK(ew,ns+1)))) then
             MASK(ew,ns) = ior(MASK(ew,ns),GLIDE_MASK_MARINE_EDGE)
          end if
       end do
    end do
  end subroutine glide_set_mask

end module glide_mask
