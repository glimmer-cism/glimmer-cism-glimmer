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
    use glimmer_global, only : dp, NaN
  !*FD masking ice thicknesses

contains
  subroutine glide_set_mask(model)
    use glide_types
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

    subroutine glide_marine_margin_normal(thck, mask, marine_bc_normal)
        !*FD This subroutine derives from the given mask the normal to an ice shelf
        !*FD each point on the marine margin.
        real(dp), dimension(:,:), intent(in) :: thck
        integer, dimension(:,:), intent(in) :: mask
        real(dp), dimension(:,:), intent(out) :: marine_bc_normal

        integer :: i, j

        do i = 1, size(mask, 1)
            do j = 1, size(mask, 2)
                if (GLIDE_IS_MARINE_ICE_EDGE(mask(i,j))) then
                    marine_bc_normal(i,j) = calc_normal_45deg(thck(i-1:i+1,j-1:j+1))
                else
                    marine_bc_normal(i,j) = NaN
                end if
            end do
        end do

    end subroutine

    function calc_normal_45deg(thck3x3)
        !*FD Computes the angle of the normal vector, in radians, for the given
        !*FD 3x3 segment of ice geometry.
        !*FD The normal is given in increments of 45 degrees (no nicer
        !*FD interpolation is currently done)
        !*FD This is based on the Payne and Price GLAM code, if/when this is
        !*FD integrated into CISM it should probably be refactored to use this.
        real(dp), dimension(3,3) :: thck3x3

        real(dp) :: calc_normal_45deg
         
        real (kind = dp), dimension(3,3) :: mask, maskcorners
        real (kind = dp), dimension(3,3) :: thckmask
        real (kind = dp), dimension(3) :: testvect
        real (kind = dp) :: phi, deg2rad

        deg2rad = 3.141592654d0 / 180.0d0
        loc_latbc = 0; phi = 0
        mask(:,1) = (/ 0.0d0, 180.0d0, 0.0d0 /)
        mask(:,2) = (/ 270.0d0, 0.0d0, 90.0d0 /)
        mask(:,3) = (/ 0.0d0, 360.0d0, 0.0d0 /)
        maskcorners(:,1) = (/ 225.0d0, 0.0d0, 135.0d0 /)
        maskcorners(:,2) = (/ 0.0d0, 0.0d0, 0.0d0 /)
        maskcorners(:,3) = (/ 315.0d0, 0.0d0, 45.0d0 /)

        ! specify new value of 'loc' vector such that fwd/bwd diffs. are set up correctly in sparse matrix
        ! when function 'fillsprsebndy' is called. Also, specify appropriate values for the vectors 'normal'
        ! and 'fwdorbwd', which specify the orientation of the boundary normal and the direction of forward or
        ! backward differencing to be done in the lateral boundary condition functions 'normhorizmainbc_lat'
        ! and 'crosshorizmainbc_lat'

        ! following is algorithm for calculating boundary normal at 45 deg. increments, based on arbitray
        ! boundary shape

        where( thck3x3 .ne. 0.0d0 )
            thckmask = 0.0_dp
        elsewhere( thck3x3 .eq. 0.0d0 )
            thckmask = 1.0d0
        endwhere

        testvect = sum( thckmask * mask, 1 )

        !if( up .eq. 3 )then ! temporary code for debugging
        !  do i = 3,1,-1
        !  print *, 'thck = ', thck(:,i)
        !  end do
        !  print *, ' '
        !
        !  do i = 3,1,-1
        !      print *, 'thckmask = ', thckmask(:,i)
        !  end do
        !  print *, ' '
        !
        !  print *, 'testvect =  ', testvect
        !  print *, ' '
        !end if

        ! calculate the angle of the normal in cart. (x,y) system w/ 0 deg. at 12 O'clock, 90 deg. at 3 O'clock, etc.
        if( sum( sum( thckmask, 1 ) ) .eq. 1.0d0 )then
            phi = sum( sum( thckmask * maskcorners, 1 ) )
        else
            if( any( testvect .eq. 360.0d0 ) )then
                if( sum( testvect ) .eq. 450.0d0 )then
                    phi = 45.0d0
                elseif( sum( testvect ) .eq. 630.0d0 )then
                    phi = 315.0d0
                else
                    phi = 0.0d0
                end if
            elseif( all( testvect .ne. 360 ) )then
                phi = sum( testvect ) / sum( testvect/testvect, testvect .ne. 0.0d0 )
            end if
        end if
        
        calc_normal_45deg = deg2rad * phi

    end function

end module glide_mask
