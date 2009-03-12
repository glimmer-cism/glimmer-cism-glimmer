!Helper module containing routines to move between staggered and
!unstaggered grids
#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

module glide_grids
    implicit none
contains
  subroutine stagvarb(ipvr,opvr,ewn,nsn,choice_arg,usrf,thklim)

    use glimmer_global, only : dp ! ewn, nsn
    use glimmer_paramets, only: thk0
    implicit none 

    real(dp), intent(out), dimension(:,:) :: opvr 
    real(dp), intent(in), dimension(:,:) :: ipvr
    
    !Optional arguments that are used for fixing the staggering algorithm for
    !nunataks
    !If choice_arg is 1, then usrf and thklim had better be present!
    integer, intent(in), optional :: choice_arg
    real(dp), intent(in), dimension(:,:), optional :: usrf
    real(dp), optional :: thklim
    
    integer :: ewn,nsn,ew,ns
    integer :: choice 

    if (present(choice_arg)) then
        choice = choice_arg !If a choice was specified, use it
    else
        choice = 0 !Default to "vanilla" staggering algorithm
    end if

    if (choice == 1) then !when the input var is thk
        !ALBnew fix:
        !If given cell has less than 100m thickness, check to see whether 
        !it is upstream of a cell that does have thickness greater than 100m,
        !then we can set stagthck to be 0.
        !we have to have this extra check (to make sure the downstream cell has > 100m of ice) -
        !because at margin cells (where only one cell out of the four is ice free), there will be 
        !two ice free neighbours, and one must be upstream of the other, so it will set 
        !stagthck to be zero at the margins if we don't check the amount of ice in the downstream cell...
    
    do ns = 1,nsn-1
          do ew = 1,ewn-1

        !ew,ns cell is ice free:
        if (ipvr(ew,ns) <= thklim/thk0 .and. ((usrf(ew,ns) >= usrf(ew+1,ns) .and. ipvr(ew+1,ns) >= thklim/thk0) &
            .or. (usrf(ew,ns) >= usrf(ew,ns+1) .and. ipvr(ew,ns+1) >= thklim/thk0))) then
                opvr(ew,ns) = 0.0

        !ew+1,ns cell is ice free:
        else if (ipvr(ew+1,ns) <= thklim/thk0 .and. ((usrf(ew+1,ns) >= usrf(ew,ns) .and. ipvr(ew,ns) >= thklim/thk0) &
            .or. (usrf(ew+1,ns) >= usrf(ew+1,ns+1) .and. ipvr(ew+1,ns+1) >= thklim/thk0))) then
                opvr(ew,ns) = 0.0
    
        !ew,ns+1 cell is ice free:
        else if (ipvr(ew,ns+1) <= thklim/thk0 .and. ((usrf(ew,ns+1) >= usrf(ew,ns) .and. ipvr(ew,ns) >= thklim/thk0) &
            .or. (usrf(ew,ns+1) >= usrf(ew+1,ns+1) .and. ipvr(ew+1,ns+1) >= thklim/thk0))) then
                opvr(ew,ns) = 0.0
    
        !ew+1,ns+1 cell is ice free:
        else if (ipvr(ew+1,ns+1) <= thklim/thk0 .and. ((usrf(ew+1,ns+1) >= usrf(ew+1,ns) .and. ipvr(ew+1,ns) >=thklim/thk0) &
            .or. (usrf(ew+1,ns+1) >= usrf(ew,ns+1) .and. ipvr(ew,ns+1) >=thklim/thk0))) then
                opvr(ew,ns) = 0.0
        else
                opvr(ew,ns) = (ipvr(ew+1,ns) + ipvr(ew,ns+1) + &
                               ipvr(ew+1,ns+1) + ipvr(ew,ns)) / 4.0d0
        end if
  
        end do
    end do
    
    else
        opvr(1:ewn-1,1:nsn-1) = (ipvr(2:ewn,1:nsn-1) + ipvr(1:ewn-1,2:nsn) + &
                                 ipvr(2:ewn,2:nsn)   + ipvr(1:ewn-1,1:nsn-1)) / 4.0d0
    end if
  end subroutine stagvarb

end module glide_grids
