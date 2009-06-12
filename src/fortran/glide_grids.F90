!Helper module containing routines to move between staggered and
!unstaggered grids
#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

module glide_grids
    implicit none

    integer, parameter :: STAGGER_CHOICE_STD = 0      !Simple linear interpolation onto staggered grid
    integer, parameter :: STAGGER_CHOICE_ALBNEW = 1   !Use the ALBnew mass conservation fix
    integer, parameter :: STAGGER_CHOICE_HARMONIC = 2 !Use the harmonic mean to bias points on shelf edges
    integer, parameter :: STAGGER_IGNORE_ZEROS = 3

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
    
    integer :: ewn,nsn,ew,ns,n
    integer :: choice
    real(dp) :: tot

    if (present(choice_arg)) then
        choice = choice_arg !If a choice was specified, use it
    else
        choice = STAGGER_CHOICE_STD !Default to "vanilla" staggering algorithm
    end if

    if (choice == STAGGER_CHOICE_ALBNEW) then !when the input var is thk
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
                if (ipvr(ew,ns) <= thklim/thk0 .and. &
                   ((usrf(ew,ns) >= usrf(ew+1,ns) .and. ipvr(ew+1,ns) >= thklim/thk0) &
                    .or. (usrf(ew,ns) >= usrf(ew,ns+1) .and. ipvr(ew,ns+1) >= thklim/thk0))) then
                        opvr(ew,ns) = 0.0

                !ew+1,ns cell is ice free:
                else if (ipvr(ew+1,ns) <= thklim/thk0 .and. &
                    ((usrf(ew+1,ns) >= usrf(ew,ns) .and. ipvr(ew,ns) >= thklim/thk0) &
                    .or. (usrf(ew+1,ns) >= usrf(ew+1,ns+1) .and. ipvr(ew+1,ns+1) >= thklim/thk0))) then
                        opvr(ew,ns) = 0.0
    
                !ew,ns+1 cell is ice free:
                else if (ipvr(ew,ns+1) <= thklim/thk0 .and. &
                    ((usrf(ew,ns+1) >= usrf(ew,ns) .and. ipvr(ew,ns) >= thklim/thk0) &
                    .or. (usrf(ew,ns+1) >= usrf(ew+1,ns+1) .and. ipvr(ew+1,ns+1) >= thklim/thk0))) then
                        opvr(ew,ns) = 0.0
    
                !ew+1,ns+1 cell is ice free:
                else if (ipvr(ew+1,ns+1) <= thklim/thk0 .and. &
                    ((usrf(ew+1,ns+1) >= usrf(ew+1,ns) .and. ipvr(ew+1,ns) >=thklim/thk0) &
                    .or. (usrf(ew+1,ns+1) >= usrf(ew,ns+1) .and. ipvr(ew,ns+1) >=thklim/thk0))) then
                        opvr(ew,ns) = 0.0
                
                !Standard Staggering
                else
                        opvr(ew,ns) = (ipvr(ew+1,ns) + ipvr(ew,ns+1) + &
                               ipvr(ew+1,ns+1) + ipvr(ew,ns)) / 4.0d0
                end if
  
        end do
    end do



    else if (choice == STAGGER_CHOICE_HARMONIC) then
        !Harmonic mean scheme from Price 
        do ns = 1, nsn-1
            do ew = 1, ewn-1
                !Only apply the harmonic mean if we are next to a point of zero
                !thickness.
                if (any(ipvr(ew:ew+1, ns:ns+1) > 0.0D0) .and. &
                    any(ipvr(ew:ew+1, ns:ns+1) > 0.0D0)) then
                        !Compute the harmonic mean
                        !1. harm. mean is biased on the LOW side, so to use it to bias thickness 
                        !HIGH, take the max thickness of the four values available and subtract 
                        !all other values from this max value
                        !
                        !2. take the harm. mean of these differences and subtract that from the max 
                        !value (the factor of 0.02 out front is used to give some non-zero 
                        !thickness at points on the normal grid where thickness is actually zero. 
                        !It's needed to avoid div. by zero when calc. the harm mean. I choose the 
                        !value by playing around w/ some representative thickness stencils in matlab)
                        !
                        !3. apply this diff. (max thick on normal stencil minus harm mean adjustment)
                        !and apply that as the value for the stag thickness
                        opvr(ew,ns) = maxval( ipvr(ew:ew+1,ns:ns+1) ) - 4.0d0 / ( &
                            sum( 1.0d0 / (-(ipvr(ew:ew+1,ns:ns+1) - &
                            maxval(ipvr(ew:ew+1,ns:ns+1))) + &
                            0.02d0*maxval(ipvr(ew:ew+1,ns:ns+1))  ) ) ) 
                else
                    opvr(ew, ns) = sum(ipvr(ew:ew+1, ns:ns+1))/4.0D0
                end if
            end do
        end do

    !Apply the standard staggering algorithm, but do not include points with
    !zero thickness in the average
    else if (choice == STAGGER_IGNORE_ZEROS) then
        do ns = 1, nsn-1
            do ew = 1, ewn-1
                n = 0
                tot = 0

                if (abs(ipvr(ew,ns)) > 1e-10) then
                    tot = tot + ipvr(ew,ns)
                    n   = n   + 1
                end if
                if (abs(ipvr(ew+1,ns)) > 1e-10) then
                    tot = tot + ipvr(ew+1,ns)
                    n   = n   + 1
                end if
                if (abs(ipvr(ew,ns+1)) > 1e-10) then
                    tot = tot + ipvr(ew,ns+1)
                    n   = n   + 1
                end if
                if (abs(ipvr(ew+1,ns+1)) > 1e-10) then
                    tot = tot + ipvr(ew+1,ns+1)
                    n   = n   + 1
                end if
                opvr(ew,ns) = tot/n
            end do
        end do

    else !Use standard staggering
        opvr(1:ewn-1,1:nsn-1) = (ipvr(2:ewn,1:nsn-1) + ipvr(1:ewn-1,2:nsn) + &
                                 ipvr(2:ewn,2:nsn)   + ipvr(1:ewn-1,1:nsn-1)) / 4.0d0
    end if
  end subroutine stagvarb

end module glide_grids
