!Helper module containing routines to move between staggered and
!unstaggered grids
#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

#include <glide_nan.inc>

module glide_grids
    use glimmer_global, only : dp, NaN
    implicit none

    integer, parameter :: STAGGER_CHOICE_STD = 0      !Simple linear interpolation onto staggered grid
    integer, parameter :: STAGGER_CHOICE_ALBNEW = 1   !Use the ALBnew mass conservation fix
    integer, parameter :: STAGGER_CHOICE_HARMONIC = 2 !Use the harmonic mean to bias points on shelf edges
    integer, parameter :: STAGGER_IGNORE_ZEROS = 3

contains
  subroutine stagvarb(ipvr,opvr,ewn,nsn,choice_arg,usrf,thklim)
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
                if (n > 0) then
                    opvr(ew,ns) = tot/n
                else
                    opvr(ew,ns) = 0
                end if
            end do
        end do

    else !Use standard staggering
        opvr(1:ewn-1,1:nsn-1) = (ipvr(2:ewn,1:nsn-1) + ipvr(1:ewn-1,2:nsn) + &
                                 ipvr(2:ewn,2:nsn)   + ipvr(1:ewn-1,1:nsn-1)) / 4.0d0
    end if
  end subroutine stagvarb

  subroutine stagvarb_3d(ipvr, opvr, ewn, nsn, upn)
    real(dp), intent(in), dimension(:,:,:) :: ipvr
    real(dp), intent(out), dimension(:,:,:) :: opvr
    integer, intent(in) :: ewn, nsn, upn
    integer :: k

    do k = 1, upn
        call stagvarb(ipvr(k,:,:), opvr(k,:,:), ewn, nsn)
    end do
  end subroutine stagvarb_3d


    !*FD Copies a staggered grid onto a nonstaggered grid.  This verion
    !*FD assumes periodic boundary conditions.
    subroutine unstagger_field_2d(f_stag, f, periodic_x, periodic_y)
        real(dp), dimension(:,:), intent(in) :: f_stag
        real(dp), dimension(:,:), intent(out) :: f
        logical, intent(in) :: periodic_x, periodic_y

        real(dp), dimension(4) :: pts

        real(dp) :: s,n

        integer :: i,j, k,i1, i2, j1, j2, ni, nj
        
        ni = size(f, 1)
        nj = size(f, 2)

        do i = 1, size(f, 1)
            do j = 1, size(f, 2)
                s = 0
                n = 0
                
                i1 = i-1
                i2 = i
                
                !If we're unstaggering with periodic boundaries, we cross over to the
                !other side of the domain when we "de-average".  Otherwise, we just ignore
                !the point that's off the domain.
                if (i1 == 0) then
                    if (periodic_y) then
                        i1 = ni - 1
                    else
                        i1 = 1
                    end if
                end if
    
                if (i2 == ni) then
                    if (periodic_y) then
                        i2 = 1
                    else
                        i2 = ni - 1
                    end if
                end if
    
                j1 = j-1
                j2 = j
    
                if (j1 == 0) then
                    if (periodic_y) then
                        j1 = nj - 1
                    else
                        j1 = 1
                    end if
                end if
    
                if (j2 == nj) then
                    if (periodic_x) then
                        j2 = 1
                    else
                        j2 = nj - 1
                    end if
                end if
                
                !Place the points into an array, loop over them, and average
                !all the points that AREN'T NaN.
                pts = (/f_stag(i1, j1), f_stag(i2, j1), f_stag(i1, j2), f_stag(i2, j2)/)
            
                do k=1,4
                    if (.not. (IS_NAN(pts(k)))) then
                        s = s + pts(k)
                        n = n + 1
                    end if
                end do
                if (n /= 0) then
                    f(i,j) = s/n
                else
                    f(i,j) = NaN
                end if
            end do
        end do
    
    end subroutine unstagger_field_2d

    subroutine unstagger_field_3d(f, f_stag, periodic_x, periodic_y)
        real(dp), dimension(:,:,:) :: f, f_stag
        logical, intent(in) :: periodic_x, periodic_y

        integer :: i

        do i = 1,size(f,1)
            call unstagger_field_2d(f(i,:,:), f_stag(i,:,:), periodic_x, periodic_y)
        end do
        
    end subroutine unstagger_field_3d


    subroutine periodic_boundaries(m, apply_to_x, apply_to_y, nlayers_arg)
        !*FD Applies periodic boundary conditions to a 2D array
        real(dp), dimension(:,:), intent(inout) :: m
        integer :: maxx, maxy
        logical :: apply_to_x, apply_to_y
        integer, optional :: nlayers_arg

        integer :: nlayers 

        if (present(nlayers_arg)) then
            nlayers = nlayers_arg
        else
            nlayers = 1
        end if

        maxx = size(m, 1)
        maxy = size(m, 2)
       
        if (apply_to_x) then
            m( 1 : nlayers, : ) = m( maxx-nlayers*2 + 1 : maxx - nlayers, :)
            m( maxx-nlayers+1 : maxx, : ) = m(nlayers + 1 : nlayers*2, : )
        end if

        if (apply_to_y) then
            m( :, 1 : nlayers ) = m( :, maxy-nlayers*2 + 1 : maxy - nlayers )
            m( :, maxy-nlayers+1 : maxy ) = m( :, nlayers + 1 : nlayers*2 )
        end if

        
        !If both directions are periodic, treat the corners specially.
        if(apply_to_x .and. apply_to_y) then
            m(1:nlayers, 1:nlayers) = m(maxx-nlayers*2+1:maxx-nlayers, maxy-nlayers*2+1:maxy-nlayers)
            m(nlayers+1:2*nlayers, nlayers+1:2*nlayers) = m(maxx-nlayers+1:maxx, maxy-nlayers+1:maxy)
            
            m(1:nlayers, maxy-nlayers+1:maxy) = m(maxx-nlayers*2+1:maxx-nlayers, nlayers+1:2*nlayers)
            m(nlayers+1:2*nlayers, maxy-nlayers*2+1:maxy-nlayers) = m(maxx-nlayers+1:maxx, 1:nlayers)
        end if
    end subroutine periodic_boundaries
    
    subroutine periodic_boundaries_3d(m, apply_to_x, apply_to_y, nlayers_arg)
        !*FD Applies periodic boundary conditions to a 3D array
        real(dp), dimension(:,:,:), intent(inout) :: m
        integer :: maxx, maxy
        logical :: apply_to_x, apply_to_y
        integer, optional :: nlayers_arg
    
        integer :: i
        
        do i = 1, size(m,1)
            call periodic_boundaries(m(i,:,:), apply_to_x, apply_to_y, nlayers_arg)
        end do
    end subroutine periodic_boundaries_3d
 
end module glide_grids
