!Contains helper functions for nonlinear iteration, both to embed in the
!iteration loop and to serialize the data into the vector format that these
!functions require.
!Currently only unstable manifold correction is implemented.
module glide_nonlin
    use glimmer_global, only: dp
    use glimmer_physcon, only: pi
    implicit none
contains

    subroutine linearize_3d(vector, start, field)
        real(dp), dimension(:) :: vector
        integer :: start
        real(dp), dimension(:,:,:) :: field
        integer :: ni, nj, nk
        integer :: i,j,k

        ni = size(field, 1)
        nj = size(field, 2)
        nk = size(field, 3)
        
        do i=1,ni
            do j=1,nj
                do k=1,nk
                    vector(start) = field(i,j,k)
                    start = start + 1
                end do
            end do
        end do
    end subroutine

    subroutine linearize_2d(vector, start, field)
        real(dp), dimension(:) :: vector
        integer :: start
        real(dp), dimension(:,:) :: field
        integer :: ni, nj
        integer :: i,j

        ni = size(field, 1)
        nj = size(field, 2)
        
        do i=1,ni
            do j=1,nj
                vector(start) = field(i,j)
                start = start + 1
            end do
        end do
    end subroutine
    
    subroutine delinearize_3d(vector, start, field)
        real(dp), dimension(:) :: vector
        integer :: start
        real(dp), dimension(:,:,:) :: field
        integer :: ni, nj, nk
        integer :: i,j,k

        ni = size(field, 1)
        nj = size(field, 2)
        nk = size(field, 3)
        
        do i=1,ni
            do j=1,nj
                do k=1,nk
                    field(i,j,k) = vector(start)
                    start = start + 1
                end do
            end do
        end do
    end subroutine

    subroutine delinearize_2d(vector, start, field)
        real(dp), dimension(:) :: vector
        integer :: start
        real(dp), dimension(:,:) :: field
        integer :: ni, nj
        integer :: i,j

        ni = size(field, 1)
        nj = size(field, 2)
        
        do i=1,ni
            do j=1,nj
                field(i,j) = vector(start)
                start = start + 1
            end do
        end do
    end subroutine

    function unstable_manifold_correction(vec_new, vec_old, vec_correction, &
                                          vec_size, toler, tot_out, theta_out)
        logical :: unstable_manifold_correction
        
        real(dp), dimension(:), intent(in) :: vec_new
        real(dp), dimension(:), intent(inout) :: vec_old
        real(dp), dimension(:), intent(inout) :: vec_correction
        integer :: vec_size
        real(dp) :: toler
        real(dp), optional, intent(out) :: tot_out
        real(dp), optional, intent(out) :: theta_out

        real(dp) :: norm1, norm2, norm3, norm4, norm5
        real(dp) :: tot
        real(dp) :: theta
        real(dp) :: alpha
        integer :: i
        real(dp), dimension(vec_size) :: vec_correction_new

        !Assume we need to iterate again until proven otherwise
        unstable_manifold_correction = .true.

        norm1 = 0
        norm2 = 0
        norm3 = 0
        norm4 = 0
        norm5 = 0

        vec_correction_new = vec_new(1:vec_size) - vec_old(1:vec_size)

        do i = 1,vec_size
            norm1 = norm1 + (vec_correction_new(i) - vec_correction(i)) ** 2
            norm2 = norm2 + vec_correction(i) ** 2
            norm3 = norm3 + vec_correction_new(i) ** 2
            norm4 = norm4 + vec_correction(i) * vec_correction_new(i)
            norm5 = norm5 + vec_new(i) ** 2
        end do

        !Compute the angle between successive correction vectors
        if ((abs(norm2) < 1d-10) .or. (abs(norm3) < 1d-10)) then
            theta=PI/2.
        else
            theta=acos(norm4/sqrt(norm2*norm3))
        endif

        if ( (theta <= (5.*PI/6.) ) ) then
            !We've requested unstable manifold correction, and the angle is
            !small (less than 5pi/6, a value identified by Hindmarsh and Payne
            !to work well).   If this is the case, we compute and apply
            !a correction vector.
            
            !Compute the error between the last two *correction vectors* (not
            !the last two iteration values!)  (See (51) in Pattyn's paper)
            if (abs(norm2) > 0.) then !We're just avoiding a divide by 0 here
                alpha=sqrt(norm1/norm2)
            else
                alpha=1.
            endif
            
            if (alpha < 1.e-6) then
                !If the correction vector didn't change much, we're done
                unstable_manifold_correction = .false.
            else
                !Update the previous guess of the velocity with the correction
                !vector.  This throws out the current iteration's computed
                !velocity, and instead uses the computed correction vector.
                vec_old = vec_old + vec_correction_new / alpha
                vec_correction = vec_correction_new

            endif
        else
            !Copy this iteration's new values to the old values
            !for the next iteration - because the angle between correction
            !vectors is large we do not want to apply a correction, so
            !we just go Picard instead
            vec_old = vec_new
            vec_correction = vec_correction_new
        endif
        
        tot=sqrt(norm3/(norm5+1d-10)) !Regularize the denominator so we don't get NAN with simple geometries

        if (present(tot_out)) then
            tot_out = tot
        end if
        
        if (present(theta_out)) then
            theta_out = theta * 180 / pi
        end if

        if (tot < toler) unstable_manifold_correction = .false. 
    end function unstable_manifold_correction
end module glide_nonlin
