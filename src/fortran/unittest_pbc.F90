!Unit test module for periodic boundary conditions
module unittest_pbc
    use glide_grids
    use glimmer_global, only: dp
    use unittest_framework
contains
    subroutine test_pbc()

       real(dp), dimension(8, 1) :: arr
       real(dp), dimension(1, 8) :: arr_t

       call unittest_begin_section("Periodic Boundary Conditions")

       arr(:,1) = (/777, 1, 2, 3, 4, 3, 2, 777/)
       call periodic_boundaries(arr, .true., .false.)
       call unittest_result("1 ghost cell, x direction", all(arr(1:2,:) == arr(7:8,:)))

       arr(:,1) = (/777, 777, 1, 2, 3, 2, 777, 777/)
       call periodic_boundaries(arr, .true., .false., 2)
       call unittest_result("2 ghost cells, x direction", all(arr(1:4,:) == arr(5:8,:)))

       arr_t(1,:) = (/777, 777, 1, 2, 3, 2, 777, 777/)
       call periodic_boundaries(arr_t, .false., .true., 2)
       call unittest_result("2 ghost cells, y direction", all(arr_t(:,1:4) == arr_t(:,5:8)))


       arr_t(1,:) = (/777, 777, 1, 2, 3, 2, 777, 777/)
       call periodic_boundaries(arr_t, .false., .true., 2)
       call unittest_result("2 ghost cells, y direction", all(arr_t(:,1:4) == arr_t(:,5:8)))


       call unittest_end_section()

    end subroutine
end module unittest_pbc
