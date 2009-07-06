module unittest_ice3d_lib
    use ice3d_lib
    implicit none
contains
    subroutine test_ice3d_lib_stencil()
        integer :: i

        do i=1,25
            write(*,*) i, stencil_y(i,0), stencil_x(i,0),stencil_z(i,0)
        end do
    end subroutine
end module unittest_ice3d_lib

