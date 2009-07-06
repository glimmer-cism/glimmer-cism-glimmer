!Contains a framework for registering test results
module unittest_framework

    character(128) :: g_section_name

    integer :: g_global_run, g_global_passed, g_section_run, g_section_passed

contains
    subroutine unittest_init()
    
        g_global_run = 0
        g_global_passed = 0

    end subroutine unittest_init

    subroutine unittest_finalise()

        write(*,*) "*********************************************"
        if (g_global_passed == g_global_run) then
            write(*,*) "All tests passed"
        else
            write(*,*) g_global_passed, "of", g_global_run, "tests passed."
        end if
    end subroutine unittest_finalise

    subroutine unittest_begin_section(name)
        character(*) :: name
        g_section_name = name

        g_section_run = 0
        g_section_passed = 0
        
    end subroutine

    subroutine unittest_end_section()

        if (g_global_passed == g_global_run) then
            write(*,*) g_section_name, ": All tests passed"
        else
            write(*,*) g_section_name, ":", g_global_passed, "of", g_global_run, "tests passed."
        end if
    end subroutine

    subroutine unittest_result(test_descr, flag)
        character(*) :: test_descr
        logical :: flag

        g_global_run = g_global_run + 1
        g_section_run = g_section_run + 1

        if (.not. flag) then
            write(*,*) g_section_name,": test FAILED:", test_descr
        else
            g_global_passed = g_global_passed + 1
            g_section_passed = g_section_passed + 1
        end if
    end subroutine

end module
