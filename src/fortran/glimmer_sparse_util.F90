!This module contains convenient cover routines for sparse matrix
!solving and error handling that are agnostic of the sparse matrix
!solver being used.  It is separate from the glimmer_sparse module
!because placing these routines in glimmer_sparse would create a 
!circular dependancy between the storage module and the sparse solver
module glimmer_sparse_util
    use glimmer_global, only: sp, dp
    use glimmer_sparse
    use glimmer_sparse_solver
    implicit none
contains
    subroutine sparse_easy_solve(matrix, rhs, answer, err, iter, error_flag, calling_file, calling_line)
        !This subroutine wraps the basic (though probably the most inefficient)
        !workflow to solve a sparse matrix using the sparse matrix solver
        !framework.  It handles errors gracefully, and reports back the
        !iterations required and the error estimate in the case of an iterative
        !solver.  At the very least it is an encapsulated example of how to
        !use the sparse solver routines, and is easy enough to drop in your
        !code if you don't care about allocating and deallocating workspace
        !every single timestep.

        type(sparse_matrix_type) :: matrix
        real(dp), dimension(:) :: rhs
        real(dp), dimension(:) :: answer
        
        real(dp), intent(out) :: err
        integer, intent(out) :: iter

        logical, intent(out) :: error_flag

        character(100), optional :: calling_file
        integer, optional :: calling_line

        type(sparse_solver_options) :: opt
        type(sparse_solver_workspace) :: wk

        integer :: ierr

        call sparse_solver_default_options(opt)
        call sparse_allocate_workspace(matrix, opt, wk)
        call sparse_solver_preprocess(matrix, opt, wk)

        ierr = sparse_solve(matrix, rhs, answer, opt, wk, err, iter, .false.)
        
        if (ierr /= 0) then
            if (present(calling_file) .and. present(calling_line)) then
                call log_sparse_error(matrix, ierr, calling_file, calling_line)
            else
                call log_sparse_error(matrix, ierr, __FILE__, __LINE__)
            end if
            error_flag = .true.
        end if

        call sparse_destroy_workspace(matrix, opt, wk)

    end subroutine sparse_easy_solve

    subroutine log_sparse_error(matrix, error, error_file, error_line, time)
        !Checks a sparse error flag and, if an error occurred, log it to
        !the GLIMMER log file.  This does not stop Glimmer, it just writes
        !to the log
        use glide_stop
        use glimmer_log
        use glimmer_filenames
        
        integer :: error, error_line
        character(*) :: error_file
        real(dp), optional :: time

        type(sparse_matrix_type) :: matrix

        integer :: isym
        integer :: lunit

        character(512) :: message
        character(128) :: errfname
        character(256) :: errdesc

        !If no error happened, this routine should be a nop
        if (error == 0) return

        !Aquire a file unit, and open the file
        lunit = get_free_unit()
        errfname = trim(process_path('sparse_dump.txt'))
        open(lunit,file=errfname)

        if (matrix%symmetric) then
            isym = 1
        else
            isym = 0
        end if

        !Output sparse matrix data to the file
        call dcpplt(matrix%order, matrix%n, matrix%row, matrix%col, matrix%val,&
                    isym, lunit)

        write(lunit,*) '***Sparse matrix structure ends.  Value listing begins'
        write(lunit,*) matrix%val

        !Close unit and finish off
        close(lunit)
        !call glide_finalise(model, .true.)

        !Grab the error message from the sparse solver
        call sparse_interpret_error(error, errdesc)

        !construct the error message and write it to the log file
        if (present(time)) then
            write(message, *)'Sparse matrix error at time: ', time, &
                             'Error description: ', errdesc, &
                             'Data dumped to ', trim(errfname)
        else
            write(message, *)'Sparse matrix error. Error description: ', errdesc, &
                             'Data dumped to ', trim(errfname)
        end if

        call write_log(trim(message), GM_ERROR, error_file, error_line)
    end subroutine log_sparse_error

end module glimmer_sparse_util
