module glimmer_sparse_solver
    !*FD This module builds on the glimmer_sparse module to provide an easy
    !*FD interface to SLAP.  The SLAP interface is intended to be both
    !*FD usable and a guide to implementing other interfaces
    
    use glimmer_sparse
    use glimmer_global, only: dp
    implicit none

#ifdef PARDISO_64BIT
        integer, parameter :: parint = 8       
#else
        integer, parameter :: parint = 4
#endif


    type sparse_solver_workspace
        !*FD This type contains any working memory needed for the sparse solver.
        !*FD It is used to store states between calls to the solver
        !*FD In the SLAP implementation, it is used to store the SLAP workspace
        !*FD This module must have this type, but its contents should be opaque
        !*FD to the user (e.g. client code should only manipulate the
        !*FD sparse_solver_workspace as a whole and should never touch its members)
        integer(kind=parint), dimension(64) :: pt
        integer(kind=parint) :: mytpe
        integer(kind=parint), dimension(64) :: iparam
        logical :: alloc
    end type sparse_solver_workspace

    type sparse_solver_options
        !*FD This type holds options that are passed to the sparse solver, such
        !*FD as preconditioner type, error tolerances, etc.  At a minimum, it
        !*FD must define the tolerance and maxiters field, as these will be
        !*FD common to any iterative sparse linear solver.  Other options
        !*FD can be defined as necessary.
        !*FD
        !*FD Design note: the options are seperated from the workspace because
        !*FD one set of options could apply to multiple matrices, and the
        !*FD lifecycles for each could be different (a workspace need only
        !*FD exist as long as the matrix does, the options could persist
        !*FD throughout the entire program)
        real(kind=dp) :: tolerance !*FD Error tolerance
        integer :: maxiters !*FD Max iterations before giving up
    end type sparse_solver_options

contains
    subroutine sparse_solver_default_options(opt)
        !*FD Populates a sparse_solver_options (defined above) with default
        !*FD options.  This is necessary because different solvers may define
        !*FD different options beyond the required fields defined above.
        !*FD Filling them in this function allows client code to pick "good"
        !*FD values in a generic way.
        type(sparse_solver_options), intent(out) :: opt
        opt%tolerance  = 1e-5
        opt%maxiters = 2000
    end subroutine sparse_solver_default_options

    subroutine sparse_allocate_workspace(matrix, options, workspace, max_nonzeros_arg)
        !*FD Allocate solver workspace.  This needs to be done once
        !*FD (when the maximum number of nonzero entries is first known)
        
        !*FD Note that the max_nonzeros argument must be optional, and if
        !*FD the current number of nonzeroes must be used.
        type(sparse_matrix_type) :: matrix
        type(sparse_solver_options) :: options
        type(sparse_solver_workspace) :: workspace
        integer, optional :: max_nonzeros_arg
   
        call pardisoinit(workspace%pt, workspace%mytpe, workspace%iparam)
    
    end subroutine sparse_allocate_workspace

    subroutine sparse_solver_preprocess(matrix, options, workspace)
        !*FD Performs any preprocessing needed to be performed on the sparse
        !*FD matrix.  Workspace must have already been allocated. 
        !*FD This function should be safe to call more than once.
        !*FD 
        !*FD It is an error to call this function on a workspace without already
        !*FD allocated memory.
        !*FD
        !*FD In general sparse_allocate_workspace should perform any actions
        !*FD that depend on the *size* of the sparse matrix, and
        !*FD sprase_solver_preprocess should perform any actions that depend
        !*FD upon the *contents* of the sparse matrix.
        type(sparse_matrix_type) :: matrix
        type(sparse_solver_options) :: options
        type(sparse_solver_workspace) :: workspace
        !Convert to the UMFPACK format.
        !Move from triad format to column format
        call to_column_format(matrix)

        !Sort by row within each column in the column format.
        !to_column_format does not necessarily do this
        call sort_column_format(matrix)

        !Move from 1-based to 0-based indexing
        matrix%row = matrix%row - 1
        matrix%col = matrix%col - 1

    end subroutine sparse_solver_preprocess

    function sparse_solve(matrix, rhs, solution, options, workspace,err,niters, verbose)
        !*FD Solves the sparse linear system, and reports status information.
        !*FD This function returns an error code that should be zero if the
        !*FD call succeeded and nonzero if it failed.  No additional error codes
        !*FD are defined.  Although this function reports back the final error
        !*FD and the number of iterations needed to converge, these should *not*
        !*FD be relied upon as not every sparse linear solver may report them.
        type(sparse_matrix_type), intent(inout) :: matrix 
        !*FD Sparse matrix to solve.  This is inout because the sparse solver
        !*FD may have to do some re-arranging of the matrix.
        
        real(kind=dp), dimension(:), intent(in) :: rhs 
        !*FD Right hand side of the solution vector
        
        real(kind=dp), dimension(:), intent(inout) :: solution 
        !*FD Solution vector, containing an initial guess.

        type(sparse_solver_options), intent(in) :: options
        !*FD Options such as convergence criteria
        
        type(sparse_solver_workspace), intent(inout) :: workspace
        !*FD Internal solver workspace
        
        real(kind=dp), intent(out) :: err
        !*FD Final solution error
        
        integer, intent(out) :: niters
        !*FD Number of iterations required to reach the solution

        logical, optional, intent(in) :: verbose
        !*FD If present and true, this argument may cause diagnostic information
        !*FD to be printed by the solver (not every solver may implement this).
        
        integer :: sparse_solve

        integer :: iunit !Unit number to print verbose output to (6=stdout, 0=no output)
        integer :: mtype

        iunit=0
        if (present(verbose)) then
            if(verbose) then
                iunit=1
                write(*,*),"Tolerance=",options%tolerance
            end if
        end if


        !Detect symmetric matrix
        if (matrix%symmetric) then
            mtype = 1
        else
            mtype = 11
        end if
        
        !PARDISO arguments
        call pardiso(workspace%pt, 1, 1, mtype, 33, matrix%order, matrix%val, &
                     matrix%row, matrix%col, 0, matrix%order, workspace%iparam, &
                     iunit, rhs, solution, sparse_solve)
        
        err = 0
        niters = 0
    end function sparse_solve

    subroutine sparse_destroy_workspace(matrix, options, workspace)
        !*FD Deallocates all working memory for the sparse linear solver.
        !*FD This need *not* be safe to call of an unallocated workspace
        !*FD No sparse solver should call this automatically.
        type(sparse_matrix_type) :: matrix
        type(sparse_solver_options) :: options
        type(sparse_solver_workspace) :: workspace

        !TODO: Call pardiso with flag to delete

    end subroutine sparse_destroy_workspace

    subroutine sparse_interpret_error(error_code, error_string)
        !*FD takes an error code output from sparse_solve and interprets it.
        !*FD error_string must be an optional argument.
        !*FD If it is not provided, the error is printed to standard out
        !*FD instead of being put in the string
        integer :: error_code
        character(*), optional, intent(out) :: error_string
        character(256) :: tmp_error_string
        
        select case (error_code)
            case (0)
                tmp_error_string = "No error"
            case (-1)
                tmp_error_string = "Input inconsistent"
            case (-2)
                tmp_error_string = "Not enough memory"
            case (-3)
                tmp_error_string = "Reordering problem"
            case (-4)
                tmp_error_string = "Zero pivot, numerical fact. or iterative refinement problem"
            case (-5)
                tmp_error_string = "Unclassified (internal) errror"
            case (-6)
                tmp_error_string = "Preordering failed (matrix types 11, 13 only)"
            case (-7)
                tmp_error_string = "Diagonal matrix problem"
            case (-8)
                tmp_error_string = "32-bit integer overflow problem"
        end select

        if (present(error_string)) then
            error_string = tmp_error_string
        else
            write(*,*) tmp_error_string
        endif
    end subroutine sparse_interpret_error
end module glimmer_sparse_solver
