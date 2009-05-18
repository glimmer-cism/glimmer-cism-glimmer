! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_sparse.f90 - part of the GLIMMER ice model       + 
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

module glimmer_sparse_type
  use glimmer_global, only:dp
  type sparse_matrix_type
     !*FD sparse matrix type
     integer :: nonzeros                                    !*FD number of nonzero elements currently stored
     integer :: order                                       !*FD order of the matrix (e.g. number of rows)
     logical :: symmetric                                   !*FD True if only one triangle of the symmetric matrix is stored
     integer, dimension(:), pointer :: col => NULL()        !*FD column index
     integer, dimension(:), pointer :: row => NULL()        !*FD row index
     real(kind=dp), dimension(:), pointer :: val => NULL()  !*FD values

  end type sparse_matrix_type

  type sparse_solver_options_base
        real(kind=dp) :: tolerance !*FD Error tolerance
        integer :: maxiters !*FD Max iterations before giving up
        integer :: method
  end type

  ! size of sparse matrix 
  integer, parameter, private :: chunksize=1000

  !MAKE_RESTART
#ifdef RESTARTS
#define RST_GLIMMER_SPARSE
#include "glimmer_rst_head.inc"
#undef RST_GLIMMER_SPARSE
#endif

contains

#ifdef RESTARTS
#define RST_GLIMMER_SPARSE
#include "glimmer_rst_body.inc"
#undef RST_GLIMMER_SPARSE
#endif

  subroutine new_sparse_matrix(order,n,mat)
    !*FD create a new sparse matrix
    implicit none
    integer, intent(in) :: n          !*FD initial size of matrix
    type(sparse_matrix_type) :: mat   !*FD matrix
    integer, intent(in) :: order      !*FD Order (number of rows and columns) of the matrix
    
    if (.not.associated(mat%col)) then
       allocate(mat%row(n))
       !SLAP's sparse column scheme looks past the assumed bounds of col to see
       !what sparse storage format we're in.  To avoid array bounds problems, we
       !add 2 to the column size.  See mailing list discussion at:
       !http://forge.nesc.ac.uk/pipermail/glimmer-discuss/2005-February/000078.html
       allocate(mat%col(n+2))
       allocate(mat%val(n))
    else
       if (size(mat%row).lt.n) then
          call del_sparse_matrix(mat)
          allocate(mat%row(n))
          allocate(mat%col(n+2))
          allocate(mat%val(n))
       end if
    end if
    mat%nonzeros = 0
    mat%order = order
    mat%symmetric = .false.
  end subroutine new_sparse_matrix

  subroutine copy_sparse_matrix(inmat,outmat)
    !*FD copy a sparse matrix.
    !*FD Slap workspace allocation on the new
    !*FD matrix is *not* done.
    implicit none
    type(sparse_matrix_type) :: inmat  !*FD matrix to be copied
    type(sparse_matrix_type) :: outmat !*FD result matrix

    call new_sparse_matrix(inmat%order,inmat%nonzeros,outmat)
    outmat%row(:) = inmat%row(:)
    outmat%col(:) = inmat%col(:)
    outmat%val(:) = inmat%val(:)
    outmat%nonzeros = inmat%nonzeros
    outmat%symmetric = inmat%symmetric
  end subroutine copy_sparse_matrix

  subroutine grow_sparse_matrix(matrix)
    !*FD grow sparse matrix
    implicit none
    type(sparse_matrix_type) :: matrix !*FD matrix

    integer, dimension(:), pointer :: newrow,newcol
    real(kind=dp), dimension(:), pointer :: newval
    integer oldsize

    oldsize = size(matrix%val)
    
    allocate(newrow(chunksize+oldsize))
    allocate(newcol(chunksize+oldsize))
    allocate(newval(chunksize+oldsize))
    write(*,*)size(matrix%col), size(matrix%row), size(matrix%val), size(newcol), size(newrow), size(newval)
    newcol(1:oldsize) = matrix%col(:)
    newrow(1:oldsize) = matrix%row(:)
    newval(1:oldsize) = matrix%val(:)

    deallocate(matrix%col)
    deallocate(matrix%row)
    deallocate(matrix%val)

    matrix%col => newcol
    matrix%row => newrow
    matrix%val => newval

  end subroutine grow_sparse_matrix

  subroutine del_sparse_matrix(matrix)
    !*FD delete sparse matrix
    implicit none
    type(sparse_matrix_type) :: matrix !*FD matrix

    deallocate(matrix%col)
    deallocate(matrix%row)
    deallocate(matrix%val)

  end subroutine del_sparse_matrix

  subroutine print_sparse(matrix, unit)
    !*FD print sparse matrix
    implicit none
    type(sparse_matrix_type) :: matrix !*FD matrix
    integer, intent(in) :: unit        !*FD unit to be printed to

    integer i
    do i = 1, matrix%nonzeros
       write(unit,*) matrix%col(i), matrix%row(i), matrix%val(i)
    end do
  end subroutine print_sparse

  subroutine sparse_matrix_vec_prod(matrix, vec, res)
    !*FD sparse matrix vector product
    implicit none
    type(sparse_matrix_type) :: matrix                !*FD matrix
    real(kind=dp), intent(in), dimension(:) :: vec    !*FD input vector
    real(kind=dp), intent(out), dimension(:) :: res   !*FD result vector

    integer i

    res = 0.
    do i=1,matrix%nonzeros
       res(matrix%col(i)) = res(matrix%col(i)) + vec(matrix%row(i))*matrix%val(i)
    end do
  end subroutine sparse_matrix_vec_prod

  subroutine sparse_insert_val(matrix, i, j, val)
    !*FD insert value into sparse matrix.  This is safe to call even if val=0
    implicit none
    type(sparse_matrix_type) :: matrix !*FD matrix
    integer, intent(in) :: i,j         !*FD column and row
    real(kind=dp), intent(in) :: val   !*FD value
    if (val /= 0.0 .and. i > 0 .and. j > 0 .and. i <= matrix%order .and. j <= matrix%order) then
        matrix%nonzeros =  matrix%nonzeros + 1
        matrix%row(matrix%nonzeros) = i
        matrix%col(matrix%nonzeros) = j
        matrix%val(matrix%nonzeros) = val

        if (matrix%nonzeros .eq. size(matrix%val)) then
            call grow_sparse_matrix(matrix)
        end if
    end if
  end subroutine sparse_insert_val

  subroutine sparse_clear(matrix)
    !*FD Clears the sparse matrix, without deallocating any of the
    !*FD previously used memory
    type(sparse_matrix_type) :: matrix
    
    matrix%nonzeros = 0
    !Clearing these shouldn't be strictly necessary, but SLAP barfs if we don't
    matrix%row = 0
    matrix%col = 0
    matrix%val = 0
  end subroutine

  function is_triad_format(matrix)
    type(sparse_matrix_type) :: matrix
    logical :: is_triad_format

    is_triad_format = .not. is_column_format(matrix)

  end function

  function is_column_format(matrix)
    type(sparse_matrix_type) :: matrix
    logical :: is_column_format

    is_column_format = matrix%col(matrix%order + 1) == matrix%nonzeros + 1
  end function

  subroutine to_column_format(matrix)
    type(sparse_matrix_type) :: matrix
     
    if(is_triad_format(matrix)) then
        call ds2y(matrix%order, matrix%nonzeros, matrix%row, matrix%col, matrix%val, 0)
    end if
  end subroutine

  subroutine sort_column_format(matrix)
    !*FD Takes a column format matrix and sorts the row indices within each column
    !*FD This is not strictly needed in some compressed-column matrices
    !*FD (e.g. those used in SLAP), but it *is* necessary in some other libraries
    !*FD (e.g. UMFPACK).  For this reason, it is not done automatically in
    !*FD to_column_format.
    type(sparse_matrix_type) :: matrix
    integer :: i
    do i=1,matrix%order !Loop through each column index
      call sort_column(matrix%val, matrix%row, matrix%col(i), matrix%col(i+1)-1)
    end do
  end subroutine

  subroutine sort_column(values, row_indices, startindex, endindex)
    real(dp),dimension(:) :: values
    integer,dimension(:) :: row_indices
    integer :: startindex
    integer :: endindex
    
    !Insertion Sort
    !TODO: something faster?
    do i=startindex+1,endindex
        currentrowindex = row_indices(i)
        currentvalue = values(i)

        j = i-1
        do while (j >= startindex .and. row_indices(j) > currentrowindex)
            row_indices(j+1) = row_indices(j)
            values(j+1) = values(j)
            j = j - 1
        end do
        row_indices(j+1) = currentrowindex
        values(j+1) = currentvalue
    end do
    
  end subroutine

end module glimmer_sparse_type
