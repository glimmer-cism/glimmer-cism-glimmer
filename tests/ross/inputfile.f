      subroutine INPUTFILE()

c     ---------------------------------------------------------
c     Entree de la geometrie et des conditions aux limites
c     ---------------------------------------------------------

          IMPLICIT NONE
          include 'shelfy.h'
	  integer i0,j0
	  DOUBLE PRECISION azi,magv

	  DX = 6822.
	  DY = 6822.

          DO I=1,NX
             DO J=1,NY
             CLIMMX(I,J) = .FALSE.
             CLIMMY(I,J) = .FALSE.
	     H(i,j)=H(i,j)*mask(i,j)
             if (icefront(i,j).eq.1) H(i,j)=1.
             END DO
          END DO
          
	  open(11,file='kbc.dat',err=12)
	  do i=1,nx
	     do j=1,ny
	     read(11,*,err=12) j0,i0
	     climmx(i0,j0) = .true.
	     climmy(i0,j0) = .true.
	     uxbar(i0,j0) = magvel(i0,j0)
     &       * sin(3.1415926/180.*azvel(i0,j0))
	     uybar(i0,j0) = magvel(i0,j0)
     &       * cos(3.1415926/180.*azvel(i0,j0))
	     end do
          end do
 12       close(11)

	  open(13,file='inlets.dat',err=14)
	  do i=1,nx
	     do j=1,ny
	     read(13,*,err=14) j0,i0,azi,magv
	     climmx(i0,j0) = .true.
	     climmy(i0,j0) = .true.
	     uxbar(i0,j0) = magv
     &       * sin(3.1415926/180.*azi)
	     uybar(i0,j0) = magv
     &       * cos(3.1415926/180.*azi)
	     end do
          end do
 14       close(13)

          RETURN
          END

