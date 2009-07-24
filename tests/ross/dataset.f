     
      Subroutine ReadDataset()

      implicit none
      include 'shelfy.h'

      integer aux(NX)
      Double precision auxr(NX)

      open (10,file='111by147Grid.dat')
    
      do i=1,nx+ny+10
         read(10,*) 
      end do

      do i=1,NY
         read(10,*) aux
	 do j=1,NX
	 mask(j,i) = aux(j)
	 end do
      end do

      read(10,*) 
      read(10,*)

      do i=1,NY
         read(10,*) auxr
	 do j=1,NX
	 azvel(j,i) = auxr(j)
	 end do
      end do

      read(10,*)
      read(10,*)

      do i=1,NY
         read(10,*) auxr
         do j=1,NX
         magvel(j,i) = auxr(j)
         end do
      end do

      read(10,*)
      read(10,*)

      do i=1,NY
         read(10,*) auxr
         do j=1,NX
         H(j,i) = auxr(j)
         end do
      end do

      read(10,*)
      read(10,*)

      do i=1,NY
         read(10,*) auxr
         do j=1,NX
         velobs(j,i) = auxr(j)
         end do
      end do

      read(10,*)
      read(10,*)

      do i=1,NY
         read(10,*) auxr
         do j=1,NX
         seabed(j,i) = auxr(j)
         end do
      end do

      read(10,*)
      read(10,*)

      do i=1,NY
         read(10,*) auxr
         do j=1,NX
         icefront(j,i) = auxr(j)
         end do
      end do

      read(10,*)
      read(10,*)

      do i=1,NY
         read(10,*) auxr
         do j=1,NX
         accumulation(j,i) = auxr(j)
	 auxr(j)=0.
         end do
      end do

      read(10,*)
      read(10,*)

      do i=1,NY
         read(10,*) auxr
         do j=1,NX
         bbar(j,i) = auxr(j)/((365.*3600.*24.)**(1./3.))
c        bbar(j,i) = 1.9E8/((365.*3600.*24.)**(1./3.))
         end do
      end do

      read(10,*)
      read(10,*)

      do i=1,NY
         read(10,*) auxr
         do j=1,NX
         temperature(j,i) = auxr(j)*mask(j,i)
         end do
      end do

      close(10)
      
      return
      end
