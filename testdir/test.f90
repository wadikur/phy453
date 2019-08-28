program test
      implicit none

      character(len=32)::filename

      integer::n=10,i

      do i=1,n

      write(filename,*) n

      print*,filename

      end do
      end program
