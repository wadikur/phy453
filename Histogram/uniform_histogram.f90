program Uniform
      implicit none

      real*8:: r
      integer,allocatable,dimension(:):: bin
      integer:: n,i,nbin

      write(*,*) "Enter the value of n"
      read(*,*) n
      write(*,*) "Enter number of bin"
      read(*,*) nbin

      open(21,file="uniform_histogram.dat")
      open(22,file="uniform_data.dat")

      allocate(bin(nbin))

      bin=0
      do i= 1,n
        call random_number(r)
        write(22,*) r
        r=nbin*r
        bin(int(r)+1)=bin(int(r)+1)+1
      end do
      write(*,*) bin
      
      do i=1,nbin
        
        write(21,*) i,bin(i)
      end do

      call system('gnuplot uniform_histogram.plt')
      call system('feh uniform_histogram.png')

end program Uniform


