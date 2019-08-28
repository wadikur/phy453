program Exponential
      implicit none

      real*8:: r
      integer, allocatable,dimension(:) ::bin
      integer:: n,i,nbin

      write(*,*) "Enter the value of n"
      read(*,*) n
      write(*,*) "Enter the number of bin"
      read(*,*) nbin

      open(21,file="exponential_histogram.dat")
      open(22,file="exponential_data.dat")

      allocate(bin(0:nbin))


      bin=0
      do i=1,n
        call random_number(r)
        r=-log(1-r)                 ! pdf=e**(-x) from (0,infinity)

        if (abs(int(r))<=nbin) then
                bin(int(r))=bin(int(r))+1
                write(22,*) r

        end if
      end do
      write(*,*) bin

      do i= 0,nbin
        write(21,*) i, bin(i)
      end do
      call system('gnuplot exponential_histogram.plt')
      call system('feh exponential_histogram.png')

end program Exponential

        
        
