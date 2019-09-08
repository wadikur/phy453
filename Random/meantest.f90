program MeanTest
      implicit none

      real*8:: r, mean, sums,abs_error

      integer:: n,i
      open(21,file="meantest.dat")

      n=1000
      sums=0.d0

      do while (n<1e6)
        do i=1,n
                call random_number(r)
                sums=sums+r
        end do
        mean=sums/real(n)
        abs_error=abs(mean-0.5)
        write(21,*) 1/sqrt(real(n)),abs_error
        n=n+1000

      end do
      
end program MeanTest
        

