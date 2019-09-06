program Crude
      implicit none
      real*8:: sums,sq_sums,r,sd,results,actual=exp(3.0d0)-1.0d0,error
      integer:: n,i

      write(*,*) "Enter the value of n"
      read(*,*) n
      open(21,file="montecarlo_crude.dat")

1     sums=0.d0
      sq_sums=0.0d0

      do i=1,n 
        call random_number(r)
        r=3.0d0*r

        sums=sums+exp(r)
        sq_sums= sq_sums + exp(r)*exp(r)
      end do

      sums=sums/real(n)
      sq_sums=sq_sums/real(n)
      sd=3.0d0*sqrt((sq_sums-sums**2)/n)
      results= 3.0d0*sums
      error=actual-results


      write(*,*) n, results, sd
      write(21,*) 1.0d0/sqrt(real(n)),error

      n=n*10
      if (n<=1e8) goto 1
end program Crude
