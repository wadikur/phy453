program HitMiss
      !Program to calculate the integration 0_3 exp(x)dx
      implicit none
      real*8:: realvalue, r(2), results, error
      integer:: n, i, counter


      realvalue=exp(3.0d0)-1
      write(*,*) "Enter the value of n"
      read(*,*) n
      open(21,file="HitMiss_error.dat")

      
1      counter=0
      do i=1, n
        call random_number(r)
        r(1)=3.0d0*r(1)
        r(2)=exp(3.0d0)*r(2)

        if (r(2)<exp(r(1))) then
                counter=counter+1
        end if
      end do

      results=3.0d0*exp(3.0d0)*(real(counter)/real(n))

      error=realvalue-results

      write(*,*) n,results,error
!      write(21,*) 1.0d0/sqrt(real(n)), error
      write(21,*) n, error
      n=n*10
      if (n<=1e8) goto 1
end program HitMiss
