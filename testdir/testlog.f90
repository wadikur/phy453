program testlog
      implicit none

      real*8:: r,n=10.0d0,pi=4*atan(1.0d0)


1      r=log(real(n))
      write(*,*) n,r
      n=n*10
      if (n<1e15) goto 1
      write(*,*) pi
end program testlog
