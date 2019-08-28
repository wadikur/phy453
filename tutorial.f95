program Fibonacci
      implicit none
      integer::x,y,z,n
      real::phi_est

      y=1
      z=1
      n=3
      do while(n<46)
        x=y+z
        phi_est=x/real(y)
        print*,n,"th number is",x
        print '(f10.6)',phi_est
        n=n+1
        z=y
        y=x
      end do
end program Fibonacci
