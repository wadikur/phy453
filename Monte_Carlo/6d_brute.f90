program MultiBrute
      !g calculate gaussian integration in the range [-5,5]
      implicit none
      real*8::r(6),func,length,volume,sums,sq_sums,mean,sigma,results
      integer:: n,i
      

      length=5.0d0
      volume=(2.0d0*length)**6
      n=100
      open(21,file="6d_brute.dat")



1     sums=0.0d0
      sq_sums=0.0d0
      do i=1,n
        call random_number(r)
        r=2.0d0*length*r-length

        sums=sums+func(r)
        sq_sums=sq_sums+func(r)*func(r)
      end do

      mean=sums/real(n)
      sq_sums=sq_sums/real(n)
      sigma=sqrt((sq_sums-mean*mean)/real(n))
      
      results=volume*mean
      sigma=volume*sigma
      write(*,*)n, results,sigma
      write(21,*) n, results, sigma
      
      n=n*10

      if (n .lt. 1e9) goto 1



end program MultiBrute

function func(x) result(y)
      real*8 :: x(6),y,xx,yy,xy,b
      b=0.5d0

      xx=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
      yy=x(4)*x(4)+x(5)*x(5)+x(6)*x(6)
      xy=(x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2
      y=exp(-xx-yy-b*xy)
end function func

