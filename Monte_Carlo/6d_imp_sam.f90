program Imp6D
      implicit none

      real*8:: func,gaussrandom,mu=0.0d0,sd=1.0d0/sqrt(2.0d0),sums,sq_sums,mean,sigma,r(6),volume

      real*8:: results
      integer:: i,n,j

      n=10 
      volume=acos(-1.0d0)**3


      sums=0d0
      sq_sums=0d0
      

1     do i=1,n
        do j=1,6
                r(j)=gaussrandom(mu,sd)
        end do
        sums=sums+func(r)
        sq_sums=sq_sums+func(r)*func(r)
      end do
      mean=sums/real(n)
      sq_sums=sq_sums/real(n)
      sigma=sqrt((sq_sums-mean*mean)/real(n))

      results=volume*mean
      sigma=volume*sigma

      write(*,*) n,results,sigma

      n=n*10
      if (n .lt. 1e9) goto 1



end program Imp6D


function func(x) result(y)
        real*8:: x(6),xy,y,b=0.5d0
        xy=(x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2

        y=exp(-b*xy)
end function func


function gaussrandom(mu,sd) result(g)
        real*8:: mu,sd,g,u,v

3       call random_number(u)
        call random_number(v)
        u=2*u-1.0d0
        v=2*v-1.0d0
        g=u*u+v*v
        if (g==0.0d0 .or. g>=1.0d0) goto 3
        g=sqrt(-2*log(g)/g)
        g=u*g
        g=sd*g+mu
end function gaussrandom


      
