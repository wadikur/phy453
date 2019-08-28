program Trapazoidal
        !This program finds the value of pi using Trapezoidal integration.
        !The function used is 4/(1+x**2)
      implicit none
      real*8:: h,x0,xn,Area,x,error,f0,fn,func
      integer:: n,i

      x0=0.0d0
      xn=1.0d0
      f0=func(0.0d0)
      fn=func(1.0d0)
      n=1000
      error=100
     open(unit=21,file="error.dat") 
     open(unit=22,file="error_vs_n.dat")
      do while (abs(error)>1e-10)
        n=n+100
        h=(xn-x0)/real(n)
        Area=(h*(f0+fn))/2.0d0
        
        do i=1,n-1
                x=x0+h*i
                Area=Area+func(x)*h
        end do
        
        error=4.0d0*atan(1.0d0)-Area
        print*, n, Area, error

        write(21,*) 1.0d0/(n**2), error
        write(22,*) n,error
      
      end do
        


end program Trapazoidal


function func(x) result(y)
      real*8::x,y
      y=(4.0d0/(1+x**2))
end function func
