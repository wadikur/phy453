program RK4
      implicit none
      real*8::xi,yi,x,y,fxy,y1,f1,y2,f2,y3,f3,h
      integer:: i,num_iteration


      xi=0.0d0
      yi=0.0d0
      h=0.02d0
      num_iteration=int(1.55/h)
      open(21,file="RK4_0_02.dat")
      x=xi
      y=yi
      write(21,*)x,y

      do i=1,num_iteration 
        fxy=1+y*y

        y1=y+(h*fxy)/2.0d0

        f1=1+y1*y1

        y2=y+(h*f1)/2.0d0

        f2=1+y2*y2

        y3=y+h*f2

        f3=1+y3*y3

        y=y+h*(fxy+2*f1+2*f2+f3)/6.0d0

        x=x+h

        write(21,*)x,y
      end do
      close(21)
end program RK4

