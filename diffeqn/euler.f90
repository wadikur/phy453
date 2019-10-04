program Euler
      implicit none
      
      real*8::x_i,y_i,x,y,fxy,h
      integer::i,num_iteration



      x_i=0.0d0
      y_i=0.0d0

      h=0.02

      num_iteration=int(1.55/h)

      open(unit=21,file="Euler_h_0_02.dat")
      y=y_i
      x=x_i
      write(21,*)x,y
      do i=1,num_iteration
        x=x+h
        fxy=1+y*y
        y=y+h*fxy

        write(21,*)x,y
      end do
      close(21)
end program Euler
