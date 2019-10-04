program impEuler
      implicit none

      real*8::x_i,y_i,x,y,fxy,h,f_i,y_temp,f_e
      integer::i,num_iteration



      x_i=0.0d0
      y_i=0.0d0

      h=0.02

      num_iteration=int(1.55/h)

      open(unit=21,file="improvedEuler_h_0_02.dat")
      y=y_i
      x=x_i
      write(21,*)x,y

      do i=1,num_iteration
        f_i=1+y*y
        y_temp=y+h*f_i
        f_e=1+y_temp*y_temp

        y=y+h*(f_i+f_e)/2.0d0
        x=x+h

        write(21,*)x,y
      end do
      close(21)
end program impEuler
