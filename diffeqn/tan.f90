program tanx
      implicit none
      real*8::x,h
      integer:: num_iteration,i

      h=0.0001
      num_iteration=int(1.55/h)
      x=0.0d0
      open(21,file="tan.dat")
      do i=0,num_iteration
        
        write(21,*)x,tan(x)
        x=x+h
      end do

      close(21)

end program tanx
