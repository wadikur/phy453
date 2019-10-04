program Sinkx
      implicit none
      real::h,k,t

      k=1.0d0
      h=0.0001
      open(21,file='Sin.dat')
       
      do while (t<=20.0d0)
        write(21,*)t,0.1d0*sin(k*t) 
        t=t+h

      end do
      close(21)
end program Sinkx
