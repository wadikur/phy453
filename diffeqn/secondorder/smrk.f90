program SMRK4
      implicit none
      
      real*8::h,k,time_range,t,x,v,fx0,x1,fv0,v1,fx1,x2,fv1,v2,fx2,x3,fv2,v3,fx3,fv3,E



      open(21,file='energy.dat')
      open(22,file='x.dat')
      open(23,file='v.dat')
     

      h=0.02
      k=1.0d0
      time_range=20.0d0
      t=0.0d0
      x=0.0d0
      v=0.1d0 

      do while (t<=time_range)

        fx0=v; x1=x+(h*fx0)/2.0d0

        fv0=-k*x; v1=v+(h*fv0)/2.0d0

        fx1=v1; x2=x+(h*fx1)/2.0d0

        fv1=-k*x1; v2=v+(h*fv1)/2.0d0

        fx2=v2; x3=x+h*fx2

        fv2=-k*x2; v3=v+h*fv2

        fx3=v3;
        fv3=-k*x3  

        x=x+h*(fx0+2*fx1+2*fx2+fx3)/6.0d0
        v=v+h*(fv0+2*fv1+2*fv2+fv3)/6.0d0

        E=(v*v+k*x*x)/2.0d0
        
        t=t+h
        write(21,*)t,E
        write(22,*)t,x
        write(23,*)t,v
      end do

      close(21)
end program SMRK4
