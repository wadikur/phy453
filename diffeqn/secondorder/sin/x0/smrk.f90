program SMRK4
      implicit none
      
      real*8::h,t,x,v,fx0,x1,fv0,v1,fx1,x2,fv1,v2,fx2,x3,fv2,v3,fx3,fv3,E
      integer:: i,num_iteration



      open(21,file='energy.dat')
      open(22,file='x.dat')
      open(23,file='v.dat')
     

      h=0.01
      num_iteration=5000
      t=0.0d0
      x=0.0d0
      v=1.9990d0 

      do i=1,num_iteration 

        fx0=v; x1=x+(h*fx0)/2.0d0

        fv0=-sin(x); v1=v+(h*fv0)/2.0d0

        fx1=v1; x2=x+(h*fx1)/2.0d0

        fv1=-sin(x1); v2=v+(h*fv1)/2.0d0

        fx2=v2; x3=x+h*fx2

        fv2=-sin(x2); v3=v+h*fv2

        fx3=v3;
        fv3=-sin(x3)  

        x=x+h*(fx0+2*fx1+2*fx2+fx3)/6.0d0
        v=v+h*(fv0+2*fv1+2*fv2+fv3)/6.0d0

        E=((v*v)/2.0d0)-cos(x)
        
        t=t+h
        write(21,*)t,E
        write(22,*)t,x
        write(23,*)t,v
        if (i==num_iteration) then
                print*,i,'t=',t,'x=',x,'E=',E
        end if
      end do

      close(21)
end program SMRK4
