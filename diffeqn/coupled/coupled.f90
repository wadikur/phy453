program Coupled
      implicit none

      integer,parameter:: nop=50
      
      real*8::y(nop),yl(nop),yr(nop),y1(nop),y1l(nop),y1r(nop),y2(nop),y2l(nop),y2r(nop),y3(nop),y3l(nop),y3r(nop)
      real*8::f0(nop),f1(nop),f2(nop),f3(nop),f0v(nop),f1v(nop),f2v(nop),f3v(nop)
      real*8::v(nop),v1(nop),v2(nop),v3(nop)
      real*8::h,t,k,hb2,twopi=8.0d0*atan(1.0d0)
      integer:: num_iteration,i,j,r=5

      open(21,file="jmol_data.xyz")
      open(22,file="postiton_1st.dat")


      h=0.02d0
      t=0.0d0
      num_iteration=5000
      k=1.0d0
      hb2=h/2.0d0

      !Initial Condition
      y=0.0d0
      y(1)=0.8d0;y(26)=0.8d0
      v=0.0d0

      write(21,*)"50"
      write(21,*)""
      do j=1,nop
       if (j==1) then
       write(21,*)"O",r*cos(twopi/50),r*sin(twopi/50),y(1)
       elseif (j==26) then
               write(21,*) "O",r*cos(26*twopi/50.0d0),sin(26*twopi/50.0d0),y(26)
       else
               write(21,*) "N",r*cos(j*twopi/50.0d0),r*sin(j*twopi/50.0d0),y(j)
       end if
      end do
      !LOOP
      do i=1,num_iteration
      yl=cshift(y,-1)
      yr=cshift(y,1)

      f0=v; y1=y+f0*hb2
      y1l=cshift(y1,-1);y1r=cshift(y1,1)

      f0v=k*(yl+yr-2.0d0*y); v1=v+f0v*hb2

      f1=v1; y2=y+f1*hb2;
      y2l=cshift(y2,-1);y2r=cshift(y2,1)

      f1v=k*(y1l+y1r-2.0d0*y1); v2=v+f1v*hb2

      f2=v2; y3=y+f2*h
      y3l=cshift(y3,-1);y3r=cshift(y3,1)

      f2v=k*(y2l+y2r-2.0d0*y2); v3=v+f2v*h

      f3=v3
      f3v=k*(y3l+y3r-2.0d0*y3)

      y=y+h*(f0+2.0d0*f1+2.0d0*f2+f3)/6.0d0
      v=v+h*(f0v+2.0d0*f1v+2.0d0*f2v+f3v)/6.0d0

      t=t+h

      write(22,*)t,y(1)
      
      if (i==2000) then
              print*, i,"t=",t,"position of the first particle=",y(1)
      end if
      
      if (mod(i,10)==0) then
      write(21,*)"50"
      write(21,*)""

      do j=1,nop
       if (j==1) then
       write(21,*)"O",r*cos(twopi/50),r*sin(twopi/50),y(1)
       elseif (j==26) then
               write(21,*) "O",r*cos(26*twopi/50.0d0),r*sin(26*twopi/50.0d0),y(26)
       else
               write(21,*) "N",r*cos(j*twopi/50.0d0),r*sin(j*twopi/50.0d0),y(j)
       end if
       end do
       end if
      end do

end program Coupled
