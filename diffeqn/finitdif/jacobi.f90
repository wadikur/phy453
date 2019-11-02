program FinitDifference
      implicit none
      real*8::f_y,end_y
      real*8,parameter::dx=0.01d0,f_x=0.0d0,end_x=0.9d0
      integer,parameter::nop=int((end_x-f_x)/dx)+1
      integer::i,ll,cond
      real*8::x(nop),y(nop),y_old(nop),limit
     real*8,parameter::d1=1.0d0/(2.0d0-10.0d0*dx*dx),d2=(1.0d0-2.5d0*dx)
      real*8,parameter::d3=(1.0d0+2.5d0*dx),d4=-10.0d0*dx*dx
      
      limit=0.0001d0
      ll=0
      cond=0
      open(33,file='initial_value.dat')
      open(30,file='dx_0.01_limit_10_4.dat')
      
      x(1)=f_x;x(nop)=end_x
      y(1)=0.0d0;y(nop)=50
      f_y=y(1);end_y=y(nop)

      do i=2,nop-1
        x(i)=x(i-1)+dx
        y(i)=(end_y-f_y)*x(i)/(end_x-f_x)

        write(33,*)x(i),y(i)
      end do
      
      do
        ll=ll+1
        if(cond==1)exit

        y_old=y

        do i=2,nop-1
           y(i)=d1*(d2*y_old(i+1)+d3*y_old(i-1)+d4*x(i))
        end do

        cond=1

        do i=2,nop-1
        if(abs(y_old(i)-y(i)) .ge. limit)cond=0
        end do

      end do

        print*, 'no of iteration',ll
        print*, "At x=",x(81),"y=",y(81)
        do i=1,nop
        write(30,*)x(i),y(i)
        end do
        close(30)
end program FinitDifference
