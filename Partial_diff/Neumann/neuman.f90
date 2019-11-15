program Neumann
      implicit none
      integer,parameter::Lx=34,Ly=34
      real*8::old_Temp(1:Lx,1:Ly),Temp(1:Lx,1:Ly)
      integer::i,j,counter,test
      real*8::dx,dy,d1,reference
      real*8::A(Ly),B(Ly),C(Lx),D(Lx)

      !  open(33,file='initial_neuman_.001.dat')
      open(35,file='final_neuman.dat')

      A=-70.0d0;B=-40.0d0;C=20.0d0;D=-10.0d0

     old_Temp=0.0d0;dx=1.0d0;dy=1.0d0
      ! D1=(0.5d0*dx*dx*dy*dy)/(dx*dx+dy*dy)

     counter=0
      do
      counter=counter+1
      test=0

      do j=2,Ly-1
         Temp(1,j)=0.25d0*(2.0d0*old_Temp(2,j)-2.0d0*dx*A(j)+old_Temp(1,j+1)+old_Temp(1,j-1))
         Temp(Lx,j)=0.25d0*(2.0d0*old_Temp(Lx-1,j)+2.0d0*dx*B(j)+old_Temp(Lx,j+1)+old_Temp(Lx,j-1))
      end do


      do i=2,Lx-1
         Temp(i,1)=0.25d0*(old_Temp(i+1,1)+old_Temp(i-1,1)+2.0d0*old_Temp(i,2)-2.0d0*dx*C(i))
         Temp(i,Ly)=0.25d0*(old_Temp(i+1,Ly)+old_Temp(i-1,Ly)+2.0d0*old_Temp(i,Ly-1)+2.0d0*dx*D(i))
      end do




      Temp(1,1)=0.5d0*(old_Temp(1,2)-dx*C(1)+old_Temp(2,1)-dx*A(1))
      Temp(1,Ly)=0.5d0*(old_Temp(1,Ly-1)+dx*D(1)+old_Temp(2,Ly)-dx*A(Ly))
      Temp(Lx,1)=0.5d0*(old_Temp(Lx-1,1)+dx*B(1)+old_Temp(Lx,2)-dx*C(Lx))
      Temp(Lx,Ly)=0.5d0*(old_Temp(Lx-1,Ly)+dx*B(Ly)+old_Temp(Lx,Ly-1)+dx*D(Lx))

      do i=2,Lx-1
        do j=2,Ly-1
         Temp(i,j)=0.25d0*(old_Temp(i-1,j)+old_Temp(i+1,j)+old_Temp(i,j-1)+old_Temp(i,j+1))
         end do
      end do

      do i=1,Lx
        do j=1,Ly
          if(abs(Temp(i,j)-old_Temp(i,j)) .gt. 0.001d0)test=1
          if (abs(Temp(i,j)-old_Temp(i,j)) .gt. 0.1) print*,"greater than 0.1"
        end do
      end do

      print*,counter
      if(test==0)exit


      old_Temp=Temp

     ! do n=1,Ly
      !  do m=1,Lx
       !  write(35,*) m,n, T(m,n)
        ! end do
     ! end do

  end do
reference=2000-Temp(1,1)
     do i=1,Lx
        do j=1,Ly
         write(35,*) i,j, Temp(i,j)+reference
         end do
      end do

      write(*,*) 'number of iteration=',counter
end program Neumann
