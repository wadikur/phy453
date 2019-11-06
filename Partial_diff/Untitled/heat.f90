program Heat_Equation
  implicit none
  integer::i,j,k,counter,Test
  real*8::dt,tolerance=0.0001d0
  integer,parameter:: Lx=34,Ly=34
  real*8::temp_old(1:Lx,1:Ly),temp(1:Lx,1:Ly)

  temp_old=0.0d0;dt=0.1d0
  open(unit=21,file='initial_condition.dat')
  open(unit=22,file='final_temp_steady.dat')

  temp_old(1,1)=3.7d0;temp(Lx,1)=0.4
  do i=1,Lx
     temp_old(i,1)=temp_old(1,1)-(i-1)*dt
     temp_old(i,Ly)=temp_old(i,1)
     do j=1,Ly
        temp_old(1,j)=temp_old(1,1)
        temp_old(Lx,j)=temp_old(Lx,1)
        write(21,*)i,j,temp_old(i,j)
     end do

  end do

  temp=temp_old



  counter=0

  do
     counter=counter+1


     Test=0


     do i=2,Lx-1
        do j=2,Ly-1

           temp(i,j)=0.25d0*(temp_old(i+1,j)+temp_old(i-1,j)+temp_old(i,j+1)+temp_old(i,j-1))

        end do
     end do

     do i=2,Lx-1
        do j=2,Ly-1
           if(abs(temp(i,j)-temp_old(i,j)) .gt. tolerance)Test=1
        end do
     end do

     if(Test==0)exit
     temp_old=temp

  end do

  print*,'no of steps for convergence=',counter

  do i=1,Lx
     do j=1,Ly
        write(22,*)i,j,temp(i,j)

     end do
  end do

end program Heat_equation

