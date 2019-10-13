program test
      implicit none

      integer::i,j

      open(21,file="jmol_test.xyz")
      do j=1,100
      write(21,*)"20"
      do i=1,20
      write(21,*) "N", 5*cos(i*8.0d0*atan(1.0d0)/20), 5*sin(i*8.0d0*atan(1.0d0)/20),mod(j*i,100) 
      end do
      end do

      end program



