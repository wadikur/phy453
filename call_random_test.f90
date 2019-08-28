! Integration of the function e^x from 0 to 3 using Monte Carlo method

Program Montecarlo
  implicit none
  real*8 :: x(10001),randDivs(10001),f(10001), sum, area(10001), error(10001)
  integer*8 :: i,j,k,n


  call random_number(x)
  
  randDivs = (x - 0.5)*3/0.5

  
!!$  do i= 1, 1000
!!$     do j= 2, 1001
!!$        if(randDivs(i)>randDivs(j)) then
!!$           swap( randDivs(i), randDivs(j))
!!$        end if
!!$     end do
!!$  end do

  do j=1,100
     n =j*100
        do k=1, n
        write(24,*) x(k), randDivs(k)
         end do

  sum =0;
  do i=1, n
     f(i) = exp(randDivs(i))
     sum = sum +f(i)
  end do
  
  
  area(j) = (sum/n)*10
  error(j)=exp(3.0d0)-1.0d0-area(j)

  write(13,*) n, area(j),error(j)
  end do
end program Montecarlo
