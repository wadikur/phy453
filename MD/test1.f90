program test
  implicit none
  integer:: n
  do
     n=n+1
     if (mod(n,2) == 0) then
        write(*,100,advance='no') "Even"
     else
        ! write(unit=6,fmt='(a)',advance='no' ) "odd"
     end if
  end do
end program test
