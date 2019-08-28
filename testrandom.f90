program Test
      implicit none
      real:: r
      integer :: n=6,i

      do i=1,n
        call random_number(r)
        print*, r
      end do
end program Test
