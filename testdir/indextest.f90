program test
      implicit none
      real*8::r(3)

      call random_number(r)

      print*,r(3),r(1),r(2)
      end program test
