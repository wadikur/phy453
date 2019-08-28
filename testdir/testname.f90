program foo
      implicit none
      character(len=1024):: filename
      write(filename, "(A5,I2)") "hello", 10
      print*, trim(filename)
end program
