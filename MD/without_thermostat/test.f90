program Test
  implicit none
  integer::i
  real*8::npartf=dfloat(1200)
  real*8::kinetic(5000),potential(5000),time

  open(21,file="kin.dat",status='old')
  open(22,file='pot.dat',status='old')
  open(23,file='kin_mod.dat')
  open(24,file='pot_mod.dat')

  do i=1,5000
     read(21,*)time, kinetic(i)
  end do
  do i=1,5000
     read(22,*)time,potential(i)
  end do

  do i=1,5000
     write(23,*) i*0.005d0,kinetic(i)
  end do
  do i=1,5000
     write(24,*) i*0.005d0,potential(i)
  end do
end program Test
