module md_module
  implicit none
  integer,parameter :: lx=20,ly=20,lz=20
  real*8,parameter:: llx=dfloat(lx), lly=dfloat(ly), llz=dfloat(lz)
  real*8,parameter:: llxby2=llx/2.0d0,llyby2=lly/2.0d0,llzby2=llz/2.0d0
  integer, parameter:: npart= 1200, niter_eq=5000, niter=50000, n_calc_av=10
  real*8, parameter:: mass=1.0d0, temp= 1.0d0, sigma=1.0d0,eps=1.0d0
  real*8, parameter:: rc=2.50d0*sigma
  real*8, parameter:: sigma6=sigma**6,sigma12=sigma**12,epst4=4*eps
  real*8, parameter:: fc=epst4*((12.0d0*sigma12/(rc**13))-(6.0d0*sigma6/(rc**7)))
  real*8, parameter:: vfc=fc*rc+epst4*(((sigma/rc)**12)-((sigma/rc)**6))

  real*8:: pos(3*npart),vel(3*npart),force(3*npart),old_force(3*npart)
  integer:: i,j

contains
  
end module md_module

program MD
  use md_module
  implicit none

  open(21,file="test.dat")
  call init_pos()
  do i=1,3*npart
     write(21,*) pos(i)
  end do

  close(21)
end program MD
