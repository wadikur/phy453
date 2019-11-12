PROGRAM VARIATIONAL

  IMPLICIT NONE
  INTEGER,PARAMETER::dp = selected_real_kind(14,200)
  REAL,PARAMETER:: pi = 3.14159265359_dp
  REAL(dp),ALLOCATABLE::H(:,:),kn(:),Energy(:),Work(:)
  REAL(dp)::a,b,V0 ,dx,x,norm,prob
  COMPLEX(dp)::psi
  INTEGER::i,j,n,nn,lwork,info,n_a,n_b,w

open(38,file="Energy_vs_a.dat")
do w=1,10
  N = 10
  nn = 2*n+1
  a = 5.0*w
  V0 = 2.0
  b = 4.0
  allocate(kn(nn),H(nn,nn),Energy(nn),work(3*nn))
  !allowed values of momenta are n*2*pi/2a
  kn(1) = 0.0_dp
  do i = 2,nn-1,2
     kn(i) = (i/2)*2.0_dp*pi/a
     kn(i+1) = -(i/2)*2.0_dp*pi/a
  end do

  ! Hamiltonian Matrix
  do i = 1,nn
     do j = 1,nn
        If(i==j) then
           H(i,j) = kn(i)**2 - V0*b/a
        else
           H(i,j) = -(V0/a)*2.0_dp*sin( (kn(j)-kn(i))*b/2.0_dp )/ (kn(j)-kn(i))
        endif
     enddo
  enddo

  lwork = 3*nn
  call dsyev('V','U',nn,H,nn,Energy,work,lwork,info)!Output to Energy array ascending order
  if (info/=0)stop

  print*,"Enerygy",Energy(1)
  write(38,*) a,Energy(1)




  dx = 0.01_dp!grid size
  n_a = nint(a*0.50_dp/dx)!# points on x>0 space
  n_b = nint(b*0.50_dp/dx)
  norm = 0.0d0
  x = 0.0_dp

  open(1,file='variational_wave3.dat',form='formatted',status='unknown',action='write')
  open(2,file='variational_potential3.dat',form='formatted',status='unknown',action='write')

  do i = -n_a,-n_b
     x = i*dx
     write(2,*)x,V0
  enddo
  do i = -n_b,n_b
     x = i*dx
     write(2,*)x,0.0
  enddo
  do i = n_b,n_a
     x = i*dx
     write(2,*)x,V0
  enddo

  x = 0.0
  do i = -n_a,n_a
     x = dx*i
     psi = 0.0
     do j =1,nn
        Psi = Psi + H(j,1)*exp((0.0,1.0)*kn(j)*x)/sqrt(a)!ground state
     enddo
     prob = Psi*conjg(Psi)
     norm = norm + prob*dx
     !print*,i,x,prob
     write(1,*)x,prob,Psi
  enddo
  deallocate(H,kn,Energy,work)




end do



END PROGRAM VARIATIONAL
