module md_module
  implicit none
  integer,parameter:: npart= 1200, niter_eq=10000, niter=100000, n_calc_av=10
  real*8, parameter:: mass=1.0d0, temp= 1.0d0, sigma=1.0d0,eps=1.0d0,dt=0.005d0
  integer,parameter :: lx=20,ly=20,lz=20
  real*8,parameter:: llx=dfloat(lx), lly=dfloat(ly), llz=dfloat(lz),npartf=dfloat(npart)
  real*8,parameter:: llxby2=llx/2.0d0,llyby2=lly/2.0d0,llzby2=llz/2.0d0
  real*8,parameter:: dt2by2=0.50d0*dt**2,dtby2=0.50d0*dt
  real*8,parameter:: rc=2.50d0*sigma
  real*8,parameter:: sigma6=sigma**6,sigma12=sigma**12,epst4=4*eps
  real*8,parameter:: fc=epst4*((12.0d0*sigma12/(rc**13))-(6.0d0*sigma6/(rc**7)))
  real*8,parameter:: vfc=fc*rc+epst4*(((sigma/rc)**12)-((sigma/rc)**6))

  real*8:: vel_const=dsqrt(12*temp/mass),av_vx=0.0d0,av_vy=0.0d0,av_vz=0.0d0
  real*8:: pos(3*npart),vel(3*npart),force(3*npart),old_force(3*npart)
  real*8:: potential_energy,kinetic_energy
  integer:: i,j

end module md_module


program MD
  use md_module
  implicit none
  integer:: time

  open(31,file="energy.dat")
  open(32,file="kinetic_energy.dat")
  open(33,file="potential_energy.dat")
  print*,"Initializing positions"
  call pos_init
  print*, "Positions initialized"
  print*,"Initializing velocity"
  call vel_init
  print*, "Velocity initialized"
  call calc_force
  do time=1,niter_eq
     call update_pos
     call calc_force
     call update_vel
  end do
  do time=1,niter
     write(31,*) kinetic_energy+potential_energy
     write(32,*) kinetic_energy
     write(33,*) potential_energy
     print*, "No_of_iteration>",time
     print*, "momentum along x",av_vx
     print*, "momentum along y",av_vy
     print*, "momentum along z",av_vz
     call update_pos
     call update_force
     call update_vel
     ! if (mod(time,n_calc_av)==0) call calc_thermodyn
  end do

end program MD
subroutine init_pos
  use md_module
  implicit none
  real*8:: rno(3)
  real*8:: r,x,y,z
  integer:: condition

  call random_number(rno)
  pos(1)=rno(1);pos(2)=rno(2);pos(3)=rno(3)

  do i=2,npart
     do
        call random_number(rno)
        pos(3*i-2)= llx*rno(1)
        pos(3*i-1)= lly*rno(2)
        pos(3*i)= llz*rno(3)

        condition=1
        do j=1,i-1
           x=pos(3*i-2)-pos(3*j-2)
           y=pos(3*i-1)-pos(3*j-1)
           z=pos(3*i)-pos(3*j)

           if (abs(x).gt.llxby2)x=(llx-abs(x))*(-1.0d0*x/abs(x))
           if (abs(y).gt.llyby2)y=(lly-abs(y))*(-1.0d0*y/abs(y))
           if (abs(z).gt.llzby2)z=(llx-abs(z))*(-1.0d0*z/abs(z))
           r=dsqrt(x*x+y*y+z*z)
           if (r .lt. 1.50d0*sigma) condition=0
           if(condition==0)exit

        end do
        if (condition==1)exit
     end do
  end do

end subroutine init_pos
subroutine init_vel
  use md_module
  implicit none
  real*8:: rno(3)
  kinetic_energy=0.0d0
  av_vx=0.0d0
  av_vy=0.0d0
  av_vz=0.0d0

  ! velocity initializing through random number
  do i=1,npart
     call random_number(rno)
     rno=vel_const*(rno-0.50d0)
     vel(3*i-2)=rno(1)
     vel(3*i-1)=rno(2)
     vel(3*i)=rno(3)
     av_vx=vel(3*i-2)+av_vx
     av_vy=vel(3*i-1)+av_vy
     av_vz=vel(3*i)+av_vz
  end do
  av_vx=av_vx/npartf
  av_vy=av_vy/npartf
  av_vz=av_vz/npartf
  do i=1,npart
     vel(3*i-2)=vel(3*i-2)-av_vx
     vel(3*i-1)=vel(3*i-1)-av_vy
     vel(3*i)=vel(3*i)-av_vz
     kinetic_energy=kinetic_energy+vel(3*i-2)*vel(3*i-2)+vel(3*i-1)*vel(3*i-1)+vel(3*i)*vel(3*i)
  end do
  kinetic_energy=mass*kinetic_energy

end subroutine init_vel
subroutine update_pos
  use md_module
  do i=1,npart
     pos(3*i-2)=pos(3*i-2)+vel(3*i-2)*dt+dt2by2*force(3*i-2)/mass
     pos(3*i-1)=pos(3*i-1)+vel(3*i-1)*dt+dt2by2*force(3*i-1)/mass
     pos(3*i)=pos(3*i)+vel(3*i)*dt+dt2by2*force(3*i)/mass
     !PBC
     pos(3*i-2)=modulo(pos(3*i-2),llx)
     pos(3*i-1)=modulo(pos(3*i-1),llx)
     pos(3*i)=modulo(pos(3*i),llx)
  end do
end subroutine update_pos
subroutine update_vel
  use md_module
  kinetic_energy=0.0d0
  av_vx=0.0d0
  av_vy=0.0d0
  av_vz=0.0d0

  do i=1,npart
     vel(3*i-2)=vel(3*i-2)+dtby2*(old_force(3*i-2)+force(3*i-2))/mass
     vel(3*i-1)=vel(3*i-1)+dtby2*(old_force(3*i-1)+force(3*i-1))/mass
     vel(3*i)=vel(3*i)+dtby2*(old_force(3*i)+force(3*i))/mass
     av_vx=vel(3*i-2)+av_vx
     av_vy=vel(3*i-1)+av_vy
     av_vz=vel(3*i)+av_vz
     kinetic_energy=kinetic_energy+vel(3*i-2)*vel(3*i-2)+vel(3*i-1)*vel(3*i-1)+vel(3*i)*vel(3*i)
  end do
  av_vx=av_vx/npartf
  av_vy=av_vy/npartf
  av_vz=av_vz/npartf
  kinetic_energy=mass*kinetic_energy

end subroutine update_vel
subroutine calc_force
  use md_module
  implicit none
  real*8:: x=0.0d0,y=0.0d0,z=0.0d0,r=0.0d0
  real*8:: pair_pot=0.0d0,pair_force=0.0d0
  old_force=force
  force=0.0d0;potential_energy=0.0d0
  do i=1,npart-1
     do j=i,npart
        x=pos(3*i-2)-pos(3*j-2)
        y=pos(3*i-1)-pos(3*j-1)
        z=pos(3*i)-pos(3*j)

        if (abs(x).gt.llxby2)x=(llx-abs(x))*(-1.0d0*x/abs(x))
        if (abs(y).gt.llyby2)y=(lly-abs(y))*(-1.0d0*y/abs(y))
        if (abs(z).gt.llzby2)z=(llx-abs(z))*(-1.0d0*z/abs(z))
        r=dsqrt(x*x+y*y+z*z)
        if (r<=rc) then
           pair_pot= epst4*(((sigma/r)**12)-((sigma/r)**6))-vfc+fc*r
           potential_energy=potential_energy+pair_pot

           pair_force=epst4*((12.0d0*sigma12/(r)**13)-(6.0d0*sigma6/(r)**7))-fc
           force(3*i-2)=force(3*i-2)+pair_force*(x/r)
           force(3*i-1)=force(3*i-1)+pair_force*(y/r)
           force(3*i)=force(3*i)+pair_force*(z/r)
           force(3*j-2)=force(3*j-2)-pair_force*(x/r)
           force(3*j-1)=force(3*j-1)-pair_force*(y/r)
           force(3*j)=force(3*j)-pair_force*(z/r)
        end if
     end do
  end do
end subroutine calc_force
