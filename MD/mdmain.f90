module md_module
  implicit none
  integer,parameter :: lx=20,ly=20,lz=20
  real*8,parameter:: llx=dfloat(lx), lly=dfloat(ly), llz=dfloat(lz)
  real*8,parameter:: llxby2=llx/2.0d0,llyby2=lly/2.0d0,llzby2=llz/2.0d0
  integer,parameter:: npart= 1200, niter_eq=5000, niter=50000, n_calc_av=10
  real*8, parameter:: mass=1.0d0, temp= 1.0d0, sigma=1.0d0,eps=1.0d0
  real*8, parameter:: rc=2.50d0*sigma
  real*8, parameter:: sigma6=sigma**6,sigma12=sigma**12,epst4=4*eps
  real*8, parameter:: fc=epst4*((12.0d0*sigma12/(rc**13))-(6.0d0*sigma6/(rc**7)))
  real*8, parameter:: vfc=fc*rc+epst4*(((sigma/rc)**12)-((sigma/rc)**6))

  real*8:: pos(3*npart),vel(3*npart),force(3*npart),old_force(3*npart)
  integer:: i,j


contains
  ! Initialize position randomly
  subroutine init_pos()
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
             write(*,*) i,j, r,condition
             if(condition==0)exit

          end do
          if (condition==1)exit
       end do
    end do

  end subroutine init_pos


  subroutine init_vel()
    implicit none
    real*8,private:: vel_const,av_vx=0.0d0,av_vy=0.0d0,av_vz=0.0d0
    real*8,private:: rno(3)

    vel_const=dsqrt(12.0d0)*dfloat(temp)

    ! velocity initializing through random number
    do i=1,npart
       call random_number(rno)
       rno=vel_const*(rno-0.5d0)
       vel(3*i-2)=rno(1)
       vel(3*i-1)=rno(2)
       vel(3*i)=rno(3)
       av_vx=vel(3*i-2)+av_vx
       av_vy=vel(3*i-1)+av_vy
       av_vz=vel(3*i)+av_vz
    end do
    av_vx=av_vx/dfloat(npart)
    av_vy=av_vy/dfloat(npart)
    av_vz=av_vz/dfloat(npart)
    do i=1,npart
       vel(3*i-2)=vel(3*i-2)-av_vx
       vel(3*i-1)=vel(3*i-1)-av_vy
       vel(3*i)=vel(3*i)-av_vz
    end do

  end subroutine init_vel()


  subroutine update_pos(delta_t)
    real*8,intent(in):: delta_t,dt2by2,r
    real*8
    dt2by2=0.50d0*delta_t**2
    !updating the position
    do i=1,npart
       pos(3*i-2)=pos(3*i-2)+vel(3*i-2)*delta_t+dt2by2*force(3*i-2)
       pos(3*i-1)=pos(3*i-1)+vel(3*i-1)*delta_t+dt2by2*force(3*i-1)
       pos(3*i)=pos(3*i)+vel(3*i)*delta_t+dt2by2*force(3*i)
       !PBC
       pos(3*i-2)=modulo(pos(3*i-2),llx)
       pos(3*i-1)=modulo(pos(3*i-1),llx)
       pos(3*i)=modulo(pos(3*i),llx)
    end do
  end subroutine update_pos()
  !calculate the force
  subroutine calc_force()
    implicit none
    real*8,private:: x,y,z
    force=0.0d0;new_pot_energy=0.0d0
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
             new_pot_energy=new_pot_energy+pair_pot

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
  end subroutine calc_force()

end subroutine init_pos()
end module md_module


program MD
use md_module
implicit none
integer:: time

call pos_init
call vel_init
call calc_force
do time=1,niter_eq
   call update_pos
end do
do time=1,niter
   call update_pos
   call update_force
   call update_vel
   if (mod(time,n_calc_av)==0) call calc_thermodyn
end do
