program Ising3D
   implicit none

   integer :: i,j,L,a,b,c,d,f,g,ni,tm,N,k,mm,n_equil=10000,n_stat=10,T_temp
   real*8:: r,E,M,h,Mag,Ei,Ef,dE,u,absm,absm_m
   real*8:: T, j_ising=1.0d0
   real*8:: av_m, av_e, cv, av_e2,av_m4,binc, chi, av_m2,av_m_n, av_e_n

   integer, dimension(:,:,:), allocatable::spin
   character(32)::filename
   integer :: munit,eunit,iunit,avmunit,aveunit,abmunit,cvunit,chiunit,bcunit
   print*, 'enter number of iterations'
   read*, ni


   munit=100
   eunit=200
   iunit=300
   avmunit=400
   aveunit=500
   abmunit=600
   cvunit=700
   chiunit=800
   bcunit=900
   !Automated over L=6,12

   do L=6,12
   N=L*L*L
   print*,"Value of L",L

   allocate(spin(L,L,L))
   E=0.0d0
   M=0.0d0

   write(filename,'("initial_lattice_L",I2.2,".dat")')L
   open(iunit,file=filename)
   !lattice initialization

   do i=1,L
     do j=1,L
        do k=1,L
           call random_number(r)


           if (r<0.5d0)then
                spin(i,j,k)=-1
           else
                spin(i,j,k)=1
           end if
           write(iunit,*) i,j,k,(spin(j,i,k))
        end do
     end do
   end do
   close(iunit)
! Calculating the Initial Energy and Magnetization
   do i=1,L
     do j=1,L
       do k=1,L
          a=i+1 ; b=i-1 ;c=j+1 ;d=j-1 ;f=k+1; g=k-1
          !Periodic Boundary condition
          if(i==L) a=1
          if(i==1) b=L
          if(j==1) d=L
          if(j==L) c=1
          if(k==L) f=1
          if(k==1) g=L

       M=M+spin(i,j,k)
       E=E-J_ising*dfloat(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,f)+spin(i,j,g)))
       end do
     end do
   end do

   Mag=M/(dfloat(N))
   E=E*0.5d0

   print*, 'initial energy E, E per spin=', E, E/dfloat(N)
   print*, 'initial magnetization M, M per spin=', M, Mag

   write(filename,'("av_mag_L",I2.2,".dat")')L
   open(avmunit,file=filename)
   write(filename,'("av_en_L",I2.2,".dat")')L
   open(aveunit,file=filename)
   write(filename,'("abs_mag_L",I2.2,".dat")')L
   open(abmunit,file=filename)
   write(filename,'("cv_L",I2.2,".dat")')L
   open(cvunit,file=filename)
   write(filename,'("chi_L",I2.2,".dat")')L
   open(chiunit,file=filename)

   write(filename,'("bc_L",I2.2,".dat")')L
   open(bcunit,file=filename)
   !Temperature loop
   do T_temp=56,25,-1
        T=dfloat(T_temp)/10.0d0
        print*,"temperaturei:",T


   av_m=0.0d0;av_e=0.0d0;av_m_n=0.0d0;absm=0.0d0; absm_m=0.0d0
   av_e_n=0.0d0; av_m2=0.0d0; av_e2=0.0d0;av_m4=0.0d0; binc=0.0d0

   write(filename,'("mag_T",F4.2,"_L",I2.2,".dat")') T,L
   open(munit,file=filename)
   write(filename,'("energy_T",F4.2,"_L",I2.2,".dat")') T,L
   open(eunit,file=filename)
!Monte carlo iteration
   do tm=1,ni
         do mm=1,N
                  call random_number(r) ; i=int(r*dfloat(L))+1
                  call random_number(r) ; j=int(r*dfloat(L))+1
                  call random_number(r) ; k=int(r*dfloat(L))+1

                   a=i+1 ; b=i-1 ; c=j+1 ; d=j-1 ; f=k+1 ; g=k-1
                   if(i==L) a=1 ; if(i==1) b=L ; if(j==L) c=1 ; if(j==1) d=L ; if(k==1) g=L; if(k==L) f=1

                   Ei= -J_ising*dfloat(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,f)+spin(i,j,g)))
                   ! Trial flip
                   spin(i,j,k)=-spin(i,j,k) 

                   Ef= -J_ising*dfloat(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,f)+spin(i,j,g)))


                   dE=Ef-Ei
                   if(dE<0.0d0)then
                      E=E+dE
                      M=M+(2.0d0*dfloat(spin(i,j,k)))
                   else
                      u=exp(-dE/(T))
                      call random_number(h)
                      if(h<u) then
                          E=E+dE
                          M=M+(2.0d0*dfloat(spin(i,j,k)))
                      else
                          spin(i,j,k)=-spin(i,j,k)
                      end if
                   end if
         end do
         mag=M/(dfloat(N))
         if(tm .gt. n_equil) then
            !if(mod(tm,n_stat) .eq. 0) then
                absm=abs(M)/dfloat(N)

                av_m=av_m+mag;
                absm_m=absm_m+absm; av_e=av_e+E/dfloat(N)
                av_m_n=av_m_n+M; av_e_n=av_e_n+E
                av_m2=av_m2+M*M; av_e2=av_e2+(E*E)
                av_m4=av_m4+M*M*M*M
            !end if
         end if
       write(munit,*)tm,mag
       write(eunit,*)tm,E/dfloat(N)
    end do !monte carlo iteration loop end

!Calucating average quantities
    av_m_n=av_m_n/dfloat(ni-n_equil)
    av_e_n=av_e_n/dfloat(ni-n_equil)
    av_m2=av_m2/dfloat(ni-n_equil)
    av_e2=av_e2/dfloat(ni-n_equil)
    cv=(av_e2-av_e_n*av_e_n)/(T*T)
    chi=(av_m2-av_m_n*av_m_n)/T
    binc=binc+(1-(av_m4/(3*av_m2*av_m2)))


    write(avmunit,*)T,av_m/dfloat(ni-n_equil)
    write(aveunit,*)T,av_e/dfloat(ni-n_equil)
    write(abmunit,*)T,absm_m/dfloat(ni-n_equil) 
    write(chiunit,*)T,chi
    write(cvunit,*)T,cv
    write(bcunit,*)T,binc


    close(munit)
    close(eunit)
    end do !Temperature loop end

    deallocate(spin)
    close(avmunit)
    close(aveunit)
    close(abmunit)
    close(chiunit)
    close(cvunit)
    close(bcunit)
    end do !Lattice size loop end
end program Ising3D
