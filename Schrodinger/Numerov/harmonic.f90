Program Harmonic
      implicit none




      print '(s,$)', 'Max value for x: '
      read (*,*),xmax
      print '(s,$)', 'Number of grid points: '
      read (*,*),mesh

      allocate ( x(0:mesh),y(0:mesh),p(0:mesh)

      dx=xmax/mesh
      ddx12=dx*dx/12.0_dp

      ! set up the potential

      do i=0,mesh
          x(i)=float(i)*dx
          vpot(i)=0.5_dp*x(i)*x(i)
      end do

      print '(a,$)','output file name= '
      read (*,'(a)') fileout
      open (21,file=fileout,status='unknown',form='formatted')


      ! Eigenvalue Search

      do

      ! Set initial lower and upper bounds to the eigenvalue

      eup=maxval(vpot(:))
      elw=minval(vpot(:))

      ! Set trial energy

      print '(a,$)', 'Trial energy ('

      e= 0.5_dp*(elw+eup)
      n_iter=1000


      !Iteration
      do=kkk=1,n_inter

      ! set up the f-function
      ! f<0 classically allowed
      ! f>0 classically not allowd

      f(0)=ddx12*(2.0_dp*(vpot(0)-e))
      icl=-1

      do i=1,mesh
      f(i)=ddx12*2.0_dp(vpot(i)-e)
      if (f(i)==0.0_dp) f(i)=1.d-20

      ! icl stores where the last change of sign observed
      if (f(i)/=sign(f(i),f(i-1))) icl=i
      end do


