program Gaussian
      implicit none

      real*8:: r,gauss
      integer, allocatable,dimension(:) ::bin
      integer:: n,i,nbin

      write(*,*) "Enter the value of n"
      read(*,*) n
      write(*,*) "Enter the number of bin"
      read(*,*) nbin

      open(21,file="gaussian_histogram.dat")
      open(22,file="gaussian_data.dat")

      allocate(bin(-nbin:nbin))


      bin=0
      do i=1,n
        r=gauss(0.0d0,1.0d0)

        if (abs(int(r))<=nbin) then
                bin(int(r))=bin(int(r))+1
                write(22,*) r

        end if
      end do
      write(*,*) bin

      do i= -nbin,nbin
        write(21,*) i, bin(i)
      end do
      call system('gnuplot gaussian_histogram.plt')
      call system('feh gaussian_histogram.png')

end program Gaussian

        
       
function gauss(mu,sigma) result(y)
        ! This function generate a random number from a uniform distribution [-1,1]
        ! Uses Box Mullar polar form
        real*8:: mu,sigma,y,u,v
3       call random_number(u)
        call random_number(v)
        u=2*u-1
        v=2*v-1
        y=u*u+v*v
        if (y==0.d0 .or. y>=1.0d0) goto 3
        y=sqrt(-2*log(y)/y)
        y=u*y
        y=sigma*y+mu
end function gauss
        

        
