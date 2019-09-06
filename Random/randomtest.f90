program randomtest
        implicit none

        real*8:: sums,sq_sums,mean,sq_mean,sigma,c_k,sum2
        integer:: n,i,j,k,m
        real*8, allocatable, dimension(:):: r
        character(len=32):: filename

        open(21,file="meandata.dat")
        open(22,file="sd.dat")


        do j=2,8 
                n=10**j
        
                allocate(r(n))

                call random_number(r)
                ! Calculating the mean and standard deviation
                sums=0.0d0
                sq_sums=0.0d0
                do i=1, n
                        sums=sums+r(i)
                        sq_sums=sq_sums+r(i)**2
                end do

                mean=sums/real(n)
                sq_mean=sq_sums/real(n)
                sigma=sqrt(sq_mean-mean**2)
                write(*,*) "number of random number:",n
                write(*,*) "The mean is:", mean
                write(*,*) "The standard deviation is:", sigma
                write(21,*) 1.0d0/sqrt(real(n)),abs(mean-0.5d0)
                write(22,*) n,sigma

                
                if (n==10000) then
                        open(23,file="scatterdata.dat")
                ! Scatter plot
                        do i=1,n-1
                                write(23,*) r(i),r(i+1)
                        end do
                        
                       
                end if
                !Correlation plot
                if (n==10000) then
                        open(24,file="correlationdata.dat")
                        do m=1,n/2
                                k=m-1
                                sum2=0.0d0
                                do i=1,n-k
                                        sum2=sum2+r(i)*r(i+k)
                                end do
                                sum2=sum2/real(n-k)
                                c_k=(sum2-mean**2)/sigma**2

                                write(24,*) k,c_k
                       end do
                end if
                if (n==10000) then

                        open(25,file="correlationdata_small.dat")
                        do m=1,50
                                k=m-1
                                sum2=0.0d0
                                do i=1,n-k
                                        sum2=sum2+r(i)*r(i+k)
                                end do
                                sum2=sum2/real(n-k)
                                c_k=(sum2-mean**2)/sigma**2

                                write(25,*) k,c_k
                       end do
                end if


                                        
                                
                deallocate(r)
       end do
                
end program randomtest

