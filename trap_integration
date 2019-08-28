program Trapezoidal
	implicit none
	real*8 :: xi, xf,width,x(100000),fx(100000),Area,error
	integer :: i,n 
	

	xi=0.0
	xf=1.0
	n=100
	
	error=100.0
	do while (abs(error)>1e-10)
		n=n+100
		width=(xf-xi)/n
		Area=0.0
	
		do i=1, n+1
			x(i)=xi+(i-1)*width
			fx(i)=4*(1.0d0/(1+x(i)**2))
		end do
		
		do i=1,n
			Area=Area+((fx(i)+fx(i+1))/2)*width
	
		end do

		error=4.0d0*atan(1.0d0)-Area
		print*, n
		print*,Area
	
		print*,(4.0d0*atan(1.0d0))
	
		print*,error
	end do
end program Trapezoidal
