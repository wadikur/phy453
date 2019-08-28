program trapezoidal
implicit none
real*8::h,a,b,area,sums,f,results,x,error
integer::n,i
a=0.0
b=1.0
n=1000
h=(b-a)/n
area=(f(a)+f(b))/2

error=100.0

do while(abs(error)>1e-6)
n=n+1000
h=(b-a)/n
sums=0

do i=0,n-1
x=a+h*i
sums=sums+f(x)
end do

results=h*(area+sums)
error=4*atan(1.0d0)-results

write(*,*)n,"integral value",results,error
end do

end program trapezoidal

function f(x) result(y)
real*8::x,y
y=4*(1/(1+x**2))
end function f



