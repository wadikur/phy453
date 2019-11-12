program Neumerov
implicit none
!integer, parameter :: dp = selected_real_kind(14,200)
integer :: mesh, i, icl
integer :: nodes, hnodes, ncross, j, n_iter
real*8 :: xmax, dx, ddx12, norm, arg, djump, fac
real*8 :: eup, elw, e
real*8, allocatable :: x(:), y(:), p(:), vpot(:), f(:)
character (len=80) :: fileout

print '(a,$)','Max value for x? '
read (*,*) xmax
print '(a,$)', 'Number of grid points? '
read (*,*) mesh

allocate ( x(0:mesh), y(0:mesh), p(0:mesh), vpot(0:mesh), f(0:mesh) )
dx =  xmax/mesh 
ddx12=dx*dx/12.0d0

do i = 0, mesh
   x(i) = float(i) * dx
   vpot(i) = 0.5d0 * x(i)*x(i)
end do

print '(a,$)','Output file name? '
read (*,'(a)') fileout

if ( fileout /= ' ' ) &
   open (7,file=fileout,status='unknown',form='formatted')


do 
  print '(a,$)', 'nodes (type -1 to stop) > '
  read (*,*) nodes
  if (nodes < 0) then
     close(7)
     deallocate ( f, vpot, p, y, x )
     stop 
  end if

  eup=maxval (vpot(:))
  elw=minval (vpot(:))

  print '(a,$)','Trial energy (0=search with bisection) '
  read (*,*) e
  if ( e == 0.0d0 ) then

     e = 0.5d0 * (elw + eup)
     n_iter = 1000
  else
     n_iter = 1
  endif
  
  do j = 1, n_iter
     f(0)=ddx12*(2.0d0*(vpot(0)-e))
     icl=-1
     do i=1,mesh
        f(i)=ddx12*2.0d0*(vpot(i)-e)
        if ( f(i) == 0.0d0) f(i)=1.d-20
        if ( f(i) /= sign(f(i),f(i-1)) ) icl=i
     end do
  
     if (icl >= mesh-2) then
        deallocate ( f, vpot, p, y, x )
        stop 'last change of sign too far'
     else if (icl < 1) then
        deallocate ( f, vpot, p, y, x )
        stop 'no classical turning point?'
     end if

     f = 1.0d0 - f
     y = 0.0d0

     hnodes = nodes/2

     if (2*hnodes == nodes) then
        y(0) = 1.0d0
        y(1) = 0.5d0*(12.0d0-10.0d0*f(0))*y(0)/f(1)
     else
        y(0) = 0.0d0
        y(1) = dx
     end if

     ncross=0
     do i =1,icl-1
        y(i+1)=((12.0d0-10.0d0*f(i))*y(i)-f(i-1)*y(i-1))/f(i+1)
        if ( y(i) /= sign(y(i),y(i+1)) ) ncross=ncross+1
     end do
     fac = y(icl)

     if (2*hnodes == nodes) then
        ncross = 2*ncross
     else
        ncross = 2*ncross+1
     end if

     if ( n_iter > 1 ) then
        if (ncross /= nodes) then
           if ( j == 1) &
               print '("Bisection         Energy       Nodes  Discontinuity")'
           print '(i5,f25.15,i5)', j, e, ncross
           if (ncross > nodes) then
              eup = e
           else 
              elw = e
           end if
           e = 0.5d0 * (eup+elw)
           cycle
        end if
     else
        print *, e, ncross, nodes
     end if

     y(mesh) = dx
     y(mesh-1) = (12.0d0-10.0d0*f(mesh))*y(mesh)/f(mesh-1)
     do i = mesh-1,icl+1,-1
        y(i-1)=((12.0d0-10.0d0*f(i))*y(i)-f(i+1)*y(i+1))/f(i-1)
     end do

     fac = fac/y(icl)
     y(icl:) = y(icl:)*fac

     norm = (2.0d0*dot_product (y, y) - y(0)*y(0)) * dx 
     y = y / sqrt(norm)
     !
     if ( n_iter > 1 ) then
        djump = ( y(icl+1) + y(icl-1) - (14.0d0-12.0d0*f(icl))*y(icl) ) / dx
        print '(i5,f25.15,i5,f14.8)', j, e, nodes, djump
        if (djump*y(icl) > 0.0d0) then
           eup = e
        else
           elw = e
        endif
        e = 0.5d0 * (eup+elw)
        if ( eup-elw < 1.d-10) exit 
     endif
  end do 

  norm = 0.0d0
  p(icl:) = 0.0d0
  do i=0,icl
     arg = (e - x(i)**2/2.0d0)
     if ( arg > 0.0d0) then
        p(i) = 1.0d0/sqrt(arg)
     else
        p(i) = 0.0d0
     end if
     norm = norm + 2.0d0*dx*p(i)
  enddo
  norm = norm - dx*p(0)
  p(:icl-1) = p(:icl-1)/norm
  write (7,'("#   x       y(x)         y(x)^2         classical p(x)      V")')
  do i=mesh,1,-1
     write (7,'(f7.3,3e16.8,f12.6)') & 
          -x(i), (-1)**nodes*y(i), y(i)*y(i), p(i), vpot(i)
  enddo
  do i=0,mesh
     write (7,'(f7.3,3e16.8,f12.6)') &
          x(i), y(i), y(i)*y(i), p(i), vpot(i)
  enddo
  write (7,'(/)')

end do 

end program Neumerov
