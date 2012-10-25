module globals
	implicit none 
	integer, parameter :: wp = selected_real_kind(p=14)
end module globals


Program RootFinding
	use globals
	implicit none
	integer :: i 
	real(wp) :: x

!	do i = -628, 628,1
!		x = i/100._wp
!		write(*,*) x,tester(x)
!	end do 

	write (*,*) NewtonsMethod(tester,6._wp)

	contains
	
	function diffy(dyn_func, c, h,accuracy)
		use globals
		implicit none
	  	real(wp) :: c,h
	  	real(wp), external :: dyn_func
	 	real(wp) :: diffy
		integer :: accuracy
	  !
	  	if(accuracy.eq.2) then
			diffy = (dyn_func(c+h)-dyn_func(c-h))/(2._wp*h)
		else if(accuracy.eq.4) then
			diffy = (-dyn_func(c+2._wp*h)+8._wp*dyn_func(c+h) &
			&       -8._wp*dyn_func(c-h)+dyn_func(c-2._wp*h))/(12._wp*h)
		else if(accuracy.eq.6) then
			diffy = (-1._wp/60._wp*dyn_func(c-3._wp*h)+3._wp/20._wp*dyn_func(c-2._wp*h)-3._wp/4._wp*dyn_func(c-h)&
			&			+1._wp/60._wp*dyn_func(c+3._wp*h)-3._wp/20._wp*dyn_func(c+2._wp*h)+3._wp/4._wp*dyn_func(c+h))/h
		else if(accuracy.eq.8) then
			diffy = (1._wp/280._wp*dyn_func(c-4._wp*h)-4._wp/105._wp*dyn_func(c-3._wp*h)&
			&			+1._wp/5._wp*dyn_func(c-2._wp*h)-4._wp/5._wp*dyn_func(c-h)&
			&			-1._wp/280._wp*dyn_func(c+4._wp*h)+4._wp/105._wp*dyn_func(c+3._wp*h)&
			&			-1._wp/5._wp*dyn_func(c+2._wp*h)+4._wp/5._wp*dyn_func(c+h))/h
		end if
	  !
	end function diffy
	
	function tester(x)
		use globals
		implicit none	
		real(wp) :: x,tester
		!
			tester = Sin(x)*x**2-1._wp
		!
	end function tester

	function NewtonsMethod(Zfunc,x0) result (x1)
		use globals
		implicit none
		real(wp), external :: Zfunc
		real(wp) :: x0, x0i, x1, ZfuncValAbs
		integer :: itts = 0 
		
		! Accuracy of the zero.
		real(wp) :: eps = 0.000000001_wp

		! Settings of the numerical differentiation function.
		real(wp) :: h = 0.0000001_wp
		integer :: accuracy = 8

		x0i = x0
		x1 = x0i - Zfunc(x0i)/diffy(Zfunc,x0i,h,accuracy)		
		x0i = x1

		ZfuncValAbs = abs(Zfunc(x1))
		do while((ZfuncValAbs.gt.eps).and.(itts.lt.10))
			x1 = x0i - Zfunc(x0i)/diffy(Zfunc,x0i,h,accuracy)
			x0i = x1
			itts = itts + 1
			ZfuncValAbs = abs(Zfunc(x1))
			write(0,*) itts, Zfunc(x1)
		end do

	end function NewtonsMethod



end Program RootFinding