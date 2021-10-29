program integral
	implicit none
	
	integer n
	double precision, parameter:: a = -1.0d0
	double precision, parameter:: b = 1.0d0
	double precision s_kubun, s_daikei, s_simpson
	
	n = 10
	
	call kubun(n, a, b, s_kubun)
	call daikei(n, a, b, s_daikei)
	call simpson(n, a, b, s_simpson)
	
	write(*, *) 'n = 10', s_kubun, s_daikei, s_simpson
	
	n = 20
	
	call kubun(n, a, b, s_kubun)
	call daikei(n, a, b, s_daikei)
	call simpson(n, a, b, s_simpson)
	
	write(*, *) 'n = 20', s_kubun, s_daikei, s_simpson
	
	n = 40
	
	call kubun(n, a, b, s_kubun)
	call daikei(n, a, b, s_daikei)
	call simpson(n, a, b, s_simpson)
	
	write(*, *) 'n = 40', s_kubun, s_daikei, s_simpson
   
end program integral

function f(xin)
	implicit none
	double precision, intent(in):: xin
	double precision f
	
	f = 2 / (1 + xin ** 2)
	
end function f

subroutine kubun(n, a, b, s_kubun)
	implicit none
	integer, intent(in):: n
	double precision, intent(in):: a, b
	double precision, intent(out):: s_kubun
	double precision h
	double precision f
	integer i
	
	h = (b - a) / dble(n)
	s_kubun = 0.0d0
	
	do i = 0, n-1
		s_kubun = s_kubun + f(a + i * h) * h
	end do
   
end subroutine kubun

subroutine daikei(n, a, b, s_daikei)
	implicit none
	integer,intent(in):: n
	double precision, intent(in):: a, b
	double precision, intent(out):: s_daikei
	double precision h
	double precision f
	integer i
	double precision f1, f2
	
	h = (b - a) / dble(n)
	s_daikei = 0.0d0
	
	do i = 0, n-1
		f1 = f(a + i * h)
		f2 = f(a + (i + 1.0d0) * h)
		s_daikei = s_daikei + (f1 + f2) * h / 2.0d0
	end do

end subroutine daikei

subroutine simpson(n, a, b, s_simpson)
	implicit none
	integer, intent(in):: n
	double precision, intent(in):: a, b
	double precision, intent(out):: s_simpson
	double precision h
	double precision f
	integer i
	double precision f1, f2, f3
	
	h = (b - a) / dble(n)
	s_simpson = 0.0d0
	
	do i = 0, n-1
		f1 = f(a + i * h)
		f2 = f(a + (i + 0.5d0) * h)
		f3 = f(a + (i + 1.0d0) * h)
		s_simpson = s_simpson + (f1 + 4 * f2 + f3) * h / 6.0d0
	end do

end subroutine simpson