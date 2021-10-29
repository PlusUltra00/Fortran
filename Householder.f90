program householder
	implicit none
	integer , parameter::nd = 3
	double precision alp(nd), bet(nd-1), eigen(nd)
	integer i, ind

	do i = 1, 3
		alp(i) = 2
	end do
	do i = 1, 2
		bet(i) = -1
	end do

	call bisection(alp, bet, nd, eigen, ind)
	do i = 1, nd
		write(*, *) sqrt(eigen(i))
	end do

end program householder

subroutine bisection(alp, bet, nd, eigen, ind)
	implicit none
	integer, intent(in)::nd
	integer, intent(out)::ind
	double precision, intent(in)::alp(nd), bet(nd-1)
	double precision, intent(out)::eigen(nd)
	integer i, k, it, itm, nc, ip
	double precision r1, r2, a, b, c, g
	integer, parameter::itmax = 1000
	double precision, parameter::eps = 1.0d-8


	r2 = abs(alp(1)) + abs(bet(1))
	do i = 2, nd-1
		r2 = max(r2, abs(bet(i-1)) + abs(alp(i)) + abs(bet(i)))
	end do
	r2 = max(r2, abs(bet(nd-1)) + abs(alp(nd)))
	r1 = -r2
	itm = 0

	do k = 1, nd
		a = r1
		b = r2
		do it = 1, itmax
			c = (a + b) / 2.0d0
			g = c - alp(1)
			nc = 0
			ip = 0
			if(g <= 0.0d0) then
				nc = 1
			end if
			do i = 2, nd
				if(ip == 1) then
					g = c - alp(i)
				else if(abs(g) < 0.10d-12) then
					ip = 1
					cycle
				else
					g = c - alp(i) - bet(i - 1) ** 2.0d0 / g
				end if
				if(g <= 0.0d0) then
					nc = nc + 1
				end if
			end do
			if(nc >= nd .and. r1 < a) then
				r1 = a
			end if
			if(nc >= k) then
				a = c
			else
				b = c
			end if
			if (b - a < eps) exit
		end do
		itm = max(it, itm)
		eigen(k) = c
		r2 = c
	end do
	if (itm > itmax) then
		ind = -itmax
	else
		ind = itm
	end if

end subroutine bisection