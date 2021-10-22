program gauss
    implicit none

    integer i, j, k, n
    double precision a(100, 100)
    double precision b(100), x(100)
    double precision tmp
    double precision amax
    integer ip
    double precision:: eps = 1.0d-14
    
    open(10, file = 'input.dat', action = 'read')
    read(10, *) n
    do i = 1, n
        read(10, *) (a(i, j), j = 1, n)
        write(*, *) (a(i, j), j = 1, n)
    end do
    
    do i = 1, n
    read(10, *) b(i)
    write(*, *) b(i)
    end do
    
    close(10)
    
    do k = 1, n-1
        amax = abs(a(k, k))
        ip = k
        do i = k+1, n
            if (abs(a(i, k)) > amax) then
                amax = abs(a(i, k))
                ip = i
            end if
        end do
        if (amax < eps) stop
        if (k /= ip) then
            do j = k, n
                tmp = a(k, j)
                a(k, j) = a(ip, j)
                a(ip, j) = tmp
            end do
            tmp = b(k)
            b(k) = b(ip)
            b(ip) = tmp
        end if
        
        do i = k+1, n
            tmp = a(i, k) / a(k, k)
            b(i) = b(i) - tmp * b(k)
            do j = k+1, n
                a(i, j) = a(i, j) - tmp * a(k, j)
            end do
        end do
    end do

    x(n) = b(n) / a(n, n)
    do k = n-1, 1, -1
        x(k) = b(k)
        do j = k+1, n
            x(k) = x(k) - a(k, j) * x(j)
        end do
        x(k) = x(k) / a(k, k)
    end do
    
    open(20, file = 'output.dat', action = 'write')
    
    write(*, *) "Ax=bの解xは以下の通りです．"
    write(20, *) "Ax=bの解xは以下の通りです．"
    
    do i = 1, n
        write(*, '(F12.4)') x(i)
        write(20, '(F12.4)') x(i)
    end do
    
    close(20)
    
end program gauss