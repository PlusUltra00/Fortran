program jacobi
    implicit none
    
    integer i, j, k, n
    double precision, allocatable :: a(:, :)
    double precision b(100), x(100)
    double precision x_1(100)
    double precision dr
    double precision :: eps = 1.0d-8
    integer :: kmax = 1000
    allocate(a(100, 100))
    
    open(10, file = 'input.dat', action = 'read')
    read(10, *) n
    
    do i = 1, n
        read(10, *) (a(i, j), j = 1, n)
        write(*, *) (a(i, j), j = 1, n)
    end do
    
    do i = 1, n
        read(10, *) b(i)
        write(*, *) b(i)
        x_1(i) = 1.0d0
    end do
    
    close(10)
    
    do k = 1, kmax
        do i = 1, n
            x(i) = b(i)
            do j = 1, n
                if(i /= j) then
                x(i) = x(i) - a(i, j) * x_1(j)
                end if
            end do
            x(i) = x(i) / a(i, i)
        end do
        dr = 0.0d0
        do i = 1, n
            dr = dr + (x(i) - x_1(i)) ** 2
            x_1(i) = x(i)
        end do
        dr = sqrt(dr)
        if (dr < eps) exit
    end do
    
    open(20, file = 'output.dat', action = 'write')
    
    if(k >= kmax) then
        write(*, *) "解が収束しませんでした"
        write(20, *) "解が収束しませんでした"
    else
        write(*, *) "Ax=bの解xは以下の通りです．"
        write(20, *) "Ax=bの解xは以下の通りです．"
        do i = 1, n
            write(*, '(F12.4)') x(i)
            write(20, '(F12.4)') x(i)
        end do
        write(*, *) "なお反復回数k=", k
        write(20, *) "なお反復回数k=", k
    end if
    
    close(20)

end program jacobi