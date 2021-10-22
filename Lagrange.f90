program ex7_3
    implicit none
    
    integer k
    integer, parameter :: n = 10
    double precision xp(1:n), yp(1:n)
    integer, parameter :: out = 100
    double precision x, y
    
    open(10, file = 'input7_3.txt', action = 'read')
    do k = 1, n
        read(10, *) xp(k), yp(k)
    end do
    close(10)
    
    open(20, file = 'output7_3.csv', status = 'replace')
    do k = 0, out
        x = xp(1) + k * (xp(n) - xp(1)) / dble(out)
        y = lagrange(x)
        write(20, *) x, ',', y
    end do
    close(20)
    
contains

function lagrange(xin)
    implicit none
    
    integer in, jn
    double precision, intent(in) :: xin
    double precision lagrange
    double precision L
    
    lagrange = 0.0d0
    
    do in = 1, n
        L = 1.0d0
    
        do jn = 1, n
            if(in /= jn) then
                L = L * (xin - xp(jn)) / (xp(in) - xp(jn))
            end if
        end do
    
        lagrange = lagrange + L * yp(in)
    end do

end function

end program