program newton
    implicit none
    
    integer i, n
    double precision x, d
    double precision :: eps = 1.0d-8
    double precision f, df
    
    n = 1000
    
    write(*,*)"初期値x0を入力してください．"
    read(*,*) x
    
    do i = 1, n
        d = - f(x) / df(x)
        x = x + d
        write(*, *) x
        if( abs(d) < eps ) exit
    end do
    
    if( i < n ) then
        write(*, '(A12,F12.8)') "求める解は", x
    else
        write(*, *) "解が収束しませんでした．"
        write(*, *) "なおxの値は", x
        
    end if
   
end program newton

function f(x)
    implicit none
    double precision, intent(in) :: x
    double precision f
    
    f = x ** 3 - 6 * x ** 2 + 7 * x + 2.0d0

end function f

function df(x)
    implicit none
    double precision, intent(in) :: x
    double precision df
    
    df = 3 * x ** 2 - 12 * x + 7.0d0

end function df
