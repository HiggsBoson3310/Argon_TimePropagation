program testing
    use matmod
    implicit none
    complex*16, parameter :: zone = (1.0d0,0.0d0), zima = (0.0d0,1.0d0)
    complex*16 :: xmat(2,2)
    real*8 :: tst(2,2)
    integer :: n
    real*8 :: eigs(2)
    complex*16 :: U(2,2)
    n=2
    xmat = reshape((/z(1.0d0,0.0d0),z(2.0d0,-1.0d0),z(2.0d0,1.0d0),z(2.0d0,0.0d0)/), shape(xmat))
    write(6,*) n
    call zdiag(n, xmat, eigs, U)
    tst = 0.d0
    write(6,*) xmat
    !write(6,*) U
    !write(6,*) eigs
    tst = reshape((/eigs(1), 0.d0, 0.d0, eigs(2)/),shape(tst))
    write(6,*) matmul(U,matmul(tst, TRANSPOSE(CONJG(U))))
    contains 
    
    function Z(x,y)
        real*8 :: x,y
        complex*16 :: z
        z = x*zone+zima*y
    end function
end program