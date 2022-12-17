module matmod
    implicit none
    contains
    !
    subroutine zdiag(n, A, eigs, U)
        integer LWORK, info
        integer :: n, aux
        complex*16 :: A(n,n), U(n,n)
        real*8 :: eigs(n)
        complex*16 :: Acopy(n,n), dummy(10)
        real*8, allocatable :: rwork(:)
        complex*16, allocatable :: work(:)
        external :: zheev
        write(6,*) "The almighty zdiag has been called"
        Acopy = A
        aux = max(1, 3*n-1)
        allocate(rwork(aux))
        write(6,*) n
        LWORK=-1
        call zheev('V', 'U', n, Acopy, n, eigs, dummy, LWORK, rwork, info)
        
        aux = int(real(dummy(1)))
        allocate(work(aux))
        LWORK = max(1, 2*N-1)+1
        call zheev('V', 'U', n, Acopy, n, eigs, work, LWORK, rwork, info)
        U = Acopy
        deallocate(rwork)
        deallocate(work)
        return
    end subroutine zdiag
    !
    subroutine zdiag_save(n, A, eigs, U)
        integer LWORK, info
        integer :: n, aux
        complex*16 :: A(n,n), U(n,n)
        real*8 :: eigs(n)
        complex*16 :: Acopy(n,n), dummy(10)
        real*8, allocatable :: rwork(:)
        complex*16, allocatable :: work(:)
        external :: zheev
        write(6,*) "The almighty zdiag has been called"
        Acopy = A
        aux = max(1, 3*n-1)
        allocate(rwork(aux))
        write(6,*) n
        LWORK=-1
        call zheev('V', 'U', n, Acopy, n, eigs, dummy, LWORK, rwork, info)
        
        aux = int(real(dummy(1)))
        allocate(work(aux))
        LWORK = max(1, 2*N-1)+1
        call zheev('V', 'U', n, Acopy, n, eigs, work, LWORK, rwork, info)
        U = Acopy
        
        open(412,file='DiagH.db',STATUS='REPLACE', &
        ACCESS='STREAM',FORM='UNFORMATTED')
        open(512,file='RUmat.db',STATUS='REPLACE', &
        ACCESS='STREAM',FORM='UNFORMATTED')
        open(612,file='IUmat.db',STATUS='REPLACE', &
        ACCESS='STREAM',FORM='UNFORMATTED')
        open(712,file='HMatR.db',STATUS='REPLACE', &
        ACCESS='STREAM',FORM='UNFORMATTED')
        open(812,file='HMatI.db',STATUS='REPLACE', &
        ACCESS='STREAM',FORM='UNFORMATTED')

        write(412) eigs
        write(512) REAL(U)
        write(612) AIMAG(U)
        write(712) REAL(A)
        write(812) AIMAG(A)

        close(412)
        close(512)
        close(612)
        CLOSE(712)
        CLOSE(812)

        deallocate(rwork)
        deallocate(work)
        return
    end subroutine zdiag_save
    !
    subroutine zmatvec(n,m, A, X, Y)
        ! A is n by m matrix
        ! B is m vector
        ! C is n vector
        complex*16, parameter :: Uno=(1.0d0,0.d0), Cer=(0.0d0,0.0d0)
        integer, intent(IN) :: n, m
        complex*16, intent(IN) :: A(n,m), X(m)
        complex*16, intent(OUT) :: Y(n)
        external :: zgemm
        Y = Cer
        call zgemv('N',n,m,Uno,A,n,X,1,Cer,Y,1)
        return
    end subroutine
    !
    subroutine zmatTvec(n,m, A, X, Y)
        ! A is m by n matrix
        ! B is m vector
        ! C is n vector
        complex*16, parameter :: Uno=(1.0d0,0.d0), Cer=(0.0d0,0.0d0)
        integer, intent(IN) :: n, m
        complex*16, intent(IN) :: A(m,n), X(m)
        complex*16, intent(OUT) :: Y(n)
        external :: zgemm
        Y = Cer
        call zgemv('T',n,m,Uno,A,n,X,1,Cer,Y,1)
        return
    end subroutine
    !
    subroutine zmatHvec(n,m, A, X, Y)
        ! A is m by n matrix
        ! B is m vector
        ! C is n vector
        complex*16, parameter :: Uno=(1.0d0,0.d0), Cer=(0.0d0,0.0d0)
        integer, intent(IN) :: n, m
        complex*16, intent(IN) :: A(m,n), X(m)
        complex*16, intent(OUT) :: Y(n)
        external :: zgemm
        Y = Cer
        call zgemv('C',n,m,Uno,A,n,X,1,Cer,Y,1)
        return
    end subroutine
    !
    subroutine zInvMat(n,A,Ai)
        integer :: n, ipiv(n), info, lwork
        complex*16, allocatable :: work(:)
        complex*16 :: A(n,n), Ai(n,n)
        external :: zgetrf, zgetri
            Ai = A
            call zgetrf(n,n,Ai,n,ipiv,info)
            allocate(work(1))
            call zgetri(n,Ai,n,ipiv,work,-1,info)
            lwork = int(work(1))
            deallocate(work)
            allocate(work(lwork))
            call zgetri(n,Ai,n,ipiv,work,lwork,info)
    end subroutine ZInvMat
    !
end module matmod
