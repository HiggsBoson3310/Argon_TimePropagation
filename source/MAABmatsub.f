      subroutine eig_r(n, M, eigr, eigi, vecr)
          implicit none
          integer, parameter :: LWMAX = 1000
          integer :: n
          real*8 :: M(n,n), eigr(n), eigi(n)
          real*8 :: vecr(n,n), vecl(n,n)
          integer :: LDA, INFO, LWORK
          real*8 :: work(LWMAX)
          LDA = n
          LWORk = -1
          work = 0.0d0
          call DGEEV('N','V', n, M, LDA, eigr, eigi, vecl, 
     x      n, vecr, n, work, LWORK, INFO)
          LWORK = min(LWMAX,int(work(1)))
          call DGEEV('N','V', n, M, LDA, eigr, eigi, vecl,
     x     n, vecr, n, work, LWORK, INFO)
          return
      end subroutine     
      
      subroutine eig_rs(n, M, eig)
          implicit none
          integer, parameter :: LWMAX = 1000
          integer :: n
          real*8 :: M(n,n), eig(n)
          integer :: LDA, INFO, LWORK
          real*8 :: work(LWMAX)
          LDA = n
          LWORk = -1
          call DSYEV('V','U', n, M, LDA, eig, work, LWORK, INFO)
          LWORK = min(LWMAX,int(work(1)))
          call DSYEV('V','U', n, M, LDA, eig, work, LWORK, INFO)
          return
      end subroutine

      subroutine find_r_rough(funk, xl, xg, dx, nr, root)
      implicit none
      logical :: concave, convex
      real*8, external :: funk
      real*8 :: xl, xg, dx, fl, fg, ff
      integer :: N, nr, i
      real*8 :: root(700), linroot
      N = int(abs(xl-xg)/dx)
      nr = 0
      DO i=1,N
        ff = funk(xl+(i+1)*dx)
        fl = funk(xl+i*dx)
        fg = funk(xl+(i-1)*dx)
        concave = (ff>fl).and.(fl<fg)
        convex = (ff<fl).and.(fl>fg)
        if((fl*fg).lt.0) then
            nr = nr+1
            root(nr) = xl+i*dx-0.5d0*dx
        endif 
        if(nr.gt.700) then
            write(6,*) "More than 700 resonances."
            return 
        endif
      END DO
      return 
      end subroutine

      function linroot(xl,xh,yl,yh)
        real*8 :: xl,xh,yl,yh
        real*8 :: m,b, linroot
        m = (yl-yh)/(xl-xh)
        b = 0.5d0*((yl+yh)-m*(xl+xh))
        linroot = -b/m
      end function

