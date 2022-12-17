module j_recoup
    implicit none
    contains
    subroutine bJcK_kLS(n, lc, sc, jc, le, K, se, L, S, J, UFT)
         
        integer, intent(in) :: n, J
        integer, intent(in) :: lc(n), le(n), L(n), S(n)
        real, intent(in) :: sc(n), jc(n), K(n), se(n)
        real*8, intent(out) :: UFT(n,n)
        real*8 :: sj1, sj2
        integer :: ii, jj, deltas
        real*8, external :: SIXJ
        !integer :: kij
        UFT = 0.0d0
        do ii=1,n
            do jj=1,n
                deltas = kij(lc(ii),lc(jj))*kij(le(ii),le(jj))
                if(deltas.eq.0) then
                    UFT(ii,jj) = 0
                else
                    !print *, "Values of ang momentum ", lc(ii), le(ii), L(jj), K(ii), sc(ii), jc(ii)
                    sj1 = SIXJ(int(2*lc(ii)), int(2*le(ii)), int(2*L(jj)), int(2*K(ii)), int(2*sc(ii)), int(2*jc(ii)))
                    sj2 = SIXJ(int(2*sc(ii)), int(2*se(ii)), int(2*S(jj)), int(2*J), int(2*L(jj)), int(2*K(ii))) 
                    UFT(ii,jj) = (-1)**(le(ii) + L(jj) + lc(ii) + se(ii) + S(jj) + 3.*sc(ii))*&
                    sqrt((1 + 2*jc(ii))*(1 + 2*K(ii))*(1 + 2*L(jj))*(1 + 2*S(jj)))*sj1*sj2
                endif
            end do
        end do
        return
    end subroutine bJcK_kLS
    !
    subroutine bJcK_kJcs(n, lc, sc, jc, Jcs, le, K, se, J, UFT)
         
        integer, intent(in) :: n, J
        integer, intent(in) :: lc(n), le(n), Jcs(n)
        real, intent(in) :: sc(n), jc(n), se(n), K(n)
        real*8, intent(out) :: UFT(n,n)
        real*8 :: sj1
        integer :: ii, jj, deltas
        real*8, external :: SIXJ
        !integer :: kij
        UFT = 0.0d0
        do ii=1,n
            do jj=1,n
                deltas = kij(lc(ii),lc(jj))*kij(le(ii),le(jj))*kij(int(2*jc(ii)),int(2*jc(jj)))
                if(deltas.eq.0) then
                    UFT(jj,ii) = 0
                else
                    !print *, "Values of ang momentum ", lc(ii), le(ii), L(jj), K(ii), sc(ii), jc(ii)
                    sj1 = SIXJ(int(2*jc(ii)), int(2*le(ii)), int(2*K(ii)), int(2*J), int(2*se(ii)), int(2*Jcs(jj)))
                    UFT(jj,ii) = (-1)**(le(ii) + K(jj) + Jcs(jj) + se(ii))*&
                    sqrt((1 + 2*K(ii))*(1 + 2*Jcs(jj)))*sj1
                endif
            end do
        end do
        return
    end subroutine bJcK_kJcs
    !
    subroutine bjj_kJcK(n, lc, sc, jc, le, K, se, je, J, UFT)
         
        integer, intent(in) :: n, J
        integer, intent(in) :: lc(n), le(n)
        real, intent(in) :: sc(n), jc(n), se(n), K(n), je(n)
        real*8, intent(out) :: UFT(n,n)
        real*8 :: sj1
        integer :: ii, jj, deltas
        real*8, external :: SIXJ
        !integer :: kij
        UFT = 0.0d0
        do ii=1,n
            do jj=1,n
                deltas = kij(lc(ii),lc(jj))*kij(le(ii),le(jj))*kij(int(2*jc(ii)),int(2*jc(jj)))
                if(deltas.eq.0) then
                    UFT(ii,jj) = 0
                else
                    !print *, "Values of ang momentum ", lc(ii), le(ii), L(jj), K(ii), sc(ii), jc(ii)
                    sj1 = SIXJ(int(2*le(ii)), int(2*se(ii)), int(2*je(ii)), int(2*J), int(2*jc(ii)), int(2*K(jj)))
                    UFT(ii,jj) = (-1)**(je(ii) + jc(ii) - J)*&
                    sqrt((1 + 2*K(jj))*(1 + 2*je(ii)))*sj1
                endif
            end do
        end do
        return
    end subroutine bjj_kJcK

    subroutine bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, UFT)
         
        integer, intent(in) :: n, J
        integer, intent(in) :: lc(n), le(n), L(n), S(n), Jcs(n)
        real, intent(in) :: sc(n), jc(n), se(n)
        real*8, intent(out) :: UFT(n,n)
        real*8 :: sj1, sj2
        integer :: ii, jj, deltas
        real*8, external :: SIXJ
        !integer :: kij
        UFT = 0.0d0
        do ii=1,n
            do jj=1,n
                deltas = kij(lc(ii),lc(jj))*kij(le(ii),le(jj))
                if(deltas.eq.0) then
                    UFT(ii,jj) = 0
                else
                    !print *, "Values of ang momentum ", lc(ii), le(ii), L(jj), Jcs(ii), sc(ii), jc(ii)
                    sj1 = SIXJ(int(2*sc(ii)), int(2*se(ii)), int(2*S(jj)), int(2*Jcs(ii)), int(2*lc(ii)), int(2*jc(ii)))
                    sj2 = SIXJ(int(2*lc(ii)), int(2*le(ii)), int(2*L(jj)), int(2*J), int(2*S(jj)), int(2*Jcs(ii)))
                    UFT(ii,jj) = (-1)**(se(ii) + jc(ii) + le(ii) -J)*&
                    sqrt((1 + 2*jc(ii))*(1 + 2*Jcs(ii))*(1 + 2*L(jj))*(1 + 2*S(jj)))*sj1*sj2
                endif
            end do
        end do
        return
    end subroutine bJcs_kLS

    subroutine bjj_kJcs(n, lc, sc, jc, le, Jcs, se, je, J, UFT)
         
        integer, intent(in) :: n, J
        integer, intent(in) :: lc(n), le(n)
        real, intent(in) :: sc(n), jc(n), se(n), je(n), Jcs(n)
        real*8, intent(out) :: UFT(n,n)
        real*8 :: sj1
        integer :: ii, jj, deltas
        real*8, external :: SIXJ
        !integer :: kij
        UFT = 0.0d0
        do ii=1,n
            do jj=1,n
                deltas = kij(lc(ii),lc(jj))*kij(le(ii),le(jj))*kij(int(2*jc(ii)),int(2*jc(jj)))
                if(deltas.eq.0) then
                    UFT(ii,jj) = 0
                else
                    !print *, "Values of ang momentum ", lc(ii), le(ii), L(jj), K(ii), sc(ii), jc(ii)
                    sj1 = SIXJ(int(2*se(ii)), int(2*jc(ii)), int(2*Jcs(jj)), int(2*J), int(2*le(ii)), int(2*je(ii)))
                    UFT(ii,jj) = (-1)**(2*(je(ii) + jc(ii))+se(ii)+jc(ii)+le(ii) - J)*&
                    sqrt((1 + 2*Jcs(jj))*(1 + 2*je(ii)))*sj1
                endif
            end do
        end do
        return
    end subroutine bjj_kJcs

    subroutine bjj_kLS(n, lc, sc, jc, le, se, je, L, S, J, UFT)
         
        integer, intent(in) :: n, J
        integer, intent(in) :: lc(n), le(n), L(n), S(n)
        real, intent(in) :: sc(n), jc(n), se(n), je(n)
        real*8, intent(out) :: UFT(n,n)
        real*8 :: nj1
        integer :: ii, jj, deltas
        real*8, external :: XNINEJ
        !integer, external :: kij
        UFT = 0.0d0
        do ii=1,n
            do jj=1,n
                deltas = kij(lc(ii),lc(jj))*kij(le(ii),le(jj))
                if(deltas.eq.0) then
                    UFT(ii,jj) = 0
                else
                    !print *, "Values of ang momentum ", lc(ii), le(ii), L(jj), K(ii), sc(ii), jc(ii)
                    nj1 = XNINEJ(int(2*lc(ii)),int(2*sc(ii)),int(2*jc(ii)), int(2*le(ii)), int(2*se(ii)),&
                            int(2*je(ii)), int(2*L(jj)), int(2*S(jj)), int(2*J))
                    UFT(ii,jj) =sqrt((1+2*je(ii))*(1+2*jc(ii))*(1 + 2*L(jj))*(1 + 2*S(jj)))*nj1
                endif
            end do
        end do
        return
    end subroutine bjj_kLS

    function kij(ip,jp)
        implicit none
        integer, intent(in) :: ip,jp
        integer :: kij
        if(ip.eq.0.and.jp.eq.0) then
            kij = 1
        else
            kij = int(float((ip+jp)-abs(ip-jp))/float((ip+jp)+abs(ip-jp)))
        endif
        return
    end function kij
    
end module