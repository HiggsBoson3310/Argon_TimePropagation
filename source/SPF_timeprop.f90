module spf_prop
    use errfun
    use arj1j2
    use matmod
    implicit none
    integer, parameter :: nos1=4, nos3=3, nos=nos1+nos3, nosgs=nos+1, nbs2f = 24, nbs2p=30, nbs4 = 26, nbs0=17, &
    nbs=nbs2p+nbs2f+nbs0+nbs4, nbs1 = nbs2p+nbs2f+nbs0, &
    ncs=5+6+6, ncs1 = 5+6, ne=300, nau=36, ntot = nos+nbs+ncs*ne, ntot1 = nos1+nbs1+ncs1*ne, &
    ntotau=ntot+5*nau, nosau=nosgs
    real*8, private  :: Mb(nbs, nos)
    complex*16, private :: Mc(nbs, ncs)
    integer,private :: spi, spi1, spi2, spi3, j, k 
    !TODO: initialize the damn matrix elements.
    contains
    !
    subroutine SPF_step(n, HD, V, dt, cn, cnp1, Es)
        integer, INTENT(IN) :: n
        complex*16, INTENT(IN) :: V(n,n)
        complex*16, INTENT(INOUT) :: cn(n)
        real*8, INTENT(IN) :: dt, Es(n),  HD(n)
        complex*16, INTENT(OUT) :: cnp1(n)

        do spi=1, n
            cnp1(spi) =  EXP(-Im*Es(spi)*dt/2.0d0)*cn(spi)
        end do
        
        call zmatHvec(n,n,V,cnp1,cn)

        do spi=1,n
            cnp1(spi) = EXP(-Im*HD(spi))*cn(spi)
        end do

        call zmatvec(n,n,V,cnp1,cn)

        do spi=1, n
            cnp1(spi) = EXP(-Im*Es(spi)*dt/2.0d0)*cn(spi)
        end do

        RETURN
    end subroutine SPf_step
    !
    subroutine SPF_step_slow(n, HD, V, dt, cn, cnp1, Es)
        integer, INTENT(IN) :: n
        complex*16, INTENT(IN) :: V(n,n)
        complex*16, INTENT(INOUT) :: cn(n)
        real*8, INTENT(IN) :: dt, Es(n),  HD(n)
        complex*16, INTENT(OUT) :: cnp1(n)

        do spi=1, n
            cnp1(spi) =  EXP(-Im*Es(spi)*dt/2.0d0)*cn(spi)
        end do
        
        cnp1 = matmul(CONJG(transpose(V)),cnp1)

        do spi=1,n
            cnp1(spi) = EXP(-Im*HD(spi))*cnp1(spi)
        end do

        cnp1 = matmul(V,cnp1)

        do spi=1, n
            cnp1(spi) = EXP(-Im*Es(spi)*dt/2.0d0)*cnp1(spi)
        end do

        RETURN
    end subroutine SPF_step_slow
    !
    function Mt(t,to,gam, w, Eo, Ebs, Exi)
        real*8, intent(IN) :: t, to, gam, w, Eo, Ebs(nos), Exi(nbs)
        complex*16 :: exnt(nos), erfarg(nos), sume
        complex*16 :: exntc, sumc
        complex*16 :: Mt(nbs,nbs)
        integer :: iloc
        Mt = zer
        do iloc=1,nbs
            do j=1, nbs
                exnt = IM*(Exi(iloc)-Ebs-w)*t-(t-to)**2/gam**2 + IM*(Ebs+w-Exi(j))*to
                exnt = exnt-gam**2/4 * (Ebs+w-Exi(j))**2
                erfarg = (t-to)/gam + Im*(gam/2.0d0*(Exi(j)-Ebs-w))
                sume = Zer
                do k=1,nos
                    sume = sume + Mb(iloc,k)*Mb(j,k)*EXP(exnt(k))*(Zone+erf_pop(erfarg(k)))
                end do
                Mt(iloc,j) = SQRT(PI)*gam/2.0d0 * sume
                exntc = IM*(Exi(iloc)-Exi(j))*t-2*(t-to)**2/gam**2
                sumc = Zer
                do k=1,ncs
                    sumc = sumc + CONJG(Mc(iloc,k))*Mc(j,k)*EXP(exntc)
                end do
                Mt(iloc,j) = Mt(iloc,j) + 2*PI*sumc

                Mt(iloc,j) = -1.*Eo**2 * 1/4. * Mt(iloc,j)
            end do
        end do
        return
    end function Mt
    !
    function Dtf(t,to,Ao,gam,w,Eo, Ebs, Exi)
        real*8, INTENT(IN) :: t,to,gam,w,Eo, Ebs(nos), Exi(nbs)
        complex*16, INTENT(IN) :: Ao(nos)
        complex*16 :: sume
        complex*16 :: Dtf(nbs)
        integer :: iloc
        do iloc=1,nbs
            sume = Zer
            do j =1, nos
                sume = sume+Ao(j)*EXP(Im*(Exi(iloc)-Ebs(k)-w)*t+Im*w*to-(t-to)**2/gam**2)*Mb(iloc,j)*Eo/2.0d0
            end do
            Dtf(iloc) = -Im*sume
        end do
        return
    end function Dtf
    !
    subroutine RKsolution(nt,t,to,Bo,gam,w,Eo,Ao)
        integer :: nt
        real*8 :: t(nt), to, gam, w, Eo, dt, Ebs(nos), Exi(nbs)
        complex*16 :: Ao(nbs), Bo(nbs), Bn(nbs)
        complex*16 :: sols(nt, nbs), Mtt(3, nbs, nbs)
        complex*16 :: Dtt(3, nbs)
        dt = (t(nt)-t(1))/nt
        sols=Zer
        sols(1,:) = Bo
        Bn = Bo
        do spi=2, nt
            Mtt(1,:,:) = Mt(t(spi), to, gam, w, Eo, Ebs, Exi)
            Dtt(1,:) = Dtf(t(spi), to, Ao, gam, w, Eo, Ebs, Exi)

        end do
    end subroutine
    !
    function ztrapz(n, y, x)
        integer :: n
        real*8 :: x(n), dx, fmin, fmax
        complex*16 :: y(n)
        complex*16 :: ztrapz
        dx = (x(n)-x(1))/n
        ztrapz = Zer
        do spi=1, n-1
            fmin = min(REAL(y(spi)), REAL(y(spi+1)))
            fmax = max(REAL(y(spi)), REAL(y(spi+1)))
            ztrapz = ztrapz + (fmin*dx + (fmax-fmin)*dx/2)*Zone
            fmin = min(AIMAG(y(spi)), AIMAG(y(spi+1)))
            fmax = max(AIMAG(y(spi)), AIMAG(y(spi+1)))
            ztrapz = ztrapz + (fmin*dx + (fmax-fmin)*dx/2)*IM
        end do
        return
    end function ztrapz
    !
    subroutine TPMat(evec, Elo, dEcm,&
        bbJ2f, bbJ2p, bbJ0p, Zs, ZJ2f, ZJ2p, ZJ0p, &
        regJ2fJ3, iregJ2fJ3, regJ2fJ1, iregJ2fJ1, regJ2pJ3, iregJ2pJ3, regJ2pJ1, iregJ2pJ1,&
        regJ0pJ3, iregJ0pJ3, regJ0pJ1, iregJ0pJ1, H)
        real*8, INTENT(IN) :: evec(3)
        real*8, INTENT(IN) :: Elo, dEcm
        real*8, INTENT(IN) :: Zs(nos,5), ZJ2f(nbs2f, 3), Zj2p(nbs2p, 3), ZJ0p(nbs0, 2)
        real*8, INTENT(IN) :: bbJ2f(nbs2f, nos1, 3, 5), bbJ2p(nbs2p, nos1, 3, 5), bbJ0p(nbs0, nos1, 2, 5)
        real*8, INTENT(IN) :: regJ2fJ3(ne,nbs2f, 6, 3), iregJ2fJ3(ne, nbs2p, 6, 3), &
        regJ2fJ1(ne, nbs2f, 5, 3), iregJ2fJ1(ne, nbs2f, 5, 3)
        real*8, INTENT(IN) :: regJ2pJ3(ne,nbs2p, 6, 3), iregJ2pJ3(ne, nbs2p, 6, 3), &
        regJ2pJ1(ne, nbs2p, 5, 3), iregJ2pJ1(ne, nbs2p, 5, 3)
        real*8, INTENT(IN) :: regJ0pJ3(ne,nbs0,6, 2), iregJ0pJ3(ne, nbs0,6, 2), &
        regJ0pJ1(ne, nbs0,5, 2), iregJ0pJ1(ne, nbs0,5, 2)
        real*8 :: Econt
        integer :: l2f(3), l2p(3), l0p(2)
        integer :: lJ1(5), lJ3(6), jcs2(3), jcs0(2), jcsJ1(5), jcsJ3(6)
        real :: jc2(3), jc0(2), jcJ1(5), jcJ3(6) 
        real*8, allocatable :: Zx(:), Zy(:)
        complex*16, allocatable :: Sdmat(:,:)
        real*8, allocatable :: Kmat(:,:)
        complex*16, intent(out) :: H(ntot, ntot)
        integer :: counterx, countery
        counterx = 0
        countery = 0
        l2f = (/ 3, 3, 3 /)
        l2p = (/ 1, 1, 1/)
        jc2 = (/ 3./2., 3./2., 1./2./)
        jcs2 =  (/ 2, 1, 1 /)
        l0p = (/1, 1/)
        jc0 = (/3/2., 1/2./)
        jcs0 = (/1, 1/)
        jcJ1 = jcb
        lJ1 = lb
        jcsJ1 = jcsb
        jcJ3 = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        lJ3 =  (/ 2, 2, 2, 4, 4, 4/)
        jcsJ3 = (/ 1, 2, 1, 1, 2, 1/)
        counterx = 0
        countery = 0
        H = Zer

        ! The matrix is really sparse, so we start filling by parts.
        ! Start with the coupling between the initial states to the intermediates:
        do counterx=1,nos
            allocate(Zx(5))
            Zx = Zs(counterx,:)
            do countery =  nos+1, nos+nbs
                if(countery.le.nos+nbs2f) then
                    allocate(Zy(3))
                    spi1 = counterx
                    spi2 = countery-nos
                    Zy = ZJ2f(spi2,:) 
                    H(counterx, countery) = RedElem_b(3, Zy, 5, Zx, jcb, jcsb, lb, 1, jc2, jcs2, l2f, 2, bbJ2f(spi2,spi1,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum(2, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nos+nbs2f+nbs2p) then
                    allocate(Zy(3))
                    spi1 = counterx
                    spi2 = countery-nos-nbs2f
                    Zy = ZJ2p(spi2,:)
                    H(counterx, countery) = RedElem_b(3, Zy, 5, Zx, jcb, jcsb, lb, 1, jc2, jcs2, l2p, 2, bbJ2p(spi2,spi1,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum(2, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                else
                    allocate(Zy(2))
                    spi1 = counterx
                    spi2 = countery-nos-nbs2f-nbs2p
                    Zy = ZJ0p(spi2,:)
                    H(counterx, countery) = RedElem_b(2, Zy, 5, Zx, jcb, jcsb, lb, 1, jc0, jcs0, l0p, 0, bbJ0p(spi2,spi1,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(0-0) * qsum(0,0,evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                endif
                deallocate(Zy)
                H(countery, counterx) = H(counterx, countery)
            end do
            deallocate(Zx)
        end do

        ! Then we just have to complete the coupling to the continuum. 
        ! We have to consider only the continuum states
        write(6,*) 'Combining with the continuum'
        write(6,*) 'Density in atomic units', dEcm/CMperAU
        do counterx=nos+1, nos+nbs
            if(counterx.le.nos+nbs2f) then
                allocate(Zx(3))
                spi1 = counterx-nos
                Zx = ZJ2f(spi1,:)
            elseif(counterx.le.nos+nbs2p+nbs2f) then
                allocate(Zx(3))
                spi1=counterx-nos-nbs2f
                Zx = Zj2p(spi1,:)
            else
                allocate(Zx(2))
                spi1 = counterx-nos-nbs2f-nbs2p
                Zx = ZJ0p(spi1,:)
            endif
            do spi2=1, ne
                Econt = Elo+spi2*dEcm
                do spi3=1, 5
                    allocate(Sdmat(5,5))
                    allocate(Kmat(5,5))
                    call J1_MQDT_at(Econt,Sdmat,Kmat,0)
                    countery = nos+nbs+(spi3-1)*ne+spi2
                    if(counterx.le.nos+nbs2f) then
                        H(counterx, countery) = CONJG(RedElem_c(5, 3, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc2,&
                                                 l2f, jcs2, 2, regJ2fJ1(spi2,spi1,:,:), iregJ2fJ1(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(1.0d0)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nos+nbs2f+nbs2p) then
                        H(counterx, countery) = CONJG(RedElem_c(5, 3, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc2,&
                                                 l2p, jcs2, 2, regJ2pJ1(spi2,spi1,:,:), iregJ2pJ1(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(1.0d0)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    else
                        H(counterx, countery) = CONJG(RedElem_c(5, 2, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc0,&
                                                 l0p, jcs0, 0, regJ0pJ1(spi2,spi1,:,:), iregJ0pJ1(spi2,spi1,:,:)))
                        !write(6,*) "before the sum ", spi1, spi2, counterx, countery, H(counterx, countery) 
                        H(counterx, countery) = SQRT(1.0d0)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,0,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    endif
                    deallocate(Sdmat)
                    deallocate(Kmat)
                enddo
                do spi3=1, 6
                    allocate(Sdmat(6,6))
                    allocate(Kmat(6,6))
                    call J3_MQDT_at(Econt,Sdmat,Kmat,0)
                    countery = nos+nbs+5*ne+(spi3-1)*ne+spi2
                    !write(6,*) "Channel ",spi3, countery
                    if(counterx.le.nos+nbs2f) then
                        H(counterx, countery) = CONJG(RedElem_c(6, 3, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc2,&
                                                 l2f, jcs2, 2, regJ2fJ3(spi2,spi1,:,:), iregJ2fJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(1.0d0)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nos+nbs2f+nbs2p) then
                        H(counterx, countery) = CONJG(RedElem_c(6, 3, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc2,&
                                                 l2p, jcs2, 2, regJ2pJ3(spi2,spi1,:,:), iregJ2pJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(1.0d0)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    else
                        H(counterx, countery) = CONJG(RedElem_c(6, 2, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc0,&
                                                 l0p, jcs0, 0, regJ0pJ3(spi2,spi1,:,:), iregJ0pJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(1.0d0)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,0,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    endif
                    deallocate(Sdmat)
                    deallocate(Kmat)
                enddo
            end do
            deallocate(Zx)
        end do

        open(38, file='SPF_matrices/SPF_re.dat', access='sequential')
        open(39, file='SPF_matrices/SPF_im.dat', access='sequential')
        do spi=1, ntot
            write(38, 1985) ( REAL(H(spi,j)),j=spi,ntot)
            write(39, 1985) (AIMAG(H(spi,j)),j=spi,ntot)
        end do
        close(38)
        close(39)
        return 
        1985 format(2500(E15.8,1x))
    end subroutine TPMat
    !
    subroutine TPMatED(evec, Elo, dEcm,&
        bbJ1J0p, bbJ1J2f, bbJ1J2p, bbJ3J2f, bbJ3J2p, bbJ3J4l, &
        Z1, Z3, ZJ2f, ZJ2p, ZJ0p, ZJ4l, &
        regJ0pJ1, iregJ0pJ1, regJ2fJ1, iregJ2fJ1, regJ2pJ1, iregJ2pJ1,&
        regJ2pJ3, iregJ2pJ3, regJ2fJ3, iregJ2fJ3, regJ4lJ3, iregJ4lJ3, &
        regJ4lJ5, iregJ4lJ5, H)
        ! Polarization vectors
        real*8, INTENT(IN) :: evec(3)
        ! Continuum energy parameters
        real*8, INTENT(IN) :: Elo, dEcm
        real*8 :: Econt
        ! Z coeffs
        real*8, INTENT(IN) :: Z1(nos1,5), Z3(nos3,6), ZJ2f(nbs2f, 3), ZJ2p(nbs2p, 3), ZJ0p(nbs0, 2), &
        ZJ4l(nbs4,6)
        ! Bound bound radial integrals.
        real*8, INTENT(IN) :: bbJ1J0p(nos1, nbs0, 2, 5), bbJ1J2f(nos1, nbs2f, 3, 5), bbJ1J2p(nos1, nbs2p, 3, 5)
        real*8, INTENT(IN) :: bbJ3J2f(nos3, nbs2f, 3, 6), bbJ3J2p(nos3, nbs2p, 3, 6), bbJ3J4l(nos3, nbs4, 6, 6)
        ! Bound Continuum radial integrals.
        real*8, INTENT(IN) :: regJ0pJ1(ne, nbs0, 5, 2), iregJ0pJ1(ne, nbs0, 5, 2)
        real*8, INTENT(IN) :: regJ2fJ3(ne,nbs2f, 6, 3), iregJ2fJ3(ne, nbs2f, 6, 3), &
        regJ2fJ1(ne, nbs2f, 5, 3), iregJ2fJ1(ne, nbs2f, 5, 3)
        real*8, INTENT(IN) :: regJ2pJ3(ne,nbs2p, 6, 3), iregJ2pJ3(ne, nbs2p, 6, 3), &
        regJ2pJ1(ne, nbs2p, 5, 3), iregJ2pJ1(ne, nbs2p, 5, 3)
        real*8, INTENT(IN) :: regJ4lJ3(ne, nbs4, 6, 6), iregJ4lJ3(ne, nbs4, 6, 6), &
        regJ4lJ5(ne, nbs4, 6, 6), iregJ4lJ5(ne, nbs4, 6, 6)
        ! Angular quantum numbers 
        integer :: l2f(3), l2p(3), l0p(2), l4l(6), lJ1(5), lJ3(6), lJ5(6)
        integer :: jcs2(3), jcs0(2), jcs4(6), jcsJ1(5), jcsJ3(6), jcsJ5(6)
        real :: jc2(3), jc0(2), jc4(6), jcJ1(5), jcJ3(6), jcJ5(6) 
        real*8, allocatable :: Zx(:), Zy(:)
        complex*16, allocatable :: Sdmat(:,:)
        real*8, allocatable :: Kmat(:,:)
        complex*16, intent(out) :: H(ntot, ntot)
        integer :: counterx, countery
        counterx = 0
        countery = 0

        l2f = (/ 3, 3, 3 /)
        l2p = (/ 1, 1, 1/)
        jc2 = (/ 3./2., 3./2., 1./2./)
        jcs2 =  (/ 2, 1, 1 /)
        
        l0p = (/1, 1/)
        jc0 = (/3/2., 1/2./)
        jcs0 = (/1, 1/)

        l4l = (/3, 3, 3, 5, 5, 5/)
        jc4 = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        jcs4 = (/1, 2, 1, 1, 2, 1/)

        lJ1 = lb
        jcJ1 = jcb
        jcsJ1 = jcsb

        lJ3 = (/ 2, 2, 2, 4, 4, 4/)
        jcJ3 = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        jcsJ3 = (/ 1, 2, 1, 1, 2, 1/)

        lJ5 = (/4,4,4,6,6,6/)
        jcJ5 = (/3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        jcsJ5 = (/1, 2, 1, 1, 2, 1/)

        H = Zer

        ! The matrix is really sparse, so we start filling by parts.
        ! Start with the coupling between the initial states to the intermediates:
        do counterx=1,nos1
            allocate(Zx(5))
            Zx = Z1(counterx,:)
            do countery =  nos+1, nos+nbs
                if(countery.le.nos+nbs2f) then
                    allocate(Zy(3))
                    spi1 = counterx
                    spi2 = countery-nos
                    Zy = ZJ2f(spi2,:) 
                    H(counterx, countery) = RedElem_b(3, Zy, 5, Zx, jcJ1, jcsJ1, lJ1, 1, jc2, jcs2, l2f, 2, bbJ1J2f(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum(2, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nos+nbs2f+nbs2p) then
                    allocate(Zy(3))
                    spi1 = counterx
                    spi2 = countery-nos-nbs2f
                    Zy = ZJ2p(spi2,:)
                    H(counterx, countery) = RedElem_b(3, Zy, 5, Zx, jcJ1, jcsJ1, lJ1, 1, jc2, jcs2, l2p, 2, bbJ1J2p(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum(2, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nos+nbs2p+nbs2f+nbs0) then
                    allocate(Zy(2))
                    spi1 = counterx
                    spi2 = countery-nos-nbs2p-nbs2f
                    Zy = ZJ0p(spi2,:)
                    H(counterx, countery) = RedElem_b(2, Zy, 5, Zx, jcJ1, jcsJ1, lJ1, 1, jc0, jcs0, l0p, 0, bbJ1J0p(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(0-0) * qsum(0,0,evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                else
                    allocate(Zy(6))
                    Zy = 0.d0
                    H(counterx, countery) = Zer
                endif
                deallocate(Zy)
                H(countery, counterx) = H(counterx, countery)
            end do
            deallocate(Zx)
        end do
        ! Now the initial J=3 states
        do counterx=nos1+1,nos
            allocate(Zx(6))
            spi1 = counterx-nos1
            Zx = Z3(spi1,:)
            do countery =  nos+1, nos+nbs
                if(countery.le.nos+nbs2f) then
                    allocate(Zy(3))
                    spi2 = countery-nos
                    Zy = ZJ2f(spi2,:) 
                    H(counterx, countery) = RedElem_b(3, Zy, 6, Zx, jcJ3, jcsJ3, lJ3, 3, jc2, jcs2, l2f, 2, bbJ3J2f(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum_c(2, 0, 3, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nos+nbs2f+nbs2p) then
                    allocate(Zy(3))
                    spi2 = countery-nos-nbs2f
                    Zy = ZJ2p(spi2,:)
                    H(counterx, countery) = RedElem_b(3, Zy, 6, Zx, jcJ3, jcsJ3, lJ3, 3, jc2, jcs2, l2p, 2, bbJ3J2p(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum_c(2, 0, 3, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nos+nbs2f+nbs2p+nbs0) then
                    allocate(Zy(2))
                    Zy = 0.d0
                    H(counterx, countery) = Zer
                else
                    allocate(Zy(6))
                    spi2 = countery-nos-nbs2p-nbs2f-nbs0
                    Zy = ZJ4l(spi2,:)
                    H(counterx, countery) = RedElem_b(6, Zy, 6, Zx, jcJ3, jcsJ3, lJ3, 4, jc4, jcs4, l4l, 4, bbJ3J4l(spi1,spi2,:,:))
                !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                H(counterx, countery) = H(counterx, countery)*(-1)**(4-0) * qsum_c(4,0,3,0,evec)
                !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                endif
                deallocate(Zy)
                H(countery, counterx) = H(counterx, countery)
            end do
            deallocate(Zx)
        end do

        ! Then we just have to complete the coupling to the continuum. 
        ! We have to consider only the continuum states
        write(6,*) 'Combining with the continuum'
        write(6,*) 'Density in atomic units', dEcm/CMperAU
        ! J2f states
        do counterx=nos+1, nos+nbs
            if(counterx.le.nos+nbs2f) then
                allocate(Zx(3))
                spi1 = counterx-nos
                Zx = ZJ2f(spi1,:)
            elseif(counterx.le.nos+nbs2p+nbs2f) then
                allocate(Zx(3))
                spi1=counterx-nos-nbs2f
                Zx = Zj2p(spi1,:)
            elseif(counterx.le.nos+nbs2f+nbs2p+nbs0) then
                allocate(Zx(2))
                spi1 = counterx-nos-nbs2f-nbs2p
                Zx = ZJ0p(spi1,:)
            else
                allocate(Zx(6))
                spi1 = counterx-nos-nbs2p-nbs2f-nbs0
                Zx = ZJ4l(spi1,:)
            endif
            do spi2=1, ne
                Econt = Elo+spi2*dEcm
                do spi3=1, 5
                    allocate(Sdmat(5,5))
                    allocate(Kmat(5,5))
                    call J1_MQDT_at(Econt,Sdmat,Kmat,0)
                    countery = nos+nbs+(spi3-1)*ne+spi2
                    if(counterx.le.nos+nbs2f) then
                        H(counterx, countery) = CONJG(RedElem_c(5, 3, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc2,&
                                                 l2f, jcs2, 2, regJ2fJ1(spi2,spi1,:,:), iregJ2fJ1(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nos+nbs2f+nbs2p) then
                        H(counterx, countery) = CONJG(RedElem_c(5, 3, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc2,&
                                                 l2p, jcs2, 2, regJ2pJ1(spi2,spi1,:,:), iregJ2pJ1(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nos+nbs2f+nbs2p+nbs0) then
                        H(counterx, countery) = CONJG(RedElem_c(5, 2, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc0,&
                                                 l0p, jcs0, 0, regJ0pJ1(spi2,spi1,:,:), iregJ0pJ1(spi2,spi1,:,:)))
                        !write(6,*) "before the sum ", spi1, spi2, counterx, countery, H(counterx, countery) 
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,0,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    else 
                        H(counterx, countery) = Zer
                        H(countery, counterx) = Zer
                    endif
                    deallocate(Sdmat)
                    deallocate(Kmat)
                enddo
                do spi3=1, 6
                    allocate(Sdmat(6,6))
                    allocate(Kmat(6,6))
                    call J3_MQDT_at(Econt,Sdmat,Kmat,0)
                    countery = nos+nbs+5*ne+(spi3-1)*ne+spi2
                    !write(6,*) "Channel ",spi3, countery
                    if(counterx.le.nos+nbs2f) then
                        H(counterx, countery) = CONJG(RedElem_c(6, 3, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc2,&
                                                 l2f, jcs2, 2, regJ2fJ3(spi2,spi1,:,:), iregJ2fJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nos+nbs2f+nbs2p) then
                        H(counterx, countery) = CONJG(RedElem_c(6, 3, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc2,&
                                                 l2p, jcs2, 2, regJ2pJ3(spi2,spi1,:,:), iregJ2pJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nos+nbs2f+nbs2p+nbs0) then
                        H(counterx, countery) = Zer
                        H(countery, counterx) = Zer
                    else
                        H(counterx, countery) = CONJG(RedElem_c(6, 6, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc4,&
                                                 l4l, jcs4, 4, regJ4lJ3(spi2,spi1,:,:), iregJ4lJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,4,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    endif
                    deallocate(Sdmat)
                    deallocate(Kmat)
                enddo
                do spi3=1,6
                    allocate(Sdmat(6,6))
                    allocate(Kmat(6,6))
                    call J5_MQDT_at(Econt,Sdmat,Kmat,0)
                    countery = nos+nbs+(5+6)*ne+(spi3-1)*ne+spi2
                    if(counterx.gt.nos+nbs2f+nbs2p+nbs0) then
                        H(counterx, countery) = CONJG(RedElem_c(6, 6, Zx, Sdmat, spi3, jcJ5, lJ5, jcsJ5, 5, jc4,&
                                                 l4l, jcs4, 4, regJ4lJ5(spi2,spi1,:,:), iregJ4lJ5(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(5-(0)) * qsum_c(5,0,4,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    else 
                        H(counterx, countery) = Zer
                        H(countery, counterx) = Zer
                    endif
                    deallocate(Sdmat)
                    deallocate(Kmat)
                end do
            end do
            deallocate(Zx)
        end do
    end subroutine TPMatED
    !
    subroutine TPMatEDwGS(evec, Elo, dEcm,&
        bbJ1J0p, bbJ1J2f, bbJ1J2p, bbJ3J2f, bbJ3J2p, bbJ3J4l, &
        Z1, Z3, ZJ2f, ZJ2p, ZJ0p, ZJ4l, &
        regJ0pJ1, iregJ0pJ1, regJ2fJ1, iregJ2fJ1, regJ2pJ1, iregJ2pJ1,&
        regJ2pJ3, iregJ2pJ3, regJ2fJ3, iregJ2fJ3, regJ4lJ3, iregJ4lJ3, &
        regJ4lJ5, iregJ4lJ5, H)
        ! Polarization vectors
        real*8, INTENT(IN) :: evec(3)
        ! Continuum energy parameters
        real*8, INTENT(IN) :: Elo, dEcm
        real*8 :: Econt
        ! Z coeffs
        real*8, INTENT(IN) :: Z1(nos1,5), Z3(nos3,6), ZJ2f(nbs2f, 3), ZJ2p(nbs2p, 3), ZJ0p(nbs0, 2), &
        ZJ4l(nbs4,6)
        ! Dipoles with the gs
        real*8 :: Dwgs(4)
        ! Bound bound radial integrals.
        real*8, INTENT(IN) :: bbJ1J0p(nos1, nbs0, 2, 5), bbJ1J2f(nos1, nbs2f, 3, 5), bbJ1J2p(nos1, nbs2p, 3, 5)
        real*8, INTENT(IN) :: bbJ3J2f(nos3, nbs2f, 3, 6), bbJ3J2p(nos3, nbs2p, 3, 6), bbJ3J4l(nos3, nbs4, 6, 6)
        ! Bound Continuum radial integrals.
        real*8, INTENT(IN) :: regJ0pJ1(ne, nbs0, 5, 2), iregJ0pJ1(ne, nbs0, 5, 2)
        real*8, INTENT(IN) :: regJ2fJ3(ne,nbs2f, 6, 3), iregJ2fJ3(ne, nbs2f, 6, 3), &
        regJ2fJ1(ne, nbs2f, 5, 3), iregJ2fJ1(ne, nbs2f, 5, 3)
        real*8, INTENT(IN) :: regJ2pJ3(ne,nbs2p, 6, 3), iregJ2pJ3(ne, nbs2p, 6, 3), &
        regJ2pJ1(ne, nbs2p, 5, 3), iregJ2pJ1(ne, nbs2p, 5, 3)
        real*8, INTENT(IN) :: regJ4lJ3(ne, nbs4, 6, 6), iregJ4lJ3(ne, nbs4, 6, 6), &
        regJ4lJ5(ne, nbs4, 6, 6), iregJ4lJ5(ne, nbs4, 6, 6)
        ! Angular quantum numbers 
        integer :: l2f(3), l2p(3), l0p(2), l4l(6), lJ1(5), lJ3(6), lJ5(6)
        integer :: jcs2(3), jcs0(2), jcs4(6), jcsJ1(5), jcsJ3(6), jcsJ5(6)
        real :: jc2(3), jc0(2), jc4(6), jcJ1(5), jcJ3(6), jcJ5(6) 
        real*8, allocatable :: Zx(:), Zy(:)
        complex*16, allocatable :: Sdmat(:,:)
        real*8, allocatable :: Kmat(:,:)
        complex*16, intent(out) :: H(ntot+1, ntot+1)
        integer :: counterx, countery
        counterx = 0
        countery = 0
    
        l2f = (/ 3, 3, 3 /)
        l2p = (/ 1, 1, 1/)
        jc2 = (/ 3./2., 3./2., 1./2./)
        jcs2 =  (/ 2, 1, 1 /)
        
        l0p = (/1, 1/)
        jc0 = (/3/2., 1/2./)
        jcs0 = (/1, 1/)
    
        l4l = (/3, 3, 3, 5, 5, 5/)
        jc4 = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        jcs4 = (/1, 2, 1, 1, 2, 1/)
    
        lJ1 = lb
        jcJ1 = jcb
        jcsJ1 = jcsb
    
        lJ3 = (/ 2, 2, 2, 4, 4, 4/)
        jcJ3 = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        jcsJ3 = (/ 1, 2, 1, 1, 2, 1/)
    
        lJ5 = (/4,4,4,6,6,6/)
        jcJ5 = (/3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        jcsJ5 = (/1, 2, 1, 1, 2, 1/)
    
        H = Zer
        Dwgs = (/0.278603d0, 0.518247d0, 0.18458d0, 0.54917d0/)
        ! Fill the coupling between the ground state and the J=1 bound states
        do countery=1,nos1
            H(1,countery+1) = Dwgs(countery)
            H(countery+1,1) = H(1,countery+1)
        end do
        
    
        ! The matrix is really sparse, so we start filling by parts.
        ! Start with the coupling between the initial states to the intermediates:
        do counterx=2,nos1+1
            allocate(Zx(5))
            spi1 = counterx-1
            Zx = Z1(spi1,:)
            do countery =  nosgs+1, nosgs+nbs
                if(countery.le.nosgs+nbs2f) then
                    allocate(Zy(3))
                    spi2 = countery-nosgs
                    Zy = ZJ2f(spi2,:) 
                    H(counterx, countery) = RedElem_b(3, Zy, 5, Zx, jcJ1, jcsJ1, lJ1, 1, jc2, jcs2, l2f, 2, bbJ1J2f(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum(2, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nosgs+nbs2f+nbs2p) then
                    allocate(Zy(3))
                    spi2 = countery-nosgs-nbs2f
                    Zy = ZJ2p(spi2,:)
                    H(counterx, countery) = RedElem_b(3, Zy, 5, Zx, jcJ1, jcsJ1, lJ1, 1, jc2, jcs2, l2p, 2, bbJ1J2p(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum(2, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nosgs+nbs2f+nbs2p+nbs0) then
                    allocate(Zy(2))
                    spi2 = countery-nosgs-nbs2f-nbs2p
                    Zy = ZJ0p(spi2,:)
                    H(counterx, countery) = RedElem_b(2, Zy, 5, Zx, jcJ1, jcsJ1, lJ1, 1, jc0, jcs0, l0p, 0, bbJ1J0p(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(0-0) * qsum(0,0,evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                else
                    allocate(Zy(6))
                    Zy = 0.d0
                    H(counterx, countery) = Zer
                endif
                deallocate(Zy)
                H(countery, counterx) = H(counterx, countery)
            end do
            deallocate(Zx)
        end do
        ! Now the initial J=3 states
        do counterx=nos1+2,nosgs
            allocate(Zx(6))
            spi1 = counterx-nos1-1
            Zx = Z3(spi1,:)
            do countery =  nosgs+1, nosgs+nbs
                if(countery.le.nosgs+nbs2f) then
                    allocate(Zy(3))
                    spi2 = countery-nosgs
                    Zy = ZJ2f(spi2,:) 
                    H(counterx, countery) = RedElem_b(3, Zy, 6, Zx, jcJ3, jcsJ3, lJ3, 3, jc2, jcs2, l2f, 2, bbJ3J2f(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum_c(2, 0, 3, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nosgs+nbs2f+nbs2p) then
                    allocate(Zy(3))
                    spi2 = countery-nosgs-nbs2f
                    Zy = ZJ2p(spi2,:)
                    H(counterx, countery) = RedElem_b(3, Zy, 6, Zx, jcJ3, jcsJ3, lJ3, 3, jc2, jcs2, l2p, 2, bbJ3J2p(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum_c(2, 0, 3, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nosgs+nbs2f+nbs2p+nbs0) then
                    allocate(Zy(2))
                    Zy = 0.d0
                    H(counterx, countery) = Zer
                else
                    allocate(Zy(6))
                    spi2 = countery-nosgs-nbs2f-nbs2p-nbs0
                    Zy = ZJ4l(spi2,:)
                    H(counterx, countery) = RedElem_b(6, Zy, 6, Zx, jcJ3, jcsJ3, lJ3, 4, jc4, jcs4, l4l, 4, bbJ3J4l(spi1,spi2,:,:))
                !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                H(counterx, countery) = H(counterx, countery)*(-1)**(4-0) * qsum_c(4,0,3,0,evec)
                !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                endif
                deallocate(Zy)
                H(countery, counterx) = H(counterx, countery)
            end do
            deallocate(Zx)
        end do
    
        ! Then we just have to complete the coupling to the continuum. 
        ! We have to consider only the continuum states
        write(6,*) 'Combining with the continuum'
        write(6,*) 'Density in atomic units', dEcm/CMperAU
        ! Couple the intermediate states with the continuum
        do counterx=nosgs+1, nosgs+nbs
            if(counterx.le.nosgs+nbs2f) then
                allocate(Zx(3))
                spi1 = counterx-nosgs
                Zx = ZJ2f(spi1,:)
            elseif(counterx.le.nosgs+nbs2f+nbs2p) then
                allocate(Zx(3))
                spi1=counterx-nosgs-nbs2f
                Zx = Zj2p(spi1,:)
            elseif(counterx.le.nosgs+nbs2f+nbs2p+nbs0) then
                allocate(Zx(2))
                spi1 = counterx-nosgs-nbs2f-nbs2p
                Zx = ZJ0p(spi1,:)
            else
                allocate(Zx(6))
                spi1 = counterx-nosgs-nbs2f-nbs2p-nbs0
                Zx = ZJ4l(spi1,:)
            endif
            do spi2=1, ne
                Econt = Elo+spi2*dEcm
                do spi3=1, 5
                    allocate(Sdmat(5,5))
                    allocate(Kmat(5,5))
                    call J1_MQDT_at(Econt,Sdmat,Kmat,0)
                    countery = nosgs+nbs+(spi3-1)*ne+spi2
                    if(counterx.le.nosgs+nbs2f) then
                        H(counterx, countery) = CONJG(RedElem_c(5, 3, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc2,&
                                                 l2f, jcs2, 2, regJ2fJ1(spi2,spi1,:,:), iregJ2fJ1(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nosgs+nbs2f+nbs2f) then
                        H(counterx, countery) = CONJG(RedElem_c(5, 3, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc2,&
                                                 l2p, jcs2, 2, regJ2pJ1(spi2,spi1,:,:), iregJ2pJ1(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nosgs+nbs2f+nbs2p+nbs0) then
                        H(counterx, countery) = CONJG(RedElem_c(5, 2, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc0,&
                                                 l0p, jcs0, 0, regJ0pJ1(spi2,spi1,:,:), iregJ0pJ1(spi2,spi1,:,:)))
                        !write(6,*) "before the sum ", spi1, spi2, counterx, countery, H(counterx, countery) 
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,0,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    else 
                        H(counterx, countery) = Zer
                        H(countery, counterx) = Zer
                    endif
                    deallocate(Sdmat)
                    deallocate(Kmat)
                enddo
                do spi3=1, 6
                    allocate(Sdmat(6,6))
                    allocate(Kmat(6,6))
                    call J3_MQDT_at(Econt,Sdmat,Kmat,0)
                    countery = nosgs+nbs+5*ne+(spi3-1)*ne+spi2
                    !write(6,*) "Channel ",spi3, countery
                    if(counterx.le.nosgs+nbs2f) then
                        H(counterx, countery) = CONJG(RedElem_c(6, 3, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc2,&
                                                 l2f, jcs2, 2, regJ2fJ3(spi2,spi1,:,:), iregJ2fJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nosgs+nbs2f+nbs2p) then
                        H(counterx, countery) = CONJG(RedElem_c(6, 3, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc2,&
                                                 l2p, jcs2, 2, regJ2pJ3(spi2,spi1,:,:), iregJ2pJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nosgs+nbs2f+nbs2p+nbs0) then
                        H(counterx, countery) = Zer
                        H(countery, counterx) = Zer
                    else
                        H(counterx, countery) = CONJG(RedElem_c(6, 6, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc4,&
                                                 l4l, jcs4, 4, regJ4lJ3(spi2,spi1,:,:), iregJ4lJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,4,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    endif
                    deallocate(Sdmat)
                    deallocate(Kmat)
                enddo
                do spi3=1,6
                    allocate(Sdmat(6,6))
                    allocate(Kmat(6,6))
                    call J5_MQDT_at(Econt,Sdmat,Kmat,0)
                    countery = nosgs+nbs+(5+6)*ne+(spi3-1)*ne+spi2
                    if(counterx.gt.nosgs+nbs2f+nbs2p+nbs0) then
                        H(counterx, countery) = CONJG(RedElem_c(6, 6, Zx, Sdmat, spi3, jcJ5, lJ5, jcsJ5, 5, jc4,&
                                                 l4l, jcs4, 4, regJ4lJ5(spi2,spi1,:,:), iregJ4lJ5(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(5-(0)) * qsum_c(5,0,4,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    else 
                        H(counterx, countery) = Zer
                        H(countery, counterx) = Zer
                    endif
                    deallocate(Sdmat)
                    deallocate(Kmat)
                end do
            end do
            deallocate(Zx)
        end do
    end subroutine TPMatEDwGs
    !
    subroutine TPMatEDwAu(evec, Elo, dEcm,&
        bbJ1J0p, bbJ1J2f, bbJ1J2p, bbJ3J2f, bbJ3J2p, bbJ3J4l, &
        Z1, Z3, ZJ2f, ZJ2p, ZJ0p, ZJ4l, &
        EauJ2f, EauJ2p, EauJ0p,&
        regJ1J2f, iregJ1J2f, regJ1J2p, iregJ1J2p, regJ1J0p, iregJ1J0p,&
        regJ3J2f, iregJ3J2f, regJ3J2p, iregJ3J2p, &
        regJ0pJ1, iregJ0pJ1, regJ2fJ1, iregJ2fJ1, regJ2pJ1, iregJ2pJ1,&
        regJ2pJ3, iregJ2pJ3, regJ2fJ3, iregJ2fJ3, regJ4lJ3, iregJ4lJ3, &
        regJ4lJ5, iregJ4lJ5, H)
        ! Polarization vectors
        real*8, INTENT(IN) :: evec(3)
        ! Continuum energy parameters
        real*8, INTENT(IN) :: Elo, dEcm
        ! Continuum energy of autoionizing states
        real*8, INTENT(IN) :: EauJ2f(nau,2), EauJ2p(nau,2), EauJ0p(nau,2)
        ! Variable to evalute stuff
        real*8 :: Econt
        ! Z coeffs
        real*8, INTENT(IN) :: Z1(nos1,5), Z3(nos3,6), ZJ2f(nbs2f, 3), ZJ2p(nbs2p, 3), ZJ0p(nbs0, 2), &
        ZJ4l(nbs4,6)
        ! Dipoles with the gs
        real*8 :: Dwgs(4)
        ! Bound bound radial integrals.
        real*8, INTENT(IN) :: bbJ1J0p(nos1, nbs0, 2, 5), bbJ1J2f(nos1, nbs2f, 3, 5), bbJ1J2p(nos1, nbs2p, 3, 5)
        real*8, INTENT(IN) :: bbJ3J2f(nos3, nbs2f, 3, 6), bbJ3J2p(nos3, nbs2p, 3, 6), bbJ3J4l(nos3, nbs4, 6, 6)
        ! Bound Continuum radial integrals.
        ! Autoionizing
        real*8, INTENT(IN) :: regJ1J2f(nau, nos1, 3, 5), iregJ1J2f(nau, nos1, 3, 5)
        real*8, INTENT(IN) :: regJ1J2p(nau, nos1, 3, 5), iregJ1J2p(nau, nos1, 3, 5)
        real*8, INTENT(IN) :: regJ1J0p(nau, nos1, 2, 5), iregJ1J0p(nau, nos1, 2, 5)
        real*8, INTENT(IN) :: regJ3J2f(nau, nos3, 3, 6), iregJ3J2f(nau, nos3, 3, 6)
        real*8, INTENT(IN) :: regJ3J2p(nau, nos3, 3, 6), iregJ3J2p(nau, nos3, 3, 6)
        ! True continuum
        real*8, INTENT(IN) :: regJ0pJ1(ne, nbs0, 5, 2), iregJ0pJ1(ne, nbs0, 5, 2)
        real*8, INTENT(IN) :: regJ2fJ3(ne,nbs2f, 6, 3), iregJ2fJ3(ne, nbs2f, 6, 3), &
        regJ2fJ1(ne, nbs2f, 5, 3), iregJ2fJ1(ne, nbs2f, 5, 3)
        real*8, INTENT(IN) :: regJ2pJ3(ne,nbs2p, 6, 3), iregJ2pJ3(ne, nbs2p, 6, 3), &
        regJ2pJ1(ne, nbs2p, 5, 3), iregJ2pJ1(ne, nbs2p, 5, 3)
        real*8, INTENT(IN) :: regJ4lJ3(ne, nbs4, 6, 6), iregJ4lJ3(ne, nbs4, 6, 6), &
        regJ4lJ5(ne, nbs4, 6, 6), iregJ4lJ5(ne, nbs4, 6, 6)
        ! Angular quantum numbers 
        integer :: l2f(3), l2p(3), l0p(2), l4l(6), lJ1(5), lJ3(6), lJ5(6)
        integer :: jcs2(3), jcs0(2), jcs4(6), jcsJ1(5), jcsJ3(6), jcsJ5(6)
        real :: jc2(3), jc0(2), jc4(6), jcJ1(5), jcJ3(6), jcJ5(6) 
        ! allocatable arrays
        real*8, allocatable :: Zx(:), Zy(:)
        complex*16, allocatable :: Sdmat(:,:)
        real*8, allocatable :: Kmat(:,:)
        complex*16, allocatable :: Sdp(:,:), Zco(:,:)
        integer, allocatable :: iopen(:), iclosed(:)
        real*8, allocatable :: beta(:)
        ! Hamiltonian
        complex*16, intent(out) :: H(ntotau+1, ntotau+1)
        integer :: counterx, countery
        counterx = 0
        countery = 0
    
        l2f = (/ 3, 3, 3 /)
        l2p = (/ 1, 1, 1/)
        jc2 = (/ 3./2., 3./2., 1./2./)
        jcs2 =  (/ 2, 1, 1 /)
        
        l0p = (/1, 1/)
        jc0 = (/3/2., 1/2./)
        jcs0 = (/1, 1/)
    
        l4l = (/3, 3, 3, 5, 5, 5/)
        jc4 = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        jcs4 = (/1, 2, 1, 1, 2, 1/)
    
        lJ1 = lb
        jcJ1 = jcb
        jcsJ1 = jcsb
    
        lJ3 = (/ 2, 2, 2, 4, 4, 4/)
        jcJ3 = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        jcsJ3 = (/ 1, 2, 1, 1, 2, 1/)
    
        lJ5 = (/4,4,4,6,6,6/)
        jcJ5 = (/3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        jcsJ5 = (/1, 2, 1, 1, 2, 1/)
    
        H = Zer
        Dwgs = (/0.278603d0, 0.518247d0, 0.18458d0, 0.54917d0/)
        ! Fill the coupling between the ground state and the J=1 bound states
        do countery=1,nos1
            H(1,countery+1) = Dwgs(countery)
            H(countery+1,1) = H(1,countery+1)
        end do
        
    
        ! The matrix is really sparse, so we start filling by parts.
        ! Start with the coupling between the initial states to the intermediates:
        do counterx=2,nos1+1
            allocate(Zx(5))
            spi1 = counterx-1
            Zx = Z1(spi1,:)
            do countery =  nosgs+1, nosgs+nbs
                if(countery.le.nosgs+nbs2f) then
                    allocate(Zy(3))
                    spi2 = countery-nosgs
                    Zy = ZJ2f(spi2,:) 
                    H(counterx, countery) = RedElem_b(3, Zy, 5, Zx, jcJ1, jcsJ1, lJ1, 1, jc2, jcs2, l2f, 2, bbJ1J2f(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum(2, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nosgs+nbs2f+nbs2p) then
                    allocate(Zy(3))
                    spi2 = countery-nosgs-nbs2f
                    Zy = ZJ2p(spi2,:)
                    H(counterx, countery) = RedElem_b(3, Zy, 5, Zx, jcJ1, jcsJ1, lJ1, 1, jc2, jcs2, l2p, 2, bbJ1J2p(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum(2, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nosgs+nbs2p+nbs2f+nbs0) then
                    allocate(Zy(2))
                    spi2 = countery-nosgs-nbs2p-nbs2f
                    Zy = ZJ0p(spi2,:)
                    H(counterx, countery) = RedElem_b(2, Zy, 5, Zx, jcJ1, jcsJ1, lJ1, 1, jc0, jcs0, l0p, 0, bbJ1J0p(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(0-0) * qsum(0,0,evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                else
                    allocate(Zy(6))
                    Zy = 0.d0
                    H(counterx, countery) = Zer
                endif
                deallocate(Zy)
                H(countery, counterx) = H(counterx, countery)
            end do
            deallocate(Zx)
        end do
        ! Now the initial J=3 states
        do counterx=nos1+2,nosgs
            allocate(Zx(6))
            spi1 = counterx-nos1-1
            Zx = Z3(spi1,:)
            do countery =  nosgs+1, nosgs+nbs
                if(countery.le.nosgs+nbs2f) then
                    allocate(Zy(3))
                    spi2 = countery-nosgs
                    Zy = ZJ2f(spi2,:) 
                    H(counterx, countery) = RedElem_b(3, Zy, 6, Zx, jcJ3, jcsJ3, lJ3, 3, jc2, jcs2, l2f, 2, bbJ3J2f(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum_c(2, 0, 3, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nosgs+nbs2f+nbs2p) then
                    allocate(Zy(3))
                    spi2 = countery-nosgs-nbs2f
                    Zy = ZJ2p(spi2,:)
                    H(counterx, countery) = RedElem_b(3, Zy, 6, Zx, jcJ3, jcsJ3, lJ3, 3, jc2, jcs2, l2p, 2, bbJ3J2p(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum_c(2, 0, 3, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nosgs+nbs2p+nbs2f+nbs0) then
                    allocate(Zy(2))
                    Zy = 0.d0
                    H(counterx, countery) = Zer
                else
                    allocate(Zy(6))
                    spi2 = countery-nosgs-nbs2p-nbs2f-nbs0
                    Zy = ZJ4l(spi2,:)
                    H(counterx, countery) = RedElem_b(6, Zy, 6, Zx, jcJ3, jcsJ3, lJ3, 4, jc4, jcs4, l4l, 4, bbJ3J4l(spi1,spi2,:,:))
                !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                H(counterx, countery) = H(counterx, countery)*(-1)**(4-0) * qsum_c(4,0,3,0,evec)
                !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                endif
                deallocate(Zy)
                H(countery, counterx) = H(counterx, countery)
            end do
            deallocate(Zx)
        end do
        
        ! Now we add the coupling with the autoionizing states
        ! J=1
        do counterx=2,nos1+1
            allocate(Zx(5))
            spi1 = counterx-1
            Zx = Z1(spi1,:)
            !
            do spi2=1,2
                do spi3=1,nau
                countery = nosgs+nbs+(spi2-1)*nau+spi3
                call J2f_MQDT_bet(EauJ2f(spi3,1), Sdp,Zco,iopen,iclosed,beta)
                H(counterx, countery) = CONJG(RedElem_bet_closed(2, 1, 5, Zx, Sdp, Zco, beta,spi2, iopen, iclosed, &
                jc2, l2f, jcs2, 2, jcJ1, lJ1, jcsJ1, 1, regJ1J2f(spi3,spi1,:,:), iregJ1J2f(spi3,spi1,:,:)))
                H(counterx, countery) = SQRT(1.d0)*H(counterx, countery) * (-1)**(2-0) * qsum(2,0,evec)
                H(countery, counterx) = CONJG(H(counterx,countery))
                !write(6, *) counterx, countery, H(counterx,countery)
                deallocate(Sdp,Zco,iopen,iclosed,beta)
                end do
            end do
            !
            do spi2=1,2
                do spi3=1,nau
                countery = nosgs+nbs+(spi2-1)*nau+spi3+2*nau
                call J2p_MQDT_bet(EauJ2p(spi3,1), Sdp,Zco,iopen,iclosed,beta)
                H(counterx, countery) = CONJG(RedElem_bet_closed(2, 1, 5, Zx, Sdp, Zco, beta,spi2, iopen, iclosed, &
                jc2, l2p, jcs2, 2, jcJ1, lJ1, jcsJ1, 1, regJ1J2p(spi3,spi1,:,:), iregJ1J2p(spi3,spi1,:,:)))
                H(counterx, countery) = SQRT(1.d0)*H(counterx, countery) * (-1)**(2-0) * qsum(2,0,evec)
                H(countery, counterx) = CONJG(H(counterx,countery))
                deallocate(Sdp,Zco,iopen,iclosed,beta)
                end do
            end do
            !
            do spi2=1,1
                do spi3=1,nau
                countery = nosgs+nbs+(spi2-1)*nau+spi3+4*nau
                call J0_MQDT_bet(EauJ0p(spi3,1), Sdp,Zco,iopen,iclosed,beta)
                H(counterx, countery) = CONJG(RedElem_bet_closed(1, 1, 5, Zx, Sdp, Zco, beta,spi2, iopen, iclosed, &
                jc0, l0p, jcs0, 0, jcJ1, lJ1, jcsJ1, 1, regJ1J0p(spi3,spi1,:,:), iregJ1J0p(spi3,spi1,:,:)))
                H(counterx, countery) = SQRT(1.d0)*H(counterx, countery) * (-1)**(0-0) * qsum(0,0,evec)
                H(countery, counterx) = CONJG(H(counterx,countery))
                end do
            end do
            DEALLOCATE(Zx)
        end do
        ! J=3
        do counterx=nos1+2, nosgs
            allocate(Zx(6))
            spi1 = counterx-1-nos1
            Zx = Z3(spi1,:)
            !
            do spi2=1,2
                do spi3=1,nau
                countery = nosgs+nbs+(spi2-1)*nau+spi3
                call J2f_MQDT_bet(EauJ2f(spi3,1), Sdp,Zco,iopen,iclosed,beta)
                H(counterx, countery) = CONJG(RedElem_bet_closed(2, 1, 6, Zx, Sdp, Zco, beta,spi2, iopen, iclosed, &
                jc2, l2f, jcs2, 2, jcJ3, lJ3, jcsJ3, 3, regJ3J2f(spi3,spi1,:,:), iregJ3J2f(spi3,spi1,:,:)))
                H(counterx, countery) = SQRT(1.d0)*H(counterx, countery) * (-1)**(3-0) * qsum_c(2,0,3,0,evec)
                H(countery, counterx) = CONJG(H(counterx,countery))
                deallocate(Sdp,Zco,iopen,iclosed,beta)
                end do
            end do
            !
            do spi2=1,2
                do spi3=1,nau
                countery = nosgs+nbs+(spi2-1)*nau+spi3+2*nau
                call J2p_MQDT_bet(EauJ2p(spi3,1), Sdp,Zco,iopen,iclosed,beta)
                H(counterx, countery) = CONJG(RedElem_bet_closed(2, 1, 6, Zx, Sdp, Zco, beta,spi2, iopen, iclosed, &
                jc2, l2p, jcs2, 2, jcJ3, lJ3, jcsJ3, 3, regJ3J2p(spi3,spi1,:,:), iregJ3J2p(spi3,spi1,:,:)))
                H(counterx, countery) = SQRT(1.d0)*H(counterx, countery) * (-1)**(2-0) * qsum_c(2,0,3,0,evec)
                H(countery, counterx) = CONJG(H(counterx,countery))
                deallocate(Sdp,Zco,iopen,iclosed,beta)
                end do
            end do
            !
            DEALLOCATE(Zx)
        end do
        ! Then we just have to complete the coupling to the continuum. 
        ! We have to consider only the continuum states
        write(6,*) 'Combining with the continuum'
        write(6,*) 'Density in atomic units', dEcm/CMperAU
    
        do counterx=nosau+1, nosau+nbs
            if(counterx.le.nosau+nbs2f) then
                allocate(Zx(3))
                spi1 = counterx-nosau
                Zx = ZJ2f(spi1,:)
            elseif(counterx.le.nosau+nbs2f+nbs2p) then
                allocate(Zx(3))
                spi1=counterx-nosau-nbs2f
                Zx = Zj2p(spi1,:)
            elseif(counterx.le.nosau+nbs2f+nbs2p+nbs0) then
                allocate(Zx(2))
                spi1 = counterx-nosau-nbs2f-nbs2p
                Zx = ZJ0p(spi1,:)
                write(6,*) "Getting the 0 Zx @", counterx, spi1 
            else
                allocate(Zx(6))
                spi1 = counterx-nosau-nbs2f-nbs2p-nbs0
                Zx = ZJ4l(spi1,:)
            endif
            do spi2=1, ne
                Econt = Elo+spi2*dEcm
                do spi3=1, 5
                    allocate(Sdmat(5,5))
                    allocate(Kmat(5,5))
                    call J1_MQDT_at(Econt,Sdmat,Kmat,0)
                    countery = nosau+nbs+5*nau+(spi3-1)*ne+spi2
                    if(counterx.le.nosau+nbs2f) then
                        H(counterx, countery) = CONJG(RedElem_c(5, 3, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc2,&
                                                 l2f, jcs2, 2, regJ2fJ1(spi2,spi1,:,:), iregJ2fJ1(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nosau+nbs2f+nbs2p) then
                        H(counterx, countery) = CONJG(RedElem_c(5, 3, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc2,&
                                                 l2p, jcs2, 2, regJ2pJ1(spi2,spi1,:,:), iregJ2pJ1(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nosau+nbs2f+nbs2p+nbs0) then
                        H(counterx, countery) = CONJG(RedElem_c(5, 2, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc0,&
                                                 l0p, jcs0, 0, regJ0pJ1(spi2,spi1,:,:), iregJ0pJ1(spi2,spi1,:,:)))
                        !write(6,*) "before the sum ", spi1, spi2, counterx, countery, H(counterx, countery) 
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,0,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    else 
                        H(counterx, countery) = Zer
                        H(countery, counterx) = Zer
                    endif
                    deallocate(Sdmat)
                    deallocate(Kmat)
                enddo
                do spi3=1, 6
                    allocate(Sdmat(6,6))
                    allocate(Kmat(6,6))
                    call J3_MQDT_at(Econt,Sdmat,Kmat,0)
                    countery = nosau+nbs+5*nau+5*ne+(spi3-1)*ne+spi2
                    !write(6,*) "Channel ",spi3, countery
                    if(counterx.le.nosau+nbs2f) then
                        H(counterx, countery) = CONJG(RedElem_c(6, 3, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc2,&
                                                 l2f, jcs2, 2, regJ2fJ3(spi2,spi1,:,:), iregJ2fJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nosau+nbs2f+nbs2p) then
                        H(counterx, countery) = CONJG(RedElem_c(6, 3, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc2,&
                                                 l2p, jcs2, 2, regJ2pJ3(spi2,spi1,:,:), iregJ2pJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nosau+nbs2f+nbs2p+nbs0) then
                        H(counterx, countery) = Zer
                        H(countery, counterx) = Zer
                    else
                        H(counterx, countery) = CONJG(RedElem_c(6, 6, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc4,&
                                                 l4l, jcs4, 4, regJ4lJ3(spi2,spi1,:,:), iregJ4lJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,4,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    endif
                    deallocate(Sdmat)
                    deallocate(Kmat)
                enddo
                do spi3=1,6
                    allocate(Sdmat(6,6))
                    allocate(Kmat(6,6))
                    call J5_MQDT_at(Econt,Sdmat,Kmat,0)
                    countery = nosau+nbs+5*nau+(5+6)*ne+(spi3-1)*ne+spi2
                    if(counterx.gt.nosau+nbs2f+nbs2p+nbs0) then
                        H(counterx, countery) = CONJG(RedElem_c(6, 6, Zx, Sdmat, spi3, jcJ5, lJ5, jcsJ5, 5, jc4,&
                                                 l4l, jcs4, 4, regJ4lJ5(spi2,spi1,:,:), iregJ4lJ5(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(5-(0)) * qsum_c(5,0,4,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    else 
                        H(counterx, countery) = Zer
                        H(countery, counterx) = Zer
                    endif
                    deallocate(Sdmat)
                    deallocate(Kmat)
                end do
            end do
            deallocate(Zx)
        end do
    end subroutine TPMatEDwAu
    !
    subroutine TPMatNED(evec, Elo, dEcm,&
        bbJ1J0p, bbJ1J2f, bbJ1J2p, bbJ3J2f, bbJ3J2p, bbJ3J4l, &
        Z1, Z3, ZJ2f, ZJ2p, ZJ0p, ZJ4l, &
        regJ0pJ1, iregJ0pJ1, regJ2fJ1, iregJ2fJ1, regJ2pJ1, iregJ2pJ1,&
        regJ2pJ3, iregJ2pJ3, regJ2fJ3, iregJ2fJ3, regJ4lJ3, iregJ4lJ3, &
        regJ4lJ5, iregJ4lJ5, H)
        ! Polarization vectors
        real*8, INTENT(IN) :: evec(3)
        ! Continuum energy parameters
        real*8, INTENT(IN) :: Elo, dEcm
        real*8 :: Econt
        ! Z coeffs
        real*8, INTENT(IN) :: Z1(nos1,5), Z3(nos3,6), ZJ2f(nbs2f, 3), ZJ2p(nbs2p, 3), ZJ0p(nbs0, 2), &
        ZJ4l(nbs4,6)
        ! Bound bound radial integrals.
        real*8, INTENT(IN) :: bbJ1J0p(nos1, nbs0, 2, 5), bbJ1J2f(nos1, nbs2f, 3, 5), bbJ1J2p(nos1, nbs2p, 3, 5)
        real*8, INTENT(IN) :: bbJ3J2f(nos3, nbs2f, 3, 6), bbJ3J2p(nos3, nbs2p, 3, 6), bbJ3J4l(nos3, nbs4, 6, 6)
        ! Bound Continuum radial integrals.
        real*8, INTENT(IN) :: regJ0pJ1(ne, nbs0, 5, 2), iregJ0pJ1(ne, nbs0, 5, 2)
        real*8, INTENT(IN) :: regJ2fJ3(ne,nbs2f, 6, 3), iregJ2fJ3(ne, nbs2f, 6, 3), &
        regJ2fJ1(ne, nbs2f, 5, 3), iregJ2fJ1(ne, nbs2f, 5, 3)
        real*8, INTENT(IN) :: regJ2pJ3(ne,nbs2p, 6, 3), iregJ2pJ3(ne, nbs2p, 6, 3), &
        regJ2pJ1(ne, nbs2p, 5, 3), iregJ2pJ1(ne, nbs2p, 5, 3)
        real*8, INTENT(IN) :: regJ4lJ3(ne, nbs4, 6, 6), iregJ4lJ3(ne, nbs4, 6, 6), &
        regJ4lJ5(ne, nbs4, 6, 6), iregJ4lJ5(ne, nbs4, 6, 6)
        ! Angular quantum numbers 
        integer :: l2f(3), l2p(3), l0p(2), l4l(6), lJ1(5), lJ3(6), lJ5(6)
        integer :: jcs2(3), jcs0(2), jcs4(6), jcsJ1(5), jcsJ3(6), jcsJ5(6)
        real :: jc2(3), jc0(2), jc4(6), jcJ1(5), jcJ3(6), jcJ5(6) 
        real*8, allocatable :: Zx(:), Zy(:)
        complex*16, allocatable :: Sdmat(:,:)
        real*8, allocatable :: Kmat(:,:)
        complex*16, intent(out) :: H(ntot, ntot)
        integer :: counterx, countery
        counterx = 0
        countery = 0

        l2f = (/ 3, 3, 3 /)
        l2p = (/ 1, 1, 1/)
        jc2 = (/ 3./2., 3./2., 1./2./)
        jcs2 =  (/ 2, 1, 1 /)
        
        l0p = (/1, 1/)
        jc0 = (/3/2., 1/2./)
        jcs0 = (/1, 1/)

        l4l = (/3, 3, 3, 5, 5, 5/)
        jc4 = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        jcs4 = (/1, 2, 1, 1, 2, 1/)

        lJ1 = lb
        jcJ1 = jcb
        jcsJ1 = jcsb

        lJ3 = (/ 2, 2, 2, 4, 4, 4/)
        jcJ3 = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        jcsJ3 = (/ 1, 2, 1, 1, 2, 1/)

        lJ5 = (/4,4,4,6,6,6/)
        jcJ5 = (/3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        jcsJ5 = (/1, 2, 1, 1, 2, 1/)

        H = Zer

        ! The matrix is really sparse, so we start filling by parts.
        ! Start with the coupling between the initial states to the intermediates:
        do counterx=1,nos1
            allocate(Zx(5))
            Zx = Z1(counterx,:)
            do countery =  nos+1, nos+nbs
                if(countery.le.nos+nbs2f) then
                    allocate(Zy(3))
                    spi1 = counterx
                    spi2 = countery-nos
                    Zy = ZJ2f(spi2,:) 
                    H(counterx, countery) = RedElem_b(3, Zy, 5, Zx, jcJ1, jcsJ1, lJ1, 1, jc2, jcs2, l2f, 2, bbJ1J2f(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum(2, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nos+nbs2f+nbs2p) then
                    allocate(Zy(3))
                    spi1 = counterx
                    spi2 = countery-nos-nbs2f
                    Zy = ZJ2p(spi2,:)
                    H(counterx, countery) = RedElem_b(3, Zy, 5, Zx, jcJ1, jcsJ1, lJ1, 1, jc2, jcs2, l2p, 2, bbJ1J2p(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum(2, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nos+nbs2f+nbs2p+nbs0) then
                    allocate(Zy(2))
                    spi1 = counterx
                    spi2 = countery-nos-nbs2f-nbs2p
                    Zy = ZJ0p(spi2,:)
                    H(counterx, countery) = RedElem_b(2, Zy, 5, Zx, jcJ1, jcsJ1, lJ1, 1, jc0, jcs0, l0p, 0, bbJ1J0p(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(0-0) * qsum(0,0,evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                else
                    allocate(Zy(6))
                    Zy = 0.d0
                    H(counterx, countery) = Zer
                endif
                deallocate(Zy)
                H(countery, counterx) = H(counterx, countery)
            end do
            deallocate(Zx)
        end do
        ! Now the initial J=3 states
        do counterx=nos1+1,nos
            allocate(Zx(6))
            spi1 = counterx-nos1
            Zx = Z3(spi1,:)
            do countery =  nos+1, nos+nbs
                if(countery.le.nos+nbs2f) then
                    allocate(Zy(3))
                    spi2 = countery-nos
                    Zy = ZJ2f(spi2,:) 
                    H(counterx, countery) = RedElem_b(3, Zy, 6, Zx, jcJ3, jcsJ3, lJ3, 3, jc2, jcs2, l2f, 2, bbJ3J2f(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum_c(2, 0, 3, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nos+nbs2f+nbs2p) then
                    allocate(Zy(3))
                    spi1 = counterx
                    spi2 = countery-nos-nbs2f
                    Zy = ZJ2p(spi2,:)
                    H(counterx, countery) = RedElem_b(3, Zy, 6, Zx, jcJ3, jcsJ3, lJ3, 3, jc2, jcs2, l2p, 2, bbJ3J2p(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum_c(2, 0, 3, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nos+nbs2f+nbs2p+nbs0) then
                    allocate(Zy(2))
                    Zy = 0.d0
                    H(counterx, countery) = Zer
                else
                    allocate(Zy(6))
                    spi1 = counterx
                    spi2 = countery-nos-nbs2f-nbs2p-nbs0
                    Zy = ZJ4l(spi2,:)
                    H(counterx, countery) = RedElem_b(6, Zy, 6, Zx, jcJ3, jcsJ3, lJ3, 4, jc4, jcs4, l4l, 4, bbJ3J4l(spi1,spi2,:,:))
                !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                H(counterx, countery) = H(counterx, countery)*(-1)**(4-0) * qsum_c(4,0,3,0,evec)
                !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                endif
                deallocate(Zy)
                H(countery, counterx) = H(counterx, countery)
            end do
            deallocate(Zx)
        end do

        ! Then we just have to complete the coupling to the continuum. 
        ! We have to consider only the continuum states
        write(6,*) 'Combining with the continuum'
        write(6,*) 'Density in atomic units', dEcm/CMperAU
        ! J2f states
        do counterx=nos+1, nos+nbs
            if(counterx.le.nos+nbs2f) then
                allocate(Zx(3))
                spi1 = counterx-nos
                Zx = ZJ2f(spi1,:)
            elseif(counterx.le.nos+nbs2f+nbs2p) then
                allocate(Zx(3))
                spi1=counterx-nos-nbs2f
                Zx = Zj2p(spi1,:)
            elseif(counterx.le.nos+nbs2f+nbs2p+nbs0) then
                allocate(Zx(2))
                spi1 = counterx-nos-nbs2f-nbs2p
                Zx = ZJ0p(spi1,:)
            else
                allocate(Zx(6))
                spi1 = counterx-nos-nbs2f-nbs2p-nbs0
                Zx = ZJ4l(spi1,:)
            endif
            do spi2=1, ne
                Econt = Elo+spi2*dEcm
                do spi3=1, 5
                    allocate(Sdmat(5,5))
                    allocate(Kmat(5,5))
                    call J1_MQDT_at(Econt,Sdmat,Kmat,0)
                    countery = nos+nbs+(spi3-1)*ne+spi2
                    if(counterx.le.nos+nbs2f) then
                        H(counterx, countery) = CONJG(RedElem_c(5, 3, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc2,&
                                                 l2f, jcs2, 2, regJ2fJ1(spi2,spi1,:,:), iregJ2fJ1(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(1.0d0)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nos+nbs2f+nbs2p) then
                        H(counterx, countery) = CONJG(RedElem_c(5, 3, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc2,&
                                                 l2p, jcs2, 2, regJ2pJ1(spi2,spi1,:,:), iregJ2pJ1(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(1.0d0)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nos+nbs2f+nbs2p+nbs0) then
                        H(counterx, countery) = CONJG(RedElem_c(5, 2, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc0,&
                                                 l0p, jcs0, 0, regJ0pJ1(spi2,spi1,:,:), iregJ0pJ1(spi2,spi1,:,:)))
                        !write(6,*) "before the sum ", spi1, spi2, counterx, countery, H(counterx, countery) 
                        H(counterx, countery) = SQRT(1.0d0)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,0,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    else 
                        H(counterx, countery) = Zer
                        H(countery, counterx) = Zer
                    endif
                    deallocate(Sdmat)
                    deallocate(Kmat)
                enddo
                do spi3=1, 6
                    allocate(Sdmat(6,6))
                    allocate(Kmat(6,6))
                    call J3_MQDT_at(Econt,Sdmat,Kmat,0)
                    countery = nos+nbs+5*ne+(spi3-1)*ne+spi2
                    !write(6,*) "Channel ",spi3, countery
                    if(counterx.le.nos+nbs2f) then
                        H(counterx, countery) = CONJG(RedElem_c(6, 3, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc2,&
                                                 l2f, jcs2, 2, regJ2fJ3(spi2,spi1,:,:), iregJ2fJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(1.0d0)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nos+nbs2f+nbs2p) then
                        H(counterx, countery) = CONJG(RedElem_c(6, 3, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc2,&
                                                 l2p, jcs2, 2, regJ2pJ3(spi2,spi1,:,:), iregJ2pJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(1.0d0)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nos+nbs2f+nbs2p+nbs0) then
                        H(counterx, countery) = Zer
                        H(countery, counterx) = Zer
                    else
                        H(counterx, countery) = CONJG(RedElem_c(6, 6, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc4,&
                                                 l4l, jcs4, 4, regJ4lJ3(spi2,spi1,:,:), iregJ4lJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(1.0d0)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,4,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    endif
                    deallocate(Sdmat)
                    deallocate(Kmat)
                enddo
                do spi3=1,6
                    allocate(Sdmat(6,6))
                    allocate(Kmat(6,6))
                    call J5_MQDT_at(Econt,Sdmat,Kmat,0)
                    countery = nos+nbs+(5+6)*ne+(spi3-1)*ne+spi2
                    if(counterx.gt.nos+nbs2f+nbs2p+nbs0) then
                        H(counterx, countery) = CONJG(RedElem_c(6, 6, Zx, Sdmat, spi3, jcJ5, lJ5, jcsJ5, 5, jc4,&
                                                 l4l, jcs4, 4, regJ4lJ5(spi2,spi1,:,:), iregJ4lJ5(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(1.0d0)*H(counterx, countery)*(-1)**(5-(0)) * qsum_c(5,0,4,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    else 
                        H(counterx, countery) = Zer
                        H(countery, counterx) = Zer
                    endif
                    deallocate(Sdmat)
                    deallocate(Kmat)
                end do
            end do
            deallocate(Zx)
        end do
        return 
        !1985 format(2500(E15.8,1x))
    end subroutine TPMatNED
    !
    subroutine TPMatED_J1(evec, Elo, dEcm,&
        bbJ2f, bbJ2p, bbJ0p, Z1, ZJ2f, ZJ2p, ZJ0p, &
        regJ2fJ3, iregJ2fJ3, regJ2fJ1, iregJ2fJ1, regJ2pJ3, iregJ2pJ3, regJ2pJ1, iregJ2pJ1,&
        regJ0pJ3, iregJ0pJ3, regJ0pJ1, iregJ0pJ1, H)
        real*8, INTENT(IN) :: evec(3)
        real*8, INTENT(IN) :: Elo, dEcm
        real*8, INTENT(IN) :: Z1(nos1,5), ZJ2f(nbs2f, 3), Zj2p(nbs2p, 3), ZJ0p(nbs0, 2)
        real*8, INTENT(IN) :: bbJ2f(nos1, nbs2f, 3, 5), bbJ2p(nos1, nbs2p, 3, 5), bbJ0p(nos1, nbs0, 2, 5)
        real*8, INTENT(IN) :: regJ2fJ3(ne,nbs2f, 6, 3), iregJ2fJ3(ne, nbs2f, 6, 3), &
        regJ2fJ1(ne, nbs2f, 5, 3), iregJ2fJ1(ne, nbs2f, 5, 3)
        real*8, INTENT(IN) :: regJ2pJ3(ne,nbs2p, 6, 3), iregJ2pJ3(ne, nbs2p, 6, 3), &
        regJ2pJ1(ne, nbs2p, 5, 3), iregJ2pJ1(ne, nbs2p, 5, 3)
        real*8, INTENT(IN) :: regJ0pJ3(ne,nbs0,6, 2), iregJ0pJ3(ne, nbs0,6, 2), &
        regJ0pJ1(ne, nbs0,5, 2), iregJ0pJ1(ne, nbs0,5, 2)
        real*8 :: Econt
        integer :: l2f(3), l2p(3), l0p(2)
        integer :: lJ1(5), lJ3(6), jcs2(3), jcs0(2), jcsJ1(5), jcsJ3(6)
        real :: jc2(3), jc0(2), jcJ1(5), jcJ3(6) 
        real*8, allocatable :: Zx(:), Zy(:)
        complex*16, allocatable :: Sdmat(:,:)
        real*8, allocatable :: Kmat(:,:)
        complex*16, intent(out) :: H(ntot1, ntot1)
        integer :: counterx, countery
        counterx = 0
        countery = 0
        l2f = (/ 3, 3, 3 /)
        l2p = (/ 1, 1, 1/)
        jc2 = (/ 3./2., 3./2., 1./2./)
        jcs2 =  (/ 2, 1, 1 /)
        l0p = (/1, 1/)
        jc0 = (/3/2., 1/2./)
        jcs0 = (/1, 1/)
        jcJ1 = jcb
        lJ1 = lb
        jcsJ1 = jcsb
        jcJ3 = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        lJ3 =  (/ 2, 2, 2, 4, 4, 4/)
        jcsJ3 = (/ 1, 2, 1, 1, 2, 1/)
        counterx = 0
        countery = 0
        H = Zer

        ! The matrix is really sparse, so we start filling by parts.
        ! Start with the coupling between the initial states to the intermediates:
        do counterx=1,nos1
            allocate(Zx(5))
            Zx = Z1(counterx,:)
            do countery =  nos1+1, nos1+nbs1
                if(countery.le.nos1+nbs2f) then
                    allocate(Zy(3))
                    spi1 = counterx
                    spi2 = countery-nos1
                    Zy = ZJ2f(spi2,:) 
                    H(counterx, countery) = RedElem_b(3, Zy, 5, Zx, jcb, jcsb, lb, 1, jc2, jcs2, l2f, 2, bbJ2f(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum(2, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nos1+nbs2f+nbs2p) then
                    allocate(Zy(3))
                    spi1 = counterx
                    spi2 = countery-nos1-nbs2f
                    Zy = ZJ2p(spi2,:)
                    H(counterx, countery) = RedElem_b(3, Zy, 5, Zx, jcb, jcsb, lb, 1, jc2, jcs2, l2p, 2, bbJ2p(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum(2, 0, evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                else
                    allocate(Zy(2))
                    spi1 = counterx
                    spi2 = countery-nos1-nbs2f-nbs2p
                    Zy = ZJ0p(spi2,:)
                    H(counterx, countery) = RedElem_b(2, Zy, 5, Zx, jcb, jcsb, lb, 1, jc0, jcs0, l0p, 0, bbJ0p(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(0-0) * qsum(0,0,evec)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                endif
                deallocate(Zy)
                H(countery, counterx) = H(counterx, countery)
            end do
            deallocate(Zx)
        end do

        ! Then we just have to complete the coupling to the continuum. 
        ! We have to consider only the continuum states
        write(6,*) 'Combining with the continuum'
        write(6,*) 'Density in atomic units', dEcm/CMperAU
        ! J2f states
        do counterx=nos1+1, nos1+nbs1
            if(counterx.le.nos1+nbs2f) then
                allocate(Zx(3))
                spi1 = counterx-nos1
                Zx = ZJ2f(spi1,:)
            elseif(counterx.le.nos1+nbs2f+nbs2p) then
                allocate(Zx(3))
                spi1=counterx-nos1-nbs2f
                Zx = Zj2p(spi1,:)
            else
                allocate(Zx(2))
                spi1 = counterx-nos1-nbs2f-nbs2p
                Zx = ZJ0p(spi1,:)
            endif
            do spi2=1, ne
                Econt = Elo+spi2*dEcm
                do spi3=1, 5
                    allocate(Sdmat(5,5))
                    allocate(Kmat(5,5))
                    call J1_MQDT_at(Econt,Sdmat,Kmat,0)
                    countery = nos1+nbs1+(spi3-1)*ne+spi2
                    if(counterx.le.nos1+nbs2f) then
                        H(counterx, countery) = CONJG(RedElem_c(5, 3, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc2,&
                                                 l2f, jcs2, 2, regJ2fJ1(spi2,spi1,:,:), iregJ2fJ1(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nos1+nbs2f+nbs2p) then
                        H(counterx, countery) = CONJG(RedElem_c(5, 3, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc2,&
                                                 l2p, jcs2, 2, regJ2pJ1(spi2,spi1,:,:), iregJ2pJ1(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    else
                        H(counterx, countery) = CONJG(RedElem_c(5, 2, Zx, Sdmat, spi3, jcJ1, lJ1, jcsJ1, 1, jc0,&
                                                 l0p, jcs0, 0, regJ0pJ1(spi2,spi1,:,:), iregJ0pJ1(spi2,spi1,:,:)))
                        !write(6,*) "before the sum ", spi1, spi2, counterx, countery, H(counterx, countery) 
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(1-(0)) * qsum_c(1,0,0,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    endif
                    deallocate(Sdmat)
                    deallocate(Kmat)
                enddo
                do spi3=1, 6
                    allocate(Sdmat(6,6))
                    allocate(Kmat(6,6))
                    call J3_MQDT_at(Econt,Sdmat,Kmat,0)
                    countery = nos1+nbs1+5*ne+(spi3-1)*ne+spi2
                    !write(6,*) "Channel ",spi3, countery
                    if(counterx.le.nos1+nbs2f) then
                        H(counterx, countery) = CONJG(RedElem_c(6, 3, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc2,&
                                                 l2f, jcs2, 2, regJ2fJ3(spi2,spi1,:,:), iregJ2fJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    elseif(counterx.le.nos1+nbs2f+nbs2p) then
                        H(counterx, countery) = CONJG(RedElem_c(6, 3, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc2,&
                                                 l2p, jcs2, 2, regJ2pJ3(spi2,spi1,:,:), iregJ2pJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,2,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    else
                        H(counterx, countery) = CONJG(RedElem_c(6, 2, Zx, Sdmat, spi3, jcJ3, lJ3, jcsJ3, 3, jc0,&
                                                 l0p, jcs0, 0, regJ0pJ3(spi2,spi1,:,:), iregJ0pJ3(spi2,spi1,:,:)))
                        H(counterx, countery) = SQRT(dEcm/CMperAU)*H(counterx, countery)*(-1)**(3-(0)) * qsum_c(3,0,0,0,evec)
                        H(countery, counterx) = CONJG(H(counterx, countery))
                        !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                    endif
                    deallocate(Sdmat)
                    deallocate(Kmat)
                enddo
            end do
            deallocate(Zx)
        end do
        return 
        !1985 format(2500(E15.8,1x))
    end subroutine TPMatED_J1
    !
    subroutine TPmatNC(evec, bbJ2f1, bbJ2p1, bbJ0p1, bbJ2f3, bbJ2p3, bbJ4l3, ZJ1, ZJ3, ZJ2f, ZJ2p, ZJ0p, ZJ4, H)
        real*8, INTENT(IN) :: evec(3)
        real*8, INTENT(IN) :: ZJ1(nos1,5), ZJ3(nos3, 6), ZJ2f(nbs2f, 3), Zj2p(nbs2p, 3), ZJ0p(nbs0, 2), ZJ4(nbs4,6)
        real*8, INTENT(IN) :: bbJ2f1(nos1, nbs2f, 3, 5), bbJ2p1(nos1, nbs2p, 3, 5), bbJ0p1(nos1, nbs0, 2, 5)
        real*8, INTENT(IN) :: bbJ2f3(nos3, nbs2f, 3, 6), bbJ2p3(nos3, nbs2p, 3, 6), bbJ4l3(nos3, nbs4, 6, 6)
        integer :: l2f(3), l2p(3), l0p(2), l4l(6)
        integer :: lJ1(5), lJ3(6), jcs2(3), jcs0(2), jcs4(6), jcsJ1(5), jcsJ3(6)
        real :: jc2(3), jc0(2), jc4(6), jcJ1(5), jcJ3(6) 
        real*8, allocatable :: Zx(:), Zy(:)
        complex*16, intent(out) :: H(nos+nbs, nos+nbs)
        integer :: counterx, countery
        counterx = 0
        countery = 0
        
        l2f = (/ 3, 3, 3 /)
        l2p = (/ 1, 1, 1/)
        jc2 = (/ 3./2., 3./2., 1./2./)
        jcs2 =  (/ 2, 1, 1 /)
        
        l0p = (/1, 1/)
        jc0 = (/3/2., 1/2./)
        jcs0 = (/1, 1/)

        l4l = (/3, 3, 3, 5, 5, 5/)
        jc4 = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        jcs4 = (/1, 2, 1, 1, 2, 1/)


        jcJ1 = jcb
        lJ1 = lb
        jcsJ1 = jcsb

        jcJ3 = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        lJ3 =  (/ 2, 2, 2, 4, 4, 4/)
        jcsJ3 = (/ 1, 2, 1, 1, 2, 1/)
        counterx = 0
        countery = 0
        H = Zer
        ! The matrix is really sparse, so we start filling by parts.
        ! Start with the coupling between the initial states to the intermediates:
        do counterx=1,nos1
            allocate(Zx(5))
            Zx = ZJ1(counterx,:)
            do countery =  nos+1, nos+nbs
                if(countery.le.nos+nbs2f) then
                    allocate(Zy(3))
                    spi1 = counterx
                    spi2 = countery-nos
                    Zy = ZJ2f(spi2,:) 
                    H(counterx, countery) = RedElem_b(3, Zy, 5, Zx, jcJ1, jcsJ1, lJ1, 1, jc2, jcs2, l2f, 2, bbJ2f1(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum_c(2, 0, 1 , 0, evec)
                    !write(6,*) "For J2f state ", spi2, " in ", spi1, H(counterx, countery)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nos+nbs2f+nbs2p) then
                    allocate(Zy(3))
                    spi1 = counterx
                    spi2 = countery-nos-nbs2f
                    Zy = ZJ2p(spi2,:)
                    H(counterx, countery) = RedElem_b(3, Zy, 5, Zx, jcJ1, jcsJ1, lJ1, 1, jc2, jcs2, l2p, 2, bbJ2p1(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum_c(2, 0, 1, 0, evec)
                    !write(6,*) "For J2p state ", spi2, " in ", spi1, H(counterx, countery)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nos+nbs2f+nbs2p+nbs0) then
                    allocate(Zy(2))
                    spi1 = counterx
                    spi2 = countery-nos-nbs2f-nbs2p
                    Zy = ZJ0p(spi2,:)
                    H(counterx, countery) = RedElem_b(2, Zy, 5, Zx, jcJ1, jcsJ1, lJ1, 1, jc0, jcs0, l0p, 0, bbJ0p1(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(0-0) * qsum_c(0,0,1,0,evec)
                    !write(6,*) "For J0p state ", spi2, " in ", spi1, H(counterx, countery)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                else
                    allocate(Zy(6))
                    Zy = 0
                endif
            
                deallocate(Zy)
                H(countery, counterx) = H(counterx, countery)
            end do
            deallocate(Zx)
        end do

        do counterx=nos1+1,nos
            allocate(Zx(6))
            Zx = ZJ3(counterx-nos1,:)
            do countery =  nos+1, nos+nbs
                if(countery.le.nos+nbs2f) then
                    allocate(Zy(3))
                    spi1 = counterx-nos1
                    spi2 = countery-nos
                    Zy = ZJ2f(spi2,:) 
                    H(counterx, countery) = RedElem_b(3, Zy, 6, Zx, jcJ3, jcsJ3, lJ3, 3, jc2, jcs2, l2f, 2, bbJ2f3(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum_c(2, 0, 3, 0, evec)
                    !write(6,*) "For J2f state ", spi2, " in ", spi1, H(counterx, countery)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nos+nbs2f+nbs2p) then
                    allocate(Zy(3))
                    spi1 = counterx-nos1
                    spi2 = countery-nos-nbs2f
                    Zy = ZJ2p(spi2,:)
                    H(counterx, countery) = RedElem_b(3, Zy, 6, Zx, jcJ3, jcsJ3, lJ3, 3, jc2, jcs2, l2p, 2, bbJ2p3(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum_c(2, 0, 3, 0, evec)
                    !write(6,*) "For J2p state ", spi2, " in ", spi1, H(counterx, countery)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.gt.nos+nbs2f+nbs2p+nbs0) then
                    allocate(Zy(6))
                    spi1 = counterx-nos1
                    spi2 = countery-nos-nbs2f-nbs2p-nbs0
                    Zy = ZJ4(spi2,:)
                    H(counterx, countery) = RedElem_b(6, Zy, 6, Zx, jcJ3, jcsJ3, lJ3, 3, jc4, jcs4, l4l, 4, bbJ4l3(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(4-0) * qsum_c(4, 0, 3, 0, evec)
                    !write(6,*) "For J4l state ", spi2, " in ", spi1, H(counterx, countery)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                else 
                    allocate(Zy(2))
                    Zy = 0.d0
                endif
                deallocate(Zy)
                H(countery, counterx) = H(counterx, countery)
            end do
            deallocate(Zx)
        end do
    end subroutine TPmatNC
    !
    subroutine TPmatNC_J1(evec, bbJ2f, bbJ2p, bbJ0p, Zs, ZJ2f, ZJ2p, ZJ0p, H)
        real*8, INTENT(IN) :: evec(3)
        real*8, INTENT(IN) :: Zs(nos1,5), ZJ2f(nbs2f, 3), Zj2p(nbs2p, 3), ZJ0p(nbs0, 2)
        real*8, INTENT(IN) :: bbJ2f(nos1, nbs2f, 3, 5), bbJ2p(nos1, nbs2p, 3, 5), bbJ0p(nos1, nbs0, 2, 5)
        integer :: l2f(3), l2p(3), l0p(2)
        integer :: lJ1(5), lJ3(6), jcs2(3), jcs0(2), jcsJ1(5), jcsJ3(6)
        real :: jc2(3), jc0(2), jcJ1(5), jcJ3(6) 
        real*8, allocatable :: Zx(:), Zy(:)
        complex*16, intent(out) :: H(nos1+nbs2f+nbs2p+nbs0, nos1+nbs2f+nbs2p+nbs0)
        integer :: counterx, countery
        counterx = 0
        countery = 0
        l2f = (/ 3, 3, 3 /)
        l2p = (/ 1, 1, 1/)
        jc2 = (/ 3./2., 3./2., 1./2./)
        jcs2 =  (/ 2, 1, 1 /)
        l0p = (/1, 1/)
        jc0 = (/3/2., 1/2./)
        jcs0 = (/1, 1/)
        jcJ1 = jcb
        lJ1 = lb
        jcsJ1 = jcsb
        jcJ3 = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        lJ3 =  (/ 2, 2, 2, 4, 4, 4/)
        jcsJ3 = (/ 1, 2, 1, 1, 2, 1/)
        counterx = 0
        countery = 0
        H = Zer
        ! The matrix is really sparse, so we start filling by parts.
        ! Start with the coupling between the initial states to the intermediates:
        do counterx=1,nos1
            allocate(Zx(5))
            Zx = Zs(counterx,:)
            do countery =  nos1+1, nos1+nbs2f+nbs2p+nbs0
                if(countery.le.nos1+nbs2f) then
                    allocate(Zy(3))
                    spi1 = counterx
                    spi2 = countery-nos1
                    Zy = ZJ2f(spi2,:) 
                    H(counterx, countery) = RedElem_b(3, Zy, 5, Zx, jcb, jcsb, lb, 1, jc2, jcs2, l2f, 2, bbJ2f(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum(2, 0, evec)
                    !write(6,*) "For J2f state ", spi2, " in ", spi1, H(counterx, countery)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                elseif(countery.le.nos1+nbs2f+nbs2p) then
                    allocate(Zy(3))
                    spi1 = counterx
                    spi2 = countery-nos1-nbs2f
                    Zy = ZJ2p(spi2,:)
                    H(counterx, countery) = RedElem_b(3, Zy, 5, Zx, jcb, jcsb, lb, 1, jc2, jcs2, l2p, 2, bbJ2p(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(2-(0)) * qsum(2, 0, evec)
                    !write(6,*) "For J2p state ", spi2, " in ", spi1, H(counterx, countery)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                else
                    allocate(Zy(2))
                    spi1 = counterx
                    spi2 = countery-nos1-nbs2f-nbs2p
                    Zy = ZJ0p(spi2,:)
                    H(counterx, countery) = RedElem_b(2, Zy, 5, Zx, jcb, jcsb, lb, 1, jc0, jcs0, l0p, 0, bbJ0p(spi1,spi2,:,:))
                    !write(6,*) "before q sum ", spi1, spi2, counterx, countery, H(counterx, countery)
                    H(counterx, countery) = H(counterx, countery)*(-1)**(0-0) * qsum(0,0,evec)
                    !write(6,*) "For J0p state ", spi2, " in ", spi1, H(counterx, countery)
                    !write(6,*) "Currently at ", spi1, spi2, counterx, countery, H(counterx, countery)
                endif
                deallocate(Zy)
                H(countery, counterx) = H(counterx, countery)
            end do
            deallocate(Zx)
        end do
    end subroutine TPmatNC_J1
    !
    function nmod(x, n)
        integer :: n, x, nmod
        if(mod(x,n).eq.0) then
            nmod=n
        else
            nmod = mod(x,n)
        endif
    end function nmod
end module