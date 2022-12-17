module ArJ1J2
    use j_recoup
    use errfun
    use matmod
    implicit none
    real*8, parameter :: RydAr = 109735.81d0, CMperAU = 2.1947463136320d5
    real*8, parameter :: eVperAU = 27.211
    real*8, parameter :: PI=4.0d0*ATAN(1.0d0)
    real*8, parameter :: auoftime = 2.418884326d-17, auofI = 3.50944758e16
    integer, parameter :: nb=5
    integer, parameter :: jcsb(5)=(/1, 2, 1, 1, 1/), lb(5)=(/2, 2, 2, 0, 0/)
    real, parameter :: jcb(5)=(/1/2., 3/2., 3/2., 1/2., 3/2./)
    complex*16, parameter :: Zer = (0.d0, 0.0d0), Im = (0.0d0, 1.0d0), Zone = (1.0d0, 0.0d0)
    integer, private :: i, i1, i2, i3
    contains
    !
    subroutine MQDT_bt(n, Ecm, Uiba, mu0, mu1, nus, Vbaa, Ilcm, A, Z, Zn, det, ryd, idet, Uiap)
        integer, intent(in) :: n, ryd, idet
        real*8, intent(in) :: Ecm, Uiba(n,n), Vbaa(n,n), mu0(n), mu1(n), nus(n), Ilcm, Uiap(n,n)
        real*8, intent(out) :: A(n), Zn(n), Z(n)
        real*8 :: Uia(n,n), Mia(n,n), Zmat(n,n),  eigr(n), eigi(n), vec(n,n), nrm, mun(n)
        integer :: mine, ii, jj
        complex*16 :: det
        external eig_rs
        Mia = 0.0d0
        Zmat = 0.0d0
        mine = 0
        vec = 0.0d0
        Uia = matmul(Uiba, Vbaa)
        if(ryd.eq.3) then
            mun = mu0-mu1*((Ecm-Ilcm)/RydAr)-(/-0.023601,0.016847/)
            call Matij(n, mu0-mu1*((Ecm-Ilcm)/RydAr)-(/-0.023601,0.016847/), nus, Uia, Mia, Zmat)
        elseif(ryd.eq.0) then
            mun = mu0-mu1*0.5d0/fnu(Ecm,Ilcm)**2
            call Matij(n, mu0-mu1*0.5d0/fnu(Ecm,Ilcm)**2, nus, Uia, Mia, Zmat)
        else
            mun = mu0-mu1*1.0d0/fnu(Ecm,Ilcm)**2
            call Matij(n, mu0-mu1*1.0d0/fnu(Ecm,Ilcm)**2, nus, Uia, Mia, Zmat)
        endif

        call eig_r(n, Mia, eigr, eigi, vec)
        det = 1.0d0
        do i=1,n
            det = det * cmplx( eigr(i), eigi(i), 8 )
        end do
        mine = minloc(eigr**2+eigi**2, DIM=n)
        if(.not.(eigi(mine).eq.0.0d0).and.(idet.eq.0)) then
            print *, "Eigen value of minimum norm non unique"
            print *, 'List of eigen values: '
            print *, eigr(:)
            print *, eigi(:)
            print *, "Smallest norm ", mine
        else
        endif
        A = vec(:,mine)
        Zn = matmul(Zmat,A)
        Z = Zn
        if(ryd.eq.0) then
            nrm = sum(nus**3 * Zn**2)+sum(mu1 * A**2)
        elseif(ryd.eq.1) then
            nrm = sum(nus**3 * Zn**2)+sum(2.0d0*mu1 * A**2)
        elseif(ryd.eq.3) then
            nrm = sum(nus**3 * Zn**2)+sum(-2.0d0*mu1 * A**2)
        endif
        
        do i1=1,n
        do ii=1, n
            do jj=1,n
                nrm = nrm+1/PI * Uiap(i1,ii)*Uia(i1,jj)&
                               *SIN(PI*(mun(ii)-mun(jj)))&
                               *Zn(ii)*Zn(jj)
            enddo
        enddo
        enddo
        Zn = Zn/SQRT(nrm)
        return 
    end subroutine MQDT_bt
    !
    subroutine MQDT_bet(n,Ecm,Uiba,mu0,mu1,Vbaa,Is,ryd,Sdp,Zco, iopen, iclosed, nus)
        integer :: n, nc, no, ryd, ii1, ii2
        real*8 :: Ecm, Uiba(n,n), Vbaa(n,n), Is(n), mu0(n), mu1(n), Ilcm
        real*8 :: Uia(n,n), etas(n)
        complex*16 :: Smatd(n,n)
        integer, allocatable,INTENT(OUT) :: iopen(:), iclosed(:)
        real*8, allocatable, INTENT(OUT) :: nus(:)
        complex*16, allocatable :: Soo(:,:), Sco(:,:), Soc(:,:), Scc(:,:), Ph(:,:)
        complex*16, allocatable, INTENT(OUT) :: Sdp(:,:),Zco(:,:)
        complex*16, allocatable :: Inv(:,:), diff(:,:)
        
        Ilcm = minval(Is)
        Uia = 0.0d0
        etas = 0.d0
        Uia = MATMUL(Uiba,Vbaa)
        
        if(ryd.eq.0) then
            call Smatij(n, mu0-mu1*0.5d0*((Ecm-Ilcm)/RydAr), etas, Uia, Smatd)
        else
            write(65,*) Ecm,mu0-mu1*1.0d0*((Ecm-Ilcm)/RydAr) 
            call Smatij(n, mu0-mu1*1.0d0*((Ecm-Ilcm)/RydAr), etas, Uia, Smatd)
        endif
        
        nc = 0
        
        do ii1=1,n
            if(Ecm.lt.Is(ii1)) nc=nc+1
        end do
        
        no = n-nc
        allocate(Soo(no,no),Sco(nc,no),Soc(no,nc),Scc(nc,nc))
        allocate(Sdp(no,no), Zco(nc,no),Inv(nc,nc),nus(nc))
        allocate(iclosed(nc),iopen(no))
        allocate(diff(nc,nc),Ph(nc,nc))
        

        nc = 0
        no = 0
        do ii1=1,n
            if(Ecm.lt.Is(ii1)) then
                nc=nc+1
                iclosed(nc) = ii1 
            else
                no = no+1
                iopen(no) = ii1
            endif
        end do
        !write(6,*) "Filling the matrices for closed channels"
        ! Fill in the matrix
        do ii1=1,nc
            nus(ii1) = fnu(Ecm,Is(iclosed(ii1)))
            do ii2=1,nc
                Scc(ii1,ii2) = Smatd(iclosed(ii1),iclosed(ii2))
                if(ii1.eq.ii2) then
                    diff(ii1,ii2) = Scc(ii1,ii2)-EXP(2*Im*PI*nus(ii1))
                    Ph(ii1,ii2) = EXP(Im*PI*nus(ii1))
                else
                    Ph(ii1,ii2) = Zer
                    diff(ii1,ii2) = Scc(ii1,ii2)
                endif
            end do
            do ii2=1,no
                Sco(ii1,ii2) = Smatd(iclosed(ii1),iopen(ii2))
            end do
        end do

        !write(6,*) "Filling the matrices for open channels"
        ! Keep filling them
        do ii1=1,no
            do ii2=1,nc
                Soc(ii1,ii2) = Smatd(iopen(ii1),iclosed(ii2))
            end do
            do ii2=1,no
                Soo(ii1,ii2) = Smatd(iopen(ii1),iopen(ii2))
            end do
        end do
        !write(6,*) diff(:,:)
        !write(6,*) Ph(:,:)
        call zInvMat(nc,diff,Inv)
        Sdp = Soo-MATMUL(Soc,MATMUL(Inv,Sco))
        Zco = MATMUL(Ph,MATMUL(Inv,Sco))
        deallocate(Soo,Sco,Soc,Scc)
        deallocate(Inv,diff,Ph)
    end subroutine MQDT_bet
    !
    subroutine MQDT_at(n, Ecm, Uiba, mu0, mu1, etas, Vbaa, Ilcm, ryd, Smatd, Kmat)
        integer, intent(in) :: n, ryd
        real*8, intent(in) :: Ecm, Uiba(n,n), Vbaa(n,n), mu0(n), mu1(n), Ilcm, etas(n)
        real*8 :: Uia(n,n)
        complex*16, intent(out) :: Smatd(n,n)
        real*8, intent(out) :: Kmat(n,n)
        Uia = 0.0d0
        Uia = MATMUL(Uiba,Vbaa)
        if(ryd.eq.0) then
            call Smatij(n, mu0+mu1*0.5d0*((Ecm-Ilcm)/RydAr), etas, Uia, Smatd)
            call Kmatij(n, mu0+mu1*0.5d0*((Ecm-Ilcm)/RydAr), Uia, Kmat)
        else
            call Smatij(n, mu0+mu1*1.0d0*((Ecm-Ilcm)/RydAr), etas, Uia, Smatd)
            call Kmatij(n, mu0+mu1*1.0d0*((Ecm-Ilcm)/RydAr), Uia, Kmat)
        endif        
    end subroutine MQDT_at
    !
    subroutine cont_Amp(ne, Econt, in, lc, Jin, Min, A, Zs, Ebs, evec, B, regint, iregint, gam, Eo, w, to)
        integer :: in, Jin, Min
        integer :: ne, lc
        real :: jce(ne)
        real*8:: Eo, gam, w, to, DEau
        integer:: jcse(ne), le(ne)
        real*8 :: Econt, Ebs(4), evec(3), regint(4, nb, ne), iregint(4, nb, ne)
        real*8, intent(in):: Zs(4,5)
        real*8 :: Zi(5)
        complex*16, intent(in) :: A(4)
        complex*16 :: matelem, Sdmat(ne,ne)
        real*8 :: Kmat(ne,ne)
        complex*16, intent(out) :: B
        if(Jin.eq.0) then
            call J0_MQDT_at(Econt, Sdmat, Kmat, 0)
            le = (/ 1, 1/)
            jce = (/ 3./2, 1./2./)
            jcse = (/ 1, 1/)
        else
            if(lc.eq.3) then
                call J2f_MQDT_at(Econt, Sdmat, Kmat, 0)
                le = (/ 3, 3, 3 /)
                jce = (/ 3./2., 3./2., 1./2. /)
                jcse =  (/ 2, 1, 1 /)
            else
                call J2p_MQDT_at(Econt, Sdmat, Kmat, 0)
                le = (/ 1, 1, 1 /)
                jce = (/ 3./2., 3./2., 1./2. /)
                jcse = (/ 2, 1, 1 /)
            endif
        endif 

        B = Zer
        
        do i=1,4
            Zi = Zs(i,:)
            matelem = (-1)**(Jin-Min) * qsum(Jin, Min, evec)
            matelem = matelem * RedElem(ne, Zi, Sdmat, in, jce, le, jcse, Jin, regint(i,:,:), iregint(i,:,:))
            DEau = Ebs(i)/CMperAU+w-Econt/CMperAU
            B = B - Im/2.0d0 * Eo * gam * SQRT(PI) * A(i) * matelem * FourF(gam, DEau, to)
        end do
    end subroutine cont_Amp
    !
    subroutine cont_Amp_Re(ne, Econt, in, lc, Jin, Min, A, Zs, Ebs, evec, B, regint, iregint, gam, Eo, w, to)
        integer :: in, Jin, Min
        integer :: ne, lc
        real :: jce(ne)
        real*8:: Eo, gam, w, to, DEau
        integer:: jcse(ne), le(ne)
        real*8 :: Econt, Ebs(4), evec(3), regint(4, nb, ne), iregint(4, nb, ne)
        real*8, intent(in):: Zs(4,5)
        real*8 :: Zi(5)
        complex*16, intent(in) :: A(4)
        complex*16 :: matelem
        real*8 :: Uia(ne,ne), mua(ne)
        complex*16 :: Bas(ne)
        complex*16, intent(out) :: B(4)
        if(Jin.eq.0) then
            call J0_MQDT_params(Econt, Uia, mua)
            le = (/ 1, 1/)
            jce = (/ 3./2, 1./2./)
            jcse = (/ 1, 1/)
        else
            if(lc.eq.3) then
                call J2f_MQDT_params(Econt, Uia, mua)
                le = (/ 3, 3, 3 /)
                jce = (/ 3./2., 3./2., 1./2. /)
                jcse =  (/ 2, 1, 1 /)
            else
                call J2p_MQDT_params(Econt, Uia, mua)
                le = (/ 1, 1, 1 /)
                jce = (/ 3./2., 3./2., 1./2. /)
                jcse = (/ 2, 1, 1 /)
            endif
        endif 

        Bas = 0.0d0
        B = Zer
        do i=1,4
            do i1=1, ne
                Zi = Zs(i,:)
                matelem = (-1)**(Jin-Min) * qsum(Jin, Min, evec)              
                matelem = matelem * RedElem_c_SW(ne, nb, Zs(i,:), Uia, mua, i1, jce, le, jcse,&
                                Jin, jcb, lb, jcsb, 1, regint(i,:,:), iregint(i,:,:))
                Bas(i1) = matelem
            end do

            do i2=1,ne
                B(i) = B(i) + BAs(i2)*EXP(-Im*PI*mua(i2))*Uia(in,i2)
            end do
            DEau = Ebs(i)/CMperAU+w-Econt/CMperAU
            B(i) = - Im/2.0d0 * Eo * gam * SQRT(PI) * B(i) * FourF(gam, DEau, to)
        end do
    end subroutine cont_Amp_Re
    !
    subroutine cont_Amp_Z(ne, Econt, in, lc, Jin, Min, A, Zs, Ebs, evec, B, regint, iregint, gam, Eo, w, to)
        integer :: in, Jin, Min
        integer :: ne, lc
        real :: jce(ne)
        real*8:: Eo, gam, w, to, DEau
        integer:: jcse(ne), le(ne)
        real*8 :: Econt, Ebs(4), evec(3), regint(4, nb, ne), iregint(4, nb, ne)
        real*8, intent(in):: Zs(4,5)
        real*8 :: Zi(5)
        complex*16, intent(in) :: A(4)
        complex*16 :: matelem, Sdmat(ne,ne)
        real*8 :: Kmat(ne,ne)
        complex*16, intent(out) :: B(4)
        if(Jin.eq.0) then
            call J0_MQDT_at(Econt, Sdmat, Kmat, 0)
            le = (/ 1, 1/)
            jce = (/ 3./2, 1./2./)
            jcse = (/ 1, 1/)
        else
            if(lc.eq.3) then
                call J2f_MQDT_at(Econt, Sdmat, Kmat, 0)
                le = (/ 3, 3, 3 /)
                jce = (/ 3./2., 3./2., 1./2. /)
                jcse =  (/ 2, 1, 1 /)
            else
                call J2p_MQDT_at(Econt, Sdmat, Kmat, 0)
                le = (/ 1, 1, 1 /)
                jce = (/ 3./2., 3./2., 1./2. /)
                jcse = (/ 2, 1, 1 /)
            endif
        endif 

        B = Zer
        
        do i=1,4
            Zi = Zs(i,:)
            matelem = (-1)**(Jin-Min) * qsum(Jin, Min, evec)
            matelem = matelem * RedElem_c(ne, nb, Zs(i,:), Sdmat, in, jce, le, jcse,&
                             Jin, jcb, lb, jcsb, 1, regint(i,:,:), iregint(i,:,:))
            DEau = Ebs(i)/CMperAU+w-Econt/CMperAU
            B(i) = B(i) - Im/2.0d0 * Eo * gam * SQRT(PI) * matelem * FourF(gam, DEau, to)
        end do
    end subroutine cont_Amp_Z
    !
    subroutine cont_Amp_TP(ne, Econt, in, Jin, Min, B,& !Continuum state info
        A, Zs, Ebs,& !Initial bound state info
        nin, EintsJ2f, ZintsJ2f, EintsJ2p, ZintsJ2p, & !Intermediate state information
        nin0, EintsJ0p, ZintsJ0p, &
        bbJ2f, bbJ2p, bbJ0p, regJ2fJin, iregJ2fJin, regJ2pJin, iregJ2pJin,&
        iregJ0pJin, regJ0pJin, & !Radial integral tensors
        evec, gam, Eo, w, to)
        integer :: nin, ne, nin0
        integer :: in, Jin, Min
        real*8 :: Econt, Ebs(4), EintsJ2f(nin), EintsJ2p(nin), EintsJ0p(nin0)
        real*8 :: Zs(4,5), ZintsJ2f(nin,3), ZintsJ2p(nin,3), ZintsJ0p(nin0, 2)
        real*8 :: bbJ2f(nin, 4, 3, 5), bbJ2p(nin, 4, 3, 5), bbJ0p(nin0, 4, 2, 5)
        real*8 :: regJ2fJin(nin, ne, 3), iregJ2fJin(nin, ne, 3)
        real*8 :: regJ2pJin(nin, ne, 3), iregJ2pJin(nin, ne, 3)
        real*8 :: regJ0pJin(nin0,ne, 2), iregJ0pJin(nin0,ne, 2)
        real*8 :: evec(3), gam, Eo, w, to
        real :: jc2(3), jc0(2), jce(ne)
        integer :: l2f(3), l2p(3), l0p(2), le(ne), jcs2(3), jcs0(2), jcse(ne)
        complex*16 :: A(4), B, xisum
        complex*16 :: matelem1, matelem2, Sdmat(ne,ne), summ
        real*8 :: Kmat(ne,ne)

        l2f = (/ 3, 3, 3 /)
        l2p = (/ 1, 1, 1/)
        jc2 = (/ 3./2., 3./2., 1./2./)
        jcs2 =  (/ 2, 1, 1 /)
        
        l0p = (/1, 1/)
        jc0 = (/3/2., 1/2./)
        jcs0 = (/1, 1/)
        if(Jin.eq.1) then
            jce = jcb
            le = lb
            jcse = jcsb
            call J1_MQDT_at(Econt, Sdmat, Kmat, 0)
        elseif(Jin.eq.3) then
            jce = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
            le =  (/ 2, 2, 2, 4, 4, 4/)
            jcse = (/ 1, 2, 1, 1, 2, 1/)
            call J3_MQDT_at(Econt, Sdmat, Kmat, 0)
        else 
            write(6,*) "You got the J wrong in the continuum, correct."
            stop
        endif
        B = Zer
        do i1=1, 4
        !write(6,*) " Bound State", i1
        xisum = Zer
        !Sum over J=2 states
        do i2=1, nin
            !write(6,*) "Doing J=2 with the outer f electron state ", i2 
            summ = Zer
            matelem1 = RedElem_b(3, ZintsJ2f(i2,:), 5, Zs(i1,:), jcb, jcsb, lb, 1, jc2, jcs2, l2f, 2, bbJ2f(i2,i1,:,:)) 
            do i3=-2,2
                summ=summ+(-1)**(2-(i3)) * qsum(2, i3, evec)
            end do
            !write(6,*) "Bound J1 to J2f", matelem1, summ
            matelem1 = matelem1 * summ

            summ = Zer
            matelem2 = RedElem_c(ne, 3, ZintsJ2f(i2,:), Sdmat, in, jce, le, jcse, Jin, jc2,&
             l2f, jcs2, 2, regJ2fJin(i2,:,:), iregJ2fJin(i2,:,:))
             do i3=-2,2
                summ=summ+(-1)**(Jin-Min) * qsum_c(Jin, Min, 2, i3, evec)
            end do
            !write(6,*) "Bound to c J2f", "to ", Jin, " ", matelem2, summ
            matelem2 = matelem2*summ 
            
            
            xisum = xisum + matelem1*matelem2*twoPhotonAmp(Ebs(i1)/CMperAU,Econt/CMperAU,EintsJ2f(i2)/CMperAU,gam,to,w)

            !write(6,*) "Doing J=2 with p electron " 

            summ = Zer
            matelem1 = RedElem_b(3, ZintsJ2p(i2,:), 5, Zs(i1,:), jcb, jcsb, lb, 1, jc2, jcs2, l2p, 2, bbJ2p(i2,i1,:,:)) 
            do i3=-2,2
                summ=summ+(-1)**(2-(i3)) * qsum(2, i3, evec)
            end do
            !write(6,*) "Bound J1 to J2p", matelem1, summ
            matelem1 = matelem1 * summ
            
            summ = Zer
            matelem2 = RedElem_c(ne, 3, ZintsJ2p(i2,:), Sdmat, in, jce, le, jcse, Jin, jc2,&
             l2p, jcs2, 2, regJ2pJin(i2,:,:), iregJ2pJin(i2,:,:))
            do i3=-2,2
                summ=summ+(-1)**(Jin-Min) * qsum_c(Jin, Min, 2, i3, evec)
            end do
            !write(6,*) "Bound to c J2p to Jin(1)", matelem2, summ
            matelem2 = matelem2* summ

            xisum = xisum + matelem1*matelem2*twoPhotonAmp(Ebs(i1)/CMperAU,Econt/CMperAU,EintsJ2p(i2)/CMperAU,gam,to,w)
            !read(*,*)
        end do
        !Sum over J=0 states
        do i2=1, nin0
            !write(6,*) "Now adding the J=0 for state ", i2
            summ = Zer
            matelem1 = RedElem_b(2, ZintsJ0p(i2,:), 5, Zs(i1,:), jcb, jcsb, lb, 1, jc0, jcs0, l0p, 0, bbJ0p(i2,i1,:,:)) 
            summ=summ+(-1)**(0-(0)) * qsum(0, 0, evec)
            !write(6,*) "Bound J1 to J0p", matelem1, summ
            matelem1 = matelem1 * summ

            summ = Zer
            matelem2 = RedElem_c(ne, 2, ZintsJ0p(i2,:), Sdmat, in, jce, le, jcse, Jin, jc0,&
             l0p, jcs0, 0, regJ0pJin(i2,:,:), iregJ0pJin(i2,:,:))
            summ=summ+(-1)**(Jin-Min) * qsum_c(Jin, Min, 0, 0, evec)
            !write(6,*) "Bound to c J0p", "to ", Jin, " ", matelem2, summ
            matelem2 = matelem2*summ 

            xisum = xisum + matelem1*matelem2*twoPhotonAmp(Ebs(i1)/CMperAU,Econt/CMperAU,EintsJ0p(i2)/CMperAU,gam,to,w)
            !read(*,*)
        end do
        !write(6,*) "For ", i1, "with amplitude", A(i1)," finalized xi sum this is B ", xisum
        !read(*,*)
        B = B + A(i1)*xisum
        end do
        B =  B*(-1* Eo**2 * gam**2 * PI * 1.0d0/8.0d0)
        !write(6,*) B
        return 
    end subroutine cont_Amp_TP
    !
    subroutine cont_Amp_TP_vec(ne, Econt, in, Jin, Min, B,& !Continuum state info
        Zs, Ebs,& !Initial bound state info
        nin, EintsJ2f, ZintsJ2f, EintsJ2p, ZintsJ2p, & !Intermediate state information
        nin0, EintsJ0p, ZintsJ0p, &
        bbJ2f, bbJ2p, bbJ0p, regJ2fJin, iregJ2fJin, regJ2pJin, iregJ2pJin,&
        iregJ0pJin, regJ0pJin, & !Radial integral tensors
        evec, gam, Eo, w)
        integer :: nin, ne, nin0
        integer :: in, Jin, Min
        real*8 :: Econt, Ebs(4), EintsJ2f(nin), EintsJ2p(nin), EintsJ0p(nin0)
        real*8 :: Zs(4,5), ZintsJ2f(nin,3), ZintsJ2p(nin,3), ZintsJ0p(nin0, 2)
        real*8 :: bbJ2f(nin, 4, 3, 5), bbJ2p(nin, 4, 3, 5), bbJ0p(nin0, 4, 2, 5)
        real*8 :: regJ2fJin(nin, ne, 3), iregJ2fJin(nin, ne, 3)
        real*8 :: regJ2pJin(nin, ne, 3), iregJ2pJin(nin, ne, 3)
        real*8 :: regJ0pJin(nin0,ne, 2), iregJ0pJin(nin0,ne, 2)
        real*8 :: evec(3), gam, Eo, w, to
        real :: jc2(3), jc0(2), jce(ne)
        integer :: l2f(3), l2p(3), l0p(2), le(ne), jcs2(3), jcs0(2), jcse(ne)
        complex*16 :: B(4), xisum
        complex*16 :: matelem1, matelem2, Sdmat(ne,ne), summ
        real*8 :: Kmat(ne,ne)
        to = 0.0d0
        l2f = (/ 3, 3, 3 /)
        l2p = (/ 1, 1, 1/)
        jc2 = (/ 3./2., 3./2., 1./2./)
        jcs2 =  (/ 2, 1, 1 /)
        
        l0p = (/1, 1/)
        jc0 = (/3/2., 1/2./)
        jcs0 = (/1, 1/)
        if(Jin.eq.1) then
            jce = jcb
            le = lb
            jcse = jcsb
            call J1_MQDT_at(Econt, Sdmat, Kmat, 0)
        elseif(Jin.eq.3) then
            jce = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
            le =  (/ 2, 2, 2, 4, 4, 4/)
            jcse = (/ 1, 2, 1, 1, 2, 1/)
            call J3_MQDT_at(Econt, Sdmat, Kmat, 0)
        else 
            write(6,*) "You got the J wrong in the continuum, correct."
            stop
        endif
        
        do i1=1, 4
        !write(6,*) "Initial bound state", i1
        xisum = Zer
        !Sum over J=2 states
        do i2=1, nin
            !write(6,*) "Doing J=2 with the outer f electron. State ", i2
            summ = Zer
            matelem1 = RedElem_b(3, ZintsJ2f(i2,:), 5, Zs(i1,:), jcb, jcsb, lb, 1, jc2, jcs2, l2f, 2, bbJ2f(i2,i1,:,:)) 
            do i3=-2,2
                summ=summ+(-1)**(2-(i3)) * qsum(2, i3, evec)
            end do
            !write(6,*) "Bound J1 to J2f", matelem1, summ
            matelem1 = matelem1 * summ

            summ = Zer
            matelem2 = RedElem_c(ne, 3, ZintsJ2f(i2,:), Sdmat, in, jce, le, jcse, Jin, jc2,&
             l2f, jcs2, 2, regJ2fJin(i2,:,:), iregJ2fJin(i2,:,:))
             do i3=-2,2
                summ=summ+(-1)**(Jin-Min) * qsum_c(Jin, Min, 2, i3, evec)
            end do
            !write(6,*) "Bound to c J2f", "to ", Jin, " ", matelem2, summ
            matelem2 = matelem2*summ 
            
            
            xisum = xisum + matelem1*matelem2*twoPhotonAmp(Ebs(i1)/CMperAU,Econt/CMperAU,EintsJ2f(i2)/CMperAU,gam,to,w)

            !Doing J=2 with p electron 

            summ = Zer
            matelem1 = RedElem_b(3, ZintsJ2p(i2,:), 5, Zs(i1,:), jcb, jcsb, lb, 1, jc2, jcs2, l2p, 2, bbJ2p(i2,i1,:,:)) 
            do i3=-2,2
                summ=summ+(-1)**(2-(i3)) * qsum(2, i3, evec)
            end do
            !write(6,*) "Bound J1 to J2p", matelem1, summ
            matelem1 = matelem1 * summ
            
            summ = Zer
            matelem2 = RedElem_c(ne, 3, ZintsJ2p(i2,:), Sdmat, in, jce, le, jcse, Jin, jc2,&
             l2p, jcs2, 2, regJ2pJin(i2,:,:), iregJ2pJin(i2,:,:))
            do i3=-2,2
                summ=summ+(-1)**(Jin-Min) * qsum_c(Jin, Min, 2, i3, evec)
            end do
            !write(6,*) "Bound to c J2f to Jin(1)", matelem2, summ
            matelem2 = matelem2* summ

            xisum = xisum + matelem1*matelem2*twoPhotonAmp(Ebs(i1)/CMperAU,Econt/CMperAU,EintsJ2p(i2)/CMperAU,gam,to,w)
            !read(*,*)
        end do
        !Sum over J=0 states
        do i2=1, nin0
            !Now adding the J=0
            !write(6,*) "Now going to the J0 state ", i2
            summ = Zer
            matelem1 = RedElem_b(2, ZintsJ0p(i2,:), 5, Zs(i1,:), jcb, jcsb, lb, 1, jc0, jcs0, l0p, 0, bbJ0p(i2,i1,:,:)) 
            summ=summ+(-1)**(0-(0)) * qsum(0, 0, evec)
            !write(6,*) "Bound J1 to J2f", matelem1, summ
            matelem1 = matelem1 * summ

            summ = Zer
            matelem2 = RedElem_c(ne, 2, ZintsJ0p(i2,:), Sdmat, in, jce, le, jcse, Jin, jc0,&
             l0p, jcs0, 0, regJ0pJin(i2,:,:), iregJ0pJin(i2,:,:))
            summ=summ+(-1)**(Jin-Min) * qsum_c(Jin, Min, 0, 0, evec)
            !write(6,*) "Bound to c J2f", "to ", Jin, " ", matelem2, summ
            matelem2 = matelem2*summ 

            xisum = xisum + matelem1*matelem2*twoPhotonAmp(Ebs(i1)/CMperAU,Econt/CMperAU,EintsJ0p(i2)/CMperAU,gam,to,w)
            !read(*,*)
        end do
        B(i1) = xisum*(-1* Eo**2 * gam**2 * PI * 1.0d0/8.0d0)
        end do
        return 
    end subroutine cont_Amp_TP_vec
    !
    subroutine cont_Amp_TP_vec_SI(ne, Econt, in, Jin, Min, B,& !Continuum state info
        Zs, Ebs,& !Initial bound state info
        nin2f, SelJ2f, EintsJ2f, ZintsJ2f,&
        nin2p, SelJ2p, EintsJ2p, ZintsJ2p, & !Intermediate state information
        nin0, SelJ0, EintsJ0p, ZintsJ0p, &
        bbJ2f, bbJ2p, bbJ0p, regJ2fJin, iregJ2fJin, regJ2pJin, iregJ2pJin,&
        iregJ0pJin, regJ0pJin, & !Radial integral tensors
        evec, gam, Eo, w)
        integer :: ne, nin2f, nin2p, nin0
        integer :: SelJ2f(2), SelJ2p(2), SelJ0(2)
        integer :: in, Jin, Min
        real*8 :: Econt, Ebs(4), EintsJ2f(nin2f), EintsJ2p(nin2p), EintsJ0p(nin0)
        real*8 :: Zs(4,5), ZintsJ2f(nin2f,3), ZintsJ2p(nin2p,3), ZintsJ0p(nin0, 2)
        real*8 :: bbJ2f(nin2f, 4, 3, 5), bbJ2p(nin2p, 4, 3, 5), bbJ0p(nin0, 4, 2, 5)
        real*8 :: regJ2fJin(nin2f, ne, 3), iregJ2fJin(nin2f, ne, 3)
        real*8 :: regJ2pJin(nin2p, ne, 3), iregJ2pJin(nin2p, ne, 3)
        real*8 :: regJ0pJin(nin0,ne, 2), iregJ0pJin(nin0,ne, 2)
        real*8 :: evec(3), gam, Eo, w, to
        real :: jc2(3), jc0(2), jce(ne)
        integer :: l2f(3), l2p(3), l0p(2), le(ne), jcs2(3), jcs0(2), jcse(ne)
        complex*16 :: B(4), xisum
        complex*16 :: matelem1, matelem2, Sdmat(ne,ne), summ
        real*8 :: Kmat(ne,ne)
        to = 0.0d0
        l2f = (/ 3, 3, 3 /)
        l2p = (/ 1, 1, 1/)
        jc2 = (/ 3./2., 3./2., 1./2./)
        jcs2 =  (/ 2, 1, 1 /)
        
        l0p = (/1, 1/)
        jc0 = (/3/2., 1/2./)
        jcs0 = (/1, 1/)
        if(Jin.eq.1) then
            jce = jcb
            le = lb
            jcse = jcsb
            call J1_MQDT_at(Econt, Sdmat, Kmat, 0)
        elseif(Jin.eq.3) then
            jce = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
            le =  (/ 2, 2, 2, 4, 4, 4/)
            jcse = (/ 1, 2, 1, 1, 2, 1/)
            call J3_MQDT_at(Econt, Sdmat, Kmat, 0)
        else 
            write(6,*) "You got the J wrong in the continuum, correct."
            stop
        endif
        
        do i1=1, 4
        !write(6,*) "Initial bound state", i1
        xisum = Zer
        
        !Sum over J=2 states with f electrons
        if(SelJ2f(1).gt.0) then
            !write(6,*) "Inclusion of J2f interm states from ", SelJ2f(1), "to ", SelJ2f(2)
            do i2=SelJ2f(1), SelJ2f(2)
                !write(6,*) "Looping over f", i2
                !write(6,*) "Doing J=2 with the outer f electron. State ", i2
                summ = Zer
                matelem1 = RedElem_b(3, ZintsJ2f(i2,:), 5, Zs(i1,:), jcb, jcsb, lb, 1, jc2, jcs2, l2f, 2, bbJ2f(i2,i1,:,:)) 
                do i3=-2,2
                    summ=summ+(-1)**(2-(i3)) * qsum(2, i3, evec)
                end do
                !write(6,*) "Bound J1 to J2f", matelem1, summ
                matelem1 = matelem1 * summ

                summ = Zer
                matelem2 = RedElem_c(ne, 3, ZintsJ2f(i2,:), Sdmat, in, jce, le, jcse, Jin, jc2,&
                l2f, jcs2, 2, regJ2fJin(i2,:,:), iregJ2fJin(i2,:,:))
                do i3=-2,2
                    summ=summ+(-1)**(Jin-Min) * qsum_c(Jin, Min, 2, i3, evec)
                end do
                write(6,*) "Bound to c J2f", "to ", Jin, " ", matelem2, summ
                matelem2 = matelem2*summ 
                
                
                xisum = xisum + matelem1*matelem2*twoPhotonAmp(Ebs(i1)/CMperAU,Econt/CMperAU,EintsJ2f(i2)/CMperAU,gam,to,w)
            end do
        endif
        !Doing J=2 with p electron 
        if(SelJ2p(1).gt.0) then
            !write(6,*) "Inclusion of J2p interm states from ", SelJ2p(1), "to ", SelJ2p(2)
            do i2=SelJ2p(1), SelJ2p(2)
                !write(6,*) "Looping over p", i2
                summ = Zer
                matelem1 = RedElem_b(3, ZintsJ2p(i2,:), 5, Zs(i1,:), jcb, jcsb, lb, 1, jc2, jcs2, l2p, 2, bbJ2p(i2,i1,:,:)) 
                do i3=-2,2
                    summ=summ+(-1)**(2-(i3)) * qsum(2, i3, evec)
                end do
                !write(6,*) "Bound J1 to J2p", matelem1, summ
                matelem1 = matelem1 * summ
                
                summ = Zer
                matelem2 = RedElem_c(ne, 3, ZintsJ2p(i2,:), Sdmat, in, jce, le, jcse, Jin, jc2,&
                l2p, jcs2, 2, regJ2pJin(i2,:,:), iregJ2pJin(i2,:,:))
                do i3=-2,2
                    summ=summ+(-1)**(Jin-Min) * qsum_c(Jin, Min, 2, i3, evec)
                end do
                !write(6,*) "Bound to c J2f to Jin(1)", matelem2, summ
                matelem2 = matelem2* summ

                xisum = xisum + matelem1*matelem2*twoPhotonAmp(Ebs(i1)/CMperAU,Econt/CMperAU,EintsJ2p(i2)/CMperAU,gam,to,w)
                !read(*,*)
            end do
        end if

        !Sum over J=0 states
        if(SelJ0(1).gt.0) then
            !write(6,*) "Inclusion of J0p interm states from ", SelJ0(1), "to ", SelJ0(2)
            do i2=SelJ0(1), SelJ0(2)
                !Now adding the J=0
                !write(6,*) "Now going to the J0 state ", i2
                summ = Zer
                matelem1 = RedElem_b(2, ZintsJ0p(i2,:), 5, Zs(i1,:), jcb, jcsb, lb, 1, jc0, jcs0, l0p, 0, bbJ0p(i2,i1,:,:)) 
                summ=summ+(-1)**(0-(0)) * qsum(0, 0, evec)
                !write(6,*) "Bound J1 to J2f", matelem1, summ
                matelem1 = matelem1 * summ

                summ = Zer
                matelem2 = RedElem_c(ne, 2, ZintsJ0p(i2,:), Sdmat, in, jce, le, jcse, Jin, jc0,&
                l0p, jcs0, 0, regJ0pJin(i2,:,:), iregJ0pJin(i2,:,:))
                summ=summ+(-1)**(Jin-Min) * qsum_c(Jin, Min, 0, 0, evec)
                !write(6,*) "Bound to c J2f", "to ", Jin, " ", matelem2, summ
                matelem2 = matelem2*summ 
                !write(6,*) i2, EintsJ0p(i2)
                xisum = xisum + matelem1*matelem2*twoPhotonAmp(Ebs(i1)/CMperAU,Econt/CMperAU,EintsJ0p(i2)/CMperAU,gam,to,w)
                !read(*,*)
            end do
        end if
        B(i1) = xisum*(-1* Eo**2 * gam**2 * PI * 1.0d0/8.0d0)
        end do
        return 
    end subroutine cont_Amp_TP_vec_SI
    !
    subroutine cont_Amp_TP_vec_SI_Re(ne, Econt, in, Jin, Min, B,& !Continuum state info
        Zs, Ebs,& !Initial bound state info
        nin2f, SelJ2f, EintsJ2f, ZintsJ2f,&
        nin2p, SelJ2p, EintsJ2p, ZintsJ2p, & !Intermediate state information
        nin0, SelJ0, EintsJ0p, ZintsJ0p, &
        bbJ2f, bbJ2p, bbJ0p, regJ2fJin, iregJ2fJin, regJ2pJin, iregJ2pJin,&
        iregJ0pJin, regJ0pJin, & !Radial integral tensors
        evec, gam, Eo, w)
        integer :: ne, nin2f, nin2p, nin0
        integer :: SelJ2f(2), SelJ2p(2), SelJ0(2)
        integer :: in, Jin, Min
        real*8 :: Econt, Ebs(4), EintsJ2f(nin2f), EintsJ2p(nin2p), EintsJ0p(nin0)
        real*8 :: Zs(4,5), ZintsJ2f(nin2f,3), ZintsJ2p(nin2p,3), ZintsJ0p(nin0, 2)
        real*8 :: bbJ2f(nin2f, 4, 3, 5), bbJ2p(nin2p, 4, 3, 5), bbJ0p(nin0, 4, 2, 5)
        real*8 :: regJ2fJin(nin2f, ne, 3), iregJ2fJin(nin2f, ne, 3)
        real*8 :: regJ2pJin(nin2p, ne, 3), iregJ2pJin(nin2p, ne, 3)
        real*8 :: regJ0pJin(nin0,ne, 2), iregJ0pJin(nin0,ne, 2)
        real*8 :: evec(3), gam, Eo, w, to
        real :: jc2(3), jc0(2), jce(ne)
        integer :: l2f(3), l2p(3), l0p(2), le(ne), jcs2(3), jcs0(2), jcse(ne)
        complex*16 :: B(4), xisum
        complex*16 :: matelem1, matelem2, summ, matelema
        real*8 :: UIa(ne,ne), mua(ne)
        to = 0.0d0
        l2f = (/ 3, 3, 3 /)
        l2p = (/ 1, 1, 1/)
        jc2 = (/ 3./2., 3./2., 1./2./)
        jcs2 =  (/ 2, 1, 1 /)
        l0p = (/1, 1/)
        jc0 = (/3/2., 1/2./)
        jcs0 = (/1, 1/)
        if(Jin.eq.1) then
            jce = jcb
            le = lb
            jcse = jcsb
            call J1_MQDT_params(Econt, Uia, mua)
        elseif(Jin.eq.3) then
            jce = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
            le =  (/ 2, 2, 2, 4, 4, 4/)
            jcse = (/ 1, 2, 1, 1, 2, 1/)
            call J3_MQDT_params(Econt, Uia, mua)
        else 
            write(6,*) "You got the J wrong in the continuum, correct."
            stop
        endif
        !write(6,*)Jin, Econt, ((AIMAG(Sdmat(ik,ij)),ik=ij,ne),ij=1,ne)
        do i1=1, 4
        !write(6,*) "Initial bound state", i1
        xisum = Zer
        
        !Sum over J=2 states with f electrons
        if(SelJ2f(1).gt.0) then
            !write(6,*) "Inclusion of J2f interm states from ", SelJ2f(1), "to ", SelJ2f(2)
            do i2=SelJ2f(1), SelJ2f(2)
                !write(6,*) "Looping over f", i2
                !write(6,*) "Doing J=2 with the outer f electron. State ", i2
                summ = Zer
                matelem1 = RedElem_b(3, ZintsJ2f(i2,:), 5, Zs(i1,:), jcb, jcsb, lb, 1, jc2, jcs2, l2f, 2, bbJ2f(i2,i1,:,:)) 
                do i3=-2,2
                    summ=summ+(-1)**(2-(i3)) * qsum(2, i3, evec)
                end do
                !write(6,*) "Bound J1 to J2f", matelem1, summ
                matelem1 = matelem1 * summ

                summ = Zer
                matelem2 = Zer
                do i3 = 1, ne
                    matelema = RedElem_c_SW(ne, 3, ZintsJ2f(i2,:), Uia, mua, i3, jce, le, jcse, Jin, jc2,&
                    l2f, jcs2, 2, regJ2fJin(i2,:,:), iregJ2fJin(i2,:,:))
                    matelem2 = matelem2 + Uia(in,i3)*EXP(Im*PI*mua(i3))*matelema
                end do
                do i3=-2,2
                    summ=summ+(-1)**(Jin-Min) * qsum_c(Jin, Min, 2, i3, evec)
                end do
                !write(6,*) "Bound to c J2f", "to ", Jin, " ", matelem2, summ
                matelem2 = matelem2*summ 
                
                
                xisum = xisum + matelem1*matelem2*twoPhotonAmp(Ebs(i1)/CMperAU,Econt/CMperAU,EintsJ2f(i2)/CMperAU,gam,to,w)
            end do
        endif
        !Doing J=2 with p electron 
        if(SelJ2p(1).gt.0) then
            !write(6,*) "Inclusion of J2p interm states from ", SelJ2p(1), "to ", SelJ2p(2)
            do i2=SelJ2p(1), SelJ2p(2)
                !write(6,*) "Looping over p", i2
                summ = Zer
                matelem1 = RedElem_b(3, ZintsJ2p(i2,:), 5, Zs(i1,:), jcb, jcsb, lb, 1, jc2, jcs2, l2p, 2, bbJ2p(i2,i1,:,:)) 
                do i3=-2,2
                    summ=summ+(-1)**(2-(i3)) * qsum(2, i3, evec)
                end do
                !write(6,*) "Bound J1 to J2p", matelem1, summ
                matelem1 = matelem1 * summ
                
                summ = Zer
                matelem2 = Zer
                do i3=1, ne
                    matelema = RedElem_c_SW(ne, 3, ZintsJ2p(i2,:), Uia, mua, i3, jce, le, jcse, Jin, jc2,&
                    l2p, jcs2, 2, regJ2pJin(i2,:,:), iregJ2pJin(i2,:,:))
                    matelem2 = matelem2 + Uia(in,i3)*EXP(Im*PI*mua(i3))*matelema
                end do 
                do i3=-2,2
                    summ=summ+(-1)**(Jin-Min) * qsum_c(Jin, Min, 2, i3, evec)
                end do
                !write(6,*) "Bound to c J2f to Jin(1)", matelem2, summ
                matelem2 = matelem2* summ

                xisum = xisum + matelem1*matelem2*twoPhotonAmp(Ebs(i1)/CMperAU,Econt/CMperAU,EintsJ2p(i2)/CMperAU,gam,to,w)
                !read(*,*)
            end do
        end if

        !Sum over J=0 states
        if(SelJ0(1).gt.0) then
            !write(6,*) "Inclusion of J0p interm states from ", SelJ0(1), "to ", SelJ0(2)
            do i2=SelJ0(1), SelJ0(2)
                !Now adding the J=0
                !write(6,*) "Now going to the J0 state ", i2
                summ = Zer
                matelem1 = RedElem_b(2, ZintsJ0p(i2,:), 5, Zs(i1,:), jcb, jcsb, lb, 1, jc0, jcs0, l0p, 0, bbJ0p(i2,i1,:,:)) 
                summ=summ+(-1)**(0-0) * qsum(0, 0, evec)
                !write(6,*) "Bound J1 to J2f", matelem1, summ
                matelem1 = matelem1 * summ

                summ = Zer
                matelem2 = Zer
                do i3=1,ne
                    matelema = RedElem_c_SW(ne, 2, ZintsJ0p(i2,:), Uia, mua, i3, jce, le, jcse, Jin, jc0,&
                    l0p, jcs0, 0, regJ0pJin(i2,:,:), iregJ0pJin(i2,:,:))
                    matelem2 = matelem2 + Uia(in,i3)*EXP(Im*PI*mua(i3))*matelema
                end do
                summ=summ+(-1)**(Jin-Min) * qsum_c(Jin, Min, 0, 0, evec)
                !write(6,*) "Bound to c J2f", "to ", Jin, " ", matelem2, summ
                matelem2 = matelem2*summ 
                !write(6,*) i2, EintsJ0p(i2)
                xisum = xisum + matelem1*matelem2*twoPhotonAmp(Ebs(i1)/CMperAU,Econt/CMperAU,EintsJ0p(i2)/CMperAU,gam,to,w)
                !read(*,*)
            end do
        end if
        B(i1) = xisum*(-1* Eo**2 * gam**2 * PI * 1.0d0/8.0d0)
        end do
        return 
    end subroutine cont_Amp_TP_vec_SI_Re
    !
    subroutine Rabi_freq(nin, ni, Zint, Jin, Min, l, W,& 
                         nos, no,  Zos, Res, evec)
        integer :: nin, ni, nos, no, l
        real*8 :: Zint(nin, ni), Zos(nos, no)
        integer :: Jin, Min
        real*8 :: W(nin, nos, ni, no), evec(3)
        real*8 :: matelem1
        complex*16 :: Res(nos, nin)
        complex*16 :: summ
        integer :: jcsi(ni), li(ni)
        real :: jci(ni)
        li = l 
        if(Jin.eq.2) then
            jcsi = (/2,1,1/)
            jci  = (/3/2.,3/2.,1/2./)
        else if(Jin.eq.0) then
            jcsi = (/1,1/)
            jci = (/3/2.,1/2./)
        else
            write(6,*) "Wrong Jin value"
        endif
        do i1=1, nos
            do i2=1, nin
                matelem1 = RedElem_b(ni, Zint(i2,:), no, Zos(i1,:), jcb, jcsb, lb, 1, jci, jcsi, li, Jin, W(i2,i1,:,:)) 
                summ=(-1)**(Jin-Min) * qsum(Jin, Min, evec)
                write(6,*) "Inside rabi for i1 ",i1," and ", i2," we get: ", matelem1, summ
                Res(i1,i2) = matelem1 * summ
            end do
        end do

    end subroutine Rabi_freq
    !
    subroutine J0_MQDT(Ecm, Zn, det, ip, idet,tst)
        integer, parameter :: n=2, J=0
        real*8, parameter :: IP32 = 127109.9d0, IP12= IP32+1432.0d0
        integer :: nn
        integer, intent(in) :: ip, idet
        integer, intent(out) :: tst
        real*8, intent(in) :: Ecm
        real*8, intent(out) :: Zn(n)
        real*8 :: Zjn(n), Z(n), Zj(n)
        complex*16, intent(out) :: det
        complex*16 :: det1
        real*8 :: Ujba(n,n), Ujap(n,n), Uiba(n,n), Uiap(n,n), Vbaa(n,n), A(n)
        integer :: L(n), S(n), lc(n), le(n), Jcs(n)
        real :: je(n), jc(n), se(n), sc(n)
        real*8 :: mu0(n), mu1(n), nus(n), Il
        real*8 :: cang
        nn = n
        tst = 1
        Il = 124000.0d0 !min(IP32, IP12)
        L = (/1, 0/)
        S = (/1, 0/)
        lc = (/1, 1/);
        le = (/1, 1/);
        je = (/3/2., 1/2./)
        jc = (/3/2., 1/2./)
        nus =(/ fnu(Ecm, IP32), fnu(Ecm, IP12) /)
        Jcs = (/1, 1/)
        se = (/1/2., 1/2./)
        sc = (/1/2., 1/2./) 
        mu0 = (/ 0.70d0, 0.52d0 /)
        mu1 = (/ 0.103d0, 0.1832d0 /)
        cang = 0.02076d0!-((Ecm-Il)/RydAr)*(-0.1185d0)
        Vbaa = reshape((/COS(cang), SIN(cang),&
                  -1.0d0*SIN(cang), COS(cang) /),shape(Vbaa), order = (/2,1/))
        call bjcs_kLS(nn, lc, sc, jc, le, jcs, se, L, S, J, Uiba)
        Uiap = 2*0.118*reshape((/-SIN(cang), COS(cang),&
        -1.0d0*COS(cang), -SIN(cang) /),shape(Vbaa), order = (/2,1/))
        Uiap = matmul(Uiba,Uiap)
        call bjj_kLS(nn, lc, sc, jc, le, se, je, L, S, J, Ujba)
        Ujap = 2*0.118*reshape((/-SIN(cang), COS(cang),&
        -1.0d0*COS(cang), -SIN(cang) /),shape(Vbaa), order = (/2,1/))
        Ujap = matmul(Ujba,Ujap)
        call MQDT_bt(nn,Ecm, Uiba, mu0, mu1, nus, Vbaa, Il, A, Z, Zn, det, 3, idet, Uiap)
        call MQDT_bt(nn,Ecm, Ujba, mu0, mu1, nus, Vbaa, Il, A, Zj, Zjn, det1, 3, idet, Ujap)
        do i=1, n
            if((nus(i).lt.le(i)).and.(ABS(Zn(i)).ge.1d-8)) then 
                tst = 0
            endif
        end do
        if(ip.eq.1.and.tst.eq.1) then
            print('((A,e25.9))'), "For energy ", Ecm
            print('(5(A,f12.9))'), "nu12 ", fnu(Ecm, IP12), ', nu32 ', fnu(Ecm, IP32)
            print ('(5(A,f12.9))'), "p1/2 p1/2: ", Zj(1), ", p3/2 p1/2: ", Zj(2), "sum sqrd ", sum(Zj**2)
            print ('(5(A,f12.9))'), "p1/2[1]p: ", Z(1), ", p3/2[1]p: ", Z(2), "sum sqrd ", sum(Z**2)
            print ('(A,e22.9)'), "Value of the determinant: ", abs(det)
        endif
    end subroutine J0_MQDT
    !
    subroutine J1_MQDT(Ecm, Zn, det,ip, idet,tst)
        integer, parameter :: n=5, J=1
        real*8, parameter :: IP32 = 127109.9d0, IP12= IP32+1432.0d0
        integer :: nn
        integer, intent(in) :: ip, idet 
        integer, intent(out) :: tst
        real*8, intent(in) :: Ecm
        real*8, intent(out) :: Zn(n)
        real*8 :: Zjn(n), Z(n), Zj(n)
        complex*16, intent(out) :: det
        complex*16 :: det1
        real*8 :: Uiba(n,n), Ujba(n,n), Uip(n,n), Vbaa(n,n), A(n)
        integer :: L(n), S(n), lc(n), le(n), Jcs(n)
        real :: je(n), jc(n), se(n), sc(n)
        real*8 :: mu0(n), mu1(n), nus(n), Il
        nn = n
        tst = 1
        Uip = 0.0d0
        Il = min(IP12, IP32)
        L  = (/2, 1, 1, 1, 1/)
        S = (/1, 0, 1, 1, 0/)
        lc = (/1, 1, 1, 1, 1/);
        le = (/2, 2, 2, 0, 0/);
        je = (/3/2., 5/2., 3/2., 1/2., 1/2./)
        jc = (/1/2., 3/2., 3/2., 1/2., 3/2./)
        nus =(/ fnu(Ecm, IP12), fnu(Ecm, IP32), fnu(Ecm, IP32), fnu(Ecm, IP12), fnu(Ecm, IP32) /)
        Jcs = (/1, 2, 1, 1, 1/)
        se = (/1/2., 1/2., 1/2., 1/2., 1/2./)
        sc = (/1/2., 1/2., 1/2., 1/2., 1/2./)
        mu0 = (/ 0.214d0, 0.070d0, 0.500d0, 0.154d0, 0.109d0 /)
        mu1 = (/ 1.22d0, 0.66d0, 1.84d0, -0.28d0, -0.40d0 /)
        Vbaa = reshape((/  1.000d0, 0.001d0,  0.000d0, 0.000d0, 0.001d0, &
                        -0.001d0, 1.000d0, -0.005d0, 0.000d0, 0.000d0, &
                        0.000d0, 0.005d0,  1.000d0, 0.000d0, 0.003d0, &
                        0.000d0, 0.000d0,  0.000d0, 1.000d0, 0.000d0, &
                        -0.001d0, 0.000d0, -0.003d0, 0.000d0, 1.000d0 /),shape(Vbaa), order = (/2,1/))
        call bJcs_kLS(nn, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        call MQDT_bt(nn,Ecm, Uiba,mu0, mu1, nus, Vbaa, Il, A, Z, Zn, det, 0, idet, Uip)
        call bjj_kLS(nn, lc, se, jc, le, se, je, L, S, J, Ujba)
        call MQDT_bt(nn,Ecm, Ujba,mu0, mu1, nus, Vbaa, Il, A, Zj, Zjn, det1, 0, idet, Uip)
        do i=1, n
            if((nus(i).lt.le(i)).and.(ABS(Zn(i)).ge.1d-8)) then 
                tst = 0
            endif
        end do
        if(ip.eq.1.and.tst.eq.1) then
            print('((A,f25.9))'), "For energy ", Ecm
            print('(5(A,f12.9))'), "nu12 ", fnu(Ecm, IP12), ', nu32 ', fnu(Ecm, IP32)
            print ('(5(A,f12.9))'), "p1/2 d3/2: ", Zj(1), ", p3/2 d5/2: ", Zj(2), ", p3/2 d3/2: ", Zj(3), &
            ", p1/2 s1/2: ", Zj(4), ", p3/2 s1/2: ", Zj(5), "sum sqrd ", sum(Zj**2)
            print ('(5(A,f12.9))'), "p1/2[1]d: ", Z(1), ", p3/2[2]d: ", Z(2), ", p3/2[1]d: ", Z(3), &
            ", p1/2[1]s: ", Z(4), ", p3/2[1]s: ", Z(5), "sum sqrd ", sum(Z**2)
            print ('(A,e22.9)'), "Value of the determinant: ", abs(det)
        endif
        return
    end subroutine J1_MQDT
    !
    subroutine J2f_MQDT(Ecm, Zn, det, ip, idet,tst)
        integer, parameter :: n=3 , J=2
        real*8, parameter :: IP32 = 127109.88d0, IP12= IP32+1431.41d0
        integer, intent(in) :: ip, idet
        integer, intent(out) :: tst
        real*8, intent(in) :: Ecm
        real*8, intent(out) :: Zn(n)
        complex*16, intent(out) :: det
        integer :: nn, lc(n), le(n), L(n), S(n),  Jcs(n)
        real :: sc(n), jc(n), je(n), K(n), se(n)
        real*8 :: Uiba(n,n), Uip(n,n), Vbaa(n,n), Ujji(n,n)
        real*8 :: mu0(n), mu1(n), nus(n), Ilcm, A(n)
        real*8 :: Zjj(n), Z(n)
        A = 0.0d0
        Ilcm = min(IP32, IP12)
        nn = n
        tst = 1
        lc= (/ 1, 1, 1 /)
        le =(/ 3, 3, 3 /)
        jc = (/ 3./2., 3./2., 1./2. /)
        nus = (/ fnu(Ecm, IP32), fnu(Ecm,IP32), fnu(Ecm,IP12)/)
        mu0 = (/ 0.022d0,-0.001d0,0.020d0 /)
        mu1 = (/-0.045d0,-0.043d0,-0.060d0/)
        je =(/ 5./2., 7./2., 5./2. /)
        K = (/ 3./2., 5./2., 5./2. /)
        Jcs = (/ 2, 1, 1 /)
        L =(/ 2, 3, 2 /)
        S = (/ 1, 1, 0 /)
        sc = (/ 1./2., 1./2., 1./2. /)
        se = (/ 1./2., 1./2., 1./2. /)
        Vbaa = reshape((/1.0d0, 0.d0, 0.d0, & 
                        0.0d0, 1.d0, 0.d0, &
                        0.0d0, 0.d0, 1.d0 /),shape(Vbaa))
        call bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        !call bJcK_kLS(n, lc, sc, jc, le, K, se, L, S, J, Uiba)
        call MQDT_bt(n, Ecm, Uiba, mu0, -1.0d0*mu1, nus, Vbaa, Ilcm, A, Z, Zn, det, 1, idet, Uip)
        call bjj_kJcK(n, lc, sc, jc, le, K, se, je, J, Ujji)
        do i=1, n
            if((nus(i).lt.le(i)).and.(ABS(Zn(i)).ge.1d-8)) then 
                tst = 0
            endif
        end do
        if(ip.eq.1.and.tst.eq.1) then
            Zjj = matmul(Ujji, Z)
            print('(A,f25.9)'), "For energy ", Ecm
            print('(5(A,f12.9))'), "nu12 ", fnu(Ecm, IP12), ', nu32 ', fnu(Ecm, IP32)
            print ('(5(A,f12.9))'), "p3/2 f5/2: ", Zjj(1), ", p3/2 f7/2: ", Zjj(2), ", p1/2 f5/2: ", Zjj(3),&
            "sum sqrd ", sum(Zjj**2)
            print ('(5(A,f12.9))'), "p3/2[2]f: ", Z(1), ", p3/2[1]f: ", Z(2), ", p1/2[1]f: ", Z(3),  "sum sqrd ", sum(Z**2)
            print ('(A,e22.9)'), "Value of the determinant: ", abs(det)
        endif
    end subroutine J2f_MQDT
    !
    subroutine J2p_MQDT(Ecm, Zn, det, ip, idet,tst)
        integer, parameter :: n=3 , J=2
        real*8, parameter :: IP32 = 127109.88d0, IP12= IP32+1431.41d0
        integer, intent(in) :: ip, idet
        integer, intent(out) :: tst
        real*8, intent(in) :: Ecm
        real*8, intent(out) :: Zn(n)
        complex*16, intent(out) :: det
        integer :: nn, lc(n), le(n), L(n), S(n),  Jcs(n)
        real :: sc(n), jc(n), je(n), K(n), se(n)
        real*8 :: Uiba(n,n), Uip(n,n), Vbaa(n,n), Ujji(n,n)
        real*8 :: mu0(n), mu1(n), nus(n), Ilcm, A(n)
        real*8 :: Zjj(n),Z(n)
        Uip = 0.d0
        A = 0.0d0
        Ilcm = min(IP32, IP12)
        nn = n
        tst = 1
        lc= (/ 1, 1, 1 /)
        le =(/ 1, 1, 1 /)
        jc = (/ 3./2., 3./2., 1./2. /)
        nus = (/ fnu(Ecm, IP32), fnu(Ecm,IP32), fnu(Ecm,IP12)/)
        mu0 = (/ 0.649d0,0.668d0,0.694d0 /)
        mu1 = (/0.239d0,0.243d0,0.301d0/)
        je =(/ 3./2., 1./2., 3./2. /)
        K = (/ 3./2., 5./2., 3./2. /)
        Jcs = (/ 2, 1, 1 /)
        L =(/ 1, 2, 2 /)
        S = (/ 1, 0, 1 /)
        sc = (/ 1./2., 1./2., 1./2. /)
        se = (/ 1./2., 1./2., 1./2. /)
        Vbaa = reshape((/0.9985079636885027d0,-0.060044423946200465d0,-0.00016432773856367433d0, & 
                        0.06040476384223767d0,0.9981651117253054d0,0.00023239451655082677d0, &
                        -0.0003331002905561653d0,0.0002772857325016309d0,0.999564364605493d0 /),shape(Vbaa))
        call bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        call MQDT_bt(n, Ecm, Uiba, mu0, -1.0d0*mu1, nus, Vbaa, Ilcm, A, Z, Zn, det, 1, idet, Uip)
        call bjj_kJcs(n,lc,sc,jc,le,real(JCs),se,je,J,Ujji)
        do i=1, n
            if((nus(i).lt.le(i)).and.(ABS(Zn(i)).ge.1d-8)) then 
                tst = 0
            endif
        end do
        if(ip.eq.1.and.tst.eq.1) then
            Zjj = matmul(Ujji, Z)
            print('(A, f25.9)'), "For energy ", Ecm
            print('(5(A,f12.9))'), "nu12 ", fnu(Ecm, IP12), ', nu32 ', fnu(Ecm, IP32)
            print ('(5(A,f12.9))'), "p3/2 p3/2: ", Zjj(1), ", p3/2 p1/2: ", Zjj(2), ", p1/2 p3/2: ", Zjj(3), & 
            "sum sqrd ", sum(Zjj**2)
            print ('(5(A,f12.9))'), "p3/2[2]p: ", Z(1), ", p3/2[1]p: ", Z(2), ", p1/2[1]: ", Z(3),  "sum sqrd ", sum(Z**2)
            print ('(A,e22.9)'), "Value of the determinant: ", abs(det)
        endif
    end subroutine J2p_MQDT
    !
    subroutine J3_MQDT(Ecm, Zn, det, ip, idet,tst)
        integer, parameter :: n=6 , J=3
        real*8, parameter :: IP32 = 127109.88d0, IP12= IP32+1431.41d0
        integer, intent(in) :: ip, idet
        integer, intent(out) :: tst
        real*8, intent(in) :: Ecm
        real*8, intent(out) :: Zn(n)
        complex*16, intent(out) :: det
        integer :: nn, lc(n), le(n), L(n), S(n),  Jcs(n)
        real :: sc(n), jc(n), je(n), K(n), se(n)
        real*8 :: Uiba(n,n), Uip(n,n), th12(n,n), th13(n,n), th23(n,n), Vbaa(n,n), Ujji(n,n)
        real*8 :: thp(n,n)
        real*8 :: mu0(n), mu1(n), nus(n), Ilcm, A(n)
        real*8 :: thsd(3), thsd1(3)
        real*8 :: Zjj(n),Z(n)
        nn = n
        tst = 1.0
        Ilcm = min(IP12, IP32)
        Uip = 0.d0
        A = 0.0d0
        L  = (/3, 3, 2, 3, 3, 4/)
        S =  (/0, 1, 1, 0, 1, 1/)
        lc = (/1, 1, 1, 1, 1, 1/)
        le = (/2, 2, 2, 4, 4, 4/)
        je = (/ 5/2., 3/2., 5/2., 7/2., 9/2., 7/2./)
        jc = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        nus = (/fnu(Ecm, IP32), fnu(Ecm, IP32), fnu(Ecm, IP12), &
                fnu(Ecm, IP32), fnu(Ecm, IP32), fnu(Ecm, IP12) /) !(/ fnu(Ecm, IP12), fnu(Ecm, IP32), fnu(Ecm, IP32), fnu(Ecm, IP12), fnu(Ecm, IP32) /)
        Jcs = (/ 1, 2, 1, 1, 2, 1/)
        se = (/1/2., 1/2., 1/2., 1/2., 1/2., 1/2./)
        sc = (/1/2., 1/2., 1/2., 1/2., 1/2., 1/2./)
        mu0 = (/ 0.23923729d0, 0.40119589d0,  0.27013775d0, 0.0d0, 0.0d0, 0.0d0 /)
        mu1 = (/ 1.01300869d0, 1.32975167d0,  0.89723505d0, 0.0d0, 0.0d0, 0.0d0 /)
        thsd =  (/ 0.01373472d0,  0.2175326d0,   0.0137986d0/)&
        -1.d0/fnu(Ecm,IP32)**2 * (/-1.03428852d0, -0.01790728d0,  0.10383079d0/)
        thsd1 = (/-1.03428852d0, -0.01790728d0,  0.10383079d0/)

        th12 = reshape((/ COS(thsd(1)), -SIN(thsd(1)),  0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                          SIN(thsd(1)),  COS(thsd(1)),  0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                                0.0d0 ,         0.0d0,  1.0d0,  0.0d0,  0.0d0,  0.0d0, &
                                0.0d0,          0.0d0,  0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                                0.0d0,          0.0d0,  0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                                0.0d0,          0.0d0,  0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th12), order = (/2,1/))

        th13 = reshape((/ COS(thsd(2)), 0.0d0, -SIN(thsd(2)),  0.0d0,  0.0d0,  0.0d0, &
                                0.00d0, 1.0d0,        0.00d0,  0.0d0,  0.0d0,  0.0d0, &
                          SIN(thsd(2)), 0.0d0,  COS(thsd(2)),  0.0d0,  0.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th13), order = (/2,1/))
        
        th23 = reshape((/1.00d0,        0.0d0,         0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                         0.00d0, COS(thsd(3)), -SIN(thsd(3)),  0.0d0,  0.0d0,  0.0d0, &
                         0.00d0, SIN(thsd(3)),  COS(thsd(3)),  0.0d0,  0.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th23), order = (/2,1/))
        Vbaa = matmul(th12,matmul(th13, th23))
        
        Uip = 0.d0
        thp = reshape((/ SIN(thsd(1)),  COS(thsd(1)),  0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                        -COS(thsd(1)),  SIN(thsd(1)),  0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                               0.0d0 ,         0.0d0,  1.0d0,  0.0d0,  0.0d0,  0.0d0, &
                               0.0d0,          0.0d0,  0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                               0.0d0,          0.0d0,  0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                               0.0d0,          0.0d0,  0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th12), order = (/2,1/))
        Uip = Uip - 2.0d0*thsd1(1)*matmul(thp, matmul(th13, th23))
        
        thp = reshape((/  SIN(thsd(2)), 0.0d0,  COS(thsd(2)),  0.0d0,  0.0d0,  0.0d0, &
                                0.00d0, 1.0d0,        0.00d0,  0.0d0,  0.0d0,  0.0d0, &
                         -COS(thsd(2)), 0.0d0,  SIN(thsd(2)),  0.0d0,  0.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th13), order = (/2,1/))
        Uip = Uip - 2.0d0*thsd1(2)*matmul(th12, matmul(thp, th23))

        thp = reshape((/ 1.00d0,        0.0d0,         0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                         0.00d0, SIN(thsd(3)),  COS(thsd(3)),  0.0d0,  0.0d0,  0.0d0, &
                         0.00d0,-COS(thsd(3)),  SIN(thsd(3)),  0.0d0,  0.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th23), order = (/2,1/))
        Uip = Uip -2.0d0*thsd1(3)*matmul(th12, matmul(th13, thp))
        call bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        call MQDT_bt(n, Ecm, Uiba, mu0, mu1, nus, Vbaa, Ilcm, A, Z, Zn, det, 1, idet, Uip)
        call bjj_kJcs(n, lc, sc, jc, le, real(Jcs), se, je, J, Ujji)
        do i=1, n
            if((nus(i).lt.le(i)).and.(ABS(Zn(i)).ge.1d-8)) then 
                tst = 0
            endif
        end do
        if(ip.eq.1.and.tst.eq.1) then
            Zjj = matmul(Ujji, Z)
            print('(A, f25.9)'), "For energy ", Ecm
            print('(5(A,f12.9))'), "nu12 ", fnu(Ecm, IP12), ', nu32 ', fnu(Ecm, IP32)
            print ('(5(A,f12.9))'), "p3/2 d3/2: ", Zjj(1), ", p3/2 d3/2: ", Zjj(2), ", p1/2 d1/2: ", Zjj(3), &
            ", p3/2 g7/2: ", Zjj(4), ", p3/2 g9/2: ", Zjj(5), ", p1/2 g7/2: ", Zjj(6), & 
            "  sum sqrd ", sum(Zjj**2)
            print ('(5(A,f12.9))'), "p3/2[1]d: ", Z(1), ", p3/2[2]d: ", Z(2), ", p1/2[1]d: ", Z(3), &
            ", p3/2[1]g: ", Z(4), ", p3/2[2]g: ", Z(5), ", p1/2[1]g: ", Z(6),& 
            "  sum sqrd ", sum(Z**2)
            print ('(A,e22.9)'), "Value of the determinant: ", abs(det)
        endif
    end subroutine J3_MQDT
    !
    subroutine J4_MQDT(Ecm, Zn,det, ip, idet,tst)
        integer, parameter :: n=6 , J=4
        real*8, parameter :: IP32 = 127109.88d0, IP12= IP32+1431.41d0
        integer, intent(in) :: ip, idet
        integer, intent(out) :: tst
        real*8, intent(in) :: Ecm
        real*8, intent(out) :: Zn(n)
        complex*16, intent(out) :: det
        integer :: nn, lc(n), le(n), L(n), S(n),  Jcs(n)
        real :: sc(n), jc(n), je(n), K(n), se(n)
        real*8 :: Uiba(n,n), Uip(n,n), th12(n,n), th13(n,n), th23(n,n), Vbaa(n,n), Ujji(n,n)
        real*8 :: thp(n,n)
        real*8 :: mu0(n), mu1(n), nus(n), Ilcm, A(n)
        real*8 :: thsf(3), thsf1(3)
        real*8 :: Zjj(n),Z(n)
        nn = n
        tst = 1
        Ilcm = min(IP12, IP32)
        A = 0.0d0
        Uip = 0.d0
        L  = (/4, 4, 3, 4, 4, 5/)
        S =  (/0, 1, 1, 0, 1, 1/)
        lc = (/1, 1, 1, 1, 1, 1/)
        le = (/3, 3, 3, 5, 5, 5/)
        je = (/ 5/2., 7/2., 7/2., 9/2., 11/2., 9/2./)
        jc = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        nus = (/fnu(Ecm, IP32), fnu(Ecm, IP32), fnu(Ecm, IP12), &
                fnu(Ecm, IP32), fnu(Ecm, IP32), fnu(Ecm, IP12) /) !(/ fnu(Ecm, IP12), fnu(Ecm, IP32), fnu(Ecm, IP32), fnu(Ecm, IP12), fnu(Ecm, IP32) /)
        Jcs = (/ 1, 2, 1, 1, 2, 1/)
        se = (/1/2., 1/2., 1/2., 1/2., 1/2., 1/2./)
        sc = (/1/2., 1/2., 1/2., 1/2., 1/2., 1/2./)
        mu0 = (/ 0.00992992d0,  0.02038345d0, -0.00285065d0, 0.0d0, 0.0d0, 0.0d0 /)
        mu1 = 0.d0
        thsf =  (/-0.03702163d0,  0.03709483d0,  0.03758697d0/)
        thsf1 = 0.d0

        th12 = reshape((/ COS(thsf(1)), -SIN(thsf(1)),  0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                          SIN(thsf(1)),  COS(thsf(1)),  0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                                0.0d0 ,         0.0d0,  1.0d0,  0.0d0,  0.0d0,  0.0d0, &
                                0.0d0,          0.0d0,  0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                                0.0d0,          0.0d0,  0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                                0.0d0,          0.0d0,  0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th12), order = (/2,1/))

        th13 = reshape((/ COS(thsf(2)), 0.0d0, -SIN(thsf(2)),  0.0d0,  0.0d0,  0.0d0, &
                                0.00d0, 1.0d0,        0.00d0,  0.0d0,  0.0d0,  0.0d0, &
                          SIN(thsf(2)), 0.0d0,  COS(thsf(2)),  0.0d0,  0.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th13), order = (/2,1/))
        
        th23 = reshape((/1.00d0,        0.0d0,         0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                         0.00d0, COS(thsf(3)), -SIN(thsf(3)),  0.0d0,  0.0d0,  0.0d0, &
                         0.00d0, SIN(thsf(3)),  COS(thsf(3)),  0.0d0,  0.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th23), order = (/2,1/))
        Vbaa = matmul(th12,matmul(th13, th23))
        
        Uip = 0.d0
        thp = reshape((/ SIN(thsf(1)),  COS(thsf(1)),  0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                        -COS(thsf(1)),  SIN(thsf(1)),  0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                               0.0d0 ,         0.0d0,  1.0d0,  0.0d0,  0.0d0,  0.0d0, &
                               0.0d0,          0.0d0,  0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                               0.0d0,          0.0d0,  0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                               0.0d0,          0.0d0,  0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th12), order = (/2,1/))
        Uip = Uip - 2.0d0*thsf1(1)*matmul(thp, matmul(th13, th23))
        
        thp = reshape((/  SIN(thsf(2)), 0.0d0,  COS(thsf(2)),  0.0d0,  0.0d0,  0.0d0, &
                                0.00d0, 1.0d0,        0.00d0,  0.0d0,  0.0d0,  0.0d0, &
                         -COS(thsf(2)), 0.0d0,  SIN(thsf(2)),  0.0d0,  0.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th13), order = (/2,1/))
        Uip = Uip - 2.0d0*thsf1(2)*matmul(th12, matmul(thp, th23))

        thp = reshape((/ 1.00d0,        0.0d0,         0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                         0.00d0, SIN(thsf(3)),  COS(thsf(3)),  0.0d0,  0.0d0,  0.0d0, &
                         0.00d0,-COS(thsf(3)),  SIN(thsf(3)),  0.0d0,  0.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th23), order = (/2,1/))
        Uip = Uip -2.0d0*thsf1(3)*matmul(th12, matmul(th13, thp))
        call bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        call MQDT_bt(n, Ecm, Uiba, mu0, mu1, nus, Vbaa, Ilcm, A, Z, Zn, det, 1, idet, Uip)
        call bjj_kJcs(n, lc, sc, jc, le, real(Jcs), se, je, J, Ujji)
        do i=1, n
            if((nus(i).lt.le(i)).and.(ABS(Zn(i)).ge.1d-8)) then 
                tst = 0
            endif
        end do
        if(ip.eq.1.and.tst.eq.1) then
            Zjj = matmul(Ujji, Z)
            print('(A, f25.9)'), "For energy ", Ecm
            print('(5(A,f12.9))'), "nu12 ", fnu(Ecm, IP12), ', nu32 ', fnu(Ecm, IP32)
            print ('(5(A,f12.9))'), "p3/2 f5/2: ", Zjj(1), ", p3/2 f7/2: ", Zjj(2), ", p1/2 f7/2: ", Zjj(3), &
            ", p3/2 h9/2: ", Zjj(4), ", p3/2 h11/2: ", Zjj(5), ", p1/2 h9/2: ", Zjj(6), & 
            "  sum sqrd ", sum(Zjj**2)
            print ('(5(A,f12.9))'), "p3/2[1]f: ", Z(1), ", p3/2[2]f: ", Z(2), ", p1/2[1]f: ", Z(3), &
            ", p3/2[1]h: ", Z(4), ", p3/2[2]h: ", Z(5), ", p1/2[1]h: ", Z(6),& 
            "  sum sqrd ", sum(Z**2)
            print ('(A,e22.9)'), "Value of the determinant: ", abs(det)
        endif
    end subroutine
    !
    subroutine J0_MQDT_bet(Ecm, sdp, Zco, iopen, iclosed, beta)
        integer, parameter :: n=2, J=0
        real*8, parameter :: IP32 = 127109.9d0, IP12= IP32+1432.0d0
        integer :: nn
        real*8, intent(in) :: Ecm
        complex*16, allocatable, intent(out) :: sdp(:,:), Zco(:,:)
        integer, ALLOCATABLE, INTENT(OUT) :: iopen(:), iclosed(:)
        real*8, ALLOCATABLE, INTENT(OUT) :: beta(:)
        real*8 :: Uiba(n,n), Vbaa(n,n)
        integer :: L(n), S(n), lc(n), le(n), Jcs(n)
        real ::  jc(n), se(n), sc(n)
        real*8 :: mu0(n), mu1(n), etas(n), Il
        real*8 :: cang
        nn = n
        Il = 124000.0d0 !min(IP32, IP12)
        L  = (/1, 0/)
        S = (/1, 0/)
        lc = (/1, 1/);
        le = (/1,1/);
        !je = (/3/2., 5/2., 3/2., 1/2., 1/2./)
        jc = (/3/2., 1/2./)
        etas = 0.0d0 !(/ fnu(Ecm, IP32), fnu(Ecm, IP12) /)
        Jcs = (/1, 1/)
        se = (/1/2., 1/2./)
        sc = (/1/2., 1/2./) 
        mu0 = (/ 0.70d0, 0.52d0 /)
        mu1 = (/ 0.103d0, 0.1832d0 /)
        cang = 0.02076d0-((Ecm-124.0d3)/RydAr)*(-0.1185d0)
        Vbaa = reshape((/COS(cang), SIN(cang),&
                  -1.0d0*SIN(cang), COS(cang) /),shape(Vbaa), order = (/2,1/))
        call bjcs_kLS(nn, lc, sc, jc, le, jcs, se, L, S, J, Uiba)
        call MQDT_bet(nn,Ecm,Uiba,mu0,mu1,Vbaa,(/IP32,IP12/),0,Sdp,Zco,iopen, iclosed, beta)

        do i=1,size(beta)
            beta(i) = PI*(beta(i)-le(iclosed(i)))
        end do
    end subroutine
    !
    subroutine J1_MQDT_bet(Ecm, sdp, Zco,iopen, iclosed, beta)
        integer, parameter :: n=5, J=1
        real*8, parameter :: IP32 = 127109.9d0, IP12= IP32+1432.0d0
        integer :: nn
        real*8, intent(in) :: Ecm
        complex*16, allocatable, intent(out) :: sdp(:,:), Zco(:,:)
        integer, allocatable, intent(out) :: iclosed(:), iopen(:)
        real*8, allocatable, intent(out) :: beta(:)
        real*8 :: Uiba(n,n), Ujji(n,n), Vbaa(n,n)
        integer :: L(n), S(n), lc(n), le(n), Jcs(n)
        real :: je(n), jc(n), se(n), sc(n)
        real*8 :: mu0(n), mu1(n), etas(n), Ilcm
        complex*16 :: sdmatjj(n,n)
        nn = n
        Ilcm = min(IP12, IP32)
        L  = (/2, 1, 1, 1, 1/)
        S = (/1, 0, 1, 1, 0/)
        lc = (/1, 1, 1, 1, 1/);
        le = (/2, 2, 2, 0, 0/);
        je = (/3/2., 5/2., 3/2., 1/2., 1/2./)
        jc = (/1/2., 3/2., 3/2., 1/2., 3/2./)
        etas = 0.0d0 !(/ fnu(Ecm, IP12), fnu(Ecm, IP32), fnu(Ecm, IP32), fnu(Ecm, IP12), fnu(Ecm, IP32) /)
        Jcs = (/1, 2, 1, 1, 1/)
        se = (/1/2., 1/2., 1/2., 1/2., 1/2./)
        sc = (/1/2., 1/2., 1/2., 1/2., 1/2./)
        mu0 = (/ 0.214d0, 0.070d0, 0.500d0, 0.154d0, 0.109d0 /)
        mu1 = (/ 1.22d0, 0.66d0, 1.84d0, -0.28d0, -0.40d0 /)
        Vbaa = reshape((/  1.000d0, 0.001d0,  0.000d0, 0.000d0, 0.001d0, &
                        -0.001d0, 1.000d0, -0.005d0, 0.000d0, 0.000d0, &
                        0.000d0, 0.005d0,  1.000d0, 0.000d0, 0.003d0, &
                        0.000d0, 0.000d0,  0.000d0, 1.000d0, 0.000d0, &
                        -0.001d0, 0.000d0, -0.003d0, 0.000d0, 1.000d0 /),shape(Vbaa), order = (/2,1/))
        call bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        call MQDT_bet(n,Ecm,Uiba,mu0,mu1,Vbaa,(/IP12,IP32,IP32,IP12,IP32/),1,Sdp,Zco,iopen,iclosed,beta)
        do i=1,size(beta)
            beta(i) = PI*(beta(i)-le(iclosed(i)))
        end do
    end subroutine J1_MQDT_bet
    !
    subroutine J2f_MQDT_bet(Ecm, sdp, Zco, iopen, iclosed, beta)
        integer, parameter :: n=3 , J=2
        real*8, parameter :: IP32 = 127109.88d0, IP12= IP32+1432.0d0
        real*8, intent(in) :: Ecm
        integer :: nn, lc(n), le(n), L(n), S(n),  Jcs(n)
        real :: sc(n), jc(n), je(n), K(n), se(n)
        real*8 :: Uiba(n,n), Vbaa(n,n)
        real*8 :: mu0(n), mu1(n), etas(n), Ilcm
        complex*16, allocatable,intent(out) :: sdp(:,:), Zco(:,:)
        integer, allocatable, intent(out) :: iopen(:), iclosed(:)
        real*8, ALLOCATABLE, INTENT(OUT) :: beta(:)

        Ilcm = min(IP32, IP12)
        nn = n
        lc= (/ 1, 1, 1 /)
        le =(/ 3, 3, 3 /)
        jc = (/ 3./2., 3./2., 1./2. /)
        etas = 0.0d0 !(/ fnu(Ecm, IP32), fnu(Ecm,IP32), fnu(Ecm,IP12)/)
        mu0 = (/ 0.022d0,-0.001d0,0.020d0 /)
        mu1 = (/-0.045d0,-0.043d0,-0.060d0/)
        je =(/ 5./2., 7./2., 5./2. /)
        K = (/ 3./2., 5./2., 5./2. /)
        Jcs = (/ 2, 1, 1 /)
        L =(/ 2, 3, 2 /)
        S = (/ 1, 1, 0 /)
        sc = (/ 1./2., 1./2., 1./2. /)
        se = (/ 1./2., 1./2., 1./2. /)
        Vbaa = reshape((/1.0d0, 0.d0, 0.d0, & 
                        0.0d0, 1.d0, 0.d0, &
                        0.0d0, 0.d0, 1.d0 /),shape(Vbaa))
        call bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        call MQDT_bet(n,Ecm,Uiba,mu0,mu1,Vbaa,(/IP32,IP32,IP12/),1,Sdp,Zco,iopen,iclosed,beta)
        do i=1,size(beta)
            beta(i) = PI*(beta(i)-le(iclosed(i)))
        end do
    end subroutine J2f_MQDT_bet
    !
    subroutine J2p_MQDT_bet(Ecm, sdp, Zco, iopen, iclosed, beta)
        integer, parameter :: n=3 , J=2
        real*8, parameter :: IP32 = 127109.88d0, IP12= IP32+1431.41d0
        real*8, intent(in) :: Ecm
        integer :: nn, lc(n), le(n), L(n), S(n),  Jcs(n)
        real :: sc(n), jc(n), je(n), K(n), se(n)
        real*8 :: Uiba(n,n), Vbaa(n,n)
        real*8 :: mu0(n), mu1(n), etas(n), Ilcm
        complex*16, allocatable, intent(out) :: sdp(:,:), Zco(:,:)
        integer, allocatable, intent(out) :: iopen(:), iclosed(:)
        real*8, allocatable, intent(out) :: beta(:)
        Ilcm = min(IP32, IP12)
        nn = n
        lc= (/ 1, 1, 1 /)
        le =(/ 1, 1, 1 /)
        jc = (/ 3./2., 3./2., 1./2. /)
        etas = 0.0d0 !(/ fnu(Ecm, IP32), fnu(Ecm,IP32), fnu(Ecm,IP12)/)
        mu0 = (/ 0.649d0,0.668d0,0.694d0 /)
        mu1 = (/0.239d0,0.243d0,0.301d0/)
        je =(/ 3./2., 1./2., 3./2. /)
        K = (/ 3./2., 5./2., 3./2. /)
        Jcs = (/ 2, 1, 1 /)
        L =(/ 1, 2, 2 /)
        S = (/ 1, 0, 1 /)
        sc = (/ 1./2., 1./2., 1./2. /)
        se = (/ 1./2., 1./2., 1./2. /)
        Vbaa = reshape((/0.9985079636885027d0,-0.060044423946200465d0,-0.00016432773856367433d0, & 
                        0.06040476384223767d0,0.9981651117253054d0,0.00023239451655082677d0, &
                        -0.0003331002905561653d0,0.0002772857325016309d0,0.999564364605493d0 /),shape(Vbaa))
        call bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        call MQDT_bet(n,Ecm,Uiba,mu0,mu1,Vbaa,(/IP32,IP32,IP12/),1,Sdp,Zco,iopen, iclosed, beta)
        do i=1,size(beta)
            beta(i) = PI*(beta(i)-le(iclosed(i)))
        end do
    end subroutine J2p_MQDT_bet
    !
    subroutine J3_MQDT_bet(Ecm, sdp, Zco, iopen, iclosed, beta)
        integer, parameter :: n=6, J=3
        real*8, parameter :: IP32 = 127109.9d0, IP12= IP32+1432.0d0
        integer :: nn
        real*8, intent(in) :: Ecm
        complex*16, allocatable,intent(out) :: sdp(:,:), Zco(:,:)
        integer, allocatable, intent(out) :: iopen(:), iclosed(:)
        real*8,allocatable, INTENT(OUT) :: beta(:)
        real*8 :: Uiba(n,n), Ujji(n,n), th13(n,n), th23(n,n), Vbaa(n,n)
        integer :: L(n), S(n), lc(n), le(n), Jcs(n)
        real :: je(n), jc(n), se(n), sc(n)
        real*8 :: mu0(n), mu1(n), etas(n), Ilcm

        nn = n
        Ilcm = min(IP12, IP32)
        L  = (/3, 3, 2, 3, 3, 4/)
        S = (/0, 1, 1, 0, 1, 1/)
        lc = (/1, 1, 1, 1, 1, 1/)
        le = (/2, 2, 2, 4, 4, 4/)
        je = (/ 5/2., 3/2., 5/2., 7/2., 9/2., 7/2./)
        jc = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        etas = 0.0d0 !(/ fnu(Ecm, IP12), fnu(Ecm, IP32), fnu(Ecm, IP32), fnu(Ecm, IP12), fnu(Ecm, IP32) /)
        Jcs = (/ 1, 2, 1, 1, 2, 1/)
        se = (/1/2., 1/2., 1/2., 1/2., 1/2., 1/2./)
        sc = (/1/2., 1/2., 1/2., 1/2., 1/2., 1/2./)
        mu0 = (/ 0.2663d0, 0.3924d0, 0.2325d0, 0.0d0, 0.0d0, 0.0d0 /)
        mu1 = 0.0d0
       th13 = reshape((/ COS(0.071d0), 0.0d0, SIN(0.071d0),  0.0d0,  0.0d0,  0.0d0, &
                                0.00d0, 1.0d0,       0.00d0,  0.0d0,  0.0d0,  0.0d0, &
                         -SIN(0.071d0), 0.0d0,  COS(0.071d0), 0.0d0,  0.0d0,  0.0d0, &
                                0.0d0,  0.0d0,        0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                                0.0d0,  0.0d0,        0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                                0.0d0,  0.0d0,        0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th13), order = (/2,1/))
        
        th23 = reshape((/       1.00d0,        0.0d0,         1.0d0,  0.0d0,  0.0d0,  0.0d0, &
                                0.00d0, COS(0.047d0), -SIN(0.047d0),  0.0d0,  0.0d0,  0.0d0, &
                                0.00d0, SIN(0.047d0),  COS(0.047d0),  0.0d0,  0.0d0,  0.0d0, &
                                0.0d0,         0.0d0,         0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                                0.0d0,         0.0d0,         0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                                0.0d0,         0.0d0,         0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th23), order = (/2,1/))
        Vbaa = matmul(th13, th23)
        call bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        call MQDT_bet(n,Ecm,Uiba,mu0,mu1,Vbaa,(/IP32,IP32,IP12,IP32,IP12/),1,Sdp,Zco, iopen, iclosed,beta)
        do i=1,size(beta)
            beta(i) = PI*(beta(i)-le(iclosed(i)))
        end do
    end subroutine J3_MQDT_bet
    !
    subroutine J4_MQDT_bet(Ecm, sdp, Zco, iopen, iclosed, beta)
        integer, parameter :: n=6 , J=4
        real*8, parameter :: IP32 = 127109.88d0, IP12= IP32+1431.41d0
        real*8, intent(in) :: Ecm
        complex*16, allocatable, intent(out) :: sdp(:,:), Zco(:,:)
        integer, allocatable, intent(out) :: iopen(:), iclosed(:)
        real*8, allocatable, intent(out) :: beta(:)
        integer :: nn, lc(n), le(n), L(n), S(n),  Jcs(n)
        real :: sc(n), jc(n), je(n), se(n)
        real*8 :: Uiba(n,n), Uip(n,n), th12(n,n), th13(n,n), th23(n,n), Vbaa(n,n)
        real*8 :: thp(n,n)
        real*8 :: mu0(n), mu1(n), Ilcm, A(n)
        real*8 :: thsf(3), thsf1(3)
        nn = n
        Ilcm = min(IP12, IP32)
        Uip = 0.d0
        L  = (/4, 4, 3, 4, 4, 5/)
        S =  (/0, 1, 1, 0, 1, 1/)
        lc = (/1, 1, 1, 1, 1, 1/)
        le = (/3, 3, 3, 5, 5, 5/)
        je = (/ 5/2., 7/2., 7/2., 9/2., 11/2., 9/2./)
        jc = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        Jcs = (/ 1, 2, 1, 1, 2, 1/)
        se = (/1/2., 1/2., 1/2., 1/2., 1/2., 1/2./)
        sc = (/1/2., 1/2., 1/2., 1/2., 1/2., 1/2./)
        mu0 = (/ 0.00992992d0,  0.02038345d0, -0.00285065d0, 0.0d0, 0.0d0, 0.0d0 /)
        mu1 = 0.d0
        thsf =  (/-0.03702163d0,  0.03709483d0,  0.03758697d0/)
        thsf1 = 0.d0

        th12 = reshape((/ COS(thsf(1)), -SIN(thsf(1)),  0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                          SIN(thsf(1)),  COS(thsf(1)),  0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                                0.0d0 ,         0.0d0,  1.0d0,  0.0d0,  0.0d0,  0.0d0, &
                                0.0d0,          0.0d0,  0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                                0.0d0,          0.0d0,  0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                                0.0d0,          0.0d0,  0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th12), order = (/2,1/))

        th13 = reshape((/ COS(thsf(2)), 0.0d0, -SIN(thsf(2)),  0.0d0,  0.0d0,  0.0d0, &
                                0.00d0, 1.0d0,        0.00d0,  0.0d0,  0.0d0,  0.0d0, &
                          SIN(thsf(2)), 0.0d0,  COS(thsf(2)),  0.0d0,  0.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th13), order = (/2,1/))
        
        th23 = reshape((/1.00d0,        0.0d0,         0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                         0.00d0, COS(thsf(3)), -SIN(thsf(3)),  0.0d0,  0.0d0,  0.0d0, &
                         0.00d0, SIN(thsf(3)),  COS(thsf(3)),  0.0d0,  0.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th23), order = (/2,1/))
        Vbaa = matmul(th12,matmul(th13, th23))
        
        Uip = 0.d0
        thp = reshape((/ SIN(thsf(1)),  COS(thsf(1)),  0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                        -COS(thsf(1)),  SIN(thsf(1)),  0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                               0.0d0 ,         0.0d0,  1.0d0,  0.0d0,  0.0d0,  0.0d0, &
                               0.0d0,          0.0d0,  0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                               0.0d0,          0.0d0,  0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                               0.0d0,          0.0d0,  0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th12), order = (/2,1/))
        Uip = Uip - 2.0d0*thsf1(1)*matmul(thp, matmul(th13, th23))
        
        thp = reshape((/  SIN(thsf(2)), 0.0d0,  COS(thsf(2)),  0.0d0,  0.0d0,  0.0d0, &
                                0.00d0, 1.0d0,        0.00d0,  0.0d0,  0.0d0,  0.0d0, &
                         -COS(thsf(2)), 0.0d0,  SIN(thsf(2)),  0.0d0,  0.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                                 0.0d0, 0.0d0,         0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th13), order = (/2,1/))
        Uip = Uip - 2.0d0*thsf1(2)*matmul(th12, matmul(thp, th23))

        thp = reshape((/ 1.00d0,        0.0d0,         0.0d0,  0.0d0,  0.0d0,  0.0d0, &
                         0.00d0, SIN(thsf(3)),  COS(thsf(3)),  0.0d0,  0.0d0,  0.0d0, &
                         0.00d0,-COS(thsf(3)),  SIN(thsf(3)),  0.0d0,  0.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                         0.0d0,         0.0d0,         0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th23), order = (/2,1/))
        Uip = Uip -2.0d0*thsf1(3)*matmul(th12, matmul(th13, thp))
        call MQDT_bet(n,Ecm,Uiba,mu0,mu1,Vbaa,(/IP32, IP32, IP12, IP32, IP32, IP12/),1,Sdp,Zco,iopen,iclosed,beta)
        do i=1,size(beta)
            beta(i) = PI*(beta(i)-le(iclosed(i)))
        end do
    end subroutine
    !
    subroutine J0_MQDT_at(Ecm, sdmat, kmat, ip)
        integer, parameter :: n=2, J=0
        real*8, parameter :: IP32 = 127109.9d0, IP12= IP32+1432.0d0
        integer :: nn
        integer, intent(in) :: ip
        real*8, intent(in) :: Ecm
        complex*16, intent(out) :: sdmat(n,n)
        complex*16 :: sdmatjj(n,n)
        real*8 :: kmat(n,n)
        real*8 :: Ujji(n,n),Uiba(n,n), Vbaa(n,n)
        integer :: L(n), S(n), lc(n), le(n), Jcs(n)
        real :: je(n), jc(n), se(n), sc(n)
        real*8 :: mu0(n), mu1(n), etas(n), Il
        real*8 :: cang
        nn = n
        Il = 124000.0d0 !min(IP32, IP12)
        L  = (/1, 0/)
        S = (/1, 0/)
        lc = (/1, 1/);
        le = (/1,1/);
        !je = (/3/2., 5/2., 3/2., 1/2., 1/2./)
        jc = (/3/2., 1/2./)
        etas = 0.0d0 !(/ fnu(Ecm, IP32), fnu(Ecm, IP12) /)
        Jcs = (/1, 1/)
        se = (/1/2., 1/2./)
        sc = (/1/2., 1/2./) 
        mu0 = (/ 0.70d0, 0.52d0 /)
        mu1 = (/ 0.103d0, 0.1832d0 /)
        cang = 0.02076d0-((Ecm-124.0d3)/RydAr)*(-0.1185d0)
        Vbaa = reshape((/COS(cang), SIN(cang),&
                  -1.0d0*SIN(cang), COS(cang) /),shape(Vbaa), order = (/2,1/))
        call bjcs_kLS(nn, lc, sc, jc, le, jcs, se, L, S, J, Uiba)
        call MQDT_at(nn, Ecm, UIba, mu0, mu1, etas, Vbaa, Il,0,Sdmat,kmat)
        call bjj_kJcs(nn, lc, sc, jc, le, real(Jcs), se, je, J, Ujji)
        if(ip.eq.1) then
            sdmatjj = matmul(matmul(Ujji, Sdmat),TRANSPOSE(Ujji))
            print ('(A,f22.9)'), "Sdmatrix for energy E(cm-1): ", Ecm
            call zwrite_matrix(n,n,sdmatjj)
        end if
    end subroutine
    !
    subroutine J0_MQDT_params(Ecm, Uia, mua)
        integer, parameter :: n=2, J=0
        real*8, parameter :: IP32 = 127109.9d0, IP12= IP32+1432.0d0
        integer :: nn
        real*8, intent(in) :: Ecm
        real*8, intent(out) :: uia(n,n), mua(n)
        real*8 :: Uiba(n,n), Vbaa(n,n)
        integer :: L(n), S(n), lc(n), le(n), Jcs(n)
        real :: je(n), jc(n), se(n), sc(n)
        real*8 :: mu0(n), mu1(n), etas(n), Il
        real*8 :: cang
        nn = n
        Il = 124000.0d0 !min(IP32, IP12)
        L  = (/1, 0/)
        S = (/1, 0/)
        lc = (/1, 1/);
        le = (/1,1/);
        je = (/3/2., 1/2./)
        jc = (/3/2., 1/2./)
        etas = 0.0d0 !(/ fnu(Ecm, IP32), fnu(Ecm, IP12) /)
        Jcs = (/1, 1/)
        se = (/1/2., 1/2./)
        sc = (/1/2., 1/2./) 
        mu0 = (/ 0.70d0, 0.52d0 /)
        mu1 = (/ 0.103d0, 0.1832d0 /)
        cang = 0.02076d0-((Ecm-124.0d3)/RydAr)*(-0.1185d0)
        Vbaa = reshape((/COS(cang), SIN(cang),&
                  -1.0d0*SIN(cang), COS(cang) /),shape(Vbaa), order = (/2,1/))
        call bjcs_kLS(nn, lc, sc, jc, le, jcs, se, L, S, J, Uiba)
        Uia = matmul(Uiba, Vbaa)
        mua = mu0 + mu0+mu1*0.5d0*((Ecm-Il)/RydAr)
    end subroutine J0_MQDT_params
    !
    subroutine J1_MQDT_at(Ecm, sdmat, kmat, ip)
        integer, parameter :: n=5, J=1
        real*8, parameter :: IP32 = 127109.9d0, IP12= IP32+1432.0d0
        integer :: nn
        integer, intent(in) :: ip
        real*8, intent(in) :: Ecm
        complex*16, intent(out) :: sdmat(n,n)
        real*8, intent(out) :: kmat(n,n)
        real*8 :: Uiba(n,n), Ujji(n,n), Vbaa(n,n)
        integer :: L(n), S(n), lc(n), le(n), Jcs(n)
        real :: je(n), jc(n), se(n), sc(n)
        real*8 :: mu0(n), mu1(n), etas(n), Ilcm
        complex*16 :: sdmatjj(n,n)
        nn = n
        Ilcm = min(IP12, IP32)
        L  = (/2, 1, 1, 1, 1/)
        S = (/1, 0, 1, 1, 0/)
        lc = (/1, 1, 1, 1, 1/);
        le = (/2, 2, 2, 0, 0/);
        je = (/3/2., 5/2., 3/2., 1/2., 1/2./)
        jc = (/1/2., 3/2., 3/2., 1/2., 3/2./)
        etas = 0.0d0 !(/ fnu(Ecm, IP12), fnu(Ecm, IP32), fnu(Ecm, IP32), fnu(Ecm, IP12), fnu(Ecm, IP32) /)
        Jcs = (/1, 2, 1, 1, 1/)
        se = (/1/2., 1/2., 1/2., 1/2., 1/2./)
        sc = (/1/2., 1/2., 1/2., 1/2., 1/2./)
        mu0 = (/ 0.214d0, 0.070d0, 0.500d0, 0.154d0, 0.109d0 /)
        mu1 = (/ 1.22d0, 0.66d0, 1.84d0, -0.28d0, -0.40d0 /)
        Vbaa = reshape((/  1.000d0, 0.001d0,  0.000d0, 0.000d0, 0.001d0, &
                        -0.001d0, 1.000d0, -0.005d0, 0.000d0, 0.000d0, &
                        0.000d0, 0.005d0,  1.000d0, 0.000d0, 0.003d0, &
                        0.000d0, 0.000d0,  0.000d0, 1.000d0, 0.000d0, &
                        -0.001d0, 0.000d0, -0.003d0, 0.000d0, 1.000d0 /),shape(Vbaa), order = (/2,1/))
        call bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        call MQDT_at(n, Ecm, Uiba, mu0, mu1, etas, Vbaa, Ilcm, 1, Sdmat, kmat)
        call bjj_kJcs(n, lc, sc, jc, le, real(Jcs), se, je, J, Ujji)
        if(ip.eq.1) then
            sdmatjj = matmul(matmul(Ujji, Sdmat),TRANSPOSE(Ujji))
            print ('(A,f22.9)'), "Sdmatrix for energy E(cm-1): ", Ecm
            call zwrite_matrix(n,n,sdmatjj)
        end if
    end subroutine J1_MQDT_at
    !
    subroutine J1_MQDT_params(Ecm, Uia, mua)
        integer, parameter :: n=5, J=1
        real*8, parameter :: IP32 = 127109.9d0, IP12= IP32+1432.0d0
        integer :: nn
        real*8, intent(in) :: Ecm
        real*8, intent(out) :: Uia(n,n), mua(n)
        real*8 :: Uiba(n,n), Vbaa(n,n)
        integer :: L(n), S(n), lc(n), le(n), Jcs(n)
        real :: je(n), jc(n), se(n), sc(n)
        real*8 :: mu0(n), mu1(n), etas(n), Ilcm
        nn = n
        Ilcm = min(IP12, IP32)
        L  = (/2, 1, 1, 1, 1/)
        S = (/1, 0, 1, 1, 0/)
        lc = (/1, 1, 1, 1, 1/);
        le = (/2, 2, 2, 0, 0/);
        je = (/3/2., 5/2., 3/2., 1/2., 1/2./)
        jc = (/1/2., 3/2., 3/2., 1/2., 3/2./)
        etas = 0.0d0 !(/ fnu(Ecm, IP12), fnu(Ecm, IP32), fnu(Ecm, IP32), fnu(Ecm, IP12), fnu(Ecm, IP32) /)
        Jcs = (/1, 2, 1, 1, 1/)
        se = (/1/2., 1/2., 1/2., 1/2., 1/2./)
        sc = (/1/2., 1/2., 1/2., 1/2., 1/2./)
        mu0 = (/ 0.214d0, 0.070d0, 0.500d0, 0.154d0, 0.109d0 /)
        mu1 = (/ 1.22d0, 0.66d0, 1.84d0, -0.28d0, -0.40d0 /)
        Vbaa = reshape((/  1.000d0, 0.001d0,  0.000d0, 0.000d0, 0.001d0, &
                        -0.001d0, 1.000d0, -0.005d0, 0.000d0, 0.000d0, &
                        0.000d0, 0.005d0,  1.000d0, 0.000d0, 0.003d0, &
                        0.000d0, 0.000d0,  0.000d0, 1.000d0, 0.000d0, &
                        -0.001d0, 0.000d0, -0.003d0, 0.000d0, 1.000d0 /),shape(Vbaa), order = (/2,1/))
        call bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        Uia = matmul(Uiba, Vbaa)
        mua =  mu0+mu1*1.0d0*((Ecm-Ilcm)/RydAr)
    end subroutine J1_MQDT_params
    !
    subroutine J2f_MQDT_at(Ecm, sdmat, kmat, ip)
        integer, parameter :: n=3 , J=2
        real*8, parameter :: IP32 = 127109.88d0, IP12= IP32+1431.41d0
        real*8, intent(in) :: Ecm
        integer :: ip
        complex*16 :: sdmatjj(n,n)
        integer :: nn, lc(n), le(n), L(n), S(n),  Jcs(n)
        real :: sc(n), jc(n), je(n), K(n), se(n)
        real*8 :: Uiba(n,n), Vbaa(n,n), Ujji(n,n)
        real*8 :: mu0(n), mu1(n), etas(n), Ilcm
        complex*16, intent(out) :: sdmat(n,n)
        real*8, intent(out) :: kmat(n,n)
        Ilcm = min(IP32, IP12)
        nn = n
        lc= (/ 1, 1, 1 /)
        le =(/ 3, 3, 3 /)
        jc = (/ 3./2., 3./2., 1./2. /)
        etas = 0.0d0 !(/ fnu(Ecm, IP32), fnu(Ecm,IP32), fnu(Ecm,IP12)/)
        mu0 = (/ 0.022d0,-0.001d0,0.020d0 /)
        mu1 = (/-0.045d0,-0.043d0,-0.060d0/)
        je =(/ 5./2., 7./2., 5./2. /)
        K = (/ 3./2., 5./2., 5./2. /)
        Jcs = (/ 2, 1, 1 /)
        L =(/ 2, 3, 2 /)
        S = (/ 1, 1, 0 /)
        sc = (/ 1./2., 1./2., 1./2. /)
        se = (/ 1./2., 1./2., 1./2. /)
        Vbaa = reshape((/1.0d0, 0.d0, 0.d0, & 
                        0.0d0, 1.d0, 0.d0, &
                        0.0d0, 0.d0, 1.d0 /),shape(Vbaa))
        call bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        call MQDT_at(n, Ecm, Uiba, mu0, mu1, etas, Vbaa, Ilcm, 1, Sdmat, Kmat)
        call bjj_kJcs(n, lc, sc, jc, le, real(Jcs), se, je, J, Ujji)
        if(ip.eq.1) then
            sdmatjj = matmul(matmul(Ujji, Sdmat),TRANSPOSE(Ujji))
            print ('(A,f22.9)'), "Sdmatrix for energy E(cm-1): ", Ecm
            call zwrite_matrix(n,n,sdmatjj)
        end if
    end subroutine J2f_MQDT_at
    !
    subroutine J2f_MQDT_params(Ecm, Uia, mua)
        integer, parameter :: n=3 , J=2
        real*8, parameter :: IP32 = 127109.88d0, IP12= IP32+1431.41d0
        real*8, intent(in) :: Ecm
        integer :: nn, lc(n), le(n), L(n), S(n),  Jcs(n)
        real :: sc(n), jc(n), je(n), K(n), se(n)
        real*8 :: Uiba(n,n), Vbaa(n,n)
        real*8 :: mu0(n), mu1(n), etas(n), Ilcm
        real*8, intent(out) :: Uia(n,n), mua(n)
        Ilcm = min(IP32, IP12)
        nn = n
        lc= (/ 1, 1, 1 /)
        le =(/ 3, 3, 3 /)
        jc = (/ 3./2., 3./2., 1./2. /)
        etas = 0.0d0 !(/ fnu(Ecm, IP32), fnu(Ecm,IP32), fnu(Ecm,IP12)/)
        mu0 = (/ 0.022d0,-0.001d0,0.020d0 /)
        mu1 = (/-0.045d0,-0.043d0,-0.060d0/)
        je =(/ 5./2., 7./2., 5./2. /)
        K = (/ 3./2., 5./2., 5./2. /)
        Jcs = (/ 2, 1, 1 /)
        L =(/ 2, 3, 2 /)
        S = (/ 1, 1, 0 /)
        sc = (/ 1./2., 1./2., 1./2. /)
        se = (/ 1./2., 1./2., 1./2. /)
        Vbaa = reshape((/1.0d0, 0.d0, 0.d0, & 
                        0.0d0, 1.d0, 0.d0, &
                        0.0d0, 0.d0, 1.d0 /),shape(Vbaa))
        call bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        Uia = matmul(Uiba, Vbaa)
        mua =mu0+mu1*1.0d0*((Ecm-Ilcm)/RydAr)
    end subroutine J2f_MQDT_params 
    !
    subroutine J2p_MQDT_at(Ecm, sdmat, kmat, ip)
        integer, parameter :: n=3 , J=2
        real*8, parameter :: IP32 = 127109.88d0, IP12= IP32+1431.41d0
        real*8, intent(in) :: Ecm
        integer :: ip
        complex*16 :: sdmatjj(n,n)
        integer :: nn, lc(n), le(n), L(n), S(n),  Jcs(n)
        real :: sc(n), jc(n), je(n), K(n), se(n)
        real*8 :: Uiba(n,n), Vbaa(n,n), Ujji(n,n)
        real*8 :: mu0(n), mu1(n), etas(n), Ilcm
        complex*16, intent(out) :: sdmat(n,n)
        real*8, intent(out) :: kmat(n,n)
        Ilcm = min(IP32, IP12)
        nn = n
        lc= (/ 1, 1, 1 /)
        le =(/ 1, 1, 1 /)
        jc = (/ 3./2., 3./2., 1./2. /)
        etas = 0.0d0 !(/ fnu(Ecm, IP32), fnu(Ecm,IP32), fnu(Ecm,IP12)/)
        mu0 = (/ 0.649d0,0.668d0,0.694d0 /)
        mu1 = (/0.239d0,0.243d0,0.301d0/)
        je =(/ 3./2., 1./2., 3./2. /)
        K = (/ 3./2., 5./2., 3./2. /)
        Jcs = (/ 2, 1, 1 /)
        L =(/ 1, 2, 2 /)
        S = (/ 1, 0, 1 /)
        sc = (/ 1./2., 1./2., 1./2. /)
        se = (/ 1./2., 1./2., 1./2. /)
        Vbaa = reshape((/0.9985079636885027d0,-0.060044423946200465d0,-0.00016432773856367433d0, & 
                        0.06040476384223767d0,0.9981651117253054d0,0.00023239451655082677d0, &
                        -0.0003331002905561653d0,0.0002772857325016309d0,0.999564364605493d0 /),shape(Vbaa))
        call bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        call MQDT_at(n, Ecm, Uiba, mu0, mu1, etas, Vbaa, Ilcm, 1, Sdmat, Kmat)
        call bjj_kJcs(n, lc, sc, jc, le, real(Jcs), se, je, J, Ujji)
        if(ip.eq.1) then
            sdmatjj = matmul(matmul(Ujji, Sdmat),TRANSPOSE(Ujji))
            print ('(A,f22.9)'), "Sdmatrix for energy E(cm-1): ", Ecm
            call zwrite_matrix(n,n,sdmatjj)
        end if
    end subroutine J2p_MQDT_at
    !
    subroutine J2p_MQDT_params(Ecm, Uia, mua)
        integer, parameter :: n=3 , J=2
        real*8, parameter :: IP32 = 127109.88d0, IP12= IP32+1431.41d0
        real*8, intent(in) :: Ecm
        integer :: nn, lc(n), le(n), L(n), S(n),  Jcs(n)
        real :: sc(n), jc(n), je(n), K(n), se(n)
        real*8 :: Uiba(n,n), Vbaa(n,n)
        real*8 :: mu0(n), mu1(n), etas(n), Ilcm
        real*8, intent(out) :: Uia(n,n), mua(n)
        Ilcm = min(IP32, IP12)
        nn = n
        lc= (/ 1, 1, 1 /)
        le =(/ 1, 1, 1 /)
        jc = (/ 3./2., 3./2., 1./2. /)
        etas = 0.0d0 !(/ fnu(Ecm, IP32), fnu(Ecm,IP32), fnu(Ecm,IP12)/)
        mu0 = (/ 0.649d0,0.668d0,0.694d0 /)
        mu1 = (/0.239d0,0.243d0,0.301d0/)
        je =(/ 3./2., 1./2., 3./2. /)
        K = (/ 3./2., 5./2., 3./2. /)
        Jcs = (/ 2, 1, 1 /)
        L =(/ 1, 2, 2 /)
        S = (/ 1, 0, 1 /)
        sc = (/ 1./2., 1./2., 1./2. /)
        se = (/ 1./2., 1./2., 1./2. /)
        Vbaa = reshape((/0.9985079636885027d0,-0.060044423946200465d0,-0.00016432773856367433d0, & 
                        0.06040476384223767d0,0.9981651117253054d0,0.00023239451655082677d0, &
                        -0.0003331002905561653d0,0.0002772857325016309d0,0.999564364605493d0 /),shape(Vbaa))
        call bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        Uia = matmul(Uiba, Vbaa)
        mua = mu0+mu1*1.0d0*((Ecm-Ilcm)/RydAr)
    end subroutine J2p_MQDT_params
    !
    subroutine J3_MQDT_at(Ecm, sdmat, kmat, ip)
        integer, parameter :: n=6, J=3
        real*8, parameter :: IP32 = 127109.9d0, IP12= IP32+1432.0d0
        integer :: nn
        integer, intent(in) :: ip
        real*8, intent(in) :: Ecm
        complex*16, intent(out) :: sdmat(n,n)
        real*8 :: Uiba(n,n), Ujji(n,n), th13(n,n), th23(n,n), Vbaa(n,n)
        integer :: L(n), S(n), lc(n), le(n), Jcs(n)
        real :: je(n), jc(n), se(n), sc(n)
        real*8 :: mu0(n), mu1(n), etas(n), Ilcm
        complex*16 :: sdmatjj(n,n)
        real*8 :: kmat(n,n)
        nn = n
        Ilcm = min(IP12, IP32)
        L  = (/3, 3, 2, 3, 3, 4/)
        S = (/0, 1, 1, 0, 1, 1/)
        lc = (/1, 1, 1, 1, 1, 1/)
        le = (/2, 2, 2, 4, 4, 4/)
        je = (/ 5/2., 3/2., 5/2., 7/2., 9/2., 7/2./)
        jc = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        etas = 0.0d0 !(/ fnu(Ecm, IP12), fnu(Ecm, IP32), fnu(Ecm, IP32), fnu(Ecm, IP12), fnu(Ecm, IP32) /)
        Jcs = (/ 1, 2, 1, 1, 2, 1/)
        se = (/1/2., 1/2., 1/2., 1/2., 1/2., 1/2./)
        sc = (/1/2., 1/2., 1/2., 1/2., 1/2., 1/2./)
        mu0 = (/ 0.2663d0, 0.3924d0, 0.2325d0, 0.0d0, 0.0d0, 0.0d0 /)
        mu1 = 0.0d0
       th13 = reshape((/ COS(0.071d0), 0.0d0, SIN(0.071d0),  0.0d0,  0.0d0,  0.0d0, &
                                0.00d0, 1.0d0,       0.00d0,  0.0d0,  0.0d0,  0.0d0, &
                         -SIN(0.071d0), 0.0d0,  COS(0.071d0), 0.0d0,  0.0d0,  0.0d0, &
                                0.0d0,  0.0d0,        0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                                0.0d0,  0.0d0,        0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                                0.0d0,  0.0d0,        0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th13), order = (/2,1/))
        
        th23 = reshape((/       1.00d0,        0.0d0,         1.0d0,  0.0d0,  0.0d0,  0.0d0, &
                                0.00d0, COS(0.047d0), -SIN(0.047d0),  0.0d0,  0.0d0,  0.0d0, &
                                0.00d0, SIN(0.047d0),  COS(0.047d0),  0.0d0,  0.0d0,  0.0d0, &
                                0.0d0,         0.0d0,         0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                                0.0d0,         0.0d0,         0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                                0.0d0,         0.0d0,         0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th23), order = (/2,1/))
        Vbaa = matmul(th13, th23)
        call bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        call MQDT_at(n, Ecm, Uiba, mu0, mu1, etas, Vbaa, Ilcm, 1, Sdmat, Kmat)
        call bjj_kJcs(n, lc, sc, jc, le, real(Jcs), se, je, J, Ujji)
        if(ip.eq.1) then
            sdmatjj = matmul(matmul(Ujji, Sdmat),TRANSPOSE(Ujji))
            print ('(A,f22.9)'), "Sdmatrix for energy E(cm-1): ", Ecm
            call zwrite_matrix(n,n,sdmatjj)
        end if
    end subroutine J3_MQDT_at
    !
    subroutine J3_MQDT_params(Ecm, Uia, mua)
        integer, parameter :: n=6, J=3
        real*8, parameter :: IP32 = 127109.9d0, IP12= IP32+1432.0d0
        integer :: nn
        real*8, intent(in) :: Ecm
        real*8, intent(out) :: Uia(n,n), mua(n)
        real*8 :: Uiba(n,n), th13(n,n), th23(n,n), Vbaa(n,n)
        integer :: L(n), S(n), lc(n), le(n), Jcs(n)
        real :: je(n), jc(n), se(n), sc(n)
        real*8 :: mu0(n), mu1(n), etas(n), Ilcm
        nn = n
        Ilcm = min(IP12, IP32)
        L  = (/3, 3, 2, 3, 3, 4/)
        S = (/0, 1, 1, 0, 1, 1/)
        lc = (/1, 1, 1, 1, 1, 1/)
        le = (/2, 2, 2, 4, 4, 4/)
        je = (/ 5/2., 3/2., 5/2., 7/2., 9/2., 7/2./)
        jc = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
        etas = 0.0d0 !(/ fnu(Ecm, IP12), fnu(Ecm, IP32), fnu(Ecm, IP32), fnu(Ecm, IP12), fnu(Ecm, IP32) /)
        Jcs = (/ 1, 2, 1, 1, 2, 1/)
        se = (/1/2., 1/2., 1/2., 1/2., 1/2., 1/2./)
        sc = (/1/2., 1/2., 1/2., 1/2., 1/2., 1/2./)
        mu0 = (/ 0.2663d0, 0.3924d0, 0.2325d0, 0.0d0, 0.0d0, 0.0d0 /)
        mu1 = 0.00d0
        th13 = reshape((/ COS(0.071d0), 0.0d0, SIN(0.071d0),  0.0d0,  0.0d0,  0.0d0, &
                                0.00d0, 1.0d0,       0.00d0,  0.0d0,  0.0d0,  0.0d0, &
                         -SIN(0.071d0), 0.0d0,  COS(0.071d0), 0.0d0,  0.0d0,  0.0d0, &
                                0.0d0,  0.0d0,        0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                                0.0d0,  0.0d0,        0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                                0.0d0,  0.0d0,        0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th13), order = (/2,1/))
        
        th23 = reshape((/       1.00d0,        0.0d0,         1.0d0,  0.0d0,  0.0d0,  0.0d0, &
                                0.00d0, COS(0.047d0), -SIN(0.047d0),  0.0d0,  0.0d0,  0.0d0, &
                                0.00d0, SIN(0.047d0),  COS(0.047d0),  0.0d0,  0.0d0,  0.0d0, &
                                0.0d0,         0.0d0,         0.0d0,  1.0d0,  0.0d0,  0.0d0, &
                                0.0d0,         0.0d0,         0.0d0,  0.0d0,  1.0d0,  0.0d0, &
                                0.0d0,         0.0d0,         0.0d0,  0.0d0,  0.0d0,  1.0d0 /),shape(th23), order = (/2,1/))
        Vbaa = matmul(th13, th23)
        call bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        Uia = matmul(Uiba, Vbaa)
        mua = mu0+mu1*1.0d0*((Ecm-Ilcm)/RydAr)
    end subroutine J3_MQDT_params
    !
    subroutine J5_MQDT_at(Ecm, sdmat, kmat, ip)
        integer, parameter :: n=6, J=3
        real*8, parameter :: IP32 = 127109.9d0, IP12= IP32+1432.0d0
        integer :: nn
        integer, intent(in) :: ip
        real*8, intent(in) :: Ecm
        complex*16, intent(out) :: sdmat(n,n)
        real*8 :: Uiba(n,n), Ujji(n,n), Vbaa(n,n)
        integer :: L(n), S(n), lc(n), le(n), Jcs(n)
        real :: je(n), jc(n), se(n), sc(n)
        real*8 :: mu0(n), mu1(n), etas(n), Ilcm
        complex*16 :: sdmatjj(n,n)
        real*8 :: kmat(n,n)
        nn = n
        Ilcm = min(IP12, IP32)
        L  = (/5, 5, 4, 5, 5, 6/)
        S = (/0, 1, 1, 0, 1, 1/)
        lc = (/1, 1, 1, 1, 1, 1/)
        le = (/4, 4, 4, 6, 6, 6/)
        je = (/ 9/2., 7/2., 9/2., 13/2., 11/2., 11/2./)
        jc = (/ 3/2., 3/2., 1/2.,  3/2.,  3/2.,  1/2./)
        etas = 0.0d0 !(/ fnu(Ecm, IP12), fnu(Ecm, IP32), fnu(Ecm, IP32), fnu(Ecm, IP12), fnu(Ecm, IP32) /)
        Jcs = (/ 1, 2, 1, 1, 2, 1/)
        se = (/1/2., 1/2., 1/2., 1/2., 1/2., 1/2./)
        sc = (/1/2., 1/2., 1/2., 1/2., 1/2., 1/2./)
        mu0 = 0.0d0
        mu1 = 0.0d0
        Vbaa = 0.0d0

        do i=1, n
            Vbaa(i,i) = 1.0d0
        end do

        call bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J, Uiba)
        call MQDT_at(n, Ecm, Uiba, mu0, mu1, etas, Vbaa, Ilcm, 1, Sdmat, Kmat)
        call bjj_kJcs(n, lc, sc, jc, le, real(Jcs), se, je, J, Ujji)
        if(ip.eq.1) then
            sdmatjj = matmul(matmul(Ujji, Sdmat),TRANSPOSE(Ujji))
            print ('(A,f22.9)'), "Sdmatrix for energy E(cm-1): ", Ecm
            call zwrite_matrix(n,n,sdmatjj)
        end if
    end subroutine
    !
    subroutine Matij(n, mus, nus, U, mat, Zmat)
        integer, intent(in) :: n
        real*8, intent(in) :: mus(n), nus(n), U(n,n)
        real*8, intent(out) :: mat(n,n), Zmat(n,n)
        real*8 :: d_s_m(n,n), d_c_m(n,n), d_s_n(n,n), d_c_n(n,n)
        mat = 0.0d0
        Zmat = 0.0d0
        d_s_m = 0.0d0
        d_c_m = 0.0d0
        d_s_n = 0.0d0
        d_c_n = 0.0d0
        do i=1,n
            d_s_m(i,i) = SIN(PI*mus(i))
            d_c_m(i,i) = COS(PI*mus(i))
            d_s_n(i,i) = SIN(PI*nus(i))
            d_c_n(i,i) = COS(PI*nus(i))
        end do
        mat = matmul(d_s_n,matmul(U,d_c_m))+matmul(d_c_n,matmul(U,d_s_m))
        Zmat = matmul(d_c_n, matmul(U, d_c_m))-matmul(d_s_n,matmul(U,d_s_m))
        return
    end subroutine Matij
    !
    subroutine Smatij(n,mus,etas, U, Smatd)
        integer, intent(in) :: n
        real*8, intent(in) :: mus(n), U(n,n)
        real*8, intent(in) :: etas(n)
        complex*16 :: DSmatd(n,n), etasm(n,n)
        complex*16, intent(out) :: Smatd(n,n)
        DSmatd = Zer
        Smatd = Zer
        etasm = Zer
        do i=1,n
            DSmatd(i,i) = EXP(-2*PI*Im*mus(i))
            etasm(i,i) = EXP(-Im*etas(i))
        enddo
        Smatd = matmul(U,matmul(DSmatd,TRANSPOSE(U)))
        !TODO: We need to figure out the complex gamma function.
    end subroutine Smatij
    !
    subroutine Kmatij(n, mus, U, Kmat)
        integer, intent(in) :: n
        real*8, intent(in) :: mus(n), U(n,n)
        real*8 :: Tanmu(n,n)
        real*8, intent(out) :: Kmat(n,n)
        Kmat = 0.d0
        Tanmu = 0.d0
        do i=1,n
            Tanmu(i,i) = TAN(PI*mus(i))
        enddo
        Kmat = matmul(U,matmul(Tanmu,TRANSPOSE(U)))
    end subroutine 
    !
    function qsum(J,M,evec) result(y)
        integer :: J,M
        real*8 :: evec(3)
        complex*16 ::  ep, em, eo, y
        real*8, external :: THRJ
        ep = -1.0d0/SQRT(2.0d0) * (evec(1)+Im*evec(2))
        em =  1.0d0/SQRT(2.0d0) * (evec(1)-Im*evec(2))
        eo = Zone*evec(3)
        y = (-1)**(1) * em * THRJ(2*J, 2*1, 2*1, -2*M, 2*1, 0)
        y = (-1)**(0) * eo * THRJ(2*J, 2*1, 2*1, -2*M, 2*0, 0) + y
        y = (-1)**(-1) * ep * THRJ(2*J, 2*1, 2*1, -2*M, 2*(-1), 0) + y
    end function
    !
    function qsum_c(J,M,Jo,Mo,evec) result(y)
        integer :: J,M,Jo,Mo
        real*8 :: evec(3)
        complex*16 ::  ep, em, eo, y
        real*8, external :: THRJ
        ep = -1.0d0/SQRT(2.0d0) * (evec(1)+Im*evec(2))
        em =  1.0d0/SQRT(2.0d0) * (evec(1)-Im*evec(2))
        eo = Zone*evec(3)
        y = (-1)**(1) * em * THRJ(2*J, 2*1, 2*Jo, -2*M, 2*1, 2*Mo)
        y = (-1)**(0) * eo * THRJ(2*J, 2*1, 2*Jo, -2*M, 2*0, 2*Mo) + y
        y = (-1)**(-1) * ep * THRJ(2*J, 2*1, 2*Jo, -2*M, 2*(-1), 2*Mo) + y
    end function
    !
    function fnu(Ecm, Icm)
        real*8, intent(in) :: Ecm, Icm
        real*8 :: fnu
        if(Icm.lt.Ecm) then
            print *,  "Energy above threshold"
            stop
        endif
        fnu = sqrt(RydAr/(Icm-Ecm))
        return 
    end function fnu
    !
    function FourF(gam, w, to) result(y)
        real*8 :: gam, w, to
        complex*16 :: y
        y = EXP(-(gam/2)**2 * w**2)*EXP(Im*w*to)
    end function
    !
    function twoPhotonAmp(ea, e, exi, gam, to, w) result(y)
        real*8 :: ea, e, exi, gam, to, w
        real*8 :: exps
        complex*16 :: y, error
        exps = EXP(-1.0d0/8.0d0 * (gam)**2 * (ea+2.0d0*w-e)**2)
        exps = exps*EXP(-1/2.0 * (gam)**2 * (exi-(ea+w)+1.0d0/4.0d0 * (ea+2*w-e))**2)
        if(exps.lt.1.0d-285) then
            error =  0
        else
            error = (1.0d0-erf_pop(Im*(gam/sqrt(2.0d0) * (exi-(ea+w)+1.0d0/4.0d0*(ea+2.0d0*w-e)))))
        endif
        !print*, error, exps
        if(exps /= exps) then
            print*, "Nan error in exponentials", exps
        endif
        if(error /= error) then
            print*, "Nan error in erf calculations, arg ", Im*(gam/sqrt(2.0d0) * (exi-(ea+w)+1.0d0/4.0d0*(ea+2.0d0*w-e))),&
             "gaussians ", exps
             stop
        endif
        y = exps*error*EXP(Im*(e-ea)*to)
        return 
    end function twoPhotonAmp
    !
    function RetwoPhotonAmp(ea, e, exi, gam, to, w) result(y)
        real*8 :: ea, e, exi, gam, to, w
        real*8 :: exps
        complex*16 :: y
        exps = EXP(-1.0d0/8.0d0 * (gam)**2 * (ea+2.0d0*w-e)**2)
        exps = exps*EXP(-1/2.0 * (gam)**2 * (exi-(ea+w)+1.0d0/4.0d0 * (ea+2*w-e))**2)
        y = exps*EXP(Im*(e-ea)*to)
        return 
    end function RetwoPhotonAmp
    !
    function RedElem(ne, Za, Smatd, in, jce, le, jcse, Jj, regints, iregints) result(y)
        integer :: ne, in, k, j, Jj
        integer :: le(ne), jcse(ne)
        real*8 :: regints(nb,ne), iregints(nb,ne)
        real :: jce(ne)
        real*8 :: Za(nb)
        complex*16 :: Smatd(ne,ne)
        complex*16 :: yt
        complex*16 :: y
        real*8, external :: SIXJ
        y = Zer
        do k=1, ne
            do j=1, nb
                yt = Za(j) * 0.5d0 * kronij(int(2*jce(k)), int(2*jcb(j))) * kronij(int(jcsb(j)), int(jcse(k)))
                yt = yt * (-1)**(jcsb(j)+lb(j)+Jj+1) * sqrt(3.0d0*(2.0d0*Jj+1.0d0))
                yt = yt * SIXJ(2*le(k),2*Jj,int(2*jcsb(j)),2*1,2*lb(j),2*1)
                yt = yt * (-1)**(le(k)+0.5d0*(le(k)+lb(j)+1)) * sqrt(real(max(le(k),lb(j))))
                yt = yt * CONJG((kronij(in,k)+Smatd(in,k))*regints(j,k)+Im*(kronij(in,k)-Smatd(in,k))*iregints(j,k))
                y = y + yt
            end do
        end do
    end function
    !
    function RedElem_b(ni, Zint, no, Zo, jco, jcso, leo, Jo, jci, jcsi, lei, Ji, Wints) result(y)
        integer :: ni, no, ii, io
        real :: jco(no), jci(ni)
        integer :: leo(no), lei(ni), jcso(no), jcsi(ni), Jo, Ji
        real*8 :: Wints(ni,no)
        real*8 :: Zint(ni), Zo(no)
        real*8 :: y, yt
        real*8, EXTERNAL :: SIXJ
        y = 0.0d0
        do ii=1,ni
            do io=1,no
                !write(6,*) "Intermediate channel ", ii, "initial channel ", io
                yt = Zint(ii)*Zo(io)*Wints(ii,io)
                !write(6,*) "Z times radial integral", ii, io, yt
                yt = yt*kronij(int(2*jco(io)),int(2*jci(ii)))*kronij(jcso(io), jcsi(ii))
                !write(6,*) "Kroneckers", yt, "for jco ", jco(io), "and jci ", jci(ii),", jcso ",jcso(io), "and jcsi ", jcsi(ii)
                yt = yt*(-1)**(jcsi(ii)+leo(io)+Ji+1)*SQRT(real(2*Ji+1)*real(2*Jo+1))
                !write(6,*) "Phase", yt, (jci(ii)+leo(io)+Ji+1)
                yt = yt*SIXJ(2*lei(ii), 2*Ji, int(2*jcsi(ii)),2*Jo,2*leo(io),2*1)
                !write(6,*) "Six j", yt
                yt = yt * (-1)**(lei(ii)+0.5d0*(leo(io)+lei(ii)+1)) * sqrt(real(max(leo(io),lei(ii))))
                !write(6,*) "Orbital", yt
                y = y + yt 
            end do
        end do 
        !write(6,*) "Final result is ", y
    end function
    !
    function RedElem_bet(no, nc, ni, Zi, Smatdp, Zco, beta,in, iopen, iclosed, &
    jce, le, jcse, Je, jci, li, jcsi, Ji, regints, iregints) result(y)
        integer :: no, nc, ni, iclosed(nc), iopen(no)
        integer :: in, k, kk, j, Je, Ji
        integer :: le(no+nc), jcse(no+nc), li(ni), jcsi(ni)
        real*8 :: regints(no+nc,ni), iregints(no+nc,ni)
        real :: jce(no+nc), jci(ni)
        real*8 :: Zi(ni)
        complex*16 :: Smatdp(no,no), Zco(nc,no)
        real*8 :: beta(nc)
        complex*16 :: yt
        complex*16 :: y
        real*8, external :: SIXJ
        
        ! The radial element with the open part follows the same as with the continuum but with the right s matrix
        y = Zer
        do kk=1, no
            do j=1, ni
                k = iopen(kk)
                yt = Zi(j) * 0.5d0 * kronij(int(2*jce(k)), int(2*jci(j))) * kronij(int(jcsi(j)), int(jcse(k)))
                !write(6,*) "krons", k, j, yt
                yt = yt * (-1)**(jcsi(j)+li(j)+Je+1) * sqrt((2.0d0*Ji+1.0d0)*(2.0d0*Je+1.0d0))
                !write(6,*) "first term", yt
                yt = yt * SIXJ(2*le(k),2*Je,int(2*jcsi(j)),2*Ji,2*li(j),2*1)
                !write(6,*) "Six j", yt
                yt = yt * (-1)**(le(k)+0.5d0*(le(k)+li(j)+1)) * sqrt(real(max(le(k),li(j))))
                !write(6,*) "Orbital", le(k), li(j), le(k)+0.5d0*(le(k)+li(j)+1), yt
                yt = yt * CONJG((kronij(in,kk)+Smatdp(in,kk))*regints(k,j)+Im*(kronij(in,kk)-Smatdp(in,kk))*iregints(k,j))
                !write(6,*) "Radials", yt
                y = y + yt
            end do
        end do

        ! And now the closed channel part:
        
        do kk=1, nc
            do j=1, ni
                k = iclosed(kk)
                yt = Zi(j) * 0.5d0 * kronij(int(2*jce(k)), int(2*jci(j))) * kronij(int(jcsi(j)), int(jcse(k)))
                !write(6,*) "krons", k, j, yt
                yt = yt * (-1)**(jcsi(j)+li(j)+Je+1) * sqrt((2.0d0*Ji+1.0d0)*(2.0d0*Je+1.0d0))
                !write(6,*) "first term", yt
                yt = yt * SIXJ(2*le(k),2*Je,int(2*jcsi(j)),2*Ji,2*li(j),2*1)
                !write(6,*) "Six j", yt
                yt = yt * (-1)**(le(k)+0.5d0*(le(k)+li(j)+1)) * sqrt(real(max(le(k),li(j))))
                !write(6,*) "Orbital", le(k), li(j), le(k)+0.5d0*(le(k)+li(j)+1), yt
                yt = yt * CONJG(-1.d0*Zco(kk,in)*(COS(beta(kk))*regints(k,j)+SIN(beta(kk))*iregints(k,j)))
                !write(6,*) "Radials", yt
                y = y+yt
            end do
        end do
    end function
    !
    function RedElem_bet_closed(no, nc, ni, Zi, Smatdp, Zco, beta,in, iopen, iclosed, &
        jce, le, jcse, Je, jci, li, jcsi, Ji, regints, iregints) result(y)
            integer :: no, nc, ni, iclosed(nc), iopen(no)
            integer :: in, k, kk, j, Je, Ji
            integer :: le(no+nc), jcse(no+nc), li(ni), jcsi(ni)
            real*8 :: regints(no+nc,ni), iregints(no+nc,ni)
            real :: jce(no+nc), jci(ni)
            real*8 :: Zi(ni)
            complex*16 :: Smatdp(no,no), Zco(nc,no)
            real*8 :: beta(nc)
            complex*16 :: yt
            complex*16 :: y
            real*8 :: norm
            real*8, external :: SIXJ
            
            ! The radial element with the open part follows the same as with the continuum but with the right s matrix
            y = Zer
            ! Bound state normalization
            norm = 0.d0
            do kk=1, nc
                norm = norm + Abs(Zco(kk,in))**2 * (beta(kk)/PI-le(iclosed(kk)))**3
            end do
            ! And now the closed channel part:
            
            do kk=1, nc
                do j=1, ni
                    k = iclosed(kk)
                    yt = Zi(j) * 0.5d0 * kronij(int(2*jce(k)), int(2*jci(j))) * kronij(int(jcsi(j)), int(jcse(k)))
                    !write(6,*) "krons", k, j, yt
                    yt = yt * (-1)**(jcsi(j)+li(j)+Je+1) * sqrt((2.0d0*Ji+1.0d0)*(2.0d0*Je+1.0d0))
                    !write(6,*) "first term", yt
                    yt = yt * SIXJ(2*le(k),2*Je,int(2*jcsi(j)),2*Ji,2*li(j),2*1)
                    !write(6,*) "Six j", yt
                    yt = yt * (-1)**(le(k)+0.5d0*(le(k)+li(j)+1)) * sqrt(real(max(le(k),li(j))))
                    !write(6,*) "Orbital", le(k), li(j), le(k)+0.5d0*(le(k)+li(j)+1), yt
                    yt = yt * CONJG(-1.d0*Zco(kk,in)*(COS(beta(kk))*regints(k,j)+SIN(beta(kk))*iregints(k,j)))/sqrt(norm)
                    !write(6,*) "Radials", yt
                    y = y+yt
                end do
            end do
        end function
    !
    function RedElem_c(ne, ni, Zi, Smatd, in, jce, le, jcse, Je, jci, li, jcsi, Ji, regints, iregints) result(y)
        integer :: ne, ni 
        integer :: in, k, j, Je, Ji
        integer :: le(ne), jcse(ne), li(ni), jcsi(ni)
        real*8 :: regints(ne,ni), iregints(ne,ni)
        real :: jce(ne), jci(ni)
        real*8 :: Zi(ni)
        complex*16 :: Smatd(ne,ne)
        complex*16 :: yt
        complex*16 :: y
        real*8, external :: SIXJ
        y = Zer
        do k=1, ne
            do j=1, ni
                yt = Zi(j) * 0.5d0 * kronij(int(2*jce(k)), int(2*jci(j))) * kronij(int(jcsi(j)), int(jcse(k)))
                !write(6,*) "krons", k, j, yt
                yt = yt * (-1)**(jcsi(j)+li(j)+Je+1) * sqrt((2.0d0*Ji+1.0d0)*(2.0d0*Je+1.0d0))
                !write(6,*) "first term", yt
                yt = yt * SIXJ(2*le(k),2*Je,int(2*jcsi(j)),2*Ji,2*li(j),2*1)
                !write(6,*) "Six j", yt
                yt = yt * (-1)**(le(k)+0.5d0*(le(k)+li(j)+1)) * sqrt(real(max(le(k),li(j))))
                !write(6,*) "Orbital", le(k), li(j), le(k)+0.5d0*(le(k)+li(j)+1), yt
                yt = yt * CONJG((kronij(in,k)+Smatd(in,k))*regints(k,j)+Im*(kronij(in,k)-Smatd(in,k))*iregints(k,j))
                !write(6,*) "Radials", yt
                y = y + yt
            end do
        end do
    end function
    !
    function RedElem_c_SW(ne, ni, Zi, Uia, mua, in, jce, le, jcse, Je, jci, li, jcsi, Ji, regints, iregints) result(y)
        integer :: ne, ni 
        integer :: in, k, j, Je, Ji
        integer :: le(ne), jcse(ne), li(ni), jcsi(ni)
        real*8 :: regints(ne,ni), iregints(ne,ni)
        real :: jce(ne), jci(ni)
        real*8 :: Zi(ni)
        real*8 :: Uia(ne,ne), mua(ne)
        real*8 :: yt
        real*8 :: y
        real*8, external :: SIXJ
        y = Zer
        do k=1, ne
            do j=1, ni
                yt = Zi(j) * kronij(int(2*jce(k)), int(2*jci(j))) * kronij(int(jcsi(j)), int(jcse(k)))
                !write(6,*) "krons", k, j, yt
                yt = yt * (-1)**(jcsi(j)+li(j)+Je+1) * sqrt((2.0d0*Ji+1.0d0)*(2.0d0*Je+1.0d0))
                !write(6,*) "first term", yt
                yt = yt * SIXJ(2*le(k),2*Je,int(2*jcsi(j)),2*Ji,2*li(j),2*1)
                !write(6,*) "Six j", yt
                yt = yt * (-1)**(le(k)+0.5d0*(le(k)+li(j)+1)) * sqrt(real(max(le(k),li(j))))
                !write(6,*) "Orbital", le(k), li(j), le(k)+0.5d0*(le(k)+li(j)+1), yt
                yt = yt * Uia(k,in) * (regints(k,j)*COS(PI*mua(in))-SIN(PI*mua(in))*iregints(k,j))
                !write(6,*) "Radials", yt
                y = y + yt
            end do
        end do
    end function
    !
    function kronij(ip,j)
        implicit none
        integer, intent(in) :: ip,j
        integer :: kronij
        kronij = int(float((ip+j)-abs(ip-j))/float((ip+j)+abs(ip-j)))
        return
    end function kronij
    !
    function J0p_det(Ecml)
        implicit none 
        real*8 :: Ecml, J0p_Det, Zc(2)
        complex*16 :: J0_Det_c
        integer :: tst
        call J0_MQDT(Ecml, Zc, J0_Det_c,0,1,tst)
        if(Abs(AIMAG(J0_Det_c)).ge.1.0d-10) then
            print *, "Determinant complex, something is not right"
            print *, "Giving real part only"
        endif
        J0p_Det= real(J0_Det_c)
    end function J0p_det
    !
    function J1_Det(Ecml)
        implicit none 
        real*8 :: Ecml, J1_Det, Zc(5)
        complex*16 :: J1_Det_c
        integer :: tst
        call J1_MQDT(Ecml, Zc, J1_Det_c,0,1,tst)
        if(Abs(AIMAG(J1_Det_c)).ge.1.0d-10) then
            print *, "Determinant complex, something is not right"
            print *, "Giving real part only"
        endif
        J1_Det= real(J1_Det_c)
    end function J1_Det
    !
    function J2f_Det(Ecml)
        implicit none 
        real*8 :: Ecml, J2f_Det, Zc(3)
        complex*16 :: J2f_Det_c
        integer :: tst
        !print*, Ecm
        call J2f_MQDT(Ecml, Zc, J2f_Det_c,0,1,tst)
        if(Abs(AIMAG(J2f_Det_c)).ge.1.0d-10) then
            print *, "Determinant complex"
            print *, "Giving real part only"
        endif
        J2f_Det= real(J2f_Det_c)
        !print*, J2f_Det
    end function J2f_Det
    !
    function J2p_Det(Ecml)
        implicit none 
        real*8 :: Ecml, J2p_Det, Zc(3)
        complex*16 :: J2p_Det_c
        integer :: tst
        !print*, Ecm
        call J2p_MQDT(Ecml, Zc, J2p_Det_c,0,1,tst)
        if(Abs(AIMAG(J2p_Det_c)).ge.1.0d-10) then
            print *, "Determinant complex"
            print *, "Giving real part only"
        endif
        J2p_Det= real(J2p_Det_c)
        !print*, J2f_Det
    end function J2p_Det
    !
    function J3_Det(Ecml)
        implicit none 
        real*8 :: Ecml, J3_Det, Zc(6)
        complex*16 :: J3_Det_c
        integer :: tst
        !print*, Ecm
        call J3_MQDT(Ecml, Zc, J3_Det_c,0,1,tst)
        if(Abs(AIMAG(J3_Det_c)).ge.1.0d-10) then
            print *, "Determinant complex"
            print *, "Giving real part only"
        endif
        J3_Det= real(J3_Det_c)
        !print*, J2f_Det
    end function J3_Det
    !
    function J4_Det(Ecml)
        implicit none 
        real*8 :: Ecml, J4_Det, Zc(6)
        complex*16 :: J4_Det_c
        integer :: tst
        !print*, Ecm
        call J4_MQDT(Ecml, Zc, J4_Det_c,0,1,tst)
        if(Abs(AIMAG(J4_Det_c)).ge.1.0d-10) then
            print *, "Determinant complex"
            print *, "Giving real part only"
        endif
        J4_Det= real(J4_Det_c)
    end function J4_Det
    !
    subroutine InitialStates()
        implicit none
        real*8 :: Z(5)
        real*8 :: Z2f(3), Z2p(3), Z0p(2), Z3(6), Z4l(6)
        real*8 :: Ecm, hcm, ep, Ecm_t, det_t
        real*8 :: Ecmlo, Ecmhi
        real*8 :: J1res(100), J3res(100)
        real*8 :: J2f_res(700), J2p_res(700), J0p_res(700), J4l_res(700)
        integer :: nJ2fr, nJ2pr, nJ0pr, nJ4lr, nJ1, nJ3
        complex*16 :: det
        external :: root, find_r_rough
        integer :: ii, jj, true_state
        ! file to store the bound states
        open(25, file='mixingCoeff.txt',status='replace')
        open(15, file='interm_bounds_J2f.dat', status='replace')
        open(16, file='interm_bounds_J2p.dat', status='replace')
        open(18, file='interm_bounds_J0p.dat', status='replace')
        open(20, file='interm_bounds_J4l.dat', status='replace')
        open(17, file='init_bounds_J1.dat', status='replace')
        open(19, file='init_bounds_J3.dat', status='replace')
        !
        Ecm_t = 0.0d0
        det_t = 0.0d0
        Ecmlo = 119398.556d0
        Ecmhi = 126525.032d0
        print *, "----------------------------------------------------------------------------------------------------------------"
        print *, "                                                  J=1 states                                                 "
        print *, "----------------------------------------------------------------------------------------------------------------"
     !
        call find_r_rough(J1_Det, Ecmlo-1/1200.0d-7+1000d0, Ecmhi-1/1200.d-7-1.0d3, 0.5d0, nJ1, J1res )
        print *, "Found ", nJ1, " possible states. Looking between ", Ecmlo-1/1200.0d-7+1.0d3, "and ", Ecmhi-1/1200.0d-7-1.0d3
        print *, J1res(1:4)
        print *, "One state"
        hcm = 0.1d0
        ep = 1.d-06
        det = 1.0d0
        Ecm = J1res(1)
        call root(Ecm,hcm,1,ep,Ecm_t,det_t,J1_Det)
        call J1_MQDT(Ecm_t, Z, det,1,0, true_state)
        if(true_state.eq.1) write(17,1787) Ecm_t, (Z(jj), jj=1, 5)
    
        print*, "Two states"
        hcm = 0.1d0
        ep = 1.d-06
        det = 1.0d0
        Ecm = J1res(2)
        call root(Ecm,hcm,1,ep,Ecm_t,det_t,J1_Det)
        call J1_MQDT(Ecm_t, Z, det,1,0, true_state)
        if(true_state.eq.1) write(17,1787) Ecm_t, (Z(jj), jj=1, 5)
    
        print *, "Three states"
        hcm = 0.1d0
        ep = 1.d-06
        det = 1.0d0
        Ecm = J1res(3)
        call root(Ecm,hcm,1,ep,Ecm_t,det_t,J1_Det)
        call J1_MQDT(Ecm_t, Z, det,1,0,true_state)
        if(true_state.eq.1) write(17,1787) Ecm_t, (Z(jj), jj=1, 5)
    
        print *, "Four states"
        hcm = 0.1d0
        ep = 1.d-06
        det = 1.0d0
        Ecm = J1res(4)
        call root(Ecm,hcm,1,ep,Ecm_t,det_t,J1_Det)
        call J1_MQDT(Ecm_t, Z, det,1,0, true_state)
        if(true_state.eq.1) write(17,1787) Ecm_t, (Z(jj), jj=1, 5)
     !
        print *, "----------------------------------------------------------------------------------------------------------------"
        print *, "                                                  J=3 states                                                 "
        print *, "----------------------------------------------------------------------------------------------------------------"
     !
        call find_r_rough(J3_Det, Ecmlo-1/1200.0d-7, Ecmhi-1/1200.d-7, 0.5d0, nJ3, J3res )
        print *, "Found ", nJ3, " possible states. Looking between ", Ecmlo-1/1200.0d-7+1.0d3, "and ", Ecmhi-1/1200.0d-7-1.0d3
        do ii=1,nJ3
            hcm = 0.1d0
            ep = 1.d-06
            det = 1.0d0
            Ecm = J3res(ii)
            call root(Ecm,hcm,1,ep,Ecm_t,det_t,J3_Det)
            call J3_MQDT(Ecm_t, Z3, det,1,0,true_state)
            if(true_state.eq.1) write(19,1789) Ecm_t, (Z3(jj), jj=1, 6)
        end do
        
        Ecm_t = 0.0d0
        det_t = 0.0d0
        Ecmlo = 105000.556d0
        Ecmhi = 126525.032d0
     !
        print *, "-----------------------------------------------------------------------------------------------------------------"
        print *, "                                                 J=2 f states                                                 "
        print *, "-----------------------------------------------------------------------------------------------------------------"
        call find_r_rough(J2f_Det, Ecmlo, Ecmhi, 0.5d0, nJ2fr, J2f_res)
        do ii =1, nJ2fr
            Z2f= 0.0d0
            hcm = 0.01d0
            ep = 1.d-06
            det = 1.0d0
            Ecm = J2f_res(ii)
            call root(Ecm,hcm,1,ep,Ecm_t,det_t,J2f_Det)
            call J2f_MQDT(Ecm_t, Z2f, det,1,0,true_state)
            if(true_state.eq.1) write(15,1687) Ecm_t, (Z2f(jj), jj=1, 3)
        end do
        print *, "-----------------------------------------------------------------------------------------------------------------"
        print *, "                                                 J=2 p states                                                    "
        print *, "-----------------------------------------------------------------------------------------------------------------"
        call find_r_rough(J2p_Det, Ecmlo, Ecmhi, 0.5d0, nJ2pr, J2p_res)
        do ii =1, nJ2pr
            Z2p= 0.0d0
            hcm = 0.01d0
            ep = 1.d-06
            det = 1.0d0
            Ecm = J2p_res(ii)
            call root(Ecm,hcm,1,ep,Ecm_t,det_t,J2p_Det)
            call J2p_MQDT(Ecm_t, Z2p, det,1,0,true_state)
            if(true_state.eq.1) write(16,1687) Ecm_t, (Z2p(jj), jj=1,3)
        end do
        print *, "-----------------------------------------------------------------------------------------------------------------"
        print *, "                                                 J=0 p states                                                 "
        print *, "-----------------------------------------------------------------------------------------------------------------"
        call find_r_rough(J0p_Det, Ecmlo, Ecmhi, 0.5d0, nJ0pr, J0p_res)
        do ii =1, nJ0pr
            Z0p= 0.0d0
            hcm = 1.0d-03
            ep = 1.d-06
            det = 1.0d0
            Ecm = J0p_res(ii)
            call root(Ecm,hcm,2,ep,Ecm_t,det_t,J0p_Det)
            call J0_MQDT(Ecm_t, Z0p, det,1,0, true_state)
            if(true_state.eq.1) write(18,1887) Ecm_t, (Z0p(jj), jj=1, 2)
        end do

        print *, "-----------------------------------------------------------------------------------------------------------------"
        print *, "                                                 J=4 l states                                                 "
        print *, "-----------------------------------------------------------------------------------------------------------------"
        call find_r_rough(J4_Det, Ecmlo, Ecmhi, 0.5d0, nJ4lr, J4l_res)
        print *, "Found ", nJ4lr, " states"
        do ii =1, nJ4lr
            Z4l= 0.0d0
            hcm = 1.0d-03
            ep = 1.d-06
            det = 1.0d0
            Ecm = J4l_res(ii)
            call root(Ecm,hcm,2,ep,Ecm_t,det_t,J4_Det)
            call J4_MQDT(Ecm_t, Z4l, det,1,0, true_state)
            if(true_state.eq.1) write(20,1789) Ecm_t, (Z4l(jj), jj=1, 6)
        end do

        print*, nJ2pr+nJ2fr+nJ0pr+nJ4lr
        close(15)
        close(16)
        close(17)
        close(18)
        close(19)
        close(20)
     1787 format(6(e25.10e4, 3x))
     1789 format(7(e25.10e4,3x))
     1687 format(4(e25.10e4, 3x))
     1887 format(3(e25.10e4, 3x))

    end subroutine InitialStates
    !
    subroutine BoundBoundMatrices(nso, no, Zos, nsi, ni, Zis, jin, jci, jcsi, lei, maxWints)
        integer nso, no, jo, nsi, ni, jin
        integer :: lei(ni), jcsi(ni)
        real :: jci(ni)
        real*8 :: Zos(nso, no), Zis(nsi, ni)
        real*8 :: maxWints(nso, nsi, ni, no)
        real*8 :: Res(nso, nsi)
        jo = 1
        open(10, file = "BoundRedElemsJin"//char(jin+48)//"Jo"//char(jo+48)//"outerl"//char(lei(1)+48)//".dat")
        write(6,*) "Working on interms ", Jin, "with orbital ", lei(:)
        do i1=1, nso
            do i2=1, nsi
                write(6,*) "Currently at initial bound ", i1, "and interm ", i2
                Res(i1, i2) = RedElem_b(ni, Zis(i2,:),no, Zos(i1, :), jcb, jcsb, lb, jo, jci, jcsi, lei, jin, maxWints(i1,i2,:,:))
            end do
            write(10,1985) Res(i1,:)
        end do
        1985 format(24(e18.10e2, 3x))
        close(10)
    end subroutine BoundBoundMatrices
    !
    subroutine write_matrix(nr, nc, a)
        integer :: nr, nc
        integer :: j
        real*8, dimension(nr,nc) :: a
        write(*,*)
        do i = 1, nr
           write(*,*) (a(i,j), j = 1, nc)
        end do
    end subroutine write_matrix
    !
    subroutine zwrite_matrix(nr, nc, a)
        integer :: nr, nc
        integer :: j
        complex*16, dimension(nr,nc) :: a
        write(*,*)
        do i = 1, nr
           write(*,*) (a(i,j), j = 1, nc)
        end do
    end subroutine zwrite_matrix
    !
    function i2a(ii) result(res)
        character(:),allocatable :: res
        integer,intent(in) :: ii
        character(range(ii)+2) :: tmp
        write(tmp,'(i0)') ii
        res = trim(tmp)
    end function
    !
    function Stark_Shift_b(np, ncp, Zp, Ep, jcp, jcsp, lep, Jp, Wp, &
                           nm, ncm, Zm, Em, jcm, jcsm, lem, Jm, Wm, &
                               nco, Zo, Eo, jco, jcso, leo, Jo, Fo, ww, gam) result(ss)
        integer :: np, ncp, Jp, nm, ncm, Jm, nco, Jo, iloc
        real*8 :: Zp(np,ncp), Ep(np), Zm(nm,ncm), Em(nm), Zo(nco), Eo
        real :: jcp(ncp), jcm(ncm), jco(nco)
        integer :: jcsp(ncp),  lep(ncp), jcsm(ncm),  lem(ncm), jcso(nco),  leo(nco)
        real*8 :: Wp(np,ncp,nco), Wm(nm,ncm,nco)
        complex*16 :: dipole
        real*8 :: ss(2), Fo, ww, gam
        
        ss = 0.d0
        do iloc=1, np
            dipole = Zone*RedElem_b(ncp, Zp(iloc,:), nco, Zo(:), jco, jcso, leo, Jo, jcp, jcsp, lep, Jp, Wp(iloc,:,:))
            dipole = dipole * (-1)**(Jp-0) * qsum_c(Jp,0,Jo,0,(/0.d0,0.d0,1.d0/))
            ss(1) = ss(2) + (-0.5d0 * ABS(dipole)**2 * 1/(Ep(iloc)-Eo-ww))
        end do

        do iloc=1, nm
            dipole = Zone*RedElem_b(ncm, Zm(iloc,:), nco, Zo(:), jco, jcso, leo, Jo, jcm, jcsm, lem, Jm, Wm(iloc,:,:))
            dipole = dipole * (-1)**(Jm-0) * qsum_c(Jm,0,Jo,0,(/0.d0,0.d0,1.d0/))
            ss(2) = ss(2) + (-0.5d0 * ABS(dipole)**2 * 1/(Em(iloc)-Eo-ww))
        end do
        ss = ss * Fo**2 * 0.25d0 * SQRT(PI/2) * gam * (EXP(-0.5d0*gam**2 * ww**2)+1.d0)
        return
    end function

    function Stark_Shift_c(np, ncp, Ep, jcp, jcsp, lep, Jp, rp, irp,&
                           nm, ncm, Em, jcm, jcsm, lem, Jm, rm, irm,&
                           dE, nco, Zo, Eo, jco, jcso, leo, Jo, Fo, ww, gam) result(ss)
        integer :: np,ncp,Jp,nm,ncm,Jm,nco,Jo, iloc, ic
        real*8 :: Zo(nco)
        real*8 :: Ep(np), Em(nm), Eo, dE
        real*8 :: rp(np,ncp,nco), irp(np,ncp,nco), rm(nm,ncm,nco), irm(nm,ncm,nco)
        real :: jcp(ncp), jcm(ncm), jco(nco)
        integer :: jcsp(ncp),  lep(ncp), jcsm(ncm),  lem(ncm), jcso(nco),  leo(nco)
        complex*16 :: Sdp(ncp,ncp), Sdm(ncm,ncm)
        real*8 :: Kp(ncp,ncp), Km(ncm,ncm)
        complex*16 ::dipole
        real*8 :: ss(2), Fo, ww, gam
        
        ss = 0.d0

        do iloc=1,np
            if(Jo.eq.0) then
                call J1_MQDT_at(Ep(iloc), Sdp, Kp, 0)
            else if(Jo.eq.2) then
                call J3_MQDT_at(Ep(iloc), Sdp, Kp, 0)
            else if(Jo.eq.4) then
                call J5_MQDT_at(Ep(iloc), Sdp, Kp, 0)
            else 
                stop "Not programed to handle that state coupled with the continuum"
            endif
            do ic=1,ncp
            dipole = RedElem_c(ncp, nco, Zo, Sdp, ic, jcp, lep, jcsp, Jp, jco, leo, jcso,&
             Jo, rp(iloc,:,:), irp(iloc,:,:))
            dipole = dipole * (-1)**(Jp-0) * qsum_c(Jp,0,Jo,0,(/0.d0,0.d0,1.d0/))
            ss(1) = ss(1) + dE*ABS(dipole)**2/(Ep(iloc)-Eo-ww)
            end do
        end do

        do iloc=1,nm
            if(Jo.eq.0) then
                go to 11
            else if(Jo.eq.2) then
                call J1_MQDT_at(Em(iloc), Sdm, Km, 0)
            else if(Jo.eq.4) then
                call J3_MQDT_at(Em(iloc), Sdm, Km, 0)
            else 
                stop "Not programed to handle that state coupled with the continuum"
            endif
            do ic=1,ncm
            dipole = RedElem_c(ncm, nco, Zo, Sdm, ic, jcm, lem, jcsm, Jm, jco, leo, jcso,&
             Jo, rm(iloc,:,:), irm(iloc,:,:))
            dipole = dipole * (-1)**(Jm-0) * qsum_c(Jm,0,Jo,0,(/0.d0,0.d0,1.d0/))
            ss(2) = ss(2) + dE*ABS(dipole)**2/(Em(iloc)-Eo-ww)
            end do
    11      continue
        end do
        ss = ss * Fo**2 * 0.25d0 * SQRT(PI/2) * gam * (EXP(-0.5d0*gam**2 * ww**2)+1.d0)
        return
    end function
end module

