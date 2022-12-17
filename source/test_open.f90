program main
    use matmod
    use ArJ1J2
    implicit none
    real*8, parameter :: IP32 = 127109.88d0, IP12= IP32+1432.0d0
    integer :: gsamps, ii, jj
    real*8 :: ognu
    real*8, external :: rint

    call setup()
    !
    gsamps = 100
    ognu = fnu(IP32,IP12)+0.05
    !
    open(66,file='Z_profiles_J2f.dat')
    call profiles(IP32+10.,IP12-RydAr/15**2, 5.d-3, 66, 23)
    close(66)
    open(66,file='Z_profiles_J2p.dat')
    call profiles(IP32+10.,IP12-RydAr/15**2, 5.d-3, 66, 21)
    close(66)
    open(66,file='Z_profiles_J0p.dat')
    call profiles(IP32+10.,IP12-RydAr/15**2, 5.d-3, 66, 1)
    close(66)
    !
    contains
    !
    subroutine find_resonances(n, nuo, unit)
        integer :: n, fails, samples, unit
        real*8 :: nuo, nump, num, dEcm, Ecmb, Zb, Ecmf, Zf
        complex*16, allocatable :: Sdp(:,:), Zco(:,:)
        integer, allocatable :: iopen(:), iclosed(:)
        real*8, allocatable :: beta(:)
        do ii=0,n
            fails = 0
            num = nuo+ii+0.5
            nump = 0.d0
        11  dEcm = (1000/num**3)/500
            !write(6,*) nump, num
            samples = 250
            Ecmb = IP12-RydAR/num**2-samples*dEcm
            call J2f_MQDT_bet(Ecmb, Sdp,Zco,iopen,iclosed,beta)
            Zb = ABS(Zco(1,1))**2
            deallocate(sdp,Zco,iopen, iclosed,beta)
            Ecmf = IP12-RydAR/num**2+samples*dEcm
            call J2f_MQDT_bet(Ecmf, Sdp,Zco,iopen,iclosed,beta)
            Zf = ABS(Zco(1,1))**2
            deallocate(sdp,Zco,iopen, iclosed,beta)
            if((nump.ne.num).and.(fails.lt.2000)) then
                nump = num
                call selector(num,dEcm,samples,Ecmf)
                num = fnu(Ecmf,IP12)
                fails = fails+1
                if(fails.eq.2000) stop "Excedeed number of iterations"
                go to 11
            else
                call writer(num,dEcm,samples,unit)
                go to 12
            endif
        12 continue
        end do
    end subroutine 
    !
    subroutine selector(numc, dE, sams, Ef)
        real*8 :: numc, dE, Ef, Ec, Zb
        integer :: sams
        complex*16, allocatable :: Sdp(:,:), Zco(:,:)
        integer, allocatable :: iopen(:), iclosed(:)
        real*8, allocatable :: beta(:)
        zb = 0.0d0
        do jj=-sams,sams
            Ec = IP12-RydAR/numc**2+jj*dE
            call J2f_MQDT_bet(Ec, Sdp,Zco,iopen,iclosed,beta)
            if(zb.lt.ABS(Zco(1,1))**2) then
                Ef = Ec
                Zb = ABS(Zco(1,1))**2
            endif
            deallocate(sdp,Zco,iopen, iclosed,beta)
        end do
    end subroutine 
    !
    subroutine writer(numc, dE, sams, unit)
        real*8 :: numc, dE, Ec, zb
        integer :: sams, unit
        real*8 :: es(2*sams),zs(2*sams)
        complex*16 :: zoo
        complex*16, allocatable :: Sdp(:,:), Zco(:,:)
        integer, allocatable :: iopen(:), iclosed(:)
        real*8, allocatable :: beta(:)
        zb = 0.0d0
        do jj=1,2*sams
            Ec = IP12-RydAR/numc**2-sams*dE+jj*dE
            call J2f_MQDT_bet(Ec, Sdp,Zco,iopen,iclosed,beta)
            if(jj.eq.sams) zoo=Zco(1,1)
            es(jj) = Ec/CMperAU
            zs(jj) = ABS(Zco(1,1))**2
            deallocate(sdp,Zco,iopen,iclosed,beta)
        end do
        write(unit,*) es(sams)*CMperAU,CMperAU*rint(zs,1,2*sams,10,dE/CmperAU)*2/(PI*zs(sams))
    end subroutine 
    !
    subroutine profiles(Elo, Ehi, dEcm, unit, J)
        integer :: samples, unit,J
        real*8 :: Elo, Ehi, dEcm, Ecmb
        complex*16, allocatable :: Sdp(:,:), Zco(:,:)
        integer, allocatable :: iopen(:), iclosed(:)
        real*8, allocatable :: beta(:)
        samples = int((Ehi-Elo)/dEcm)
        do ii=0,samples
            Ecmb = Elo + ii * dEcm
            if(J.eq.23) then 
                call J2f_MQDT_bet(Ecmb, Sdp,Zco,iopen,iclosed,beta)
            else if(J.eq.21) then
                call J2p_MQDT_bet(Ecmb, Sdp,Zco,iopen,iclosed,beta)
            else 
                call J0_MQDT_bet(Ecmb, Sdp,Zco,iopen,iclosed,beta)
            endif 
            write(unit,*) Ecmb, ABS(Zco(:,:))**2
        end do
    end subroutine
end program 
