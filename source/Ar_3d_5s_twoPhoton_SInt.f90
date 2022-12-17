program main
    use j_recoup
    use ArJ1J2
    use ERRFUN
    implicit none
    integer, parameter :: nbs=24, nbs0=13, nconts=200, nos = 4
    integer :: nbst
    integer :: nbs0t, Na
    integer :: jj, mm, kk, ii1
    integer :: nl, ncon
    real*8 :: Ar, Ai
    complex*16 :: A(4)
    character(len=:), allocatable :: file_root
    real*8 :: Eo, w, gam, evec(3)
    real*8 :: Zo(4, 5), Ebs(4)
    real*8 :: ZJ0(nbs0,2), EnJ0(nbs0)
    real*8 :: ZJ2f(nbs,3), EnJ2f(nbs)
    real*8 :: ZJ2p(nbs,3), EnJ2p(nbs)
    real*8 :: rRI_CI_J0p_J1(nconts,nbs0,5,2), rRI_CI_J0p_J3(nconts,nbs0,5,2)
    real*8 :: iRI_CI_J0p_J1(nconts,nbs0,5,2), iRI_CI_J0p_J3(nconts,nbs0,6,2)
    real*8 :: rRI_CI_J2f_J1(nconts,nbs,5,3), rRI_CI_J2f_J3(nconts,nbs,6,3)
    real*8 :: iRI_CI_J2f_J1(nconts,nbs,5,3), iRI_CI_J2f_J3(nconts,nbs,6,3)
    real*8 :: rRI_CI_J2p_J1(nconts,nbs,5,3), rRI_CI_J2p_J3(nconts,nbs,6,3)
    real*8 :: iRI_CI_J2p_J1(nconts,nbs,5,3), iRI_CI_J2p_J3(nconts,nbs,6,3)
    real*8 :: RI_OI_J1_J2f(nos, nbs, 3, 5), RI_OI_J1_J2p(nos, nbs, 3,5)
    real*8 :: RI_OI_J1_J0p(nos, nbs0,2, 5)
    real*8 ::  Elo, dEcont
    nbst = 24
    nbs0t = 13
    Na=-5
    !-------------------------------- calling set up from gensub -----------------------
    call SETUP
    !-----------------------------------------------------------------------------------
    !call InitialStates
    !stop
    ! Call if the we don't hacve data on the initial bound state and the initial 
    !call InitialStates
    !Retrieving from precomputed files.
    !Loading the Zcoefficients
        open(10, file='init_bounds_J1.dat')
        do mm=1,4
            read(10,*) Ebs(mm), (Zo(mm,kk), kk=1,5)
            !write(6,*) Ebs(mm), (Zo(mm,kk), kk=1,5)
        end do
        close(10)

        open(10, file='interm_bounds_J2f.dat')
        do mm=1,24
            read(10,*) EnJ2f(mm), (ZJ2f(mm,kk), kk=1,3)
            !write(6,*) EnJ2f(mm), (ZJ2f(mm,kk), kk=1,3)
        end do
        close(10)

        open(10, file='interm_bounds_J2p.dat')
        do mm=1,24
            read(10,*) EnJ2p(mm), (ZJ2p(mm,kk), kk=1,3)
            !write(6,*) EnJ2p(mm), (ZJ2p(mm,kk), kk=1,3)
        end do
        close(10)

        open(10, file='interm_bounds_J0p.dat')
        do mm=1,13
            read(10,*) EnJ0(mm), (ZJ0(mm,kk), kk=1,2)
            !write(6,*) EnJ0(mm), (ZJ0(mm,kk), kk=1,2)
        end do
        close(10)
    !
    !Loading the radial elements
        write(6,*) "Loading Radial Elements into four dimensional arrays"
        ! Radial elements from initial J=1 state, to intermediate J=2 state
        ! f states
        file_root = 'Python_Code_RadInts/Rad Ints bb for J=2 lint f/interm_rad_ints_'
        call Init_big_matsW(nos, nbs, RI_OI_J1_J2f, 3, 5, len(file_root), file_root)
        ! p states
        !J=2
        file_root = 'Python_Code_RadInts/Rad Ints bb for J=2 lint p/interm_rad_ints_'
        call Init_big_matsW(nos, nbs, RI_OI_J1_J2p, 3, 5, len(file_root), file_root)
        !J=0
        file_root = 'Python_Code_RadInts/Rad Ints bb for J=0 lint p/interm_rad_ints_'
        call Init_big_matsW(nos, nbs0, RI_OI_J1_J0p, 2, 5, len(file_root), file_root)
        ! Radial elements to J=1 continuum, from an intermediate J=2 state
        ! f states
        file_root = 'Python_Code_RadInts/RI b(J=2,l=f) c(J=1)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs, nconts, rRI_CI_J2f_J1, 5, 3, len(file_root), file_root)
        file_root = 'Python_Code_RadInts/RI b(J=2,l=f) c(J=1)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs, nconts, iRI_CI_J2f_J1, 5, 3, len(file_root), file_root)
        ! p states
        !J=2
        file_root = 'Python_Code_RadInts/RI b(J=2,l=p) c(J=1)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs, nconts, rRI_CI_J2p_J1, 5, 3, len(file_root), file_root)
        file_root = 'Python_Code_RadInts/RI b(J=2,l=p) c(J=1)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs, nconts, iRI_CI_J2p_J1, 5, 3, len(file_root), file_root)
        !J=0
        file_root = 'Python_Code_RadInts/RI b(J=0,l=p) c(J=1)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs0, nconts, rRI_CI_J0p_J1, 5, 2, len(file_root), file_root)
        file_root = 'Python_Code_RadInts/RI b(J=0,l=p) c(J=1)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs0, nconts, iRI_CI_J0p_J1, 5, 2, len(file_root), file_root)

        ! Radial elements to J=3 continuum, from an intermediate J=2 state
        ! f states
        file_root = 'Python_Code_RadInts/RI b(J=2,l=f) c(J=3)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs, nconts, rRI_CI_J2f_J3, 6, 3, len(file_root), file_root)
        file_root = 'Python_Code_RadInts/RI b(J=2,l=f) c(J=3)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs, nconts, iRI_CI_J2f_J3, 6, 3, len(file_root), file_root)
        ! p states
        !J=2
        file_root = 'Python_Code_RadInts/RI b(J=2,l=p) c(J=3)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs, nconts, rRI_CI_J2p_J3, 6, 3, len(file_root), file_root)
        file_root = 'Python_Code_RadInts/RI b(J=2,l=p) c(J=3)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs, nconts, iRI_CI_J2p_J3, 6, 3, len(file_root), file_root)
        !J=0
        rRI_CI_J0p_J3=0.0d0
        iRI_CI_J0p_J3=0.0d0

    !

    !Laser parameters
    open(98, file='SingleZs_laser_params_two_photon.dat', status='old')
    read(98,*) Eo, w, gam
    read(98,*) evec(:)
    read(98,*) nl
    close(98)
    print *, gam, Eo
    !Initial state amplitudes
    open(98, file='init_amplitudes.dat', status='old')
    !write(6,*) "don't mind me pritin"
    do ii1=1,4
        read(98,*) Ar, Ai
        write(6,*) Ar, Ai
        A(ii1) = Ar*Zone+Im*Ai
    !    write(6,*) A(ii1)
    enddo
    close(98)
    !Continuum state energies.
    open(08,file='Python_Code_RadInts/EcontEV_det.dat')
    read(08,*) ncon, Elo, dEcont
    close(08)

    !Rabi
    call CoreAllIs()
    stop
    !
    do jj=1,nbs
        call CoreZs((/NA,NA/),(/Na,Na/),(/jj,jj/))
        call CoreZs((/NA,NA/),(/jj,jj/), (/NA,NA/))
    end do
    
    do jj=1,nbs0
        call CoreZs((/jj,jj/),(/NA,NA/),(/NA,NA/))
    end do
    
    contains
    subroutine Core(SelJ0,SelJ2p, SelJ2f)
        integer :: ie, ntd
        integer :: Jt, Mt, nt
        integer :: iit, ii2
        integer :: SelJ0(2), SelJ2p(2), SelJ2f(2)
        real*8 :: ProbsJ12(200), ProbsJ32(200)
        real*8 ::  Econt
        real*8 :: tdlo, tdf, dtd
        real*8 :: to
        !real*8 :: RI_CI_J0_J1(nos, nbs, 2, 5)
        complex*16 :: BJ1(5), BJ3(6), BtmpJ1(5, 4), BtmpJ3(6,4), Phases(4)
        complex*16 :: BwJ1(nconts,5,4), BwJ3(nconts, 6,4)


        open(78, file="IndividualInterms/'"//itoa(nl)//"SI"//filepreface(SelJ0, SelJ2p, SelJ2f)//&
        "RFinalContJ1States.dat", status='replace')
        open(178, file="IndividualInterms/'"//itoa(nl)//"SI"//filepreface(SelJ0, SelJ2p, SelJ2f)//&
        "IFinalContJ1States.dat", status='replace')
        open(79, file="IndividualInterms/'"//itoa(nl)//"SI"//filepreface(SelJ0, SelJ2p, SelJ2f)//&
        "RFinalContJ3States.dat", status='replace')
        open(179, file="IndividualInterms/'"//itoa(nl)//"SI"//filepreface(SelJ0, SelJ2p, SelJ2f)//&
        "IFinalContJ3States.dat", status='replace')
        Jt = 1
        Mt = 0
        nt = 5
        do ie=1, ncon
            BtmpJ1 = Zer
            do ii1=1,nt
                call cont_Amp_TP_vec_SI(nt, Elo+ie*dEcont, ii1, Jt, Mt, BtmpJ1(ii1,:),&
                Zo, Ebs,& !Initial bound state info
                nbs, SelJ2f, EnJ2f, ZJ2f, &
                nbs, SelJ2p, EnJ2p, ZJ2p, & !Intermediate state information
                nbs0, SelJ0, EnJ0, ZJ0, &
                RI_OI_J1_J2f, RI_OI_J1_J2f, RI_OI_J1_J0p,&
                rRI_CI_J2f_J1(ie,:,:,:), iRI_CI_J2f_J1(ie,:,:,:),&
                rRI_CI_J2p_J1(ie,:,:,:), iRI_CI_J2p_J1(ie,:,:,:), &
                rRI_CI_J0p_J1(ie,:,:,:), iRI_CI_J0p_J1(ie,:,:,:), &  !Radial integral tensors
                evec, gam, Eo, w)
            end do
            write(78,1978) ((REAL(BtmpJ1(ii1,ii2)), ii2=1,4),ii1=1,nt)
            write(178,1978) ((AIMAG(BtmpJ1(ii1,ii2)), ii2=1,4),ii1=1,nt)
            BwJ1(ie,:,:)=BtmpJ1(:,:)
        end do
        close(78)
        !Computations for J=3 
        Jt = 3
        Mt = 0
        nt = 6
        do ie=1,nconts
            BJ3 = Zer
            do ii1=1,nt
                call cont_Amp_TP_vec_SI(nt, Elo+ie*dEcont, ii1, Jt, Mt, BtmpJ3(ii1,:),&
                Zo, Ebs,& !Initial bound state info
                nbs, SelJ2f,EnJ2f, ZJ2f,&
                nbs, SelJ2p, EnJ2p, ZJ2p, &
                nbs0, SelJ0, EnJ0, ZJ0,& !Intermediate state information
                RI_OI_J1_J2f, RI_OI_J1_J2f, RI_OI_J1_J0p, &
                rRI_CI_J2f_J3(ie,:,:,:), iRI_CI_J2f_J3(ie,:,:,:),&
                rRI_CI_J2p_J3(ie,:,:,:), iRI_CI_J2p_J3(ie,:,:,:), & 
                rRI_CI_J0p_J3(ie,:,:,:), iRI_CI_J0p_J3(ie,:,:,:), &!Radial integral tensors
                evec, gam, Eo, w)
            end do
            write(79,1979) ((REAL(BtmpJ3(ii1,ii2)), ii2=1,4),ii1=1,nt)
            write(179,1979) ((AIMAG(BtmpJ3(ii1,ii2)), ii2=1,4),ii1=1,nt)
            BwJ3(ie, :, :) = BtmpJ3(:,:)
        end do
        close(79)

        ntd = 1800
        tdlo = 0.00d0
        tdf = 7.5d4 
        dtd = (tdf-tdlo)/ntd
        open(10, file="IndividualInterms/'"//itoa(nl)//"SI"//filepreface(SelJ0, SelJ2p, SelJ2f)//"tdel_probsJ12.dat")
        open(11, file="IndividualInterms/'"//itoa(nl)//"SI"//filepreface(SelJ0, SelJ2p, SelJ2f)//"tdel_probsJ32.dat")
        
        do iit=0,ntd
            to = tdlo+iit*dtd
            ProbsJ12 = 0.0d0
            ProbsJ32 = 0.0d0
            do ie=1, ncon
                Econt = Elo+ie*dEcont
                Phases = EXP(Im*to/CMperAU*(/Econt-Ebs/))
                BtmpJ1 = BwJ1(ie,:,:)
                BtmpJ3 = BwJ3(ie,:,:)
                !Go over the channels of J1
                Bj1 = Zer
                do ii1=1,5
                    do ii2=1,4
                        BJ1(ii1) = BJ1(ii1) + A(ii2)*BtmpJ1(ii1,ii2)*Phases(ii2)
                    end do
                    if((ii1.eq.1).or.(ii1.eq.4)) then
                        ProbsJ12(ie) = ProbsJ12(ie)+ABS(BJ1(ii1))**2
                    else
                        ProbsJ32(ie) = ProbsJ32(ie)+ABS(BJ1(ii1))**2
                    endif
                end do

                !Go over the channels of J3
                Bj3 = Zer
                do ii1=1,6
                    do ii2=1,4
                        BJ3(ii1) = BJ3(ii1) + A(ii2)*BtmpJ3(ii1,ii2)*Phases(ii2)
                    end do
                    if((ii1.eq.3).or.(ii1.eq.6)) then
                        ProbsJ12(ie) = ProbsJ12(ie)+ABS(BJ3(ii1))**2
                    else
                        ProbsJ32(ie) = ProbsJ32(ie)+ABS(BJ3(ii1))**2
                    endif
                end do 
            end do
            write(10,1910) to, ProbsJ12(:)
            write(11,1910) to, ProbsJ32(:)
        end do
        close(10)
        close(11)
        1978 format(5(4(E20.9E3,1x),3x))
        1979 format(6(4(E20.9E3,1x),3x))
        1910 format(201(E20.9E3,1x))
    end subroutine
    !
    subroutine CoreZs(SelJ0, SelJ2p, SelJ2f)
        integer :: ie
        integer :: Jt, Mt, nt
        integer ::  ii2
        integer :: SelJ0(2), SelJ2p(2), SelJ2f(2)
        !real*8 :: RI_CI_J0_J1(nos, nbs, 2, 5)
        complex*16 :: BtmpJ1(5, 4), BtmpJ3(6,4)
        complex*16 :: BwJ1(nconts,5,4), BwJ3(nconts, 6,4)


        open(78, file="ZetaMatrices/"//itoa(nl)//"SI"//fileZpreface(SelJ0(1), SelJ2p(1), SelJ2f(1))//&
        "RFinalContJ1States.dat", status='replace')
        open(178, file="ZetaMatrices/"//itoa(nl)//"SI"//fileZpreface(SelJ0(1), SelJ2p(1), SelJ2f(1))//&
        "IFinalContJ1States.dat", status='replace')
        open(79, file="ZetaMatrices/"//itoa(nl)//"SI"//fileZpreface(SelJ0(1), SelJ2p(1), SelJ2f(1))//&
        "RFinalContJ3States.dat", status='replace')
        open(179, file="ZetaMatrices/"//itoa(nl)//"SI"//fileZpreface(SelJ0(1), SelJ2p(1), SelJ2f(1))//&
        "IFinalContJ3States.dat", status='replace')
        Jt = 1
        Mt = 0
        nt = 5
        do ie=1, ncon
            BtmpJ1 = Zer
            do ii1=1,nt
                call cont_Amp_TP_vec_SI(nt, Elo+ie*dEcont, ii1, Jt, Mt, BtmpJ1(ii1,:),&
                Zo, Ebs,& !Initial bound state info
                nbs, SelJ2f, EnJ2f, ZJ2f, &
                nbs, SelJ2p, EnJ2p, ZJ2p, & !Intermediate state information
                nbs0, SelJ0, EnJ0, ZJ0, &
                RI_OI_J1_J2f, RI_OI_J1_J2f, RI_OI_J1_J0p,&
                rRI_CI_J2f_J1(ie,:,:,:), iRI_CI_J2f_J1(ie,:,:,:),&
                rRI_CI_J2p_J1(ie,:,:,:), iRI_CI_J2p_J1(ie,:,:,:), &
                rRI_CI_J0p_J1(ie,:,:,:), iRI_CI_J0p_J1(ie,:,:,:), &  !Radial integral tensors
                evec, gam, Eo, w)
            end do
            write(78,1978) ((REAL(BtmpJ1(ii1,ii2)), ii2=1,4),ii1=1,nt)
            write(178,1978) ((AIMAG(BtmpJ1(ii1,ii2)), ii2=1,4),ii1=1,nt)
            BwJ1(ie,:,:)=BtmpJ1(:,:)
        end do
        close(78)
        !Computations for J=3 
        Jt = 3
        Mt = 0
        nt = 6
        do ie=1,nconts
            BtmpJ3 = Zer
            do ii1=1,nt
                call cont_Amp_TP_vec_SI(nt, Elo+ie*dEcont, ii1, Jt, Mt, BtmpJ3(ii1,:),&
                Zo, Ebs,& !Initial bound state info
                nbs, SelJ2f,EnJ2f, ZJ2f,&
                nbs, SelJ2p, EnJ2p, ZJ2p, &
                nbs0, SelJ0, EnJ0, ZJ0,& !Intermediate state information
                RI_OI_J1_J2f, RI_OI_J1_J2f, RI_OI_J1_J0p, &
                rRI_CI_J2f_J3(ie,:,:,:), iRI_CI_J2f_J3(ie,:,:,:),&
                rRI_CI_J2p_J3(ie,:,:,:), iRI_CI_J2p_J3(ie,:,:,:), & 
                rRI_CI_J0p_J3(ie,:,:,:), iRI_CI_J0p_J3(ie,:,:,:), &!Radial integral tensors
                evec, gam, Eo, w)
            end do
            write(79,1979) ((REAL(BtmpJ3(ii1,ii2)), ii2=1,4),ii1=1,nt)
            write(179,1979) ((AIMAG(BtmpJ3(ii1,ii2)), ii2=1,4),ii1=1,nt)
            BwJ3(ie, :, :) = BtmpJ3(:,:)
        end do
        close(79)
        close(179)
        close(78)
        close(178)
        1978 format(5(4(E20.9E3,1x),3x))
        1979 format(6(4(E20.9E3,1x),3x))
    end subroutine 
    !
    subroutine CoreAllIs()
        integer :: ie
        integer :: Jt, Mt, nt
        integer ::  ii2
        integer :: SelJ0(2), SelJ2p(2), SelJ2f(2)
        !real*8 :: RI_CI_J0_J1(nos, nbs, 2, 5)
        complex*16 :: BtmpJ1(5, 4), BtmpJ3(6,4)
        complex*16 :: BwJ1(nconts,5,4), BwJ3(nconts, 6,4)
        SelJ0  = (/1,nbs0/)
        SelJ2p = (/1,nbs/)
        SelJ2f = (/1,nbs/)
        write(6,*) "Run with using Incoming wave wave function"
        open(78, file="ZetaMatrices/"//itoa(nl)//"AIn"//filepreface(SelJ0, SelJ2p, SelJ2f)//&
        "RFinalContJ1States.dat", status='replace')
        open(178, file="ZetaMatrices/"//itoa(nl)//"AIn"//filepreface(SelJ0, SelJ2p, SelJ2f)//&
        "IFinalContJ1States.dat", status='replace')
        open(79, file="ZetaMatrices/"//itoa(nl)//"AIn"//filepreface(SelJ0, SelJ2p, SelJ2f)//&
        "RFinalContJ3States.dat", status='replace')
        open(179, file="ZetaMatrices/"//itoa(nl)//"AIn"//filepreface(SelJ0, SelJ2p, SelJ2f)//&
        "IFinalContJ3States.dat", status='replace')
        Jt = 1
        Mt = 0
        nt = 5
        do ie=1, ncon
            BtmpJ1 = Zer
            do ii1=1,nt
                call cont_Amp_TP_vec_SI(nt, Elo+ie*dEcont, ii1, Jt, Mt, BtmpJ1(ii1,:),&
                Zo, Ebs,& !Initial bound state info
                nbs, SelJ2f, EnJ2f, ZJ2f, &
                nbs, SelJ2p, EnJ2p, ZJ2p, & !Intermediate state information
                nbs0, SelJ0, EnJ0, ZJ0, &
                RI_OI_J1_J2f, RI_OI_J1_J2f, RI_OI_J1_J0p,&
                rRI_CI_J2f_J1(ie,:,:,:), iRI_CI_J2f_J1(ie,:,:,:),&
                rRI_CI_J2p_J1(ie,:,:,:), iRI_CI_J2p_J1(ie,:,:,:), &
                rRI_CI_J0p_J1(ie,:,:,:), iRI_CI_J0p_J1(ie,:,:,:), &  !Radial integral tensors
                evec, gam, Eo, w)
            end do
            write(78,1978) ((REAL(BtmpJ1(ii1,ii2)), ii2=1,4),ii1=1,nt)
            write(178,1978) ((AIMAG(BtmpJ1(ii1,ii2)), ii2=1,4),ii1=1,nt)
            BwJ1(ie,:,:)=BtmpJ1(:,:)
        end do
        close(78)
        !Computations for J=3 
        Jt = 3
        Mt = 0
        nt = 6
        do ie=1,nconts
            BtmpJ3 = Zer
            do ii1=1,nt
                call cont_Amp_TP_vec_SI(nt, Elo+ie*dEcont, ii1, Jt, Mt, BtmpJ3(ii1,:),&
                Zo, Ebs,& !Initial bound state info
                nbs, SelJ2f,EnJ2f, ZJ2f,&
                nbs, SelJ2p, EnJ2p, ZJ2p, &
                nbs0, SelJ0, EnJ0, ZJ0,& !Intermediate state information
                RI_OI_J1_J2f, RI_OI_J1_J2f, RI_OI_J1_J0p, &
                rRI_CI_J2f_J3(ie,:,:,:), iRI_CI_J2f_J3(ie,:,:,:),&
                rRI_CI_J2p_J3(ie,:,:,:), iRI_CI_J2p_J3(ie,:,:,:), & 
                rRI_CI_J0p_J3(ie,:,:,:), iRI_CI_J0p_J3(ie,:,:,:), &!Radial integral tensors
                evec, gam, Eo, w)
            end do
            write(79,1979) ((REAL(BtmpJ3(ii1,ii2)), ii2=1,4),ii1=1,nt)
            write(179,1979) ((AIMAG(BtmpJ3(ii1,ii2)), ii2=1,4),ii1=1,nt)
            BwJ3(ie, :, :) = BtmpJ3(:,:)
        end do
        close(79)
        close(179)
        close(78)
        close(178)
        1978 format(5(4(E20.9E3,1x),3x))
        1979 format(6(4(E20.9E3,1x),3x))
    end subroutine
    !
    subroutine CoreAllIs_Re()
        integer :: ie
        integer :: Jt, Mt, nt
        integer ::  ii2
        integer :: SelJ0(2), SelJ2p(2), SelJ2f(2)
        !real*8 :: RI_CI_J0_J1(nos, nbs, 2, 5)
        complex*16 :: BtmpJ1(5, 4), BtmpJ3(6,4)
        complex*16 :: BwJ1(nconts,5,4), BwJ3(nconts, 6,4)
        SelJ0  = (/1,nbs0/)
        SelJ2p = (/1,nbs/)
        SelJ2f = (/1,nbs/)
        write(6,*) "Run with the real Eige Phaseshift wave function."
        open(78, file="ZetaMatrices/"//itoa(nl)//"AIn"//filepreface(SelJ0, SelJ2p, SelJ2f)//&
        "RFinalContJ1States.dat", status='replace')
        open(178, file="ZetaMatrices/"//itoa(nl)//"AIn"//filepreface(SelJ0, SelJ2p, SelJ2f)//&
        "IFinalContJ1States.dat", status='replace')
        open(79, file="ZetaMatrices/"//itoa(nl)//"AIn"//filepreface(SelJ0, SelJ2p, SelJ2f)//&
        "RFinalContJ3States.dat", status='replace')
        open(179, file="ZetaMatrices/"//itoa(nl)//"AIn"//filepreface(SelJ0, SelJ2p, SelJ2f)//&
        "IFinalContJ3States.dat", status='replace')
        Jt = 1
        Mt = 0
        nt = 5
        do ie=1, ncon
            BtmpJ1 = Zer
            do ii1=1,nt
                call cont_Amp_TP_vec_SI_Re(nt, Elo+ie*dEcont, ii1, Jt, Mt, BtmpJ1(ii1,:),&
                Zo, Ebs,& !Initial bound state info
                nbs, SelJ2f, EnJ2f, ZJ2f, &
                nbs, SelJ2p, EnJ2p, ZJ2p, & !Intermediate state information
                nbs0, SelJ0, EnJ0, ZJ0, &
                RI_OI_J1_J2f, RI_OI_J1_J2f, RI_OI_J1_J0p,&
                rRI_CI_J2f_J1(ie,:,:,:), iRI_CI_J2f_J1(ie,:,:,:),&
                rRI_CI_J2p_J1(ie,:,:,:), iRI_CI_J2p_J1(ie,:,:,:), &
                rRI_CI_J0p_J1(ie,:,:,:), iRI_CI_J0p_J1(ie,:,:,:), &  !Radial integral tensors
                evec, gam, Eo, w)
            end do
            write(78,1978) ((REAL(BtmpJ1(ii1,ii2)), ii2=1,4),ii1=1,nt)
            write(178,1978) ((AIMAG(BtmpJ1(ii1,ii2)), ii2=1,4),ii1=1,nt)
            BwJ1(ie,:,:)=BtmpJ1(:,:)
        end do
        close(78)
        !Computations for J=3 
        Jt = 3
        Mt = 0
        nt = 6
        do ie=1,nconts
            BtmpJ3 = Zer
            do ii1=1,nt
                call cont_Amp_TP_vec_SI_Re(nt, Elo+ie*dEcont, ii1, Jt, Mt, BtmpJ3(ii1,:),&
                Zo, Ebs,& !Initial bound state info
                nbs, SelJ2f,EnJ2f, ZJ2f,&
                nbs, SelJ2p, EnJ2p, ZJ2p, &
                nbs0, SelJ0, EnJ0, ZJ0,& !Intermediate state information
                RI_OI_J1_J2f, RI_OI_J1_J2f, RI_OI_J1_J0p, &
                rRI_CI_J2f_J3(ie,:,:,:), iRI_CI_J2f_J3(ie,:,:,:),&
                rRI_CI_J2p_J3(ie,:,:,:), iRI_CI_J2p_J3(ie,:,:,:), & 
                rRI_CI_J0p_J3(ie,:,:,:), iRI_CI_J0p_J3(ie,:,:,:), &!Radial integral tensors
                evec, gam, Eo, w)
            end do
            write(79,1979) ((REAL(BtmpJ3(ii1,ii2)), ii2=1,4),ii1=1,nt)
            write(179,1979) ((AIMAG(BtmpJ3(ii1,ii2)), ii2=1,4),ii1=1,nt)
            BwJ3(ie, :, :) = BtmpJ3(:,:)
        end do
        close(79)
        close(179)
        close(78)
        close(178)
        1978 format(5(4(E20.9E3,1x),3x))
        1979 format(6(4(E20.9E3,1x),3x))
    end subroutine
    !
    subroutine Rabi()
        complex*16 :: RJ2f(nos, nbs), RJ2p(nos, nbs), RJ0(nos,nbs0)
        open(301, file="Rabi_Freqs/Rabi_J2f.dat")
        open(302, file="Rabi_Freqs/Rabi_J2p.dat")
        open(303, file="Rabi_Freqs/Rabi_J0p.dat")
        call Rabi_freq(nbs, 3, ZJ2f, 2, 0, 3, RI_OI_J1_J2f,& 
        nos, 5,  Zo, RJ2f, evec)
        call Rabi_freq(nbs, 3, ZJ2p, 2, 0, 1, RI_OI_J1_J2p,& 
        nos, 5,  Zo, RJ2p, evec)
        call Rabi_freq(nbs0, 2, ZJ0, 0, 0, 1, RI_OI_J1_J0p,& 
        nos, 5,  Zo, RJ0, evec)
        do i1=1,4
            write(301,1301) (REAL(RJ2f(i1,i2)),i2=1,24)
            write(302,1301) (REAL(RJ2p(i1,i2)),i2=1,24)
            write(303,1303) (REAL(RJ0(i1,i2)) ,i2=1,13)
        enddo
        close(301)
        close(302)
        close(303)
        1301 format(24(E20.9E3,2x))
        1303 format(13(E20.9E3,2x))
    end subroutine 
    !
    function fileZpreface(nJ0, nJ2p, nJ2f) result(prefac)
        integer :: nJ0, nJ2p, nJ2f
        character(LEN=:), ALLOCATABLE :: prefac
        prefac = ''
        if(nJ0.gt.0) then
            prefac = prefac//"_J0p"//itoa(nJ0)//"_"
        endif
        if(nJ2p.gt.0) then
            prefac = prefac//"_J2p"//itoa(nJ2p)//"_"
        endif
        if(nJ2f.gt.0) then
            prefac = prefac//"_J2f"//itoa(nJ2f)//"_"
        endif
    end function 
    !
    function filepreface(SelJ0, SelJ2p, SelJ2f) result(prefac)
        integer:: SelJ0(2), SelJ2p(2), SelJ2f(2)
        character(LEN=:), allocatable :: prefac
        prefac = ''
        if(SelJ0(1).gt.0) then
            prefac = prefac//"_J0p"//itoa(SelJ0(1))//"-"//itoa(SelJ0(2))//"_"
        endif
        if(SelJ2p(1).gt.0) then
            prefac = prefac//"_J2p"//itoa(SelJ2p(1))//"-"//itoa(SelJ2p(2))//"_"
        endif
        if(SelJ2f(1).gt.0) then
            prefac = prefac//"_J2f"//itoa(SelJ2f(1))//"-"//itoa(SelJ2f(2))//"_"
        endif
        !allocate(prefac(len(PJ2f)+len(PJ2p)+len(PJ0)))
        !write(prefac,*) PJ0, PJ2p, PJ2f
        !prefac = adjustr(prefac))
    end function
    !
    function itoa(ii) result(res)
        character(:),allocatable :: res
        integer,intent(in) :: ii
        character(range(ii)+2) :: tmp
        write(tmp,'(i0)') ii
        res = trim(tmp)
    end function
end program main
!
function kij(i,j)
    implicit none
    integer, intent(in) :: i,j
    integer :: kij
    kij = int(float((i+j)-abs(i-j))/float((i+j)+abs(i-j)))
    return
end function kij
!
subroutine Init_big_mats(nnbs, nnconts, Rs, chanc, chanb, nn,file_name)
    implicit none
    integer :: nnbs, nnconts, chanc, chanb, ii, jj, kk, ll,nn
    real*8 :: Rs(nnconts, nnbs, chanc, chanb), Rtmp(chanc*chanb)
    character(len=nn) :: file_name
    character(len=:), ALLOCATABLE :: fi
    !character(len=2) :: indx
    !write(6,*) shape(Rs)
    do jj=1, nnbs
        !write(6,*) TRIM(ADJUSTL(indx)), ADJUSTL(indx)
        fi = file_name//itoa(jj-1)//'.dat'
        !write(6,*) fi
        open(15, file=fi, status='old', access='sequential')
    do ii=1, nnconts        
        read(15,*) Rtmp
        do ll=1,chanc
            do kk=1, chanb
                !if((ii.eq.1).and.(jj.eq.1)) write(6,*) ll, kk, kk+(ll-1)*chanb
                Rs(ii, jj, ll, kk) = Rtmp(kk+(ll-1)*chanb)
            enddo
        enddo
    end do
        close(15)
    end do
    contains
    !
    function itoa(ii1) result(res)
        character(:),allocatable :: res
        integer,intent(in) :: ii1
        character(range(ii1)+2) :: tmp
        write(tmp,'(i0)') ii1
        res = trim(tmp)
    end function
end subroutine Init_big_mats

subroutine Init_big_matsW(nnbs, nnconts, Rs, chanc, chanb, nn,file_name)
    implicit none
    integer :: nnbs, nnconts, chanc, chanb, ii, jj, kk, ll,nn
    real*8 :: Rs(nnbs, nnconts, chanc, chanb), Rtmp(chanc*chanb)
    character(len=nn) :: file_name
    character(len=:), allocatable :: fi
    !character(len=2) :: indx
    write(6,*) "INside init big matsW ", nnbs, nnconts, chanc, chanb
    write(6,*) shape(Rs)
    do jj=1, nnbs
        !write(indx,'(I2)') jj-1
        !write(6,*) TRIM(ADJUSTL(indx)), ADJUSTL(indx)
        fi = file_name//itoa(jj-1)//'.dat'
        write(6,*) fi
        open(15, file=fi, status='old', access='sequential')
    do ii=1, nnconts        
        read(15,*) Rtmp
        do ll=1,chanc
            do kk=1, chanb
                !if((ii.eq.1).and.(jj.eq.1)) write(6,*) ll, kk, kk+(ll-1)*chanb
                Rs(jj, ii, ll, kk) = Rtmp(kk+(ll-1)*chanb)
            enddo
        enddo
    end do
        close(15)
    end do
    contains
    !
    function itoa(ii1) result(res)
        character(:),allocatable :: res
        integer,intent(in) :: ii1
        character(range(ii1)+2) :: tmp
        write(tmp,'(i0)') ii1
        res = trim(tmp)
    end function
end subroutine Init_big_matsW
