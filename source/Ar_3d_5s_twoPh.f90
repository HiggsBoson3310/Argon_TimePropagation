program main
    use j_recoup
    use ArJ1J2
    use ERRFUN
    implicit none
    integer, parameter :: nbs=24, nbs0=13, nconts=200, nos = 4
    integer :: ncon, mm, kk, ie, ntd
    integer :: nl
    integer :: Jt, Mt, nt
    integer :: ii1, iit, ii2
    real*8 :: ProbsJ12(200), ProbsJ32(200)
    real*8 :: Elo, dEcont, Ebs(4)
    real*8 :: tdlo, tdf, dtd
    real*8 :: Zo(4, 5)
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
    real*8 :: Eo, w, gam, evec(3), to
    real*8 :: Ar, Ai
    !real*8 :: RI_CI_J0_J1(nos, nbs, 2, 5)
    character(len=143) :: file_root
    complex*16 :: BJ1(5), BJ3(6), A(4), BtmpJ1(5, 4), BtmpJ3(6,4)

    !call InitialStates
    !-------------------------------- calling set up from gensub -----------------------
    call SETUP
    !-----------------------------------------------------------------------------------
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
        call Init_big_matsW(nos, nbs, RI_OI_J1_J2f, 3, 5, file_root)
        ! p states
        !J=2
        file_root = 'Python_Code_RadInts/Rad Ints bb for J=2 lint p/interm_rad_ints_'
        call Init_big_mats(nos, nbs, RI_OI_J1_J2p, 3, 5, file_root)
        !J=0
        file_root = 'Python_Code_RadInts/Rad Ints bb for J=0 lint p/interm_rad_ints_'
        call Init_big_mats(nos, nbs0, RI_OI_J1_J0p, 2, 5, file_root)
        ! Radial elements to J=1 continuum, from an intermediate J=2 state
        ! f states
        file_root = 'Python_Code_RadInts/RI b(J=2,l=f) c(J=1)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs, nconts, rRI_CI_J2f_J1, 5, 3, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=2,l=f) c(J=1)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs, nconts, iRI_CI_J2f_J1, 5, 3, file_root)
        ! p states
        !J=2
        file_root = 'Python_Code_RadInts/RI b(J=2,l=p) c(J=1)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs, nconts, rRI_CI_J2p_J1, 5, 3, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=2,l=p) c(J=1)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs, nconts, iRI_CI_J2p_J1, 5, 3, file_root)
        !J=0
        file_root = 'Python_Code_RadInts/RI b(J=0,l=p) c(J=1)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs0, nconts, rRI_CI_J0p_J1, 5, 2, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=0,l=p) c(J=1)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs0, nconts, iRI_CI_J0p_J1, 5, 2, file_root)

        ! Radial elements to J=3 continuum, from an intermediate J=2 state
        ! f states
        file_root = 'Python_Code_RadInts/RI b(J=2,l=f) c(J=3)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs, nconts, rRI_CI_J2f_J3, 6, 3, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=2,l=f) c(J=3)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs, nconts, iRI_CI_J2f_J3, 6, 3, file_root)
        ! p states
        !J=2
        file_root = 'Python_Code_RadInts/RI b(J=2,l=p) c(J=3)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs, nconts, rRI_CI_J2p_J3, 6, 3, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=2,l=p) c(J=3)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs, nconts, iRI_CI_J2p_J3, 6, 3, file_root)
        !J=0
        rRI_CI_J0p_J3=0.0d0
        iRI_CI_J0p_J3=0.0d0

    !

    !Laser parameters
    open(98, file='laser_params_two_photon.dat', status='old')
    read(98,*) Eo, w, gam
    read(98,*) evec(:)
    read(98,*) nl
    close(98)
    print *, gam, Eo
    !Initial state amplitudes
    open(98, file='init_amplitudes.dat', status='old')
    write(6,*) "don't mind me pritin"
    do ii1=1,4
        read(98,*) Ar, Ai
        write(6,*) Ar, Ai
        A(ii1) = Ar*Zone+Im*Ai
        write(6,*) A(ii1)
    enddo
    close(98)
    !Continuum state energies.
    open(08,file='Python_Code_RadInts/EcontEV_det.dat')
    read(08,*) ncon, Elo, dEcont
    close(08)
    !call BoundBoundMatrices(4, 5, Zo, 24, 3, ZJ2f, 2, (/3/2.,3/2.,1/2./), (/2,1,1/), (/3,3,3/), RI_OI_J1_J2f)
    !call BoundBoundMatrices(4, 5, Zo, 24, 3, ZJ2p, 2, (/3/2.,3/2.,1/2./), (/2,1,1/), (/1,1,1/), RI_OI_J1_J2p)
    !call BoundBoundMatrices(4, 5, Zo, 13, 2, ZJ0,  0, (/1/2.,3/2./), (/1,1/), (/1,1/), RI_OI_J1_J0p)
    !Timde delay grid


    
    ntd = 1800
    tdlo = 0.00d0
    tdf = 7.5d4 
    dtd = (tdf-tdlo)/ntd
    open(10, file=""//char(nl+48)//"tdel_probsJ12_twophoton.dat")
    open(11, file=""//char(nl+48)//"tdel_probsJ32_twophoton.dat")
    
    do iit=0,ntd
        write(6,*) iit
        ProbsJ12 = 0.0d0
        ProbsJ32 = 0.0d0
        to = tdlo + iit*dtd
        if(iit.eq.0) then
            open(81, file=""//char(nl+48)//"TwoPhoton_J1_Channel_probabilities.dat") 
            open(83, file=""//char(nl+48)//"TwoPhoton_J3_Channel_probabilities.dat")
        endif
        !Computation of J=1
        Jt = 1
        Mt = 0
        nt = 5
        do ie=1,ncon
            BJ1 = Zer
            do ii1=1,5
                call cont_Amp_TP(nt, Elo+ie*dEcont, ii1, Jt, Mt, BJ1(ii1),&
                A, Zo, Ebs,& !Initial bound state info
                nbs, EnJ2f, ZJ2f, EnJ2p, ZJ2p, & !Intermediate state information
                nbs0, EnJ0, ZJ0, &
                RI_OI_J1_J2f, RI_OI_J1_J2f, RI_OI_J1_J0p,&
                rRI_CI_J2f_J1(ie,:,:,:), iRI_CI_J2f_J1(ie,:,:,:),&
                rRI_CI_J2p_J1(ie,:,:,:), iRI_CI_J2p_J1(ie,:,:,:), &
                rRI_CI_J0p_J1(ie,:,:,:), iRI_CI_J0p_J1(ie,:,:,:), &  !Radial integral tensors
                evec, gam, Eo, w, to)
                if((ii1.eq.1).or.(ii1.eq.4)) then
                    ProbsJ12(ie) = ProbsJ12(ie)+ABS(BJ1(ii1))**2
                else
                    ProbsJ32(ie) = ProbsJ32(ie)+ABS(BJ1(ii1))**2
                endif
            end do
            !write(6,*) Elo+ie*dEcont, Bj1
            if(iit.eq.0) write(81, 1981) Elo+ie*dEcont, (ABS(BJ1(mm))**2, mm=1, nt)
        end do
        !Computations for J=3 
        Jt = 3
        Mt = 0
        nt = 6
        do ie=1,nconts
            BJ3 = Zer
            do ii1=1,nt
                call cont_Amp_TP(nt, Elo+ie*dEcont, ii1, Jt, Mt, BJ3(ii1),&
                A, Zo, Ebs,& !Initial bound state info
                nbs, EnJ2f, ZJ2f, EnJ2p, ZJ2p, &
                nbs0, EnJ0, ZJ0,& !Intermediate state information
                RI_OI_J1_J2f, RI_OI_J1_J2f, RI_OI_J1_J0p, &
                rRI_CI_J2f_J3(ie,:,:,:), iRI_CI_J2f_J3(ie,:,:,:),&
                rRI_CI_J2p_J3(ie,:,:,:), iRI_CI_J2p_J3(ie,:,:,:), & 
                rRI_CI_J0p_J3(ie,:,:,:), iRI_CI_J0p_J3(ie,:,:,:), &!Radial integral tensors
                evec, gam, Eo, w, to)
                if((ii1.eq.3).or.(ii1.eq.6)) then
                    ProbsJ12(ie) = ProbsJ12(ie)+ABS(BJ3(ii1))**2
                else
                    ProbsJ32(ie) = ProbsJ32(ie)+ABS(BJ3(ii1))**2
                endif
            end do
            if(iit.eq.0) write(83, 1983) Elo+ie*dEcont, (ABS(BJ3(mm))**2, mm=1, nt)
        end do
        if(iit.eq.-12) then
            close(81)
            close(83)
        endif
        write(10, 1910) to, (ProbsJ12(mm), mm=1,nconts)
        write(11, 1910) to, (ProbsJ32(mm), mm=1,nconts)
    end do
    close(10)
    close(11)
 1981 format(6(E20.9E3, 3x))
 1983 format(7(E20.9E3, 3x))
 1910 format(201(E20.9E3, 3x))
    !Testing the beats with the new amplitude.
        !open(45, file = 'test_TP.dat')
        !do iit=0,168
        !    ff = Zer
        !    do ii1=1,4
        !        ff = ff+twoPhotonAmp(Ebs(ii1)/CMperAU, (Elo+60*dEcont)/CMperAU,EnJ2f(10)/CMperAU,gam,tdlo+iit*dtd,w)
        !    enddo
        !    write(45,*) tdlo+iit*dtd, ABS(ff)**2
        !end do
        !close(45)
    !stop
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
subroutine Init_big_mats(nnbs, nnconts, Rs, chanc, chanb, file_name)
    implicit none
    integer :: nnbs, nnconts, chanc, chanb, ii, jj, kk, ll
    real*8 :: Rs(nnconts, nnbs, chanc, chanb), Rtmp(chanc*chanb)
    character(len=143) :: file_name
    character(len=149) :: fi
    character(len=2) :: indx
    !write(6,*) shape(Rs)
    do jj=1, nnbs
        write(indx,'(I2)') jj-1
        !write(6,*) TRIM(ADJUSTL(indx)), ADJUSTL(indx)
        fi = TRIM(file_name)//TRIM(ADJUSTL(indx))//'.dat'
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
end subroutine Init_big_mats

subroutine Init_big_matsW(nnbs, nnconts, Rs, chanc, chanb, file_name)
    implicit none
    integer :: nnbs, nnconts, chanc, chanb, ii, jj, kk, ll
    real*8 :: Rs(nnbs, nnconts, chanc, chanb), Rtmp(chanc*chanb)
    character(len=143) :: file_name
    character(len=149) :: fi
    character(len=2) :: indx
    !write(6,*) shape(Rs)
    do jj=1, nnbs
        write(indx,'(I2)') jj-1
        !write(6,*) TRIM(ADJUSTL(indx)), ADJUSTL(indx)
        fi = TRIM(file_name)//TRIM(ADJUSTL(indx))//'.dat'
        !write(6,*) fi
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
end subroutine Init_big_matsW
