program main
    use j_recoup
    use ArJ1J2
    use ERRFUN
    implicit none
    integer, parameter :: nconts=200, nos = 4
    integer :: Na
    integer :: mm, kk, ii1
    integer :: nl, ncon
    real*8 :: Ar, Ai
    complex*16 :: A(4)
    character(len=:), allocatable :: file_root
    real*8 :: Eo, w, gam, evec(3)
    real*8 :: Zo(4, 5), Ebs(4)
    real*8 :: rRI_CI_J0p_J1(nconts, 4, 5, 2)
    real*8 :: iRI_CI_J0p_J1(nconts, 4, 5, 2)
    real*8 :: rRI_CI_J2f_J1(nconts, 4, 5, 3)
    real*8 :: iRI_CI_J2f_J1(nconts, 4, 5, 3)
    real*8 :: rRI_CI_J2p_J1(nconts, 4, 5, 3)
    real*8 :: iRI_CI_J2p_J1(nconts, 4, 5, 3)
    real*8 ::  Elo, dEcont
    Na=-5
    !-------------------------------- calling set up from gensub -----------------------
    call SETUP
    !-----------------------------------------------------------------------------------
    call InitialStates
    stop
    ! Call if the we don't hacve data on the initial bound state and the initial 
    !call InitialStates
    !Retrieving from precomputed files.
    !Loading the Zcoefficients
    open(10, file='init_bounds_J1.dat')
    do mm=1,4
        read(10,*) Ebs(mm), (Zo(mm,kk), kk=1,5)
    end do
    close(10)
    !
    !Loading the radial elements
        write(6,*) "Loading Radial Elements into four dimensional arrays"
        !
        ! Continuum J=2 states
        ! f states
        file_root = 'Python_Code_RadInts/RI b(J=1) c(J2f)/cont_rad_ints_regb'
        call Init_big_mats(nos, nconts, rRI_CI_J2f_J1, 5, 3, len(file_root), file_root)
        call Init_big_mats(nos, nconts, iRI_CI_J2f_J1, 5, 3, len(file_root), file_root)
        ! p states
        file_root = 'Python_Code_RadInts/RI b(J=1) c(J2p)/cont_rad_ints_regb'
        call Init_big_mats(nos, nconts, rRI_CI_J2p_J1, 5, 3, len(file_root), file_root)
        call Init_big_mats(nos, nconts, iRI_CI_J2p_J1, 5, 3, len(file_root), file_root)
        !
        !Continuum J=0 states.
        ! p states
        file_root = 'Python_Code_RadInts/RI b(J=1) c(J0p)/cont_rad_ints_regb'
        call Init_big_mats(nos, nconts, rRI_CI_J0p_J1, 5, 2, len(file_root), file_root)
        call Init_big_mats(nos, nconts, iRI_CI_J0p_J1, 5, 2, len(file_root), file_root)
    !

    !Laser parameters
    open(98, file='laser_params_one_photon.dat', status='old')
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

    call Core_Re()
    
    contains
    subroutine Core_Re()
        integer :: ie
        integer :: Jt, Mt, nt
        integer ::  ii2
        !real*8 :: RI_CI_J0_J1(nos, nbs, 2, 5)
        complex*16 :: BtmpJ2f(3, 4), BtmpJ2p(3,4), BtmpJ0(2,4)
        !complex*16 :: BwJ2f(nconts,3,4),  BwJ2p(nconts,3,4), BwJ0(nconts, 2,4)
        real*8 :: to
        write(6,*) "Running calculation with Eigenphase wave functions"
        open(77,  file="OnePhotonAmps/"//itoa(nl)//"SI_RFinalContJ2fStates.dat", status='replace')
        open(177, file="OnePhotonAmps/"//itoa(nl)//"SI_IFinalContJ2fStates.dat", status='replace')
        open(78,  file="OnePhotonAmps/"//itoa(nl)//"SI_RFinalContJ2pStates.dat", status='replace')
        open(178, file="OnePhotonAmps/"//itoa(nl)//"SI_IFinalContJ2pStates.dat", status='replace')
        open(79,  file="OnePhotonAmps/"//itoa(nl)//"SI_RFinalContJ0States.dat", status='replace')
        open(179, file="OnePhotonAmps/"//itoa(nl)//"SI_IFinalContJ0States.dat", status='replace')
        
        Jt = 2
        Mt = 0
        nt = 3
        do ie=1, ncon
            BtmpJ2f = Zer
            BtmpJ2p = Zer
            do ii1=1,nt
                call cont_Amp_Re(nt, Elo+ie*dEcont, ii1, 3, Jt, 0, A, Zo, Ebs, evec, BtmpJ2f(ii1,:),&
                rRI_CI_J2f_J1(ie,:,:,:), iRI_CI_J2f_J1(ie,:,:,:), gam, Eo, w, to)

                call cont_Amp_Re(nt, Elo+ie*dEcont, ii1, 1, Jt, 0, A, Zo, Ebs, evec, BtmpJ2p(ii1,:),&
                rRI_CI_J2p_J1(ie,:,:,:), iRI_CI_J2p_J1(ie,:,:,:), gam, Eo, w, to)
            end do
            write(77,1977) ((REAL(BtmpJ2f(ii1,ii2)), ii2=1,4),ii1=1,nt)
            write(177,1977) ((AIMAG(BtmpJ2f(ii1,ii2)), ii2=1,4),ii1=1,nt)
            write(78,1977) ((REAL(BtmpJ2f(ii1,ii2)), ii2=1,4),ii1=1,nt)
            write(178,1977) ((AIMAG(BtmpJ2f(ii1,ii2)), ii2=1,4),ii1=1,nt)
        end do
        !Computations for J=0
        Jt = 0
        Mt = 0
        nt = 2
        do ie=1,nconts
            BtmpJ0 = Zer
            do ii1=1,nt
                call cont_Amp_Re(nt, Elo+ie*dEcont, ii1, 3, Jt, 0, A, Zo, Ebs, evec, BtmpJ0(ii1,:),&
                rRI_CI_J2f_J1(ie,:,:,:), iRI_CI_J2f_J1(ie,:,:,:), gam, Eo, w, to)
            end do
            write(79,1979) ((REAL(BtmpJ0(ii1,ii2)), ii2=1,4),ii1=1,nt)
            write(179,1979) ((AIMAG(BtmpJ0(ii1,ii2)), ii2=1,4),ii1=1,nt)
        end do
        close(77)
        close(177)
        close(78)
        close(178)
        close(79)
        close(179)
        1977 format(3(4(E20.9E3,1x),3x))
        1979 format(2(4(E20.9E3,1x),3x))
    end subroutine 
    !
    subroutine Core_Z()
        integer :: ie
        integer :: Jt, Mt, nt
        integer ::  ii2
        !real*8 :: RI_CI_J0_J1(nos, nbs, 2, 5)
        complex*16 :: BtmpJ2f(3, 4), BtmpJ2p(3,4), BtmpJ0(2,4)
        !complex*16 :: BwJ2f(nconts,3,4),  BwJ2p(nconts,3,4), BwJ0(nconts, 2,4)
        real*8 :: to
        write(6,*) "Running calculation with Incoming wave functions"
        open(77,  file="OnePhotonAmps/"//itoa(nl)//"SI_RFinalContJ2fStates.dat", status='replace')
        open(177, file="OnePhotonAmps/"//itoa(nl)//"SI_IFinalContJ2fStates.dat", status='replace')
        open(78,  file="OnePhotonAmps/"//itoa(nl)//"SI_RFinalContJ2pStates.dat", status='replace')
        open(178, file="OnePhotonAmps/"//itoa(nl)//"SI_IFinalContJ2pStates.dat", status='replace')
        open(79,  file="OnePhotonAmps/"//itoa(nl)//"SI_RFinalContJ0States.dat", status='replace')
        open(179, file="OnePhotonAmps/"//itoa(nl)//"SI_IFinalContJ0States.dat", status='replace')
        
        Jt = 2
        Mt = 0
        nt = 3
        do ie=1, ncon
            BtmpJ2f = Zer
            BtmpJ2p = Zer
            do ii1=1,nt
                call cont_Amp_Z(nt, Elo+ie*dEcont, ii1, 3, Jt, 0, A, Zo, Ebs, evec, BtmpJ2f(ii1,:),&
                rRI_CI_J2f_J1(ie,:,:,:), iRI_CI_J2f_J1(ie,:,:,:), gam, Eo, w, to)

                call cont_Amp_Z(nt, Elo+ie*dEcont, ii1, 1, Jt, 0, A, Zo, Ebs, evec, BtmpJ2p(ii1,:),&
                rRI_CI_J2p_J1(ie,:,:,:), iRI_CI_J2p_J1(ie,:,:,:), gam, Eo, w, to)
            end do
            write(77,1977) ((REAL(BtmpJ2f(ii1,ii2)), ii2=1,4),ii1=1,nt)
            write(177,1977) ((AIMAG(BtmpJ2f(ii1,ii2)), ii2=1,4),ii1=1,nt)
            write(78,1977) ((REAL(BtmpJ2f(ii1,ii2)), ii2=1,4),ii1=1,nt)
            write(178,1977) ((AIMAG(BtmpJ2f(ii1,ii2)), ii2=1,4),ii1=1,nt)
        end do
        !Computations for J=0
        Jt = 0
        Mt = 0
        nt = 2
        do ie=1,nconts
            BtmpJ0 = Zer
            do ii1=1,nt
                call cont_Amp_Re(nt, Elo+ie*dEcont, ii1, 3, Jt, 0, A, Zo, Ebs, evec, BtmpJ0(ii1,:),&
                rRI_CI_J2f_J1(ie,:,:,:), iRI_CI_J2f_J1(ie,:,:,:), gam, Eo, w, to)
            end do
            write(79,1979) ((REAL(BtmpJ0(ii1,ii2)), ii2=1,4),ii1=1,nt)
            write(179,1979) ((AIMAG(BtmpJ0(ii1,ii2)), ii2=1,4),ii1=1,nt)
        end do
        close(77)
        close(177)
        close(78)
        close(178)
        close(79)
        close(179)
        1977 format(3(4(E20.9E3,1x),3x))
        1979 format(2(4(E20.9E3,1x),3x))
    end subroutine  
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