program main
    use j_recoup
    use ArJ1J2
    use ERRFUN
    implicit none
    integer, parameter :: nbs=24, nbs0=13, nconts=200, nos = 4
    logical :: exp_flag
    integer :: ncon, mm, kk, ie, ntd
    integer :: nl
    integer :: Jt, Mt, nt
    integer :: ii1, iit, ii2
    real*8 :: ProbsJ12(200), ProbsJ32(200)
    real*8 :: Elo, dEcont, Ebs(4), Econt
    real*8 :: tdlo, tdf, dtd
    real*8 :: Zo(4, 5)
    real*8 :: ZJ0(nbs0,2), EnJ0(nbs0)
    real*8 :: ZJ2f(nbs,3), EnJ2f(nbs)
    real*8 :: ZJ2p(nbs,3), EnJ2p(nbs)
    real*8 :: rRI_CI_J0p_J1(nconts,nbs0,5,2), rRI_CI_J0p_J3(nconts,nbs0,6,2)
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
    complex*16 :: BJ1(5), BJ3(6), A(4), BtmpJ1(5, 4), BtmpJ3(6,4), Phases(4)
    complex*16 :: BwJ1(nconts,5,4), BwJ3(nconts, 6,4)
    !-------------------------------- calling set up from gensub -----------------------
    call SETUP
    !-----------------------------------------------------------------------------------
    
    !---------------------- Call if initial states are not calculated ------------------
    !call InitialStates
    !-----------------------------------------------------------------------------------
    
    !-------------------------- Define flag for the use of exp energies ----------------
    exp_flag = .true.
    !-----------------------------------------------------------------------------------

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
        call Init_big_matsW(nos, nbs, RI_OI_J1_J2p, 3, 5, file_root)
        !J=0
        file_root = 'Python_Code_RadInts/Rad Ints bb for J=0 lint p/interm_rad_ints_'
        call Init_big_matsW(nos, nbs0, RI_OI_J1_J0p, 2, 5, file_root)
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
    
    if(exp_flag) then 
        open(98, file='init_amplitudes_expE.dat', status='old')
        write(6,*) "The following amplitudes and energies will be used in each state 5s, 3d, 5s'  and 3d'"
        do ii1=1,4
            read(98,*) Ar, Ai, Ebs(ii1)
            write(6,*) Ar, Ai, Ebs(ii1)
            A(ii1) = Ar*Zone+Im*Ai
            write(6,*) A(ii1)
        enddo
    else
        open(98, file='init_amplitudes.dat', status='old')
        write(6,*) "The following amplitudes will be used in each state 5s, 3d, 5s'  and 3d'"
        do ii1=1,4
            read(98,*) Ar, Ai
            write(6,*) Ar, Ai
            A(ii1) = Ar*Zone+Im*Ai
            write(6,*) A(ii1)
        enddo
    endif
    close(98)
    
    !Continuum state energies.
    open(08,file='Python_Code_RadInts/EcontEV_det.dat')
    read(08,*) ncon, Elo, dEcont
    close(08)
    open(78, file="FinalContJ1States_Re.dat", status='replace')
    open(68, file="FinalContJ1States_Im.dat", status='replace')
    open(79, file="FinalContJ3States_Re.dat", status='replace')
    open(69, file="FinalContJ3States_Im.dat", status='replace')
    Jt = 1
    Mt = 0
    nt = 5
    do ie=1,ncon
        BtmpJ1 = Zer
        do ii1=1,nt
            call cont_Amp_TP_vec(nt, Elo+ie*dEcont, ii1, Jt, Mt, BtmpJ1(ii1,:),&
            Zo, Ebs,& !Initial bound state info
            nbs, EnJ2f, ZJ2f, EnJ2p, ZJ2p, & !Intermediate state information
            nbs0, EnJ0, ZJ0, &
            RI_OI_J1_J2f, RI_OI_J1_J2f, RI_OI_J1_J0p,&
            rRI_CI_J2f_J1(ie,:,:,:), iRI_CI_J2f_J1(ie,:,:,:),&
            rRI_CI_J2p_J1(ie,:,:,:), iRI_CI_J2p_J1(ie,:,:,:), &
            rRI_CI_J0p_J1(ie,:,:,:), iRI_CI_J0p_J1(ie,:,:,:), &  !Radial integral tensors
            evec, gam, Eo, w)
        end do
        write(78,1978)  ((REAL(BtmpJ1(ii1,ii2)), ii2=1,4),ii1=1,nt)
        write(68,1978) ((AIMAG(BtmpJ1(ii1,ii2)), ii2=1,4),ii1=1,nt)
        BwJ1(ie,:,:)=BtmpJ1(:,:)
    end do
    close(78)
    close(68)

    !Computations for J=3 
    Jt = 3
    Mt = 0
    nt = 6
    do ie=1,nconts
        BJ3 = Zer
        do ii1=1,nt
            call cont_Amp_TP_vec(nt, Elo+ie*dEcont, ii1, Jt, Mt, BtmpJ3(ii1,:),&
            Zo, Ebs,& !Initial bound state info
            nbs, EnJ2f, ZJ2f, EnJ2p, ZJ2p, &
            nbs0, EnJ0, ZJ0,& !Intermediate state information
            RI_OI_J1_J2f, RI_OI_J1_J2f, RI_OI_J1_J0p, &
            rRI_CI_J2f_J3(ie,:,:,:), iRI_CI_J2f_J3(ie,:,:,:),&
            rRI_CI_J2p_J3(ie,:,:,:), iRI_CI_J2p_J3(ie,:,:,:), & 
            rRI_CI_J0p_J3(ie,:,:,:), iRI_CI_J0p_J3(ie,:,:,:), &!Radial integral tensors
            evec, gam, Eo, w)
        end do
        write(79,1979)  ((REAL(BtmpJ3(ii1,ii2)), ii2=1,4),ii1=1,nt)
        write(69,1979) ((AIMAG(BtmpJ3(ii1,ii2)), ii2=1,4),ii1=1,nt)
        BwJ3(ie, :, :) = BtmpJ3(:,:)
    end do
    close(79)
    close(69)
    
    write(6,*) "Time delay calculation started"
    ntd = 2000
    tdlo = 0.00d0
    tdf = 4.0d-12/auoftime
    dtd = (tdf-tdlo)/ntd
    write(6,*) "with ", ntd, "points at a density of ", dtd*auoftime*1.0d15, "fs"
    open(10, file=""//itoa(nl)//"tdel_probsJ12.dat")
    open(11, file=""//itoa(nl)//"tdel_probsJ32.dat")
    do iit=0,ntd
        to = tdlo+iit*dtd
        ProbsJ12 = 0.0d0
        ProbsJ32 = 0.0d0
        do ie=0, ncon
            Econt = Elo+ie*dEcont
            Phases = EXP(Im*to/CMperAU*(/Econt-Ebs/))
            BtmpJ1 = BwJ1(ie,:,:)
            BtmpJ3 = BwJ3(ie,:,:)
            !Go over the channels of J1
            Bj1 = Zer
            do ii1=1,5
                do ii2=1,4
                    if(BtmpJ1(ii1,ii2) /= BtmpJ1(ii1,ii2)) then 
                        print*, "Nan in Btmp" 
                        stop
                    endif
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
    contains
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
   