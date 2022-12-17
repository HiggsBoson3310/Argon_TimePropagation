program main
    use j_recoup
    use ArJ1J2
    use ERRFUN
    use spf_prop
    use matmod
    implicit none
    integer, parameter :: nconts=200, nred=nos+nbs, nred1 = nos1+nbs2p+nbs2f+nbs0
    logical :: exp_flag
    integer :: ncon, mm, kk
    integer :: nl
    integer :: mi, mi1, mi2, ii1
    integer :: ieo, igam
    real*8 :: Elo, dEcont, EJ1(nos1), EJ3(nos3)
    real*8 :: ZJ1(nos1, 5), ZJ3(nos3,6)
    real*8 :: ZJ0(nbs0,2), EnJ0(nbs0)
    real*8 :: ZJ2f(nbs2f,3), EnJ2f(nbs2f)
    real*8 :: ZJ2p(nbs2p,3), EnJ2p(nbs2p)
    real*8 :: ZJ4l(nbs4,6), EnJ4l(nbs4)
    real*8 :: RI_OI_J1_J2f(nos1, nbs2f, 3, 5), RI_OI_J1_J2p(nos1, nbs2p, 3,5)
    real*8 :: RI_OI_J1_J0p(nos1, nbs0, 2, 5)
    real*8 :: RI_OI_J3_J2f(nos3, nbs2f, 3, 6), RI_OI_J3_J2p(nos3, nbs2p, 3,6)
    real*8 :: RI_OI_J3_J4l(nos3, nbs4, 6, 6)
    real*8 :: Eo, w, gam, evec(3)
    real*8 :: Ar, Ai
    complex*16 :: A(4)
    complex*16 :: H(nred, nred), V(nred,nred), c(nred), caux(nred)
    complex*16 :: HJ1(nred1,nred1), V1(nred1, nred1), c1(nred1), caux1(nred1)
    real*8 :: Es(nred), HD(nred)
    real*8 :: Es1(nred1), HD1(nred1)
    integer :: nt
    real*8 :: toolo, toofo, dtoo
    integer :: noo
    !real*8 :: RI_CI_J0_J1(nos, nbs2, 2, 5)
    character(len=143) :: file_root
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
        do mm=1,nos1
            read(10,*) EJ1(mm), (ZJ1(mm,kk), kk=1,5)
            !write(6,*) EJ1(mm), (ZJ1(mm,kk), kk=1,5)
        end do
        close(10)

        open(10, file='init_bounds_J3.dat')
        do mm=1,nos3
            read(10,*) EJ3(mm), (ZJ3(mm,kk), kk=1,6)
            !write(6,*) EJ3(mm), (ZJ3(mm,kk), kk=1,6)
        end do
        close(10)

        open(10, file='interm_bounds_J2f.dat')
        do mm=1,nbs2f
            read(10,*) EnJ2f(mm), (ZJ2f(mm,kk), kk=1,3)
            !write(6,*) EnJ2f(mm), (ZJ2f(mm,kk), kk=1,3)
        end do
        close(10)

        open(10, file='interm_bounds_J2p.dat')
        do mm=1,nbs2p
            read(10,*) EnJ2p(mm), (ZJ2p(mm,kk), kk=1,3)
            !write(6,*) EnJ2p(mm), (ZJ2p(mm,kk), kk=1,3)
        end do
        close(10)

        open(10, file='interm_bounds_J0p.dat')
        do mm=1,nbs0
            read(10,*) EnJ0(mm), (ZJ0(mm,kk), kk=1,2)
            !write(6,*) EnJ0(mm), (ZJ0(mm,kk), kk=1,2)
        end do
        close(10)

        open(10, file='interm_bounds_J4l.dat')
        do mm=1,nbs4
            read(10,*) EnJ4l(mm), (ZJ4l(mm,kk), kk=1,6)
            !write(6,*) EnJ4l(mm), (ZJ4l(mm,kk), kk=1,6)
        end do
        close(10)
    !
    
    ! Initializing the matrix elements to a consitent manageble value
        RI_OI_J1_J2f = 2.d3
        RI_OI_J1_J0p = 2.d3
        RI_OI_J3_J2f = 2.d3
        RI_OI_J3_J4l = 2.d3
        RI_OI_J1_J2p = 2.d3
        RI_OI_J3_J2p = 2.d3
    !Loading the radial elements
        write(6,*) "Loading Radial Elements into four dimensional arrays"
        ! Radial elements from initial J=1 state, to intermediate J=2 state
        ! f states
        file_root = 'Python_Code_RadInts/RI b(J=1)  b(J=2 lint f)/interm_rad_ints_'
        call Init_big_matsW(nos1, nbs2f, RI_OI_J1_J2f, 3, 5, file_root)
        ! p states
        !J=2
        file_root = 'Python_Code_RadInts/RI b(J=1)  b(J=2 lint p)/interm_rad_ints_'
        call Init_big_matsW(nos1, nbs2p, RI_OI_J1_J2p, 3, 5, file_root)
        !J=0
        file_root = 'Python_Code_RadInts/RI b(J=1)  b(J=0 lint p)/interm_rad_ints_'
        call Init_big_matsW(nos1, nbs0, RI_OI_J1_J0p, 2, 5, file_root)

        ! Radial elements from initial J=3 state, to intermediate J=2 state
        ! f states
        file_root = 'Python_Code_RadInts/RI b(J=3)  b(J=2 lint f)/interm_rad_ints_'
        call Init_big_matsW(nos3, nbs2f, RI_OI_J3_J2f, 3, 6, file_root)
        ! p states
        !J=2
        file_root = 'Python_Code_RadInts/RI b(J=3)  b(J=2 lint p)/interm_rad_ints_'
        call Init_big_matsW(nos3, nbs2p, RI_OI_J3_J2p, 3, 6, file_root)
        !J=0
        file_root = 'Python_Code_RadInts/RI b(J=3)  b(J=4 lint f)/interm_rad_ints_'
        call Init_big_matsW(nos3, nbs4, RI_OI_J3_J4l, 6, 6, file_root)

    !

    !Laser parameters
        open(98, file='laser_params_two_photon.dat', status='old')
        read(98,*) Eo, w, gam
        read(98,*) evec(:)
        read(98,*) nl
        close(98)
        print *, gam, Eo

    !
        
    ! Continuum paramerers
    open(98, file='Python_Code_RadInts/EcontEV_det.dat',status='old')
    read(98,*) ncon, Elo, dEcont
    close(98)
    ! Initial state amplitudes
    if(exp_flag) then 
        open(98, file='init_amplitudes_expE.dat', status='old')
        write(6,*) "The following amplitudes and energies will be used in each state 5s, 3d, 5s'  and 3d'"
        do ii1=1,4
            read(98,*) Ar, Ai, EJ1(ii1)
            write(6,*) Ar, Ai, EJ1(ii1)
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

    ! Vector with all he energies in atomic units
    Es(1:nos1) = EJ1/CMperAU
    Es(nos1+1:nos) = EJ3/CMperAU
    Es(nos+1:nbs2f+nos) = EnJ2f/CMperAU
    Es(nbs2f+nos+1: nbs2f+nbs2p+nos) = EnJ2p/CMperAU
    Es(nbs2f+nbs2p+nos+1 : nos+nbs2f+nbs2p+nbs0) = EnJ0/CMperAU
    Es(nos+nbs2f+nbs2p+nbs0+1:) = Enj4l/CMperAU

    open(15, file='SPF_no_cont/energiesJ3.dat')
    do ii1=1, nred
        write(15,*) Es(ii1)
    end do
    close(15)

    Es1(1:nos1) = EJ1/CmperAU
    Es1(nos1+1:nbs2f+nos1) = EnJ2f/CMperAU
    Es1(nbs2f+nos1+1:nbs2f+nbs2p+nos1) = EnJ2p/CMperAU
    Es1(nbs2f+nbs2p+nos1+1:) = EnJ0/CMperAU
    open(15, file='SPF_no_cont/energies.dat')
    do ii1=1, nred1
        write(15,*) Es1(ii1)
    end do
    close(15)
    
    call no_delay_J3()
    
    stop
    ! Open a couple files for the storage of the probabilities we are interested in
    !open(45, file='SPF_no_cont/SPF_prob32.dat', access='sequential')
    !open(46, file='SPF_no_cont/SPF_prob12.dat', access='sequential')
    !do mi2 = 1, noo
    !    call single_to(toolo+mi2*dtoo,gam,46,45)
    !end do
    ! Defining a test mesh 
    ! Call if you want to see an example calculation and see the evolution in time.
    !open(45, file='SPF_no_cont/testing_cs.dat', access='sequential',status='replace') 
    !call example_to(toolo+50*dtoo, gam, 45)!single_to(toolo+itoo*dtoo, gam, 45)
    !close(45)
    !close(46)
    contains
    !
    function itoa(ii) result(res)
        character(:),allocatable :: res
        integer,intent(in) :: ii
        character(range(ii)+2) :: tmp
        write(tmp,'(i0)') ii
        res = trim(tmp)
    end function
    !
    subroutine profile(thi, tlo, dtt, too, gaam)
        integer :: iloc
        real*8 :: thi, tlo, dtt, too, gaam
        open(12, file="SPF_no_cont/laser_profile.dat", status='replace', access='sequential')
        do iloc=1, int((thi-tlo)/dtt)
            write(12,*) tlo+iloc*dtt, Eo * EXP(-(tlo+iloc*dtt-too)**2/gaam**2) * COS(w*(tlo+iloc*dtt-too))
        end do
        close(12)
    end subroutine
    !
    function Las(t, dtt, too, Eoo, ww,gaam)
        real*8 :: t, dtt, tt(30), too, Eoo, ww, gaam
        real*8 :: Las
        real*8 :: fss(30)
        real*8, external :: rint
        do mi1=1,30
            tt(mi1) = t+mi1*dtt/30
            fss(mi1) =  Eoo * EXP(-(tt(mi1)-too)**2/gaam**2) * COS(ww*(tt(mi1)-too))
        end do
        Las =  rint(fss,1,30,10,dtt/30)
        return 
    end function Las
    !
    subroutine single_to(too, gaam, nunit12, nunit32)
        integer :: offset, nunit12, nunit32, itt, iloc
        real*8 :: too, gaam
        real*8 :: probJ32(ncon), probJ12(ncon)
        real*8 :: tlo, dt, thi

        tlo = too-5*gaam*SQRT(LOG(2.0d0))
        if(tlo.lt.0) stop 'Use a larger time delay, the initial time is becoming negative'
        thi = too+5*gaam*SQRT(LOG(2.0d0))
        nt = 700
        dt = (thi-tlo)/nt
        ! Reinitialize the vector 
        c = Zer
        do iloc=1, 4
            c(iloc) = A(iloc)*EXP(-Im*Es(iloc)*tlo)
        end do
        write(6,*) 'Time loop started', too
        do itt=1, nt
            call SPf_step(nred, HD*Las(tlo+itt*dt, dt, too, Eo, w, gaam), V, dt, c, caux, Es)
            c = caux
        end do
        probJ12 = 0.0d0
        probJ32 = 0.0d0
        do iloc=1,5
            offset = nos+nbs+(iloc-1)*ncon
            if((iloc.eq.1).or.(iloc.eq.4)) then
                probJ12 = probJ12 + (/(ABS(c(mi2))**2, mi2=offset+1,offset+ncon)/)
            else
                probJ32 = probJ32 + (/(ABS(c(mi2))**2, mi2=offset+1,offset+ncon)/)
            endif
        end do
        do iloc=1,6
            offset = nos+nbs+(iloc-1)*ncon+5*ncon
            if((iloc.eq.3).or.(iloc.eq.6)) then
                probJ12 = probJ12 + (/(ABS(c(mi2))**2, mi2=offset+1,offset+ncon)/)
            else
                probJ32 = probJ32 + (/(ABS(c(mi2))**2, mi2=offset+1,offset+ncon)/)
            endif
        end do
        write(nunit12, 1958) too, probJ12(:)
        write(nunit32, 1958) too, probJ32(:)
        1958 format(201(E20.8E4, 2x))
    end subroutine 
    !
    subroutine example_to(too, Eoo, ww, gaam, nunit)
        integer :: nunit, ntt, iloc
        real*8 :: too, Eoo, ww, gaam
        real*8 :: ttlo, tthi, dtt

        ttlo = too-5*gaam*SQRT(LOG(2.0d0))
        tthi = too+5*gaam*SQRT(LOG(2.0d0))
        dtt = 2*PI/w*1/50
        ntt = int((tthi-ttlo)/dtt)

        c = Zer
        do iloc=1, 4
            c(iloc) = A(iloc)*EXP(-Im*Es(iloc)*ttlo)
        end do
        write(6,*) "Time loop has started"
        do iloc=1, ntt
            call SPf_step(nred, HD*Las(ttlo+iloc*dtt, dtt, too, Eoo, ww, gaam), V, dtt, c, caux, Es)
            c = caux
            write(nunit, 1986) ttlo+iloc*dtt, c(:)
        end do
        call profile(tthi, ttlo, dtt, too, gaam)
        1986 format((E25.8E4,5x,2500(E20.8E4,1x,E20.8E4,3x)))
    end subroutine 
    !
    subroutine example_to_J1(too, Eoo, ww, gaam, nunit)
        integer :: nunit, ntt, iloc
        real*8 :: too, Eoo, ww, gaam
        real*8 :: ttlo, tthi, dtt

        ttlo = too-5*gaam*SQRT(LOG(2.0d0))
        tthi = too+5*gaam*SQRT(LOG(2.0d0))
        dtt = 2*PI/w*1/100
        ntt = int((tthi-ttlo)/dtt)

        c1 = Zer
        do iloc=1, nos1
            c1(iloc) = A(iloc)*EXP(-Im*Es1(iloc)*ttlo)
        end do
        write(6,*) "Time loop has started"
        do iloc=1, ntt
            call SPf_step(nred1, HD1*Las(ttlo+iloc*dtt, dtt, too, Eoo, ww, gaam), V1, dtt, c1, caux1, Es1)
            c1 = caux1
            if(mod(iloc,5).eq.0) write(nunit, 1986) ttlo+iloc*dtt, c1(:)
        end do
        call profile(tthi, ttlo, dtt, too, gaam)
        1986 format((E25.8E4,5x,2500(E20.8E4,1x,E20.8E4,3x)))
    end subroutine
    !
    subroutine exploring_Eo_and_gam()
        do ieo=0, 20
            c = Zer
            call example_to(5992.27734375d0, Eo+ieo*(1.0D-1-Eo)/20, w, gam, 15)
            write(15, 1514) Eo+ieo*(1.0D-1-Eo)/20, c(:)
        end do
        stop
        do igam=0, 20
            c = Zer
            call example_to(5992.27734375d0, Eo, w, gam+igam*gam/20, 15)
            write(15, 1514) Eo+ieo*(1.0D-1-Eo)/20, c(:)
        end do
    1514 format(136(E25.8E4,1x))
    end subroutine 
    !
    subroutine comparing_wJ3_nJ3()
        write(6,*) "Creating the matrix ignoring the continuum"
        ! Get the interaction hamiltonian matrix
        call TPMatNC(evec, RI_OI_J1_J2f, RI_OI_J1_J2p, RI_OI_J1_J0p, RI_OI_J3_J2f, RI_OI_J3_J2p, RI_OI_J3_J4l &
        , ZJ1, ZJ3, ZJ2f, ZJ2p, ZJ0, ZJ4l, H)
        call TPmatNC_J1(evec, RI_OI_J1_J2f, RI_OI_J1_J2p, RI_OI_J1_J0p, ZJ1, ZJ2f, ZJ2p, ZJ0, HJ1)
        write(6,*) "Diagonalizing the large matrix."
        call zdiag(nred, H, Hd, V)
        call zdiag(nred1, HJ1, Hd1, V1)
        ! Initialize the vector 
        c = Zer
        c1 = Zer
        do mi=1, nos1
            c(mi) = A(mi)
            c1(mi) = A(mi)
        end do
        toolo = 5*gam*SQRT(LOG(2.0d0))
        toofo = toolo + 2.0d-12/auoftime
        noo = 100
        dtoo = 2.0d-12/auoftime * 1/noo
        write(6,*) "Evaluating at delay: ", 5992.27734375d0
        write(6,*) "Eo loop"
        open(15, file='SPF_no_cont/wJ3_raman.dat', status='replace')
        open(16, file='SPF_no_cont/nJ3_raman.dat', status='replace')
        call example_to(toolo+20*dtoo, Eo, w, gam, 15)
        call example_to_J1(toolo+20*dtoo, Eo, w, gam,16)
        close(15)
        close(16)
    end subroutine
    !
    subroutine delay_independent()
        call TPmatNC_J1(evec, RI_OI_J1_J2f, RI_OI_J1_J2p, RI_OI_J1_J0p, ZJ1, ZJ2f, ZJ2p, ZJ0, HJ1)
        write(6,*) "Diagonalizing the large matrix."
        call zdiag(nred1, HJ1, Hd1, V1)
        ! Initialize the vector on a single J=1 state.
        A = Zer
        A(1) = Zone
        ! Propagate with zero time delay
        open(64, file='SPF_no_cont/no_delay_c.dat')
        call example_to_J1(0.d0, Eo, w, gam, 64)
        close(64)
        ! Now add a time delay.
        open(65)
        toolo = 5*gam*SQRT(LOG(2.0d0))
        dtoo = 2.0d-12/auoftime * 1/noo
        write(6,*) toolo+20*dtoo
        call example_to_J1(toolo+20*dtoo, Eo, w, gam, 65)
        close(65)
        return
    end subroutine
    !
    subroutine no_delay_J3()
        integer :: ist
        call TPMatNC(evec, RI_OI_J1_J2f, RI_OI_J1_J2p, RI_OI_J1_J0p, RI_OI_J3_J2f, RI_OI_J3_J2p, RI_OI_J3_J4l &
        , ZJ1, ZJ3, ZJ2f, ZJ2p, ZJ0, ZJ4l, H)
        call zdiag(nred, H, Hd, V)
        write(6,*) "Diagonalizing the large matrix."
        call zdiag(nred, H, HD, V)
        ! Now do the things for each one with one amplitude on each
        do ist=1,4
            A = Zer
            A(ist) = Zone
            open(66, file='SPF_no_cont/no_delay_c'//itoa(ist)//'.dat')
            call example_to(0.d0, Eo, w, gam, 66)
            close(66)
        end do
        return
    end subroutine
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