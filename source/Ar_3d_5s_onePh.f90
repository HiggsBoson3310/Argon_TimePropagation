program two_photon
    use ArJ1J2
    use j_recoup
    implicit none
    logical :: exp_flag
    real*8, parameter :: IP32 = 127109.88d0, IP12= IP32+1431.41d0
    !real*8, parameter :: CMperAU = 2.1947463136320e+5
    integer :: nl
    integer :: unit, ii, jj, kk
    integer :: tdelc, ncon
    real*8 :: tdelo, tdelf, dtd
    real*8 :: Ecmi, Elo, dEcont
    real*8, allocatable :: Ecmil(:)
    real*8 :: r1, r2, r3, r4, r5, r6, r7 ,r8
    real*8 :: Eo, w, gam, evecD(3)
    complex*16 :: A(4)
    real*8 :: Ebs(4), Deau
    real*8 :: Zs(4,5)
    real*8, ALLOCATABLE :: regintJ2f(:,:,:,:), iregintJ2f(:,:,:,:)
    real*8, ALLOCATABLE :: regintJ2p(:,:,:,:), iregintJ2p(:,:,:,:)
    real*8, ALLOCATABLE :: regintJ0p(:,:,:,:), iregintJ0p(:,:,:,:)
    real*8 :: dE, Emin
    real*8 :: Ar, Ai
    !integer :: ekinl, line(180), ekin
    !real*8 :: ekinr, frac, frac1, test
    complex*16 :: J2fB(3), J2pB(3), J0pB(2)
    real*8 :: ProbsJ12(200), ProbsJ32(200)

    !--------------------------- CALLING SET UP ROUTINE --------------------------------------
    call setup 
    !-----------------------------------------------------------------------------------------
    
    !---------------------- Call if initial states are not calculated ------------------
    !call InitialStates
    !-----------------------------------------------------------------------------------

    !-------------------------- Define flag for the use of exp energies ----------------
    exp_flag = .true.
    !-----------------------------------------------------------------------------------
    
    !
    open(98, file='init_bounds_J1.dat',status='old')
    do ii=1,4
        read(98,*) Ebs(ii), Zs(ii,:)
    end do
    close(98)

    !Read laser parameters of the laser.
    open(98, file='laser_params_one_photon.dat', status='old')
    read(98,*) Eo, w, gam
    read(98,*) evecD(:)
    read(98, *) nl
    close(98)
    print *, gam, nl, itoa(nl)
    !Intial amplitudes
    !We will divide it into real and imaginary part 
    !We expect them to be normalized and in order of  increasing energy
     if(exp_flag) then 
        open(98, file='init_amplitudes_expE.dat', status='old')
        write(6,*) "The following amplitudes and energies will be used in each state 5s, 3d, 5s'  and 3d'"
        do ii=1,4
            read(98,*) Ar, Ai, Ebs(ii)
            write(6,*) Ar, Ai, Ebs(ii)
            A(ii) = Ar*Zone+Im*Ai
            write(6,*) A(ii)
        enddo
    else
        open(98, file='init_amplitudes.dat', status='old')
        write(6,*) "The following amplitudes will be used in each state 5s, 3d, 5s'  and 3d'"
        do ii=1,4
            read(98,*) Ar, Ai
            write(6,*) Ar, Ai
            A(ii) = Ar*Zone+Im*Ai
            write(6,*) A(ii)
        enddo
    endif
    close(98)
    
    !Initial distribution of states
        !Order of the bound states: 5s, 3d, 5s' and 3d'
        !Amplitudes.
        !Equally excitation in all states.
        !A = (/Zone/SQRT(4.0d0), Zone/SQRT(4.0d0), Zone/SQRT(4.0d0), Zone/SQRT(4.0d0) /)
        !Equally excited 5s and 3d
        !A = (/SQRT(Zone/2.0), SQRT(Zone/2.0d0), Zer, Zer /)
        !Equally excited 5s' and 3d'
        !A = (/Zer, Zer, Zone/SQRT(2.0d0), Zone/SQRT(2.0d0) /)
        !Z coefficients of the bound states
        !read file with information about the Econt grid
    !
    open(08,file='Python_Code_RadInts/EcontEV_det.dat', status='old')
    read(08,*) ncon, Elo, dEcont
    close(08)
    allocate(Ecmil(ncon))
    allocate(regintJ2f(ncon,4, nb, 3), iregintJ2f(ncon,4,nb,3))
    allocate(regintJ2p(ncon,4, nb, 3), iregintJ2p(ncon,4,nb,3))
    allocate(regintJ0p(ncon,4, nb, 2), iregintJ0p(ncon,4,nb,2))
    do jj=0, ncon-1
        Ecmil(jj+1) = Elo+jj*dEcont
    end do
    ! kinetic energy axis and time delay axis:
    dE = ((Ecmil(ncon)-IP32)-(Ecmil(1)-IP12))/180
    Emin = (Ecmil(1)-IP12)
    write(6,*) "Emin: ", Emin
    open(66, file='timde_delays_photoionization.dat', status='replace')
    !open(67, file='testing_the_beats.dat', status='replace')
    tdelc = 1800
    tdelo = 0.00d0
    tdelf = 7.5d4 
    dtd = (tdelf-tdelo)/tdelc
    !Single time delay sanity check list
    open(650, file=""//itoa(nl)//"profile_las.dat")
    open(651, file=""//itoa(nl)//"tdel_probsJ12_oneph.dat")
    open(652, file=""//itoa(nl)//"tdel_probsJ32_oneph.dat")

    do kk=0, tdelc
        ProbsJ12 = 0.0d0
        ProbsJ32 = 0.0d0
        if(kk.eq.1) then
            write(6,*) "File opening"
            open(10, file="Python_Code_RadInts/Rads f/Radial_ef_bd_3d.dat", status="old")
                open(11, file="Python_Code_RadInts/Rads f/Radial_ef_bs_3d.dat", status="old")    
                open(12, file="Python_Code_RadInts/Rads f/Radial_Ir_ef_bd_3d.dat", status="old")    
                open(13, file="Python_Code_RadInts/Rads f/Radial_Ir_ef_bs_3d.dat", status="old")
                

            open(20, file="Python_Code_RadInts/Rads f/Radial_ef_bd_5s.dat", status="old")
                open(21, file="Python_Code_RadInts/Rads f/Radial_ef_bs_5s.dat", status="old")    
                open(22, file="Python_Code_RadInts/Rads f/Radial_Ir_ef_bd_5s.dat", status="old")    
                open(23, file="Python_Code_RadInts/Rads f/Radial_Ir_ef_bs_5s.dat", status="old")
                

            open(30, file="Python_Code_RadInts/Rads f/Radial_ef_bd_3dp.dat", status='old')
                open(31, file="Python_Code_RadInts/Rads f/Radial_ef_bs_3dp.dat", status='old')    
                open(32, file="Python_Code_RadInts/Rads f/Radial_Ir_ef_bd_3dp.dat", status='old')    
                open(33, file="Python_Code_RadInts/Rads f/Radial_Ir_ef_bs_3dp.dat", status='old')
                

            open(40, file="Python_Code_RadInts/Rads f/Radial_ef_bd_5sp.dat", status='old')
                open(41, file="Python_Code_RadInts/Rads f/Radial_ef_bs_5sp.dat", status='old')    
                open(42, file="Python_Code_RadInts/Rads f/Radial_Ir_ef_bd_5sp.dat", status='old')    
                open(43, file="Python_Code_RadInts/Rads f/Radial_Ir_ef_bs_5sp.dat", status='old')
                
            !
            open(50, file="Python_Code_RadInts/Rads p/Radial_ep_bd_3d.dat", status='old')
                open(51, file="Python_Code_RadInts/Rads p/Radial_ep_bs_3d.dat", status='old')    
                open(52, file="Python_Code_RadInts/Rads p/Radial_Ir_ep_bd_3d.dat", status='old')    
                open(53, file="Python_Code_RadInts/Rads p/Radial_Ir_ep_bs_3d.dat", status='old')
                
            open(60, file="Python_Code_RadInts/Rads p/Radial_ep_bd_5s.dat", status='old')
                open(61, file="Python_Code_RadInts/Rads p/Radial_ep_bs_5s.dat", status='old')    
                open(62, file="Python_Code_RadInts/Rads p/Radial_Ir_ep_bd_5s.dat", status='old')    
                open(63, file="Python_Code_RadInts/Rads p/Radial_Ir_ep_bs_5s.dat", status='old')
                
            open(70, file="Python_Code_RadInts/Rads p/Radial_ep_bd_3dp.dat", status='old')
                open(71, file="Python_Code_RadInts/Rads p/Radial_ep_bs_3dp.dat", status='old')    
                open(72, file="Python_Code_RadInts/Rads p/Radial_Ir_ep_bd_3dp.dat", status='old')    
                open(73, file="Python_Code_RadInts/Rads p/Radial_Ir_ep_bs_3dp.dat", status='old')
                
            open(80, file="Python_Code_RadInts/Rads p/Radial_ep_bd_5sp.dat", status='old')
                open(81, file="Python_Code_RadInts/Rads p/Radial_ep_bs_5sp.dat", status='old')    
                open(82, file="Python_Code_RadInts/Rads p/Radial_Ir_ep_bd_5sp.dat", status='old')    
                open(83, file="Python_Code_RadInts/Rads p/Radial_Ir_ep_bs_5sp.dat", status='old')
            
            open(90, file="Python_Code_RadInts/Rads p/Radial_ep_bd_3d.dat", status='old')
                open(91, file="Python_Code_RadInts/Rads p/Radial_ep_bs_3d.dat", status='old')    
                open(92, file="Python_Code_RadInts/Rads p/Radial_Ir_ep_bd_3d.dat", status='old')    
                open(93, file="Python_Code_RadInts/Rads p/Radial_Ir_ep_bs_3d.dat", status='old')
                
            open(100, file="Python_Code_RadInts/Rads p/Radial_ep_bd_5s.dat", status='old')
                open(101, file="Python_Code_RadInts/Rads p/Radial_ep_bs_5s.dat", status='old')    
                open(102, file="Python_Code_RadInts/Rads p/Radial_Ir_ep_bd_5s.dat", status='old')    
                open(103, file="Python_Code_RadInts/Rads p/Radial_Ir_ep_bs_5s.dat", status='old')
                
            open(110, file="Python_Code_RadInts/Rads p/Radial_ep_bd_3dp.dat", status='old')
                open(111, file="Python_Code_RadInts/Rads p/Radial_ep_bs_3dp.dat", status='old')    
                open(112, file="Python_Code_RadInts/Rads p/Radial_Ir_ep_bd_3dp.dat", status='old')    
                open(113, file="Python_Code_RadInts/Rads p/Radial_Ir_ep_bs_3dp.dat", status='old')
                
            open(120, file="Python_Code_RadInts/Rads p/Radial_ep_bd_5sp.dat", status='old')
                open(121, file="Python_Code_RadInts/Rads p/Radial_ep_bs_5sp.dat", status='old')    
                open(122, file="Python_Code_RadInts/Rads p/Radial_Ir_ep_bd_5sp.dat", status='old')    
                open(123, file="Python_Code_RadInts/Rads p/Radial_Ir_ep_bs_5sp.dat", status='old')
                
       
        endif
            !    
        do jj=1,ncon
            Ecmi = Ecmil(jj)
            !Fill the radial integrals for all open channels and bound states
            if(kk.eq.1) then
                do ii=1,4
                    unit = 10*ii
                    call fillradintsJ2(unit, unit+1, regintJ2f(jj,ii,:,:))
                    call fillradintsJ2(unit+2, unit+3, iregintJ2f(jj,ii,:,:))
                    call fillradintsJ2(unit+40, unit+41, regintJ2p(jj,ii,:,:))
                    call fillradintsJ2(unit+42, unit+43, iregintJ2p(jj,ii,:,:))
                    call fillradintsJ0(unit+80, unit+81, regintJ0p(jj,ii,:,:))
                    call fillradintsJ0(unit+82, unit+83, iregintJ0p(jj,ii,:,:))
                end do
            endif
            do ii=1,3
                call cont_Amp(3, Ecmi, ii, 3, 2, 0, A, Zs, Ebs, evecD, J2fB(ii), regintJ2f(jj,:,:,:), iregintJ2f(jj,:,:,:), gam,&
                 Eo, w, tdelo+kk*dtd)
                call cont_Amp(3, Ecmi, ii, 1, 2, 0, A, Zs, Ebs, evecD, J2pB(ii), regintJ2p(jj,:,:,:), iregintJ2p(jj,:,:,:), gam,&
                 Eo, w, tdelo+kk*dtd)
                if(ii.eq.3) then
                    ProbsJ12(jj) = ProbsJ12(jj) + ABS(J2pB(ii))**2 + ABS(J2fB(ii))**2
                else
                    ProbsJ32(jj) = ProbsJ32(jj) + ABS(J2pB(ii))**2 + ABS(J2fB(ii))**2
                endif
            end do

            do ii=1,2
                call cont_Amp(2, Ecmi, ii, 1, 0, 0, A, Zs, Ebs, evecD, J0pB(ii), regintJ0p(jj,:,:,:), iregintJ0p(jj,:,:,:), gam,&
                 Eo, w, tdelo+kk*dtd)
                if(ii.eq.2) then
                    ProbsJ12(jj) = ProbsJ12(jj) + ABS(J0pB(ii))**2
                else
                    ProbsJ32(jj) = ProbsJ32(jj) + ABS(J0pB(ii))**2
                endif
            end do
            if(kk.eq.int(1200/dtd)) then
                write(650,1985) Ecmi, ABS(FourF(gam, Deau, 0.d0))**2
            end if
        end do
        if(kk.eq.1) then
            write(6,*) "File closing"
            do ii=1,12
                unit = 10*ii
                write(6,*) "Closing", unit
                close(unit)
                close(unit+1)
                close(unit+2)
                close(unit+3)
            enddo
        endif
        !write(6,*) kk, ncon
        write(651,6831) tdelo+kk*dtd, ProbsJ12(:)
        write(652,6831) tdelo+kk*dtd, ProbsJ32(:)
    end do
    close(650)
    close(651)
    close(652)
1985 format(2(e18.10E3,2x))
!1986 format(4(ES18.10E3,2x))
6831 format(201(ES18.10E3,2x))
    
    contains 
    subroutine fillradintsJ2(unit1, unit2, matrix)
        integer :: unit1, unit2
        real*8 :: matrix(5,3)
        !Order of this files is:
        !Ecmi  rad_ef(12)_nu(12)  rad_ef(32)_nu(12)  rad_ef(12)_nu(32)  rad_ef(32)_nu(32)
        read(unit1,*) Ecmi, r1, r2, r3, r4 !d outer electron for J=1
        read(unit2,*) Ecmi, r5, r6, r7, r8 !s outer electron for J=1
        matrix(1,1) = r2
        matrix(1,2) = r2
        matrix(1,3) = r1
        matrix(2,1) = r4
        matrix(2,2) = r4
        matrix(2,3) = r3
        matrix(3,1) = r4
        matrix(3,2) = r4
        matrix(3,3) = r3
        matrix(4,1) = r6
        matrix(4,2) = r6
        matrix(4,3) = r5
        matrix(5,1) = r8
        matrix(5,2) = r8
        matrix(5,3) = r7
    end subroutine
    subroutine fillradintsJ0(unit1, unit2, matrix)
        integer :: unit1, unit2
        real*8 :: matrix(5,2)
        !Order of this files is:
        !Ecmi  rad_ep(12)_nu(12)  rad_ep32)_nu(12)  rad_ep(12)_nu(32)  rad_ep(32)_nu(32)
        read(unit1,*) Ecmi, r1, r2, r3, r4 !d outer electron for J=1
        read(unit2,*) Ecmi, r5, r6, r7, r8 !s outer electron for J=1
        matrix(1,1) = r2
        matrix(1,2) = r1
        matrix(2,1) = r4
        matrix(2,2) = r3
        matrix(3,1) = r4
        matrix(3,2) = r3
        matrix(4,1) = r6
        matrix(4,2) = r5
        matrix(5,1) = r8
        matrix(5,2) = r7
    end subroutine

    function close_int(aa) result(y)
        real*8 :: aa
        integer :: y, yl, yh
        yl = int(aa)
        yh = yl+1
        if(abs(yl-aa).le.abs(yh-aa)) then
            y = yl
        else 
            y = yh
        endif
    end function

    function itoa(iii) result(res)
        character(:),allocatable :: res
        integer,intent(in) :: iii
        character(range(iii)+2) :: tmp
        write(tmp,'(i0)') iii
        res = trim(tmp)
    end function
end program two_photon




