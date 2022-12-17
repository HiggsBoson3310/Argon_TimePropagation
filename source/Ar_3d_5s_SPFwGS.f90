program main
    use j_recoup
    use arj1j2
    use ERRFUN
    use spf_prop
    use matmod
    use omp_lib
    implicit none
    integer, parameter :: nconts=300
    logical :: exp_flag
    integer :: ncon, mm, kk
    integer :: nl, nlf
    integer :: mi, mi1, ii1
    real*8 :: Elo, dEcont
    real*8 :: Z1(nos1, 5), EJ1(nos1)
    real*8 :: Z3(nos3,6), EJ3(nos3)
    real*8 :: ZJ0(nbs0,2), EnJ0(nbs0)
    real*8 :: ZJ2f(nbs2f,3), EnJ2f(nbs2f)
    real*8 :: ZJ2p(nbs2p,3), EnJ2p(nbs2p)
    real*8 :: ZJ0au(nau,2), ZJ2fau(nau,2), ZJ2pau(nau,2)
    real*8 :: ZJ4l(nbs4,6), EnJ4l(nbs4)
    ! Bound state core selector 0 jc=1/2, 1 jc=3/2
    integer :: C_J1(nos1), C_J3(nos3), C_J0(nbs0), C_J2f(nbs2f), C_J2p(nbs2p), C_J4(nbs4)
    ! Continuum radial elements
    real*8 :: rRI_CI_J0p_J1(nconts,nbs0,5,2)
    real*8 :: iRI_CI_J0p_J1(nconts,nbs0,5,2)
    real*8 :: rRI_CI_J2f_J1(nconts,nbs2f,5,3), rRI_CI_J2f_J3(nconts,nbs2f,6,3)
    real*8 :: iRI_CI_J2f_J1(nconts,nbs2f,5,3), iRI_CI_J2f_J3(nconts,nbs2f,6,3)
    real*8 :: rRI_CI_J2p_J1(nconts,nbs2p,5,3), rRI_CI_J2p_J3(nconts,nbs2p,6,3)
    real*8 :: iRI_CI_J2p_J1(nconts,nbs2p,5,3), iRI_CI_J2p_J3(nconts,nbs2p,6,3)
    real*8 :: rRI_CI_J4l_J3(nconts,nbs4,6,6), rRI_CI_J4l_J5(nconts,nbs4,6,6)
    real*8 :: iRI_CI_J4l_J3(nconts,nbs4,6,6), iRI_CI_J4l_J5(nconts,nbs4,6,6)
    ! Bound to autoionizing states
    real*8 :: rRI_CI_J1_J0p(nau,nos1,2,5), iRI_CI_J1_J0p(nau,nos1,2,5)
    real*8 :: rRI_CI_J1_J2p(nau,nos1,3,5), iRI_CI_J1_J2p(nau,nos1,3,5)
    real*8 :: rRI_CI_J1_J2f(nau,nos1,3,5), iRI_CI_J1_J2f(nau,nos1,3,5)
    real*8 :: rRI_CI_J3_J2p(nau,nos3,3,6), iRI_CI_J3_J2p(nau,nos3,3,6)
    real*8 :: rRI_CI_J3_J2f(nau,nos3,3,6), iRI_CI_J3_J2f(nau,nos3,3,6)
    ! bouns to bound radial integrals
    real*8 :: RI_OI_J1_J2f(nos1, nbs2f, 3, 5), RI_OI_J1_J2p(nos1, nbs2p, 3,5)
    real*8 :: RI_OI_J1_J0p(nos1, nbs0,2, 5)
    real*8 :: RI_OI_J3_J2f(nos3, nbs2f, 3, 6), RI_OI_J3_J2p(nos3, nbs2p, 3, 6)
    real*8 :: RI_OI_J3_J4l(nos3, nbs4, 6, 6)
    real*8 :: Eo, w, gam, evec(3)
    real*8 :: Epa, wopa, gamopa, tdelay
    real*8 :: Ar, Ai
    complex*16 :: A(nos)
    complex*16 :: H(ntotau+1, ntotau+1), V(ntotau+1,ntotau+1),&
     c(ntotau+1), caux(ntotau+1), csf(nosgs, ntotau+1)
    real*8 :: Es(ntotau+1), HD(ntotau+1)
    !integer :: nt
    real :: finish, start
    character(len=143) :: file_root
    character(len=10) :: nts
    !real*8 :: tttloo, ttthoi, dttt
    !integer :: nttt

    ! Commands for setting up what variable to loop over:

    integer :: which_loop
    integer :: name_displace, loop_end
    real*8 :: lower_lim, differential, scale_low, scale_differential
    logical :: reference_exist


    call cpu_time(start)
    !-------------------------------- calling set up from gensub -----------------------
    call SETUP
    !-----------------------------------------------------------------------------------
    
    !---------------------- Call if initial states are not calculated ------------------
    !call InitialStates
    !-----------------------------------------------------------------------------------
    
    !-------------------------- Define flag for the use of exp energies ----------------
    exp_flag = .true.
    !-----------------------------------------------------------------------------------

    !-------------------------- Get the folder with the laser inputs -------------------
    call get_command_argument(1, nts)
    read(nts,'(I10)') nl
    !-------------------------- Indicators for looping over param ----------------------
    call get_command_argument(2,nts)
    read(nts,'(I10)') which_loop
    call get_command_argument(3,nts)
    read(nts,'(I10)') loop_end
    call get_command_argument(4, nts)
    read(nts,'(f4.6)') lower_lim
    call get_command_argument(5, nts)
    read(nts,'(f4.6)') differential
    call get_command_argument(6,nts)
    read(nts,'(I10)') name_displace

    !-----------------------------------------------------------------------------------

    ! Retrieving from precomputed files.
    
    ! Loading the Zcoefficients
    call load_Zs()
    !

    ! Loading the radial elements
    call load_rads()
    !

    ! Laser parameters
    open(98, file='./SPF_DRIVER_'//itoa(nl)//'/laser_params_driver.dat', status='old')
    read(98,*) Eo, w, gam
    read(98,*) evec(:)
    read(98,*) nlf
    close(98)
    print *, "For the driver/harmonics", gam, Eo
    open(98, file='./SPF_DRIVER_'//itoa(nl)//'/laser_params_opa.dat', status='old')
    read(98,*) Epa, wopa, gamopa
    read(98,*) tdelay
    close(98)
    print *, "For the OPA", gamopa, Epa
    !
    if(.not.(nlf.eq.nlf)) stop "Incorrect file-folder combination."
    
    ! Continuum paramerers
    open(98, file='Python_Code_RadInts/EcontEV_det.dat',status='old')
    read(98,*) ncon, Elo, dEcont
    close(98)
    write(6,*) "Found ", ncon, " continuum states, and density ", dEcont, "cm(-1)"
    if(ncon.ne.nconts) then
        write(6,*) "Not the same number of continuum states in the file as in the program"
        STOP
    endif

    ! Initial state amplitudes
    call A_init()

    write(6,*) "The following amplitudes and energies will be used in each state 5s, 3d, 5s'  and 3d'"
    write(6,*) "if the option to not start everything in the grounds state is chosen"
    do ii1=1, nos1
        write(6,*) REAL(A(ii1)), AIMAG(A(ii1)), EJ1(ii1)
    end do

    do ii1=1, nos3
        write(6,*) REAL(A(ii1+nos1)), AIMAG(A(ii1)+nos1), EJ3(ii1)
    end do

    ! Put all energies in a single array
    call all_Energiesau()

    ! Initialize and diagonalize the dipole op.

    !
    !call set_and_diagHau()
    call set_and_diagHau
    !

    !------ Determine if a reference calculation is needed -----------
    inquire(file='SPF_DRIVER_'//itoa(nl)//'/reference_spec.dat',EXIST=reference_exist)
    if(.not.(reference_exist)) then
        write(6,*) 'There was no reference file calculated.' 
        write(6,*) 'Running one to have the driver information.'
        open(18,file='./SPF_DRIVER_'//itoa(nl)//'/reference_spec.dat',status='replace')
        call OPA_time_delay_phase(Eo,gam, w, Epa,gamopa,wopa,18)
        close(18)
    endif
    ! ----- Determine which loops are to be made ---------------------
    if(which_loop.eq.0) then
        scale_low = (lower_lim/eVperAU)
        scale_differential = differential/eVperAU
        call Spec_freqs_loop(Eo,gam,w,Epa,gamopa,18, &
        loop_end,scale_low,scale_differential,name_displace)
    else if(which_loop.eq.1) then
        scale_low = lower_lim*1.0d12/auofI
        scale_differential = differential*1.0d12/auofI
        call Spec_Intensity_loop(Eo,gam,w,gamopa,wopa,18,&
        loop_end,scale_low,scale_differential,name_displace)
    else if(which_loop.eq.2) then
        scale_low = lower_lim*1.0d-15/auoftime*1/SQRT(2*LOG(2.d0))
        scale_differential = differential*1.0d-15/auoftime*1/SQRT(2*LOG(2.d0))
        call Spec_duration_loop(Eo,gam,w,Epa,wopa,18,&
        loop_end,scale_low,scale_differential,name_displace)
    else
        STOP 'Invalid loop requested. Options are only 0 for frequency, 1 for intensity and 2 for duration'
    endif

    
    !call dur_in_table(Epa,gamopa, wopa)

    call cpu_time(finish)
    write(6,*) "Total calculation took ", finish-start, " seconds."
    stop
!    1986 format((E25.8E4,5x,50000(E20.8E4,1x,E20.8E4,3x)))
!    1988 format((E25.8E4,5x,E20.8E4,1x,E20.8E4))
    contains
    !
    subroutine Load_Zs()
        integer :: nnmax
        open(10, file='init_bounds_J1.dat')
        do mm=1,nos1
            read(10,*) EJ1(mm), (Z1(mm,kk), kk=1,5)
            nnmax = maxloc(Z1(mm,:),1)
            if((nnmax.eq.2).or.(nnmax.eq.3).or.(nnmax.eq.5)) then
                C_J1(mm) = 1
            else
                C_J1(mm) = 0
            endif
            !write(6,*) Ebs(mm), (Zo(mm,kk), kk=1,5)
        end do
        close(10)

        open(10, file='init_bounds_J3.dat')
        do mm=1,nos3
            read(10,*) EJ3(mm), (Z3(mm,kk), kk=1,6)
            nnmax = maxloc(Z3(mm,:),1)
            if((nnmax.eq.1).or.(nnmax.eq.2).or.(nnmax.eq.4).or.(nnmax.eq.5)) then
                C_J3(mm) = 1
            else
                C_J1(mm) = 0
            endif
            !write(6,*) Ebs(mm), (Zo(mm,kk), kk=1,5)
        end do
        close(10)

        open(10, file='interm_bounds_J2f.dat')
        do mm=1,nbs2f
            read(10,*) EnJ2f(mm), (ZJ2f(mm,kk), kk=1,3)
            nnmax = maxloc(ZJ2f(mm,:),1)
            if((nnmax.eq.1).or.(nnmax.eq.2)) then
                C_J2f(mm) = 1
            else
                C_J2f(mm) = 0
            endif
            !write(6,*) EnJ2f(mm), (ZJ2f(mm,kk), kk=1,3)
        end do
        close(10)

        open(10, file='test_of_zs_J2f.dat')
        do mm=1,nau
            read(10,*) (ZJ2fau(mm,kk), kk=1,2)
        end do
        close(10)

        open(10, file='interm_bounds_J2p.dat')
        do mm=1,nbs2p
            read(10,*) EnJ2p(mm), (ZJ2p(mm,kk), kk=1,3)
            nnmax = maxloc(ZJ2p(mm,:),1)
            if((nnmax.eq.1).or.(nnmax.eq.2)) then
                C_J2p(mm) = 1
            else
                C_J2p(mm) = 0
            endif
            !write(6,*) EnJ2p(mm), (ZJ2p(mm,kk), kk=1,3)
        end do
        close(10)

        open(10, file='test_of_zs_J2p.dat')
        do mm=1,nau
            read(10,*) (ZJ2pau(mm,kk), kk=1,2)
        end do
        close(10)

        open(10, file='interm_bounds_J0p.dat')
        do mm=1,nbs0
            read(10,*) EnJ0(mm), (ZJ0(mm,kk), kk=1,2)
            nnmax = maxloc(ZJ0(mm,:),1)
            if((nnmax.eq.1)) then
                C_J0(mm) = 1
            else
                C_J0(mm) = 0
            endif
            !write(6,*) EnJ0(mm), (ZJ0(mm,kk), kk=1,2)
        end do
        close(10)

        open(10, file='test_of_zs_J0.dat')
        do mm=1,nau
            read(10,*) (ZJ0au(mm,kk), kk=1,2)
        end do
        close(10)

        open(10, file='interm_bounds_J4l.dat')
        do mm=1,nbs4
            read(10,*) EnJ4l(mm), (ZJ4l(mm,kk), kk=1,6)
            nnmax = maxloc(ZJ4l(mm,:),1)
            if((nnmax.eq.1).or.(nnmax.eq.2).or.(nnmax.eq.4).or.(nnmax.eq.5)) then
                C_J4(mm) = 1
            else
                C_J4(mm) = 0
            endif
            !write(6,*) EnJ4l(mm), (ZJ4l(mm,kk), kk=1,2)
        end do
        close(10)
    end subroutine
    !
    subroutine load_rads()
        write(6,*) "Loading Radial Elements into four dimensional arrays"

        ! Radial elements from initial J=1 state, to intermediate J=2 state
        ! f states
        file_root = 'Python_Code_RadInts/RI b(J=1)  b(J=2 lint f)/interm_rad_ints_'
        call Init_big_matsW(nos1, nbs2f, RI_OI_J1_J2f, 3, 5, file_root)
        ! p states
        file_root = 'Python_Code_RadInts/RI b(J=1)  b(J=2 lint p)/interm_rad_ints_'
        call Init_big_matsW(nos1, nbs2p, RI_OI_J1_J2p, 3, 5, file_root)
        !J=0
        file_root = 'Python_Code_RadInts/RI b(J=1)  b(J=0 lint p)/interm_rad_ints_'
        call Init_big_matsW(nos1, nbs0, RI_OI_J1_J0p, 2, 5, file_root)
        ! And to the autoionizing states
        ! Init_big_mats(nnbs, nnconts, Rs, chanc, chanb, file_name)
        file_root = 'Python_Code_RadInts/RI b(J=1) c(J=0,l=p)/cont_rad_ints_regb'
        call Init_big_mats(nos1, nau, rRI_CI_J1_J0p, 2, 5, file_root)
        ! autoionizing
        file_root = 'Python_Code_RadInts/RI b(J=1) c(J=0,l=p)/cont_rad_ints_iregb'
        call Init_big_mats(nos1, nau, iRI_CI_J1_J0p, 2, 5, file_root)
        
        file_root = 'Python_Code_RadInts/RI b(J=1) c(J=2,l=p)/cont_rad_ints_regb'
        call Init_big_mats(nos1, nau, rRI_CI_J1_J2p, 3, 5, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=1) c(J=2,l=p)/cont_rad_ints_iregb'
        call Init_big_mats(nos1, nau, iRI_CI_J1_J2p, 3, 5, file_root)

        file_root = 'Python_Code_RadInts/RI b(J=1) c(J=2,l=f)/cont_rad_ints_regb'
        call Init_big_mats(nos1, nau, rRI_CI_J1_J2f, 3, 5, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=1) c(J=2,l=f)/cont_rad_ints_iregb'
        call Init_big_mats(nos1, nau, iRI_CI_J1_J2f, 3, 5, file_root)

        ! Radial elements from initial J=3 state, to intermediate J=2 state
        ! f states
        file_root = 'Python_Code_RadInts/RI b(J=3)  b(J=2 lint f)/interm_rad_ints_'
        call Init_big_matsW(nos3, nbs2f, RI_OI_J3_J2f, 3, 6, file_root)
        ! p states
        file_root = 'Python_Code_RadInts/RI b(J=3)  b(J=2 lint p)/interm_rad_ints_'
        call Init_big_matsW(nos3, nbs2p, RI_OI_J3_J2p, 3, 6, file_root)
        !J=4
        file_root = 'Python_Code_RadInts/RI b(J=3)  b(J=4 lint f)/interm_rad_ints_'
        call Init_big_matsW(nos3, nbs4, RI_OI_J3_J4l, 6, 6, file_root)
        !  Autoionizing states
        file_root = 'Python_Code_RadInts/RI b(J=3) c(J=2,l=p)/cont_rad_ints_reg_b'
        call Init_big_mats(nos3, nau, rRI_CI_J3_J2p, 3, 6, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=3) c(J=2,l=p)/cont_rad_ints_ireg_b'
        call Init_big_mats(nos3, nau, iRI_CI_J3_J2p, 3, 6, file_root)

        file_root = 'Python_Code_RadInts/RI b(J=3) c(J=2,l=f)/cont_rad_ints_reg_b'
        call Init_big_mats(nos3, nau, rRI_CI_J3_J2f, 3, 6, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=3) c(J=2,l=f)/cont_rad_ints_ireg_b'
        call Init_big_mats(nos3, nau, iRI_CI_J3_J2f, 3, 6, file_root)
        ! Radial elements to J=1 continuum, from an intermediate J=2 state
        ! f states
        file_root = 'Python_Code_RadInts/RI b(J=2,l=f) c(J=1)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs2f, nconts, rRI_CI_J2f_J1, 5, 3, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=2,l=f) c(J=1)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs2f, nconts, iRI_CI_J2f_J1, 5, 3, file_root)
        ! p states
        file_root = 'Python_Code_RadInts/RI b(J=2,l=p) c(J=1)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs2p, nconts, rRI_CI_J2p_J1, 5, 3, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=2,l=p) c(J=1)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs2p, nconts, iRI_CI_J2p_J1, 5, 3, file_root)
        ! Radial elements to J=1 continuum, from an intermediate J=0
        file_root = 'Python_Code_RadInts/RI b(J=0,l=p) c(J=1)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs0, nconts, rRI_CI_J0p_J1, 5, 2, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=0,l=p) c(J=1)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs0, nconts, iRI_CI_J0p_J1, 5, 2, file_root)

        ! Radial elements to J=3 continuum, from an intermediate J=2 state
        ! f states
        file_root = 'Python_Code_RadInts/RI b(J=2,l=f) c(J=3)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs2f, nconts, rRI_CI_J2f_J3, 6, 3, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=2,l=f) c(J=3)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs2f, nconts, iRI_CI_J2f_J3, 6, 3, file_root)
        ! p states
        !J=2
        file_root = 'Python_Code_RadInts/RI b(J=2,l=p) c(J=3)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs2p, nconts, rRI_CI_J2p_J3, 6, 3, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=2,l=p) c(J=3)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs2p, nconts, iRI_CI_J2p_J3, 6, 3, file_root)
        ! Radial elements to J=3 continuum, from an intermediate J=4 state
        file_root = 'Python_Code_RadInts/RI b(J=4,l=f) c(J=3)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs4, nconts, rRI_CI_J4l_J3, 6, 6, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=4,l=f) c(J=3)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs4, nconts, iRI_CI_J4l_J3, 6, 6, file_root)

        ! Radial elements to J=5 continuum, from an intermediate J=4 state
        file_root = 'Python_Code_RadInts/RI b(J=4,l=f) c(J=5)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs4, nconts, rRI_CI_J4l_J5, 6, 6, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=4,l=f) c(J=5)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs4, nconts, iRI_CI_J4l_J5, 6, 6, file_root)
    end subroutine
    !
    subroutine A_init()
        if(exp_flag) then 
            open(98, file='./SPF_DRIVER_'//itoa(nl)//'/init_amplitudes_expE.dat', status='old')
            do ii1=1,nos1
                read(98,*) Ar, Ai, EJ1(ii1)
                A(ii1) = Ar*Zone+Im*Ai
            enddo
            do ii1=1,nos3
                read(98,*) Ar, Ai, EJ3(ii1)
                A(ii1+nos1) = Ar*Zone+Im*Ai
            enddo
            close(98)
        else
            open(98, file='./SPF_DRIVER_'//itoa(nl)//'/init_amplitudes.dat', status='old')
            do ii1=1,nos1
                read(98,*) Ar, Ai
                A(ii1) = Ar*Zone+Im*Ai
            enddo
            do ii1=1,nos3
                read(98,*) Ar, Ai
                A(ii1+nos1) = Ar*Zone+Im*Ai
            enddo
            close(98)
        endif
    end subroutine A_init
    !
    subroutine all_Energies()
        ! Vector with all he energies in atomic units
        Es(1) = 0.d0
        Es(2:nos1+1) = EJ1/CMperAU
        Es(nos1+2:nosgs) = EJ3/CMperAU
        Es(nosgs+1:nbs2f+nosgs) = EnJ2f/CMperAU
        Es(nbs2f+nosgs+1: nbs2f+nbs2p+nosgs) = EnJ2p/CMperAU
        Es(nbs2p+nbs2f+nosgs+1:nosgs+nbs2f+nbs2p+nbs0) = EnJ0/CMperAU
        Es(nosgs+nbs2f+nbs2p+nbs0+1:nosgs+nbs) = Enj4l/CMperAU
        do mi1=1,5
        Es(nbs+nosgs+(mi1-1)*ncon+1:nbs+nosgs+(mi1-1)*ncon+ncon) = (/ ((Elo+mi*dEcont)/CMperAU, mi=1, ncon) /)
        end do
        do mi1=1,6
        Es(nbs+nosgs+(mi1-1)*ncon+5*ncon+1:nbs+nosgs+(mi1-1)*ncon+ncon+5*ncon) = (/ ((Elo+mi*dEcont)/CMperAU, mi=1, ncon) /)
        end do
        do mi1=1,6
        Es(nbs+nosgs+(mi1-1)*ncon+5*ncon+6*ncon+1:nbs+nosgs+(mi1-1)*ncon+ncon+5*ncon+6*ncon) &
        = (/ ((Elo+mi*dEcont)/CMperAU, mi=1, ncon) /)
        end do
        !
        open(41, file='SPF_DRIVER_'//itoa(nl)//'/energies.dat')
        write(41,1981) Es(:)
        close(41)
        1981 format(50000(E20.7E4,1x))
    end subroutine all_Energies
    !
    subroutine all_Energiesau()
        ! Vector with all the energies in atomic units
        Es(1) = 0.d0
        Es(2:nos1+1) = EJ1/CMperAU
        Es(nos1+2:nosgs) = EJ3/CMperAU
        Es(nosgs+1:nbs2f+nosgs) = EnJ2f/CMperAU
        Es(nbs2f+nosgs+1: nbs2f+nbs2p+nosgs) = EnJ2p/CMperAU
        Es(nbs2f+nbs2p+nosgs+1:nosgs+nbs2f+nbs2p+nbs0) = EnJ0/CMperAU
        Es(nosgs+nbs2f+nbs2p+nbs0+1:nosgs+nbs) = Enj4l/CMperAU
        do mi1=1,2
            Es(nbs+nosgs+(mi1-1)*nau+1:nbs+nosgs+(mi1-1)*nau+nau) = (/ (ZJ2fau(mi,1)/CMperAU, mi=1, nau) /)
        end do
        do mi1=1,2
            Es(nbs+nosgs+(mi1-1)*nau+1+2*nau:nbs+nosgs+(mi1-1)*nau+nau+2*nau) = (/ (ZJ2pau(mi,1)/CMperAU, mi=1, nau) /)
        end do
        do mi1=1,1
            Es(nbs+nosgs+(mi1-1)*nau+1+4*nau:nbs+nosgs+(mi1-1)*nau+nau+4*nau) = (/ (ZJ0au(mi,1)/CMperAU, mi=1, nau) /)
        end do
    
        do mi1=1,5
        Es(nbs+nosgs+5*nau+(mi1-1)*ncon+1:nbs+nosgs+5*nau+(mi1-1)*ncon+ncon) = (/ ((Elo+mi*dEcont)/CMperAU, mi=1, ncon) /)
        end do
        do mi1=1,6
        Es(nbs+nosgs+5*nau+(mi1-1)*ncon+5*ncon+1:nbs+nosgs+5*nau+(mi1-1)*ncon+ncon+5*ncon) &
        = (/ ((Elo+mi*dEcont)/CMperAU, mi=1, ncon) /)
        end do
        do mi1=1,6
        Es(nbs+nosgs+5*nau+(mi1-1)*ncon+5*ncon+6*ncon+1:nbs+nosgs+5*nau+(mi1-1)*ncon+ncon+5*ncon+6*ncon) &
        = (/ ((Elo+mi*dEcont)/CMperAU, mi=1, ncon) /)
        end do
        !
        open(41, file='SPF_DRIVER_'//itoa(nl)//'/energies.dat')
        write(41,1981) Es(:)
        close(41)
        1981 format(50000(E20.7E4,1x))
    end subroutine all_Energiesau
    !
    subroutine set_and_diagH()
        real :: tst, tnd
        call TPMatEDwGS(evec, Elo, dEcont,&
        RI_OI_J1_J0p, RI_OI_J1_J2f, RI_OI_J1_J2p, RI_OI_J3_J2f,RI_OI_J3_J2p, RI_OI_J3_J4l,&
        Z1, Z3, ZJ2f, ZJ2p, ZJ0, ZJ4l, &
        rRI_CI_J0p_J1, iRI_CI_J0p_J1, rRI_CI_J2f_J1, iRI_CI_J2f_J1, &
        rRI_CI_J2p_J1, iRI_CI_J2p_J1,&
        rRI_CI_J2p_J3, iRI_CI_J2p_J3, rRI_CI_J2f_J3, iRI_CI_J2f_J3, rRI_CI_J4l_J3, iRI_CI_J4l_J3, &
        rRI_CI_J4l_J5, iRI_CI_J4l_J5, H)
        call cpu_time(tst)
        call zdiag(ntot+1, H, HD, V)
        call cpu_time(tnd)
        write(6,*) "The diagionalization took", tnd-tst," seconds"
        
    end subroutine set_and_diagH
    !
    subroutine set_and_diagHau()
        real :: tst, tnd
        call TPMatEDwAu(evec, Elo, dEcont,&
        RI_OI_J1_J0p, RI_OI_J1_J2f, RI_OI_J1_J2p, RI_OI_J3_J2f,RI_OI_J3_J2p, RI_OI_J3_J4l,&
        Z1, Z3, ZJ2f, ZJ2p, ZJ0, ZJ4l, &
        ZJ2fau,ZJ2pau,ZJ0au,&
        rRI_CI_J1_J2f, iRI_CI_J1_J2f, rRI_CI_J1_J2p, iRI_CI_J1_J2p, rRI_CI_J1_J0p, iRI_CI_J1_J0p,&
        rRI_CI_J3_J2f, iRI_CI_J3_J2f, rRI_CI_J3_J2p, iRI_CI_J3_J2p,&
        rRI_CI_J0p_J1, iRI_CI_J0p_J1, rRI_CI_J2f_J1, iRI_CI_J2f_J1, &
        rRI_CI_J2p_J1, iRI_CI_J2p_J1,&
        rRI_CI_J2p_J3, iRI_CI_J2p_J3, rRI_CI_J2f_J3, iRI_CI_J2f_J3, rRI_CI_J4l_J3, iRI_CI_J4l_J3, &
        rRI_CI_J4l_J5, iRI_CI_J4l_J5, H)
        call cpu_time(tst)
        call zdiag_save(ntotau+1, H, HD, V)
        call cpu_time(tnd)
        write(6,*) "The diagionalization took", tnd-tst," seconds"
    end subroutine set_and_diagHau
    !
    subroutine set_and_diagHau_seg()
        integer :: ii, jj, Jist,Jjst
        real :: tst, tnd
        call TPMatEDwAu(evec, Elo, dEcont,&
        RI_OI_J1_J0p, RI_OI_J1_J2f, RI_OI_J1_J2p, RI_OI_J3_J2f,RI_OI_J3_J2p, RI_OI_J3_J4l,&
        Z1, Z3, ZJ2f, ZJ2p, ZJ0, ZJ4l, &
        ZJ2fau,ZJ2pau,ZJ0au,&
        rRI_CI_J1_J2f, iRI_CI_J1_J2f, rRI_CI_J1_J2p, iRI_CI_J1_J2p, rRI_CI_J1_J0p, iRI_CI_J1_J0p,&
        rRI_CI_J3_J2f, iRI_CI_J3_J2f, rRI_CI_J3_J2p, iRI_CI_J3_J2p,&
        rRI_CI_J0p_J1, iRI_CI_J0p_J1, rRI_CI_J2f_J1, iRI_CI_J2f_J1, &
        rRI_CI_J2p_J1, iRI_CI_J2p_J1,&
        rRI_CI_J2p_J3, iRI_CI_J2p_J3, rRI_CI_J2f_J3, iRI_CI_J2f_J3, rRI_CI_J4l_J3, iRI_CI_J4l_J3, &
        rRI_CI_J4l_J5, iRI_CI_J4l_J5, H)
        Jist = 1
        do ii=1, nos1
            Jjst = nosgs
            do jj=1, nbs2f
            if(C_J1(ii).ne.C_J2f(jj)) then
                H(Jist+ii,Jjst+jj) = Zer
                H(Jjst+jj,Jist+ii) = Zer
            endif
            end do
            Jjst = nosgs+nbs2f
            do jj=1,nbs2p
            if(C_J1(ii).ne.C_J2p(jj)) then
                H(Jist+ii,Jjst+jj) = Zer
                H(Jjst+jj,Jist+ii) = Zer
            endif
            end do
            Jjst = nosgs+nbs2f+nbs2p
            do jj=1,nbs0
            if(C_J1(ii).ne.C_J0(jj)) then
                H(Jist+ii,Jjst+jj) = Zer
                H(Jjst+jj,Jist+ii) = Zer
            endif
            end do
        end do
        Jist = nos1+1
        do ii=1, nos3
            Jjst = nosgs
            do jj=1, nbs2f
            if(C_J3(ii).ne.C_J2f(jj)) then
                H(Jist+ii,Jjst+jj) = Zer
                H(Jjst+jj,Jist+ii) = Zer
            endif
            end do
            Jjst = nosgs+nbs2f
            do jj=1,nbs2p
            if(C_J3(ii).ne.C_J2p(jj)) then
                H(Jist+ii,Jjst+jj) = Zer
                H(Jjst+jj,Jist+ii) = Zer
            endif
            end do
            Jjst = nosgs+nbs2f+nbs2p+nbs0
            do jj=1,nbs4
            if(C_J3(ii).ne.C_J4(jj)) then
                H(Jist+ii,Jjst+jj) = Zer
                H(Jjst+jj,Jist+ii) = Zer
            endif
            end do
        end do

        open(712,file='H_segMatR.db',STATUS='REPLACE', &
        ACCESS='STREAM',FORM='UNFORMATTED')
        open(812,file='H_segMatI.db',STATUS='REPLACE', &
        ACCESS='STREAM',FORM='UNFORMATTED')

        write(712) REAL(H)
        write(812) AIMAG(H)

        CLOSE(712)
        CLOSE(812)

        call cpu_time(tst)
        call zdiag(ntotau+1, H, HD, V)
        

        call cpu_time(tnd)
        write(6,*) "The diagionalization took", tnd-tst," seconds"
    end subroutine set_and_diagHau_seg
    !
    subroutine set_and_diagHau_from_file()
        real*8 :: aux(ntotau+1,ntotau+1)
        open(412,file='DiagH.db',STATUS='OLD', &
        ACCESS='STREAM',FORM='UNFORMATTED')
        open(512,file='RUmat.db',STATUS='OLD', &
        ACCESS='STREAM',FORM='UNFORMATTED')
        open(612,file='IUmat.db',STATUS='OLD', &
        ACCESS='STREAM',FORM='UNFORMATTED')
        open(712,file='HMatR.db',STATUS='OLD', &
        ACCESS='STREAM',FORM='UNFORMATTED')
        open(812,file='HMatI.db',STATUS='OLD', &
        ACCESS='STREAM',FORM='UNFORMATTED')
        !
        read(412) HD
        !
        read(512) aux
        V = reshape(aux,(/ntotau+1,ntotau+1/))*Zone
        read(612) aux
        V = reshape(aux,(/ntotau+1,ntotau+1/))*Im + V
        !
        read(712) aux
        H = reshape(aux,(/ntotau+1,ntotau+1/))*Zone
        read(812) aux
        H = reshape(aux,(/ntotau+1,ntotau+1/))*Im + H
        !
        close(412)
        close(512)
        close(612)
        CLOSE(712)
        CLOSE(812)
    end subroutine set_and_diagHau_from_file
    !
    subroutine set_and_diagHnaungs()
        integer :: ix, iy
        real*8 :: sparse, sparseH
        real :: tst, tnd
        call TPMatEDwAu(evec, Elo, dEcont,&
        RI_OI_J1_J0p, RI_OI_J1_J2f, RI_OI_J1_J2p, RI_OI_J3_J2f,RI_OI_J3_J2p, RI_OI_J3_J4l,&
        Z1, Z3, ZJ2f, ZJ2p, ZJ0, ZJ4l, &
        ZJ2fau,ZJ2pau,ZJ0au,&
        rRI_CI_J1_J2f, iRI_CI_J1_J2f, rRI_CI_J1_J2p, iRI_CI_J1_J2p, rRI_CI_J1_J0p, iRI_CI_J1_J0p,&
        rRI_CI_J3_J2f, iRI_CI_J3_J2f, rRI_CI_J3_J2p, iRI_CI_J3_J2p,&
        rRI_CI_J0p_J1, iRI_CI_J0p_J1, rRI_CI_J2f_J1, iRI_CI_J2f_J1, &
        rRI_CI_J2p_J1, iRI_CI_J2p_J1,&
        rRI_CI_J2p_J3, iRI_CI_J2p_J3, rRI_CI_J2f_J3, iRI_CI_J2f_J3, rRI_CI_J4l_J3, iRI_CI_J4l_J3, &
        rRI_CI_J4l_J5, iRI_CI_J4l_J5, H)
        H(1,:) = Zer
        H(:,1) = Zer
        do ix=nosgs+nbs+1, nosgs+nbs+nau*5 
            write(6,*) " Setting column and row ", ix, " to zero"
            H(ix,:) = Zer
            H(:,ix) = Zer 
        end do
        call cpu_time(tst)
        call zdiag(ntotau+1, H, HD, V)
        call cpu_time(tnd)
        write(6,*) "The diagionalization took", tnd-tst," seconds"
        sparse = 0.0
        sparseH = 0.0
        do ix=1,ntotau+1
            do iy=1, ntotau+1
                if(ABS(V(ix,iy)).le.1.0d-10) then 
                    sparse = sparse + 1
                end if
                if(ABS(H(ix,iy)).le.1.0d-10) then 
                    sparseH = sparseH + 1
                end if
            end do
        end do
        open(6666,file="SPFwGSDiag_noaunogs.dat")
        write(6666,*) HD(:)
        close(6666)
        write(6,*) "V matrix is", sparse/((ntotau+1)**2) * 100, "sparse"
        write(6,*) "H matrix is", sparseH/((ntotau+1)**2) * 100, "sparse"
    end subroutine set_and_diagHnaungs
    !
    function driver(t,Eoo, ww,gaam)
        real*8 :: t, Eoo, ww, gaam 
        real*8 :: driver, nine_ham
        nine_ham = 12.1661676d-15/auoftime * 1/SQRT(2.d0*LOG(2.d0)) ! gaam/9.0d0
        driver = Eoo * EXP(-(t)**2/(gaam/1.0d0)**2) * COS(1.d0*ww*(t)) + &
        SQRT(1.0d-5)*(Eoo * EXP(-(t)**2/(9.0d0/9.0d0 * nine_ham)**2) * COS(3.0d0*ww*(t)) + &
                      Eoo * EXP(-(t)**2/(9.0d0/9.0d0 * nine_ham)**2) * COS(5.0d0*ww*(t)) + &
                      Eoo * EXP(-(t)**2/(9.0d0/9.0d0 * nine_ham)**2) * COS(7.0d0*ww*(t)) + &
                      Eoo * EXP(-(t)**2/(9.0d0/9.0d0 * nine_ham)**2) * COS(9.0d0*ww*(t)) )
    end function driver
    !
    function opa(t, too,Eopa,wwopa,gaamopa)
        real*8 :: t,too,Eopa,wwopa,gaamopa
        real*8 :: opa
        opa = Eopa * EXP(-(t-too)**2/gaamopa**2) * COS(wwopa*(t-too))
    end function
    !
    function driver_opa(t,too,Eoo,ww,gaam,Eopa,wwopa,gaamopa)
        real*8 :: t,too,Eoo,ww,gaam,Eopa,wwopa,gaamopa
        real*8 :: driver_opa
        driver_opa = driver(t,Eoo,ww,gaam) + opa(t,too,Eopa,wwopa,gaamopa)
    end function
    !
    subroutine profiledriv(thi, tlo, dtt, Eoo, ww, gaam)
        integer :: iloc
        real*8 :: thi, tlo, dtt, gaam, Eoo, ww
        open(122, file="SPF_DRIVER_"//itoa(nl)//"/driver_profile.dat", status='replace', access='sequential')
        do iloc=1, int((thi-tlo)/dtt)
            write(122,*) tlo+iloc*dtt, driver(tlo+iloc*dtt,Eoo,ww,gaam)
        end do
        close(122)
    end subroutine
    !
    subroutine profileopa(thi, tlo, dtt, too, Eoo, ww, gaam)
        integer :: iloc
        real*8 :: thi, tlo, dtt, too, gaam, Eoo, ww
        open(122, file="SPF_DRIVER_"//itoa(nl)//"/opa_profile.dat", status='replace', access='sequential')
        do iloc=1, int((thi-tlo)/dtt)
            write(122,*) tlo+iloc*dtt, opa(tlo+iloc*dtt, too, Eoo,ww,gaam)
        end do
        close(122)
    end subroutine
    !
    subroutine profiledrivopa(thi, tlo, dtt, too, Eoo, ww, gaam, Eopa, wwopa,gaamopa)
        integer :: iloc
        real*8 :: thi, tlo, dtt, too, Eoo, ww, gaam
        real*8 :: Eopa, wwopa, gaamopa
        real*8 :: dur
        dur = 0.0d0!gaam*SQRT(2*LOG(2.0d0))/2.0d0
        open(122, file="SPF_DRIVER_"//itoa(nl)//"/opa_profile.dat", status='replace', access='sequential')
        do iloc=1, int((thi-tlo)/dtt)
            write(122,*) tlo+iloc*dtt, driver_opa(tlo+iloc*dtt,too,eoo,ww,gaam,Eopa,wwopa,gaamopa)
        end do
        close(122)
    end subroutine
    !
    function Driver_int(t, dtt, Eoo, ww,gaam)
        integer :: iloc
        real*8 :: t, dtt, tt(31), Eoo, ww, gaam
        real*8 :: Driver_int
        real*8 :: fss(31)
        real*8 :: dur
        real*8, external :: rint
        dur = 0.0d0!gaam*SQRT(2.d0*LOG(2.0d0))/2.0d0
        do iloc=0,30
            tt(iloc+1) = t+iloc*dtt/30
            fss(iloc+1) = driver(tt(iloc+1),Eoo,ww,gaam)
        end do
        Driver_int =  rint(fss,1,31,10,dtt/30)
        return 
    end function Driver_int
    !
    function Opa_int(t, dtt, too, Eoo, ww,gaam)
        integer :: iloc
        real*8 :: t, dtt, too, tt(31), Eoo, ww, gaam
        real*8 :: Opa_int
        real*8 :: fss(31)
        real*8 :: dur
        real*8, external :: rint
        dur = 0.0d0!gaam*SQRT(2.d0*LOG(2.0d0))/2.0d0
        do iloc=0,30
            tt(iloc+1) = t+iloc*dtt/30
            fss(iloc+1) = opa(tt(iloc+1),too,Eoo,ww,gaam)
        end do
        Opa_int =  rint(fss,1,31,10,dtt/30)
        return 
    end function Opa_int
    !
    function DriverOpa_int(t, dtt, too, Eoo, ww,gaam,Eopa,wwopa,gaamopa)
        integer :: iloc
        real*8 :: t, dtt, tt(31), too, Eoo, ww, gaam
        real*8 :: Eopa,wwopa,gaamopa
        real*8 :: DriverOpa_int
        real*8 :: fss(31)
        real*8, external :: rint
        do iloc=0,30
            tt(iloc+1) = t+iloc*dtt/30
            fss(iloc+1) = driver_opa(tt(iloc+1),too,Eoo,ww,gaam,Eopa,wwopa,gaamopa)
        end do
        DriverOpa_int =  rint(fss,1,31,10,dtt/30)
        return 
    end function DriverOpa_int
    !
    subroutine driver_prop( Eoo, ww, gaam, nunit,ip, cinit)
        integer :: nunit, ntt, iloc
        real*8 ::  Eoo, ww, gaam
        real*8 :: ttlo, tthi, dtt
        complex*16 :: cinit(size(c))
        logical :: ip
        real :: tst, tnd
        ttlo = -5*gaam*SQRT(LOG(2.0d0))
        tthi = +5*gaam*SQRT(LOG(2.0d0))
        dtt = 2*PI/(ww) * 1/100
        ntt = int((tthi-ttlo)/dtt)
        write(6,*) "A driver propagation on ", ntt, " time steps will be done"
        c = cinit
        write(6,*) "Time loop has started"
        call cpu_time(tst)
        do iloc=1, ntt
            call SPf_step(ntotau+1, HD*Driver_int(ttlo+iloc*dtt, dtt, Eoo, ww, gaam), &
            V, dtt, c, caux, Es)
            c = caux
            if(ip.and.(nmod(iloc,20).eq.20)) then
                 write(nunit, 1986) ttlo+iloc*dtt, c(:)
            endif
        end do
        call cpu_time(tnd)
        if(.not.(ip)) write(nunit,1986) Eoo, ww, gaam, c(:)
        write(6,*) "Propagation time ",tnd-tst," seconds"
        call profiledriv(tthi, ttlo, dtt,Eoo,ww, gaam)
        1986 format((E25.8E4,5x,50000(E20.8E4,1x,E20.8E4,3x)))
    end subroutine
    !
    subroutine opa_prop( Eoopa, wwopa, gaamopa, too, nunit, cinit, tlst, ip)
        integer :: nunit, ntt, iloc
        real*8 ::  Eoopa, wwopa, gaamopa
        real*8 :: ttlo, tthi, dtt, too, tlst
        complex*16 :: cinit(size(c))
        logical :: ip
        real :: tst, tnd

        ttlo = too-5*gaamopa*SQRT(LOG(2.0d0))
        tthi = too+5*gaamopa*SQRT(LOG(2.0d0))
        if(tlst.gt.ttlo) then
            write(6,*) "In OPA only propagation there is overlap between driver and OPA"
            STOP
        end if
        dtt = 2*PI/(wwopa) * 1/100
        ntt = int((tthi-ttlo)/dtt)
        write(6,*) "An OPA propagation on ", ntt, " time steps will be done"
        do iloc=1,size(c)
            c(iloc) = cinit(iloc)*EXP(-IM*(ttlo-tlst)*Es(iloc))
        end do
        write(6,*) "Time loop has started"
        call cpu_time(tst)
        do iloc=1, ntt
            call SPF_step(ntotau+1, HD*OPA_int(ttlo+iloc*dtt, dtt, too, Eoopa, wwopa, gaamopa), &
            V, dtt, c, caux, Es)
            c = caux
            if(ip.and.(nmod(iloc,20).eq.20)) then
                 write(nunit, 1987) ttlo+iloc*dtt, c(:)
            endif
        end do
        call cpu_time(tnd)
        if(.not.(ip)) write(nunit,1986) Eoopa, wwopa, gaamopa, c(:)
        write(6,*) "Propagation time ",tnd-tst," seconds"
        call profileopa(tthi, ttlo, dtt, too, Eoopa,wwopa, gaamopa)
        1986 format((3E25.8E4,5x,50000(E20.8E4,1x,E20.8E4,3x)))
        1987 format((1E25.8E4,5x,50000(E20.8E4,1x,E20.8E4,3x)))
    end subroutine
    !
    subroutine driver_opa_prop(too, Eoo, ww, gaam, Eoopa, wwopa,gaamopa, nunit,ip,ags)
        integer :: nunit, ntt, iloc
        real*8 :: too, Eoo, ww, gaam
        real*8 :: Eoopa, wwopa, gaamopa
        real*8 :: ttlo, tthi, dtt, tdel, nrm
        logical :: ip, ags
        real :: tst, tnd
        if(5*gaam*SQRT(LOG(2.0d0)).lt.too-5*gaamopa*SQRT(LOG(2.0d0))) then
            write(6,*) "Non overlapping driver and opa, doing separate propagations"
            ttlo = -5*gaam*SQRT(LOG(2.0d0))
            tthi = +5*gaam*SQRT(LOG(2.0d0))
            tdel = too-5*gaamopa*SQRT(LOG(2.0d0))-tthi
            dtt = 2*PI/(ww) * 1/100
            ntt = int((tthi-ttlo)/dtt)
            write(6,*) "When propagating fundamental and harmonics ", ntt, " time steps will be done"
            if(ags) then
                write(6,*) "Initializing all in the GS"
                c = Zer
                c(1) = Zone
            else 
                write(6,*) "Initializing the J=1 states to the given amplitudes and adjusting the phase to the past."
                c = Zer
                nrm = 0.d0
                do iloc=2,8
                    c(iloc) = A(iloc-1)*EXP(+IM*Es(iloc)*(tdel+tthi-ttlo))
                    nrm = 1.0d0 !nrm + ABS(A(iloc-1))**2
                end do
                c = c/sqrt(nrm)
                write(6,*) "Initia States Amplitude"
                do iloc=2,8
                    write(6,*) c(iloc) 
                end do
            endif
            write(6,*) "Time loop has started"
            call cpu_time(tst)
            do iloc=1, ntt
                call SPf_step(ntotau+1, HD*DriverOpa_int(ttlo+iloc*dtt, dtt, too, Eoo, ww, gaam, &
                Eoopa,wwopa,gaamopa), V, dtt, c, caux, Es)
                c = caux
                if(ip.and.(nmod(iloc,20).eq.20)) then
                    !write(6,*) "Writing to unit in line 679"
                    write(nunit, 1986) ttlo+iloc*dtt, c(:)
                endif
            end do
            call cpu_time(tnd)
            write(6,*) "Time of the driver ",tnd-tst," seconds"
            write(6,*) "----------------- fin of driver and harmonics -----------------"
            ttlo = too-5*gaamopa*SQRT(LOG(2.0d0))
            tdel = ttlo-tthi
            do iloc=1, size(c)
                c(iloc) = EXP(-Im*Es(iloc)*tdel)*c(iloc)
            end do
            tthi = too+5*gaamopa*SQRT(LOG(2.0d0))
            dtt = 2*PI/(wwopa) * 1/100
            ntt = int((tthi-ttlo)/dtt)
            write(6,*) "When propagating OPA ", ntt, " time steps will be done"
            write(6,*) "Time loop has started"
            call cpu_time(tst)
            do iloc=1, ntt
                call SPf_step(ntotau+1, HD*DriverOpa_int(ttlo+iloc*dtt, dtt, too, Eoo, ww, gaam, &
                Eoopa,wwopa,gaamopa),&
                V, dtt, c, caux, Es)
                c = caux
                if(ip.and.(nmod(iloc,20).eq.20)) then
                    !write(6,*) "Writing to unit in line 702"
                    write(nunit, 1986) ttlo+iloc*dtt, c(:)
                endif
            end do
            call cpu_time(tnd)
            if(.not.(ip)) then
                !write(6,*) "Writing to unit in line 707"
                write(nunit,1986) too, c(:)
                !write(nunit,*) "----------------- fin of OPA -----------------"
            end if
            write(6,*) "Time of the OPA ",tnd-tst," seconds"
            write(6,*) "----------------- fin of OPA -----------------"
            write(6,*) "writing profile"
            ttlo = -5*gaam*SQRT(LOG(2.0d0))
            tthi = too+5*gaamopa*SQRT(LOG(2.0d0))
            dtt = 2*PI/(max(ww,wwopa)) * 1/100
            call profiledrivopa(tthi, ttlo, dtt, too,Eoo,ww,gaam,Eoopa,wwopa,gaamopa)
        else
            write(6,*) "Overlapping driver and opa, doing joint propagation"
            ttlo = -5*gaam*SQRT(LOG(2.0d0))
            tthi = too+5*gamopa*SQRT(LOG(2.0d0))
            dtt = 2*PI/(max(ww,wwopa)) * 1/100
            ntt = int((tthi-ttlo)/dtt)
            write(6,*) "When propagating ", ntt, " time steps will be done"
            if(ags) then
                write(6,*) "Initializing all in the GS"
                c = Zer
                c(1) = Zone
            else 
                write(6,*) "Initializing the J=1 states to the given amplitudes and adjusting the phase to the past."
                c = Zer
                nrm = 0.d0
                do iloc=2,8
                    c(iloc) = A(iloc-1)*EXP(+IM*Es(iloc)*(tdel+tthi-ttlo))
                    nrm = 1.0d0 !nrm + ABS(A(iloc-1))**2
                end do
                c = c/sqrt(nrm)
                write(6,*) "Normalized vector"
                do iloc=2,8
                    write(6,*) c(iloc) 
                end do
            endif
            write(6,*) "Time loop has started"
            call cpu_time(tst)
            do iloc=1, ntt
                call SPf_step(ntotau+1, HD*DriverOpa_int(ttlo+iloc*dtt, dtt, too, Eoo, ww, gaam,&
                Eoopa,wwopa,gaamopa), V, dtt, c, caux, Es)
                c = caux
                if(ip.and.(nmod(iloc,20).eq.20)) then
                        !write(6,*) "Writing to unit in line 750"
                        write(nunit, 1986) ttlo+iloc*dtt, c(:)
                endif
            end do
            call cpu_time(tnd)
            if(.not.(ip)) then
                    !write(6,*) "Writing to unit in line 755"
                    write(nunit,1986) too, c(:)
            end if
            write(6,*) "Time of conjoined propagation ",tnd-tst," seconds"
            write(6,*) "----------- fin of driver and harmonics and OPA -----------------"
            call profiledrivopa(tthi, ttlo, dtt, too,Eoo,ww,gaam,Eoopa,wwopa,gaamopa)
        endif
        1986 format((E25.8E4,5x,50000(E20.8E4,1x,E20.8E4,3x)))
    end subroutine
    !
    subroutine driver_and_opa_prop(too, Eoo, ww, gaam, Eoopa, wwopa,gaamopa, nunit,ip,ags)
        integer :: nunit, iloc
        real*8 :: too, Eoo, ww, gaam
        real*8 :: Eoopa, wwopa, gaamopa
        real*8 :: ttlo, tthi, tdel, nrm
        logical :: ip, ags
        if(+5*gaam*SQRT(LOG(2.d0)).gt.too-5*gaamopa*SQRT(LOG(2.d0))) then 
            write(6,*) "Overlapping lasers, call driver_opa_prop"
            STOP
        end if
        if(ags) then
            write(6,*) "Initializing all in the GS"
            c = Zer
            c(1) = Zone
        else 
            write(6,*) "Initializing the J=1 states to the given amplitudes and adjusting the phase to the past."
            c = Zer
            nrm = 0.d0
            tthi = 5*gaam*SQRT(LOG(2.d0))
            tdel = too-5*gaamopa*SQRT(LOG(2.d0))-5*gaam*SQRT(LOG(2.d0))
            do iloc=2,8
                c(iloc) = A(iloc-1)*EXP(+IM*Es(iloc)*(tdel+tthi-ttlo))
                nrm = 1.0d0 !nrm + ABS(A(iloc-1))**2
            end do
            c = c/sqrt(nrm)
            write(6,*) "Initia States Amplitude"
            do iloc=2,8
                write(6,*) c(iloc) 
            end do
        endif
        call driver_prop(Eoo,ww,gaam,nunit,ip,c)
        call opa_prop(Eoopa,wwopa,gaamopa,too,nunit,c,+5*gaam*SQRT(LOG(2.d0)),ip)
    end subroutine
    !
    subroutine driver_excitation(Eoo, ww, gaam)
        real*8 :: Eoo, ww, gaam
        logical :: pf 
        pf = .true.
        write(6,*) "Computing 0 delay coeffiecients"
        write(6,*) "State ", "gs"
        if(pf) open(44,file='SPF_DRIVER_'//itoa(nl)//'/only_driver_csgs.dat')
        c = Zer
        c(1) = Zone
        call driver_prop( Eoo, ww, gaam, 44, pf, c)
        close(44)
    end subroutine driver_excitation
    !
    subroutine driver_opa_excitation(Eoo, ww, gaam, Eoopa, wwopa,gaamopa, too)
        real*8 :: Eoo, ww, gaam
        real*8 :: Eoopa, wwopa, gaamopa, too
        logical :: pf, ags 
        pf = .true.
        ags = .true.
        if(pf) open(44,file='SPF_DRIVER_'//itoa(nl)//'/joint_prop_csgs.dat')
        call driver_opa_prop(too, Eoo, ww, gaam,Eoopa,wwopa,gaamopa, 44, pf,ags)
        close(44)
    end subroutine driver_opa_excitation
    !
    subroutine driver_and_opa_excitation(Eoo, ww, gaam, Eoopa, wwopa,gaamopa, too)
        real*8 :: Eoo, ww, gaam
        real*8 :: Eoopa, wwopa, gaamopa, too
        logical :: pf, ags 
        pf = .true.
        ags = .true.
        if(pf) open(44,file='SPF_DRIVER_'//itoa(nl)//'/sep_prop_csgs_sep.dat')
        call driver_and_opa_prop(too, Eoo, ww, gaam,Eoopa,wwopa,gaamopa, 44, pf,ags)
        close(44)
    end subroutine driver_and_opa_excitation
    !
    subroutine opa_single_state_excitation(Eoopa, wwopa,gaamopa)
        integer :: iloc
        real*8 :: wwopa,Eoopa,gaamopa
        logical :: ip
        ip = .false.
        do iloc=2,8
            if(ip) then
                open(15,file="SPF_DRIVER_"//itoa(nl)//"/t_coeffs_state"//itoa(iloc-1)//".dat")
            else
                open(15,file="SPF_DRIVER_"//itoa(nl)//"/tend_coeffs_state"//itoa(iloc-1)//".dat")
            endif
            c = Zer
            c(iloc) = Zone
            write(6,*) "Propagating State ", iloc
            call opa_prop(Eoopa,wwopa,gaamopa,0.d0,15,c,-5*gaamopa*SQRT(LOG(2.d0)),ip)
            csf(iloc-1,:) = c(:)
            close(15)
        end do
    end subroutine
    !
    subroutine opa_delay_loop(Eoo,ww,gaam,Eoopa,wwopa,gaamopa,tooM)
        real*8 :: Eoo, ww, gaam
        real*8 :: Eoopa, wwopa, gaamopa, tooM, dtoop, too
        integer :: itoloc, ntoot
        logical :: pf, ags
        pf = .false.
        ags = .true.
        open(18, file='SPF_DRIVER_'//itoa(nl)//'/tdel_loop_cs.dat')
        ntoot = 100
        dtoop = (tooM-10.d-15/auoftime)/ntoot
        do itoloc=0,ntoot
            too = 10.d-15/auoftime+itoloc*dtoop
            call driver_opa_prop(too,Eoo,ww,gaam,Eoopa,wwopa,gaamopa,18,pf,ags)
        end do
        close(18)
    end subroutine opa_delay_loop
    !
    subroutine driver_params_loop(Eoo, ww, gaam)
        integer :: ist
        real*8 :: Eoo, ww, gaam
        logical :: pf 
        pf = .false.
        write(6,*) "Computing the final population with variable parameters 33 propagations will be carried."
        write(6,*) "State ", "gs"
        
        open(44,file='SPF_DRIVER_'//itoa(nl)//'/Driver_pulse_gamma_loop_amplitudes.dat')
        do ist=0,6
            c = Zer
            c(1) = Zone
            call driver_prop( Eoo, ww, 0.98*gaam+ist*(1.02-0.98)*gaam/6.0d0, 44, pf,c)
        end do
        close(44)
        open(44,file='SPF_DRIVER_'//itoa(nl)//'/Driver_pulse_intensity_loop_amplitudes.dat')
        do ist=0,6
            c = Zer
            c(1) = Zone
            call driver_prop( sqrt(0.001d0)*Eoo+ist*(SQRT(4.d0)-SQRT(0.001d0))*Eo/6.0d0, ww, gaam,&
                               44, pf,c)
        end do
        close(44)
        open(44,file='SPF_DRIVER_'//itoa(nl)//'/Driver_frequency_loop_amplitudes.dat')
        do ist=0,6
            c = Zer
            c(1) = Zone
            call driver_prop(Eoo, 0.9862*ww+ist*(1.002-0.9862)*ww/6.0d0, gaam, 44, pf, c)
        end do
    end subroutine driver_params_loop
    !
    subroutine computed_cs()
        integer :: ist
        do ist=1,nos
            open(68, file="SPF_DRIVER_"//itoa(nl)//"/tend_coeffs_state"//itoa(ist)//".dat")
            write(6,*) "Reading state ", ist 
            read(68,1987) Epa, wopa, gamopa, csf(ist,:)
            write(6,*) csf(ist,1)
            close(68)
        end do
        1987 format((3E25.8E4,5x,50000(E20.8E4,1x,E20.8E4,3x)))
    end subroutine computed_cs
    !
    subroutine time_delay_phase(gaamopa,tlast, dtoo,tofin, nunit)
        real*8 :: tlast, dtoo, tofin, tbel, too, gaamopa
        integer :: nunit
        integer :: iloc, jloc, nend
        complex*16 :: cloc(size(c))
        tbel = tlast+5*gaamopa*SQRT(LOG(2.0))
        nend = int((tofin)/dtoo)
        write(6,*) nend
        do iloc=0, nend
            cloc = Zer
            too = tbel+iloc*dtoo
            do jloc=1,nosgs-1
                cloc = cloc + A(jloc)*EXP(-Im*(too-5*gaamopa*SQRT(LOG(2.0))-tlast)*Es(jloc+1))*csf(jloc,:)
            end do
            write(nunit,1986) too, cloc
        end do
        1986 format((1E25.8E4,5x,50000(E20.8E4,1x,E20.8E4,3x)))
    end subroutine time_delay_phase
    !
    subroutine OPA_time_delay_phase(Eoo,gaam, ww, Eoopa,gaamopa,wwopa,nunit)
        integer :: nunit
        logical :: dri_prop
        real*8 :: Eoo, gaam, ww
        real*8 :: Eoopa, gaamopa, wwopa, thoi
        real*8 :: toof, dtoo
        INQUIRE(FILE='SPF_DRIVER_'//itoa(nl)//"/init_amps_driver.dat", EXIST=dri_prop)
        if(dri_prop.and.(Eoo.eq.Eo).and.(ww.eq.w).and.(gaam.eq.gam)) then 
            write(6,*) "Existing driver propagation with parameters"
            write(6,*) "Gamma: ", gam, " Eo:", Eo, " w:",w
            write(6,*) "Reading into csf."
            open(45,file='SPF_DRIVER_'//itoa(nl)//"/init_amps_driver.dat")
            do mi1=1,nos
                read(45,1988) thoi, A(mi1)
            end do
            write(6,*) "End time of driver ", thoi*auoftime/1.0d-15
            close(45)
        else
            write(6,*) "No driver propagation found. Propagating it and saving to file init_amps_driver."
            c= Zer
            c(1) = Zone
            call driver_prop(Eoo,ww,gaam,44,.false.,c)
            open(45,file='SPF_DRIVER_'//itoa(nl)//"/init_amps_driver.dat")
            do mi1=2,8
                A(mi1-1) = c(mi1)
                write(45,1988) 5*gam*SQRT(LOG(2.d0)),c(mi1)
            end do
            close(45)
        end if

        write(6,*) "Progating an OPA for each state and saving in a single file."
        call opa_single_state_excitation(Eoopa,wwopa,gaamopa)
        dtoo = 243.914d0
        toof = 5d-12/auoftime
        write(6,*) "Setting phases for delays up to 5 ps with spacing ", dtoo*auoftime/1e-15, " fs"
        call time_delay_phase(gaamopa,thoi, dtoo, toof, nunit)
        1988 format((E25.8E4,5x,E20.8E4,1x,E20.8E4))
    end subroutine
    !
    subroutine Spec_freqs_loop(Eoo,gaam,ww,Eoopa,gaamopa,nunit, &
    Tot_Iter_w,lower_scale,dw,name_offset)
        real*8 :: Eoo, gaam, ww, Eoopa, gaamopa
        integer :: nunit
        integer :: Tot_Iter_w
        integer :: iter_w, name_offset
        real*8 :: wi, dw
        real*8 :: lower_scale
        open(12,file='SPF_DRIVER_'//itoa(nl)//'/OPA_Frequency_loop_params'//itoa(name_offset)//'.dat')
        do iter_w=0,Tot_Iter_w
            open(nunit,file='SPF_DRIVER_'//itoa(nl)//'/OPA_Frequency_Variation'//itoa(iter_w+name_offset)//'.dat',status='replace')
            wi = lower_scale + iter_w*dw
            write(12,*) Eoopa, wi, gaamopa
            write(6,*)'Running a spectrogram at frequency', wi*eVperAU
            call OPA_time_delay_phase(Eoo,gaam,ww,Eoopa,gaamopa,wi,nunit)
            close(nunit)
        end do
        close(12)
    end subroutine
    !
    subroutine Spec_Intensity_loop(Eoo,gaam,ww,gaamopa,wwopa,nunit, &
    Tot_Iter_e, lower_scale,dI,name_offset)
        real*8 :: Eoo, gaam, ww, gaamopa, wwopa
        integer :: nunit
        integer :: Tot_Iter_e
        integer :: iter_e, name_offset
        real*8 :: Ei,dI
        real*8 :: lower_scale
        open(12,file='SPF_DRIVER_'//itoa(nl)//'/OPA_Itensity_loop_params'//itoa(name_offset)//'.dat')
        do iter_e=0,Tot_Iter_e
            open(nunit,file='SPF_DRIVER_'//itoa(nl)//'/OPA_Intensity_Variation'//itoa(iter_e+name_offset)//'.dat',status='replace')
            Ei = SQRT(lower_scale+dI*iter_e)
            write(12,*) Ei, wwopa, gaamopa
            write(6,*)'Running a spectrogram at Intensity', Ei**2 * auofI/1d12
            call OPA_time_delay_phase(Eoo,gaam,ww,Ei,gaamopa,wwopa,nunit)
            close(nunit)
        end do
        close(12)
    end subroutine 
    !
    subroutine Spec_duration_loop(Eoo,gaam,ww,Eoopa,wwopa,nunit, &
    Tot_Iter_g,lower_scale,dg,name_offset)
        real*8 :: Eoo, gaam, ww, Eoopa, wwopa
        integer :: nunit
        integer :: Tot_Iter_g
        integer :: iter_g, name_offset
        real*8 :: gi,dg
        real*8 :: lower_scale
        open(12,file='SPF_DRIVER_'//itoa(nl)//'/OPA_Duration_loop_params'//itoa(name_offset)//'.dat')
        do iter_g=0, Tot_Iter_g
            open(nunit,file='SPF_DRIVER_'//itoa(nl)//'/OPA_Duration_Variation'//itoa(iter_g+name_offset)//'.dat',status='replace')
            gi = lower_scale + dg*iter_g
            write(12,*) Eoopa, wwopa, gi
            write(6,*) 'Running a spectrogram at duration ',gi*SQRT(2.d0*LOG(2.d0))*auoftime/1d-15
            call OPA_time_delay_phase(Eoo,gaam,ww,Eoopa,gi,wwopa,nunit)
            close(nunit)
        end do
        close(12)
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
subroutine Init_big_mats(nnbs, nnconts, Rs, chanc, chanb, file_name)
    implicit none
    integer :: nnbs, nnconts, chanc, chanb, ii, jj, kk, ll
    real*8 :: Rs(nnconts, nnbs, chanc, chanb), Rtmp(chanc*chanb)
    character(len=143) :: file_name
    character(len=149) :: fi
    character(len=2) :: indx
    !write(6,*) shape(Rs)
    Rs = 0.d0
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
!
subroutine Init_big_matsW(nnbs, nnconts, Rs, chanc, chanb, file_name)
    implicit none
    integer :: nnbs, nnconts, chanc, chanb, ii, jj, kk, ll
    real*8 :: Rs(nnbs, nnconts, chanc, chanb), Rtmp(chanc*chanb)
    character(len=143) :: file_name
    character(len=149) :: fi
    character(len=2) :: indx, indx1
    !write(6,*) shape(Rs)
    Rs = 0.0d0
    do jj=1, nnbs
        write(indx,'(I2)') jj-1
        !write(6,*) TRIM(ADJUSTL(indx)), ADJUSTL(indx)
        fi = TRIM(file_name)//TRIM(ADJUSTL(indx))//'.dat'
        !write(6,*) fi
        open(15, file=fi, status='old', access='sequential')
    do ii=1, nnconts    
        read(15,*) Rtmp
        write(indx1,'(I2)') ii
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
!
! call set_and_diagHnaungs()
    ! 
    !open(667,file="J1Hr_gs.dat")
    !open(778,file="J1Hi_gs.dat")
    !write(667,*) (REAL(H(1,mi2)),mi2=1,ntotau)
    !write(778,*) (AIMAG(H(1,mi2)),mi2=1,ntotau)
    !end do
    !close(667)
    !close(778)
    !open(889,file="J3Hr_SPFwGSnaungs.dat")
    !open(990,file="J3Hi_SPFwGSnaungs.dat")
    !do ii1=9+24,8+24+30
    !    write(889,*) (REAL(H(ii1,mi2)),mi2=1,ntotau)
    !    write(990,*) (AIMAG(H(ii1,mi2)),mi2=1,ntotau)
    !end do
    !close(889)
    !close(990)
    !stop
    
    ! 

    !call driver_opa_excitation(Eo,w,gam,Epa,wopa,gamopa,tdelay)