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
    integer :: mi, mi1, mi2, ii1
    real*8 :: Elo, dEcont
    real*8 :: Z1(nos1, 5), EJ1(nos1)
    real*8 :: Z3(nos3,6), EJ3(nos3)
    real*8 :: ZJ0(nbs0,2), EnJ0(nbs0)
    real*8 :: ZJ2f(nbs2f,3), EnJ2f(nbs2f)
    real*8 :: ZJ2p(nbs2p,3), EnJ2p(nbs2p)
    real*8 :: ZJ4l(nbs4,6), EnJ4l(nbs4)
    real*8 :: rRI_CI_J0p_J1(nconts,nbs0,5,2)
    real*8 :: iRI_CI_J0p_J1(nconts,nbs0,5,2)
    real*8 :: rRI_CI_J2f_J1(nconts,nbs2f,5,3), rRI_CI_J2f_J3(nconts,nbs2f,6,3)
    real*8 :: iRI_CI_J2f_J1(nconts,nbs2f,5,3), iRI_CI_J2f_J3(nconts,nbs2f,6,3)
    real*8 :: rRI_CI_J2p_J1(nconts,nbs2p,5,3), rRI_CI_J2p_J3(nconts,nbs2p,6,3)
    real*8 :: iRI_CI_J2p_J1(nconts,nbs2p,5,3), iRI_CI_J2p_J3(nconts,nbs2p,6,3)
    real*8 :: rRI_CI_J4l_J3(nconts,nbs4,6,6), rRI_CI_J4l_J5(nconts,nbs4,6,6)
    real*8 :: iRI_CI_J4l_J3(nconts,nbs4,6,6), iRI_CI_J4l_J5(nconts,nbs4,6,6)
    real*8 :: RI_OI_J1_J2f(nos1, nbs2f, 3, 5), RI_OI_J1_J2p(nos1, nbs2p, 3,5)
    real*8 :: RI_OI_J1_J0p(nos1, nbs0,2, 5)
    real*8 :: RI_OI_J3_J2f(nos3, nbs2f, 3, 6), RI_OI_J3_J2p(nos3, nbs2p, 3, 6)
    real*8 :: RI_OI_J3_J4l(nos3, nbs4, 6, 6)
    real*8 :: Eo, w, gam, evec(3)
    real*8 :: Ar, Ai
    complex*16 :: A(nos)
    complex*16 :: H(ntot, ntot), V(ntot,ntot), c(ntot), caux(ntot), csf(nos, ntot)
    real*8 :: Es(ntot), HD(ntot)
    complex*16 :: H1(ntot1, ntot1), V1(ntot1,ntot1), c1(ntot1), caux1(ntot1), cfs1(4, ntot1)
    real*8 :: Es1(ntot1), HD1(ntot1)
    integer :: nt
    !real*8 :: toolo, toofo, dtoo, too
    !integer :: noo
    real :: start, finish
    character(len=143) :: file_root
    character(len=4) :: nts
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
    !-----------------------------------------------------------------------------------
    

    ! For timing purposes
    call cpu_time(start)

    ! Retrieving from precomputed files.
    
    ! Loading the Zcoefficients
    call load_Zs()
    !
    
    ! Loading the radial elements
    call load_rads()
    !

    ! Laser parameters
    open(98, file='./SPF_J3_'//itoa(nl)//'/laser_params_two_photon.dat', status='old')
    read(98,*) Eo, w, gam
    read(98,*) evec(:)
    read(98,*) nlf
    close(98)
    print *, gam, Eo
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
    do ii1=1, nos1
        write(6,*) REAL(A(ii1)), AIMAG(A(ii1)), EJ1(ii1)
    end do

    do ii1=1, nos3
        write(6,*) REAL(A(ii1+nos1)), AIMAG(A(ii1)+nos1), EJ3(ii1)
    end do

    ! Put all energies in a single array
    call all_Energies()

    ! Initialize and diagonalize the dipole op.
    call set_and_diagH()
    
    !open(667,file="J1Hr_SPF.dat")
    !open(778,file="J1Hi_SPF.dat")
    !do ii1=8,7+24
    !    write(667,*) (REAL(H(ii1,mi2)),mi2=1,ntot)
    !    write(778,*) (AIMAG(H(ii1,mi2)),mi2=1,ntot)
    !end do
    !close(667)
    !close(778)

    !open(889,file="J3Hr_SPF.dat")
    !open(990,file="J3Hi_SPF.dat")
    !do ii1=24+7+1,24+7+30
    !    write(889,*) (REAL(H(ii1,mi2)),mi2=1,ntot)
    !    write(990,*) (AIMAG(H(ii1,mi2)),mi2=1,ntot)
    !end do
    !close(889)
    !close(990)
    !stop
    
    ! call set_and_diagH_ignorelower()
    ! No energy density dipole
    ! call set_and_diagH_nED()
    ! Call the subroutine for the zero delay coefficients
    call no_delay_J3(Eo, w, gam)
    
    call cpu_time(finish)
    write(6,*) "Total calculation took ", finish-start, " seconds."
    stop
    
    contains
    !
    subroutine Load_Zs()
        open(10, file='init_bounds_J1.dat')
        do mm=1,nos1
            read(10,*) EJ1(mm), (Z1(mm,kk), kk=1,5)
            !write(6,*) Ebs(mm), (Zo(mm,kk), kk=1,5)
        end do
        close(10)

        open(10, file='init_bounds_J3.dat')
        do mm=1,nos3
            read(10,*) EJ3(mm), (Z3(mm,kk), kk=1,6)
            !write(6,*) Ebs(mm), (Zo(mm,kk), kk=1,5)
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
            open(98, file='./SPF_J3_'//itoa(nl)//'/init_amplitudes_expE.dat', status='old')
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
            open(98, file='./SPF_J3_'//itoa(nl)//'/init_amplitudes.dat', status='old')
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
        Es(1:nos1) = EJ1/CMperAU
        Es1(1:nos1) = EJ1/CMperAU
        Es(nos1+1:nos) = EJ3/CMperAU
        Es(nos+1:nbs2f+nos) = EnJ2f/CMperAU
        Es1(nos1+1:nbs2f+nos1) = EnJ2f/CMperAU
        Es(nbs2f+nos+1: nbs2f+nbs2p+nos) = EnJ2p/CMperAU
        Es1(nbs2f+nos1+1:nbs2f+nbs2p+nos1) = EnJ2p/CMperAU
        Es(nbs2f+nbs2p+nos+1:nos+nbs2f+nbs2p+nbs0) = EnJ0/CMperAU
        Es1(nbs2f+nbs2p+nos1+1:nos1+nbs2f+nbs2p+nbs0) = EnJ0/CMperAU
        Es(nos+nbs2f+nbs2p+nbs0+1:nos+nbs) = Enj4l/CMperAU
        do mi1=1,5
        Es(nbs+nos+(mi1-1)*ncon+1:nbs+nos+(mi1-1)*ncon+ncon) = (/ ((Elo+mi*dEcont)/CMperAU, mi=1, ncon) /)
        Es1(nbs1+nos1+(mi1-1)*ncon+1:nbs1+nos1+(mi1-1)*ncon+ncon) = (/ ((Elo+mi*dEcont)/CMperAU, mi=1, ncon) /)
        end do
        do mi1=1,6
        Es(nbs+nos+(mi1-1)*ncon+5*ncon+1:nbs+nos+(mi1-1)*ncon+ncon+5*ncon) = (/ ((Elo+mi*dEcont)/CMperAU, mi=1, ncon) /)
        Es1(nbs1+nos1+(mi1-1)*ncon+5*ncon+1:nbs1+nos1+(mi1-1)*ncon+ncon+5*ncon) = (/ ((Elo+mi*dEcont)/CMperAU, mi=1, ncon) /)
        end do
        do mi1=1,6
        Es(nbs+nos+(mi1-1)*ncon+5*ncon+6*ncon+1:nbs+nos+(mi1-1)*ncon+ncon+5*ncon+6*ncon) &
        = (/ ((Elo+mi*dEcont)/CMperAU, mi=1, ncon) /)
        end do
        !
        open(41, file='SPF_J3_'//itoa(nl)//'/energies.dat')
        write(41,1981) Es(:)
        close(41)
        1981 format(3500(E20.7E4,1x))
    end subroutine all_Energies
    !
    subroutine set_and_diagH()
        integer :: ix, iy
        real*8 :: sparse, sparseH
        real :: tst, tnd
        call TPMatED(evec, Elo, dEcont,&
        RI_OI_J1_J0p, RI_OI_J1_J2f, RI_OI_J1_J2p, RI_OI_J3_J2f,RI_OI_J3_J2p, RI_OI_J3_J4l,&
        Z1, Z3, ZJ2f, ZJ2p, ZJ0, ZJ4l, &
        rRI_CI_J0p_J1, iRI_CI_J0p_J1, rRI_CI_J2f_J1, iRI_CI_J2f_J1, &
        rRI_CI_J2p_J1, iRI_CI_J2p_J1,&
        rRI_CI_J2p_J3, iRI_CI_J2p_J3, rRI_CI_J2f_J3, iRI_CI_J2f_J3, rRI_CI_J4l_J3, iRI_CI_J4l_J3, &
        rRI_CI_J4l_J5, iRI_CI_J4l_J5, H)
        call cpu_time(tst)
        call zdiag(ntot, H, HD, V)
        call cpu_time(tnd)
        write(6,*) "The diagionalization took", tnd-tst," seconds"
        sparse = 0.0
        sparseH = 0.0
        do ix=1,ntot
            do iy=1, ntot
                if(ABS(V(ix,iy)).le.1.0d-10) then 
                    sparse = sparse + 1
                end if
                if(ABS(H(ix,iy)).le.1.0d-10) then 
                    sparseH = sparseH + 1
                end if
            end do
        end do
        open(6666,file="SPF_diag.dat")
        write(6666,*) HD(:)
        close(6666)
        write(6,*) "V matrix is", sparse/(ntot**2) * 100, "sparse"
        write(6,*) "H matrix is", sparseH/(ntot**2) * 100, "sparse"
    end subroutine set_and_diagH
    !
    subroutine set_and_diagH_ignorelower()
        integer :: ix, iy
        real*8 :: sparse, sparseH
        real :: tst, tnd
        call TPMatED(evec, Elo, dEcont,&
        RI_OI_J1_J0p, RI_OI_J1_J2f, RI_OI_J1_J2p, RI_OI_J3_J2f,RI_OI_J3_J2p, RI_OI_J3_J4l,&
        Z1, Z3, ZJ2f, ZJ2p, ZJ0, ZJ4l, &
        rRI_CI_J0p_J1, iRI_CI_J0p_J1, rRI_CI_J2f_J1, iRI_CI_J2f_J1, &
        rRI_CI_J2p_J1, iRI_CI_J2p_J1,&
        rRI_CI_J2p_J3, iRI_CI_J2p_J3, rRI_CI_J2f_J3, iRI_CI_J2f_J3, rRI_CI_J4l_J3, iRI_CI_J4l_J3, &
        rRI_CI_J4l_J5, iRI_CI_J4l_J5, H)
        do ix=1, ntot
            if(Es(ix).lt.MINVAL(Es(1:7))) then 
                H(ix,:) = Zer
                H(:,ix) = Zer
                write(6,*) "Ignoring ",ix
            end if
        end do
        call cpu_time(tst)
        call zdiag(ntot, H, HD, V)
        call cpu_time(tnd)
        write(6,*) "The diagionalization took", tnd-tst," seconds"
        sparse = 0.0
        sparseH = 0.0
        do ix=1,ntot
            do iy=1, ntot
                if(ABS(V(ix,iy)).le.1.0d-10) then 
                    sparse = sparse + 1
                end if
                if(ABS(H(ix,iy)).le.1.0d-10) then 
                    sparseH = sparseH + 1
                end if
            end do
        end do
        write(6,*) "V matrix is", sparse/(ntot**2) * 100, "sparse"
        write(6,*) "H matrix is", sparseH/(ntot**2) * 100, "sparse"
    end subroutine set_and_diagH_ignorelower
    !
    subroutine set_and_diagH_nED()
        integer :: ix, iy
        real*8 :: sparse, sparseH, tst, tnd
        write(6,*) "Initializing the H matrix"
        call TPMatNED(evec, Elo, dEcont,&
        RI_OI_J1_J0p, RI_OI_J1_J2f, RI_OI_J1_J2p, RI_OI_J3_J2f,RI_OI_J3_J2p, RI_OI_J3_J4l,&
        Z1, Z3, ZJ2f, ZJ2p, ZJ0, ZJ4l, &
        rRI_CI_J0p_J1, iRI_CI_J0p_J1, rRI_CI_J2f_J1, iRI_CI_J2f_J1, &
        rRI_CI_J2p_J1, iRI_CI_J2p_J1,&
        rRI_CI_J2p_J3, iRI_CI_J2p_J3, rRI_CI_J2f_J3, iRI_CI_J2f_J3, rRI_CI_J4l_J3, iRI_CI_J4l_J3, &
        rRI_CI_J4l_J5, iRI_CI_J4l_J5, H)
        write(6,*) "Diagonalizing matrix"
        call cpu_time(tst)
        call zdiag(ntot, H, HD, V)
        call cpu_time(tnd)
        write(6,*) "Took ", tnd-tst, " seconds to diagonalize."
        sparse = 0.0
        sparseH = 0.0
        do ix=1,ntot
            do iy=1, ntot
                if(ABS(V(ix,iy)).le.1.0d-10) then 
                    sparse = sparse + 1
                end if
                if(ABS(H(ix,iy)).le.1.0d-10) then 
                    sparseH = sparseH + 1
                end if
            end do
        end do
        write(6,*) "V matrix is", sparse/(ntot**2) * 100, "sparse"
        write(6,*) "H matrix is", sparseH/(ntot**2) * 100, "sparse"
    end subroutine set_and_diagH_nED
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
        open(12, file="SPF_J3_"//itoa(nl)//"/laser_profile.dat", status='replace', access='sequential')
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
    !-----OLD ROUTINE NEEDS TO BE CHANGED TO INCLUDE THE J=5, but this is old-----
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
        do iloc=1, nos
            c(iloc) = A(iloc)*EXP(-Im*Es(iloc)*tlo)
        end do
        write(6,*) 'Time loop started', too
        do itt=1, nt
            call SPf_step(ntot, HD*Las(tlo+itt*dt, dt, too, Eo, w, gaam), V, dt, c, caux, Es)
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
        1958 format(201(E17.8E4, 2x))
    end subroutine 
    !-----------------------------------------------------------------------------
    !
    subroutine probs_from_cs(cs, too, nunit12, nunit32)
        integer :: iloc, nunit12, nunit32, offset
        real*8 :: probJ12(ncon), probJ32(ncon)
        complex*16 :: cs(ntot)
        real*8 :: too
        probJ12 = 0.0d0
        probJ32 = 0.0d0
        do iloc=1,5
            offset = nos+nbs+(iloc-1)*ncon
            if((iloc.eq.1).or.(iloc.eq.4)) then
                probJ12 = probJ12 + (/(ABS(cs(mi2))**2, mi2=offset+1,offset+ncon)/)
            else
                probJ32 = probJ32 + (/(ABS(cs(mi2))**2, mi2=offset+1,offset+ncon)/)
            endif
        end do
        !
        do iloc=1,6
            offset = nos+nbs+(iloc-1)*ncon+5*ncon
            if((iloc.eq.3).or.(iloc.eq.6)) then
                probJ12 = probJ12 + (/(ABS(cs(mi2))**2, mi2=offset+1,offset+ncon)/)
            else
                probJ32 = probJ32 + (/(ABS(cs(mi2))**2, mi2=offset+1,offset+ncon)/)
            endif
        end do
        !
        do iloc=1,6
            offset = nos+nbs+(iloc-1)*ncon+11*ncon
            if((iloc.eq.3).or.(iloc.eq.6)) then
                probJ12 = probJ12 + (/(ABS(cs(mi2))**2, mi2=offset+1,offset+ncon)/)
            else
                probJ32 = probJ32 + (/(ABS(cs(mi2))**2, mi2=offset+1,offset+ncon)/)
            endif
        end do
        !
        write(nunit12, 1958) too, probJ12(:)
        write(nunit32, 1958) too, probJ32(:)
        1958 format(201(E17.8E4, 2x))
    end subroutine probs_from_cs
    !
    subroutine example_to(too, Eoo, ww, gaam, nunit,ip)
        integer :: nunit, ntt, iloc
        real*8 :: too, Eoo, ww, gaam
        real*8 :: ttlo, tthi, dtt
        logical :: ip
        real :: tst, tnd

        ttlo = too-5*gaam*SQRT(LOG(2.0d0))
        tthi = too+5*gaam*SQRT(LOG(2.0d0))
        dtt = 2*PI/ww * 1/100
        ntt = int((tthi-ttlo)/dtt)
        write(6,*) "A propagation on ", ntt, " time steps will be done"
        c = Zer
        do iloc=1, nos
            c(iloc) = A(iloc)*EXP(-Im*Es(iloc)*ttlo)
        end do
        write(6,*) "Time loop has started"
        call cpu_time(tst)
        do iloc=1, ntt
            call SPf_step(ntot, HD*Las(ttlo+iloc*dtt, dtt, too, Eoo, ww, gaam), V, dtt, c, caux, Es)
            c = caux
            if(ip.and.(nmod(iloc,10).eq.10)) then
                 write(nunit, 1986) ttlo+iloc*dtt, c(:)
            endif
        end do
        call cpu_time(tnd)
        write(6,*) "Propagation time ",tnd-tst," seconds"
        !call profile(tthi, ttlo, dtt, too, gaam)
        1986 format((E25.8E4,5x,5000(E20.8E4,1x,E20.8E4,3x)))
    end subroutine 
    !
    subroutine example_toJ1(too, Eoo, ww, gaam, nunit, ip)
        integer :: nunit, ntt, iloc
        real*8 :: too, Eoo, ww, gaam
        real*8 :: ttlo, tthi, dtt
        logical :: ip

        ttlo = too-5*gaam*SQRT(LOG(2.0d0))
        tthi = too+5*gaam*SQRT(LOG(2.0d0))
        dtt = 2*PI/ww * 1/100
        ntt = int((tthi-ttlo)/dtt)

        c1 = Zer
        do iloc=1, 4
            c1(iloc) = A(iloc)*EXP(-Im*Es1(iloc)*ttlo)
        end do
        write(6,*) "Time loop has started, with ", ntt, " steps."
        do iloc=1, ntt
            call SPf_step(ntot1, HD1*Las(ttlo+iloc*dtt, dtt, too, Eoo, ww, gaam), V1, dtt, c1, caux1, Es1)
            c1 = caux1
            if(ip.eqv..true..and.mod(iloc,10).eq.0) write(nunit, 1986) ttlo+iloc*dtt, c1(:)
        end do
        call profile(tthi, ttlo, dtt, too, gaam)
        1986 format((E25.8E4,5x,3500(E20.8E4,1x,E20.8E4,3x)))
    end subroutine 
    !
    subroutine no_delay_J3(Eoo, ww, gaam)
        integer :: ist
        real*8 :: Eoo, ww, gaam
        logical :: pf 
        pf = .true.
        open(17, file='SPF_J3_'//itoa(nl)//'/no_del_cs.dat')
        write(6,*) "Computing 0 delay coeffiecients"
        do ist=1,nos
            write(6,*) "State ", ist
            A = Zer
            A(ist) = Zone
            if(pf) open(44,file='SPF_J3_'//itoa(nl)//'/no_del_cs'//itoa(ist)//'.dat')
            call example_to(0.d0, Eoo, ww, gaam, 44, pf)
            close(44)
            csf(ist,:) = c(:)
            if(pf) write(17,1987) csf(ist,:)
        end do
        1987 format(5000(E20.8E4,1x,E20.8E4,3x))
    end subroutine no_delay_J3
    !
    subroutine computed_cs()
        integer :: ist, osq
        open(68, file="SPF_J3_"//itoa(nl)//"no_del_cs.dat")
        do ist=1,nos
            write(6,*) "Reading state ", ist 
            read(68,FMT=1987,IOSTAT=osq) csf(ist,:)
            write(6,*) csf(ist,1)
        end do
        close(68)
        1987 format(3500(E20.8E4,1x,E20.8E4,3x))
    end subroutine computed_cs
    !
    subroutine time_delay_loop()
        ! Run the time delay loop
        real*8 :: toolo, toofo, dtoo, too
        integer :: noo
        write(6, *) "Time delay loop"
        toolo = 5*gam*SQRT(LOG(2.0d0))
        toofo = toolo + 2.0d-12/auoftime
        noo = 700
        dtoo = 2.0d-12/auoftime * 1/noo

        open(78, file='SPF_J3_'//itoa(nl)//'/probJ12.dat')
        open(79, file='SPF_J3_'//itoa(nl)//'/probJ32.dat')

        do ii1=0, noo
            too = toolo+ii1*dtoo
            c = Zer
            do mi1=1,nos
                c = c + A(mi1)*EXP(-IM*Es(mi1)*too)*csf(mi1,:)
            end do
            ! Compute the probabilities and dump them in the files
            call probs_from_cs(c, too, 78, 79)
        end do

        close(78)
        close(79)
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
        ! write(6,*) ii        
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
   