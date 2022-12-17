program raman_p2up
    use j_recoup
    use errfun
    use arj1j2
    use matmod
    implicit none
    integer, parameter :: nos=4, nbs2=24, nbs0=13
    logical :: exp_flag
    integer :: mm, kk, nl, ii1, ii2, ii3
    real*8 :: Delt, Deltp
    real*8 :: Eo, w, to, gam, evec(3)
    real*8 :: Ar, Ai
    complex*16 :: A(nos)
    real*8 :: Zo(nos,5), Ebs(nos)
    real*8 :: ZJ0(nbs0,2), EnJ0(nbs0)
    real*8 :: ZJ2f(nbs2,3), EnJ2f(nbs2)
    real*8 :: ZJ2p(nbs2,3), EnJ2p(nbs2)
    real*8 :: RI_OI_J1_J2f(nos, nbs2, 3, 5), RI_OI_J1_J2p(nos, nbs2, 3,5)
    real*8 :: RI_OI_J1_J0p(nos, nbs0,2, 5)
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
        do mm=1,nos
            read(10,*) Ebs(mm), (Zo(mm,kk), kk=1,5)
            !write(6,*) Ebs(mm), (Zo(mm,kk), kk=1,5)
        end do
        close(10)

        open(10, file='interm_bounds_J2f.dat')
        do mm=1,nbs2
            read(10,*) EnJ2f(mm), (ZJ2f(mm,kk), kk=1,3)
            !write(6,*) EnJ2f(mm), (ZJ2f(mm,kk), kk=1,3)
        end do
        close(10)

        open(10, file='interm_bounds_J2p.dat')
        do mm=1,nbs2
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
    !
    
    !Loading the radial elements
        write(6,*) "Loading Radial Elements into four dimensional arrays"
        ! Radial elements from initial J=1 state, to intermediate J=2 state
        ! f states
        file_root = 'Python_Code_RadInts/Rad Ints bb for J=2 lint f/interm_rad_ints_'
        call Init_big_matsW(nos, nbs2, RI_OI_J1_J2f, 3, 5, file_root)
        ! p states
        !J=2
        file_root = 'Python_Code_RadInts/Rad Ints bb for J=2 lint p/interm_rad_ints_'
        call Init_big_matsW(nos, nbs2, RI_OI_J1_J2p, 3, 5, file_root)
        !J=0
        file_root = 'Python_Code_RadInts/Rad Ints bb for J=0 lint p/interm_rad_ints_'
        call Init_big_matsW(nos, nbs0, RI_OI_J1_J0p, 2, 5, file_root)
    
    !
    
            !Laser parameters
    open(98, file='laser_params_two_photon.dat', status='old')
    read(98,*) Eo, w, gam
    read(98,*) evec(:)
    read(98,*) nl
    close(98)
    print *, gam, Eo

    ! Initial state amplitudes
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
    ! Let's test the integrals:
    open(18,file='Raman_p2up/kernel_integrals.dat')
    ii1 = 3
    ii2 = 15
    Delt =  Deltxia(EnJ2f(ii2)/CMperAU,w,Ebs(ii1)/CMperAU,gam)
    Deltp = Deltxia(EnJ2f(ii2)/CMperAU,w,Ebs(1)/CMperAU,gam)
    write(18,*) Deltp, Delt, kernel(Deltp, Delt) 
    close(18)
    
    contains

    function Deltxia(Exi, ww, Ea, gaam)
        real*8 :: Exi, ww, Ea, gaam
        real*8 :: Deltxia
        Deltxia = 0.5d0*(Exi-Ea-ww)*gaam
    end function

    function kernel(delp, del)
        integer :: nx, ix, lnchf, ip
        real*8 :: dx, xlim, delp, del, ki, kr, xx
        real*8, allocatable :: xl(:), funci(:), funcr(:)
        complex*16 :: fi
        complex*16 :: kernel
        real*8, external :: rint
        complex*16, external :: CONHYP
        lnchf = 0
        ip = 0
        if(.not.(delp.eq.0.0d0)) then
            dx = 2*PI/(10*32*delp)
        else
            dx = 4.0d0/30.0d0
        endif

        xlim = 2
        nx = max(30,int(2*xlim/dx))
        allocate(xl(nx+1), funci(nx+1), funcr(nx+1))
        open(19,file='Raman_p2up/integrand.dat', status='replace')
        do ix=0, nx
            xl(ix+1) = -xlim+ix*dx
            xx = xl(ix+1)
            fi = EXP(-xx**2-2*IM*delp*xx)*(EXP(-del**2)+(2*(xx-IM*del))/SQRT(PI)*EXP(-xx**2+2*IM*del*xx)*&
            CONHYP(1.d0*ZONE,3.d0/2.d0*ZONE,(xx-IM*del)**2,LNCHF,IP))
            funcr(ix+1) = REAL(fi)
            funci(ix+1) = AIMAG(fi)
            write(19,*) (xx-IM*del)**2, CONHYP(1.d0*CMPLX(1.),3.d0/2.d0*CMPLX(1.),(xx-IM*del)**2,LNCHF,IP)
        end do
        close(19)
        ki = rint(funci, 1, nx, 10, dx)
        kr = rint(funcr, 1, nx, 20, dx)

        kernel = Zone*kr+Im*ki

    end function
    !
end program
!
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
   
