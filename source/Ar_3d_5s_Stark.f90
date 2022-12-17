program main
    use arj1j2
    implicit none
    integer, parameter :: nconts=200
    integer, parameter :: nos1=4, nos3=3, nos=nos1+nos3, nbs2 = 24, nbs4 = 26, nbs0=13, &
    nbs=nbs2+nbs2+nbs0+nbs4, nbs1 = 2*nbs2+nbs0, &
    ncs=5+6+6, ncs1 = 5+6, ne=200, ntot = nos+nbs+ncs*ne, ntot1 = nos1+nbs1+ncs1*ne
    logical :: exp_flag
    integer :: ncon, mm, kk
    integer :: nl, nlf
    integer :: mi, mi1, mi2, ii1
    real*8 :: Elo, dEcont, Econ(nconts)
    real*8 :: Z1(nos1, 5), EJ1(nos1)
    real*8 :: Z3(nos3,6), EJ3(nos3)
    real*8 :: ZJ0(nbs0,2), EnJ0(nbs0)
    real*8 :: ZJ2f(nbs2,3), EnJ2f(nbs2)
    real*8 :: ZJ2p(nbs2,3), EnJ2p(nbs2)
    real*8 :: ZJ4l(nbs4,6), EnJ4l(nbs4)

    real*8 :: rRI_CI_J0p_J1(nconts,nbs0,5,2)
    real*8 :: iRI_CI_J0p_J1(nconts,nbs0,5,2)
    real*8 :: rRI_CI_J2f_J1(nconts,nbs2,5,3), rRI_CI_J2f_J3(nconts,nbs2,6,3)
    real*8 :: iRI_CI_J2f_J1(nconts,nbs2,5,3), iRI_CI_J2f_J3(nconts,nbs2,6,3)
    real*8 :: rRI_CI_J2p_J1(nconts,nbs2,5,3), rRI_CI_J2p_J3(nconts,nbs2,6,3)
    real*8 :: iRI_CI_J2p_J1(nconts,nbs2,5,3), iRI_CI_J2p_J3(nconts,nbs2,6,3)
    real*8 :: rRI_CI_J4l_J3(nconts,nbs4,6,6), rRI_CI_J4l_J5(nconts,nbs4,6,6)
    real*8 :: iRI_CI_J4l_J3(nconts,nbs4,6,6), iRI_CI_J4l_J5(nconts,nbs4,6,6)
    
    real*8 :: RI_OI_J1_J2f(nos1, nbs2, 3, 5), RI_OI_J1_J2p(nos1, nbs2, 3,5)
    real*8 :: RI_OI_J2f_J1(nbs2, nos1, 5, 3), RI_OI_J2p_J1(nbs2, nos1, 5,3)
    
    real*8 :: RI_OI_J1_J0p(nos1, nbs0,2, 5)
    real*8 :: RI_OI_J0p_J1(nbs0, nos1,5, 2)
    
    real*8 :: RI_OI_J3_J2f(nos3, nbs2, 3, 6), RI_OI_J3_J2p(nos3, nbs2, 3, 6)
    real*8 :: RI_OI_J2f_J3(nbs2, nos3, 6, 3), RI_OI_J2p_J3(nbs2, nos3, 6, 3)
    
    real*8 :: RI_OI_J3_J4l(nos3, nbs4, 6, 6)
    real*8 :: RI_OI_J4l_J3(nbs4, nos3, 6, 6)

    real*8 :: Eo, w, gam, evec(3), ss(2), ssw
    real, allocatable :: jcp(:), jcm(:), jco(:)
    integer, allocatable :: jcsp(:),  lep(:), jcsm(:),  lem(:), jcso(:),  leo(:)
    character(10) :: nts
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
    do mi1=1,ncon
        Econ = (Elo+mi1*dEcont)/CMperAU
    end do
    !
    ! Compute the stark shifts for the initial states. This explicitely excludes the copuling with the ground state
    
    open(97, file="J1_starkshifts.dat")
    allocate(jcp(3), jcm(2), jco(5))
    allocate(jcsp(3),  lep(3), jcsm(2),  lem(2), jcso(5),  leo(5))
    jcp = (/3/2.,3/2.,1/2./)
    jcsp = (/2, 1, 1/)

    jcm = (/3/2., 1/2./)
    jcsm = (/1, 1/)
    lem = (/1, 1/)

    leo = lb
    jco = jcb
    jcso = jcsb

    do mi1=1,nos1
        ss = 0.d0
        lep = (/3,3,3/)
        ss = Stark_Shift_b(nbs2,3,ZJ2f,EnJ2f,jcp,jcsp,lep,2,RI_OI_J1_J2f(mi1,:,:,:),&
                            nbs0,2,ZJ0 ,EnJ0 ,jcm,jcsm,lem,0,RI_OI_J1_J0p(mi1,:,:,:),&
                            5,Z1(mi1,:),EJ1(mi1),jco, jcso, leo, 1,Eo,-1*w,gam)
        ssw = ss(1)+0.5d0*ss(2)                   
        lep = (/1,1,1/)
        ss = ss + Stark_Shift_b(nbs2,3,ZJ2p,EnJ2p,jcp,jcsp,lep,2,RI_OI_J1_J2p(mi1,:,:,:),&
                            nbs0,2,ZJ0 ,EnJ0 ,jcm,jcsm,lem,0,RI_OI_J1_J0p(mi1,:,:,:),&
                            5,Z1(mi1,:),EJ1(mi1),jco, jcso, leo, 1,Eo,-1*w,gam)
        ssw = ssw + ss(1)+0.5d0*ss(2)
        write(97,*) mi1, ssw
    enddo
    deallocate(jcp, jcm, jco)
    deallocate(jcsp,  lep, jcsm,  lem, jcso,  leo)
    close(97)

    open(97, file='J3_starkshifts.dat')
    allocate(jcp(6), jcm(3), jco(6))
    allocate(jcsp(6),  lep(6), jcsm(3),  lem(3), jcso(6),  leo(6))
    lep = (/3, 3, 3, 5, 5, 5/)
    jcp = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
    jcsp = (/1, 2, 1, 1, 2, 1/)

    jcm = (/ 3./2., 3./2., 1./2./)
    jcsm =  (/ 2, 1, 1 /)

    leo = (/ 2, 2, 2, 4, 4, 4/)
    jco = (/ 3/2., 3/2., 1/2., 3/2., 3/2., 1/2./)
    jcso = (/ 1, 2, 1, 1, 2, 1/)

    do mi1=1,nos3
        ss = 0.d0
        lem = (/3,3,3/)
        ss = Stark_Shift_b(nbs4,6,ZJ4l,EnJ4l,jcp,jcsp,lep,4,RI_OI_J3_J4l(mi1,:,:,:),&
                           nbs2,3,ZJ2f,EnJ2f,jcm,jcsm,lem,2,RI_OI_J3_J2f(mi1,:,:,:),&
                            6,Z3(mi1,:),EJ3(mi1),jco, jcso, leo, 3,Eo,-1*w,gam)
        ssw = 0.5d0*ss(1)+ss(2)                   
        lem = (/1,1,1/)
        ss = ss + Stark_Shift_b(nbs4,6,ZJ4l,EnJ4l,jcp,jcsp,lep,4,RI_OI_J3_J4l(mi1,:,:,:),&
                                nbs2,3,ZJ2p,EnJ2p,jcm,jcsm,lem,2,RI_OI_J3_J2p(mi1,:,:,:),&
                                6,Z3(mi1,:),EJ3(mi1),jco, jcso, leo, 3,Eo,-1*w,gam)
        ssw = ssw + 0.5d0*ss(1)+ss(2)
        write(97,*) mi1, ssw
    enddo
    deallocate(jcp, jcm, jco)
    deallocate(jcsp,  lep, jcsm,  lem, jcso,  leo)
    close(97)

    ! Now move on to the intermediate

    open(97, file='J0_starkshifts.dat')
    allocate(jcp(5), jcm(5), jco(2))
    allocate(jcsp(5),  lep(5), jcsm(5),  lem(5), jcso(2),  leo(2))
    
    leo = (/1, 1/)
    jco = (/3/2., 1/2./)
    jcso = (/1, 1/)

    lep = lb
    jcp = jcb
    jcsp = jcsb

    lem = lb
    jcm = jcb
    jcsm = jcsb
    
    ! Reorganize the radial integrals
    do mi1=1, nos1
        do mi2=1, nbs0
            RI_OI_J0p_J1(mi2,mi1,:,:) = TRANSPOSE(RI_OI_J1_J0p(mi1,mi2,:,:))
        end do
    end do

    do mi1=1,nbs0
        ss = 0.d0
        ss = Stark_Shift_b(nos1,5,Z1,EJ1,jcp,jcsp,lep,1,RI_OI_J0p_J1(mi1,:,:,:),&
                           nos1,5,Z1,EJ1,jcm,jcsm,lem,1,RI_OI_J0p_J1(mi1,:,:,:),&
                            2, ZJ0(mi1,:),EnJ0(mi1),jco, jcso, leo, 0,Eo, w, gam)
        ssw = 0.5d0*(ss(1)+ss(2))                   
        ss = Stark_Shift_c(ncon, 5, Econ, jcp, jcsp, lep, 1, rRI_CI_J0p_J1(:,mi1,:,:), iRI_CI_J0p_J1(:,mi1,:,:),&
                           ncon, 5, Econ, jcm, jcsm, lem, 1, rRI_CI_J0p_J1(:,mi1,:,:), iRI_CI_J0p_J1(:,mi1,:,:),&
                            dEcont/CMperAU, 2, ZJ0(mi1,:), EnJ0(mi1), jco, jcso, leo, 0, Eo, w, gam)
        write(97,*) mi1, ssw
    enddo
    deallocate(jcp, jcm, jco)
    deallocate(jcsp,  lep, jcsm,  lem, jcso,  leo)
    close(97)
    
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

        open(10, file='interm_bounds_J4l.dat')
        do mm=1,nbs4
            read(10,*) EnJ4l(mm), (ZJ4l(mm,kk), kk=1,6)
            !write(6,*) EnJ4l(mm), (ZJ4l(mm,kk), kk=1,2)
        end do
        close(10)
    end subroutine
    !
    subroutine load_rads()
        character(len=143) :: file_root
        write(6,*) "Loading Radial Elements into four dimensional arrays"
        ! Radial elements from initial J=1 state, to intermediate J=2 state
        ! f states
        file_root = 'Python_Code_RadInts/RI b(J=1)  b(J=2 lint f)/interm_rad_ints_'
        call Init_big_matsW(nos1, nbs2, RI_OI_J1_J2f, 3, 5, file_root)
        ! p states
        file_root = 'Python_Code_RadInts/RI b(J=1)  b(J=2 lint p)/interm_rad_ints_'
        call Init_big_matsW(nos1, nbs2, RI_OI_J1_J2p, 3, 5, file_root)
        !J=0
        file_root = 'Python_Code_RadInts/RI b(J=1)  b(J=0 lint p)/interm_rad_ints_'
        call Init_big_matsW(nos1, nbs0, RI_OI_J1_J0p, 2, 5, file_root)

        ! Radial elements from initial J=3 state, to intermediate J=2 state
        ! f states
        file_root = 'Python_Code_RadInts/RI b(J=3)  b(J=2 lint f)/interm_rad_ints_'
        call Init_big_matsW(nos3, nbs2, RI_OI_J3_J2f, 3, 6, file_root)
        ! p states
        file_root = 'Python_Code_RadInts/RI b(J=3)  b(J=2 lint p)/interm_rad_ints_'
        call Init_big_matsW(nos3, nbs2, RI_OI_J3_J2p, 3, 6, file_root)
        !J=4
        file_root = 'Python_Code_RadInts/RI b(J=3)  b(J=4 lint f)/interm_rad_ints_'
        call Init_big_matsW(nos3, nbs4, RI_OI_J3_J4l, 6, 6, file_root)

        ! Radial elements to J=1 continuum, from an intermediate J=2 state
        ! f states
        file_root = 'Python_Code_RadInts/RI b(J=2,l=f) c(J=1)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs2, nconts, rRI_CI_J2f_J1, 5, 3, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=2,l=f) c(J=1)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs2, nconts, iRI_CI_J2f_J1, 5, 3, file_root)
        ! p states
        file_root = 'Python_Code_RadInts/RI b(J=2,l=p) c(J=1)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs2, nconts, rRI_CI_J2p_J1, 5, 3, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=2,l=p) c(J=1)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs2, nconts, iRI_CI_J2p_J1, 5, 3, file_root)
        ! Radial elements to J=1 continuum, from an intermediate J=0
        file_root = 'Python_Code_RadInts/RI b(J=0,l=p) c(J=1)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs0, nconts, rRI_CI_J0p_J1, 5, 2, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=0,l=p) c(J=1)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs0, nconts, iRI_CI_J0p_J1, 5, 2, file_root)

        ! Radial elements to J=3 continuum, from an intermediate J=2 state
        ! f states
        file_root = 'Python_Code_RadInts/RI b(J=2,l=f) c(J=3)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs2, nconts, rRI_CI_J2f_J3, 6, 3, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=2,l=f) c(J=3)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs2, nconts, iRI_CI_J2f_J3, 6, 3, file_root)
        ! p states
        !J=2
        file_root = 'Python_Code_RadInts/RI b(J=2,l=p) c(J=3)/cont_rad_ints_reg_b'
        call Init_big_mats(nbs2, nconts, rRI_CI_J2p_J3, 6, 3, file_root)
        file_root = 'Python_Code_RadInts/RI b(J=2,l=p) c(J=3)/cont_rad_ints_ireg_b'
        call Init_big_mats(nbs2, nconts, iRI_CI_J2p_J3, 6, 3, file_root)
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