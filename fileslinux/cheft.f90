! Alex, November 14, 2013

! Corrected version of the chiral potentials, 
! also containing cutoffs at r0 = 0.9 and 1.1 fm.
! The CIB and CSB terms which were previously mixing
! Evgeny and Pudliner conventions, have been changed.

! SFR cutoff = 1000 MeV

! To be determined: 
! if we should bother including a (1-sigma1.sigma2)/4 term in CIB
! and CSB, the latter of which would require more than 18 operators

! Note that:
! * this contains the same (modified due to the Goldberger-Treiman relation)
!   one-pion exchange potential at all orders, so that only intermediate 
!   and short-range behavior changes from order to order
! * Now CIB and CSB are present at all orders (with different parameters
!   at each order) in order to always produce decent scattering lengths

module cheft

    implicit none

    integer, private, parameter :: i4=selected_int_kind(9)
    integer, private, parameter :: r8=selected_real_kind(15,9)

    real(kind=r8), private, parameter :: m_pi=138.03_r8
    real(kind=r8), private, parameter :: m_pi_not=134.98_r8
    real(kind=r8), private, parameter :: m_pi_plus=139.570_r8
    real(kind=r8), private, parameter :: m_n=939.565_r8
    real(kind=r8), private, parameter :: g_a=1.267_r8
    real(kind=r8), private, parameter :: g_a_Tre=1.29_r8
    real(kind=r8), private, parameter :: f_pi=92.4_r8
    real(kind=r8), private, parameter :: pi=3.1415926535897_r8
    real(kind=r8), private, parameter :: hbarc=197.327_r8
    real(kind=r8), private, parameter :: g34 = 1.22541670246517_r8

    real(kind=r8), private, parameter :: EGM_C1=-0.00081_r8
    real(kind=r8), private, parameter :: EGM_C3=-0.0034_r8
    real(kind=r8), private, parameter :: EGM_C4=0.0034_r8

    real(kind=r8), private, parameter :: Lam=1000.0_r8

    !r0 = 0.9
!   real(kind=r8), private, parameter :: r_0p=0.9_r8
!   real(kind=r8), private, parameter :: LO_C_s=-0.128977_r8
!   real(kind=r8), private, parameter :: LO_C_t=0.518088_r8
!   real(kind=r8), private, parameter :: LO_CIB=-0.026542_r8
!   real(kind=r8), private, parameter :: LO_CSB=-0.019226_r8
!   real(kind=r8), private, parameter :: NLO_C_s=9.550477_r8
!   real(kind=r8), private, parameter :: NLO_C_t=3.148778_r8
!   real(kind=r8), private, dimension(7) :: NLO_C=(/0.370853_r8,&
!   0.363613_r8,-0.184321_r8,0.218210_r8,-3.003296_r8,0.543729_r8,&
!   -0.502584_r8/)
!   real(kind=r8), private, parameter :: NLO_CIB=0.058030_r8
!   real(kind=r8), private, parameter :: NLO_CSB=0.011477_r8
!   real(kind=r8), private, parameter :: N2LO_C_s=7.747830_r8
!   real(kind=r8), private, parameter :: N2LO_C_t=0.452466_r8
!   real(kind=r8), private, dimension(7) :: N2LO_C=(/-0.217154_r8,&
!   0.034570_r8,-0.115352_r8,0.118178_r8,-2.416024_r8,0.154628_r8,&
!   -0.267084_r8/)
!   real(kind=r8), private, parameter :: N2LO_CIB=0.052466_r8
!   real(kind=r8), private, parameter :: N2LO_CSB=0.010824_r8

     !r0 = 1.0
!    real(kind=r8), private, parameter :: r_0p=1.0_r8
!    real(kind=r8), private, parameter :: LO_C_s=-0.751120_r8
!    real(kind=r8), private, parameter :: LO_C_t=0.374089_r8
!    real(kind=r8), private, parameter :: LO_CIB=-0.023605_r8
!    real(kind=r8), private, parameter :: LO_CSB=-0.019876_r8
!    real(kind=r8), private, parameter :: NLO_C_s=3.168026_r8
!    real(kind=r8), private, parameter :: NLO_C_t=1.413952_r8
!    real(kind=r8), private, dimension(7) :: NLO_C=(/0.314202_r8,&
!    0.257857_r8,-0.131344_r8,0.118613_r8,-2.385514_r8,0.373188_r8,&
!    -0.356684_r8/)
!    real(kind=r8), private, parameter :: NLO_CIB=0.050944_r8
!    real(kind=r8), private, parameter :: NLO_CSB=0.008231_r8
!    real(kind=r8), private, parameter :: N2LO_C_s=5.438495_r8
!    real(kind=r8), private, parameter :: N2LO_C_t=0.276718_r8
!    real(kind=r8), private, dimension(7) :: N2LO_C=(/-0.140842_r8,&
!    0.042427_r8,-0.123378_r8,0.110184_r8,-2.112533_r8,0.158979_r8,&
!    -0.269935_r8/)
!    real(kind=r8), private, parameter :: N2LO_CIB=0.053204_r8
!    real(kind=r8), private, parameter :: N2LO_CSB=0.009759_r8

     !r0 = 1.1
!    real(kind=r8), private, parameter :: r_0p=1.1_r8
!    real(kind=r8), private, parameter :: LO_C_s=-1.296304_r8
!    real(kind=r8), private, parameter :: LO_C_t=0.256476_r8
!    real(kind=r8), private, parameter :: LO_CIB=-0.019217_r8
!    real(kind=r8), private, parameter :: LO_CSB=-0.020007_r8
!    real(kind=r8), private, parameter :: NLO_C_s=1.030748_r8
!    real(kind=r8), private, parameter :: NLO_C_t=0.906986_r8
!    real(kind=r8), private, dimension(7) :: NLO_C=(/0.272389_r8,&
!    0.220322_r8,-0.136407_r8,0.094205_r8,-2.162377_r8,0.330646_r8,&
!    -0.335697_r8/)
!    real(kind=r8), private, parameter :: NLO_CIB=0.051529_r8
!    real(kind=r8), private, parameter :: NLO_CSB=0.007041_r8
!    real(kind=r8), private, parameter :: N2LO_C_s=3.886979_r8
!    real(kind=r8), private, parameter :: N2LO_C_t=0.244154_r8
!    real(kind=r8), private, dimension(7) :: N2LO_C=(/-0.096496_r8,&
!    0.059465_r8,-0.141827_r8,0.111460_r8,-2.008237_r8,0.183174_r8,&
!    -0.301046_r8/)
!    real(kind=r8), private, parameter :: N2LO_CIB=0.055377_r8
!    real(kind=r8), private, parameter :: N2LO_CSB=0.009018_r8

     !r0 = 1.2
!    real(kind=r8), private, parameter :: r_0p=1.2_r8
!    real(kind=r8), private, parameter :: LO_C_s=-1.796928_r8
!    real(kind=r8), private, parameter :: LO_C_t=0.154414_r8
!    real(kind=r8), private, parameter :: LO_CIB=-0.013348_r8
!    real(kind=r8), private, parameter :: LO_CSB=-0.019589_r8
!    real(kind=r8), private, parameter :: NLO_C_s=0.035507_r8
!    real(kind=r8), private, parameter :: NLO_C_t=0.717286_r8
!    real(kind=r8), private, dimension(7) :: NLO_C=(/0.222888_r8,&
!    0.228779_r8,-0.150425_r8,0.089285_r8,-2.029313_r8,0.340112_r8,&
!    -0.362474_r8/)
!    real(kind=r8), private, parameter :: NLO_CIB=0.054768_r8
!    real(kind=r8), private, parameter :: NLO_CSB=0.006596_r8
!    real(kind=r8), private, parameter :: N2LO_C_s=2.687641_r8
!    real(kind=r8), private, parameter :: N2LO_C_t=0.233817_r8
!    real(kind=r8), private, dimension(7) :: N2LO_C=(/-0.079513_r8,&
!    0.076102_r8,-0.169261_r8,0.123588_r8,-1.942800_r8,0.214206_r8,&
!    -0.341926_r8/)
!    real(kind=r8), private, parameter :: N2LO_CIB=0.056482_r8
!    real(kind=r8), private, parameter :: N2LO_CSB=0.007710_r8

     real(kind=r8), private, save :: r_0p
     real(kind=r8), private, save :: LO_C_s
     real(kind=r8), private, save :: LO_C_t
     real(kind=r8), private, save :: LO_CIB
     real(kind=r8), private, save :: LO_CSB
     real(kind=r8), private, save :: NLO_C_s
     real(kind=r8), private, save :: NLO_C_t
     real(kind=r8), private, save, dimension(7) :: NLO_C
     real(kind=r8), private, save :: NLO_CIB
     real(kind=r8), private, save :: NLO_CSB
     real(kind=r8), private, save :: N2LO_C_s
     real(kind=r8), private, save :: N2LO_C_t
     real(kind=r8), private, save, dimension(7) :: N2LO_C
     real(kind=r8), private, save :: N2LO_CIB
     real(kind=r8), private, save :: N2LO_CSB

contains

    subroutine cheft_pot(lpot,rr,vnn)
    integer :: lpot
    real(kind=r8) :: rr
    real(kind=r8), dimension(18) :: vnn

    real(kind=r8) :: r
    real(kind=r8) :: x_not, x_plus, g_r0p, del_r
    real(kind=r8) :: pre_not, pre_plus
    real(kind=r8) :: W_C, V_S, V_T
    real(kind=r8) :: V_C, W_S, W_T
    real(kind=r8) :: lap_res, other_res
    real(kind=r8) :: sp_not, sp_plus
    real(kind=r8) :: ten_not, ten_plus

    ! Hi-jacking Bob's lpot to describe chiral EFT orders:
    ! lpot=112 means LO, R0=1.0
    ! lpot=113 means NLO, R0=1,0
    ! lpot=114 means N2LO, R0=1.0
    ! lpot=122 means LO, R0=1.2
    ! lpot=123 means NLO, R0=1,2
    ! lpot=124 means N2LO, R0=1.2

    select case(lpot)
       case(112:114)
          r_0p=1.0_r8
          LO_C_s=-0.751120_r8
          LO_C_t=0.374089_r8
          LO_CIB=-0.023605_r8
          LO_CSB=-0.019876_r8
          NLO_C_s=3.168026_r8
          NLO_C_t=1.413952_r8
          NLO_C=(/0.314202_r8,0.257857_r8,-0.131344_r8,0.118613_r8,-2.385514_r8,0.373188_r8,-0.356684_r8/)
          NLO_CIB=0.050944_r8
          NLO_CSB=0.008231_r8
          N2LO_C_s=5.438495_r8
          N2LO_C_t=0.276718_r8
          N2LO_C=(/-0.140842_r8,0.042427_r8,-0.123378_r8,0.110184_r8,-2.112533_r8,0.158979_r8,-0.269935_r8/)
          N2LO_CIB=0.053204_r8
          N2LO_CSB=0.009759_r8
       case(122:124)
          r_0p=1.2_r8
          LO_C_s=-1.796928_r8
          LO_C_t=0.154414_r8
          LO_CIB=-0.013348_r8
          LO_CSB=-0.019589_r8
          NLO_C_s=0.035507_r8
          NLO_C_t=0.717286_r8
          NLO_C=(/0.222888_r8,0.228779_r8,-0.150425_r8,0.089285_r8,-2.029313_r8,0.340112_r8,-0.362474_r8/)
          NLO_CIB=0.054768_r8
          NLO_CSB=0.006596_r8
          N2LO_C_s=2.687641_r8
          N2LO_C_t=0.233817_r8
          N2LO_C=(/-0.079513_r8,0.076102_r8,-0.169261_r8,0.123588_r8,-1.942800_r8,0.214206_r8,-0.341926_r8/)
          N2LO_CIB=0.056482_r8
          N2LO_CSB=0.007710_r8
    end select

    !Sometimes the input rr is 0, so I need to make sure
    !this routine handles that properly:
    if (abs(rr).LT.0.001_r8) then
        r = 0.001_r8
    else
        r = rr
    end if

    vnn(1:18) = 0.0_r8

    g_r0p = 1.0_r8/(pi*g34*r_0p**3)

    !Note that the following is used both for contacts and to regulate OPE
    del_r = exp(-(r/r_0p)**4)

    !Also note that contact parameters are always multiplied by hbarc

    !We take OPE as Goldberger-Treiman modified, i.e. the same, at all orders
    !x = r*mu = r*m_pi*c/hbar
    x_not = m_pi_not*r/hbarc
    pre_not = (m_pi_not**3/(12.0_r8*pi))*(g_a_Tre/(2*f_pi))**2*&
                  exp(-x_not)/x_not 
    x_plus = m_pi_plus*r/hbarc
    pre_plus = (m_pi_plus**3/(12.0_r8*pi))*(g_a_Tre/(2*f_pi))**2*&
                   exp(-x_plus)/x_plus 

    !The pion masses are according to Eq. (2.7) in the 1997 Pudliner PRC
    sp_not = (1.0_r8-del_r)*pre_not
    sp_plus = (1.0_r8-del_r)*pre_plus
    vnn(4) = (sp_not + 2.0_r8*sp_plus)/3.0_r8
    ten_not = (1.0_r8-del_r)*pre_not*(1.0_r8 + 3.0_r8/x_not + &
              3.0_r8/x_not**2)
    ten_plus = (1.0_r8-del_r)*pre_plus*(1.0_r8 + 3.0_r8/x_plus + &
               3.0_r8/x_plus**2)
    vnn(6) = (ten_not + 2.0_r8*ten_plus)/3.0_r8
    vnn(16) = (sp_not - sp_plus)/3.0_r8
    vnn(17) = (ten_not - ten_plus)/3.0_r8


    if (lpot.EQ.112.or.lpot.eq.122) then
        !We are including charge-breaking at each order
        vnn(1) = vnn(1) + hbarc*(LO_CIB/2.0_r8)*del_r*g_r0p
        vnn(2) = vnn(2) + hbarc*(LO_CIB/6.0_r8)*del_r*g_r0p
        vnn(15) = vnn(15) + hbarc*(LO_CIB/6.0_r8)*del_r*g_r0p
        vnn(18) = vnn(18) + hbarc*(LO_CSB/2.0_r8)*del_r*g_r0p

        !LO has 2 contacts
        vnn(1) = vnn(1) + hbarc*LO_C_s*del_r*g_r0p
        vnn(3) = vnn(3) + hbarc*LO_C_t*del_r*g_r0p
    end if

    if (lpot.EQ.113.or.lpot.eq.123) then
        !We are including charge-breaking at each order
        vnn(1) = vnn(1) + hbarc*(NLO_CIB/2.0_r8)*del_r*g_r0p
        vnn(2) = vnn(2) + hbarc*(NLO_CIB/6.0_r8)*del_r*g_r0p
        vnn(15) = vnn(15) + hbarc*(NLO_CIB/6.0_r8)*del_r*g_r0p
        vnn(18) = vnn(18) + hbarc*(NLO_CSB/2.0_r8)*del_r*g_r0p

        !At NLO the LO 2 contacts change
        vnn(1) = vnn(1) + hbarc*NLO_C_s*del_r*g_r0p
        vnn(3) = vnn(3) + hbarc*NLO_C_t*del_r*g_r0p

        !NLO also has a TPE
        call spectralNLO(r,W_C,V_S,V_T,Lam)
        vnn(2) = vnn(2) + (1.0_r8-del_r)*W_C
        vnn(3) = vnn(3) + (1.0_r8-del_r)*V_S 
        vnn(5) = vnn(5) + (1.0_r8-del_r)*V_T 

        !NLO also has 7 new contacts
        lap_res = 20.0_r8*(r**2/r_0p**4) - 16.0_r8*(r**6/r_0p**8)
        vnn(1) = vnn(1) + hbarc*NLO_C(1)*lap_res*del_r*g_r0p
        vnn(2) = vnn(2) + hbarc*NLO_C(2)*lap_res*del_r*g_r0p
        vnn(3) = vnn(3) + hbarc*NLO_C(3)*lap_res*del_r*g_r0p 
        vnn(4) = vnn(4) + hbarc*NLO_C(4)*lap_res*del_r*g_r0p 
        vnn(3) = vnn(3) + hbarc*NLO_C(6)*(lap_res/3.0_r8)*del_r*g_r0p 
        vnn(4) = vnn(4) + hbarc*NLO_C(7)*(lap_res/3.0_r8)*del_r*g_r0p 
        other_res = (8.0_r8*(r**2/r_0p**4) - 16.0_r8*(r**6/r_0p**8))/3.0_r8
        vnn(5) = vnn(5) + hbarc*NLO_C(6)*other_res*del_r*g_r0p 
        vnn(6) = vnn(6) + hbarc*NLO_C(7)*other_res*del_r*g_r0p 
        vnn(7) = hbarc*NLO_C(5)*2.0_r8*(r**2/r_0p**4)*del_r*g_r0p 
    end if

    if (lpot.EQ.114.or.lpot.eq.124) then
        !We are including charge-breaking at each order
        vnn(1) = vnn(1) + hbarc*(N2LO_CIB/2.0_r8)*del_r*g_r0p
        vnn(2) = vnn(2) + hbarc*(N2LO_CIB/6.0_r8)*del_r*g_r0p
        vnn(15) = vnn(15) + hbarc*(N2LO_CIB/6.0_r8)*del_r*g_r0p
        vnn(18) = vnn(18) + hbarc*(N2LO_CSB/2.0_r8)*del_r*g_r0p

        !At N2LO the LO 2 contacts change
        vnn(1) = vnn(1) + hbarc*N2LO_C_s*del_r*g_r0p
        vnn(3) = vnn(3) + hbarc*N2LO_C_t*del_r*g_r0p

        !N2LO also has NLO's TPE
        call spectralNLO(r,W_C,V_S,V_T,Lam)
        vnn(2) = vnn(2) + (1.0_r8-del_r)*W_C
        vnn(3) = vnn(3) + (1.0_r8-del_r)*V_S 
        vnn(5) = vnn(5) + (1.0_r8-del_r)*V_T 

        !N2LO also has a different TPE
        call spectralN2LO(r,V_C,W_S,W_T,Lam)
        vnn(1) = vnn(1) + (1.0_r8-del_r)*V_C
        vnn(4) = vnn(4) + (1.0_r8-del_r)*W_S 
        vnn(6) = vnn(6) + (1.0_r8-del_r)*W_T 
        
        !N2LO also has new parameters for the 7 NLO contacts
        lap_res = 20.0_r8*(r**2/r_0p**4) - 16.0_r8*(r**6/r_0p**8)
        vnn(1) = vnn(1) + hbarc*N2LO_C(1)*lap_res*del_r*g_r0p
        vnn(2) = vnn(2) + hbarc*N2LO_C(2)*lap_res*del_r*g_r0p
        vnn(3) = vnn(3) + hbarc*N2LO_C(3)*lap_res*del_r*g_r0p 
        vnn(4) = vnn(4) + hbarc*N2LO_C(4)*lap_res*del_r*g_r0p 
        vnn(3) = vnn(3) + hbarc*N2LO_C(6)*(lap_res/3.0_r8)*del_r*g_r0p 
        vnn(4) = vnn(4) + hbarc*N2LO_C(7)*(lap_res/3.0_r8)*del_r*g_r0p 
        other_res = (8.0_r8*(r**2/r_0p**4) - 16.0_r8*(r**6/r_0p**8))/3.0_r8
        vnn(5) = vnn(5) + hbarc*N2LO_C(6)*other_res*del_r*g_r0p 
        vnn(6) = vnn(6) + hbarc*N2LO_C(7)*other_res*del_r*g_r0p 
        vnn(7) = hbarc*N2LO_C(5)*2.0_r8*(r**2/r_0p**4)*del_r*g_r0p 
    end if

    return
    end subroutine cheft_pot

!----------------------------------------------------------------------!

    SUBROUTINE gauleg(x1,x2,x,w,n)
    INTEGER n
    real(kind=r8) ::  x1,x2,x(n),w(n)
    real(kind=r8) ::  EPS
    PARAMETER (EPS=1.d-14)
    INTEGER i,j,m
    real(kind=r8) ::  p1,p2,p3,pp,xl,xm,z,z1
    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
    do 12 i=1,m
      z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1     continue
        p1=1.d0
        p2=0.d0
        do 11 j=1,n
          p3=p2
          p2=p1
          p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11      continue
        pp=n*(z*p1-p2)/(z*z-1.d0)
        z1=z
        z=z1-p1/pp
      if(abs(z-z1).gt.EPS)goto 1
      x(i)=xm-xl*z
      x(n+1-i)=xm+xl*z
      w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
      w(n+1-i)=w(i)
12  continue
    return
    END SUBROUTINE gauleg

!----------------------------------------------------------------------!

    subroutine spectralNLO(r,W_C,V_S,V_T,L_cut)
    
    real(kind=r8) :: r, W_C, V_S, V_T, L_cut
    integer, parameter :: n_int=50

    integer :: i
    real(kind=r8), dimension(n_int) :: Xgauleg, Wgauleg
    real(kind=r8) :: f_pi_o_fm, m_pi_o_fm, L_cut_o_fm
    real(kind=r8) :: fmpisq, overall, co1, co2, co3, res
    real(kind=r8) :: a, b, x

    f_pi_o_fm=f_pi/hbarc
    m_pi_o_fm=m_pi/hbarc
    fmpisq = 4.0*(m_pi_o_fm**2)
    L_cut_o_fm=L_cut/hbarc

    a = 2.0*m_pi_o_fm
    b = L_cut_o_fm

    call gauleg(a,b,Xgauleg,Wgauleg,n_int)

    W_C = 0.0_r8
    V_S = 0.0_r8
    V_T = 0.0_r8

    do i = 1, n_int
        x = Xgauleg(i)

        overall = 1.0/(768.0*pi*(f_pi_o_fm**4))
        co1 = 5.0*g_a**4 - 4.0*g_a**2 - 1.0
        co2 = 23.0*g_a**4 - 10.0*g_a**2 - 1.0
        co3 = 48.0*(g_a**4)*(m_pi_o_fm**4)
        res = exp(-x*r)*sqrt(x**2 - fmpisq)
        res = res*(fmpisq*co1 - (x**2)*co2 + co3/(fmpisq - x**2))
        res = res*overall
        res = (hbarc*res)/(2.0*(pi**2)*r)
        W_C = W_C + Wgauleg(i)*res

        overall = 3.0*(g_a**4)/(128.0*pi*(f_pi_o_fm**4))
        res = exp(-x*r)*sqrt(x**2 - fmpisq)
        res = res*(x**2)*overall
        res = (hbarc*res)/(3.0*(pi**2)*r)
        V_S = V_S + Wgauleg(i)*res

        overall = 3.0*(g_a**4)/(128.0*pi*(f_pi_o_fm**4))
        res = exp(-x*r)*sqrt(x**2 - fmpisq)
        res = res*(3.0 + 3.0*x*r + (x**2)*(r**2))*overall
        res = -(hbarc*res)/(6.0*(pi**2)*(r**3))
        V_T = V_T + Wgauleg(i)*res

    end do

    end subroutine spectralNLO

!----------------------------------------------------------------------!

    subroutine spectralN2LO(r,V_C,W_S,W_T,L_cut)
    
    real(kind=r8) :: r, V_C, W_S, W_T, L_cut

    real(kind=r8) :: f_pi_o_fm, m_pi_o_fm, L_cut_o_fm
    real(kind=r8) :: fpifour, overall, co1, co2
    real(kind=r8) :: x, y

    f_pi_o_fm=f_pi/hbarc
    m_pi_o_fm=m_pi/hbarc
    fpifour = f_pi_o_fm**4
    L_cut_o_fm=L_cut/hbarc

    x = m_pi_o_fm*r
    y = L_cut_o_fm*r

    overall = (3.0*g_a**2)/(32.0*(pi**2)*fpifour)
    co1 = 2.0*EGM_C1*x**2*(1.0+x)**2 &
          + EGM_C3*(6.0+12.0*x+10.0*x**2+4.0*x**3+x**4)
    co1 = co1*exp(-2.0*x)/(r**6)
    co2 = 4.0*EGM_C1*x**2*(2.0+y*(2.0+y)-2.0*x**2) &
          + EGM_C3*(24.0 + y*(24.0+12.0*y+4.0*y**2+y**3) &
                - 4.0*x**2*(2.0+2.0*y+y**2) + 4.0*x**4)
    co2 = -co2*exp(-y)/(4.0*r**6)
    V_C = overall*(co1 + co2)
    V_C = V_C*(hbarc**2)

    overall = (g_a**2)/(48.0*(pi**2)*fpifour)
    co1 = EGM_C4*(1.0+x)*(3.0+3.0*x+2.0*x**2)
    co1 = co1*exp(-2.0*x)/(r**6)
    co2 = EGM_C4*(24.0+24.0*y+12.0*y**2+4.0*y**3+y**4 &
          - 4.0*x**2*(2.0+2.0*y+y**2))
    co2 = -co2*exp(-y)/(8.0*r**6)
    W_S = overall*(co1 + co2)
    W_S = W_S*(hbarc**2)

    overall = (g_a**2)/(48.0*(pi**2)*fpifour)
    co1 = EGM_C4*(1.0+x)*(3.0+3.0*x+x**2)
    co1 = -co1*exp(-2.0*x)/(r**6)
    co2 = EGM_C4*(48.0+48.0*y+24.0*y**2+7.0*y**3+y**4 &
          - 4.0*x**2*(8.0+5.0*y+y**2))
    co2 = co2*exp(-y)/(16.0*r**6)
    W_T = overall*(co1 + co2)
    W_T = W_T*(hbarc**2)

    end subroutine spectralN2LO

!----------------------------------------------------------------------!
!                                                                      !
!  A partial-wave projection for V8-like potentials (very much like    !
!  Bob's); now modified to include partial-wave projections for        !
!  operators 15 through 18 (in the AV18-classification scheme).        !
!  Calls subroutines cheft_pot.                                        !
!----------------------------------------------------------------------!
!  arguments for cheft_fullpw                                          !
!  lpot: As in Alex's code above.                                      !
!  l:    orbital angular momentum of pair (0,1,2,...)                  !
!  s:    total spin of pair (0 or 1)                                   !
!  j:    total angular momentum of pair (0,1,2,...)                    !
!  t:    total isospin of pair (0 or 1)                                !
!  t1z:  isospin of particle 1 (1 for p, -1 for n)                     !
!  t2z:  isospin of particle 2 (1 for p, -1 for n)                     !
!  rr:   separation in fm                                              !
!  vpw:  returned potential in MeV (2x2 array)                         !
!        (includes all strong and NO em terms)                         !
!----------------------------------------------------------------------!
!  order of terms in v(l,m):                                           !
!       single channel                 coupled channel (l=j-1,s=1)     !
!       v(1,1) = v(l,s,j,t,t1z,t2z)    v(1,1) = v(l,s,j,t,t1z,t2z)     !
!       v(2,1) = 0                     v(2,1) = v(l<->l+2)             !
!       v(1,2) = 0                     v(1,2) = v(l<->l+2)             !
!       v(2,2) = 0                     v(2,2) = v(l+2,s,j,t,t1z,t2z)   !
!----------------------------------------------------------------------!
   subroutine cheft_fullpw(lpot,l,s,j,t,t1z,t2z,rr,vpw)

   implicit none

   integer(kind=i4),intent(in) :: lpot,l,s,j,t,t1z,t2z
   real(kind=r8),intent(in) :: rr
   real(kind=r8),intent(out) :: vpw(2,2)

   integer(kind=i4) :: ncc,ls,lsm,lsp
   real(kind=r8) :: vnn(18),vc,vt,vls,s12,s12m,s12p,s1ds2,t1dt2,t12
                                                                                                
   call cheft_pot(lpot,rr,vnn)

   s1ds2=4*s-3
   t1dt2=4*t-3
   t12=3*t1z*t2z-t1dt2
   vc=vnn(1)+t1dt2*vnn(2)+s1ds2*vnn(3)+s1ds2*t1dt2*vnn(4)+&
   &t12*vnn(15)+s1ds2*t12*vnn(16)+(t1z+t2z)*vnn(18)
   vt=vnn(5)+t1dt2*vnn(6)+t12*vnn(17)
   vls=vnn(7)+t1dt2*vnn(8)

   ncc=1

   if(s.eq.1.and.j.gt.l)then
      ncc=2
   endif

   if(ncc.eq.1)then
      s12=0.0_r8

      if(s.eq.1.and.l.eq.j)then
         s12=2.0_r8
      endif

      if(l.eq.(j+1))then
         s12=-2.0_r8*(j+2)/(2.0_r8*j+1.0_r8)
      endif

      ls=(j*(j+1)-l*(l+1)-s*(s+1))/2

      vpw(1,1)=vc+s12*vt+ls*vls
      vpw(2,1)=0
      vpw(1,2)=0
      vpw(2,2)=0

   elseif(ncc.eq.2)then

      s12m=-2.0_r8*(j-1.0_r8)/(2.0_r8*j+1.0_r8)
      s12=sqrt(36.0_r8*j*(j+1))/(2.0_r8*j+1.0_r8)
      s12p=-2.0_r8*(j+2.0_r8)/(2.0_r8*j+1.0_r8)
      lsm=j-1
      lsp=-(j+2)
      vpw(1,1)=vc+s12m*vt+lsm*vls
      vpw(2,1)=s12*vt
      vpw(1,2)=s12*vt
      vpw(2,2)=vc+s12p*vt+lsp*vls

   endif

   return

   endsubroutine cheft_fullpw

!----------------------------------------------------------------------!
!                                                                      !
!  A partial-wave projection for V8-like potentials (very much like    !
!  Bob's).                                                             !
!  Calls subroutines cheft_pot.                                        !
!----------------------------------------------------------------------!
!  arguments for cheft_pw                                              !
!  lpot: As in Alex's code above.                                      !
!  l:    orbital angular momentum of pair (0,1,2,...)                  !
!  s:    total spin of pair (0 or 1)                                   !
!  j:    total angular momentum of pair (0,1,2,...)                    !
!  t:    total isospin of pair (0 or 1)                                !
!  rr:   separation in fm                                              !
!  vpw:  returned potential in MeV (2x2 array)                         !
!        (includes all strong and NO em terms)                         !
!----------------------------------------------------------------------!
!  order of terms in v(l,m):                                           !
!       single channel                 coupled channel (l=j-1,s=1)     !
!       v(1,1) = v(l,s,j,t)            v(1,1) = v(l,s,j,t)             !
!       v(2,1) = 0                     v(2,1) = v(l<->l+2)             !
!       v(1,2) = 0                     v(1,2) = v(l<->l+2)             !
!       v(2,2) = 0                     v(2,2) = v(l+2,s,j,t)           !
!----------------------------------------------------------------------!
   subroutine cheft_pw(lpot,l,s,j,t,rr,vpw)

   implicit none

   integer(kind=i4),intent(in) :: lpot,l,s,j,t
   real(kind=r8),intent(in) :: rr
   real(kind=r8),intent(out) :: vpw(2,2)

   integer(kind=i4) :: ncc,ls,lsm,lsp
   real(kind=r8) :: vnn(18),vc,vt,vls,s12,s12m,s12p,s1ds2,t1dt2

   call cheft_pot(lpot,rr,vnn)

   s1ds2=4*s-3
   t1dt2=4*t-3
   vc=vnn(1)+t1dt2*vnn(2)+s1ds2*vnn(3)+s1ds2*t1dt2*vnn(4)
   vt=vnn(5)+t1dt2*vnn(6)
   vls=vnn(7)+t1dt2*vnn(8)

   ncc=1

   if(s.eq.1.and.j.gt.l)then
      ncc=2
   endif

   if(ncc.eq.1)then
      s12=0.0_r8

      if(s.eq.1.and.l.eq.j)then
         s12=2.0_r8
      endif

      if(l.eq.(j+1))then
         s12=-2.0_r8*(j+2)/(2.0_r8*j+1.0_r8)
      endif

      ls=(j*(j+1)-l*(l+1)-s*(s+1))/2

      vpw(1,1)=vc+s12*vt+ls*vls
      vpw(2,1)=0
      vpw(1,2)=0
      vpw(2,2)=0

   elseif(ncc.eq.2)then

      s12m=-2.0_r8*(j-1.0_r8)/(2.0_r8*j+1.0_r8)
      s12=sqrt(36.0_r8*j*(j+1))/(2.0_r8*j+1.0_r8)
      s12p=-2.0_r8*(j+2.0_r8)/(2.0_r8*j+1.0_r8)
      lsm=j-1
      lsp=-(j+2)
      vpw(1,1)=vc+s12m*vt+lsm*vls
      vpw(2,1)=s12*vt
      vpw(1,2)=s12*vt
      vpw(2,2)=vc+s12p*vt+lsp*vls

   endif

   return

   endsubroutine cheft_pw

!----------------------------------------------------------------------!
end module cheft
