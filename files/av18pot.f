c **********************************************************************
c Argonne v18 and vn' potential package 
c
c prepared 1 Sep 1994 by R. B. Wiringa, Physics Division, 
c Argonne National Laboratory, Argonne, IL 60439
c e-mail: wiringa@theory.phy.anl.gov 
c 
c reference: 
c  "Accurate nucleon-nucleon potential with charge-independence breaking"
c   R. B. Wiringa, V. G. J. Stoks, and R. Schiavilla, 
c   Physical Review C51, 38 (1995) - WSS95.
c
c option for v8' reprojection of v18 potential added 10 Jan 2001
c
c reference:
c  "Quantum Monte Carlo calculations of nuclei with A<=7"
c   B. S. Pudliner, V. R. Pandharipande, J. Carlson, Steven C. Pieper,
c   and R. B. Wiringa, Physical Review C56, 1720 (1997) - PPCPW97
c
c option for v6', v4', vx', v2', v1' potentials added 16 Jul 2002
c
c reference:
c  "Evolution of Nuclear Spectra with Nuclear Forces"
c   R. B. Wiringa and Steven C. Pieper,
c   Physical Review Letters 89, 182501 (2002) - WP02
c
c this file contains 4 subroutines:
c   subroutine av18pw(lpot,l,s,j,t,t1z,t2z,r,vpw)
c   subroutine av18op(lpot,r,vnn)
c   subroutine empot(lpot,r,vem)
c   subroutine consts(hc,mpi0,mpic,mp,mn,alpha,mup,mun)
c
c the variable lpot selects between v18, v8' and other options
c
c av18pw gives the full potential in a particular partial wave
c av18op gives the strong interaction part in operator format
c empot  gives the electromagnetic part in operator format
c consts gives values of fundamental constants and masses used
c
c notes: 
c 1) av18pw includes full EM interaction for lpot=1;
c    for lpot>1 it includes only C1(pp), i.e.,
c    Coulomb with a form factor for pp channels.
c
c 2) empot does not include the energy-dependence of the Coulomb
c    interaction used in eq.(4) of WSS95, i.e., it uses alpha, 
c    not alpha'.
c
c 3) the vacuum polarization in empot is a short-range approximation 
c    to eq.(7) suitable for bound states, but not for scattering.
c    it is taken from eq.(3.13) of Rev. Mod. Phys. 44, 48 (1972)
c
c    8/28/97 error in this formula detected and corrected:
c    should be -(gamma+5/6) instead of printed (-gamma+5/6)
c
c 4) these subroutines should be compiled with a compiler option
c    that forces all floating point constants to be evaluated at
c    real*8 significance, e.g., on an IBM RS6000 the xlf compiler
c    option qdpc=e should be used; on SGI machines, the -r8 option
c    should be used; on a Cray no action is needed.
c    if such an option is not available and the default precision is 
c    real*4 (32 bits), then all constants should be explicitly 
c    converted to double precision by appending a D0.
c
c *id* av18pw **********************************************************
c subroutine for partial-wave projection of argonne v18 potential
c or reprojected vn' potential
c calls subroutines av18op, empot
c ----------------------------------------------------------------------
c arguments for av18pw
c lpot: switch for potential choice
c       = 1 : v18
c       = 2 : v8'
c       = 3 : v6'
c       = 4 : v4'
c       = 5 : vx'
c       = 6 : v2'
c       = 7 : v1'
c l:    orbital angular momentum of pair (0,1,2,...)
c s:    total spin of pair (0 or 1)
c j:    total angular momentum of pair (0,1,2,...)
c t:    total isospin of pair (0 or 1)
c t1z:  isospin of particle 1 (1 for p, -1 for n)
c t2z:     "    "     "     2 (1 for p, -1 for n)
c r:    separation in fm
c v:    returned potential in MeV (2x2 array)
c       (includes all strong and em terms)
c ----------------------------------------------------------------------
c order of terms in v(l,m):
c      single channel                 coupled channel (l=j-1,s=1)
c      v(1,1) = v(l,s,j,t,t1z,t2z)    v(1,1) = v(l,s,j,t,t1z,t2z)
c      v(2,1) = 0                     v(2,1) = v(l<->l+2)
c      v(1,2) = 0                     v(1,2) = v(l<->l+2)
c      v(2,2) = 0                     v(2,2) = v(l+2,s,j,t,t1z,t2z)
c ----------------------------------------------------------------------
      subroutine av18pw(lpot,l,s,j,t,t1z,t2z,r,vpw)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer l,s,j,t,t1z,t2z,s1ds2,t1dt2,t12
      dimension vnn(18),vem(14),vpw(2,2)
c ------------------------
c strong interaction terms
c ------------------------
      call av18op(lpot,r,vnn)
      s1ds2=4*s-3
      t1dt2=4*t-3
      t12=3*t1z*t2z-t1dt2
      vc=vnn(1)+t1dt2*vnn(2)+s1ds2*vnn(3)+s1ds2*t1dt2*vnn(4)
     &  +t12*vnn(15)+s1ds2*t12*vnn(16)+(t1z+t2z)*vnn(18)
      vt=vnn(5)+t1dt2*vnn(6)+t12*vnn(17)
      vls=vnn(7)+t1dt2*vnn(8)
      vl2=vnn(9)+t1dt2*vnn(10)+s1ds2*vnn(11)+s1ds2*t1dt2*vnn(12)
      vls2=vnn(13)+t1dt2*vnn(14)
c ---------------------
c electromagnetic terms
c ---------------------
      call empot(lpot,r,vem)
      if (t1z+t2z) 10,20,30
   10 vc=vc+s1ds2*vem(7)
      vt=vt+vem(10)
      go to 40
   20 vc=vc+vem(5)+s1ds2*vem(8)
      vt=vt+vem(11)
      vls=vls+vem(14)
      go to 40
   30 vc=vc+vem(1)+vem(2)+vem(3)+vem(4)+s1ds2*vem(6)
      vt=vt+vem(9)
      vls=vls+vem(12)
   40 continue
c ---------------------
   50 ncc=1
      if (s.eq.1.and.j.gt.l) ncc=2
      if (ncc.eq.1) then
        s12=0.
        if (s.eq.1.and.l.eq.j) s12=2.
        if (l.eq.(j+1)) s12=-2.*(j+2.)/(2.*j+1.)
        ls=(j*(j+1)-l*(l+1)-s*(s+1))/2
        vpw(1,1)=vc+s12*vt+ls*vls+l*(l+1)*vl2+ls**2*vls2
        vpw(2,1)=0
        vpw(1,2)=0
        vpw(2,2)=0
      else if (ncc.eq.2) then
        s12m=-2.*(j-1.)/(2.*j+1.)
        s12=sqrt(36.*j*(j+1))/(2.*j+1.)
        s12p=-2.*(j+2.)/(2.*j+1.)
        lsm=j-1
        lsp=-(j+2)
        vpw(1,1)=vc+s12m*vt+lsm*vls+l*(l+1)*vl2+lsm**2*vls2
        vpw(2,1)=s12*vt
        vpw(1,2)=s12*vt
        vpw(2,2)=vc+s12p*vt+lsp*vls+(l+2)*(l+3)*vl2+lsp**2*vls2
      end if
      return
      end
c *id* av18op **********************************************************
c subroutine for strong interaction part of argonne v18 potential 
c or reprojected vn' potential in operator format
c calls subroutine consts
c ----------------------------------------------------------------------
c arguments for av18pot
c lpot: switch for potential choice
c       = 1 : v18
c       = 2 : v8'
c       = 3 : v6'
c       = 4 : v4'
c       = 5 : vx'
c       = 6 : v2'
c       = 7 : v1'
c r:    separation in fm
c vnn:  output potential in MeV (18 component array)
c ----------------------------------------------------------------------
c order of operators l in vnn(l):
c l:    1=1                              2=t1.t2
c       3=s1.s2                          4=(s1.s2)(t1.t2)
c       5=S12 [=3(s1.r)(s2.r)-s1.s2]     6=S12(t1.t2)
c       7=L.S                            8=L.S(t1.t2)
c       9=L**2                          10=L**2(t1.t2)
c      11=L**2(s1.s2)                   12=L**2(s1.s2)(t1.t2)
c      13=(L.S)**2                      14=(L.S)**2(t1.t2)
c      15=T12 [=3*t1z*t2z-t1.t2]        16=(s1.s2)T12
c      17=S12*T12                       18=t1z+t2z
c where s1=sigma_1, t1=tau_1, t1z=tau_1(z), etc.
c ----------------------------------------------------------------------
      subroutine av18op(lpot,r,vnn)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      dimension vnn(18)
      real*8 mpi0,mpic,mp,mn,mup,mun
      real*8 mpi,mu0,muc,mu
      data small/1.d-4/
      do 5 l=1,18
        vnn(l)=0
    5 continue
      call consts(hc,mpi0,mpic,mp,mn,alpha,mup,mun)
      mpi=(mpi0+2.*mpic)/3.
      mu0=mpi0/hc
      muc=mpic/hc
      mu=mpi/hc
      fsq=.075
      cpi=2.1
      rws=.5
      aiws=5.
      x=mu*r
      x0=mu0*r
      xc=muc*r
      if (r.le.small) then
        tpi=3*cpi**2*r/mu**3
        ypi0=(mpi0/mpic)**2*(mpi0/3)*cpi*r/mu0
        tpi0=3*cpi*ypi0/mu0**2
        ypic=(mpic/3)*cpi*r/muc
        tpic=3*cpi*ypic/muc**2
      else
        rcut=1-exp(-cpi*r*r)
        ypi=exp(-x)*rcut/x
        tpi=(1+(3+3/x)/x)*ypi*rcut
        ypi0=(mpi0/mpic)**2*(mpi0/3)*exp(-x0)*rcut/x0
        tpi0=(1+(3+3/x0)/x0)*ypi0*rcut
        ypic=(mpic/3)*exp(-xc)*rcut/xc
        tpic=(1+(3+3/xc)/xc)*ypic*rcut
      end if
      ypi0=fsq*ypi0
      ypic=fsq*ypic
      tpi0=fsq*tpi0
      tpic=fsq*tpic
      tpi2=tpi*tpi
      ws=1/(1+exp((r-rws)*aiws))
      ws0=1/(1+exp(-rws*aiws))
      wsp=ws*(1+aiws*exp(-rws*aiws)*ws0*r)
      wsx=ws*x
      wsx2=wsx*x
      dypi00=(mpi0/mpic)**2*(mpi0/3)*cpi/mu0
      dypic0=(mpic/3)*cpi/muc
      ypi0p=ypi0-fsq*dypi00*ws*r/ws0
      ypicp=ypic-fsq*dypic0*ws*r/ws0
      ypi=(ypi0+2*ypic)/3
      tpi=(tpi0+2*tpic)/3
      p11pp=  -7.62701*tpi2+1815.4920*wsp+1847.8059*wsx2+ypi0p
      p11np=  -7.62701*tpi2+1813.5315*wsp+1847.8059*wsx2-ypi0p+2*ypicp
      p11nn=  -7.62701*tpi2+1811.5710*wsp+1847.8059*wsx2+ypi0p
      pt1pp=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2+tpi0
      pt1np=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2-tpi0+2*tpic
      pt1nn=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2+tpi0
      pls1=    -.62697*tpi2 -570.5571*wsp +819.1222*wsx2
      pl211=    .06709*tpi2 +342.0669*wsp -615.2339*wsx2
      pls21=    .74129*tpi2   +9.3418*wsp -376.4384*wsx2
      p10=    -8.62770*tpi2+2605.2682*wsp +441.9733*wsx2-ypi0p-2*ypicp
      pt0=    1.485601*tpi2-1126.8359*wsx +370.1324*wsx2-tpi0-2*tpic
      pls0=     .10180*tpi2  +86.0658*wsp -356.5175*wsx2
      pl210=   -.13201*tpi2 +253.4350*wsp   -1.0076*wsx2
      pls20=    .07357*tpi2 -217.5791*wsp  +18.3935*wsx2
      p01pp= -11.27028*tpi2+3346.6874*wsp-3*ypi0p
      p01np= -10.66788*tpi2+3126.5542*wsp-3*(-ypi0p+2*ypicp)
      p01nn= -11.27028*tpi2+3342.7664*wsp-3*ypi0p
      pl201=    .12472*tpi2  +16.7780*wsp
      p00=    -2.09971*tpi2+1204.4301*wsp-3*(-ypi0p-2*ypicp)
      pl200=   -.31452*tpi2 +217.4559*wsp
      p11=(p11pp+p11nn+p11np)/3
      p11cd=(.5*(p11pp+p11nn)-p11np)/6
      p11cs=(p11pp-p11nn)/4
      pt1=(pt1pp+pt1nn+pt1np)/3
      pt1cd=(.5*(pt1pp+pt1nn)-pt1np)/6
      pt1cs=(pt1pp-pt1nn)/4
      p01=(p01pp+p01nn+p01np)/3
      p01cd=(.5*(p01pp+p01nn)-p01np)/6
      p01cs=(p01pp-p01nn)/4
c ------------------------
c option for v8' reduction
c ------------------------
      if (lpot.ge.2) then
        p00=p00+2*pl200
        p11=p11+2*pl211+4*pls21/3
        pt1=pt1-5*pls21/12
        pls1=pls1-.5*pls21
        pls0=pls0-2*pl210-3*pls20
      end if
c ------------------------
c option for v6' redutcion
c ------------------------
      if (lpot.ge.3) p10=p10-.3*pls0
c ------------------------
c option for v4' redutcion
c ------------------------
      if (lpot.ge.4) p10=p10+.8735*pt0
c ------------------------
c option for vx' reduction
c ------------------------
      if (lpot.eq.5) then
        vnn(1)=.0625*(9*p11+3*p10+3*p01+p00)
        vnn(2)=.0125*(9*p11-5*p10-5*p01+p00)
        vnn(3)=vnn(2)
        vnn(4)=vnn(2)
        return
c ------------------------
c option for v2' reduction
c ------------------------
      else if (lpot.eq.6) then
        vnn(1)=.25*(3*p01+p10)
        vnn(2)=.25*(  p01-p10)
        return
c ------------------------
c option for v1' reduction
c ------------------------
      else if (lpot.eq.7) then
        vnn(1)=.5*(p01+p10)
        return
      end if
c ------------------------
      vnn(1)=.0625*(9*p11+3*p10+3*p01+p00)
      vnn(2)=.0625*(3*p11-3*p10  +p01-p00)
      vnn(3)=.0625*(3*p11  +p10-3*p01-p00)
      vnn(4)=.0625*(  p11  -p10  -p01+p00)
c ------------------------
      if (lpot.eq.4) return
c ------------------------
      vnn(5)=.25*(3*pt1+pt0)
      vnn(6)=.25*(  pt1-pt0)
c ------------------------
      if (lpot.eq.3) return
c ------------------------
      vnn(7)=.25*(3*pls1+pls0)
      vnn(8)=.25*(  pls1-pls0)
c ------------------------
      if (lpot.eq.2) return
c ------------------------
      vnn(9)= .0625*(9*pl211+3*pl210+3*pl201+pl200)
      vnn(10)=.0625*(3*pl211-3*pl210+  pl201-pl200)
      vnn(11)=.0625*(3*pl211+  pl210-3*pl201-pl200)
      vnn(12)=.0625*(  pl211-  pl210-  pl201+pl200)
      vnn(13)=.25*(3*pls21+pls20)
      vnn(14)=.25*(  pls21-pls20)
      vnn(15)=.25*(3*p11cd+p01cd)
      vnn(16)=.25*(  p11cd-p01cd)
      vnn(17)=pt1cd
      vnn(18)=p01cs
      return
      end
c *id* empot ***********************************************************
c subroutine for electromagnetic part of Argonne v18 potential
c for avn' models returns pp Coulomb only
c calls subroutine consts
c ----------------------------------------------------------------------
c arguments for empot
c lpot: switch for potential choice
c       = 1 : full EM potential
c       > 1 : C1(pp) only
c r:    input separation in fm
c vem:  output potential in MeV (14 component array)
c ----------------------------------------------------------------------
c order of operators in vem(l)
c l:    1=C1    (pp)          2=DF    (pp)          3=C2      (pp)
c       4=VP    (pp)                                5=C1      (np)
c       6=s1.s2 (pp)          7=s1.s2 (nn)          8=s1.s2   (np)
c       9=S12   (pp)         10=S12   (nn)         11=S12     (np)
c      12=L.S   (pp)         13=L.S   (nn)         14=L.S     (np)
c C1 = one-photon-exchange Coulomb with form factor
c C2 = two-photon-exchange Coulomb
c DF = Darwin-Foldy
c VP = vacuum polarization (short-range approximation)
c all other terms from magnetic moment (MM) interactions
c ----------------------------------------------------------------------
      subroutine empot(lpot,r,vem)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      dimension vem(14)
      real*8 mpi0,mpic,mp,mn,mup,mun
      real*8 kr,me,mr
      data small/1.d-5/
      do 5 l=1,14
        vem(l)=0
    5 continue
      call consts(hc,mpi0,mpic,mp,mn,alpha,mup,mun)
      b=4.27
      br=b*r
      pi=acos(-1.)
      me=0.510999
      mr=mp*mn/(mp+mn)
      gamma=0.577216
      beta=.0189
      if (r.lt.small) then
        fcoulr=5*b/16
        ftr3=b**3*br**2/720
        flsr3=b**3/48
        kr=me*small/hc
      else
        fcoulr=(1-(1+11*br/16+3*br**2/16+br**3/48)*exp(-br))/r
        ftr3=(1-(1+br+br**2/2+br**3/6+br**4/24+br**5/144)*exp(-br))/r**3
        flsr3=(1-(1+br+br**2/2+7*br**3/48+br**4/48)*exp(-br))/r**3
        kr=me*r/hc
      end if
      fivp=-gamma-5./6.+abs(log(kr))+6*pi*kr/8
      fdelta=b**3*(1+br+br**2/3)*exp(-br)/16
      fnpr=b**3*(15+15*br+6*br**2+br**3)*exp(-br)/384
      vem(1)=alpha*hc*fcoulr
c ------------------------
      if (lpot.ge.2) return
c ------------------------
      vem(2)=-alpha*hc**3*fdelta/(4*mp**2)
      vem(3)=-vem(1)**2/mp
      vem(4)=2*alpha*vem(1)*fivp/(3*pi)
      vem(5)=alpha*hc*beta*fnpr
      vem(6)=-alpha*hc**3*mup**2*fdelta/(6*mp**2)
      vem(7)=-alpha*hc**3*mun**2*fdelta/(6*mn**2)
      vem(8)=-alpha*hc**3*mup*mun*fdelta/(6*mn*mp)
      vem(9)=-alpha*hc**3*mup**2*ftr3/(4*mp**2)
      vem(10)=-alpha*hc**3*mun**2*ftr3/(4*mn**2)
      vem(11)=-alpha*hc**3*mup*mun*ftr3/(4*mp*mn)
      vem(12)=-alpha*hc**3*(4*mup-1)*flsr3/(2*mp**2)
      vem(13)=0
      vem(14)=-alpha*hc**3*mun*flsr3/(2*mn*mr)
      return
      end
c *id* consts **********************************************************
c subroutine for constants in av18 potential
c ----------------------------------------------------------------------
c arguments for consts
c hc:    output value for hbar*c (MeV-fm)
c mpi0:    "      "    "  neutral pion mass (MeV)
c mpic:    "      "    "  charged pion mass (MeV)
c mp:      "      "    "  proton mass (MeV)
c mn:      "      "    "  neutron mass (MeV)
c alpha:   "      "    "  electromagnetic constant alpha
c mup:     "      "    "  proton magnetic moment (nm)
c mun:     "      "    "  neutron magnetic moment (nm)
c ----------------------------------------------------------------------
      subroutine consts(hc,mpi0,mpic,mp,mn,alpha,mup,mun)
      real*8 hc,mpi0,mpic,mp,mn,alpha,mup,mun
      hc=197.327053d0
      mpi0=134.9739d0
      mpic=139.5675d0
      mp=938.27231d0
      mn=939.56563d0
      alpha=1.0d0/137.035989
      mup= 2.7928474d0
      mun=-1.9130427d0
      return
      end
