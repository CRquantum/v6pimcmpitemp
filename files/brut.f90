!$Id: wavefunction.f90,v 1.2 2003/02/25 18:11:59 nuclear Exp $
    
! si, sj change sign. CR      
module brut ! this module mainly deal with correlation functions
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), private, parameter :: zero=0.0_r8,one=1.0_r8,two=2.0_r8,three=3.0_r8,four=4.0_r8 &
      ,five=5.0_r8,six=6.0_r8,seven=7.0_r8,eight=8.0_r8,anine=9.0_r8 &
      ,ten=10.0_r8,tenth=.1_r8,half=.5_r8,third=1.0_r8/3,sixth=1.0_r8/6 &
      ,pi=4.0_r8*atan(1.0_r8) 
! others
	integer(kind=i4), private, save :: nscale,irn(4)
	character*3 part(4)	
	data part /'pup','pdn','nup','ndn'/  
    integer(kind=i4), private, save :: nprot,n2p
	integer(kind=i4), private, save :: lpot 
	logical, private, save :: iem ! EM switch
! /fops/   
   real(kind=r8), private, allocatable, save :: f(:,:)
   real(kind=r8), private, save :: drop
   integer(kind=i4), private, save :: ng 
!
   data ng /2000/
   data nscale /2/
! /channel/   
   real(kind=r8), private, save :: dr,hbi,h12,range,akappa,alpha,omeg,dxke
   integer(kind=i4), private, save :: ngrid   

! /varpar/
   real(kind=r8), private, save :: ebind,capc,capr,amu,beta,varran
! /const/ 
   real(kind=r8), private, save ::  hb,delta,deltai
   integer(kind=i4), private, save :: nspin,niso,npair,n,np1,nisospin
!! /scrach/ 
   real(kind=r8), private, allocatable, save :: fa(:),fb(:) 
! /spins/
   complex(kind=r8), private, allocatable, save :: cwtgnd(:,:)
   integer(kind=i4), private, allocatable, save :: liso(:),lisoi(:)   
! /rnyucm/
   integer(kind=i4), private, save :: m1ran,m2ran,m3ran,m4ran,l1ran,l2ran,l3ran,l4ran
! /ground/ 
    integer(kind=i4), private, allocatable, save :: ignd(:)  
!	/estcum/ 
	real(kind=r8), private, save :: ecum,pecum,tpbcum
	integer(kind=i4), private, save :: jump
! /contrl/ 
	integer(kind=i4), private, save :: nstep,nblk,nblkeq,idump,irstar
! /config/ 
	real(kind=r8), private, allocatable, save :: xold(:,:),xnew(:,:),fold(:,:),fnew(:,:)
    real(kind=r8), private, save ::	psi2o,psi2n,eold,enew,peo,pen
	integer(kind=i4), private, allocatable, save :: iplo(:),jplo(:),ipro(:),jpro(:)	
! /contrl/	

! definition
! 0: n, down. 
! 1: p, up. 
! nisospin 1: 0011  =3  
!          2: 0101  =5
!          3: 0110  =6
!          4: 1001  =9
!          5: 1010  =10
!          6: 1100  =12
! nspin 0:15
	  
!
! common block variables
! /const/
!      n=# of particles
!      pi=3.14159...
!      hb=hbar**2/2m
!      delta=2*step size
!      deltai=1/delta
!      nspin=# of spin states
!      niso=#of isospin states =nisospin
!      npair=# of pairs of particles
! /spins/
!      cwtgnd(0:nspin-1,niso) the weight of the noninteracting
!                              ground state
!      liso(niso)=pointer from the allowed isospin states to their
!                 binary representation
!      lisoi(0:nspin-1)=inverse pointer
! /varpar/
!      ebind=guess for binding energy of last particle
!      capc=large r cutoff where wave function should
!           fall off exponentially
!      capr=wood-saxon width
!      amu=wood-saxon smoothness
!      beta=tensor-isospin strength
!      varran=range of variational functions
! /contrl/
!      nstep=# of steps/block
!      nblk=# of blocks of steps
!      nblkeq=# of steps before equilibrium
!      idump=flag to dump partial sums for a restart
!      irstar=flag to attempt a restart
! /config/
!      xold(3,n)=current position of particles
!      xnew(3,n)=new        "      "    "
!      f(3,n)=psi star grad psi/psi**2
!      psi2=psi**2
!      e=psi star h psi/psi**2
!      pe=psi star v psi/psi**2
!      iplo(npairs),jplo(npairs)=current pair order on left
!      ipro(npairs),jpro(npairs)=   "      "    "    " right
! /ground/
!      ignd(n)=ground state configuration read in
!              1=pup
!              2=pdn
!              3=nup
!              4=ndn
!
!
!
! read in data
! quantities not given above
!      title=6 character title for run
!      date=6 character date
!      eunit=units of run
!      rn=20 octal digit random # seed
!
!
! calculate constants and set up random # generator
!	  
!
! write out data
!
	  
contains	   
	  subroutine brutinit(npartin,nprotin,hbarin,iemin) 
      use mympi
      use math
	  real(kind=r8) :: hbarin
	  integer(kind=i4) :: npartin,nprotin,lpotin
	  integer(kind=i4) :: i,irn(4)
	  logical :: iemin

	  iem=iemin
	  n=npartin
      nprot=nprotin
	  hb=hbarin
      irn=1
	  
	  delta=3.0_r8
	  ebind=16.0_r8
	  capc=1.0_r8
	  capr=1.0_r8
	  amu=0.5_r8
	  beta=1.0_r8
	  varran=4.0_r8	  
	  !lpot=lpotin ! 11 is v6' in pot.f and 3 in av18pot.f
      
      ngrid=nscale*ng
      range=varran
      alpha=one/(n-one)
      akappa=alpha*sqrt(ebind*(n-1)/(hb*float(n)))
      hbi=one/(two*hb)
      dr=range/ngrid ! do not forget to set dxx.
      h12=dr**2/(six*two)
      drop=dr*nscale ! drop=range/ng
      dxke=dr ! the dxx for KE
      n2p=6*n+1 ! deal with KE

      ! call setrn0(irn) I don't use rannyu
      deltai=one/delta
 

      allocate(f(6,ng))
      allocate(ignd(n))
      allocate(xold(3,n),xnew(3,n),fold(3,n),fnew(3,n))      
      
	  do i=1,n
		  ignd(i)=i   ! the ignd in v6old. 1234 = pu pd nu nd
      enddo
	  
      if (myrank().eq.0) then
          write (6,'(t10,''trial binding energy ='',t40,f10.5)') ebind
          write (6,'(t10,''large r cut off ='',t40,f10.5)') capc
          write (6,'(t10,''wood saxon width ='',t40,f10.5)') capr
          write (6,'(t10,''wood saxon smoothness ='',t40,f10.5)') amu
          write (6,'(t10,''t tau strength ='',t40,f10.5)') beta
          write (6,'(t10,''range of correlations ='',t40,f10.5)') varran
          write (6,'(t10,''ground state configuration '',(t50,a3))') (part(ignd(i)),i=1,n)	 
      endif
       
      !call setpot(lpot,one,one,one,one,one,one,h2m,h2mcsb)	  
	  
      nspin=2**n
      call combin(n,nprot,nisospin) ! give value for nisospin, n=4, it is 4C2= 6.
      !write(6,*) 'nisospin=',nisospin
      allocate(liso(nisospin),lisoi(0:nspin-1),cwtgnd(0:nspin-1,nisospin))      
	  call spinit0(ignd) ! get nspin, niso
	  npair=n*(n-1)/2
      np1=npair-1 
      allocate(iplo(npair),jplo(npair),ipro(npair),jpro(npair))
	  call getord(iplo,jplo) ! first set npair in it.
      call getord(ipro,jpro)
         !write(6,*) 'iplo=',iplo
         !write(6,*) 'jplo=',jplo  
         !write(6,*) 'ipro=',ipro
         !write(6,*) 'jpro=',jpro
      allocate(fa(ngrid),fb(ngrid))  
	  call fcorop
	   
	  return
	  end subroutine brutinit		  

      subroutine coropl(x,cwtn)
      complex(kind=r8) :: cwtn(0:nspin-1,nisospin)
      real(kind=r8) :: x(3,n)	  
	  call corop(x,iplo,jplo,cwtn)
         !write(6,*) 'iplo=',iplo
         !write(6,*) 'jplo=',jplo
	  end
	  
      subroutine coropr(x,cwtn)
      complex(kind=r8) :: cwtn(0:nspin-1,nisospin)
      real(kind=r8) :: x(3,n)
	  call corop(x,ipro,jpro,cwtn)
	  end

      subroutine corop(x,ip,jp,cwtn)
!
! subroutine to calculate the weight of the wave function
! in each spin isospin state given the correlation functions
!
      complex(kind=r8) :: xif,xjf,wtx
      complex(kind=r8) :: cwtn(0:nspin-1,nisospin),wt(0:nspin-1,nisospin,2)
      real(kind=r8) :: x(3,n),dx(3),dxri(3),fop(6)
      integer(kind=i4) :: ip(npair),jp(npair)
      real(kind=r8) :: fc,fs,fi,fsi,ft,fti
      equivalence (fop(1),fc),(fop(2),fs),(fop(3),fi),(fop(4),fsi),(fop(5),ft),(fop(6),fti)
      integer(kind=i4) :: k,l,mn,mo,ii,i,j,ispin,ic,jspin,kspin,lspin,ifsgni,jfsgni,maskex,jiexi &
		                 ,iflps,jflps,jiflps,ifsgn,jfsgn,jiexs,m1
      real(kind=r8) :: rij,riji,si,sj,zi,zj
      real(kind=r8) :: fceff,fsex,fiex,fsiex,fteff,ftieff
!     integer(kind=i4) shiftl,shiftr,xor,and,or
!
! vms vax and pc booleans
!
!     shiftl(ishift,nbit)=ishft(ishift,nbit)
!     shiftr(ishift,nbit)=ishft(ishift,-nbit)
!     xor(i1,i2)=ieor(i1,i2)
!     and(i1,i2)=iand(i1,i2)
!     or(i1,i2)=ior(i1,i2)
!
! sun booleans
!
!     shiftl(ishift,nbit)=lshift(ishift,nbit)
!     shiftr(ishift,nbit)=rshift(ishift,nbit)
!
! set initial weight to ground state
!
      do 10 k=1,nisospin
      do 10 l=1,nspin
   10 wt(l-1,k,1)=cwtgnd(l-1,k)
      mn=1
      mo=2
!
! loop over pairs
!     
    !open(11,FILE='rij.txt',FORM='FORMATTED')
      do 20 ii=1,npair

      i=ip(ii)
      j=jp(ii)
      ispin=shiftl(1,i-1)
!
! calculate pair distances
!
!      rij=zero
!      do 40 ic=1,3
!      dx(ic)=x(ic,i)-x(ic,j)
!   40 rij=rij+dx(ic)**2
!      rij=sqrt(rij)
!	  !write(6,*) 'rij=',rij
!      riji=one/rij
!      do 50 ic=1,3
!50    dxri(ic)=dx(ic)*riji
!      
!      write(6,*) 'drxi old=', dxri
         
      dx(:)=x(:,i)-x(:,j)
      rij=sqrt(sum(dx(:)**2))
      dxri(:)=dx(:)/rij      
      
      mn=mo
      mo=3-mn
!
! zero new weight
!
!      do 60 k=1,nisospin
!      do 60 l=1,nspin
!60    wt(l-1,k,mn)=(zero,zero)
      
      wt(:,:,mn)=0
!
! calculate pair correlations
! loop over isospin states
!
      call calf(rij,fop)
	  !write (6,'(2x,7(f15.7,2x))') rij,fop
		!fop=0.
		!fop(6)=1. 
	  !write(6,*) 'fop=',fop 
      fceff=fc-fs-fi+fsi+ft-fti
      fsex=two*(fs-fsi-ft+fti)
      fiex=two*(fi-fsi+fti)
      fsiex=four*(fsi-fti)
      fteff=three*(ft-fti)
      ftieff=six*fti
      
      !write (11,'(2x,6(f15.7,2x))') fceff,fsex,fiex,fsiex,fteff,ftieff
        
      jspin=shiftl(1,j-1)
      do 70 l=1,nisospin
      lspin=liso(l)
      ifsgni=shiftr(lspin,i-1)
      jfsgni=shiftr(lspin,j-1)
      m1=and(1,xor(ifsgni,jfsgni))
      maskex=or(shiftl(m1,i-1),shiftl(m1,j-1))
      jiexi=xor(maskex,lspin)
      jiexi=lisoi(jiexi)
!
! jiexi=isospin exchange index
!
!
! loop over spin states
!
      do 80 k=1,nspin
      kspin=k-1
      iflps=xor(kspin,ispin)
      jflps=xor(kspin,jspin)
      jiflps=xor(iflps,jspin)
      ifsgn=and(1,shiftr(kspin,i-1))
      jfsgn=and(1,shiftr(kspin,j-1))
      m1=xor(ifsgn,jfsgn)
      maskex=or(shiftl(m1,i-1),shiftl(m1,j-1))
      jiexs=xor(maskex,kspin)
         
      !write (11,*) 'others1='
      !write (11,'(2x,10(f15.7,2x))') iflps,jflps,jiflps,ifsgn,jfsgn,m1,maskex,jiexs
         
!
! iflps=index for i th spin flipped
! jflps=  "    "  j  "  "    "
! jiflps= "    " j+i   spins  "
! jiexs=  "    "  "     "    exchanged
! si=eigenvalue of sigma z for i th spin
! sj=   "       "    "   "  "  j th "
!
      si=-(-2*ifsgn+1)
      
      !write(11,*) 'jfsgn,sj=',jfsgn,-(-2*jfsgn+1)
      
      sj=-(-2*jfsgn+1)
      
      !write(11,*) 'jfsgn,sj again=',jfsgn,-(-2*jfsgn+1),sj
      
      zi=si*dxri(3)
      zj=sj*dxri(3)
      xif=cmplx(dxri(1),dxri(2)*si)
      xjf=cmplx(dxri(1),dxri(2)*sj)
      wtx=wt(kspin,l,mo)
          
      !write (11,*) 'wtx='
      !write (11,'(2x,2(f15.7,2x))') wtx
      !write (11,*) 'others='
      !write (11,'(2x,10(f15.7,2x))') si,sj,zi,zj,xif,xjf     
            
!
! calculate new weight
!
      wt(kspin,l,mn)=wt(kspin,l,mn)+wtx*(fceff+fteff*zi*zj)
      wt(iflps,l,mn)=wt(iflps,l,mn)+wtx*fteff*xif*zj
      wt(jflps,l,mn)=wt(jflps,l,mn)+wtx*fteff*xjf*zi
      wt(jiflps,l,mn)=wt(jiflps,l,mn)+wtx*fteff*xif*xjf
      wt(jiexs,l,mn)=wt(jiexs,l,mn)+wtx*fsex
      wt(kspin,jiexi,mn)=wt(kspin,jiexi,mn)+wtx*(fiex+zi*zj*ftieff)
      wt(iflps,jiexi,mn)=wt(iflps,jiexi,mn)+wtx*ftieff*zj*xif
      wt(jflps,jiexi,mn)=wt(jflps,jiexi,mn)+wtx*ftieff*zi*xjf
      wt(jiflps,jiexi,mn)=wt(jiflps,jiexi,mn)+wtx*ftieff*xif*xjf
      wt(jiexs,jiexi,mn)=wt(jiexs,jiexi,mn)+wtx*fsiex
      
      !write (11,*) 'wt='
      !write (11,'(2x,2(f15.7,2x))') wt(kspin,l,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(iflps,l,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(jflps,l,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(jiflps,l,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(jiexs,l,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(kspin,jiexi,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(iflps,jiexi,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(jflps,jiexi,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(jiflps,jiexi,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(jiexs,jiexi,mn)
            
   80 continue
   70 continue
   20 continue
!
! put new weight in cwtn
!
      do 90 k=1,nisospin
      do 90 l=1,nspin
      cwtn(l-1,k)=wt(l-1,k,mn)
      
      !write (11,'(2x,2(f15.7,2x))') wt(l-1,k,mn)
      
90    continue      
      !write (11,'(2x,2(f30.20,2x))') cwtn     
      !close(11)
      return
      end
      

      subroutine coropv(x,ip,jp,dxx,cwta)  
      use mympi
! subroutine to calculate the weight of the wave function
! in each spin isospin state given the correlation functions
!
      complex(kind=r8) :: xif(n2p),xjf(n2p)
      complex(kind=r8) :: cwta(n2p,0:nspin-1,nisospin),wt(n2p,0:nspin-1,nisospin,2),wtx(n2p)
      real(kind=r8) :: x(3,n),dx(n2p,3),dxri(n2p,3) &
       ,rij(n2p),zi(n2p),zj(n2p),riji(n2p) &
       ,xa(n2p,3,4),fop(6,n2p),fceff(n2p),fsex(n2p),fiex(n2p) &
       ,fsiex(n2p),fteff(n2p),ftieff(n2p)
      integer(kind=i4) :: ip(npair),jp(npair)
      integer(kind=i4) :: k,l,la,mn,mo,i,ic,ic1,i1,ii,j,ispin,jspin,kspin,lspin &
		                 ,ifsgni,jfsgni,maskex,jiexi,iflps,jflps,jiflps,ifsgn,jfsgn,jiexs,m1
      real(kind=r8) :: dxx,si,sj
!     integer(kind=i4) shiftl,shiftr,xor,and,or
!
! vms vax and pc booleans
!
!     shiftl(ishift,nbit)=ishft(ishift,nbit)
!     shiftr(ishift,nbit)=ishft(ishift,-nbit)
!     xor(i1,i2)=ieor(i1,i2)
!     and(i1,i2)=iand(i1,i2)
!     or(i1,i2)=ior(i1,i2)
!
! sun booleans
!
!     shiftl(ishift,nbit)=lshift(ishift,nbit)
!     shiftr(ishift,nbit)=rshift(ishift,nbit)
!
! set initial weight to ground state
!

!      do 10 k=1,nisospin
!      do 10 l=1,nspin
!      do 10 la=1,n2p
!10    wt(la,l-1,k,1)=cwtgnd(l-1,k)
      
      do la=1,n2p
         wt(la,:,:,1)=cwtgnd(:,:)    
      enddo
 
      mn=1
      mo=2
!
! loop over pairs
!
!      do 100 ic1=1,3
!      do 100 i1=1,n
!100   xa(1,ic1,i1)=x(ic1,i1)
      xa(1,:,:)=x(:,:)  
      l=1  
!      do 110 i=1,n
!      do 110 ic=1,3
!      x(ic,i)=x(ic,i)+dxx
!      l=l+1
!      do 120 ic1=1,3
!      do 120 i1=1,n
!  120 xa(l,ic1,i1)=x(ic1,i1)
!      x(ic,i)=x(ic,i)-two*dxx
!      l=l+1
!      do 130 ic1=1,3
!      do 130 i1=1,n
!  130 xa(l,ic1,i1)=x(ic1,i1)
!      x(ic,i)=x(ic,i)+dxx
!110   continue
      
      do i=1,n
       do ic=1,3
          x(ic,i)=x(ic,i)+dxx 
          l=l+1
          xa(l,:,:)=x(:,:)
          x(ic,i)=x(ic,i)-two*dxx
          l=l+1
          xa(l,:,:)=x(:,:)
          x(ic,i)=x(ic,i)+dxx 
       enddo
      enddo
     
      !open(11,FILE='rija.txt',FORM='FORMATTED')
   
      do 20 ii=1,npair
      i=ip(ii)
      j=jp(ii)
      ispin=shiftl(1,i-1)
!
! calculate pair distances
!
          
!      do 140 la=1,n2p
!140   rij(la)=zero
        
      rij(:)=zero
      
!      do 40 ic=1,3
!      do 40 la=1,n2p
!      dx(la,ic)=xa(la,ic,i)-xa(la,ic,j)
!40    rij(la)=rij(la)+dx(la,ic)**2
      
      do ic=1,3
         dx(:,ic)=xa(:,ic,i)-xa(:,ic,j) 
         rij(:)=rij(:)+dx(:,ic)**2 
      enddo
         
!      do 150 la=1,n2p
!150   rij(la)=sqrt(rij(la))
      
      rij(:)=sqrt(rij(:))
        
!      do 160 la=1,n2p
!160   riji(la)=one/rij(la)  
!      do 50 ic=1,3
!      do 50 la=1,n2p
!50    dxri(la,ic)=dx(la,ic)*riji(la)
      
      do ic=1,3
        dxri(:,ic)=dx(:,ic)/rij(:)
      enddo
 
      mn=mo
      mo=3-mn
!
! zero new weight
!
      
!      do 60 k=1,nisospin
!      do 60 l=1,nspin
!      do 60 la=1,n2p
!60    wt(la,l-1,k,mn)=(zero,zero)
      
      wt(:,:,:,mn)=(zero,zero)
      
!      do 170 la=1,n2p
!170   call calf(rij(la),fop(:,la))
      
      do la=1,n2p
         call calf(rij(la),fop(:,la)) 
      enddo
   
      !write (11,'(2x,7(f15.7,2x))') rij(1),fop(:,1)
!
! calculate pair correlations
!
!
! loop over isospin states
!
      do 180 la=1,n2p
      fceff(la)=fop(1,la)-fop(2,la)-fop(3,la)+fop(4,la)+fop(5,la)-fop(6,la)
      fsex(la)=two*(fop(2,la)-fop(4,la)-fop(5,la)+fop(6,la))
      fiex(la)=two*(fop(3,la)-fop(4,la)+fop(6,la))
      fsiex(la)=four*(fop(4,la)-fop(6,la))
      fteff(la)=three*(fop(5,la)-fop(6,la))
180   ftieff(la)=six*fop(6,la)
       
      !write (11,'(2x,6(f15.7,2x))') fceff(1),fsex(1),fiex(1),fsiex(1),fteff(1),ftieff(1)
       
      jspin=shiftl(1,j-1)
      do 70 l=1,nisospin
      lspin=liso(l)
      ifsgni=shiftr(lspin,i-1)
      jfsgni=shiftr(lspin,j-1)
      m1=and(1,xor(ifsgni,jfsgni))
      maskex=or(shiftl(m1,i-1),shiftl(m1,j-1))
      jiexi=xor(maskex,lspin)
      jiexi=lisoi(jiexi)
!
! jiexi=isospin exchange index
!
!
! loop over spin states
!
      do 80 k=1,nspin
      kspin=k-1
      iflps=xor(kspin,ispin)
      jflps=xor(kspin,jspin)
      jiflps=xor(iflps,jspin)
      ifsgn=and(1,shiftr(kspin,i-1))
      jfsgn=and(1,shiftr(kspin,j-1))
      m1=xor(ifsgn,jfsgn)
      maskex=or(shiftl(m1,i-1),shiftl(m1,j-1))
      jiexs=xor(maskex,kspin)
      
      !write (11,*) 'others1='
      !write (11,'(2x,10(f15.7,2x))') iflps,jflps,jiflps,ifsgn,jfsgn,m1,maskex,jiexs
      
      
!
! iflps=index for i th spin flipped
! jflps=  "    "  j  "  "    "
! jiflps= "    " j+i   spins  "
! jiexs=  "    "  "     "    exchanged
! si=eigenvalue of sigma z for i th spin
! sj=   "       "    "   "  "  j th "
!
      si=-(-2*ifsgn+1)
      
      !write(11,*) 'jfsgn,sj=',jfsgn,-(-2*jfsgn+1)
      
      sj=-(-2*jfsgn+1)
      
      !write(11,*) 'jfsgn,sj again=',jfsgn,-(-2*jfsgn+1),sj
      
      do 190 la=1,n2p
      zi(la)=si*dxri(la,3)
      zj(la)=sj*dxri(la,3)
      xif(la)=cmplx(dxri(la,1),dxri(la,2)*si)
      xjf(la)=cmplx(dxri(la,1),dxri(la,2)*sj)
190   wtx(la)=wt(la,kspin,l,mo)
      
       !write (11,*) 'wtx='      
       !write (11,'(2x,2(f15.7,2x))') wtx(1)
       !write (11,*) 'others='
       !write (11,'(2x,10(f15.7,2x))') si,sj,zi(1),zj(1),xif(1),xjf(1) 
!
! calculate new weight
!
      do 200 la=1,n2p
      wt(la,kspin,l,mn)=wt(la,kspin,l,mn)+wtx(la)*(fceff(la)+fteff(la)*zi(la)*zj(la))
      wt(la,iflps,l,mn)=wt(la,iflps,l,mn)+wtx(la)*fteff(la)*xif(la)*zj(la)
      wt(la,jflps,l,mn)=wt(la,jflps,l,mn)+wtx(la)*fteff(la)*xjf(la)*zi(la)
      wt(la,jiflps,l,mn)=wt(la,jiflps,l,mn)+wtx(la)*fteff(la)*xif(la)*xjf(la)
      wt(la,jiexs,l,mn)=wt(la,jiexs,l,mn)+wtx(la)*fsex(la)
      wt(la,kspin,jiexi,mn)=wt(la,kspin,jiexi,mn)+wtx(la)*(fiex(la)+zi(la)*zj(la)*ftieff(la))
      wt(la,iflps,jiexi,mn)=wt(la,iflps,jiexi,mn)+wtx(la)*ftieff(la)*zj(la)*xif(la)
      wt(la,jflps,jiexi,mn)=wt(la,jflps,jiexi,mn)+wtx(la)*ftieff(la)*zi(la)*xjf(la)
      wt(la,jiflps,jiexi,mn)=wt(la,jiflps,jiexi,mn)+wtx(la)*ftieff(la)*xif(la)*xjf(la)
      wt(la,jiexs,jiexi,mn)=wt(la,jiexs,jiexi,mn)+wtx(la)*fsiex(la)
      
      !if (la.eq.1) then
      !write (11,*) 'wt='
      !write (11,'(2x,2(f15.7,2x))') wt(1,kspin,l,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(1,iflps,l,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(1,jflps,l,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(1,jiflps,l,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(1,jiexs,l,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(1,kspin,jiexi,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(1,iflps,jiexi,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(1,jflps,jiexi,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(1,jiflps,jiexi,mn)
      !write (11,'(2x,2(f15.7,2x))') wt(1,jiexs,jiexi,mn)   
      !endif
      
   200 continue     
   80 continue
   70 continue
   20 continue
      
   
!
! put new weight in cwtn
!
      do 90 k=1,nisospin
      do 90 l=1,nspin
      do 90 la=1,n2p
      cwta(la,l-1,k)=wt(la,l-1,k,mn)
      !if (la.eq.1) write (11,'(2x,2(f15.7,2x))') wt(la,l-1,k,mn)  
   90 continue
      !write (11,'(2x,2(f30.20,2x))') cwta(1,:,:)
      !close(11)
      return
	  end
	  
      subroutine pecal(x,cwtl,cwtr,petot)
!
! subroutine to calculate potential energy
! compare with corop
! cwtl and cwtr are the weights for the spin isospin
! states for the wave function on the right and left
!
      complex(kind=r8) xif,xjf,pex,wtx
	  complex(kind=r8) :: cwtl(0:nspin-1,nisospin),cwtr(0:nspin-1,nisospin)
      real(kind=r8) :: x(3,n),dx(3),dxri(3)
      integer(kind=i4) :: i,im1,ispin,j,ic,jspin,l,lspin,k,kspin &
		                 ,ifsgni,jfsgni,maskex,jiexi,iflps,jflps,jiflps,ifsgn,jfsgn,jiexs,m1
      real(kind=r8) :: rij,rij2,riji,si,sj,zi,zj,petot
	  real(kind=r8) :: vc,vs,vi,vsi,vt,vti,vb,vbi,vceff,vsex,viex,vsiex,vteff,vtieff
!
! vms vax and pc booleans
!
!     shiftl(ishift,nbit)=ishft(ishift,nbit)
!     shiftr(ishift,nbit)=ishft(ishift,-nbit)
!     xor(i1,i2)=ieor(i1,i2)
!     and(i1,i2)=iand(i1,i2)
!     or(i1,i2)=ior(i1,i2)
!
! sun booleans
!
!     shiftl(ishift,nbit)=lshift(ishift,nbit)
!     shiftr(ishift,nbit)=rshift(ishift,nbit)
      pex=(zero,zero)      
      do 20 i=2,n
      im1=i-1
      ispin=shiftl(1,im1)
      do 30 j=1,im1
!      rij=zero
!      rij2=zero
!      do 40 ic=1,3
!      dx(ic)=x(ic,i)-x(ic,j)
!40    rij=rij+dx(ic)**2      
      dx(:)=x(:,i)-x(:,j)
      rij=sqrt(sum(dx(:)**2))
      dxri(:)=dx(:)/rij
      !riji=one/rij
      !rij2=sum(dx(:)**2)
      !dxri(:)= sign(1.0,dx(:))*sqrt( dx(:)**2/rij2 )
!      do 50 ic=1,3
!50    dxri(ic)=dx(ic)/rij   !*riji
      call potmod(vc,vs,vi,vsi,vt,vti,vb,vbi,rij)
      vceff=vc-vs-vi+vsi+vt-vti
      vsex=two*(vs-vsi-vt+vti)
      viex=two*(vi-vsi+vti)
      vsiex=four*(vsi-vti)
      vteff=three*(vt-vti)
      vtieff=six*vti
      jspin=shiftl(1,j-1)
      do 70 l=1,nisospin
      lspin=liso(l)
      ifsgni=shiftr(lspin,i-1)
      jfsgni=shiftr(lspin,j-1)
      m1=and(1,xor(ifsgni,jfsgni))
      maskex=or(shiftl(m1,i-1),shiftl(m1,j-1))
      jiexi=xor(maskex,lspin)
      jiexi=lisoi(jiexi)
      do 80 k=1,nspin
      kspin=k-1
      iflps=xor(kspin,ispin)
      jflps=xor(kspin,jspin)
      jiflps=xor(iflps,jspin)
      ifsgn=and(1,shiftr(kspin,i-1))
      jfsgn=and(1,shiftr(kspin,j-1))
      m1=xor(ifsgn,jfsgn)
      maskex=or(shiftl(m1,i-1),shiftl(m1,j-1))
      jiexs=xor(maskex,kspin)
      si=-(-2*ifsgn+1)
      sj=-(-2*jfsgn+1)
      zi=si*dxri(3)
      zj=sj*dxri(3)
      xif=dcmplx(dxri(1),-dxri(2)*si)
      xjf=dcmplx(dxri(1),-dxri(2)*sj)
      wtx=conjg(cwtl(kspin,l))
      
      !pex=pex+cwtr(kspin,l)*wtx*(vceff+vteff*zi*zj)
      !pex=pex+cwtr(iflps,l)*wtx*vteff*xif*zj
      !pex=pex+cwtr(jflps,l)*wtx*vteff*xjf*zi
      !pex=pex+cwtr(jiflps,l)*wtx*vteff*xif*xjf
      !pex=pex+cwtr(jiexs,l)*wtx*vsex
      !pex=pex+cwtr(kspin,jiexi)*wtx*(viex+zi*zj*vtieff)
      !pex=pex+cwtr(iflps,jiexi)*wtx*vtieff*zj*xif
      !pex=pex+cwtr(jflps,jiexi)*wtx*vtieff*zi*xjf
      !pex=pex+cwtr(jiflps,jiexi)*wtx*vtieff*xif*xjf
      !pex=pex+cwtr(jiexs,jiexi)*wtx*vsiex

      pex=pex+wtx*(  cwtr(kspin,l)*(vceff+vteff*zi*zj) &
                    +cwtr(iflps,l)*vteff*xif*zj &
                    +cwtr(jflps,l)*vteff*xjf*zi &
                    +cwtr(jiflps,l)*vteff*xif*xjf &
                    +cwtr(jiexs,l)*vsex &
                    +cwtr(kspin,jiexi)*(viex+zi*zj*vtieff) &
                    +cwtr(iflps,jiexi)*vtieff*zj*xif &
                    +cwtr(jflps,jiexi)*vtieff*zi*xjf &
                    +cwtr(jiflps,jiexi)*vtieff*xif*xjf &
                    +cwtr(jiexs,jiexi)*vsiex  )
      
   80 continue
   70 continue
   30 continue
   20 continue
      petot=pex      
      return
      end	  

      subroutine empecal(x,cwtl,cwtr,empetot)
      complex(kind=r8) :: cwtl(0:nspin-1,nisospin),cwtr(0:nspin-1,nisospin)
      real(kind=r8) :: x(3,n),dx(3),dxri(3),vem(14),vem2(14),emtot,r,empetot
	  integer(kind=i4) :: l,lspin,i,j,ipi,ipj
         
      do l=1,nisospin
         lspin=liso(l)       
         emtot=0.
         do i=1,n-1               
          do j=i+1,n                
             dx(:)=x(:,i)-x(:,j)		   
		     r=sqrt(dot_product(dx,dx)) !r=sqrt(sum(dx(:)**2))         
             ipi=(2*(and(1,shiftr(lspin,i-1))))/2
             ipj=(2*(and(1,shiftr(lspin,j-1))))/2  
             if ((ipi*ipj).eq.1) then
              call empot(2,r,vem) ! only pp simpliest coulomb here.
			  
			  !call empotmod(0,one,one,r,vem2)
			  !write(6,*) 'vem, vem2 are', vem(1),vem2(1)
			  !if ( abs(vem(1)-vem2(1)).gt.1.d-5 ) stop
			  
              emtot=emtot+vem(1) 
             endif          
          enddo       
         enddo  
         cwtr(:,l)=cwtr(:,l)*emtot        
       enddo
      empetot=sum(conjg(cwtl(:,:))*cwtr(:,:))
      
      return
      end      
      
      
      
      subroutine calf(r,fop)
!
! subroutine to interpolate value of correlation function
! from tables
!
      real(kind=r8) :: fop(6)
	  real(kind=r8) :: dri,r,rx,drx,xm1,xp1,xp2,w1,w2,w3,w4
	  integer(kind=i4) :: ix,k
	  
	  dri=one/drop  !!!!!!!!!! either do it again here or make it common in this module.
	  
      if (r.gt.range) go to 20
      rx=r*dri
      ix=rx
      !write(6,*) 'ix=', ix,r,dri,rx,ng-2,max(2,ix)
	  
      ix=min(ng-2,max(2,ix))
	
      drx=ix-rx
      xm1=drx-one
      xp1=drx+one
      xp2=xp1+one
      w1=sixth*drx*xp1*xp2
      w2=-half*xm1*xp1*xp2
      w3=half*xm1*drx*xp2
      w4=-sixth*xm1*drx*xp1
      do 30 k=1,6
		 ! write(6,*) 'ix=', ix,r,drop,jump,dri,rx
   30 fop(k)=w1*f(k,ix-1)+w2*f(k,ix)+w3*f(k,ix+1)+w4*f(k,ix+2)
      return
   20 fop(1)=exp(-akappa*r)*r**(-alpha)
      do 40 k=2,6
   40 fop(k)=zero
      fop(6)=-beta*omeg*fop(1)*third
      return
	  end
	    
      subroutine fcorop
      use mympi
!
! subroutine to calculate tables of correlation functions
!
      real(kind=r8) :: f01(ngrid),f10(ngrid),ft0(ngrid)
	  integer(kind=i4) :: i,j,ic
      !data ng /100/
      !data nscale /2/
!
! set up parameters
!
!
! solve equations in coupled s=1 t=0 channels
!
      call cupchn(f10,ft0)
!
! solve equations in s=0 t=1 channel
!
      call sngchn(f01)
!
! set up tables
!
      do 10 j=1,ng
      i=j*nscale
      f(1,j)=half*half*(f01(i)+three*f10(i))
      f(2,j)=half*half*(f10(i)-f01(i))
      f(3,j)=zero
      f(4,j)=zero
      f(5,j)=zero
10    f(6,j)=-ft0(i)*beta*third   
      if (myrank().eq.0) then
      write (6,'(///,1x,''correlation functions calculated'')')
      open(unit=77,file='correlations.unf',form='unformatted')
	  open(unit=11,file='correlations.dat',form='formatted')
      write (77) ng
      write (77) varran
      write (77) akappa
      write (77) alpha
      write (77) beta
      write (77) omeg
      do 20 j=1,ng
      write (77) (f(ic,j),ic=1,6)
      write (11,'(6f15.7)') (f(ic,j),ic=1,6)
   20 continue
      close(77)
	  close(11)
      endif   
      return
	  end
	  	  
      subroutine cupchn(f10,ft0)
      real(kind=r8), parameter :: small=1.d-8
!
! subroutine to solve coupled channel equations in s=1 t=0
! channels
!
      real(kind=r8) :: f10(ngrid),ft0(ngrid)
      real(kind=r8) :: gam1,gam2,gam3,beta1,fabc11,fbbc11,beta2,fabc12,fbbc12,beta3,fabc2,fax,fbx
	  integer(kind=i4) :: k,kk,i
	  real(kind=r8) :: r
      gam1=ten
      do 10 k=1,20
      if (k.eq.1) beta1=tenth*tenth
      call cupbc0(gam1,beta1,fabc11,fbbc11)
      beta2=beta1+tenth**3
      do 20 kk=1,20
      call cupbc0(gam1,beta2,fabc12,fbbc12)
      if (abs(fbbc12-fbbc11).lt.small) go to 30
      beta3=beta1-fbbc11*(beta2-beta1)/(fbbc12-fbbc11)
      fbbc11=fbbc12
      beta1=beta2
   20 beta2=beta3
   30 continue
      if (k.ne.1) go to 40
      gam2=gam1
      gam1=gam1+one
      fabc2=fabc12
      go to 10
   40 if (abs(fabc12-fabc2).lt.small) go to 50
      gam3=gam2-fabc2*(gam1-gam2)/(fabc12-fabc2)
      fabc2=fabc12
      gam2=gam1
      gam1=gam3
   10 continue
   50 continue
      do 100 i=1,ngrid
      r=dr*i
      fax=fa(i)/(r)
      fbx=fb(i)/(r)
      f10(i)=fax
      ft0(i)=fbx
  100 continue
      omeg=beta2
      return
	  end
	  
	  
      subroutine cupbc0(gam,omegx,fabc,fbbc)
!
! subroutine to integrate coupled channel equations in s=1 t=0
! channels
!
      real(kind=r8) :: gam,omegx,fabc,fbbc
	  integer(kind=i4) :: i,j
	  real(kind=r8) :: r,v00,v01,v10,v11,vt0,vt1,vb0,vb1,ri,fa0,fb0,fa1,fb1,xa,ya,yb,deti,xb
      real(kind=r8), allocatable, save :: vaa(:),vab(:),vba(:),vbb(:),vgam(:)
      logical :: precaldone=.false.
      
      if (.not.precaldone) then
          allocate(vaa(ngrid),vab(ngrid),vba(ngrid),vbb(ngrid),vgam(ngrid))  
          do i=1,ngrid
              r=dr*i
              call potchn(v00,v01,v10,v11,vt0,vt1,vb0,vb1,r)
              ri=one/r
              vaa(i)=hbi*v10+(akappa**2-two*(one-alpha) &
              *akappa*ri+alpha*(alpha-one)*ri**2)* &
              (one-exp(-(r/capc)**2))
              vab(i)=hbi*eight*vt0
              vbb(i)=six*ri**2+vaa(i)-two*hbi*vt0-three*hbi*vb0
              vba(i)=hbi*vt0
              vgam(i)=hbi/(one+exp((r-capr)/amu))      
          enddo    
          precaldone=.true.
      endif

      fa(ngrid)=exp(-akappa*range)*range**(one-alpha)
      fb(ngrid)=omegx*fa(ngrid)
      fa(ngrid-1)=exp(-akappa*(range-dr))*(range-dr)**(one-alpha)
      fb(ngrid-1)=omegx*fa(ngrid-1)
      fa0=fa(ngrid)
      fb0=fb(ngrid)
      fa1=fa(ngrid-1)
      fb1=fb(ngrid-1)
      do 30 j=3,ngrid
      i=ngrid-j+1
      xa=(-fa0*(one-h12*(vaa(i+2)+gam*vgam(i+2)))+fa1*(two+ten* &
      h12*(vaa(i+1)+gam*vgam(i+1))) &
        +h12*(fb0*vab(i+2)+ten*fb1*vab(i+1)))/(one-h12*(vaa(i) &
      +gam*vgam(i)))
      xb=(-fb0*(one-h12*(vbb(i+2)+gam*vgam(i+2)))+fb1*(two+ten* &
      h12*(vbb(i+1)+gam*vgam(i+1))) &
        +h12*(fa0*vba(i+2)+ten*fa1*vba(i+1)))/(one-h12*(vbb(i) &
      +gam*vgam(i)))
      ya=h12*vab(i)
      yb=h12*vba(i)
      deti=one/(one-ya*yb)
      fa(i)=(xa+ya*xb)*deti
      fb(i)=(xb+yb*xa)*deti
      fa0=fa1
      fa1=fa(i)
      fb0=fb1
      fb1=fb(i)
   30 continue
      fabc=two*fa(1)-fa(2)
      fbbc=two*fb(1)-fb(2)
      return
	  end
	    
      subroutine sngchn(f01)
      real(kind=r8), parameter :: small=1.d-8      
      real(kind=r8) :: f00(ngrid),f01(ngrid)
!
! subroutine to solve single channel equation for s=0 t=1
!
      real(kind=r8) :: gam1,fbc1,gam2,fbc2,gam3
	  integer(kind=i4) :: i,k
	  real(kind=r8) :: r
      gam1=ten
      call sngbc1(gam1,fbc1)
      gam2=gam1+one
      do 110 k=1,20
      call sngbc1(gam2,fbc2)
      if (abs(fbc2-fbc1).lt.small) go to 120
      gam3=gam1-fbc1*(gam2-gam1)/(fbc2-fbc1)
      fbc1=fbc2
      gam1=gam2
  110 gam2=gam3	  
  120 continue
      do 200 i=1,ngrid
      r=dr*i
      f01(i)=fa(i)/(r)
  200 continue
      return
	  end
	  	  
      subroutine sngbc1(gam,fbc)
!
! subroutine to integrate single channel equation in s=0 t=1
!
      real(kind=r8) :: gam,fbc
	  integer(kind=i4) :: i,j
	  real(kind=r8) :: r,v00,v01,v10,v11,vt0,vt1,vb0,vb1,ri,fa0,fa1
      real(kind=r8), allocatable, save :: vaa(:),vab(:),vba(:),vbb(:),vgam(:)
      logical :: precaldone=.false.
      
      if (.not.precaldone) then
          allocate(vaa(ngrid),vab(ngrid),vba(ngrid),vbb(ngrid),vgam(ngrid))       
          do i=1,ngrid
              r=dr*i
              call potchn(v00,v01,v10,v11,vt0,vt1,vb0,vb1,r)
              ri=one/r
              vaa(i)=hbi*v01+(akappa**2-two*(one-alpha) &
              *akappa*ri+alpha*(alpha-one)*ri**2)* &
              (one-exp(-(r/capc)**2))
              vgam(i)=hbi/(one+exp((r-capr)/amu))
          enddo
          precaldone=.true.
      endif      
      
      fa(ngrid)=exp(-akappa*range)*range**(one-alpha)
      fa(ngrid-1)=exp(-akappa*(range-dr))*(range-dr)**(one-alpha)
      fa0=fa(ngrid)
      fa1=fa(ngrid-1)
      do 30 j=3,ngrid
      i=ngrid-j+1
      fa(i)=(-fa0*(one-h12*(vaa(i+2)+gam*vgam(i+2)))+fa1*(two &
      +ten*h12*(vaa(i+1)+gam*vgam(i+1)))) &
        /(one-h12*(vaa(i)+gam*vgam(i)))
      fa0=fa1
      fa1=fa(i)
   30 continue
      fbc=two*fa(1)-fa(2)
      return
	  end
	  	  
      subroutine potchn(v00,v01,v10,v11,vt0,vt1,vb0,vb1,r)
      real(kind=r8) :: v00,v01,v10,v11,vt0,vt1,vb0,vb1
	  integer(kind=i4) :: i,j
	  real(kind=r8) :: r,rx,vc,vs,vi,vsi,vt,vti,vb,vbi  
      rx=r
!
! subroutine to convert the 8 operator potentials calculated
! by pot to the 8 channel potentials used in the constrained varia
! equations. In v00,v01,v10,v11 the first no. is the total
! spin, and the second the total isospin. Vt0,vt1 are the tensor
! potentials in the isospin=0,1 channels.
!
      call potmod(vc,vs,vi,vsi,vt,vti,vb,vbi,rx)
      v00=vc-three*vs-three*vi+anine*vsi
      v01=vc-three*vs+vi-three*vsi
      v10=vc+vs-three*vi-three*vsi
      v11=vc+vs+vi+vsi
      vt0=vt-three*vti
      vt1=vt+vti
      vb1=vb+vbi
      vb0=vb-three*vbi
      return
	  end
	  
      subroutine potmod (vc,vs,vi,vsi,vt,vti,vb,vbi,r)
      real(kind=r8) :: vv(18),ww(14),vp(12)
      real(kind=r8) :: vc,vs,vi,vsi,vt,vti,vb,vbi,r
      !call pot(0,r,vv,vp,ww)
	  call av18op(3,r,vv)  ! 3 is v6'
      vc=vv(1)
      vi=vv(2)
      vs=vv(3)
      vsi=vv(4)
      vt=vv(5)
      vti=vv(6)
      vb=vv(7)
      vbi=vv(8)
      vb=zero
      vbi=zero
      return
      end

    
      subroutine hpsitcwt(x,ip,jp,cwt,vj,e,ke,pe,empe,cwta) ! it changes x.
      use mympi
      real(kind=r8), parameter :: dxx=0.001_r8  ! roughly dr
!
! subroutine to calculate the values of psi e pe and grad psi
! given positions and order of correlations
!
      real(kind=r8) :: x(3,n),d2(3,n),ff(3,n),fd(3,n),d2d(3,n),fj(3,n),d2j(3,n),x0(3,n)
      complex(kind=r8) :: cwta(n2p,0:nspin-1,nisospin),cwtl(0:nspin-1,nisospin),cwtr(0:nspin-1,nisospin),cwt(0:nspin-1,nisospin)
      integer(kind=i4) :: ip(npair),jp(npair),ipl(npair),ipr(npair),jpl(npair),jpr(npair) ! should be 6 i think.
	  real(kind=r8) :: xcm,psi2,e,ke,pe,vj,empe,vjp,vjm
	  integer(kind=i4) ::  ic,i,l,k,kspin,la      
!
! calculate wave functions on left and right

      do ic=1,3 
        xcm=sum(x(ic,:))/n
        x(ic,:)=x(ic,:)-xcm
      enddo
      cwtl=cwt
      call coropv(x,ip,jp,dxke,cwta)
      
      !call corop(x,ip,jp,cwtl) ! need to comment.  
      !write(6,'( ''cwtl ='', 12(g15.7,1x) )') cwtl
      !write(6,'( ''cwta ='', 12(g15.7,1x) )') cwta(1,:,:)
      !write(6,'( ''cwta - cwtl ='', 12(g15.7,1x) )') cwta(1,:,:)-cwtl
       
      vj=sum(conjg(cwtl(:,:))*cwta(1,:,:))  ! this is psi2 which is conjg(cwt)*psitcwt, real part
!
! calculate the potential energy
!   
      cwtr(:,:)=cwta(1,:,:)
      call pecal(x,cwtl,cwtr,pe)
      
      empe=0
	  if (iem) call empecal(x,cwtl,cwtr,empe) ! EM potential on
 
! new ke     
      la=1
      do i=1,n
       do ic=1,3
          la=la+1
          vjp=sum(conjg(cwtl(:,:))*cwta(la,:,:))
          la=la+1
          vjm=sum(conjg(cwtl(:,:))*cwta(la,:,:))
          d2(ic,i)=vjp+vjm-2*vj
       enddo
      enddo  
      ke=-hb*sum(d2(:,:))/dxke**2    
      !write(6,*) 'new num ke=', ke 
      e=ke+pe+empe !  e is just the numerator, without the denominator conjg(cwt)*rs0hpsit
     
      return
      end subroutine hpsitcwt 
	  

      subroutine vpsitcwt(x,ip,jp,cwt,pe) ! it changes x.
      real(kind=r8), parameter :: dxx=.02d0 
      real(kind=r8) :: x(3,n),d2(3,1),ff(3,1),fd(3,4),d2d(3,4),fj(3,4),d2j(3,4)
      complex(kind=r8) :: cwta(n2p,0:nspin-1,nisospin),cwtl(0:nspin-1,nisospin),cwtr(0:nspin-1,nisospin),cwt(0:nspin-1,nisospin)
      integer(kind=i4) :: ip(npair),jp(npair),ipl(npair),ipr(npair),jpl(npair),jpr(npair) ! should be 6 i think.
	  real(kind=r8) :: xcm,psi2,e,pe,vj,empe,vjp,vjm
	  integer(kind=i4) ::  ic,i,l,k,kspin,la
   
!      do 200 ic=1,3
!      xcm=zero
!      do 210 i=1,n
!  210 xcm=xcm+x(ic,i)
!      xcm=xcm/n
!      do 220 i=1,n
!  220 x(ic,i)=x(ic,i)-xcm
!200   continue
  
      
     do ic=1,3 
        xcm=sum(x(ic,:))/n
        x(ic,:)=x(ic,:)-xcm
     enddo      
   
      !call corop(x,ip,jp,cwtl)
	  cwtl=cwt
      call corop(x,ip,jp,cwta(1,:,:))  ! the cwta(la=a) in coropv should be the same as the cwt for corop
!
! calculate the potential energy
!
      do 100 l=1,nisospin
      do 100 k=1,nspin
  100 cwtr(k-1,l)=cwta(1,k-1,l)
      call pecal(x,cwtl,cwtr,pe)
      return
      end subroutine vpsitcwt    
      
      subroutine psi02(x,ip,jp,cwt,psi2) 
! it changes x.
! given positions and order of correlations
!
      real(kind=r8) :: x(3,n)
      complex(kind=r8) :: cwta(n2p,0:nspin-1,nisospin),cwtl(0:nspin-1,nisospin),cwtr(0:nspin-1,nisospin),cwt(0:nspin-1,nisospin)
      integer(kind=i4) :: ip(npair),jp(npair),ipl(npair),ipr(npair),jpl(npair),jpr(npair) ! should be 6 i think.
	  real(kind=r8) :: xcm,psi2,vj
	  integer(kind=i4) ::  ic,i,l,k,kspin,la
    
!      do 200 ic=1,3
!      xcm=zero
!      do 210 i=1,n
!  210 xcm=xcm+x(ic,i)
!      xcm=xcm/n
!      do 220 i=1,n
!  220 x(ic,i)=x(ic,i)-xcm
!200   continue
 
      do ic=1,3 
        xcm=sum(x(ic,:))/n
        x(ic,:)=x(ic,:)-xcm
      enddo      
 
      !call corop(x,ip,jp,cwtl)
	  cwtl=cwt
      call corop(x,ip,jp,cwta(1,:,:))
  
      vj=zero
      do 10 l=1,nisospin
      do 20 k=1,nspin
      kspin=k-1
   20 vj=vj+conjg(cwtl(kspin,l))*cwta(1,kspin,l)
   10 continue
      psi2=vj ! psi2 is conjg(cwt)*psitcwt, real part
      return
      end subroutine psi02     
      
      subroutine xconvert(x) ! it changes x to x-xcm
      real(kind=r8) :: x(3,n)
	  real(kind=r8) :: xcm
	  integer(kind=i4) ::  ic  
      do ic=1,3 
        xcm=sum(x(ic,:))/n
        x(ic,:)=x(ic,:)-xcm
      enddo
      return
      end subroutine xconvert
      	   
      subroutine getord(ip,jp)
      use math
	  use random
	  real(kind=r8) :: rn(1),xtemp(npair)
      
!
! subroutine to get a random order of pairs
!
      integer(kind=i4) :: ip(npair),jp(npair)
      integer(kind=i4), allocatable, save :: ijpair(:,:) 
! if write like ijpair(npair), it will not remember the previous ijpair array. but if write like ijpair(2,6), it will remember. be careful. 
      integer(kind=i4) :: l,i,im1,j,ip1,iswap
	  real(kind=r8) :: swap
      logical :: skip=.false.
      
! set up initial array containing valid pairs
!      l=0
!      do 20 i=2,n
!      im1=i-1
!      do 20 j=1,im1
!      l=l+1
!      ijpair(1,l)=i
!20    ijpair(2,l)=j
           
 	  !write(6,*) 'skip=', skip    
      if (.not.skip) then
          if (.not. allocated(ijpair)) allocate(ijpair(2,npair))
          l=0
          do i=2,n
             im1=i-1 
             do j=1,im1
                l=l+1
                ijpair(1,l)=i
                ijpair(2,l)=j    
             enddo  
          enddo
          skip=.true.    
      endif
    
! check, try not to repeat the same calculations again.      
!
! get an array of random # s
!
!      write(6,*) 'npair=', npair
  
!      do 30 i=1,npair
!	  rn=randn(1)
!      xtemp(i)=rn(1) ! I don't use rannyu
!	  !write(6,*) 'xtemp=',xtemp(i) 
!      ip(i)=ijpair(1,i)
!30    jp(i)=ijpair(2,i) 
 
         
      do i=1,npair
         rn=randn(1)
         xtemp(i)=rn(1) 
         ip(i)=ijpair(1,i)
         jp(i)=ijpair(2,i)
      enddo   
!
! sort random numbers to get a random permutation
!
      do 40 i=1,np1
      ip1=i+1
      do 40 j=ip1,npair
      if (xtemp(i).gt.xtemp(j)) go to 40
      swap=xtemp(i)
      xtemp(i)=xtemp(j)
      xtemp(j)=swap
      iswap=ip(i)
      ip(i)=ip(j)
      ip(j)=iswap
      iswap=jp(i)
      jp(i)=jp(j)
      jp(j)=iswap
40	  continue
	  !write(6,*) 'getord finish!'
      return
	  end subroutine getord
	  
      subroutine getordlr(ipl,jpl,ipr,jpr)
	  integer(kind=i4) :: ipl(npair),jpl(npair),ipr(npair),jpr(npair)	 
	   !write(6,*) 'getordlr start!'	  
	  call getord(ipl,jpl) ! first set npair in it.
      call getord(ipr,jpr)  
      return
	  end
	  
      subroutine outputordlr(ipl,jpl,ipr,jpr)
	  integer(kind=i4) :: ipl(npair),jpl(npair),ipr(npair),jpr(npair)	
	  ipl=iplo
	  jpl=jplo
	  ipr=ipro
	  jpr=jpro	
      return
      end	  
	  
      subroutine inputordlr(ipl,jpl,ipr,jpr)
	  integer(kind=i4) :: ipl(npair),jpl(npair),ipr(npair),jpr(npair)	
	  iplo=ipl
	  jplo=jpl
	  ipro=ipr
	  jpro=jpr	
      return
      end
      
      
      subroutine spinit0(ignd) ! this is the original spinit
!
! subroutine to initialize spin variables and the non
! interacting ground state
!
	  integer(kind=i4) :: in
      integer(kind=i4) :: itemp1(4),itemp2(4),ignd(4)
	  integer(kind=i4) :: ineut,i,l,lspin,j,lx,n1,i1,iswap,k,kspin
	  real(kind=r8) :: si
!     integer(kind=i4) shiftl,shiftr,xor,and,or
!
! vms vax and pc booleans
!
!     shiftl(ishift,nbit)=ishft(ishift,nbit)
!     shiftr(ishift,nbit)=ishft(ishift,-nbit)
!     xor(i1,i2)=ieor(i1,i2)
!     and(i1,i2)=iand(i1,i2)
!     or(i1,i2)=ior(i1,i2)
!
! sun booleans
!
!     shiftl(ishift,nbit)=lshift(ishift,nbit)
!     shiftr(ishift,nbit)=rshift(ishift,nbit)
	  integer(kind=i4) :: m1=1,m2=2
!
! count # of neutrons
!
      ineut=0
      do 10 i=1,n
      if (ignd(i).gt.2) ineut=ineut+1
   10 continue
!      nspin=2**n
      niso=0
!
! loop through isospin states and record the ones
! satisfying charge conservation
!
      do 20 l=1,nspin
      lspin=l-1
      in=0
      do 30 j=1,n
      lx=and(m1,shiftr(lspin,j-1))
      if (lx.eq.1) in=in+1
   30 continue
      if (in.ne.ineut) go to 20
      niso=niso+1
      liso(niso)=lspin
      lisoi(lspin)=niso
   20 continue
!
! put the read in state in array
!
      do 40 i=1,n
   40 itemp1(i)=ignd(i)
      n1=n-1
!
! sort array
!
      do 50 i=1,n1
      i1=i+1
      do 60 j=i1,n
      if (itemp1(i).gt.itemp1(j)) go to 60
      iswap=itemp1(i)
      itemp1(i)=itemp1(j)
      itemp1(j)=iswap
   60 continue
   50 continue
!
! go through allowed spin isospin states and determine
! which are permutations of th ground state
!
      do 70 l=1,niso
      lspin=shiftl(liso(l),1)
      do 80 k=1,nspin
      kspin=k-1
      cwtgnd(kspin,l)=(zero,zero)
      do 90 i=1,n
   90 itemp2(i)=and(m1,shiftr(kspin,i-1))+and(m2,shiftr(lspin,i-1))+1
!
! sort states and record sign of exchanges
!
      si=one
      do 100 i=1,n1
      i1=i+1
      do 110 j=i1,n
      if (itemp2(i).gt.itemp2(j)) go to 110
      si=-si
      iswap=itemp2(i)
      itemp2(i)=itemp2(j)
      itemp2(j)=iswap
  110 continue
  100 continue
!
! if same as ground state then fix weight
!
      do 120 i=1,n
      if (itemp1(i).ne.itemp2(i)) go to 130
  120 continue
      cwtgnd(kspin,l)=cmplx(si,zero)
  130 continue
   80 continue
   70 continue
      return
	  end subroutine spinit0	  

!      function rannyu(i)
!      real(kind=r8), parameter :: twom12=0.000244140625d0
!	  real(kind=r8) :: rannyu
!	  integer(kind=i4) :: i
!!
!! generic statement functions
!!
!!     ishft12(ii)=ii/4096
!!     mask12(ii)=mod(ii,4096)
!!
!! unix f77 statement functions
!!
!      ishft12(ii)=rshift(ii,12)
!      mask12(ii)=and(ii,4095)
!!
!! fps statement functions
!!
!!     ishft12(ii)=shift(ii,-12)
!!     mask12(ii)=and(ii,4095)
!!
!! cray cft statement functions
!!
!!     ishft12(ii)=shiftr(ii.12)
!!     mask12(ii)=and(ii,4095)
!!
!! vms fortran and convex fc statement functions
!!
!!     ishft12(ii)=ishft(ii,-12)
!!     mask12(ii)=iand(ii,4095)
!      i1=l1ran*m4ran+l2ran*m3ran+l3ran*m2ran+l4ran*m1ran
!      i2=l2ran*m4ran+l3ran*m3ran+l4ran*m2ran
!      i3=l3ran*m4ran+l4ran*m3ran
!      i4=l4ran*m4ran
!      l4ran=mask12(i4)
!      i3=i3+ishft12(i4)
!      l3ran=mask12(i3)
!      i2=i2+ishft12(i3)
!      l2ran=mask12(i2)
!      l1ran=mask12(i1+ishft12(i2))
!      rannyu=twom12*(l1ran+twom12*(l2ran+twom12*(l3ran+twom12*(l4ran))))
!      return
!	  end
!      subroutine setrn0(iseed)
!      integer(kind=i4) iseed(4),m(4),l(4)
!      data m / 502,1521,4071,2107/
!      data l /   0,   0,   0,   1/  
!	  m1ran=m(1)
!	  m2ran=m(2)
!	  m3ran=m(3)
!	  m4ran=m(4)
!	  l1ran=l(1)
!	  l2ran=l(2)
!	  l3ran=l(3)
!	  l4ran=l(4) 
!	  do 10 i=1,4
!         l(i)=iseed(i)
!   10    continue
!      l(4)=2*(l(4)/2)+1
!      return
!	  end
!      subroutine savern0(iseed)
!      integer(kind=i4) iseed(4),m(4),l(4)
!	  m(1)=m1ran
!	  m(2)=m2ran
!	  m(3)=m3ran
!	  m(4)=m4ran
!	  l(1)=l1ran
!	  l(2)=l2ran
!	  l(3)=l3ran
!	  l(4)=l4ran	  
!	  do 10 i=1,4
!         iseed(i)=l(i)
!   10    continue
!      return
!      end	
	  
!      subroutine hpsivmc(x,ipl,jpl,ipr,jpr,psi2,e,d2,pe,f)
!      
!      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
!      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,anine=9.d0)
!      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)
!      parameter (dxx=.02d0)
!!
!! subroutine to calculate the values of psi e pe and grad psi
!! given positions and order of correlations
!!
!      complex(kind=r8) cwtl,cwta,cwtr
!      dimension x(3,1),d2(3,1),f(3,1)
!      dimension cwta(n2p,0:nspin-1,nisospin),cwtl(0:nspin-1,nisospin),cwtr(0:nspin-1,nisospin)
!      dimension fd(3,4),d2d(3,4),fj(3,4),d2j(3,4)
!      dimension ipl(npair),ipr(npair),jpl(npair),jpr(npair)
!!
!! calculate wave functions on left and right
!!
!      do 200 ic=1,3
!      xcm=zero
!      do 210 i=1,n
!  210 xcm=xcm+x(ic,i)
!      xcm=xcm/n
!      do 220 i=1,n
!  220 x(ic,i)=x(ic,i)-xcm
!  200 continue
!      call corop(x,ipl,jpl,cwtl)
!      call coropv(x,ipr,jpr,dxx,cwta)
!!
!! calculate value of psi**2
!!
!      vj=zero
!      do 10 l=1,nisospin
!      do 20 k=1,nspin
!      kspin=k-1
!   20 vj=vj+dconjg(cwtl(kspin,l))*cwta(1,kspin,l)
!   10 continue
!      psi2=vj
!!
!! calculate the potential energy
!!
!      do 100 l=1,nisospin
!      do 100 k=1,nspin
!  100 cwtr(k-1,l)=cwta(1,k-1,l)
!      call pecal(x,cwtl,cwtr,pe)
!      pe=pe/vj
!!
!! calculate the kinetic energy numerically
!!
!      la=1
!      do 30 i=1,n
!      do 30 ic=1,3
!      la=la+1
!      vjp=zero
!      do 40 l=1,nisospin
!      do 50 k=1,nspin
!      kspin=k-1
!   50 vjp=vjp+dconjg(cwtl(kspin,l))*cwta(la,kspin,l)
!   40 continue
!      la=la+1
!      vjm=zero
!      do 60 l=1,nisospin
!      do 70 k=1,nspin
!      kspin=k-1
!   70 vjm=vjm+dconjg(cwtl(kspin,l))*cwta(la,kspin,l)
!   60 continue
!      fj(ic,i)=(vjp-vjm)/(two*vj*dxx)
!      d2j(ic,i)=(vjp+vjm-two*vj)/(dxx**2*vj)-fj(ic,i)**2
!   30 continue
!!
!! sum up stuff
!!
!      e=zero
!      do 80 i=1,n
!      do 80 ic=1,3
!      d2(ic,i)=d2j(ic,i)
!   80 f(ic,i)=fj(ic,i)
!      do 90 i=1,n
!      do 90 ic=1,3
!   90 e=e-hb*(d2(ic,i)+f(ic,i)**2)
!      e=e+pe
!      return
!	  end	
	
	
	
	
	
	
 
end module brut
