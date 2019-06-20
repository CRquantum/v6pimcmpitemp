!$Id: wavefunction.f90,v 1.2 2003/02/25 18:11:59 nuclear Exp $
module v6stepcalc
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   complex(kind=r8), private, parameter :: czero=(0.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: ci=(0.0_r8,1.0_r8)
   complex(kind=r8), private, parameter :: cone=(1.0_r8,0.0_r8)
   integer(kind=i4), private, save :: npart,nprot,nneut
   integer(kind=i4), private, save :: nisospin,nspin,nbasis,niso
   integer(kind=i4), private, save :: mmax,nbisect,ntab,nrepmax
   real(kind=r8), private, save :: scalep,rangev,rangevsafe
   real(kind=r8), private, save :: hbar,dt
   complex(kind=r8), private, save, allocatable :: cwt(:,:),invcwt(:,:) ! invcwt given nspin,niso, give nstate accordingly.
   real(kind=r8), private, allocatable, save :: utab(:,:,:),vtab(:,:),evdttab(:,:,:),vemtab(:) &
                                               ,utabpr(:,:,:),evdttabpr(:,:,:)
   real(kind=r8), private,save, dimension(6) :: u0tab
! 
   integer(kind=i4), private, allocatable, save :: invspin(:),invispin(:),liso(:),lisoi(:)
   complex(kind=r8), private, save, allocatable :: cwtgnd(:,:)    
   logical, private, save :: iem
   integer(kind=i4), private, save :: irep
! 0: n, down. 
! 1: p, up. 
! nisospin 1: 0011  =3  
!          2: 0101  =5
!          3: 0110  =6
!          4: 1001  =9
!          5: 1010  =10
!          6: 1100  =12
! nspin 0:15
contains	   
	  subroutine v6stepcalcinit(npartin,nprotin,dtin,hbarin,ntabin,rangevin,mmaxin,nrepmaxin,iemin,irepin)
	  use v6pot
      use mympi
	  real(kind=r8) :: dtin,hbarin,rangevin,dr
	  integer(kind=i4) :: npartin,nprotin,mmaxin,nbisectin,ntabin,nrepmaxin
	  real(kind=r8) :: v(6),e(6) ! v(6) should already absorbed the stime step dt. 
	  integer(kind=i4) :: i,j,k,ibox,lpot,lpotpr,irepin
      real(kind=r8) :: vfact(8)   
	  logical :: iemin
	  irep=irepin
	  iem=iemin ! em force switch
	  call spinit(npartin,nprotin)
	  dt=dtin
	  hbar=hbarin
	  ntab=ntabin
	  rangev=rangevin
	  mmax=mmaxin
      nrepmax=nrepmaxin
	  nbisect=2**mmax
	  dr=rangev/ntab 
      scalep=1._r8/dr 
! make a table for u for different r and dt.  time step=2**j dt.
	  !write(6,*) 'u table start!'
	  vfact=1.0_r8
	  ibox=0
	  lpot=11  ! 11 for v6
	  lpotpr=11
      call v6potinit(npart,ibox,-rangev,lpot,ntab,vfact,lpotpr,iem) !lpot=11 for v6'.
	  if (.not.allocated(vtab)) allocate(vtab(6,0:ntab)) 
	  
      !call getvtab(vtab)
      call getvtabschmidt(vtab)

	  if (iem) then
		allocate(vemtab(0:ntab)) 
		call getvemtab(vemtab) 
      endif
	  
	!vtab=0.
	!vtab(6,:)=1.	 
	    
	  !write(6,*) 'get v table done!'  
      if (.not.allocated(utab)) allocate(utab(6,0:ntab,-mmax:mmax))
      if (.not.allocated(utabpr)) allocate(utabpr(6,0:ntab,-mmax:mmax)) ! utabpr deal with +delta t
	  if (.not.allocated(evdttab)) allocate(evdttab(6,0:ntab,-mmax:mmax))
	  if (.not.allocated(evdttabpr)) allocate(evdttabpr(6,0:ntab,-mmax:mmax))	  
	  do j=-mmax,mmax ! dt=2**j
	   do i=0,ntab
       ! we use Schmidt order, note that it is different from the order in v18pot. so need to change it.    
	   v(:)=vtab(:,i)*dt*2.0_r8**j  
	   e(1)=exp(-(v(1)+v(2)+2*v(3)+v(4)+v(5)+2*v(6))) !c
	   e(2)=exp(-(v(1)+v(2)-4*v(3)+v(4)+v(5)-4*v(6))) !sigma
	   e(3)=exp(-(v(1)+v(2)+2*v(3)-3*v(4)-3*v(5)-6*v(6))) !t	
	   e(4)=exp(-(v(1)+v(2)-4*v(3)-3*v(4)-3*v(5)+12*v(6))) !tau
	   e(5)=exp(-(v(1)-3*v(2)+v(4)-3*v(5))) !sigma tau
	   e(6)=exp(-(v(1)-3*v(2)-3*v(4)+9*v(5))) ! t tau	
	   evdttab(:,i,j)=e(:)	   
	   utab(1,i,j)=(6*e(1)+3*e(2)+2*e(3)+e(4)+3*e(5)+e(6))/16
	   utab(2,i,j)=(6*e(1)+3*e(2)+2*e(3)+e(4)-9*e(5)-3*e(6))/48
	   utab(3,i,j)=(3*e(1)-3*e(2)+e(3)-e(4))/24
	   utab(4,i,j)=(2*e(1)+e(2)-2*e(3)-e(4)+e(5)-e(6))/16
	   utab(5,i,j)=(2*e(1)+e(2)-2*e(3)-e(4)-3*e(5)+3*e(6))/48
	   utab(6,i,j)=(e(1)-e(2)-e(3)+e(4))/24	
       
       v(:)=vtab(:,i)*(-dt)*2.0_r8**j  
	   e(1)=exp(-(v(1)+v(2)+2*v(3)+v(4)+v(5)+2*v(6))) !c
	   e(2)=exp(-(v(1)+v(2)-4*v(3)+v(4)+v(5)-4*v(6))) !sigma
	   e(3)=exp(-(v(1)+v(2)+2*v(3)-3*v(4)-3*v(5)-6*v(6))) !t	
	   e(4)=exp(-(v(1)+v(2)-4*v(3)-3*v(4)-3*v(5)+12*v(6))) !tau
	   e(5)=exp(-(v(1)-3*v(2)+v(4)-3*v(5))) !sigma tau
	   e(6)=exp(-(v(1)-3*v(2)-3*v(4)+9*v(5))) ! t tau	
	   evdttabpr(:,i,j)=e(:)	
	   utabpr(1,i,j)=(6*e(1)+3*e(2)+2*e(3)+e(4)+3*e(5)+e(6))/16
	   utabpr(2,i,j)=(6*e(1)+3*e(2)+2*e(3)+e(4)-9*e(5)-3*e(6))/48
	   utabpr(3,i,j)=(3*e(1)-3*e(2)+e(3)-e(4))/24
	   utabpr(4,i,j)=(2*e(1)+e(2)-2*e(3)-e(4)+e(5)-e(6))/16
	   utabpr(5,i,j)=(2*e(1)+e(2)-2*e(3)-e(4)-3*e(5)+3*e(6))/48
	   utabpr(6,i,j)=(e(1)-e(2)-e(3)+e(4))/24	
       
       enddo
	  enddo
	  u0tab(1)=1
	  u0tab(2:6)=0	

      
     if (myrank().eq.0) then 
      open(unit=9,form='formatted',file='vtableout.dat')
      rewind 9
      do i=0,ntab
	   write(9,'(6(e15.7,2x))') vtab(:,i)  
      enddo
      close(9)    
      
	  open(unit=9,form='formatted',file='dtutab.dat')
      rewind 9
      do i=0,ntab
	  write(9,'(6(e15.7,2x))') utab(1:6,i,0)    
      enddo
      close(9)  
	  open(unit=9,form='formatted',file='dtutabpr.dat')
      rewind 9
      do i=0,ntab
	  write(9,'(6(e15.7,2x))') utabpr(1:6,i,0)    
      enddo
      close(9)       
	  write(6,*) 'u table and u pr table done!'
     endif 
      
	  return
	  end subroutine v6stepcalcinit	

	  subroutine getutab(utabout)
	  real(kind=r8) :: utabout(6,0:ntab,-mmax:mmax)
      
      
	  utabout=utab

      
	  return
	  end subroutine getutab
	  
	  subroutine getevdttab(evdttabout)
	  real(kind=r8) :: evdttabout(6,0:ntab,-mmax:mmax)
	  evdttabout=evdttab
	  return
      end subroutine getevdttab
      
      subroutine vtableout(vtabout)
	  real(kind=r8) :: vtabout(6,0:ntab)
	  vtabout=vtab
	  return
	  end subroutine vtableout
	  
  !    subroutine readvtab ! call v6potinit first
  !    integer(kind=i4) :: i,j  
  !    open(unit=9,form='formatted',file='vtable.dat')
  !    rewind 9
  !    do i=0,ntab
	 !  do j=1,6
	 !   read(9,*) vtab(j,i)
		!!write(13,*) j,i,vtab(j,i)
	 !  enddo      
  !    enddo
	 ! close(9)   
  !    return
	 ! end subroutine readvtab
	  
	    	
      subroutine spinit(npartin,nprotin) ! give initial weight of 96 states by cwtgnd
	  use math
      real(kind=r8) :: si
	  integer(kind=i4) :: npartin,nprotin,nisospinin,ineut
	  integer(kind=i4) :: iispin,ispin,nstate
      real(kind=r8), parameter :: zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0 &
                                 ,five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,nine=9.d0 &
                                 ,ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0	
      integer(kind=i4), allocatable :: ignd(:),itemp1(:),itemp2(:)
	  integer(kind=i4) :: m1=1,m2=2
	  integer(kind=i4) :: i,j,k,l,lspin,kspin,in,lx,niso,n,n1,i1,iswap 
	  npart=npartin
	  nprot=nprotin
	  nneut=npart-nprot
      nspin=2**npart ! n=4, it is 16.
      call combin(npart,nprot,nisospin) ! give value for nisospin, n=4, it is 4C2= 6.
      nbasis=nspin*nisospin ! n=4, it is 96.
	  if (.not.allocated(invspin)) allocate(invspin(nbasis))	
      if (.not.allocated(invispin)) allocate(invispin(nbasis))
      if (.not.allocated(cwt)) allocate(cwt(0:nspin-1,nisospin))
	  if (.not.allocated(invcwt)) allocate(invcwt(0:nspin-1,nisospin))
      if (.not.allocated(liso)) allocate(liso(nisospin))
	  if (.not.allocated(lisoi)) allocate(lisoi(0:nspin-1))
	  if (.not.allocated(cwtgnd)) allocate(cwtgnd(0:nspin-1,nisospin))
      if (.not.allocated(ignd)) allocate(ignd(npart))
	  if (.not.allocated(itemp1)) allocate(itemp1(npart))
	  if (.not.allocated(itemp2)) allocate(itemp2(npart))	
      
	  do iispin=1,nisospin
	   do ispin=0,nspin-1
		   nstate=(iispin-1)*nspin+ispin+1 ! 1:96
		   invcwt(ispin,iispin)=nstate
		   invspin(nstate)=ispin ! 0:15
		   invispin(nstate)=iispin ! 1:6
	   enddo
	  enddo
	  do i=1,npart
		  ignd(i)=i   ! the ignd in v6old. 1234 = pu pd nu nd
	  enddo
! spinit begin	  
	  n=npart
	  ineut=nneut
	  niso=0
! loop through isospin states and record the ones satisfying charge conservation	  
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
! put the read in state in array
      do 40 i=1,n
   40 itemp1(i)=ignd(i)
      n1=n-1
! sort array   4321 
      do 50 i=1,n1
      i1=i+1
      do 60 j=i1,n
      if (itemp1(i).gt.itemp1(j)) go to 60
      iswap=itemp1(i)
      itemp1(i)=itemp1(j)
      itemp1(j)=iswap
   60 continue
   50 continue	  
! go through allowed spin isospin states and determine which are permutations of the ground state
      do 70 l=1,niso
      lspin=shiftl(liso(l),1)
      do 80 k=1,nspin
      kspin=k-1
      cwtgnd(kspin,l)=(zero,zero)
      do 90 i=1,n
   90 itemp2(i)=and(m1,shiftr(kspin,i-1)) &
               +and(m2,shiftr(lspin,i-1))+1
! sort states and record sign of exchanges
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
! if same as ground state then fix weight
      do 120 i=1,n
      if (itemp1(i).ne.itemp2(i)) go to 130
  120 continue
      cwtgnd(kspin,l)=cmplx(si,zero)
  130 continue
   80 continue
   70 continue 	  
! test here:	  	    
	return	    
	end subroutine spinit  

	subroutine getcwtgnd(cwttmp,invspintmp,invispintmp)
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin)	
	integer(kind=i4) :: invspintmp(nbasis),invispintmp(nbasis)
	cwttmp=cwtgnd
	invspintmp=invspin
	invispintmp=invispin
	return
	end subroutine getcwtgnd
	
	subroutine getcwtgnd2(cwttmp)
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin)	
	cwttmp=cwtgnd
	return
	end subroutine getcwtgnd2
	  
	subroutine op1(i,j,nstate,cwtnew) ! op1=1
	integer(kind=i4) :: i,j,nstate,nsp,niso
	complex(kind=r8) :: cwto
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin)
	cwtnew=0
	cwtnew(invspin(nstate),invispin(nstate))=1
	end subroutine op1
		
	subroutine op2(i,j,nstate,cwtnew) ! op2 = sigma i dot sigma j
	integer(kind=i4) :: i,j,nstate,nsp,niso,sign,spn,newspn,signi,signj
	complex(kind=r8) :: cwto
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin)
	cwtnew=0
	spn=invspin(nstate) ! spn, 0:15.
	signi=2*(and(1,shiftr(spn,i-1)))-1
	signj=2*(and(1,shiftr(spn,j-1)))-1

	newspn=xor(xor(spn,shiftl(1,i-1)),shiftl(1,j-1)) ! flip i, j spin.
	cwtnew(newspn,invispin(nstate))=1 & ! six sjx	
	                               -signi*signj & ! siy sjy, -1 bc i*i
	                               +cwtnew(newspn,invispin(nstate))  
! siz sjz
	cwtnew(spn,invispin(nstate))=signi*signj+cwtnew(spn,invispin(nstate))	
	return
	end subroutine op2
! can check op2 with 2opexsigma-1. that should be slightly faster i think.
	
	subroutine op3(x,i,j,nstate,cwtnew) ! op3 = tensor op, make sure dx neq 0, otherwise need to add protection.
	integer(kind=i4) :: i,j,nstate,nsp,niso,sign,spn,newspn,newspni,newspnj,signi,signj
	complex(kind=r8) :: cwto
    real(kind=r8) :: r(3),x(3,npart),dx(3)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin)
	cwtnew=0
	spn=invspin(nstate) ! spn, 0:15.
	signi=2*(and(1,shiftr(spn,i-1)))-1
	signj=2*(and(1,shiftr(spn,j-1)))-1
	newspni=xor(spn,shiftl(1,i-1))
	newspnj=xor(spn,shiftl(1,j-1))
	newspn=xor(newspni,shiftl(1,j-1))
	dx(:)=x(:,i)-x(:,j)
	r(:)=dx(:)/sqrt(dot_product(dx,dx)) !unit vector
	!write(6,*) r
	cwtnew(newspn,invispin(nstate))=(r(1)**2 &	! six sjx
	                               -(signi*signj)*r(2)**2 & ! siy sjy
	                               +(signi+signj)*ci*r(1)*r(2))*3 & ! six sjy,  siy sjx
	                               +cwtnew(newspn,invispin(nstate))
	cwtnew(newspni,invispin(nstate))=(signj*r(1)*r(3) & ! six sjz
									+signi*ci*signj*r(2)*r(3))*3 & ! siy sjz
									+cwtnew(newspni,invispin(nstate))
	cwtnew(newspnj,invispin(nstate))=(signi*r(1)*r(3) & ! six sjz
									+signj*ci*signi*r(2)*r(3))*3 & ! siy sjz	
	                                +cwtnew(newspnj,invispin(nstate))
	cwtnew(spn,invispin(nstate))=signi*signj*r(3)**2*3+cwtnew(spn,invispin(nstate)) ! siz sjz
! -op2	
	cwtnew(newspn,invispin(nstate))=cwtnew(newspn,invispin(nstate)) &
								   -1 & ! six sjx
	                               +signi*signj  ! siy sjy 
	cwtnew(spn,invispin(nstate))=cwtnew(spn,invispin(nstate)) &
	                            -signi*signj ! siz sjz
	return
	end subroutine op3	
			
	subroutine op5(i,j,nstate,cwtnew) ! sigma i dot sigma j tau i dot tau j
	integer(kind=i4) :: i,j,nstate,nstatenew,nstatenew1,nstatenew2,nstatenew3,nsp,niso,sign,spn,newspn,signi,signj
	complex(kind=r8) :: cwto
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin)
	cwtnew=0	
	call opextau(i,j,nstate,nstatenew1)	
	call opexsigma(i,j,nstatenew1,nstatenew2)
	call opexsigma(i,j,nstate,nstatenew3)
	cwtnew(invspin(nstatenew2),invispin(nstatenew2))=4+cwtnew(invspin(nstatenew2),invispin(nstatenew2))
	cwtnew(invspin(nstatenew3),invispin(nstatenew3))=-2+cwtnew(invspin(nstatenew3),invispin(nstatenew3))
	cwtnew(invspin(nstatenew1),invispin(nstatenew1))=-2+cwtnew(invspin(nstatenew1),invispin(nstatenew1))
	cwtnew(invspin(nstate),invispin(nstate))=1+cwtnew(invspin(nstate),invispin(nstate))
	return
	end subroutine op5	  
	    
	subroutine op4(i,j,nstate,cwtnew) ! tau i dot tau j
	integer(kind=i4) :: i,j,nstate,nstatenew,nstatenew1,nstatenew2,nstatenew3,nsp,niso,sign,spn,newspn,signi,signj
	complex(kind=r8) :: cwto
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin)
	cwtnew=0	
	call opextau(i,j,nstate,nstatenew1)	
	cwtnew(invspin(nstatenew1),invispin(nstatenew1))=cwtnew(invspin(nstatenew1),invispin(nstatenew1))+2
	cwtnew(invspin(nstate),invispin(nstate))=cwtnew(invspin(nstate),invispin(nstate))-1
	return
	end subroutine op4	  
	
	subroutine op6(x,i,j,nstate,cwtnew) ! op6 = tensor  tau i dot tau j
	integer(kind=i4) :: i,j,nstate,nstatenew,nstatenew1,nstatenew2,nstatenew3,nsp,niso,sign,spn,newspn,signi,signj
	complex(kind=r8) :: cwto
    real(kind=r8) :: r(3),x(3,npart)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin)
	cwtnew=0
	call opextau(i,j,nstate,nstatenew1)	
	call op3(x,i,j,nstatenew1,cwtnew1)
	call op3(x,i,j,nstate,cwtnew2)
	cwtnew=2*cwtnew1-cwtnew2
	return
	end subroutine op6		
!   
	subroutine opexsigma(i,j,nstate,nstatenew) ! spin exchange, I think this is the fastest.
	integer(kind=i4) :: i,j,iv,jv,nstate,nstatenew,nsp,niso,sign,spn,newspn,signi,signj
	complex(kind=r8) :: cwto
	!cwtnew=0.
	spn=invspin(nstate) ! spn, 0:15.
	iv=shiftr(and(spn,shiftl(1,i-1)),i-1)
	jv=shiftr(and(spn,shiftl(1,j-1)),j-1)
	newspn=spn+(jv-iv)*2**(i-1)+(iv-jv)*2**(j-1)
	!write(6,*) and(spn,shiftl(1,i-1)),and(spn,shiftl(1,j-1)),i,j,iv,jv,newspn,spn
	!cwtnew(newspn,invispin(nstate))=1.
	nstatenew=invcwt(newspn,invispin(nstate))
	return
	end subroutine opexsigma

	subroutine opexsigma1sumbasis(i,j,cwtold,cwtnew)
	integer(kind=i4) :: i,j,is	
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin)
	cwtnew=0.
	do is=1,nbasis
		 call opexsigma1(i,j,is,cwtnew1) 
		 cwtnew=cwtnew+cwtnew1*cwtold(invspin(is),invispin(is))
	enddo 
	return
	end subroutine opexsigma1sumbasis
	
	subroutine opexsigmasumbasis(i,j,cwtold,cwtnew)
	integer(kind=i4) :: i,j,is,in
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin)
	cwtnew=0
	do is=1,nbasis
		 call opexsigma(i,j,is,in) 
		 cwtnew2=0
		 cwtnew2(invspin(in),invispin(in))=1
		 cwtnew=cwtnew+cwtnew2*cwtold(invspin(is),invispin(is))
	enddo 
	return
	end subroutine opexsigmasumbasis
	
	
	
	subroutine opexsigma0(i,j,nstate,nstatenew) ! spin exchange with if then else 
	integer(kind=i4) :: i,j,iv,jv,nstate,nstatenew,nsp,niso,sign,spn,newspn,signi,signj
	complex(kind=r8) :: cwto
	!cwtnew=0.
	spn=invspin(nstate) ! spn, 0:15.
	iv=shiftr(and(spn,shiftl(1,i-1)),i-1)
	jv=shiftr(and(spn,shiftl(1,j-1)),j-1)	
	if (iv.eq.jv) then
	!cwtnew(spn,invispin(nstate))=1.	
	nstatenew=nstate
	else
	newspn=xor(spn,shiftl(1,i-1))
	newspn=xor(newspn,shiftl(1,j-1))	
	!cwtnew(newspn,invispin(nstate))=1.	
	nstatenew=invcwt(newspn,invispin(nstate))
	endif
	return
	end subroutine opexsigma0
	
	subroutine opexsigma1(i,j,nstate,cwtnew) ! spin exchange, use op2, haven't check. may be slow.
	integer(kind=i4) :: i,j,nstate,nstatenew,nsp,niso,sign,spn,newspn,signi,signj
	complex(kind=r8) :: cwto
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin)
	!call op2(i,j,nstate,1.,cwtnew)
	cwtnew=0
	spn=invspin(nstate) ! spn, 0:15.
	signi=2*(and(1,shiftr(spn,i-1)))-1
	signj=2*(and(1,shiftr(spn,j-1)))-1
	newspn=xor(spn,shiftl(1,i-1))
	newspn=xor(newspn,shiftl(1,j-1))	
	cwtnew(newspn,invispin(nstate))=cwtnew(newspn,invispin(nstate))+1 & ! six sjx	
	                               -signi*signj ! siy sjy
! siz sjz
	cwtnew(spn,invispin(nstate))=cwtnew(spn,invispin(nstate))+signi*signj
!			
	cwtnew(spn,invispin(nstate))=cwtnew(spn,invispin(nstate))+1
	cwtnew=cwtnew/2
	return
	end subroutine opexsigma1
	
	subroutine opextau(i,j,nstate,nstatenew) ! tau exchange 
	integer(kind=i4) :: i,j,iv,jv,nstate,nstatenew,nsp,niso,sign,ispn,newispn,signi,signj
	complex(kind=r8) :: cwto
	!cwtnew=0.
	ispn=liso(invispin(nstate)) 
	iv=shiftr(and(ispn,shiftl(1,i-1)),i-1)
	jv=shiftr(and(ispn,shiftl(1,j-1)),j-1)	
	newispn=ispn+(jv-iv)*2**(i-1)+(iv-jv)*2**(j-1)
	!cwtnew(invspin(nstate),lisoi(newispn))=1.
	nstatenew=invcwt(invspin(nstate),lisoi(newispn))
	!write(6,*) invspin(nstate),lisoi(newispn),ispn,newispn
	return
	end subroutine opextau  
	  
	subroutine opextau0(i,j,nstate,nstatenew) ! tau exchange with if then else
	integer(kind=i4) :: i,j,iv,jv,nstate,nstatenew,nsp,niso,sign,ispn,newispn,signi,signj
	complex(kind=r8) :: cwto
	!cwtnew=0.
	ispn=liso(invispin(nstate)) 
	iv=shiftr(and(ispn,shiftl(1,i-1)),i-1)
	jv=shiftr(and(ispn,shiftl(1,j-1)),j-1)	
	if (iv.eq.jv) then
	!cwtnew(invspin(nstate),lisoi(ispn))=1.	
	nstatenew=nstate
	else
	newispn=xor(ispn,shiftl(1,i-1))
	newispn=xor(newispn,shiftl(1,j-1))	
	!cwtnew(invspin(nstate),lisoi(newispn))=1.
	nstatenew=invcwt(invspin(nstate),lisoi(newispn))
	endif
	return
	end subroutine opextau0	
!	
	subroutine op1sumbasis(i,j,cwtold,cwtnew)
	integer(kind=i4) :: i,j,is	
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin)
	cwtnew=0
	do is=1,nbasis
		 call op1(i,j,is,cwtnew1) 
		 cwtnew=cwtnew+cwtnew1*cwtold(invspin(is),invispin(is))
	enddo 
	return
	end subroutine op1sumbasis

	subroutine op2sumbasis(i,j,cwtold,cwtnew)
	integer(kind=i4) :: i,j,is	
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin)
	cwtnew=0
	do is=1,nbasis
		 call op2(i,j,is,cwtnew1) 
		 cwtnew=cwtnew+cwtnew1*cwtold(invspin(is),invispin(is))
	enddo 
	return
	end subroutine op2sumbasis
	
	subroutine op3sumbasis(x,i,j,cwtold,cwtnew)
	integer(kind=i4) :: i,j,is,k
    real(kind=r8) :: r(3),x(3,npart)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin) &
		               ,cwtnews(0:nspin-1,nisospin,nbasis)
	cwtnew=0
	do is=1,nbasis
		call op3(x,i,j,is,cwtnew1) 
		cwtnew=cwtnew+cwtnew1*cwtold(invspin(is),invispin(is))
		!cwtnews(:,:,is)=cwtnew1(:,:)*cwtold(invspin(is),invispin(is))
		!write(11,*) is,cwtnew(invspin(is),invispin(is)),cwtnew1(invspin(is),invispin(is)),cwtold(invspin(is),invispin(is))
	enddo 
	!do is=1,nbasis		
	!   cwtnew(:,:)=cwtnew(:,:)+cwtnews(:,:,is)
	!   write(11,*) is 
	!   k=11
	!   write(11,*) k,cwtnews(invspin(k),invispin(k),is),cwtnew(invspin(k),invispin(k))
	!enddo
	return
	end subroutine op3sumbasis	
	
	subroutine op4sumbasis(i,j,cwtold,cwtnew)
	integer(kind=i4) :: i,j,is	
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin)
	cwtnew=0
	do is=1,nbasis
		 call op4(i,j,is,cwtnew1) 
		 cwtnew=cwtnew+cwtnew1*cwtold(invspin(is),invispin(is))
	enddo 
	return
	end subroutine op4sumbasis
	
	subroutine op5sumbasis(i,j,cwtold,cwtnew)
	integer(kind=i4) :: i,j,is	
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin)
	cwtnew=0.
	do is=1,nbasis
		 call op5(i,j,is,cwtnew1) 
		 cwtnew=cwtnew+cwtnew1*cwtold(invspin(is),invispin(is))
	enddo 
	return
	end subroutine op5sumbasis
	
	subroutine op6sumbasis(x,i,j,cwtold,cwtnew)
	integer(kind=i4) :: i,j,is	
    real(kind=r8) :: r(3),x(3,npart)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin)
	cwtnew=0
	do is=1,nbasis
		 call op6(x,i,j,is,cwtnew1) 
		 cwtnew=cwtnew+cwtnew1*cwtold(invspin(is),invispin(is))
	enddo
	return
	end subroutine op6sumbasis
	
	subroutine v6prop1(x,i,j,u,cwtold,cwtnew)
	integer(kind=i4) :: i,j	
    real(kind=r8) :: r(3),x(3,npart)
	real(kind=r8) :: u(6) ! v6 so 6, coefficients after conversion. check if real or complex, I think real.
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin)
	cwtnew=0
	call op1sumbasis(i,j,cwtold,cwtnew1)
	cwtnew=cwtnew+cwtnew1*u(1)
	call op2sumbasis(i,j,cwtold,cwtnew1)
	cwtnew=cwtnew+cwtnew1*u(2)
	call op3sumbasis(x,i,j,cwtold,cwtnew1)
	cwtnew=cwtnew+cwtnew1*u(3)
	call op4sumbasis(i,j,cwtold,cwtnew1)
	cwtnew=cwtnew+cwtnew1*u(4)
	call op5sumbasis(i,j,cwtold,cwtnew1)
	cwtnew=cwtnew+cwtnew1*u(5)
	call op6sumbasis(x,i,j,cwtold,cwtnew1)
	cwtnew=cwtnew+cwtnew1*u(6)
	return
    end subroutine v6prop1

    subroutine v6prop1schmidt(x,i,j,u,cwtold,cwtnew)
	integer(kind=i4) :: i,j	
    real(kind=r8) :: r(3),x(3,npart)
	real(kind=r8) :: u(6) ! v6 so 6, coefficients after conversion. check if real or complex, I think real.
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin)
	cwtnew=0
	call op1sumbasis(i,j,cwtold,cwtnew1)
	cwtnew=cwtnew+cwtnew1*u(1)
	call op2sumbasis(i,j,cwtold,cwtnew1)
	cwtnew=cwtnew+cwtnew1*u(2)
	call op3sumbasis(x,i,j,cwtold,cwtnew1)
	cwtnew=cwtnew+cwtnew1*u(5)
	call op4sumbasis(i,j,cwtold,cwtnew1)
	cwtnew=cwtnew+cwtnew1*u(3)
	call op5sumbasis(i,j,cwtold,cwtnew1)
	cwtnew=cwtnew+cwtnew1*u(4)
	call op6sumbasis(x,i,j,cwtold,cwtnew1)
	cwtnew=cwtnew+cwtnew1*u(6)
	return
	end subroutine v6prop1schmidt
    
	subroutine vemprop1(r,i,j,it,cwtold,cwtnew)
	integer(kind=i4) :: i,j,l,it,lspin,ipi,ipj
	real(kind=r8) :: vem(14),emtot,r
	complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin)
	cwtnew=cwtold
	
	if (irep.eq.4) return
	
	do l=1,nisospin
         lspin=liso(l)       
		   ipi=(2*(and(1,shiftr(lspin,i-1))))/2
		   ipj=(2*(and(1,shiftr(lspin,j-1))))/2  
		   if ((ipi*ipj).eq.1) then
              call empot(2,r,vem) ! only pp simpliest coulomb here.
              emtot=vem(1) 
              cwtnew(:,l)=cwtold(:,l)*exp(-emtot*2.0_r8**it*dt)              
		   endif               
	enddo
	return
	end subroutine vemprop1
	
	subroutine v6propr(x,it,cwtold,cwtnew)
	integer(kind=i4) :: i,j,it,index	
	real(kind=r8) :: x(3,npart) ! the configuration for the system
    real(kind=r8), dimension(3) :: dx
    real(kind=r8) :: r,c1,c2,c3,c4,dr,tstep ! tstep=dt*2**it
	real(kind=r8) :: u(6)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin) &
		               ,cwt1(0:nspin-1,nisospin)
	
	if (irep.eq.4) then
	   cwtnew=cwtold ! because for vmc, dt=0.	
	   return	
	endif
	
	cwt1=cwtold
	
	!write(6,*) 'x=',x
	
	do i=1,npart-1
	    do j=i+1,npart 		
           dx(:)=x(:,i)-x(:,j)		   
		   r=sqrt(dot_product(dx,dx)) !r=sqrt(sum(dx(:)**2))
		   
		   if (r.eq.0) then
		     write(6,*) 'two particles have the same position, stop'
			 stop
		   endif

           if (r.lt.rangev) then
			  dr=scalep*r
			  index=dr
			  index=max(1,min(index,ntab-2))
			  dr=dr-index
              c1=-dr*(dr-1)*(dr-2)/6
              c2=(dr+1)*(dr-1)*(dr-2)/2
              c3=-(dr+1)*dr*(dr-2)/2
              c4=(dr+1)*dr*(dr-1)/6
              u(:)=c1*utab(:,index-1,it)+c2*utab(:,index,it) &
			      +c3*utab(:,index+1,it)+c4*utab(:,index+2,it)
			  
           !write(6,*) 'r,dr=',r,dr
           !write(6,*) 'dx=',dx	
           !write(6,*) 'i,j,u=',i,j,u
           !write(6,*) 'c1,2,3,4=',c1,c2,c3,c4
           !write(6,*) 'utab=',utab(:,index-1,it),utab(:,index,it) &
           !                  ,utab(:,index+1,it),utab(:,index+2,it)  
           
		   else
			  u=u0tab
		   endif
		   
		   
		   if (iem) then
               
			   !call vemprop1(r,i,j,it-1,cwt1,cwtnew1) ! EM is on.   
			   !call v6prop1(x,i,j,u,cwtnew1,cwtnew2)
			   !call vemprop1(r,i,j,it-1,cwtnew2,cwtnew) ! EM is on. 
               
               ! since vem pp commute with v6, the following is enough.
               call vemprop1(r,i,j,it,cwt1,cwtnew1) ! EM is on.   
			   call v6prop1(x,i,j,u,cwtnew1,cwtnew)
               
		   else
			   call v6prop1(x,i,j,u,cwt1,cwtnew)  
     
               ! write(6,*) 'vpropr! i=,j=,cwtnew=', i, j, sum(cwtnew)
 	   
		   endif
		   cwt1=cwtnew
		   
		   
		enddo
	enddo
    !cwtnew=cwt1
	
	
	return
    end subroutine v6propr
    
    subroutine v6proprpr(x,it,cwtold,cwtnew)
	integer(kind=i4) :: i,j,it,index	
	real(kind=r8) :: x(3,npart) ! the configuration for the system
    real(kind=r8), dimension(3) :: dx
    real(kind=r8) :: r,c1,c2,c3,c4,dr,tstep ! tstep=dt*2**it
	real(kind=r8) :: u(6)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin) &
		               ,cwt1(0:nspin-1,nisospin)
	
	if (irep.eq.4) then
	   cwtnew=cwtold ! because for vmc, dt=0.	
	   return	
	endif
	
	cwt1=cwtold
	
	!write(6,*) 'x=',x
	
	do i=1,npart-1
	    do j=i+1,npart 		
           dx(:)=x(:,i)-x(:,j)		   
		   r=sqrt(dot_product(dx,dx)) !r=sqrt(sum(dx(:)**2))
		   
		   if (r.eq.0) then
		     write(6,*) 'two particles have the same position, stop'
			 stop
		   endif

           if (r.lt.rangev) then
			  dr=scalep*r
			  index=dr
			  index=max(1,min(index,ntab-2))
			  dr=dr-index
              c1=-dr*(dr-1)*(dr-2)/6
              c2=(dr+1)*(dr-1)*(dr-2)/2
              c3=-(dr+1)*dr*(dr-2)/2
              c4=(dr+1)*dr*(dr-1)/6
              u(:)=c1*utabpr(:,index-1,it)+c2*utabpr(:,index,it) &
			      +c3*utabpr(:,index+1,it)+c4*utabpr(:,index+2,it)
			  
           !write(6,*) 'r,dr=',r,dr
           !write(6,*) 'dx=',dx	
           !write(6,*) 'i,j,u=',i,j,u
           !write(6,*) 'c1,2,3,4=',c1,c2,c3,c4
           !write(6,*) 'utab=',utab(:,index-1,it),utab(:,index,it) &
           !                  ,utab(:,index+1,it),utab(:,index+2,it)  
           
		   else
			  u=u0tab
		   endif
		   
		   
		   if (iem) then
			   !call vemprop1(r,i,j,it-1,cwt1,cwtnew1) ! EM is on.   
			   !call v6prop1(x,i,j,u,cwtnew1,cwtnew2)
			   !call vemprop1(r,i,j,it-1,cwtnew2,cwtnew) ! EM is on.   
               
			   call vemprop1(r,i,j,it,cwt1,cwtnew1) ! EM is on.   
			   call v6prop1(x,i,j,u,cwtnew1,cwtnew)
               
               
		   else
			   call v6prop1(x,i,j,u,cwt1,cwtnew)  
			   
		   endif
		   cwt1=cwtnew
		   
		   
		enddo
	enddo
    !cwtnew=cwt1
	
	
	return
	end subroutine v6proprpr
    
    
	
	subroutine v6propl(x,it,cwtold,cwtnew)
	integer(kind=i4) :: i,j,it,index	
	real(kind=r8) :: x(3,npart)
    real(kind=r8), dimension(3) :: dx
    real(kind=r8) :: r,c1,c2,c3,c4,dr,tstep ! tstep=dt*2**it
	real(kind=r8) :: u(6)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin) &
		               ,cwt1(0:nspin-1,nisospin)

	if (irep.eq.4) then
	   cwtnew=cwtold ! because for vmc, dt=0.	
	   return	
	endif	
	
	cwt1=cwtold
	!write(6,*) 'x=',x	
	
	do i=npart-1,1,-1
		do j=npart,i+1,-1 
		   dx(:)=x(:,i)-x(:,j)
		   r=sqrt(dot_product(dx,dx)) 
           if (r.lt.rangev) then
			  dr=scalep*r
			  index=dr
			  index=max(1,min(index,ntab-2))
			  dr=dr-index
              c1=-dr*(dr-1)*(dr-2)/6
              c2=(dr+1)*(dr-1)*(dr-2)/2
              c3=-(dr+1)*dr*(dr-2)/2
              c4=(dr+1)*dr*(dr-1)/6
              u(:)=c1*utab(:,index-1,it)+c2*utab(:,index,it) &
			      +c3*utab(:,index+1,it)+c4*utab(:,index+2,it)
			 		  			  
			  !write(6,*) 'r,dr=',r,dr  
			  !write(6,*) 'i,j,u=',i,j,u
     !         write(6,*) 'c1,2,3,4=',c1,c2,c3,c4
     !         write(6,*) 'dx=',dx		  
			  		  		  
		   else
			  u=u0tab
		   endif
		   
		   
		   if (iem) then
			   !call vemprop1(r,i,j,it-1,cwt1,cwtnew1) ! EM is on.   
			   !call v6prop1(x,i,j,u,cwtnew1,cwtnew2)
			   !call vemprop1(r,i,j,it-1,cwtnew2,cwtnew) ! EM is on.   
               
			   call vemprop1(r,i,j,it,cwt1,cwtnew1) ! EM is on.   
			   call v6prop1(x,i,j,u,cwtnew1,cwtnew)
		   else
			   call v6prop1(x,i,j,u,cwt1,cwtnew)  ! need to use schmidt order. 
			   
		   endif
		   cwt1=cwtnew
		   
		   
		enddo
	enddo

	
	
	return
    end subroutine v6propl
    
    
    subroutine v6proplpr(x,it,cwtold,cwtnew)
	integer(kind=i4) :: i,j,it,index	
	real(kind=r8) :: x(3,npart)
    real(kind=r8), dimension(3) :: dx
    real(kind=r8) :: r,c1,c2,c3,c4,dr,tstep ! tstep=dt*2**it
	real(kind=r8) :: u(6)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin) &
		               ,cwt1(0:nspin-1,nisospin)
	
	if (irep.eq.4) then
	   cwtnew=cwtold ! because for vmc, dt=0.	
	   return	
	endif
	
	cwt1=cwtold
	!write(6,*) 'x=',x	
	
	do i=npart-1,1,-1
		do j=npart,i+1,-1 
		   dx(:)=x(:,i)-x(:,j)
		   r=sqrt(dot_product(dx,dx)) 
           if (r.lt.rangev) then
			  dr=scalep*r
			  index=dr
			  index=max(1,min(index,ntab-2))
			  dr=dr-index
              c1=-dr*(dr-1)*(dr-2)/6
              c2=(dr+1)*(dr-1)*(dr-2)/2
              c3=-(dr+1)*dr*(dr-2)/2
              c4=(dr+1)*dr*(dr-1)/6
              u(:)=c1*utabpr(:,index-1,it)+c2*utabpr(:,index,it) &
			      +c3*utabpr(:,index+1,it)+c4*utabpr(:,index+2,it)
			 		  			  
			  !write(6,*) 'r,dr=',r,dr  
			  !write(6,*) 'i,j,u=',i,j,u
     !         write(6,*) 'c1,2,3,4=',c1,c2,c3,c4
     !         write(6,*) 'dx=',dx		  
			  		  		  
		   else
			  u=u0tab
		   endif
		   
		   
		   if (iem) then
			   !call vemprop1(r,i,j,it-1,cwt1,cwtnew1) ! EM is on.   
			   !call v6prop1(x,i,j,u,cwtnew1,cwtnew2)
			   !call vemprop1(r,i,j,it-1,cwtnew2,cwtnew) ! EM is on.   
               
			   call vemprop1(r,i,j,it,cwt1,cwtnew1) ! EM is on.   
			   call v6prop1(x,i,j,u,cwtnew1,cwtnew)
		   else
			   call v6prop1(x,i,j,u,cwt1,cwtnew)  
			   
		   endif
		   cwt1=cwtnew
		   
		   
		enddo
	enddo


	return
	end subroutine v6proplpr
    
    
	
	subroutine v6proprl(xl,xr,it,cwtold,cwtnew) ! timestep=2**it dt
	integer(kind=i4) :: it
	real(kind=r8) :: xl(3,npart),xr(3,npart) 
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin)
	if (irep.eq.4) then
	   cwtnew=cwtold ! because for vmc, dt=0.	
	   return	
	endif	
    call v6propr(xr,it,cwtold,cwtnew1)
	call v6propl(xl,it,cwtnew1,cwtnew)
	return
	end subroutine v6proprl
	
	subroutine v6proplr(xl,xr,it,cwtold,cwtnew) ! timestep=2**it dt
	integer(kind=i4) :: it
	real(kind=r8) :: xl(3,npart),xr(3,npart) 
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin)
	if (irep.eq.4) then
	   cwtnew=cwtold ! because for vmc, dt=0.	
	   return	
	endif
	call v6propl(xl,it,cwtold,cwtnew1)
    call v6propr(xr,it,cwtnew1,cwtnew)  
	return
    end subroutine v6proplr    	
    
   subroutine v6mat(x,vij,vijeigvec,vijeigval)
   use matrixmod
   use brut
   real(kind=r8) :: x(3,npart),x0(3,npart)
   complex(kind=r8) :: cwt(0:nspin-1,nisospin),cwtnew(0:nspin-1,nisospin),cwtgnd(0:nspin-1,nisospin)
   integer(kind=i4) ::  i,j,k,is,iiso 
   complex(kind=r8) :: vij(nbasis,nbasis),vijeigvec(nbasis,nbasis)
   real(kind=r8) :: vijeigval(nbasis),t 
! <i|V6|j> matrix first.  
!    invcwt 
! maybe do a x convert. 
! maybe make it 16 * 16 eigenstates. see if it can be faster.
    x0=x
    call xconvert(x0)
    call getcwtgnd2(cwtgnd) ! no need perhaps.
    do j=1,nbasis   
       call v6val1(x0,j,cwtnew) 
       do i=1,nbasis  
          vij(i,j)=cwtnew(invspin(i),invispin(i)) ! vij should be hermitian. can check.
       enddo
    enddo
! can check hermitian.
    
! diagnalize the hermitian vij
        
    call eigen(vij,vijeigvec,vijeigval,nbasis)
 
    return
    end subroutine v6mat   
    
    
    
    
	subroutine v6val(x,cwtold,cwtnew)
	integer(kind=i4) :: i,j,it,index	
	real(kind=r8) :: x(3,npart) ! the configuration for the system
    real(kind=r8), dimension(3) :: dx
    real(kind=r8) :: r,c1,c2,c3,c4,dr,tstep ! tstep=dt*2**it
	real(kind=r8) :: u(6),v(6)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin) &
		               ,cwt1(0:nspin-1,nisospin)	
	cwtnew=0
	do i=1,npart-1
		do j=i+1,npart
           dx(:)=x(:,i)-x(:,j)
		   r=sqrt(dot_product(dx,dx)) !r=sqrt(sum(dx(:)**2))
           if (r.lt.rangev) then
			  dr=scalep*r
			  index=dr
			  index=max(1,min(index,ntab-2))
			  dr=dr-index
              c1=-dr*(dr-1)*(dr-2)/6
              c2=(dr+1)*(dr-1)*(dr-2)/2
              c3=-(dr+1)*dr*(dr-2)/2
              c4=(dr+1)*dr*(dr-1)/6
              v(:)=c1*vtab(:,index-1)+c2*vtab(:,index)+c3*vtab(:,index+1)+c4*vtab(:,index+2)
		   else
			  v=0
		   endif
		   call v6prop1(x,i,j,v,cwtold,cwt1)
		   cwtnew=cwt1+cwtnew				
		enddo
	enddo
	return
	end subroutine v6val	
	
    subroutine v6val1(x,n,cwtnew) ! can be optimized to calculate for a particular state n.
	integer(kind=i4) :: i,j,n,index	
	real(kind=r8) :: x(3,npart) ! the configuration for the system
    real(kind=r8), dimension(3) :: dx
    real(kind=r8) :: r,c1,c2,c3,c4,dr,tstep ! tstep=dt*2**it
	real(kind=r8) :: u(6),v(6)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin) &
		               ,cwt1(0:nspin-1,nisospin)	
	cwtold=0
    cwtold(invspin(n),invispin(n))=1
	cwtnew=0.
	do i=1,npart-1
		do j=i+1,npart
           dx(:)=x(:,i)-x(:,j)
		   r=sqrt(dot_product(dx,dx)) !r=sqrt(sum(dx(:)**2))
           if (r.lt.rangev) then
			  dr=scalep*r
			  index=dr
			  index=max(1,min(index,ntab-2))
			  dr=dr-index
              c1=-dr*(dr-1)*(dr-2)/6
              c2=(dr+1)*(dr-1)*(dr-2)/2
              c3=-(dr+1)*dr*(dr-2)/2
              c4=(dr+1)*dr*(dr-1)/6
              v(:)=c1*vtab(:,index-1)+c2*vtab(:,index)+c3*vtab(:,index+1)+c4*vtab(:,index+2)
		   else
			  v=0.
		   endif
		   call v6prop1(x,i,j,v,cwtold,cwt1)
		   cwtnew=cwt1+cwtnew				
		enddo
	enddo
	return
	end subroutine v6val1		
	

	
	
	subroutine psitprop(x,cwtold,cwtnew) ! Schmidt order.
	use brut
	integer(kind=i4) :: i,j,it,index,k	
	real(kind=r8) :: x(3,npart) ! the configuration for the system
    real(kind=r8), dimension(3) :: dx
    real(kind=r8) :: r,c1,c2,c3,c4,dr,tstep ! tstep=dt*2**it
	real(kind=r8) :: u(6),fop(6)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		               ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin) &
		               ,cwt1(0:nspin-1,nisospin)
	integer :: ipl(6),jpl(6),ipr(6),jpr(6)
	
	cwt1=cwtold
	
	call outputordlr(ipl,jpl,ipr,jpr)
	
	do k=1,6 ! 6 is npair
		i=ipl(k)
		j=jpl(k)
	  !write(6,*) 'i=', i
	  !write(6,*) 'j=', j
		dx(:)=x(:,i)-x(:,j)		   
		r=sqrt(dot_product(dx,dx))
		call calf(r,fop)
		
		!fop=0.
		!fop(6)=1.
		
	    write(6,*) 'r=',r
		write(6,*) 'fop=',fop
		call v6prop1schmidt(x,i,j,fop,cwt1,cwtnew) 
		cwt1=cwtnew	
	enddo
	
	return
	end subroutine psitprop	
	
	
	
	
	
	
	

	
	
	
	
	
	!subroutine v6cof(r,v) ! can be defined as a precise way of calculating v(r) and u(r) if needed.
	!integer(kind=i4) :: i,j		
	!real(kind=r8) :: x(3,npart),v(6)
 !   real(kind=r8), dimension(3) :: dx
 !   real(kind=r8) :: r
	!
	!	  v(1)=
	!	  
	!	  v(2)=
	!	  
	!	  v(3)=
	!	  
	!	  v(4)=
	!	  
	!	  v(5)=
	!	  
	!	  v(6)=
	!	
	!
	!return
	!end subroutine v6cof	
	
 
	
	

	

	
	
	
	
	
	
	
 
end module v6stepcalc
