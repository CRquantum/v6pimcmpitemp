module v6stepcalc
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   complex(kind=r8), private, parameter :: czero=(0.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: ci=(0.0_r8,1.0_r8),forthci=(0.0_r8,-0.25_r8) ! forthci=1/(4i)
   complex(kind=r8), private, parameter :: cone=(1.0_r8,0.0_r8)
   integer(kind=i4), private, save :: npart,nprot,nneut,n2p
   integer(kind=i4), private, save :: lpot,lpotpr
   integer(kind=i4), private, save :: nisospin,nspin,nbasis,niso
   integer(kind=i4), private, save :: mmax,nbisect,ntab,ntabrsmall,ntabrbig,ntablag4,nrepmax
   real(kind=r8), private, save :: scalep,rangev,rangevsafe,rangevsmall,rangevbig,drlag4
   real(kind=r8), private, save :: hbar,dt
   complex(kind=r8), private, save, allocatable :: cwt(:,:),cwtgnd(:,:)  
   real(kind=r8), private, allocatable, save :: utab(:,:,:),vtab(:,:),vtabls(:,:),evdttab(:,:,:),vemtab(:) &
                                               ,utabpr(:,:,:),evdttabpr(:,:,:) &
                                               ,utabrsmall(:,:,:),utabrsmallpr(:,:,:),utabrbig(:,:,:),utabrbigpr(:,:,:) &
                                               ,vtabrsmall(:,:),vtabrbig(:,:) &
                                               ,laginterpc(:,:),v6taborigin(:,:)
   real(kind=r8), private,save, dimension(6) :: u0tab
   complex(kind=r8), private, save :: cls1 ! coefficient for LS term.
   integer(kind=i4), private, allocatable, save :: ignd(:),invspin(:),invispin(:),liso(:),lisoi(:) &
                                                  ,invcwt(:,:) ! invcwt given nspin,niso, give nstate accordingly.  
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
    subroutine v6stepcalcinit(npartin,nprotin,lpotin,dtin,hbarin,ntabin,rangevin,mmaxin,nrepmaxin,iemin,irepin)
    use brut
    use math
    use v6pot
    use mympi
    real(kind=r8) :: dtin,hbarin,rangevin,dr,drc
    integer(kind=i4) :: npartin,nprotin,mmaxin,ntabin,nrepmaxin,lpotin
    real(kind=r8) :: v(6),e(6) ! v(6) should already absorbed the stime step dt. 
    integer(kind=i4) :: i,j,ibox,irepin
    integer(kind=i4) :: myunit
    real(kind=r8) :: vfact(8)   
    logical :: iemin
    integer(kind=i4) :: iispin,ispin,nstate
    
    irep=irepin
    iem=iemin ! em force switch

    npart=npartin
    n2p=6*npart+1
    nprot=nprotin
    nneut=npart-nprot
    nspin=2**npart ! n=4, it is 16.
    call combin(npart,nprot,nisospin) ! give value for nisospin, n=4, it is 4C2= 6.
    nbasis=nspin*nisospin ! n=4, it is 96.
    allocate(invspin(nbasis))	
    allocate(invispin(nbasis))
    allocate(cwt(0:nspin-1,nisospin))
    allocate(invcwt(0:nspin-1,nisospin))
    allocate(liso(nisospin))
    allocate(lisoi(0:nspin-1))
    allocate(cwtgnd(0:nspin-1,nisospin))
    allocate(ignd(npart))
      
    !do iispin=1,nisospin
    !    do ispin=0,nspin-1
	   !     nstate=(iispin-1)*nspin+ispin+1 ! 1:96
	   !     invcwt(ispin,iispin)=nstate
	   !     invspin(nstate)=ispin ! 0:15
	   !     invispin(nstate)=iispin ! 1:6
    !    enddo
    !enddo    

    !call spinit(npartin,nprotin)
    call getinvstarrays(nspin,nisospin,invcwt,invspin,invispin)
    call spinittest(npart,nprot,nspin,nisospin,liso,lisoi,ignd,cwtgnd) 
    
    dt=dtin
    hbar=hbarin
    ntab=ntabin
    rangev=rangevin
    mmax=mmaxin
    nrepmax=nrepmaxin
    nbisect=2**mmax
    dr=rangev/ntab 
    scalep=1/dr 
    ! make a table for u for different r and dt.  time step=2**j dt.
    !write(6,*) 'u table start!'
    vfact=1
    ibox=0
    lpot=lpotin ! 3 for v6
    lpotpr=lpot

    ntabrsmall=5000 ! small v range table
    ntabrbig=5000  ! big v range table    
    
    rangevsmall=2*dr
    rangevbig=4*rangev
    
    call v6potinit(npart,ibox,-rangev,lpot,ntab,vfact,lpotpr,iem &
                   ,ntabrsmall,ntabrbig,rangevsmall,rangevbig) 
    
    if (.not.allocated(vtab)) allocate(vtab(6,0:ntab)) 
    if (.not.allocated(vtabls)) allocate(vtabls(2,0:ntab)) 
    if (.not.allocated(v6taborigin)) allocate(v6taborigin(6,0:ntab)) 

    if (.not.allocated(vtabrsmall)) allocate(vtabrsmall(6,0:ntabrsmall)) 
    if (.not.allocated(vtabrbig)) allocate(vtabrbig(6,0:ntabrbig)) 
	  
    
    call getvtab(v6taborigin)
    call getvtabschmidt(vtab,vtabrsmall,vtabrbig) ! vtab is v6 schmidt order, this is what we use in the code.

    call getvtabls(vtabls)
    
    if (iem) then
    allocate(vemtab(0:ntab)) 
    call getvemtab(vemtab) 
    endif
	  
    !vtab=0.
    !vtab(6,:)=1.	 
	    
    !write(6,*) 'get v table done!'  
    if (.not.allocated(utab)) allocate(utab(6,0:ntab,-mmax:mmax))
    if (.not.allocated(utabpr)) allocate(utabpr(6,0:ntab,-mmax:mmax)) ! utabpr deal with +delta t
    !if (.not.allocated(evdttab)) allocate(evdttab(6,0:ntab,-mmax:mmax))
    !if (.not.allocated(evdttabpr)) allocate(evdttabpr(6,0:ntab,-mmax:mmax))	
    
    if (.not.allocated(utabrsmall)) allocate(utabrsmall(6,0:ntabrsmall,-mmax:mmax))
    if (.not.allocated(utabrsmallpr)) allocate(utabrsmallpr(6,0:ntabrsmall,-mmax:mmax)) ! utabpr deal with +delta t
    if (.not.allocated(utabrbig)) allocate(utabrbig(6,0:ntabrbig,-mmax:mmax))
    if (.not.allocated(utabrbigpr)) allocate(utabrbigpr(6,0:ntabrbig,-mmax:mmax)) ! utabpr deal with +delta t    
    
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
            !evdttab(:,i,j)=e(:)	  
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
            !evdttabpr(:,i,j)=e(:)	
            utabpr(1,i,j)=(6*e(1)+3*e(2)+2*e(3)+e(4)+3*e(5)+e(6))/16
            utabpr(2,i,j)=(6*e(1)+3*e(2)+2*e(3)+e(4)-9*e(5)-3*e(6))/48
            utabpr(3,i,j)=(3*e(1)-3*e(2)+e(3)-e(4))/24
            utabpr(4,i,j)=(2*e(1)+e(2)-2*e(3)-e(4)+e(5)-e(6))/16
            utabpr(5,i,j)=(2*e(1)+e(2)-2*e(3)-e(4)-3*e(5)+3*e(6))/48
            utabpr(6,i,j)=(e(1)-e(2)-e(3)+e(4))/24	    
        enddo     
        
        do i=0,ntabrsmall
            ! we use Schmidt order, note that it is different from the order in v18pot. so need to change it.    
            v(:)=vtabrsmall(:,i)*dt*2.0_r8**j  
            e(1)=exp(-(v(1)+v(2)+2*v(3)+v(4)+v(5)+2*v(6))) !c
            e(2)=exp(-(v(1)+v(2)-4*v(3)+v(4)+v(5)-4*v(6))) !sigma
            e(3)=exp(-(v(1)+v(2)+2*v(3)-3*v(4)-3*v(5)-6*v(6))) !t	
            e(4)=exp(-(v(1)+v(2)-4*v(3)-3*v(4)-3*v(5)+12*v(6))) !tau
            e(5)=exp(-(v(1)-3*v(2)+v(4)-3*v(5))) !sigma tau
            e(6)=exp(-(v(1)-3*v(2)-3*v(4)+9*v(5))) ! t tau	
            !evdttab(:,i,j)=e(:)	  
            utabrsmall(1,i,j)=(6*e(1)+3*e(2)+2*e(3)+e(4)+3*e(5)+e(6))/16
            utabrsmall(2,i,j)=(6*e(1)+3*e(2)+2*e(3)+e(4)-9*e(5)-3*e(6))/48
            utabrsmall(3,i,j)=(3*e(1)-3*e(2)+e(3)-e(4))/24
            utabrsmall(4,i,j)=(2*e(1)+e(2)-2*e(3)-e(4)+e(5)-e(6))/16
            utabrsmall(5,i,j)=(2*e(1)+e(2)-2*e(3)-e(4)-3*e(5)+3*e(6))/48
            utabrsmall(6,i,j)=(e(1)-e(2)-e(3)+e(4))/24	
       
            v(:)=vtabrsmall(:,i)*(-dt)*2.0_r8**j  
            e(1)=exp(-(v(1)+v(2)+2*v(3)+v(4)+v(5)+2*v(6))) !c
            e(2)=exp(-(v(1)+v(2)-4*v(3)+v(4)+v(5)-4*v(6))) !sigma
            e(3)=exp(-(v(1)+v(2)+2*v(3)-3*v(4)-3*v(5)-6*v(6))) !t	
            e(4)=exp(-(v(1)+v(2)-4*v(3)-3*v(4)-3*v(5)+12*v(6))) !tau
            e(5)=exp(-(v(1)-3*v(2)+v(4)-3*v(5))) !sigma tau
            e(6)=exp(-(v(1)-3*v(2)-3*v(4)+9*v(5))) ! t tau	
            !evdttabpr(:,i,j)=e(:)	
            utabrsmallpr(1,i,j)=(6*e(1)+3*e(2)+2*e(3)+e(4)+3*e(5)+e(6))/16
            utabrsmallpr(2,i,j)=(6*e(1)+3*e(2)+2*e(3)+e(4)-9*e(5)-3*e(6))/48
            utabrsmallpr(3,i,j)=(3*e(1)-3*e(2)+e(3)-e(4))/24
            utabrsmallpr(4,i,j)=(2*e(1)+e(2)-2*e(3)-e(4)+e(5)-e(6))/16
            utabrsmallpr(5,i,j)=(2*e(1)+e(2)-2*e(3)-e(4)-3*e(5)+3*e(6))/48
            utabrsmallpr(6,i,j)=(e(1)-e(2)-e(3)+e(4))/24	    
        enddo   
        
        do i=0,ntabrbig
            ! we use Schmidt order, note that it is different from the order in v18pot. so need to change it.    
            v(:)=vtabrbig(:,i)*dt*2.0_r8**j  
            e(1)=exp(-(v(1)+v(2)+2*v(3)+v(4)+v(5)+2*v(6))) !c
            e(2)=exp(-(v(1)+v(2)-4*v(3)+v(4)+v(5)-4*v(6))) !sigma
            e(3)=exp(-(v(1)+v(2)+2*v(3)-3*v(4)-3*v(5)-6*v(6))) !t	
            e(4)=exp(-(v(1)+v(2)-4*v(3)-3*v(4)-3*v(5)+12*v(6))) !tau
            e(5)=exp(-(v(1)-3*v(2)+v(4)-3*v(5))) !sigma tau
            e(6)=exp(-(v(1)-3*v(2)-3*v(4)+9*v(5))) ! t tau	
            !evdttab(:,i,j)=e(:)	  
            utabrbig(1,i,j)=(6*e(1)+3*e(2)+2*e(3)+e(4)+3*e(5)+e(6))/16
            utabrbig(2,i,j)=(6*e(1)+3*e(2)+2*e(3)+e(4)-9*e(5)-3*e(6))/48
            utabrbig(3,i,j)=(3*e(1)-3*e(2)+e(3)-e(4))/24
            utabrbig(4,i,j)=(2*e(1)+e(2)-2*e(3)-e(4)+e(5)-e(6))/16
            utabrbig(5,i,j)=(2*e(1)+e(2)-2*e(3)-e(4)-3*e(5)+3*e(6))/48
            utabrbig(6,i,j)=(e(1)-e(2)-e(3)+e(4))/24	
       
            v(:)=vtabrbig(:,i)*(-dt)*2.0_r8**j  
            e(1)=exp(-(v(1)+v(2)+2*v(3)+v(4)+v(5)+2*v(6))) !c
            e(2)=exp(-(v(1)+v(2)-4*v(3)+v(4)+v(5)-4*v(6))) !sigma
            e(3)=exp(-(v(1)+v(2)+2*v(3)-3*v(4)-3*v(5)-6*v(6))) !t	
            e(4)=exp(-(v(1)+v(2)-4*v(3)-3*v(4)-3*v(5)+12*v(6))) !tau
            e(5)=exp(-(v(1)-3*v(2)+v(4)-3*v(5))) !sigma tau
            e(6)=exp(-(v(1)-3*v(2)-3*v(4)+9*v(5))) ! t tau	
            !evdttabpr(:,i,j)=e(:)	
            utabrbigpr(1,i,j)=(6*e(1)+3*e(2)+2*e(3)+e(4)+3*e(5)+e(6))/16
            utabrbigpr(2,i,j)=(6*e(1)+3*e(2)+2*e(3)+e(4)-9*e(5)-3*e(6))/48
            utabrbigpr(3,i,j)=(3*e(1)-3*e(2)+e(3)-e(4))/24
            utabrbigpr(4,i,j)=(2*e(1)+e(2)-2*e(3)-e(4)+e(5)-e(6))/16
            utabrbigpr(5,i,j)=(2*e(1)+e(2)-2*e(3)-e(4)-3*e(5)+3*e(6))/48
            utabrbigpr(6,i,j)=(e(1)-e(2)-e(3)+e(4))/24	    
        enddo        
 
    enddo
    u0tab(1)=1
    u0tab(2:6)=0	
    
! calculate lagrange c table. however this does not have any performance increase.
    
    ntablag4=10000
    if (.not.allocated(laginterpc)) allocate(laginterpc(4,0:ntablag4)) 
    drlag4=1.0_r8/ntablag4  
    do i=0,ntablag4    
        drc=i*drlag4
        laginterpc(1,i)=-drc*(drc-1)*(drc-2)/6
        laginterpc(2,i)=(drc+1)*(drc-1)*(drc-2)/2
        laginterpc(3,i)=-(drc+1)*drc*(drc-2)/2
        laginterpc(4,i)=(drc+1)*drc*(drc-1)/6            
    enddo
    
    if (myrank().eq.0) then 
        open(unit=9,form='formatted',file='v6tableorigin.dat')
        rewind 9
        do i=0,ntab
            write(9,'(7(e15.7,2x))') i*dr,v6taborigin(:,i)  
        enddo
        close(9) 
        
        open(unit=9,form='formatted',file='vlstableorigin.dat')
        rewind 9
        do i=0,ntab
            write(9,'(7(e15.7,2x))') i*dr,vtabls(:,i)  
        enddo
        close(9)
        
        open(newunit=myunit,form='formatted',file='vtableSchmidt.dat')
        rewind myunit
        do i=0,ntab
	        write(myunit,'(7(e15.7,2x))') i*dr,vtab(:,i)
        enddo
        close(myunit) 
       
        open(newunit=myunit,form='formatted',file='vrsmalltableSchmidt.dat')
        rewind myunit
        do i=0,ntabrsmall
	        write(myunit,'(6(e15.7,2x))') vtabrsmall(:,i)
        enddo
        close(myunit)        
       
        open(newunit=myunit,form='formatted',file='vrbigtableSchmidt.dat')
        rewind myunit
        do i=0,ntabrbig
	        write(myunit,'(6(e15.7,2x))') vtabrbig(:,i)
        enddo
        close(myunit)          

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
        
        open(unit=9,form='formatted',file='lag4pctable.dat')
        rewind 9
            do i=0,ntablag4
        write(9,'(4(e15.7,2x))') laginterpc(:,i)
        enddo
        close(9)         
        
        write(6,*) 'u, v tables and pr tables done!'
    endif     
    
    
    

    cls1=(-1.0_r8/(2*hbar))/(4*ci)
    
    
    return
    end subroutine v6stepcalcinit	

    
    subroutine v2ucalc(r,sgn,it,u) ! convert v(6) to u(6) if sgn=1. v(6) to upr(6) if sgn=-1. 
    use v6pot
    use cheft
    real(kind=r8) :: r,u(6),v(6),vv(18),e(6) 
    integer(kind=i4) :: it,sgn
    
    select case (lpot)
        case (1:7)
            call av18op(lpot,r,vv) ! vv is v6'
        case (112:114,122:124)
            call cheft_pot(lpot,r,vv)
        case default
            u(:)=u0tab(:)
            return
    end select  ! can just add switch in the input file.    
    
    call v6schmidtord(vv,v)
    v(1:6)=v(1:6)*sgn*dt*2.0_r8**it 
    e(1)=exp(-(v(1)+v(2)+2*v(3)+v(4)+v(5)+2*v(6))) !c
    e(2)=exp(-(v(1)+v(2)-4*v(3)+v(4)+v(5)-4*v(6))) !sigma
    e(3)=exp(-(v(1)+v(2)+2*v(3)-3*v(4)-3*v(5)-6*v(6))) !t	
    e(4)=exp(-(v(1)+v(2)-4*v(3)-3*v(4)-3*v(5)+12*v(6))) !tau
    e(5)=exp(-(v(1)-3*v(2)+v(4)-3*v(5))) !sigma tau
    e(6)=exp(-(v(1)-3*v(2)-3*v(4)+9*v(5))) ! t tau	    
    
    u(1)=(6*e(1)+3*e(2)+2*e(3)+e(4)+3*e(5)+e(6))/16
    u(2)=(6*e(1)+3*e(2)+2*e(3)+e(4)-9*e(5)-3*e(6))/48
    u(3)=(3*e(1)-3*e(2)+e(3)-e(4))/24
    u(4)=(2*e(1)+e(2)-2*e(3)-e(4)+e(5)-e(6))/16
    u(5)=(2*e(1)+e(2)-2*e(3)-e(4)-3*e(5)+3*e(6))/48
    u(6)=(e(1)-e(2)-e(3)+e(4))/24    
    
    return
    end subroutine

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
	  
	    	
    !subroutine spinit(npartin,nprotin) ! give initial weight of 96 states by cwtgnd
    !use math
    !real(kind=r8) :: si
    !integer(kind=i4) :: npartin,nprotin,ineut
    !integer(kind=i4) :: iispin,ispin,nstate
    !real(kind=r8), parameter :: zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0 &
    !                            ,five=5.d0,six=6.d0,seven=7.d0,eight=8.d0,nine=9.d0 &
    !                            ,ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0	
    !integer(kind=i4), allocatable :: ignd(npartin),itemp1(npartin),itemp2(npartin)
    !integer(kind=i4) :: m1=1,m2=2
    !integer(kind=i4) :: i,j,k,l,lspin,kspin,in,lx,niso,n,n1,i1,iswap 
    !
    !do i=1,npart
	   ! ignd(i)=i   ! the ignd in v6old. 1234 = pu pd nu nd
    !enddo
    !! spinit begin	  
    !n=npart
    !ineut=nneut
    !niso=0
    !! loop through isospin states and record the ones satisfying charge conservation	  
    !do 20 l=1,nspin
    !lspin=l-1
    !in=0
    !do 30 j=1,n
    !lx=and(m1,shiftr(lspin,j-1))
    !if (lx.eq.1) in=in+1
    !30 continue
    !if (in.ne.ineut) go to 20
    !niso=niso+1
    !liso(niso)=lspin
    !lisoi(lspin)=niso
    !20 continue
    !
    !
    !! put the read in state in array
    !do 40 i=1,n
    !40 itemp1(i)=ignd(i)
    !n1=n-1
    !! sort array   4321 
    !do 50 i=1,n1
    !i1=i+1
    !do 60 j=i1,n
    !if (itemp1(i).gt.itemp1(j)) go to 60
    !iswap=itemp1(i)
    !itemp1(i)=itemp1(j)
    !itemp1(j)=iswap
    !60 continue
    !50 continue	  
    !! go through allowed spin isospin states and determine which are permutations of the ground state
    !do 70 l=1,niso
    !lspin=shiftl(liso(l),1)
    !do 80 k=1,nspin
    !kspin=k-1
    !cwtgnd(kspin,l)=(zero,zero)
    !do 90 i=1,n
    !90 itemp2(i)=and(m1,shiftr(kspin,i-1)) &
    !        +and(m2,shiftr(lspin,i-1))+1
    !! sort states and record sign of exchanges
    !si=one
    !do 100 i=1,n1
    !i1=i+1
    !do 110 j=i1,n
    !if (itemp2(i).gt.itemp2(j)) go to 110
    !si=-si
    !iswap=itemp2(i)
    !itemp2(i)=itemp2(j)
    !itemp2(j)=iswap
    !110 continue
    !100 continue
    !! if same as ground state then fix weight
    !do 120 i=1,n
    !if (itemp1(i).ne.itemp2(i)) go to 130
    !120 continue
    !cwtgnd(kspin,l)=cmplx(si,zero)
    !130 continue
    !80 continue
    !70 continue 	  
    !! test here:	
    !return	    
    !end subroutine spinit  

    subroutine getcwtgnd(cwttmp,invspintmp,invispintmp)
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin)	
    integer(kind=i4) :: invspintmp(nbasis),invispintmp(nbasis)
    cwttmp(:,:)=cwtgnd(:,:)
    invspintmp(:)=invspin(:)
    invispintmp(:)=invispin(:)
    return
    end subroutine getcwtgnd
	
    subroutine getcwtgnd2(cwttmp)
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin)	
    cwttmp(:,:)=cwtgnd(:,:)
    return
    end subroutine getcwtgnd2
    
    subroutine getliso(lisotmp,lisoitmp)
    integer(kind=i4) :: lisotmp(nisospin),lisoitmp(0:nspin-1)
    lisotmp(:)=liso(:)
    lisoitmp(:)=lisoi(:)
    return
    end subroutine getliso   
	  
    subroutine op3pure(r,i,j,nstate,feff,cwtnew) ! op3pure: tensor witout op2 and without factor 3.
    integer(kind=i4) :: i,j,nstate,spn,ispn,newspn,newspni,newspnj,signi,signj,signij
    real(kind=r8) :: r(3)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),feff,cir2
    spn=invspin(nstate) ! spn, 0:15.
    ispn=invispin(nstate)
    !signi=2*(and(1,shiftr(spn,i-1)))-1
    !signj=2*(and(1,shiftr(spn,j-1)))-1  
    signi=2*ibits(spn,i-1,1)-1
    signj=2*ibits(spn,j-1,1)-1    
    newspni=xor(spn,shiftl(1,i-1))
    newspnj=xor(spn,shiftl(1,j-1))
    newspn=xor(newspni,shiftl(1,j-1))
    signij=signi*signj
    cir2=ci*r(2)
    cwtnew(newspn,ispn)=feff*(r(1)**2-r(2)*(signij*r(2)-(signi+signj)*ci*r(1)))+cwtnew(newspn,ispn)    
    cwtnew(newspni,ispn)=feff*(signj*r(3)*(r(1)+signi*cir2))+cwtnew(newspni,ispn)
    cwtnew(newspnj,ispn)=feff*(signi*r(3)*(r(1)+signj*cir2))+cwtnew(newspnj,ispn)    
    cwtnew(spn,ispn)=feff*signij*r(3)**2+cwtnew(spn,ispn)       
    !cwtnew(newspn,ispn)=feff*(r(1)**2-r(2)*((signi*signj)*r(2)-(signi+signj)*ci*r(1)))+cwtnew(newspn,ispn)
    !cwtnew(newspni,ispn)=feff*(signj*r(3)*(r(1)+signi*ci*r(2)))+cwtnew(newspni,ispn)
    !cwtnew(newspnj,ispn)=feff*(signi*r(3)*(r(1)+signj*ci*r(2)))+cwtnew(newspnj,ispn)
    !cwtnew(spn,ispn)=feff*signi*signj*r(3)**2+cwtnew(spn,ispn) 
    return
    end subroutine op3pure    
  
    subroutine opexsigma(i,j,nstate,nstatenew) ! spin exchange, I think this is the fastest.
    integer(kind=i4) :: i,j,nstate,nstatenew,spn,newspn
    spn=invspin(nstate) ! spn, 0:15.
    
    !iv=ibits(spn,i-1,1)
    !jv=ibits(spn,j-1,1)
    
    !iv=shiftr(and(spn,shiftl(1,i-1)),i-1)
    !jv=shiftr(and(spn,shiftl(1,j-1)),j-1)
    
    !newspn=spn+(jv-iv)*2**(i-1)+(iv-jv)*2**(j-1)

    newspn=spn
    call mvbits(spn,i-1,1,newspn,j-1)
    call mvbits(spn,j-1,1,newspn,i-1)    
    nstatenew=invcwt(newspn,invispin(nstate))
    return
    end subroutine opexsigma
 
    subroutine opextau(i,j,nstate,nstatenew) ! tau exchange 
    integer(kind=i4) :: i,j,nstate,nstatenew,ispn,newispn
    ispn=liso(invispin(nstate)) 
    
    !iv=ibits(ispn,i-1,1)
    !jv=ibits(ispn,j-1,1)    

    !iv=shiftr(and(ispn,shiftl(1,i-1)),i-1)
    !jv=shiftr(and(ispn,shiftl(1,j-1)),j-1)
    
    !newispn=ispn+(jv-iv)*2**(i-1)+(iv-jv)*2**(j-1)

    newispn=ispn
    call mvbits(ispn,i-1,1,newispn,j-1)
    call mvbits(ispn,j-1,1,newispn,i-1)
    nstatenew=invcwt(invspin(nstate),lisoi(newispn))
    return
    end subroutine opextau     
    
    subroutine opexsigmatau(i,j,nstate,nstatenew) ! sigma, tau exchange at the same time. 
    integer(kind=i4) :: i,j,nstate,nstatenew,ispn,newispn,spn,newspn
    spn=invspin(nstate) ! spn, 0:15.
    newspn=spn
    call mvbits(spn,i-1,1,newspn,j-1)
    call mvbits(spn,j-1,1,newspn,i-1)     
    ispn=liso(invispin(nstate)) 
    newispn=ispn
    call mvbits(ispn,i-1,1,newispn,j-1)
    call mvbits(ispn,j-1,1,newispn,i-1)   
    nstatenew=invcwt(newspn,lisoi(newispn))
    return
    end subroutine opexsigmatau    
  
    subroutine v6prop1(x,i,j,u,cwtold,cwtnew) ! optimized version, 30+ times faster than old verison.
    integer(kind=i4) :: i,j,nstate,nstatenew,nstatenew1,nstatenew2
    
    real(kind=r8) :: x(3,npart),r(3),rsq(3),dx(3),r3r1
    real(kind=r8) :: u(6) ! v6 so 6, coefficients after conversion. check if real or complex, I think real.
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
                       ,fceff,fsex,fiex,fsiex,fteff,fteffiex,f(6) &
                       ,cir2,cir1,c(4),cfteff(4),cfteffiex(4),r2cir1,r3cir2
    integer(kind=i4) :: spn,ispn,newspn,newispn,newspni,newspnj,signi,signj,signij  
    
    if (irep.eq.4) then
        cwtnew(:,:)=cwtold(:,:)
        return  
    endif
     
    dx(:)=x(:,i)-x(:,j)
    r(:)=dx(:)/sqrt(dot_product(dx,dx)) !unit vector
    rsq(:)=r(:)**2
    r3r1=r(3)*r(1)
    cir1=ci*r(1)
    cir2=ci*r(2)
    r2cir1=r(2)*cir1
    r3cir2=r(3)*cir2
    cwtnew(:,:)=0 
    do nstate=1,nbasis   
        f(:)=u(:)*cwtold(invspin(nstate),invispin(nstate))
        fceff=f(1)-f(2)-f(4)+f(5)+f(3)-f(6)
        fsex=2*(f(2)-f(5)-f(3)+f(6))
        fiex=2*(f(4)-f(5)+f(6))
        fsiex=4*(f(5)-f(6))
        fteff=3*(f(3)-f(6))
        fteffiex=6*f(6)
        ! fceff,op1
        cwtnew(invspin(nstate),invispin(nstate))=cwtnew(invspin(nstate),invispin(nstate))+fceff
        ! fsex,opexsigma     
        call opexsigma(i,j,nstate,nstatenew1)
        cwtnew(invspin(nstatenew1),invispin(nstatenew1))=cwtnew(invspin(nstatenew1),invispin(nstatenew1))+fsex
        ! fiex,opextau        
        call opextau(i,j,nstate,nstatenew2)  
        cwtnew(invspin(nstatenew2),invispin(nstatenew2))=cwtnew(invspin(nstatenew2),invispin(nstatenew2))+fiex        
        ! fsiex
        call opextau(i,j,nstatenew1,nstatenew)  
        cwtnew(invspin(nstatenew),invispin(nstatenew))=cwtnew(invspin(nstatenew),invispin(nstatenew))+fsiex      
        ! fteff,op3 
        !call op3pure(r,i,j,nstate,fteff,cwtnew) 
        spn=invspin(nstate) ! spn, 0:15.
        ispn=invispin(nstate) 
        signi=2*ibits(spn,i-1,1)-1
        signj=2*ibits(spn,j-1,1)-1    
        newspni=xor(spn,shiftl(1,i-1))
        newspnj=xor(spn,shiftl(1,j-1))
        newspn=xor(newspni,shiftl(1,j-1))
        signij=signi*signj
        !c(1)=rsq(1)-r(2)*(signij*r(2)-(signi+signj)*cir1) 
        c(1)=rsq(1)-rsq(2)*signij+(signi+signj)*r2cir1
        !c(2)=signj*r(3)*(r(1)+signi*cir2)
        c(2)=signj*r3r1+signij*r3cir2
        !c(3)=signi*r(3)*(r(1)+signj*cir2)
        c(3)=signi*r3r1+signij*r3cir2
        c(4)=signij*rsq(3) 
        cfteff(:)=fteff*c(:)
        cfteffiex(:)=fteffiex*c(:)   
        cwtnew(newspn,ispn)=cfteff(1)+cwtnew(newspn,ispn)    
        cwtnew(newspni,ispn)=cfteff(2)+cwtnew(newspni,ispn)
        cwtnew(newspnj,ispn)=cfteff(3)+cwtnew(newspnj,ispn)    
        cwtnew(spn,ispn)=cfteff(4)+cwtnew(spn,ispn)        
        ! fteffiex,opextau,op3
        !call opextau(i,j,nstate,nstatenew)
        !call op3pure(r,i,j,nstatenew2,fteffiex,cwtnew)    
        ispn=invispin(nstatenew2)      
        !spn=invspin(nstatenew2) ! spn, 0:15.
        !signi=2*ibits(spn,i-1,1)-1
        !signj=2*ibits(spn,j-1,1)-1  
        !newspni=xor(spn,shiftl(1,i-1))
        !newspnj=xor(spn,shiftl(1,j-1))
        !newspn=xor(newspni,shiftl(1,j-1))  
        cwtnew(newspn,ispn)=cfteffiex(1)+cwtnew(newspn,ispn)    
        cwtnew(newspni,ispn)=cfteffiex(2)+cwtnew(newspni,ispn)
        cwtnew(newspnj,ispn)=cfteffiex(3)+cwtnew(newspnj,ispn)    
        cwtnew(spn,ispn)=cfteffiex(4)+cwtnew(spn,ispn)         
    enddo 
    return
    end subroutine v6prop1         
   
    subroutine vemprop1(expvemcdt,i,j,cwtold,cwtnew)
    integer(kind=i4) :: i,j,l,it,lspin,ipi,ipj,index
    real(kind=r8) :: expvemcdt  
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin)
   
    cwtnew(:,:)=cwtold(:,:)
    if (irep.eq.4) return    
    
    do l=1,nisospin
        lspin=liso(l)       
	    ipi=and(1,shiftr(lspin,i-1))  !(2*(and(1,shiftr(lspin,i-1))))/2
	    ipj=and(1,shiftr(lspin,j-1))    !(2*(and(1,shiftr(lspin,j-1))))/2  
	    if ((ipi*ipj).eq.1) then
            cwtnew(:,l)=cwtold(:,l)*expvemcdt        
	    endif               
    enddo
    return
    end subroutine vemprop1
    
    subroutine vlsprop1(x0,x1,i,j,vls,sign,cwtold,cwtnew) 
    integer(kind=i4) :: i,j,nstate,nstatenew,nstatenew1,nstatenew2,sign
    real(kind=r8) :: x0(3,npart),x1(3,npart),rij(3),dr01i(3),dr01j(3),dr01ij(3)
    real(kind=r8) :: vls 
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
                       ,fx,fyi,fz,fold,cls1vls
    integer(kind=i4) :: spn,ispn,newspn,newispn,newspni,newspnj,signi,signj,signij  
    
    if (irep.eq.4) then
        cwtnew(:,:)=cwtold(:,:)
        return  
    endif
    
    rij(:)=x0(:,i)-x0(:,j)
    
    dr01i(:)=x0(:,i)-x1(:,i)
    dr01j(:)=x0(:,j)-x1(:,j) 
    dr01ij(:)=dr01i(:)-dr01j(:)
    
    cls1vls=cls1*vls*sign
    fx=(rij(2)*dr01ij(3)-rij(3)*dr01ij(2))*cls1vls
    fyi=-(rij(1)*dr01ij(3)-rij(3)*dr01ij(1))*cls1vls*ci
    fz=(rij(1)*dr01ij(2)-rij(2)*dr01ij(1))*cls1vls
      
    cwtnew(:,:)=0 
    do nstate=1,nbasis 
        spn=invspin(nstate) ! spn, 0:15.
        ispn=invispin(nstate) 
        fold=cwtold(spn,ispn)
       
        signi=2*ibits(spn,i-1,1)-1
        signj=2*ibits(spn,j-1,1)-1 
       
        newspni=xor(spn,shiftl(1,i-1)) ! flip spin i
        newspnj=xor(spn,shiftl(1,j-1)) ! flip spin j
        
        cwtnew(newspni,ispn)=cwtnew(newspni,ispn)+(fx+signi*fyi)*fold 
        cwtnew(newspnj,ispn)=cwtnew(newspnj,ispn)+(fx+signj*fyi)*fold  
        cwtnew(spn,ispn)=cwtnew(spn,ispn)+fz*(signi+signj)*fold                   
    enddo 
    
    return
    end subroutine vlsprop1   
       
    subroutine vlsprop0(spinorbit,sign,x0,x1,cwtold,cwtnew) ! without no LS*t1*t2. old version and works. 2020/5/20
! only LS for now. no LS*t1*t2. Linear terms. If cwtold and cwtl2r, then sign=-1. If cwtr2l, sign=1.
    use cheft
    integer(kind=i4) :: i,j,it,index,indexlag4,sign,nstate
    real(kind=r8) :: x0(3,npart),x1(3,npart) ! the configuration for the system
    real(kind=r8) :: r,c1,c2,c3,c4,dr,rij(3),dr01i(3),dr01j(3),dr01ij(3)
    real(kind=r8) :: vls
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		                ,cwtnew1(0:nspin-1,nisospin),cwt1(0:nspin-1,nisospin)
    complex(kind=r8) :: fx,fyi,fz,fold,cls1vls
    integer(kind=i4) :: spn,ispn,newspn,newispn,newspni,newspnj,signi,signj,signij  
    logical :: spinorbit
    if (irep.eq.4) then
	    cwtnew(:,:)=cwtold(:,:) ! because for vmc, dt=0.	
	    return	
    endif
    
    cwtnew(:,:)=cwtold(:,:)
    if (spinorbit) then    
        do i=1,npart-1
	        do j=i+1,npart 		
                rij(:)=x0(:,i)-x0(:,j)		   
		        r=sqrt(dot_product(rij,rij)) !r=sqrt(sum(dx(:)**2))
                if (r.lt.rangev) then
			        dr=scalep*r
			        index=dr
			        index=max(1,min(index,ntab-2)) 
			        dr=dr-index
                    c1=-dr*(dr-1)*(dr-2)/6
                    c2=(dr+1)*(dr-1)*(dr-2)/2
                    c3=-(dr+1)*dr*(dr-2)/2
                    c4=(dr+1)*dr*(dr-1)/6
                    vls=c1*vtabls(1,index-1)+c2*vtabls(1,index)+c3*vtabls(1,index+1)+c4*vtabls(1,index+2)      
                else
                    !vls=0
                    call cheft_pot_ls(lpot,r,vls)
                endif	
                dr01i(:)=x0(:,i)-x1(:,i)
                dr01j(:)=x0(:,j)-x1(:,j) 
                dr01ij(:)=dr01i(:)-dr01j(:)  
                cls1vls=cls1*vls*sign
                fx=(rij(2)*dr01ij(3)-rij(3)*dr01ij(2))*cls1vls
                fyi=-(rij(1)*dr01ij(3)-rij(3)*dr01ij(1))*cls1vls*ci
                fz=(rij(1)*dr01ij(2)-rij(2)*dr01ij(1))*cls1vls    
                do spn=0,nspin-1      
                    signi=2*ibits(spn,i-1,1)-1
                    signj=2*ibits(spn,j-1,1)-1     
                    newspni=xor(spn,shiftl(1,i-1)) ! flip spin i
                    newspnj=xor(spn,shiftl(1,j-1)) ! flip spin j
                    cwtnew(newspni,:)=cwtnew(newspni,:)-(fx+signi*fyi)*cwtold(spn,:)
                    cwtnew(newspnj,:)=cwtnew(newspnj,:)-(fx+signj*fyi)*cwtold(spn,:)
                    cwtnew(spn,:)=cwtnew(spn,:)-fz*(signi+signj)*cwtold(spn,:)           
                enddo                 
                !call vlsprop1(x0,x1,i,j,vls,sign,cwtold,cwtnew1)                   
		        !cwtnew(:,:)=cwtnew(:,:)+cwtnew1(:,:)	
	        enddo
        enddo
    endif
    return
    end subroutine vlsprop0
    
    subroutine vlsprop(spinorbit,sign,x0,x1,cwtold,cwtnew) ! added LS*t1*t2. 2020/5/20
! only LS for now. and LS*t1*t2. Linear terms. If cwtold and cwtl2r, then sign=-1. If cwtr2l, sign=1.
    use cheft
    integer(kind=i4) :: i,j,it,index,indexlag4,sign,nstate
    real(kind=r8) :: x0(3,npart),x1(3,npart) ! the configuration for the system
    real(kind=r8) :: r,c1,c2,c3,c4,dr,rij(3),dr01i(3),dr01j(3),dr01ij(3),rdx,rdy,rdz
    real(kind=r8) :: vls(2),vv(18)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		                ,cwtnew1(0:nspin-1,nisospin),cwt1(0:nspin-1,nisospin)
    complex(kind=r8) :: fx,fyi,fz,f2x,f2yi,f2z,fold,cls1vls,cls1vls2
    integer(kind=i4) :: spn,ispn,newspn,newispn,newspni,newspnj,signi,signj,signij,nsnew
    logical :: spinorbit
    if (irep.eq.4) then
	    cwtnew(:,:)=cwtold(:,:) ! because for vmc, dt=0.	
	    return	
    endif
       
    if (spinorbit) then
        if (lpot>2) then
            call vlsprop0(spinorbit,sign,x0,x1,cwtold,cwtnew)
            return
        endif     
        ! AV8' added LS*t1*t2
        cwtnew(:,:)=cwtold(:,:)
        do i=1,npart-1 ! all the rij thing for a bead can be optimized, no need to repeatedly calculate them.
	        do j=i+1,npart 		
                rij(:)=x0(:,i)-x0(:,j)		   
		        r=sqrt(dot_product(rij,rij)) !r=sqrt(sum(dx(:)**2))
                if (r.lt.rangev) then
			        dr=scalep*r
			        index=dr
			        index=max(1,min(index,ntab-2)) 
			        dr=dr-index
                    c1=-dr*(dr-1)*(dr-2)/6
                    c2=(dr+1)*(dr-1)*(dr-2)/2
                    c3=-(dr+1)*dr*(dr-2)/2
                    c4=(dr+1)*dr*(dr-1)/6
                    vls(:)=c1*vtabls(:,index-1)+c2*vtabls(:,index)+c3*vtabls(:,index+1)+c4*vtabls(:,index+2)       
                else
                    !vls(:)=0              
                    call av18op(lpot,r,vv)
                    vls(1:2)=vv(7:8)
                endif	
                         
                !write(6,*) 'vls in prop=',r, vls(1), vls(2)
                
                dr01i(:)=x0(:,i)-x1(:,i)
                dr01j(:)=x0(:,j)-x1(:,j) 
                dr01ij(:)=dr01i(:)-dr01j(:)  

                rdx=rij(2)*dr01ij(3)-rij(3)*dr01ij(2)
                rdy=-(rij(1)*dr01ij(3)-rij(3)*dr01ij(1))
                rdz=rij(1)*dr01ij(2)-rij(2)*dr01ij(1)                
                
                cls1vls=cls1*(vls(1)-vls(2))*sign
                cls1vls2=cls1*vls(2)*2*sign
                         
                fx=rdx*cls1vls
                fyi=rdy*cls1vls*ci
                fz=rdz*cls1vls 
                
                f2x=rdx*cls1vls2
                f2yi=rdy*cls1vls2*ci
                f2z=rdz*cls1vls2
                             
                !!!!! check
                !write (6,'( ''rdx, rdy, rdz = '', 3(g15.7,1x)  )')rdx,rdy,rdz
                !write (6,'( ''vls1,2 = '', 3(g15.7,1x)  )') vls(1:2),vls(1)-vls(2)  
                !write (6,'( ''cls1vls = '', 2g15.7  )') cls1vls
                !write (6,'( ''cls1vls2 = '', 2g15.7  )') cls1vls2 
                !!!!!
                          
                do ispn=1,nisospin
                    call opextau(i,j,invcwt(0,ispn),nsnew)        ! take spn as 0 here bc we only care about ispin exchange.             
                    newispn=invispin(nsnew)                      
                    do spn=0,nspin-1                      
                        fold=cwtold(spn,ispn)
                        signi=2*ibits(spn,i-1,1)-1
                        signj=2*ibits(spn,j-1,1)-1     
                        newspni=xor(spn,shiftl(1,i-1)) ! flip spin i
                        newspnj=xor(spn,shiftl(1,j-1)) ! flip spin j
                        cwtnew(newspni,ispn)=cwtnew(newspni,ispn)-(fx+signi*fyi)*fold
                        cwtnew(newspnj,ispn)=cwtnew(newspnj,ispn)-(fx+signj*fyi)*fold
                        cwtnew(spn,ispn)=cwtnew(spn,ispn)-fz*(signi+signj)*fold
                                     
                    !call opextau(i,j,invcwt(spn,ispn),nsnew)                    
                    !newispn=invispin(nsnew) 
                        
                        !write(6,*) 'isospin in prop:',i,j,ispn,newispn,nisospin
                        
                        cwtnew(newspni,newispn)=cwtnew(newspni,newispn)-(f2x+signi*f2yi)*fold
                        cwtnew(newspnj,newispn)=cwtnew(newspnj,newispn)-(f2x+signj*f2yi)*fold
                        cwtnew(spn,newispn)=cwtnew(spn,newispn)-f2z*(signi+signj)*fold                                
                    enddo   
                enddo                 
                !call vlsprop1(x0,x1,i,j,vls,sign,cwtold,cwtnew1)                   
		        !cwtnew(:,:)=cwtnew(:,:)+cwtnew1(:,:)	
	        enddo
        enddo
    else
        cwtnew(:,:)=cwtold(:,:)
    endif
    return
    end subroutine vlsprop
  
    subroutine projectornp(cwtold,cwtnnew,cwtpnew) ! (1 +/- tau_{iz})/2 the projection operator
    integer(kind=i4) :: i,j,l,lspin,ipi,ipj
    complex(kind=r8) :: cwtnnew(0:nspin-1,nisospin,npart),cwtpnew(0:nspin-1,nisospin,npart) &
                       ,cwtold(0:nspin-1,nisospin)
    cwtnnew=0
    cwtpnew=0
    do l=1,nisospin
        lspin=liso(l)     
        do i=1,npart           
		    ipi=and(1,shiftr(lspin,i-1))    
            if (ipi.eq.1) then ! proton
                cwtpnew(:,l,i)=cwtold(:,l)
            else ! neutron
                cwtnnew(:,l,i)=cwtold(:,l)
            endif
        enddo       
    enddo
    return
    end subroutine projectornp
    
    subroutine v6propr(x,it,cwtold,cwtnew)
    integer(kind=i4) :: i,j,it,index,indexlag4	
    real(kind=r8) :: x(3,npart) ! the configuration for the system
    real(kind=r8), dimension(3) :: dx
    real(kind=r8) :: r,c1,c2,c3,c4,dr
    real(kind=r8) :: u(6),utmp(6)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		                ,cwtnew1(0:nspin-1,nisospin),cwt1(0:nspin-1,nisospin)
    real(kind=r8) :: expvemcdt,vem(14)	
    if (irep.eq.4) then
	    cwtnew(:,:)=cwtold(:,:) ! because for vmc, dt=0.	
	    return	
    endif
	
    cwt1(:,:)=cwtold(:,:)
	
    !write(6,*) 'x=',x
	
    do i=1,npart-1
	    do j=i+1,npart 		
            dx(:)=x(:,i)-x(:,j)		   
		    r=sqrt(dot_product(dx,dx)) !r=sqrt(sum(dx(:)**2))
		   
		    !if (r.eq.0) then
		    !    write(6,*) 'two particles have the same position, stop'
			   ! stop
		    !endif

            if (r.lt.rangev) then
			    dr=scalep*r
			    index=dr
			    index=max(1,min(index,ntab-2)) 
			    dr=dr-index

                c1=-dr*(dr-1)*(dr-2)/6
                c2=(dr+1)*(dr-1)*(dr-2)/2
                c3=-(dr+1)*dr*(dr-2)/2
                c4=(dr+1)*dr*(dr-1)/6
 
                !indexlag4=nint(dr/drlag4)
                !c1=laginterpc(1,indexlag4)
                !c2=laginterpc(2,indexlag4)
                !c3=laginterpc(3,indexlag4)
                !c4=laginterpc(4,indexlag4) 

                u(:)=c1*utab(:,index-1,it)+c2*utab(:,index,it) &
			        +c3*utab(:,index+1,it)+c4*utab(:,index+2,it)
			  
            !write(6,*) 'r,dr=',r,dr
            !write(6,*) 'dx=',dx	
            !write(6,*) 'i,j,u=',i,j,u
            !write(6,*) 'c1,2,3,4=',c1,c2,c3,c4
            !write(6,*) 'utab=',utab(:,index-1,it),utab(:,index,it) &
            !                  ,utab(:,index+1,it),utab(:,index+2,it)  
           
            else
                call v2ucalc(r,1,it,u)
			    !u(:)=u0tab(:)
            endif
		   
            !!!!!!!!!!!!!!!! check
            !call v2ucalc(r,1,it,utmp)
            !write(6,'(''u6 ='', t20, 6(g15.7,1x) )') u(1:6)
            !write(6,'(''utmp6 ='', t20,  6(g15.7,1x) )') utmp(1:6)
            !write(6,'(''r utmp6 - u6 ='', t20,  6(g15.7,1x) )') utmp(1:6) - u(1:6) 
            !if ( sum((u(:)-utmp(:))**2)>=1.d-5   ) then
            !  write(6,'(''utmp6 ='', t20, 6(g15.7,1x) )') utmp(1:6)  
            !  write(6,*) '$$$$$$$$$  warning V6 part has issue in v6propr $$$$$$$$$$$$'
            !endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            
            
            
            
            
		   
		    if (iem) then
               
			    !call vemprop1(r,i,j,it-1,cwt1,cwtnew1) ! EM is on.   
			    !call v6prop1(x,i,j,u,cwtnew1,cwtnew2)
			    !call vemprop1(r,i,j,it-1,cwtnew2,cwtnew) ! EM is on. 
               
                ! since vem pp commute with v6, the following is enough.
                if (r.lt.rangev) then
                    vem(1)=c1*vemtab(index-1)+c2*vemtab(index)+c3*vemtab(index+1)+c4*vemtab(index+2)      
                else
                    call empot(2,r,vem) ! only pp simpliest coulomb here.
                endif 
                expvemcdt=exp(-vem(1)*2.0_r8**it*dt)      
			    call vemprop1(expvemcdt,i,j,cwt1,cwtnew1) ! EM is on.   
			    call v6prop1(x,i,j,u,cwtnew1,cwtnew)
               
		    else
			    call v6prop1(x,i,j,u,cwt1,cwtnew)  
     
                ! write(6,*) 'vpropr! i=,j=,cwtnew=', i, j, sum(cwtnew)
 	   
		    endif
		    cwt1(:,:)=cwtnew(:,:)
		   
		   
	    enddo
    enddo
    !cwtnew=cwt1
	
    return
    end subroutine v6propr
    
    subroutine v6proprpr(x,it,cwtold,cwtnew)
    integer(kind=i4) :: i,j,it,index,indexlag4	
    real(kind=r8) :: x(3,npart) ! the configuration for the system
    real(kind=r8), dimension(3) :: dx
    real(kind=r8) :: r,c1,c2,c3,c4,dr
    real(kind=r8) :: u(6)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		                ,cwtnew1(0:nspin-1,nisospin),cwt1(0:nspin-1,nisospin)
    real(kind=r8) :: expvemcdt,vem(14)	
    if (irep.eq.4) then
	    cwtnew(:,:)=cwtold(:,:) ! because for vmc, dt=0.	
	    return	
    endif
	
    cwt1(:,:)=cwtold(:,:)
	
    !write(6,*) 'x=',x
	
    do i=1,npart-1
	    do j=i+1,npart 		
            dx(:)=x(:,i)-x(:,j)		   
		    r=sqrt(dot_product(dx,dx)) !r=sqrt(sum(dx(:)**2))
		   
		    !if (r.eq.0) then
		    !    write(6,*) 'two particles have the same position, stop'
			   ! stop
		    !endif

            if (r.lt.rangev) then
			    dr=scalep*r
			    index=dr
			    index=max(1,min(index,ntab-2))
			    dr=dr-index

                c1=-dr*(dr-1)*(dr-2)/6
                c2=(dr+1)*(dr-1)*(dr-2)/2
                c3=-(dr+1)*dr*(dr-2)/2
                c4=(dr+1)*dr*(dr-1)/6

                !indexlag4=nint(dr/drlag4)
                !c1=laginterpc(1,indexlag4)
                !c2=laginterpc(2,indexlag4)
                !c3=laginterpc(3,indexlag4)
                !c4=laginterpc(4,indexlag4)              

                u(:)=c1*utabpr(:,index-1,it)+c2*utabpr(:,index,it) &
			        +c3*utabpr(:,index+1,it)+c4*utabpr(:,index+2,it)
			  
            !write(6,*) 'r,dr=',r,dr
            !write(6,*) 'dx=',dx	
            !write(6,*) 'i,j,u=',i,j,u
            !write(6,*) 'c1,2,3,4=',c1,c2,c3,c4
            !write(6,*) 'utab=',utab(:,index-1,it),utab(:,index,it) &
            !                  ,utab(:,index+1,it),utab(:,index+2,it)  
           
            else
                call v2ucalc(r,-1,it,u)
			    !u(:)=u0tab(:)
		    endif
		   
		   
		    if (iem) then
			    !call vemprop1(r,i,j,it-1,cwt1,cwtnew1) ! EM is on.   
			    !call v6prop1(x,i,j,u,cwtnew1,cwtnew2)
			    !call vemprop1(r,i,j,it-1,cwtnew2,cwtnew) ! EM is on.   
                if (r.lt.rangev) then
                    vem(1)=c1*vemtab(index-1)+c2*vemtab(index)+c3*vemtab(index+1)+c4*vemtab(index+2)      
                else
                    call empot(2,r,vem) ! only pp simpliest coulomb here.
                endif 
                expvemcdt=exp(-vem(1)*2.0_r8**it*dt)      
			    call vemprop1(expvemcdt,i,j,cwt1,cwtnew1) ! EM is on.   
			    call v6prop1(x,i,j,u,cwtnew1,cwtnew)
               
               
		    else
			    call v6prop1(x,i,j,u,cwt1,cwtnew)  
			   
		    endif
		    cwt1(:,:)=cwtnew(:,:)
		   
		   
	    enddo
    enddo
    !cwtnew=cwt1
	
	
    return
    end subroutine v6proprpr
    
    subroutine v6propl(x,it,cwtold,cwtnew)
    use mympi  
    integer(kind=i4) :: i,j,it,index,indexlag4	
    real(kind=r8) :: x(3,npart)
    real(kind=r8), dimension(3) :: dx
    real(kind=r8) :: r,c1,c2,c3,c4,dr
    real(kind=r8) :: u(6),utmp(6)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		                ,cwtnew1(0:nspin-1,nisospin),cwt1(0:nspin-1,nisospin)
    real(kind=r8) :: expvemcdt,vem(14)
    if (irep.eq.4) then
	    cwtnew(:,:)=cwtold(:,:) ! because for vmc, dt=0.	
	    return	
    endif	
	
    cwt1(:,:)=cwtold(:,:)
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

                !indexlag4=nint(dr/drlag4)
                !c1=laginterpc(1,indexlag4)
                !c2=laginterpc(2,indexlag4)
                !c3=laginterpc(3,indexlag4)
                !c4=laginterpc(4,indexlag4)                
          
                u(:)=c1*utab(:,index-1,it)+c2*utab(:,index,it) &
			        +c3*utab(:,index+1,it)+c4*utab(:,index+2,it)
			 		  			  
			    !write(6,*) 'r,dr=',r,dr  
			    !write(6,*) 'i,j,u=',i,j,u
        !         write(6,*) 'c1,2,3,4=',c1,c2,c3,c4
        !         write(6,*) 'dx=',dx		  
			  		  		  
            else
                call v2ucalc(r,1,it,u)
			    !u(:)=u0tab(:)
            endif
	
            !!!!!!!!!!!!!! check !!!!!!!!!!!!!!!!!!!         
            !call v2ucalc(r,1,it,utmp)
            !!write(6,'(''u6 ='', t20, 6(g15.7,1x) )') u(1:6)
            !!write(6,'(''utmp6 ='', t20,  6(g15.7,1x) )') utmp(1:6)
            !write(6,'(''l utmp6 - u6 ='', t20,  6(g15.7,1x) )') utmp(1:6) - u(1:6) 
            !if ( sum((u(:)-utmp(:))**2)>=1.d-5   ) then
            !  write(6,'(''utmp6 ='', t20, 6(g15.7,1x) )') utmp(1:6)  
            !  write(6,*) '$$$$$$$$$  warning V6 part has issue in v6propl $$$$$$$$$$$$'
            !endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
		   
		    if (iem) then
			    !call vemprop1(r,i,j,it-1,cwt1,cwtnew1) ! EM is on.   
			    !call v6prop1(x,i,j,u,cwtnew1,cwtnew2)
			    !call vemprop1(r,i,j,it-1,cwtnew2,cwtnew) ! EM is on.   
                if (r.lt.rangev) then
                    vem(1)=c1*vemtab(index-1)+c2*vemtab(index)+c3*vemtab(index+1)+c4*vemtab(index+2)      
                else
                    call empot(2,r,vem) ! only pp simpliest coulomb here.
                endif 
                expvemcdt=exp(-vem(1)*2.0_r8**it*dt)      
			    call vemprop1(expvemcdt,i,j,cwt1,cwtnew1) ! EM is on.   
			    call v6prop1(x,i,j,u,cwtnew1,cwtnew)
		    else
			    call v6prop1(x,i,j,u,cwt1,cwtnew)  ! need to use schmidt order. 
			   
		    endif
		    cwt1(:,:)=cwtnew(:,:)
		   
		   
	    enddo
    enddo

    return
    end subroutine v6propl
    
    subroutine v6proplpr(x,it,cwtold,cwtnew)
    integer(kind=i4) :: i,j,it,index,indexlag4	
    real(kind=r8) :: x(3,npart)
    real(kind=r8), dimension(3) :: dx
    real(kind=r8) :: r,c1,c2,c3,c4,dr
    real(kind=r8) :: u(6)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		                ,cwtnew1(0:nspin-1,nisospin),cwt1(0:nspin-1,nisospin)
    real(kind=r8) :: expvemcdt,vem(14)	
    if (irep.eq.4) then
	    cwtnew(:,:)=cwtold(:,:) ! because for vmc, dt=0.	
	    return	
    endif
	
    cwt1(:,:)=cwtold(:,:)
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
                
                !indexlag4=nint(dr/drlag4)
                !c1=laginterpc(1,indexlag4)
                !c2=laginterpc(2,indexlag4)
                !c3=laginterpc(3,indexlag4)
                !c4=laginterpc(4,indexlag4)                
                
                u(:)=c1*utabpr(:,index-1,it)+c2*utabpr(:,index,it) &
			        +c3*utabpr(:,index+1,it)+c4*utabpr(:,index+2,it)
			 		  			  
			    !write(6,*) 'r,dr=',r,dr  
			    !write(6,*) 'i,j,u=',i,j,u
        !         write(6,*) 'c1,2,3,4=',c1,c2,c3,c4
        !         write(6,*) 'dx=',dx		  
			  		  		  
            else
                call v2ucalc(r,-1,it,u)
			    !u(:)=u0tab(:)
		    endif
		   
		   
		    if (iem) then
			    !call vemprop1(r,i,j,it-1,cwt1,cwtnew1) ! EM is on.   
			    !call v6prop1(x,i,j,u,cwtnew1,cwtnew2)
			    !call vemprop1(r,i,j,it-1,cwtnew2,cwtnew) ! EM is on.   
                if (r.lt.rangev) then
                    vem(1)=c1*vemtab(index-1)+c2*vemtab(index)+c3*vemtab(index+1)+c4*vemtab(index+2)      
                else
                    call empot(2,r,vem) ! only pp simpliest coulomb here.
                endif 
                expvemcdt=exp(-vem(1)*2.0_r8**it*dt)                
			    call vemprop1(expvemcdt,i,j,cwt1,cwtnew1) ! EM is on.   
			    call v6prop1(x,i,j,u,cwtnew1,cwtnew)
		    else
			    call v6prop1(x,i,j,u,cwt1,cwtnew)  
			   
		    endif
		    cwt1(:,:)=cwtnew(:,:)
		   	   
	    enddo
    enddo

    return
    end subroutine v6proplpr  
	
    subroutine v6proprl(xl,xr,it,cwtold,cwtnew) ! timestep=2**it dt
    integer(kind=i4) :: it
    real(kind=r8) :: xl(3,npart),xr(3,npart) 
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		                ,cwtnew1(0:nspin-1,nisospin)
    if (irep.eq.4) then
	    cwtnew(:,:)=cwtold(:,:) ! because for vmc, dt=0.	
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
		                ,cwtnew1(0:nspin-1,nisospin)
    if (irep.eq.4) then
	    cwtnew(:,:)=cwtold(:,:) ! because for vmc, dt=0.	
	    return	
    endif
    call v6propl(xl,it,cwtold,cwtnew1)
    call v6propr(xr,it,cwtnew1,cwtnew)  
    return
    end subroutine v6proplr    	
       
    subroutine v6val(x,cwtold,cwtnew)
    use v6pot
    integer(kind=i4) :: i,j,index,indexlag4	
    real(kind=r8) :: x(3,npart) ! the configuration for the system
    real(kind=r8), dimension(3) :: dx
    real(kind=r8) :: r,c1,c2,c3,c4,dr
    real(kind=r8) :: v(6),vv(18)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin),cwt1(0:nspin-1,nisospin)	
    cwtnew(:,:)=0
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
   
                !indexlag4=nint(dr/drlag4)
                !c1=laginterpc(1,indexlag4)
                !c2=laginterpc(2,indexlag4)
                !c3=laginterpc(3,indexlag4)
                !c4=laginterpc(4,indexlag4)                

                v(:)=c1*vtab(:,index-1)+c2*vtab(:,index)+c3*vtab(:,index+1)+c4*vtab(:,index+2)
            else
                call av18op(3,r,vv)
			    call v6schmidtord(vv,v)
                !v(:)=0
		    endif
		    call v6prop1(x,i,j,v,cwtold,cwt1)
		    cwtnew=cwt1+cwtnew				
	    enddo
    enddo
    return
    end subroutine v6val	
	
    subroutine v6val1(x,n,cwtnew) ! can be optimized to calculate for a particular state n.
    use v6pot
    integer(kind=i4) :: i,j,n,index,indexlag4	
    real(kind=r8) :: x(3,npart) ! the configuration for the system
    real(kind=r8), dimension(3) :: dx
    real(kind=r8) :: r,c1,c2,c3,c4,dr
    real(kind=r8) :: v(6),vv(18)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin),cwt1(0:nspin-1,nisospin)	
    cwtold=0
    cwtold(invspin(n),invispin(n))=1
    cwtnew(:,:)=0
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

                !indexlag4=nint(dr/drlag4)
                !c1=laginterpc(1,indexlag4)
                !c2=laginterpc(2,indexlag4)
                !c3=laginterpc(3,indexlag4)
                !c4=laginterpc(4,indexlag4)              

                v(:)=c1*vtab(:,index-1)+c2*vtab(:,index)+c3*vtab(:,index+1)+c4*vtab(:,index+2)
            else
                call av18op(3,r,vv)
                call v6schmidtord(vv,v)                  
			    !v(:)=0
		    endif
		    call v6prop1(x,i,j,v,cwtold,cwt1)
		    cwtnew=cwt1+cwtnew				
	    enddo
    enddo
    return
    end subroutine v6val1	
    

    !subroutine v6prop1b(x,i,j,u,cwtold,cwtnew) ! optimized version 20190703, 25 times faster than old verison.
    !integer(kind=i4) :: i,j,nstate,nstatenew,nstatenew1,nstatenew2
    !real(kind=r8) :: x(3,npart),r(3),dx(3)
    !real(kind=r8) :: u(6) ! v6 so 6, coefficients after conversion. check if real or complex, I think real.
    !complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
    !                   ,fceff,fsex,fiex,fsiex,fteff,fteffiex,f(6)
    !cwtnew(:,:)=0 
    !dx(:)=x(:,i)-x(:,j)
    !r(:)=dx(:)/sqrt(dot_product(dx,dx)) !unit vector
    !do nstate=1,nbasis   
    !    f(:)=u(:)*cwtold(invspin(nstate),invispin(nstate))
    !    fceff=f(1)-f(2)-f(4)+f(5)+f(3)-f(6)
    !    fsex=2*(f(2)-f(5)-f(3)+f(6))
    !    fiex=2*(f(4)-f(5)+f(6))
    !    fsiex=4*(f(5)-f(6))
    !    fteff=3*(f(3)-f(6))
    !    fteffiex=6*f(6)
    !    ! fceff,op1
    !    cwtnew(invspin(nstate),invispin(nstate))=cwtnew(invspin(nstate),invispin(nstate))+fceff
    !    ! fsex,opexsigma
    !    call opexsigma(i,j,nstate,nstatenew1)
    !    cwtnew(invspin(nstatenew1),invispin(nstatenew1))=cwtnew(invspin(nstatenew1),invispin(nstatenew1))+fsex
    !    ! fiex,opextau        
    !    call opextau(i,j,nstate,nstatenew2)
    !    cwtnew(invspin(nstatenew2),invispin(nstatenew2))=cwtnew(invspin(nstatenew2),invispin(nstatenew2))+fiex        
    !    ! fsiex
    !    call opextau(i,j,nstatenew1,nstatenew) 
    !    !call opexsigma(i,j,nstatenew1,nstatenew)   
    !    !call opexsigmatau(i,j,nstate,nstatenew)
    !    cwtnew(invspin(nstatenew),invispin(nstatenew))=cwtnew(invspin(nstatenew),invispin(nstatenew))+fsiex      
    !    ! fteff,op3
    !    call op3pure(r,i,j,nstate,fteff,cwtnew)
    !    ! fteffiex,opextau,op3
    !    !call opextau(i,j,nstate,nstatenew) 
    !    call op3pure(r,i,j,nstatenew2,fteffiex,cwtnew)
    !enddo 
    !return
    !end subroutine v6prop1b     
    
    !subroutine v6prop1a(x,i,j,u,cwtold,cwtnew) ! optimized version 20190703, 15-20 times faster than old verison.
    !integer(kind=i4) :: i,j,nstate,nstatenew,nstatenew1
    !real(kind=r8) :: x(3,npart),r(3),dx(3)
    !real(kind=r8) :: u(6) ! v6 so 6, coefficients after conversion. check if real or complex, I think real.
    !complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
    !                   ,fceff,fsex,fiex,fsiex,fteff,fteffiex,f(6)
    !
    !cwtnew(:,:)=0
    !
    !dx(:)=x(:,i)-x(:,j)
    !r(:)=dx(:)/sqrt(dot_product(dx,dx)) !unit vector
    !
    !do nstate=1,nbasis
    !      
    !    f(:)=u(:)*cwtold(invspin(nstate),invispin(nstate))
    !
    !    fceff=f(1)-f(2)-f(4)+f(5)+f(3)-f(6)
    !    fsex=2*(f(2)-f(5)-f(3)+f(6))
    !    fiex=2*(f(4)-f(5)+f(6))
    !    fsiex=4*(f(5)-f(6))
    !    fteff=3*(f(3)-f(6))
    !    fteffiex=6*f(6)
    ! 
    !    ! fceff,op1
    !    cwtnew(invspin(nstate),invispin(nstate))=cwtnew(invspin(nstate),invispin(nstate)) &
    !                                            +fceff
    !    ! fsex,opexsigma
    !    call opexsigma(i,j,nstate,nstatenew)
    !    cwtnew(invspin(nstatenew),invispin(nstatenew))=cwtnew(invspin(nstatenew),invispin(nstatenew)) &
    !                                                +fsex
    !    ! fiex,opextau        
    !    call opextau(i,j,nstate,nstatenew)
    !    cwtnew(invspin(nstatenew),invispin(nstatenew))=cwtnew(invspin(nstatenew),invispin(nstatenew)) &
    !                                                +fiex        
    !    ! fisex
    !    !call opextau(i,j,nstate,nstatenew1) 
    !    !call opexsigma(i,j,nstatenew1,nstatenew)
    !    
    !    call opexsigmatau(i,j,nstate,nstatenew)
    !    cwtnew(invspin(nstatenew),invispin(nstatenew))=cwtnew(invspin(nstatenew),invispin(nstatenew)) &
    !                                                +fsiex      
    !    ! fteff,op3
    !    call op3a(r,i,j,nstate,fteff,cwtnew)
    !    ! fteffiex,opextau,op3
    !    call opextau(i,j,nstate,nstatenew)
    !    call op3a(r,i,j,nstatenew,fteffiex,cwtnew)
    !
    !enddo
    !
    !return
    !end subroutine v6prop1a  	
    !
    !subroutine op1(nstate,cwtnew) ! op1=1
    !integer(kind=i4) :: nstate
    !complex(kind=r8) :: cwtnew(0:nspin-1,nisospin)
    !cwtnew(:,:)=0
    !cwtnew(invspin(nstate),invispin(nstate))=1
    !return
    !end subroutine op1
		  !
    !subroutine op2(i,j,nstate,cwtnew) ! op2 = sigma i dot sigma j
    !integer(kind=i4) :: i,j,nstate,spn,newspn,signi,signj
    !complex(kind=r8) :: cwtnew(0:nspin-1,nisospin)
    !cwtnew(:,:)=0
    !spn=invspin(nstate) ! spn, 0:15.
    !signi=2*(and(1,shiftr(spn,i-1)))-1
    !signj=2*(and(1,shiftr(spn,j-1)))-1
    !
    !newspn=xor(xor(spn,shiftl(1,i-1)),shiftl(1,j-1)) ! flip i, j spin.
    !cwtnew(newspn,invispin(nstate))=1 & ! six sjx	
	   !                             -signi*signj & ! siy sjy, -1 bc i*i
	   !                             +cwtnew(newspn,invispin(nstate))  
    !! siz sjz
    !cwtnew(spn,invispin(nstate))=signi*signj+cwtnew(spn,invispin(nstate))	
    !return
    !end subroutine op2
    !! can check op2 with 2opexsigma-1. that should be slightly faster i think.
    !
    !subroutine op2ex(i,j,nstate,cwtnew) ! op2 = 2*opexsigma-1. should be slightly faster  subroutine op2.
    !integer(kind=i4) :: i,j,nstate,nstatenew
    !complex(kind=r8) :: cwtnew(0:nspin-1,nisospin)
    !cwtnew(:,:)=0
    !cwtnew(invspin(nstate),invispin(nstate))=-1
    !call opexsigma(i,j,nstate,nstatenew)
    !cwtnew(invspin(nstatenew),invispin(nstatenew))=cwtnew(invspin(nstatenew),invispin(nstatenew))+2
    !return
    !end subroutine op2ex    
    !
    !subroutine op3(x,i,j,nstate,cwtnew) ! op3 = tensor op, make sure dx neq 0, otherwise need to add protection.
    !integer(kind=i4) :: i,j,nstate,spn,newspn,newspni,newspnj,signi,signj
    !real(kind=r8) :: r(3),x(3,npart),dx(3)
    !complex(kind=r8) :: cwtnew(0:nspin-1,nisospin)
    !cwtnew(:,:)=0
    !spn=invspin(nstate) ! spn, 0:15.
    !signi=2*(and(1,shiftr(spn,i-1)))-1
    !signj=2*(and(1,shiftr(spn,j-1)))-1
    !newspni=xor(spn,shiftl(1,i-1))
    !newspnj=xor(spn,shiftl(1,j-1))
    !newspn=xor(newspni,shiftl(1,j-1))
    !dx(:)=x(:,i)-x(:,j)
    !r(:)=dx(:)/sqrt(dot_product(dx,dx)) !unit vector
    !!write(6,*) r
    !cwtnew(newspn,invispin(nstate))=(r(1)**2 &	! six sjx
	   !                             -(signi*signj)*r(2)**2 & ! siy sjy
	   !                             +(signi+signj)*ci*r(1)*r(2))*3 & ! six sjy,  siy sjx
	   !                             +cwtnew(newspn,invispin(nstate))
    !cwtnew(newspni,invispin(nstate))=(signj*r(1)*r(3) & ! six sjz
				!				    +signi*ci*signj*r(2)*r(3))*3 & ! siy sjz
				!				    +cwtnew(newspni,invispin(nstate))
    !cwtnew(newspnj,invispin(nstate))=(signi*r(1)*r(3) & ! six sjz
				!				    +signj*ci*signi*r(2)*r(3))*3 & ! siy sjz	
	   !                             +cwtnew(newspnj,invispin(nstate))
    !cwtnew(spn,invispin(nstate))=signi*signj*r(3)**2*3+cwtnew(spn,invispin(nstate)) ! siz sjz
    !! -op2	
    !cwtnew(newspn,invispin(nstate))=cwtnew(newspn,invispin(nstate)) &
				!				    -1 & ! six sjx
	   !                             +signi*signj  ! siy sjy 
    !cwtnew(spn,invispin(nstate))=cwtnew(spn,invispin(nstate)) &
	   !                         -signi*signj ! siz sjz
    !return
    !end subroutine op3	
    !
    !subroutine op4(i,j,nstate,cwtnew) ! tau i dot tau j
    !integer(kind=i4) :: i,j,nstate,nstatenew
    !complex(kind=r8) :: cwtnew(0:nspin-1,nisospin)
    !cwtnew(:,:)=0	
    !cwtnew(invspin(nstate),invispin(nstate))=-1
    !call opextau(i,j,nstate,nstatenew)	
    !cwtnew(invspin(nstatenew),invispin(nstatenew))=cwtnew(invspin(nstatenew),invispin(nstatenew))+2
    !return
    !end subroutine op4    
    !
    !subroutine op5(i,j,nstate,cwtnew) ! sigma i dot sigma j tau i dot tau j
    !integer(kind=i4) :: i,j,nstate,nstatenew1,nstatenew2,nstatenew3
    !complex(kind=r8) :: cwtnew(0:nspin-1,nisospin)
    !cwtnew(:,:)=0	
    !call opextau(i,j,nstate,nstatenew1)	
    !call opexsigma(i,j,nstatenew1,nstatenew2)
    !call opexsigma(i,j,nstate,nstatenew3)
    !cwtnew(invspin(nstatenew2),invispin(nstatenew2))=4+cwtnew(invspin(nstatenew2),invispin(nstatenew2))
    !cwtnew(invspin(nstatenew3),invispin(nstatenew3))=-2+cwtnew(invspin(nstatenew3),invispin(nstatenew3))
    !cwtnew(invspin(nstatenew1),invispin(nstatenew1))=-2+cwtnew(invspin(nstatenew1),invispin(nstatenew1))
    !cwtnew(invspin(nstate),invispin(nstate))=1+cwtnew(invspin(nstate),invispin(nstate))
    !return
    !end subroutine op5	  
	   ! 
    !subroutine op6(x,i,j,nstate,cwtnew) ! op6 = tensor  tau i dot tau j
    !integer(kind=i4) :: i,j,nstate,nstatenew1
    !real(kind=r8) :: x(3,npart)
    !complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin)
    !cwtnew(:,:)=0
    !call opextau(i,j,nstate,nstatenew1)	
    !call op3(x,i,j,nstatenew1,cwtnew1)
    !call op3(x,i,j,nstate,cwtnew2)
    !cwtnew=2*cwtnew1-cwtnew2
    !return
    !end subroutine op6	    
    !
  !  subroutine op1sumbasis(cwtold,cwtnew)
  !  complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin)
  !  cwtnew(:,:)=cwtold(:,:)
  !  return
  !  end subroutine op1sumbasis
  !
  !  subroutine op2sumbasis(i,j,cwtold,cwtnew)
  !  integer(kind=i4) :: i,j,is	
  !  complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		!                ,cwtnew1(0:nspin-1,nisospin)
  !  cwtnew(:,:)=0
  !  do is=1,nbasis
		!    call op2ex(i,j,is,cwtnew1) ! call op2(i,j,is,cwtnew1) 
		!    cwtnew=cwtnew+cwtnew1*cwtold(invspin(is),invispin(is))
  !  enddo 
  !  return
  !  end subroutine op2sumbasis
	 !
  !  subroutine op3sumbasis(x,i,j,cwtold,cwtnew)
  !  integer(kind=i4) :: i,j,is
  !  real(kind=r8) :: x(3,npart)
  !  complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		!                ,cwtnew1(0:nspin-1,nisospin)
  !  cwtnew(:,:)=0
  !  do is=1,nbasis
	 !   call op3(x,i,j,is,cwtnew1) 
	 !   cwtnew=cwtnew+cwtnew1*cwtold(invspin(is),invispin(is))
	 !   !cwtnews(:,:,is)=cwtnew1(:,:)*cwtold(invspin(is),invispin(is))
	 !   !write(11,*) is,cwtnew(invspin(is),invispin(is)),cwtnew1(invspin(is),invispin(is)),cwtold(invspin(is),invispin(is))
  !  enddo 
  !  !do is=1,nbasis		
  !  !   cwtnew(:,:)=cwtnew(:,:)+cwtnews(:,:,is)
  !  !   write(11,*) is 
  !  !   k=11
  !  !   write(11,*) k,cwtnews(invspin(k),invispin(k),is),cwtnew(invspin(k),invispin(k))
  !  !enddo
  !  return
  !  end subroutine op3sumbasis	
	 !
  !  subroutine op4sumbasis(i,j,cwtold,cwtnew)
  !  integer(kind=i4) :: i,j,is	
  !  complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		!                ,cwtnew1(0:nspin-1,nisospin)
  !  cwtnew(:,:)=0
  !  do is=1,nbasis
		!call op4(i,j,is,cwtnew1) 
		!cwtnew=cwtnew+cwtnew1*cwtold(invspin(is),invispin(is))
  !  enddo 
  !  return
  !  end subroutine op4sumbasis
	 !
  !  subroutine op5sumbasis(i,j,cwtold,cwtnew)
  !  integer(kind=i4) :: i,j,is	
  !  complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		!                ,cwtnew1(0:nspin-1,nisospin)
  !  cwtnew(:,:)=0
  !  do is=1,nbasis
		!call op5(i,j,is,cwtnew1) 
		!cwtnew=cwtnew+cwtnew1*cwtold(invspin(is),invispin(is))
  !  enddo 
  !  return
  !  end subroutine op5sumbasis
	 !
  !  subroutine op6sumbasis(x,i,j,cwtold,cwtnew)
  !  integer(kind=i4) :: i,j,is	
  !  real(kind=r8) :: x(3,npart)
  !  complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		!                ,cwtnew1(0:nspin-1,nisospin)
  !  cwtnew(:,:)=0
  !  do is=1,nbasis
		!call op6(x,i,j,is,cwtnew1) 
		!cwtnew=cwtnew+cwtnew1*cwtold(invspin(is),invispin(is))
  !  enddo
  !  return
  !  end subroutine op6sumbasis    

    

    !subroutine opexsigma1sumbasis(i,j,cwtold,cwtnew)
    !integer(kind=i4) :: i,j,is	
    !complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		  !              ,cwtnew1(0:nspin-1,nisospin)
    !cwtnew(:,:)=0
    !do is=1,nbasis
		  !  call opexsigma1(i,j,is,cwtnew1) 
		  !  cwtnew=cwtnew+cwtnew1*cwtold(invspin(is),invispin(is))
    !enddo 
    !return
    !end subroutine opexsigma1sumbasis
	   !
    !subroutine opexsigmasumbasis(i,j,cwtold,cwtnew)
    !integer(kind=i4) :: i,j,is,in
    !complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) 
    !cwtnew(:,:)=0
    !do is=1,nbasis
		  !  call opexsigma(i,j,is,in) 
		  !  cwtnew2=0
		  !  cwtnew2(invspin(in),invispin(in))=1
		  !  cwtnew=cwtnew+cwtnew2*cwtold(invspin(is),invispin(is))
    !enddo 
    !return
    !end subroutine opexsigmasumbasis
	   !
    !subroutine opexsigma0(i,j,nstate,nstatenew) ! spin exchange with if then else 
    !integer(kind=i4) :: i,j,iv,jv,nstate,nstatenew,spn,newspn
    !!cwtnew(:,:)=0
    !spn=invspin(nstate) ! spn, 0:15.
    !iv=shiftr(and(spn,shiftl(1,i-1)),i-1)
    !jv=shiftr(and(spn,shiftl(1,j-1)),j-1)	
    !if (iv.eq.jv) then
    !!cwtnew(spn,invispin(nstate))=1.	
    !nstatenew=nstate
    !else
    !newspn=xor(spn,shiftl(1,i-1))
    !newspn=xor(newspn,shiftl(1,j-1))	
    !!cwtnew(newspn,invispin(nstate))=1.	
    !nstatenew=invcwt(newspn,invispin(nstate))
    !endif
    !return
    !end subroutine opexsigma0
	   !
    !subroutine opexsigma1(i,j,nstate,cwtnew) ! spin exchange, use op2, haven't check. may be slow.
    !integer(kind=i4) :: i,j,nstate,spn,newspn,signi,signj
    !complex(kind=r8) :: cwtnew(0:nspin-1,nisospin)
    !!call op2(i,j,nstate,1.,cwtnew)
    !cwtnew(:,:)=0
    !spn=invspin(nstate) ! spn, 0:15.
    !signi=2*(and(1,shiftr(spn,i-1)))-1
    !signj=2*(and(1,shiftr(spn,j-1)))-1
    !newspn=xor(spn,shiftl(1,i-1))
    !newspn=xor(newspn,shiftl(1,j-1))	
    !cwtnew(newspn,invispin(nstate))=cwtnew(newspn,invispin(nstate))+1 & ! six sjx	
	   !                             -signi*signj ! siy sjy
    !! siz sjz
    !cwtnew(spn,invispin(nstate))=cwtnew(spn,invispin(nstate))+signi*signj
    !!			
    !cwtnew(spn,invispin(nstate))=cwtnew(spn,invispin(nstate))+1
    !cwtnew=cwtnew/2
    !return
    !end subroutine opexsigma1
	   !
    !
	   !
    !subroutine opextau0(i,j,nstate,nstatenew) ! tau exchange with if then else
    !integer(kind=i4) :: i,j,iv,jv,nstate,nstatenew,ispn,newispn
    !!cwtnew(:,:)=0
    !ispn=liso(invispin(nstate)) 
    !iv=shiftr(and(ispn,shiftl(1,i-1)),i-1)
    !jv=shiftr(and(ispn,shiftl(1,j-1)),j-1)	
    !if (iv.eq.jv) then
    !!cwtnew(invspin(nstate),lisoi(ispn))=1.	
    !nstatenew=nstate
    !else
    !newispn=xor(ispn,shiftl(1,i-1))
    !newispn=xor(newispn,shiftl(1,j-1))	
    !!cwtnew(invspin(nstate),lisoi(newispn))=1.
    !nstatenew=invcwt(invspin(nstate),lisoi(newispn))
    !endif
    !return
    !end subroutine opextau0	
    	
    
    
    !subroutine v6prop1old(x,i,j,u,cwtold,cwtnew) ! old and slow   nama: v6prop1
    !integer(kind=i4) :: i,j	
    !real(kind=r8) :: x(3,npart)
    !real(kind=r8) :: u(6) ! v6 so 6, coefficients after conversion. check if real or complex, I think real.
    !complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		  !              ,cwtnew1(0:nspin-1,nisospin)
    !cwtnew(:,:)=0
    !call op1sumbasis(cwtold,cwtnew1)
    !cwtnew=cwtnew+cwtnew1*u(1)
    !call op2sumbasis(i,j,cwtold,cwtnew1)
    !cwtnew=cwtnew+cwtnew1*u(2)
    !call op3sumbasis(x,i,j,cwtold,cwtnew1)
    !cwtnew=cwtnew+cwtnew1*u(3)
    !call op4sumbasis(i,j,cwtold,cwtnew1)
    !cwtnew=cwtnew+cwtnew1*u(4)
    !call op5sumbasis(i,j,cwtold,cwtnew1)
    !cwtnew=cwtnew+cwtnew1*u(5)
    !call op6sumbasis(x,i,j,cwtold,cwtnew1)
    !cwtnew=cwtnew+cwtnew1*u(6)
    !return
    !end subroutine v6prop1old
    !
    !
    !subroutine v6prop1schmidt(x,i,j,u,cwtold,cwtnew)
    !integer(kind=i4) :: i,j	
    !real(kind=r8) :: x(3,npart)
    !real(kind=r8) :: u(6) ! v6 so 6, coefficients after conversion. check if real or complex, I think real.
    !complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		  !              ,cwtnew1(0:nspin-1,nisospin)
    !cwtnew(:,:)=0
    !call op1sumbasis(cwtold,cwtnew1)
    !cwtnew=cwtnew+cwtnew1*u(1)
    !call op2sumbasis(i,j,cwtold,cwtnew1)
    !cwtnew=cwtnew+cwtnew1*u(2)
    !call op3sumbasis(x,i,j,cwtold,cwtnew1)
    !cwtnew=cwtnew+cwtnew1*u(5)
    !call op4sumbasis(i,j,cwtold,cwtnew1)
    !cwtnew=cwtnew+cwtnew1*u(3)
    !call op5sumbasis(i,j,cwtold,cwtnew1)
    !cwtnew=cwtnew+cwtnew1*u(4)
    !call op6sumbasis(x,i,j,cwtold,cwtnew1)
    !cwtnew=cwtnew+cwtnew1*u(6)
    !return
    !end subroutine v6prop1schmidt	
	
    !subroutine psitprop(x,cwtold,cwtnew) ! Schmidt order.
    !use brut
    !integer(kind=i4) :: i,j,k	
    !real(kind=r8) :: x(3,npart) ! the configuration for the system
    !real(kind=r8), dimension(3) :: dx
    !real(kind=r8) :: r
    !real(kind=r8) :: fop(6)
    !complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin),cwt1(0:nspin-1,nisospin)
    !integer :: ipl(6),jpl(6),ipr(6),jpr(6)
	   !
    !cwt1(:,:)=cwtold(:,:)
	   !
    !call outputordlr(ipl,jpl,ipr,jpr)
	   !
    !do k=1,6 ! 6 is npair
	   ! i=ipl(k)
	   ! j=jpl(k)
	   ! !write(6,*) 'i=', i
	   ! !write(6,*) 'j=', j
	   ! dx(:)=x(:,i)-x(:,j)		   
	   ! r=sqrt(dot_product(dx,dx))
	   ! call calf(r,fop)
		  !
	   ! !fop=0.
	   ! !fop(6)=1.
		  !
	   ! write(6,*) 'r=',r
	   ! write(6,*) 'fop=',fop
	   ! call v6prop1schmidt(x,i,j,fop,cwt1,cwtnew) 
	   ! cwt1(:,:)=cwtnew(:,:)	
    !enddo
	   !
    !return
    !end subroutine psitprop	
	
    !subroutine v6mat(x,vij,vijeigvec,vijeigval)
    !use matrixmod
    !use brut
    !real(kind=r8) :: x(3,npart),x0(3,npart)
    !complex(kind=r8) :: cwt(0:nspin-1,nisospin),cwtnew(0:nspin-1,nisospin),cwtgnd(0:nspin-1,nisospin)
    !integer(kind=i4) ::  i,j,k,is,iiso 
    !complex(kind=r8) :: vij(nbasis,nbasis),vijeigvec(nbasis,nbasis)
    !real(kind=r8) :: vijeigval(nbasis),t 
    !! <i|V6|j> matrix first.  
    !!    invcwt 
    !! maybe do a x convert. 
    !! maybe make it 16 * 16 eigenstates. see if it can be faster.
    !x0=x
    !call xconvert(x0)
    !call getcwtgnd2(cwtgnd) ! no need perhaps.
    !do j=1,nbasis   
    !    call v6val1(x0,j,cwtnew) 
    !    do i=1,nbasis  
    !        vij(i,j)=cwtnew(invspin(i),invispin(i)) ! vij should be hermitian. can check.
    !    enddo
    !enddo
    !! can check hermitian.
    !
    !! diagnalize the hermitian vij
    !    
    !call eigen(vij,vijeigvec,vijeigval,nbasis)
    !
    !return
    !end subroutine v6mat 	
	
	
	
	
	

	
	
	
	
	
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
