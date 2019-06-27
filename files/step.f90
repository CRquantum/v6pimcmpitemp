module step
    use chorizos
    implicit none
    integer, private, parameter :: i4=selected_int_kind(9)
    integer, private, parameter :: i8=selected_int_kind(15)
    integer, private, parameter :: r8=selected_real_kind(15,9)
    type (chorizo), private, allocatable, save :: ristra(:),ristraold(:),ristranew(:),ristranew1(:)
    type (chorizo), private, allocatable, save :: newristra(:),oldristra(:)
    integer(kind=i4), private, save :: npart,nprot,nspin,nisospin,nbasis
    integer(kind=i4), private, save :: nav,neq,nstep,nstepnow,nstepstuck,nblocknow,nblockstuck
    integer(kind=i4), private, save :: nchorizo,nchorizomid,nchorizomod,nrepmax,irep,nstepdecor,nrhobin
    integer(kind=i4), private, save :: mmax,nbisect,nbisecthalf,mmaxhalf,nreparrow,icrephistlast,icrephist
    real(kind=r8), private, save ::  hbar,dt,sigma,driftm,el,repprob
    real(kind=r8), private, save ::  etrial,ecut,eloctot,rhobinsize
    real(kind=r8), private, save ::  x0step,mov1step,mov2step,shftstep
    integer(kind=i4), private, save :: ic0,ic0tot,icn,icntot,ibisect,ibisecttot,ieloc &
                                        ,icmov1,icmov1tot,icstdmove1tot,icshft,icshfttot,ic23,ic23tot,ice,icbisecttot &
	                                    ,ibisectimptot,icstepstd &
                                        ,ibisectl,ibisecttotl,ibisectr,ibisecttotr &
                                    ,iclr,iclrtot,icrep,icreptot &
                                    ,ibisecthalf,ibisecthalftot,ibisectdouble,ibisectdoubletot &
                                    ,icstdmove2tot,icstdmove2 &
                                    ,irepstepct,irepstepmonmax 
    integer(kind=i4), private, save, allocatable :: repstepct1(:),repstepct2(:),repstepct1sum(:),repstepct2sum(:) &
                                                   ,ireparrowhistory(:) &
                                                   ,ibisectstuck(:)
    real(kind=r8), private, save, allocatable ::  bisectrate(:),bisectcount(:),bisecttot(:) &
                                                ,bisectcountl(:),bisecttotl(:) &
                                                ,bisectcountr(:),bisecttotr(:) &
                                                ,bisecthalfcount(:),bisecthalftot(:)
    complex(kind=r8), private, save, allocatable :: cwtbegin(:,:),cwtend(:,:),cwtl2r(:,:,:),cwtr2l(:,:,:) &
                                                    ,cwtbegin1(:,:),cwtend1(:,:),cwtl2r1(:,:,:),cwtr2l1(:,:,:) &
                                                    ,cwtbeginnew(:,:),cwtendnew(:,:),cwtl2rnew(:,:,:),cwtr2lnew(:,:,:) &
                                                    ,cwtbeginnew1(:,:),cwtendnew1(:,:),cwtl2rnew1(:,:,:),cwtr2lnew1(:,:,:)  &
                                                    ,cwtground(:,:)
    real(kind=r8), private, save, allocatable :: xin(:),yin(:),corr(:)
    integer(kind=i4), private, allocatable, save :: invspin(:),invispin(:)
    logical, private, save :: icorrchk,isite
    character(len=70), private, save :: infile,outfile
    integer(kind=i4), private, save :: iplo(6),jplo(6),ipro(6),jpro(6) &
                            ,ipl(6),jpl(6),ipr(6),jpr(6)	
    logical, private, save :: rejectbisectv6improve,rejectbisectv6improvehalf
    character(len=120), private, save, allocatable :: answer(:)
    integer(kind=i4), private, save :: ipath,ipathtmp,npo,ixtempo    

contains
    subroutine stepinit(npartin,nchorizoin,hbarin,dtin &
        ,mmaxin,nrepmaxin,irepin,nstepdecorin &
	    ,x0stepin,mov1stepin,mov2stepin,shftstepin,nrhobinin,icorrchkin &
	    ,nprotin,navin,nstepin,neqin)
    use math
    use v6stepcalc
    real(kind=r8) :: hbarin,dtin &
                    ,x0stepin,mov1stepin,mov2stepin,shftstepin
    integer(kind=i4) :: npartin,nchorizoin,nrepmaxin,irepin,mmaxin,nstepdecorin,nrhobinin,nprotin &
                        ,navin,nstepin,neqin
    integer(kind=i4) :: i,nristra
    logical :: icorrchkin
    !repprob=repprobin
    npart=npartin    
    x0step=x0stepin
    !nrepmax=nrepmaxin
    mmax=mmaxin
    nrepmax=nrepmaxin
    irep=irepin
    nreparrow=1 ! initialize reptation direction. 1 is right, -1 is to the left.  
    irepstepct=0
    ipath=1 ! initialize the number of path as 1. count the path number for equilibriated path.
    ipathtmp=1 ! count the path number for unequilibriated path.
    hbar=hbarin
    dt=dtin
    sigma=sqrt(2.0_r8*hbar*dt) ! free particle propagator's sigma.
    !driftm=sigma*dcut
    nstepdecor=nstepdecorin
    mov1step=mov1stepin
    mov2step=mov2stepin
    shftstep=shftstepin
    nrhobin=nrhobinin
    icorrchk=icorrchkin
    call chorizoinit(npart)
    nchorizo=nchorizoin
    nchorizomod=nchorizo+1
    nchorizomid=nint(dble(nchorizo)/2.)
    nbisect=2**mmax
    nbisecthalf=nbisect/2
    mmaxhalf=mmax-1
    nristra=nbisect ! max(nrepmax,nbisect)
    allocate(ristra(0:nchorizo),ristraold(0:nchorizo),ristranew(0:nchorizo),ristranew1(0:nchorizo))
   
    do i=0,nchorizo
        allocate(ristra(i)%x(3,npart),ristra(i)%dpsi(3,npart))
	    allocate(ristraold(i)%x(3,npart),ristraold(i)%dpsi(3,npart))
	    allocate(ristranew(i)%x(3,npart),ristranew(i)%dpsi(3,npart))
        allocate(ristranew1(i)%x(3,npart),ristranew1(i)%dpsi(3,npart))
    enddo
    allocate( bisectrate(0:mmax),bisectcount(0:mmax),bisecttot(0:mmax) &
            ,bisectcountl(0:mmax),bisecttotl(0:mmax) &
            ,bisectcountr(0:mmax),bisecttotr(0:mmax) &
            ,bisecthalfcount(0:mmaxhalf),bisecthalftot(0:mmaxhalf) )
   
    allocate(newristra(0:nristra),oldristra(0:nristra))
    do i=0,nristra
        allocate(oldristra(i)%x(3,npart),oldristra(i)%dpsi(3,npart))
        allocate(newristra(i)%x(3,npart),newristra(i)%dpsi(3,npart))
    enddo
    !etrial=etrialin*npart
    !ecut=-log(ecutin)/dt
    ! spin part:
    nprot=nprotin
    nspin=2**npart ! n=4, it is 16.
    call combin(npart,nprot,nisospin) ! give value for nisospin, n=4, it is 4C2= 6.
    nbasis=nspin*nisospin ! n=4, it is 96.   
    allocate(cwtbegin(0:nspin-1,nisospin),cwtend(0:nspin-1,nisospin) &
    ,cwtl2r(0:nchorizo,0:nspin-1,nisospin),cwtr2l(0:nchorizo,0:nspin-1,nisospin))
    allocate(cwtbegin1(0:nspin-1,nisospin),cwtend1(0:nspin-1,nisospin) &
    ,cwtl2r1(0:nchorizo,0:nspin-1,nisospin),cwtr2l1(0:nchorizo,0:nspin-1,nisospin))
    allocate(cwtbeginnew(0:nspin-1,nisospin),cwtendnew(0:nspin-1,nisospin) &
    ,cwtl2rnew(0:nchorizo,0:nspin-1,nisospin),cwtr2lnew(0:nchorizo,0:nspin-1,nisospin))
    allocate(cwtbeginnew1(0:nspin-1,nisospin),cwtendnew1(0:nspin-1,nisospin) &
    ,cwtl2rnew1(0:nchorizo,0:nspin-1,nisospin),cwtr2lnew1(0:nchorizo,0:nspin-1,nisospin))   
    allocate(cwtground(0:nspin-1,nisospin))
    allocate(invspin(nbasis),invispin(nbasis))  
    call getcwtgnd(cwtground,invspin,invispin)
   
    irepstepmonmax=10*nchorizo
    allocate(repstepct1(0:irepstepmonmax),repstepct2(0:irepstepmonmax))
    allocate(repstepct1sum(0:irepstepmonmax),repstepct2sum(0:irepstepmonmax))
    repstepct1=0
    repstepct2=0
    repstepct1sum=0
    repstepct2sum=0  
   
    nav=navin
    nstep=nstepin
    neq=neqin
    
    allocate(xin(nstep),yin(nstep),corr(0:nstep)) ! for correlation check.
    allocate(ibisectstuck(nav+neq))
    ibisectstuck=0 ! check how many times bisection stuck.    
    
    if (nrepmax /= 0) then
        allocate(ireparrowhistory(0:nav*nstep))
        ireparrowhistory=0
    endif
   
    return
    end subroutine stepinit

    subroutine chorizoin(x,i)
    use wavefunction
    integer(kind=i4) :: i
    real(kind=r8) :: x(:,:)
    ristra(i)%x=x
    !call hpsi(ristra(i))
    return
    end subroutine chorizoin

    subroutine chorizoout(x,i)
    integer(kind=i4) :: i
    real(kind=r8) :: x(:,:)
    x=ristra(i)%x
    return
    end subroutine chorizoout

    subroutine chorizoallin(xallin)
    use wavefunction
    integer(kind=i4) :: i
    real(kind=r8) :: xallin(3,npart,0:nchorizo)
    do i=0,nchorizo
    ristra(i)%x=xallin(:,:,i)
    enddo
    return
    end subroutine chorizoallin   
   
    subroutine chorizoallout(xout,ipl,jpl,ipr,jpr)
    integer(kind=i4) :: i
    real(kind=r8) :: xout(3,npart,0:nchorizo)
    integer(kind=i4) :: ipl(6),jpl(6),ipr(6),jpr(6)	
    do i=0,nchorizo
        xout(:,:,i)=ristra(i)%x
    enddo
    call stepoutputordlr(ipl,jpl,ipr,jpr)
    return
    end subroutine chorizoallout   
   
    subroutine stepchorizoinit(isitein,infilein,outfilein) ! initialize the position, calculate the weights for each bead. l2r,r2l
    use wavefunction
    use v6stepcalc
    use random
    use brut
    use mympi
    integer(kind=i4) :: i,j,k,nptot,ixtemp,icheck,l
    real(kind=r8) :: xr(3,npart),xl(3,npart),x(3,npart),x0(3,npart) &
                    ,x0tot(3,npart,0:nchorizo),x0tot1d(3*npart*(nchorizo+1))
    real(kind=r8) ::  rangex
    character(len=70) :: infilein,outfilein  
    logical :: isitein
    real(kind=r8) :: psi20,psi2n
    real(kind=r8) :: vn,vd,rf,vnke,vnpev6,vnpeem
    real(kind=r8), parameter :: dxx=.02d0  
    integer(kind=i4), allocatable :: icheck1d(:)
    integer(kind=i4), allocatable :: iplall(:),jplall(:),iprall(:),jprall(:)	! for mpi_gather
    real(kind=r8), allocatable :: xall(:) ! for mpi_gather
   
    ! call stepinit first   
    isite=isitein
    infile=infilein
    outfile=outfilein
    
    if (myrank()==0) then
       allocate(xall(3*npart*(nchorizo+1)*nproc()))
       allocate(iplall(6*nproc()),jplall(6*nproc()),iprall(6*nproc()),jprall(6*nproc()))    
    endif
    ixtemp=3*npart*(nchorizo+1)  
   
   
    if (isite) then
        if (myrank().eq.0) write (6,'(''input from sites'')')
	    rangex=x0step
	    x0=rangex*reshape((randn(size(x0))-0.5_r8),shape(x0)) ! let all the beads begin from the same x0. good for vmc.
        do i=0,nchorizo
            call chorizoin(x0,i) 
        enddo
        call outputordlr(iplo,jplo,ipro,jpro)
 
        call chorizoallout(x0tot,iplo,jplo,ipro,jpro)      
        call gather(reshape(x0tot,(/3*npart*(nchorizo+1)/)),xall)
        call gather(iplo,iplall)
        call gather(jplo,jplall)
        call gather(ipro,iprall)
        call gather(jpro,jprall)     
       
    else 
           
    if (myrank()==0) then      
        write (6,'(''input from file'',t40,a20)') infile
        open(unit=9,form='formatted',file=trim(infile),position='rewind')
        read(9,'(i10)') nptot
        if (nproc()/=nptot) then  
            if (nproc()<nptot) then
                write (6,'(''# of cores now < # of cores in infile.'',i10, '' vs '',i10  )'),nproc(),nptot   
                write (12,'(''# of cores now < # of cores in infile.'',i10, '' vs '',i10  )'),nproc(),nptot  
                do l=0,nproc()-1
                    read (9,'(6i10)') iplall((6*l+1):(6*(l+1))) 
                    read (9,'(6i10)') jplall((6*l+1):(6*(l+1))) 
                    read (9,'(6i10)') iprall((6*l+1):(6*(l+1))) 
                    read (9,'(6i10)') jprall((6*l+1):(6*(l+1))) 
                    read (9,'(3e15.7)') xall((ixtemp*l+1):(ixtemp*(l+1)))    
                enddo                    
            else
                write (6,'(''# of cores now > # of cores in infile.'',i10, '' vs '',i10  )'),nproc(),nptot 
                write (12,'(''# of cores now > # of cores in infile.'',i10, '' vs '',i10  )'),nproc(),nptot 
                do l=0,nptot-1
                    read (9,'(6i10)') iplall((6*l+1):(6*(l+1))) 
                    read (9,'(6i10)') jplall((6*l+1):(6*(l+1))) 
                    read (9,'(6i10)') iprall((6*l+1):(6*(l+1))) 
                    read (9,'(6i10)') jprall((6*l+1):(6*(l+1))) 
                    read (9,'(3e15.7)') xall((ixtemp*l+1):(ixtemp*(l+1)))    
                enddo                 
                ! add more replica
                do l=nptot,nproc()-1
                   k=mod(l-nptot,nptot)                                    
                    iplall((6*l+1):(6*(l+1)))=iplall((6*k+1):(6*(k+1))) 
                    jplall((6*l+1):(6*(l+1)))=jplall((6*k+1):(6*(k+1)))  
                    iprall((6*l+1):(6*(l+1)))=iprall((6*k+1):(6*(k+1)))  
                    jprall((6*l+1):(6*(l+1)))=jprall((6*k+1):(6*(k+1)))  
                    xall((ixtemp*l+1):(ixtemp*(l+1)))=xall((ixtemp*k+1):(ixtemp*(k+1)))                        
                enddo 
            endif       
        else    
            do l=0,nproc()-1
                read (9,'(6i10)') iplall((6*l+1):(6*(l+1))) 
                read (9,'(6i10)') jplall((6*l+1):(6*(l+1))) 
                read (9,'(6i10)') iprall((6*l+1):(6*(l+1))) 
                read (9,'(6i10)') jprall((6*l+1):(6*(l+1))) 
                read (9,'(3e15.7)') xall((ixtemp*l+1):(ixtemp*(l+1)))    
            enddo               
        endif
        close(9)
    endif   
    
    ! I think, when all the process see this scatter command, they will wait here, until process 0 reach here and give this scatter order. Just like bcast.    
        
        call scatter(iplall,iplo)   
        call scatter(jplall,jplo)
        call scatter(iprall,ipro)
        call scatter(jprall,jpro)
        call scatter(xall,x0tot1d)
        
        call barrier ! just wait until receive scatter info from process 0.
        x0tot=reshape(x0tot1d,shape(x0tot))
        call chorizoallin(x0tot)
        call inputordlr(iplo,jplo,ipro,jpro)     
   
    endif
 
    if (myrank()==0) write (6,'(''output to file'',t40,a20)') outfile    
    xr=ristra(nchorizo)%x
    xl=ristra(0)%x 

    call corop(xl,iplo,jplo,cwtbegin) ! psitcwt    
    call v6propr(xl,-1,cwtbegin,cwtl2r(0,:,:))
    do i=1,nchorizo-1
	    x=ristra(i)%x
	    call v6proplr(x,x,-1,cwtl2r(i-1,:,:),cwtl2r(i,:,:))
    enddo 

    call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))
    ! right to left  

    call corop(xr,ipro,jpro,cwtend)
    call v6propr(xr,-1,cwtend,cwtr2l(nchorizo,:,:))
    do j=nchorizo-1,1,-1
	    x=ristra(j)%x
	    call v6proplr(x,x,-1,cwtr2l(j+1,:,:),cwtr2l(j,:,:))
    enddo
    call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))
    
    call hpsi(ristra(0),ristra(nchorizo),vn,vd,vnke,vnpev6,vnpeem,psi20,psi2n,rf,icheck)  
   
    if (myrank()==0) allocate(icheck1d(0:nproc()-1))   
    call gather(icheck,icheck1d)
    if (myrank()==0) then
     if (sum(icheck1d)/=0) then
        open(unit=9,form='formatted',file=trim(outfile),position='rewind')
        do l=0,nproc()-1
        if ( icheck1d(l) /= 0 ) then
            write(6,'( ''process= '', i10 , '' psi20 and psi2 error! '' )') l
            write (9,'(i10)') l
            write (9,'(6i10)') iplall((6*l+1):(6*(l+1))) 
            write (9,'(6i10)') jplall((6*l+1):(6*(l+1))) 
            write (9,'(6i10)') iprall((6*l+1):(6*(l+1))) 
            write (9,'(6i10)') jprall((6*l+1):(6*(l+1))) 
            write (9,'(3e15.7)') xall((ixtemp*l+1):(ixtemp*(l+1)))               
        endif    
        enddo 
        close(9)  
        write(6,*) 'check hpsi fail! Bug x recorded in outfile! use hpsi to check!'
        call abort         
     else
        write(6,*) 'hpsi pass!'       
     endif  
     deallocate(icheck1d)
     deallocate(xall,iplall,jplall,iprall,jprall)     
    endif
    
    return
    end subroutine stepchorizoinit 
   
    subroutine stepchorizoinitlr 
    !  initialize the position, calculate the weights for each bead. l2r,r2l
    use wavefunction
    use v6stepcalc
    use brut 
    use random
    integer(kind=i4) :: i,j
    real(kind=r8) :: xr(3,npart),xl(3,npart),x(3,npart)
    real(kind=r8) ::  newp,oldp,rn(1)
   
    iclrtot=iclrtot+1
   
    call getordlr(ipl,jpl,ipr,jpr) ! pick the new l and r order. 
  
    ! need to input cwtbegin and cwtend!
    ! left to right
    xl=ristra(0)%x
    call corop(xl,ipl,jpl,cwtbeginnew) ! psit cwt
    !write(6,*) 'coropl finish'
    ! right to left   
    xr=ristra(nchorizo)%x 
    call corop(xr,ipr,jpr,cwtendnew)
    call v6propr(xr,-1,cwtendnew,cwtr2lnew(nchorizo,:,:))
    do j=nchorizo-1,1,-1
	    x=ristra(j)%x
	    call v6proplr(x,x,-1,cwtr2lnew(j+1,:,:),cwtr2lnew(j,:,:))
    enddo
    call v6propl(xl,-1,cwtr2lnew(1,:,:),cwtr2lnew(0,:,:))
   
    newp=abs(real( sum(conjg(cwtbeginnew(:,:))*cwtr2lnew(0,:,:)) )) 	  
    oldp=abs(real( sum(conjg(cwtbegin(:,:))*cwtr2l(0,:,:)) ))   
   
    rn=randn(1)
    
    if (rn(1).lt.(newp/oldp)) then
        
        iclr=iclr+1
           
    ! l2r        
    call v6propr(xl,-1,cwtbeginnew,cwtl2rnew(0,:,:))
    do i=1,nchorizo-1
	    x=ristra(i)%x
	    call v6proplr(x,x,-1,cwtl2rnew(i-1,:,:),cwtl2rnew(i,:,:))
    enddo   
    call v6propl(xr,-1,cwtl2rnew(nchorizo-1,:,:),cwtl2rnew(nchorizo,:,:))

        cwtbegin=cwtbeginnew
        cwtend=cwtendnew
        cwtl2r=cwtl2rnew
        cwtr2l=cwtr2lnew
        iplo=ipl
        jplo=jpl
        ipro=ipr
        jpro=jpr
        call inputordlr(iplo,jplo,ipro,jpro)    
    endif

	
    return
    end subroutine stepchorizoinitlr   
   
    subroutine stepstd(nblocknowin) ! for non-constant trial wavefunction
    use estimator
    use estimatorristra
    use wavefunction
    use random
    use brut
    use mympi
    real(kind=r8) :: rn(1),rm,rmsq,r1,r2,r3,r23,rb
    real(kind=r8) :: v(0:nchorizo),xtmp(3,npart) &
	                ,x0tot(3,npart,0:nchorizo)
    real(kind=r8), dimension(0:nrhobin) :: rhodistout
    integer(kind=i4) :: i,l,ixd,lmax,icheck,ixtemp,nblocknowin
    real(kind=r8) :: vn,vd,vnke,vnpev6,vnpeem,rf,vnum,vdenom &
                    ,sum1,sum2
    real(kind=r8) :: psi20,psi2n
    integer(kind=i4), allocatable :: icheck1d(:)
    real(kind=r8), allocatable :: psi201d(:),psi2n1d(:)
    integer(kind=i4), allocatable :: iplall(:),jplall(:),iprall(:),jprall(:)	! for mpi_gather
    real(kind=r8), allocatable :: xall(:) ! for mpi_gather
   
    nblocknow=nblocknowin
    
    icstepstd=icstepstd+1
   
    if (irep.ne.4) then
 	   
        !call getordlr(iplo,jplo,ipro,jpro)  
       

	    !rn=randn(1)
	    !ipick=min(int(2.*rn(1))+1,2)	   
	    !if (ipick.eq.1) call stdmove1old
	    !if (ipick.eq.2) then
	    ! call stdmove2
	    ! call stdmove3
	    !endif
	    !call stdmove1gauss
	    !call stdmove1
	    !call stdmove1sym
	    !call stdmove1old
	    !call stdmove2
	    !call stdmove3
	    !call stdmove23
	    !call stdmove
	  
       
        rb=0.9
        r1=0.95
        r2=0.967
        r3=0.984
        r23=1
     
        !rb=0.
        !r1=0.
        !r2=1./3.
        !r3=2./3.
        !r23=1.     
     

        rn=randn(1)
        if (rn(1)<=rb) then
            call bisect
        else
    ! different moves are important. especially it seems move 1 by 1 is needed. Otherwise may stuck at bad configs.
        if ((rn(1)>rb).and.(rn(1)<=r1))  call stdmove1fastv6(0,nchorizo)   ! move 1 by 1     
        if ((rn(1)>r1).and.(rn(1)<=r2))  call stdmove2v6(0,nchorizo)   ! move all beads  
        if ((rn(1)>r2).and.(rn(1)<=r3))  call stdmove3v6(0,nchorizo)  ! shift beads
        if ((rn(1)>r3).and.(rn(1)<=r23))  call stdmove23v6(0,nchorizo)  ! move all + shift     
        endif
     
        !call bisectv6improvew11    
 
    else
	    call stdvmcmove2v6 
    endif
        call stepchorizoinitlr       
      
        !if (mod( icstepstd,(nchorizo/nbisect+2) ).eq.0) call stepchorizoinitlr
      
        !call stdmove1rangefastv6(0,nchorizo)
  
	    !call bisectimprove2
	    !call bisect
	    !call bisectold
	 
	    !rn=randn(1)
	    !ipick=min(int(2.*rn(1))+1,2)	   
	    !if (ipick.eq.1) call stdmove1c
	    !if (ipick.eq.2) then
	    ! call stdmove2c
	    ! call stdmove3c
	    !endif 

    ! check correlation steps
    if (icorrchk) then    
        if (myrank().eq.0) then  
            write(6,*) 'icstepstd=', icstepstd,nstep
            if (icstepstd <= nstep) then   ! here nstepdecor usually set the same as nstep. 
   	            xtmp=ristra(nchorizomid)%x
	            call funcrmsq(xtmp,rm,rmsq)
	            xin(icstepstd)=rm	   
                !call hpsi(ristra(0),ristra(nchorizo),vn,vd,vnke,vnpev6,vnpeem,psi20,psi2n,rf,icheck)
                !yin(icstepstd)=vn
	            if (icstepstd.eq.nstep) then
	                ixd=icstepstd
	                lmax=nstep                   
                    call corrchk(xin,ixd,corr,lmax)       
                endif   
            endif
        endif
    endif
     
    ! calculations	  
    if (mod(icstepstd,nstepdecor).eq.0) then   
    
        !write(6,*) 'icstepstd=',icstepstd    
		   
        ! do not need to calculate elocal too often.  
        ice=ice+1
        call hpsi(ristra(0),ristra(nchorizo),vn,vd,vnke,vnpev6,vnpeem,psi20,psi2n,rf,icheck)
     
        !write(6,'( ''stepstd hpsicheck'', 8(g15.7,1x) )') vn,vd,vnke,vnpev6,vnpeem,psi20,psi2n,rf
   

        ! check	hpsi begin. can check to see if just let rank 0 allocate those arrays works or not.
        if (myrank()==0) then
            allocate(xall(3*npart*(nchorizo+1)*nproc()))
            allocate(iplall(6*nproc()),jplall(6*nproc()),iprall(6*nproc()),jprall(6*nproc()))       
            allocate(icheck1d(0:nproc()-1)) 
            allocate(psi201d(0:nproc()-1),psi2n1d(0:nproc()-1))
        endif
        ! when we do this nproc allocation, better check how much memory it consumes at this moment.   
   
        ixtemp=3*npart*(nchorizo+1)  
        call chorizoallout(x0tot,iplo,jplo,ipro,jpro)      
        call gather(reshape(x0tot,(/3*npart*(nchorizo+1)/)),xall)
        call gather(iplo,iplall)
        call gather(jplo,jplall)
        call gather(ipro,iprall)
        call gather(jpro,jprall)    
        call gather(icheck,icheck1d)
        call gather(psi20,psi201d)
        call gather(psi2n,psi2n1d)  
   
        call addall(repstepct1,repstepct1sum)
        call addall(repstepct2,repstepct2sum)   
        sum1=sum(repstepct1sum)
        sum2=sum(repstepct2sum)   
   
        if (myrank().eq.0) then
            if (sum(icheck1d) /= 0) then
                open(unit=19,form='formatted',file='bugbeads.out',position='rewind')
                do l=0,nproc()-1
                    if ( icheck1d(l) /= 0 ) then
                     write(6,'( ''process= '', i10 , '' psi20 and psi2 error! '' )') l
                     write (6,'(3e15.7)') psi201d(l),psi2n1d(l),abs((psi201d(l)-psi2n1d(l))/((psi201d(l)+psi2n1d(l))/2.))
                     write (19,'(i10)') l
                     write (19,'(3e15.7)') psi201d(l),psi2n1d(l),abs((psi201d(l)-psi2n1d(l))/((psi201d(l)+psi2n1d(l))/2.))
                     write (19,'(6i10)') iplall((6*l+1):(6*(l+1))) 
                     write (19,'(6i10)') jplall((6*l+1):(6*(l+1))) 
                     write (19,'(6i10)') iprall((6*l+1):(6*(l+1))) 
                     write (19,'(6i10)') jprall((6*l+1):(6*(l+1))) 
                     write (19,'(3e15.7)') xall((ixtemp*l+1):(ixtemp*(l+1)))               
                    endif    
                enddo 
                close(19)  
                write(6,*) 'check hpsi fail! Bug x recorded in outfile! use hpsi to check!'
                call abort 
            else              
                if (nblocknow > neq) then             
                ! output all the x.  call it path.out
                    open(unit=19,form='unformatted',file='path.unf',position='append')
                    write(19) ipath
                    do l=0,nproc()-1    
                        write (19) l
                        write (19) iplall((6*l+1):(6*(l+1))) 
                        write (19) jplall((6*l+1):(6*(l+1))) 
                        write (19) iprall((6*l+1):(6*(l+1))) 
                        write (19) jprall((6*l+1):(6*(l+1))) 
                        write (19) xall((ixtemp*l+1):(ixtemp*(l+1)))               
                    enddo 
                    close(19)  
                    open(unit=19,form='formatted',file='pathnumbercounts.txt',position='rewind')
                    write (19,'(3i10)') nproc(),ipath,ixtemp
                    close(19)
                    ipath=ipath+1              
                else                
                    open(unit=19,form='unformatted',file='path0.unf',position='append')
                    write(19) ipathtmp
                    do l=0,nproc()-1    
                        write (19) l
                        write (19) iplall((6*l+1):(6*(l+1))) 
                        write (19) jplall((6*l+1):(6*(l+1))) 
                        write (19) iprall((6*l+1):(6*(l+1))) 
                        write (19) jprall((6*l+1):(6*(l+1))) 
                        write (19) xall((ixtemp*l+1):(ixtemp*(l+1)))               
                    enddo 
                    close(19)  
                    open(unit=19,form='formatted',file='path0numbercounts.txt',position='rewind')
                    write (19,'(3i10)') nproc(),ipathtmp,ixtemp
                    close(19)
                    ipathtmp=ipathtmp+1                                             
                endif         
                if (nrepmax /= 0) then
                    open(unit=19,form='formatted',file='reptationMON.txt',position='rewind')
                ! check reptation number of steps distributions.
                    do i=0,irepstepmonmax
                        write(19,'( i10 , i10, f15.7, i10, f15.7 )') i,repstepct1sum(i),repstepct1sum(i)/sum1 &
                                                                    ,repstepct2sum(i),repstepct2sum(i)/sum2 
                    enddo 
                    close(19)     
                    open(unit=19,form='formatted',file='reparrowhist.txt',position='append')
                ! check reptation direction history.
                    do i=icrephistlast,icrephist
                        write(19,'( i10,1x,i10 )') i,ireparrowhistory(i)
                    enddo 
                    icrephistlast=icrephist+1
                    close(19) 
                endif               
            endif
            deallocate(xall,iplall,jplall,iprall,jprall,icheck1d,psi201d,psi2n1d)   
        endif
    
        ! check hpsi end 
        vnum=vn
        vdenom=vd
        ! update the v beads   
        do i=0,nchorizo
        call vtot(i,v(i))
        enddo
        v=v/rf
        call addvalristra(1,v) ! do not need this part   
   
        call addval(3,vnum,1.0_r8)
        call addval(4,vdenom,1.0_r8)
        call addval(5,rf,1.0_r8)
 
        !do i=0,nproc()-1
        !   if (myrank().eq.i) then
        !   write(6,'( ''myrank='', i10, 1x, 3(g15.7,1x) )')  myrank(),vnum,vdenom,rf
        !   endif
        !   call barrier
        !enddo       
        !!! irn=84, bug bead occurs. focus on this.     
      
        call funcrmsq(ristra(nchorizomid)%x,rm,rmsq)
        call addval(6,rmsq,1.0_r8)
        call addval(7,rm,1.0_r8)
        !write(6,*) 'nchorizomid=',nchorizomid
        call samplerhodistrbt(ristra(nchorizomid)%x,nchorizomid,rhodistout)	  
        call addval(8,rhodistout(0),1.0_r8)

        if (myrank().eq.0)  then     
            open (unit=11,FILE='values.txt',FORM='FORMATTED',position='append')    
            write(11,'(10(G15.7,2x))') vnum,vdenom,vnke,vnpev6,vnpeem,rm,rhodistout(0) !11 means write to 'values.txt'
            close (11)   
        endif  

    endif
   
    return
    end subroutine stepstd 

    subroutine reptation 
    use wavefunction
    use v6stepcalc
    use brut 
    use random
    use mympi
    use estimator
    integer(kind=i4) :: i,j
    real(kind=r8) :: xr(3,npart),xl(3,npart),x(3,npart),gauss(3,npart)
    real(kind=r8) ::  newp,oldp,rn(1)
    
    icrephist=icrephist+1
    icreptot=icreptot+1
    ristraold(0:nchorizo)=ristra(0:nchorizo)
    ! this is the more general case, nrepmax do not have to be 1.    
    ! nreparrow: -1 means to the left, +1 means to the right. 
    if (nreparrow == 1) then
    ! move right    
       
        ristranew(0:nchorizo-nrepmax)=ristraold(nrepmax:nchorizo)   
           
    !do i=nrepmax,nchorizo
    !   ristranew(i-nrepmax)=ristraold(i)    
    !enddo

    do i=nchorizo-nrepmax+1,nchorizo
        gauss=sigma*reshape(gaussian(3*npart),(/3,npart/))
        ristranew(i)%x=ristranew(i-1)%x+gauss
    enddo
    else
    if (nreparrow /= -1) then
        write(6,*) 'reptation direction error! Stop!',nreparrow
        call abort
    endif
    ! move left

        ristranew(nrepmax:nchorizo)=ristraold(0:nchorizo-nrepmax)   
    
    !if (myrank() == 0) then     
    !  open(unit=19,form='formatted',file='checkbeadsristrawritten1.txt',position='rewind')  
    !  do i=0,nchorizo-nrepmax
    !     write (19,'(6i10)') i
    !     write (19,'(3e15.7)') ristraold(i)%x-ristranew(i+nrepmax)%x         
    !  enddo
    !  close(19)    
    !  
    !  call abort
    !endif
    
    
    !do i=0,nchorizo-nrepmax
    !   ristranew(i+nrepmax)=ristraold(i)   
    !enddo

    !if (myrank() == 0) then     
    !  open(unit=19,form='formatted',file='checkbeadsristrawritten2.txt',position='rewind')  
    !  do i=0,nchorizo-nrepmax
    !     write (19,'(6i10)') i
    !     write (19,'(3e15.7)') ristraold(i)%x-ristranew(i+nrepmax)%x         
    !  enddo
    !  close(19)     
    !  
    !  call abort
    !  
    !endif    
 
    do i=nrepmax-1,0,-1
        gauss=sigma*reshape(gaussian(3*npart),(/3,npart/))
        ristranew(i)%x=ristranew(i+1)%x+gauss
    enddo       
    endif
   
    ! need to input cwtbegin and cwtend!
    ! left to right
    xl=ristranew(0)%x
    call corop(xl,iplo,jplo,cwtbeginnew) ! psit cwt
    !write(6,*) 'coropl finish'
    ! right to left   
    xr=ristranew(nchorizo)%x 
    call corop(xr,ipro,jpro,cwtendnew)
    call v6propr(xr,-1,cwtendnew,cwtr2lnew(nchorizo,:,:))
    do j=nchorizo-1,1,-1
	    x=ristranew(j)%x
	    call v6proplr(x,x,-1,cwtr2lnew(j+1,:,:),cwtr2lnew(j,:,:))
    enddo
    call v6propl(xl,-1,cwtr2lnew(1,:,:),cwtr2lnew(0,:,:))
   
    newp=abs(real( sum(conjg(cwtbeginnew(:,:))*cwtr2lnew(0,:,:)) )) 	  
    oldp=abs(real( sum(conjg(cwtbegin(:,:))*cwtr2l(0,:,:)) ))   
   
    rn=randn(1)
    
    !write (6,*) 'rn=, newp/oldp=',rn(1), newp/oldp
    
    if (rn(1).lt.(newp/oldp)) then
        icrep=icrep+1
	    call addval(1,1.0_r8,1.0_r8)   
        ristra(0:nchorizo)=ristranew(0:nchorizo) 
     
        irepstepct=irepstepct+nrepmax
         
        if ( nreparrow ==1 ) then 
        ireparrowhistory(icrephist)=ireparrowhistory(icrephist-1)+nrepmax   
        else
        ireparrowhistory(icrephist)=ireparrowhistory(icrephist-1)-nrepmax   
        endif
   
    ! l2r        
        call v6propr(xl,-1,cwtbeginnew,cwtl2rnew(0,:,:))
        do i=1,nchorizo-1
	    x=ristra(i)%x
	    call v6proplr(x,x,-1,cwtl2rnew(i-1,:,:),cwtl2rnew(i,:,:))
        enddo   
        call v6propl(xr,-1,cwtl2rnew(nchorizo-1,:,:),cwtl2rnew(nchorizo,:,:))
        cwtbegin=cwtbeginnew
        cwtend=cwtendnew
        cwtl2r=cwtl2rnew
        cwtr2l=cwtr2lnew
                
    else
	    call addval(1,0.0_r8,1.0_r8) 
     
        if (irepstepct > irepstepmonmax) then
            write(6,*) 'stop for now bc irepstepct is very big, it is =', irepstepct
            call abort
        endif     
        if (nreparrow == 1) then 
        repstepct1(irepstepct)=repstepct1(irepstepct)+1     
        else
        repstepct2(irepstepct)=repstepct2(irepstepct)+1    
        endif
        irepstepct=0
        nreparrow=-nreparrow  ! change reptation direction. 
        !write(6,*) 'reptation reject once on proc ', myrank()
        ireparrowhistory(icrephist)=ireparrowhistory(icrephist-1)
     

    endif   
   
    
    return
    end subroutine reptation
   
    subroutine reptation1 
    ! move 1 bead, small importance sampling is added. nrepmax=1 case.
    use wavefunction
    use v6stepcalc
    use brut 
    use random
    use mympi
    use estimator
    integer(kind=i4) :: i,j
    real(kind=r8) :: xr(3,npart),xl(3,npart),x(3,npart),gauss(3,npart),xr1(3,npart),xl1(3,npart)
    real(kind=r8) ::  newp,oldp,rn(1),newp1,oldp1 
   
    icrephist=icrephist+1
    icreptot=icreptot+1
    ristraold(0:nchorizo)=ristra(0:nchorizo)
        
    ! nreparrow: -1 means to the left, +1 means to the right. 
    if (nreparrow == 1) then
    ! move right        
    ristranew(0:nchorizo-1)=ristraold(1:nchorizo)   
    ristranew1(0:nchorizo-1)=ristranew(0:nchorizo-1)      
    !do i=nrepmax,nchorizo
    !   ristranew(i-nrepmax)=ristraold(i)    
    !enddo
    gauss=sigma*reshape(gaussian(3*npart),(/3,npart/))
    ristranew(nchorizo)%x=ristranew(nchorizo-1)%x+gauss ! + direction.
    ristranew1(nchorizo)%x=ristranew1(nchorizo-1)%x-gauss ! ristranew1 is the last bead - direction.

    ! l2r    
    xl=ristranew(0)%x
    call corop(xl,iplo,jplo,cwtbeginnew) ! psit cwt    
    call v6propr(xl,-1,cwtbeginnew,cwtl2rnew(0,:,:))
    do i=1,nchorizo-1
	    x=ristranew(i)%x
	    call v6proplr(x,x,-1,cwtl2rnew(i-1,:,:),cwtl2rnew(i,:,:))
    enddo 
    cwtl2rnew1(0:nchorizo-1,:,:)=cwtl2rnew(0:nchorizo-1,:,:)
    xr=ristranew(nchorizo)%x 
    xr1=ristranew1(nchorizo)%x 
    call v6propl(xr,-1,cwtl2rnew(nchorizo-1,:,:),cwtl2rnew(nchorizo,:,:)) 
    call v6propl(xr1,-1,cwtl2rnew1(nchorizo-1,:,:),cwtl2rnew1(nchorizo,:,:))     
    call corop(xr,ipro,jpro,cwtendnew)  
    call corop(xr1,ipro,jpro,cwtendnew1)  

    newp=abs(real( sum(conjg(cwtl2rnew(nchorizo,:,:))*cwtendnew(:,:)) )) 	 
    newp1=abs(real( sum(conjg(cwtl2rnew1(nchorizo,:,:))*cwtendnew1(:,:)) )) 
    
    oldp=abs(real( sum(conjg(cwtbegin(:,:))*cwtr2l(0,:,:)) ))   
    
    xl1=2.*ristra(1)%x-ristra(0)%x
    call corop(xl1,iplo,jplo,cwtbegin1) 
    call v6propl(xl1,-1,cwtr2l(1,:,:),cwtr2l1(0,:,:))    
    
    oldp1=abs(real( sum(conjg(cwtbegin1(:,:))*cwtr2l1(0,:,:)) ))   
    
    rn=randn(1)
    if ( rn(1) < ((newp+newp1)/(oldp+oldp1)) ) then    
        icrep=icrep+1
	    call addval(1,1.0_r8,1.0_r8)  
        irepstepct=irepstepct+1
        ireparrowhistory(icrephist)=ireparrowhistory(icrephist-1)+1
     
     
        rn=randn(1)
        if ( rn(1) < (newp/(newp+newp1)) ) then
        ristra(0:nchorizo)=ristranew(0:nchorizo)       
    ! r2l        
        call v6propr(xr,-1,cwtendnew,cwtr2lnew(nchorizo,:,:))
        do j=nchorizo-1,1,-1
	        x=ristra(j)%x
	        call v6proplr(x,x,-1,cwtr2lnew(j+1,:,:),cwtr2lnew(j,:,:))
        enddo
        call v6propl(xl,-1,cwtr2lnew(1,:,:),cwtr2lnew(0,:,:))  
        cwtbegin=cwtbeginnew
        cwtend=cwtendnew
        cwtl2r=cwtl2rnew 
        cwtr2l=cwtr2lnew      
        else
        ristra(0:nchorizo)=ristranew1(0:nchorizo)   
    ! r2l 
        call v6propr(xr1,-1,cwtendnew1,cwtr2lnew1(nchorizo,:,:))
        do j=nchorizo-1,1,-1
	        x=ristra(j)%x
	        call v6proplr(x,x,-1,cwtr2lnew1(j+1,:,:),cwtr2lnew1(j,:,:))
        enddo
        call v6propl(xl,-1,cwtr2lnew1(1,:,:),cwtr2lnew1(0,:,:))     
        cwtbegin=cwtbeginnew
        cwtend=cwtendnew1
        cwtl2r=cwtl2rnew1  
        cwtr2l=cwtr2lnew1         
        endif       
    else   
        call addval(1,0.0_r8,1.0_r8)   
     
        if (irepstepct > irepstepmonmax) then
            write(6,*) 'stop for now bc irepstepct is very big, it is =', irepstepct
            call abort
        endif 
     
        if (nreparrow == 1) then
        repstepct1(irepstepct)=repstepct1(irepstepct)+1     
        else
        repstepct2(irepstepct)=repstepct2(irepstepct)+1    
        endif
        irepstepct=0
        nreparrow=-nreparrow     
        ireparrowhistory(icrephist)=ireparrowhistory(icrephist-1)
    endif
  
    else
       
    if (nreparrow /= -1) then
        write(6,*) 'reptation direction error! Stop!',nreparrow
        call abort
    endif
    ! move left, the 0th bead changes.
  
        ristranew(1:nchorizo)=ristraold(0:nchorizo-1)     
        ristranew1(1:nchorizo)=ristranew(1:nchorizo)
        gauss=sigma*reshape(gaussian(3*npart),(/3,npart/))
        ristranew(0)%x=ristranew(1)%x+gauss
        ristranew1(0)%x=ristranew1(1)%x-gauss
      
    ! right to left.    
    xr=ristranew(nchorizo)%x
    xl=ristranew(0)%x
    xl1=ristranew1(0)%x
    call corop(xr,ipro,jpro,cwtendnew)     
    call v6propr(xr,-1,cwtendnew,cwtr2lnew(nchorizo,:,:))
    do j=nchorizo-1,1,-1
        x=ristranew(j)%x
        call v6proplr(x,x,-1,cwtr2lnew(j+1,:,:),cwtr2lnew(j,:,:))
    enddo
    cwtr2lnew1(1:nchorizo,:,:)=cwtr2lnew(1:nchorizo,:,:)
    call v6propl(xl,-1,cwtr2lnew(1,:,:),cwtr2lnew(0,:,:))      
    call v6propl(xl1,-1,cwtr2lnew(1,:,:),cwtr2lnew1(0,:,:))    

    call corop(xl,iplo,jplo,cwtbeginnew)
    call corop(xl1,iplo,jplo,cwtbeginnew1)   
    
    newp=abs(real( sum(conjg(cwtbeginnew(:,:))*cwtr2lnew(0,:,:)) ))  
    newp1=abs(real( sum(conjg(cwtbeginnew1(:,:))*cwtr2lnew1(0,:,:)) )) 
    
    oldp=abs(real( sum(conjg(cwtbegin(:,:))*cwtr2l(0,:,:)) ))   
    
    xr1=2.*ristra(nchorizo-1)%x-ristra(nchorizo)%x
    call corop(xr1,ipro,jpro,cwtend1)   
    call v6propl(xr1,-1,cwtl2r(nchorizo-1,:,:),cwtl2r1(nchorizo,:,:))   
    
    oldp1=abs(real( sum(conjg(cwtl2r1(nchorizo,:,:))*cwtend1(:,:)) ))   
    
    rn=randn(1)
    if ( rn(1) < ((newp+newp1)/(oldp+oldp1)) ) then    
        icrep=icrep+1
	    call addval(1,1.0_r8,1.0_r8)    
     
        irepstepct=irepstepct+1     
        ireparrowhistory(icrephist)=ireparrowhistory(icrephist-1)-1
     
        rn=randn(1)
        if ( rn(1) < (newp/(newp+newp1)) ) then
        ristra(0:nchorizo)=ristranew(0:nchorizo)              
    ! l2r
        call v6propr(xl,-1,cwtbeginnew,cwtl2rnew(0,:,:))
        do i=1,nchorizo-1
	        x=ristra(i)%x
	        call v6proplr(x,x,-1,cwtl2rnew(i-1,:,:),cwtl2rnew(i,:,:))
        enddo   
        call v6propl(xr,-1,cwtl2rnew(nchorizo-1,:,:),cwtl2rnew(nchorizo,:,:))       
        cwtbegin=cwtbeginnew
        cwtend=cwtendnew
        cwtl2r=cwtl2rnew 
        cwtr2l=cwtr2lnew   
  
        else
     
        ristra(0:nchorizo)=ristranew1(0:nchorizo)   
    ! l2r
        call v6propr(xl1,-1,cwtbeginnew1,cwtl2rnew1(0,:,:))
        do i=1,nchorizo-1
	        x=ristra(i)%x
	        call v6proplr(x,x,-1,cwtl2rnew1(i-1,:,:),cwtl2rnew1(i,:,:))
        enddo   
        call v6propl(xr,-1,cwtl2rnew1(nchorizo-1,:,:),cwtl2rnew1(nchorizo,:,:))       
        cwtbegin=cwtbeginnew1
        cwtend=cwtendnew
        cwtl2r=cwtl2rnew1 
        cwtr2l=cwtr2lnew1   
  
        endif       
    else   
        call addval(1,0.0_r8,1.0_r8)   
     
        if (irepstepct > irepstepmonmax) then
            write(6,*) 'stop for now bc irepstepct is very big, it is =', irepstepct
            call abort
        endif 
     
        if (nreparrow == 1) then
        repstepct1(irepstepct)=repstepct1(irepstepct)+1     
        else
        repstepct2(irepstepct)=repstepct2(irepstepct)+1    
        endif
        irepstepct=0     
        nreparrow=-nreparrow    
        ireparrowhistory(icrephist)=ireparrowhistory(icrephist-1)
    endif    
     
    endif
 
    return
    end subroutine reptation1   
   
   
    subroutine rqmcmonitor
    end subroutine rqmcmonitor
   
   
    subroutine pathprocess
    end subroutine pathprocess
   
   
   
   
   
   
   
    subroutine bisect
    use random
    integer(kind=i4) :: ileftbisect,irightbisect
 
    icbisecttot=icbisecttot+1   
   
    if (nbisect.lt.nchorizo) then
    call bisectpicksliceswide(ileftbisect,irightbisect) ! pick the ileftbisect here
    
        !write(12,*) 'ileft,iright=', ileftbisect,irightbisect,nbisect,nbisecthalf
        !write(6,*) 'ileft,iright=', ileftbisect,irightbisect,nbisect,nbisecthalf
    
    if (ileftbisect.le.(nchorizo-nbisect)) then
     
        !call bisectv6improve(ileftbisect)
        !call bisectv6a(ileftbisect)
        call bisectv6b(ileftbisect)
        !write(6,*) 'bisectv6improve called =', ileftbisect
     
    else

        !rn=randn(1)
        !if (rn(1) <= 0.5) then
            call bisectv6lb   
        !else
            call bisectv6rb  
        !endif
          
    !!! reptation        
    !    if (nrepmax == 1) then
    !     !call reptation1 
    !        call reptation
    !    else 
    !        if (nrepmax /= 0) then
    !            call reptation  
    !        endif
    !    endif
    !!! reptation, consider a bisection-ike reptation.
   

    ! below are mixed move, turn off now.  
    !----------------------------------------------      
        !rn=randn(1)
        !if (rn(1)<=0.99999) then  
        
    ! new improved  before 2019-4-7        
        !if (ileftbisect.le.(nchorizo-nbisecthalf)) then    
        ! !write(6,*) 'ileftbisect -', ileftbisect    
        ! !call bisectv6improve(nchorizo-nbisect)
        ! call bisectv6a(nchorizo-nbisect)
        ! !write(6,*) 'bisectv6r called =', ileftbisect
        ! call bisectv6ra
        ! !write(6,*) 'bisectv6l called =', ileftbisect 
        ! call bisectv6la      
        !else
        ! !write(6,*) 'ileftbisect +', ileftbisect
        ! !write(6,*) 'bisectv6r called =', ileftbisect   
        ! call bisectv6ra
        ! !write(6,*) 'bisectv6l called =', ileftbisect 
        ! call bisectv6la 
        ! !call bisectv6improve(0) 
        ! call bisectv6a(0)
        !endif           
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
         
    ! half lr bisection        
        !if (ileftbisect.le.(nchorizo-nbisecthalf)) then    
        ! !write(6,*) 'ileftbisect -', ileftbisect    
        ! call bisectv6improve(nchorizo-nbisect)
        ! !write(6,*) 'bisectv6r called =', ileftbisect
        ! call bisectv6rhalf(nchorizo)
        ! call bisectv6lhalf(0)      
        !else
        ! !write(6,*) 'ileftbisect +', ileftbisect 
        ! call bisectv6rhalf(nchorizo)
        !! write(6,*) 'bisectv6l called =', ileftbisect 
        ! call bisectv6lhalf(0) 
        ! call bisectv6improve(0) 
        !endif  
    !     

    ! old lr bisection        
        !if (ileftbisect.le.(nchorizo-nbisecthalf)) then    
        ! !write(6,*) 'ileftbisect -', ileftbisect    
        ! call stdmove1rangefastv6(nchorizo,nchorizo)  
        ! !write(6,*) 'bisectv6r called =', ileftbisect
        ! call bisectv6r
        ! call bisectv6l      
        !else
        ! !write(6,*) 'ileftbisect +', ileftbisect
        ! call stdmove1rangefastv6(0,0)
        ! call bisectv6r
        ! !write(6,*) 'bisectv6l called =', ileftbisect 
        ! call bisectv6l 
        !endif  
    !        
      
      
        !else   

    ! move 1 by 1   
	    ! if (ileftbisect.lt.nchorizo) then ! nbisect should be >= 2.
	    !n1=nchorizo-ileftbisect
	    !n2=nbisect-1-n1 
    !      !write (12,*) 'n1, n2=', n1,n2
    !      !write (6,*) 'n1, n2=', n1,n2
    !      !call stdmove1rangefastv6(ileftbisect+1,nchorizo)	
    !      call stdmove1fastv6(ileftbisect+1,nchorizo)     
	    !if (n2.gt.0) then    
	    !!call stdmove1rangefastv6(0,n2-1)
    !      call stdmove1fastv6(0,n2-1)
    !      endif
    !    else ! ileftbisect.eq.nchorizo	 
    !      !call stdmove1rangefastv6(0,nbisect-2)	
    !      call stdmove1fastv6(0,nbisect-2)
    !    endif	         
    !     
        !endif 
    !------------------------------------------------   
    endif
    
    else
	    stop 'nbisect >= nchorizo, stop!'
    endif

    return 
    end subroutine bisect

    subroutine bisectv6b(ileftbisect) 
    ! based on v6a version. move beads first, then do an overall judge.
    use wavefunction
    use estimator
    use random
    use v6stepcalc
    integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
					    ,k,bilvl &
	                    ,dtexpt
    real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    real(kind=r8) :: xmid(3,npart),xr(3,npart),xl(3,npart)
    real(kind=r8) :: x(3,npart),xnew(3,npart)
    real(kind=r8) :: tmp 
    complex(kind=r8) :: oldlv(0:mmax),newlv(0:mmax)
    complex(kind=r8) :: cwtl2rtmp1(0:nbisect,0:nspin-1,nisospin),cwtr2ltmp1(0:nbisect,0:nspin-1,nisospin) &
	                    ,cwtl2rtmp2(0:nbisect,0:nspin-1,nisospin),cwtr2ltmp2(0:nbisect,0:nspin-1,nisospin) &
	                    ,cwtbegintmp(0:nspin-1,nisospin),cwtendtmp(0:nspin-1,nisospin)
    logical :: reject	
 
    ibisecttot=ibisecttot+1  
   
    !do i=0,nbisect
	    !oldristra(i)=ristra(ileftbisect+i)
    !enddo
   
    oldristra(0:nbisect)=ristra(ileftbisect:ileftbisect+nbisect)
   
    newristra(0)=oldristra(0)
    newristra(nbisect)=oldristra(nbisect)

   
    !do i=0,nbisect
    !xoldtmp(i,:,:)=oldristra(i)%x
    !enddo
     
    cwtbegintmp(:,:)=cwtl2r(ileftbisect,:,:)
    cwtendtmp(:,:)=cwtr2l(ileftbisect+nbisect,:,:)  
   
    ! tmp1 deal with oldlv
    cwtl2rtmp1(0:nbisect,:,:)=cwtl2r(ileftbisect:ileftbisect+nbisect,:,:)
    cwtr2ltmp1(0:nbisect,:,:)=cwtr2l(ileftbisect:ileftbisect+nbisect,:,:)
    ! tmp2 deal with newlv   
    cwtl2rtmp2(0:nbisect,:,:)=cwtl2r(ileftbisect:ileftbisect+nbisect,:,:)
    cwtr2ltmp2(0:nbisect,:,:)=cwtr2l(ileftbisect:ileftbisect+nbisect,:,:)
   
    oldlv=1.
    newlv=1.   
    
    tau=dt*nbisect ! total time for the bisection slice.

    ! up to here, newristra has been loaded.
   
    ! judge if we accept the newristra or not.
    ! in this v6 case, v is not justpotential v, it should be the value of < ... >. lv.
        
    do bilvl=1,mmax  ! besection.  level 1 to N. mmax = the total level N.	 
        ! move the current level   
        sigmamid=sqrt(hbar*tau*2.0_r8**(-bilvl))  
	    !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	    do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		    imid=2**(mmax-bilvl)+i*2**(mmax-bilvl+1)
		    il=imid-2**(mmax-bilvl)
		    ir=imid+2**(mmax-bilvl)
		    gauss=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		    xmid=(newristra(il)%x+newristra(ir)%x)/2
		    newristra(imid)%x=xmid+gauss  
            !write(12,*) 'bilvl il ir imid sigmamid =',bilvl, il,ir,imid,sigmamid
            !write(12,*) 'gauss =', gauss
        enddo         
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!    
    ! then judge the last level
	    bisecttot(mmax)=bisecttot(mmax)+1	   
	    lmax=2**mmax-1 
	    jd=1 ! interval
	    dtexpt=-1 ! for the R L part, that's why the extra -1. dt=2**dtexp
	    do l=1,lmax
		    j=l*jd
		    jp=(l-1)*jd
		    !xold=oldristra(j)%x
            xnew=newristra(j)%x
		    ! call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(jp,:,:),cwtl2rtmp1(j,:,:))    ! should be lr, not rl I believe.
		    call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(jp,:,:),cwtl2rtmp2(j,:,:))  	   
        enddo
	    jmax=lmax*jd
	    newlv(mmax)=abs(real( sum(conjg(cwtl2rtmp2(jmax,:,:))*cwtendtmp(:,:)) )) 	
	    oldlv(mmax)=abs(real( sum(conjg(cwtl2rtmp1(jmax,:,:))*cwtendtmp(:,:)) ))	   
	    tmp=log(newlv(mmax))-log(oldlv(mmax)) 
  	    !tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
	    rn=randn(1)	   
    ! can check if l2r and r2l give	the same lv.		    
            !write(12,*) 'tmp ln rn=',tmp,log(rn(1))
	    if (log(rn(1)).lt.tmp) then
	        reject=.false. 
	        bisectcount(mmax)=bisectcount(mmax)+1	 
	    else
	        reject=.true. 
            !write(6,*) 'rn=',rn(1),log(rn(1))
            !write(6,*) 'tmp=',tmp,newlv(mmax),oldlv(mmax),log(newlv(mmax)),log(oldlv(mmax)) 	
            
        endif 
		    !write(6,*) 'imid=',imid	
    !!!!!!!!!!!!!!!!!!!!!!!!!!!        
        
    if (reject) then
	    call addval(2,0.0_r8,1.0_r8)
    else
	    ibisect=ibisect+1 
	    call addval(2,1.0_r8,1.0_r8)  
 
	    !do i=0,nbisect
    !       ristra(ileftbisect+i)=newristra(i)
    !   enddo	 
        ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good. 04/10/2019£¬ remove pointer, problem fixed I tihnk.	 

        xl=ristra(0)%x     
	    xr=ristra(nchorizo)%x     
    ! update cwtl2r, cwtr2l.	 
    ! update cwtl2r for all beads. this can be optimized in the future.
     
	    cwtl2r(ileftbisect+1:ileftbisect+nbisect-1,:,:)=cwtl2rtmp2(1:jmax,:,:) ! here jmax should already be 2**mmax-1 
	    if ((ileftbisect+nbisect).le.nchorizo-1) then
        do k=ileftbisect+nbisect,nchorizo-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtl2r(k-1,:,:),cwtl2r(k,:,:))
        enddo
	    call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))
        else
        !ileftbisect+nbisect=nchorizo  
	    call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))   
        endif
   
    ! update cwtr2l for all the beads
	    do k=ileftbisect+nbisect-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
        call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))  	

    endif  

    rejectbisectv6improve=reject
    
    return  
    end subroutine bisectv6b            
   
    subroutine bisectv6lb 
    ! b version, move the whole nbiect beads instead of bisecthalf.
    ! should be better than bisectv6l,rhalf     
    ! move the left end, 0th bead,ileftbiect=0; if ileftbisect /= 0 then the small modification on 0th level is needed.
    use wavefunction
    use estimator
    use random
    use v6stepcalc
    use brut
    integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
					    ,irightbisect,k,bilvl &
	                    ,dtexpt &
                        ,nbisectnow,mmaxnow
    real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    real(kind=r8) :: xmid(3,npart),xr(3,npart),xl(3,npart)
    real(kind=r8) :: x(3,npart),xnew(3,npart),x0new(3,npart)
    real(kind=r8) :: tmp 
    complex(kind=r8) :: oldlv(0:mmax),newlv(0:mmax)
    complex(kind=r8) :: cwtl2rtmp1(0:nbisect,0:nspin-1,nisospin) &
	                    ,cwtl2rtmp2(0:nbisect,0:nspin-1,nisospin) &
	                    ,cwtendtmp(0:nspin-1,nisospin) &
                        ,cwtbegintmp2(0:nspin-1,nisospin)        
    logical :: reject	

    ibisecttotl=ibisecttotl+1  
   
    ileftbisect=0
    irightbisect=ileftbisect+nbisect
   
    !nbisectnow=irightbisect-ileftbisect ! has to be 2^n.
    !mmaxnow=log(real(nbisectnow))/log(real(2)) ! mmaxnow = mmax-1
   
    nbisectnow=nbisect
    mmaxnow=mmax
   
    tau=dt*nbisectnow ! total time for the bisection slice. nisectnow is usually nbisecthalf
    sigmamid=sqrt(2.*hbar*tau)  

    oldristra(0:nbisectnow)=ristra(ileftbisect:ileftbisect+nbisectnow)
 
    ! deal with the two ends first, then do the classic bisection.   
   
    newristra(nbisectnow)=oldristra(nbisectnow)   
    cwtendtmp(:,:)=cwtr2l(ileftbisect+nbisectnow,:,:)
 
    newristra(0)%x=oldristra(nbisectnow)%x+sigmamid*reshape(gaussian(3*npart),(/3,npart/))
   
    x0new=newristra(0)%x
   
    ! i,jpro,i,jplo, the o means old are actually the current order.
   
    call corop(x0new,iplo,jplo,cwtbegintmp2)     
   
    oldlv(0)=1.
    newlv(0)=1. ! it seems this does not need to be set to 1. check this.

    cwtl2rtmp1(0:nbisect,:,:)=cwtl2r(ileftbisect:ileftbisect+nbisect,:,:)
    call v6propr(x0new,-1,cwtbegintmp2,cwtl2rtmp2(0,:,:)) 
   
    ! prepare the further new bisection positions
    do bilvl=1,mmaxnow   ! level 1 to N. mmax = the total level N.
	    sigmamid=sqrt(hbar*tau*2.0_r8**(-bilvl))  
	    !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	    do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		    imid=2**(mmaxnow-bilvl)+i*2**(mmaxnow-bilvl+1)
		    il=imid-2**(mmaxnow-bilvl)
		    ir=imid+2**(mmaxnow-bilvl)
		    !write(6,*) 'bilvl il ir imid sigmamid=',bilvl, il,ir,imid,sigmamid
		    gauss=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		    xmid=(newristra(il)%x+newristra(ir)%x)/2.
		    newristra(imid)%x=xmid+gauss
	    enddo 
    enddo 

        bilvl=mmaxnow  ! besection.  level 1 to N. mmax = the total level N.	   
	    bisecttotl(bilvl)=bisecttotl(bilvl)+1.	   
	    lmax=2**bilvl-1 
	    jd=2**(mmaxnow-bilvl) ! interval
	    dtexpt=mmaxnow-bilvl-1 ! for the R L part, that's why the extra -1. dt=2**dtexp 
	    do l=1,lmax
		    j=l*jd
		    jp=(l-1)*jd
		    !xold=oldristra(j)%x
            xnew=newristra(j)%x
		    !call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(jp,:,:),cwtl2rtmp1(j,:,:))    ! should be lr, not rl I believe.
		    call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(jp,:,:),cwtl2rtmp2(j,:,:))  	   
		    !if ( sum((xnew-xold)**2).eq.0) then
			    !write(6,*) 'bilvl=',bilvl	
			    !write(6,*) 'l=',l,j,jp
			    !write(6,*) 'sigmamid=',sigmamid
		    !   write(6,*) 'xnew-xold=',xnew-xold
		    !endif   
	    enddo
	    jmax=lmax*jd  
	    oldlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp1(jmax,:,:))*cwtendtmp(:,:)) ))	
	    newlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp2(jmax,:,:))*cwtendtmp(:,:)) )) 
	    tmp=log(newlv(bilvl))-log(oldlv(bilvl))
  	    !tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
	    rn=randn(1)	   
    ! can check if l2r and r2l give	the same lv.		    
	        !write(6,*) 'tmp ln rn=',tmp,log(rn(1))
	    if (log(rn(1)).lt.tmp) then
	        reject=.false. 
            bisectcountl(bilvl)=bisectcountl(bilvl)+1.
	    else
	        reject=.true.  
	    endif   
 
    if (reject) then
        call addval(10,0.0_r8,1.0_r8)   
    else
        call addval(10,1.0_r8,1.0_r8)  
	    ibisectl=ibisectl+1 
	    do i=0,nbisectnow
            ristra(ileftbisect+i)=newristra(i)
        enddo	 
    !ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good.	 
        xl=ristra(0)%x     
	    xr=ristra(nchorizo)%x 
    ! update cwtl2r, cwtr2l.	 
        cwtbegin=cwtbegintmp2
	    cwtl2r(ileftbisect:ileftbisect+nbisectnow-1,:,:)=cwtl2rtmp2(0:jmax,:,:) ! here jmax should already be 2**mmax-1    
	    if ((ileftbisect+nbisectnow).le.nchorizo-1) then
        do k=ileftbisect+nbisectnow,nchorizo-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtl2r(k-1,:,:),cwtl2r(k,:,:))
        enddo
	    call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))
        else
        !ileftbisect+nbisect=nchorizo  
	    call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))   
        endif       
    ! update cwtr2l for all the beads
	    do k=ileftbisect+nbisectnow-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
        call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))  	   
    endif    
  
    return  
    end subroutine bisectv6lb  
   
   
    subroutine bisectv6rb ! move the right end. irightbisect=nchorizo
    ! move beads, then do an overall judge.
    use wavefunction
    use estimator
    use random
    use v6stepcalc
    use brut
    integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
					    ,irightbisect,k,bilvl &
	                    ,dtexpt &
                        ,nbisectnow,mmaxnow
    real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    real(kind=r8) :: xmid(3,npart),xl(3,npart)
    real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart) 
    real(kind=r8) :: tmp 
    complex(kind=r8) :: oldlv(0:mmax),newlv(0:mmax)
    complex(kind=r8) :: cwtl2rtmp1(0:nbisect,0:nspin-1,nisospin),cwtl2rtmp2(0:nbisect,0:nspin-1,nisospin) &
	                    ,cwtbegintmp(0:nspin-1,nisospin) &
                        ,cwtendtmp1now(0:nspin-1,nisospin),cwtendtmp2now(0:nspin-1,nisospin)   	              
                      
    logical :: reject	
   
    ibisecttotr=ibisecttotr+1    
   
    ileftbisect=nchorizo-nbisect
    irightbisect=nchorizo
   
    !nbisectnow=irightbisect-ileftbisect
    !mmaxnow=log(real(nbisectnow))/log(real(2))
 
    nbisectnow=nbisect
    mmaxnow=mmax  
  
    tau=dt*nbisectnow ! total time for the bisection slice.
    sigmamid=sqrt(2.*hbar*tau)  
   
    oldristra(0:nbisectnow)=ristra(ileftbisect:ileftbisect+nbisectnow)
 
    ! deal with the two ends first, then do the classic bisection.   
    newristra(0)=oldristra(0) 
    cwtbegintmp(:,:)=cwtl2r(ileftbisect,:,:)
   
    newristra(nbisectnow)%x=oldristra(0)%x+sigmamid*reshape(gaussian(3*npart),(/3,npart/))
   
    xold=oldristra(nbisectnow)%x
    xnew=newristra(nbisectnow)%x   
   
    call corop(xnew,ipro,jpro,cwtendnew)     

    oldlv(0)=1.
    newlv(0)=1. ! why 1. check again.
   
    cwtl2rtmp1(0:nbisect,:,:)=cwtl2r(ileftbisect:ileftbisect+nbisect,:,:)
    cwtl2rtmp2(0,:,:)=cwtl2rtmp1(0,:,:)
    call v6propr(oldristra(nbisectnow)%x,-1,cwtend,cwtendtmp1now)  
    call v6propr(newristra(nbisectnow)%x,-1,cwtendnew,cwtendtmp2now)    
     
    do bilvl=1,mmaxnow   ! level 1 to N. mmax = the total level N.
	    sigmamid=sqrt(hbar*tau*2.0_r8**(-bilvl))  
	    !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	    do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		    imid=2**(mmaxnow-bilvl)+i*2**(mmaxnow-bilvl+1)
		    il=imid-2**(mmaxnow-bilvl)
		    ir=imid+2**(mmaxnow-bilvl)
		    !write(6,*) 'bilvl il ir imid sigmamid=',bilvl, il,ir,imid,sigmamid
		    gauss=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		    xmid=(newristra(il)%x+newristra(ir)%x)/2.
		    newristra(imid)%x=xmid+gauss
	    enddo 
    enddo     
    ! judge
        bilvl=mmaxnow  ! besection.  level 1 to N. mmax = the total level N.	   
	    bisecttotr(bilvl)=bisecttotr(bilvl)+1.	   
	    lmax=2**bilvl-1 
	    jd=2**(mmaxnow-bilvl) ! interval
	    dtexpt=mmaxnow-bilvl-1 ! for the R L part, that's why the extra -1. dt=2**dtexp         
	    do l=1,lmax
		    j=l*jd
		    jp=(l-1)*jd
		    !xold=oldristra(j)%x
            xnew=newristra(j)%x
		    !call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(jp,:,:),cwtl2rtmp1(j,:,:))    ! should be lr, not rl I believe.
		    call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(jp,:,:),cwtl2rtmp2(j,:,:))  	   
		    !if ( sum((xnew-xold)**2).eq.0) then
			    !write(6,*) 'bilvl=',bilvl	
			    !write(6,*) 'l=',l,j,jp
			    !write(6,*) 'sigmamid=',sigmamid
		    !   write(6,*) 'xnew-xold=',xnew-xold
		    !endif   
	    enddo
	    jmax=lmax*jd  
	    oldlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp1(jmax,:,:))*cwtendtmp1now(:,:)) ))	
	    newlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp2(jmax,:,:))*cwtendtmp2now(:,:)) )) 
	    tmp=log(newlv(bilvl))-log(oldlv(bilvl))
  	    !tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
	    rn=randn(1)	   
    ! can check if l2r and r2l give	the same lv.		    
	        !write(6,*) 'tmp ln rn=',tmp,log(rn(1))
	    if (log(rn(1)).lt.tmp) then
	        reject=.false. 
            bisectcountr(bilvl)=bisectcountr(bilvl)+1.
	    else
	        reject=.true.  
	    endif 
   
    if (reject) then
        call addval(11,0.0_r8,1.0_r8)  
    else
        call addval(11,1.0_r8,1.0_r8)  
	    ibisectr=ibisectr+1 
 
        ristra(ileftbisect:ileftbisect+nbisectnow)=newristra(0:nbisectnow)     
    !ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good.	 
	 	 
    ! update cwtl2r, cwtr2l.	 
    ! update cwtl2r for all beads. this can be optimized in the future.
        cwtl2r(ileftbisect+1:ileftbisect+nbisectnow-1,:,:)=cwtl2rtmp2(1:jmax,:,:) ! here jmax should already be 2**mmax-1  
  
        call v6propl(ristra(nchorizo)%x,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:)) ! ileftbisect+nbisectnow=nchorizo  

     
        ! update cwtr2l for all the beads
        cwtend=cwtendnew         
        cwtr2l(nchorizo,:,:)=cwtendtmp2now
	    do k=nchorizo-1,1,-1
        x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
        xl=ristra(0)%x
        call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))  	

    endif    
                     
    return  
    end subroutine bisectv6rb      
   
   
    subroutine bisectv6a(ileftbisect) 
    ! slighly faster than the bisectv6improve verison. only move the beads at the currently level then judge.
    use wavefunction
    use estimator
    use random
    use v6stepcalc
    integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
					    ,k,bilvl &
	                    ,dtexpt
    real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    real(kind=r8) :: xmid(3,npart),xr(3,npart),xl(3,npart)
    real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart)
    real(kind=r8) :: tmp
    complex(kind=r8) :: oldlv(0:mmax),newlv(0:mmax)
    complex(kind=r8) :: cwtl2rtmp1(0:nbisect,0:nspin-1,nisospin),cwtr2ltmp1(0:nbisect,0:nspin-1,nisospin) &
	                    ,cwtl2rtmp2(0:nbisect,0:nspin-1,nisospin),cwtr2ltmp2(0:nbisect,0:nspin-1,nisospin) &
	                    ,cwtbegintmp(0:nspin-1,nisospin),cwtendtmp(0:nspin-1,nisospin) 
    logical :: reject	
 
    ibisecttot=ibisecttot+1  
   
    !do i=0,nbisect
	    !oldristra(i)=ristra(ileftbisect+i)
    !enddo
   
    oldristra(0:nbisect)=ristra(ileftbisect:ileftbisect+nbisect)
   
    newristra(0)=oldristra(0)
    newristra(nbisect)=oldristra(nbisect)

   
    !do i=0,nbisect
    !xoldtmp(i,:,:)=oldristra(i)%x
    !enddo
     
    cwtbegintmp(:,:)=cwtl2r(ileftbisect,:,:)
    cwtendtmp(:,:)=cwtr2l(ileftbisect+nbisect,:,:)  
   
    ! tmp1 deal with oldlv
    cwtl2rtmp1(0:nbisect,:,:)=cwtl2r(ileftbisect:ileftbisect+nbisect,:,:)
    cwtr2ltmp1(0:nbisect,:,:)=cwtr2l(ileftbisect:ileftbisect+nbisect,:,:)
    ! tmp2 deal with newlv   
    cwtl2rtmp2(0:nbisect,:,:)=cwtl2r(ileftbisect:ileftbisect+nbisect,:,:)
    cwtr2ltmp2(0:nbisect,:,:)=cwtr2l(ileftbisect:ileftbisect+nbisect,:,:)
   
    oldlv=1.
    newlv=1.   
    
    tau=dt*nbisect ! total time for the bisection slice.

    ! up to here, newristra has been loaded.
   
    ! judge if we accept the newristra or not.
    ! in this v6 case, v is not justpotential v, it should be the value of < ... >. lv.
        
    all: do bilvl=1,mmax  ! besection.  level 1 to N. mmax = the total level N.	 
        ! move the current level   
        sigmamid=sqrt(hbar*tau*2.0_r8**(-bilvl))  
	    !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	    do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		    imid=2**(mmax-bilvl)+i*2**(mmax-bilvl+1)
		    il=imid-2**(mmax-bilvl)
		    ir=imid+2**(mmax-bilvl)
		    gauss=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		    xmid=(newristra(il)%x+newristra(ir)%x)/2.
		    newristra(imid)%x=xmid+gauss
           
            !write(12,*) 'bilvl il ir imid sigmamid =',bilvl, il,ir,imid,sigmamid
            !write(12,*) 'gauss =', gauss
        enddo         
        ! then judge  
	    bisecttot(bilvl)=bisecttot(bilvl)+1.	   
	    lmax=2**bilvl-1 
	    jd=2**(mmax-bilvl) ! interval
	    dtexpt=mmax-bilvl-1 ! for the R L part, that's why the extra -1. dt=2**dtexp
	    do l=1,lmax
		    j=l*jd
		    jp=(l-1)*jd
		    xold=oldristra(j)%x
            xnew=newristra(j)%x
		    call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(jp,:,:),cwtl2rtmp1(j,:,:))    ! should be lr, not rl I believe.
    ! in the last bisection level, no need to calculate this because it old is already there. check.		   
            call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(jp,:,:),cwtl2rtmp2(j,:,:))  	   

	    enddo
	    jmax=lmax*jd
	    newlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp2(jmax,:,:))*cwtendtmp(:,:)) )) 	   
	    oldlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp1(jmax,:,:))*cwtendtmp(:,:)) ))
    
	    tmp=log(newlv(bilvl))-log(oldlv(bilvl))+log(oldlv(bilvl-1))-log(newlv(bilvl-1))   
  	    !tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
	    rn=randn(1)	   
    ! can check if l2r and r2l give	the same lv.		    
            !write(12,*) 'tmp ln rn=',tmp,log(rn(1))
	    if (log(rn(1)).lt.tmp) then
	        reject=.false. 
	        bisectcount(bilvl)=bisectcount(bilvl)+1.	 
	    else
	        reject=.true.  
	        exit all 
	    endif 
		    !write(6,*) 'imid=',imid	   
    enddo all
    if (reject) then
	    call addval(2,0.0_r8,1.0_r8)
    else
	    ibisect=ibisect+1 
	    call addval(2,1.0_r8,1.0_r8)  
 
	    !do i=0,nbisect
    !       ristra(ileftbisect+i)=newristra(i)
    !   enddo	 
        ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good. 04/10/2019£¬ remove pointer, problem fixed I tihnk.	 

        xl=ristra(0)%x     
	    xr=ristra(nchorizo)%x     
    ! update cwtl2r, cwtr2l.	 
    ! update cwtl2r for all beads. this can be optimized in the future.
     
	    cwtl2r(ileftbisect+1:ileftbisect+nbisect-1,:,:)=cwtl2rtmp2(1:jmax,:,:) ! here jmax should already be 2**mmax-1 
	    if ((ileftbisect+nbisect).le.nchorizo-1) then
        do k=ileftbisect+nbisect,nchorizo-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtl2r(k-1,:,:),cwtl2r(k,:,:))
        enddo
	    call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))
        else
        !ileftbisect+nbisect=nchorizo  
	    call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))   
        endif
   
    ! update cwtr2l for all the beads
	    do k=ileftbisect+nbisect-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
        call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))  	
  
     
    endif  

    rejectbisectv6improve=reject
    return  
    end subroutine bisectv6a         
   
   
    subroutine bisectv6improve(ileftbisect) ! for v6 interaction, nbisect beads
    use wavefunction
    use estimator
    use random
    use v6stepcalc
    integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
					    ,irightbisect,k,bilvl &
	                    ,dtexpt
    real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    real(kind=r8) :: xmid(3,npart),xr(3,npart),xl(3,npart)
    real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart)
    real(kind=r8) :: tmp 
    complex(kind=r8) :: oldlv(0:mmax),newlv(0:mmax)
    complex(kind=r8) :: cwtl2rtmp1(0:nbisect,0:nspin-1,nisospin),cwtr2ltmp1(0:nbisect,0:nspin-1,nisospin) &
	                    ,cwtl2rtmp2(0:nbisect,0:nspin-1,nisospin),cwtr2ltmp2(0:nbisect,0:nspin-1,nisospin) &
	                    ,cwtbegintmp(0:nspin-1,nisospin),cwtendtmp(0:nspin-1,nisospin)
    logical :: reject	
      
    irightbisect=ileftbisect+nbisect
    if ((nbisect<2).or.(ileftbisect < 0).or.(irightbisect>nchorizo)) stop ' bisectv6improve does not work! '
   
    ibisecttot=ibisecttot+1    
    do i=0,nbisect
	    oldristra(i)=ristra(ileftbisect+i)
    enddo
   
    newristra(0)=oldristra(0)
    newristra(nbisect)=oldristra(nbisect)
    
    !do i=0,nbisect
    !xoldtmp(i,:,:)=oldristra(i)%x
    !enddo
     
    cwtbegintmp(:,:)=cwtl2r(ileftbisect,:,:)
    cwtendtmp(:,:)=cwtr2l(ileftbisect+nbisect,:,:)  
   
    ! tmp1 deal with oldlv
    cwtl2rtmp1(0:nbisect,:,:)=cwtl2r(ileftbisect:ileftbisect+nbisect,:,:)
    cwtr2ltmp1(0:nbisect,:,:)=cwtr2l(ileftbisect:ileftbisect+nbisect,:,:)
    ! tmp2 deal with newlv   
    cwtl2rtmp2(0:nbisect,:,:)=cwtl2r(ileftbisect:ileftbisect+nbisect,:,:)
    cwtr2ltmp2(0:nbisect,:,:)=cwtr2l(ileftbisect:ileftbisect+nbisect,:,:)
   
    oldlv=1.
    newlv=1.   
    
    tau=dt*nbisect ! total time for the bisection slice.
    do bilvl=1,mmax   ! level 1 to N. mmax = the total level N.
	    sigmamid=sqrt(hbar*tau*2.0_r8**(-bilvl))  
	    !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	    do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		    imid=2**(mmax-bilvl)+i*2**(mmax-bilvl+1)
		    il=imid-2**(mmax-bilvl)
		    ir=imid+2**(mmax-bilvl)
		    !write(6,*) 'bilvl il ir imid sigmamid=',bilvl, il,ir,imid,sigmamid
		    gauss=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		    xmid=(newristra(il)%x+newristra(ir)%x)/2.
		    newristra(imid)%x=xmid+gauss
		   
		    !xold=oldristra(imid)%x
		    !xnew=newristra(imid)%x

    !       if ( sum((newristra(imid)%x-oldristra(imid)%x)**2).eq.0) then
		    !  write(6,*) '-------------------------'
		    !  write(6,*) 'bilvl=',bilvl
		    !  write(6,*) 'imid,il,ir=',imid,il,ir
		    !  write(6,*) 'sigmamid=',sigmamid
		    !  write(6,*) 'gauss=',gauss
		    !     write(6,*) 'xnew=',newristra(imid)%x
		    !  write(6,*) 'xold=',oldristra(imid)%x
    !           write(6,*) 'oldristra imid x=',oldristra(imid)%x
		    !  write(6,*) '-------------------------'
		    !
		    !  
    !do j=0,nbisect
    !write(6,*) 'j=',j  
    !write(6,*) 'xold ristra diff=', oldristra(j)%x-xoldtmp(j,:,:)
    !enddo			  
		    !  
		    !  stop
		    !  endif		   	   
		    !call hpsi(newristra(imid))
		    !call hpsi(oldristra(imid)) ! put here for save. But is not really needed if do not move beads before or after bisection
	    enddo 
    enddo  
    ! up to here, newristra has been loaded.
   
    ! judge if we accept the newristra or not.
    ! in this v6 case, v is not justpotential v, it should be the value of < ... >. lv.
		    !call hpsi0(newristra(imid))
		    !tmp=(0.5d0*(newristra(il)%v+newristra(ir)%v)-newristra(imid)%v &
        !          -0.5d0*(oldristra(il)%v+oldristra(ir)%v)+oldristra(imid)%v)*tau/(2.0_r8**bilvl) 	      

    all: do bilvl=1,mmax  ! besection.  level 1 to N. mmax = the total level N.	   
	    bisecttot(bilvl)=bisecttot(bilvl)+1.	   
	    lmax=2**bilvl-1 
	    jd=2**(mmax-bilvl) ! interval
	    dtexpt=mmax-bilvl-1 ! for the R L part, that's why the extra -1. dt=2**dtexp
	    do l=1,lmax
		    j=l*jd
		    jp=(l-1)*jd
		    xold=oldristra(j)%x
            xnew=newristra(j)%x
		    call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(jp,:,:),cwtl2rtmp1(j,:,:))    ! should be lr, not rl I believe.
		    call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(jp,:,:),cwtl2rtmp2(j,:,:))  	   
		    !if ( sum((xnew-xold)**2).eq.0) then
			    !write(6,*) 'bilvl=',bilvl	
			    !write(6,*) 'l=',l,j,jp
			    !write(6,*) 'sigmamid=',sigmamid
		    !   write(6,*) 'xnew-xold=',xnew-xold
		    !endif   
	    enddo
	    jmax=lmax*jd
	    newlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp2(jmax,:,:))*cwtendtmp(:,:)) )) 	   
	    oldlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp1(jmax,:,:))*cwtendtmp(:,:)) ))	   
	    tmp=log(newlv(bilvl))-log(oldlv(bilvl))+log(oldlv(bilvl-1))-log(newlv(bilvl-1))   
  	    !tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
	    rn=randn(1)	   
    ! can check if l2r and r2l give	the same lv.		    
	        !write(6,*) 'tmp ln rn=',tmp,log(rn(1))
	    if (log(rn(1)).lt.tmp) then
	        reject=.false. 
			 
			    !if ( sum((xnew-xold)**2).eq.0) then
			    !
			    !write(6,*) 'bilvl=',bilvl	
			    !write(6,*) 'sigmamid=',sigmamid
		    !   write(6,*) 'xnew-xold=',xnew-xold
        !         write(6,*) 'tmp ln(rn)=',tmp,log(rn(1))		   
		    !   write (6,*) 'lvnew-lvold=',newlv(bilvl)-oldlv(bilvl)
			    !write (6,*) 'previous lvnew-lvold=',newlv(bilvl-1)-oldlv(bilvl-1)		  
				    ! 
		    !
			    !endif  
				
		    !if (mod(icstepstd,5*nstepdecor).eq.0) then   
		    ! !if (tmp.gt.300.) then
		    !  !write(6,*) 'bug bead appear! exit!'
		    !     write(6,*) 'ileft bilvl il ir imid sigmamid=',ileftbisect,bilvl,il,ir,imid,sigmamid	   
	    !        write(6,*) 'tmp ln(rn)=',tmp,log(rn(1))		   
		    !     write (6,*) 'lv=',oldlv(bilvl),newlv(bilvl),oldlv(bilvl-1),newlv(bilvl-1)
		    !  write(6,*) 'cwtbegintmp+*cwtendtmp=',sum(conjg(cwtbegintmp(:,:))*cwtendtmp(:,:))
		    !  !write(6,*) 'cwtl2rtmp1=',cwtl2rtmp1(jmax,:,:),jmax
		    !  !write(6,*) 'xold=',xold	
    !           !write(6,*) 'xnew=',xnew	
		    !  !write(6,*) 'cwtl2rtmp2=',cwtl2rtmp2(jmax,:,:)	
		    !  !write(6,*) 'cwtendtmp=',cwtendtmp
		    !
		    !     !stop
		    !    endif
	    else
	        reject=.true.  
	        exit all 
	    endif 
		    !write(6,*) 'imid=',imid

	    bisectcount(bilvl)=bisectcount(bilvl)+1.
    enddo all
    if (reject) then
	    call addval(2,0.0_r8,1.0_r8)
    else
	    ibisect=ibisect+1 
	    call addval(2,1.0_r8,1.0_r8)  
 
	    do i=0,nbisect
            ristra(ileftbisect+i)=newristra(i)
        enddo	 
    !ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good.	 

        xl=ristra(0)%x     
	    xr=ristra(nchorizo)%x     
    ! update cwtl2r, cwtr2l.	 
    ! update cwtl2r for all beads. this can be optimized in the future.
     
	    cwtl2r(ileftbisect+1:ileftbisect+nbisect-1,:,:)=cwtl2rtmp2(1:jmax,:,:) ! here jmax should already be 2**mmax-1 
	    if ((ileftbisect+nbisect).le.nchorizo-1) then
        do k=ileftbisect+nbisect,nchorizo-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtl2r(k-1,:,:),cwtl2r(k,:,:))
        enddo
	    call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))
        else
        !ileftbisect+nbisect=nchorizo  
	    call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))   
        endif
   
    ! update cwtr2l for all the beads
	    do k=ileftbisect+nbisect-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
        call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))  	
      
    endif     
   
    rejectbisectv6improve=reject
    return  
    end subroutine bisectv6improve 
   
    subroutine bisectv6improvehalf(ileftbisect) ! for v6 interaction, nbisect/2 beads
    use wavefunction
    use estimator
    use random
    use v6stepcalc
    integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
					    ,irightbisect,k,bilvl &
	                    ,dtexpt
    real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    real(kind=r8) :: xmid(3,npart),xr(3,npart),xl(3,npart)
    real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart)
    real(kind=r8) :: tmp 
    complex(kind=r8) :: oldlv(0:mmaxhalf),newlv(0:mmaxhalf)
    complex(kind=r8) :: cwtl2rtmp1(0:nbisecthalf,0:nspin-1,nisospin),cwtr2ltmp1(0:nbisecthalf,0:nspin-1,nisospin) &
	                    ,cwtl2rtmp2(0:nbisecthalf,0:nspin-1,nisospin),cwtr2ltmp2(0:nbisecthalf,0:nspin-1,nisospin) &
	                    ,cwtbegintmp(0:nspin-1,nisospin),cwtendtmp(0:nspin-1,nisospin)
    logical :: reject	
   
    irightbisect=ileftbisect+nbisecthalf
    if ((nbisecthalf<2).or.(ileftbisect < 0).or.(irightbisect>nchorizo))  stop ' bisectv6improve does not work! '   
   
   
    ibisecthalftot=ibisecthalftot+1    
    do i=0,nbisecthalf
	    oldristra(i)=ristra(ileftbisect+i)
    enddo
   
    newristra(0)=oldristra(0)
    newristra(nbisecthalf)=oldristra(nbisecthalf)
    
    !do i=0,nbisect
    !xoldtmp(i,:,:)=oldristra(i)%x
    !enddo
     
    cwtbegintmp(:,:)=cwtl2r(ileftbisect,:,:)
    cwtendtmp(:,:)=cwtr2l(ileftbisect+nbisecthalf,:,:)  
   
    ! tmp1 deal with oldlv
    cwtl2rtmp1(0:nbisecthalf,:,:)=cwtl2r(ileftbisect:ileftbisect+nbisecthalf,:,:)
    cwtr2ltmp1(0:nbisecthalf,:,:)=cwtr2l(ileftbisect:ileftbisect+nbisecthalf,:,:)
    ! tmp2 deal with newlv   
    cwtl2rtmp2(0:nbisecthalf,:,:)=cwtl2r(ileftbisect:ileftbisect+nbisecthalf,:,:)
    cwtr2ltmp2(0:nbisecthalf,:,:)=cwtr2l(ileftbisect:ileftbisect+nbisecthalf,:,:)
   
    oldlv=1.
    newlv=1.   
    
    tau=dt*nbisecthalf ! total time for the bisection slice.
    do bilvl=1,mmaxhalf   ! level 1 to N. mmax = the total level N.
	    sigmamid=sqrt(hbar*tau*2.0_r8**(-bilvl))  
	    !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	    do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		    imid=2**(mmaxhalf-bilvl)+i*2**(mmaxhalf-bilvl+1)
		    il=imid-2**(mmaxhalf-bilvl)
		    ir=imid+2**(mmaxhalf-bilvl)
		    !write(6,*) 'bilvl il ir imid sigmamid=',bilvl, il,ir,imid,sigmamid
		    gauss=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		    xmid=(newristra(il)%x+newristra(ir)%x)/2.
		    newristra(imid)%x=xmid+gauss
	    enddo 
    enddo  
    ! up to here, newristra has been loaded.
   
    ! judge if we accept the newristra or not.
	      
    all: do bilvl=1,mmaxhalf  ! besection.  level 1 to N. mmax = the total level N.	   
	    bisecthalftot(bilvl)=bisecthalftot(bilvl)+1.	   
	    lmax=2**bilvl-1 
	    jd=2**(mmaxhalf-bilvl) ! interval
	    dtexpt=mmaxhalf-bilvl-1 ! for the R L part, that's why the extra -1. dt=2**dtexp
	    do l=1,lmax
		    j=l*jd
		    jp=(l-1)*jd
		    xold=oldristra(j)%x
            xnew=newristra(j)%x
		    call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(jp,:,:),cwtl2rtmp1(j,:,:))    ! should be lr, not rl I believe.
		    call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(jp,:,:),cwtl2rtmp2(j,:,:))  	   
 
	    enddo
	    jmax=lmax*jd
	    newlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp2(jmax,:,:))*cwtendtmp(:,:)) )) 	   
	    oldlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp1(jmax,:,:))*cwtendtmp(:,:)) ))	   
	    tmp=log(newlv(bilvl))-log(oldlv(bilvl))+log(oldlv(bilvl-1))-log(newlv(bilvl-1))   
  	    !tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
	    rn=randn(1)	   
    ! can check if l2r and r2l give	the same lv.		    
	        !write(6,*) 'tmp ln rn=',tmp,log(rn(1))
	    if (log(rn(1)).lt.tmp) then
	        reject=.false. 
            bisecthalfcount(bilvl)=bisecthalfcount(bilvl)+1.
	    else
	        reject=.true.  
	        exit all 
	    endif 
	   
    enddo all
    if (reject) then
	    call addval(9,0.0_r8,1.0_r8)
    else
	    ibisecthalf=ibisecthalf+1 
	    call addval(9,1.0_r8,1.0_r8)  
 
	    do i=0,nbisecthalf
            ristra(ileftbisect+i)=newristra(i)
        enddo	 
    !ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good.	 

        xl=ristra(0)%x     
	    xr=ristra(nchorizo)%x     
 
	    cwtl2r(ileftbisect+1:ileftbisect+nbisecthalf-1,:,:)=cwtl2rtmp2(1:jmax,:,:) ! here jmax should already be 2**mmax-1 
	    if ((ileftbisect+nbisecthalf).le.nchorizo-1) then
        do k=ileftbisect+nbisecthalf,nchorizo-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtl2r(k-1,:,:),cwtl2r(k,:,:))
        enddo
	    call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))
        else
        !ileftbisect+nbisect=nchorizo  
	    call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))   
        endif
    
    ! update cwtr2l for all the beads
	    do k=ileftbisect+nbisecthalf-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
        call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))  	
    
    endif     
   
    rejectbisectv6improvehalf=reject
    return  
    end subroutine bisectv6improvehalf 
   
    
    subroutine bisectv6rhalf(irightbisect) ! nbisecthalf, this is not efficient perhaps.
    use wavefunction
    use estimator
    use random
    use v6stepcalc
    use brut
    integer(kind=i4) :: ileftbisect,i &
					    ,irightbisect 
    real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    real(kind=r8) :: xnew(3,npart),xold(3,npart),xref(3,npart) &
                    ,xd2old(3,npart),xd2new(3,npart)
    real(kind=r8) :: tmp1,tmp2,logratio &
	                ,rt2,rt3
    complex(kind=r8) :: cwtold(0:nspin-1,nisospin),cwtnew(0:nspin-1,nisospin)
    complex(kind=r8) :: cwtl2rtmp(0:nspin-1,nisospin),cwtr2ltmp(0:nspin-1,nisospin)
	
    ileftbisect=irightbisect-nbisecthalf
   
    if ((nbisecthalf<2).or.(ileftbisect < 0).or.(irightbisect>nchorizo)) stop ' bisectv6rhalf does not work! '
   
    ibisecttotr=ibisecttotr+1   
   
    ristraold(0:nchorizo)=ristra(0:nchorizo)
   
    i=irightbisect
   
    xold=ristraold(i)%x 
    xref=ristraold(i-nbisecthalf)%x 
   
        tau=dt*nbisect
     
        sigmamid=sqrt(hbar*tau)  
     
        gauss=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
     
        !gauss=mov1step*reshape((randn(3*npart)-0.5),(/3,npart/)) ! this is not gauss move. for vmc(dt=0) only
     
	    xnew=xref+gauss
    	 
	    if (i.eq.nchorizo) then   
	    !call psitcwt(xold,cwtold)
      
	    cwtold=cwtend
        cwtr2ltmp(:,:)=cwtr2l(i,:,:)
         
	    tmp1=log(abs(real(sum(conjg(cwtl2r(i,:,:))*cwtold(:,:)))))
	    !call psitcwt(xnew,cwtnew)
	    call corop(xnew,ipro,jpro,cwtnew)
	    call v6propl(xnew,-1,cwtl2r(i-1,:,:),cwtl2rtmp(:,:))
	    tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtnew(:,:)))))
	    xd2new=(ristra(i-1)%x-xnew)**2
	    xd2old=(ristra(i-1)%x-xold)**2 
	    rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	  
        !rt2=0. ! vmc(dt=0) only.
	  
	    xd2new=(xref-xnew)**2
	    xd2old=(xref-xold)**2       
        rt3=1./(2.*tau*hbar)*(sum(xd2new)-sum(xd2old))
      
	    logratio=tmp2-tmp1+rt2+rt3
	    rn=randn(1)
	  
	    if (log(rn(1)).lt.logratio) then

	    ristra(i)%x=xnew
        cwtend=cwtnew
        call v6propr(xnew,-1,cwtend(:,:),cwtr2l(i,:,:))
        call bisectv6improvehalf(ileftbisect)
        if (rejectbisectv6improvehalf) then 
        ! set the values back. nothing changed.   
        cwtend=cwtold   
        ristra(i)%x=xold
        cwtr2l(i,:,:)=cwtr2ltmp(:,:) 
        else
            ibisectr=ibisectr+1   
        endif 
        endif 		     
  
        else
	   
        cwtr2ltmp(:,:)=cwtr2l(i,:,:)
        cwtold=cwtl2r(i,:,:)
	    call v6proplr(xnew,xnew,-1,cwtl2r(i-1,:,:),cwtl2rtmp(:,:))
	   
	    tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtr2l(i+1,:,:)))))
	    tmp1=log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(i+1,:,:)))))
	    
	    xd2new=(ristra(i-1)%x-xnew)**2+(ristra(i+1)%x-xnew)**2
	    xd2old=(ristra(i-1)%x-xold)**2+(ristra(i+1)%x-xold)**2
	   
	    rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	
        !rt2=0.  ! vmc(dt=0) only.
	       
        xd2new=(xref-xnew)**2
	    xd2old=(xref-xold)**2       
        rt3=1./(2.*tau*hbar)*(sum(xd2new)-sum(xd2old))
       
	    logratio=tmp2-tmp1+rt2+rt3	 
	    rn=randn(1)

	    if (log(rn(1)).lt.logratio) then
	    ristra(i)%x=xnew 
        call v6proplr(xnew,xnew,-1,cwtr2l(i+1,:,:),cwtr2l(i,:,:))
        call bisectv6improvehalf(ileftbisect)
        if (rejectbisectv6improvehalf) then 
        ! set the values back. nothing changed.   
        cwtr2l(i,:,:)=cwtr2ltmp(:,:)   
        ristra(i)%x=xold
        else
            ibisectr=ibisectr+1   
        endif 
        endif
         
	    endif 
   
    return  
    end subroutine bisectv6rhalf    
   

    subroutine bisectv6lhalf(ileftbisect) ! move the left end, 0th bead,ileftbiect=0
    use wavefunction
    use estimator
    use random
    use v6stepcalc
    use brut
    integer(kind=i4) :: ileftbisect,i &
					    ,irightbisect 
    real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    real(kind=r8) :: xnew(3,npart),xold(3,npart),xref(3,npart) &
                    ,xd2old(3,npart),xd2new(3,npart)
    real(kind=r8) :: tmp1,tmp2,logratio &
	                ,rt2,rt3
    complex(kind=r8) :: cwtold(0:nspin-1,nisospin),cwtnew(0:nspin-1,nisospin)
    complex(kind=r8) :: cwtl2rtmp(0:nspin-1,nisospin),cwtr2ltmp(0:nspin-1,nisospin) 	              
   
    irightbisect=ileftbisect+nbisecthalf
    if ((nbisecthalf<2).or.(ileftbisect < 0).or.(irightbisect>nchorizo)) stop ' bisectv6lhalf does not work! '
   
    ibisecttotl=ibisecttotl+1   
   
    ristraold(0:nchorizo)=ristra(0:nchorizo)

    i=ileftbisect
   
    xold=ristraold(i)%x 
    xref=ristraold(i+nbisecthalf)%x 
   
        tau=dt*nbisect
     
        sigmamid=sqrt(hbar*tau)  
     
        gauss=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
     
        !gauss=mov1step*reshape((randn(3*npart)-0.5),(/3,npart/)) ! this is not gauss move. for vmc(dt=0) only
     
	    xnew=xref+gauss
    	 
	    if (i.eq.0) then   
    
	    cwtold=cwtbegin
        cwtl2rtmp(:,:)=cwtl2r(i,:,:)
        
	    tmp1=log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(i,:,:)))))
	    !call psitcwt(xnew,cwtnew)  
	    call corop(xnew,iplo,jplo,cwtnew)
	    call v6propl(xnew,-1,cwtr2l(i+1,:,:),cwtr2ltmp(:,:))
	    tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
	    xd2new=(ristra(i+1)%x-xnew)**2
	    xd2old=(ristra(i+1)%x-xold)**2
	  
	    rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	  
        !rt2=0. ! vmc(dt=0) only.
	  
	    xd2new=(xref-xnew)**2
	    xd2old=(xref-xold)**2       
        rt3=1./(2.*tau*hbar)*(sum(xd2new)-sum(xd2old))      
      
	    logratio=tmp2-tmp1+rt2+rt3
	    rn=randn(1)
  
	    if (log(rn(1)).lt.logratio) then ! dt ne 0
	    ristra(i)%x=xnew
    ! update cwtbegin
        cwtbegin=cwtnew	   
    ! update cwtl2r	   
	    call v6propr(xnew,-1,cwtnew(:,:),cwtl2r(i,:,:)) ! just update the i. 
    
        call bisectv6improvehalf(ileftbisect)
        if (rejectbisectv6improvehalf) then 
        ! set the values back. nothing changed.   
        cwtbegin=cwtold   
        ristra(i)%x=xold
        cwtl2r(i,:,:)=cwtl2rtmp(:,:) 
        else
            ibisectl=ibisectl+1   
        endif   
	    endif       
      
        else
	        
        cwtold=cwtl2r(i,:,:)
	    call v6proplr(xnew,xnew,-1,cwtl2r(i-1,:,:),cwtl2rtmp(:,:))
	   
	    tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtr2l(i+1,:,:)))))
	    tmp1=log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(i+1,:,:)))))
	    
	    xd2new=(ristra(i-1)%x-xnew)**2+(ristra(i+1)%x-xnew)**2
	    xd2old=(ristra(i-1)%x-xold)**2+(ristra(i+1)%x-xold)**2
	   
	    rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	
        !rt2=0.  ! vmc(dt=0) only.
	       
        xd2new=(xref-xnew)**2
	    xd2old=(xref-xold)**2       
        rt3=1./(2.*tau*hbar)*(sum(xd2new)-sum(xd2old))
       
	    logratio=tmp2-tmp1+rt2+rt3	 
	    rn=randn(1)

	    if (log(rn(1)).lt.logratio) then
	    ristra(i)%x=xnew 
        call v6proplr(xnew,xnew,-1,cwtl2r(i-1,:,:),cwtl2r(i,:,:))
        call bisectv6improvehalf(ileftbisect)
        if (rejectbisectv6improvehalf) then 
        ! set the values back. nothing changed.   
        cwtl2r(i,:,:)=cwtold
        ristra(i)%x=xold
        else
            ibisectl=ibisectl+1   
        endif 
        endif
         
	    endif 
 
    return  
    end subroutine bisectv6lhalf         
   
    subroutine bisectv6la ! should be better than bisectv6l,rhalf     
    ! move the left end, 0th bead,ileftbiect=0; if ileftbisect /= 0 then the small modification on 0th level is needed.
    use wavefunction
    use estimator
    use random
    use v6stepcalc
    use brut
    integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
					    ,irightbisect,k,bilvl &
	                    ,dtexpt &
                        ,nbisectnow,mmaxnow
    real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    real(kind=r8) :: xmid(3,npart),xr(3,npart),xl(3,npart)
    real(kind=r8) :: x(3,npart),xnew(3,npart),x0new(3,npart),xold(3,npart) 
    real(kind=r8) :: tmp 
    complex(kind=r8) :: oldlv(0:mmaxhalf),newlv(0:mmaxhalf)
    complex(kind=r8) :: cwtl2rtmp1(0:nbisecthalf,0:nspin-1,nisospin) &
	                    ,cwtl2rtmp2(0:nbisecthalf,0:nspin-1,nisospin) &
	                    ,cwtendtmp(0:nspin-1,nisospin) &
                        ,cwtbegintmp2(0:nspin-1,nisospin)        
    logical :: reject	

    ibisecttotl=ibisecttotl+1  
   
    ileftbisect=0
    irightbisect=ileftbisect+nbisecthalf
   
    !nbisectnow=irightbisect-ileftbisect ! has to be 2^n.
    !mmaxnow=log(real(nbisectnow))/log(real(2)) ! mmaxnow = mmax-1
   
    nbisectnow=nbisecthalf
    mmaxnow=mmaxhalf
   
    tau=dt*nbisectnow ! total time for the bisection slice. nisectnow is usually nbisecthalf
    sigmamid=sqrt(2.*hbar*tau)  
   
    do i=0,nbisectnow
	    oldristra(i)=ristra(ileftbisect+i)
    enddo   
 
    ! deal with the two ends first, then do the classic bisection.   
   
    newristra(nbisectnow)=oldristra(nbisectnow)   
    cwtendtmp(:,:)=cwtr2l(ileftbisect+nbisectnow,:,:)
 
    newristra(0)%x=oldristra(nbisectnow)%x+sigmamid*reshape(gaussian(3*npart),(/3,npart/))
   
    x0new=newristra(0)%x
   
    ! i,jpro,i,jplo, the o means old are actually the current order.
   
    call corop(x0new,iplo,jplo,cwtbegintmp2)     
   
    oldlv(0)=1.
    newlv(0)=1.

    cwtl2rtmp1(0,:,:)=cwtl2r(ileftbisect,:,:)
    call v6propr(x0new,-1,cwtbegintmp2,cwtl2rtmp2(0,:,:)) 
   
    ! prepare the further new bisection positions
    do bilvl=1,mmaxnow   ! level 1 to N. mmax = the total level N.
	    sigmamid=sqrt(hbar*tau*2.0_r8**(-bilvl))  
	    !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	    do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		    imid=2**(mmaxnow-bilvl)+i*2**(mmaxnow-bilvl+1)
		    il=imid-2**(mmaxnow-bilvl)
		    ir=imid+2**(mmaxnow-bilvl)
		    !write(6,*) 'bilvl il ir imid sigmamid=',bilvl, il,ir,imid,sigmamid
		    gauss=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		    xmid=(newristra(il)%x+newristra(ir)%x)/2.
		    newristra(imid)%x=xmid+gauss
	    enddo 
    enddo      
    ! further bisection    
    all: do bilvl=1,mmaxnow  ! besection.  level 1 to N. mmax = the total level N.	   
	    bisecttotl(bilvl)=bisecttotl(bilvl)+1.	   
	    lmax=2**bilvl-1 
	    jd=2**(mmaxnow-bilvl) ! interval
	    dtexpt=mmaxnow-bilvl-1 ! for the R L part, that's why the extra -1. dt=2**dtexp 
	    do l=1,lmax
		    j=l*jd
		    jp=(l-1)*jd
		    xold=oldristra(j)%x
            xnew=newristra(j)%x
		    call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(jp,:,:),cwtl2rtmp1(j,:,:))    ! should be lr, not rl I believe.
		    call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(jp,:,:),cwtl2rtmp2(j,:,:))  	   
		    !if ( sum((xnew-xold)**2).eq.0) then
			    !write(6,*) 'bilvl=',bilvl	
			    !write(6,*) 'l=',l,j,jp
			    !write(6,*) 'sigmamid=',sigmamid
		    !   write(6,*) 'xnew-xold=',xnew-xold
		    !endif   
	    enddo
	    jmax=lmax*jd  
	    oldlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp1(jmax,:,:))*cwtendtmp(:,:)) ))	
	    newlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp2(jmax,:,:))*cwtendtmp(:,:)) )) 
	    tmp=log(newlv(bilvl))-log(oldlv(bilvl))+log(oldlv(bilvl-1))-log(newlv(bilvl-1))   
  	    !tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
	    rn=randn(1)	   
    ! can check if l2r and r2l give	the same lv.		    
	        !write(6,*) 'tmp ln rn=',tmp,log(rn(1))
	    if (log(rn(1)).lt.tmp) then
	        reject=.false. 
            bisectcountl(bilvl)=bisectcountl(bilvl)+1.
	    else
	        reject=.true.  
	        exit all 
	    endif   
    enddo all   
    
    if (reject) then
    else
	    ibisectl=ibisectl+1 
	    do i=0,nbisectnow
            ristra(ileftbisect+i)=newristra(i)
        enddo	 
    !ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good.	 
        xl=ristra(0)%x     
	    xr=ristra(nchorizo)%x 
    ! update cwtl2r, cwtr2l.	 
        cwtbegin=cwtbegintmp2
	    cwtl2r(ileftbisect:ileftbisect+nbisectnow-1,:,:)=cwtl2rtmp2(0:jmax,:,:) ! here jmax should already be 2**mmax-1    
	    if ((ileftbisect+nbisectnow).le.nchorizo-1) then
        do k=ileftbisect+nbisectnow,nchorizo-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtl2r(k-1,:,:),cwtl2r(k,:,:))
        enddo
	    call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))
        else
        !ileftbisect+nbisect=nchorizo  
	    call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))   
        endif       
    ! update cwtr2l for all the beads
	    do k=ileftbisect+nbisectnow-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
        call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))  	   
    endif    
  
    return  
    end subroutine bisectv6la  
   
   
    subroutine bisectv6ra ! move the right end. irightbisect=nchorizo
    use wavefunction
    use estimator
    use random
    use v6stepcalc
    use brut
    integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
					    ,irightbisect,k,bilvl &
	                    ,dtexpt &
                        ,nbisectnow,mmaxnow
    real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    real(kind=r8) :: xmid(3,npart),xl(3,npart)
    real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart)
    real(kind=r8) :: tmp 
    complex(kind=r8) :: oldlv(0:mmaxhalf),newlv(0:mmaxhalf)
    complex(kind=r8) :: cwtl2rtmp1(0:nbisecthalf,0:nspin-1,nisospin),cwtl2rtmp2(0:nbisecthalf,0:nspin-1,nisospin) &
	                    ,cwtbegintmp(0:nspin-1,nisospin) &
                        ,cwtendtmp1now(0:nspin-1,nisospin),cwtendtmp2now(0:nspin-1,nisospin)   	              
                      
    logical :: reject	
   
    ibisecttotr=ibisecttotr+1    
   
    ileftbisect=nchorizo-nbisecthalf
    irightbisect=nchorizo
   
    !nbisectnow=irightbisect-ileftbisect
    !mmaxnow=log(real(nbisectnow))/log(real(2))
 
    nbisectnow=nbisecthalf
    mmaxnow=mmaxhalf   
  
    tau=dt*nbisectnow ! total time for the bisection slice.
    sigmamid=sqrt(2.*hbar*tau)  
   
    do i=0,nbisectnow
	    oldristra(i)=ristra(ileftbisect+i)
    enddo   
 
    ! deal with the two ends first, then do the classic bisection.   
    newristra(0)=oldristra(0) 
    cwtbegintmp(:,:)=cwtl2r(ileftbisect,:,:)
   
    newristra(nbisectnow)%x=oldristra(0)%x+sigmamid*reshape(gaussian(3*npart),(/3,npart/))
   
    xold=oldristra(nbisectnow)%x
    xnew=newristra(nbisectnow)%x   
   
    call corop(xnew,ipro,jpro,cwtendnew)     

    oldlv(0)=1.
    newlv(0)=1.
   
    cwtl2rtmp1(0,:,:)=cwtl2r(ileftbisect,:,:)
    cwtl2rtmp2(0,:,:)=cwtl2rtmp1(0,:,:)
    call v6propr(oldristra(nbisectnow)%x,-1,cwtend,cwtendtmp1now)  
    call v6propr(newristra(nbisectnow)%x,-1,cwtendnew,cwtendtmp2now)    
     
    do bilvl=1,mmaxnow   ! level 1 to N. mmax = the total level N.
	    sigmamid=sqrt(hbar*tau*2.0_r8**(-bilvl))  
	    !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	    do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		    imid=2**(mmaxnow-bilvl)+i*2**(mmaxnow-bilvl+1)
		    il=imid-2**(mmaxnow-bilvl)
		    ir=imid+2**(mmaxnow-bilvl)
		    !write(6,*) 'bilvl il ir imid sigmamid=',bilvl, il,ir,imid,sigmamid
		    gauss=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		    xmid=(newristra(il)%x+newristra(ir)%x)/2.
		    newristra(imid)%x=xmid+gauss
	    enddo 
    enddo     
    ! further bisection
    all: do bilvl=1,mmaxnow  ! besection.  level 1 to N. mmax = the total level N.	   
	    bisecttotr(bilvl)=bisecttotr(bilvl)+1.	   
	    lmax=2**bilvl-1 
	    jd=2**(mmaxnow-bilvl) ! interval
	    dtexpt=mmaxnow-bilvl-1 ! for the R L part, that's why the extra -1. dt=2**dtexp         
	    do l=1,lmax
		    j=l*jd
		    jp=(l-1)*jd
		    xold=oldristra(j)%x
            xnew=newristra(j)%x
		    call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(jp,:,:),cwtl2rtmp1(j,:,:))    ! should be lr, not rl I believe.
		    call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(jp,:,:),cwtl2rtmp2(j,:,:))  	   
		    !if ( sum((xnew-xold)**2).eq.0) then
			    !write(6,*) 'bilvl=',bilvl	
			    !write(6,*) 'l=',l,j,jp
			    !write(6,*) 'sigmamid=',sigmamid
		    !   write(6,*) 'xnew-xold=',xnew-xold
		    !endif   
	    enddo
	    jmax=lmax*jd  
	    oldlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp1(jmax,:,:))*cwtendtmp1now(:,:)) ))	
	    newlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp2(jmax,:,:))*cwtendtmp2now(:,:)) )) 
	    tmp=log(newlv(bilvl))-log(oldlv(bilvl))+log(oldlv(bilvl-1))-log(newlv(bilvl-1))   
  	    !tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
	    rn=randn(1)	   
    ! can check if l2r and r2l give	the same lv.		    
	        !write(6,*) 'tmp ln rn=',tmp,log(rn(1))
	    if (log(rn(1)).lt.tmp) then
	        reject=.false. 
            bisectcountr(bilvl)=bisectcountr(bilvl)+1.
	    else
	        reject=.true.  
	        exit all 
	    endif 
    enddo all   
    
    if (reject) then
    else
	    ibisectr=ibisectr+1 
 
	    do i=0,nbisectnow
            ristra(ileftbisect+i)=newristra(i)
        enddo	 

    !ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good.	 
	 	 
    ! update cwtl2r, cwtr2l.	 
    ! update cwtl2r for all beads. this can be optimized in the future.
        cwtl2r(ileftbisect+1:ileftbisect+nbisectnow-1,:,:)=cwtl2rtmp2(1:jmax,:,:) ! here jmax should already be 2**mmax-1  
  
        call v6propl(ristra(nchorizo)%x,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:)) ! ileftbisect+nbisectnow=nchorizo  

     
        ! update cwtr2l for all the beads
        cwtend=cwtendnew         
        cwtr2l(nchorizo,:,:)=cwtendtmp2now
	    do k=nchorizo-1,1,-1
        x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
        xl=ristra(0)%x
        call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))  	

    endif    
                     
    return  
    end subroutine bisectv6ra   
   
    !   this is not the best bisectv6l and r, but is ok.
    subroutine bisectv6r ! move the right end. irightbisect=nchorizo
    ibisecttotr=ibisecttotr+1 
    call bisectv6improve(nchorizo-nbisect)
    if (.not.rejectbisectv6improve) ibisectr=ibisectr+1 
    return  
    end subroutine bisectv6r    
   
    subroutine bisectv6l ! move the left end, 0th bead,ileftbiect=0
    ibisecttotl=ibisecttotl+1  
    call bisectv6improve(0)
    if (.not.rejectbisectv6improve) ibisectl=ibisectl+1 
    return  
    end subroutine bisectv6l    


    subroutine bisectv6improvew11 ! for v6 interaction with move 1 by 1.
    use wavefunction
    use estimator
    use random
    use v6stepcalc
    integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
					    ,irightbisect,k,bilvl,n1,n2 &
	                    ,dtexpt
    real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    real(kind=r8) :: xmid(3,npart),xr(3,npart),xl(3,npart)
    real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart)
    real(kind=r8) :: tmp 
    complex(kind=r8) :: oldlv(0:mmax),newlv(0:mmax)
    complex(kind=r8) :: cwtl2rtmp1(0:nbisect,0:nspin-1,nisospin),cwtr2ltmp1(0:nbisect,0:nspin-1,nisospin) &
	                    ,cwtl2rtmp2(0:nbisect,0:nspin-1,nisospin),cwtr2ltmp2(0:nbisect,0:nspin-1,nisospin) &
	                    ,cwtbegintmp(0:nspin-1,nisospin),cwtendtmp(0:nspin-1,nisospin)
    logical :: reject	
    icbisecttot=icbisecttot+1
    if (nbisect.lt.nchorizo) then
    call bisectpicksliceswide(ileftbisect,irightbisect) ! pick the ileftbisect here
    
    !write(6,*) 'ileftbisect=', ileftbisect
    !write(12,*) 'ileftbisect=', ileftbisect  
    
    if (ileftbisect.le.(nchorizo-nbisect)) then
    ibisecttot=ibisecttot+1 
   
    !write(6,*) 'bisection is running'
   
    !oldristra(0:nbisect)=ristra(ileftbisect:ileftbisect+nbisect)
   
    do i=0,nbisect
	    oldristra(i)=ristra(ileftbisect+i)
    enddo
   
    newristra(0)=oldristra(0)
    newristra(nbisect)=oldristra(nbisect)
    
    !do i=0,nbisect
    !xoldtmp(i,:,:)=oldristra(i)%x
    !enddo
     
    cwtbegintmp(:,:)=cwtl2r(ileftbisect,:,:)
    cwtendtmp(:,:)=cwtr2l(ileftbisect+nbisect,:,:)
   
    ! tmp1 deal with oldlv
    cwtl2rtmp1(0:nbisect,:,:)=cwtl2r(ileftbisect:ileftbisect+nbisect,:,:)
    cwtr2ltmp1(0:nbisect,:,:)=cwtr2l(ileftbisect:ileftbisect+nbisect,:,:)
    ! tmp2 deal with newlv   
    cwtl2rtmp2(0:nbisect,:,:)=cwtl2r(ileftbisect:ileftbisect+nbisect,:,:)
    cwtr2ltmp2(0:nbisect,:,:)=cwtr2l(ileftbisect:ileftbisect+nbisect,:,:)
   
    oldlv=1.
    newlv=1.
    
    tau=dt*nbisect ! total time for the bisection slice.
    do bilvl=1,mmax   ! level 1 to N. mmax = the total level N.
	    sigmamid=sqrt(hbar*tau*2.0d0**(-bilvl))  
	    !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	    do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		    imid=2**(mmax-bilvl)+i*2**(mmax-bilvl+1)
		    il=imid-2**(mmax-bilvl)
		    ir=imid+2**(mmax-bilvl)
		    !write(6,*) 'bilvl il ir imid sigmamid=',bilvl, il,ir,imid,sigmamid
		    gauss=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		    xmid=(newristra(il)%x+newristra(ir)%x)/2.
		    newristra(imid)%x=xmid+gauss
		   
		    !xold=oldristra(imid)%x
		    !xnew=newristra(imid)%x

    !       if ( sum((newristra(imid)%x-oldristra(imid)%x)**2).eq.0) then
		    !  write(6,*) '-------------------------'
		    !  write(6,*) 'bilvl=',bilvl
		    !  write(6,*) 'imid,il,ir=',imid,il,ir
		    !  write(6,*) 'sigmamid=',sigmamid
		    !  write(6,*) 'gauss=',gauss
		    !     write(6,*) 'xnew=',newristra(imid)%x
		    !  write(6,*) 'xold=',oldristra(imid)%x
    !           write(6,*) 'oldristra imid x=',oldristra(imid)%x
		    !  write(6,*) '-------------------------'
		    !
		    !  
    !do j=0,nbisect
    !write(6,*) 'j=',j  
    !write(6,*) 'xold ristra diff=', oldristra(j)%x-xoldtmp(j,:,:)
    !enddo			  
		    !  
		    !  stop
		    !  endif		   	   
		    !call hpsi(newristra(imid))
		    !call hpsi(oldristra(imid)) ! put here for save. But is not really needed if do not move beads before or after bisection
	    enddo 
    enddo  
    ! up to here, newristra has been loaded.
   
    ! judge if we accept the newristra or not.
    ! in this v6 case, v is not justpotential v, it should be the value of < ... >. lv.
		    !call hpsi0(newristra(imid))
		    !tmp=(0.5d0*(newristra(il)%v+newristra(ir)%v)-newristra(imid)%v &
        !          -0.5d0*(oldristra(il)%v+oldristra(ir)%v)+oldristra(imid)%v)*tau/(2.0d0**bilvl) 	      

    all: do bilvl=1,mmax  ! besection.  level 1 to N. mmax = the total level N.	   
	    bisecttot(bilvl)=bisecttot(bilvl)+1.	   
	    lmax=2**bilvl-1 
	    jd=2**(mmax-bilvl) ! interval
	    dtexpt=mmax-bilvl-1 ! for the R L part, that's why the extra -1. dt=2**dtexp
	    do l=1,lmax
		    j=l*jd
		    jp=(l-1)*jd
		    xold=oldristra(j)%x
            xnew=newristra(j)%x
		    call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(jp,:,:),cwtl2rtmp1(j,:,:))    ! should be lr, not rl I believe.
		    call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(jp,:,:),cwtl2rtmp2(j,:,:))  	   
		    !if ( sum((xnew-xold)**2).eq.0) then
			    !write(6,*) 'bilvl=',bilvl	
			    !write(6,*) 'l=',l,j,jp
			    !write(6,*) 'sigmamid=',sigmamid
		    !   write(6,*) 'xnew-xold=',xnew-xold
		    !endif   
	    enddo
	    jmax=lmax*jd
	    newlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp2(jmax,:,:))*cwtendtmp(:,:)) )) 	   
	    oldlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp1(jmax,:,:))*cwtendtmp(:,:)) ))	   
	    tmp=log(newlv(bilvl))-log(oldlv(bilvl))+log(oldlv(bilvl-1))-log(newlv(bilvl-1))   
  	    !tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
	    rn=randn(1)	   
    ! can check if l2r and r2l give	the same lv.		    
	        !write(6,*) 'tmp ln rn=',tmp,log(rn(1))
	    if (log(rn(1)).lt.tmp) then
	        reject=.false. 
			 
			    !if ( sum((xnew-xold)**2).eq.0) then
			    !
			    !write(6,*) 'bilvl=',bilvl	
			    !write(6,*) 'sigmamid=',sigmamid
		    !   write(6,*) 'xnew-xold=',xnew-xold
        !         write(6,*) 'tmp ln(rn)=',tmp,log(rn(1))		   
		    !   write (6,*) 'lvnew-lvold=',newlv(bilvl)-oldlv(bilvl)
			    !write (6,*) 'previous lvnew-lvold=',newlv(bilvl-1)-oldlv(bilvl-1)		  
				    ! 
		    !
			    !endif  
				
		    !if (mod(icstepstd,5*nstepdecor).eq.0) then   
		    ! !if (tmp.gt.300.) then
		    !  !write(6,*) 'bug bead appear! exit!'
		    !     write(6,*) 'ileft bilvl il ir imid sigmamid=',ileftbisect,bilvl,il,ir,imid,sigmamid	   
	    !        write(6,*) 'tmp ln(rn)=',tmp,log(rn(1))		   
		    !     write (6,*) 'lv=',oldlv(bilvl),newlv(bilvl),oldlv(bilvl-1),newlv(bilvl-1)
		    !  write(6,*) 'cwtbegintmp+*cwtendtmp=',sum(conjg(cwtbegintmp(:,:))*cwtendtmp(:,:))
		    !  !write(6,*) 'cwtl2rtmp1=',cwtl2rtmp1(jmax,:,:),jmax
		    !  !write(6,*) 'xold=',xold	
    !           !write(6,*) 'xnew=',xnew	
		    !  !write(6,*) 'cwtl2rtmp2=',cwtl2rtmp2(jmax,:,:)	
		    !  !write(6,*) 'cwtendtmp=',cwtendtmp
		    !
		    !     !stop
		    !    endif
	    else
	        reject=.true.  
	        exit all 
	    endif 
		    !write(6,*) 'imid=',imid

	    bisectcount(bilvl)=bisectcount(bilvl)+1.
    enddo all
    if (reject) then
	    call addval(2,0.0_r8,1.0_r8)
    else
	    ibisect=ibisect+1 
	    call addval(2,1.0_r8,1.0_r8)  
 
	    do i=0,nbisect
            ristra(ileftbisect+i)=newristra(i)
        enddo	 
    !ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good.	 
	 	 
    ! update cwtl2r, cwtr2l.	 
    ! update cwtl2r for all beads. this can be optimized in the future.
	 
        xr=ristra(nchorizo)%x
        xl=ristra(0)%x
     
        cwtl2r(ileftbisect+1:ileftbisect+nbisect-1,:,:)=cwtl2rtmp2(1:jmax,:,:) ! here jmax should already be 2**mmax-1 
        if ((ileftbisect+nbisect).le.nchorizo-1) then
        do k=ileftbisect+nbisect,nchorizo-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtl2r(k-1,:,:),cwtl2r(k,:,:))
	    enddo
	    call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))
        else
        call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))    
        endif
     
     
    ! update cwtr2l for all the beads
     
	    do k=ileftbisect+nbisect-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
        call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))  	
    endif  
    else
   
            !write(6,*) 'bisection is not running'
       
	    ! this part is somewhat like reptate perhaps.
    ! need to modify here   
	
       
    !!!--------------------    
       
	    if (ileftbisect.lt.nchorizo) then ! nbisect should be >= 2.
	    n1=nchorizo-ileftbisect
	    n2=nbisect-1-n1
	    !write (6,*) 'bisect improve', ileftbisect+1,nchorizo,n1,n2,n2-1
	    !call stdmove1range(ileftbisect+1,nchorizo)
        
        !call stdmove1rangefastv6(ileftbisect+1,nchorizo)	
        call stdmove1fastv6(ileftbisect+1,nchorizo)
        
	    if (n2.gt.0) then
	    !call stdmove1range(0,n2-1)
            
	    !call stdmove1rangefastv6(0,n2-1)
        call stdmove1fastv6(0,n2-1)
        
        endif
        else ! ileftbisect.eq.nchorizo
		    !write (6,*) 'bisect improve ileftbisect=', ileftbisect+1,nchorizo,n1,n2,n2-1
	    !call stdmove1range(0,nbisect-2)	
         
        !call stdmove1rangefastv6(0,nbisect-2)	
        call stdmove1fastv6(0,nbisect-2)
        endif	
       
    !!!--------------------     
     
    endif
    else
    !call bisectnchorizo  
	    write (6,*) 'nbisect >= nchorizo, stop!'
	    stop 
    endif
  
  
    ! do not need to calculate elocal too often.  
    !if (mod(icbisecttot,nstepdecor).eq.0) then
    !ice=ice+1
    !call hpsi(ristra(0))
    !eleft=ristra(0)%elocal  
    !call hpsi(ristra(nchorizo))
    !eright=ristra(nchorizo)%elocal 
    !select case (irep)
    !case (1:2)
    !   call addval(3,eleft,1.0d0)
    !   call addval(4,eright,1.0d0)
	    !call addval(5,(eleft+eright)/2.,1.0d0)
	    !call funcrmsq(ristra(nchorizomid)%x,rm,rmsq)
	    !call addval(6,rmsq,1.0d0)
	    !call addval(7,rm,1.0d0)
	    !write(11,'(2x,3(G12.5,2x))') eleft,eright,(eleft+eright)/2.,rm!11 means write to 'values.txt'
    !case (3:4)
    !   call addval(4,eleft,1.0d0)
    !   call addval(5,eright,1.0d0)
	    !write(11,'(2x,3(G12.5,2x))') eleft,eright,(eleft+eright)/2.!11 means write to 'values.txt'
    !end select 
    !endif
    return  
    end subroutine bisectv6improvew11      

   
    subroutine stdmove1fastv6a(ileft,iright) ! need to modify.
    ! ileft < iright; update from left to right. 
    ! This is move 1 by 1, but different way of move, this is like a 2-beads bisect.
    use wavefunction
    use estimator
    use v6stepcalc
    use random
    use brut
    integer(kind=i4) :: i,k,ileft,iright,icnow
    real(kind=r8) :: rn(1),gauss(3,npart)
    real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart),xoldnxt(3,npart),xd2old(3,npart),xd2new(3,npart) &
	                ,xr(3,npart),xl(3,npart)
    real(kind=r8) :: rt2,tmp1,tmp2,logratio
    complex(kind=r8) :: cwtold(0:nspin-1,nisospin),cwtnew(0:nspin-1,nisospin)
    complex(kind=r8) :: cwtl2rtmp(0:nspin-1,nisospin),cwtr2ltmp(0:nspin-1,nisospin)

    if ((ileft.lt.0).or.(iright.gt.nchorizo).or.(iright.lt.0).or.(ileft.gt.nchorizo).or.(ileft.gt.iright)) then
	    write (6,*) 'stdmove1fastv6 range error', ileft,iright
	    stop
    endif
       
    ristraold(0:nchorizo)=ristra(0:nchorizo)
    icstdmove1tot=icstdmove1tot+1
    icnow=0
   
    do i=ileft,iright
      
        !write(6,*) 'i=',i
       
	    icmov1tot=icmov1tot+1  
	    xold=ristraold(i)%x 
    ! ran move	 
	    !rn3=reshape(randn(3*npart),shape(rn3)) 
	    !xnew=xold+mov1step*(rn3-0.5_r8) 
    ! gaussian move	 
	 
        gauss=mov1step*reshape(gaussian(3*npart),(/3,npart/))
     
        !gauss=mov1step*reshape((randn(3*npart)-0.5),(/3,npart/)) ! this is not gauss move. for vmc(dt=0) only
     
	    xnew=xold+gauss
	 
        if (i.eq.0) then
		 
	    !call psitcwt(xold,cwtold)
	  
	    cwtold=cwtbegin
	  
	    tmp1=log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(i,:,:)))))
	    !call psitcwt(xnew,cwtnew)  
	    call corop(xnew,iplo,jplo,cwtnew)
	    call v6propl(xnew,-1,cwtr2l(i+1,:,:),cwtr2ltmp(:,:))
	    tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
	    xd2new=(ristra(i+1)%x-xnew)**2
	    xd2old=(ristra(i+1)%x-xold)**2
	  
	    rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	  
        !rt2=0. ! vmc(dt=0) only.
	  
	    logratio=tmp2-tmp1+rt2
	    rn=randn(1)
  
	    if (log(rn(1)).lt.logratio) then ! dt ne 0
        icnow=icnow+1 
	    icmov1=icmov1+1	
	    ristra(i)%x=xnew
	    cwtr2l(i,:,:)=cwtr2ltmp(:,:)   
    ! update cwtbegin
        cwtbegin=cwtnew	   
    ! update cwtl2r	   
	    call v6propr(xnew,-1,cwtnew(:,:),cwtl2r(i,:,:)) ! just update the i.
    ! prepare for the next cwtl2r
        xoldnxt=ristra(i+1)%x
        call v6proplr(xoldnxt,xoldnxt,-1,cwtl2r(i,:,:),cwtl2r(i+1,:,:))            
	    endif 
	  
	    else
		 
	    if (i.eq.nchorizo) then   
	    !call psitcwt(xold,cwtold)
	    cwtold=cwtend
	    tmp1=log(abs(real(sum(conjg(cwtl2r(i,:,:))*cwtold(:,:)))))
	    !call psitcwt(xnew,cwtnew)
	    call corop(xnew,ipro,jpro,cwtnew)
	    call v6propl(xnew,-1,cwtl2r(i-1,:,:),cwtl2rtmp(:,:))
	    tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtnew(:,:)))))
	    xd2new=(ristra(i-1)%x-xnew)**2
	    xd2old=(ristra(i-1)%x-xold)**2
	  
	    rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	  
        !rt2=0. ! vmc(dt=0) only.
	  
	    logratio=tmp2-tmp1+rt2
	    rn=randn(1)
	  
	    if (log(rn(1)).lt.logratio) then
        icnow=icnow+1
	    icmov1=icmov1+1	
	    ristra(i)%x=xnew
	    cwtl2r(i,:,:)=cwtl2rtmp(:,:)
    ! update cwtend
        cwtend=cwtnew
	    endif 		     
	   
	    else
	   
        cwtold=cwtl2r(i,:,:)
	    call v6proplr(xnew,xnew,-1,cwtl2r(i-1,:,:),cwtl2rtmp(:,:))
	   
	    tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtr2l(i+1,:,:)))))
	    tmp1=log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(i+1,:,:)))))
	    
	    xd2new=(ristra(i-1)%x-xnew)**2+(ristra(i+1)%x-xnew)**2
	    xd2old=(ristra(i-1)%x-xold)**2+(ristra(i+1)%x-xold)**2
	   
	    rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	
        !rt2=0.  ! vmc(dt=0) only.
	   
	    logratio=tmp2-tmp1+rt2	 
	    rn=randn(1)

	    if (log(rn(1)).lt.logratio) then
            !write(6,*) 'chorizo i updated!'
        icnow=icnow+1           
	    icmov1=icmov1+1	
	    ristra(i)%x=xnew
	    cwtl2r(i,:,:)=cwtl2rtmp(:,:)	
        endif
       
        if ( icnow /= 0 ) then
    ! update l2r
            if (i.le.(nchorizo-2)) then
		    x=ristra(i+1)%x
		    call v6proplr(x,x,-1,cwtl2r(i,:,:),cwtl2r(i+1,:,:))
            else ! i=nchorizo-1
            xr=ristra(i+1)%x
            call v6propl(xr,-1,cwtl2r(i,:,:),cwtl2r(i+1,:,:))	    
            endif               
        endif 
    
	    endif 
	    endif 	 
    enddo 

    if (icnow /= 0) then 
        xl=ristra(0)%x
        xr=ristra(nchorizo)%x
    ! update cwtl2r 
    if (iright <= nchorizo-2) then 
        do k=iright,nchorizo-2
            x=ristra(k+1)%x
            call v6proplr(x,x,-1,cwtl2r(k,:,:),cwtl2r(k+1,:,:))   
        enddo
        call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))	      
    else  
        if (iright == nchorizo-1) then   
            call v6propl(xr,-1,cwtl2r(iright,:,:),cwtl2r(nchorizo,:,:))	          
        endif 
    endif     
    ! update cwtr2l from iright th bead until the 0th.
    if (iright == nchorizo) then    
        !x=ristra(iright)%x   
	    call v6propr(xr,-1,cwtend(:,:),cwtr2l(iright,:,:))
	    do k=iright-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	    
    else      
        if (iright >= 2) then  
        x=ristra(iright)%x
        call v6proplr(x,x,-1,cwtr2l(iright+1,:,:),cwtr2l(iright,:,:))
        do k=iright-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	   
        else ! iright =1 or 0
          
        if (iright == 1) then   
        x=ristra(iright)%x
        call v6proplr(x,x,-1,cwtr2l(iright+1,:,:),cwtr2l(iright,:,:))   
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(iright,:,:),cwtr2l(0,:,:))	
        else ! iright =0
        call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	   
        endif
       
        endif  
    endif
    endif
   

    return
    end subroutine stdmove1fastv6a                  
   
   
   
    subroutine stdmove1fastv6(ileft,iright) 
    ! ileft < iright; update from left to right. this subroutine can replace stdmove1rangefast.
    use wavefunction
    use estimator
    use v6stepcalc
    use random
    use brut
    integer(kind=i4) :: i,k,ileft,iright,icnow
    real(kind=r8) :: rn(1),gauss(3,npart)
    real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart),xoldnxt(3,npart),xd2old(3,npart),xd2new(3,npart) &
	                ,xr(3,npart),xl(3,npart)
    real(kind=r8) :: rt2,tmp1,tmp2,logratio
    complex(kind=r8) :: cwtold(0:nspin-1,nisospin),cwtnew(0:nspin-1,nisospin)
    complex(kind=r8) :: cwtl2rtmp(0:nspin-1,nisospin),cwtr2ltmp(0:nspin-1,nisospin)

    if ((ileft.lt.0).or.(iright.gt.nchorizo).or.(iright.lt.0).or.(ileft.gt.nchorizo).or.(ileft.gt.iright)) then
	    write (6,*) 'stdmove1fastv6 range error', ileft,iright
	    stop
    endif
       
    ristraold(0:nchorizo)=ristra(0:nchorizo)
    icstdmove1tot=icstdmove1tot+1
    icnow=0
   
    do i=ileft,iright
      
        !write(6,*) 'i=',i
       
	    icmov1tot=icmov1tot+1  
	    xold=ristraold(i)%x 
    ! ran move	 
	    !rn3=reshape(randn(3*npart),shape(rn3)) 
	    !xnew=xold+mov1step*(rn3-0.5_r8) 
    ! gaussian move	 
	 
        gauss=mov1step*reshape(gaussian(3*npart),(/3,npart/))
     
        !gauss=mov1step*reshape((randn(3*npart)-0.5),(/3,npart/)) ! this is not gauss move. for vmc(dt=0) only
     
	    xnew=xold+gauss
	 
        if (i.eq.0) then
		 
	    !call psitcwt(xold,cwtold)
	  
	    cwtold=cwtbegin
	  
	    tmp1=log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(i,:,:)))))
	    !call psitcwt(xnew,cwtnew)  
	    call corop(xnew,iplo,jplo,cwtnew)
	    call v6propl(xnew,-1,cwtr2l(i+1,:,:),cwtr2ltmp(:,:))
	    tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
	    xd2new=(ristra(i+1)%x-xnew)**2
	    xd2old=(ristra(i+1)%x-xold)**2
	  
	    rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	  
        !rt2=0. ! vmc(dt=0) only.
	  
	    logratio=tmp2-tmp1+rt2
	    rn=randn(1)
  
	    if (log(rn(1)).lt.logratio) then ! dt ne 0
        icnow=icnow+1 
	    icmov1=icmov1+1	
	    ristra(i)%x=xnew
	    cwtr2l(i,:,:)=cwtr2ltmp(:,:)   
    ! update cwtbegin
        cwtbegin=cwtnew	   
    ! update cwtl2r	   
	    call v6propr(xnew,-1,cwtnew(:,:),cwtl2r(i,:,:)) ! just update the i.
    ! prepare for the next cwtl2r
        xoldnxt=ristra(i+1)%x
        call v6proplr(xoldnxt,xoldnxt,-1,cwtl2r(i,:,:),cwtl2r(i+1,:,:))            
	    endif 
	  
	    else
		 
	    if (i.eq.nchorizo) then   
	    !call psitcwt(xold,cwtold)
	    cwtold=cwtend
	    tmp1=log(abs(real(sum(conjg(cwtl2r(i,:,:))*cwtold(:,:)))))
	    !call psitcwt(xnew,cwtnew)
	    call corop(xnew,ipro,jpro,cwtnew)
	    call v6propl(xnew,-1,cwtl2r(i-1,:,:),cwtl2rtmp(:,:))
	    tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtnew(:,:)))))
	    xd2new=(ristra(i-1)%x-xnew)**2
	    xd2old=(ristra(i-1)%x-xold)**2
	  
	    rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	  
        !rt2=0. ! vmc(dt=0) only.
	  
	    logratio=tmp2-tmp1+rt2
	    rn=randn(1)
	  
	    if (log(rn(1)).lt.logratio) then
        icnow=icnow+1
	    icmov1=icmov1+1	
	    ristra(i)%x=xnew
	    cwtl2r(i,:,:)=cwtl2rtmp(:,:)
    ! update cwtend
        cwtend=cwtnew
	    endif 		     
	   
	    else
	   
        cwtold=cwtl2r(i,:,:)
	    call v6proplr(xnew,xnew,-1,cwtl2r(i-1,:,:),cwtl2rtmp(:,:))
	   
	    tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtr2l(i+1,:,:)))))
	    tmp1=log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(i+1,:,:)))))
	    
	    xd2new=(ristra(i-1)%x-xnew)**2+(ristra(i+1)%x-xnew)**2
	    xd2old=(ristra(i-1)%x-xold)**2+(ristra(i+1)%x-xold)**2
	   
	    rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	
        !rt2=0.  ! vmc(dt=0) only.
	   
	    logratio=tmp2-tmp1+rt2	 
	    rn=randn(1)

	    if (log(rn(1)).lt.logratio) then
            !write(6,*) 'chorizo i updated!'
        icnow=icnow+1           
	    icmov1=icmov1+1	
	    ristra(i)%x=xnew
	    cwtl2r(i,:,:)=cwtl2rtmp(:,:)	
        endif
       
        if ( icnow /= 0 ) then
    ! update l2r
            if (i.le.(nchorizo-2)) then
		    x=ristra(i+1)%x
		    call v6proplr(x,x,-1,cwtl2r(i,:,:),cwtl2r(i+1,:,:))
            else ! i=nchorizo-1
            xr=ristra(i+1)%x
            call v6propl(xr,-1,cwtl2r(i,:,:),cwtl2r(i+1,:,:))	    
            endif               
        endif 
    
	    endif 
	    endif 	 
    enddo 

    if (icnow /= 0) then 
        xl=ristra(0)%x
        xr=ristra(nchorizo)%x
    ! update cwtl2r 
    if (iright <= nchorizo-2) then 
        do k=iright,nchorizo-2
            x=ristra(k+1)%x
            call v6proplr(x,x,-1,cwtl2r(k,:,:),cwtl2r(k+1,:,:))   
        enddo
        call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))	      
    else  
        if (iright == nchorizo-1) then   
            call v6propl(xr,-1,cwtl2r(iright,:,:),cwtl2r(nchorizo,:,:))	          
        endif 
    endif     
    ! update cwtr2l from iright th bead until the 0th.
    if (iright == nchorizo) then    
        !x=ristra(iright)%x   
	    call v6propr(xr,-1,cwtend(:,:),cwtr2l(iright,:,:))
	    do k=iright-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	    
    else      
        if (iright >= 2) then  
        x=ristra(iright)%x
        call v6proplr(x,x,-1,cwtr2l(iright+1,:,:),cwtr2l(iright,:,:))
        do k=iright-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	   
        else ! iright =1 or 0
          
        if (iright == 1) then   
        x=ristra(iright)%x
        call v6proplr(x,x,-1,cwtr2l(iright+1,:,:),cwtr2l(iright,:,:))   
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(iright,:,:),cwtr2l(0,:,:))	
        else ! iright =0
        call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	   
        endif
       
        endif  
    endif
    endif
   

    return
    end subroutine stdmove1fastv6               
   
 
    subroutine stdmove1rangefastv6(ileft,iright) ! this is not bisection, just propose a move, and not efficient.
    use wavefunction
    use estimator
    use v6stepcalc
    use random
    use brut
    integer(kind=i4) :: i,k,ileft,iright,istep
    real(kind=r8) :: rn(1),gauss(3,npart)
    real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart),xd2old(3,npart),xd2new(3,npart) &
	                ,xr(3,npart),xl(3,npart)
    real(kind=r8) :: rt2,tmp1,tmp2,logratio
    complex(kind=r8) :: cwtold(0:nspin-1,nisospin),cwtnew(0:nspin-1,nisospin)
    complex(kind=r8) :: cwtl2rtmp(0:nspin-1,nisospin),cwtr2ltmp(0:nspin-1,nisospin)
   
    ristraold(0:nchorizo)=ristra(0:nchorizo)
    icstdmove1tot=icstdmove1tot+1
    if ((ileft.lt.0).or.(iright.gt.nchorizo).or.(iright.lt.0).or.(ileft.gt.nchorizo)) then
	    write (6,*) 'stdmove1 range error', ileft,iright
	    stop
    else
	    if ( ileft.le.iright ) then 
		    istep=1
	    else
		    istep=-1
	    endif 
    endif
    do i=ileft,iright,istep
      
        !write(6,*) 'i=',i
       
	    icmov1tot=icmov1tot+1  
	    xold=ristraold(i)%x 
    ! ran move	 
	    !rn3=reshape(randn(3*npart),shape(rn3)) 
	    !xnew=xold+mov1step*(rn3-0.5_r8) 
    ! gaussian move	 
	 
        gauss=mov1step*reshape(gaussian(3*npart),(/3,npart/))
     
        !gauss=mov1step*reshape((randn(3*npart)-0.5),(/3,npart/)) ! this is not gauss move. for vmc(dt=0) only
     
	    xnew=xold+gauss
	 
        if (i.eq.0) then
		 
	    !call psitcwt(xold,cwtold)
	  
	    cwtold=cwtbegin
	  
	    tmp1=log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(i,:,:)))))
	    !call psitcwt(xnew,cwtnew)  
	    call corop(xnew,iplo,jplo,cwtnew)
	    call v6propl(xnew,-1,cwtr2l(i+1,:,:),cwtr2ltmp(:,:))
	    tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
	    xd2new=(ristra(i+1)%x-xnew)**2
	    xd2old=(ristra(i+1)%x-xold)**2
	  
	    rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	  
        !rt2=0. ! vmc(dt=0) only.
	  
	    logratio=tmp2-tmp1+rt2
	    rn=randn(1)
  
	    if (log(rn(1)).lt.logratio) then ! dt ne 0
         
	    icmov1=icmov1+1	
	    ristra(i)%x=xnew
	    cwtr2l(i,:,:)=cwtr2ltmp(:,:)
    ! update cwtbegin
        cwtbegin=cwtnew	   
    ! update cwtl2r	   
	    call v6propr(xnew,-1,cwtnew(:,:),cwtl2r(i,:,:))
	    do k=i+1,nchorizo-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtl2r(k-1,:,:),cwtl2r(k,:,:))
	    enddo
	    xr=ristra(nchorizo)%x
	    call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))	     
	    !write(6,*) 'move succeed 0',i
     
	    endif 
	  
	    else
		 
	    if (i.eq.nchorizo) then   
	    !call psitcwt(xold,cwtold)
	    cwtold=cwtend
	    tmp1=log(abs(real(sum(conjg(cwtl2r(i,:,:))*cwtold(:,:)))))
	    !call psitcwt(xnew,cwtnew)
	    call corop(xnew,ipro,jpro,cwtnew)
	    call v6propl(xnew,-1,cwtl2r(i-1,:,:),cwtl2rtmp(:,:))
	    tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtnew(:,:)))))
	    xd2new=(ristra(i-1)%x-xnew)**2
	    xd2old=(ristra(i-1)%x-xold)**2
	  
	    rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	  
        !rt2=0. ! vmc(dt=0) only.
	  
	    logratio=tmp2-tmp1+rt2
	    rn=randn(1)
	  
	    if (log(rn(1)).lt.logratio) then
	    icmov1=icmov1+1	
	    ristra(i)%x=xnew
	    cwtl2r(i,:,:)=cwtl2rtmp(:,:)
    ! update cwtend
        cwtend=cwtnew	   
    ! update cwtr2l	   
	    call v6propr(xnew,-1,cwtnew(:,:),cwtr2l(i,:,:))
	    do k=nchorizo-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
	    xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	    
	    !write(6,*) 'move succeed 0',i
	    endif 		     
	   
	    else
	   
        cwtold=cwtl2r(i,:,:)

	    call v6proplr(xnew,xnew,-1,cwtl2r(i-1,:,:),cwtl2rtmp(:,:))
	   
	    tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtr2l(i+1,:,:)))))
	    tmp1=log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(i+1,:,:)))))
	    
	    xd2new=(ristra(i-1)%x-xnew)**2+(ristra(i+1)%x-xnew)**2
	    xd2old=(ristra(i-1)%x-xold)**2+(ristra(i+1)%x-xold)**2
	   
	    rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	
        !rt2=0.  ! vmc(dt=0) only.
	   
	    logratio=tmp2-tmp1+rt2	 
	    rn=randn(1)

	    if (log(rn(1)).lt.logratio) then
           
            !write(6,*) 'chorizo i updated!'
           
	    icmov1=icmov1+1	
	    ristra(i)%x=xnew
	    cwtl2r(i,:,:)=cwtl2rtmp(:,:)	

        xr=ristra(nchorizo)%x
	    xl=ristra(0)%x
    ! update l2r
        if (i.le.(nchorizo-2)) then
	    do k=i+1,nchorizo-1
		    x=ristra(k)%x
		    call v6proplr(x,x,-1,cwtl2r(k-1,:,:),cwtl2r(k,:,:))
	    enddo
	    call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))	
        else ! i=nchorizo-1
        call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))	    
        endif     
    ! update r2l      
	    do k=i,1,-1
            x=ristra(k)%x
		    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))			
	    enddo
	    call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	
	    endif    
	    endif 
	    endif 	 
    enddo 
    !if (mod(icstdmove1tot,nstepdecor).eq.0) then
    !!write (6,*) 'icmov1tot=',icmov1tot
    !ice=ice+1
    !call hpsi(ristra(0))
    !eleft=ristra(0)%elocal  
    !call hpsi(ristra(nchorizo))
    !eright=ristra(nchorizo)%elocal 
    !select case (irep)
    !case (1:2)
    !   call addval(3,eleft,1.0d0)
    !   call addval(4,eright,1.0d0)
	    !call funcrmsq(ristra(nchorizomid)%x,rm,rmsq)
	    !call addval(5,rmsq,1.0d0)
	    !write(11,'(2x,3(G12.5,2x))') eleft,eright,(eleft+eright)/2.,rm!11 means write to 'values.txt'
    !case (3:4)
    !   call addval(4,eleft,1.0d0)
    !   call addval(5,eright,1.0d0)
	    !write(11,'(2x,3(G12.5,2x))') eleft,eright,(eleft+eright)/2.!11 means write to 'values.txt'
    !end select
    !endif 
    return
    end subroutine stdmove1rangefastv6         
    
    subroutine stdmove2v6(ileft,iright) ! move2 means moving all the beads at the same time 
    ! ileft < iright; update from left to right. this subroutine can replace stdmove1rangefast.
    use wavefunction
    use estimator
    use v6stepcalc
    use random
    use brut
    integer(kind=i4) :: i,k,ileft,iright
    real(kind=r8) :: rn(1),gauss(3,npart)
    real(kind=r8) :: x(3,npart),xd2old(3,npart),xd2new(3,npart) &
	                ,xr(3,npart),xl(3,npart)
    real(kind=r8) :: rt2(0:nchorizo),tmp1,tmp2,logratio

    if ((ileft.lt.0).or.(iright.gt.nchorizo).or.(iright.lt.0).or.(ileft.gt.nchorizo).or.(ileft.gt.iright)) then
	    write (6,*) 'stdmove1fastv6 range error', ileft,iright
	    stop
    else
    if (ileft.eq.iright) then
        call stdmove1fastv6(ileft,iright) 
        return
    endif          
    endif

    ! now, ileft<iright   
   
    ristraold(0:nchorizo)=ristra(0:nchorizo)
    icstdmove2tot=icstdmove2tot+1   
   
    do i=ileft,iright
    ! ran move	 
	    !rn3=reshape(randn(3*npart),shape(rn3)) 
	    !xnew=xold+mov1step*(rn3-0.5_r8) 
    ! gaussian move	 
        gauss=mov2step*reshape(gaussian(3*npart),(/3,npart/))
        !gauss=mov1step*reshape((randn(3*npart)-0.5),(/3,npart/)) ! this is not gauss move. for vmc(dt=0) only
	    ristranew(i)%x=ristraold(i)%x+gauss    
    enddo
   
    tmp1=log(abs(real(sum(conjg(cwtl2r(0,:,:))*cwtr2l(1,:,:)))))
    tmp2=tmp1   
    rt2=0.
    cwtl2rnew=cwtl2r
    cwtr2lnew=cwtr2l
   
    do i=ileft,iright  
   
        if (i.eq.0) then
		 
	    !call psitcwt(xold,cwtold)  
	    call corop(ristranew(i)%x,iplo,jplo,cwtbeginnew)
	    call v6propr(ristranew(i)%x,-1,cwtbeginnew(:,:),cwtl2rnew(i,:,:)) ! just update the i.
      
	    xd2new=(  ristranew(i)%x - ristranew(i+1)%x  )**2
	    xd2old=(  ristraold(i)%x - ristraold(i+1)%x  )**2
	    rt2(i)=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	  
        !rt2=0. ! vmc(dt=0) only.
	     
	    else
		 
	    if (i.eq.nchorizo) then   
	    !call psitcwt(xold,cwtold)

	    !call psitcwt(xnew,cwtnew)
	    call corop(ristranew(i)%x,ipro,jpro,cwtendnew)
	    call v6propl(ristranew(i)%x,-1,cwtl2rnew(i-1,:,:),cwtl2rnew(i,:,:))
        call v6propr(ristranew(i)%x,-1,cwtendnew,cwtr2lnew(i,:,:))
         
	    else
	   
	    call v6proplr(ristranew(i)%x,ristranew(i)%x,-1,cwtl2rnew(i-1,:,:),cwtl2rnew(i,:,:))
	    
        xd2new=(  ristranew(i)%x - ristranew(i+1)%x  )**2
	    xd2old=(  ristraold(i)%x - ristraold(i+1)%x  )**2
	    rt2(i)=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	
        !rt2=0.  ! vmc(dt=0) only.

	    endif 
	    endif 	 
    enddo 

    if (iright /= nchorizo) then
    tmp2=log(abs(real(sum(conjg(cwtl2rnew(iright,:,:))*cwtr2lnew(iright+1,:,:)))))
    else
    tmp2=log(abs(real(sum(conjg(cwtl2rnew(iright,:,:))*cwtendnew(:,:)))))   
    endif
   
    if (ileft/=0) then
        xd2new=(  ristranew(ileft)%x - ristraold(ileft-1)%x  )**2
	    xd2old=(  ristraold(ileft)%x - ristraold(ileft-1)%x  )**2       
        rt2(ileft-1)=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	
    endif
   
    logratio=tmp2-tmp1+sum(rt2)
    rn=randn(1)
   
    if (log(rn(1)).lt.logratio) then 
       
        icstdmove2=icstdmove2+1
       
        do i=ileft,iright
        ristra(i)=ristranew(i)
        cwtl2r(i,:,:)=cwtl2rnew(i,:,:)
        enddo 
      
        if (ileft.eq.0) cwtbegin=cwtbeginnew
        if (iright.eq.nchorizo) cwtend=cwtendnew 
       
    ! update the chain.
        xl=ristra(0)%x
        xr=ristra(nchorizo)%x
    ! update cwtl2r 
    if (iright <= nchorizo-2) then 
        do k=iright,nchorizo-2
            x=ristra(k+1)%x
            call v6proplr(x,x,-1,cwtl2r(k,:,:),cwtl2r(k+1,:,:))   
        enddo
        call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))	      
    else  
        if (iright == nchorizo-1) then   
            call v6propl(xr,-1,cwtl2r(iright,:,:),cwtl2r(nchorizo,:,:))	          
        endif 
    endif     
    ! update cwtr2l from iright th bead until the 0th.
    if (iright == nchorizo) then    
        !x=ristra(iright)%x   
	    call v6propr(xr,-1,cwtend(:,:),cwtr2l(iright,:,:))
	    do k=iright-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	    
    else      
        if (iright >= 2) then  
        x=ristra(iright)%x
        call v6proplr(x,x,-1,cwtr2l(iright+1,:,:),cwtr2l(iright,:,:))
        do k=iright-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	   
        else ! iright =1 or 0
          
        if (iright == 1) then   
        x=ristra(iright)%x
        call v6proplr(x,x,-1,cwtr2l(iright+1,:,:),cwtr2l(iright,:,:))   
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(iright,:,:),cwtr2l(0,:,:))	
        else ! iright =0
        call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	   
        endif
       
        endif  
    endif
   
       
    endif
   

    return
    end subroutine stdmove2v6             
   
    subroutine stdmove3v6(ileft,iright) ! move3 means shift all the beads at the same time 
    ! ileft < iright; update from left to right. this subroutine can replace stdmove1rangefast.
    use wavefunction
    use estimator
    use v6stepcalc
    use random
    use brut
    integer(kind=i4) :: i,k,ileft,iright
    real(kind=r8) :: rn(1),gauss(3,npart)
    real(kind=r8) :: x(3,npart),xd2old(3,npart),xd2new(3,npart) &
	                ,xr(3,npart),xl(3,npart)
    real(kind=r8) :: rt2(0:nchorizo),tmp1,tmp2,logratio

    if ((ileft.lt.0).or.(iright.gt.nchorizo).or.(iright.lt.0).or.(ileft.gt.nchorizo).or.(ileft.gt.iright)) then
	    write (6,*) 'stdmove1fastv6 range error', ileft,iright
	    stop
    else
    if (ileft.eq.iright) then
        call stdmove1fastv6(ileft,iright) 
        return
    endif          
    endif

    ! now, ileft<iright   
   
    ristraold(0:nchorizo)=ristra(0:nchorizo)
    icshfttot=icshfttot+1   
    gauss=shftstep*reshape(gaussian(3*npart),(/3,npart/))
   
    do i=ileft,iright
    ! ran move	 
	    !rn3=reshape(randn(3*npart),shape(rn3)) 
	    !xnew=xold+mov1step*(rn3-0.5_r8) 
    ! gaussian move	 
     
        !gauss=mov1step*reshape((randn(3*npart)-0.5),(/3,npart/)) ! this is not gauss move. for vmc(dt=0) only
	    ristranew(i)%x=ristraold(i)%x+gauss    
    enddo ! just the shift. all beads shifted by the same amount.
   
    tmp1=log(abs(real(sum(conjg(cwtl2r(0,:,:))*cwtr2l(1,:,:)))))
    tmp2=tmp1   
    rt2=0.
    cwtl2rnew=cwtl2r
    cwtr2lnew=cwtr2l
   
    do i=ileft,iright  
   
        if (i.eq.0) then
		 
	    !call psitcwt(xold,cwtold)  
	    call corop(ristranew(i)%x,iplo,jplo,cwtbeginnew)
	    call v6propr(ristranew(i)%x,-1,cwtbeginnew(:,:),cwtl2rnew(i,:,:)) ! just update the i.
      
	    xd2new=(  ristranew(i)%x - ristranew(i+1)%x  )**2
	    xd2old=(  ristraold(i)%x - ristraold(i+1)%x  )**2
	    rt2(i)=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	  
        !rt2=0. ! vmc(dt=0) only.
	     
	    else
		 
	    if (i.eq.nchorizo) then   
	    !call psitcwt(xold,cwtold)

	    !call psitcwt(xnew,cwtnew)
	    call corop(ristranew(i)%x,ipro,jpro,cwtendnew)
	    call v6propl(ristranew(i)%x,-1,cwtl2rnew(i-1,:,:),cwtl2rnew(i,:,:))
        call v6propr(ristranew(i)%x,-1,cwtendnew,cwtr2lnew(i,:,:))
         
	    else
	   
	    call v6proplr(ristranew(i)%x,ristranew(i)%x,-1,cwtl2rnew(i-1,:,:),cwtl2rnew(i,:,:))
	    
        xd2new=(  ristranew(i)%x - ristranew(i+1)%x  )**2
	    xd2old=(  ristraold(i)%x - ristraold(i+1)%x  )**2
	    rt2(i)=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	
        !rt2=0.  ! vmc(dt=0) only.

	    endif 
	    endif 	 
    enddo 

    if (iright /= nchorizo) then
    tmp2=log(abs(real(sum(conjg(cwtl2rnew(iright,:,:))*cwtr2lnew(iright+1,:,:)))))
    else
    tmp2=log(abs(real(sum(conjg(cwtl2rnew(iright,:,:))*cwtendnew(:,:)))))   
    endif
   
    if (ileft/=0) then
        xd2new=(  ristranew(i)%x - ristranew(i-1)%x  )**2
	    xd2old=(  ristraold(i)%x - ristraold(i-1)%x  )**2       
        rt2(ileft-1)=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	
    endif
   
    logratio=tmp2-tmp1+sum(rt2)
    rn=randn(1)
   
    if (log(rn(1)).lt.logratio) then 
 
        icshft=icshft+1   
       
        do i=ileft,iright
        ristra(i)=ristranew(i)
        cwtl2r(i,:,:)=cwtl2rnew(i,:,:)
        enddo 
      
        if (ileft.eq.0) cwtbegin=cwtbeginnew
        if (iright.eq.nchorizo) cwtend=cwtendnew 
       
    ! update the chain.
        xl=ristra(0)%x
        xr=ristra(nchorizo)%x
    ! update cwtl2r 
    if (iright <= nchorizo-2) then 
        do k=iright,nchorizo-2
            x=ristra(k+1)%x
            call v6proplr(x,x,-1,cwtl2r(k,:,:),cwtl2r(k+1,:,:))   
        enddo
        call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))	      
    else  
        if (iright == nchorizo-1) then   
            call v6propl(xr,-1,cwtl2r(iright,:,:),cwtl2r(nchorizo,:,:))	          
        endif 
    endif     
    ! update cwtr2l from iright th bead until the 0th.
    if (iright == nchorizo) then    
        !x=ristra(iright)%x   
	    call v6propr(xr,-1,cwtend(:,:),cwtr2l(iright,:,:))
	    do k=iright-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	    
    else      
        if (iright >= 2) then  
        x=ristra(iright)%x
        call v6proplr(x,x,-1,cwtr2l(iright+1,:,:),cwtr2l(iright,:,:))
        do k=iright-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	   
        else ! iright =1 or 0
          
        if (iright == 1) then   
        x=ristra(iright)%x
        call v6proplr(x,x,-1,cwtr2l(iright+1,:,:),cwtr2l(iright,:,:))   
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(iright,:,:),cwtr2l(0,:,:))	
        else ! iright =0
        call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	   
        endif
       
        endif  
    endif
   
       
    endif
   

    return
    end subroutine stdmove3v6                
   
    subroutine stdmove23v6(ileft,iright) ! move3 means shift all the beads at the same time 
    ! ileft < iright; update from left to right. this subroutine can replace stdmove1rangefast.
    use wavefunction
    use estimator
    use v6stepcalc
    use random
    use brut
    integer(kind=i4) :: i,k,ileft,iright
    real(kind=r8) :: rn(1),gauss(3,npart),gauss1(3,npart)
    real(kind=r8) :: x(3,npart),xd2old(3,npart),xd2new(3,npart) &
	                ,xr(3,npart),xl(3,npart)
    real(kind=r8) :: rt2(0:nchorizo),tmp1,tmp2,logratio

    if ((ileft.lt.0).or.(iright.gt.nchorizo).or.(iright.lt.0).or.(ileft.gt.nchorizo).or.(ileft.gt.iright)) then
	    write (6,*) 'stdmove1fastv6 range error', ileft,iright
	    stop
    else
    if (ileft.eq.iright) then
        call stdmove1fastv6(ileft,iright) 
        return
    endif          
    endif

    ! now, ileft<iright   
   
    ristraold(0:nchorizo)=ristra(0:nchorizo)
    ic23tot=ic23tot+1   
   
    gauss=shftstep*reshape(gaussian(3*npart),(/3,npart/))
    do i=ileft,iright
    ! ran move	 
	    !rn3=reshape(randn(3*npart),shape(rn3)) 
	    !xnew=xold+mov1step*(rn3-0.5_r8) 
    ! gaussian move	 
        gauss1=mov2step*reshape(gaussian(3*npart),(/3,npart/))
        !gauss=mov1step*reshape((randn(3*npart)-0.5),(/3,npart/)) ! this is not gauss move. for vmc(dt=0) only
	    ristranew(i)%x=ristraold(i)%x+gauss1+gauss    
    enddo ! just the shift. all beads shifted by the same amount.
   
    tmp1=log(abs(real(sum(conjg(cwtl2r(0,:,:))*cwtr2l(1,:,:)))))
    tmp2=tmp1   
    rt2=0.
    cwtl2rnew=cwtl2r
    cwtr2lnew=cwtr2l
   
    do i=ileft,iright  
   
        if (i.eq.0) then
		 
	    !call psitcwt(xold,cwtold)  
	    call corop(ristranew(i)%x,iplo,jplo,cwtbeginnew)
	    call v6propr(ristranew(i)%x,-1,cwtbeginnew(:,:),cwtl2rnew(i,:,:)) ! just update the i.
      
	    xd2new=(  ristranew(i)%x - ristranew(i+1)%x  )**2
	    xd2old=(  ristraold(i)%x - ristraold(i+1)%x  )**2
	    rt2(i)=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	  
        !rt2=0. ! vmc(dt=0) only.
	     
	    else
		 
	    if (i.eq.nchorizo) then   
	    !call psitcwt(xold,cwtold)

	    !call psitcwt(xnew,cwtnew)
	    call corop(ristranew(i)%x,ipro,jpro,cwtendnew)
	    call v6propl(ristranew(i)%x,-1,cwtl2rnew(i-1,:,:),cwtl2rnew(i,:,:))
        call v6propr(ristranew(i)%x,-1,cwtendnew,cwtr2lnew(i,:,:))
         
	    else
	   
	    call v6proplr(ristranew(i)%x,ristranew(i)%x,-1,cwtl2rnew(i-1,:,:),cwtl2rnew(i,:,:))
	    
        xd2new=(  ristranew(i)%x - ristranew(i+1)%x  )**2
	    xd2old=(  ristraold(i)%x - ristraold(i+1)%x  )**2
	    rt2(i)=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	
        !rt2=0.  ! vmc(dt=0) only.

	    endif 
	    endif 	 
    enddo 

    if (iright /= nchorizo) then
    tmp2=log(abs(real(sum(conjg(cwtl2rnew(iright,:,:))*cwtr2lnew(iright+1,:,:)))))
    else
    tmp2=log(abs(real(sum(conjg(cwtl2rnew(iright,:,:))*cwtendnew(:,:)))))   
    endif
   
    if (ileft/=0) then
        xd2new=(  ristranew(i)%x - ristranew(i-1)%x  )**2
	    xd2old=(  ristraold(i)%x - ristraold(i-1)%x  )**2       
        rt2(ileft-1)=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	
    endif
   
    logratio=tmp2-tmp1+sum(rt2)
    rn=randn(1)
   
    if (log(rn(1)).lt.logratio) then 
 
        ic23=ic23+1   
       
        do i=ileft,iright
        ristra(i)=ristranew(i)
        cwtl2r(i,:,:)=cwtl2rnew(i,:,:)
        enddo 
      
        if (ileft.eq.0) cwtbegin=cwtbeginnew
        if (iright.eq.nchorizo) cwtend=cwtendnew 
       
    ! update the chain.
        xl=ristra(0)%x
        xr=ristra(nchorizo)%x
    ! update cwtl2r 
    if (iright <= nchorizo-2) then 
        do k=iright,nchorizo-2
            x=ristra(k+1)%x
            call v6proplr(x,x,-1,cwtl2r(k,:,:),cwtl2r(k+1,:,:))   
        enddo
        call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))	      
    else  
        if (iright == nchorizo-1) then   
            call v6propl(xr,-1,cwtl2r(iright,:,:),cwtl2r(nchorizo,:,:))	          
        endif 
    endif     
    ! update cwtr2l from iright th bead until the 0th.
    if (iright == nchorizo) then    
        !x=ristra(iright)%x   
	    call v6propr(xr,-1,cwtend(:,:),cwtr2l(iright,:,:))
	    do k=iright-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	    
    else      
        if (iright >= 2) then  
        x=ristra(iright)%x
        call v6proplr(x,x,-1,cwtr2l(iright+1,:,:),cwtr2l(iright,:,:))
        do k=iright-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	   
        else ! iright =1 or 0
          
        if (iright == 1) then   
        x=ristra(iright)%x
        call v6proplr(x,x,-1,cwtr2l(iright+1,:,:),cwtr2l(iright,:,:))   
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(iright,:,:),cwtr2l(0,:,:))	
        else ! iright =0
        call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	   
        endif
       
        endif  
    endif
   
       
    endif
   

    return
    end subroutine stdmove23v6                   
   

    subroutine stdvmcmove2v6 ! move2 means moving all the beads at the same time, vmc.
    ! ileft < iright; update from left to right. this subroutine can replace stdmove1rangefast.  
    use wavefunction
    use estimator
    use v6stepcalc
    use random
    use brut
    integer(kind=i4) :: i,k,ileft,iright
    real(kind=r8) :: rn(1),gauss(3,npart)
    real(kind=r8) :: x(3,npart),xr(3,npart),xl(3,npart)
    real(kind=r8) :: rt2(0:nchorizo),tmp1,tmp2,logratio

   
   
    ileft=0
    iright=nchorizo ! we do not set nchorizo as 0.
   
    if ((ileft.lt.0).or.(iright.gt.nchorizo).or.(iright.lt.0).or.(ileft.gt.nchorizo).or.(ileft.gt.iright)) then
	    write (6,*) 'stdmove1fastv6 range error', ileft,iright
	    stop      
    endif

    ! now, ileft<iright   
   
    ristraold(0:nchorizo)=ristra(0:nchorizo)
    icstdmove2tot=icstdmove2tot+1   
   
    gauss=mov2step*reshape((randn(3*npart)-0.5),(/3,npart/)) ! this is not gauss move. for vmc(dt=0) only
    do i=ileft,iright   
	    ristranew(i)%x=ristraold(i)%x+gauss    
    enddo
   
    tmp1=log(abs(real(sum(conjg(cwtl2r(0,:,:))*cwtr2l(1,:,:)))))
    tmp2=tmp1   
    rt2=0.
    cwtl2rnew=cwtl2r
    cwtr2lnew=cwtr2l
   
    do i=ileft,iright  
   
        if (i.eq.0) then
		 
	    !call psitcwt(xold,cwtold)  
	    call corop(ristranew(i)%x,iplo,jplo,cwtbeginnew)
	    call v6propr(ristranew(i)%x,-1,cwtbeginnew(:,:),cwtl2rnew(i,:,:)) ! just update the i.
      
	    !xd2new=(  ristranew(i)%x - ristranew(i+1)%x  )**2
	    !xd2old=(  ristraold(i)%x - ristraold(i+1)%x  )**2
	    !rt2(i)=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	  
        !rt2=0. ! vmc(dt=0) only.
	     
	    else
		 
	    if (i.eq.nchorizo) then   
	    !call psitcwt(xold,cwtold)

	    !call psitcwt(xnew,cwtnew)
	    call corop(ristranew(i)%x,ipro,jpro,cwtendnew)
	    call v6propl(ristranew(i)%x,-1,cwtl2rnew(i-1,:,:),cwtl2rnew(i,:,:))
        call v6propr(ristranew(i)%x,-1,cwtendnew,cwtr2lnew(i,:,:))
         
	    else
	   
	    call v6proplr(ristranew(i)%x,ristranew(i)%x,-1,cwtl2rnew(i-1,:,:),cwtl2rnew(i,:,:))
	    
        !xd2new=(  ristranew(i)%x - ristranew(i+1)%x  )**2
	    !xd2old=(  ristraold(i)%x - ristraold(i+1)%x  )**2
	    !rt2(i)=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	
        !rt2=0.  ! vmc(dt=0) only.

	    endif 
	    endif 	 
    enddo 

    tmp2=log(abs(real(sum(conjg(cwtl2rnew(iright,:,:))*cwtendnew(:,:)))))   

   
    logratio=tmp2-tmp1+sum(rt2)
    rn=randn(1)
   
    if (log(rn(1)).lt.logratio) then 
       
        icstdmove2=icstdmove2+1
       
        do i=ileft,iright
        ristra(i)=ristranew(i)
        cwtl2r(i,:,:)=cwtl2rnew(i,:,:)
        enddo 
      
        if (ileft.eq.0) cwtbegin=cwtbeginnew
        if (iright.eq.nchorizo) cwtend=cwtendnew 
       
    ! update the chain.
        xl=ristra(0)%x
        xr=ristra(nchorizo)%x
    ! update cwtl2r 
    if (iright <= nchorizo-2) then 
        do k=iright,nchorizo-2
            x=ristra(k+1)%x
            call v6proplr(x,x,-1,cwtl2r(k,:,:),cwtl2r(k+1,:,:))   
        enddo
        call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))	      
    else  
        if (iright == nchorizo-1) then   
            call v6propl(xr,-1,cwtl2r(iright,:,:),cwtl2r(nchorizo,:,:))	          
        endif 
    endif     
    ! update cwtr2l from iright th bead until the 0th.
    if (iright == nchorizo) then    
        !x=ristra(iright)%x   
	    call v6propr(xr,-1,cwtend(:,:),cwtr2l(iright,:,:))
	    do k=iright-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	    
    else      
        if (iright >= 2) then  
        x=ristra(iright)%x
        call v6proplr(x,x,-1,cwtr2l(iright+1,:,:),cwtr2l(iright,:,:))
        do k=iright-1,1,-1
	    x=ristra(k)%x
	    call v6proplr(x,x,-1,cwtr2l(k+1,:,:),cwtr2l(k,:,:))
	    enddo
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	   
        else ! iright =1 or 0
          
        if (iright == 1) then   
        x=ristra(iright)%x
        call v6proplr(x,x,-1,cwtr2l(iright+1,:,:),cwtr2l(iright,:,:))   
	    !xl=ristra(0)%x
	    call v6propl(xl,-1,cwtr2l(iright,:,:),cwtr2l(0,:,:))	
        else ! iright =0
        call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))	   
        endif
       
        endif  
    endif
    
    call addval(2,1.0_r8,1.0_r8)  
    else
    call addval(2,0.0_r8,1.0_r8)   
    endif
   

    return
    end subroutine stdvmcmove2v6             
   
   
    !   subroutine bisectnchorizov6 ! * to be continued if needed...
    !! can always set nbisect < nchorizo. don't really it but can be added later. 
    !! just the nbisect=nchorizo case, and just the move   
    !   use wavefunction
    !   use estimator
    !   use random
    !   integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,irightbisect,imidbisect,k,bilvl,ipick,na,nb &
    !   ,nbisect1,mmax1
    !   real(kind=r8) :: rn(1),t,tau,sigmamid,gauss(3,npart),prob,eleft,eright,rm,rmsq
    !   real(kind=r8) :: xmid(3,npart)
    !   real(kind=r8) :: rn3(3,npart),xnew(3,npart),xold(3,npart),xoldnxt(3,npart)
    !   real(kind=r8) :: tmp,tmp1,tmp2,tmpq,logratio,logval,logval2 &
    !	               ,logno,p1n,p1,rt2,p2n,p2,rt1,p3n,p3
    !   real(kind=r8) :: vdiff(0:mmax),vsumold(0:mmax),vsumnew(0:mmax)
    !   logical :: reject
    !   if (nbisect.ne.nchorizo) then
    !	 write (6,*) 'nbisect ne nchorizo, cannot call bisectnchorizo here!', nbisect,nchorizo
    !	 stop 
    !   else
    !	ibisecttot=ibisecttot+1
    !    ileftbisect=0
    !    irightbisect=nchorizo
    !! need to show the successful rate for the 0 and nchorizo movement. better make it 50%.   
    !! bisection   
    !   oldristra(0:nbisect)=ristra(ileftbisect:ileftbisect+nbisect)
    !   newristra(0)=oldristra(0)
    !   newristra(nbisect)=oldristra(nbisect) 
    !   tau=dt*nbisect ! total time for the bisection slice.
    !   all: do bilvl=1,mmax   ! level 1 to N. mmax = the total level N.
    !	  sigmamid=sqrt(hbar*tau*2.0d0**(-bilvl))  
    !	  bisecttot(bilvl)=bisecttot(bilvl)+1.
    !	   !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
    !	   do i=0,2**(bilvl-1)-1
    !		   imid=2**(mmax-bilvl)+i*2**(mmax-bilvl+1)
    !		   il=imid-2**(mmax-bilvl)
    !		   ir=imid+2**(mmax-bilvl)
    !		   !write(6,*) 'bilvl il ir imid sigmamid=',bilvl, il,ir,imid,sigmamid
    !		   gauss=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
    !		   xmid=(newristra(il)%x+newristra(ir)%x)/2.
    !		   newristra(imid)%x=xmid+gauss
    !		   call hpsi(newristra(imid))
    !		   call hpsi(oldristra(imid)) ! put here for save. But is not really needed if do not move beads before or after bisection
    !		   tmp=(0.5d0*(newristra(il)%v+newristra(ir)%v)-newristra(imid)%v &
    !               -0.5d0*(oldristra(il)%v+oldristra(ir)%v)+oldristra(imid)%v)*tau/(2.0d0**bilvl)
    !		   rn=randn(1)
    !	       !write(6,*) 'tmp ln rn=',tmp,log(rn(1))
    !	       if (log(rn(1)).lt.tmp) then
    !	     	reject=.false.  
    !			 if (tmp.gt.1.0d10) then
    !			  write(6,*) 'bug bead appear! exit!'
    !		      write(6,*) 'ileft bilvl il ir imid sigmamid=',ileftbisect,bilvl, il,ir,imid,sigmamid	   
    !	          write(6,*) 'tmp ln(rn)=',tmp,log(rn(1))		   
    !		      write (6,*) 'v=',newristra(il)%v,newristra(ir)%v,newristra(imid)%v &
    !			                ,oldristra(il)%v,oldristra(ir)%v,oldristra(imid)%v
    !		      stop
    !		     endif
    !	       else
    !	    	reject=.true.  
    !	    	exit all 
    !	       endif 
    !		   !write(6,*) 'imid=',imid
    !	   enddo
    !	   bisectcount(bilvl)=bisectcount(bilvl)+1.
    !   enddo all
    !   if (reject) then
    !	 call addval(2,0.0d0,1.0_r8)
    !   else
    !	 ibisect=ibisect+1 
    !	 do i=0,nbisect
    !         ristra(ileftbisect+i)=newristra(i)
    !	 enddo
    !	 call addval(2,1.0_r8,1.0_r8)   
    !   endif 
    !	 call stdmove1 ! added the gaussian move	
    !   endif	   
    !   return  
    !   end subroutine bisectnchorizov6   
    !   
   

    subroutine corrchk(x,ixd,corrin,lmax)
! Interacting Electrons, Ceperley et al, p585. MC in ab initio quantum chemstry, Hammond et al, p59.    
    use wavefunction
    use estimator
    use random
    integer(kind=i4) :: ixd,lmax,lcorr,icorr,j! lmax <= ixd, usally set ixd=lmax=nstepdecor
    real(kind=r8) :: xave,xvar
    real(kind=r8) :: x(:),corrin(:),kappa(ixd-1)
    integer(kind=i4) :: iounit
    logical :: findstepdecorr 
    ! nstepdecor usually set the same as nstep when do the correlation check.
    xave=sum(x(1:ixd))/ixd
    xvar=sum(x(1:ixd)**2)/ixd-xave**2      
    findstepdecorr=.false.
    corrin=0.
    kappa=0
    OPEN(newunit=iounit,FILE='xcorrelationcheck.txt',FORM='FORMATTED')
    WRITE(iounit,'(2x,5(G15.8,2x))') 'k','corr(k)','kappa' 
    do lcorr=1,ixd-1 ! k     
	    do icorr=1,ixd-lcorr ! n
		    corrin(lcorr)=corrin(lcorr)+ (x(icorr)-xave)*(x(icorr+lcorr)-xave)
	    enddo
	    corrin(lcorr)=corrin(lcorr)/(ixd-lcorr)/xvar 
        kappa(lcorr)=1+2*sum(corrin(1:lcorr))
        WRITE(iounit,'(2x,20(G15.8,2x))') lcorr,corrin(lcorr),kappa(lcorr)
    enddo    
    close(iounit) 
	OPEN(newunit=iounit,FILE='xcorrchk.txt',FORM='FORMATTED')
    WRITE(iounit,'(2x,2(G15.8,2x))') 'k','x' 
    do j=1,lmax 
	    WRITE(iounit,'(2x,2(G15.8,2x))') j,x(j)	
    enddo
	close(iounit)
    return  
    end subroutine corrchk       
   
    subroutine counterinitVMC
    ic0=0
    ic0tot=0
    return
    end subroutine counterinitVMC
   
    subroutine showstatVMC
    write (6,'(''VMC move sucess ratio ='',t40,f10.5)') dble(ic0)/dble(ic0tot)
    write (6,'(''E Local ave ='',t40,f15.8)') dble(eloctot)/dble(ieloc)
    return
    end subroutine showstatVMC
    
    subroutine zerepstepmon ! monitor the history of rqmc steps.
    repstepct1=0
    repstepct2=0
    repstepct1sum=0
    repstepct2sum=0  
    icrephist=0
    icrephistlast=0
    end subroutine zerepstepmon
   
    subroutine counterinit ! zero the stats after finishing one block.
    icstepstd=0
    ice=0
    ic0=0
    ic0tot=0
    icn=0
    icntot=0
    icmov1=0
    icmov1tot=0
    icstdmove1tot=0
    icstdmove2=0
    icstdmove2tot=0
    icshft=0
    icshfttot=0
    ic23=0
    ic23tot=0
    ibisect=0
    ibisecttot=0
    icbisecttot=0
    bisectrate=0
    bisectcount=0
    bisecttot=0
    bisectcountl=0
    bisecttotl=0
    bisectcountr=0
    bisecttotr=0
    ibisectl=0
    ibisecttotl=0
    ibisectr=0
    ibisecttotr=0
    iclr=0
    iclrtot=0
    icrep=0
    icreptot=0
    return
    end subroutine counterinit
   
    subroutine zeroeloc
    ieloc=0
    eloctot=0.
    return
    end subroutine zeroeloc
    
    subroutine checkblkstat ! the check stuck config clock.
    use mympi
    integer(kind=i4) :: myunit
    real(kind=r8) :: r(10)
    character(len=30) :: filename    
    real(kind=r8) :: xtot(3,npart,0:nchorizo)
    integer(kind=i4) :: ipl(6),jpl(6),ipr(6),jpr(6)
      
    r(1)=dble(icmov1)/icmov1tot
    r(2)=dble(icstdmove2)/icstdmove2tot
    r(3)=dble(icshft)/icshfttot
    r(4)=dble(ic23)/ic23tot
    r(5)=dble(icn)/icntot
    r(6)=dble(iclr)/iclrtot
    r(7)=dble(icrep)/icreptot
    r(8)=dble(ibisect)/ibisecttot
    r(9)=dble(ibisectl)/ibisecttotl
    r(10)=dble(ibisectr)/ibisecttotr
     
    if ((ibisect.eq.0).or.(ibisectl.eq.0).or.(ibisectr.eq.0)) then
        nblockstuck=nblocknow
        ibisectstuck(nblockstuck)=1        
        if ( nblockstuck >= 5 ) then 
            if (sum(ibisectstuck(nblockstuck-4:nblockstuck)).eq.5) then
                call chorizoallout(xtot,ipl,jpl,ipr,jpr) 
                write(filename,'("myrank",i10,".stuck")') myrank()
                open(newunit=myunit,form='formatted',file=trim(filename),position='rewind')  
                    write (myunit,*) 'nsteps really stuck are ',nblockstuck            
                    write (myunit,'(6i10)') ibisect,ibisecttot,ibisectl,ibisecttotl,ibisectr,ibisecttotr
                    write (myunit,'(6i10)') ipl
                    write (myunit,'(6i10)') jpl
                    write (myunit,'(6i10)') ipr
                    write (myunit,'(6i10)') jpr
                    write (myunit,'(3e15.7)') reshape(xtot,(/3*npart*(nchorizo+1)/))      
                close(myunit)
                write(6,*) 'stuck configuration, abort! myrank= ', myrank()
                call abort                         
            endif       
        endif
    endif

    return
    end subroutine checkblkstat
     
    subroutine showstat
    use mympi
    integer(kind=i4) :: k
    real(kind=r8) :: r(10)
      
    r(1)=dble(icmov1)/icmov1tot
    r(2)=dble(icstdmove2)/icstdmove2tot
    r(3)=dble(icshft)/icshfttot
    r(4)=dble(ic23)/ic23tot
    r(5)=dble(icn)/icntot
    r(6)=dble(iclr)/iclrtot
    r(7)=dble(icrep)/icreptot
    r(8)=dble(ibisect)/ibisecttot
    r(9)=dble(ibisectl)/ibisecttotl
    r(10)=dble(ibisectr)/ibisecttotr
    
    bisectrate(:)=bisectcount(:)/bisecttot(:)
   
    write (6,'(''dt ='',t40,f10.5)') dt
    write (6,'(''nchorizo ='',t40,i10)') nchorizo 
   
    write (6,'(''acceptance ratio (mov1) ='',t40,f10.5,3(i10))') r(1),icmov1,icmov1tot,icstdmove1tot
    write (6,'(''acceptance ratio (mov2) ='',t40,f10.5,2(i10))') r(2),icstdmove2,icstdmove2tot
    write (6,'(''acceptance ratio (mov3 shift) ='',t40,f10.5,2(i10))') r(3),icshft,icshfttot
    write (6,'(''ic23 ratio (stdmov23) ='',t40,f10.5,2(i10))') r(4),ic23,ic23tot
    write (6,'(''icn ratio ='',t40,f10.5,2(i10))') r(5),icn,icntot
    write (6,'(''lr ratio ='',t40,f10.5,2(i10))') r(6),iclr,iclrtot
    write (6,'(''reptation ratio ='',t40,f10.5,3(i10))') r(7),icrep,icreptot
    write (6,'(''ibisect ratio ='',t40,f10.5,3(i10))') r(8),ibisect,ibisecttot,icbisecttot
    write (6,'(''bisection level 1 to mmax ratio ='',t40,f10.5)') (bisectrate(k),k=1,mmax)
    write (6,'(''ibisect lend ratio ='',t40,f10.5,3(i10))') r(9),ibisectl,ibisecttotl
    !write (6,'(''bisection lend level 0 to mmaxnow ratio ='',t40,f10.5)') (bisectcountl(k)/bisecttotl(k),k=0,mmax-1)
    write (6,'(''ibisect rend ratio ='',t40,f10.5,3(i10))') r(10),ibisectr,ibisecttotr
    !write (6,'(''bisection rend level 0 to mmaxnow ratio ='',t40,f10.5)') (bisectcountr(k)/bisecttotr(k),k=0,mmax-1)
    write (6,'(''# of energy samples ='',t40,i10)') ice
    write (6,*) 
   
    write (12,'(''acceptance ratio (mov1) ='',t40,f10.5,3(i10))') r(1),icmov1,icmov1tot,icstdmove1tot
    write (12,'(''acceptance ratio (mov2) ='',t40,f10.5,2(i10))') r(2),icstdmove2,icstdmove2tot
    write (12,'(''acceptance ratio (mov3 shift) ='',t40,f10.5,2(i10))') r(3),icshft,icshfttot
    write (12,'(''ic23 ratio (stdmov23) ='',t40,f10.5,2(i10))') r(4),ic23,ic23tot
    write (12,'(''icn ratio ='',t40,f10.5,2(i10))') r(5),icn,icntot
    write (12,'(''lr ratio ='',t40,f10.5,2(i10))') r(6),iclr,iclrtot
    write (12,'(''reptation ratio ='',t40,f10.5,3(i10))') r(7),icrep,icreptot
    write (12,'(''ibisect ratio ='',t40,f10.5,3(i10))') r(8),ibisect,ibisecttot,icbisecttot
    write (12,'(''bisection level 1 to mmax ratio ='',t40,f10.5)') (bisectrate(k),k=1,mmax)
    write (12,'(''ibisect lend ratio ='',t40,f10.5,3(i10))') r(9),ibisectl,ibisecttotl
    !write (12,'(''bisection lend level 0 to mmaxnow ratio ='',t40,f10.5)') (bisectcountl(k)/bisecttotl(k),k=0,mmax-1)
    write (12,'(''ibisect rend ratio ='',t40,f10.5,3(i10))') r(10),ibisectr,ibisecttotr
    !write (12,'(''bisection rend level 0 to mmaxnow ratio ='',t40,f10.5)') (bisectcountr(k)/bisecttotr(k),k=0,mmax-1)
    write (12,'(''# of energy samples ='',t40,i10)') ice
    write (12,*) 	  
    
    return
    end subroutine showstat 
   
    subroutine bisectpickslices(ileftbisect,irightbisect,ipick) 
    use random
    integer(kind=i4) :: ileftbisect,irightbisect,ipick,na,nb
    real(kind=r8) :: rn(1)
    rn=randn(1)
    ! need to initialize nbisect and nchorizo first.
    na=nbisect/2
    nb=nchorizo-na
    ipick=na+int(dble(nb-na+1)*rn(1)) ! I believe this picking strategy is good. 
    ileftbisect=ipick-na
    irightbisect=ipick+na 
    return    
    end subroutine bisectpickslices
   
    subroutine bisectpicksliceswide(ileftbisect,irightbisect) 
    use random
    integer(kind=i4) :: ileftbisect,irightbisect
    real(kind=r8) :: rn(1)
    rn=randn(1)
    ! need to initialize nbisect and nchorizo first.
    ileftbisect=min(int(dble(nchorizo+1)*rn(1)),nchorizo) ! min just to make sure, not really need.
    irightbisect=int(dble(nchorizo+nbisect-1)*rn(1))+1
    return    
    end subroutine bisectpicksliceswide    
      
 
   
    subroutine xoutput(x) ! initialize the position, calculate the weights for each bead. l2r,r2l
    real(kind=r8) :: x(3,npart)
    x=ristra(0)%x
    return
    end subroutine xoutput   
   
   
    subroutine cwtoutr(i,cwttmp)
    use wavefunction
    use estimator
    use random
    use v6stepcalc
    integer(kind=i4) :: i
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin),cwt(0:nspin-1,nisospin)
    call getcwtgnd2(cwt)
    call v6propr(ristra(0)%x,i,cwt,cwttmp)
    return
    end subroutine
   
    subroutine cwtoutl(i,cwttmp)
    use wavefunction
    use estimator
    use random
    use v6stepcalc
    integer(kind=i4) :: i
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin),cwt(0:nspin-1,nisospin)
   
    call getcwtgnd2(cwt)
    call v6propl(ristra(0)%x,i,cwt,cwttmp)
    return
    end subroutine   
   
    subroutine cwtr(k,i,cwttmp,cwt)
    use v6stepcalc
    integer(kind=i4) :: i,k
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin),cwt(0:nspin-1,nisospin)
   
    call v6propr(ristra(k)%x,i,cwttmp,cwt)    
    return
    end subroutine 
   
    subroutine cwtrpr(k,i,cwttmp,cwt)
    use v6stepcalc
    integer(kind=i4) :: i,k
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin),cwt(0:nspin-1,nisospin)
   
    call v6proprpr(ristra(k)%x,i,cwttmp,cwt)    
    return
    end subroutine     
   
   
    subroutine cwtl(k,i,cwttmp,cwt)
    use v6stepcalc
    integer(kind=i4) :: i,k
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin),cwt(0:nspin-1,nisospin)
   
    call v6propl(ristra(k)%x,i,cwttmp,cwt)    
    return
    end subroutine   
   
    subroutine cwtlpr(k,i,cwttmp,cwt) ! pr means use -dt.
    use v6stepcalc
    integer(kind=i4) :: i,k
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin),cwt(0:nspin-1,nisospin)
   
    call v6proplpr(ristra(k)%x,i,cwttmp,cwt)    
    return
    end subroutine    
   
    subroutine cwt1test(ntab,rangev,it,i,j,cwtold,cwtnew,val,u,evdt) ! just for a particular pair i,j
    use v6stepcalc
    integer(kind=i4) :: i,j,it,index,ntab	
    real(kind=r8) :: x(3,npart) ! the configuration for the system
    real(kind=r8) :: dx(3)
    real(kind=r8) :: r,c1,c2,c3,c4,dr,val,e1,e5! tstep=dt*2**it
    real(kind=r8) :: u(6),evdt(6),scalep,rangev,u1(6),u2(6),evdtinv(6)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin),cwt1(0:nspin-1,nisospin)
    real(kind=r8) :: utab(6,0:ntab,-mmax:mmax),evdttab(6,0:ntab,-mmax:mmax) 
    real(kind=r8) :: vtab(6,0:ntab)

    

    call getutab(utab)
    
    
    call getevdttab(evdttab)  
   

    x=ristra(0)%x
   
    cwt1=cwtold
	
    !write(6,*) 'x=',x

            dx(:)=x(:,i)-x(:,j)		   
		    r=sqrt(dot_product(dx,dx)) !r=sqrt(sum(dx(:)**2))
		   
		    if (r.eq.0) then
		        write(6,*) 'two particles have the same position, stop'
			    stop
		    endif
	   
		    scalep=ntab/rangev
		   
            if (r.lt.rangev) then
			    dr=scalep*r
			    index=dr
			    index=max(1,min(index,ntab-2))
			    dr=dr-index
                c1=-dr*(dr-1.0_r8)*(dr-2.0_r8)/6.0_r8
                c2=(dr+1.0_r8)*(dr-1.0_r8)*(dr-2.0_r8)/2.0_r8
                c3=-(dr+1.0_r8)*dr*(dr-2.0_r8)/2.0_r8
                c4=(dr+1.0_r8)*dr*(dr-1.0_r8)/6.0_r8
                u(:)=c1*utab(:,index-1,it)+c2*utab(:,index,it) &
			        +c3*utab(:,index+1,it)+c4*utab(:,index+2,it)
			    evdt(:)=c1*evdttab(:,index-1,it)+c2*evdttab(:,index,it) &
			        +c3*evdttab(:,index+1,it)+c4*evdttab(:,index+2,it)
			
              
                write(6,'(''index r dr='',t50, i20, 2(f15.7,1x) )') index,r,dr
                write(6,'(''c1-4 ='',t50,4(e27.20,1x) )') c1,c2,c3,c4
              
                call vtableout(vtab)
              
                write(6,'(''vtab='',t50,6(e27.20,1x) )') vtab(:,index)
                write(6,'(''evdt='',t50,6(e27.20,1x) )') evdt
                write(6,'(''1/evdt='',t50,6(e27.20,1x) )') 1./evdt 
			    write(6,'(''lagrange u='',t50,6(e27.20,1x) )') u
              
        u1(1)=(6.*evdt(1)+3.*evdt(2)+2.*evdt(3)+evdt(4)+3.*evdt(5)+evdt(6))/16.
	    u1(2)=(6.*evdt(1)+3.*evdt(2)+2.*evdt(3)+evdt(4)-9.*evdt(5)-3.*evdt(6))/48.
	    u1(3)=(3.*evdt(1)-3.*evdt(2)+evdt(3)-evdt(4))/24.
	    u1(4)=(2.*evdt(1)+evdt(2)-2.*evdt(3)-evdt(4)+evdt(5)-evdt(6))/16.
	    u1(5)=(2.*evdt(1)+evdt(2)-2.*evdt(3)-evdt(4)-3.*evdt(5)+3.*evdt(6))/48.
	    u1(6)=(evdt(1)-evdt(2)-evdt(3)+evdt(4))/24.
       
        evdtinv=1./evdt 
        u2(1)=(6.*evdtinv(1)+3.*evdtinv(2)+2.*evdtinv(3)+evdtinv(4)+3.*evdtinv(5)+evdtinv(6))/16.
	    u2(2)=(6.*evdtinv(1)+3.*evdtinv(2)+2.*evdtinv(3)+evdtinv(4)-9.*evdtinv(5)-3.*evdtinv(6))/48.
	    u2(3)=(3.*evdtinv(1)-3.*evdtinv(2)+evdtinv(3)-evdtinv(4))/24.
	    u2(4)=(2.*evdtinv(1)+evdtinv(2)-2.*evdtinv(3)-evdtinv(4)+evdtinv(5)-evdtinv(6))/16.
	    u2(5)=(2.*evdtinv(1)+evdtinv(2)-2.*evdtinv(3)-evdtinv(4)-3.*evdtinv(5)+3.*evdtinv(6))/48.
	    u2(6)=(evdtinv(1)-evdtinv(2)-evdtinv(3)+evdtinv(4))/24.       
              
                write(6,'(''assume evdt is correct true u='',t50,6(e27.20,1x) )') u1
                write(6,'(''difference btw true u and u is='',t50,6(e27.20,1x) )') u1-u
                write(6,'(''true u for -dt the u2 ='',t50,6(e27.20,1x) )') u2
              
               
               
                e1=evdt(1)
                e5=evdt(5)
                val=(3./e1+1./e5)*(3.*e1+e5)/16.+(3./e1-3./e5)*(e1-e5)/16.
              
                write(6,*) 'check 1= ',val
               
              
		    else
			    write(6,*) 'Damn!!!!!!!!!!!!!!!!!!' 
			    u=0.
                u1=0.
			    stop
		    endif
		    call v6prop1(x,i,j,u1,cwt1,cwtnew)
    return
    end subroutine cwt1test    
   
    subroutine vtot(i,v) ! v is the numerator, in the end need the denominator.
    use brut
    use v6stepcalc   
    use chorizos
    complex(kind=r8) :: cwtr2lm(0:nchorizo,0:nspin-1,nisospin),cwtl2rm(0:nchorizo,0:nspin-1,nisospin) ! m means middle.
    real(kind=r8) :: v
    integer(kind=i4) :: i
    real(kind=r8) :: x(3,npart) 

    x=ristra(i)%x 
    call xconvert(x)  
    if (i.eq.0) then             
        call pecal(x,cwtbegin,cwtr2l(i,:,:),v)     
    else
        if (i.eq.nchorizo) then           
            call pecal(x,cwtl2r(i,:,:),cwtend,v)   
        else
    ! in the somewhere middle  
        call v6proplpr(x,-1,cwtr2l(i,:,:),cwtr2lm(i,:,:))
        call v6proplpr(x,-1,cwtl2r(i,:,:),cwtl2rm(i,:,:))
        call pecal(x,cwtl2rm(i,:,:),cwtr2lm(i,:,:),v)        
        endif   
    endif
   
    return
    end subroutine vtot

    subroutine hpsi(c0,cn,valnum,valdenom,valnumke,valnumpev6,valnumpeem,psi20,psi2n,rf,icheck) ! calculate numerator and denominator.
    use brut
    use v6stepcalc   
    use chorizos
    use mympi
    type (chorizo) :: c0,cn
    complex(kind=r8) ::  cwta1(25,0:nspin-1,nisospin),cwta2(25,0:nspin-1,nisospin) 
    complex(kind=r8) :: f
    real(kind=r8) :: rf,absrf,valnum,valdenom,alr,blr,valnumpev6,valnumpeem,valnumke
    real(kind=r8) :: psi20,psi2n,e0,en,psi2,pe0,pen,empe0,empen,psi,ke0,ken
    integer(kind=i4) :: icheck ! 0=pass,1=fail.
    real(kind=r8) :: x0(3,npart),xn(3,npart)
 
    !cwtr2lout=cwtr2l
    !cwtl2rout=cwtl2r
    x0=c0%x
    xn=cn%x
   
    !write(6,*) 'hpsi start!'
     
    ! add the cwt r2l l2r checking.

    call hpsitcwt(x0,iplo,jplo,cwtr2l(0,:,:),psi20,e0,ke0,pe0,empe0,cwta1)
     
    !call hpsitcwt(x0,iplo,jplo,cwtr2l(0,:,:),psi20,e0,d20,pe0,f0,cwta0) 
   
    !write(6,*) 'cwta1 - cwta0 =', sum(abs(cwta1(1,:,:)-cwta0(1,:,:)))!,(cwta1(1,:,:)-cwta0(1,:,:))   
  
    !call corop(x0,iplo,jplo,cwt)   
    !write(6,*) 'cwta1 - cwtbegin =', sum(abs(cwta1(1,:,:)-cwt(:,:)))!,(cwta1(1,:,:)-cwt(:,:))
    !write(6,*) 'cwta1 =', cwta1(1,:,:)
    !write(6,*) 'cwtbegin =', cwt
 
    call hpsitcwt(xn,ipro,jpro,cwtl2r(nchorizo,:,:),psi2n,en,ken,pen,empen,cwta2)
    !call corop(xn,ipro,jpro,cwt)  
    !write(6,*) 'cwta2 - cwtbegin =', sum(abs(cwta2(1,:,:)-cwt(:,:)))!,(cwta2(1,:,:)-cwt(:,:))       
    !write(6,*) 'cwta2 =', cwta2(1,:,:)
    !write(6,*) 'cwtend =', cwt  
   
    psi2=psi20+psi2n
    psi=psi2/2  ! psi is real f. ! psi20 and psi2n should be the same.
   
    alr=0.5*(e0+en)/abs(psi)
    blr=psi/abs(psi)
   
    
    !open(unit=19,form='formatted',file='check1.txt')
    !write(6,'( ''hpsitcwt check'', 8(g15.7,1x) )') e0,en,psi20,psi2n
    !write(19,'( ''hpsitcwt check'', 8(g15.7,1x) )') e0,en,psi20,psi2n
    !
    !do i=1,npart-1
    !    do j=i+1,npart
    ! 
    !    !write(6,'( ''i,j='', 2(i10,1x) )') i,j    
    !    write(19,'( ''i,j='', 2(i10,1x) )') i,j   
    ! 
	   ! dx(:)=x0(:,i)-x0(:,j)
	   ! r(:)=dx(:)/sqrt(dot_product(dx,dx)) !unit vector   
    !
    !    !write(6,'( ''rij='', 3(g15.7,1x) )') sqrt(dot_product(dx,dx)) 
    !    !write(6,'( ''x0 check'', 3(g15.7,1x) )') r
    !    write(19,'( ''rij='', 3(g15.7,1x) )') sqrt(dot_product(dx,dx)) 
    !    write(19,'( ''x0 check'', 3(g15.7,1x) )') r
    !    write(19,'( ''x0 check x0'', 3(g15.7,1x) )') x0
    ! 
    ! 
	   ! dx(:)=xn(:,i)-xn(:,j)
	   ! r(:)=dx(:)/sqrt(dot_product(dx,dx)) !unit vector   
    !
    !    !write(6,'( ''xn check'', 3(g15.7,1x) )') r   
    !    write(19,'( ''xn check'', 3(g15.7,1x) )') r   
    !    write(19,'( ''xn check xn'', 3(g15.7,1x) )') xn     
    ! 
    !enddo
    !enddo
    !
    !
    !do i=0,nchorizo
    !!write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'   
    !!write(6,*) 'check cwtr2l,l2r'
    !!write(6,*) 'i=',i
    !!write(6,'( ''cwtl2r ='', 12(g15.7,1x) )') cwtl2r(i,:,:) 
    !!write(6,*) 'nchorizo-i=',nchorizo-i    
    !!write(6,'( ''cwtr2l ='', 12(g15.7,1x) )') cwtr2l(nchorizo-i,:,:)
    !!write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    !
    !
    !write(19,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'   
    !write(19,*) 'check cwtr2l,l2r'
    !write(19,*) 'i=',i
    !write(19,'( ''cwtl2r ='', 12(g15.7,1x) )') cwtl2r(i,:,:) 
    !write(19,*) 'nchorizo-i=',nchorizo-i    
    !write(19,'( ''cwtr2l ='', 12(g15.7,1x) )') cwtr2l(nchorizo-i,:,:)
    !write(19,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    !
    !enddo
    !
    !
    !close(19)
   
    valnum=alr
    valdenom=blr
   
    valnumpev6=0.5*(pe0+pen)/abs(psi)
    valnumpeem=0.5*(empe0+empen)/abs(psi)
    valnumke=alr-valnumpev6-valnumpeem 
   
    !write(6,*) 'e0=',e0
    !write(6,*) 'en=',en,(e0+en)
    !write(6,*) 'psi20, psi2n=',psi20,psi2n,psi,abs(psi)
    !write(6,*) 'valnum=',valnum,valnumke,alr,valnumpev6,valnumpeem
    !write(6,*) 'valdenom=',valdenom

    !rgl=e0
    !rgr=en
    !rg=(rgl+rgr)/2.
    !
    f=cmplx(psi) ! neglect the imaginary part actually.
    rf=real(f)
    absrf=abs(real(f))
  
    !write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'  
    !write(6,*) 'psi20 and psi2n=',psi20,psi2n  
    !write(12,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'  
    !write(12,*) 'psi20 and psi2n=',psi20,psi2n   
   
   
   
    if ( abs((psi20-psi2n)/((psi20+psi2n)/2.)).gt.(1.d-3)   ) then
       
    icheck=1 ! fail.
       
    !write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'  
    !write(6,*) 'somthing wrong with psi20 and psi2n=',psi20,psi2n  
    !write(12,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'  
    !write(12,*) 'somthing wrong with psi20 and psi2n=',psi20,psi2n  	
    !write(12,*) 'iplo=',iplo
    !write(12,*) 'jplo=',jplo   
    !write(12,*) 'ipro=',ipro
    !write(12,*) 'jpro=',jpro      
    !write(12,*) 'cwtl2r(nchorizo,:,:)=',cwtl2r(nchorizo,:,:)     
    !write(12,*) 'cwtr2l(0,:,:)=',cwtr2l(0,:,:)

    ! scan and compare.
    
    ! write out x.      
	    !open(unit=9,form='formatted',file=trim(outfile))
    !   rewind 9
    !   do l=0,nchorizo
	    ! call chorizoout(x,l)
	    ! write (9,'(6e15.7)') x
    !   enddo
    !   write (9,'(6i10)') iplo 
    !   write (9,'(6i10)') jplo 
    !   write (9,'(6i10)') ipro 
    !   write (9,'(6i10)') jpro 
    !   close(9)  
    !   write(6,*) 'bug x recorded in outfile! use hpsi to check!'
    !  stop
    
    
    else
    icheck=0 ! pass
    endif
   
    ! calculate V at all the beads.

    return
    end subroutine hpsi     
   
    subroutine stepoutputordlr(ipl,jpl,ipr,jpr)
    integer(kind=i4) :: ipl(6),jpl(6),ipr(6),jpr(6)	
	    ipl=iplo
	    jpl=jplo
	    ipr=ipro
	    jpr=jpro	
    return
    end subroutine stepoutputordlr 
   
    
    subroutine cleanfolder(irepin)
    integer(kind=i4) :: irepin
    ! clean those file with position=append
    
    select case (irepin)     
    case (2)
        
        open(unit=19,form='unformatted',file='path.unf',status='replace',position='rewind') 
        close(19)
        open(unit=19,form='unformatted',file='path0.unf',status='replace',position='rewind') 
        close(19)
        open(unit=19,form='formatted',file='reparrowhist.txt',status='replace',position='rewind')
        close(19)
        open(unit=11,form='formatted',file='values.txt',status='replace',position='rewind')
        close(11)  
        open(unit=13,form='formatted',file='energy_evolution.txt',status='replace',position='rewind')
        close(13)          
        
    case (4,5:6)
        
        open(unit=19,form='formatted',file='reparrowhist.txt',status='replace',position='rewind')
        close(19)
        open(unit=11,form='formatted',file='values.txt',status='replace',position='rewind')
        close(11)  
        open(unit=13,form='formatted',file='energy_evolution.txt',status='replace',position='rewind')
        close(13)  
        
    case default
        write (6,'(''irep value has not been decided in subroutine cleanfolder yet '')') irep
        call abort         
           
    end select

    return
    end subroutine cleanfolder
    
    subroutine updatecwtchain ! update the whole cwtr2l, cwtl2r, cwtbegin and cwtend.
    use wavefunction
    use v6stepcalc
    use brut 
    use random
    integer(kind=i4) :: i,j
    real(kind=r8) :: xr(3,npart),xl(3,npart),x(3,npart)
   
    xl=ristra(0)%x 
    xr=ristra(nchorizo)%x
    ! l2r   
    call corop(xl,iplo,jplo,cwtbegin) ! psitcwt
    call v6propr(xl,-1,cwtbegin,cwtl2r(0,:,:))
    do i=1,nchorizo-1
    x=ristra(i)%x
    call v6proplr(x,x,-1,cwtl2r(i-1,:,:),cwtl2r(i,:,:))
    enddo   
    call v6propl(xr,-1,cwtl2r(nchorizo-1,:,:),cwtl2r(nchorizo,:,:))
    ! r2l 
    call corop(xr,ipro,jpro,cwtend)
    call v6propr(xr,-1,cwtend,cwtr2l(nchorizo,:,:))
    do j=nchorizo-1,1,-1
    x=ristra(j)%x
    call v6proplr(x,x,-1,cwtr2l(j+1,:,:),cwtr2l(j,:,:))
    enddo
    call v6propl(xl,-1,cwtr2l(1,:,:),cwtr2l(0,:,:))

    return
    end subroutine updatecwtchain   

    
    subroutine compute ! call stepinit first so that get irep first.
    use estimator
    use estimatorristra
    use wavefunction
    use random
    use brut
    use mympi
    use math
    real(kind=r8) :: time0,time1,timeall,second
    real(kind=r8) :: x0tot(3,npart,0:nchorizo),x0tot1d(3*npart*(nchorizo+1))
    integer(kind=i4) :: i,j,k,k0,l,ixtemp,nmax,nmax1,nmax2
    integer(kind=i4) :: np,n1,m1,n2,m2,n3,n4,n4c,m4,n5,n51,n5c,n6,n61,n6c,lnow &
                      ,npatheach0,npatheach &
                      ,day,hour,minute
    integer(kind=i4) :: ipathtot,ipathnow
    integer(kind=i4), allocatable :: iplall(:,:),jplall(:,:),iprall(:,:),jprall(:,:)	! for mpi_gather
    real(kind=r8), allocatable :: xall(:,:) ! for mpi_gather
    integer(kind=i4), allocatable :: iplallo(:),jplallo(:),iprallo(:),jprallo(:)	
    real(kind=r8), allocatable :: xallo(:) 
    integer(kind=i4), allocatable :: iplalln(:),jplalln(:),ipralln(:),jpralln(:)	
    real(kind=r8), allocatable :: xalln(:)  
    logical :: loaded,loadmore,loadextra,loadedextra
    
    ixtemp=3*npart*(nchorizo+1)    
    if (ixtemp /= ixtempo) then
       if (myrank().eq.0) then     
           write (6,*) 'ixtemp and ixtemp old does not match! check!'
           call abort
       endif  
    endif
            
    ipathtot=ipath*npo
    np=nproc()        
    n1=npo/np
    m1=mod(npo,np)
        
    n2=ipathtot/np
    m2=mod(ipathtot,np)
        
    npatheach0=10 ! it can change to a proper value according to the speed.
    npatheach=min(npatheach0,n2) ! # of paths per core for each step of calculation.
    n3=npatheach/npatheach0
        
    n4=np/npo 
    m4=mod(np,npo)

    loaded=.false.
    loadmore=.false.
    loadextra=.false.
    loadedextra=.false.        
       
    select case (irep)
           
    case (5)  ! the usual way, read in a path then calcualte. 
        
        if (myrank().eq.0) then   
            allocate(xallo(ixtemp*npo))
            allocate(iplallo(6*npo),jplallo(6*npo),iprallo(6*npo),jprallo(6*npo))   
            allocate(xalln(ixtemp*np))
            allocate(iplalln(6*np),jplalln(6*np),ipralln(6*np),jpralln(6*np))   
            open(unit=19,form='unformatted',file='path.unf',position='rewind')   
        endif                             
        if (n4 == 0) then ! npo # > np         
           if (m1>0) n6=np/m1  ! np<npo             
           if (myrank().eq.0) then   
                if (n6>1) then       
                    allocate(xall( (ixtemp*n1*np+1):ixtemp*npo,n6))
                    allocate(iplall((6*n1*np+1):6*npo,n6),jplall((6*n1*np+1):6*npo,n6) &
                            ,iprall((6*n1*np+1):6*npo,n6),jprall((6*n1*np+1):6*npo,n6))                         
                endif
           endif                   
            n6c=0          
            do i=1,ipath  
                if (myrank()==0) then  ! read in the old path file.  
                    if (i==1) time0=mpi_wtime()
                    read (19) ipathnow
                 do l=0,npo-1    
                    read (19) lnow
                    read (19) iplallo((6*l+1):(6*(l+1))) 
                    read (19) jplallo((6*l+1):(6*(l+1))) 
                    read (19) iprallo((6*l+1):(6*(l+1))) 
                    read (19) jprallo((6*l+1):(6*(l+1))) 
                    read (19) xallo((ixtemp*l+1):(ixtemp*(l+1)))               
                 enddo                          
                endif  
                do j=1,n1   
                    if (myrank().eq.0) then
                        do k=0,np-1                
                            l=np*(j-1)+k ! convert the current core # from new to old. 
                            if (k /= 0) then
                                call send(iplallo((6*l+1):(6*(l+1))),k,1) 
                                call send(jplallo((6*l+1):(6*(l+1))),k,2)
                                call send(iprallo((6*l+1):(6*(l+1))),k,3) 
                                call send(jprallo((6*l+1):(6*(l+1))),k,4) 
                                call send(xallo((ixtemp*l+1):(ixtemp*(l+1))),k,5)      
                            else ! send to process 0.
                                iplo=iplallo((6*l+1):(6*(l+1)))
                                jplo=jplallo((6*l+1):(6*(l+1)))
                                ipro=iprallo((6*l+1):(6*(l+1)))
                                jpro=jprallo((6*l+1):(6*(l+1)))
                                x0tot1d=xallo((ixtemp*l+1):(ixtemp*(l+1)))
                            endif          
                        enddo
                    else                                        
                        call recv(iplo,0,1)
                        call recv(jplo,0,2)
                        call recv(ipro,0,3)
                        call recv(jpro,0,4)
                        call recv(x0tot1d,0,5)                  
                    endif            
                    ! calculate 
                    x0tot=reshape(x0tot1d,shape(x0tot)) 
                    call computecore(x0tot,iplo,jplo,ipro,jpro,0,np-1)    
                    ! update                       
                    call barrier
                    call computeupdate(i) ! require all the cores.   
                    if ((j.eq.n1).and.(myrank().eq.0).and.(i==1)) then
                        time1=mpi_wtime()
                        timeall=(time1-time0)*dble(ipathtot)/dble(n1*np)
                        call timedhms(timeall,day,hour,minute,second)                        
    write (6,'(/,''Estimated total time ='',i10,'' days'',i10,'' hours'',i10,'' minutes'',f10.3,'' seconds'')') &
           day,hour,minute,second
    write (12,'(/,''Estimated total time ='',i10,'' days'',i10,'' hours'',i10,'' minutes'',f10.3,'' seconds'')') &
           day,hour,minute,second         
                    endif  
                enddo   
                if (m1 > 0) then  !        deal with the rest and calculate                              
                 if (n6 <= 1) then                   
                    if ( myrank()<=m1-1 ) then  
                        if (myrank().eq.0) then            
                            do k=0,m1-1             
                                l=np*n1+k
                                if (k /= 0) then
                                    call send(iplallo((6*l+1):(6*(l+1))),k,1) 
                                    call send(jplallo((6*l+1):(6*(l+1))),k,2)
                                    call send(iprallo((6*l+1):(6*(l+1))),k,3) 
                                    call send(jprallo((6*l+1):(6*(l+1))),k,4) 
                                    call send(xallo((ixtemp*l+1):(ixtemp*(l+1))),k,5) 
                                else
                                    iplo=iplallo((6*l+1):(6*(l+1)))
                                    jplo=jplallo((6*l+1):(6*(l+1)))
                                    ipro=iprallo((6*l+1):(6*(l+1)))
                                    jpro=jprallo((6*l+1):(6*(l+1)))
                                    x0tot1d=xallo((ixtemp*l+1):(ixtemp*(l+1)))
                                endif        
                            enddo                                                         
                        else                    
                            call recv(iplo,0,1)
                            call recv(jplo,0,2)
                            call recv(ipro,0,3)
                            call recv(jpro,0,4)
                            call recv(x0tot1d,0,5)                    
                        endif                   
                    endif  
                    if ( myrank()<=m1-1 ) then                  
                        x0tot=reshape(x0tot1d,shape(x0tot))             
                        call computecore(x0tot,iplo,jplo,ipro,jpro,0,m1-1) 
                    endif
                    call barrier 
                    call computeupdate(i)                     
                 else ! n6 >1             
                     n6c=n6c+1                    
                     if (myrank().eq.0) then                         
                        do l=n1*np,npo-1 
                            iplall((6*l+1):(6*(l+1)),n6c)=iplallo((6*l+1):(6*(l+1)))
                            jplall((6*l+1):(6*(l+1)),n6c)=jplallo((6*l+1):(6*(l+1)))
                            iprall((6*l+1):(6*(l+1)),n6c)=iprallo((6*l+1):(6*(l+1)))
                            jprall((6*l+1):(6*(l+1)),n6c)=jprallo((6*l+1):(6*(l+1)))
                            xall((ixtemp*l+1):(ixtemp*(l+1)),n6c)=xallo((ixtemp*l+1):(ixtemp*(l+1)))                        
                        enddo                  
                     endif                   
                     if ( (n6c.eq.n6).or.((i.eq.ipath).and.(n6c<n6).and.(n6c>0)) ) then                        
                         if (myrank()<=m1*n6c-1) then                            
                             if (myrank().eq.0) then                                
                                 do k=0,m1*n6c-1                   
                                   n61=k/m1+1    
                                   l=n1*np+mod(k,m1)
                                   if (k/=0) then
                                       call send( iplall((6*l+1):(6*(l+1)),n61),k,1   )
                                       call send( jplall((6*l+1):(6*(l+1)),n61),k,2   )
                                       call send( iprall((6*l+1):(6*(l+1)),n61),k,3   )
                                       call send( jprall((6*l+1):(6*(l+1)),n61),k,4   )
                                       call send( xall((ixtemp*l+1):(ixtemp*(l+1)),n61),k,5 )
                                   else
                                       iplo=iplall((6*l+1):(6*(l+1)),n61)
                                       jplo=jplall((6*l+1):(6*(l+1)),n61)
                                       ipro=iprall((6*l+1):(6*(l+1)),n61)
                                       jpro=jprall((6*l+1):(6*(l+1)),n61)
                                       x0tot1d=xall((ixtemp*l+1):(ixtemp*(l+1)),n61)
                                   endif                                     
                                 enddo
                             else                
                                call recv(iplo,0,1)
                                call recv(jplo,0,2)
                                call recv(ipro,0,3)
                                call recv(jpro,0,4)
                                call recv(x0tot1d,0,5)                
                             endif                   
                          x0tot=reshape(x0tot1d,shape(x0tot))             
                          call computecore(x0tot,iplo,jplo,ipro,jpro,0,m1*n6c-1)                        
                         endif   
                       call barrier
                       call computeupdate(i) ! require all the cores.                                            
                         n6c=0                                      
                     endif                 
                 endif                                
                endif                
            enddo                      
        else ! np >= npo                 
            n5=np/(npo-m4)       ! save the rest ones.     
            if (myrank().eq.0) then     
                allocate(xall((ixtemp*m4+1):ixtemp*npo,n5))
                allocate(iplall((6*m4+1):6*npo,n5),jplall((6*m4+1):6*npo,n5) &
                        ,iprall((6*m4+1):6*npo,n5),jprall((6*m4+1):6*npo,n5))         
            endif        
            n4c=0 ! count n4
            n5c=0   
            do i=1,ipath
                if (myrank()==0) then
                    time0=mpi_wtime()
                    read (19) ipathnow
                    do l=0,npo-1    
                        read (19) lnow
                        read (19) iplallo((6*l+1):(6*(l+1))) 
                        read (19) jplallo((6*l+1):(6*(l+1))) 
                        read (19) iprallo((6*l+1):(6*(l+1))) 
                        read (19) jprallo((6*l+1):(6*(l+1))) 
                        read (19) xallo((ixtemp*l+1):(ixtemp*(l+1)))               
                    enddo  
                endif                                     
                n4c=n4c+1 ! n4                   
                if (npo*n4c>np) then
                   loadextra =.true. 
                else     
                   loadextra =.false. 
                endif  
                if (loadextra) then
                    nmax1=m4-1
                else
                    nmax1=npo-1
                endif
                k0=npo*(n4c-1) 
                nmax=k0+nmax1            
                if ( myrank().eq.0 ) then        
                    do l=0,nmax1
                        k=k0+l
                        if (k/=0) then
                            call send(iplallo((6*l+1):(6*(l+1))),k,1) 
                            call send(jplallo((6*l+1):(6*(l+1))),k,2)
                            call send(iprallo((6*l+1):(6*(l+1))),k,3) 
                            call send(jprallo((6*l+1):(6*(l+1))),k,4) 
                            call send(xallo((ixtemp*l+1):(ixtemp*(l+1))),k,5)
                        else ! send to rank0 itself.
                            iplo=iplallo((6*l+1):(6*(l+1)))
                            jplo=jplallo((6*l+1):(6*(l+1)))
                            ipro=iprallo((6*l+1):(6*(l+1)))
                            jpro=jprallo((6*l+1):(6*(l+1)))
                            x0tot1d=xallo((ixtemp*l+1):(ixtemp*(l+1)))                        
                        endif 
                    enddo                                   
                    if ( loadextra ) then
                        n5c=n5c+1
                        do l=m4,npo-1
                            iplall((6*l+1):(6*(l+1)),n5c)=iplallo((6*l+1):(6*(l+1)))
                            jplall((6*l+1):(6*(l+1)),n5c)=jplallo((6*l+1):(6*(l+1)))
                            iprall((6*l+1):(6*(l+1)),n5c)=iprallo((6*l+1):(6*(l+1)))
                            jprall((6*l+1):(6*(l+1)),n5c)=jprallo((6*l+1):(6*(l+1)))
                            xall((ixtemp*l+1):(ixtemp*(l+1)),n5c)=xallo((ixtemp*l+1):(ixtemp*(l+1)))
                        enddo
                    endif     
                endif                             
                if ( (myrank() /= 0).and.(myrank() >= k0).and.( myrank() <= nmax ) ) then
                    call recv(iplo,0,1)
                    call recv(jplo,0,2)
                    call recv(ipro,0,3)
                    call recv(jpro,0,4)
                    call recv(x0tot1d,0,5)            
                endif                                      
               if (loadextra) loaded=.true.  
               call bcast(n5c)
               if (n5c.eq.n5) loadedextra=.true.                                 
                if (m4==0) then              
                   if (  (n4c==n4).or.(i==ipath)   ) then              
                    loaded=.true.
                   else                  
                    loaded=.false.                 
                   endif              
                else  ! load extra may needed. loadmore              
                   if (i<ipath) then 
                       if (.not.loadextra  ) then
                            loaded=.false.                              
                       else
                            loaded=.true.
                       endif                        
                   else                    
                       loaded=.true.                    
                   endif
                endif                                          
               if (loaded) then
                  if (myrank()<=nmax) then
                      x0tot=reshape(x0tot1d,shape(x0tot))             
                      call computecore(x0tot,iplo,jplo,ipro,jpro,0,nmax)       
                  endif
                  call barrier
                  call computeupdate(i) ! require all the cores.                                                 
                   ! reset
                   n4c=0
                   loaded=.false.                     
                   if (((i.eq.n4).or.(i.eq.n4+1)).and.myrank().eq.0) then
                        time1=mpi_wtime()
                        timeall=(time1-time0)*dble(ipathtot)/dble(np)
                        call timedhms(timeall,day,hour,minute,second)                        
    write (6,'(/,''Estimated total time ='',i10,'' days'',i10,'' hours'',i10,'' minutes'',f10.3,'' seconds'')') &
           day,hour,minute,second
    write (12,'(/,''Estimated total time ='',i10,'' days'',i10,'' hours'',i10,'' minutes'',f10.3,'' seconds'')') &
           day,hour,minute,second                                          
                   endif
                   if ( (loadedextra).or.((i.eq.ipath).and.(.not.loadedextra).and.(n5c/=0)) ) then                                                 
                       nmax2=n5c*(npo-m4)-1           
                       if (myrank().eq.0) then                 
                           do  k=0,nmax2
                               n51=k/(npo-m4)+1
                               l=m4+mod(k,npo-m4)  !m4+k-(n51-1)*(npo-m4)   ! check!!!!!!!!!!!!!!!!!!                                   
                               if (k/=0) then
                                   call send( iplall((6*l+1):(6*(l+1)),n51),k,1   )
                                   call send( jplall((6*l+1):(6*(l+1)),n51),k,2   )
                                   call send( iprall((6*l+1):(6*(l+1)),n51),k,3   )
                                   call send( jprall((6*l+1):(6*(l+1)),n51),k,4   )
                                   call send( xall((ixtemp*l+1):(ixtemp*(l+1)),n51),k,5 )
                               else
                                   iplo=iplall((6*l+1):(6*(l+1)),n51)
                                   jplo=jplall((6*l+1):(6*(l+1)),n51)
                                   ipro=iprall((6*l+1):(6*(l+1)),n51)
                                   jpro=jprall((6*l+1):(6*(l+1)),n51)
                                   x0tot1d=xall((ixtemp*l+1):(ixtemp*(l+1)),n51)
                               endif                                          
                           enddo
                       else if ( myrank() <= nmax2 ) then
                            call recv(iplo,0,1)
                            call recv(jplo,0,2)
                            call recv(ipro,0,3)
                            call recv(jpro,0,4)
                            call recv(x0tot1d,0,5)                    
                       endif
                       if (myrank() <= nmax2) then
                          x0tot=reshape(x0tot1d,shape(x0tot))             
                          call computecore(x0tot,iplo,jplo,ipro,jpro,0,nmax2)    
                       endif   
                       call barrier
                       call computeupdate(i) ! require all the cores.
                           ! reset
                           n5c=0
                           loadedextra=.false.      
                   endif                                                 
               endif                     
            call barrier
            enddo
        endif       
        call barrier
        if (myrank().eq.0) close (19)        
            
    case (6)
        
        ! aother mode, change the way the pimc output the path.unf,
        ! not decided and checked yet.
        ! or may besupermode, do not need to consider ram resource problem, big ram is ready. bug need to fix
               
          if (myrank().eq.0) then
              
              write(6,*) 'irep=6 case has not been decided yet, abort!'
              call abort
          endif
                 
    case default
        write (6,'(''Illegal irep value in compute, irep= '')') irep
        call abort      
        
    end select

    return
    end subroutine compute
    
    subroutine computecore(x0totin,iploin,jploin,iproin,jproin,rankl,rankr) 
! read in path, calculate all the results. and write to a file
! rank1 must be zero   
    use estimator
    use estimatorristra
    use wavefunction
    use random
    use brut
    use mympi
    real(kind=r8) :: rm,rmsq
    real(kind=r8) :: v(0:nchorizo) &
	                ,x0totin(3,npart,0:nchorizo)
    real(kind=r8), dimension(0:nrhobin) :: rhodistout
    integer(kind=i4) :: i,rankl,rankr,icheck
    real(kind=r8) :: vn,vd,vnke,vnpev6,vnpeem,rf,vnum,vdenom
    real(kind=r8) :: psi20,psi2n
    integer(kind=i4) :: iploin(6),jploin(6),iproin(6),jproin(6)
 
! init   
    call chorizoallin(x0totin)
    call inputordlr(iploin,jploin,iproin,jproin)   
    call updatecwtchain
    
! calculate    
    call hpsi(ristra(0),ristra(nchorizo),vn,vd,vnke,vnpev6,vnpeem,psi20,psi2n,rf,icheck)
    vnum=vn
    vdenom=vd
    ! update the v beads   
    do i=0,nchorizo
        call vtot(i,v(i))
    enddo
    v=v/rf
    call addvalristra(1,v) ! do not need this part   
   
    call addval(3,vnum,1.0_r8)
    call addval(4,vdenom,1.0_r8)
	call addval(5,rf,1.0_r8)    

	call funcrmsq(ristra(nchorizomid)%x,rm,rmsq)
	call addval(6,rmsq,1.0_r8)
	call addval(7,rm,1.0_r8)
    
	call samplerhodistrbt(ristra(nchorizomid)%x,nchorizomid,rhodistout)	  
	call addval(8,rhodistout(0),1.0_r8)
  
    ! write out all the num and denom of all the cores

    if (myrank().eq.0) open (unit=11,FILE='values.txt',FORM='FORMATTED',position='append')    
    do i=rankl,rankr ! rank1 must be zero     
        if ( myrank().eq.i  ) then         
            if (myrank() /= 0) then         
                call send(vnum,0,1)
                call send(vdenom,0,2)
                call send(vnke,0,3)
                call send(vnpev6,0,4)
                call send(vnpeem,0,5)
                call send(rm,0,6)
                call send(rhodistout(0),0,7)               
            else                
                ! rank0 write out.                   
                write(11,'(10(G15.7,2x))') vnum,vdenom,vnke,vnpev6,vnpeem,rm,rhodistout(0) !11 means write to 'values.txt'                          
            endif          
        endif
        if ( (myrank()==0).and.(i /= 0) ) then      
            call recv(vnum,i,1)
            call recv(vdenom,i,2)
            call recv(vnke,i,3)
            call recv(vnpev6,i,4)
            call recv(vnpeem,i,5)
            call recv(rm,i,6)
            call recv(rhodistout(0),i,7)         
      ! rank i write out.        
        write(11,'(10(G15.7,2x))') vnum,vdenom,vnke,vnpev6,vnpeem,rm,rhodistout(0)                     
        endif 
    enddo
    if (myrank().eq.0) close (11)
    return
    end subroutine computecore   
    
    subroutine computeupdate(iin) ! require all the cores.
    use estimator
    use estimatorristra
    use wavefunction
    use random
    use brut
    use mympi
    integer(kind=i4) :: j,k,iin
    real(kind=r8) :: val2,err2,val1,err1,error
    real(kind=r8), allocatable :: rhodist(:),rhodisterr(:)
    real(kind=r8), allocatable :: valnow(:),valaverage(:),valerr(:)
    character(len=30) :: vpropchorizo

    call update   !collect block averages
    call updateristra
      
    if (myrank().eq.0) then   
        write (6,*) 'path #: ', iin 
	    write (12,*) 'path #: ', iin  
        answer=resstring()
        write (6,'(a120)') (answer(k),k=1,size(answer))
	    write (12,'(a120)') (answer(k),k=1,size(answer))  
    endif    
     
    ! show rho distribution on the fly
	allocate(rhodist(0:nrhobin),rhodisterr(0:nrhobin))
	call updaterhodistrbt(rhodist,rhodisterr,nchorizomid) ! choose the middle bead
    if (myrank().eq.0) then
	    open(unit=32,FILE='rhodistribution.txt',FORM='FORMATTED')	 
	    do k=0,nrhobin
	    write (32,'(i10,1x,3(f10.5,1x))') k,(k+0.5)*rhobinsize,rhodist(k),rhodisterr(k)
	    enddo
	    close(32) 
    endif
    deallocate(rhodist,rhodisterr)
     
    ! write out a file with the Vtot values at different chorizos
    if (myrank().eq.0) then  
    allocate(valnow(0:nchorizo),valaverage(0:nchorizo),valerr(0:nchorizo))
    call resultristra(1,valnow,valaverage,valerr,vpropchorizo)
    open(unit=33,form='formatted',file=vpropchorizo)
    call resultest(4,val2,err2) ! #4 is the denominator 
    do j=0,nchorizo
        val1=valaverage(j)
        err1=valerr(j)  
        error=sqrt((1./val2*err1)**2+(val1/(val2**2)*err2)**2) 
        ! val1/val2 and error are the Vtot and error we want.
        write (33,'(i10,f10.5,6g15.6)') j,j*dt,val1/val2,error,valaverage(j),valerr(j),val2,err2
    enddo
    close(33) 
    deallocate(valnow,valaverage,valerr)
    endif     
    

    return
    end subroutine computeupdate  
    
    subroutine computeinit(answerin,rhobinsizein) ! call it before call compute
    use mympi
    character(len=120) :: answerin(:)
    real(kind=r8) :: rhobinsizein

    allocate(answer(size(answerin)))
    rhobinsize=rhobinsizein 
                                                                                   
    if (myrank().eq.0) then 
        open(unit=19,form='formatted',file='pathnumbercounts.txt',position='rewind')
        read(19,'(3i10)') npo,ipath,ixtempo
        close(19) 

        write (6,'(''npo ='',t30,i20)') npo
        write (6,'(''ipath ='',t30,i20)') ipath        
        
        write (12,'(''npo ='',t30,i20)') npo
        write (12,'(''ipath ='',t30,i20)') ipath
        
        if (ipath < 1) then
            write(6,*) 'path is empty, stop!'
            call abort
        endif  
    endif
    
    call bcast(npo)
    call bcast(ipath)  
    call bcast(ixtempo)
       
    call barrier
    return
    end subroutine computeinit      
    
    
    

!!!!!!!!!!!! this part is just to check v6 matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!    subroutine checkv6mat(x,psi2,psi2a)
!!  make sure this subroutine is called after some initializaiotn.
!    use matrixmod
!    use brut
!    use v6stepcalc
!    real(kind=r8) :: x(3,npart),x0(3,npart)
!    complex(kind=r8) :: cwt(0:nspin-1,nisospin),cwtnew(0:nspin-1,nisospin),cwtgnd(0:nspin-1,nisospin)
!    integer(kind=i4) ::  i,j,k,is,iiso 
!    complex(kind=r8) :: vij(nbasis,nbasis),vijeigvec(nbasis,nbasis),npsitr(nbasis),npsitl(nbasis) &
!                       ,evijdt(nbasis,nbasis),evijdt1(nbasis,nbasis),evijtt(nbasis,nbasis),evijtt1(nbasis,nbasis) &
!                       ,cwtbegin1(nbasis),cwtend1(nbasis),cwtend1t(nbasis,1),cwtbegin1t(1,nbasis) &
!                       ,c1(nbasis,nbasis)
!    real(kind=r8) :: elambdadt(nbasis,nbasis)
!    real(kind=r8) :: vijeigval(nbasis),psi2,tt,psi2a(1,1),psi2b,psi2tmp
!    
!    tt=dt*nchorizo    
!    
!    call v6mat(x,vij,vijeigvec,vijeigval)
! 
!    CALL PRINT_MATRIX( 'Vijmatrix', nbasis, nbasis, vij, nbasis ) 
!    CALL PRINT_MATRIX( 'Vijdagmatrix', nbasis, nbasis, conjg((transpose(vij))), nbasis ) 
!    CALL PRINT_RMATRIX( 'Eigenvalues', 1, nbasis, vijeigval, 1 )
!    CALL PRINT_MATRIX( 'Eigenvectors (stored columnwise)', nbasis, nbasis, vijeigvec, nbasis )  
!    CALL PRINT_RMATRIX( 'exp Eigenvalues totol time', 1, nbasis, exp(-vijeigval*tt), 1 ) 
!      
!!   get <n|psi_t(r)>  
!        
!    do i=1,nbasis    
!       cwtend1(i)=cwtend(invspin(i),invispin(i))
!       cwtbegin1(i)=cwtbegin(invspin(i),invispin(i))
!    enddo
!
!    do i=1,nbasis
!       npsitr(i)=sum(conjg(vijeigvec(:,i))*cwtend1(:))
!       npsitl(i)=sum(conjg(vijeigvec(:,i))*cwtbegin1(:)) 
!       !npsitl(i)=sum(conjg(cwtbegin1(:))*vijeigvec(:,i)) ! already conjugated
!    enddo
!
!    !CALL PRINT_CMATRIX( 'npsitr', 1, nbasis, npsitr, 1 )
!    !CALL PRINT_CMATRIX( 'npsitl', 1, nbasis, npsitl, 1 )
!    
!    psi2=sum(conjg(npsitl(:))*exp(-vijeigval(:)*tt)*npsitr(:)) 
!!psi2 is real part of f. if npsitl is already conjugated then no need to conjg it here again.
!    
!    !psi2tmp=0
!    !do i=1,nbasis
!    !   psi2tmp=psi2tmp+(npsitl(i))*exp(-vijeigval(i)*tt)*npsitr(i)  
!    !   write (6,'(10g15.7)') i,(npsitl(i)),exp(-vijeigval(i)*tt),npsitr(i) &
!    !                             ,(npsitl(i))*exp(-vijeigval(i)*tt)*npsitr(i)  
!    !   write (6,'(2g15.7)') i, psi2tmp
!    !enddo
!    
! 
!! matmul
!    
!    elambdadt=0
!    do i=1,nbasis
!       elambdadt(i,i)=exp(-dt*vijeigval(i))
!    enddo
!   
!    evijdt=matmul(vijeigvec,matmul(elambdadt,conjg(transpose(vijeigvec))))
!  
!    evijtt=0
!    do i=1,nbasis
!       evijtt(i,i)=1 
!    enddo 
!    do i=1,nchorizo    
!        evijtt=matmul(evijtt,evijdt)   
!    enddo
!
!    cwtbegin1t(1,:)=cwtbegin1(:)
!    cwtend1t(:,1)=cwtend1(:)
!    psi2a=matmul(conjg(cwtbegin1t),matmul(evijtt,cwtend1t))
!     
!! regular
!    
!    do j=1,nbasis
!        do i=1,nbasis
!           evijdt1(i,j)=sum(vijeigvec(i,:)*exp(-dt*vijeigval(:))*conjg(vijeigvec(j,:)))       
!        enddo  
!    enddo
!    
!    CALL PRINT_MATRIX( 'evdt', nbasis, nbasis, evijdt, nbasis ) 
!    call PRINT_MATRIX( 'evdt1', nbasis, nbasis, evijdt1, nbasis ) 
!      
!    evijtt1=0
!    do i=1,nbasis
!       evijtt1(i,i)=1 
!    enddo     
!    do i=1,nchorizo    
!        evijtt1=matmul(evijtt1,evijdt1)   
!    enddo   
!    
!    CALL PRINT_MATRIX( 'evtt', nbasis, nbasis, evijtt, nbasis ) 
!    call PRINT_MATRIX( 'evtt1', nbasis, nbasis, evijtt1, nbasis )     
!    
!    c1=matmul(vijeigvec,conjg(transpose(vijeigvec)))
!    CALL PRINT_MATRIX( 'nnh', nbasis, nbasis, c1, nbasis ) ! eigenvec*eigenvec^dagger=1
!    write(6,*) 'c1=', sum(c1)
!    
!    return
!    end subroutine checkv6mat 
!
!    subroutine sdv6(x)
!!  Method of steepest descent
!    use matrixmod
!    use brut
!    use v6stepcalc
!    real(kind=r8) :: x(3,npart),xo(3,npart),xn(3,npart),xrvs(3,npart) &
!                    ,lambda,t,error,errorx,g(3,npart),go(3,npart),gn(3,npart),f,fo,fn,ta,tb,tc &
!                    ,ga(3,npart),gb(3,npart),gc(3,npart) &
!                    ,step1,fa,fb,fc,step,step0,stepnow
!    real(kind=r8), allocatable :: trial1(:)
!    real(kind=r8), parameter :: small=1.0d-10,tiny=1.0d-10,epsilon=1.0d-8 &
!                               ,w=(3-sqrt(5.0_r8))/2,phi=(1+sqrt(5.0_r8))/2
!    integer(kind=i4) :: find1(1),pieces,i,j,k,l,lmax,ct,ia,ib,ic
!    logical :: abcdone,adone,bdone,cdone
!    
!    lmax=1000
!    step0=1.0d-5
!    step=step0
!    abcdone=.false.
!    adone=.true.
!    bdone=.false.
!    cdone=.false.
!    pieces=33
!    allocate(trial1(pieces))
!     
!    ta=0
!    xo=x
!
!    call ftest(xo,fo)
!    write(6,*) 'f initial =', fo
!    fa=fo
!    !write(6,*) 'x initial =', xo 
!    call gtest(xo,go)
!    write(6,*) 'g initial =', go 
!    g=go
!    
!    l=0
!    ct=20
!    open(unit=9,form='formatted',file='sdv6',position='rewind')
!    write (9,'(13g15.7)') xo,fo
!    
!    all: do while (ct>=0)
!        
!!!!!!!!!!!!!!!!!!!!!!!!! 
!        
!! for a given x, find the min t.           
!           ib=0
!           sub1: do while (.not.bdone)          
!               t=ta+step 
!               call ftabc(x,t,g,fn)
!               if (fn<fa) then
!                   tb=t
!                   fb=fn
!                   bdone=.true.
!                   stepnow=step
!               else
!                   step=step/10  ! do this to ensure ta=0, fa>fb.
!                   write(6,*) 'try b again',step,fn,fo
!                   ib=ib+1
!                   if (ib >= 10) then 
!                      step=stepnow/(10.0_r8**(ib-10)) ! this 10._r8, the _r8 is important.
!                      x=xrvs
!                      call gtest(x,g) 
!                   endif   
! 
!                  ! can add a random monte carlo search here perhaps.
!                   
!                   if (ib >= 300) then
!                       write(6,*) 'it looks like minimal is found and recorded already.'
!                       exit all
!                   endif    
!               endif   
!           enddo sub1
!           sub2: do while (.not.cdone)
!               t=tb+phi*(tb-ta) 
!               call ftabc(x,t,g,fn)
!               if (fn>fb) then
!                   tc=t
!                   fc=fn
!                   cdone=.true.
!                   write(6,*) 'ta,tb,tc: ', ta,tb,tc
!                   write(6,*) 'fa,fb,fc: ', fa,fb,fc
!               else
!                   ta=tb
!                   tb=t 
!                   fa=fb
!                   fb=fn
!                   write(6,*) 'try c again',fn,fb,fa
!               endif   
!           enddo sub2          
!     ! until here we find ta,tb,tc. next need to fine the t to minimize f
!                
!           step1=(tc-ta)/pieces
!           do i=1,pieces
!               t=ta+i*step1
!               call ftabc(x,t,g,fn)
!               trial1(i)=fn
!           enddo
!           
!           find1=minloc(trial1)
!           
!           t=ta+find1(1)*step1
!           xrvs=xo ! save the old x.
!           !xrvs=xo+(ta+min(pieces,max(0,find1(1)-1))*step1)*g
!           x=x+t*g   
!  
!           call ftest(x,fn)
!           xn=x
!           write(6,*) 'fn=', fn
!           !write(6,*) 'xnew=', xn          
!           call gtest(x,g) ! update g.
!           write (9,'(25g15.7)') xn,fn,g
!           error=abs(fn-fo)
!           errorx=abs(sum((xn(:,:)-xo(:,:))**2))
!           write(6,*) 'error,errorx=', error,errorx
!
!! update           
!           fa=fn
!           bdone=.false.
!           cdone=.false.           
!           ta=0           
!           xo=xn 
!           fo=fb
!           go=g
!           
!!!!!!!!!!!!!!!!!!!!                     
!   
!    l=l+1
!             
!           if ( (error <= small).and.(errorx <= tiny) ) then  
!               write(6,*) 'almost there, countdown =', ct
!               ct=ct-1
!           endif
!  
!    if (l > lmax) then
!        write(6,*) 'too many loops, exit do while: all'
!        exit all  
!    endif
!        
!    enddo all
!
!    close(9)
!    
!    return
!    end subroutine sdv6      
!
!    subroutine cgv6(x)
!!  Method of conjugate gradient
!    use matrixmod
!    use brut
!    use v6stepcalc
!    real(kind=r8) :: x(3,npart),xo(3,npart),xn(3,npart),xrvs(3,npart),xnew(3,npart) &
!                    ,lambda,t,error,errorx,ptgdiff,g(3,npart),go(3,npart),gn(3,npart) &
!                     ,f,fo,fn,ta,tb,tc,tx &
!                    ,ga(3,npart),gb(3,npart),gc(3,npart) &
!                    ,alpha,p(3,npart),po(3,npart),pn(3,npart) &
!                    ,step1,fa,fb,fc,step,step0,stepnow
!    real(kind=r8), allocatable :: trial1(:)
!    real(kind=r8), parameter :: small=1.0d-10,tiny=1.0d-10,epsilon=1.0d-10 &
!                               ,w=(3-sqrt(5.0_r8))/2,wg=(sqrt(5.0_r8)-1)/2,phi=(1+sqrt(5.0_r8))/2
!    integer(kind=i4) :: find1(1),pieces,i,j,k,l,lmax,ct,ia,ib,ic,ig
!    logical :: abcdone,adone,bdone,cdone
!
!    lmax=1000
!    step0=1.0d-5
!    step=step0
!    abcdone=.false.
!    adone=.true.
!    bdone=.false.
!    cdone=.false.
!    pieces=33
!    allocate(trial1(pieces))
!     
!    ta=0
!    xo=x
!
!    call ftest(xo,fo)
!    write(6,*) 'f initial =', fo
!    fa=fo
!    !write(6,*) 'x initial =', xo 
!    call gtest(xo,go)
!    write(6,*) 'g initial =', go 
!    g=go
!    
!    po=go
!    p=go
!    
!    l=0
!    ct=20
!    open(unit=9,form='formatted',file='cgv6',position='rewind')
!    write (9,'(13g15.7)') xo,fo
!    
!    all: do while (ct>=0)
!        
!!!!!!!!!!!!!!!!!!!!!!!!! 
!        
!! for a given x, find the min t.           
!           ib=0
!           sub1: do while (.not.bdone)          
!               t=ta+step 
!               call ftabc(x,t,p,fn)
!               if (fn<fa) then
!                   tb=t
!                   fb=fn
!                   bdone=.true.
!                   stepnow=step
!               else
!                   step=step/10  ! do this to ensure ta=0, fa>fb.
!                   write(6,*) 'try b again',step,fn,fo
!                   ib=ib+1
!                   if (ib >= 10) then 
!                      step=stepnow/(10.0_r8**(ib-10)) ! this 10._r8, the _r8 is important.
!                      x=xrvs
!                      call gtest(x,g) 
!                   endif   
! 
!                  ! can add a random monte carlo search here perhaps.
!                   
!                   if (ib >= 300) then
!                       write(6,*) 'it looks like minimal is found and recorded already.'
!                       exit all
!                   endif    
!               endif   
!           enddo sub1
!           sub2: do while (.not.cdone)
!               t=tb+phi*(tb-ta) ! phi is golden number 1.618 
!               call ftabc(x,t,p,fn)
!               if (fn>fb) then
!                   tc=t
!                   fc=fn
!                   cdone=.true.
!                   write(6,*) 'ta,tb,tc: ', ta,tb,tc
!                   write(6,*) 'fa,fb,fc: ', fa,fb,fc
!               else
!                   ta=tb
!                   tb=t 
!                   fa=fb
!                   fb=fn
!                   write(6,*) 'try c again',fn,fb,fa
!               endif   
!           enddo sub2          
!     ! until here we find ta,tb,tc. next need to fine the t to minimize f
!           
!!!!!!!!!  stupid search.              
!           !step1=(tc-ta)/pieces
!           !do i=1,pieces
!           !    t=ta+i*step1
!           !    call ftabc(x,t,p,fn)
!           !    trial1(i)=fn
!           !enddo
!           !
!           !find1=minloc(trial1)
!           !
!           !t=ta+find1(1)*step1
!           !xrvs=xo ! save the old x.
!           !!xrvs=xo+(ta+min(pieces,max(0,find1(1)-1))*step1)*g
!           !x=x+t*p    
!           !call ftest(x,fn)
!           !call gtest(x,g) ! update g.
!!!!!!!!!!!!!!!!
!           
!!!! add complete golden-section search to find x.
!           
!           ptgdiff=abs(2*(fb-fo)/(fb+fo))
!           ig=0
!           sub3: do while (ptgdiff>epsilon)
!             fo=fb
!             if ((tc-tb)>(tb-ta)) then 
!               tx=ta+wg*(tc-ta) 
!               call ftabc(x,tx,p,fn)  
!               if (fn>fb) then
!                tc=tx
!                fc=fn
!               else
!                   ta=tb
!                   fa=fb
!                   tb=tx
!                   fb=fn
!                   ptgdiff=abs(2*(fn-fo)/(fn+fo))
!               endif    
!             else 
!               tx=ta+w*(tc-ta) 
!               call ftabc(x,tx,p,fn) 
!               if (fn>fb) then
!                ta=tx
!                fa=fn
!               else  
!                   tc=tb
!                   fc=fb
!                   fb=fn
!                   tb=tx        
!                   ptgdiff=abs(2*(fn-fo)/(fn+fo))      
!               endif                         
!             endif 
!             ig=ig+1
!             !write(6,*) 'ig=',ig
!             !write(6,'(8g15.7)') ta,tb,tc,(tb-ta)/(tc-ta),(tc-tb)/(tc-ta),fo,fb,tc-ta
!             if (ig>100) then
!                 write(6,*) 'looks like golden-section cannot go further, stop sub3 do loop.',ptgdiff,fn,fo
!                 exit sub3
!             endif       
!           enddo sub3
!           xrvs=xo
!           x=x+tb*p  
!           call ftest(x,fn)
!           call gtest(x,g) ! update g.
!           
!     
!! new p.
!           alpha=sum(g(:,:)**2)/sum(go(:,:)**2)
!           p=g+alpha*po
!           xn=x
!           write(6,*) 'fn=', fn
!           !write(6,*) 'xnew=', xn          
!           write (9,'(38g15.7)') xn,fn,p,g,alpha
!           error=abs(fn-fo)
!           errorx=abs(sum((xn(:,:)-xo(:,:))**2))
!           write(6,*) 'error,errorx=', error,errorx
!
!! update           
!           fa=fn
!           bdone=.false.
!           cdone=.false.           
!           ta=0           
!           xo=xn 
!           fo=fb
!           go=g
!           po=p
!           
!!!!!!!!!!!!!!!!!!!!                     
!   
!    l=l+1 
!           if (mod(l,size(x)).eq.0) then  ! reset. this is important for cg.
!               p=g
!               po=g  
!               write(6,*) 'reset',l
!           endif
!           
!           
!           if ( (error <= small).and.(errorx <= tiny) ) then  
!               write(6,*) 'almost there, countdown =', ct
!               ct=ct-1
!           endif
!  
!    if (l > lmax) then
!        write(6,*) 'too many loops, exit do while: all'
!        exit all  
!    endif
!        
!    enddo all
!
!    close(9)
!    
!    return
!    end subroutine cgv6  
!!!!!!!!!!!! this part is just to check m6matrix end !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

    
    !subroutine ftabc(xin,tin,gin,f)
    !real(kind=r8) :: xin(3,npart),xnew(3,npart),f,gin(3,npart),tin
    !xnew=xin+tin*gin
    !call ftest(xnew,f)
    !return
    !end subroutine ftabc    
    !
    !subroutine ftest(xin,f)
    !use v6stepcalc 
    !real(kind=r8) :: xin(3,npart),f,vijeigval(nbasis)
    !complex(kind=r8) :: vij(nbasis,nbasis),vijeigvec(nbasis,nbasis)
    !call v6mat(xin,vij,vijeigvec,vijeigval)
    !f=vijeigval(1)
    !return
    !end subroutine ftest
    !
    !subroutine gtest(xin,g) ! g=-gradiant
    !use v6stepcalc
    !real(kind=r8) :: g(3,npart)  
    !real(kind=r8) :: dx=0.001_r8
    !real(kind=r8) :: xin(3,npart),f,vijeigval(nbasis),gp,g0,xnew(3,npart)
    !complex(kind=r8) :: vij(nbasis,nbasis),vijeigvec(nbasis,nbasis)
    !integer(kind=i4) :: i,j,k
    !call v6mat(xin,vij,vijeigvec,vijeigval)
    !g0=vijeigval(1)    
    !do i=1,npart 
    !    do j=1,3
    !        xnew=xin
    !        xnew(j,i)=xin(j,i)+dx
    !        call v6mat(xnew,vij,vijeigvec,vijeigval)
    !        gp=vijeigval(1)
    !        g(j,i)=-(gp-g0)/dx
    !    enddo
    !enddo
    !return
    !end subroutine gtest    

    
    
    !subroutine test1229a
    !use brut
    !use v6stepcalc
    !real(kind=r8) :: x(3,npart) ! the configuration for the system
    !complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
		  !              ,cwtnew1(0:nspin-1,nisospin),cwtnew2(0:nspin-1,nisospin) &
		  !              ,cwt1(0:nspin-1,nisospin),cwtgnd(0:nspin-1,nisospin)
    !integer :: ipl(6),jpl(6),ipr(6),jpr(6)
    !call xoutput(x)
    !call getcwtgnd2(cwtgnd)
    !call outputordlr(ipl,jpl,ipr,jpr)
	   !
    !call psitprop(x,cwtgnd,cwtnew)
	   !
    !call corop(x,ipl,jpl,cwtnew1)
	   !
	   !
    !!write(6,*) 'cwtnew1-cwtnew=', cwtnew1-cwtnew
    !!
    !!write(6,*) 'ipl=', ipl
    !!write(6,*) 'jpl=', jpl
    !!
    !!write(6,*) 'cwtnew1=', cwtnew1
    !!write(6,*) 'cwtnew=', cwtnew
	   !
    !open(11,FILE='cwtmine.txt',FORM='FORMATTED')
    !write (11,'(2x,2(f15.7,2x))') cwtnew
    !close(11)
    !
    !open(11,FILE='cwtSchmidt.txt',FORM='FORMATTED')
    !write (11,'(2x,2(f15.7,2x))') cwtnew1
    !close(11)
	   !
	   !
    !return
    !end subroutine test1229a        


    
    !!!------------------------------------------------   
    !!! test      
    !   call coropr(xr,cwt)
    !   val1=sum(conjg(cwtl2r(nchorizo,:,:))*cwt(:,:))
    !   
    !   !write(6,*) 'cwt=', cwt(:,:)
    !   !write(6,*) 'cwtl2r=', cwtl2r(nchorizo,:,:)  
    !   !write(6,*) 'xl before=', xl
    !   
    !   call coropl(xl,cwt)
    !   
    !   !write(6,*) 'xl after=', xl
    !    
    !   val2=sum(conjg(cwtr2l(0,:,:))*cwt(:,:))
    !
    !   !write(6,*) 'cwt=', cwt(:,:)
    !   !write(6,*) 'cwtr2l=', cwtr2l(0,:,:)
    !   
    !   write(6,*) 'val1, val2=',val1,val2
    !!!! result is correct, val1=val2.
    !   
    !   call outputordlr(iplo,jplo,ipro,jpro)
    ! 
    !   write(6,*) '------------------------------------------------------'
    !   write(6,*) 'iplo=',iplo
    !   write(6,*) 'jplo=',jplo   
    !   
    !   write(6,*) 'ipro=',ipro
    !   write(6,*) 'jpro=',jpro
    !   
    !   xtmp=xl
    !   call hpsitcwt(xtmp,iplo,jplo,cwtr2l(0,:,:),psi20,e0,d20,pe0,f0,cwta1)
    !   
    !  ! write(6,*) 'xl after=', xl
    !   
    !   !write(6,*) 'cwtr2l=', cwtr2l(0,:,:)
    !   !write(6,*) 'cwta1=', cwta1(1,:,:)
    !   
    !   !write(6,*) 'x before=', x  
    !   !call coropv(xl,iplo,jplo,dxx,cwta)
    !   !write(6,*) 'x after=', x  
    !   !write(6,*) 'cwta=', cwta(1,:,:)  
    !   !write(6,*) 'cwta1-cwta=', cwta1(1,:,:)-cwta(1,:,:)  
    !   
    !   !write(6,*) 'xl before=', xl
    !   call coropl(xl,cwt)
    !   !write(6,*) 'xl after=', xl   
    !   write(6,*) 'cwt-cwta1=', cwt-cwta1(1,:,:)
    !   !write(6,*) 'x=', x
    !   !write(6,*) 'x1=', ristra(1)%x
    !   !write(6,*) 'xl=', xl 
    !   !write(6,*) 'x0=', ristra(0)%x   
    !   
    !   xtmp=xr
    !   call hpsitcwt(xtmp,ipro,jpro,cwtl2r(nchorizo,:,:),psi2n,en,d2n,pen,fn,cwta)
    !   call coropr(xr,cwt)
    !   write(6,*) 'cwt-cwta=', cwt-cwta(1,:,:)   
    !   
    !! Schmidt order and my order are different, check and fix. it doesn't matter I think.  
    !   
    !   write(6,*) 'psi20, psi2n (should be the same)=',psi20,psi2n    
    
    !subroutine hpsi(c0,cn,valnum,valdenom,rgl,rgr,rg,f,rf,absrf) ! not really hpsi anymore.
    !use v6stepcalc   
    !use chorizos
    !type (chorizo) :: c0,cn
    !complex(kind=r8) :: cwtr2lout(0:nchorizo,0:nspin-1,nisospin),cwtl2rout(0:nchorizo,0:nspin-1,nisospin) &
	    !                ,psith0cwt(0:nspin-1,nisospin),nhpsitcwt(0:nspin-1,nisospin) &
	    !                ,cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin) &
	    !                ,cwtnew1(0:nspin-1,nisospin),cwtold1(0:nspin-1,nisospin) &
	    !                ,cwtnew2(0:nspin-1,nisospin),cwtold2(0:nspin-1,nisospin)	   
    !complex(kind=r8) :: gl,gr,f
    !real(kind=r8) :: rg,rgl,rgr,rf,absrf,valnum,valdenom
    !real(kind=r8) :: x0(3,npart),xn(3,npart)
    !integer(kind=i4) :: i,j,k
    !integer(kind=i4) :: invspin(nbasis),invispin(nbasis)
    !
    !cwtr2lout=cwtr2l
    !cwtl2rout=cwtl2r
    !x0=c0%x
    !xn=cn%x
    !call getcwtgnd(cwtold,invspin,invispin) ! cwtold is cwtgnd   
    !
    !do k=1,nbasis
	    !call v6val1(x0,k,cwtnew1)
	    !psith0cwt(invspin(k),invispin(k))=sum(cwtnew1*cwtold) 
	    !call v6val1(xn,k,cwtnew2)
	    !nhpsitcwt(invspin(k),invispin(k))=sum(conjg(cwtnew2)*cwtold) 
    !enddo
    !
    !gl=sum(psith0cwt(:,:)*cwtr2l(0,:,:))
    !rgl=real(gl)
    !gr=sum(conjg(cwtl2r(nchorizo,:,:))*nhpsitcwt(:,:))
    !rgr=real(gr)   
    !rg=0.5*(rgl+rgr)
    !
    !f=sum(cwtold(:,:)*cwtr2l(0,:,:))
    !rf=real(f)
    !absrf=abs(rf)
    !
    !valnum=rg/absrf
    !valdenom=rf/absrf
    !
    !!write(6,*) 'gl=',gl
    !!write(6,*) 'gr=',gr   
    !!write(6,*) 'f=',f  
    !!
    !!write(6,*) 'cwtr2l(0)=',cwtr2l(0,:,:)
    !!write(6,*) 'psith0cwt=',psith0cwt(:,:)
    !return
    !end subroutine hpsi       
    
   
   
    end module step
