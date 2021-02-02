module step
    use chorizos
    implicit none
    integer, private, parameter :: i4=selected_int_kind(9)
    integer, private, parameter :: i8=selected_int_kind(15)
    integer, private, parameter :: r8=selected_real_kind(15,9)
    type (chorizo), private, allocatable, save :: ristra(:),ristraold(:),ristranew(:),ristranew1(:)
    type (chorizo), private, allocatable, save :: newristra(:),oldristra(:)
    real(kind=r8), private, parameter ::  hc=197.327053d0,mpi0=134.9739d0,mpic=139.5675d0,mp=938.27231d0,mn=939.56563d0 &
                                         ,alpha=1.0d0/137.035989,mup= 2.7928474d0,mun=-1.9130427d0,zero=0.0d0 &
                                         ,tiny=1.0d-8
    complex(kind=r8), private, parameter :: czero=(0.0_r8,0.0_r8)
    complex(kind=r8), private, parameter :: ci=(0.0_r8,1.0_r8)
    complex(kind=r8), private, parameter :: cone=(1.0_r8,0.0_r8)
    integer(kind=i4), private, save :: npart,nprot,npair,nspin,nisospin,nbasis
    integer(kind=i4), private, save :: nav,neq,nstep,nstepcorrchk,nstepnow,nstepstuck,nblocknow,nblockstuck
    integer(kind=i4), private, save :: nchorizo,nchorizomid,nchorizomod,nrepmax,irep,nstepdecor,nrhobin
    integer(kind=i4), private, save :: nest,nestristra
    integer(kind=i4), private, save :: mmax,nbisect,nbisecthalf,mmaxhalf,nreparrow,icrephistlast,icrephist
    real(kind=r8), private, save ::  hbar,dt,sigma,driftm,el,repprob
    real(kind=r8), private, save ::  etrial,ecut,eloctot,rhobinsize
    real(kind=r8), private, save ::  x0step,mov1step,mov2step,shftstep
    real(kind=r8), private, save ::  rbisect,rmov1,rmov2,rmov3,rmov23
    integer(kind=i4), private, save :: ic0,ic0tot,icn,icntot,ibisect,ibisecttot,ibisectextra,ibisectextratot,ieloc &
                                        ,icmov1,icmov1rvs,icmov1tot,icmov1rvstot,icstdmove1tot,icstdmove1rvstot &
                                        ,icshft,icshfttot,ic23,ic23tot,ice,icbisecttot &
	                                    ,ibisectimptot,icstepstd &
                                        ,ibisectl,ibisecttotl,ibisectr,ibisecttotr &
                                        ,iclr,iclrtot,icrep,icreptot &
                                        ,ibisecthalf,ibisecthalftot,ibisectdouble,ibisectdoubletot &
                                        ,icstdmove2tot,icstdmove2 &
                                        ,irepstepct,irepstepmonmax &
                                        ,icstdmove1bisect1rvstot,icmov1bisect1rvstot,icstdmove1bisect1tot &
                                        ,icmov1bisect1tot,icmov1bisect1rvs,icmov1bisect1
    integer(kind=i4), private, save, allocatable :: repstepct1(:),repstepct2(:),repstepct1sum(:),repstepct2sum(:) &
                                                   ,ireparrowhistory(:) &
                                                   ,ibisectstuck(:)
    real(kind=r8), private, save, allocatable ::  bisectrate(:),bisectextrarate(:),bisectcount(:),bisectextracount(:) &
                                                 ,bisectextratot(:),bisecttot(:) &
                                                 ,bisectcountl(:),bisecttotl(:) &
                                                 ,bisectcountr(:),bisecttotr(:) &
                                                 ,bisecthalfcount(:),bisecthalftot(:)
    complex(kind=r8), private, save, allocatable :: cwtbegin(:,:),cwtend(:,:),cwtl2r(:,:,:),cwtr2l(:,:,:) &
                                                    ,cwtbegin1(:,:),cwtend1(:,:),cwtl2r1(:,:,:),cwtr2l1(:,:,:) &
                                                    ,cwtbeginnew(:,:),cwtendnew(:,:),cwtl2rnew(:,:,:),cwtr2lnew(:,:,:) &
                                                    ,cwtbeginnew1(:,:),cwtendnew1(:,:),cwtl2rnew1(:,:,:),cwtr2lnew1(:,:,:)  &
                                                    ,cwtground(:,:)
    real(kind=r8), private, save, allocatable :: xin(:),x1in(:),x2in(:),corr(:),corr1(:),corr2(:) ! for correlation check.
    integer(kind=i4), private, allocatable, save :: invspin(:),invispin(:),liso(:),lisoi(:)
    logical, private, save :: icorrchk,isite
    character(len=70), private, save :: infile,outfile
    integer(kind=i4), private, save, allocatable :: iplo(:),jplo(:),ipro(:),jpro(:),ipl(:),jpl(:),ipr(:),jpr(:)	
    logical, private, save :: rejectbisectv6improve,rejectbisectv6improvehalf
    character(len=120), private, save, allocatable :: answer(:)
    integer(kind=i4), private, save :: ipath,ipathtmp,npo,ixtempo,ipathcountold 
    integer(kind=i4), private, save :: p3map(3,6),p4map(4,24),p5map(5,120),p6map(6,720),p7map(7,5040) ! permutation map. 
    real(kind=r8), private, save :: pathcorevalnew,pathcoreval ! the current abs real f without gaussian part and the normalization, the thing.
    logical, private, save :: iem,spinorbit
    real(kind=r8), private, save :: time0,time1,time2,time3,time4,time5,time00 &
                                   ,time0tot,time1tot,time2tot,time3tot,time4tot,time5tot
    integer(kind=i4), private, save :: irlintervalmax,irlintervalmax1,irlstep1,irlstep2,ne0respbeads,nrespest  
    real(kind=r8), private, save, allocatable :: ebeads(:,:,:),ewbeads(:,:,:),e0beads(:,:),ew0beads(:,:) &
                                                ,e0beadscrs(:,:),ew0beadscrs(:,:)
    logical, private, save :: resume
    integer(kind=i4), parameter :: ipathmark=2147483647 !2**31-1
    integer(kind=i4), parameter :: pathunfstep1=6 ! if the structure of path.unf changed then pathunfstep1 need to be changed accordingly.
    integer(kind=i4) :: ipathmarkchk
    integer(kind=i4), private, save :: nprocess
    real(kind=r8), private, save, allocatable :: xgatherall(:) ! for mpi_gather, this will be bigger when beads * cores are big. need mpi_write
    real(kind=r8), private, save :: sizexgatherall
    integer(kind=i4), private, save, allocatable :: icheckhpsi1d(:) ! for mpi_gather 
    integer(kind=i4), allocatable :: iplgatherall(:),jplgatherall(:),iprgatherall(:),jprgatherall(:)	! for mpi_gather
    integer(kind=i4), allocatable :: iplocache(:),jplpcache(:),iprcache(:),jprcache(:)	! cache
    real(kind=r8), private, save, allocatable :: x0tot1dcache(:)    
    logical :: enablewritepathgs
    integer(kind=i4), private, save :: cvalstrange
    character(len=120), private, save :: file_path,file_path0,file_pathnumbercount,file_path0numbercount &
                                        ,file_checkpoint,file_pimc_checkpoint,file_values
contains
    subroutine stepinit(npartin,nchorizoin,hbarin,dtin &
        ,mmaxin,nrepmaxin,iemin,spinorbitin,irepin,nstepdecorin &
	    ,x0stepin,mov1stepin,mov2stepin,shftstepin,nrhobinin,icorrchkin &
	    ,nprotin,navin,nstepin,neqin,rbin,r1in,r2in,r3in,r23in &
        ,nestin,nestristrain,resumein)
    use math
    use v6stepcalc
    use mympi
    real(kind=r8) :: hbarin,dtin &
                    ,x0stepin,mov1stepin,mov2stepin,shftstepin &
                    ,rbin,r1in,r2in,r3in,r23in
    integer(kind=i4) :: npartin,nchorizoin,nrepmaxin,irepin,mmaxin,nstepdecorin,nrhobinin,nprotin &
                        ,navin,nstepin,neqin
    integer(kind=i4) :: nestin,nestristrain
    integer(kind=i4) :: i,nristra
    logical :: icorrchkin,iemin,resumein,spinorbitin
    !repprob=repprobin
    
    npart=npartin  
    npair=npart*(npart-1)/2
    allocate(iplo(npair),jplo(npair),ipro(npair),jpro(npair),ipl(npair),jpl(npair),ipr(npair),jpr(npair))    

    x0step=x0stepin
    !nrepmax=nrepmaxin
    mmax=mmaxin
    nrepmax=nrepmaxin
    iem=iemin
    spinorbit=spinorbitin
    irep=irepin
    nreparrow=1 ! initialize reptation direction. 1 is right, -1 is to the left.  
    irepstepct=0
     
    resume=resumein
    ! ipath series is only useful for rank0.
    if (.not.resume) then
        ipath=1 ! initialize the number of path as 1. count the path number for equilibriated path.
        ipathtmp=1 ! count the path number for unequilibriated path.
    endif

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
    allocate( bisectrate(0:mmax),bisectextrarate(0:mmax),bisectcount(0:mmax),bisectextracount(0:mmax) &
            ,bisectextratot(0:mmax),bisecttot(0:mmax) &
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
    ,cwtl2r(0:nspin-1,nisospin,0:nchorizo),cwtr2l(0:nspin-1,nisospin,0:nchorizo))
    allocate(cwtbegin1(0:nspin-1,nisospin),cwtend1(0:nspin-1,nisospin) &
    ,cwtl2r1(0:nspin-1,nisospin,0:nchorizo),cwtr2l1(0:nspin-1,nisospin,0:nchorizo))
    allocate(cwtbeginnew(0:nspin-1,nisospin),cwtendnew(0:nspin-1,nisospin) &
    ,cwtl2rnew(0:nspin-1,nisospin,0:nchorizo),cwtr2lnew(0:nspin-1,nisospin,0:nchorizo))
    allocate(cwtbeginnew1(0:nspin-1,nisospin),cwtendnew1(0:nspin-1,nisospin) &
    ,cwtl2rnew1(0:nspin-1,nisospin,0:nchorizo),cwtr2lnew1(0:nspin-1,nisospin,0:nchorizo))   
    allocate(cwtground(0:nspin-1,nisospin))
    allocate(invspin(nbasis),invispin(nbasis))  
    call getcwtgnd(cwtground,invspin,invispin)
    
    nav=navin
    nstep=nstepin
    neq=neqin
    
    if (icorrchk) then
        allocate(xin(nstep),x1in(nstep),x2in(nstep),corr(0:nstep),corr1(0:nstep),corr2(0:nstep)) ! for correlation check.
    endif
    
    allocate(ibisectstuck(nav+neq))
    ibisectstuck=0 ! initialize it. check how many times bisection stuck.    
    
    if (nrepmax /= 0) then ! for reptation.
        irepstepmonmax=10*nchorizo
        allocate(repstepct1(0:irepstepmonmax),repstepct2(0:irepstepmonmax))
        allocate(repstepct1sum(0:irepstepmonmax),repstepct2sum(0:irepstepmonmax))
        repstepct1=0
        repstepct2=0
        repstepct1sum=0
        repstepct2sum=0         
        allocate(ireparrowhistory(0:nav*nstep))
        ireparrowhistory=0 
    endif
   
    rbisect=rbin
    rmov1=r1in
    rmov2=r2in
    rmov3=r3in
    rmov23=r23in
    
    call mathinit
    call getpermutmap(3,p3map)
    call getpermutmap(4,p4map)
    call getpermutmap(5,p5map)
    call getpermutmap(6,p6map)
    call getpermutmap(7,p7map)
    
    allocate(lisoi(0:nspin-1),liso(nisospin))
    call getliso(liso,lisoi)

 
    nest=nestin
    nestristra=nestristrain
      
    nprocess=nproc()
    allocate(icheckhpsi1d(0:nprocess-1))
    
    !sizexgatherall=dble(3*npart*(nchorizo+1)*nproc()*32)/dble(1024**3) ! size unit is GB.
    !if (sizexgatherall>3.5) then
    !    enablewritepathgs=.false.
    !    if (myrank().eq.0) write (6,*) 'size of xgatherall array is too large, writepathgs is Disabled.', sizexgatherall   
    !else
    !    enablewritepathgs=.true.
    !    allocate(xgatherall(3*npart*(nchorizo+1)*nproc()))
    !    allocate(iplgatherall(npair*nproc()))
    !    allocate(jplgatherall(npair*nproc()))
    !    allocate(iprgatherall(npair*nproc()))
    !    allocate(jprgatherall(npair*nproc()))        
    !endif
    enablewritepathgs=.false.
    
    cvalstrange=0 ! initially there is no strange value, the counting set to 0. for checkhpsi.

    return
    end subroutine stepinit
        
    subroutine stepfilesinit(irepin)
    integer(kind=i4) :: irepin
    character(len=120) :: file_id  
    
    select case (irepin) ! check cleanfolder subroutine match the name.
    case (4) ! VMC
        write(file_path,'("path_VMC.unf")')  
        write(file_path0,'("path0_VMC.unf")')
        write(file_pathnumbercount,'("pathnumbercount_VMC.txt")') 
        write(file_path0numbercount,'("path0numbercount_VMC.txt")') 
        write(file_checkpoint,'("checkpoint_VMC.unf")')        
        write(file_values,'("values_VMC.txt")')         
       
    case (5:6) ! Calc Mode
        write(file_path,'("path.unf")')  
        write(file_path0,'("path0.unf")')
        write(file_pathnumbercount,'("pathnumbercount.txt")') 
        write(file_path0numbercount,'("path0numbercount.txt")') 
        
        write(file_pimc_checkpoint,'("checkpoint.unf")') 
        write(file_id,'(i5)') irepin
        file_checkpoint='checkpoint_CalcMode' // trim(adjustl(file_id)) // '.unf'    
        file_values='values_CalcMode' // trim(adjustl(file_id)) // '.txt'   
       
    case default ! PIMC
        
        write(file_path,'("path.unf")')  
        write(file_path0,'("path0.unf")')
        write(file_pathnumbercount,'("pathnumbercount.txt")') 
        write(file_path0numbercount,'("path0numbercount.txt")') 
        write(file_checkpoint,'("checkpoint.unf")')        
        write(file_values,'("values.txt")')  

    end select     
    return
    end subroutine stepfilesinit

    subroutine checkresume(irepin,resumenow,outfilenow,isitenow) ! call stepfilesinit(irep) first
    use mympi
    character(len=70) :: outfilenow
    logical :: resumenow,alive,alive1,alive2,alive3,alive1a,isitenow
    integer(kind=i4) :: irepin
    if (resumenow) then
        alive=.false.
        if ((irepin.eq.2).or.(irepin.eq.4)) then
                inquire(file=file_checkpoint,exist=alive1)
                inquire(file=file_path,exist=alive2)
                inquire(file=outfilenow,exist=alive3)
                if (.not.alive1) write (6,*) 'Missing: ', file_checkpoint
                if (.not.alive2) write (6,*) 'Missing: ', file_path
                if (.not.alive3) write (6,*) 'Missing: ', outfilenow
                if (alive1.and.alive2.and.alive3) alive=.true.
        else if ((irepin.eq.5).or.(irepin.eq.6)) then
                inquire(file=file_checkpoint,exist=alive1)               
                inquire(file=file_path,exist=alive2)              
                inquire(file=file_pathnumbercount,exist=alive3)  
                if (.not.alive1) write (6,*) 'Missing: ', file_checkpoint
                if (.not.alive2) write (6,*) 'Missing: ', file_path
                if (.not.alive3) write (6,*) 'Missing: ', file_pathnumbercount
                if (alive1.and.alive2.and.alive3) alive=.true.             
        else
            write (6,*) 'checkresume Stop because irep is not defined, ', irepin
            call abort
        endif
        if (alive) then
            write (6,'(''============== Resume Mode On =============='')') 
            isitenow=.false. ! if resume then isite has to be false no matter what.
        else
            write (6,'(''=== No checkpoint => Resume Mode Off => Normal Mode On ==='')') 
            resumenow=.false.
            isitenow=.true.
        endif     
    endif     
    return
    end subroutine checkresume 
        
    subroutine cleanfolder(irepin,resumein)
    use estimator
    integer(kind=i4) :: irepin
    logical :: resumein
    character(len=120) :: file_pathold,file_path0old,file_checkpointold
    ! clean those file with position=append

    file_pathold=trim(file_path) // '.old'     
    file_path0old=trim(file_path0) // '.old'  
    file_checkpointold=trim(file_checkpoint) // '.old'     
    
    select case (irepin)     
    case (2,4) ! PIMC, VMC
        if (.not.resumein) then
            call rename(file_path,trim(file_pathold))  ! do a protection here.
            call rename(file_path0,trim(file_path0old))
            call rename(file_checkpoint,trim(file_checkpointold))
            
            open(unit=19,form='unformatted',file=file_path,status='replace',position='rewind') 
            close(19)
            open(unit=19,form='unformatted',file=file_path0,status='replace',position='rewind') 
            close(19)
            open(unit=19,form='unformatted',file=file_checkpoint,status='replace',position='rewind')
            close(19,status='delete')
            open(unit=11,form='formatted',file=file_values,status='replace',position='rewind')
            close(11) 
            if (irepin.eq.2) then
                open(unit=19,form='formatted',file='reparrowhist.txt',status='replace',position='rewind')
                close(19)   
            endif
            
            call estimatorfilesclean
              
        endif 
    case (5:6) ! Calc Modes 5, 6.
        if (.not.resumein) then
            !open(unit=19,form='formatted',file='reparrowhist.txt',status='replace',position='rewind')
            !close(19)
            !open(unit=11,form='formatted',file='values.txt',status='replace',position='rewind')
            !close(11)  
            call estimatorfilesclean
        endif       
        
    case default
        write (6,'(''irep value has not been decided in subroutine cleanfolder yet '')') irepin
        call abort         
           
    end select

    return
    end subroutine cleanfolder         
   
    subroutine stepchorizoinit(isitein,infilein,outfilein) ! call stepinit first
! initialize the position, calculate the weights for each bead. l2r,r2l
    use wavefunction
    use v6stepcalc
    use estimator
    use estimatorristra
    use random
    use brut
    use mympi
    integer(kind=i4) :: i,j,k,nptot,ixtemp,icheck,l
    real(kind=r8) :: xr(3,npart),xl(3,npart),x(3,npart),x0(3,npart) &
                    ,x0tot(3,npart,0:nchorizo),x0vmc(3,npart)
    real(kind=r8) ::  rangex,tmp,tmp1,tmp2,tmp3,r
    character(len=70) :: infilein,outfilein  
    logical :: isitein
    real(kind=r8) :: psi20,psi2n
    real(kind=r8) :: vn,vd,rf,signrf,vnke,vnpev6,vnpeem,vnpels,vnum,vdenom,rm,rmsq,rmn,rmnsq,rmp,rmpsq
    real(kind=r8), dimension(0:nrhobin) :: rhodistout
    character(len=120), dimension(:), allocatable :: answer
    real(kind=r8) :: x0tot1d(3*npart*(nchorizo+1)),x0tot1dnow(3*npart*(nchorizo+1))
    character(len=120) :: filename
    integer(kind=i4) :: myunit,ndummy
    logical :: hpsipass
    complex(kind=r8) :: cwtl2rm(0:nspin-1,nisospin),cwtr2lm(0:nspin-1,nisospin),cwttmp(0:nspin-1,nisospin)
    
    ! call stepinit first   
    isite=isitein
    infile=infilein
    outfile=outfilein

    !ixtemp=3*npart*(nchorizo+1)  
   
    if (isite) then
        if (myrank().eq.0) write (6,'(''input from sites'')')
	    rangex=1.12*dble(npart)**(1.0_r8/3) ! +0.54*log(9.0_r8)   ! decrease to 1/10 rho0, Fermi charge distribution.  x0step
	    if (irep.eq.4) then
            x0vmc=rangex*reshape((randn(size(x0))-0.5_r8),shape(x0)) ! let all the beads begin from the same x0. good for vmc.
            do i=0,nchorizo
	            x0=x0vmc 
                call chorizoin(x0,i) 
            enddo             
        else
            do i=0,nchorizo
                
                
	            x0=rangex*reshape((randn(size(x0))-0.5_r8),shape(x0)) ! all the beads begin from the different x0. good for pimc.
       
                call chorizoin(x0,i) 
                
                ! check H2 rij
                !r=sqrt(sum((x0(:,1)-x0(:,2))**2))
                !write (6,*) 'i chorizo, r12=',i,r
                !
                
                
            enddo             
        endif
        
        call outputordlr(iplo,jplo,ipro,jpro) ! bc brutinit has already been called.  
        call chorizoallout(x0tot,iplo,jplo,ipro,jpro)
        !x0tot1d=reshape(x0tot,(/3*npart*(nchorizo+1)/)) 
            
    else  
    
        if (.not.resume) then
            if (myrank().eq.0) then
                write (6,'(''input from file'',t40,a20)') infile
                write (12,'(''input from file'',t40,a20)') infile
            endif
            call readconfig(infile,iplo,jplo,ipro,jpro,x0tot1d)
        else ! resume mode
            if (myrank().eq.0) then
                write (6,'(''resume mode, input from file'',t40,a20)') outfile
                write (12,'(''resume mode, input from file'',t40,a20)') outfile
            endif  
            call readconfig(outfile,iplo,jplo,ipro,jpro,x0tot1d)
        endif

        x0tot=reshape(x0tot1d,shape(x0tot))
        call chorizoallin(x0tot)
        call inputordlr(iplo,jplo,ipro,jpro)     

    endif
 
    if (myrank()==0) write (6,'(''output to file'',t40,a20)') outfile    
    xr=ristra(nchorizo)%x
    xl=ristra(0)%x 

    call corop(xl,iplo,jplo,cwtbegin) ! psitcwt    
    call v6propr(xl,-1,cwtbegin,cwtl2r(:,:,0)) 
    do i=1,nchorizo-1
	    x=ristra(i)%x
        call vlsprop(spinorbit,-1,ristra(i-1)%x,x,cwtl2r(:,:,i-1),cwttmp(:,:))
	    call v6proplr(x,x,-1,cwttmp(:,:),cwtl2r(:,:,i))    
    enddo 

    call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,nchorizo-1),cwttmp(:,:))
    call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))
    
    
    !write (6,*) '##################### above is L2R ##########################'
    
    
    ! right to left  
    call corop(xr,ipro,jpro,cwtend)
    call v6propr(xr,-1,cwtend,cwtr2l(:,:,nchorizo))
    do j=nchorizo-1,1,-1
	    x=ristra(j)%x        
        call vlsprop(spinorbit,1,x,ristra(j+1)%x,cwtr2l(:,:,j+1),cwttmp(:,:))
	    call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,j))   
    enddo
    
    call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
    call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))
    
    !write (6,*) '##################### above is R2L ##########################'
    
    call hpsi(ristra(0),ristra(nchorizo),vn,vd,vnke,vnpev6,vnpeem,vnpels,psi20,psi2n,rf,icheck)    
    pathcoreval=abs(rf) ! the f without gaussian part and the normalization, the thing.
    
    
    
     !!!!!!!!!!!! check  
        
        !write (6,*) 'icheck=',icheck
        !
        !tmp1=abs(real(sum(conjg(cwtbegin(:,:))*cwtr2l(:,:,0))))
        !
        !!call corop(xr,ipro,jpro,cwtendtmp)
        !!write (6,*) 'bisectv6bflex    cwtend = ', cwtend
        !!write (6,*) 'cwtendtmp = ',cwtendtmp
        !!write (6,*) 'cwtend diff = ',cwtendtmp-cwtend
        !
        !write (6,*) 'cwtbegin = ', cwtbegin
        !write (6,*) 'cwtend = ', cwtend
        !
        !do k=0,nchorizo
        !    write (6,*) 'k=,',k
        !write (6,*) 'cwtr2l(:,:,k) = ', cwtr2l(:,:,k)
        !!write (6,*) 'cwtr2l(:,:,nchorizo) = ', cwtr2l(:,:,nchorizo)
        !enddo
        
        
    
        !do k=0,nchorizo
        !    write (6,*) 'k=,',k
        !    if (k<=nchorizo-1) then
        !        call vlsprop(spinorbit,1,ristra(k)%x,ristra(k+1)%x,cwtr2l(:,:,k+1),cwttmp(:,:))
        !        tmp2=abs(real(sum(conjg(cwtl2r(:,:,k))*cwttmp(:,:))))    
        !    else
        !        tmp2=abs(real(sum(conjg(cwtl2r(:,:,nchorizo))*cwtend(:,:))))
        !    endif
        !    
        !    tmp=(tmp1+tmp2)/2  
        !    tmp3=abs(  (tmp1-tmp2)/( (tmp1+tmp2)/2  ) )*100
        !    write(6,'(''pathcoreval check stepchorizoinit :'', 5(g15.7,1x))') ,pathcoreval, tmp,tmp1,tmp2,tmp3
        !    if (  tmp3 >= 0.1   ) then
        !        write(6,*) '################## Caution ###################'
        !        write(6,'(''pathcoreval check stepchorizoinit failed err = '',g15.7, ''%'' )')  tmp3
        !        write(6,*) '###############################################'
        !        !stop
        !    endif
        !enddo
        !
        !write (6,*) 'check again---'
        !tmp2=abs(real(sum(conjg(cwtl2r(:,:,nchorizo))*cwtend(:,:))))
        !tmp=(tmp1+tmp2)/2  
        !    tmp3=abs(  (tmp1-tmp2)/( (tmp1+tmp2)/2  ) )*100
        !    write(6,'(''pathcoreval check stepchorizoinit :'', 5(g15.7,1x))') ,pathcoreval, tmp,tmp1,tmp2,tmp3
        !    if (  tmp3 >= 0.1   ) then
        !        write(6,*) '################## Caution ###################'
        !        write(6,'(''pathcoreval check stepchorizoinit failed err = '',g15.7, ''%'' )')  tmp3
        !        write(6,*) '###############################################'
        !        !stop
        !    endif        
        !write (6,*) 'check again over---'
    
     !!!!!!!!!!!!!!!!!!      
    
    
    
    

! chekhpsi
    call checkhpsi(vn,vd,vnke,vnpev6,vnpeem,vnpels,psi20,psi2n,rf,icheck,iplo,jplo,ipro,jpro,x0tot1d,-1,ndummy,hpsipass)
    if (hpsipass) then
        if (myrank().eq.0) then
            write (6,'(/,''hpsi check pass! eg: rank 0 result is: '',3(1x,g15.7) )') psi20,psi2n,rf 
            write (12,'(/,''hpsi check pass! eg: rank 0 result is: '',3(1x,g15.7) )') psi20,psi2n,rf 
        endif     
    endif
!
    if (.not.resume) then
        call zerest
        call zerrhodist 
        vnum=vn
        vdenom=vd
        call addval(3,vnum,1.0_r8)
        call addval(4,vdenom,1.0_r8)
        call addval(5,rf,1.0_r8)
        
        call getcwtlrm(nchorizomid,cwtl2rm,cwtr2lm)
        call samplerhodistrbt(ristra(nchorizomid)%x,nchorizomid,cwtl2rm,cwtr2lm,rf,rhodistout,rm,rmsq,rmn,rmnsq,rmp,rmpsq)
        call addval(6,rmsq,1.0_r8)
        call addval(7,rm,1.0_r8)
        call addval(8,rhodistout(0),1.0_r8) 
        call addval(13,rmnsq,1.0_r8)
        call addval(14,rmn,1.0_r8)      
        call addval(15,rmpsq,1.0_r8)
        call addval(16,rmp,1.0_r8)
 
        call addval(17,vnke,1.0_r8)
        call addval(18,vnpev6+vnpeem+vnpels,1.0_r8)
        call addval(19,vnpev6,1.0_r8)
        call addval(20,vnpeem,1.0_r8) 
        call addval(21,vnpels,1.0_r8) 
        
        call update   !collect block averages
        if (myrank().eq.0) then
            allocate(answer(nest))
            answer=resstring()
            call writeanswer(answer) 
            deallocate(answer)
            write (6,'(/,''initial config calculated. ----------------------------------------'')') 
	        write (12,'(/,''initial config calculated. ----------------------------------------'')')    
        endif  
   
    endif
    
    
    !!!!
        !stop
    !!!!
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
    real(kind=r8) ::  newp,oldp,rn(1),tmp1test1
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin)
    real(kind=r8) :: tmp,tmp1,tmp2,tmp3
   
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
    call v6propr(xr,-1,cwtendnew,cwtr2lnew(:,:,nchorizo))
    do j=nchorizo-1,1,-1
	    x=ristra(j)%x     
        call vlsprop(spinorbit,1,x,ristra(j+1)%x,cwtr2lnew(:,:,j+1),cwttmp(:,:))
	    call v6proplr(x,x,-1,cwttmp(:,:),cwtr2lnew(:,:,j))   
    enddo
    
    call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2lnew(:,:,1),cwttmp(:,:))
    call v6propl(xl,-1,cwttmp(:,:),cwtr2lnew(:,:,0))
   
    newp=abs(real( sum(conjg(cwtbeginnew(:,:))*cwtr2lnew(:,:,0)) )) 	  
    oldp=pathcoreval !abs(real( sum(conjg(cwtbegin(:,:))*cwtr2l(:,:,0)) ))   
    
    !tmp1test1=abs(real( sum(conjg(cwtbegin(:,:))*cwtr2l(:,:,0)) ))   
    !write(6,'(''lr tmp1-oldp='',t30,i5,3(f15.7,1x))') ,nblocknow, abs(tmp1test1-oldp),tmp1test1,oldp     
    !if (abs(tmp1test1-oldp)>=1.0d-7) stop 
    
    rn=randn(1)  
    if (rn(1).lt.(newp/oldp)) then  
        iclr=iclr+1 
        pathcoreval=newp
    ! l2r        
        call v6propr(xl,-1,cwtbeginnew,cwtl2rnew(:,:,0))
        do i=1,nchorizo-1
	        x=ristra(i)%x
           
            call vlsprop(spinorbit,-1,ristra(i-1)%x,x,cwtl2rnew(:,:,i-1),cwttmp(:,:))       
	        call v6proplr(x,x,-1,cwttmp(:,:),cwtl2rnew(:,:,i))
   
        enddo   

        call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2rnew(:,:,nchorizo-1),cwttmp(:,:))
        call v6propl(xr,-1,cwttmp(:,:),cwtl2rnew(:,:,nchorizo))
              
        cwtbegin(:,:)=cwtbeginnew(:,:)
        cwtend(:,:)=cwtendnew(:,:) 
        cwtl2r=cwtl2rnew
        cwtr2l=cwtr2lnew
        iplo=ipl
        jplo=jpl
        ipro=ipr
        jpro=jpr
        call inputordlr(iplo,jplo,ipro,jpro)   
        
        
        
     !!!!!!!!!!!! check   
        !tmp1=abs(real(sum(conjg(cwtbegin(:,:))*cwtr2l(:,:,0))))
        !tmp2=abs(real(sum(conjg(cwtl2r(:,:,nchorizo))*cwtend(:,:))))
        !tmp=(tmp1+tmp2)/2  
        !tmp3=abs(  (pathcoreval-tmp)/( (pathcoreval+tmp)/2  )        )*100
        !write(6,*) 'pathcoreval check L and R order:', oldp,pathcoreval, tmp,tmp1,tmp2,tmp3
        !if (  tmp3 >= 0.1   ) then
        !    write(6,*) '################## Caution ###################'
        !    write(6,*) 'pathcoreval check L and R order failed % err = ', tmp3
        !    write(6,*) '###############################################'
        !    !stop
        !endif
     !!!!!!!!!!!!!!!!!!         

        
    endif
    
    
    !write (6,*) 'stepchorizoinitlr   cwtend = ', cwtend
    
    

    return
    end subroutine stepchorizoinitlr   
   
    subroutine stepstd(nblocknowin) ! for non-constant trial wavefunction
    use estimator
    use estimatorristra
    use wavefunction
    use random
    use brut
    use mympi
    real(kind=r8) :: rn(1),rn1(1),rm,rmsq,rmn,rmnsq,rmp,rmpsq
    real(kind=r8) :: v(0:nchorizo),vnn(0:nchorizo),vem(0:nchorizo) &
                    ,xtmp(3,npart),x0tot(3,npart,0:nchorizo)
    real(kind=r8), dimension(0:nrhobin) :: rhodistout
    integer(kind=i4) :: i,l,ixd,lmax,icheck,ixtemp,nblocknowin
    real(kind=r8) :: vn,vd,vnke,vnpev6,vnpeem,vnpels,rf,signrf,vnum,vdenom &
                    ,sum1,sum2,vnumsum,vdenomsum,rfsum,vnkesum,vnpev6sum,vnpeemsum,vnpelssum
    real(kind=r8) :: psi20,psi2n
    real(kind=r8) :: rb,r1,r2,r3,r23 ! set percentage upper bar of moves for different types of moves.  
    logical :: imp
    character(len=120) :: filename
    integer(kind=i4) :: myunit,myunit1,myunit2
    integer(kind=i4) :: iplnow(npair),jplnow(npair),iprnow(npair),jprnow(npair)	
    real(kind=r8) :: x0tot1d(3*npart*(nchorizo+1)),x0tot1dnow(3*npart*(nchorizo+1))
    logical :: hpsipass,neqfin
    complex(kind=r8) :: cwtl2rm(0:nspin-1,nisospin),cwtr2lm(0:nspin-1,nisospin)

    
    rb=rbisect ! 0.9   bisection rate(reptation if enabled is in it.)
    r1=rmov1 ! 0.95 move 1 by 1    
    r2=rmov2 ! 0.967 !move all beads  
    r3=rmov3 ! 0.984 ! shift beads
    r23=rmov23 ! 1 move all + shift 
    
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
	  
       
        !rb=1 ! 0.9   bisection rate(reptation if enabled is in it.)
        !r1=0.95 ! move 1 by 1    
        !r2=0.967 ! move all beads  
        !r3=0.984 ! shift beads
        !r23=1  ! move all + shift     
     
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
            if ((rn(1)>rb).and.(rn(1)<=r1)) then !if ((rn(1)>rb).and.(rn(1)<=r1)) call stdmove1fastv6(0,nchorizo)   ! move 1 by 1    
                rn1=randn(1)
                imp=.true.          
                call stdmove1fastv6flex(0,nchorizo,2*nint(rn1(1))-1) 
                !call stdmove1fastv6flexa(0,nchorizo,1,imp) 
            
                !if (rn1(1)<0.35) call stdmove1bisect1v6flex(0,nchorizo,1)    
                !if ((rn1(1)>=0.35).and.(rn1(1)<0.7)) call stdmove1bisect1v6flex(0,nchorizo,-1)    
                !if ((rn1(1)>=0.7).and.(rn1(1)<0.85)) call stdmove1fastv6flex(0,nchorizo,1)
                !if ((rn1(1)>=0.85)) call stdmove1fastv6flex(0,nchorizo,-1)   
            endif        
            if ((rn(1)>r1).and.(rn(1)<=r2))  call stdmove2v6(0,nchorizo)   ! move all beads  
            if ((rn(1)>r2).and.(rn(1)<=r3))  call stdmove3v6(0,nchorizo)  ! shift beads
            if ((rn(1)>r3).and.(rn(1)<=r23))  call stdmove23v6(0,nchorizo)  ! move all + shift     
        endif
     
        !call bisectv6improvew11    
  
    else    
        call stdvmcmove2v6
    endif
    
    call stepchorizoinitlr ! this part may be combined into bisection part to inrease the speed by just a little bit (maynot really worth to do it). 2020-05-10, CR.      
      
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

    ! check correlation steps. 
    ! generate a equilibrated config, then read in it, and set a big nstep, let nstepdecor=nstep, then run to check correlation time.
    if (icorrchk) then    
        if (myrank().eq.0) then  
            if (  mod(icstepstd,50).eq.0  ) write(6,*) 'icstepstd=', icstepstd,nstep
            if (icstepstd <= nstep) then   ! here nstepdecor usually set the same as nstep. 
	            call funcrmsq(ristra(nchorizomid)%x,xin(icstepstd),rmsq)  
                call funcrmsq(ristra(0)%x,x1in(icstepstd),rmsq)
                call funcrmsq(ristra(nchorizo)%x,x2in(icstepstd),rmsq)
                !call hpsi(ristra(0),ristra(nchorizo),vn,vd,vnke,vnpev6,vnpeem,psi20,psi2n,rf,icheck)
                !yin(icstepstd)=vn
	            if (icstepstd.eq.nstep) then
	                ixd=icstepstd
	                lmax=nstep                   
                    call corrchk(xin,x1in,x2in,ixd,corr,corr1,corr2,lmax)       
                endif   
            endif
        endif
    endif
    
    ! calculations	  
    if (mod(icstepstd,nstepdecor).eq.0) then  
        
        call barrier ! need to add barrier make all the cores sync here. Otherwise may have racing problem.
        
        
        if ((nblocknow.eq.1).and.(myrank().eq.0)) then
            time0=mpi_wtime() 
            if (icstepstd.eq.nstepdecor) then
                time00=0
                time5tot=0
            endif   
        endif
        
        !write(6,*) 'icstepstd=',icstepstd    
		   
        ! do not need to calculate elocal too often.  
        ice=ice+1
        
        call hpsi(ristra(0),ristra(nchorizo),vn,vd,vnke,vnpev6,vnpeem,vnpels,psi20,psi2n,rf,icheck)
        
    
        if ((nblocknow.eq.1).and.(myrank().eq.0)) time1=mpi_wtime()        

        !write(6,'( ''stepstd hpsicheck'', 8(g15.7,1x) )') vn,vd,vnke,vnpev6,vnpeem,psi20,psi2n,rf
          
! check	hpsi begin. can check to see if just let rank 0 allocate those arrays works or not.
! no need to check it those way, this may be slow.        
             
        ! ixtemp=3*npart*(nchorizo+1)  
        
        call chorizoallout(x0tot,iplo,jplo,ipro,jpro)   
        x0tot1d=reshape(x0tot,(/3*npart*(nchorizo+1)/))
             
! chekhpsi
        call checkhpsi(vn,vd,vnke,vnpev6,vnpeem,vnpels,psi20,psi2n,rf,icheck,iplo,jplo,ipro,jpro,x0tot1d,nblocknow,ipath,hpsipass)        
!
                
        if ((nblocknow.eq.1).and.(myrank().eq.0)) time2=mpi_wtime()    
        
        if (irep==2) then ! Only for PIMC we need to write out path.
            if (nblocknow > neq) then   ! write out path  
                neqfin=.true.
                call writepath(ipath,iplo,jplo,ipro,jpro,x0tot1d,neqfin)
                !call writepathgs(ipath,iplo,jplo,ipro,jpro,x0tot1d,neqfin)
                ipathcountold=ipath
                ipath=ipath+1  
            else   
                neqfin=.false.
                call writepath(ipathtmp,iplo,jplo,ipro,jpro,x0tot1d,neqfin)
                !call writepathgs(ipathtmp,iplo,jplo,ipro,jpro,x0tot1d,neqfin)
                ipathtmp=ipathtmp+1 
            endif          
        endif
        
        if (nrepmax /= 0) then ! stat on reptation. only for PIMC.
            call addall(repstepct1,repstepct1sum)
            call addall(repstepct2,repstepct2sum)   
            sum1=sum(repstepct1sum)
            sum2=sum(repstepct2sum) 
            if (myrank().eq.0) then
                open(newunit=myunit,form='formatted',file='reptationMON.txt',position='rewind')
            ! check reptation number of steps distributions.
                do i=0,irepstepmonmax
                    write(myunit,'( i10 , i10, f15.7, i10, f15.7 )') i,repstepct1sum(i),repstepct1sum(i)/sum1 &
                                                                ,repstepct2sum(i),repstepct2sum(i)/sum2 
                enddo 
                close(myunit)     
                open(newunit=myunit1,form='formatted',file='reparrowhist.txt',position='append')
            ! check reptation direction history.
                do i=icrephistlast,icrephist
                    write(myunit1,'( i10,1x,i10 )') i,ireparrowhistory(i)
                enddo 
                icrephistlast=icrephist+1
                close(myunit1)       
            endif
        endif          
        
!
  
        if ((nblocknow.eq.1).and.(myrank().eq.0)) time5=mpi_wtime() 
  
        vnum=vn
        vdenom=vd
        ! update the v beads   
        
        
         !write (6,*) 'stepstd     cwtend = ', cwtend
        
        
        do i=0,nchorizo
            call vtot(i,v(i),vnn(i),vem(i))
        enddo
        v=v/abs(rf) ! in the end need to be devided by average of vdenom and consider its error.
        vnn=vnn/abs(rf)
        vem=vem/abs(rf)
        call addvalristra(1,v) 
        call addvalristra(2,vnn) 
        call addvalristra(3,vem)
        
        if ((nblocknow.eq.1).and.(myrank().eq.0)) time3=mpi_wtime() 
       
        call addval(3,vnum,1.0_r8)
        call addval(4,vdenom,1.0_r8)
        call addval(5,rf,1.0_r8)
               
        call getcwtlrm(nchorizomid,cwtl2rm,cwtr2lm)
        call samplerhodistrbt(ristra(nchorizomid)%x,nchorizomid,cwtl2rm,cwtr2lm,rf,rhodistout,rm,rmsq,rmn,rmnsq,rmp,rmpsq)
        call addval(6,rmsq,1.0_r8)
        call addval(7,rm,1.0_r8)
        call addval(8,rhodistout(0),1.0_r8) 
        
        
        
        call addval(13,rmnsq,1.0_r8)
        call addval(14,rmn,1.0_r8)      
        call addval(15,rmpsq,1.0_r8)
        call addval(16,rmp,1.0_r8)
        
        call addval(17,vnke,1.0_r8)
        call addval(18,vnpev6+vnpeem+vnpels,1.0_r8)
        call addval(19,vnpev6,1.0_r8)
        call addval(20,vnpeem,1.0_r8)  
        call addval(21,vnpels,1.0_r8) 
        
        
        call addall(vnum,vnumsum)
        call addall(vdenom,vdenomsum)
        call addall(rf,rfsum)
        call addall(vnke,vnkesum)
        call addall(vnpev6,vnpev6sum)
        call addall(vnpeem,vnpeemsum)
        call addall(vnpels,vnpelssum)
        
        
        if (myrank().eq.0) then     
            open (newunit=myunit2,FILE=file_values,FORM='FORMATTED',position='append')    ! may be used for data ayalysis. 
                write(myunit2,'(15(G15.7,2x))') vnumsum,vdenomsum,vnkesum,vnpev6sum,vnpeem &         
                     ,vnumsum/nproc(),vdenomsum/nproc(),rfsum/nproc(),vnkesum/nproc(),vnpev6sum/nproc(),vnpeemsum/nproc()
            close (myunit2)   
        endif            
        

!!! this part can be removed.          
        if ((nblocknow.eq.1).and.(myrank().eq.0)) then
            time4=mpi_wtime() 
             write (6,*) '-------------------------------------------------------------------------------'
             write (12,*) '-------------------------------------------------------------------------------'
             write (6,'(/,''Time for one calculation ='',f15.7,'' seconds'')') time4-time0
             write (12,'(/,''Time for one calculation ='',f15.7,'' seconds'')') time4-time0
             write (6,'(/,''Time for one block all calculation est ='',f15.7,'' seconds'')') (time4-time0)*(nstep/nstepdecor)
             write (12,'(/,''Time for one block all calculation est ='',f15.7,'' seconds'')') (time4-time0)*(nstep/nstepdecor)
             write (6,'(/,''Time for one hpsi ='',f15.7,'' seconds'')') time1-time0
             write (12,'(/,''Time for one hpsi ='',f15.7,'' seconds'')') time1-time0       
             write (6,'(/,''Time for one hpsi check ='',f15.7,'' seconds'')') time2-time1
             write (12,'(/,''Time for one hpsi check ='',f15.7,'' seconds'')') time2-time1 
             write (6,'(/,''Time for one vtot ='',f15.7,'' seconds'')') time3-time5
             write (12,'(/,''Time for one vtot ='',f15.7,'' seconds'')') time3-time5 
             write (6,'(/,''Time for one MPI path IO ='',f15.7,'' seconds'')') time5-time2
             write (12,'(/,''Time for one MPI path IO ='',f15.7,'' seconds'')') time5-time2  
             write (6,'(/,''Time for one the rest ='',f15.7,'' seconds'')') time4-time5
             write (12,'(/,''Time for one the rest ='',f15.7,'' seconds'')') time4-time5
             write (6,*) '-------------------------------------------------------------------------------'
             write (12,*) '-------------------------------------------------------------------------------'             
             time00=time00+time4-time0
             time5tot=time5tot+time5-time2
        endif
        
        
        
        
        
!!!   
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
    
    if (nrepmax.eq.0) return
    
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
    !     write (19,'(i10)') i
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
    !     write (19,'(i10)') i
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
    call v6propr(xr,-1,cwtendnew,cwtr2lnew(:,:,nchorizo))
    do j=nchorizo-1,1,-1
	    x=ristranew(j)%x
	    call v6proplr(x,x,-1,cwtr2lnew(:,:,j+1),cwtr2lnew(:,:,j))
    enddo
    call v6propl(xl,-1,cwtr2lnew(:,:,1),cwtr2lnew(:,:,0))
   
    newp=abs(real( sum(conjg(cwtbeginnew(:,:))*cwtr2lnew(:,:,0)) )) 	  
    oldp=pathcoreval !abs(real( sum(conjg(cwtbegin(:,:))*cwtr2l(:,:,0)) ))   
   
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
        call v6propr(xl,-1,cwtbeginnew,cwtl2rnew(:,:,0))
        do i=1,nchorizo-1
	    x=ristra(i)%x
	    call v6proplr(x,x,-1,cwtl2rnew(:,:,i-1),cwtl2rnew(:,:,i))
        enddo   
        call v6propl(xr,-1,cwtl2rnew(:,:,nchorizo-1),cwtl2rnew(:,:,nchorizo))
        cwtbegin(:,:)=cwtbeginnew(:,:)
        cwtend(:,:)=cwtendnew(:,:) 
        cwtl2r=cwtl2rnew
        cwtr2l=cwtr2lnew
        
        pathcoreval=newp
                
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
   
    subroutine reptation1 ! move 1 bead, small importance sampling is added. nrepmax=1 case.
    use wavefunction
    use v6stepcalc
    use brut 
    use random
    use mympi
    use estimator
    integer(kind=i4) :: i,j
    real(kind=r8) :: xr(3,npart),xl(3,npart),x(3,npart),gauss(3,npart),xr1(3,npart),xl1(3,npart)
    real(kind=r8) ::  newp,oldp,rn(1),newp1,oldp1 
   
    if (nrepmax.eq.0) return
    
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
        call v6propr(xl,-1,cwtbeginnew,cwtl2rnew(:,:,0))
        
        do i=1,nchorizo-1
	        x=ristranew(i)%x
	        call v6proplr(x,x,-1,cwtl2rnew(:,:,i-1),cwtl2rnew(:,:,i))
            
        enddo 
        cwtl2rnew1(:,:,0:nchorizo-1)=cwtl2rnew(:,:,0:nchorizo-1)
        xr=ristranew(nchorizo)%x 
        xr1=ristranew1(nchorizo)%x 
        call v6propl(xr,-1,cwtl2rnew(:,:,nchorizo-1),cwtl2rnew(:,:,nchorizo)) 
        call v6propl(xr1,-1,cwtl2rnew1(:,:,nchorizo-1),cwtl2rnew1(:,:,nchorizo))     
        call corop(xr,ipro,jpro,cwtendnew)  
        call corop(xr1,ipro,jpro,cwtendnew1)  

        newp=abs(real( sum(conjg(cwtl2rnew(:,:,nchorizo))*cwtendnew(:,:)) )) 	 
        newp1=abs(real( sum(conjg(cwtl2rnew1(:,:,nchorizo))*cwtendnew1(:,:)) )) 
    
        oldp=pathcoreval !abs(real( sum(conjg(cwtbegin(:,:))*cwtr2l(:,:,0)) ))   
    
        xl1=2*ristra(1)%x-ristra(0)%x
        call corop(xl1,iplo,jplo,cwtbegin1) 
        call v6propl(xl1,-1,cwtr2l(:,:,1),cwtr2l1(:,:,0))    
    
        oldp1=abs(real( sum(conjg(cwtbegin1(:,:))*cwtr2l1(:,:,0)) ))   
    
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
                call v6propr(xr,-1,cwtendnew,cwtr2lnew(:,:,nchorizo))
                do j=nchorizo-1,1,-1
	                x=ristra(j)%x
	                call v6proplr(x,x,-1,cwtr2lnew(:,:,j+1),cwtr2lnew(:,:,j))
                enddo
                call v6propl(xl,-1,cwtr2lnew(:,:,1),cwtr2lnew(:,:,0))  
                cwtbegin(:,:)=cwtbeginnew(:,:)
                cwtend(:,:)=cwtendnew(:,:) 
                cwtl2r=cwtl2rnew 
                cwtr2l=cwtr2lnew  
                
                pathcoreval=newp
                
            else
                ristra(0:nchorizo)=ristranew1(0:nchorizo)   
            ! r2l 
                call v6propr(xr1,-1,cwtendnew1,cwtr2lnew1(:,:,nchorizo))
                do j=nchorizo-1,1,-1
	                x=ristra(j)%x
	                call v6proplr(x,x,-1,cwtr2lnew1(:,:,j+1),cwtr2lnew1(:,:,j))
                enddo
                call v6propl(xl,-1,cwtr2lnew1(:,:,1),cwtr2lnew1(:,:,0))     
                cwtbegin(:,:)=cwtbeginnew(:,:)
                cwtend=cwtendnew1
                cwtl2r=cwtl2rnew1  
                cwtr2l=cwtr2lnew1  
                
                pathcoreval=newp1
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
        call v6propr(xr,-1,cwtendnew,cwtr2lnew(:,:,nchorizo))
        do j=nchorizo-1,1,-1
            x=ristranew(j)%x
            call v6proplr(x,x,-1,cwtr2lnew(:,:,j+1),cwtr2lnew(:,:,j))
        enddo
        cwtr2lnew1(:,:,1:nchorizo)=cwtr2lnew(:,:,1:nchorizo)
        call v6propl(xl,-1,cwtr2lnew(:,:,1),cwtr2lnew(:,:,0))      
        call v6propl(xl1,-1,cwtr2lnew(:,:,1),cwtr2lnew1(:,:,0))    

        call corop(xl,iplo,jplo,cwtbeginnew)
        call corop(xl1,iplo,jplo,cwtbeginnew1)   
    
        newp=abs(real( sum(conjg(cwtbeginnew(:,:))*cwtr2lnew(:,:,0)) ))  
        newp1=abs(real( sum(conjg(cwtbeginnew1(:,:))*cwtr2lnew1(:,:,0)) )) 
    
        oldp=pathcoreval !abs(real( sum(conjg(cwtbegin(:,:))*cwtr2l(:,:,0)) ))   
    
        xr1=2.*ristra(nchorizo-1)%x-ristra(nchorizo)%x
        call corop(xr1,ipro,jpro,cwtend1)   
        call v6propl(xr1,-1,cwtl2r(:,:,nchorizo-1),cwtl2r1(:,:,nchorizo))   
    
        oldp1=abs(real( sum(conjg(cwtl2r1(:,:,nchorizo))*cwtend1(:,:)) ))   
    
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
                call v6propr(xl,-1,cwtbeginnew,cwtl2rnew(:,:,0))
                do i=1,nchorizo-1
	                x=ristra(i)%x
	                call v6proplr(x,x,-1,cwtl2rnew(:,:,i-1),cwtl2rnew(:,:,i))
                enddo   
                call v6propl(xr,-1,cwtl2rnew(:,:,nchorizo-1),cwtl2rnew(:,:,nchorizo))       
                cwtbegin(:,:)=cwtbeginnew(:,:)
                cwtend(:,:)=cwtendnew(:,:) 
                cwtl2r=cwtl2rnew 
                cwtr2l=cwtr2lnew 
                
                pathcoreval=newp
  
            else
     
                ristra(0:nchorizo)=ristranew1(0:nchorizo)   
            ! l2r
                call v6propr(xl1,-1,cwtbeginnew1,cwtl2rnew1(:,:,0))
                do i=1,nchorizo-1
	                x=ristra(i)%x
	                call v6proplr(x,x,-1,cwtl2rnew1(:,:,i-1),cwtl2rnew1(:,:,i))
                enddo   
                call v6propl(xr,-1,cwtl2rnew1(:,:,nchorizo-1),cwtl2rnew1(:,:,nchorizo))       
                cwtbegin=cwtbeginnew1
                cwtend(:,:)=cwtendnew(:,:) 
                cwtl2r=cwtl2rnew1 
                cwtr2l=cwtr2lnew1 
                
                pathcoreval=newp1
  
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
   
    subroutine bisect
    use random
    use math
    real(kind=r8) :: rn(1),rn2(1)
    integer(kind=i4) :: ileftbisect,irightbisect,mmaxnow,mmaxnow2,n,p,i,ipick &
                        ,ileftbisect2,irightbisect2,ipick2,nbisectnow,nbisectnow2
    
    icbisecttot=icbisecttot+1   
   
    if (nbisect.le.nchorizo) then
           
        rn2=randn(1)
        if ( rn2(1) <= 0.9 ) then
            mmaxnow=mmax   
        else
            if (mmax<3) then
                mmaxnow=mmax 
            else
                mmaxnow=mmax-1
            endif
        endif  
        nbisectnow=2**mmaxnow
        
        call bisectpicksliceswide(ileftbisect,irightbisect,nbisectnow) ! pick the ileftbisect here
    
            !write(12,*) 'ileft,iright=', ileftbisect,irightbisect,nbisect,nbisecthalf
            !write(6,*) 'ileft,iright=', ileftbisect,irightbisect,nbisect,nbisecthalf
    
        if (ileftbisect <= (nchorizo-nbisectnow)) then
     
            !call bisectv6improve(ileftbisect)
            !call bisectv6a(ileftbisect)
        
        
            !call bisectv6b(ileftbisect)
            call bisectv6bflex(ileftbisect,mmaxnow)
        
            !write(6,*) 'bisectv6improve called =', ileftbisect
     
        else
            
            
            !call bisectv6lb
            !call bisectv6rb  ! those two pair seem to give rank611 partial stuck config.
          
            !rn=randn(1)
            !if (rn(1) < 0.5) then
            !    call bisectv6lbflex(mmaxnow)  
            !    call bisectv6rbflex(mmaxnow)    
            !else
            !    call bisectv6rbflex(mmaxnow)    
            !    call bisectv6lbflex(mmaxnow)    
            !endif
                       
            !rn=randn(1)
            !if (rn(1) < 0.5) then  
                !call bisectv6bextraflex(ileftbisect,mmaxnow)     
            !else  
                

            
            
            
            
                mmaxnow2=mmaxnow-1   ! the smaller bisection.
                n=3
                rn=randn(1)
                nbisectnow2=2**mmaxnow2
                if (ileftbisect <= nchorizo-nbisectnow2) then     
                    do i=1,n
                        p=int(rn(1)*size(p3map,DIM=2))+1
                        ipick=p3map(i,p)
                        select case (ipick)
                        case (1)
                            call bisectv6bflex(ileftbisect,mmaxnow2)
                        case (2)
                            call bisectv6rbflex(mmaxnow2)  
                        case (3)
                            call bisectv6lbflex(mmaxnow2)        
                        case default 
                            write (6,'(''ipick value does not match the move case # '')') ipick,ileftbisect
                            call abort 
                        end select
                    enddo
                else
                    do i=1,n
                        p=int(rn(1)*size(p3map,DIM=2))+1
                        ipick=p3map(i,p)
                        select case (ipick)
                        case (1)
                            call bisectv6bflex(nbisect-(nchorizo-ileftbisect)-nbisectnow2,mmaxnow2)  
                        case (2)
                            call bisectv6rbflex(mmaxnow2)  
                        case (3)
                            call bisectv6lbflex(mmaxnow2)        
                        case default 
                        write (6,'(''ipick value does not match the move case # '')') ipick,ileftbisect
                        call abort 
                        end select
                    enddo                    
                endif                 
   
                
                
                
                
                
                
                
                
                
                
                
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
        !!! reptation, consider a bisection-ike reptation maybe.
   

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
 
    subroutine bisectv6bflex(ileftbisect,mmaxnow) 
    ! based on v6a version. move beads first, then do an overall judge.
    use wavefunction
    use estimator
    use random
    use v6stepcalc
    use brut
    integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
					    ,k,bilvl &
	                    ,dtexpt,mmaxnow,nbisectnow
    real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    real(kind=r8) :: xmid(3,npart),xr(3,npart),xl(3,npart)
    real(kind=r8) :: x(3,npart),xnew(3,npart)
    real(kind=r8) :: tmp,tmp1,tmp2,tmp3,tmp1test1
    real(kind=r8) :: oldlv(0:mmaxnow),newlv(0:mmaxnow)
    complex(kind=r8) :: cwtl2rtmp1(0:nspin-1,nisospin,0:2**mmaxnow),cwtr2ltmp1(0:nspin-1,nisospin,0:2**mmaxnow) &
	                    ,cwtl2rtmp2(0:nspin-1,nisospin,0:2**mmaxnow),cwtr2ltmp2(0:nspin-1,nisospin,0:2**mmaxnow) &
	                    ,cwtbegintmp(0:nspin-1,nisospin),cwtendtmp(0:nspin-1,nisospin),cwttmp(0:nspin-1,nisospin)
    logical :: reject	
 
    if (mmaxnow.eq.mmax) ibisecttot=ibisecttot+1 
    if (mmaxnow.eq.mmax-1) ibisecthalftot=ibisecthalftot+1
   
    nbisectnow=2**mmaxnow
    oldristra(0:nbisectnow)=ristra(ileftbisect:ileftbisect+nbisectnow)
   
    newristra(0)=oldristra(0)
    newristra(nbisectnow)=oldristra(nbisectnow)

    cwtbegintmp(:,:)=cwtl2r(:,:,ileftbisect)
    cwtendtmp(:,:)=cwtr2l(:,:,ileftbisect+nbisectnow)  
   
    ! tmp1 deal with oldlv
    cwtl2rtmp1(:,:,0:nbisectnow)=cwtl2r(:,:,ileftbisect:ileftbisect+nbisectnow)
    cwtr2ltmp1(:,:,0:nbisectnow)=cwtr2l(:,:,ileftbisect:ileftbisect+nbisectnow)
    ! tmp2 deal with newlv   
    cwtl2rtmp2(:,:,0:nbisectnow)=cwtl2r(:,:,ileftbisect:ileftbisect+nbisectnow)
    cwtr2ltmp2(:,:,0:nbisectnow)=cwtr2l(:,:,ileftbisect:ileftbisect+nbisectnow)
   
    oldlv=1
    newlv=1   
    
    tau=dt*nbisectnow ! total time for the bisection slice.

    ! up to here, newristra has been loaded.
   
    ! judge if we accept the newristra or not.
    ! in this v6 case, v is not justpotential v, it should be the value of < ... >. lv.
        
    do bilvl=1,mmaxnow  ! besection.  level 1 to N. mmaxnow = the total level N.	 
        ! move the current level   
        sigmamid=sqrt(hbar*tau/2**bilvl)  
	    !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	    do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		    imid=2**(mmaxnow-bilvl)+i*2**(mmaxnow-bilvl+1)
		    il=imid-2**(mmaxnow-bilvl)
		    ir=imid+2**(mmaxnow-bilvl)
		    gauss(:,:)=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		    xmid=(newristra(il)%x+newristra(ir)%x)/2
		    newristra(imid)%x(:,:)=xmid(:,:)+gauss(:,:)  
            !write(12,*) 'bilvl il ir imid sigmamid =',bilvl, il,ir,imid,sigmamid
            !write(12,*) 'gauss =', gauss
        enddo         
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!    
    ! then judge the last level
	    bisecttot(mmaxnow)=bisecttot(mmaxnow)+1	   
	    lmax=2**mmaxnow-1 
	    jd=1 ! interval
	    dtexpt=-1 ! for the R L part, that's why the extra -1. dt=2**dtexp
        
        
	    do l=1,lmax
		    j=l*jd
		    jp=(l-1)*jd
		    !xold=oldristra(j)%x
            xnew(:,:)=newristra(j)%x(:,:)
               
            call vlsprop(spinorbit,-1,newristra(jp)%x,xnew,cwtl2rtmp2(:,:,jp),cwttmp(:,:))
		    ! call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(:,:,jp),cwtl2rtmp1(:,:,j))    ! should be lr, not rl I believe.
		    call v6proplr(xnew,xnew,dtexpt,cwttmp(:,:),cwtl2rtmp2(:,:,j)) 
    
        enddo
          
	    jmax=lmax*jd
        
        call vlsprop(spinorbit,-1,newristra(jmax)%x,newristra(jmax+1)%x,cwtl2rtmp2(:,:,jmax),cwttmp(:,:))       
	    newlv(mmaxnow)=abs(real( sum(conjg(cwttmp(:,:))*cwtendtmp(:,:)) )) 	
	    oldlv(mmaxnow)=pathcoreval !abs(real( sum(conjg(cwtl2rtmp1(:,:,jmax))*cwtendtmp(:,:)) ))	
        
        
        !tmp1test1=abs(real( sum(conjg(cwtl2rtmp1(:,:,jmax))*cwtendtmp(:,:)) ))	
        !write(6,'(''bisect tmp1-oldp='',t30,i5,3(f15.7,1x))') ,nblocknow, abs(tmp1test1-oldlv(mmaxnow)),tmp1test1,oldlv(mmaxnow) 
        !if (abs(tmp1test1-oldlv(mmaxnow))>=1.0d-3) stop 
        
        
	    tmp=log(newlv(mmaxnow))-log(oldlv(mmaxnow)) 
  	    !tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
	    rn=randn(1)	   
    ! can check if l2r and r2l give	the same lv.		    
            !write(12,*) 'tmp ln rn=',tmp,log(rn(1))
	    if (log(rn(1)).lt.tmp) then
	        reject=.false. 
	        bisectcount(mmaxnow)=bisectcount(mmaxnow)+1	 
	    else
	        reject=.true. 
            !write(6,*) 'rn=',rn(1),log(rn(1))
            !write(6,*) 'tmp=',tmp,newlv(mmaxnow),oldlv(mmaxnow),log(newlv(mmaxnow)),log(oldlv(mmaxnow)) 	        
        endif 
           
    if (reject) then
	    if (mmaxnow.eq.mmax) call addval(2,0.0_r8,1.0_r8)
        if (mmaxnow.eq.mmax-1) call addval(9,0.0_r8,1.0_r8)
    else   
        if (mmaxnow.eq.mmax) then
            ibisect=ibisect+1 
	        call addval(2,1.0_r8,1.0_r8)  
        endif
        if (mmaxnow.eq.mmax-1) then
            ibisecthalf=ibisecthalf+1 
	        call addval(9,1.0_r8,1.0_r8)  
        endif        
        ristra(ileftbisect:ileftbisect+nbisectnow)=newristra(0:nbisectnow) 
        xl=ristra(0)%x     
	    xr=ristra(nchorizo)%x     
    ! update cwtl2r, cwtr2l.	 
    ! update cwtl2r for all beads. this can be optimized in the future.
     
	    cwtl2r(:,:,ileftbisect+1:ileftbisect+nbisectnow-1)=cwtl2rtmp2(:,:,1:jmax) ! here jmax should already be 2**mmaxnow-1 
	    if ((ileftbisect+nbisectnow).le.nchorizo-1) then
            do k=ileftbisect+nbisectnow,nchorizo-1
	            x=ristra(k)%x
                call vlsprop(spinorbit,-1,ristra(k-1)%x,x,cwtl2r(:,:,k-1),cwttmp(:,:))
	            call v6proplr(x,x,-1,cwttmp(:,:),cwtl2r(:,:,k))
                
            enddo

            call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,nchorizo-1),cwttmp(:,:))
	        call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))
        else !ileftbisect+nbisectnow=nchorizo    

            call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,nchorizo-1),cwttmp(:,:))
	        call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))   
        endif
   
    ! update cwtr2l for all the beads
	    do k=ileftbisect+nbisectnow-1,1,-1
	        x=ristra(k)%x        
            call vlsprop(spinorbit,1,x,ristra(k+1)%x,cwtr2l(:,:,k+1),cwttmp(:,:))
	        call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,k))
            
        enddo
         
        call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
        call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))  	

     ! update pathcoreval   
        pathcoreval=newlv(mmaxnow)
        
        
     !!!!!!!!!!!! check  
        
        
       
        !tmp1=abs(real(sum(conjg(cwtbegin(:,:))*cwtr2l(:,:,0))))
        !
        !
        !!k=0
        !!call vlsprop(spinorbit,1,ristra(k)%x,ristra(k+1)%x,cwtr2l(:,:,k+1),cwttmp(:,:))
        !!tmp2=abs(real(sum(conjg(cwtl2r(:,:,k))*cwttmp(:,:))))
        !
        !
        !
        !!call corop(xr,ipro,jpro,cwtendtmp)
        !!write (6,*) 'bisectv6bflex    cwtend = ', cwtend
        !!write (6,*) 'cwtendtmp = ',cwtendtmp
        !!write (6,*) 'cwtend diff = ',cwtendtmp-cwtend
        !
        !tmp2=abs(real(sum(conjg(cwtl2r(:,:,nchorizo))*cwtend(:,:))))
        !tmp=(tmp1+tmp2)/2  
        !tmp3=abs(  (pathcoreval-tmp)/( (pathcoreval+tmp)/2  )        )*100
        !write(6,*) 'ileftbisect,rborder:', ileftbisect,ileftbisect+nbisectnow
        !write(6,*) 'pathcoreval check bisect:',oldlv(mmaxnow), pathcoreval, tmp,tmp1,tmp2,tmp3
        !if (  tmp3 >= 0.1   ) then
        !    write(6,*) '################## Caution ###################'
        !    write(6,*) 'pathcoreval check bisect failed % err = ', tmp3
        !    write(6,*) '###############################################'
        !    stop
        !endif
    
     !!!!!!!!!!!!!!!!!!   
        
            
        
        
        
    endif  
    !rejectbisectv6improve=reject
    return  
    end subroutine bisectv6bflex           
     
    subroutine bisectv6lbflex(mmaxnow)
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
    real(kind=r8) :: tmp,tmp1,tmp2,tmp3,tmp1test1
    real(kind=r8) :: oldlv(0:mmax),newlv(0:mmax)
    complex(kind=r8) :: cwtl2rtmp1(0:nspin-1,nisospin,0:2**mmaxnow) &
	                    ,cwtl2rtmp2(0:nspin-1,nisospin,0:2**mmaxnow) &
	                    ,cwtendtmp(0:nspin-1,nisospin) &
                        ,cwtbegintmp2(0:nspin-1,nisospin),cwttmp(0:nspin-1,nisospin)     
    logical :: reject	

    ibisecttotl=ibisecttotl+1  
    
    nbisectnow=2**mmaxnow
    ileftbisect=0
    irightbisect=ileftbisect+nbisectnow

    tau=dt*nbisectnow ! total time for the bisection slice. nisectnow is usually nbisecthalf
    sigmamid=sqrt(2*hbar*tau)  

    oldristra(0:nbisectnow)=ristra(ileftbisect:ileftbisect+nbisectnow)
 
    ! deal with the two ends first, then do the classic bisection.   
   
    newristra(nbisectnow)=oldristra(nbisectnow)   
    cwtendtmp(:,:)=cwtr2l(:,:,ileftbisect+nbisectnow)
 
    newristra(0)%x=oldristra(nbisectnow)%x+sigmamid*reshape(gaussian(3*npart),(/3,npart/))
   
    x0new(:,:)=newristra(0)%x(:,:)
   
    ! i,jpro,i,jplo, the o means old are actually the current order.
   
    call corop(x0new,iplo,jplo,cwtbegintmp2)     
   
    !oldlv(0)=1.
    !newlv(0)=1. ! dummy

    cwtl2rtmp1(:,:,0:nbisectnow)=cwtl2r(:,:,ileftbisect:ileftbisect+nbisectnow)
    call v6propr(x0new,-1,cwtbegintmp2,cwtl2rtmp2(:,:,0)) 
   
    ! prepare the further new bisection positions
    do bilvl=1,mmaxnow   ! level 1 to N. mmax = the total level N.
	    sigmamid=sqrt(hbar*tau/2**bilvl)  
	    !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	    do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		    imid=2**(mmaxnow-bilvl)+i*2**(mmaxnow-bilvl+1)
		    il=imid-2**(mmaxnow-bilvl)
		    ir=imid+2**(mmaxnow-bilvl)
		    !write(6,*) 'bilvl il ir imid sigmamid=',bilvl, il,ir,imid,sigmamid
		    gauss(:,:)=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		    xmid(:,:)=(newristra(il)%x(:,:)+newristra(ir)%x(:,:))/2
		    newristra(imid)%x(:,:)=xmid(:,:)+gauss(:,:)
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
        xnew(:,:)=newristra(j)%x(:,:)       
        call vlsprop(spinorbit,-1,newristra(jp)%x,xnew,cwtl2rtmp2(:,:,jp),cwttmp(:,:))
		call v6proplr(xnew,xnew,dtexpt,cwttmp(:,:),cwtl2rtmp2(:,:,j))  	   
		!if ( sum((xnew-xold)**2).eq.0) then
			!write(6,*) 'bilvl=',bilvl	
			!write(6,*) 'l=',l,j,jp
			!write(6,*) 'sigmamid=',sigmamid
		!   write(6,*) 'xnew-xold=',xnew-xold
		!endif   
	enddo
	jmax=lmax*jd  
	oldlv(bilvl)=pathcoreval !abs(real( sum(conjg(cwtl2rtmp1(:,:,jmax))*cwtendtmp(:,:)) ))	
    
    call vlsprop(spinorbit,-1,newristra(jmax)%x,newristra(jmax+1)%x,cwtl2rtmp2(:,:,jmax),cwttmp(:,:))  
	newlv(bilvl)=abs(real( sum(conjg(cwttmp(:,:))*cwtendtmp(:,:)) )) 
	tmp=log(newlv(bilvl))-log(oldlv(bilvl))
  	!tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
        
    !tmp1test1=abs(real( sum(conjg(cwtl2rtmp1(:,:,jmax))*cwtendtmp(:,:)) ))	
    !write(6,'(''bisect L tmp1-oldp='',t30,i5,3(f15.7,1x))') ,nblocknow, abs(tmp1test1-oldlv(mmaxnow)),tmp1test1,oldlv(mmaxnow)     
    !if (abs(tmp1test1-oldlv(mmaxnow))>=1.0d-3) stop 
   
	rn=randn(1)	   
! can check if l2r and r2l give	the same lv.		    
	    !write(6,*) 'tmp ln rn=',tmp,log(rn(1))
	if (log(rn(1)).lt.tmp) then
	    reject=.false. 
        bisectcountl(bilvl)=bisectcountl(bilvl)
        call addval(10,1.0_r8,1.0_r8)  
	    ibisectl=ibisectl+1 
        ristra(ileftbisect:ileftbisect+nbisectnow)=newristra(0:nbisectnow)
    !ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good.	 
        xl=ristra(0)%x     
	    xr=ristra(nchorizo)%x 
    ! update cwtl2r, cwtr2l.	 
        cwtbegin(:,:)=cwtbegintmp2(:,:)
	    cwtl2r(:,:,ileftbisect:ileftbisect+nbisectnow-1)=cwtl2rtmp2(:,:,0:jmax) ! here jmax should already be 2**mmax-1    
	    if ((ileftbisect+nbisectnow).le.nchorizo-1) then
            do k=ileftbisect+nbisectnow,nchorizo-1
	            x=ristra(k)%x             
                call vlsprop(spinorbit,-1,ristra(k-1)%x,x,cwtl2r(:,:,k-1),cwttmp(:,:))
	            call v6proplr(x,x,-1,cwttmp(:,:),cwtl2r(:,:,k))
            enddo
           
            call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,nchorizo-1),cwttmp(:,:))
	        call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))
            
        else !  ileftbisect+nbisect=nchorizo  

            call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,nchorizo-1),cwttmp(:,:))
	        call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))   
        endif       
    ! update cwtr2l for all the beads
        if ((ileftbisect+nbisectnow-1)>=1) then
	        do k=ileftbisect+nbisectnow-1,1,-1
	            x=ristra(k)%x           
            
                call vlsprop(spinorbit,1,x,ristra(k+1)%x,cwtr2l(:,:,k+1),cwttmp(:,:))
	            call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,k))
            enddo
            
            call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
            call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))  
        else

            call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
            call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))  
        endif
  
        pathcoreval=newlv(mmaxnow)
        
     !!!!!!!!!!!! check   
        !tmp1=abs(real(sum(conjg(cwtbegin(:,:))*cwtr2l(:,:,0))))
        !tmp2=abs(real(sum(conjg(cwtl2r(:,:,nchorizo))*cwtend(:,:))))
        !tmp=(tmp1+tmp2)/2  
        !tmp3=abs(  (pathcoreval-tmp)/( (pathcoreval+tmp)/2  )        )*100
        !write(6,*) 'ileftbisect,mmaxnow:', ileftbisect,mmaxnow
        !write(6,*) 'pathcoreval check bisect LB:', pathcoreval, tmp,tmp1,tmp2,tmp3
        !if (  tmp3 >= 0.1   ) then
        !    write(6,*) '################## Caution ###################'
        !    write(6,*) 'pathcoreval check bisect LB failed % err = ', tmp3
        !    write(6,*) '###############################################'
        !    !stop
        !endif
     !!!!!!!!!!!!!!!!!!         
        
        
        
        
        
        
        
        
        
        
	else
	    reject=.true.  
        call addval(10,0.0_r8,1.0_r8)   
	endif   

    return  
    end subroutine bisectv6lbflex  
    
    subroutine bisectv6rbflex(mmaxnow) ! mmaxnow <= mmax
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
    real(kind=r8) :: tmp,tmp1,tmp2,tmp3,tmp1test1
    real(kind=r8) :: oldlv(0:mmax),newlv(0:mmax)
    complex(kind=r8) :: cwtl2rtmp1(0:nspin-1,nisospin,0:2**mmaxnow),cwtl2rtmp2(0:nspin-1,nisospin,0:2**mmaxnow) &
	                    ,cwtbegintmp(0:nspin-1,nisospin) &
                        ,cwtendtmp1now(0:nspin-1,nisospin),cwtendtmp2now(0:nspin-1,nisospin),cwttmp(0:nspin-1,nisospin)   	                               
    logical :: reject	
   
    ibisecttotr=ibisecttotr+1    
   
    nbisectnow=2**mmaxnow
    ileftbisect=nchorizo-nbisectnow
    irightbisect=nchorizo
  
  
    tau=dt*nbisectnow ! total time for the bisection slice.
    sigmamid=sqrt(2*hbar*tau)  
   
    oldristra(0:nbisectnow)=ristra(ileftbisect:ileftbisect+nbisectnow)
 
    ! deal with the two ends first, then do the classic bisection.   
    newristra(0)=oldristra(0) 
    cwtbegintmp(:,:)=cwtl2r(:,:,ileftbisect)
   
    newristra(nbisectnow)%x=oldristra(0)%x+sigmamid*reshape(gaussian(3*npart),(/3,npart/))
   
    xold(:,:)=oldristra(nbisectnow)%x(:,:)
    xnew(:,:)=newristra(nbisectnow)%x(:,:)   
   
    call corop(xnew,ipro,jpro,cwtendnew)     

    !oldlv(0)=1
    !newlv(0)=1 ! dummy, no use.
   
    cwtl2rtmp1(:,:,0:nbisectnow)=cwtl2r(:,:,ileftbisect:ileftbisect+nbisectnow)
    cwtl2rtmp2(:,:,0)=cwtl2rtmp1(:,:,0)
    cwtendtmp1now=cwtr2l(:,:,nchorizo)
    !call v6propr(oldristra(nbisectnow)%x,-1,cwtend,cwtendtmp1now)  
    call v6propr(newristra(nbisectnow)%x,-1,cwtendnew,cwtendtmp2now)      
    
    do bilvl=1,mmaxnow   ! level 1 to N. mmax = the total level N.
	    sigmamid=sqrt(hbar*tau/2**bilvl)  
	    !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	    do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		    imid=2**(mmaxnow-bilvl)+i*2**(mmaxnow-bilvl+1)
		    il=imid-2**(mmaxnow-bilvl)
		    ir=imid+2**(mmaxnow-bilvl)
		    !write(6,*) 'bilvl il ir imid sigmamid=',bilvl, il,ir,imid,sigmamid
		    gauss(:,:)=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		    xmid=(newristra(il)%x+newristra(ir)%x)/2
		    newristra(imid)%x(:,:)=xmid(:,:)+gauss(:,:)
	    enddo 
    enddo     
    ! judge
    bilvl=mmaxnow  ! besection.  level 1 to N. mmax = the total level N.	   
	bisecttotr(bilvl)=bisecttotr(bilvl)+1   
	lmax=2**bilvl-1 
	jd=2**(mmaxnow-bilvl) ! interval
	dtexpt=mmaxnow-bilvl-1 ! for the R L part, that's why the extra -1. dt=2**dtexp   
       

	do l=1,lmax
		j=l*jd
		jp=(l-1)*jd
		!xold=oldristra(j)%x
        xnew(:,:)=newristra(j)%x(:,:)
        
        call vlsprop(spinorbit,-1,newristra(jp)%x,xnew,cwtl2rtmp2(:,:,jp),cwttmp(:,:))
		call v6proplr(xnew,xnew,dtexpt,cwttmp(:,:),cwtl2rtmp2(:,:,j)) 
        
        
        
		!if ( sum((xnew-xold)**2).eq.0) then
			!write(6,*) 'bilvl=',bilvl	
			!write(6,*) 'l=',l,j,jp
			!write(6,*) 'sigmamid=',sigmamid
		!   write(6,*) 'xnew-xold=',xnew-xold
		!endif   
	enddo
	jmax=lmax*jd    
	oldlv(bilvl)=pathcoreval !abs(real( sum(conjg(cwtl2rtmp1(:,:,jmax))*cwtendtmp1now(:,:)) ))
       
    call vlsprop(spinorbit,-1,newristra(jmax)%x,newristra(jmax+1)%x,cwtl2rtmp2(:,:,jmax),cwttmp(:,:))
	newlv(bilvl)=abs(real( sum(conjg(cwttmp(:,:))*cwtendtmp2now(:,:)) )) 
    
	tmp=log(newlv(bilvl))-log(oldlv(bilvl))
  	!tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
    
    !tmp1test1=abs(real( sum(conjg(cwtl2rtmp1(:,:,jmax))*cwtendtmp1now(:,:)) ))	
    !write(6,'(''bisect R tmp1-oldp='',t30,i5,3(f15.7,1x))') ,nblocknow, abs(tmp1test1-oldlv(mmaxnow)),tmp1test1,oldlv(mmaxnow)     
    !if (abs(tmp1test1-oldlv(mmaxnow))>=1.0d-7) stop 
  
	rn=randn(1)	   
! can check if l2r and r2l give	the same lv.		    
	    !write(6,*) 'tmp ln rn=',tmp,log(rn(1))
	if (log(rn(1)).lt.tmp) then
	    reject=.false. 
        bisectcountr(bilvl)=bisectcountr(bilvl)+1
        call addval(11,1.0_r8,1.0_r8)  
	    ibisectr=ibisectr+1 
        ristra(ileftbisect:ileftbisect+nbisectnow)=newristra(0:nbisectnow)      	 
    ! update cwtl2r, cwtr2l.	 
    ! update cwtl2r for all beads. this can be optimized in the future.
             
        cwtl2r(:,:,ileftbisect+1:ileftbisect+nbisectnow-1)=cwtl2rtmp2(:,:,1:jmax) ! here jmax should already be 2**mmax-1 
       
        call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,ristra(nchorizo)%x,cwtl2r(:,:,nchorizo-1),cwttmp(:,:))
        call v6propl(ristra(nchorizo)%x,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo)) ! ileftbisect+nbisectnow=nchorizo  
        
        ! update cwtr2l for all the beads
        cwtend(:,:)=cwtendnew(:,:)          
        cwtr2l(:,:,nchorizo)=cwtendtmp2now
	    do k=nchorizo-1,1,-1
            x=ristra(k)%x           
            call vlsprop(spinorbit,1,x,ristra(k+1)%x,cwtr2l(:,:,k+1),cwttmp(:,:))
	        call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,k))
	    enddo
        xl=ristra(0)%x
         
        call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
        call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))  	
        
        pathcoreval=newlv(mmaxnow)
        
        
     !!!!!!!!!!!! check   
        !tmp1=abs(real(sum(conjg(cwtbegin(:,:))*cwtr2l(:,:,0))))
        !tmp2=abs(real(sum(conjg(cwtl2r(:,:,nchorizo))*cwtend(:,:))))
        !tmp=(tmp1+tmp2)/2  
        !tmp3=abs(  (pathcoreval-tmp)/( (pathcoreval+tmp)/2  )        )*100
        !write(6,*) 'ileftbisect,mmaxnow:', ileftbisect,mmaxnow
        !write(6,*) 'pathcoreval check bisect RB:', pathcoreval, tmp,tmp1,tmp2,tmp3
        !if (  tmp3 >= 0.1   ) then
        !    write(6,*) '################## Caution ###################'
        !    write(6,*) 'pathcoreval check bisect RB failed % err = ', tmp3
        !    write(6,*) '###############################################'
        !    !stop
        !endif
     !!!!!!!!!!!!!!!!!!          
        
        
        
        
        
        
        
        
        
	else
	    reject=.true.  
        call addval(11,0.0_r8,1.0_r8)      
	endif    
    return  
    end subroutine bisectv6rbflex     

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
    real(kind=r8) :: rt2,tmp,tmp1,tmp2,tmp3,logratio,gc
    complex(kind=r8) :: cwtold(0:nspin-1,nisospin),cwtnew(0:nspin-1,nisospin),cwttmp(0:nspin-1,nisospin)
    complex(kind=r8) :: cwtl2rtmp(0:nspin-1,nisospin),cwtr2ltmp(0:nspin-1,nisospin)

    if ((ileft.lt.0).or.(iright.gt.nchorizo).or.(iright.lt.0).or.(ileft.gt.nchorizo).or.(ileft.gt.iright)) then
	    write (6,*) 'stdmove1fastv6 range error', ileft,iright
	    stop
    endif
       
    
    gc=-0.5_r8/(dt*(2*hbar))
    
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
	 
        gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
     
        !gauss=mov1step*reshape((randn(3*npart)-0.5),(/3,npart/)) ! this is not gauss move. for vmc(dt=0) only
     
	    xnew(:,:)=xold(:,:)+gauss(:,:)  
	 
        if (i.eq.0) then
		 
	        !call psitcwt(xold,cwtold)
	  
	        cwtold(:,:)=cwtbegin(:,:)
	  
	        tmp1=log(pathcoreval)!log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i)))))
	        !call psitcwt(xnew,cwtnew)  
	        call corop(xnew,iplo,jplo,cwtnew)
            
            call vlsprop(spinorbit,1,xnew,ristra(i+1)%x,cwtr2l(:,:,i+1),cwttmp(:,:))
	        call v6propl(xnew,-1,cwttmp(:,:),cwtr2ltmp(:,:))

	        tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
	        xd2new(:,:)=(ristra(i+1)%x(:,:)-xnew(:,:))**2
	        xd2old(:,:)=(ristra(i+1)%x(:,:)-xold(:,:))**2

            rt2=gc*(sum(xd2new)-sum(xd2old))	  
            !rt2=0. ! vmc(dt=0) only.
	  
	        logratio=tmp2-tmp1+rt2
	        rn=randn(1)
  
	        if (log(rn(1)).lt.logratio) then ! dt ne 0
                icnow=icnow+1 
	            icmov1=icmov1+1	
	            ristra(i)%x(:,:)=xnew(:,:)
	            cwtr2l(:,:,i)=cwtr2ltmp(:,:)   
            ! update cwtbegin
                cwtbegin(:,:)=cwtnew(:,:)	   
            ! update cwtl2r	
               
	            call v6propr(xnew,-1,cwtnew(:,:),cwtl2r(:,:,i)) ! just update the i.
            ! prepare for the next cwtl2r
      
                xoldnxt=ristra(i+1)%x
                                
                call vlsprop(spinorbit,-1,ristra(i)%x,xoldnxt,cwtl2r(:,:,i),cwttmp(:,:))
                call v6proplr(xoldnxt,xoldnxt,-1,cwttmp(:,:),cwtl2r(:,:,i+1)) 
            
                pathcoreval=exp(tmp2)
	        endif 
	  
	    else
		 
	        if (i.eq.nchorizo) then   
	            !call psitcwt(xold,cwtold)
	            cwtold(:,:)=cwtend(:,:)
	            tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtl2r(:,:,i))*cwtold(:,:)))))
	            !call psitcwt(xnew,cwtnew)
	            call corop(xnew,ipro,jpro,cwtnew)
                
                call vlsprop(spinorbit,-1,ristra(i-1)%x,xnew,cwtl2r(:,:,i-1),cwttmp(:,:))
	            call v6propl(xnew,-1,cwttmp(:,:),cwtl2rtmp(:,:))
                              
	            tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtnew(:,:)))))
	            xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2
	            xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2
	  
	            rt2=gc*(sum(xd2new)-sum(xd2old))	  
                !rt2=0. ! vmc(dt=0) only.
	  
	            logratio=tmp2-tmp1+rt2
	            rn=randn(1)
	  
	            if (log(rn(1)).lt.logratio) then
                    icnow=icnow+1
	                icmov1=icmov1+1	
	                ristra(i)%x(:,:)=xnew(:,:)
	                cwtl2r(:,:,i)=cwtl2rtmp(:,:)
                ! update cwtend
                    cwtend(:,:)=cwtnew(:,:)
                    pathcoreval=exp(tmp2)
	            endif 		     
	   
	        else
	   
                cwtold(:,:)=cwtl2r(:,:,i)
                                
                call vlsprop(spinorbit,-1,ristra(i-1)%x,xnew,cwtl2r(:,:,i-1),cwttmp(:,:))                
	            call v6proplr(xnew,xnew,-1,cwttmp(:,:),cwtl2rtmp(:,:))                             
                call vlsprop(spinorbit,-1,xnew,ristra(i+1)%x,cwtl2rtmp(:,:),cwttmp(:,:))
                
	            tmp2=log(abs(real(sum(conjg(cwttmp(:,:))*cwtr2l(:,:,i+1)))))        
	            tmp1=log(pathcoreval)!log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i+1)))))
	    
	            xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2+(ristra(i+1)%x(:,:)-xnew(:,:))**2
	            xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2+(ristra(i+1)%x(:,:)-xold(:,:))**2
	   
	            rt2=gc*(sum(xd2new)-sum(xd2old))	
                !rt2=0.  ! vmc(dt=0) only.
	   
	            logratio=tmp2-tmp1+rt2	 
	            rn=randn(1)

	            if (log(rn(1)).lt.logratio) then
                        !write(6,*) 'chorizo i updated!'
                    icnow=icnow+1           
	                icmov1=icmov1+1	
	                ristra(i)%x(:,:)=xnew(:,:)
	                cwtl2r(:,:,i)=cwtl2rtmp(:,:)
                    pathcoreval=exp(tmp2)
                endif
       
                if ( icnow /= 0 ) then
            ! update l2r
                    if (i.le.(nchorizo-2)) then
		                x=ristra(i+1)%x                      
                        call vlsprop(spinorbit,-1,ristra(i)%x,x,cwtl2r(:,:,i),cwttmp(:,:))
		                call v6proplr(x,x,-1,cwttmp(:,:),cwtl2r(:,:,i+1))     
                    else ! i=nchorizo-1
                        xr=ristra(i+1)%x                    
                        call vlsprop(spinorbit,-1,ristra(i)%x,xr,cwtl2r(:,:,i),cwttmp(:,:))
                        call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,i+1))	    
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
                call vlsprop(spinorbit,-1,ristra(k)%x,x,cwtl2r(:,:,k),cwttmp(:,:))
                call v6proplr(x,x,-1,cwttmp(:,:),cwtl2r(:,:,k+1))   
            enddo
                        
            call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,nchorizo-1),cwttmp(:,:))
            call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))	      
        else  
            if (iright == nchorizo-1) then                         
                call vlsprop(spinorbit,-1,ristra(iright)%x,xr,cwtl2r(:,:,iright),cwttmp(:,:))
                call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))	      
            endif 
        endif     
        ! update cwtr2l from iright th bead until the 0th.
        if (iright == nchorizo) then    
            !x=ristra(iright)%x   
	        call v6propr(xr,-1,cwtend(:,:),cwtr2l(:,:,iright))
	        do k=iright-1,1,-1
	           x=ristra(k)%x
               call vlsprop(spinorbit,1,x,ristra(k+1)%x,cwtr2l(:,:,k+1),cwttmp(:,:))
	           call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,k))
            enddo
	        !xl=ristra(0)%x
            call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
	        call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	    
        else      
            if (iright >= 2) then  
                x=ristra(iright)%x
                call vlsprop(spinorbit,1,x,ristra(iright+1)%x,cwtr2l(:,:,iright+1),cwttmp(:,:))
                call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,iright))
                do k=iright-1,1,-1
	               x=ristra(k)%x
                   call vlsprop(spinorbit,1,x,ristra(k+1)%x,cwtr2l(:,:,k+1),cwttmp(:,:))
	               call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,k))
                enddo
	            !xl=ristra(0)%x
                call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
	            call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	   
            else ! iright =1 or 0
                if (iright == 1) then   
                    x=ristra(iright)%x
                    call vlsprop(spinorbit,1,x,ristra(iright+1)%x,cwtr2l(:,:,iright+1),cwttmp(:,:))
                    call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,iright))   
	                !xl=ristra(0)%x
                    call vlsprop(spinorbit,1,xl,ristra(iright)%x,cwtr2l(:,:,iright),cwttmp(:,:))
	                call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	
                else ! iright =0
                    call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
                    call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	   
                endif
       
            endif  
        endif
    endif
   
    
    return
    end subroutine stdmove1fastv6     
    
    
    subroutine stdmove1fastv6flex(ileft,iright,arrow) 
    ! ileft does not have to be > iright
    use wavefunction
    use estimator
    use v6stepcalc
    use random
    use brut
    use mympi
    integer(kind=i4) :: i,k,ileft,iright,icnow
    real(kind=r8) :: rn(1),gauss(3,npart)
    real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart),xoldnxt(3,npart),xd2old(3,npart),xd2new(3,npart) &
	                ,xr(3,npart),xl(3,npart)
    real(kind=r8) :: rt2,tmp,tmp1,tmp2,tmp3,gc,logratio,tmp1test1,tmp1test2
    complex(kind=r8) :: cwtold(0:nspin-1,nisospin),cwtnew(0:nspin-1,nisospin)
    complex(kind=r8) :: cwtl2rtmp(0:nspin-1,nisospin),cwtr2ltmp(0:nspin-1,nisospin),cwttmp(0:nspin-1,nisospin)
    integer(kind=i4) :: arrow

    if ((ileft.lt.0).or.(iright.gt.nchorizo).or.(iright.lt.0).or.(ileft.gt.nchorizo).or.(ileft.gt.iright)) then
	    write (6,*) 'stdmove1fastv6flex range error', ileft,iright
	    stop
    endif
    if (ileft.eq.iright) then
        call stdmove1fastv6(ileft,iright) 
        return
    endif     
    
    gc=-0.5_r8/(dt*(2*hbar))
  
    select case (arrow)
        
    case (-1) ! inverse order from right to left
        
        icnow=0
        icstdmove1rvstot=icstdmove1rvstot+1
        icmov1rvstot=abs(iright-ileft+1)*icstdmove1rvstot
        
        do i=iright,ileft,-1  
	        xold=ristra(i)%x 
            gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
            !gauss=mov1step*reshape((randn(3*npart)-0.5),(/3,npart/)) ! this is not gauss move. for vmc(dt=0) only
	        xnew(:,:)=xold(:,:)+gauss(:,:)            
            
            if (i.eq.0) then
		 
	        !call psitcwt(xold,cwtold)
	  
	            cwtold(:,:)=cwtbegin(:,:)
	            tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i)))))      
                
                !tmp1test1=log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i)))))                
                !write(6,'(''tmp1,tmp1test1='',t30,i5,3(f15.7,1x))') ,i, abs(tmp1-tmp1test1),tmp1,tmp1test1  
                !if (abs(tmp1test1-tmp1)>=1.0d-7) stop           
                
	            !call psitcwt(xnew,cwtnew)  
	            call corop(xnew,iplo,jplo,cwtnew)
                
                call vlsprop(spinorbit,1,xnew,ristra(i+1)%x,cwtr2l(:,:,i+1),cwttmp(:,:))             
	            call v6propl(xnew,-1,cwttmp(:,:),cwtr2ltmp(:,:))
                
                
                
	            tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
	            xd2new(:,:)=(ristra(i+1)%x(:,:)-xnew(:,:))**2
	            xd2old(:,:)=(ristra(i+1)%x(:,:)-xold(:,:))**2
	  
	            rt2=gc*(sum(xd2new)-sum(xd2old))	  
                !rt2=0. ! vmc(dt=0) only.
	  
	            logratio=tmp2-tmp1+rt2
	            rn=randn(1)
  
	            if (log(rn(1)).lt.logratio) then ! dt ne 0
                    icnow=icnow+1 
	                icmov1rvs=icmov1rvs+1	
	                ristra(i)%x(:,:)=xnew(:,:)
	                cwtr2l(:,:,i)=cwtr2ltmp(:,:)   
                ! update cwtbegin
                    cwtbegin(:,:)=cwtnew(:,:)   
                    pathcoreval=exp(tmp2)
                endif 
                
	        else
		 
	            if (i.eq.nchorizo) then   
	                !call psitcwt(xold,cwtold)
	                cwtold(:,:)=cwtend(:,:)
	                tmp1=log(pathcoreval)  !log(abs(real(sum(conjg(cwtl2r(:,:,i))*cwtold(:,:))))) ! can put those tmp1 thing in the beginning, no need to repreatedly calculate                
                    
                !tmp1test1=log(abs(real(sum(conjg(cwtl2r(:,:,i))*cwtold(:,:)))))           
                !write(6,'(''tmp1,tmp1test1='',t30,i5,3(f15.7,1x))') ,i, abs(tmp1-tmp1test1),tmp1,tmp1test1                     
                ! if (abs(tmp1test1-tmp1)>=1.0d-7) stop    
                    
                    !call psitcwt(xnew,cwtnew)
	                call corop(xnew,ipro,jpro,cwtnew)     
                    call v6propr(xnew,-1,cwtnew(:,:),cwtr2ltmp(:,:))      
                    call vlsprop(spinorbit,1,ristra(i-1)%x,xnew,cwtl2r(:,:,i-1),cwttmp(:,:))
                    
                    tmp2=log(abs(real(sum(conjg(cwttmp(:,:))*cwtr2ltmp(:,:)))))
                
	                xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2
	                xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2
	  
	                rt2=gc*(sum(xd2new)-sum(xd2old))	  
                    !rt2=0. ! vmc(dt=0) only.
	  
	                logratio=tmp2-tmp1+rt2
	                rn=randn(1)
	  
	                if (log(rn(1)).lt.logratio) then
                        icnow=icnow+1
	                    icmov1rvs=icmov1rvs+1	
	                    ristra(i)%x(:,:)=xnew(:,:)
	                    cwtr2l(:,:,i)=cwtr2ltmp(:,:)
                    ! update cwtend and prepare for the left bead.
                        cwtend(:,:)=cwtnew(:,:)
                        x=ristra(i-1)%x                      
                        call vlsprop(spinorbit,1,x,ristra(i)%x,cwtr2l(:,:,i),cwttmp(:,:))                
		                call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,i-1))
                        pathcoreval=exp(tmp2)
	                endif 		     
	   
                else
                  
                    call vlsprop(spinorbit,1,xnew,ristra(i+1)%x,cwtr2l(:,:,i+1),cwttmp(:,:))      
	                call v6proplr(xnew,xnew,-1,cwttmp(:,:),cwtr2ltmp(:,:))
                    call vlsprop(spinorbit,1,ristra(i-1)%x,xnew,cwtr2ltmp(:,:),cwttmp(:,:))
                    
	                tmp2=log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwttmp(:,:)))))
                    tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2l(:,:,i)))))

                !tmp1test1=log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2l(:,:,i)))))
                !write(6,'(''tmp1,tmp1test1='',t30,i5,3(f15.7,1x))') ,i, abs(tmp1-tmp1test1),tmp1,tmp1test1                        
                !if (abs(tmp1test1-tmp1)>=1.0d-7) stop     
                    
                    
                    
	                xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2+(ristra(i+1)%x(:,:)-xnew(:,:))**2
	                xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2+(ristra(i+1)%x(:,:)-xold(:,:))**2
	   
	                rt2=gc*(sum(xd2new)-sum(xd2old))	
                    !rt2=0.  ! vmc(dt=0) only.
	   
	                logratio=tmp2-tmp1+rt2	 
	                rn=randn(1)

	                if (log(rn(1)).lt.logratio) then
                        icnow=icnow+1           
	                    icmov1rvs=icmov1rvs+1	
	                    ristra(i)%x(:,:)=xnew(:,:)
	                    cwtr2l(:,:,i)=cwtr2ltmp(:,:)	
                        pathcoreval=exp(tmp2)
                    endif
       
                    if ( icnow /= 0 ) then
                ! update r2l for the next r2l bead, the left bead
                        if (i >= 2) then
		                    x=ristra(i-1)%x
                            call vlsprop(spinorbit,1,x,ristra(i)%x,cwtr2l(:,:,i),cwttmp(:,:))
		                    call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,i-1))
                        else ! i=1
                            xl=ristra(i-1)%x
                            call vlsprop(spinorbit,1,xl,ristra(i)%x,cwtr2l(:,:,i),cwttmp(:,:))
                            call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,i-1))	    
                        endif               
                    endif 
    
	            endif 
	        endif                 
            
        enddo
        
        if (icnow /= 0) then 
            xl=ristra(0)%x
            xr=ristra(nchorizo)%x
            ! update r2l
            if ( ileft >= 3 ) then
                do k=ileft,3,-1 
                   x=ristra(k-2)%x
                   call vlsprop(spinorbit,1,x,ristra(k-1)%x,cwtr2l(:,:,k-1),cwttmp(:,:))
	               call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,k-2))  
                enddo
                call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
                call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	 
            else
                if (ileft == 2) then
                    call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
                    call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))
                endif     
            endif
            ! update l2r
            if (ileft >= 1) then
                if (ileft <= nchorizo-1 ) then    
                    do k=ileft,nchorizo-1
                       x=ristra(k)%x
                       call vlsprop(spinorbit,-1,ristra(k-1)%x,x,cwtl2r(:,:,k-1),cwttmp(:,:))
	                   call v6proplr(x,x,-1,cwttmp(:,:),cwtl2r(:,:,k))
                    enddo
                    call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,nchorizo-1),cwttmp(:,:))
                    call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))	   
                else ! ileft = nchorizo
                    call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,nchorizo-1),cwttmp(:,:))
                    call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))	   
                endif 
            else ! ileft == 0
                call v6propr(xl,-1,cwtbegin(:,:),cwtl2r(:,:,0))	   
                    do k=1,nchorizo-1
                       x=ristra(k)%x
                       call vlsprop(spinorbit,-1,ristra(k-1)%x,x,cwtl2r(:,:,k-1),cwttmp(:,:))
	                   call v6proplr(x,x,-1,cwttmp(:,:),cwtl2r(:,:,k))
                    enddo
                    call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,nchorizo-1),cwttmp(:,:))
                    call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))	     
            endif
        endif
           
  
    case (1) ! normal, from left to right, should be just the same as stdmove1fastv6(ileft,iright)  

        icnow=0
        icstdmove1tot=icstdmove1tot+1
        icmov1tot=abs(iright-ileft+1)*icstdmove1tot
        do i=ileft,iright
      
            !write(6,*) 'i=',i
	        xold=ristra(i)%x 
        ! ran move	 
	        !rn3=reshape(randn(3*npart),shape(rn3)) 
	        !xnew=xold+mov1step*(rn3-0.5_r8) 
        ! gaussian move	 
	 
            gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
     
            !gauss=mov1step*reshape((randn(3*npart)-0.5),(/3,npart/)) ! this is not gauss move. for vmc(dt=0) only
     
	        xnew(:,:)=xold(:,:)+gauss(:,:)  
	 
            if (i.eq.0) then
		 
	            !call psitcwt(xold,cwtold)
	  
	            cwtold(:,:)=cwtbegin(:,:)
	  
	            tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i)))))
                
                !tmp1test1=log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i)))))              
                !write(6,'(''-> tmp1,tmp1test1='',t30,i5,3(f15.7,1x))') ,i, abs(tmp1-tmp1test1),tmp1,tmp1test1                 
                
                
	            !call psitcwt(xnew,cwtnew)  
	            call corop(xnew,iplo,jplo,cwtnew)
                
                call vlsprop(spinorbit,1,xnew,ristra(i+1)%x,cwtr2l(:,:,i+1),cwttmp(:,:))
                
	            call v6propl(xnew,-1,cwttmp(:,:),cwtr2ltmp(:,:))
                
	            tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
	            xd2new(:,:)=(ristra(i+1)%x(:,:)-xnew(:,:))**2
	            xd2old(:,:)=(ristra(i+1)%x(:,:)-xold(:,:))**2
	  
	            rt2=gc*(sum(xd2new)-sum(xd2old))	  
                !rt2=0. ! vmc(dt=0) only.
	  
	            logratio=tmp2-tmp1+rt2
	            rn=randn(1)
  
	            if (log(rn(1)).lt.logratio) then ! dt ne 0
                    icnow=icnow+1 
	                icmov1=icmov1+1	
	                ristra(i)%x(:,:)=xnew(:,:)
	                cwtr2l(:,:,i)=cwtr2ltmp(:,:)   
                ! update cwtbegin
                    cwtbegin(:,:)=cwtnew(:,:)  
                ! update cwtl2r	
                    
	                call v6propr(xnew,-1,cwtnew(:,:),cwtl2r(:,:,i)) ! just update the i.
                ! prepare for the next cwtl2r
                    xoldnxt=ristra(i+1)%x
                    call vlsprop(spinorbit,-1,ristra(i)%x,xoldnxt,cwtl2r(:,:,i),cwttmp(:,:))
                    call v6proplr(xoldnxt,xoldnxt,-1,cwttmp(:,:),cwtl2r(:,:,i+1)) 
                    
                    pathcoreval=exp(tmp2)
	            endif 
	  
	        else
		 
	        if (i.eq.nchorizo) then   
	            !call psitcwt(xold,cwtold)
	            cwtold(:,:)=cwtend(:,:)
	            tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtl2r(:,:,i))*cwtold(:,:)))))
                
                !tmp1test1=log(abs(real(sum(conjg(cwtl2r(:,:,i))*cwtold(:,:)))))            
                !write(6,'(''-> tmp1,tmp1test1='',t30,i5,3(f15.7,1x))') ,i, abs(tmp1-tmp1test1),tmp1,tmp1test1 
                
	            !call psitcwt(xnew,cwtnew)
	            call corop(xnew,ipro,jpro,cwtnew)
                call vlsprop(spinorbit,-1,ristra(i-1)%x,xnew,cwtl2r(:,:,i-1),cwttmp(:,:))
	            call v6propl(xnew,-1,cwttmp(:,:),cwtl2rtmp(:,:))
	            tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtnew(:,:)))))
	            xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2
	            xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2
	  
	            rt2=gc*(sum(xd2new)-sum(xd2old))	  
                !rt2=0. ! vmc(dt=0) only.
	  
	            logratio=tmp2-tmp1+rt2
	            rn=randn(1)
	  
	            if (log(rn(1)).lt.logratio) then
                    icnow=icnow+1
	                icmov1=icmov1+1	
	                ristra(i)%x(:,:)=xnew(:,:)
	                cwtl2r(:,:,i)=cwtl2rtmp(:,:)
                ! update cwtend
                    cwtend(:,:)=cwtnew(:,:)
                    
                    pathcoreval=exp(tmp2)
	            endif 		     
	   
	        else
	   
                cwtold(:,:)=cwtl2r(:,:,i)
                
                call vlsprop(spinorbit,-1,ristra(i-1)%x,xnew,cwtl2r(:,:,i-1),cwttmp(:,:))          
	            call v6proplr(xnew,xnew,-1,cwttmp(:,:),cwtl2rtmp(:,:))
	            call vlsprop(spinorbit,-1,xnew,ristra(i+1)%x,cwtl2rtmp(:,:),cwttmp(:,:))
	            tmp2=log(abs(real(sum(conjg(cwttmp(:,:))*cwtr2l(:,:,i+1)))))
	            tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i+1)))))
                
                !tmp1test1=log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i+1)))))
                !write(6,'(''-> tmp1,tmp1test1='',t30,i5,3(f15.7,1x))') ,i, abs(tmp1-tmp1test1),tmp1,tmp1test1                
	    
	            xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2+(ristra(i+1)%x(:,:)-xnew(:,:))**2
	            xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2+(ristra(i+1)%x(:,:)-xold(:,:))**2
	   
	            rt2=gc*(sum(xd2new)-sum(xd2old))	
                !rt2=0.  ! vmc(dt=0) only.
	   
	            logratio=tmp2-tmp1+rt2	 
	            rn=randn(1)

	            if (log(rn(1)).lt.logratio) then
                        !write(6,*) 'chorizo i updated!'
                    icnow=icnow+1           
	                icmov1=icmov1+1	
	                ristra(i)%x(:,:)=xnew(:,:)
	                cwtl2r(:,:,i)=cwtl2rtmp(:,:)	
                    
                    pathcoreval=exp(tmp2)
                endif
       
                if ( icnow /= 0 ) then
            ! update l2r
                    if (i.le.(nchorizo-2)) then
		                x=ristra(i+1)%x
                        call vlsprop(spinorbit,-1,ristra(i)%x,x,cwtl2r(:,:,i),cwttmp(:,:))
		                call v6proplr(x,x,-1,cwttmp(:,:),cwtl2r(:,:,i+1))
                    else ! i=nchorizo-1
                        xr=ristra(i+1)%x
                        call vlsprop(spinorbit,-1,ristra(i)%x,xr,cwtl2r(:,:,i),cwttmp(:,:))
                        call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,i+1))	    
                    endif               
                endif 
    
	        endif 
	        endif 	 
        enddo 

        
! the below parts are slightly differnt from that in stdmove1fastv6(ileft,iright)          
        if (icnow /= 0) then 
            xl=ristra(0)%x
            xr=ristra(nchorizo)%x
            ! update cwtl2r    
            if (iright <= nchorizo-3) then 
                do k=iright,nchorizo-3
                    x=ristra(k+2)%x
                    call vlsprop(spinorbit,-1,ristra(k+1)%x,x,cwtl2r(:,:,k+1),cwttmp(:,:))
                    call v6proplr(x,x,-1,cwttmp(:,:),cwtl2r(:,:,k+2))   
                enddo
                call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,nchorizo-1),cwttmp(:,:))
                call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))	      
            else    
                if (iright == nchorizo-2) then 
                    call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,iright+1),cwttmp(:,:))
                    call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))	          
                endif  
            endif    
            ! update cwtr2l from iright th bead until the 0th.
            if (iright == nchorizo) then    
                !x=ristra(iright)%x   
	            call v6propr(xr,-1,cwtend(:,:),cwtr2l(:,:,iright))
	            do k=iright-1,1,-1
	                x=ristra(k)%x
                    call vlsprop(spinorbit,1,x,ristra(k+1)%x,cwtr2l(:,:,k+1),cwttmp(:,:))
	                call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,k))
                enddo
	            !xl=ristra(0)%x
                call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
	            call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	    
            else      
                if (iright >= 1) then  
                    do k=iright,1,-1
	                   x=ristra(k)%x
                       call vlsprop(spinorbit,1,x,ristra(k+1)%x,cwtr2l(:,:,k+1),cwttmp(:,:))
	                   call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,k))
                    enddo
	                !xl=ristra(0)%x
                    call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
	                call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	   
                else ! iright = 0
                    call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
                    call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	   
                endif  
            endif
        endif   
    case default 
        write (6,*) 'arrow for stdmove1fastv6flex has not been set, abort!', arrow
	    call abort
    end select
    
    
    
    
     !!!!!!!!!!!! check   
        !tmp1=abs(real(sum(conjg(cwtbegin(:,:))*cwtr2l(:,:,0))))
        !tmp2=abs(real(sum(conjg(cwtl2r(:,:,nchorizo))*cwtend(:,:))))
        !tmp=(tmp1+tmp2)/2  
        !tmp3=abs(  (pathcoreval-tmp)/( (pathcoreval+tmp)/2  )        )*100
        !write(6,*) 'pathcoreval check stdmove1fastv6flex  :', pathcoreval, tmp,tmp1,tmp2,tmp3
        !if (  tmp3 >= 0.1   ) then
        !    write(6,*) '################## Caution ###################'
        !    write(6,*) 'pathcoreval check stdmove1fastv6flex   failed % err = ', tmp3
        !    write(6,*) '###############################################'
        !    !stop
        !endif
     !!!!!!!!!!!!!!!!!!     
    
    
    
    
    
    
    return
    end subroutine stdmove1fastv6flex              

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
    real(kind=r8) :: rt2(0:nchorizo),tmp,tmp1,tmp2,tmp3,logratio,tmp1test1,gc
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin)   

    if ((ileft.lt.0).or.(iright.gt.nchorizo).or.(iright.lt.0).or.(ileft.gt.nchorizo).or.(ileft.gt.iright)) then
	    write (6,*) 'stdmove1fastv6 range error', ileft,iright
	    stop
    else
    if (ileft.eq.iright) then
        call stdmove1fastv6(ileft,iright) 
        return
    endif          
    endif

    gc=-0.5_r8/(dt*(2*hbar))
    
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
	    ristranew(i)%x(:,:)=ristraold(i)%x(:,:)+gauss(:,:)    
    enddo
   
    tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtl2r(:,:,0))*cwtr2l(:,:,1)))))
    
    !tmp1test1=log(abs(real(sum(conjg(cwtl2r(:,:,0))*cwtr2l(:,:,1)))))
    !write(6,'(''2v6 tmp1,tmp1test1:'',t30,i5,3(f15.7,1x))') ,nblocknow, abs(tmp1-tmp1test1),tmp1,tmp1test1         
    !if (abs(tmp1test1-tmp1)>=1.0d-7) stop 
    
    
    tmp2=tmp1   
    rt2=0.
    cwtl2rnew(:,:,:)=cwtl2r(:,:,:)
    cwtr2lnew(:,:,:)=cwtr2l(:,:,:)
   
    do i=ileft,iright  
   
        if (i.eq.0) then
		 
	        !call psitcwt(xold,cwtold)  
	        call corop(ristranew(i)%x,iplo,jplo,cwtbeginnew)
	        call v6propr(ristranew(i)%x,-1,cwtbeginnew(:,:),cwtl2rnew(:,:,i)) ! just update the i.
      
	        xd2new(:,:)=(ristranew(i)%x(:,:)-ristranew(i+1)%x(:,:))**2
	        xd2old(:,:)=(ristraold(i)%x(:,:)-ristraold(i+1)%x(:,:))**2
	        rt2(i)=gc*(sum(xd2new)-sum(xd2old))	  
            !rt2=0. ! vmc(dt=0) only.
	     
	    else
		 
	        if (i.eq.nchorizo) then   
	            !call psitcwt(xold,cwtold)

	            !call psitcwt(xnew,cwtnew)
	            call corop(ristranew(i)%x,ipro,jpro,cwtendnew)             
                call vlsprop(spinorbit,-1,ristranew(i-1)%x,ristranew(i)%x,cwtl2rnew(:,:,i-1),cwttmp(:,:))
	            call v6propl(ristranew(i)%x,-1,cwttmp(:,:),cwtl2rnew(:,:,i))
                call v6propr(ristranew(i)%x,-1,cwtendnew,cwtr2lnew(:,:,i))
         
            else
	   
                call vlsprop(spinorbit,-1,ristranew(i-1)%x,ristranew(i)%x,cwtl2rnew(:,:,i-1),cwttmp(:,:))
	            call v6proplr(ristranew(i)%x,ristranew(i)%x,-1,cwttmp(:,:),cwtl2rnew(:,:,i))
	    
                xd2new(:,:)=(ristranew(i)%x(:,:)-ristranew(i+1)%x(:,:))**2
	            xd2old(:,:)=(ristraold(i)%x(:,:)-ristraold(i+1)%x(:,:))**2
	            rt2(i)=gc*(sum(xd2new)-sum(xd2old))	
                !rt2=0.  ! vmc(dt=0) only.

	        endif 
	    endif 	 
    enddo 

    if (iright /= nchorizo) then
        call vlsprop(spinorbit,-1,ristranew(iright)%x,ristranew(iright+1)%x,cwtl2rnew(:,:,iright),cwttmp(:,:))
        tmp2=log(abs(real(sum(conjg(cwttmp(:,:))*cwtr2lnew(:,:,iright+1)))))
    else
        tmp2=log(abs(real(sum(conjg(cwtl2rnew(:,:,iright))*cwtendnew(:,:)))))   
    endif
   
    if (ileft/=0) then
        xd2new=(  ristranew(ileft)%x - ristraold(ileft-1)%x  )**2
	    xd2old=(  ristraold(ileft)%x - ristraold(ileft-1)%x  )**2       
        rt2(ileft-1)=gc*(sum(xd2new)-sum(xd2old))	
    endif
   
    logratio=tmp2-tmp1+sum(rt2)
    rn=randn(1)
   
    if (log(rn(1)).lt.logratio) then     
        icstdmove2=icstdmove2+1   
        do i=ileft,iright
            ristra(i)=ristranew(i)
            cwtl2r(:,:,i)=cwtl2rnew(:,:,i)
        enddo  
        if (ileft.eq.0) cwtbegin(:,:)=cwtbeginnew(:,:)
        if (iright.eq.nchorizo) cwtend(:,:)=cwtendnew(:,:)      
        pathcoreval=exp(tmp2)      
    ! update the chain.
        xl=ristra(0)%x
        xr=ristra(nchorizo)%x
    ! update cwtl2r 
        if (iright <= nchorizo-2) then 
            do k=iright,nchorizo-2
                x=ristra(k+1)%x
                call vlsprop(spinorbit,-1,ristra(k)%x,x,cwtl2r(:,:,k),cwttmp(:,:))
                call v6proplr(x,x,-1,cwttmp(:,:),cwtl2r(:,:,k+1))   
            enddo
            call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,nchorizo-1),cwttmp(:,:))
            call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))	      
        else  
            if (iright == nchorizo-1) then
                call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,iright),cwttmp(:,:))
                call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))	          
            endif 
        endif     
        ! update cwtr2l from iright th bead until the 0th.
        if (iright == nchorizo) then    
            !x=ristra(iright)%x   
	        call v6propr(xr,-1,cwtend(:,:),cwtr2l(:,:,iright))
	        do k=iright-1,1,-1
	            x=ristra(k)%x
                call vlsprop(spinorbit,1,x,ristra(k+1)%x,cwtr2l(:,:,k+1),cwttmp(:,:))
	            call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,k))
            enddo
	        !xl=ristra(0)%x
            call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
	        call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	    
        else      
            if (iright >= 2) then  
                x=ristra(iright)%x
                call vlsprop(spinorbit,1,x,ristra(iright+1)%x,cwtr2l(:,:,iright+1),cwttmp(:,:))
                call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,iright))
                do k=iright-1,1,-1
	                x=ristra(k)%x
                    call vlsprop(spinorbit,1,x,ristra(k+1)%x,cwtr2l(:,:,k+1),cwttmp(:,:))
	                call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,k))
                enddo
	            !xl=ristra(0)%x
                call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
	            call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	   
            else ! iright =1 or 0        
                if (iright == 1) then   
                    x=ristra(iright)%x
                    call vlsprop(spinorbit,1,x,ristra(iright+1)%x,cwtr2l(:,:,iright+1),cwttmp(:,:))
                    call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,iright))   
	                !xl=ristra(0)%x
                    call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,iright),cwttmp(:,:))
	                call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	
                else ! iright =0
                    call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
                    call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	   
                endif 
            endif  
        endif 
    endif
   
    
    
     !!!!!!!!!!!! check   
        !tmp1=abs(real(sum(conjg(cwtbegin(:,:))*cwtr2l(:,:,0))))
        !tmp2=abs(real(sum(conjg(cwtl2r(:,:,nchorizo))*cwtend(:,:))))
        !tmp=(tmp1+tmp2)/2  
        !tmp3=abs(  (pathcoreval-tmp)/( (pathcoreval+tmp)/2  )        )*100
        !write(6,*) 'pathcoreval check stdmove2v6   :', pathcoreval, tmp,tmp1,tmp2,tmp3
        !if (  tmp3 >= 0.1   ) then
        !    write(6,*) '################## Caution ###################'
        !    write(6,*) 'pathcoreval check stdmove2v6    failed % err = ', tmp3
        !    write(6,*) '###############################################'
        !    !stop
        !endif
     !!!!!!!!!!!!!!!!!!      
    
    
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
    real(kind=r8) :: rt2(0:nchorizo),tmp,tmp1,tmp2,tmp3,logratio,tmp1test1,gc
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin) 

    if ((ileft.lt.0).or.(iright.gt.nchorizo).or.(iright.lt.0).or.(ileft.gt.nchorizo).or.(ileft.gt.iright)) then
	    write (6,*) 'stdmove1fastv6 range error', ileft,iright
	    stop
    else
    if (ileft.eq.iright) then
        call stdmove1fastv6(ileft,iright) 
        return
    endif          
    endif
 
    gc=-0.5_r8/(dt*(2*hbar))
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
	    ristranew(i)%x(:,:)=ristraold(i)%x(:,:)+gauss(:,:)    
    enddo ! just the shift. all beads shifted by the same amount.
   
    tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtl2r(:,:,0))*cwtr2l(:,:,1)))))
    
    
    !tmp1test1=log(abs(real(sum(conjg(cwtl2r(:,:,0))*cwtr2l(:,:,1)))))
    !write(6,'(''3v6 tmp1,tmp1test1:'',t30,i5,3(f15.7,1x))') ,nblocknow, abs(tmp1-tmp1test1),tmp1,tmp1test1     
    !if (abs(tmp1test1-tmp1)>=1.0d-7) stop  
    
    tmp2=tmp1   
    rt2=0.
    cwtl2rnew(:,:,:)=cwtl2r(:,:,:)
    cwtr2lnew(:,:,:)=cwtr2l(:,:,:)
   
    do i=ileft,iright  
   
        if (i.eq.0) then
		 
	        !call psitcwt(xold,cwtold)  
	        call corop(ristranew(i)%x,iplo,jplo,cwtbeginnew)
	        call v6propr(ristranew(i)%x,-1,cwtbeginnew(:,:),cwtl2rnew(:,:,i)) ! just update the i.
      
	        xd2new(:,:)=(ristranew(i)%x(:,:)-ristranew(i+1)%x(:,:))**2
	        xd2old(:,:)=(ristraold(i)%x(:,:)-ristraold(i+1)%x(:,:))**2
	        rt2(i)=gc*(sum(xd2new)-sum(xd2old))	  
            !rt2=0. ! vmc(dt=0) only.
	     
	    else
		 
	        if (i.eq.nchorizo) then   
	            !call psitcwt(xold,cwtold)

	            !call psitcwt(xnew,cwtnew)
	            call corop(ristranew(i)%x,ipro,jpro,cwtendnew)
                call vlsprop(spinorbit,-1,ristranew(i-1)%x,ristranew(i)%x,cwtl2rnew(:,:,i-1),cwttmp(:,:))
	            call v6propl(ristranew(i)%x,-1,cwttmp(:,:),cwtl2rnew(:,:,i))
                call v6propr(ristranew(i)%x,-1,cwtendnew,cwtr2lnew(:,:,i))
         
            else
	   
                call vlsprop(spinorbit,-1,ristranew(i-1)%x,ristranew(i)%x,cwtl2rnew(:,:,i-1),cwttmp(:,:))
	            call v6proplr(ristranew(i)%x,ristranew(i)%x,-1,cwttmp(:,:),cwtl2rnew(:,:,i))
	    
                xd2new(:,:)=(ristranew(i)%x(:,:)-ristranew(i+1)%x(:,:))**2
	            xd2old(:,:)=(ristraold(i)%x(:,:)-ristraold(i+1)%x(:,:))**2
	            rt2(i)=gc*(sum(xd2new)-sum(xd2old))	
                !rt2=0.  ! vmc(dt=0) only.

	        endif 
	    endif 	 
    enddo 

    if (iright /= nchorizo) then
        call vlsprop(spinorbit,-1,ristranew(iright)%x,ristranew(iright+1)%x,cwtl2rnew(:,:,iright),cwttmp(:,:))
        tmp2=log(abs(real(sum(conjg(cwttmp(:,:))*cwtr2lnew(:,:,iright+1)))))
    else
        tmp2=log(abs(real(sum(conjg(cwtl2rnew(:,:,iright))*cwtendnew(:,:)))))   
    endif
   
    if (ileft /= 0) then
        xd2new(:,:)=(ristranew(i)%x(:,:)-ristranew(i-1)%x(:,:))**2
	    xd2old(:,:)=(ristraold(i)%x(:,:)-ristraold(i-1)%x(:,:))**2         
        rt2(ileft-1)=gc*(sum(xd2new)-sum(xd2old))	
    endif
   
    logratio=tmp2-tmp1+sum(rt2)
    rn=randn(1)
   
    if (log(rn(1)).lt.logratio) then 
        icshft=icshft+1       
        do i=ileft,iright
            ristra(i)=ristranew(i)
            cwtl2r(:,:,i)=cwtl2rnew(:,:,i)
        enddo    
        if (ileft.eq.0) cwtbegin(:,:)=cwtbeginnew(:,:)
        if (iright.eq.nchorizo) cwtend(:,:)=cwtendnew(:,:)      
        pathcoreval=exp(tmp2)    
    ! update the chain.
        xl=ristra(0)%x
        xr=ristra(nchorizo)%x
        ! update cwtl2r 
        if (iright <= nchorizo-2) then 
            do k=iright,nchorizo-2
                x=ristra(k+1)%x
                call vlsprop(spinorbit,-1,ristra(k)%x,x,cwtl2r(:,:,k),cwttmp(:,:))
                call v6proplr(x,x,-1,cwttmp(:,:),cwtl2r(:,:,k+1))   
            enddo
            call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,nchorizo-1),cwttmp(:,:))
            call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))	      
        else  
            if (iright == nchorizo-1) then   
                call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,iright),cwttmp(:,:))
                call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))	          
            endif 
        endif     
        ! update cwtr2l from iright th bead until the 0th.
        if (iright == nchorizo) then    
            !x=ristra(iright)%x   
	        call v6propr(xr,-1,cwtend(:,:),cwtr2l(:,:,iright))
	        do k=iright-1,1,-1
	            x=ristra(k)%x
                call vlsprop(spinorbit,1,x,ristra(k+1)%x,cwtr2l(:,:,k+1),cwttmp(:,:))
	            call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,k))
            enddo
	        !xl=ristra(0)%x
            call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
	        call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	    
        else      
            if (iright >= 2) then  
                x=ristra(iright)%x
                call vlsprop(spinorbit,1,x,ristra(iright+1)%x,cwtr2l(:,:,iright+1),cwttmp(:,:))
                call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,iright))
                do k=iright-1,1,-1
	                x=ristra(k)%x
                    call vlsprop(spinorbit,1,x,ristra(k+1)%x,cwtr2l(:,:,k+1),cwttmp(:,:))
	                call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,k))
                enddo
	            !xl=ristra(0)%x
                call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
	            call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	   
            else ! iright =1 or 0     
                if (iright == 1) then   
                    x=ristra(iright)%x
                    call vlsprop(spinorbit,1,x,ristra(iright+1)%x,cwtr2l(:,:,iright+1),cwttmp(:,:))
                    call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,iright))   
	                !xl=ristra(0)%x
                    call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,iright),cwttmp(:,:))
	                call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	
                else ! iright =0
                    call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
                    call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	   
                endif   
            endif  
        endif     
    endif

     !!!!!!!!!!!! check   
        !tmp1=abs(real(sum(conjg(cwtbegin(:,:))*cwtr2l(:,:,0))))
        !tmp2=abs(real(sum(conjg(cwtl2r(:,:,nchorizo))*cwtend(:,:))))
        !tmp=(tmp1+tmp2)/2  
        !tmp3=abs(  (pathcoreval-tmp)/( (pathcoreval+tmp)/2  )        )*100
        !write(6,*) 'pathcoreval check stdmove3v6   :', pathcoreval, tmp,tmp1,tmp2,tmp3
        !if (  tmp3 >= 0.1   ) then
        !    write(6,*) '################## Caution ###################'
        !    write(6,*) 'pathcoreval check stdmove3v6    failed % err = ', tmp3
        !    write(6,*) '###############################################'
        !    !stop
        !endif
     !!!!!!!!!!!!!!!!!!  

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
    real(kind=r8) :: rt2(0:nchorizo),tmp,tmp1,tmp2,tmp3,logratio,tmp1test1,gc
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin) 

    if ((ileft.lt.0).or.(iright.gt.nchorizo).or.(iright.lt.0).or.(ileft.gt.nchorizo).or.(ileft.gt.iright)) then
	    write (6,*) 'stdmove1fastv6 range error', ileft,iright
	    stop
    else
    if (ileft.eq.iright) then
        call stdmove1fastv6(ileft,iright) 
        return
    endif          
    endif

    gc=-0.5_r8/(dt*(2*hbar))
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
	    ristranew(i)%x(:,:)=ristraold(i)%x(:,:)+gauss1(:,:)+gauss(:,:)    
    enddo ! just the shift. all beads shifted by the same amount.
   
    tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtl2r(:,:,0))*cwtr2l(:,:,1)))))
    
    !tmp1test1=log(abs(real(sum(conjg(cwtl2r(:,:,0))*cwtr2l(:,:,1)))))
    !write(6,'(''23v6 tmp1,tmp1test1:'',t30,i5,3(f15.7,1x))') ,nblocknow, abs(tmp1-tmp1test1),tmp1,tmp1test1     
    !if (abs(tmp1test1-tmp1)>=1.0d-7) stop 
    
    
    
    tmp2=tmp1   
    rt2=0.
    cwtl2rnew(:,:,:)=cwtl2r(:,:,:)
    cwtr2lnew(:,:,:)=cwtr2l(:,:,:)
   
    do i=ileft,iright  
   
        if (i.eq.0) then
		 
	    !call psitcwt(xold,cwtold)  
	    call corop(ristranew(i)%x,iplo,jplo,cwtbeginnew)
	    call v6propr(ristranew(i)%x,-1,cwtbeginnew(:,:),cwtl2rnew(:,:,i)) ! just update the i.
      
	    xd2new(:,:)=(ristranew(i)%x(:,:)-ristranew(i+1)%x(:,:))**2
	    xd2old(:,:)=(ristraold(i)%x(:,:)-ristraold(i+1)%x(:,:))**2
	    rt2(i)=gc*(sum(xd2new)-sum(xd2old))	  
        !rt2=0. ! vmc(dt=0) only.
	     
	    else
		 
	    if (i.eq.nchorizo) then   
	        !call psitcwt(xold,cwtold)

	        !call psitcwt(xnew,cwtnew)
	        call corop(ristranew(i)%x,ipro,jpro,cwtendnew)             
            call vlsprop(spinorbit,-1,ristranew(i-1)%x,ristranew(i)%x,cwtl2rnew(:,:,i-1),cwttmp(:,:))
	        call v6propl(ristranew(i)%x,-1,cwttmp(:,:),cwtl2rnew(:,:,i))
            call v6propr(ristranew(i)%x,-1,cwtendnew,cwtr2lnew(:,:,i))
         
	    else
	   
            call vlsprop(spinorbit,-1,ristranew(i-1)%x,ristranew(i)%x,cwtl2rnew(:,:,i-1),cwttmp(:,:))
	        call v6proplr(ristranew(i)%x,ristranew(i)%x,-1,cwttmp(:,:),cwtl2rnew(:,:,i))
	    
            xd2new(:,:)=(ristranew(i)%x(:,:)-ristranew(i+1)%x(:,:))**2
	        xd2old(:,:)=(ristraold(i)%x(:,:)-ristraold(i+1)%x(:,:))**2
	        rt2(i)=gc*(sum(xd2new)-sum(xd2old))	
            !rt2=0.  ! vmc(dt=0) only.

	    endif 
	    endif 	 
    enddo 

    if (iright /= nchorizo) then
        call vlsprop(spinorbit,-1,ristranew(iright)%x,ristranew(iright+1)%x,cwtl2rnew(:,:,iright),cwttmp(:,:))
        tmp2=log(abs(real(sum(conjg(cwttmp(:,:))*cwtr2lnew(:,:,iright+1)))))
    else
        tmp2=log(abs(real(sum(conjg(cwtl2rnew(:,:,iright))*cwtendnew(:,:)))))   
    endif
   
    if (ileft/=0) then
        xd2new(:,:)=(ristranew(i)%x(:,:)-ristranew(i-1)%x(:,:))**2
	    xd2old(:,:)=(ristraold(i)%x(:,:)-ristraold(i-1)%x(:,:))**2       
        rt2(ileft-1)=gc*(sum(xd2new)-sum(xd2old))	
    endif
   
    logratio=tmp2-tmp1+sum(rt2)
    rn=randn(1)
   
    if (log(rn(1)).lt.logratio) then 
 
        ic23=ic23+1   
       
        do i=ileft,iright
            ristra(i)=ristranew(i)
            cwtl2r(:,:,i)=cwtl2rnew(:,:,i)
        enddo 
      
        if (ileft.eq.0) cwtbegin(:,:)=cwtbeginnew(:,:)
        if (iright.eq.nchorizo) cwtend(:,:)=cwtendnew(:,:)  
        
        pathcoreval=exp(tmp2)
       
    ! update the chain.
        xl=ristra(0)%x
        xr=ristra(nchorizo)%x
    ! update cwtl2r 
        if (iright <= nchorizo-2) then 
            do k=iright,nchorizo-2
                x=ristra(k+1)%x
                call vlsprop(spinorbit,-1,ristra(k)%x,x,cwtl2r(:,:,k),cwttmp(:,:))
                call v6proplr(x,x,-1,cwttmp(:,:),cwtl2r(:,:,k+1))   
            enddo
            call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,nchorizo-1),cwttmp(:,:))
            call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))	      
        else  
            if (iright == nchorizo-1) then
                call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,iright),cwttmp(:,:))
                call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))	          
            endif 
        endif     
        ! update cwtr2l from iright th bead until the 0th.
        if (iright == nchorizo) then    
            !x=ristra(iright)%x   
	        call v6propr(xr,-1,cwtend(:,:),cwtr2l(:,:,iright))
	        do k=iright-1,1,-1
	            x=ristra(k)%x
                call vlsprop(spinorbit,1,x,ristra(k+1)%x,cwtr2l(:,:,k+1),cwttmp(:,:))
	            call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,k))
            enddo
	        !xl=ristra(0)%x
            call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
	        call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	    
        else      
            if (iright >= 2) then  
                x=ristra(iright)%x
                call vlsprop(spinorbit,1,x,ristra(iright+1)%x,cwtr2l(:,:,iright+1),cwttmp(:,:))
                call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,iright))
                do k=iright-1,1,-1
	                x=ristra(k)%x
                    call vlsprop(spinorbit,1,x,ristra(k+1)%x,cwtr2l(:,:,k+1),cwttmp(:,:))
	                call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,k))
                enddo
	            !xl=ristra(0)%x
                call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
	            call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	   
            else ! iright =1 or 0        
                if (iright == 1) then   
                    x=ristra(iright)%x
                    call vlsprop(spinorbit,1,x,ristra(iright+1)%x,cwtr2l(:,:,iright+1),cwttmp(:,:))
                    call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,iright))   
	                !xl=ristra(0)%x
                    call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,iright),cwttmp(:,:))
	                call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	
                else ! iright =0
                    call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
                    call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))	   
                endif 
            endif  
        endif 
   
       
    endif

    
    
     !!!!!!!!!!!! check   
        !tmp1=abs(real(sum(conjg(cwtbegin(:,:))*cwtr2l(:,:,0))))
        !tmp2=abs(real(sum(conjg(cwtl2r(:,:,nchorizo))*cwtend(:,:))))
        !tmp=(tmp1+tmp2)/2  
        !tmp3=abs(  (pathcoreval-tmp)/( (pathcoreval+tmp)/2  )        )*100
        !write(6,*) 'pathcoreval check stdmove23v6   :', pathcoreval, tmp,tmp1,tmp2,tmp3
        !if (  tmp3 >= 0.1   ) then
        !    write(6,*) '################## Caution ###################'
        !    write(6,*) 'pathcoreval check stdmove23v6    failed % err = ', tmp3
        !    write(6,*) '###############################################'
        !    !stop
        !endif
     !!!!!!!!!!!!!!!!!!     
    
    

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
    iright=nchorizo ! nchorizo is set to 1 for VMC.
   
    if ((ileft.lt.0).or.(iright.gt.nchorizo).or.(iright.lt.0).or.(ileft.gt.nchorizo).or.(ileft.gt.iright)) then
	    write (6,*) 'stdmove1fastv6 range error', ileft,iright
	    stop      
    endif

    ! now, ileft<iright   
   
    ristraold(0:nchorizo)=ristra(0:nchorizo)
    icstdmove2tot=icstdmove2tot+1   
   
    gauss=mov2step*reshape((randn(3*npart)-0.5),(/3,npart/)) ! this is not gauss move. for vmc(dt=0) only
    do i=ileft,iright   
	    ristranew(i)%x(:,:)=ristraold(i)%x(:,:)+gauss(:,:)    
    enddo
   
    tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtl2r(:,:,0))*cwtr2l(:,:,1)))))
    tmp2=tmp1   
    rt2(:)=0
    cwtl2rnew(:,:,:)=cwtl2r(:,:,:)
    cwtr2lnew(:,:,:)=cwtr2l(:,:,:)
   
    do i=ileft,iright  
   
        if (i.eq.0) then
		 
	        !call psitcwt(xold,cwtold)  
	        call corop(ristranew(i)%x,iplo,jplo,cwtbeginnew)
	        call v6propr(ristranew(i)%x,-1,cwtbeginnew(:,:),cwtl2rnew(:,:,i)) ! just update the i. 
	        !xd2new(:,:)=(ristranew(i)%x(:,:)-ristranew(i+1)%x(:,:))**2
	        !xd2old(:,:)=(ristraold(i)%x(:,:)-ristraold(i+1)%x(:,:))**2
	        !rt2(i)=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	  
            !rt2=0. ! vmc(dt=0) only.
	     
	    else
		 
	        if (i.eq.nchorizo) then   
	            !call psitcwt(xold,cwtold)
	            !call psitcwt(xnew,cwtnew)
	            call corop(ristranew(i)%x,ipro,jpro,cwtendnew)
	            call v6propl(ristranew(i)%x,-1,cwtl2rnew(:,:,i-1),cwtl2rnew(:,:,i))
                call v6propr(ristranew(i)%x,-1,cwtendnew,cwtr2lnew(:,:,i))
         
	        else
	   
	            call v6proplr(ristranew(i)%x,ristranew(i)%x,-1,cwtl2rnew(:,:,i-1),cwtl2rnew(:,:,i))
                !xd2new(:,:)=(ristranew(i)%x(:,:)-ristranew(i+1)%x(:,:))**2
	            !xd2old(:,:)=(ristraold(i)%x(:,:)-ristraold(i+1)%x(:,:))**2
	            !rt2(i)=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	
                !rt2=0.  ! vmc(dt=0) only.

	        endif 
	    endif 	 
    enddo 

    tmp2=log(abs(real(sum(conjg(cwtl2rnew(:,:,iright))*cwtendnew(:,:)))))   

   
    logratio=tmp2-tmp1+sum(rt2)
    rn=randn(1)
   
    if (log(rn(1)).lt.logratio) then 
       
        icstdmove2=icstdmove2+1
       
        pathcoreval=exp(tmp2)
        do i=ileft,iright
        ristra(i)=ristranew(i)
        cwtl2r(:,:,i)=cwtl2rnew(:,:,i)
        enddo 
      
        if (ileft.eq.0) cwtbegin(:,:)=cwtbeginnew(:,:)
        if (iright.eq.nchorizo) cwtend(:,:)=cwtendnew(:,:)  
       
        ! update the chain.
            xl=ristra(0)%x
            xr=ristra(nchorizo)%x
        ! update cwtl2r 
        if (iright <= nchorizo-2) then 
            do k=iright,nchorizo-2
                x=ristra(k+1)%x
                call v6proplr(x,x,-1,cwtl2r(:,:,k),cwtl2r(:,:,k+1))   
            enddo
            call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))	      
        else  
            if (iright == nchorizo-1) then   
                call v6propl(xr,-1,cwtl2r(:,:,iright),cwtl2r(:,:,nchorizo))	          
            endif 
        endif     
        ! update cwtr2l from iright th bead until the 0th.
        if (iright == nchorizo) then    
            !x=ristra(iright)%x   
	        call v6propr(xr,-1,cwtend(:,:),cwtr2l(:,:,iright))
            if ((iright-1)>=1) then
	            do k=iright-1,1,-1
	                x=ristra(k)%x
	                call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
                enddo
            endif
	        !xl=ristra(0)%x
	        call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	    
        else      
            if (iright >= 2) then  
                x=ristra(iright)%x
                call v6proplr(x,x,-1,cwtr2l(:,:,iright+1),cwtr2l(:,:,iright))
                do k=iright-1,1,-1
	            x=ristra(k)%x
	            call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
	            enddo
	            !xl=ristra(0)%x
	            call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	   
            else ! iright =1 or 0
          
                if (iright == 1) then   
                x=ristra(iright)%x
                call v6proplr(x,x,-1,cwtr2l(:,:,iright+1),cwtr2l(:,:,iright))   
	            !xl=ristra(0)%x
	            call v6propl(xl,-1,cwtr2l(:,:,iright),cwtr2l(:,:,0))	
                else ! iright =0
                call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	   
                endif
       
            endif  
        endif
    
        call addval(2,1.0_r8,1.0_r8)  
    else
        call addval(2,0.0_r8,1.0_r8)   
    endif
   

    return
    end subroutine stdvmcmove2v6 
    
    subroutine corrchk(x,x1,x2,ixd,corrin,corr1in,corr2in,lmax)
! Interacting Electrons, Ceperley et al, p585. MC in ab initio quantum chemstry, Hammond et al, p59, eq(2.27).   
    use wavefunction
    use estimator
    use random
    integer(kind=i4) :: ixd,lmax,lcorr,icorr,j! lmax <= ixd, usally set ixd=lmax=nstepdecor
    real(kind=r8) :: xave,xvar,x1ave,x1var,x2ave,x2var
    real(kind=r8) :: x(:),corrin(:),x1(:),x2(:),corr1in(:),corr2in(:),kappa(ixd-1),kappa1(ixd-1),kappa2(ixd-1)
    integer(kind=i4) :: iounit
    logical :: findstepdecorr 
    ! nstepdecor usually set the same as nstep when do the correlation check.
    xave=sum(x(1:ixd))/ixd
    xvar=sum(x(1:ixd)**2)/ixd-xave**2 
    x1ave=sum(x1(1:ixd))/ixd
    x1var=sum(x1(1:ixd)**2)/ixd-x1ave**2
    x2ave=sum(x2(1:ixd))/ixd
    x2var=sum(x2(1:ixd)**2)/ixd-x2ave**2    
    findstepdecorr=.false.
    corrin=0
    corr1in=0
    corr2in=0
    kappa=0
    OPEN(newunit=iounit,FILE='xcorrelationcheck.txt',FORM='FORMATTED')
    WRITE(iounit,'(2x,8(G15.8,2x))') 'k','corr(k)','kappa','corr1(k)','kappa1','corr2(k)','kappa2'
    do lcorr=1,ixd-1 ! k     
	    do icorr=1,ixd-lcorr ! n or i.
		    corrin(lcorr)=corrin(lcorr)+ (x(icorr)-xave)*(x(icorr+lcorr)-xave)
            corr1in(lcorr)=corr1in(lcorr)+ (x1(icorr)-x1ave)*(x1(icorr+lcorr)-x1ave)
            corr2in(lcorr)=corr2in(lcorr)+ (x2(icorr)-x2ave)*(x2(icorr+lcorr)-x2ave)
	    enddo
	    corrin(lcorr)=corrin(lcorr)/(ixd-lcorr)/xvar 
        kappa(lcorr)=1+2*sum(corrin(1:lcorr))
	    corr1in(lcorr)=corr1in(lcorr)/(ixd-lcorr)/x1var 
        kappa1(lcorr)=1+2*sum(corr1in(1:lcorr))
	    corr2in(lcorr)=corr2in(lcorr)/(ixd-lcorr)/x2var 
        kappa2(lcorr)=1+2*sum(corr2in(1:lcorr))
        WRITE(iounit,'(2x,20(G15.8,2x))') lcorr,corrin(lcorr),kappa(lcorr) &
                                         ,corr1in(lcorr),kappa1(lcorr),corr2in(lcorr),kappa2(lcorr)
    enddo    
    close(iounit) 
	OPEN(newunit=iounit,FILE='xcorrchk.txt',FORM='FORMATTED')
    WRITE(iounit,'(2x,2(G15.8,2x))') 'k','x','x1','x2' 
    do j=1,lmax 
	    WRITE(iounit,'(2x,5(G15.8,2x))') j,x(j),x1(j),x2(j)	
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
    if (nrepmax.eq.0) return
    repstepct1=0
    repstepct2=0
    repstepct1sum=0
    repstepct2sum=0  
    icrephist=0
    icrephistlast=0
    return
    end subroutine zerepstepmon
   
    subroutine counterinit ! zero the stats after finishing one block.
    icstepstd=0
    ice=0
    ic0=0
    ic0tot=0
    icn=0
    icntot=0
    icmov1=0
    icmov1rvs=0
    icmov1tot=0
    icmov1rvstot=0
    icstdmove1tot=0
    icstdmove1rvstot=0
    
    icstdmove1bisect1rvstot=0
    icmov1bisect1rvstot=0
    icstdmove1bisect1tot=0
    icmov1bisect1tot=0    
    icmov1bisect1rvs=0
    icmov1bisect1=0
    
    icstdmove2=0
    icstdmove2tot=0
    icshft=0
    icshfttot=0
    ic23=0
    ic23tot=0
    ibisect=0
    ibisecttot=0
    ibisecthalf=0
    ibisecthalftot=0
    ibisectextra=0
    ibisectextratot=0
    
    
    
    icbisecttot=0
    bisectrate=0
    bisectcount=0
    
    bisectextracount=0
    bisectextrarate=0
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
    use random
    integer(kind=i4) :: myunit
    real(kind=r8) :: r(20)
    character(len=120) :: filename    
    real(kind=r8) :: xtot(3,npart,0:nchorizo)
    integer(kind=i4) :: ipl(npair),jpl(npair),ipr(npair),jpr(npair)
    integer(kind=i8) :: irnout
    
    call statcalc(r)
      
    if (rbisect>0) then
        if ((ibisect.eq.0).or.(ibisectl.eq.0).or.(ibisectr.eq.0)) then
            nblockstuck=nblocknow
            ibisectstuck(nblockstuck)=1  
    ! optional   
            if ( nproc().eq.1 ) then
                call chorizoallout(xtot,ipl,jpl,ipr,jpr) 
                call showirn(irnout)
                write(filename,'("myrank",i10,1x,i10".stuckonce")') myrank(),nblockstuck
                open(newunit=myunit,form='formatted',file=trim(filename),position='rewind')  
                    write (myunit,*) 'stuck block',nblockstuck            
                    write (myunit,'(6i10)') ibisect,ibisecttot,ibisectl,ibisecttotl,ibisectr,ibisecttotr
                    write (myunit,*) 'irn seed is ',irnout 
                    write (myunit,*)  
                    write (myunit,'(''dt ='',t40,f10.5)') dt
                    write (myunit,'(''nchorizo ='',t40,i10)') nchorizo   
                    call printstat(myunit)        
                    write (myunit,'(<npair>i10)') ipl
                    write (myunit,'(<npair>i10)') jpl
                    write (myunit,'(<npair>i10)') ipr
                    write (myunit,'(<npair>i10)') jpr
                    write (myunit,'(3e15.7)') reshape(xtot,(/3*npart*(nchorizo+1)/))      
                close(myunit)    
            endif
    ! optional end            
            if ( nblockstuck >= 5 ) then 
                if (sum(ibisectstuck(nblockstuck-4:nblockstuck)).eq.5) then
                    call chorizoallout(xtot,ipl,jpl,ipr,jpr) 
                    call showirn(irnout)
                    write(filename,'("myrank",i10,".stuck")') myrank()
                    open(newunit=myunit,form='formatted',file=trim(filename),position='rewind')  
                        write (myunit,*) 'nsteps really stuck are ',nblockstuck            
                        write (myunit,'(6i10)') ibisect,ibisecttot,ibisectl,ibisecttotl,ibisectr,ibisecttotr
                        write (myunit,*) 'irn seed is ',irnout 
                
                        call printstat(myunit)
                
                        write (myunit,'(<npair>i10)') ipl
                        write (myunit,'(<npair>i10)') jpl
                        write (myunit,'(<npair>i10)') ipr
                        write (myunit,'(<npair>i10)') jpr
                        write (myunit,'(3e15.7)') reshape(xtot,(/3*npart*(nchorizo+1)/))      
                    close(myunit)
                    write(6,*) 'stuck configuration, abort! myrank= ', myrank()
                    write(12,*) 'stuck configuration, abort! myrank= ', myrank()
                    call abort                         
                endif       
            endif
        endif           
    endif    
        
    return
    end subroutine checkblkstat
 
    subroutine showstat
    use mympi   
    write (6,'(''dt ='',t40,g15.7)') dt
    write (6,'(''nchorizo ='',t40,i10)') nchorizo 
    call printstat(6)
    call printstat(12)
    return
    end subroutine showstat 
    
    subroutine printstat(unit)
    use mympi
    integer(kind=i4) :: unit,k
    real(kind=r8) :: r(20)
    call statcalc(r)
    write (unit,*) 
    write (unit,'(''myrank = '',i10)') myrank()
    write (unit,'(''acceptance (mov1bisect1 ->) ='',t40,f10.5,3(i10))') r(12) &
                                                   ,icmov1bisect1,icmov1bisect1tot,icstdmove1bisect1tot
    write (unit,'(''acceptance (mov1bisect1 <-) ='',t40,f10.5,3(i10))') r(13) &
                                                   ,icmov1bisect1rvs,icmov1bisect1rvstot,icstdmove1bisect1rvstot 
    write (unit,'(''acceptance (mov1 ->) ='',t40,f10.5,3(i10))') r(1),icmov1,icmov1tot,icstdmove1tot
    write (unit,'(''acceptance (mov1 <-) ='',t40,f10.5,3(i10))') r(11),icmov1rvs,icmov1rvstot,icstdmove1rvstot 
    write (unit,'(''acceptance (mov2 overall) ='',t40,f10.5,2(i10))') r(2),icstdmove2,icstdmove2tot
    write (unit,'(''acceptance (mov3 shift) ='',t40,f10.5,2(i10))') r(3),icshft,icshfttot
    write (unit,'(''acceptance ratio (stdmov23 mov2+3) ='',t40,f10.5,2(i10))') r(4),ic23,ic23tot
    write (unit,'(''icn ratio ='',t40,f10.5,2(i10))') r(5),icn,icntot
    write (unit,'(''lr ratio ='',t40,f10.5,2(i10))') r(6),iclr,iclrtot
    write (unit,'(''acceptance reptation ='',t40,f10.5,3(i10))') r(7),icrep,icreptot
    write (unit,'(''acceptance bisection ='',t40,f10.5,3(i10))') r(8),ibisect,ibisecttot,icbisecttot
    write (unit,'(''acceptance bisectionhalf ='',t40,f10.5,3(i10))') r(15),ibisecthalf,ibisecthalftot
    write (unit,'(''acceptance bisection extra ='',t40,f10.5,2(i10))') r(14),ibisectextra,ibisectextratot
    
    write (unit,'(''bisection level 1 to mmax ratio ='',t40,<mmax>(f5.3,1x))') (bisectrate(k),k=1,mmax)
    write (unit,'(''acceptance bisection L end ='',t40,f10.5,3(i10))') r(9),ibisectl,ibisecttotl
    !write (unit,'(''bisection lend level 0 to mmaxnow ratio ='',t40,f10.5)') (bisectcountl(k)/bisecttotl(k),k=0,mmax-1)
    write (unit,'(''acceptance bisection R end ='',t40,f10.5,3(i10))') r(10),ibisectr,ibisecttotr
    !write (unit,'(''bisection rend level 0 to mmaxnow ratio ='',t40,f10.5)') (bisectcountr(k)/bisecttotr(k),k=0,mmax-1)
    write (unit,'(''energy samples per block per core ='',t40,i10)') ice
    write (unit,*) 	  
    return
    end subroutine printstat     
    
    subroutine statcalc(r)
    real(kind=r8) :: r(:) 
    
    bisectrate(:)=bisectcount(:)/bisecttot(:)
    bisectextrarate(:)=bisectextracount(:)/bisectextratot(:)    
    r(1)=dble(icmov1)/dble(icmov1tot)
    r(2)=dble(icstdmove2)/dble(icstdmove2tot)
    r(3)=dble(icshft)/dble(icshfttot)
    r(4)=dble(ic23)/dble(ic23tot)
    r(5)=dble(icn)/dble(icntot)
    r(6)=dble(iclr)/dble(iclrtot)
    r(7)=dble(icrep)/dble(icreptot)
    r(8)=dble(ibisect)/dble(ibisecttot)
    r(9)=dble(ibisectl)/dble(ibisecttotl)
    r(10)=dble(ibisectr)/dble(ibisecttotr)
    r(11)=dble(icmov1rvs)/dble(icmov1rvstot)
    r(12)=dble(icmov1bisect1)/dble(icmov1bisect1tot)
    r(13)=dble(icmov1bisect1rvs)/dble(icmov1bisect1rvstot)
    r(14)=dble(ibisectextra)/dble(ibisectextratot)
    r(15)=dble(ibisecthalf)/dble(ibisecthalftot)
    
    return
    end subroutine statcalc     
    
   
    subroutine bisectpickslices(ileftbisect,irightbisect,ipick,nbisectnow) 
    use random
    integer(kind=i4) :: ileftbisect,irightbisect,ipick,na,nb,nbisectnow
    real(kind=r8) :: rn(1)
    if (nbisectnow <= 1) stop ' nbisectnow <= 1 in bisectpickslices, stop! ' 
    rn=randn(1)
    ! need to initialize nbisect and nchorizo first.
    na=nbisectnow/2
    nb=nchorizo-na
    ipick=na+int(dble(nb-na+1)*rn(1)) ! I believe this picking strategy is good. 
    ileftbisect=ipick-na
    irightbisect=ipick+na 
    return    
    end subroutine bisectpickslices
   
    subroutine bisectpicksliceswide(ileftbisect,irightbisect,nbisectnow) 
    use random
    integer(kind=i4) :: ileftbisect,irightbisect,nbisectnow
    real(kind=r8) :: rn(1)
    rn=randn(1)
    ! need to initialize nbisect and nchorizo first.
    ileftbisect=min(int(dble(nchorizo+1)*rn(1)),nchorizo) ! min just to make sure, not really need.
    irightbisect=int(dble(nchorizo+nbisectnow-1)*rn(1))+1
    return    
    end subroutine bisectpicksliceswide    
      
 
   
    subroutine xoutput(x) ! initialize the position, calculate the weights for each bead. l2r,r2l
    real(kind=r8) :: x(3,npart)
    x=ristra(0)%x
    return
    end subroutine xoutput   
   
    subroutine getcwtr2l(ichorizoin,cwttmp)
    integer(kind=i4) :: ichorizoin
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin)
    cwttmp(:,:)=cwtr2l(:,:,ichorizoin)
    return
    end subroutine getcwtr2l
    
    subroutine getcwtl2r(ichorizoin,cwttmp)
    integer(kind=i4) :: ichorizoin
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin)
    cwttmp(:,:)=cwtl2r(:,:,ichorizoin)
    return
    end subroutine getcwtl2r 
    
    subroutine getcwtlrm(ichorizoin,cwtl2rm,cwtr2lm)
    use v6stepcalc
    integer(kind=i4) :: ichorizoin
    complex(kind=r8) :: cwtl2rl(0:nspin-1,nisospin),cwtr2lr(0:nspin-1,nisospin) &
                       ,cwtl2rm(0:nspin-1,nisospin),cwtr2lm(0:nspin-1,nisospin) 
    if (nchorizo<=1) then  
        call getcwtr2l(1,cwtr2lm)
        call getcwtl2r(0,cwtl2rm)
    else
        call getcwtr2l(ichorizoin,cwtr2lr)
        call getcwtl2r(ichorizoin,cwtl2rl)
        call v6proplpr(ristra(ichorizoin)%x,-1,cwtr2lr,cwtr2lm)
        call v6proplpr(ristra(ichorizoin)%x,-1,cwtl2rl,cwtl2rm)  
    endif    
    return
    end subroutine getcwtlrm
    
    
    
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

    subroutine vtot(i,v,vnn,vem) ! v is the numerator, in the end need the denominator.
    use brut
    use v6stepcalc   
    use chorizos
    complex(kind=r8) :: cwtr2lm(0:nspin-1,nisospin,0:nchorizo),cwtl2rm(0:nspin-1,nisospin,0:nchorizo) ! m means middle.
    real(kind=r8) :: v,vnn,vem
    integer(kind=i4) :: i
    real(kind=r8) :: x(3,npart) 
    if (.not.iem) vem=0
    x=ristra(i)%x 
    !call xconvert(x) 
    if (i.eq.0) then             
        call pecal(x,cwtbegin,cwtr2l(:,:,i),vnn)  
        if (iem) call empecal(x,cwtbegin,cwtr2l(:,:,i),vem)  
    else
        if (i.eq.nchorizo) then           
            call pecal(x,cwtl2r(:,:,i),cwtend,vnn)   
            if (iem) call empecal(x,cwtl2r(:,:,i),cwtend,vem)   
        else
    ! in the somewhere middle  
        call v6proplpr(x,-1,cwtr2l(:,:,i),cwtr2lm(:,:,i))
        call v6proplpr(x,-1,cwtl2r(:,:,i),cwtl2rm(:,:,i))
        call pecal(x,cwtl2rm(:,:,i),cwtr2lm(:,:,i),vnn)  
        if (iem) call empecal(x,cwtl2rm(:,:,i),cwtr2lm(:,:,i),vem)   
        endif   
    endif 
    v=vnn+vem
    return
    end subroutine vtot

    subroutine hpsi(c0,cn,valnum,valdenom,valnumke,valnumpev6,valnumpeem,valnumpels,psi20,psi2n,rf,icheck) ! calculate numerator and denominator.
    use brut
    use v6stepcalc   
    use chorizos
    use mympi
    type (chorizo) :: c0,cn
    complex(kind=r8) ::  cwta1(25,0:nspin-1,nisospin),cwta2(25,0:nspin-1,nisospin) 
    complex(kind=r8) :: f
    real(kind=r8) :: rf,valnum,valdenom,alr,blr,valnumpev6,valnumpeem,valnumke,valnumpels
    real(kind=r8) :: psi20,psi2n,e0,en,psi2,pe0,pen,empe0,empen,psi,ke0,ken,pels0,pelsn
    integer(kind=i4) :: icheck ! 0=pass,1=fail.
    real(kind=r8) :: x0(3,npart),xn(3,npart)
 
    !cwtr2lout=cwtr2l
    !cwtl2rout=cwtl2r
    x0=c0%x
    xn=cn%x
   
    !write(6,*) 'hpsi start!'
     
    ! add the cwt r2l l2r checking.

    call hpsitcwt(x0,iplo,jplo,cwtr2l(:,:,0),psi20,e0,ke0,pe0,empe0,pels0,cwta1)
     
    !call hpsitcwt(x0,iplo,jplo,cwtr2l(:,:,0),psi20,e0,d20,pe0,f0,cwta0) 
   
    !write(6,*) 'cwta1 - cwta0 =', sum(abs(cwta1(1,:,:)-cwta0(1,:,:)))!,(cwta1(1,:,:)-cwta0(1,:,:))   
  
    !call corop(x0,iplo,jplo,cwt)   
    !write(6,*) 'cwta1 - cwtbegin =', sum(abs(cwta1(1,:,:)-cwt(:,:)))!,(cwta1(1,:,:)-cwt(:,:))
    !write(6,*) 'cwta1 =', cwta1(1,:,:)
    !write(6,*) 'cwtbegin =', cwt
 
    call hpsitcwt(xn,ipro,jpro,cwtl2r(:,:,nchorizo),psi2n,en,ken,pen,empen,pelsn,cwta2)
    !call corop(xn,ipro,jpro,cwt)  
    !write(6,*) 'cwta2 - cwtbegin =', sum(abs(cwta2(1,:,:)-cwt(:,:)))!,(cwta2(1,:,:)-cwt(:,:))       
    !write(6,*) 'cwta2 =', cwta2(1,:,:)
    !write(6,*) 'cwtend =', cwt  
   
    psi2=psi20+psi2n
    psi=psi2/2  ! psi is real f. ! psi20 and psi2n should be the same.
   
    
    !write(6,'('' e0,en,psi20,psi2n '',t30,4(g15.7,1x))') e0,en,psi20,psi2n
    
    
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
    !!write(6,'( ''cwtl2r ='', 12(g15.7,1x) )') cwtl2r(:,:,i) 
    !!write(6,*) 'nchorizo-i=',nchorizo-i    
    !!write(6,'( ''cwtr2l ='', 12(g15.7,1x) )') cwtr2l(:,:,nchorizo-i)
    !!write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    !
    !
    !write(19,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'   
    !write(19,*) 'check cwtr2l,l2r'
    !write(19,*) 'i=',i
    !write(19,'( ''cwtl2r ='', 12(g15.7,1x) )') cwtl2r(:,:,i) 
    !write(19,*) 'nchorizo-i=',nchorizo-i    
    !write(19,'( ''cwtr2l ='', 12(g15.7,1x) )') cwtr2l(:,:,nchorizo-i)
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
    valnumpels=0.5*(pels0+pelsn)/abs(psi)
    valnumke=alr-valnumpev6-valnumpeem-valnumpels 
   
    !write(6,*) 'e0=',e0
    !write(6,*) 'en=',en,(e0+en)
    !write(6,*) 'psi20, psi2n=',psi20,psi2n,psi,abs(psi)
    !write(6,*) 'valnum=',valnum,valnumke,alr,valnumpev6,valnumpeem
    !write(6,*) 'valdenom=',valdenom

    !rgl=e0
    !rgr=en
    !rg=(rgl+rgr)/2.
    !
    rf=psi
    !f=cmplx(psi) ! neglect the imaginary part actually. 
    !rf=real(f) ! f is the stuff without gaussian part bc gaussians are cancelled ususaly.
    !absrf=abs(real(f)) ! perhaps call it pathcore
  
    !write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'  
    !write(6,*) 'psi20 and psi2n=',psi20,psi2n  
    !write(12,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'  
    !write(12,*) 'psi20 and psi2n=',psi20,psi2n   
   
   
   
    if ( abs((psi20-psi2n)/((psi20+psi2n)/2))>=1.d-3   ) then
       
        icheck=1 ! fail.
       
    !write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'  
    !write(6,*) 'somthing wrong with psi20 and psi2n=',psi20,psi2n  
    !write(12,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'  
    !write(12,*) 'somthing wrong with psi20 and psi2n=',psi20,psi2n  	
    !write(12,*) 'iplo=',iplo
    !write(12,*) 'jplo=',jplo   
    !write(12,*) 'ipro=',ipro
    !write(12,*) 'jpro=',jpro      
    !write(12,*) 'cwtl2r(:,:,nchorizo)=',cwtl2r(:,:,nchorizo)     
    !write(12,*) 'cwtr2l(:,:,0)=',cwtr2l(:,:,0)

    ! scan and compare.
    
    ! write out x.      
	    !open(unit=9,form='formatted',file=trim(outfile))
    !   rewind 9
    !   do l=0,nchorizo
	    ! call chorizoout(x,l)
	    ! write (9,'(6e15.7)') x
    !   enddo
    !   write (9,'(<npair>i10)') iplo 
    !   write (9,'(<npair>i10)') jplo 
    !   write (9,'(<npair>i10)') ipro 
    !   write (9,'(<npair>i10)') jpro 
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
    integer(kind=i4) :: ipl(npair),jpl(npair),ipr(npair),jpr(npair)	
	    ipl=iplo
	    jpl=jplo
	    ipr=ipro
	    jpr=jpro	
    return
    end subroutine stepoutputordlr 
    
    subroutine updatecwtchain ! update the whole cwtr2l, cwtl2r, cwtbegin and cwtend.
    use wavefunction
    use v6stepcalc
    use brut 
    use random
    integer(kind=i4) :: i,j
    real(kind=r8) :: xr(3,npart),xl(3,npart),x(3,npart)
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin) 
   
    xl=ristra(0)%x 
    xr=ristra(nchorizo)%x
    ! l2r   
    call corop(xl,iplo,jplo,cwtbegin) ! psitcwt
    call v6propr(xl,-1,cwtbegin,cwtl2r(:,:,0))
    do i=1,nchorizo-1
        x=ristra(i)%x
        call vlsprop(spinorbit,-1,ristra(i-1)%x,x,cwtl2r(:,:,i-1),cwttmp(:,:))
        call v6proplr(x,x,-1,cwttmp(:,:),cwtl2r(:,:,i))
    enddo  
    call vlsprop(spinorbit,-1,ristra(nchorizo-1)%x,xr,cwtl2r(:,:,nchorizo-1),cwttmp(:,:))
    call v6propl(xr,-1,cwttmp(:,:),cwtl2r(:,:,nchorizo))
    ! r2l 
    call corop(xr,ipro,jpro,cwtend)
    call v6propr(xr,-1,cwtend,cwtr2l(:,:,nchorizo))
    do j=nchorizo-1,1,-1
        x=ristra(j)%x
        call vlsprop(spinorbit,1,x,ristra(j+1)%x,cwtr2l(:,:,j+1),cwttmp(:,:))
        call v6proplr(x,x,-1,cwttmp(:,:),cwtr2l(:,:,j))
    enddo
    call vlsprop(spinorbit,1,xl,ristra(1)%x,cwtr2l(:,:,1),cwttmp(:,:))
    call v6propl(xl,-1,cwttmp(:,:),cwtr2l(:,:,0))

    return
    end subroutine updatecwtchain   
    
    subroutine calcpathcoreval(ristrain,cwtbeginout,cwtendout,cwtl2rout,cwtr2lout,pathcorevalout) 
    use wavefunction
    use v6stepcalc
    use brut 
    use chorizos
    integer(kind=i4) :: i,j
    real(kind=r8) :: xr(3,npart),xl(3,npart),x(3,npart),pathcorevalout
    type (chorizo) :: ristrain(0:nchorizo)
    complex(kind=r8) :: cwttmp(0:nspin-1,nisospin),cwtbeginout(0:nspin-1,nisospin),cwtendout(0:nspin-1,nisospin) &
                        ,cwtl2rout(0:nspin-1,nisospin,0:nchorizo),cwtr2lout(0:nspin-1,nisospin,0:nchorizo)
   
    xl=ristrain(0)%x 
    xr=ristrain(nchorizo)%x
    ! l2r   
    call corop(xl,iplo,jplo,cwtbeginout) ! psitcwt
    call v6propr(xl,-1,cwtbeginout,cwtl2rout(:,:,0))
    do i=1,nchorizo-1
        x=ristrain(i)%x
        call vlsprop(spinorbit,-1,ristrain(i-1)%x,x,cwtl2rout(:,:,i-1),cwttmp(:,:))
        call v6proplr(x,x,-1,cwttmp(:,:),cwtl2rout(:,:,i))
    enddo   
    call vlsprop(spinorbit,-1,ristrain(nchorizo-1)%x,xr,cwtl2rout(:,:,nchorizo-1),cwttmp(:,:))
    call v6propl(xr,-1,cwttmp(:,:),cwtl2rout(:,:,nchorizo))
    ! r2l 
    call corop(xr,ipro,jpro,cwtendout)
    call v6propr(xr,-1,cwtendout,cwtr2lout(:,:,nchorizo))
    do j=nchorizo-1,1,-1
        x=ristrain(j)%x
        call vlsprop(spinorbit,1,x,ristrain(j+1)%x,cwtr2lout(:,:,j+1),cwttmp(:,:))
        call v6proplr(x,x,-1,cwttmp(:,:),cwtr2lout(:,:,j))
    enddo
    call vlsprop(spinorbit,1,xl,ristrain(1)%x,cwtr2lout(:,:,1),cwttmp(:,:))
    call v6propl(xl,-1,cwttmp(:,:),cwtr2lout(:,:,0))
    pathcorevalout=abs(real( sum(conjg(cwtbeginout(:,:))*cwtr2lout(:,:,0)) )) 	 
    return
    end subroutine calcpathcoreval   

    
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
    integer(kind=i4) :: istart,idummy
    integer(kind=i4) :: iploc,ipnow,rknow
    logical :: ipathlocated    

    
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
    
    if (.not.resume) then
        istart=1
    else
        if (np>npo) then
            if (myrank().eq.0) write (6,*) 'np now > np old, resume mode cannot be done, program stop.'
            call abort
        else ! np must be <= npo
            call readcheckpoint(irep,idummy)
            call bcast(ipathcountold)
            istart=ipathcountold+1
            if (myrank().eq.0) write(6,*) 'istart = ',istart
        endif
    endif
         
    select case (irep)
           
    case (5:6)  
! the usual way, read in a path then calcualte. 
! 5 means calculate with angle averaged response functions. 
! 6 means calculate without angle averaged response functions. fast.
        
        if (myrank().eq.0) then   
            allocate(xallo(ixtemp*npo))
            allocate(iplallo(npair*npo),jplallo(npair*npo),iprallo(npair*npo),jprallo(npair*npo))   
            allocate(xalln(ixtemp*np))
            allocate(iplalln(npair*np),jplalln(npair*np),ipralln(npair*np),jpralln(npair*np))       
! locate the postion in path.unf first            
            if (.not.resume) then
                open(unit=19,form='unformatted',file=file_path,position='rewind')   
            else
            ! this scan mode may take long for big file. it depend on IO performace of storage system.  
            ! decide search from beginning or the end point.  
            ! call subroutine readcheckpoint first.
                
                ipathlocated=.false.    
                
                if ( ipathcountold < ipath/2 ) then ! search from beginning
                     write (6,*) 'search path.unf from the beginning, ',ipathcountold,ipath,ipath/2
                     open(unit=19,form='unformatted',file=file_path,position='rewind')   
                         do while (.not.ipathlocated)
                             do iploc=1,pathunfstep1*npo+1 
                                read (19) 
                             enddo                        
                             read (19) ipnow
                             backspace 19
                             if (ipnow.eq.istart) then
                                 ipathlocated=.true.
                                 write (6,*) ' istart successfully located for subroutine compute! ',ipnow
                                 exit
                             endif
                             if (ipnow.eq.ipathmark) then
                                write (6,*) ' reached bottom. previous compute already finished, stop ',ipnow
                                stop
                             endif
                         enddo  
                else  ! search from bottom
                    write (6,*) 'search path.unf from the bottom, ',ipathcountold,ipath,ipath/2
                    open(unit=19,form='unformatted',file=file_path,position='append')   
                        do while (.not.ipathlocated)               
                            do iploc=1,pathunfstep1*npo+1 
                                backspace 19
                            enddo                             
                            backspace 19
                            read (19) ipnow
                            if (ipnow.eq.istart) then
                                ipathlocated=.true.
                                write (6,*) ' istart successfully located for subroutine compute! ',ipnow
                                backspace 19
                                exit 
                            endif
                            if (ipnow.eq.1) then
                                write (6,*) ' reached top. previous compute already finished, stop ',ipnow
                                ! the path.unf may have been cut a little bit if resumed in PIMC mode to match ipathcountold.
                                ! So it may just stop at ipathcountold instead of ipath.
                                stop                              
                            endif
                        enddo
                endif
                
            endif
            write (6,*) 'Calculations begin now! '  
        endif   
        
        if (n4 == 0) then ! npo # > np 
            n6=0
            if (m1>0) n6=np/m1  ! np<npo 
            if (myrank().eq.0) then   
                if (n6>1) then       
                    allocate(xall( (ixtemp*n1*np+1):ixtemp*npo,n6))
                    allocate(iplall((npair*n1*np+1):npair*npo,n6),jplall((npair*n1*np+1):npair*npo,n6) &
                            ,iprall((npair*n1*np+1):npair*npo,n6),jprall((npair*n1*np+1):npair*npo,n6))                         
                endif
            endif                   
            n6c=0          
            do i=istart,ipath 
                           
                if ((np <= npo).and.((i-1) >= 1)) then    
                    ipathcountold=i-1
                    if (myrank().eq.0) call writecheckpoint(irep,idummy)  ! save the (i-1)th path information                     
                endif               

                if (myrank()==0) then  ! read in the old path file.  
                    if (i==1) time0=mpi_wtime()
                    read (19) ipathnow  
                    do l=0,npo-1    
                        read (19) lnow
                        read (19) iplallo((npair*l+1):(npair*(l+1))) 
                        read (19) jplallo((npair*l+1):(npair*(l+1))) 
                        read (19) iprallo((npair*l+1):(npair*(l+1))) 
                        read (19) jprallo((npair*l+1):(npair*(l+1))) 
                        read (19) xallo((ixtemp*l+1):(ixtemp*(l+1)))               
                    enddo                          
                endif  
                do j=1,n1   
                    if (myrank().eq.0) then
                        do k=0,np-1                
                            l=np*(j-1)+k ! convert the current core # from new to old. 
                            if (k /= 0) then
                                call send(iplallo((npair*l+1):(npair*(l+1))),k,1) 
                                call send(jplallo((npair*l+1):(npair*(l+1))),k,2)
                                call send(iprallo((npair*l+1):(npair*(l+1))),k,3) 
                                call send(jprallo((npair*l+1):(npair*(l+1))),k,4) 
                                call send(xallo((ixtemp*l+1):(ixtemp*(l+1))),k,5)      
                            else ! send to process 0.
                                iplo=iplallo((npair*l+1):(npair*(l+1)))
                                jplo=jplallo((npair*l+1):(npair*(l+1)))
                                ipro=iprallo((npair*l+1):(npair*(l+1)))
                                jpro=jprallo((npair*l+1):(npair*(l+1)))
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
                                    call send(iplallo((npair*l+1):(npair*(l+1))),k,1) 
                                    call send(jplallo((npair*l+1):(npair*(l+1))),k,2)
                                    call send(iprallo((npair*l+1):(npair*(l+1))),k,3) 
                                    call send(jprallo((npair*l+1):(npair*(l+1))),k,4) 
                                    call send(xallo((ixtemp*l+1):(ixtemp*(l+1))),k,5) 
                                else
                                    iplo=iplallo((npair*l+1):(npair*(l+1)))
                                    jplo=jplallo((npair*l+1):(npair*(l+1)))
                                    ipro=iprallo((npair*l+1):(npair*(l+1)))
                                    jpro=jprallo((npair*l+1):(npair*(l+1)))
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
                            iplall((npair*l+1):(npair*(l+1)),n6c)=iplallo((npair*l+1):(npair*(l+1)))
                            jplall((npair*l+1):(npair*(l+1)),n6c)=jplallo((npair*l+1):(npair*(l+1)))
                            iprall((npair*l+1):(npair*(l+1)),n6c)=iprallo((npair*l+1):(npair*(l+1)))
                            jprall((npair*l+1):(npair*(l+1)),n6c)=jprallo((npair*l+1):(npair*(l+1)))
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
                                       call send( iplall((npair*l+1):(npair*(l+1)),n61),k,1   )
                                       call send( jplall((npair*l+1):(npair*(l+1)),n61),k,2   )
                                       call send( iprall((npair*l+1):(npair*(l+1)),n61),k,3   )
                                       call send( jprall((npair*l+1):(npair*(l+1)),n61),k,4   )
                                       call send( xall((ixtemp*l+1):(ixtemp*(l+1)),n61),k,5 )
                                   else
                                       iplo=iplall((npair*l+1):(npair*(l+1)),n61)
                                       jplo=jplall((npair*l+1):(npair*(l+1)),n61)
                                       ipro=iprall((npair*l+1):(npair*(l+1)),n61)
                                       jpro=jprall((npair*l+1):(npair*(l+1)),n61)
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
                allocate(iplall((npair*m4+1):npair*npo,n5),jplall((npair*m4+1):npair*npo,n5) &
                        ,iprall((npair*m4+1):npair*npo,n5),jprall((npair*m4+1):npair*npo,n5))         
            endif        
            n4c=0 ! count n4
            n5c=0   
            do i=istart,ipath
                        
                if ((np <= npo).and.((i-1) >= 1)) then    
                    ipathcountold=i-1
                    if (myrank().eq.0) call writecheckpoint(irep,idummy)                       
                endif

                if (myrank()==0) then
                    time0=mpi_wtime()
                    read (19) ipathnow
                    do l=0,npo-1    
                        read (19) lnow
                        read (19) iplallo((npair*l+1):(npair*(l+1))) 
                        read (19) jplallo((npair*l+1):(npair*(l+1))) 
                        read (19) iprallo((npair*l+1):(npair*(l+1))) 
                        read (19) jprallo((npair*l+1):(npair*(l+1))) 
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
                            call send(iplallo((npair*l+1):(npair*(l+1))),k,1) 
                            call send(jplallo((npair*l+1):(npair*(l+1))),k,2)
                            call send(iprallo((npair*l+1):(npair*(l+1))),k,3) 
                            call send(jprallo((npair*l+1):(npair*(l+1))),k,4) 
                            call send(xallo((ixtemp*l+1):(ixtemp*(l+1))),k,5)
                        else ! send to rank0 itself.
                            iplo=iplallo((npair*l+1):(npair*(l+1)))
                            jplo=jplallo((npair*l+1):(npair*(l+1)))
                            ipro=iprallo((npair*l+1):(npair*(l+1)))
                            jpro=jprallo((npair*l+1):(npair*(l+1)))
                            x0tot1d=xallo((ixtemp*l+1):(ixtemp*(l+1)))                        
                        endif 
                    enddo                                   
                    if ( loadextra ) then
                        n5c=n5c+1
                        do l=m4,npo-1
                            iplall((npair*l+1):(npair*(l+1)),n5c)=iplallo((npair*l+1):(npair*(l+1)))
                            jplall((npair*l+1):(npair*(l+1)),n5c)=jplallo((npair*l+1):(npair*(l+1)))
                            iprall((npair*l+1):(npair*(l+1)),n5c)=iprallo((npair*l+1):(npair*(l+1)))
                            jprall((npair*l+1):(npair*(l+1)),n5c)=jprallo((npair*l+1):(npair*(l+1)))
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
                                   call send( iplall((npair*l+1):(npair*(l+1)),n51),k,1   )
                                   call send( jplall((npair*l+1):(npair*(l+1)),n51),k,2   )
                                   call send( iprall((npair*l+1):(npair*(l+1)),n51),k,3   )
                                   call send( jprall((npair*l+1):(npair*(l+1)),n51),k,4   )
                                   call send( xall((ixtemp*l+1):(ixtemp*(l+1)),n51),k,5 )
                               else
                                   iplo=iplall((npair*l+1):(npair*(l+1)),n51)
                                   jplo=jplall((npair*l+1):(npair*(l+1)),n51)
                                   ipro=iprall((npair*l+1):(npair*(l+1)),n51)
                                   jpro=jprall((npair*l+1):(npair*(l+1)),n51)
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
            
    case (7)
        
        ! aother mode, change the way the pimc output the path.unf,
        ! not decided and checked yet.
        ! or may besupermode, do not need to consider ram resource problem, big ram is ready. bug need to fix
        ! can use MPI IO method.
               
          if (myrank().eq.0) then
              
              write(6,*) 'irep=7 case has not been decided yet, abort!'
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
    real(kind=r8) :: rm,rmsq,rmn,rmnsq,rmp,rmpsq
    real(kind=r8) :: v(0:nchorizo),vnn(0:nchorizo),vem(0:nchorizo)  &
	                ,x0totin(3,npart,0:nchorizo)
    real(kind=r8), dimension(0:nrhobin) :: rhodistout
    integer(kind=i4) :: i,rankl,rankr,icheck
    real(kind=r8) :: vn,vd,vnke,vnpev6,vnpeem,vnpels,rf,signrf,vnum,vdenom
    real(kind=r8) :: psi20,psi2n
    integer(kind=i4) :: iploin(npair),jploin(npair),iproin(npair),jproin(npair)
    real(kind=r8) :: k(3),kabs,kx,ky,kz,phi,theta
    integer(kind=i4) :: myunit
    complex(kind=r8) :: cwtl2rm(0:nspin-1,nisospin),cwtr2lm(0:nspin-1,nisospin)
! init   
    call chorizoallin(x0totin)
    call inputordlr(iploin,jploin,iproin,jproin)   
    call updatecwtchain
    
! calculate    
    call hpsi(ristra(0),ristra(nchorizo),vn,vd,vnke,vnpev6,vnpeem,vnpels,psi20,psi2n,rf,icheck)
    vnum=vn
    vdenom=vd
    ! update the v beads   
    do i=0,nchorizo
        call vtot(i,v(i),vnn(i),vem(i))
    enddo
    v=v/abs(rf) ! in the end need to be devided by average of vdenom and consider its error.
    vnn=vnn/abs(rf)
    vem=vem/abs(rf)
    call addvalristra(1,v) 
    call addvalristra(2,vnn) 
    call addvalristra(3,vem)
    
! Spherical coordinate system. phi is angle with z, theta and k xy projection angle with x.
    kabs=350_r8
    phi=0
    theta=0
    kx=kabs*sin(phi)*cos(theta) 
    ky=kabs*sin(phi)*sin(theta)
    kz=kabs*cos(phi)
    
    k=(/kx,ky,kz/)
    call euclideanresponsecalc(k,ebeads,ewbeads,e0beads,ew0beads,e0beadscrs,ew0beadscrs)  
    e0beads=e0beads/abs(rf)
    ew0beads=ew0beads/abs(rf)  
    e0beadscrs=e0beadscrs/abs(rf)
    ew0beadscrs=ew0beadscrs/abs(rf)
    
    call addvalristrae0resp(1,e0beads(1,:))
    call addvalristrae0resp(2,e0beads(2,:))
    call addvalristrae0resp(3,e0beads(3,:))
    call addvalristrae0resp(4,e0beads(4,:)) 
    
    call addvalristrae0resp(5,ew0beads(1,:))
    call addvalristrae0resp(6,ew0beads(2,:))
    call addvalristrae0resp(7,ew0beads(3,:))
    call addvalristrae0resp(8,ew0beads(4,:))  
    
    call addvalristrae0resp(9,e0beadscrs(1,:)) 
    call addvalristrae0resp(10,e0beadscrs(2,:)) 
    call addvalristrae0resp(11,e0beadscrs(3,:))
    call addvalristrae0resp(12,e0beadscrs(4,:)) 
    call addvalristrae0resp(13,e0beadscrs(5,:)) 
    call addvalristrae0resp(14,e0beadscrs(6,:)) 
    call addvalristrae0resp(15,e0beadscrs(7,:))
    call addvalristrae0resp(16,e0beadscrs(8,:))
    call addvalristrae0resp(17,e0beadscrs(9,:)) 
    
    call addvalristrae0resp(18,ew0beadscrs(1,:)) 
    call addvalristrae0resp(19,ew0beadscrs(2,:)) 
    call addvalristrae0resp(20,ew0beadscrs(3,:))
    call addvalristrae0resp(21,ew0beadscrs(4,:)) 
    call addvalristrae0resp(22,ew0beadscrs(5,:)) 
    call addvalristrae0resp(23,ew0beadscrs(6,:)) 
    call addvalristrae0resp(24,ew0beadscrs(7,:))
    call addvalristrae0resp(25,ew0beadscrs(8,:))
    call addvalristrae0resp(26,ew0beadscrs(9,:)) 
    
    !call addvalristrae0resp(27,e0beadscrs(1,:)+e0beadscrs(5,:)+e0beadscrs(9,:)) ! xx+yy+zz
    !call addvalristrae0resp(28,ew0beadscrs(1,:)+ew0beadscrs(5,:)+ew0beadscrs(9,:)) 
    
    call addvalristrae0resp(27,e0beads(5,:))
    call addvalristrae0resp(28,ew0beads(5,:))
    
    
    !k=(/kabs,zero,zero/)
    !call euclideanresponsecalc(k,ebeads,ewbeads,e0beads,ew0beads,e0beadscrs,ew0beadscrs)  
    !e0beads=e0beads/abs(rf)
    !ew0beads=ew0beads/abs(rf) 
    !
    !call addvalristrae0resp(9,e0beads(1,:))
    !call addvalristrae0resp(10,e0beads(2,:))
    !call addvalristrae0resp(11,e0beads(3,:))
    !call addvalristrae0resp(12,e0beads(4,:)) 
    !call addvalristrae0resp(13,ew0beads(1,:))
    !call addvalristrae0resp(14,ew0beads(2,:))
    !call addvalristrae0resp(15,ew0beads(3,:))
    !call addvalristrae0resp(16,ew0beads(4,:))  
    !
    !k=(/zero,kabs,zero/)
    !call euclideanresponsecalc(k,ebeads,ewbeads,e0beads,ew0beads,e0beadscrs,ew0beadscrs)  
    !e0beads=e0beads/abs(rf)
    !ew0beads=ew0beads/abs(rf) 
    !
    !call addvalristrae0resp(17,e0beads(1,:))
    !call addvalristrae0resp(18,e0beads(2,:))
    !call addvalristrae0resp(19,e0beads(3,:))
    !call addvalristrae0resp(20,e0beads(4,:)) 
    !call addvalristrae0resp(21,ew0beads(1,:))
    !call addvalristrae0resp(22,ew0beads(2,:))
    !call addvalristrae0resp(23,ew0beads(3,:))
    !call addvalristrae0resp(24,ew0beads(4,:))     
    
! response function end    
    
    call addval(3,vnum,1.0_r8)
    call addval(4,vdenom,1.0_r8)
	call addval(5,rf,1.0_r8)    
    
    call getcwtlrm(nchorizomid,cwtl2rm,cwtr2lm)
	call samplerhodistrbt(ristra(nchorizomid)%x,nchorizomid,cwtl2rm,cwtr2lm,rf,rhodistout,rm,rmsq,rmn,rmnsq,rmp,rmpsq)
    call addval(6,rmsq,1.0_r8)
    call addval(7,rm,1.0_r8)
    call addval(8,rhodistout(0),1.0_r8) 
    
    call addval(13,rmnsq,1.0_r8)
    call addval(14,rmn,1.0_r8)      
    call addval(15,rmpsq,1.0_r8)
    call addval(16,rmp,1.0_r8)
    
    call addval(17,vnke,1.0_r8)
    call addval(18,vnpev6+vnpeem+vnpels,1.0_r8)
    call addval(19,vnpev6,1.0_r8)
    call addval(20,vnpeem,1.0_r8)
    call addval(21,vnpels,1.0_r8)
  
    ! write out all the num and denom of all the cores
    
    if (rankl/=0) then  ! rankl must be zero   
        if (myrank().eq.0) then
            write(6,*) ' rankl /= 0! computecore error! Abort! '
            call abort
        endif 
    endif  
      
    !if ((myrank()>=rankl).and.(myrank()<=rankr)) then ! rankl must be zero   
    !    if (myrank()/=0) then
    !        call send(vnum,0,1)
    !        call send(vdenom,0,2)
    !        call send(vnke,0,3)
    !        call send(vnpev6,0,4)
    !        call send(vnpeem,0,5)
    !        call send(rm,0,6)
    !        call send(rhodistout(0),0,7)  
    !        call send(rmsq,0,8)
    !        call send(rmn,0,9)
    !        call send(rmnsq,0,10)
    !        call send(rmp,0,11)
    !        call send(rmpsq,0,12)
    !    else ! myrank().eq.0 which should be rankl
    !        open (newunit=myunit,FILE=file_values,FORM='FORMATTED',position='append')    
    !            write(myunit,'(10(G15.7,2x))') vnum,vdenom,vnke,vnpev6,vnpeem,rm,rmsq,rmn,rmnsq,rmp,rmpsq,rhodistout(0)     !myunit means write to 'values.txt'    
    !            do i=1,rankr
    !                call recv(vnum,i,1)
    !                call recv(vdenom,i,2)
    !                call recv(vnke,i,3)
    !                call recv(vnpev6,i,4)
    !                call recv(vnpeem,i,5)
    !                call recv(rm,i,6)
    !                call recv(rhodistout(0),i,7)  
    !                call recv(rmsq,i,8)
    !                call recv(rmn,i,9)
    !                call recv(rmnsq,i,10)
    !                call recv(rmp,i,11)
    !                call recv(rmpsq,i,12)
    !                write(myunit,'(10(G15.7,2x))') vnum,vdenom,vnke,vnpev6,vnpeem,rm,rmsq,rmn,rmnsq,rmp,rmpsq,rhodistout(0)    
    !            enddo   
    !        close(myunit)
    !    endif    
    !endif

    return
    end subroutine computecore   
    
    subroutine computeupdate(iin) ! require all the cores.
    use estimator
    use estimatorristra
    use wavefunction
    use random
    use brut
    use mympi
    integer(kind=i4) :: i,j,k,iin,idummy1,idummy2,myunit
    real(kind=r8) :: val1,val2,valnn1,valem1
    real(kind=r8), allocatable :: valnow(:),valaverage(:),valerr(:) &
                                 ,valnnnow(:),valnnaverage(:),valnnerr(:) &
                                 ,valemnow(:),valemaverage(:),valemerr(:)
    real(kind=r8) :: err1,err2,error,errnn1,errem1,errornn,errorem
    character(len=30) :: vpropchorizo,vpropchorizotmp,ce0resptmp
    real(kind=r8) :: vale0resp1(nrespest),erre0resp1(nrespest),errore0resp(nrespest) &
    ,vale0respnow(nrespest,0:ne0respbeads),vale0respaverage(nrespest,0:ne0respbeads),vale0resperr(nrespest,0:ne0respbeads)
    character(len=120) :: file_response,file_responsedetails,file_id

    call update   !collect block averages
    call updateristra
    call updateristrae0resp
    call updaterhodistrbt(nchorizomid) ! choose the middle bead
      
    if (myrank().eq.0) then   
        write (6,*) 'path #: ', iin 
	    write (12,*) 'path #: ', iin  
        answer=resstring()
        call writeanswer(answer) 
    endif    
     
    ! show rho distribution on the fly
    call writerhodist(nchorizomid)
     
    ! write out a file with the Vtot values at different chorizos
    call writevtotbeads(dt)     
    
    ! write response function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (myrank().eq.0) then  
        
        do i=1,nrespest
            call resultristrae0resp(i,vale0respnow(i,:),vale0respaverage(i,:),vale0resperr(i,:),ce0resptmp)
        enddo
        
        !i=3
        !call resultristrae0resp(i,vale0respnow(i,:),vale0respaverage(i,:),vale0resperr(i,:),ce0resptmp)

        
        write(file_id,'(i5)') irep
        file_response='E0response_CalcMode' // trim(adjustl(file_id)) // '.txt'          
        file_responsedetails='E0responsedetails_CalcMode' // trim(adjustl(file_id)) // '.txt'  
        
        open(newunit=myunit,form='formatted',file=trim(file_response))
        call resultest(4,val2,err2) ! #4 is the denominator 
        do j=0,ne0respbeads 
            vale0resp1(:)=vale0respaverage(:,j)
            erre0resp1(:)=vale0resperr(:,j)  
            errore0resp(:)=sqrt((1/val2*erre0resp1(:))**2+(vale0resp1(:)/(val2**2)*err2)**2)
                ! val1/val2 and error are the Vtot and error we want.    
            !write (myunit,'(i10,f10.5,48g15.6)') j*irlstep2,j*irlstep2*dt &
            !                                ,vale0resp1(3)/val2,errore0resp(3)
            write (myunit,'(i10,f10.5,60g15.6)') j*irlstep2,j*irlstep2*dt &
                                            ,vale0resp1(5)/val2,vale0resp1(6)/val2,vale0resp1(7)/val2,vale0resp1(8)/val2 &
                                            ,vale0resp1(28)/val2 &
                                            ,errore0resp(5),errore0resp(6),errore0resp(7),errore0resp(8),errore0resp(28)
        enddo
        close(myunit) 
        
        open(newunit=myunit,form='formatted',file=trim(file_responsedetails))
        call resultest(4,val2,err2) ! #4 is the denominator 
        do j=0,ne0respbeads 
            vale0resp1(:)=vale0respaverage(:,j)
            erre0resp1(:)=vale0resperr(:,j)  
            errore0resp(:)=sqrt((1/val2*erre0resp1(:))**2+(vale0resp1(:)/(val2**2)*err2)**2)
            write (myunit,'(i10,f10.5,60g15.6)') j*irlstep2,j*irlstep2*dt &
                                            ,vale0resp1(1:nrespest)/val2,errore0resp(1:nrespest)
        enddo
        close(myunit)
    endif     
    
    !ipathcountold=iin
    !call writecheckpoint(irep,idummy1)  ! need to change its position.

    return
    end subroutine computeupdate    
    
    subroutine computeinit(answerin,rhobinsizein) ! call it before call compute
    use mympi
    use estimatorristra
    character(len=120) :: answerin(:)
    real(kind=r8) :: rhobinsizein
    integer(kind=i4) :: ntmp
    logical :: file_exists
    allocate(answer(size(answerin)))
    rhobinsize=rhobinsizein 
    ipathcountold=1 ! set initial value for ipathcountold
    if (myrank().eq.0) then 
        INQUIRE(FILE=file_pathnumbercount, EXIST=file_exists) 
        if (file_exists) then
            open(unit=19,form='formatted',file=file_pathnumbercount,position='rewind')
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
        else
            write (6,*) 'file_pathnumbercount does not exist, abort!',file_pathnumbercount
            call abort      
        endif
    endif
    call bcast(npo)
    call bcast(ipath)  
    call bcast(ixtempo)     
    call barrier
! response function init.    
    ntmp=nint(0.05/dt)
    if ((ntmp+1)>nchorizo) then
        write (6,*) 'nint(0.04/dt)+1 > nchorizo, stop: ',ntmp,ntmp+1,nchorizo
        call abort
    else   
        if (mod(ntmp,2).eq.0) then
            irlintervalmax=ntmp
        else
            irlintervalmax=ntmp+1  ! in the middle. make it even number.
        endif         
    endif
    irlintervalmax1=24 ! along beads.
    irlstep1=3  
    !n1=irlintervalmax1/irlstep1+1
    irlstep2=2 ! along the beads from the middle
    nrespest=28
    call responsecalcinit(nrespest,irlintervalmax,irlintervalmax1,irlstep1,irlstep2)
    return
    end subroutine computeinit 

    subroutine responsecalcinit(nrespestin,irlintervalmaxin,irlintervalmax1in,irlstep1in,irlstep2in)
! read in path, calculate all the results. and write to a file
! rank1 must be zero   
    use estimator
    use estimatorristra
    use wavefunction
    use random
    use brut
    use mympi
    integer(kind=i4) :: irlintervalmaxin,irlintervalmax1in,irlstep1in,irlstep2in,nrespestin
    irlintervalmax=irlintervalmaxin
    irlintervalmax1=irlintervalmax1in
    irlstep1=irlstep1in
    irlstep2=irlstep2in
    ne0respbeads=irlintervalmax/irlstep2
    nrespest=nrespestin
    call setestnumristraresp(nrespestin,irlintervalmaxin,irlintervalmax1in,irlstep1in,irlstep2in)
    allocate(ebeads(5,0:nchorizo,0:irlintervalmax1),ewbeads(5,0:nchorizo,0:irlintervalmax1) &
             ,e0beads(5,0:ne0respbeads),ew0beads(5,0:ne0respbeads) &
             ,e0beadscrs(9,0:ne0respbeads),ew0beadscrs(9,0:ne0respbeads))
    call zerestristrae0resp
    return
    end subroutine responsecalcinit  
    
    subroutine responseinfo(nrespestio,ne0respbeadsio) ! before call it, make sure variables has already been assigned values.
    integer(kind=i4) :: nrespestio,ne0respbeadsio
    nrespestio=nrespest
    ne0respbeadsio=ne0respbeads
    return
    end subroutine responseinfo
    
    subroutine euclideanresponsecalc(k,ebeads,ewbeads,e0beadsout,ew0beadsout,e0beadscrsout,ew0beadscrsout)  
    use mympi
    integer(kind=i4) :: i,il,ir,irlintervalnow
    real(kind=r8) :: k(3),e(5),ew(5),ecrs(9),ewcrs(9) &
                    ,ebeads(:,:,:),ewbeads(:,:,:) &
                    ,e0beadsout(5,0:ne0respbeads),ew0beadsout(5,0:ne0respbeads) &
                    ,e0beadscrsout(9,0:ne0respbeads),ew0beadscrsout(9,0:ne0respbeads)
    ! k should be vector.
    
! along the beads 
    if (nchorizo>=irlintervalmax) then  
    !    do irlintervalnow=0,irlintervalmax1,irlstep1
    !        do il=0,nchorizo-irlintervalnow
    !            ir=il+irlintervalnow
    !            call euclideanresponse(k,il,ir,e,ew)   
    !            ebeads(:,il,irlintervalnow)=e(:)
    !            ewbeads(:,il,irlintervalnow)=ew(:)   
    !        enddo
    !    enddo
      
! in the middle
        do irlintervalnow=0,irlintervalmax,irlstep2  
            il=nchorizo/2-irlintervalnow/2
            ir=nchorizo/2+irlintervalnow/2
            
            !write(6,*) 'myrank=', myrank()
            !write(6,*) 'il,ir check:', il,ir,nchorizo/2,irlintervalnow/2,irlintervalnow,irlstep2        
            
            select case (irep)
            case (5)    
                call euclideanresponseangleave(k,il,ir,e,ew) ! angle averaged mode, slow.  
            case (6)     
                call euclideanresponse(k,il,ir,e,ew)  ! fast mode
            case default
                call euclideanresponseangleave(k,il,ir,e,ew)   
            end select

            e0beadsout(:,irlintervalnow/2)=e(:)
            ew0beadsout(:,irlintervalnow/2)=ew(:) 

            !call euclideanresponsecrs(k,il,ir,ecrs,ewcrs)
            !e0beadscrsout(:,irlintervalnow/2)=ecrs(:)
            !ew0beadscrsout(:,irlintervalnow/2)=ewcrs(:)       
                 
        enddo
             
    else
        if (myrank().eq.0) then
            write(6,*) 'nchorizo < irlintervalmax , response cannot be calculated!'
            call abort
        endif           
    endif    
    
    return
    end subroutine euclideanresponsecalc
    
    subroutine euclideanresponse(k,il,ir,e,ew) ! scalar, fast mode. no angle average.
    integer(kind=i4) :: il,ir
    real(kind=r8) :: k(3),e(5),ew(5),edummy(9),ewdummy(9)
    
    !call resprhonchk(k,il,ir,e(1),ew(1))
    call resprhon(k,il,ir,e(1),ew(1))
    !call resprhopchk(k,il,ir,e(2),ew(2))
    call resprhop(k,il,ir,e(2),ew(2))
    !call resprhotchk(k,il,ir,e(3),ew(3))
    call resprhot(k,il,ir,e(3),ew(3))
    !call resprhostlchk(k,il,ir,e(4),ew(4))
    call resprhostl(k,il,ir,e(4),ew(4))
    !call resprhosttchk(k,il,ir,e(5),ew(5)) ! added
    call resprhostt(k,il,ir,e(5),ew(5),edummy,ewdummy)
    return
    end subroutine euclideanresponse
    
    subroutine euclideanresponsecrs(k,il,ir,e,ew) ! crs means with cross product, vector.
    integer(kind=i4) :: il,ir
    real(kind=r8) :: k(3),e(9),ew(9),eval,ewval
    call resprhostt(k,il,ir,eval,ewval,e,ew)
    return
    end subroutine euclideanresponsecrs 
    
    
    subroutine euclideanresponseangleave(k,il,ir,e,ew) ! angle averaged 5 response functions. Carlson 1994.
    use brut
    use v6stepcalc
    use math
    integer(kind=i4) :: lspin,ipi,il,ir,i,j,l,spni,jpi,i1,icase,icasestep 
    real(kind=r8) :: e(5),ew(5),k(3),x(3,npart),xil(3,npart),xir(3,npart),m,tau,omega,kmod,dr(3),r,rx,ry,rz,c1a,c1(9),kr,r1(3),sj0,sj2 &
                    ,drij(3,npart,npart),rij(npart,npart),rhonlr,c1tmp(9)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtnew1(0:nspin-1,nisospin) &
                       ,rhonl,rhonr,cpathcoreval(5) &
                       ,cwtr2lm(0:nspin-1,nisospin,0:nchorizo),cwtr2ltmp(0:nspin-1,nisospin,0:nchorizo) &
                       ,cpathcorevaltot(5),cwttmp(0:nspin-1,nisospin) 
    if (il<=ir) then
        m=(mp+mn)/2
        rhonl=0
        rhonr=0
        rhonlr=0
        tau=(ir-il)*dt
        xil=ristra(il)%x
        xir=ristra(ir)%x
        call xconvert(xil)
        call xconvert(xir)
        kmod=sqrt(dot_product(k,k))
        omega=dot_product(k,k)/(2*m) ! omega_el. 
        cpathcorevaltot(:)=czero
! rhon             
        do i=1,npart
            do j=1,npart               
                dr(:)=xil(:,i)-xir(:,j)
                drij(:,i,j)=dr(:)
		        r=sqrt(dot_product(dr,dr))  
                rij(i,j)=r
                !rhonl=rhonl+exp(-ci*dot_product(k(:),xil(:,i))/hc)
                !rhonr=rhonr+exp(ci*dot_product(k(:),xir(:,i))/hc) 
                if ((kmod*r/hc)>tiny) then
                    c1a=sin(kmod*r/hc)/(kmod*r/hc)
                else
                    c1a=1
                endif                
                rhonlr=rhonlr+c1a                        
! rhop                
                cwtr2ltmp=cwtr2l 
                if (il<ir) then     
                    if (ir.eq.nchorizo) then
                        do l=1,nisospin
                            lspin=liso(l)  
                            jpi=and(1,shiftr(lspin,j-1)) 
                            if (jpi.eq.1) then
                              cwtnew(:,l)=cwtend(:,l)
                            else
                              cwtnew(:,l)=0    
                            endif                                 
                        enddo 
                        call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                    else                     
                        call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                        do l=1,nisospin
                            lspin=liso(l)  
                            jpi=and(1,shiftr(lspin,j-1)) 
                            if (jpi.eq.1) then
                              cwtnew(:,l)=cwtr2lm(:,l,ir)
                            else
                              cwtnew(:,l)=0    
                            endif   
                        enddo  
                        call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                    endif  
                    if ((ir-il)>1) then
                        do i1=ir-1,il+1,-1        
                           x=ristra(i1)%x
                           call vlsprop(spinorbit,1,x,ristra(i1+1)%x,cwtr2ltmp(:,:,i1+1),cwttmp(:,:))
                           call v6proplr(x,x,-1,cwttmp(:,:),cwtr2ltmp(:,:,i1))      
                        enddo
                    endif
                    call vlsprop(spinorbit,1,xil,ristra(il+1)%x,cwtr2ltmp(:,:,il+1),cwttmp(:,:))
                    call v6propl(xil,-1,cwttmp(:,:),cwtr2lm(:,:,il))               
                    do l=1,nisospin
                        lspin=liso(l)  
                        ipi=and(1,shiftr(lspin,i-1)) 
                        if (ipi.eq.1) then
                            cwtnew(:,l)=cwtr2lm(:,l,il)
                        else
                            cwtnew(:,l)=0    
                        endif   
                    enddo                 
                    if (il/=0) then                
                        call v6propr(xil,-1,cwtnew(:,:),cwtr2ltmp(:,:,il))                     
                        call vlsprop(spinorbit,1,ristra(il-1)%x,xil,cwtr2ltmp(:,:,il),cwttmp(:,:))
                        cpathcoreval(2)=sum(conjg(cwtl2r(:,:,il-1))*cwttmp(:,:))
                    else ! il==0
                        cpathcoreval(2)=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))    
                    endif         
                else ! il==ir
                    if (ir.eq.nchorizo) then                   
                        do l=1,nisospin
                            lspin=liso(l)  
                            ipi=and(1,shiftr(lspin,i-1)) 
                            jpi=and(1,shiftr(lspin,j-1)) 
                            if ((ipi.eq.1).and.(jpi.eq.1)) then
                              cwtnew(:,l)=cwtend(:,l)
                            else
                              cwtnew(:,l)=0    
                            endif   
                        enddo                         
                        call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir)) 
                        call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                        cpathcoreval(2)=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))  
                    else     
                        if (ir.eq.0) then                     
                            do l=1,nisospin
                                lspin=liso(l)  
                                ipi=and(1,shiftr(lspin,i-1)) 
                                jpi=and(1,shiftr(lspin,j-1)) 
                                if ((ipi.eq.1).and.(jpi.eq.1)) then
                                  cwtnew(:,l)=cwtr2l(:,l,ir)   
                                else
                                  cwtnew(:,l)=0    
                                endif   
                            enddo                         
                            cpathcoreval(2)=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))     
                        else ! il=ir at somewhere middle
                            call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                            do l=1,nisospin
                                lspin=liso(l)  
                                ipi=and(1,shiftr(lspin,i-1)) 
                                jpi=and(1,shiftr(lspin,j-1)) 
                                if ((ipi.eq.1).and.(jpi.eq.1)) then
                                  cwtnew(:,l)=cwtr2lm(:,l,ir)
                                else
                                  cwtnew(:,l)=0    
                                endif   
                            enddo                          
                            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                            call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                            cpathcoreval(2)=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))   
                        endif   
                    endif   
                endif            
                cpathcorevaltot(2)=cpathcorevaltot(2)+cpathcoreval(2)*c1a                               
! rhotau  isoscalar       
                cwtr2ltmp=cwtr2l    
                if (il<ir) then     
                    if (ir.eq.nchorizo) then
                        do l=1,nisospin
                            lspin=liso(l)  
                            jpi=2*and(1,shiftr(lspin,j-1))-1  
                            cwtnew(:,l)=cwtend(:,l)*jpi
                        enddo                               
                        call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                    else
                        call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                        do l=1,nisospin
                            lspin=liso(l)  
                            jpi=2*and(1,shiftr(lspin,j-1))-1  
                            cwtnew(:,l)=cwtr2lm(:,l,ir)*jpi
                        enddo                  
                        call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                    endif  
                    if ((ir-il)>1) then
                        do i1=ir-1,il+1,-1        
                           x=ristra(i1)%x
                           call vlsprop(spinorbit,1,x,ristra(i1+1)%x,cwtr2ltmp(:,:,i1+1),cwttmp(:,:))
                           call v6proplr(x,x,-1,cwttmp(:,:),cwtr2ltmp(:,:,i1))      
                        enddo
                    endif
                    call vlsprop(spinorbit,1,xil,ristra(il+1)%x,cwtr2ltmp(:,:,il+1),cwttmp(:,:))
                    call v6propl(xil,-1,cwttmp(:,:),cwtr2lm(:,:,il))          
                    do l=1,nisospin
                        lspin=liso(l)  
                        ipi=2*and(1,shiftr(lspin,i-1))-1  
                        cwtnew(:,l)=cwtr2lm(:,l,il)*ipi
                    enddo             
                    if (il/=0) then
                        call v6propr(xil,-1,cwtnew(:,:),cwtr2ltmp(:,:,il)) 
                        call vlsprop(spinorbit,1,ristra(il-1)%x,xil,cwtr2ltmp(:,:,il),cwttmp(:,:))
                        cpathcoreval(3)=sum(conjg(cwtl2r(:,:,il-1))*cwttmp(:,:))
                    else ! il==0
                        cpathcoreval(3)=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))    
                    endif         
                else ! il==ir
                    if (ir.eq.nchorizo) then              
                            do l=1,nisospin
                                lspin=liso(l)  
                                jpi=2*and(1,shiftr(lspin,j-1))-1  
                                ipi=2*and(1,shiftr(lspin,i-1))-1
                                cwtnew(:,l)=cwtend(:,l)*ipi*jpi
                            enddo                 
                        call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))    
                        call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                        cpathcoreval(3)=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))  
                    else     
                        if (ir.eq.0) then   
                            do l=1,nisospin
                                lspin=liso(l)  
                                jpi=2*and(1,shiftr(lspin,j-1))-1  
                                ipi=2*and(1,shiftr(lspin,i-1))-1
                                cwtnew(:,l)=cwtr2l(:,l,ir)*ipi*jpi
                            enddo                     
                            cpathcoreval(3)=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))     
                        else ! il=ir at somewhere middle
                            call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                            do l=1,nisospin
                                lspin=liso(l)  
                                jpi=2*and(1,shiftr(lspin,j-1))-1  
                                ipi=2*and(1,shiftr(lspin,i-1))-1
                                cwtnew(:,l)=cwtr2lm(:,l,ir)*ipi*jpi
                            enddo                         
                            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                            call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                            cpathcoreval(3)=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))   
                        endif   
                    endif   
                endif            
                cpathcorevaltot(3)=cpathcorevaltot(3)+cpathcoreval(3)*c1a                
! rhostl            
                dr(:)=-drij(:,i,j) ! rji
		        r=rij(i,j)  
                kr=kmod*r/hc    
                if (kr > tiny) then 
                    r1(:)=dr(:)/r
                    rx=r1(1)
                    ry=r1(2)
                    rz=r1(3)
                    call sphbesselj(0,kr,sj0) 
                    call sphbesselj(2,kr,sj2)  
                    c1(1)=(sj0+sj2*(rz**2-2*rx**2+ry**2))/3.0_r8
                    c1(2)=-sj2*rx*ry
                    c1(3)=-sj2*rx*rz
                    c1(4)=c1(2)
                    c1(5)=(sj0+sj2*(rz**2-2*ry**2+rx**2))/3.0_r8
                    c1(6)=-sj2*ry*rz
                    c1(7)=c1(3)
                    c1(8)=c1(6)
                    c1(9)=(sj0+sj2*(rx**2-2*rz**2+ry**2))/3.0_r8                  
                    c1tmp(:)=c1(:)                   
                    icasestep=1
                else
                    c1(:)=0
                    c1(1)=1.0_r8/3
                    c1(5)=c1(1)
                    c1(9)=c1(1)  
                    c1tmp(:)=c1(:)                    
                    icasestep=4
                endif
                do icase=1,9,icasestep
                    cwtr2ltmp=cwtr2l 
                    if (il<ir) then     
                        if (ir.eq.nchorizo) then                  
                            call resprhostlchkprop1(icase,j,cwtend,cwtnew) 
                            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                        else
                            call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))                 
                            call resprhostlchkprop1(icase,j,cwtr2lm(:,:,ir),cwtnew(:,:))
                            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                        endif  
                        if ((ir-il)>1) then         
                            do i1=ir-1,il+1,-1        
                                x=ristra(i1)%x
                                call vlsprop(spinorbit,1,x,ristra(i1+1)%x,cwtr2ltmp(:,:,i1+1),cwttmp(:,:))
                                call v6proplr(x,x,-1,cwttmp(:,:),cwtr2ltmp(:,:,i1))      
                            enddo      
                        endif
                        call vlsprop(spinorbit,1,xil,ristra(il+1)%x,cwtr2ltmp(:,:,il+1),cwttmp(:,:))
                        call v6propl(xil,-1,cwttmp(:,:),cwtr2lm(:,:,il))    
                        call resprhostlchkprop1pr(icase,i,cwtr2lm(:,:,il),cwtnew(:,:)) ! -k for l.
                        if (il/=0) then
                            call v6propr(xil,-1,cwtnew(:,:),cwtr2ltmp(:,:,il)) 
                            call vlsprop(spinorbit,1,ristra(il-1)%x,xil,cwtr2ltmp(:,:,il),cwttmp(:,:))
                            cpathcoreval(4)=sum(conjg(cwtl2r(:,:,il-1))*cwttmp(:,:))
                        else ! il==0
                            cpathcoreval(4)=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))    
                        endif         
                    else ! il==ir
                        if (ir.eq.nchorizo) then   
                            call resprhostlchkprop1(icase,j,cwtend(:,:),cwtnew1(:,:))
                            call resprhostlchkprop1pr(icase,i,cwtnew1(:,:),cwtnew(:,:))  
                            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))  
                            call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                            cpathcoreval(4)=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))  
                        else     
                            if (ir.eq.0) then
                                call resprhostlchkprop1(icase,j,cwtr2l(:,:,ir),cwtnew1(:,:))
                                call resprhostlchkprop1pr(icase,i,cwtnew1(:,:),cwtnew(:,:))  
                                cpathcoreval(4)=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))     
                            else ! il=ir at somewhere middle
                                call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                                call resprhostlchkprop1(icase,j,cwtr2lm(:,:,ir),cwtnew1(:,:))
                                call resprhostlchkprop1pr(icase,i,cwtnew1(:,:),cwtnew(:,:))    
                                call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                                call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                                cpathcoreval(4)=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))   
                            endif   
                        endif   
                    endif 
                    cpathcorevaltot(4)=cpathcorevaltot(4)+cpathcoreval(4)*c1(icase) 
                enddo                 
! rhostt
                ! rearrange
                c1(1)=c1tmp(5)+c1tmp(9)
                c1(2)=-c1tmp(2)                
                c1(3)=-c1tmp(3)
                c1(4)=-c1tmp(2)
                c1(5)=c1tmp(1)+c1tmp(9)
                c1(6)=-c1tmp(6)
                c1(7)=-c1tmp(7)
                c1(8)=-c1tmp(8)
                c1(9)=c1tmp(1)+c1tmp(5)         
                do icase=1,9,icasestep  
                    cwtr2ltmp=cwtr2l 
                    if (il<ir) then     
                        if (ir.eq.nchorizo) then                               
                            call resprhostlchkprop1(icase,j,cwtend,cwtnew) 
                            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                        else
                            call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))                 
                            call resprhostlchkprop1(icase,j,cwtr2lm(:,:,ir),cwtnew(:,:))
                            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                        endif  
                        if ((ir-il)>1) then         
                            do i1=ir-1,il+1,-1        
                                x=ristra(i1)%x
                                call vlsprop(spinorbit,1,x,ristra(i1+1)%x,cwtr2ltmp(:,:,i1+1),cwttmp(:,:))
                                call v6proplr(x,x,-1,cwttmp(:,:),cwtr2ltmp(:,:,i1))      
                            enddo      
                        endif                       
                        call vlsprop(spinorbit,1,xil,ristra(il+1)%x,cwtr2ltmp(:,:,il+1),cwttmp(:,:))
                        call v6propl(xil,-1,cwttmp(:,:),cwtr2lm(:,:,il))    
                        call resprhostlchkprop1pr(icase,i,cwtr2lm(:,:,il),cwtnew(:,:)) ! -k for l.
                        if (il/=0) then
                            call v6propr(xil,-1,cwtnew(:,:),cwtr2ltmp(:,:,il)) 
                            call vlsprop(spinorbit,1,ristra(il-1)%x,xil,cwtr2ltmp(:,:,il),cwttmp(:,:))
                            cpathcoreval(5)=sum(conjg(cwtl2r(:,:,il-1))*cwttmp(:,:))
                        else ! il==0
                            cpathcoreval(5)=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))    
                        endif         
                    else ! il==ir
                        if (ir.eq.nchorizo) then   
                            call resprhostlchkprop1(icase,j,cwtend(:,:),cwtnew1(:,:))
                            call resprhostlchkprop1pr(icase,i,cwtnew1(:,:),cwtnew(:,:))  
                            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))   
                            call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                            cpathcoreval(5)=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))  
                        else     
                            if (ir.eq.0) then
                                call resprhostlchkprop1(icase,j,cwtr2l(:,:,ir),cwtnew1(:,:))
                                call resprhostlchkprop1pr(icase,i,cwtnew1(:,:),cwtnew(:,:))  
                                cpathcoreval(5)=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))     
                            else ! il=ir at somewhere middle
                                call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                                call resprhostlchkprop1(icase,j,cwtr2lm(:,:,ir),cwtnew1(:,:))
                                call resprhostlchkprop1pr(icase,i,cwtnew1(:,:),cwtnew(:,:))    
                                call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                                call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                                cpathcoreval(5)=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))   
                            endif   
                        endif   
                    endif 
                    cpathcorevaltot(5)=cpathcorevaltot(5)+cpathcoreval(5)*c1(icase) 
                enddo              
            enddo
        enddo
        cpathcoreval(1)=sum(conjg(cwtbegin(:,:))*cwtr2l(:,:,0))
        e(1)=real(rhonlr*cpathcoreval(1))      
        e(2:5)=real(cpathcorevaltot(2:5))
        ew(1:5)=exp(tau*omega)*e(1:5) 
        !e(3)=real(cpathcorevaltot(3))
        !ew(3)=exp(tau*omega)*e(3)          
        !e(4)=real(cpathcorevaltot(4))
        !ew(4)=exp(tau*omega)*e(4)     
        !e(5)=real(cpathcorevaltot(5))
        !ew(5)=exp(tau*omega)*e(5) 
    else
        write(6,*) 'il > ir in resprhon, stop'
        stop             
    endif 
    return    
    end subroutine euclideanresponseangleave

    subroutine resprhon(k,il,ir,e,ew) ! Carlson 1994 PRC paper.
    use brut
    use v6stepcalc
    integer(kind=i4) :: il,ir,i ! il<ir
    real(kind=r8) :: e,ew,k(3),xil(3,npart),xir(3,npart),m,tau,omega
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin) &
                       ,rhonl,rhonr,cpathcoreval 
    if (il<=ir) then
        m=(mp+mn)/2
        rhonl=0
        rhonr=0
        tau=(ir-il)*dt
        xil=ristra(il)%x
        xir=ristra(ir)%x
        call xconvert(xil)
        call xconvert(xir)
        do i=1,npart
            rhonl=rhonl+exp(-ci*dot_product(k(:),xil(:,i))/hc)
            rhonr=rhonr+exp(ci*dot_product(k(:),xir(:,i))/hc) 
        enddo
        cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtr2l(:,:,0))
        e=real(rhonl*rhonr*cpathcoreval)
        omega=dot_product(k,k)/(2*m) ! omega_el. 
        ew=exp(tau*omega)*e  
        !write(6,*) 'k,il,ir,tau,omega,e,ew',k,il,ir,tau,omega,e,ew   
    else
        write(6,*) 'il > ir in resprhon, stop'
        stop             
    endif
    return
    end subroutine resprhon
    
    subroutine resprhonchk(k,il,ir,e,ew) ! check spherical symmetry of k.
    use brut
    use v6stepcalc
    integer(kind=i4) :: il,ir,i,j ! il<ir
    real(kind=r8) :: e,ew,k(3),xil(3,npart),xir(3,npart),m,tau,omega,kmod,dx(3),r,c1,rhonlr
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin) &
                       ,rhonl,rhonr,cpathcoreval
    if (il<=ir) then
        m=(mp+mn)/2
        rhonl=0
        rhonr=0
        rhonlr=0
        tau=(ir-il)*dt
        xil=ristra(il)%x
        xir=ristra(ir)%x
        call xconvert(xil)
        call xconvert(xir)
        kmod=sqrt(dot_product(k,k))
        do i=1,npart
            do j=1,npart
                dx(:)=xil(:,i)-xir(:,j)
		        r=sqrt(dot_product(dx,dx))                
                !rhonl=rhonl+exp(-ci*dot_product(k(:),xil(:,i))/hc)
                !rhonr=rhonr+exp(ci*dot_product(k(:),xir(:,i))/hc) 
                if ((kmod*r/hc)>tiny) then
                    c1=sin(kmod*r/hc)/(kmod*r/hc)
                else
                    c1=1
                endif                
                rhonlr=rhonlr+c1
            enddo
        enddo
        cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtr2l(:,:,0))
        e=real(rhonlr*cpathcoreval)
        omega=dot_product(k,k)/(2*m) ! omega_el. 
        ew=exp(tau*omega)*e  
        !write(6,*) 'k,il,ir,tau,omega,e,ew',k,il,ir,tau,omega,e,ew   
    else
        write(6,*) 'il > ir in resprhon, stop'
        stop             
    endif
    return
    end subroutine resprhonchk    
    

    subroutine resprhop(k,il,ir,e,ew) ! Carlson 1994 PRC paper.
    use brut
    use v6stepcalc
    integer(kind=i4) :: lspin,ipi,il,ir,i,l ! il<ir
    real(kind=r8) :: e,ew,k(3),x(3,npart),xil(3,npart),xir(3,npart),m,tau,omega
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin) &
                       ,rhopl(nisospin),rhopr(nisospin),cpathcoreval &
                       ,cwtr2lm(0:nspin-1,nisospin,0:nchorizo),cwtr2ltmp(0:nspin-1,nisospin,0:nchorizo) &
                       ,cwttmp(0:nspin-1,nisospin) 
    if (il>ir) stop 'il > ir in resprhop, stop'
    ! r2l
    m=(mp+mn)/2
    omega=dot_product(k,k)/(2*m) ! omega_el. 
    tau=(ir-il)*dt
    xil=ristra(il)%x
    xir=ristra(ir)%x
    call xconvert(xil)
    call xconvert(xir)    
    do l=1,nisospin
        rhopr(l)=0
        rhopl(l)=0
        lspin=liso(l)  
        do i=1,npart
		    ipi=and(1,shiftr(lspin,i-1))  !(2*(and(1,shiftr(lspin,i-1))))/2
		    if (ipi.eq.1) then
                rhopr(l)=rhopr(l)+exp(ci*dot_product(k(:),xir(:,i))/hc) 
                rhopl(l)=rhopl(l)+exp(-ci*dot_product(k(:),xil(:,i))/hc) 
            endif    
        enddo                    
    enddo    
    cwtr2ltmp=cwtr2l 
    if (il<ir) then     
        if (ir.eq.nchorizo) then
            do l=1,nisospin
                cwtnew(:,l)=cwtend(:,l)*rhopr(l)                              
            enddo 
            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
        else
            call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
            do l=1,nisospin
                cwtnew(:,l)=cwtr2lm(:,l,ir)*rhopr(l)
            enddo  
            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
        endif  
        if ((ir-il)>1) then
            do i=ir-1,il+1,-1        
               x=ristra(i)%x
               call vlsprop(spinorbit,1,x,ristra(i+1)%x,cwtr2ltmp(:,:,i+1),cwttmp(:,:))
               call v6proplr(x,x,-1,cwttmp(:,:),cwtr2ltmp(:,:,i))      
            enddo
        endif
        call vlsprop(spinorbit,1,xil,ristra(il+1)%x,cwtr2ltmp(:,:,il+1),cwttmp(:,:))
        call v6propl(xil,-1,cwttmp(:,:),cwtr2lm(:,:,il))
        do l=1,nisospin
            cwtnew(:,l)=cwtr2lm(:,l,il)*rhopl(l) 
        enddo 
        if (il/=0) then
            call v6propr(xil,-1,cwtnew(:,:),cwtr2ltmp(:,:,il)) 
            call vlsprop(spinorbit,1,ristra(il-1)%x,xil,cwtr2ltmp(:,:,il),cwttmp(:,:))
            cpathcoreval=sum(conjg(cwtl2r(:,:,il-1))*cwttmp(:,:))
        else ! il==0
            cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))    
        endif         
    else ! il==ir
        if (ir.eq.nchorizo) then
            do l=1,nisospin
                cwtnew(:,l)=rhopl(l)*rhopr(l)*cwtend(:,l)
            enddo 
            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))  
            call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
            cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))  
        else     
            if (ir.eq.0) then
                do l=1,nisospin
                    cwtnew(:,l)=rhopl(l)*rhopr(l)*cwtr2l(:,l,ir)             
                enddo
                cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))     
            else ! il=ir at somewhere middle
                call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                do l=1,nisospin     
                    cwtnew(:,l)=rhopl(l)*rhopr(l)*cwtr2lm(:,l,ir)
                enddo  
                call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))   
            endif   
        endif   
    endif
    e=real(cpathcoreval)
    ew=exp(tau*omega)*e    
    return
    end subroutine resprhop

    subroutine resprhopchk(k,il,ir,e,ew) ! Carlson 1994 PRC paper.
    use brut
    use v6stepcalc
    integer(kind=i4) :: lspin,ipi,jpi,il,ir,i,j,l,i1 ! il<ir
    real(kind=r8) :: e,ew,k(3),x(3,npart),xil(3,npart),xir(3,npart),m,tau,omega,kmod,dx(3),r,c1
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin) &
                       ,cpathcoreval &
                       ,cwtr2lm(0:nspin-1,nisospin,0:nchorizo),cwtr2ltmp(0:nspin-1,nisospin,0:nchorizo) &
                       ,cpathcorevaltot,cwttmp(0:nspin-1,nisospin) 
    if (il>ir) stop 'il > ir in resprhop, stop'
    ! r2l
    m=(mp+mn)/2
    omega=dot_product(k,k)/(2*m) ! omega_el. 
    tau=(ir-il)*dt
    xil=ristra(il)%x
    xir=ristra(ir)%x
    call xconvert(xil)
    call xconvert(xir)    
  
    cpathcorevaltot=0
    kmod=sqrt(dot_product(k,k))
    do i=1,npart
        do j=1,npart
            dx(:)=xil(:,i)-xir(:,j)
		    r=sqrt(dot_product(dx,dx))                
            !rhonl=rhonl+exp(-ci*dot_product(k(:),xil(:,i))/hc)
            !rhonr=rhonr+exp(ci*dot_product(k(:),xir(:,i))/hc) 
            if ((kmod*r/hc)>tiny) then
                c1=sin(kmod*r/hc)/(kmod*r/hc)
            else
                c1=1
            endif        
            
            cwtr2ltmp=cwtr2l 
            if (il<ir) then     
                if (ir.eq.nchorizo) then
                    do l=1,nisospin
                        lspin=liso(l)  
                        jpi=and(1,shiftr(lspin,j-1)) 
                        if (jpi.eq.1) then
                          cwtnew(:,l)=cwtend(:,l)
                        else
                          cwtnew(:,l)=0    
                        endif                                 
                    enddo 
                    call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                else
                    call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                    do l=1,nisospin
                        lspin=liso(l)  
                        jpi=and(1,shiftr(lspin,j-1)) 
                        if (jpi.eq.1) then
                          cwtnew(:,l)=cwtr2lm(:,l,ir)
                        else
                          cwtnew(:,l)=0    
                        endif   
                    enddo  
                    call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                endif  
                if ((ir-il)>1) then
                    do i1=ir-1,il+1,-1        
                       x=ristra(i1)%x
                       call vlsprop(spinorbit,1,x,ristra(i1+1)%x,cwtr2ltmp(:,:,i1+1),cwttmp(:,:))
                       call v6proplr(x,x,-1,cwttmp(:,:),cwtr2ltmp(:,:,i1))      
                    enddo
                endif
                call vlsprop(spinorbit,1,xil,ristra(il+1)%x,cwtr2ltmp(:,:,il+1),cwttmp(:,:))
                call v6propl(xil,-1,cwttmp(:,:),cwtr2lm(:,:,il))           
                do l=1,nisospin
                    lspin=liso(l)  
                    ipi=and(1,shiftr(lspin,i-1)) 
                    if (ipi.eq.1) then
                        cwtnew(:,l)=cwtr2lm(:,l,il)
                    else
                        cwtnew(:,l)=0    
                    endif   
                enddo                 
                if (il/=0) then
                    call v6propr(xil,-1,cwtnew(:,:),cwtr2ltmp(:,:,il)) 
                    call vlsprop(spinorbit,1,ristra(il-1)%x,xil,cwtr2ltmp(:,:,il),cwttmp(:,:))
                    cpathcoreval=sum(conjg(cwtl2r(:,:,il-1))*cwttmp(:,:))
                else ! il==0
                    cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))    
                endif         
            else ! il==ir
                if (ir.eq.nchorizo) then             
                    do l=1,nisospin
                        lspin=liso(l)  
                        ipi=and(1,shiftr(lspin,i-1)) 
                        jpi=and(1,shiftr(lspin,j-1)) 
                        if ((ipi.eq.1).and.(jpi.eq.1)) then
                          cwtnew(:,l)=cwtend(:,l)
                        else
                          cwtnew(:,l)=0    
                        endif   
                    enddo                         
                    call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir)) 
                    call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                    cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))  
                else     
                    if (ir.eq.0) then                       
                        do l=1,nisospin
                            lspin=liso(l)  
                            ipi=and(1,shiftr(lspin,i-1)) 
                            jpi=and(1,shiftr(lspin,j-1)) 
                            if ((ipi.eq.1).and.(jpi.eq.1)) then
                              cwtnew(:,l)=cwtr2l(:,l,ir)   
                            else
                              cwtnew(:,l)=0    
                            endif   
                        enddo                               
                        cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))     
                    else ! il=ir at somewhere middle
                        call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))                
                        do l=1,nisospin
                            lspin=liso(l)  
                            ipi=and(1,shiftr(lspin,i-1)) 
                            jpi=and(1,shiftr(lspin,j-1)) 
                            if ((ipi.eq.1).and.(jpi.eq.1)) then
                              cwtnew(:,l)=cwtr2lm(:,l,ir)
                            else
                              cwtnew(:,l)=0    
                            endif   
                        enddo                                 
                        call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                        call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                        cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))   
                    endif   
                endif   
            endif            
        cpathcorevaltot=cpathcorevaltot+cpathcoreval*c1
        enddo
    enddo
    e=real(cpathcorevaltot)
    ew=exp(tau*omega)*e    
    return
    end subroutine resprhopchk
      
    subroutine resprhot(k,il,ir,e,ew)
    use brut
    use v6stepcalc
    integer(kind=i4) :: lspin,ipi,il,ir,i,l ! il<ir
    real(kind=r8) :: e,ew,k(3),x(3,npart),xil(3,npart),xir(3,npart),m,tau,omega
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin) &
                       ,rhotl(nisospin),rhotr(nisospin),cpathcoreval &
                       ,cwtr2lm(0:nspin-1,nisospin,0:nchorizo),cwtr2ltmp(0:nspin-1,nisospin,0:nchorizo) &
                       ,cwttmp(0:nspin-1,nisospin) 
    logical :: isoscalar
    if (il>ir) stop 'il > ir in resprhop, stop'
    ! r2l
    isoscalar=.true. ! for v6'.
    if (isoscalar) then ! tau_dag <- tau_z
        m=(mp+mn)/2
        omega=dot_product(k,k)/(2*m) ! omega_el. 
        tau=(ir-il)*dt
        xil=ristra(il)%x
        xir=ristra(ir)%x
        call xconvert(xil)
        call xconvert(xir) 
        rhotr=0
        rhotl=0  
        do l=1,nisospin
            lspin=liso(l)  
            do i=1,npart
		        ipi=2*and(1,shiftr(lspin,i-1))-1  
                rhotr(l)=rhotr(l)+ipi*exp(ci*dot_product(k(:),xir(:,i))/hc) 
                rhotl(l)=rhotl(l)+ipi*exp(-ci*dot_product(k(:),xil(:,i))/hc)   
            enddo                    
        enddo    
        cwtr2ltmp=cwtr2l 
        if (il<ir) then     
            if (ir.eq.nchorizo) then
                do l=1,nisospin
                    cwtnew(:,l)=cwtend(:,l)*rhotr(l)                              
                enddo 
                call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
            else
                call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                do l=1,nisospin
                    cwtnew(:,l)=cwtr2lm(:,l,ir)*rhotr(l)
                enddo  
                call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
            endif  
            if ((ir-il)>1) then
                do i=ir-1,il+1,-1        
                   x=ristra(i)%x
                   call vlsprop(spinorbit,1,x,ristra(i+1)%x,cwtr2ltmp(:,:,i+1),cwttmp(:,:))
                   call v6proplr(x,x,-1,cwttmp(:,:),cwtr2ltmp(:,:,i))      
                enddo
            endif
            call vlsprop(spinorbit,1,xil,ristra(il+1)%x,cwtr2ltmp(:,:,il+1),cwttmp(:,:))
            call v6propl(xil,-1,cwttmp(:,:),cwtr2lm(:,:,il))
            do l=1,nisospin
                cwtnew(:,l)=cwtr2lm(:,l,il)*rhotl(l) 
            enddo 
            if (il/=0) then
                call v6propr(xil,-1,cwtnew(:,:),cwtr2ltmp(:,:,il)) 
                call vlsprop(spinorbit,1,ristra(il-1)%x,xil,cwtr2ltmp(:,:,il),cwttmp(:,:))
                cpathcoreval=sum(conjg(cwtl2r(:,:,il-1))*cwttmp(:,:))
            else ! il==0
                cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))    
            endif         
        else ! il==ir
            if (ir.eq.nchorizo) then
                do l=1,nisospin
                    cwtnew(:,l)=rhotl(l)*rhotr(l)*cwtend(:,l)
                enddo 
                call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))   
                call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))  
            else     
                if (ir.eq.0) then
                    do l=1,nisospin
                        cwtnew(:,l)=rhotl(l)*rhotr(l)*cwtr2l(:,l,ir)             
                    enddo
                    cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))     
                else ! il=ir at somewhere middle
                    call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                    do l=1,nisospin     
                        cwtnew(:,l)=rhotl(l)*rhotr(l)*cwtr2lm(:,l,ir)
                    enddo  
                    call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                    call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                    cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))   
                endif   
            endif   
        endif
        e=real(cpathcoreval)
        ew=exp(tau*omega)*e     
    else
        ! for v6 leave it empty now.  
    endif
    return
    end subroutine resprhot
    
    subroutine resprhotchk(k,il,ir,e,ew)
    use brut
    use v6stepcalc
    integer(kind=i4) :: lspin,ipi,il,ir,i,l,j,i1,jpi ! il<ir
    real(kind=r8) :: e,ew,k(3),x(3,npart),xil(3,npart),xir(3,npart),m,tau,omega,kmod,dx(3),r,c1
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin) &
                       ,rhotl(nisospin),rhotr(nisospin),cpathcoreval &
                       ,cwtr2lm(0:nspin-1,nisospin,0:nchorizo),cwtr2ltmp(0:nspin-1,nisospin,0:nchorizo) &
                       ,cpathcorevaltot,cwttmp(0:nspin-1,nisospin) 
    logical :: isoscalar
    if (il>ir) stop 'il > ir in resprhop, stop'
    ! r2l
    isoscalar=.true. ! for v6'.
    if (isoscalar) then ! tau_dag <- tau_z
        m=(mp+mn)/2
        omega=dot_product(k,k)/(2*m) ! omega_el. 
        tau=(ir-il)*dt
        xil=ristra(il)%x
        xir=ristra(ir)%x
        call xconvert(xil)
        call xconvert(xir) 
        rhotr=0
        rhotl=0  

        cpathcorevaltot=0
        kmod=sqrt(dot_product(k,k))        
        
        do i=1,npart
            do j=1,npart
                dx(:)=xil(:,i)-xir(:,j)
		        r=sqrt(dot_product(dx,dx))                
                !rhonl=rhonl+exp(-ci*dot_product(k(:),xil(:,i))/hc)
                !rhonr=rhonr+exp(ci*dot_product(k(:),xir(:,i))/hc) 
                if ((kmod*r/hc)>tiny) then
                    c1=sin(kmod*r/hc)/(kmod*r/hc)
                else
                    c1=1
                endif       
                cwtr2ltmp=cwtr2l 
                if (il<ir) then     
                    if (ir.eq.nchorizo) then
                        do l=1,nisospin
                            lspin=liso(l)  
                            jpi=2*and(1,shiftr(lspin,j-1))-1  
                            cwtnew(:,l)=cwtend(:,l)*jpi
                        enddo                 
                        call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                    else
                        call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                        do l=1,nisospin
                            lspin=liso(l)  
                            jpi=2*and(1,shiftr(lspin,j-1))-1  
                            cwtnew(:,l)=cwtr2lm(:,l,ir)*jpi
                        enddo                  
                        call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                    endif  
                    if ((ir-il)>1) then
                        do i1=ir-1,il+1,-1        
                           x=ristra(i1)%x
                           call vlsprop(spinorbit,1,x,ristra(i1+1)%x,cwtr2ltmp(:,:,i1+1),cwttmp(:,:))
                           call v6proplr(x,x,-1,cwttmp(:,:),cwtr2ltmp(:,:,i1))      
                        enddo
                    endif
                    call vlsprop(spinorbit,1,xil,ristra(il+1)%x,cwtr2ltmp(:,:,il+1),cwttmp(:,:))
                    call v6propl(xil,-1,cwttmp(:,:),cwtr2lm(:,:,il))          
                    do l=1,nisospin
                        lspin=liso(l)  
                        ipi=2*and(1,shiftr(lspin,i-1))-1  
                        cwtnew(:,l)=cwtr2lm(:,l,il)*ipi
                    enddo             
                    if (il/=0) then
                        call v6propr(xil,-1,cwtnew(:,:),cwtr2ltmp(:,:,il)) 
                        call vlsprop(spinorbit,1,ristra(il-1)%x,xil,cwtr2ltmp(:,:,il),cwttmp(:,:))
                        cpathcoreval=sum(conjg(cwtl2r(:,:,il-1))*cwttmp(:,:))
                    else ! il==0
                        cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))    
                    endif         
                else ! il==ir
                    if (ir.eq.nchorizo) then              
                        do l=1,nisospin
                            lspin=liso(l)  
                            jpi=2*and(1,shiftr(lspin,j-1))-1  
                            ipi=2*and(1,shiftr(lspin,i-1))-1
                            cwtnew(:,l)=cwtend(:,l)*ipi*jpi
                        enddo                 
                        call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))  
                        call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                        cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))  
                    else     
                        if (ir.eq.0) then   
                            do l=1,nisospin
                                lspin=liso(l)  
                                jpi=2*and(1,shiftr(lspin,j-1))-1  
                                ipi=2*and(1,shiftr(lspin,i-1))-1
                                cwtnew(:,l)=cwtr2l(:,l,ir)*ipi*jpi
                            enddo                     
                            cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))     
                        else ! il=ir at somewhere middle
                            call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                            do l=1,nisospin
                                lspin=liso(l)  
                                jpi=2*and(1,shiftr(lspin,j-1))-1  
                                ipi=2*and(1,shiftr(lspin,i-1))-1
                                cwtnew(:,l)=cwtr2lm(:,l,ir)*ipi*jpi
                            enddo                         
                            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                            call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                            cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))   
                        endif   
                    endif   
                endif            
                cpathcorevaltot=cpathcorevaltot+cpathcoreval*c1    
            enddo 
        enddo
        e=real(cpathcorevaltot)
        ew=exp(tau*omega)*e     
    else
        ! for v6 leave it empty now.  
    endif
    return
    end subroutine resprhotchk
     
    subroutine resprhostl(k,il,ir,e,ew)
    use brut
    use v6stepcalc
    integer(kind=i4) :: lspin,ipi,il,ir,i,j,l,spni ! il<ir
    real(kind=r8) :: e,ew,k(3),kz1,x(3,npart),xil(3,npart),xir(3,npart),m,tau,omega
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtnew1(0:nspin-1,nisospin) &
                       ,rhostll(0:nspin-1,nisospin),rhostlr(0:nspin-1,nisospin),cpathcoreval &
                       ,cwtr2lm(0:nspin-1,nisospin,0:nchorizo),cwtr2ltmp(0:nspin-1,nisospin,0:nchorizo) &
                       ,cwttmp(0:nspin-1,nisospin) 
    logical :: isoscalar
    if (il>ir) stop 'il > ir in resprhop, stop'
    ! r2l
    isoscalar=.true. ! for v6'.
    if (isoscalar) then ! tau_dag <- tau_z
        m=(mp+mn)/2
        omega=dot_product(k,k)/(2*m) ! omega_el. 
        tau=(ir-il)*dt
        xil=ristra(il)%x
        xir=ristra(ir)%x
        call xconvert(xil)
        call xconvert(xir)  
        if ((k(1).eq.0).and.(k(2).eq.0)) then ! kz case
            rhostlr=0
            rhostll=0  
            kz1=k(3)/abs(k(3))
            do l=1,nisospin
                lspin=liso(l)        
                do j=0,nspin-1 
                    do i=1,npart  
                        spni=2*and(1,shiftr(j,i-1))-1
		                ipi=2*and(1,shiftr(lspin,i-1))-1  
                        rhostlr(j,l)=rhostlr(j,l)+kz1*spni*ipi*exp(ci*dot_product(k(:),xir(:,i))/hc) 
                        rhostll(j,l)=rhostll(j,l)+kz1*spni*ipi*exp(-ci*dot_product(k(:),xil(:,i))/hc)       
                    enddo  
                enddo 
            enddo  
            cwtr2ltmp=cwtr2l 
            if (il<ir) then     
                if (ir.eq.nchorizo) then
                    do l=1,nisospin
                        do j=0,nspin-1
                            cwtnew(j,l)=cwtend(j,l)*rhostlr(j,l)   
                        enddo
                    enddo 
                    call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                else
                    call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                    do l=1,nisospin
                        do j=0,nspin-1
                            cwtnew(j,l)=cwtr2lm(j,l,ir)*rhostlr(j,l)
                        enddo
                    enddo  
                    call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                endif  
                if ((ir-il)>1) then
                    do i=ir-1,il+1,-1        
                       x=ristra(i)%x
                       call vlsprop(spinorbit,1,x,ristra(i+1)%x,cwtr2ltmp(:,:,i+1),cwttmp(:,:))
                       call v6proplr(x,x,-1,cwttmp(:,:),cwtr2ltmp(:,:,i))      
                    enddo
                endif
                call vlsprop(spinorbit,1,xil,ristra(il+1)%x,cwtr2ltmp(:,:,il+1),cwttmp(:,:))
                call v6propl(xil,-1,cwttmp(:,:),cwtr2lm(:,:,il))
                do l=1,nisospin
                    do j=0,nspin-1
                        cwtnew(j,l)=cwtr2lm(j,l,il)*rhostll(j,l) 
                    enddo
                enddo 
                if (il/=0) then
                    call v6propr(xil,-1,cwtnew(:,:),cwtr2ltmp(:,:,il)) 
                    call vlsprop(spinorbit,1,ristra(il-1)%x,xil,cwtr2ltmp(:,:,il),cwttmp(:,:))
                    cpathcoreval=sum(conjg(cwtl2r(:,:,il-1))*cwttmp(:,:))
                else ! il==0
                    cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))    
                endif         
            else ! il==ir
                if (ir.eq.nchorizo) then
                    do l=1,nisospin
                        do j=0,nspin-1
                            cwtnew(j,l)=rhostll(j,l)*rhostlr(j,l)*cwtend(j,l)
                        enddo
                    enddo 
                    call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))    
                    call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                    cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))  
                else     
                    if (ir.eq.0) then
                        do l=1,nisospin
                            do j=0,nspin-1
                                cwtnew(j,l)=rhostll(j,l)*rhostlr(j,l)*cwtr2l(j,l,ir) 
                            enddo
                        enddo
                        cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))     
                    else ! il=ir at somewhere middle
                        call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                        do l=1,nisospin  
                            do j=0,nspin-1
                                cwtnew(j,l)=rhostll(j,l)*rhostlr(j,l)*cwtr2lm(j,l,ir)
                            enddo
                        enddo  
                        call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                        call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                        cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))   
                    endif   
                endif   
            endif
        else ! general case for k 
            cwtr2ltmp=cwtr2l 
            if (il<ir) then     
                if (ir.eq.nchorizo) then
                    call resprhostlprop1(xir,k,cwtend,cwtnew) 
                    call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                else
                    call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))                 
                    call resprhostlprop1(xir,k,cwtr2lm(:,:,ir),cwtnew(:,:))
                    call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                endif  
                if ((ir-il)>1) then         
                    do i=ir-1,il+1,-1        
                       x=ristra(i)%x
                       call vlsprop(spinorbit,1,x,ristra(i+1)%x,cwtr2ltmp(:,:,i+1),cwttmp(:,:))
                       call v6proplr(x,x,-1,cwttmp(:,:),cwtr2ltmp(:,:,i))      
                    enddo      
                endif
                call vlsprop(spinorbit,1,xil,ristra(il+1)%x,cwtr2ltmp(:,:,il+1),cwttmp(:,:))
                call v6propl(xil,-1,cwttmp(:,:),cwtr2lm(:,:,il))                        
                call resprhostlprop1pr(xil,k,cwtr2lm(:,:,il),cwtnew(:,:)) ! -k for l.
                if (il/=0) then
                    call v6propr(xil,-1,cwtnew(:,:),cwtr2ltmp(:,:,il))
                    call vlsprop(spinorbit,1,ristra(il-1)%x,xil,cwtr2ltmp(:,:,il),cwttmp(:,:))
                    cpathcoreval=sum(conjg(cwtl2r(:,:,il-1))*cwttmp(:,:))
                else ! il==0
                    cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))    
                endif         
            else ! il==ir
                if (ir.eq.nchorizo) then                 
                    call resprhostlprop1(xir,k,cwtend(:,:),cwtnew1(:,:))
                    call resprhostlprop1pr(xil,k,cwtnew1(:,:),cwtnew(:,:))  
                    call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))   
                    call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                    cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))  
                else     
                    if (ir.eq.0) then
                        call resprhostlprop1(xir,k,cwtr2l(:,:,ir),cwtnew1(:,:))
                        call resprhostlprop1pr(xil,k,cwtnew1(:,:),cwtnew(:,:))                        
                        cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))     
                    else ! il=ir at somewhere middle
                        call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                        call resprhostlprop1(xir,k,cwtr2lm(:,:,ir),cwtnew1(:,:))
                        call resprhostlprop1pr(xil,k,cwtnew1(:,:),cwtnew(:,:))                         
                        call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                        call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                        cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))   
                    endif   
                endif   
            endif               
        endif
        e=real(cpathcoreval)
        ew=exp(tau*omega)*e                 
    else
        ! for v6 leave it empty now.  
    endif    
    return
    end subroutine resprhostl
    
    subroutine resprhostlprop1(x,k,cwtold,cwtnew) ! 
    integer(kind=i4) :: i,j,nstate
    real(kind=r8) :: x(3,npart),k(3),kmod,k1(3)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin),f(6),expf
    integer(kind=i4) :: spn,spni,ispn,ipi,newspn,newispn,newspni,newspnj,signi,signj,signij    
    cwtnew(:,:)=0 
    kmod=sqrt(dot_product(k,k))
    k1(:)=k(:)/kmod !unit vector
    do nstate=1,nbasis      
        spn=invspin(nstate)
        ispn=liso(invispin(nstate))  
        do i=1,npart     
             ipi=2*and(1,shiftr(ispn,i-1))-1  
             spni=2*and(1,shiftr(spn,i-1))-1
             expf=exp(ci*dot_product(k(:),x(:,i))/hc)
             f(1)=ipi*spni*k1(3)*expf*cwtold(invspin(nstate),invispin(nstate))   ! sigmaz
             f(2)=ipi*(k1(1)+spni*ci*k1(2))*expf*cwtold(invspin(nstate),invispin(nstate)) ! sigmax, sigmay           
             newspn=xor(spn,shiftl(1,i-1)) 
             cwtnew(spn,invispin(nstate))=cwtnew(spn,invispin(nstate))+f(1)
             cwtnew(newspn,invispin(nstate))=cwtnew(newspn,invispin(nstate))+f(2) 
        enddo
    enddo 
    return
    end subroutine resprhostlprop1   
    
    subroutine resprhostlprop1pr(x,k,cwtold,cwtnew) ! pr means the rho dagger.
    integer(kind=i4) :: i,j,nstate
    real(kind=r8) :: x(3,npart),k(3),kmod,k1(3)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin),f(6),expf 
    integer(kind=i4) :: spn,spni,ispn,ipi,newspn,newispn,newspni,newspnj,signi,signj,signij    
    cwtnew(:,:)=0 
    kmod=sqrt(dot_product(k,k))
    k1(:)=k(:)/kmod !unit vector
    do nstate=1,nbasis      
        spn=invspin(nstate)
        ispn=liso(invispin(nstate))  
        do i=1,npart     
             ipi=2*and(1,shiftr(ispn,i-1))-1  
             spni=2*and(1,shiftr(spn,i-1))-1
             expf=exp(-ci*dot_product(k(:),x(:,i))/hc)
             f(1)=ipi*spni*k1(3)*expf*cwtold(invspin(nstate),invispin(nstate))   ! sigmaz
             f(2)=ipi*(k1(1)+spni*ci*k1(2))*expf*cwtold(invspin(nstate),invispin(nstate)) ! sigmax, sigmay           
             newspn=xor(spn,shiftl(1,i-1)) 
             cwtnew(spn,invispin(nstate))=cwtnew(spn,invispin(nstate))+f(1)
             cwtnew(newspn,invispin(nstate))=cwtnew(newspn,invispin(nstate))+f(2) 
        enddo
    enddo 
    return
    end subroutine resprhostlprop1pr   
    
    subroutine resprhostlchk(k,il,ir,e,ew)
    use brut
    use v6stepcalc
    use math
    integer(kind=i4) :: lspin,ipi,il,ir,i,j,l,spni,jpi,i1,icase,icasestep ! il<ir
    real(kind=r8) :: e,ew,k(3),x(3,npart),xil(3,npart),xir(3,npart),m,tau,omega,kmod,dr(3),r,rx,ry,rz,c1(9),kr,r1(3),sj0,sj2
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtnew1(0:nspin-1,nisospin) &
                       ,cpathcoreval &
                       ,cwtr2lm(0:nspin-1,nisospin,0:nchorizo),cwtr2ltmp(0:nspin-1,nisospin,0:nchorizo) &
                       ,cpathcorevaltot,cwttmp(0:nspin-1,nisospin) 
    logical :: isoscalar
    if (il>ir) stop 'il > ir in resprhop, stop'
    ! r2l
    isoscalar=.true. ! for v6'.
    if (isoscalar) then ! tau_dag <- tau_z
        m=(mp+mn)/2
        omega=dot_product(k,k)/(2*m) ! omega_el. 
        tau=(ir-il)*dt
        xil=ristra(il)%x
        xir=ristra(ir)%x
        call xconvert(xil)
        call xconvert(xir)  
    ! general case for k      
        cpathcorevaltot=0
        kmod=sqrt(dot_product(k,k))     
        do i=1,npart
            do j=1,npart  
                dr(:)=-xil(:,i)+xir(:,j)   ! make rji.
		        r=sqrt(dot_product(dr,dr))  
                kr=kmod*r/hc    
                if (kr > tiny) then 
                    r1(:)=dr(:)/r
                    rx=r1(1)
                    ry=r1(2)
                    rz=r1(3)
                    call sphbesselj(0,kr,sj0) 
                    call sphbesselj(2,kr,sj2)  
                    c1(1)=(sj0+sj2*(rz**2-2*rx**2+ry**2))/3.0_r8
                    c1(2)=-sj2*rx*ry
                    c1(3)=-sj2*rx*rz
                    c1(4)=c1(2)
                    c1(5)=(sj0+sj2*(rz**2-2*ry**2+rx**2))/3.0_r8
                    c1(6)=-sj2*ry*rz
                    c1(7)=c1(3)
                    c1(8)=c1(6)
                    c1(9)=(sj0+sj2*(rx**2-2*rz**2+ry**2))/3.0_r8      
                    icasestep=1
                else
                    c1(:)=0
                    c1(1)=1.0_r8/3
                    c1(5)=c1(1)
                    c1(9)=c1(1)
                    icasestep=4
                endif
                do icase=1,9,icasestep
                    cwtr2ltmp=cwtr2l 
                    if (il<ir) then     
                        if (ir.eq.nchorizo) then                  
                            call resprhostlchkprop1(icase,j,cwtend,cwtnew) 
                            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                        else
                            call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))                 
                            call resprhostlchkprop1(icase,j,cwtr2lm(:,:,ir),cwtnew(:,:))
                            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                        endif  
                        if ((ir-il)>1) then         
                            do i1=ir-1,il+1,-1        
                               x=ristra(i1)%x
                               call vlsprop(spinorbit,1,x,ristra(i1+1)%x,cwtr2ltmp(:,:,i1+1),cwttmp(:,:))
                               call v6proplr(x,x,-1,cwttmp(:,:),cwtr2ltmp(:,:,i1))      
                            enddo      
                        endif
                        call vlsprop(spinorbit,1,xil,ristra(il+1)%x,cwtr2ltmp(:,:,il+1),cwttmp(:,:))
                        call v6propl(xil,-1,cwttmp(:,:),cwtr2lm(:,:,il))    
                        call resprhostlchkprop1pr(icase,i,cwtr2lm(:,:,il),cwtnew(:,:)) ! -k for l.
                        if (il/=0) then
                            call v6propr(xil,-1,cwtnew(:,:),cwtr2ltmp(:,:,il)) 
                            call vlsprop(spinorbit,1,ristra(il-1)%x,xil,cwtr2ltmp(:,:,il),cwttmp(:,:))
                            cpathcoreval=sum(conjg(cwtl2r(:,:,il-1))*cwttmp(:,:))
                        else ! il==0
                            cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))    
                        endif         
                    else ! il==ir
                        if (ir.eq.nchorizo) then   
                            call resprhostlchkprop1(icase,j,cwtend(:,:),cwtnew1(:,:))
                            call resprhostlchkprop1pr(icase,i,cwtnew1(:,:),cwtnew(:,:))  
                            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))   
                            call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                            cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))  
                        else     
                            if (ir.eq.0) then
                                call resprhostlchkprop1(icase,j,cwtr2l(:,:,ir),cwtnew1(:,:))
                                call resprhostlchkprop1pr(icase,i,cwtnew1(:,:),cwtnew(:,:))  
                                cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))     
                            else ! il=ir at somewhere middle
                                call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                                call resprhostlchkprop1(icase,j,cwtr2lm(:,:,ir),cwtnew1(:,:))
                                call resprhostlchkprop1pr(icase,i,cwtnew1(:,:),cwtnew(:,:))    
                                call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                                call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                                cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))   
                            endif   
                        endif   
                    endif 
                    cpathcorevaltot=cpathcorevaltot+cpathcoreval*c1(icase) 
                enddo
            enddo
        enddo
        e=real(cpathcorevaltot)
        ew=exp(tau*omega)*e                 
    else
        ! for v6 leave it empty now.  
    endif    
    return
    end subroutine resprhostlchk    
    
    subroutine resprhostlchkprop1(icase,i,cwtold,cwtnew) ! for xir
    integer(kind=i4) :: i,j,icase
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin)
    integer(kind=i4) :: spn,spni,ispn,ipi,newspn,newispn,newspni,newspnj,signi,signj,signij    
    select case (icase)
    case (1,4,7)   
        call resprhostlxchkprop1(i,cwtold,cwtnew)    
    case (2,5,8)
        call resprhostlychkprop1(i,cwtold,cwtnew)    
    case (3,6,9) 
        call resprhostlzchkprop1(i,cwtold,cwtnew)  
    case default 
        write (6,'(''icase value does not match resprhostlchkprop1'')') icase
        call abort         
    end select    
    return
    end subroutine resprhostlchkprop1 
    
    subroutine resprhostlchkprop1pr(icase,i,cwtold,cwtnew) ! for xil
    integer(kind=i4) :: i,j,icase
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin)
    integer(kind=i4) :: spn,spni,ispn,ipi,newspn,newispn,newspni,newspnj,signi,signj,signij    
    select case (icase)
    case (1,2,3)   
        call resprhostlxchkprop1(i,cwtold,cwtnew)    
    case (4,5,6)
        call resprhostlychkprop1(i,cwtold,cwtnew)    
    case (7,8,9) 
        call resprhostlzchkprop1(i,cwtold,cwtnew)  
    case default 
        write (6,'(''icase value does not match resprhostlchkprop1pr'')') icase
        call abort         
    end select    
    return
    end subroutine resprhostlchkprop1pr
    
    subroutine resprhostlxchkprop1(i,cwtold,cwtnew) ! 
    integer(kind=i4) :: i,j,nstate
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin),f(6)
    integer(kind=i4) :: spn,spni,ispn,ipi,newspn,newispn,newspni,newspnj,signi,signj,signij    
    cwtnew(:,:)=0 
    do nstate=1,nbasis      
        spn=invspin(nstate)
        ispn=liso(invispin(nstate))           
        ipi=2*and(1,shiftr(ispn,i-1))-1  
        spni=2*and(1,shiftr(spn,i-1))-1      
        f(2)=ipi*cwtold(invspin(nstate),invispin(nstate)) ! sigmax         
        newspn=xor(spn,shiftl(1,i-1)) 
        cwtnew(newspn,invispin(nstate))=cwtnew(newspn,invispin(nstate))+f(2)   
    enddo           
    return
    end subroutine resprhostlxchkprop1 
    
    subroutine resprhostlychkprop1(i,cwtold,cwtnew) ! 
    integer(kind=i4) :: i,j,nstate
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin),f(6)
    integer(kind=i4) :: spn,spni,ispn,ipi,newspn,newispn,newspni,newspnj,signi,signj,signij    
    cwtnew(:,:)=0 
    do nstate=1,nbasis      
        spn=invspin(nstate)
        ispn=liso(invispin(nstate))       
        ipi=2*and(1,shiftr(ispn,i-1))-1  
        spni=2*and(1,shiftr(spn,i-1))-1      
        f(2)=ipi*(spni*ci)*cwtold(invspin(nstate),invispin(nstate)) ! sigmay           
        newspn=xor(spn,shiftl(1,i-1)) 
        cwtnew(newspn,invispin(nstate))=cwtnew(newspn,invispin(nstate))+f(2)  
    enddo           
    return
    end subroutine resprhostlychkprop1     
    
    subroutine resprhostlzchkprop1(i,cwtold,cwtnew) ! 
    integer(kind=i4) :: i,j,nstate
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin),f(6)
    integer(kind=i4) :: spn,spni,ispn,ipi,newspn,newispn,newspni,newspnj,signi,signj,signij    
    cwtnew(:,:)=0 
    do nstate=1,nbasis      
        spn=invspin(nstate)
        ispn=liso(invispin(nstate))     
        ipi=2*and(1,shiftr(ispn,i-1))-1  
        spni=2*and(1,shiftr(spn,i-1))-1   
        f(1)=ipi*spni*cwtold(invspin(nstate),invispin(nstate))   ! sigmaz
        cwtnew(spn,invispin(nstate))=cwtnew(spn,invispin(nstate))+f(1)
    enddo          
    return
    end subroutine resprhostlzchkprop1    
    
    subroutine resprhostt(k,il,ir,eval,ewval,e,ew)
    use brut
    use v6stepcalc
    integer(kind=i4) :: icase,lspin,ipi,il,ir,i,j,l,spni ! il<ir
    real(kind=r8) :: eval,ewval,e(9),ew(9),k(3),x(3,npart),xil(3,npart),xir(3,npart),m,tau,omega
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtnew1(0:nspin-1,nisospin) &
                       ,cpathcoreval &
                       ,cwtr2lm(0:nspin-1,nisospin,0:nchorizo),cwtr2ltmp(0:nspin-1,nisospin,0:nchorizo) &
                       ,cwttmp(0:nspin-1,nisospin) 
    logical :: isoscalar
    ! icase 1: xx. 2: xy. 3: xz. 4: yx. 5: yy. 6: yz. 7: zx. 8: zy. 9: zz. 
    if (il>ir) stop 'il > ir in resprhop, stop'
    ! r2l
    isoscalar=.true. ! for v6'.
    if (isoscalar) then ! tau_dag <- tau_z
        m=(mp+mn)/2
        omega=dot_product(k,k)/(2*m) ! omega_el. 
        tau=(ir-il)*dt
        xil=ristra(il)%x
        xir=ristra(ir)%x
        call xconvert(xil)
        call xconvert(xir)  
        do icase=1,9,4  ! only add xx+yy+zz
            cwtr2ltmp=cwtr2l    
            if (il<ir) then     
                if (ir.eq.nchorizo) then
                    call resprhosttprop1(icase,xir,k,cwtend,cwtnew) 
                    call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                else
                    call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))                 
                    call resprhosttprop1(icase,xir,k,cwtr2lm(:,:,ir),cwtnew(:,:))
                    call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                endif  
                if ((ir-il)>1) then         
                    do i=ir-1,il+1,-1        
                        x=ristra(i)%x
                        call vlsprop(spinorbit,1,x,ristra(i+1)%x,cwtr2ltmp(:,:,i+1),cwttmp(:,:))
                        call v6proplr(x,x,-1,cwttmp(:,:),cwtr2ltmp(:,:,i))      
                    enddo      
                endif
                call vlsprop(spinorbit,1,xil,ristra(il+1)%x,cwtr2ltmp(:,:,il+1),cwttmp(:,:))
                call v6propl(xil,-1,cwttmp(:,:),cwtr2lm(:,:,il))                        
                call resprhosttprop1pr(icase,xil,k,cwtr2lm(:,:,il),cwtnew(:,:))
                if (il/=0) then
                    call v6propr(xil,-1,cwtnew(:,:),cwtr2ltmp(:,:,il)) 
                    call vlsprop(spinorbit,1,ristra(il-1)%x,xil,cwtr2ltmp(:,:,il),cwttmp(:,:))
                    cpathcoreval=sum(conjg(cwtl2r(:,:,il-1))*cwttmp(:,:))
                else ! il==0
                    cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))    
                endif         
            else ! il==ir
                if (ir.eq.nchorizo) then                 
                    call resprhosttprop1(icase,xir,k,cwtend(:,:),cwtnew1(:,:))
                    call resprhosttprop1pr(icase,xil,k,cwtnew1(:,:),cwtnew(:,:))  
                    call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))   
                    call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                    cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))  
                else     
                    if (ir.eq.0) then
                        call resprhosttprop1(icase,xir,k,cwtr2l(:,:,ir),cwtnew1(:,:))
                        call resprhosttprop1pr(icase,xil,k,cwtnew1(:,:),cwtnew(:,:))                        
                        cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))     
                    else ! il=ir at somewhere middle
                        call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                        call resprhosttprop1(icase,xir,k,cwtr2lm(:,:,ir),cwtnew1(:,:))
                        call resprhosttprop1pr(icase,xil,k,cwtnew1(:,:),cwtnew(:,:))                         
                        call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                        call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                        cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))   
                    endif   
                endif   
            endif  
            e(icase)=real(cpathcoreval)
            ew(icase)=exp(tau*omega)*e(icase)        
        enddo   
    else
        ! for other than v6 leave it empty now.  
    endif     
    eval=e(1)+e(5)+e(9)
    ewval=exp(tau*omega)*eval
    return
    end subroutine resprhostt


    subroutine resprhosttprop1(icase,x,k,cwtold,cwtnew) ! 
    integer(kind=i4) :: icase
    real(kind=r8) :: x(3,npart),k(3)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin)
    select case (icase)
    case (1,4,7)
        call resprhosttxprop1(x,k,cwtold,cwtnew) 
    case (2,5,8)
        call resprhosttyprop1(x,k,cwtold,cwtnew) 
    case (3,6,9)    
        call resprhosttzprop1(x,k,cwtold,cwtnew) 
    case default 
        write (6,'(''icase value does not match resprhosttprop1'')') icase
        call abort         
    end select
    return
    end subroutine resprhosttprop1  
    
    subroutine resprhosttprop1pr(icase,x,k,cwtold,cwtnew) ! 
    integer(kind=i4) :: icase
    real(kind=r8) :: x(3,npart),k(3)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin)
    select case (icase)
    case (1,2,3)
        call resprhosttxprop1pr(x,k,cwtold,cwtnew) 
    case (4,5,6)
        call resprhosttyprop1pr(x,k,cwtold,cwtnew) 
    case (7,8,9)    
        call resprhosttzprop1pr(x,k,cwtold,cwtnew) 
    case default 
        write (6,'(''icase value does not match resprhosttprop1'')') icase
        call abort         
    end select
    return
    end subroutine resprhosttprop1pr    
    
    subroutine resprhosttxprop1(x,k,cwtold,cwtnew) ! 
    integer(kind=i4) :: i,j,nstate
    real(kind=r8) :: x(3,npart),k(3),kmod,k1(3)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin),f(6),expf
    integer(kind=i4) :: spn,spni,ispn,ipi,newspn,newispn,newspni,newspnj,signi,signj,signij    
    cwtnew(:,:)=0 
    kmod=sqrt(dot_product(k,k))
    k1(:)=k(:)/kmod !unit vector
    do nstate=1,nbasis      
        spn=invspin(nstate)
        ispn=liso(invispin(nstate))  
        do i=1,npart     
             ipi=2*and(1,shiftr(ispn,i-1))-1  
             spni=2*and(1,shiftr(spn,i-1))-1
             expf=exp(ci*dot_product(k(:),x(:,i))/hc)
             f(1)=-ipi*spni*k1(2)*expf*cwtold(invspin(nstate),invispin(nstate))   ! sigmaz
             f(2)=ipi*(spni*ci*k1(3))*expf*cwtold(invspin(nstate),invispin(nstate)) ! sigmax, sigmay           
             newspn=xor(spn,shiftl(1,i-1)) 
             cwtnew(spn,invispin(nstate))=cwtnew(spn,invispin(nstate))+f(1)
             cwtnew(newspn,invispin(nstate))=cwtnew(newspn,invispin(nstate))+f(2) 
        enddo
    enddo 
    return
    end subroutine resprhosttxprop1     
    
    subroutine resprhosttxprop1pr(x,k,cwtold,cwtnew) ! 
    integer(kind=i4) :: i,j,nstate
    real(kind=r8) :: x(3,npart),k(3),kmod,k1(3)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin),f(6),expf
    integer(kind=i4) :: spn,spni,ispn,ipi,newspn,newispn,newspni,newspnj,signi,signj,signij    
    cwtnew(:,:)=0 
    kmod=sqrt(dot_product(k,k))
    k1(:)=k(:)/kmod !unit vector
    do nstate=1,nbasis      
        spn=invspin(nstate)
        ispn=liso(invispin(nstate))  
        do i=1,npart     
             ipi=2*and(1,shiftr(ispn,i-1))-1  
             spni=2*and(1,shiftr(spn,i-1))-1
             expf=exp(-ci*dot_product(k(:),x(:,i))/hc)
             f(1)=-ipi*spni*k1(2)*expf*cwtold(invspin(nstate),invispin(nstate))   ! sigmaz
             f(2)=ipi*(spni*ci*k1(3))*expf*cwtold(invspin(nstate),invispin(nstate)) ! sigmax, sigmay           
             newspn=xor(spn,shiftl(1,i-1)) 
             cwtnew(spn,invispin(nstate))=cwtnew(spn,invispin(nstate))+f(1)
             cwtnew(newspn,invispin(nstate))=cwtnew(newspn,invispin(nstate))+f(2) 
        enddo
    enddo 
    return
    end subroutine resprhosttxprop1pr  
    
    subroutine resprhosttyprop1(x,k,cwtold,cwtnew) ! 
    integer(kind=i4) :: i,j,nstate
    real(kind=r8) :: x(3,npart),k(3),kmod,k1(3)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin),f(6),expf
    integer(kind=i4) :: spn,spni,ispn,ipi,newspn,newispn,newspni,newspnj,signi,signj,signij    
    cwtnew(:,:)=0 
    kmod=sqrt(dot_product(k,k))
    k1(:)=k(:)/kmod !unit vector
    do nstate=1,nbasis      
        spn=invspin(nstate)
        ispn=liso(invispin(nstate))  
        do i=1,npart     
             ipi=2*and(1,shiftr(ispn,i-1))-1  
             spni=2*and(1,shiftr(spn,i-1))-1
             expf=exp(ci*dot_product(k(:),x(:,i))/hc)
             f(1)=ipi*spni*k1(1)*expf*cwtold(invspin(nstate),invispin(nstate))   ! sigmaz
             f(2)=ipi*(-k1(3))*expf*cwtold(invspin(nstate),invispin(nstate)) ! sigmax, sigmay           
             newspn=xor(spn,shiftl(1,i-1)) 
             cwtnew(spn,invispin(nstate))=cwtnew(spn,invispin(nstate))+f(1)
             cwtnew(newspn,invispin(nstate))=cwtnew(newspn,invispin(nstate))+f(2) 
        enddo
    enddo 
    return
    end subroutine resprhosttyprop1
    
    subroutine resprhosttyprop1pr(x,k,cwtold,cwtnew) ! 
    integer(kind=i4) :: i,j,nstate
    real(kind=r8) :: x(3,npart),k(3),kmod,k1(3)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin),f(6),expf
    integer(kind=i4) :: spn,spni,ispn,ipi,newspn,newispn,newspni,newspnj,signi,signj,signij    
    cwtnew(:,:)=0 
    kmod=sqrt(dot_product(k,k))
    k1(:)=k(:)/kmod !unit vector
    do nstate=1,nbasis      
        spn=invspin(nstate)
        ispn=liso(invispin(nstate))  
        do i=1,npart     
             ipi=2*and(1,shiftr(ispn,i-1))-1  
             spni=2*and(1,shiftr(spn,i-1))-1
             expf=exp(-ci*dot_product(k(:),x(:,i))/hc)
             f(1)=ipi*spni*k1(1)*expf*cwtold(invspin(nstate),invispin(nstate))   ! sigmaz
             f(2)=ipi*(-k1(3))*expf*cwtold(invspin(nstate),invispin(nstate)) ! sigmax, sigmay           
             newspn=xor(spn,shiftl(1,i-1)) 
             cwtnew(spn,invispin(nstate))=cwtnew(spn,invispin(nstate))+f(1)
             cwtnew(newspn,invispin(nstate))=cwtnew(newspn,invispin(nstate))+f(2) 
        enddo
    enddo 
    return
    end subroutine resprhosttyprop1pr 
    
    subroutine resprhosttzprop1(x,k,cwtold,cwtnew) ! 
    integer(kind=i4) :: i,j,nstate
    real(kind=r8) :: x(3,npart),k(3),kmod,k1(3)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin),f(6),expf
    integer(kind=i4) :: spn,spni,ispn,ipi,newspn,newispn,newspni,newspnj,signi,signj,signij    
    cwtnew(:,:)=0 
    kmod=sqrt(dot_product(k,k))
    k1(:)=k(:)/kmod !unit vector
    do nstate=1,nbasis      
        spn=invspin(nstate)
        ispn=liso(invispin(nstate))  
        do i=1,npart     
             ipi=2*and(1,shiftr(ispn,i-1))-1  
             spni=2*and(1,shiftr(spn,i-1))-1
             expf=exp(ci*dot_product(k(:),x(:,i))/hc)
             f(2)=ipi*(k1(2)-spni*ci*k1(1))*expf*cwtold(invspin(nstate),invispin(nstate)) ! sigmax, sigmay           
             newspn=xor(spn,shiftl(1,i-1)) 
             cwtnew(newspn,invispin(nstate))=cwtnew(newspn,invispin(nstate))+f(2) 
        enddo
    enddo 
    return
    end subroutine resprhosttzprop1    
    
    subroutine resprhosttzprop1pr(x,k,cwtold,cwtnew) ! 
    integer(kind=i4) :: i,j,nstate
    real(kind=r8) :: x(3,npart),k(3),kmod,k1(3)
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin),f(6),expf
    integer(kind=i4) :: spn,spni,ispn,ipi,newspn,newispn,newspni,newspnj,signi,signj,signij    
    cwtnew(:,:)=0 
    kmod=sqrt(dot_product(k,k))
    k1(:)=k(:)/kmod !unit vector
    do nstate=1,nbasis      
        spn=invspin(nstate)
        ispn=liso(invispin(nstate))  
        do i=1,npart     
             ipi=2*and(1,shiftr(ispn,i-1))-1  
             spni=2*and(1,shiftr(spn,i-1))-1
             expf=exp(-ci*dot_product(k(:),x(:,i))/hc)
             f(2)=ipi*(k1(2)-spni*ci*k1(1))*expf*cwtold(invspin(nstate),invispin(nstate)) ! sigmax, sigmay           
             newspn=xor(spn,shiftl(1,i-1)) 
             cwtnew(newspn,invispin(nstate))=cwtnew(newspn,invispin(nstate))+f(2) 
        enddo
    enddo 
    return
    end subroutine resprhosttzprop1pr    

    subroutine resprhosttchk(k,il,ir,e,ew)
    use brut
    use v6stepcalc
    use math
    integer(kind=i4) :: icase,icasestep,lspin,ipi,il,ir,i,j,l,spni,jpi,i1 ! il<ir
    real(kind=r8) :: e,ew,k(3),x(3,npart),xil(3,npart),xir(3,npart),m,tau,omega,kmod &
                    ,dr(3),r,rx,ry,rz,c1(9),c1tmp(9),kr,r1(3),sj0,sj2
    complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtnew1(0:nspin-1,nisospin) &
                       ,cpathcoreval &
                       ,cwtr2lm(0:nspin-1,nisospin,0:nchorizo),cwtr2ltmp(0:nspin-1,nisospin,0:nchorizo) &
                       ,cpathcorevaltot,cwttmp(0:nspin-1,nisospin) 
    logical :: isoscalar
    ! icase 1: xx. 2: xy. 3: xz. 4: yx. 5: yy. 6: yz. 7: zx. 8: zy. 9: zz. 
    if (il>ir) stop 'il > ir in resprhop, stop'
    ! r2l
    isoscalar=.true. ! for v6'.
    if (isoscalar) then ! tau_dag <- tau_z
        m=(mp+mn)/2
        omega=dot_product(k,k)/(2*m) ! omega_el. 
        tau=(ir-il)*dt
        xil=ristra(il)%x
        xir=ristra(ir)%x
        call xconvert(xil)
        call xconvert(xir)  
        cpathcorevaltot=0
        kmod=sqrt(dot_product(k,k))          
        do i=1,npart
            do j=1,npart
                dr(:)=-xil(:,i)+xir(:,j)   ! make rji.
		        r=sqrt(dot_product(dr,dr))  
                kr=kmod*r/hc
                if (kr > tiny) then
                    r1(:)=dr(:)/r
                    rx=r1(1)
                    ry=r1(2)
                    rz=r1(3)
                    call sphbesselj(0,kr,sj0) 
                    call sphbesselj(2,kr,sj2)  
                    c1tmp(1)=(sj0+sj2*(rz**2-2*rx**2+ry**2))/3.0_r8 ! xx
                    c1tmp(2)=-sj2*rx*ry ! xy
                    c1tmp(3)=-sj2*rx*rz ! xz
                    c1tmp(4)=c1tmp(2) ! yx
                    c1tmp(5)=(sj0+sj2*(rz**2-2*ry**2+rx**2))/3.0_r8 ! yy
                    c1tmp(6)=-sj2*ry*rz ! yz
                    c1tmp(7)=c1tmp(3) ! zx
                    c1tmp(8)=c1tmp(6) ! zy
                    c1tmp(9)=(sj0+sj2*(rx**2-2*rz**2+ry**2))/3.0_r8 ! zz      
                    icasestep=1
                else
                    c1tmp(:)=0
                    c1tmp(1)=1.0_r8/3
                    c1tmp(5)=c1tmp(1)
                    c1tmp(9)=c1tmp(1)   
                    icasestep=4
                endif
              ! rearrange
                    c1(1)=c1tmp(5)+c1tmp(9)
                    c1(2)=-c1tmp(2)                
                    c1(3)=-c1tmp(3)
                    c1(4)=-c1tmp(2)
                    c1(5)=c1tmp(1)+c1tmp(9)
                    c1(6)=-c1tmp(6)
                    c1(7)=-c1tmp(7)
                    c1(8)=-c1tmp(8)
                    c1(9)=c1tmp(1)+c1tmp(5)                             
                do icase=1,9,icasestep  
                      cwtr2ltmp=cwtr2l 
                    if (il<ir) then     
                        if (ir.eq.nchorizo) then         
                            call resprhostlchkprop1(icase,j,cwtend,cwtnew) 
                            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                        else
                            call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))                 
                            call resprhostlchkprop1(icase,j,cwtr2lm(:,:,ir),cwtnew(:,:))
                            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                        endif  
                        if ((ir-il)>1) then         
                            do i1=ir-1,il+1,-1        
                               x=ristra(i1)%x
                               call vlsprop(spinorbit,1,x,ristra(i1+1)%x,cwtr2ltmp(:,:,i1+1),cwttmp(:,:))
                               call v6proplr(x,x,-1,cwttmp(:,:),cwtr2ltmp(:,:,i1))      
                            enddo      
                        endif
                        call vlsprop(spinorbit,1,xil,ristra(il+1)%x,cwtr2ltmp(:,:,il+1),cwttmp(:,:))
                        call v6propl(xil,-1,cwttmp(:,:),cwtr2lm(:,:,il))    
                        call resprhostlchkprop1pr(icase,i,cwtr2lm(:,:,il),cwtnew(:,:)) ! -k for l.
                        if (il/=0) then
                            call v6propr(xil,-1,cwtnew(:,:),cwtr2ltmp(:,:,il)) 
                            call vlsprop(spinorbit,1,ristra(il-1)%x,xil,cwtr2ltmp(:,:,il),cwttmp(:,:))
                            cpathcoreval=sum(conjg(cwtl2r(:,:,il-1))*cwttmp(:,:))
                        else ! il==0
                            cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))    
                        endif         
                    else ! il==ir
                        if (ir.eq.nchorizo) then   
                            call resprhostlchkprop1(icase,j,cwtend(:,:),cwtnew1(:,:))
                            call resprhostlchkprop1pr(icase,i,cwtnew1(:,:),cwtnew(:,:))  
                            call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))   
                            call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                            cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))  
                        else     
                            if (ir.eq.0) then
                                call resprhostlchkprop1(icase,j,cwtr2l(:,:,ir),cwtnew1(:,:))
                                call resprhostlchkprop1pr(icase,i,cwtnew1(:,:),cwtnew(:,:))  
                                cpathcoreval=sum(conjg(cwtbegin(:,:))*cwtnew(:,:))     
                            else ! il=ir at somewhere middle
                                call v6proplpr(xir,-1,cwtr2l(:,:,ir),cwtr2lm(:,:,ir))
                                call resprhostlchkprop1(icase,j,cwtr2lm(:,:,ir),cwtnew1(:,:))
                                call resprhostlchkprop1pr(icase,i,cwtnew1(:,:),cwtnew(:,:))    
                                call v6propr(xir,-1,cwtnew(:,:),cwtr2ltmp(:,:,ir))
                                call vlsprop(spinorbit,1,ristra(ir-1)%x,xir,cwtr2ltmp(:,:,ir),cwttmp(:,:))
                                cpathcoreval=sum(conjg(cwtl2r(:,:,ir-1))*cwttmp(:,:))   
                            endif   
                        endif   
                    endif 
                    cpathcorevaltot=cpathcorevaltot+cpathcoreval*c1(icase) 
                enddo                 
            enddo 
        enddo
        e=real(cpathcorevaltot)
        ew=exp(tau*omega)*e         
    else
        ! for v6 leave it empty now.  
    endif    
    return
    end subroutine resprhosttchk

    subroutine readconfig(filename,iploin,jploin,iproin,jproin,x0tot1din)
    use mympi
    character(len=70) :: filename  
    integer(kind=i4) :: i,nptot
    integer(kind=i4) :: iploin(npair),jploin(npair),iproin(npair),jproin(npair)	
    integer(kind=i4) :: iplnow(npair),jplnow(npair),iprnow(npair),jprnow(npair)	 ! just for one time use.
    real(kind=r8) :: x0tot1din(3*npart*(nchorizo+1)),x0tot1dnow(3*npart*(nchorizo+1))  
    integer(kind=i4) :: myunit
    
    call barrier
    if (myrank() /= 0) then
        call recv(iploin,0,1)
        call recv(jploin,0,2)
        call recv(iproin,0,3)
        call recv(jproin,0,4)
        call recv(x0tot1din,0,5)                   
    else
        
        open(newunit=myunit,form='formatted',file=trim(filename),position='rewind')
            read(myunit,'(i10)') nptot        
         
            if (nproc() < nptot) then
                if (myrank().eq.0) then
                    write (6,'(''# of cores now < # of cores in infile.'',i10, '' vs '',i10  )'),nproc(),nptot   
                    write (12,'(''# of cores now < # of cores in infile.'',i10, '' vs '',i10  )'),nproc(),nptot  
                endif
            endif
            if (nproc() .eq. nptot) then
                if (myrank().eq.0) then
                    write (6,'(''# of cores now = # of cores in infile.'',i10, '' vs '',i10  )'),nproc(),nptot 
                    write (12,'(''# of cores now = # of cores in infile.'',i10, '' vs '',i10  )'),nproc(),nptot 
                endif                
            endif    
            if (nproc() > nptot) then
                if (myrank().eq.0) then
                    write (6,'(''# of cores now > # of cores in infile.'',i10, '' vs '',i10  )'),nproc(),nptot 
                    write (12,'(''# of cores now > # of cores in infile.'',i10, '' vs '',i10  )'),nproc(),nptot 
                endif               
            endif         
  
            read (myunit,'(<npair>i10)') iploin
            read (myunit,'(<npair>i10)') jploin
            read (myunit,'(<npair>i10)') iproin
            read (myunit,'(<npair>i10)') jproin
            read (myunit,'(3e15.7)') x0tot1din   
                       
            do i=1,nproc()-1
                if ((nproc() > nptot).and.(mod(i,nptot).eq.0)) then
                    rewind myunit
                    read (myunit,*) ! skip the nptot
                endif
                    
                read (myunit,'(<npair>i10)') iplnow(:)
                read (myunit,'(<npair>i10)') jplnow
                read (myunit,'(<npair>i10)') iprnow
                read (myunit,'(<npair>i10)') jprnow
                read (myunit,'(3e15.7)') x0tot1dnow                     
                call send(iplnow,i,1)
                call send(jplnow,i,2)
                call send(iprnow,i,3)
                call send(jprnow,i,4)
                call send(x0tot1dnow,i,5)                       
            enddo
            
        close(myunit)    
            
    endif    
     
    call barrier
    return
    end subroutine readconfig    
    
    subroutine writeoutx(irepin,it,iblock)
    use mympi
    integer(kind=i4) :: iplnow(npair),jplnow(npair),iprnow(npair),jprnow(npair)	 ! just for one time use.
    real(kind=r8) :: xtot(3,npart,0:nchorizo),xtot1d(3*npart*(nchorizo+1))
    integer(kind=i4) :: i
    integer(kind=i4) :: irepin,it,iblock
 
    call chorizoallout(xtot,iplnow,jplnow,iprnow,jprnow)  
    !if ((i.eq.1).and.(myrank().eq.0)) time8=mpi_wtime()  
    xtot1d=reshape(xtot,(/3*npart*(nchorizo+1)/))

    call barrier
    if (myrank() /= 0) then
        call send(iplnow,0,1)
        call send(jplnow,0,2)
        call send(iprnow,0,3)
        call send(jprnow,0,4)
        call send(xtot1d,0,5)          
    else
        call rename(trim(outfile),'he4ristra.out.old')
        if (iblock>neq) call writecheckpoint(irepin,it) ! write out checkpoint.
        open(unit=9,form='formatted',file='xouttmp.dat',position='rewind')
            write(9,'(i10)') nproc()  
            
            write (9,'(<npair>i10)') iplnow
            write (9,'(<npair>i10)') jplnow
            write (9,'(<npair>i10)') iprnow
            write (9,'(<npair>i10)') jprnow
            write (9,'(3e15.7)') xtot1d   
                    
            do i=1,nproc()-1
                call recv(iplnow,i,1)
                call recv(jplnow,i,2)
                call recv(iprnow,i,3)
                call recv(jprnow,i,4)
                call recv(xtot1d,i,5)   
                ! rank i write out.        
                write (9,'(<npair>i10)') iplnow
                write (9,'(<npair>i10)') jplnow
                write (9,'(<npair>i10)') iprnow
                write (9,'(<npair>i10)') jprnow
                write (9,'(3e15.7)') xtot1d             
            enddo
            
        close(9) 
        call rename('xouttmp.dat',trim(outfile)) ! if the name is outfile then means writeoutx is successful.       
    endif
    
    return
    end subroutine writeoutx     
    
    subroutine writeoutxgs(irepin,it,iblock) ! gs means gather scatter.
    use mympi
    integer(kind=i4) :: iplnow(npair),jplnow(npair),iprnow(npair),jprnow(npair)	   
    integer(kind=i4) :: ixtemp,l,myunit
    integer(kind=i4) :: irepin,it,iblock
    real(kind=r8) :: xtot(3,npart,0:nchorizo)
    
    if (enablewritepathgs) then
        ixtemp=size(xtot)
        call chorizoallout(xtot,iplnow,jplnow,iprnow,jprnow)  
        !if ((i.eq.1).and.(myrank().eq.0)) time8=mpi_wtime()
        call barrier
        call gather(reshape(xtot,(/ixtemp/)),xgatherall)
        call gather(iplnow,iplgatherall)
        call gather(jplnow,jplgatherall)
        call gather(iprnow,iprgatherall)
        call gather(jprnow,jprgatherall)   
        !if ((i.eq.1).and.(myrank().eq.0)) time9=mpi_wtime()
        if (myrank().eq.0) then
            call rename(trim(outfile),'he4ristra.out.old')
            if (iblock>neq) call writecheckpoint(irepin,it) ! write out checkpoint.
            open(newunit=myunit,form='formatted',file='xouttmp.dat',position='rewind')
                write(myunit,'(i10)') nproc()
                do l=0,nproc()-1
                    write (myunit,'(<npair>i10)') iplgatherall((npair*l+1):(npair*(l+1))) 
                    write (myunit,'(<npair>i10)') jplgatherall((npair*l+1):(npair*(l+1))) 
                    write (myunit,'(<npair>i10)') iprgatherall((npair*l+1):(npair*(l+1))) 
                    write (myunit,'(<npair>i10)') jprgatherall((npair*l+1):(npair*(l+1))) 
                    write (myunit,'(3e15.7)') xgatherall((ixtemp*l+1):(ixtemp*(l+1)))    
                enddo 
            close(myunit) 
            call rename('xouttmp.dat',trim(outfile)) ! if the name is outfile then means writeoutx is successful.
        endif            
    else
        call writeoutx(irepin,it,iblock)  
    endif
    
    return
    end subroutine writeoutxgs
    
    subroutine checkhpsi(vn,vd,vnke,vnpev6,vnpeem,vnpels,psi20,psi2n,rf,icheck &
                        ,iplnow,jplnow,iprnow,jprnow,x0tot1d,nblocknow,ipathnow &
                        ,hpsipass)
    use mympi
    integer(kind=i4) :: icheck,i,nblocknow,ipathnow
    character(len=120) :: filename
    integer(kind=i4) :: myunit
    real(kind=r8) :: x0tot1d(3*npart*(nchorizo+1)),vn,vd,vnke,vnpev6,vnpeem,vnpels,psi20,psi2n,rf
    integer(kind=i4) :: iplnow(npair),jplnow(npair),iprnow(npair),jprnow(npair)
    logical :: hpsipass
    real(kind=r8), parameter :: valstrange=1000_r8,kestrange=2000_r8,pestrange=2000_r8
    real(kind=r8), parameter :: limcvalstrange=10000
    !call gather(icheck,icheckhpsi1d) ! gather to rank 0.  ! this gather bcast may cause bug, check why.
    !call bcast(icheckhpsi1d)     
    if (icheck /= 0) then  
        hpsipass=.false.    
        write(6,'( '' check hpsi fail! bugbeads recorded, abort! '' )') 
        write(filename,'("rank",i10,".bugbeads")') myrank()
        open(newunit=myunit,form='formatted',file=trim(filename),position='rewind')  
            write (myunit,'(i10)') myrank()
            write (myunit,'(2e15.7,1x, t40, "Percentage diff is: ", g15.7, " %")') psi20,psi2n,abs((psi20-psi2n)/((psi20+psi2n)/2))*100 ! % difference
            write (myunit,'(<npair>i10)') iplnow
            write (myunit,'(<npair>i10)') jplnow
            write (myunit,'(<npair>i10)') iprnow
            write (myunit,'(<npair>i10)') jprnow
            write (myunit,'(3e15.7)') x0tot1d
        close(myunit)             
    else   
        hpsipass=.true.
        if ((abs(vn/vd)>valstrange).or.(abs(vnke)>kestrange).or.(abs(vnpev6)>pestrange).and.(nblocknow>neq)) then
            call writemypath(vn,vd,vnke,vnpev6,vnpeem,vnpels,psi20,psi2n,rf,nblocknow,ipathnow)
            cvalstrange=cvalstrange+1
        endif      
    endif 
    call barrier
    if (icheck/=0) call abort  ! instead of just kill the process, can make it smarter, like copy a config from rank0.
    if (cvalstrange>limcvalstrange) then
        write (6,*) 'Abort! Too many abnormal values even after equilibrium, check bug, rank= ', myrank()
        call abort
    endif
    return
    end subroutine checkhpsi      
    
    subroutine writemypath(vn,vd,vnke,vnpev6,vnpeem,vnpels,psi20,psi2n,rf,nblocknow,ipathnow) ! write out the path for just one core.
    use mympi
    integer(kind=i4) :: i,ipathin,myunit1,myunit2,nblocknow,ipathnow
    real(kind=r8) :: x0tot(3,npart,0:nchorizo),x0tot1d(3*npart*(nchorizo+1)) &
                     ,vn,vd,vnke,vnpev6,vnpeem,vnpels,psi20,psi2n,rf 
    character(len=120) :: filename
    call chorizoallout(x0tot,iplo,jplo,ipro,jpro)   
    x0tot1d=reshape(x0tot,(/3*npart*(nchorizo+1)/))    
    write(filename,'("rank",i10,".path")') myrank()
    open(newunit=myunit1,form='formatted',file=trim(filename),position='append')
        write (myunit1,*) '---------'
        write (myunit1,'(2i10)') nblocknow,ipathnow
        write (myunit1,'(9e15.7,1x)') vn,vd,vnke,vnpev6,vnpeem,vnpels,psi20,psi2n,rf
        write (myunit1,*)
        write (myunit1,'(<npair>i10)') iplo
        write (myunit1,'(<npair>i10)') jplo
        write (myunit1,'(<npair>i10)') ipro
        write (myunit1,'(<npair>i10)') jpro
        write (myunit1,'(3e15.7)') x0tot1d
        write (myunit1,*) '---------'
    close(myunit1)  
    return
    end subroutine writemypath    

    subroutine writepath(ipathin,iploin,jploin,iproin,jproin,x0tot1din,neqfin)
    use mympi
    integer(kind=i4) :: i,ipathin,myunit1,myunit2
    integer(kind=i4) :: iploin(npair),jploin(npair),iproin(npair),jproin(npair)	
    integer(kind=i4) :: iplnow(npair),jplnow(npair),iprnow(npair),jprnow(npair)	 ! just for one time use.
    real(kind=r8) :: x0tot1din(3*npart*(nchorizo+1)),x0tot1dnow(3*npart*(nchorizo+1))  
    character(len=120) :: filename1,filename2
    logical :: neqfin
    real(kind=r8) time0,time1,time2,time3
    
    
    !if (myrank().eq.0) time0=mpi_wtime() 
    
    call barrier
    if (myrank() /= 0) then
        call send(iploin,0,1)
        call send(jploin,0,2)
        call send(iproin,0,3)
        call send(jproin,0,4)
        call send(x0tot1din,0,5)           
    else  
      
        if (neqfin) then
            filename1=file_path
            filename2=file_pathnumbercount
        else
            filename1=file_path0
            filename2=file_path0numbercount          
            
        endif

        open(newunit=myunit1,form='unformatted',file=trim(filename1),position='append')
            if (ipathin >= 2) backspace myunit1 ! to overwrite ipathmark  
            
            write(myunit1) ipathin 
            write (myunit1) myrank()
            
            write (myunit1) iploin
            write (myunit1) jploin
            write (myunit1) iproin
            write (myunit1) jproin
            write (myunit1) x0tot1din
            
            do i=1,nproc()-1
                call recv(iplnow,i,1)
                call recv(jplnow,i,2)
                call recv(iprnow,i,3)
                call recv(jprnow,i,4)
                call recv(x0tot1dnow,i,5)    
                
                write (myunit1) i
                write (myunit1) iplnow
                write (myunit1) jplnow
                write (myunit1) iprnow
                write (myunit1) jprnow
                write (myunit1) x0tot1dnow  
                
            enddo      
            write(myunit1) ipathmark   
        close(myunit1) 
        open(newunit=myunit2,form='formatted',file=trim(filename2),position='rewind')
            write (myunit2,'(3i10)') nproc(),ipathin,size(x0tot1din)
        close(myunit2)     
        
        !time3=mpi_wtime() 
        !write (6,*) 'write send recv takes: ', time3-time0
            
    endif   
    
    
    return
    end subroutine writepath
    
    
    subroutine writepathgs(ipathin,iploin,jploin,iproin,jproin,x0tot1din,neqfin)
    use mympi
    integer(kind=i4) :: ixtemp,l
    integer(kind=i4), allocatable :: iplall(:),jplall(:),iprall(:),jprall(:)	! for mpi_gather
    real(kind=r8), allocatable :: xall(:) ! for mpi_gather
    integer(kind=i4) :: ipathin,myunit1,myunit2
    integer(kind=i4) :: iploin(npair),jploin(npair),iproin(npair),jproin(npair)	
    real(kind=r8) :: x0tot1din(3*npart*(nchorizo+1))
    character(len=120) :: filename1,filename2
    logical :: neqfin   
    real(kind=r8) time0,time1,time2,time3

    if (enablewritepathgs) then
        
        !if (myrank().eq.0) time0=mpi_wtime() 
        
        ixtemp=size(x0tot1din)     
        call barrier
        call gather(x0tot1din,xgatherall)
        call gather(iploin,iplgatherall)
        call gather(jploin,jplgatherall)
        call gather(iproin,iprgatherall)
        call gather(jproin,jprgatherall)   
        
        
        !if (myrank().eq.0) time1=mpi_wtime() 
        
        
        if (myrank().eq.0) then
        
            if (neqfin) then
                filename1=file_path
                filename2=file_pathnumbercount
            else
                filename1=file_path0
                filename2=file_path0numbercount           
            
            endif
            
            
            !time2=mpi_wtime() 
            
            open(newunit=myunit1,form='unformatted',file=trim(filename1),position='append')
                if (ipathin >= 2) backspace myunit1 ! to overwrite ipathmark          
                write(myunit1) ipathin 
                do l=0,nproc()-1     
                    write (myunit1) l
                    write (myunit1) iplgatherall((npair*l+1):(npair*(l+1))) 
                    write (myunit1) jplgatherall((npair*l+1):(npair*(l+1))) 
                    write (myunit1) iprgatherall((npair*l+1):(npair*(l+1))) 
                    write (myunit1) jprgatherall((npair*l+1):(npair*(l+1))) 
                    write (myunit1) xgatherall((ixtemp*l+1):(ixtemp*(l+1)))        
                enddo   
                write(myunit1) ipathmark   
            close(myunit1) 
            open(newunit=myunit2,form='formatted',file=trim(filename2),position='rewind')
                write (myunit2,'(3i10)') nproc(),ipathin,size(x0tot1din)
            close(myunit2) 
        
   
            !time3=mpi_wtime() 
            !write (6,*) 'gather takes: ', time1-time0
            !write (6,*) 'write takes: ', time3-time2
            
        endif
 
    else
        call writepath(ipathin,iploin,jploin,iproin,jproin,x0tot1din,neqfin)  
    endif
    
    
    return
    end subroutine writepathgs
    
        
    subroutine writecheckpoint(irepin,it)
    use estimator
    use estimatorristra
    use wavefunction
    use mympi
    integer(kind=i4) :: it
    real(kind=r8) :: wttotio(nest),valtotio(nest),val2totio(nest)
    real(kind=r8) :: wttotristraio(nestristra),valtotristraio(0:nchorizo,nestristra),val2totristraio(0:nchorizo,nestristra)
    real(kind=r8) :: wtrhodisttotio(0:nchorizo),rhodisttotio(0:nrhobin,0:nchorizo),rhodist2totio(0:nrhobin,0:nchorizo) &
                    ,rhondisttotio(0:nrhobin,0:nchorizo),rhondist2totio(0:nrhobin,0:nchorizo) &
                    ,rhopdisttotio(0:nrhobin,0:nchorizo),rhopdist2totio(0:nrhobin,0:nchorizo)
    real(kind=r8) :: wte0resptotio(nrespest),vale0resptotio(0:ne0respbeads,nrespest),vale0resp2totio(0:ne0respbeads,nrespest)
    integer(kind=i4) :: irepin
    logical :: input
    integer(kind=i4) :: myunit
    character(len=120) :: file_checkpointold

    if (myrank().eq.0) then  
        input=.false. ! false means output mode.
        file_checkpointold=trim(file_checkpoint) // '.old'   
        call rename(file_checkpoint,trim(file_checkpointold))
        open(newunit=myunit,form='unformatted',file='checkpointtmp.unf',position='rewind') 
            write (myunit) nproc()
            if ((irepin.eq.2).or.(irepin.eq.4)) then
                write (myunit) it
                write (myunit) ipathcountold 
            endif
            if ((irepin.eq.5).or.(irepin.eq.6)) then
                write (myunit) ipathcountold    
            endif
  
    ! estimators    
            call checkpointest(input,wttotio,valtotio,val2totio)  
            write (myunit) wttotio
            write (myunit) valtotio
            write (myunit) val2totio
    ! estimatorsristra  
            call checkpointestristra(input,wttotristraio,valtotristraio,val2totristraio)
            write (myunit) wttotristraio
            write (myunit) valtotristraio
            write (myunit) val2totristraio          
    ! density distribution
            call checkpointrhodist(input,wtrhodisttotio,rhodisttotio,rhodist2totio &
                                                       ,rhondisttotio,rhondist2totio,rhopdisttotio,rhopdist2totio)
            write (myunit) wtrhodisttotio
            write (myunit) rhodisttotio
            write (myunit) rhodist2totio 
            write (myunit) rhondisttotio
            write (myunit) rhondist2totio
            write (myunit) rhopdisttotio
            write (myunit) rhopdist2totio    
    ! response function
            if ((irepin.eq.5).or.(irepin.eq.6)) then  ! need to set nrespest, ne0respbeads first.
                call checkpointresponse(input,wte0resptotio,vale0resptotio,vale0resp2totio)
                write (myunit) wte0resptotio
                write (myunit) vale0resptotio
                write (myunit) vale0resp2totio         
            endif
        close(myunit)   
        call rename('checkpointtmp.unf',file_checkpoint) ! if it is checkpoint.unf then it means writecheckpoint is successful.
    endif
    return 
    end subroutine writecheckpoint
    
    subroutine readcheckpoint(irepin,it) ! the initialization for resume mode.
    use estimator
    use estimatorristra
    use wavefunction
    use mympi
    integer(kind=i4) :: it,i,k,nprocprev
    real(kind=r8) :: wttotio(nest),valtotio(nest),val2totio(nest)
    real(kind=r8) :: wttotristraio(nestristra),valtotristraio(0:nchorizo,nestristra),val2totristraio(0:nchorizo,nestristra)
    real(kind=r8) :: wtrhodisttotio(0:nchorizo),rhodisttotio(0:nrhobin,0:nchorizo),rhodist2totio(0:nrhobin,0:nchorizo) &
                    ,rhondisttotio(0:nrhobin,0:nchorizo),rhondist2totio(0:nrhobin,0:nchorizo) &
                    ,rhopdisttotio(0:nrhobin,0:nchorizo),rhopdist2totio(0:nrhobin,0:nchorizo)
    real(kind=r8) :: wte0resptotio(nrespest),vale0resptotio(0:ne0respbeads,nrespest),vale0resp2totio(0:ne0respbeads,nrespest)
    character(len=120), dimension(:), allocatable :: answer 
    integer(kind=i4) :: irepin,myunit
    logical :: input
    
    if (myrank().eq.0) then  
        input=.true. ! false means output mode.
        open(newunit=myunit,form='unformatted',file=file_checkpoint,position='rewind')    
            read (myunit) nprocprev    
            if (nprocprev /= nproc()) then ! protection of resume mode
                write (6,*) ' # of processors does not match checkpoint.unf record, Resume Mode Stop! ',nproc(),nprocprev     
                call abort
            endif
            if ((irepin.eq.2).or.(irepin.eq.4)) then
                read (myunit) it
                write (6,*) 'it= ',it
                read (myunit) ipathcountold   
                write (6,*) 'ipathcountold = ',ipathcountold 
                ipath=ipathcountold+1
                ipathtmp=1 
            endif
            if ((irepin.eq.5).or.(irepin.eq.6)) then
                read (myunit) ipathcountold  
                write (6,*) 'ipathcountold = ',ipathcountold  
            endif  
    ! estimators    
            read (myunit) wttotio
            read (myunit) valtotio
            read (myunit) val2totio
            call checkpointest(input,wttotio,valtotio,val2totio)
    ! estimatorsristra  
            read (myunit) wttotristraio
            read (myunit) valtotristraio
            read (myunit) val2totristraio 
            call checkpointestristra(input,wttotristraio,valtotristraio,val2totristraio)
    ! density distribution
            read (myunit) wtrhodisttotio
            read (myunit) rhodisttotio
            read (myunit) rhodist2totio         
            read (myunit) rhondisttotio
            read (myunit) rhondist2totio
            read (myunit) rhopdisttotio
            read (myunit) rhopdist2totio              
            call checkpointrhodist(input,wtrhodisttotio,rhodisttotio,rhodist2totio &
                                                       ,rhondisttotio,rhondist2totio,rhopdisttotio,rhopdist2totio)
    ! response function
            if ((irepin.eq.5).or.(irepin.eq.6)) then  ! need to set nrespest, ne0respbeads first.    
                read (myunit) wte0resptotio
                read (myunit) vale0resptotio
                read (myunit) vale0resp2totio    
                call checkpointresponse(input,wte0resptotio,vale0resptotio,vale0resp2totio)
            endif
        close(myunit)    
              
! show last time results
        allocate(answer(nest))
        answer=resstring0()
        call writeanswer(answer)
        deallocate(answer)
        write (6,'(/,''last config recalculated. ----------------------------------------'')') 
	    write (12,'(/,''last config recalculated. ----------------------------------------'')')    

! need to check path.unf and locate it at correct position.   
        call checkpathunf(irepin,ipathcountold)     
        write (6,*) ' readcheckpoint finished ---------------------------------------- '    
    endif    
              
    return 
    end subroutine readcheckpoint     
  
    

    subroutine checkpathunf(irepin,ipathcountoldin) ! call readcheckpoint first. check and finalize path.unf
    use mympi
    integer(kind=i4) :: i,myunit,ixtemp
    integer(kind=i4) :: irepin,ipathcountoldin,ranknow,ipathnow,ipathnow2,ipathtarget,ipathtarget1,ipathtarget2
    logical :: ipathlocated    
    
    if (myrank().eq.0) then
        
       select case (irepin)  
       case (2,4)    
            ipathtarget=ipathcountoldin+1
            ipathlocated=.false. 
            open(unit=19,form='unformatted',file=file_path,position='append')         
                backspace 19
                read (19) ipathmarkchk
                backspace 19
                if (ipathmarkchk.eq.ipathmark) then ! means previous path.unf writing is finished.
                    write (6,*) 'Previous path.unf writing is finished, ipathcountoldin = ',ipathcountoldin  
                    do while (.not.ipathlocated)   
                        do i=1,pathunfstep1*nproc()  ! just use pathunfstep1*nproc()  instead of pathunfstep1(it is 6).
                            backspace 19  
                        enddo
                        read (19) ranknow
                        !write (6,*) 'ranknow=',ranknow
                        backspace 19
                        if (ranknow.eq.0) then
                           backspace 19
                           read (19) ipathnow
                           !write (6,*) 'ranknow, ipathnow = ',ranknow,ipathnow
                           if (ipathnow.eq.ipathtarget) then ! ipath is correctly located.
                               write (6,*) ' path.unf writing is located(ipathtarget) ',ipathtarget  ! may need to add target2.
                               ipathlocated=.true.
                               backspace 19                         
                               write(19) ipathmark
                               exit
                           else if (ipathnow.eq.ipathcountoldin) then
                               write (6,*) ' path.unf writing is located(ipathcountoldin) ',ipathcountoldin  ! may need to add target2.
                               ipathlocated=.true.
                               exit ! exit and do nothing.
                           else    
                               backspace 19
                           endif     
                        else
                            cycle                
                        endif              
                    enddo                          
                else ! means previously path.unf writing is not finished.
                    write (6,*) 'Previous path.unf writing is NOT finished',ipathcountoldin                
                    read (19) ipathnow
                    backspace 19
                    if (ipathnow.eq.ipathtarget) then
                        do i=1,nproc()*pathunfstep1
                            backspace 19
                        enddo
                        backspace 19
                        read (19) ipathnow2
                        backspace 19
                        if (ipathnow2.eq.ipathcountoldin) then
                            write (6,*) ' path.unf writing is now located(ipathtarget) ',ipathtarget
                            ipathlocated=.true.
                            do i=1,nproc()*pathunfstep1
                                read (19)
                            enddo
                            write(19) ipathmark                   
                        endif
                    endif
                    do while (.not.ipathlocated)                   
                       read (19) ranknow
                       backspace 19        
                       if (ranknow.eq.0) then   
                           backspace 19
                           read (19) ipathnow
                           if (ipathnow.eq.ipathtarget) then ! ipath is correctly located.
                               write (6,*) ' path.unf writing is now located(ipathtarget) ',ipathtarget
                               ipathlocated=.true.
                               backspace 19
                               write(19) ipathmark        
                               exit
                           else
                               backspace 19
                           endif     
                       else
                           backspace 19
                           cycle 
                       endif
                    enddo             
                endif  
            close (19) 
            
            ! update the content in 'pathnumbercounts.txt'
            ixtemp=3*npart*(nchorizo+1) 
            open(newunit=myunit,form='formatted',file=file_pathnumbercount,position='rewind')
            write (myunit,'(3i10)') nproc(),ipathcountoldin,ixtemp
            close(myunit)         
 
       case (5,6) 
           
           ! call computeinit first to get np and npo
            ipathtarget=ipath  
            ipathlocated=.false. 
            open(unit=19,form='unformatted',file=file_path,position='append')         
                backspace 19
                read (19) ipathmarkchk
                backspace 19
                if (ipathmarkchk.eq.ipathmark) then ! means previous path.unf writing is finished.
                    write (6,*) 'Previous path.unf writing is finished, ipath = ', ipath                         
                else ! means previously path.unf writing is not finished.
                    write (6,*) 'Previous path.unf writing is NOT finished, ipath = ',ipath               
                    do while (.not.ipathlocated)                   
                       read (19) ranknow
                       backspace 19        
                       if (ranknow.eq.0) then   
                           backspace 19
                           read (19) ipathnow
                           if (ipathnow.eq.ipathtarget) then ! ipath is correctly located.
                               write (6,*) ' path.unf writing is now located(ipathtarget) ',ipathtarget
                               ipathlocated=.true.  
                               do i=1,pathunfstep1*npo 
                                   read (19)
                               enddo
                               write(19) ipathmark ! finalize path.unf       
                               exit
                           else
                               backspace 19
                           endif     
                       else
                           backspace 19
                           cycle 
                       endif
                    enddo             
                endif  
            close (19)       

            !ipathtarget1=ipathcountold+1
            !ipathtarget2=ipathcountold+2
            !ipathlocated=.false. 
            !        
            !open(unit=19,form='unformatted',file='path.unf',position='append')           
            !read (19) ipathnow
            !if (ipathnow.eq.ipathtarget1) then
            !    read (19) ranknow
            !    if (ranknow.eq.0) then
            !        write (6,*) ' path.unf is located(special target1) '
            !        ipathlocated=.true.
            !        backspace 19
            !        backspace 19
            !    endif  
            !else
            !    backspace 19
            !    backspace 19   
            !endif
            !
            !do while (.not.ipathlocated)  
            !    read (19) ipathnow
            !    if ((ipathnow.eq.ipathtarget1).or.(ipathnow.eq.ipathtarget2)) then
            !        read (19) ranknow
            !        if (ranknow.eq.0) then
            !            if (ipathnow.eq.ipathtarget1) write (6,*) 'path.unf located(target1)',ipathcountold,ipathtarget1
            !            if (ipathnow.eq.ipathtarget2) write (6,*) 'path.unf located(target2)',ipathcountold,ipathtarget2
            !            ipathlocated=.true.
            !            backspace 19
            !            backspace 19
            !            exit 
            !        endif    
            !    endif   
            !enddo              
            
       case default
            write (6,'(''irep value has not been defined in checkpathunf '')') irepin
            call abort     
       end select
       
    write (6,*) ' check path.unf finished ---------------------------------------- '    
    endif
    
    return
    end subroutine checkpathunf

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
    integer(kind=i4) :: ipl(npair),jpl(npair),ipr(npair),jpr(npair)	
    do i=0,nchorizo
        xout(:,:,i)=ristra(i)%x
    enddo
    call stepoutputordlr(ipl,jpl,ipr,jpr)
    return
    end subroutine chorizoallout     
    
    subroutine getstepstdcalctime(tout)
    real(kind=r8) :: tout
    tout=time00
    return
    end subroutine getstepstdcalctime
    
    subroutine getstepstdMPIpathtime(tout)
    real(kind=r8) :: tout
    tout=time5tot
    return
    end subroutine getstepstdMPIpathtime 

    
    
! ---------------- the below propagator did not add LS propagator
    
    !subroutine bisectv6bextraflex(ileftbisect,mmaxnow) ! does not seem to work well.
    !! for ileftbisect > (nchorizo-nbisect) case.
    !! based on v6a version. move beads first, then do an overall judge. add an extra 0-nbisect beads on right. 
    !use wavefunction
    !use estimator
    !use random
    !use v6stepcalc
    !integer(kind=i4) :: ileftbisect,il,ir,i,imid,k,bilvl,m,mmaxnow,nbisectnow
    !real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    !real(kind=r8) :: xmid(3,npart),xd2old(3,npart),xd2new(3,npart)
    !real(kind=r8) :: tmp,tmp1test1,rt2
    !real(kind=r8) :: oldlv(0:mmaxnow),newlv(0:mmaxnow)
    !logical :: reject,fake
    !if (ileftbisect <= (nchorizo-nbisect)) then
    !    call bisectv6bflex(ileftbisect,mmaxnow)
    !    return
    !endif
    !if (mmaxnow.eq.mmax) ibisectextratot=ibisectextratot+1  
    !nbisectnow=2**mmaxnow
    !m=nbisectnow-(nchorizo-ileftbisect+1)
    !
    !!write(6,'(''m= '',t40,i10)') m
    !
    !fake=.true.
    !if (fake) then ! use wise fake gaussian move instead of bisection move. to be continued...
    !
    !else ! this one is dead now. linking 0th and nchorizo th beads does not work, and does not really make sense.
    !    newristra(0)=ristra(ileftbisect)
    !    newristra(nbisectnow)=ristra(m)
    !    tau=dt*nbisectnow ! total time for the bisection slice.
    !    ! judge if we accept the newristra or not.
    !    ! in this v6 case, v is not justpotential v, it should be the value of < ... >. lv.   
    !    do bilvl=1,mmaxnow  ! besection.  level 1 to N. mmaxnow = the total level N.	 
    !        ! move the current level   
    !        sigmamid=sqrt(hbar*tau/2**bilvl)  
	   !     !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	   !     do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		  !      imid=2**(mmaxnow-bilvl)+i*2**(mmaxnow-bilvl+1)
		  !      il=imid-2**(mmaxnow-bilvl)
		  !      ir=imid+2**(mmaxnow-bilvl)
		  !      gauss(:,:)=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		  !      xmid=(newristra(il)%x+newristra(ir)%x)/2
		  !      newristra(imid)%x(:,:)=xmid(:,:)+gauss(:,:)  
    !            !write(12,*) 'bilvl il ir imid sigmamid =',bilvl, il,ir,imid,sigmamid
    !            !write(12,*) 'gauss =', gauss
    !        enddo         
    !    enddo   
    !    ! then judge the last level 
    !    bisectextratot(mmaxnow)=bisectextratot(mmaxnow)+1
    !    ristranew=ristra 
    !    ristranew(ileftbisect:nchorizo)=newristra(0:nchorizo-ileftbisect) 
    !    ristranew(0:m)=newristra(nchorizo-ileftbisect+1:nbisectnow)     
    !    call calcpathcoreval(ristranew,cwtbeginnew,cwtendnew,cwtl2rnew,cwtr2lnew,pathcorevalnew)    
	   ! newlv(mmaxnow)=pathcorevalnew
	   ! oldlv(mmaxnow)=pathcoreval  
    !    xd2old=(ristra(nchorizo)%x-ristra(0)%x)**2  
    !    xd2new=(ristranew(nchorizo)%x-ristranew(0)%x)**2  
    !    rt2=-0.5_r8/(dt*(2*hbar))*(sum(xd2old)-sum(xd2new))	 ! this rt2 is different from normal 1 by 1 bc it is from T not PI, so old-new.         
    !    tmp=log(newlv(mmaxnow))-log(oldlv(mmaxnow))+rt2
    !    !write(6,'(''tmp,newlv,oldlv,rt2: '',t40,4(g15.7,1x))') tmp,log(newlv(mmaxnow)),log(oldlv(mmaxnow)),rt2  
    !    rn=randn(1)	   
	   ! if (log(rn(1)).lt.tmp) then
	   !     reject=.false. 
	   !     bisectextracount(mmaxnow)=bisectextracount(mmaxnow)+1	 
    !        if (mmaxnow.eq.mmax) then
    !            ibisectextra=ibisectextra+1 
	   !         call addval(12,1.0_r8,1.0_r8)  
    !        endif        
    !        ristra=ristranew
    !        cwtbegin(:,:)=cwtbeginnew(:,:)
    !        cwtend(:,:)=cwtendnew(:,:) 
    !        cwtl2r=cwtl2rnew
    !        cwtr2l=cwtr2lnew
    !        pathcoreval=pathcorevalnew            
	   ! else
	   !     reject=.true. 
    !        if (mmaxnow.eq.mmax) call addval(12,0.0_r8,1.0_r8)    
    !    endif   
    !endif
    !return  
    !end subroutine bisectv6bextraflex              
    !   
    
    
    !subroutine stdmove1fastv6flexa(ileft,iright,arrow,imp) 
    !! add importance sampling. need to Check. after having more table points for trial wave function it should work.
    !! move 1 by 1 with bisection 1 bead move.
    !use wavefunction
    !use estimator
    !use v6stepcalc
    !use random
    !use brut
    !use mympi
    !integer(kind=i4) :: i,k,ileft,iright,icnow
    !real(kind=r8) :: rn(1),rn2(2),gauss(3,npart),gc
    !real(kind=r8) :: x(3,npart),xnew(3,npart),xnew1(3,npart),xold(3,npart),xold1(3,npart) &
    !                ,xoldnxt(3,npart),xd2old(3,npart),xd2new(3,npart),xd2old1(3,npart),xd2new1(3,npart) &
	   !             ,xr(3,npart),xl(3,npart)
    !real(kind=r8) :: rt2,tmp1,tmp2,tmp3,tmp4,q1,q2,gnew,gold,gnew1,gold1,logratio,tmp1test1,tmp1test2
    !complex(kind=r8) :: cwtold(0:nspin-1,nisospin),cwtnew(0:nspin-1,nisospin) &
    !                   ,cwtold1(0:nspin-1,nisospin),cwtnew1(0:nspin-1,nisospin)
    !complex(kind=r8) :: cwtl2rtmp(0:nspin-1,nisospin),cwtr2ltmp(0:nspin-1,nisospin) &
    !                   ,cwtl2rtmpnew(0:nspin-1,nisospin),cwtl2rtmpnew1(0:nspin-1,nisospin),cwtl2rtmpold1(0:nspin-1,nisospin) &
    !                   ,cwtr2ltmpnew(0:nspin-1,nisospin),cwtr2ltmpnew1(0:nspin-1,nisospin),cwtr2ltmpold1(0:nspin-1,nisospin)
    !integer(kind=i4) :: arrow
    !logical :: imp
    !
    !if ((ileft.lt.0).or.(iright.gt.nchorizo).or.(iright.lt.0).or.(ileft.gt.nchorizo).or.(ileft.gt.iright)) then
	   ! write (6,*) 'stdmove1fastv6flex range error', ileft,iright
	   ! stop
    !endif
    !if (ileft.eq.iright) then
    !    call stdmove1fastv6(ileft,iright) 
    !    return
    !endif       
    !gc=-0.5_r8/(dt*2*hbar)   
    !select case (arrow)    
    !case (-1) ! inverse order from right to left  
    !    icnow=0
    !    icstdmove1bisect1rvstot=icstdmove1bisect1rvstot+1
    !    !icmov1bisect1rvstot=abs(iright-ileft+1)*icstdmove1bisect1rvstot 
    !    do i=iright,ileft,-1  
    !        if (i.eq.0) then      
    !            if (.not.imp) then
    !                icmov1rvstot=icmov1rvstot+1
    !                xold=ristra(i)%x 
    !                gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                xnew(:,:)=xold(:,:)+gauss(:,:)  
    !                tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i)))))       
    !                call corop(xnew,iplo,jplo,cwtnew)
    !                call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmp(:,:))
    !                tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
    !                xd2new(:,:)=(ristra(i+1)%x(:,:)-xnew(:,:))**2
    !                xd2old(:,:)=(ristra(i+1)%x(:,:)-xold(:,:))**2
    !                rt2=gc*(sum(xd2new)-sum(xd2old))	  
    !                logratio=tmp2-tmp1+rt2         
	   !             rn=randn(1)
	   !             if (log(rn(1)).lt.logratio) then ! dt ne 0
    !                    icnow=icnow+1               
    !                    icmov1rvs=icmov1rvs+1                             
	   !                 ristra(i)%x(:,:)=xnew(:,:)
	   !                 cwtr2l(:,:,i)=cwtr2ltmp(:,:)   
    !                ! update cwtbegin
    !                    cwtbegin(:,:)=cwtnew(:,:)  
    !                    pathcoreval=exp(tmp2)
    !                endif              
    !            else ! imp        
    !                icmov1rvstot=icmov1rvstot+1
    !                gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                xold=ristra(i)%x 
    !                xnew(:,:)=xold(:,:)+gauss(:,:)  
    !                xold1=2*xnew-xold
    !                xnew1=xold-gauss
    !                xd2new(:,:)=(ristra(i+1)%x(:,:)-xnew(:,:))**2
    !                xd2old(:,:)=(ristra(i+1)%x(:,:)-xold(:,:))**2
    !                xd2new1=(ristra(i+1)%x-xnew1)**2
    !                xd2old1=(ristra(i+1)%x-xold1)**2    
    !                gnew=exp(gc*sum(xd2new))
    !                gold=exp(gc*sum(xd2old))
    !                gnew1=exp(gc*sum(xd2new1))
    !                gold1=exp(gc*sum(xd2old1))
    !                call corop(xnew,iplo,jplo,cwtnew)
    !                call corop(xold1,iplo,jplo,cwtold1)
    !                call corop(xnew1,iplo,jplo,cwtnew1)
    !                call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew(:,:)) 
    !                call v6propl(xold1,-1,cwtr2l(:,:,i+1),cwtr2ltmpold1(:,:))
	   !             call v6propl(xnew1,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew1(:,:)) 
    !                tmp1=pathcoreval
    !                tmp2=abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmpnew(:,:))))       
    !                tmp3=abs(real(sum(conjg(cwtold1(:,:))*cwtr2ltmpold1(:,:))))                         
	   !             tmp4=abs(real(sum(conjg(cwtnew1(:,:))*cwtr2ltmpnew1(:,:))))                         
    !                q1=(gnew*tmp2+gnew1*tmp4)/(gold*tmp1+gold1*tmp3)
    !                q2=gnew*tmp2/(gnew*tmp2+gnew1*tmp4)                            
	   !             rn2=randn(2)
	   !             if (rn2(1) < q1) then ! dt ne 0
    !                    icnow=icnow+1 
    !                    icmov1rvs=icmov1rvs+1                    
    !                    if (rn2(2) < q2) then
	   !                     ristra(i)%x(:,:)=xnew(:,:)
	   !                     cwtr2l(:,:,i)=cwtr2ltmpnew(:,:)   
    !                        cwtbegin(:,:)=cwtnew(:,:)  
    !                        pathcoreval=tmp2                                    
    !                    else
	   !                     ristra(i)%x=xnew1
	   !                     cwtr2l(:,:,i)=cwtr2ltmpnew1(:,:)   
    !                        cwtbegin=cwtnew1	   
    !                        pathcoreval=tmp4                                   
    !                    endif
    !                endif                                         
    !            endif
    !        else  
	   !         if (i.eq.nchorizo) then      
    !                if (.not.imp) then           
    !                    icmov1rvstot=icmov1rvstot+1
    !                    xold=ristra(i)%x 
    !                    gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                    xnew(:,:)=xold(:,:)+gauss(:,:)                          
	   !                 tmp1=log(pathcoreval)  !log(abs(real(sum(conjg(cwtl2r(:,:,i))*cwtold(:,:))))) ! can put those tmp1 thing in the beginning, no need to repreatedly calculate                
	   !                 call corop(xnew,ipro,jpro,cwtnew)
    !                    call v6propr(xnew,-1,cwtnew(:,:),cwtr2ltmp(:,:))          
    !                    tmp2=log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmp(:,:)))))
	   !                 xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2
	   !                 xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2
	   !                 rt2=gc*(sum(xd2new)-sum(xd2old))	  
	   !                 logratio=tmp2-tmp1+rt2                        
	   !                 rn=randn(1)
	   !                 if (log(rn(1)).lt.logratio) then
    !                        icnow=icnow+1
    !                        icmov1rvs=icmov1rvs+1   
	   !                     ristra(i)%x(:,:)=xnew(:,:)
	   !                     cwtr2l(:,:,i)=cwtr2ltmp(:,:)
    !                    ! update cwtend and prepare for the left bead.
    !                        cwtend(:,:)=cwtnew(:,:)
    !                        x=ristra(i-1)%x
		  !                  call v6proplr(x,x,-1,cwtr2l(:,:,i),cwtr2l(:,:,i-1))
    !                        pathcoreval=exp(tmp2)
    !                    endif 	               
    !                else ! imp  i=nchorizo                
    !                    icmov1rvstot=icmov1rvstot+1
    !                    gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                    xold=ristra(i)%x 
    !                    xnew(:,:)=xold(:,:)+gauss(:,:)  
    !                    xold1=2*xnew-xold
    !                    xnew1=xold-gauss
    !                    xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2
    !                    xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2
    !                    xd2new1=(ristra(i-1)%x-xnew1)**2
    !                    xd2old1=(ristra(i-1)%x-xold1)**2    
    !                    gnew=exp(gc*sum(xd2new))
    !                    gold=exp(gc*sum(xd2old))
    !                    gnew1=exp(gc*sum(xd2new1))
    !                    gold1=exp(gc*sum(xd2old1)) 
    !                    call corop(xnew,ipro,jpro,cwtnew)
    !                    call corop(xold1,ipro,jpro,cwtold1)
    !                    call corop(xnew1,ipro,jpro,cwtnew1)                           
    !                    call v6propr(xnew,-1,cwtnew(:,:),cwtr2ltmpnew(:,:))                         
    !                    call v6propr(xold1,-1,cwtold1(:,:),cwtr2ltmpold1(:,:))          
	   !                 call v6propr(xnew1,-1,cwtnew1(:,:),cwtr2ltmpnew1(:,:))
    !                    tmp1=pathcoreval
    !                    tmp2=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpnew(:,:))))              
    !                    tmp3=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpold1(:,:))))                                 
	   !                 tmp4=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpnew1(:,:))))                         
    !                    q1=(gnew*tmp2+gnew1*tmp4)/(gold*tmp1+gold1*tmp3)
    !                    q2=gnew*tmp2/(gnew*tmp2+gnew1*tmp4)                               
	   !                 rn2=randn(2)
	   !                 if (rn2(1) < q1) then ! dt ne 0
    !                        icnow=icnow+1 
    !                        icmov1rvs=icmov1rvs+1                        
    !                        if (rn2(2) < q2) then
	   !                         ristra(i)%x(:,:)=xnew(:,:)
	   !                         cwtr2l(:,:,i)=cwtr2ltmpnew(:,:)   
    !                            cwtend(:,:)=cwtnew(:,:)	   
    !                            pathcoreval=tmp2                             
    !                        else
	   !                         ristra(i)%x=xnew1
	   !                         cwtr2l(:,:,i)=cwtr2ltmpnew1(:,:)   
    !                            cwtend=cwtnew1	   
    !                            pathcoreval=tmp4                                   
    !                        endif
    !                        x=ristra(i-1)%x
    !                        call v6proplr(x,x,-1,cwtr2l(:,:,i),cwtr2l(:,:,i-1)) 
    !                    endif                             
    !                endif       
    !            else                 
    !                if (.not.imp) then 
    !                    icmov1rvstot=icmov1rvstot+1
    !                    xold=ristra(i)%x 
    !                    gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                    xnew(:,:)=xold(:,:)+gauss(:,:)                           
	   !                 call v6proplr(xnew,xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmp(:,:))
	   !                 tmp2=log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmp(:,:)))))
    !                    tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2l(:,:,i)))))  
	   !                 xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2+(ristra(i+1)%x(:,:)-xnew(:,:))**2
	   !                 xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2+(ristra(i+1)%x(:,:)-xold(:,:))**2   
	   !                 rt2=gc*(sum(xd2new)-sum(xd2old))		   
	   !                 logratio=tmp2-tmp1+rt2                                                             
	   !                 rn=randn(1)
	   !                 if (log(rn(1)).lt.logratio) then
    !                        icnow=icnow+1           
    !                        icmov1rvs=icmov1rvs+1   
	   !                     ristra(i)%x(:,:)=xnew(:,:)
	   !                     cwtr2l(:,:,i)=cwtr2ltmp(:,:)	
    !                        pathcoreval=exp(tmp2)
    !                    endif                 
    !                else  ! imp                    
    !                    icmov1rvstot=icmov1rvstot+1
    !                    gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                    xold=ristra(i)%x
    !                    xnew(:,:)=xold(:,:)+gauss(:,:)  
    !                    xold1=2*xnew-xold
    !                    xnew1=xold-gauss
    !                    xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2+(ristra(i+1)%x(:,:)-xnew(:,:))**2
    !                    xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2+(ristra(i+1)%x(:,:)-xold(:,:))**2
    !                    xd2new1=(ristra(i-1)%x-xnew1)**2+(ristra(i+1)%x-xnew1)**2
    !                    xd2old1=(ristra(i-1)%x-xold1)**2+(ristra(i+1)%x-xold1)**2    
    !                    gnew=exp(gc*sum(xd2new))
    !                    gold=exp(gc*sum(xd2old))
    !                    gnew1=exp(gc*sum(xd2new1))
    !                    gold1=exp(gc*sum(xd2old1)) 
    !                    call v6proplr(xnew,xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew(:,:))
    !                    call v6proplr(xold1,xold1,-1,cwtr2l(:,:,i+1),cwtr2ltmpold1(:,:))
    !                    call v6proplr(xnew1,xnew1,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew1(:,:))         
    !                    tmp1=pathcoreval                   
    !                    tmp2=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpnew(:,:))))                          
    !                    tmp3=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpold1(:,:))))                                             
	   !                 tmp4=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpnew1(:,:))))                                     
    !                    q1=(gnew*tmp2+gnew1*tmp4)/(gold*tmp1+gold1*tmp3)
    !                    q2=gnew*tmp2/(gnew*tmp2+gnew1*tmp4)                                                          
	   !                 rn2=randn(2)
	   !                 if (rn2(1) < q1) then ! dt ne 0
    !                        icnow=icnow+1 
    !                        icmov1rvs=icmov1rvs+1                      
    !                        if (rn2(2) < q2) then
	   !                         ristra(i)%x(:,:)=xnew(:,:)
	   !                         cwtr2l(:,:,i)=cwtr2ltmpnew(:,:)   	   
    !                            pathcoreval=tmp2                             
    !                        else
	   !                         ristra(i)%x=xnew1
	   !                         cwtr2l(:,:,i)=cwtr2ltmpnew1(:,:)   	   
    !                            pathcoreval=tmp4                                   
    !                        endif
    !                    endif                                
    !                endif
    !                if ( icnow /= 0 ) then
    !            ! update r2l for the next r2l bead, the left bead
    !                    if (i >= 2) then
		  !                  x=ristra(i-1)%x
		  !                  call v6proplr(x,x,-1,cwtr2l(:,:,i),cwtr2l(:,:,i-1))
    !                    else ! i=1
    !                        xl=ristra(i-1)%x
    !                        call v6propl(xl,-1,cwtr2l(:,:,i),cwtr2l(:,:,i-1))	    
    !                    endif               
    !                endif 
	   !         endif 
	   !     endif                   
    !    enddo
    !    if (icnow /= 0) then 
    !        xl=ristra(0)%x
    !        xr=ristra(nchorizo)%x
    !        ! update r2l
    !        if ( ileft >= 3 ) then
    !            do k=ileft,3,-1 
    !               x=ristra(k-2)%x
	   !            call v6proplr(x,x,-1,cwtr2l(:,:,k-1),cwtr2l(:,:,k-2))  
    !            enddo
    !            call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	 
    !        else
    !            if (ileft == 2) call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	 
    !        endif
    !        ! update l2r
    !        if (ileft >= 1) then
    !            if (ileft <= nchorizo-1 ) then    
    !                do k=ileft,nchorizo-1
    !                   x=ristra(k)%x
	   !                call v6proplr(x,x,-1,cwtl2r(:,:,k-1),cwtl2r(:,:,k))
    !                enddo
    !                call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))	   
    !            else ! ileft = nchorizo
    !                call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))	   
    !            endif 
    !        else ! ileft == 0
    !            call v6propr(xl,-1,cwtbegin(:,:),cwtl2r(:,:,0))	   
    !                do k=1,nchorizo-1
    !                   x=ristra(k)%x
	   !                call v6proplr(x,x,-1,cwtl2r(:,:,k-1),cwtl2r(:,:,k))
    !                enddo
    !            call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))	     
    !        endif
    !    endif
    !case (1) ! normal, from left to right
    !    icnow=0
    !    icstdmove1bisect1tot=icstdmove1bisect1tot+1
    !    !icmov1bisect1tot=abs(iright-ileft+1)*icstdmove1bisect1tot       
    !    do i=ileft,iright
    !        if (i.eq.0) then    
    !            if (.not.imp) then  
    !                icmov1tot=icmov1tot+1
    !                xold=ristra(i)%x 
    !                gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                xnew(:,:)=xold(:,:)+gauss(:,:)   
	   !             cwtold(:,:)=cwtbegin(:,:)
	   !             tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i)))))
	   !             call corop(xnew,iplo,jplo,cwtnew)
	   !             call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmp(:,:))
	   !             tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
	   !             xd2new(:,:)=(ristra(i+1)%x(:,:)-xnew(:,:))**2
	   !             xd2old(:,:)=(ristra(i+1)%x(:,:)-xold(:,:))**2  
	   !             rt2=gc*(sum(xd2new)-sum(xd2old))	  
	   !             logratio=tmp2-tmp1+rt2                                    
	   !             rn=randn(1)
	   !             if (log(rn(1)).lt.logratio) then ! dt ne 0
    !                    icnow=icnow+1 
    !                    icmov1=icmov1+1	      
	   !                 ristra(i)%x(:,:)=xnew(:,:)
	   !                 cwtr2l(:,:,i)=cwtr2ltmp(:,:)   
    !                ! update cwtbegin
    !                    cwtbegin(:,:)=cwtnew(:,:)  
    !                ! update cwtl2r	   
	   !                 call v6propr(xnew,-1,cwtnew(:,:),cwtl2r(:,:,i)) ! just update the i.
    !                ! prepare for the next cwtl2r
    !                    xoldnxt=ristra(i+1)%x
    !                    call v6proplr(xoldnxt,xoldnxt,-1,cwtl2r(:,:,i),cwtl2r(:,:,i+1))   
    !                    pathcoreval=exp(tmp2)
    !                endif          
    !            else ! imp i=0    
    !                icmov1tot=icmov1tot+1
    !                gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                xold=ristra(i)%x
    !                xnew(:,:)=xold(:,:)+gauss(:,:)  
    !                xold1=2*xnew-xold
    !                xnew1=xold-gauss
    !                xd2new(:,:)=(ristra(i+1)%x(:,:)-xnew(:,:))**2
    !                xd2old(:,:)=(ristra(i+1)%x(:,:)-xold(:,:))**2
    !                xd2new1=(ristra(i+1)%x-xnew1)**2
    !                xd2old1=(ristra(i+1)%x-xold1)**2    
    !                gnew=exp(gc*sum(xd2new))
    !                gold=exp(gc*sum(xd2old))
    !                gnew1=exp(gc*sum(xd2new1))
    !                gold1=exp(gc*sum(xd2old1)) 
    !                call corop(xnew,iplo,jplo,cwtnew)
    !                call corop(xold1,iplo,jplo,cwtold1)
    !                call corop(xnew1,iplo,jplo,cwtnew1)
    !                call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew(:,:))
    !                call v6propl(xold1,-1,cwtr2l(:,:,i+1),cwtr2ltmpold1(:,:))
    !                call v6propl(xnew1,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew1(:,:))               
    !                tmp1=pathcoreval         
    !                tmp2=abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmpnew(:,:))))           
    !                tmp3=abs(real(sum(conjg(cwtold1(:,:))*cwtr2ltmpold1(:,:))))                                         
	   !             tmp4=abs(real(sum(conjg(cwtnew1(:,:))*cwtr2ltmpnew1(:,:))))                         
    !                q1=(gnew*tmp2+gnew1*tmp4)/(gold*tmp1+gold1*tmp3)
    !                q2=gnew*tmp2/(gnew*tmp2+gnew1*tmp4)                                              
	   !             rn2=randn(2)
	   !             if (rn2(1) < q1) then ! dt ne 0
    !                    icnow=icnow+1 
    !                    icmov1=icmov1+1                      
    !                    if (rn2(2) < q2) then
	   !                     ristra(i)%x(:,:)=xnew(:,:)
	   !                     cwtr2l(:,:,i)=cwtr2ltmpnew(:,:)   
    !                        cwtbegin(:,:)=cwtnew(:,:)  
    !                        pathcoreval=tmp2                        
    !                    ! update cwtl2r	   
	   !                     call v6propr(xnew,-1,cwtnew(:,:),cwtl2r(:,:,i)) ! just update the i.
    !                    ! prepare for the next cwtl2r
    !                        xoldnxt=ristra(i+1)%x
    !                        call v6proplr(xoldnxt,xoldnxt,-1,cwtl2r(:,:,i),cwtl2r(:,:,i+1))                               
    !                    else
	   !                     ristra(i)%x=xnew1
	   !                     cwtr2l(:,:,i)=cwtr2ltmpnew1(:,:)   
    !                        cwtbegin=cwtnew1	   
    !                        pathcoreval=tmp4   
    !                    ! update cwtl2r	   
	   !                     call v6propr(xnew1,-1,cwtnew1(:,:),cwtl2r(:,:,i)) ! just update the i.
    !                    ! prepare for the next cwtl2r
    !                        xoldnxt=ristra(i+1)%x
    !                        call v6proplr(xoldnxt,xoldnxt,-1,cwtl2r(:,:,i),cwtl2r(:,:,i+1))  
    !                    endif
    !                endif                        
    !            endif  
	   !     else
	   !         if (i.eq.nchorizo) then   
    !                if (.not.imp) then
    !                    icmov1tot=icmov1tot+1
    !                    xold=ristra(i)%x 
    !                    gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                    xnew(:,:)=xold(:,:)+gauss(:,:)                       
	   !                 cwtold(:,:)=cwtend(:,:)
	   !                 tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtl2r(:,:,i))*cwtold(:,:)))))
	   !                 call corop(xnew,ipro,jpro,cwtnew)
	   !                 call v6propl(xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmp(:,:))
	   !                 tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtnew(:,:)))))
	   !                 xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2
	   !                 xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2
	   !                 rt2=gc*(sum(xd2new)-sum(xd2old))	  
	   !                 logratio=tmp2-tmp1+rt2                         
	   !                 rn=randn(1)
	   !                 if (log(rn(1)).lt.logratio) then
    !                        icnow=icnow+1
    !                        icmov1=icmov1+1	     
	   !                     ristra(i)%x(:,:)=xnew(:,:)
	   !                     cwtl2r(:,:,i)=cwtl2rtmp(:,:)
    !                    ! update cwtend
    !                        cwtend(:,:)=cwtnew(:,:)                  
    !                        pathcoreval=exp(tmp2)
    !                    endif                        
    !                else ! imp                
    !                    icmov1tot=icmov1tot+1
    !                    gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                    xold=ristra(i)%x 
    !                    xnew(:,:)=xold(:,:)+gauss(:,:)  
    !                    xold1=2*xnew-xold
    !                    xnew1=xold-gauss
    !                    xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2
    !                    xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2
    !                    xd2new1=(ristra(i-1)%x-xnew1)**2
    !                    xd2old1=(ristra(i-1)%x-xold1)**2    
    !                    gnew=exp(gc*sum(xd2new))
    !                    gold=exp(gc*sum(xd2old))
    !                    gnew1=exp(gc*sum(xd2new1))
    !                    gold1=exp(gc*sum(xd2old1)) 
    !                    call corop(xnew,ipro,jpro,cwtnew)
    !                    call corop(xold1,ipro,jpro,cwtold1)
    !                    call corop(xnew1,ipro,jpro,cwtnew1)
    !                    call v6propl(xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmpnew(:,:))
    !                    call v6propl(xold1,-1,cwtl2r(:,:,i-1),cwtl2rtmpold1(:,:))
    !                    call v6propl(xnew1,-1,cwtl2r(:,:,i-1),cwtl2rtmpnew1(:,:))                         
    !                    tmp1=pathcoreval               
	   !                 tmp2=abs(real(sum(conjg(cwtl2rtmpnew(:,:))*cwtnew(:,:))))                        
	   !                 tmp3=abs(real(sum(conjg(cwtl2rtmpold1(:,:))*cwtold1(:,:))))                                              
	   !                 tmp4=abs(real(sum(conjg(cwtl2rtmpnew1(:,:))*cwtnew1(:,:))))                                             
    !                    q1=(gnew*tmp2+gnew1*tmp4)/(gold*tmp1+gold1*tmp3)
    !                    q2=gnew*tmp2/(gnew*tmp2+gnew1*tmp4)                                                                   
	   !                 rn2=randn(2)
	   !                 if (rn2(1) < q1) then ! dt ne 0
    !                        icnow=icnow+1 
    !                        icmov1=icmov1+1                         
    !                        if (rn2(2) < q2) then
	   !                         ristra(i)%x(:,:)=xnew(:,:)
	   !                         cwtl2r(:,:,i)=cwtl2rtmpnew(:,:)
    !                            cwtend(:,:)=cwtnew(:,:)	   
    !                            pathcoreval=tmp2                             
    !                        else
	   !                         ristra(i)%x=xnew1
	   !                         cwtl2r(:,:,i)=cwtl2rtmpnew1(:,:)
    !                            cwtend=cwtnew1	   
    !                            pathcoreval=tmp4                                   
    !                        endif 
    !                    endif                                
    !                endif
    !            else            
    !                if (.not.imp) then
    !                    icmov1tot=icmov1tot+1
    !                    xold=ristra(i)%x 
    !                    gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                    xnew(:,:)=xold(:,:)+gauss(:,:)   
    !                    cwtold(:,:)=cwtl2r(:,:,i)
	   !                 call v6proplr(xnew,xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmp(:,:)) 
	   !                 tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtr2l(:,:,i+1)))))
	   !                 tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i+1)))))
	   !                 xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2+(ristra(i+1)%x(:,:)-xnew(:,:))**2
	   !                 xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2+(ristra(i+1)%x(:,:)-xold(:,:))**2	   
	   !                 rt2=gc*(sum(xd2new)-sum(xd2old))	
	   !                 logratio=tmp2-tmp1+rt2	                                                
	   !                 rn=randn(1)
	   !                 if (log(rn(1)).lt.logratio) then
    !                            !write(6,*) 'chorizo i updated!'
    !                        icnow=icnow+1           
    !                        icmov1=icmov1+1	  	
	   !                     ristra(i)%x(:,:)=xnew(:,:)
	   !                     cwtl2r(:,:,i)=cwtl2rtmp(:,:)	    
    !                        pathcoreval=exp(tmp2)
    !                    endif       
    !                else ! imp                                 
    !                    icmov1tot=icmov1tot+1
    !                    gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                    xold=ristra(i)%x 
    !                    xnew(:,:)=xold(:,:)+gauss(:,:)  
    !                    xold1=2*xnew-xold
    !                    xnew1=xold-gauss
    !                    xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2+(ristra(i+1)%x(:,:)-xnew(:,:))**2
    !                    xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2+(ristra(i+1)%x(:,:)-xold(:,:))**2
    !                    xd2new1=(ristra(i-1)%x-xnew1)**2+(ristra(i+1)%x-xnew1)**2
    !                    xd2old1=(ristra(i-1)%x-xold1)**2+(ristra(i+1)%x-xold1)**2    
    !                    gnew=exp(gc*sum(xd2new))
    !                    gold=exp(gc*sum(xd2old))
    !                    gnew1=exp(gc*sum(xd2new1))
    !                    gold1=exp(gc*sum(xd2old1)) 
    !                    call v6proplr(xnew,xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmpnew(:,:)) 
    !                    call v6proplr(xold1,xold1,-1,cwtl2r(:,:,i-1),cwtl2rtmpold1(:,:)) 
    !                    call v6proplr(xnew1,xnew1,-1,cwtl2r(:,:,i-1),cwtl2rtmpnew1(:,:))                             
    !                    tmp1=pathcoreval                      
	   !                 tmp2=abs(real(sum(conjg(cwtl2rtmpnew(:,:))*cwtr2l(:,:,i+1))))                      
	   !                 tmp3=abs(real(sum(conjg(cwtl2rtmpold1(:,:))*cwtr2l(:,:,i+1))))                                            
	   !                 tmp4=abs(real(sum(conjg(cwtl2rtmpnew1(:,:))*cwtr2l(:,:,i+1))))                                                   
    !                    q1=(gnew*tmp2+gnew1*tmp4)/(gold*tmp1+gold1*tmp3)
    !                    q2=gnew*tmp2/(gnew*tmp2+gnew1*tmp4)                                                                  
	   !                 rn2=randn(2)
	   !                 if (rn2(1) < q1) then 
    !                        icnow=icnow+1 
    !                        icmov1=icmov1+1                    
    !                        if (rn2(2) < q2) then
	   !                         ristra(i)%x(:,:)=xnew(:,:)
	   !                         cwtl2r(:,:,i)=cwtl2rtmpnew(:,:)   
    !                            pathcoreval=tmp2                             
    !                        else
	   !                         ristra(i)%x=xnew1
	   !                         cwtl2r(:,:,i)=cwtl2rtmpnew1(:,:)	   
    !                            pathcoreval=tmp4                                   
    !                        endif 
    !                    endif                             
    !                endif                
    !                if ( icnow /= 0 ) then
    !            ! update l2r
    !                    if (i.le.(nchorizo-2)) then
		  !              x=ristra(i+1)%x
		  !              call v6proplr(x,x,-1,cwtl2r(:,:,i),cwtl2r(:,:,i+1))
    !                    else ! i=nchorizo-1
    !                    xr=ristra(i+1)%x
    !                    call v6propl(xr,-1,cwtl2r(:,:,i),cwtl2r(:,:,i+1))	    
    !                    endif               
    !                endif   
	   !         endif 
	   !     endif 	 
    !    enddo 
    !    if (icnow /= 0) then 
    !        xl=ristra(0)%x
    !        xr=ristra(nchorizo)%x
    !        ! update cwtl2r    
    !        if (iright <= nchorizo-3) then 
    !            do k=iright,nchorizo-3
    !                x=ristra(k+2)%x
    !                call v6proplr(x,x,-1,cwtl2r(:,:,k+1),cwtl2r(:,:,k+2))   
    !            enddo
    !            call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))	      
    !        else    
    !            if (iright == nchorizo-2) then   
    !                call v6propl(xr,-1,cwtl2r(:,:,iright+1),cwtl2r(:,:,nchorizo))	          
    !            endif  
    !        endif    
    !        ! update cwtr2l from iright th bead until the 0th.
    !        if (iright == nchorizo) then    
    !            !x=ristra(iright)%x   
	   !         call v6propr(xr,-1,cwtend(:,:),cwtr2l(:,:,iright))
	   !         do k=iright-1,1,-1
	   !         x=ristra(k)%x
	   !         call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
	   !         enddo
	   !         !xl=ristra(0)%x
	   !         call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	    
    !        else      
    !            if (iright >= 1) then  
    !                do k=iright,1,-1
	   !                x=ristra(k)%x
	   !                call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
	   !             enddo
	   !             !xl=ristra(0)%x
	   !             call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	   
    !            else ! iright = 0
    !                call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	   
    !            endif  
    !        endif
    !    endif   
    !case default 
    !    write (6,*) 'arrow for stdmove1bisect1flex has not been set, abort!', arrow
	   ! call abort
    !end select
    !return
    !end subroutine stdmove1fastv6flexa              
    !
    
    
    !subroutine bisectv6b(ileftbisect) 
    !! based on v6a version. move beads first, then do an overall judge.
    !use wavefunction
    !use estimator
    !use random
    !use v6stepcalc
    !integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
				!	    ,k,bilvl &
	   !                 ,dtexpt
    !real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    !real(kind=r8) :: xmid(3,npart),xr(3,npart),xl(3,npart)
    !real(kind=r8) :: x(3,npart),xnew(3,npart)
    !real(kind=r8) :: tmp 
    !real(kind=r8) :: oldlv(0:mmax),newlv(0:mmax)
    !complex(kind=r8) :: cwtl2rtmp1(0:nspin-1,nisospin,0:nbisect),cwtr2ltmp1(0:nspin-1,nisospin,0:nbisect) &
	   !                 ,cwtl2rtmp2(0:nspin-1,nisospin,0:nbisect),cwtr2ltmp2(0:nspin-1,nisospin,0:nbisect) &
	   !                 ,cwtbegintmp(0:nspin-1,nisospin),cwtendtmp(0:nspin-1,nisospin)
    !logical :: reject	
    !
    !ibisecttot=ibisecttot+1  
    !
    !!do i=0,nbisect
	   ! !oldristra(i)=ristra(ileftbisect+i)
    !!enddo
    !
    !oldristra(0:nbisect)=ristra(ileftbisect:ileftbisect+nbisect)
    !
    !newristra(0)=oldristra(0)
    !newristra(nbisect)=oldristra(nbisect)
    !
    !
    !!do i=0,nbisect
    !!xoldtmp(i,:,:)=oldristra(i)%x
    !!enddo
    ! 
    !cwtbegintmp(:,:)=cwtl2r(:,:,ileftbisect)
    !cwtendtmp(:,:)=cwtr2l(:,:,ileftbisect+nbisect)  
    !
    !! tmp1 deal with oldlv
    !cwtl2rtmp1(:,:,0:nbisect)=cwtl2r(:,:,ileftbisect:ileftbisect+nbisect)
    !cwtr2ltmp1(:,:,0:nbisect)=cwtr2l(:,:,ileftbisect:ileftbisect+nbisect)
    !! tmp2 deal with newlv   
    !cwtl2rtmp2(:,:,0:nbisect)=cwtl2r(:,:,ileftbisect:ileftbisect+nbisect)
    !cwtr2ltmp2(:,:,0:nbisect)=cwtr2l(:,:,ileftbisect:ileftbisect+nbisect)
    !
    !oldlv=1.
    !newlv=1.   
    !
    !tau=dt*nbisect ! total time for the bisection slice.
    !
    !! up to here, newristra has been loaded.
    !
    !! judge if we accept the newristra or not.
    !! in this v6 case, v is not justpotential v, it should be the value of < ... >. lv.
    !    
    !do bilvl=1,mmax  ! besection.  level 1 to N. mmax = the total level N.	 
    !    ! move the current level   
    !    sigmamid=sqrt(hbar*tau/2**bilvl)  
	   ! !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	   ! do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		  !  imid=2**(mmax-bilvl)+i*2**(mmax-bilvl+1)
		  !  il=imid-2**(mmax-bilvl)
		  !  ir=imid+2**(mmax-bilvl)
		  !  gauss(:,:)=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		  !  xmid=(newristra(il)%x+newristra(ir)%x)/2
		  !  newristra(imid)%x(:,:)=xmid(:,:)+gauss(:,:)  
    !        !write(12,*) 'bilvl il ir imid sigmamid =',bilvl, il,ir,imid,sigmamid
    !        !write(12,*) 'gauss =', gauss
    !    enddo         
    !enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    !! then judge the last level
	   ! bisecttot(mmax)=bisecttot(mmax)+1	   
	   ! lmax=2**mmax-1 
	   ! jd=1 ! interval
	   ! dtexpt=-1 ! for the R L part, that's why the extra -1. dt=2**dtexp
	   ! do l=1,lmax
		  !  j=l*jd
		  !  jp=(l-1)*jd
		  !  !xold=oldristra(j)%x
    !        xnew(:,:)=newristra(j)%x(:,:)
		  !  ! call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(:,:,jp),cwtl2rtmp1(:,:,j))    ! should be lr, not rl I believe.
		  !  call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(:,:,jp),cwtl2rtmp2(:,:,j))  	   
    !    enddo
	   ! jmax=lmax*jd
	   ! newlv(mmax)=abs(real( sum(conjg(cwtl2rtmp2(:,:,jmax))*cwtendtmp(:,:)) )) 	
	   ! oldlv(mmax)=pathcoreval !abs(real( sum(conjg(cwtl2rtmp1(:,:,jmax))*cwtendtmp(:,:)) ))	   
	   ! tmp=log(newlv(mmax))-log(oldlv(mmax)) 
  	 !   !tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
	   ! rn=randn(1)	   
    !! can check if l2r and r2l give	the same lv.		    
    !        !write(12,*) 'tmp ln rn=',tmp,log(rn(1))
	   ! if (log(rn(1)).lt.tmp) then
	   !     reject=.false. 
	   !     bisectcount(mmax)=bisectcount(mmax)+1	 
	   ! else
	   !     reject=.true. 
    !        !write(6,*) 'rn=',rn(1),log(rn(1))
    !        !write(6,*) 'tmp=',tmp,newlv(mmax),oldlv(mmax),log(newlv(mmax)),log(oldlv(mmax)) 	
    !        
    !    endif 
		  !  !write(6,*) 'imid=',imid	
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!        
    !    
    !if (reject) then
	   ! call addval(2,0.0_r8,1.0_r8)
    !else
	   ! ibisect=ibisect+1 
	   ! call addval(2,1.0_r8,1.0_r8)  
    !
	   ! !do i=0,nbisect
    !!       ristra(ileftbisect+i)=newristra(i)
    !!   enddo	 
    !    ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good. 04/10/2019�� remove pointer, problem fixed I tihnk.	 
    !
    !    xl=ristra(0)%x     
	   ! xr=ristra(nchorizo)%x     
    !    
    !    pathcoreval=newlv(mmax)
    !! update cwtl2r, cwtr2l.	 
    !! update cwtl2r for all beads. this can be optimized in the future.
    ! 
	   ! cwtl2r(:,:,ileftbisect+1:ileftbisect+nbisect-1)=cwtl2rtmp2(:,:,1:jmax) ! here jmax should already be 2**mmax-1 
	   ! if ((ileftbisect+nbisect).le.nchorizo-1) then
    !    do k=ileftbisect+nbisect,nchorizo-1
	   ! x=ristra(k)%x
	   ! call v6proplr(x,x,-1,cwtl2r(:,:,k-1),cwtl2r(:,:,k))
    !    enddo
	   ! call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))
    !    else
    !    !ileftbisect+nbisect=nchorizo  
	   ! call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))   
    !    endif
    !
    !! update cwtr2l for all the beads
	   ! do k=ileftbisect+nbisect-1,1,-1
	   ! x=ristra(k)%x
	   ! call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
	   ! enddo
    !    call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))  	
    !
    !endif  
    !
    !!rejectbisectv6improve=reject
    !
    !return  
    !end subroutine bisectv6b            

!    subroutine bisectv6lb 
!    ! b version, move the whole nbiect beads instead of bisecthalf.
!    ! should be better than bisectv6l,rhalf     
!    ! move the left end, 0th bead,ileftbiect=0; if ileftbisect /= 0 then the small modification on 0th level is needed.
!    use wavefunction
!    use estimator
!    use random
!    use v6stepcalc
!    use brut
!    integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
!					    ,irightbisect,k,bilvl &
!	                    ,dtexpt &
!                        ,nbisectnow,mmaxnow
!    real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
!    real(kind=r8) :: xmid(3,npart),xr(3,npart),xl(3,npart)
!    real(kind=r8) :: x(3,npart),xnew(3,npart),x0new(3,npart)
!    real(kind=r8) :: tmp 
!    real(kind=r8) :: oldlv(0:mmax),newlv(0:mmax)
!    complex(kind=r8) :: cwtl2rtmp1(0:nspin-1,nisospin,0:nbisect) &
!	                    ,cwtl2rtmp2(0:nspin-1,nisospin,0:nbisect) &
!	                    ,cwtendtmp(0:nspin-1,nisospin) &
!                        ,cwtbegintmp2(0:nspin-1,nisospin)        
!    logical :: reject	
!
!    ibisecttotl=ibisecttotl+1  
!   
!    ileftbisect=0
!    irightbisect=ileftbisect+nbisect
!   
!    !nbisectnow=irightbisect-ileftbisect ! has to be 2^n.
!    !mmaxnow=log(real(nbisectnow))/log(real(2)) ! mmaxnow = mmax-1
!   
!    nbisectnow=nbisect
!    mmaxnow=mmax
!   
!    tau=dt*nbisectnow ! total time for the bisection slice. nisectnow is usually nbisecthalf
!    sigmamid=sqrt(2*hbar*tau)  
!
!    oldristra(0:nbisectnow)=ristra(ileftbisect:ileftbisect+nbisectnow)
! 
!    ! deal with the two ends first, then do the classic bisection.   
!   
!    newristra(nbisectnow)=oldristra(nbisectnow)   
!    cwtendtmp(:,:)=cwtr2l(:,:,ileftbisect+nbisectnow)
! 
!    newristra(0)%x=oldristra(nbisectnow)%x+sigmamid*reshape(gaussian(3*npart),(/3,npart/))
!   
!    x0new(:,:)=newristra(0)%x(:,:)
!   
!    ! i,jpro,i,jplo, the o means old are actually the current order.
!   
!    call corop(x0new,iplo,jplo,cwtbegintmp2)     
!   
!    !oldlv(0)=1.
!    !newlv(0)=1. ! dummy.
!
!    cwtl2rtmp1(:,:,0:nbisect)=cwtl2r(:,:,ileftbisect:ileftbisect+nbisect)
!    call v6propr(x0new,-1,cwtbegintmp2,cwtl2rtmp2(:,:,0)) 
!   
!    ! prepare the further new bisection positions
!    do bilvl=1,mmaxnow   ! level 1 to N. mmax = the total level N.
!	    sigmamid=sqrt(hbar*tau*2.0_r8**(-bilvl))  
!	    !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
!	    do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
!		    imid=2**(mmaxnow-bilvl)+i*2**(mmaxnow-bilvl+1)
!		    il=imid-2**(mmaxnow-bilvl)
!		    ir=imid+2**(mmaxnow-bilvl)
!		    !write(6,*) 'bilvl il ir imid sigmamid=',bilvl, il,ir,imid,sigmamid
!		    gauss(:,:)=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
!		    xmid(:,:)=(newristra(il)%x(:,:)+newristra(ir)%x(:,:))/2
!		    newristra(imid)%x(:,:)=xmid(:,:)+gauss(:,:)
!	    enddo 
!    enddo 
!
!    bilvl=mmaxnow  ! besection.  level 1 to N. mmax = the total level N.	   
!	bisecttotl(bilvl)=bisecttotl(bilvl)+1.	   
!	lmax=2**bilvl-1 
!	jd=2**(mmaxnow-bilvl) ! interval
!	dtexpt=mmaxnow-bilvl-1 ! for the R L part, that's why the extra -1. dt=2**dtexp 
!	do l=1,lmax
!		j=l*jd
!		jp=(l-1)*jd
!		!xold=oldristra(j)%x
!        xnew(:,:)=newristra(j)%x(:,:)
!		!call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(:,:,jp),cwtl2rtmp1(:,:,j))    ! should be lr, not rl I believe.
!		call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(:,:,jp),cwtl2rtmp2(:,:,j))  	   
!		!if ( sum((xnew-xold)**2).eq.0) then
!			!write(6,*) 'bilvl=',bilvl	
!			!write(6,*) 'l=',l,j,jp
!			!write(6,*) 'sigmamid=',sigmamid
!		!   write(6,*) 'xnew-xold=',xnew-xold
!		!endif   
!	enddo
!	jmax=lmax*jd  
!	oldlv(bilvl)=pathcoreval !abs(real( sum(conjg(cwtl2rtmp1(:,:,jmax))*cwtendtmp(:,:)) ))	
!	newlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp2(:,:,jmax))*cwtendtmp(:,:)) )) 
!	tmp=log(newlv(bilvl))-log(oldlv(bilvl))
!  	!tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
!	rn=randn(1)	   
!! can check if l2r and r2l give	the same lv.		    
!	    !write(6,*) 'tmp ln rn=',tmp,log(rn(1))
!	if (log(rn(1)).lt.tmp) then
!	    reject=.false. 
!        bisectcountl(bilvl)=bisectcountl(bilvl)+1
!        call addval(10,1.0_r8,1.0_r8)  
!	    ibisectl=ibisectl+1 
!	    do i=0,nbisectnow
!            ristra(ileftbisect+i)=newristra(i)
!        enddo	 
!    !ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good.	 
!        xl=ristra(0)%x     
!	    xr=ristra(nchorizo)%x 
!    ! update cwtl2r, cwtr2l.	 
!        cwtbegin(:,:)=cwtbegintmp2(:,:)
!	    cwtl2r(:,:,ileftbisect:ileftbisect+nbisectnow-1)=cwtl2rtmp2(:,:,0:jmax) ! here jmax should already be 2**mmax-1    
!	    if ((ileftbisect+nbisectnow).le.nchorizo-1) then
!            do k=ileftbisect+nbisectnow,nchorizo-1
!	        x=ristra(k)%x
!	        call v6proplr(x,x,-1,cwtl2r(:,:,k-1),cwtl2r(:,:,k))
!            enddo
!	        call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))
!        else !ileftbisect+nbisect=nchorizo  
!	        call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))   
!        endif       
!    ! update cwtr2l for all the beads
!        if ((ileftbisect+nbisectnow-1)>=1) then
!	        do k=ileftbisect+nbisectnow-1,1,-1
!	           x=ristra(k)%x
!	           call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
!	        enddo
!            call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0)) 
!        else   
!            call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0)) 
!        endif
! 
!        pathcoreval=newlv(mmaxnow)        
!        
!	else
!	    reject=.true.  
!        call addval(10,0.0_r8,1.0_r8)   
!	endif   
!
!    return  
!    end subroutine bisectv6lb  
!     

!    subroutine bisectv6rb ! move the right end. irightbisect=nchorizo
!    ! move beads, then do an overall judge.
!    use wavefunction
!    use estimator
!    use random
!    use v6stepcalc
!    use brut
!    integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
!					    ,irightbisect,k,bilvl &
!	                    ,dtexpt &
!                        ,nbisectnow,mmaxnow
!    real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
!    real(kind=r8) :: xmid(3,npart),xl(3,npart)
!    real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart) 
!    real(kind=r8) :: tmp 
!    real(kind=r8) :: oldlv(0:mmax),newlv(0:mmax)
!    complex(kind=r8) :: cwtl2rtmp1(0:nspin-1,nisospin,0:nbisect),cwtl2rtmp2(0:nspin-1,nisospin,0:nbisect) &
!	                    ,cwtbegintmp(0:nspin-1,nisospin) &
!                        ,cwtendtmp1now(0:nspin-1,nisospin),cwtendtmp2now(0:nspin-1,nisospin)   	                               
!    logical :: reject	
!   
!    ibisecttotr=ibisecttotr+1    
!   
!    ileftbisect=nchorizo-nbisect
!    irightbisect=nchorizo
!   
!    !nbisectnow=irightbisect-ileftbisect
!    !mmaxnow=log(real(nbisectnow))/log(real(2))
! 
!    nbisectnow=nbisect
!    mmaxnow=mmax  
!  
!    tau=dt*nbisectnow ! total time for the bisection slice.
!    sigmamid=sqrt(2*hbar*tau)  
!   
!    oldristra(0:nbisectnow)=ristra(ileftbisect:ileftbisect+nbisectnow)
! 
!    ! deal with the two ends first, then do the classic bisection.   
!    newristra(0)=oldristra(0) 
!    cwtbegintmp(:,:)=cwtl2r(:,:,ileftbisect)
!   
!    newristra(nbisectnow)%x=oldristra(0)%x+sigmamid*reshape(gaussian(3*npart),(/3,npart/))
!   
!    xold(:,:)=oldristra(nbisectnow)%x(:,:)
!    xnew(:,:)=newristra(nbisectnow)%x(:,:)    
!   
!    call corop(xnew,ipro,jpro,cwtendnew)     
!
!    !oldlv(0)=1
!    !newlv(0)=1 ! dummy, no use.
!   
!    cwtl2rtmp1(:,:,0:nbisect)=cwtl2r(:,:,ileftbisect:ileftbisect+nbisect)
!    cwtl2rtmp2(:,:,0)=cwtl2rtmp1(:,:,0)
!    cwtendtmp1now=cwtr2l(:,:,nchorizo)
!    !call v6propr(oldristra(nbisectnow)%x,-1,cwtend,cwtendtmp1now)  
!    call v6propr(newristra(nbisectnow)%x,-1,cwtendnew,cwtendtmp2now)   ! r2l  
!     
!    do bilvl=1,mmaxnow   ! level 1 to N. mmax = the total level N.
!	    sigmamid=sqrt(hbar*tau/2**bilvl)  
!	    !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
!	    do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
!		    imid=2**(mmaxnow-bilvl)+i*2**(mmaxnow-bilvl+1)
!		    il=imid-2**(mmaxnow-bilvl)
!		    ir=imid+2**(mmaxnow-bilvl)
!		    !write(6,*) 'bilvl il ir imid sigmamid=',bilvl, il,ir,imid,sigmamid
!		    gauss(:,:)=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
!		    xmid=(newristra(il)%x+newristra(ir)%x)/2
!		    newristra(imid)%x(:,:)=xmid(:,:)+gauss(:,:)
!	    enddo 
!    enddo     
!    ! judge
!    bilvl=mmaxnow  ! besection.  level 1 to N. mmax = the total level N.	   
!	bisecttotr(bilvl)=bisecttotr(bilvl)+1   
!	lmax=2**bilvl-1 
!	jd=2**(mmaxnow-bilvl) ! interval
!	dtexpt=mmaxnow-bilvl-1 ! for the R L part, that's why the extra -1. dt=2**dtexp         
!	do l=1,lmax
!		j=l*jd
!		jp=(l-1)*jd
!		!xold=oldristra(j)%x
!        xnew(:,:)=newristra(j)%x(:,:)
!		!call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(:,:,jp),cwtl2rtmp1(:,:,j))    ! should be lr, not rl I believe.
!		call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(:,:,jp),cwtl2rtmp2(:,:,j))  	   
!		!if ( sum((xnew-xold)**2).eq.0) then
!			!write(6,*) 'bilvl=',bilvl	
!			!write(6,*) 'l=',l,j,jp
!			!write(6,*) 'sigmamid=',sigmamid
!		!   write(6,*) 'xnew-xold=',xnew-xold
!		!endif   
!	enddo
!	jmax=lmax*jd  
!	oldlv(bilvl)=pathcoreval !abs(real( sum(conjg(cwtl2rtmp1(:,:,jmax))*cwtendtmp1now(:,:)) ))	
!	newlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp2(:,:,jmax))*cwtendtmp2now(:,:)) )) 
!	tmp=log(newlv(bilvl))-log(oldlv(bilvl))
!  	!tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
!	rn=randn(1)	   
!! can check if l2r and r2l give	the same lv.		    
!	    !write(6,*) 'tmp ln rn=',tmp,log(rn(1))
!	if (log(rn(1)).lt.tmp) then
!	    reject=.false. 
!        bisectcountr(bilvl)=bisectcountr(bilvl)+1 
!        call addval(11,1.0_r8,1.0_r8)  
!	    ibisectr=ibisectr+1 
!        ristra(ileftbisect:ileftbisect+nbisectnow)=newristra(0:nbisectnow)     
!	 	 
!    ! update cwtl2r, cwtr2l.	 
!    ! update cwtl2r for all beads. this can be optimized in the future.
!        cwtl2r(:,:,ileftbisect+1:ileftbisect+nbisectnow-1)=cwtl2rtmp2(:,:,1:jmax) ! here jmax should already be 2**mmax-1  
!  
!        call v6propl(ristra(nchorizo)%x,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo)) ! ileftbisect+nbisectnow=nchorizo  
!
!        ! update cwtr2l for all the beads
!        cwtend(:,:)=cwtendnew(:,:)         
!        cwtr2l(:,:,nchorizo)=cwtendtmp2now
!	    do k=nchorizo-1,1,-1
!            x=ristra(k)%x
!	        call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
!	    enddo
!        xl=ristra(0)%x
!        call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))  	
!
!        pathcoreval=newlv(mmaxnow)
!        
!	else
!	    reject=.true.  
!        call addval(11,0.0_r8,1.0_r8)  
!	endif 
!            
!    return  
!    end subroutine bisectv6rb      
!  


    !subroutine bisectv6a(ileftbisect) 
    !! slighly faster than the bisectv6improve verison. only move the beads at the currently level then judge.
    !use wavefunction
    !use estimator
    !use random
    !use v6stepcalc
    !integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
				!	    ,k,bilvl &
	   !                 ,dtexpt
    !real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    !real(kind=r8) :: xmid(3,npart),xr(3,npart),xl(3,npart)
    !real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart)
    !real(kind=r8) :: tmp
    !real(kind=r8) :: oldlv(0:mmax),newlv(0:mmax)
    !complex(kind=r8) :: cwtl2rtmp1(0:nspin-1,nisospin,0:nbisect),cwtr2ltmp1(0:nspin-1,nisospin,0:nbisect) &
	   !                 ,cwtl2rtmp2(0:nspin-1,nisospin,0:nbisect),cwtr2ltmp2(0:nspin-1,nisospin,0:nbisect) &
	   !                 ,cwtbegintmp(0:nspin-1,nisospin),cwtendtmp(0:nspin-1,nisospin) 
    !logical :: reject	
    !
    !ibisecttot=ibisecttot+1  
    !
    !!do i=0,nbisect
	   ! !oldristra(i)=ristra(ileftbisect+i)
    !!enddo
    !
    !oldristra(0:nbisect)=ristra(ileftbisect:ileftbisect+nbisect)
    !
    !newristra(0)=oldristra(0)
    !newristra(nbisect)=oldristra(nbisect)
    !
    !
    !!do i=0,nbisect
    !!xoldtmp(i,:,:)=oldristra(i)%x
    !!enddo
    ! 
    !cwtbegintmp(:,:)=cwtl2r(:,:,ileftbisect)
    !cwtendtmp(:,:)=cwtr2l(:,:,ileftbisect+nbisect)  
    !
    !! tmp1 deal with oldlv
    !cwtl2rtmp1(:,:,0:nbisect)=cwtl2r(:,:,ileftbisect:ileftbisect+nbisect)
    !cwtr2ltmp1(:,:,0:nbisect)=cwtr2l(:,:,ileftbisect:ileftbisect+nbisect)
    !! tmp2 deal with newlv   
    !cwtl2rtmp2(:,:,0:nbisect)=cwtl2r(:,:,ileftbisect:ileftbisect+nbisect)
    !cwtr2ltmp2(:,:,0:nbisect)=cwtr2l(:,:,ileftbisect:ileftbisect+nbisect)
    !
    !oldlv=1.
    !newlv=1.   
    !
    !tau=dt*nbisect ! total time for the bisection slice.
    !
    !! up to here, newristra has been loaded.
    !
    !! judge if we accept the newristra or not.
    !! in this v6 case, v is not justpotential v, it should be the value of < ... >. lv.
    !    
    !all: do bilvl=1,mmax  ! besection.  level 1 to N. mmax = the total level N.	 
    !    ! move the current level   
    !    sigmamid=sqrt(hbar*tau*2.0_r8**(-bilvl))  
	   ! !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	   ! do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		  !  imid=2**(mmax-bilvl)+i*2**(mmax-bilvl+1)
		  !  il=imid-2**(mmax-bilvl)
		  !  ir=imid+2**(mmax-bilvl)
		  !  gauss(:,:)=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		  !  xmid(:,:)=(newristra(il)%x(:,:)+newristra(ir)%x(:,:))/2
		  !  newristra(imid)%x(:,:)=xmid(:,:)+gauss(:,:)
    !       
    !        !write(12,*) 'bilvl il ir imid sigmamid =',bilvl, il,ir,imid,sigmamid
    !        !write(12,*) 'gauss =', gauss
    !    enddo         
    !    ! then judge  
	   ! bisecttot(bilvl)=bisecttot(bilvl)+1.	   
	   ! lmax=2**bilvl-1 
	   ! jd=2**(mmax-bilvl) ! interval
	   ! dtexpt=mmax-bilvl-1 ! for the R L part, that's why the extra -1. dt=2**dtexp
	   ! do l=1,lmax
		  !  j=l*jd
		  !  jp=(l-1)*jd
		  !  xold=oldristra(j)%x
    !        xnew(:,:)=newristra(j)%x(:,:)
		  !  call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(:,:,jp),cwtl2rtmp1(:,:,j))    ! should be lr, not rl I believe.
    !! in the last bisection level, no need to calculate this because it old is already there. check.		   
    !        call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(:,:,jp),cwtl2rtmp2(:,:,j))  	   
    !
	   ! enddo
	   ! jmax=lmax*jd
	   ! newlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp2(:,:,jmax))*cwtendtmp(:,:)) )) 	   
	   ! oldlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp1(:,:,jmax))*cwtendtmp(:,:)) ))
    !
	   ! tmp=log(newlv(bilvl))-log(oldlv(bilvl))+log(oldlv(bilvl-1))-log(newlv(bilvl-1))   
  	 !   !tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
	   ! rn=randn(1)	   
    !! can check if l2r and r2l give	the same lv.		    
    !        !write(12,*) 'tmp ln rn=',tmp,log(rn(1))
	   ! if (log(rn(1)).lt.tmp) then
	   !     reject=.false. 
	   !     bisectcount(bilvl)=bisectcount(bilvl)+1.	 
	   ! else
	   !     reject=.true.  
	   !     exit all 
	   ! endif 
		  !  !write(6,*) 'imid=',imid	   
    !enddo all
    !if (reject) then
	   ! call addval(2,0.0_r8,1.0_r8)
    !else
	   ! ibisect=ibisect+1 
	   ! call addval(2,1.0_r8,1.0_r8)  
    !
	   ! !do i=0,nbisect
    !!       ristra(ileftbisect+i)=newristra(i)
    !!   enddo	 
    !    ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good. 04/10/2019�� remove pointer, problem fixed I tihnk.	 
    !
    !    xl=ristra(0)%x     
	   ! xr=ristra(nchorizo)%x     
    !! update cwtl2r, cwtr2l.	 
    !! update cwtl2r for all beads. this can be optimized in the future.
    ! 
	   ! cwtl2r(:,:,ileftbisect+1:ileftbisect+nbisect-1)=cwtl2rtmp2(:,:,1:jmax) ! here jmax should already be 2**mmax-1 
	   ! if ((ileftbisect+nbisect).le.nchorizo-1) then
    !    do k=ileftbisect+nbisect,nchorizo-1
	   ! x=ristra(k)%x
	   ! call v6proplr(x,x,-1,cwtl2r(:,:,k-1),cwtl2r(:,:,k))
    !    enddo
	   ! call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))
    !    else
    !    !ileftbisect+nbisect=nchorizo  
	   ! call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))   
    !    endif
    !
    !! update cwtr2l for all the beads
	   ! do k=ileftbisect+nbisect-1,1,-1
	   ! x=ristra(k)%x
	   ! call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
	   ! enddo
    !    call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))  	
    !
    !    pathcoreval=newlv(mmax)
    !endif  
    !
    !rejectbisectv6improve=reject
    !return  
    !end subroutine bisectv6a         
    !
   
    !subroutine bisectv6improve(ileftbisect) ! for v6 interaction, nbisect beads
    !use wavefunction
    !use estimator
    !use random
    !use v6stepcalc
    !integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
				!	    ,irightbisect,k,bilvl &
	   !                 ,dtexpt
    !real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    !real(kind=r8) :: xmid(3,npart),xr(3,npart),xl(3,npart)
    !real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart)
    !real(kind=r8) :: tmp 
    !real(kind=r8) :: oldlv(0:mmax),newlv(0:mmax)
    !complex(kind=r8) :: cwtl2rtmp1(0:nspin-1,nisospin,0:nbisect),cwtr2ltmp1(0:nspin-1,nisospin,0:nbisect) &
	   !                 ,cwtl2rtmp2(0:nspin-1,nisospin,0:nbisect),cwtr2ltmp2(0:nspin-1,nisospin,0:nbisect) &
	   !                 ,cwtbegintmp(0:nspin-1,nisospin),cwtendtmp(0:nspin-1,nisospin)
    !logical :: reject	
    !  
    !irightbisect=ileftbisect+nbisect
    !if ((nbisect<2).or.(ileftbisect < 0).or.(irightbisect>nchorizo)) stop ' bisectv6improve does not work! '
    !
    !ibisecttot=ibisecttot+1    
    !do i=0,nbisect
	   ! oldristra(i)=ristra(ileftbisect+i)
    !enddo
    !
    !newristra(0)=oldristra(0)
    !newristra(nbisect)=oldristra(nbisect)
    !
    !!do i=0,nbisect
    !!xoldtmp(i,:,:)=oldristra(i)%x
    !!enddo
    ! 
    !cwtbegintmp(:,:)=cwtl2r(:,:,ileftbisect)
    !cwtendtmp(:,:)=cwtr2l(:,:,ileftbisect+nbisect)  
    !
    !! tmp1 deal with oldlv
    !cwtl2rtmp1(:,:,0:nbisect)=cwtl2r(:,:,ileftbisect:ileftbisect+nbisect)
    !cwtr2ltmp1(:,:,0:nbisect)=cwtr2l(:,:,ileftbisect:ileftbisect+nbisect)
    !! tmp2 deal with newlv   
    !cwtl2rtmp2(:,:,0:nbisect)=cwtl2r(:,:,ileftbisect:ileftbisect+nbisect)
    !cwtr2ltmp2(:,:,0:nbisect)=cwtr2l(:,:,ileftbisect:ileftbisect+nbisect)
    !
    !oldlv=1.
    !newlv=1.   
    !
    !tau=dt*nbisect ! total time for the bisection slice.
    !do bilvl=1,mmax   ! level 1 to N. mmax = the total level N.
	   ! sigmamid=sqrt(hbar*tau*2.0_r8**(-bilvl))  
	   ! !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	   ! do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		  !  imid=2**(mmax-bilvl)+i*2**(mmax-bilvl+1)
		  !  il=imid-2**(mmax-bilvl)
		  !  ir=imid+2**(mmax-bilvl)
		  !  !write(6,*) 'bilvl il ir imid sigmamid=',bilvl, il,ir,imid,sigmamid
		  !  gauss(:,:)=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		  !  xmid(:,:)=(newristra(il)%x(:,:)+newristra(ir)%x(:,:))/2
		  !  newristra(imid)%x(:,:)=xmid(:,:)+gauss(:,:)
		  ! 
		  !  !xold=oldristra(imid)%x
		  !  !xnew=newristra(imid)%x
    !
    !!       if ( sum((newristra(imid)%x-oldristra(imid)%x)**2).eq.0) then
		  !  !  write(6,*) '-------------------------'
		  !  !  write(6,*) 'bilvl=',bilvl
		  !  !  write(6,*) 'imid,il,ir=',imid,il,ir
		  !  !  write(6,*) 'sigmamid=',sigmamid
		  !  !  write(6,*) 'gauss=',gauss
		  !  !     write(6,*) 'xnew=',newristra(imid)%x
		  !  !  write(6,*) 'xold=',oldristra(imid)%x
    !!           write(6,*) 'oldristra imid x=',oldristra(imid)%x
		  !  !  write(6,*) '-------------------------'
		  !  !
		  !  !  
    !!do j=0,nbisect
    !!write(6,*) 'j=',j  
    !!write(6,*) 'xold ristra diff=', oldristra(j)%x-xoldtmp(j,:,:)
    !!enddo			  
		  !  !  
		  !  !  stop
		  !  !  endif		   	   
		  !  !call hpsi(newristra(imid))
		  !  !call hpsi(oldristra(imid)) ! put here for save. But is not really needed if do not move beads before or after bisection
	   ! enddo 
    !enddo  
    !! up to here, newristra has been loaded.
    !
    !! judge if we accept the newristra or not.
    !! in this v6 case, v is not justpotential v, it should be the value of < ... >. lv.
		  !  !call hpsi0(newristra(imid))
		  !  !tmp=(0.5d0*(newristra(il)%v+newristra(ir)%v)-newristra(imid)%v &
    !    !          -0.5d0*(oldristra(il)%v+oldristra(ir)%v)+oldristra(imid)%v)*tau/(2.0_r8**bilvl) 	      
    !
    !all: do bilvl=1,mmax  ! besection.  level 1 to N. mmax = the total level N.	   
	   ! bisecttot(bilvl)=bisecttot(bilvl)+1.	   
	   ! lmax=2**bilvl-1 
	   ! jd=2**(mmax-bilvl) ! interval
	   ! dtexpt=mmax-bilvl-1 ! for the R L part, that's why the extra -1. dt=2**dtexp
	   ! do l=1,lmax
		  !  j=l*jd
		  !  jp=(l-1)*jd
		  !  xold=oldristra(j)%x
    !        xnew(:,:)=newristra(j)%x(:,:)
		  !  call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(:,:,jp),cwtl2rtmp1(:,:,j))    ! should be lr, not rl I believe.
		  !  call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(:,:,jp),cwtl2rtmp2(:,:,j))  	   
		  !  !if ( sum((xnew-xold)**2).eq.0) then
			 !   !write(6,*) 'bilvl=',bilvl	
			 !   !write(6,*) 'l=',l,j,jp
			 !   !write(6,*) 'sigmamid=',sigmamid
		  !  !   write(6,*) 'xnew-xold=',xnew-xold
		  !  !endif   
	   ! enddo
	   ! jmax=lmax*jd
	   ! newlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp2(:,:,jmax))*cwtendtmp(:,:)) )) 	   
	   ! oldlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp1(:,:,jmax))*cwtendtmp(:,:)) ))	   
	   ! tmp=log(newlv(bilvl))-log(oldlv(bilvl))+log(oldlv(bilvl-1))-log(newlv(bilvl-1))   
  	 !   !tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
	   ! rn=randn(1)	   
    !! can check if l2r and r2l give	the same lv.		    
	   !     !write(6,*) 'tmp ln rn=',tmp,log(rn(1))
	   ! if (log(rn(1)).lt.tmp) then
	   !     reject=.false. 
			 !
			 !   !if ( sum((xnew-xold)**2).eq.0) then
			 !   !
			 !   !write(6,*) 'bilvl=',bilvl	
			 !   !write(6,*) 'sigmamid=',sigmamid
		  !  !   write(6,*) 'xnew-xold=',xnew-xold
    !    !         write(6,*) 'tmp ln(rn)=',tmp,log(rn(1))		   
		  !  !   write (6,*) 'lvnew-lvold=',newlv(bilvl)-oldlv(bilvl)
			 !   !write (6,*) 'previous lvnew-lvold=',newlv(bilvl-1)-oldlv(bilvl-1)		  
				!    ! 
		  !  !
			 !   !endif  
				!
		  !  !if (mod(icstepstd,5*nstepdecor).eq.0) then   
		  !  ! !if (tmp.gt.300.) then
		  !  !  !write(6,*) 'bug bead appear! exit!'
		  !  !     write(6,*) 'ileft bilvl il ir imid sigmamid=',ileftbisect,bilvl,il,ir,imid,sigmamid	   
	   ! !        write(6,*) 'tmp ln(rn)=',tmp,log(rn(1))		   
		  !  !     write (6,*) 'lv=',oldlv(bilvl),newlv(bilvl),oldlv(bilvl-1),newlv(bilvl-1)
		  !  !  write(6,*) 'cwtbegintmp+*cwtendtmp=',sum(conjg(cwtbegintmp(:,:))*cwtendtmp(:,:))
		  !  !  !write(6,*) 'cwtl2rtmp1=',cwtl2rtmp1(:,:,jmax),jmax
		  !  !  !write(6,*) 'xold=',xold	
    !!           !write(6,*) 'xnew=',xnew	
		  !  !  !write(6,*) 'cwtl2rtmp2=',cwtl2rtmp2(:,:,jmax)	
		  !  !  !write(6,*) 'cwtendtmp=',cwtendtmp
		  !  !
		  !  !     !stop
		  !  !    endif
	   ! else
	   !     reject=.true.  
	   !     exit all 
	   ! endif 
		  !  !write(6,*) 'imid=',imid
    !
	   ! bisectcount(bilvl)=bisectcount(bilvl)+1.
    !enddo all
    !if (reject) then
	   ! call addval(2,0.0_r8,1.0_r8)
    !else
	   ! ibisect=ibisect+1 
	   ! call addval(2,1.0_r8,1.0_r8)  
    !
	   ! do i=0,nbisect
    !        ristra(ileftbisect+i)=newristra(i)
    !    enddo	 
    !!ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good.	 
    !
    !    xl=ristra(0)%x     
	   ! xr=ristra(nchorizo)%x     
    !! update cwtl2r, cwtr2l.	 
    !! update cwtl2r for all beads. this can be optimized in the future.
    ! 
	   ! cwtl2r(:,:,ileftbisect+1:ileftbisect+nbisect-1)=cwtl2rtmp2(:,:,1:jmax) ! here jmax should already be 2**mmax-1 
	   ! if ((ileftbisect+nbisect).le.nchorizo-1) then
    !    do k=ileftbisect+nbisect,nchorizo-1
	   ! x=ristra(k)%x
	   ! call v6proplr(x,x,-1,cwtl2r(:,:,k-1),cwtl2r(:,:,k))
    !    enddo
	   ! call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))
    !    else
    !    !ileftbisect+nbisect=nchorizo  
	   ! call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))   
    !    endif
    !
    !! update cwtr2l for all the beads
	   ! do k=ileftbisect+nbisect-1,1,-1
	   ! x=ristra(k)%x
	   ! call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
	   ! enddo
    !    call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))  	
    !    
    !    pathcoreval=newlv(mmax)
    !  
    !endif     
    !
    !rejectbisectv6improve=reject
    !return  
    !end subroutine bisectv6improve 
   
    !subroutine bisectv6improvehalf(ileftbisect) ! for v6 interaction, nbisect/2 beads
    !use wavefunction
    !use estimator
    !use random
    !use v6stepcalc
    !integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
				!	    ,irightbisect,k,bilvl &
	   !                 ,dtexpt
    !real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    !real(kind=r8) :: xmid(3,npart),xr(3,npart),xl(3,npart)
    !real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart)
    !real(kind=r8) :: tmp 
    !real(kind=r8) :: oldlv(0:mmaxhalf),newlv(0:mmaxhalf)
    !complex(kind=r8) :: cwtl2rtmp1(0:nspin-1,nisospin,0:nbisecthalf),cwtr2ltmp1(0:nspin-1,nisospin,0:nbisecthalf) &
	   !                 ,cwtl2rtmp2(0:nspin-1,nisospin,0:nbisecthalf),cwtr2ltmp2(0:nspin-1,nisospin,0:nbisecthalf) &
	   !                 ,cwtbegintmp(0:nspin-1,nisospin),cwtendtmp(0:nspin-1,nisospin)
    !logical :: reject	
    !
    !irightbisect=ileftbisect+nbisecthalf
    !if ((nbisecthalf<2).or.(ileftbisect < 0).or.(irightbisect>nchorizo))  stop ' bisectv6improve does not work! '   
    !
    !
    !ibisecthalftot=ibisecthalftot+1    
    !do i=0,nbisecthalf
	   ! oldristra(i)=ristra(ileftbisect+i)
    !enddo
    !
    !newristra(0)=oldristra(0)
    !newristra(nbisecthalf)=oldristra(nbisecthalf)
    !
    !!do i=0,nbisect
    !!xoldtmp(i,:,:)=oldristra(i)%x
    !!enddo
    ! 
    !cwtbegintmp(:,:)=cwtl2r(:,:,ileftbisect)
    !cwtendtmp(:,:)=cwtr2l(:,:,ileftbisect+nbisecthalf)  
    !
    !! tmp1 deal with oldlv
    !cwtl2rtmp1(:,:,0:nbisecthalf)=cwtl2r(:,:,ileftbisect:ileftbisect+nbisecthalf)
    !cwtr2ltmp1(:,:,0:nbisecthalf)=cwtr2l(:,:,ileftbisect:ileftbisect+nbisecthalf)
    !! tmp2 deal with newlv   
    !cwtl2rtmp2(:,:,0:nbisecthalf)=cwtl2r(:,:,ileftbisect:ileftbisect+nbisecthalf)
    !cwtr2ltmp2(:,:,0:nbisecthalf)=cwtr2l(:,:,ileftbisect:ileftbisect+nbisecthalf)
    !
    !oldlv=1.
    !newlv=1.   
    !
    !tau=dt*nbisecthalf ! total time for the bisection slice.
    !do bilvl=1,mmaxhalf   ! level 1 to N. mmax = the total level N.
	   ! sigmamid=sqrt(hbar*tau*2.0_r8**(-bilvl))  
	   ! !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	   ! do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		  !  imid=2**(mmaxhalf-bilvl)+i*2**(mmaxhalf-bilvl+1)
		  !  il=imid-2**(mmaxhalf-bilvl)
		  !  ir=imid+2**(mmaxhalf-bilvl)
		  !  !write(6,*) 'bilvl il ir imid sigmamid=',bilvl, il,ir,imid,sigmamid
		  !  gauss(:,:)=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		  !  xmid(:,:)=(newristra(il)%x(:,:)+newristra(ir)%x(:,:))/2
		  !  newristra(imid)%x(:,:)=xmid(:,:)+gauss(:,:)
	   ! enddo 
    !enddo  
    !! up to here, newristra has been loaded.
    !
    !! judge if we accept the newristra or not.
	   !   
    !all: do bilvl=1,mmaxhalf  ! besection.  level 1 to N. mmax = the total level N.	   
	   ! bisecthalftot(bilvl)=bisecthalftot(bilvl)+1.	   
	   ! lmax=2**bilvl-1 
	   ! jd=2**(mmaxhalf-bilvl) ! interval
	   ! dtexpt=mmaxhalf-bilvl-1 ! for the R L part, that's why the extra -1. dt=2**dtexp
	   ! do l=1,lmax
		  !  j=l*jd
		  !  jp=(l-1)*jd
		  !  xold=oldristra(j)%x
    !        xnew(:,:)=newristra(j)%x(:,:)
		  !  call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(:,:,jp),cwtl2rtmp1(:,:,j))    ! should be lr, not rl I believe.
		  !  call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(:,:,jp),cwtl2rtmp2(:,:,j))  	   
    !
	   ! enddo
	   ! jmax=lmax*jd
	   ! newlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp2(:,:,jmax))*cwtendtmp(:,:)) )) 	   
	   ! oldlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp1(:,:,jmax))*cwtendtmp(:,:)) ))	   
	   ! tmp=log(newlv(bilvl))-log(oldlv(bilvl))+log(oldlv(bilvl-1))-log(newlv(bilvl-1))   
  	 !   !tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
	   ! rn=randn(1)	   
    !! can check if l2r and r2l give	the same lv.		    
	   !     !write(6,*) 'tmp ln rn=',tmp,log(rn(1))
	   ! if (log(rn(1)).lt.tmp) then
	   !     reject=.false. 
    !        bisecthalfcount(bilvl)=bisecthalfcount(bilvl)+1.
	   ! else
	   !     reject=.true.  
	   !     exit all 
	   ! endif 
	   !
    !enddo all
    !if (reject) then
	   ! call addval(9,0.0_r8,1.0_r8)
    !else
	   ! ibisecthalf=ibisecthalf+1 
	   ! call addval(9,1.0_r8,1.0_r8)  
    !
	   ! do i=0,nbisecthalf
    !        ristra(ileftbisect+i)=newristra(i)
    !    enddo	 
    !!ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good.	 
    !
    !    xl=ristra(0)%x     
	   ! xr=ristra(nchorizo)%x     
    !
	   ! cwtl2r(:,:,ileftbisect+1:ileftbisect+nbisecthalf-1)=cwtl2rtmp2(:,:,1:jmax) ! here jmax should already be 2**mmax-1 
	   ! if ((ileftbisect+nbisecthalf).le.nchorizo-1) then
    !    do k=ileftbisect+nbisecthalf,nchorizo-1
	   ! x=ristra(k)%x
	   ! call v6proplr(x,x,-1,cwtl2r(:,:,k-1),cwtl2r(:,:,k))
    !    enddo
	   ! call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))
    !    else
    !    !ileftbisect+nbisect=nchorizo  
	   ! call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))   
    !    endif
    !
    !! update cwtr2l for all the beads
	   ! do k=ileftbisect+nbisecthalf-1,1,-1
	   ! x=ristra(k)%x
	   ! call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
	   ! enddo
    !    call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))  	
    !    
    !    pathcoreval=newlv(mmaxhalf)
    !endif     
    !
    !rejectbisectv6improvehalf=reject
    !return  
    !end subroutine bisectv6improvehalf 
    !
    
    !subroutine bisectv6rhalf(irightbisect) ! nbisecthalf, this is not efficient perhaps.
    !use wavefunction
    !use estimator
    !use random
    !use v6stepcalc
    !use brut
    !integer(kind=i4) :: ileftbisect,i &
				!	    ,irightbisect 
    !real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    !real(kind=r8) :: xnew(3,npart),xold(3,npart),xref(3,npart) &
    !                ,xd2old(3,npart),xd2new(3,npart)
    !real(kind=r8) :: tmp1,tmp2,logratio &
	   !             ,rt2,rt3
    !complex(kind=r8) :: cwtold(0:nspin-1,nisospin),cwtnew(0:nspin-1,nisospin)
    !complex(kind=r8) :: cwtl2rtmp(0:nspin-1,nisospin),cwtr2ltmp(0:nspin-1,nisospin)
	   !
    !ileftbisect=irightbisect-nbisecthalf
    !
    !if ((nbisecthalf<2).or.(ileftbisect < 0).or.(irightbisect>nchorizo)) stop ' bisectv6rhalf does not work! '
    !
    !ibisecttotr=ibisecttotr+1   
    !
    !ristraold(0:nchorizo)=ristra(0:nchorizo)
    !
    !i=irightbisect
    !
    !xold=ristraold(i)%x 
    !xref=ristraold(i-nbisecthalf)%x 
    !
    !    tau=dt*nbisect
    ! 
    !    sigmamid=sqrt(hbar*tau)  
    ! 
    !    gauss(:,:)=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
    ! 
    !    !gauss=mov1step*reshape((randn(3*npart)-0.5),(/3,npart/)) ! this is not gauss move. for vmc(dt=0) only
    ! 
	   ! xnew=xref+gauss
    !	 
	   ! if (i.eq.nchorizo) then   
	   ! !call psitcwt(xold,cwtold)
    !  
	   ! cwtold(:,:)=cwtend(:,:)
    !    cwtr2ltmp(:,:)=cwtr2l(:,:,i)
    !     
	   ! tmp1=log(abs(real(sum(conjg(cwtl2r(:,:,i))*cwtold(:,:)))))
	   ! !call psitcwt(xnew,cwtnew)
	   ! call corop(xnew,ipro,jpro,cwtnew)
	   ! call v6propl(xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmp(:,:))
	   ! tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtnew(:,:)))))
	   ! xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2
	   ! xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2 
	   ! rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	  
    !    !rt2=0. ! vmc(dt=0) only.
	   !
	   ! xd2new=(xref-xnew)**2
	   ! xd2old=(xref-xold)**2       
    !    rt3=1./(2.*tau*hbar)*(sum(xd2new)-sum(xd2old))
    !  
	   ! logratio=tmp2-tmp1+rt2+rt3
	   ! rn=randn(1)
	   !
	   ! if (log(rn(1)).lt.logratio) then
    !
	   ! ristra(i)%x(:,:)=xnew(:,:)
    !    cwtend(:,:)=cwtnew(:,:)
    !    call v6propr(xnew,-1,cwtend(:,:),cwtr2l(:,:,i))
    !    call bisectv6improvehalf(ileftbisect)
    !    if (rejectbisectv6improvehalf) then 
    !    ! set the values back. nothing changed.   
    !    cwtend=cwtold   
    !    ristra(i)%x=xold
    !    cwtr2l(:,:,i)=cwtr2ltmp(:,:) 
    !    else
    !        ibisectr=ibisectr+1   
    !    endif 
    !    endif 		     
    !
    !    else
	   !
    !    cwtr2ltmp(:,:)=cwtr2l(:,:,i)
    !    cwtold(:,:)=cwtl2r(:,:,i)
	   ! call v6proplr(xnew,xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmp(:,:))
	   !
	   ! tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtr2l(:,:,i+1)))))
	   ! tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i+1)))))
	   ! 
	   ! xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2+(ristra(i+1)%x(:,:)-xnew(:,:))**2
	   ! xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2+(ristra(i+1)%x(:,:)-xold(:,:))**2
	   !
	   ! rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	
    !    !rt2=0.  ! vmc(dt=0) only.
	   !    
    !    xd2new=(xref-xnew)**2
	   ! xd2old=(xref-xold)**2       
    !    rt3=1./(2.*tau*hbar)*(sum(xd2new)-sum(xd2old))
    !   
	   ! logratio=tmp2-tmp1+rt2+rt3	 
	   ! rn=randn(1)
    !
	   ! if (log(rn(1)).lt.logratio) then
	   ! ristra(i)%x(:,:)=xnew(:,:) 
    !    call v6proplr(xnew,xnew,-1,cwtr2l(:,:,i+1),cwtr2l(:,:,i))
    !    call bisectv6improvehalf(ileftbisect)
    !    if (rejectbisectv6improvehalf) then 
    !    ! set the values back. nothing changed.   
    !    cwtr2l(:,:,i)=cwtr2ltmp(:,:)   
    !    ristra(i)%x=xold
    !    
    !    pathcoreval=exp(tmp2)
    !    
    !    
    !    else
    !        ibisectr=ibisectr+1   
    !    endif 
    !    endif
    !     
	   ! endif 
    !
    !return  
    !end subroutine bisectv6rhalf    
    !

    !subroutine bisectv6lhalf(ileftbisect) ! move the left end, 0th bead,ileftbiect=0
    !use wavefunction
    !use estimator
    !use random
    !use v6stepcalc
    !use brut
    !integer(kind=i4) :: ileftbisect,i &
				!	    ,irightbisect 
    !real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    !real(kind=r8) :: xnew(3,npart),xold(3,npart),xref(3,npart) &
    !                ,xd2old(3,npart),xd2new(3,npart)
    !real(kind=r8) :: tmp1,tmp2,logratio &
	   !             ,rt2,rt3
    !complex(kind=r8) :: cwtold(0:nspin-1,nisospin),cwtnew(0:nspin-1,nisospin)
    !complex(kind=r8) :: cwtl2rtmp(0:nspin-1,nisospin),cwtr2ltmp(0:nspin-1,nisospin) 	              
    !
    !irightbisect=ileftbisect+nbisecthalf
    !if ((nbisecthalf<2).or.(ileftbisect < 0).or.(irightbisect>nchorizo)) stop ' bisectv6lhalf does not work! '
    !
    !ibisecttotl=ibisecttotl+1   
    !
    !ristraold(0:nchorizo)=ristra(0:nchorizo)
    !
    !i=ileftbisect
    !
    !xold=ristraold(i)%x 
    !xref=ristraold(i+nbisecthalf)%x 
    !
    !    tau=dt*nbisect
    ! 
    !    sigmamid=sqrt(hbar*tau)  
    ! 
    !    gauss(:,:)=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
    ! 
    !    !gauss=mov1step*reshape((randn(3*npart)-0.5),(/3,npart/)) ! this is not gauss move. for vmc(dt=0) only
    ! 
	   ! xnew=xref+gauss
    !	 
	   ! if (i.eq.0) then   
    !
	   ! cwtold(:,:)=cwtbegin(:,:)
    !    cwtl2rtmp(:,:)=cwtl2r(:,:,i)
    !    
	   ! tmp1=log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i)))))
	   ! !call psitcwt(xnew,cwtnew)  
	   ! call corop(xnew,iplo,jplo,cwtnew)
	   ! call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmp(:,:))
	   ! tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
	   ! xd2new(:,:)=(ristra(i+1)%x(:,:)-xnew(:,:))**2
	   ! xd2old(:,:)=(ristra(i+1)%x(:,:)-xold(:,:))**2
	   !
	   ! rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	  
    !    !rt2=0. ! vmc(dt=0) only.
	   !
	   ! xd2new=(xref-xnew)**2
	   ! xd2old=(xref-xold)**2       
    !    rt3=1./(2.*tau*hbar)*(sum(xd2new)-sum(xd2old))      
    !  
	   ! logratio=tmp2-tmp1+rt2+rt3
	   ! rn=randn(1)
    !
	   ! if (log(rn(1)).lt.logratio) then ! dt ne 0
	   ! ristra(i)%x(:,:)=xnew(:,:)
    !! update cwtbegin
    !    cwtbegin(:,:)=cwtnew(:,:)  
    !! update cwtl2r	   
	   ! call v6propr(xnew,-1,cwtnew(:,:),cwtl2r(:,:,i)) ! just update the i. 
    !
    !    call bisectv6improvehalf(ileftbisect)
    !    if (rejectbisectv6improvehalf) then 
    !    ! set the values back. nothing changed.   
    !    cwtbegin=cwtold   
    !    ristra(i)%x=xold
    !    cwtl2r(:,:,i)=cwtl2rtmp(:,:) 
    !    else
    !        ibisectl=ibisectl+1   
    !    endif   
	   ! endif       
    !  
    !    else
	   !     
    !    cwtold(:,:)=cwtl2r(:,:,i)
	   ! call v6proplr(xnew,xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmp(:,:))
	   !
	   ! tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtr2l(:,:,i+1)))))
	   ! tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i+1)))))
	   ! 
	   ! xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2+(ristra(i+1)%x(:,:)-xnew(:,:))**2
	   ! xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2+(ristra(i+1)%x(:,:)-xold(:,:))**2
	   !
	   ! rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	
    !    !rt2=0.  ! vmc(dt=0) only.
	   !    
    !    xd2new=(xref-xnew)**2
	   ! xd2old=(xref-xold)**2       
    !    rt3=1./(2.*tau*hbar)*(sum(xd2new)-sum(xd2old))
    !   
	   ! logratio=tmp2-tmp1+rt2+rt3	 
	   ! rn=randn(1)
    !
	   ! if (log(rn(1)).lt.logratio) then
	   ! ristra(i)%x(:,:)=xnew(:,:) 
    !    call v6proplr(xnew,xnew,-1,cwtl2r(:,:,i-1),cwtl2r(:,:,i))
    !    call bisectv6improvehalf(ileftbisect)
    !    if (rejectbisectv6improvehalf) then 
    !    ! set the values back. nothing changed.   
    !    cwtl2r(:,:,i)=cwtold
    !    ristra(i)%x=xold
    !    
    !    pathcoreval=exp(tmp2)
    !    else
    !        ibisectl=ibisectl+1   
    !    endif 
    !    endif
    !     
	   ! endif 
    !
    !return  
    !end subroutine bisectv6lhalf         
   
    !subroutine bisectv6la ! should be better than bisectv6l,rhalf     
    !! move the left end, 0th bead,ileftbiect=0; if ileftbisect /= 0 then the small modification on 0th level is needed.
    !use wavefunction
    !use estimator
    !use random
    !use v6stepcalc
    !use brut
    !integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
				!	    ,irightbisect,k,bilvl &
	   !                 ,dtexpt &
    !                    ,nbisectnow,mmaxnow
    !real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    !real(kind=r8) :: xmid(3,npart),xr(3,npart),xl(3,npart)
    !real(kind=r8) :: x(3,npart),xnew(3,npart),x0new(3,npart),xold(3,npart) 
    !real(kind=r8) :: tmp 
    !real(kind=r8) :: oldlv(0:mmaxhalf),newlv(0:mmaxhalf)
    !complex(kind=r8) :: cwtl2rtmp1(0:nspin-1,nisospin,0:nbisecthalf) &
	   !                 ,cwtl2rtmp2(0:nspin-1,nisospin,0:nbisecthalf) &
	   !                 ,cwtendtmp(0:nspin-1,nisospin) &
    !                    ,cwtbegintmp2(0:nspin-1,nisospin)        
    !logical :: reject	
    !
    !ibisecttotl=ibisecttotl+1  
    !
    !ileftbisect=0
    !irightbisect=ileftbisect+nbisecthalf
    !
    !!nbisectnow=irightbisect-ileftbisect ! has to be 2^n.
    !!mmaxnow=log(real(nbisectnow))/log(real(2)) ! mmaxnow = mmax-1
    !
    !nbisectnow=nbisecthalf
    !mmaxnow=mmaxhalf
    !
    !tau=dt*nbisectnow ! total time for the bisection slice. nisectnow is usually nbisecthalf
    !sigmamid=sqrt(2.*hbar*tau)  
    !
    !do i=0,nbisectnow
	   ! oldristra(i)=ristra(ileftbisect+i)
    !enddo   
    !
    !! deal with the two ends first, then do the classic bisection.   
    !
    !newristra(nbisectnow)=oldristra(nbisectnow)   
    !cwtendtmp(:,:)=cwtr2l(:,:,ileftbisect+nbisectnow)
    !
    !newristra(0)%x=oldristra(nbisectnow)%x+sigmamid*reshape(gaussian(3*npart),(/3,npart/))
    !
    !x0new(:,:)=newristra(0)%x(:,:)
    !
    !! i,jpro,i,jplo, the o means old are actually the current order.
    !
    !call corop(x0new,iplo,jplo,cwtbegintmp2)     
    !
    !oldlv(0)=1.
    !newlv(0)=1.
    !
    !cwtl2rtmp1(:,:,0)=cwtl2r(:,:,ileftbisect)
    !call v6propr(x0new,-1,cwtbegintmp2,cwtl2rtmp2(:,:,0)) 
    !
    !! prepare the further new bisection positions
    !do bilvl=1,mmaxnow   ! level 1 to N. mmax = the total level N.
	   ! sigmamid=sqrt(hbar*tau*2.0_r8**(-bilvl))  
	   ! !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	   ! do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		  !  imid=2**(mmaxnow-bilvl)+i*2**(mmaxnow-bilvl+1)
		  !  il=imid-2**(mmaxnow-bilvl)
		  !  ir=imid+2**(mmaxnow-bilvl)
		  !  !write(6,*) 'bilvl il ir imid sigmamid=',bilvl, il,ir,imid,sigmamid
		  !  gauss(:,:)=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		  !  xmid(:,:)=(newristra(il)%x(:,:)+newristra(ir)%x(:,:))/2
		  !  newristra(imid)%x(:,:)=xmid(:,:)+gauss(:,:)
	   ! enddo 
    !enddo      
    !! further bisection    
    !all: do bilvl=1,mmaxnow  ! besection.  level 1 to N. mmax = the total level N.	   
	   ! bisecttotl(bilvl)=bisecttotl(bilvl)+1.	   
	   ! lmax=2**bilvl-1 
	   ! jd=2**(mmaxnow-bilvl) ! interval
	   ! dtexpt=mmaxnow-bilvl-1 ! for the R L part, that's why the extra -1. dt=2**dtexp 
	   ! do l=1,lmax
		  !  j=l*jd
		  !  jp=(l-1)*jd
		  !  xold=oldristra(j)%x
    !        xnew(:,:)=newristra(j)%x(:,:)
		  !  call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(:,:,jp),cwtl2rtmp1(:,:,j))    ! should be lr, not rl I believe.
		  !  call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(:,:,jp),cwtl2rtmp2(:,:,j))  	   
		  !  !if ( sum((xnew-xold)**2).eq.0) then
			 !   !write(6,*) 'bilvl=',bilvl	
			 !   !write(6,*) 'l=',l,j,jp
			 !   !write(6,*) 'sigmamid=',sigmamid
		  !  !   write(6,*) 'xnew-xold=',xnew-xold
		  !  !endif   
	   ! enddo
	   ! jmax=lmax*jd  
	   ! oldlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp1(:,:,jmax))*cwtendtmp(:,:)) ))	
	   ! newlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp2(:,:,jmax))*cwtendtmp(:,:)) )) 
	   ! tmp=log(newlv(bilvl))-log(oldlv(bilvl))+log(oldlv(bilvl-1))-log(newlv(bilvl-1))   
  	 !   !tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
	   ! rn=randn(1)	   
    !! can check if l2r and r2l give	the same lv.		    
	   !     !write(6,*) 'tmp ln rn=',tmp,log(rn(1))
	   ! if (log(rn(1)).lt.tmp) then
	   !     reject=.false. 
    !        bisectcountl(bilvl)=bisectcountl(bilvl)+1.
	   ! else
	   !     reject=.true.  
	   !     exit all 
	   ! endif   
    !enddo all   
    !
    !if (.not.reject) then
	   ! ibisectl=ibisectl+1 
	   ! do i=0,nbisectnow
    !        ristra(ileftbisect+i)=newristra(i)
    !    enddo	 
    !!ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good.	 
    !    xl=ristra(0)%x     
	   ! xr=ristra(nchorizo)%x 
    !! update cwtl2r, cwtr2l.	 
    !    cwtbegin(:,:)=cwtbegintmp2(:,:)
	   ! cwtl2r(:,:,ileftbisect:ileftbisect+nbisectnow-1)=cwtl2rtmp2(:,:,0:jmax) ! here jmax should already be 2**mmax-1    
	   ! if ((ileftbisect+nbisectnow).le.nchorizo-1) then
    !    do k=ileftbisect+nbisectnow,nchorizo-1
	   ! x=ristra(k)%x
	   ! call v6proplr(x,x,-1,cwtl2r(:,:,k-1),cwtl2r(:,:,k))
    !    enddo
	   ! call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))
    !    else
    !    !ileftbisect+nbisect=nchorizo  
	   ! call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))   
    !    endif       
    !! update cwtr2l for all the beads
	   ! do k=ileftbisect+nbisectnow-1,1,-1
	   ! x=ristra(k)%x
	   ! call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
	   ! enddo
    !    call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))  
    !    pathcoreval=newlv(mmaxnow)
    !endif    
    !
    !return  
    !end subroutine bisectv6la  
    !
   
    !subroutine bisectv6ra ! move the right end. irightbisect=nchorizo
    !use wavefunction
    !use estimator
    !use random
    !use v6stepcalc
    !use brut
    !integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
				!	    ,irightbisect,k,bilvl &
	   !                 ,dtexpt &
    !                    ,nbisectnow,mmaxnow
    !real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    !real(kind=r8) :: xmid(3,npart),xl(3,npart)
    !real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart)
    !real(kind=r8) :: tmp 
    !real(kind=r8) :: oldlv(0:mmaxhalf),newlv(0:mmaxhalf)
    !complex(kind=r8) :: cwtl2rtmp1(0:nspin-1,nisospin,0:nbisecthalf),cwtl2rtmp2(0:nspin-1,nisospin,0:nbisecthalf) &
	   !                 ,cwtbegintmp(0:nspin-1,nisospin) &
    !                    ,cwtendtmp1now(0:nspin-1,nisospin),cwtendtmp2now(0:nspin-1,nisospin)   	              
    !                  
    !logical :: reject	
    !
    !ibisecttotr=ibisecttotr+1    
    !
    !ileftbisect=nchorizo-nbisecthalf
    !irightbisect=nchorizo
    !
    !!nbisectnow=irightbisect-ileftbisect
    !!mmaxnow=log(real(nbisectnow))/log(real(2))
    !
    !nbisectnow=nbisecthalf
    !mmaxnow=mmaxhalf   
    !
    !tau=dt*nbisectnow ! total time for the bisection slice.
    !sigmamid=sqrt(2.*hbar*tau)  
    !
    !do i=0,nbisectnow
	   ! oldristra(i)=ristra(ileftbisect+i)
    !enddo   
    !
    !! deal with the two ends first, then do the classic bisection.   
    !newristra(0)=oldristra(0) 
    !cwtbegintmp(:,:)=cwtl2r(:,:,ileftbisect)
    !
    !newristra(nbisectnow)%x=oldristra(0)%x+sigmamid*reshape(gaussian(3*npart),(/3,npart/))
    !
    !xold(:,:)=oldristra(nbisectnow)%x(:,:)
    !xnew(:,:)=newristra(nbisectnow)%x(:,:)    
    !
    !call corop(xnew,ipro,jpro,cwtendnew)     
    !
    !oldlv(0)=1.
    !newlv(0)=1.
    !
    !cwtl2rtmp1(:,:,0)=cwtl2r(:,:,ileftbisect)
    !cwtl2rtmp2(:,:,0)=cwtl2rtmp1(:,:,0)
    !call v6propr(oldristra(nbisectnow)%x,-1,cwtend,cwtendtmp1now)  
    !call v6propr(newristra(nbisectnow)%x,-1,cwtendnew,cwtendtmp2now)    
    ! 
    !do bilvl=1,mmaxnow   ! level 1 to N. mmax = the total level N.
	   ! sigmamid=sqrt(hbar*tau*2.0_r8**(-bilvl))  
	   ! !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	   ! do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		  !  imid=2**(mmaxnow-bilvl)+i*2**(mmaxnow-bilvl+1)
		  !  il=imid-2**(mmaxnow-bilvl)
		  !  ir=imid+2**(mmaxnow-bilvl)
		  !  !write(6,*) 'bilvl il ir imid sigmamid=',bilvl, il,ir,imid,sigmamid
		  !  gauss(:,:)=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		  !  xmid(:,:)=(newristra(il)%x(:,:)+newristra(ir)%x(:,:))/2
		  !  newristra(imid)%x(:,:)=xmid(:,:)+gauss(:,:)
	   ! enddo 
    !enddo     
    !! further bisection
    !all: do bilvl=1,mmaxnow  ! besection.  level 1 to N. mmax = the total level N.	   
	   ! bisecttotr(bilvl)=bisecttotr(bilvl)+1.	   
	   ! lmax=2**bilvl-1 
	   ! jd=2**(mmaxnow-bilvl) ! interval
	   ! dtexpt=mmaxnow-bilvl-1 ! for the R L part, that's why the extra -1. dt=2**dtexp         
	   ! do l=1,lmax
		  !  j=l*jd
		  !  jp=(l-1)*jd
		  !  xold=oldristra(j)%x
    !        xnew(:,:)=newristra(j)%x(:,:)
		  !  call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(:,:,jp),cwtl2rtmp1(:,:,j))    ! should be lr, not rl I believe.
		  !  call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(:,:,jp),cwtl2rtmp2(:,:,j))  	   
		  !  !if ( sum((xnew-xold)**2).eq.0) then
			 !   !write(6,*) 'bilvl=',bilvl	
			 !   !write(6,*) 'l=',l,j,jp
			 !   !write(6,*) 'sigmamid=',sigmamid
		  !  !   write(6,*) 'xnew-xold=',xnew-xold
		  !  !endif   
	   ! enddo
	   ! jmax=lmax*jd  
	   ! oldlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp1(:,:,jmax))*cwtendtmp1now(:,:)) ))	
	   ! newlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp2(:,:,jmax))*cwtendtmp2now(:,:)) )) 
	   ! tmp=log(newlv(bilvl))-log(oldlv(bilvl))+log(oldlv(bilvl-1))-log(newlv(bilvl-1))   
  	 !   !tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
	   ! rn=randn(1)	   
    !! can check if l2r and r2l give	the same lv.		    
	   !     !write(6,*) 'tmp ln rn=',tmp,log(rn(1))
	   ! if (log(rn(1)).lt.tmp) then
	   !     reject=.false. 
    !        bisectcountr(bilvl)=bisectcountr(bilvl)+1.
	   ! else
	   !     reject=.true.  
	   !     exit all 
	   ! endif 
    !enddo all   
    !
    !if (.not.reject) then
	   ! ibisectr=ibisectr+1 
	   ! do i=0,nbisectnow
    !        ristra(ileftbisect+i)=newristra(i)
    !    enddo	 
    !
    !!ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good.	 
	 	 !
    !! update cwtl2r, cwtr2l.	 
    !! update cwtl2r for all beads. this can be optimized in the future.
    !    cwtl2r(:,:,ileftbisect+1:ileftbisect+nbisectnow-1)=cwtl2rtmp2(:,:,1:jmax) ! here jmax should already be 2**mmax-1  
    !
    !    call v6propl(ristra(nchorizo)%x,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo)) ! ileftbisect+nbisectnow=nchorizo  
    ! 
    !    ! update cwtr2l for all the beads
    !    cwtend(:,:)=cwtendnew(:,:)          
    !    cwtr2l(:,:,nchorizo)=cwtendtmp2now
	   ! do k=nchorizo-1,1,-1
    !    x=ristra(k)%x
	   ! call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
	   ! enddo
    !    xl=ristra(0)%x
    !    call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))  	
    !    pathcoreval=newlv(mmaxnow)
    !endif    
    !                 
    !return  
    !end subroutine bisectv6ra   
    
    !   this is not the best bisectv6l and r, but is ok.
    !subroutine bisectv6r ! move the right end. irightbisect=nchorizo
    !ibisecttotr=ibisecttotr+1 
    !call bisectv6improve(nchorizo-nbisect)
    !if (.not.rejectbisectv6improve) ibisectr=ibisectr+1 
    !return  
    !end subroutine bisectv6r    
   
    !subroutine bisectv6l ! move the left end, 0th bead,ileftbiect=0
    !ibisecttotl=ibisecttotl+1  
    !call bisectv6improve(0)
    !if (.not.rejectbisectv6improve) ibisectl=ibisectl+1 
    !return  
    !end subroutine bisectv6l    
    !  
    
! ------------------------------------------------------------------------------    
    !subroutine writepath0(ipathtmpin,iploin,jploin,iproin,jproin,x0tot1din)
    !use mympi
    !integer(kind=i4) :: i,ipathtmpin,myunit1,myunit2
    !integer(kind=i4) :: iploin(npair),jploin(npair),iproin(npair),jproin(npair)	
    !integer(kind=i4) :: iplnow(npair),jplnow(npair),iprnow(npair),jprnow(npair)	 ! just for one time use.
    !real(kind=r8) :: x0tot1din(3*npart*(nchorizo+1)),x0tot1dnow(3*npart*(nchorizo+1))    
    !
    !if (myrank() /= 0) then
    !    call send(iploin,0,1)
    !    call send(jploin,0,2)
    !    call send(iproin,0,3)
    !    call send(jproin,0,4)
    !    call send(x0tot1din,0,5)        
    !else
    !    open(newunit=myunit1,form='unformatted',file='path0.unf',position='append')
    !        write(myunit1) ipathtmp  
    !        
    !        write (myunit1) i
    !        write (myunit1) iploin
    !        write (myunit1) jploin
    !        write (myunit1) iproin
    !        write (myunit1) jproin
    !        write (myunit1) x0tot1din
    !        
    !        do i=1,nproc()-1
    !            call recv(iplnow,i,1)
    !            call recv(jplnow,i,2)
    !            call recv(iprnow,i,3)
    !            call recv(jprnow,i,4)
    !            call recv(x0tot1dnow,i,5)   
    !            ! rank i write out. 
    !            write (myunit1) i
    !            write (myunit1) iplnow
    !            write (myunit1) jplnow
    !            write (myunit1) iprnow
    !            write (myunit1) jprnow
    !            write (myunit1) x0tot1dnow             
    !        enddo
    !    close (myunit1)   
    !        
    !    open(newunit=myunit2,form='formatted',file='path0numbercounts.txt',position='rewind')
    !        write (myunit2,'(3i10)') nproc(),ipathtmpin,size(x0tot1din)
    !    close(myunit2)            
    !        
    !        
    !endif
    !
    !return
    !end subroutine writepath0
    
    !subroutine stdmove1bisect1v6flex(ileft,iright,arrow) ! try to add importance sampling.
    !! move 1 by 1 with bisection 1 bead move.
    !use wavefunction
    !use estimator
    !use v6stepcalc
    !use random
    !use brut
    !use mympi
    !integer(kind=i4) :: i,k,ileft,iright,icnow
    !real(kind=r8) :: rn(1),gauss(3,npart),sigmamid,gc,sigmamid0,sigmamid1,mov1ratio,rn1(1)
    !real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart),xoldnxt(3,npart),xd2old(3,npart),xd2new(3,npart) &
	   !             ,xr(3,npart),xl(3,npart)
    !real(kind=r8) :: rt2,tmp1,tmp2,logratio,tmp1test1,tmp1test2
    !complex(kind=r8) :: cwtold(0:nspin-1,nisospin),cwtnew(0:nspin-1,nisospin)
    !complex(kind=r8) :: cwtl2rtmp(0:nspin-1,nisospin),cwtr2ltmp(0:nspin-1,nisospin)
    !integer(kind=i4) :: arrow
    !
    !if ((ileft.lt.0).or.(iright.gt.nchorizo).or.(iright.lt.0).or.(ileft.gt.nchorizo).or.(ileft.gt.iright)) then
	   ! write (6,*) 'stdmove1fastv6flex range error', ileft,iright
	   ! stop
    !endif
    !if (ileft.eq.iright) then
    !    call stdmove1fastv6(ileft,iright) 
    !    return
    !endif     
    !gc=-0.5_r8/(dt*(2*hbar))
    !sigmamid0=sqrt(2*hbar*dt)  
    !sigmamid1=sqrt(hbar*dt)
    !
    !mov1ratio=0.5 ! the ratio of bisect 1 move.
    !select case (arrow)   
    !case (-1) ! inverse order from right to left  
    !    icnow=0
    !    icstdmove1bisect1rvstot=icstdmove1bisect1rvstot+1
    !    !icmov1bisect1rvstot=abs(iright-ileft+1)*icstdmove1bisect1rvstot 
    !    do i=iright,ileft,-1  
    !        rn1=randn(1)
    !        if (i.eq.0) then       
    !            if (rn1(1)<mov1ratio) then 
    !                icmov1bisect1rvstot=icmov1bisect1rvstot+1
	   !             xold=ristra(i+1)%x 
    !                gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                sigmamid=sigmamid0  
	   !             xnew=xold+sigmamid*gauss
	   !             cwtold(:,:)=cwtbegin(:,:)
	   !             tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i)))))      
	   !             call corop(xnew,iplo,jplo,cwtnew)
	   !             call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmp(:,:))
	   !             tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
	   !             logratio=tmp2-tmp1             
    !            else
    !                icmov1rvstot=icmov1rvstot+1
    !                xold=ristra(i)%x 
    !                gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                xnew(:,:)=xold(:,:)+gauss(:,:)  
    !                cwtold(:,:)=cwtbegin(:,:)
    !                tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i)))))       
    !                call corop(xnew,iplo,jplo,cwtnew)
    !                call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmp(:,:))
    !                tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
    !                xd2new(:,:)=(ristra(i+1)%x(:,:)-xnew(:,:))**2
    !                xd2old(:,:)=(ristra(i+1)%x(:,:)-xold(:,:))**2
    !                rt2=gc*(sum(xd2new)-sum(xd2old))	  
    !                logratio=tmp2-tmp1+rt2         
    !            endif            
	   !         rn=randn(1)
	   !         if (log(rn(1)).lt.logratio) then ! dt ne 0
    !                icnow=icnow+1 
    !                if (rn1(1)<mov1ratio) then 
	   !                 icmov1bisect1rvs=icmov1bisect1rvs+1	
    !                else
    !                    icmov1rvs=icmov1rvs+1   
    !                endif                    
	   !             ristra(i)%x(:,:)=xnew(:,:)
	   !             cwtr2l(:,:,i)=cwtr2ltmp(:,:)   
    !            ! update cwtbegin
    !                cwtbegin(:,:)=cwtnew(:,:)  
    !                pathcoreval=exp(tmp2)
    !            endif 
	   !     else
	   !         if (i.eq.nchorizo) then                  
    !                if (rn1(1)<mov1ratio) then
    !                    icmov1bisect1rvstot=icmov1bisect1rvstot+1
	   !                 xold=ristra(i-1)%x 
    !                    gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                    sigmamid=sigmamid0  
	   !                 xnew=xold+sigmamid*gauss 
	   !                 cwtold(:,:)=cwtend(:,:)
	   !                 tmp1=log(pathcoreval)  !log(abs(real(sum(conjg(cwtl2r(:,:,i))*cwtold(:,:))))) ! can put those tmp1 thing in the beginning, no need to repreatedly calculate                
	   !                 call corop(xnew,ipro,jpro,cwtnew)
    !                    call v6propr(xnew,-1,cwtnew(:,:),cwtr2ltmp(:,:))          
    !                    tmp2=log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmp(:,:)))))     
	   !                 logratio=tmp2-tmp1
    !                else
    !                    icmov1rvstot=icmov1rvstot+1
    !                    xold=ristra(i)%x 
    !                    gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                    xnew(:,:)=xold(:,:)+gauss(:,:)                          
	   !                 cwtold(:,:)=cwtend(:,:)
	   !                 tmp1=log(pathcoreval)  !log(abs(real(sum(conjg(cwtl2r(:,:,i))*cwtold(:,:))))) ! can put those tmp1 thing in the beginning, no need to repreatedly calculate                
	   !                 call corop(xnew,ipro,jpro,cwtnew)
    !                    call v6propr(xnew,-1,cwtnew(:,:),cwtr2ltmp(:,:))          
    !                    tmp2=log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmp(:,:)))))
	   !                 xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2
	   !                 xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2
	   !                 rt2=gc*(sum(xd2new)-sum(xd2old))	  
	   !                 logratio=tmp2-tmp1+rt2                    
    !                endif                
	   !             rn=randn(1)
	   !             if (log(rn(1)).lt.logratio) then
    !                    icnow=icnow+1
    !                    if (rn1(1)<mov1ratio) then 
	   !                     icmov1bisect1rvs=icmov1bisect1rvs+1	
    !                    else
    !                        icmov1rvs=icmov1rvs+1   
    !                    endif
	   !                 ristra(i)%x(:,:)=xnew(:,:)
	   !                 cwtr2l(:,:,i)=cwtr2ltmp(:,:)
    !                ! update cwtend and prepare for the left bead.
    !                    cwtend(:,:)=cwtnew(:,:)
    !                    x=ristra(i-1)%x
		  !              call v6proplr(x,x,-1,cwtr2l(:,:,i),cwtr2l(:,:,i-1))
    !                    pathcoreval=exp(tmp2)
	   !             endif 		     
    !            else 
    !                if (rn1(1)<mov1ratio) then   
    !                    icmov1bisect1rvstot=icmov1bisect1rvstot+1
	   !                 xold=(ristra(i-1)%x+ristra(i+1)%x)/2 
    !                    gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                    sigmamid=sigmamid1  
	   !                 xnew=xold+sigmamid*gauss                     
	   !                 call v6proplr(xnew,xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmp(:,:))
	   !                 tmp2=log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmp(:,:)))))
    !                    tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2l(:,:,i)))))
	   !                 logratio=tmp2-tmp1   
    !                else
    !                    icmov1rvstot=icmov1rvstot+1
    !                    xold=ristra(i)%x 
    !                    gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                    xnew(:,:)=xold(:,:)+gauss(:,:)                           
	   !                 call v6proplr(xnew,xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmp(:,:))
	   !                 tmp2=log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmp(:,:)))))
    !                    tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2l(:,:,i)))))  
	   !                 xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2+(ristra(i+1)%x(:,:)-xnew(:,:))**2
	   !                 xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2+(ristra(i+1)%x(:,:)-xold(:,:))**2   
	   !                 rt2=gc*(sum(xd2new)-sum(xd2old))		   
	   !                 logratio=tmp2-tmp1+rt2                                            
    !                endif           
	   !             rn=randn(1)
	   !             if (log(rn(1)).lt.logratio) then
    !                    icnow=icnow+1           
    !                    if (rn1(1)<mov1ratio) then 
	   !                     icmov1bisect1rvs=icmov1bisect1rvs+1	
    !                    else
    !                        icmov1rvs=icmov1rvs+1   
    !                    endif	
	   !                 ristra(i)%x(:,:)=xnew(:,:)
	   !                 cwtr2l(:,:,i)=cwtr2ltmp(:,:)	
    !                    pathcoreval=exp(tmp2)
    !                endif
    !                if ( icnow /= 0 ) then
    !            ! update r2l for the next r2l bead, the left bead
    !                    if (i >= 2) then
		  !                  x=ristra(i-1)%x
		  !                  call v6proplr(x,x,-1,cwtr2l(:,:,i),cwtr2l(:,:,i-1))
    !                    else ! i=1
    !                        xl=ristra(i-1)%x
    !                        call v6propl(xl,-1,cwtr2l(:,:,i),cwtr2l(:,:,i-1))	    
    !                    endif               
    !                endif 
	   !         endif 
	   !     endif                   
    !    enddo
    !    if (icnow /= 0) then 
    !        xl=ristra(0)%x
    !        xr=ristra(nchorizo)%x
    !        ! update r2l
    !        if ( ileft >= 3 ) then
    !            do k=ileft,3,-1 
    !               x=ristra(k-2)%x
	   !            call v6proplr(x,x,-1,cwtr2l(:,:,k-1),cwtr2l(:,:,k-2))  
    !            enddo
    !            call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	 
    !        else
    !            if (ileft == 2) call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	 
    !        endif
    !        ! update l2r
    !        if (ileft >= 1) then
    !            if (ileft <= nchorizo-1 ) then    
    !                do k=ileft,nchorizo-1
    !                   x=ristra(k)%x
	   !                call v6proplr(x,x,-1,cwtl2r(:,:,k-1),cwtl2r(:,:,k))
    !                enddo
    !                call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))	   
    !            else ! ileft = nchorizo
    !                call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))	   
    !            endif 
    !        else ! ileft == 0
    !            call v6propr(xl,-1,cwtbegin(:,:),cwtl2r(:,:,0))	   
    !                do k=1,nchorizo-1
    !                   x=ristra(k)%x
	   !                call v6proplr(x,x,-1,cwtl2r(:,:,k-1),cwtl2r(:,:,k))
    !                enddo
    !            call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))	     
    !        endif
    !    endif
    !case (1) ! normal, from left to right
    !    icnow=0
    !    icstdmove1bisect1tot=icstdmove1bisect1tot+1
    !    !icmov1bisect1tot=abs(iright-ileft+1)*icstdmove1bisect1tot       
    !    do i=ileft,iright
    !        rn1=randn(1)
    !        if (i.eq.0) then  
    !            if (rn1(1)<mov1ratio) then 
    !                icmov1bisect1tot=icmov1bisect1tot+1
	   !             xold=ristra(i+1)%x 
    !                gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                sigmamid=sigmamid0  
	   !             xnew=xold+sigmamid*gauss               
	   !             cwtold(:,:)=cwtbegin(:,:)
	   !             tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i))))) 
	   !             call corop(xnew,iplo,jplo,cwtnew)
	   !             call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmp(:,:))
	   !             tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
	   !             logratio=tmp2-tmp1  
    !            else   
    !                icmov1tot=icmov1tot+1
    !                xold=ristra(i)%x 
    !                gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                xnew(:,:)=xold(:,:)+gauss(:,:)   
	   !             cwtold(:,:)=cwtbegin(:,:)
	   !             tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i)))))
	   !             call corop(xnew,iplo,jplo,cwtnew)
	   !             call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmp(:,:))
	   !             tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
	   !             xd2new(:,:)=(ristra(i+1)%x(:,:)-xnew(:,:))**2
	   !             xd2old(:,:)=(ristra(i+1)%x(:,:)-xold(:,:))**2  
	   !             rt2=gc*(sum(xd2new)-sum(xd2old))	  
	   !             logratio=tmp2-tmp1+rt2                             
    !            endif         
	   !         rn=randn(1)
	   !         if (log(rn(1)).lt.logratio) then ! dt ne 0
    !                icnow=icnow+1 
    !                if (rn1(1)<mov1ratio) then  
	   !                 icmov1bisect1=icmov1bisect1+1	
    !                else
    !                    icmov1=icmov1+1	  
    !                endif         
	   !             ristra(i)%x(:,:)=xnew(:,:)
	   !             cwtr2l(:,:,i)=cwtr2ltmp(:,:)   
    !            ! update cwtbegin
    !                cwtbegin(:,:)=cwtnew(:,:)  
    !            ! update cwtl2r	   
	   !             call v6propr(xnew,-1,cwtnew(:,:),cwtl2r(:,:,i)) ! just update the i.
    !            ! prepare for the next cwtl2r
    !                xoldnxt=ristra(i+1)%x
    !                call v6proplr(xoldnxt,xoldnxt,-1,cwtl2r(:,:,i),cwtl2r(:,:,i+1))   
    !                pathcoreval=exp(tmp2)
	   !         endif 
	   !     else
	   !         if (i.eq.nchorizo) then              
    !                if (rn1(1)<mov1ratio) then   
    !                    icmov1bisect1tot=icmov1bisect1tot+1
	   !                 xold=ristra(i-1)%x 
    !                    gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                    sigmamid=sigmamid0  
	   !                 xnew=xold+sigmamid*gauss 
	   !                 cwtold(:,:)=cwtend(:,:)
	   !                 tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtl2r(:,:,i))*cwtold(:,:)))))
	   !                 call corop(xnew,ipro,jpro,cwtnew)
	   !                 call v6propl(xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmp(:,:))
	   !                 tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtnew(:,:)))))	  
	   !                 logratio=tmp2-tmp1            
    !                else
    !                    icmov1tot=icmov1tot+1
    !                    xold=ristra(i)%x 
    !                    gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                    xnew(:,:)=xold(:,:)+gauss(:,:)                       
	   !                 cwtold(:,:)=cwtend(:,:)
	   !                 tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtl2r(:,:,i))*cwtold(:,:)))))
	   !                 call corop(xnew,ipro,jpro,cwtnew)
	   !                 call v6propl(xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmp(:,:))
	   !                 tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtnew(:,:)))))
	   !                 xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2
	   !                 xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2
	   !                 rt2=gc*(sum(xd2new)-sum(xd2old))	  
	   !                 logratio=tmp2-tmp1+rt2                    
    !                endif 
	   !             rn=randn(1)
	   !             if (log(rn(1)).lt.logratio) then
    !                    icnow=icnow+1
    !                    if (rn1(1)<mov1ratio) then  
	   !                     icmov1bisect1=icmov1bisect1+1	
    !                    else
    !                        icmov1=icmov1+1	  
    !                    endif
	   !                 ristra(i)%x(:,:)=xnew(:,:)
	   !                 cwtl2r(:,:,i)=cwtl2rtmp(:,:)
    !                ! update cwtend
    !                    cwtend(:,:)=cwtnew(:,:)                  
    !                    pathcoreval=exp(tmp2)
	   !             endif 		      
    !            else
    !                if (rn1(1)<mov1ratio) then 
    !                    icmov1bisect1tot=icmov1bisect1tot+1
	   !                 xold=(ristra(i-1)%x+ristra(i+1)%x)/2 
    !                    gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                    sigmamid=sigmamid1  
	   !                 xnew=xold+sigmamid*gauss 
    !                    cwtold(:,:)=cwtl2r(:,:,i)
	   !                 call v6proplr(xnew,xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmp(:,:))
	   !                 tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtr2l(:,:,i+1)))))
	   !                 tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i+1)))))
	   !                 logratio=tmp2-tmp1                  
    !                else
    !                    icmov1tot=icmov1tot+1
    !                    xold=ristra(i)%x 
    !                    gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                    xnew(:,:)=xold(:,:)+gauss(:,:)   
    !                    cwtold(:,:)=cwtl2r(:,:,i)
	   !                 call v6proplr(xnew,xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmp(:,:)) 
	   !                 tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtr2l(:,:,i+1)))))
	   !                 tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i+1)))))
	   !                 xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2+(ristra(i+1)%x(:,:)-xnew(:,:))**2
	   !                 xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2+(ristra(i+1)%x(:,:)-xold(:,:))**2	   
	   !                 rt2=gc*(sum(xd2new)-sum(xd2old))	
	   !                 logratio=tmp2-tmp1+rt2	                              
    !                endif      
	   !             rn=randn(1)
	   !             if (log(rn(1)).lt.logratio) then
    !                        !write(6,*) 'chorizo i updated!'
    !                    icnow=icnow+1           
    !                    if (rn1(1)<mov1ratio) then  
	   !                     icmov1bisect1=icmov1bisect1+1	
    !                    else
    !                        icmov1=icmov1+1	  
    !                    endif	
	   !                 ristra(i)%x(:,:)=xnew(:,:)
	   !                 cwtl2r(:,:,i)=cwtl2rtmp(:,:)	    
    !                    pathcoreval=exp(tmp2)
    !                endif
    !                if ( icnow /= 0 ) then
    !            ! update l2r
    !                    if (i.le.(nchorizo-2)) then
		  !              x=ristra(i+1)%x
		  !              call v6proplr(x,x,-1,cwtl2r(:,:,i),cwtl2r(:,:,i+1))
    !                    else ! i=nchorizo-1
    !                    xr=ristra(i+1)%x
    !                    call v6propl(xr,-1,cwtl2r(:,:,i),cwtl2r(:,:,i+1))	    
    !                    endif               
    !                endif 
	   !         endif 
	   !     endif 	 
    !    enddo 
    !    if (icnow /= 0) then 
    !        xl=ristra(0)%x
    !        xr=ristra(nchorizo)%x
    !        ! update cwtl2r    
    !        if (iright <= nchorizo-3) then 
    !            do k=iright,nchorizo-3
    !                x=ristra(k+2)%x
    !                call v6proplr(x,x,-1,cwtl2r(:,:,k+1),cwtl2r(:,:,k+2))   
    !            enddo
    !            call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))	      
    !        else    
    !            if (iright == nchorizo-2) then   
    !                call v6propl(xr,-1,cwtl2r(:,:,iright+1),cwtl2r(:,:,nchorizo))	          
    !            endif  
    !        endif    
    !        ! update cwtr2l from iright th bead until the 0th.
    !        if (iright == nchorizo) then    
    !            !x=ristra(iright)%x   
	   !         call v6propr(xr,-1,cwtend(:,:),cwtr2l(:,:,iright))
	   !         do k=iright-1,1,-1
	   !         x=ristra(k)%x
	   !         call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
	   !         enddo
	   !         !xl=ristra(0)%x
	   !         call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	    
    !        else      
    !            if (iright >= 1) then  
    !                do k=iright,1,-1
	   !                x=ristra(k)%x
	   !                call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
	   !             enddo
	   !             !xl=ristra(0)%x
	   !             call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	   
    !            else ! iright = 0
    !                call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	   
    !            endif  
    !        endif
    !    endif   
    !case default 
    !    write (6,*) 'arrow for stdmove1bisect1flex has not been set, abort!', arrow
	   ! call abort
    !end select
    !return
    !end subroutine stdmove1bisect1v6flex              
    !  
    
    !subroutine stdmove1bisect1v6flexa(ileft,iright,arrow,imp) ! add importance sampling. perhaps there is no need to add imp for bisect 1.
    !! move 1 by 1 with bisection 1 bead move.
    !use wavefunction
    !use estimator
    !use v6stepcalc
    !use random
    !use brut
    !use mympi
    !integer(kind=i4) :: i,k,ileft,iright,icnow
    !real(kind=r8) :: rn(1),gauss(3,npart),sigmamid,gc,sigmamid0,sigmamid1,mov1ratio,rn1(1)
    !real(kind=r8) :: x(3,npart),xnew(3,npart),xnew1(3,npart),xold(3,npart),xold1(3,npart) &
    !                ,xoldnxt(3,npart),xd2old(3,npart),xd2new(3,npart),xd2old1(3,npart),xd2new1(3,npart) &
	   !             ,xr(3,npart),xl(3,npart)
    !real(kind=r8) :: rt2,tmp1,tmp2,tmp3,tmp4,q1,q2,gnew,gold,gnew1,gold1,logratio,tmp1test1,tmp1test2
    !complex(kind=r8) :: cwtold(0:nspin-1,nisospin),cwtnew(0:nspin-1,nisospin) &
    !                   ,cwtold1(0:nspin-1,nisospin),cwtnew1(0:nspin-1,nisospin)
    !complex(kind=r8) :: cwtl2rtmp(0:nspin-1,nisospin),cwtr2ltmp(0:nspin-1,nisospin) &
    !                   ,cwtl2rtmpnew(0:nspin-1,nisospin),cwtl2rtmpnew1(0:nspin-1,nisospin),cwtl2rtmpold1(0:nspin-1,nisospin) &
    !                   ,cwtr2ltmpnew(0:nspin-1,nisospin),cwtr2ltmpnew1(0:nspin-1,nisospin),cwtr2ltmpold1(0:nspin-1,nisospin)
    !integer(kind=i4) :: arrow
    !logical :: imp
    !
    !if ((ileft.lt.0).or.(iright.gt.nchorizo).or.(iright.lt.0).or.(ileft.gt.nchorizo).or.(ileft.gt.iright)) then
	   ! write (6,*) 'stdmove1fastv6flex range error', ileft,iright
	   ! stop
    !endif
    !if (ileft.eq.iright) then
    !    call stdmove1fastv6(ileft,iright) 
    !    return
    !endif       
    !gc=-0.5_r8/(dt*(2*hbar))
    !sigmamid0=sqrt(2*hbar*dt)  
    !sigmamid1=sqrt(hbar*dt)    
    !
    !mov1ratio=0.1 ! the ratio of bisect 1 move.
    !select case (arrow)    
    !case (-1) ! inverse order from right to left  
    !    icnow=0
    !    icstdmove1bisect1rvstot=icstdmove1bisect1rvstot+1
    !    !icmov1bisect1rvstot=abs(iright-ileft+1)*icstdmove1bisect1rvstot 
    !    do i=iright,ileft,-1  
    !        rn1=randn(1)
    !        if (i.eq.0) then      
    !            if (.not.imp) then
    !                if (rn1(1)<mov1ratio) then 
    !                    icmov1bisect1rvstot=icmov1bisect1rvstot+1
	   !                 xold=ristra(i+1)%x 
    !                    gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                    sigmamid=sigmamid0  
	   !                 xnew=xold+sigmamid*gauss
	   !                 tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i)))))      
	   !                 call corop(xnew,iplo,jplo,cwtnew)
	   !                 call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmp(:,:))
	   !                 tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
	   !                 logratio=tmp2-tmp1             
    !                else
    !                    icmov1rvstot=icmov1rvstot+1
    !                    xold=ristra(i)%x 
    !                    gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                    xnew(:,:)=xold(:,:)+gauss(:,:)  
    !                    tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i)))))       
    !                    call corop(xnew,iplo,jplo,cwtnew)
    !                    call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmp(:,:))
    !                    tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
    !                    xd2new(:,:)=(ristra(i+1)%x(:,:)-xnew(:,:))**2
    !                    xd2old(:,:)=(ristra(i+1)%x(:,:)-xold(:,:))**2
    !                    rt2=gc*(sum(xd2new)-sum(xd2old))	  
    !                    logratio=tmp2-tmp1+rt2         
    !                endif            
	   !             rn=randn(1)
	   !             if (log(rn(1)).lt.logratio) then ! dt ne 0
    !                    icnow=icnow+1 
    !                    if (rn1(1)<mov1ratio) then 
	   !                     icmov1bisect1rvs=icmov1bisect1rvs+1	
    !                    else
    !                        icmov1rvs=icmov1rvs+1   
    !                    endif                    
	   !                 ristra(i)%x(:,:)=xnew(:,:)
	   !                 cwtr2l(:,:,i)=cwtr2ltmp(:,:)   
    !                ! update cwtbegin
    !                    cwtbegin(:,:)=cwtnew(:,:)  
    !                    pathcoreval=exp(tmp2)
    !                endif              
    !            else ! imp   
    !                if (rn1(1)<mov1ratio) then        ! i=0
    !                    icmov1bisect1rvstot=icmov1bisect1rvstot+1
    !                    gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                    sigmamid=sigmamid0 
    !                    xold=ristra(i+1)%x  ! ref point.
    !                    xnew=xold+sigmamid*gauss     
    !                    xold1=2*xold-ristra(i)%x 
    !                    xnew1=xold-sigmamid*gauss
    !                    call corop(xnew,iplo,jplo,cwtnew)
    !                    call corop(xold1,iplo,jplo,cwtold1)
    !                    call corop(xnew1,iplo,jplo,cwtnew1)
	   !                 call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew(:,:))        
	   !                 call v6propl(xold1,-1,cwtr2l(:,:,i+1),cwtr2ltmpold1(:,:))  
	   !                 call v6propl(xnew1,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew1(:,:))  
    !                    tmp1=pathcoreval 
	   !                 tmp2=abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmpnew(:,:))))                            
    !                    tmp3=abs(real(sum(conjg(cwtold1(:,:))*cwtr2ltmpold1(:,:))))                              
	   !                 tmp4=abs(real(sum(conjg(cwtnew1(:,:))*cwtr2ltmpnew1(:,:))))                                       
    !                    q1=(tmp2+tmp4)/(tmp1+tmp3)
    !                    q2=tmp2/(tmp2+tmp4)  
    !                else     
    !                    icmov1rvstot=icmov1rvstot+1
    !                    gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                    xold=ristra(i)%x 
    !                    xnew(:,:)=xold(:,:)+gauss(:,:)  
    !                    xold1=2*xnew-xold
    !                    xnew1=xold-gauss
    !                    xd2new(:,:)=(ristra(i+1)%x(:,:)-xnew(:,:))**2
    !                    xd2old(:,:)=(ristra(i+1)%x(:,:)-xold(:,:))**2
    !                    xd2new1=(ristra(i+1)%x-xnew1)**2
    !                    xd2old1=(ristra(i+1)%x-xold1)**2    
    !                    gnew=exp(gc*sum(xd2new))
    !                    gold=exp(gc*sum(xd2old))
    !                    gnew1=exp(gc*sum(xd2new1))
    !                    gold1=exp(gc*sum(xd2old1))
    !                    call corop(xnew,iplo,jplo,cwtnew)
    !                    call corop(xold1,iplo,jplo,cwtold1)
    !                    call corop(xnew1,iplo,jplo,cwtnew1)
    !                    call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew(:,:)) 
    !                    call v6propl(xold1,-1,cwtr2l(:,:,i+1),cwtr2ltmpold1(:,:))
	   !                 call v6propl(xnew1,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew1(:,:)) 
    !                    tmp1=pathcoreval
    !                    tmp2=abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmpnew(:,:))))       
    !                    tmp3=abs(real(sum(conjg(cwtold1(:,:))*cwtr2ltmpold1(:,:))))                         
	   !                 tmp4=abs(real(sum(conjg(cwtnew1(:,:))*cwtr2ltmpnew1(:,:))))                         
    !                    q1=(gnew*tmp2+gnew1*tmp4)/(gold*tmp1+gold1*tmp3)
    !                    q2=gnew*tmp2/(gnew*tmp2+gnew1*tmp4)                            
    !                endif
	   !             rn=randn(1)
	   !             if (rn(1) < q1) then ! dt ne 0
    !                    icnow=icnow+1 
    !                    if (rn1(1)<mov1ratio) then 
	   !                     icmov1bisect1rvs=icmov1bisect1rvs+1	
    !                    else
    !                        icmov1rvs=icmov1rvs+1   
    !                    endif                     
    !                    rn=randn(1)
    !                    if (rn(1) < q2) then
	   !                     ristra(i)%x(:,:)=xnew(:,:)
	   !                     cwtr2l(:,:,i)=cwtr2ltmpnew(:,:)   
    !                        cwtbegin(:,:)=cwtnew(:,:)  
    !                        pathcoreval=tmp2                                    
    !                    else
	   !                     ristra(i)%x=xnew1
	   !                     cwtr2l(:,:,i)=cwtr2ltmpnew1(:,:)   
    !                        cwtbegin=cwtnew1	   
    !                        pathcoreval=tmp4                                   
    !                    endif
    !                endif                                         
    !            endif
    !        else  
	   !         if (i.eq.nchorizo) then      
    !                if (.not.imp) then           
    !                    if (rn1(1)<mov1ratio) then
    !                        icmov1bisect1rvstot=icmov1bisect1rvstot+1
	   !                     xold=ristra(i-1)%x 
    !                        gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                        sigmamid=sigmamid0  
	   !                     xnew=xold+sigmamid*gauss 
	   !                     tmp1=log(pathcoreval)  !log(abs(real(sum(conjg(cwtl2r(:,:,i))*cwtold(:,:))))) ! can put those tmp1 thing in the beginning, no need to repreatedly calculate                
	   !                     call corop(xnew,ipro,jpro,cwtnew)
    !                        call v6propr(xnew,-1,cwtnew(:,:),cwtr2ltmp(:,:))          
    !                        tmp2=log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmp(:,:)))))     
	   !                     logratio=tmp2-tmp1
    !                    else
    !                        icmov1rvstot=icmov1rvstot+1
    !                        xold=ristra(i)%x 
    !                        gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                        xnew(:,:)=xold(:,:)+gauss(:,:)                          
	   !                     tmp1=log(pathcoreval)  !log(abs(real(sum(conjg(cwtl2r(:,:,i))*cwtold(:,:))))) ! can put those tmp1 thing in the beginning, no need to repreatedly calculate                
	   !                     call corop(xnew,ipro,jpro,cwtnew)
    !                        call v6propr(xnew,-1,cwtnew(:,:),cwtr2ltmp(:,:))          
    !                        tmp2=log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmp(:,:)))))
	   !                     xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2
	   !                     xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2
	   !                     rt2=gc*(sum(xd2new)-sum(xd2old))	  
	   !                     logratio=tmp2-tmp1+rt2                    
    !                    endif                
	   !                 rn=randn(1)
	   !                 if (log(rn(1)).lt.logratio) then
    !                        icnow=icnow+1
    !                        if (rn1(1)<mov1ratio) then 
	   !                         icmov1bisect1rvs=icmov1bisect1rvs+1	
    !                        else
    !                            icmov1rvs=icmov1rvs+1   
    !                        endif
	   !                     ristra(i)%x(:,:)=xnew(:,:)
	   !                     cwtr2l(:,:,i)=cwtr2ltmp(:,:)
    !                    ! update cwtend and prepare for the left bead.
    !                        cwtend(:,:)=cwtnew(:,:)
    !                        x=ristra(i-1)%x
		  !                  call v6proplr(x,x,-1,cwtr2l(:,:,i),cwtr2l(:,:,i-1))
    !                        pathcoreval=exp(tmp2)
    !                    endif 	               
    !                else ! imp  i=nchorizo                
    !                    if (rn1(1)<mov1ratio) then     
    !                        icmov1bisect1rvstot=icmov1bisect1rvstot+1
    !                        gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                        sigmamid=sigmamid0
    !                        xold=ristra(i-1)%x  ! ref point
    !                        xnew=xold+sigmamid*gauss     
    !                        xold1=2*xold-ristra(i)%x 
    !                        xnew1=xold-sigmamid*gauss                             
    !                        call corop(xnew,ipro,jpro,cwtnew)
    !                        call corop(xold1,ipro,jpro,cwtold1)
    !                        call corop(xnew1,ipro,jpro,cwtnew1)
	   !                     call v6propr(xnew,-1,cwtnew(:,:),cwtr2ltmpnew(:,:))                            
	   !                     call v6propr(xold1,-1,cwtold1(:,:),cwtr2ltmpold1(:,:))                
	   !                     call v6propr(xnew1,-1,cwtnew1(:,:),cwtr2ltmpnew1(:,:))
	   !                     tmp1=pathcoreval
    !                        tmp2=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpnew(:,:))))
    !                        tmp3=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpold1(:,:)))) 
    !                        tmp4=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpnew1(:,:))))                                        
    !                        q1=(tmp2+tmp4)/(tmp1+tmp3)
    !                        q2=tmp2/(tmp2+tmp4)  
    !                    else  
    !                        icmov1rvstot=icmov1rvstot+1
    !                        gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                        xold=ristra(i)%x 
    !                        xnew(:,:)=xold(:,:)+gauss(:,:)  
    !                        xold1=2*xnew-xold
    !                        xnew1=xold-gauss
    !                        xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2
    !                        xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2
    !                        xd2new1=(ristra(i-1)%x-xnew1)**2
    !                        xd2old1=(ristra(i-1)%x-xold1)**2    
    !                        gnew=exp(gc*sum(xd2new))
    !                        gold=exp(gc*sum(xd2old))
    !                        gnew1=exp(gc*sum(xd2new1))
    !                        gold1=exp(gc*sum(xd2old1)) 
    !                        call corop(xnew,ipro,jpro,cwtnew)
    !                        call corop(xold1,ipro,jpro,cwtold1)
    !                        call corop(xnew1,ipro,jpro,cwtnew1)                           
    !                        call v6propr(xnew,-1,cwtnew(:,:),cwtr2ltmpnew(:,:))                         
    !                        call v6propr(xold1,-1,cwtold1(:,:),cwtr2ltmpold1(:,:))          
	   !                     call v6propr(xnew1,-1,cwtnew1(:,:),cwtr2ltmpnew1(:,:))
    !                        tmp1=pathcoreval
    !                        tmp2=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpnew(:,:))))              
    !                        tmp3=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpold1(:,:))))                                 
	   !                     tmp4=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpnew1(:,:))))                         
    !                        q1=(gnew*tmp2+gnew1*tmp4)/(gold*tmp1+gold1*tmp3)
    !                        q2=gnew*tmp2/(gnew*tmp2+gnew1*tmp4)                            
    !                    endif     
	   !                 rn=randn(1)
	   !                 if (rn(1) < q1) then ! dt ne 0
    !                        icnow=icnow+1 
    !                        if (rn1(1)<mov1ratio) then 
	   !                         icmov1bisect1rvs=icmov1bisect1rvs+1	
    !                        else
    !                            icmov1rvs=icmov1rvs+1   
    !                        endif                     
    !                        rn=randn(1)
    !                        if (rn(1) < q2) then
	   !                         ristra(i)%x(:,:)=xnew(:,:)
	   !                         cwtr2l(:,:,i)=cwtr2ltmpnew(:,:)   
    !                            cwtend(:,:)=cwtnew(:,:)	   
    !                            pathcoreval=tmp2                             
    !                        else
	   !                         ristra(i)%x=xnew1
	   !                         cwtr2l(:,:,i)=cwtr2ltmpnew1(:,:)   
    !                            cwtend=cwtnew1	   
    !                            pathcoreval=tmp4                                   
    !                        endif
    !                        x=ristra(i-1)%x
    !                        call v6proplr(x,x,-1,cwtr2l(:,:,i),cwtr2l(:,:,i-1)) 
    !                    endif                             
    !                endif       
    !            else                 
    !                if (.not.imp) then 
    !                    if (rn1(1)<mov1ratio) then   
    !                        icmov1bisect1rvstot=icmov1bisect1rvstot+1
	   !                     xold=(ristra(i-1)%x+ristra(i+1)%x)/2 
    !                        gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                        sigmamid=sigmamid1  
	   !                     xnew=xold+sigmamid*gauss                     
	   !                     call v6proplr(xnew,xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmp(:,:))
	   !                     tmp2=log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmp(:,:)))))
    !                        tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2l(:,:,i)))))
	   !                     logratio=tmp2-tmp1   
    !                    else
    !                        icmov1rvstot=icmov1rvstot+1
    !                        xold=ristra(i)%x 
    !                        gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                        xnew(:,:)=xold(:,:)+gauss(:,:)                           
	   !                     call v6proplr(xnew,xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmp(:,:))
	   !                     tmp2=log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmp(:,:)))))
    !                        tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2l(:,:,i)))))  
	   !                     xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2+(ristra(i+1)%x(:,:)-xnew(:,:))**2
	   !                     xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2+(ristra(i+1)%x(:,:)-xold(:,:))**2   
	   !                     rt2=gc*(sum(xd2new)-sum(xd2old))		   
	   !                     logratio=tmp2-tmp1+rt2                                            
    !                    endif           
	   !                 rn=randn(1)
	   !                 if (log(rn(1)).lt.logratio) then
    !                        icnow=icnow+1           
    !                        if (rn1(1)<mov1ratio) then 
	   !                         icmov1bisect1rvs=icmov1bisect1rvs+1	
    !                        else
    !                            icmov1rvs=icmov1rvs+1   
    !                        endif	
	   !                     ristra(i)%x(:,:)=xnew(:,:)
	   !                     cwtr2l(:,:,i)=cwtr2ltmp(:,:)	
    !                        pathcoreval=exp(tmp2)
    !                    endif                 
    !                else  ! imp        
    !                    if (rn1(1)<mov1ratio) then 
    !                        icmov1bisect1rvstot=icmov1bisect1rvstot+1
    !                        gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                        sigmamid=sigmamid1
    !                        xold=(ristra(i-1)%x+ristra(i+1)%x)/2  ! ref point
    !                        xnew=xold+sigmamid*gauss 
    !                        xnew1=xold-sigmamid*gauss
    !                        xold1=2*xold-ristra(i)%x 
    !                        call v6proplr(xnew,xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew(:,:))
    !                        call v6proplr(xnew1,xnew1,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew1(:,:))
    !                        call v6proplr(xold1,xold1,-1,cwtr2l(:,:,i+1),cwtr2ltmpold1(:,:))    
    !                        tmp1=pathcoreval              
	   !                     tmp2=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpnew(:,:))))                         
    !                        tmp3=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpold1(:,:))))                              
	   !                     tmp4=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpnew1(:,:))))                                       
    !                        q1=(tmp2+tmp4)/(tmp1+tmp3)
    !                        q2=tmp2/(tmp2+tmp4)                                
    !                    else            
    !                        icmov1rvstot=icmov1rvstot+1
    !                        gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                        xold=ristra(i)%x
    !                        xnew(:,:)=xold(:,:)+gauss(:,:)  
    !                        xold1=2*xnew-xold
    !                        xnew1=xold-gauss
    !                        xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2+(ristra(i+1)%x(:,:)-xnew(:,:))**2
    !                        xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2+(ristra(i+1)%x(:,:)-xold(:,:))**2
    !                        xd2new1=(ristra(i-1)%x-xnew1)**2+(ristra(i+1)%x-xnew1)**2
    !                        xd2old1=(ristra(i-1)%x-xold1)**2+(ristra(i+1)%x-xold1)**2    
    !                        gnew=exp(gc*sum(xd2new))
    !                        gold=exp(gc*sum(xd2old))
    !                        gnew1=exp(gc*sum(xd2new1))
    !                        gold1=exp(gc*sum(xd2old1)) 
    !                        call v6proplr(xnew,xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew(:,:))
    !                        call v6proplr(xold1,xold1,-1,cwtr2l(:,:,i+1),cwtr2ltmpold1(:,:))
    !                        call v6proplr(xnew1,xnew1,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew1(:,:))         
    !                        tmp1=pathcoreval                   
    !                        tmp2=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpnew(:,:))))                          
    !                        tmp3=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpold1(:,:))))                                             
	   !                     tmp4=abs(real(sum(conjg(cwtl2r(:,:,i-1))*cwtr2ltmpnew1(:,:))))                                     
    !                        q1=(gnew*tmp2+gnew1*tmp4)/(gold*tmp1+gold1*tmp3)
    !                        q2=gnew*tmp2/(gnew*tmp2+gnew1*tmp4)                                  
    !                    endif                         
	   !                 rn=randn(1)
	   !                 if (rn(1) < q1) then ! dt ne 0
    !                        icnow=icnow+1 
    !                        if (rn1(1)<mov1ratio) then 
	   !                         icmov1bisect1rvs=icmov1bisect1rvs+1	
    !                        else
    !                            icmov1rvs=icmov1rvs+1   
    !                        endif                     
    !                        rn=randn(1)
    !                        if (rn(1) < q2) then
	   !                         ristra(i)%x(:,:)=xnew(:,:)
	   !                         cwtr2l(:,:,i)=cwtr2ltmpnew(:,:)   	   
    !                            pathcoreval=tmp2                             
    !                        else
	   !                         ristra(i)%x=xnew1
	   !                         cwtr2l(:,:,i)=cwtr2ltmpnew1(:,:)   	   
    !                            pathcoreval=tmp4                                   
    !                        endif
    !                    endif                                
    !                endif
    !                if ( icnow /= 0 ) then
    !            ! update r2l for the next r2l bead, the left bead
    !                    if (i >= 2) then
		  !                  x=ristra(i-1)%x
		  !                  call v6proplr(x,x,-1,cwtr2l(:,:,i),cwtr2l(:,:,i-1))
    !                    else ! i=1
    !                        xl=ristra(i-1)%x
    !                        call v6propl(xl,-1,cwtr2l(:,:,i),cwtr2l(:,:,i-1))	    
    !                    endif               
    !                endif 
	   !         endif 
	   !     endif                   
    !    enddo
    !    if (icnow /= 0) then 
    !        xl=ristra(0)%x
    !        xr=ristra(nchorizo)%x
    !        ! update r2l
    !        if ( ileft >= 3 ) then
    !            do k=ileft,3,-1 
    !               x=ristra(k-2)%x
	   !            call v6proplr(x,x,-1,cwtr2l(:,:,k-1),cwtr2l(:,:,k-2))  
    !            enddo
    !            call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	 
    !        else
    !            if (ileft == 2) call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	 
    !        endif
    !        ! update l2r
    !        if (ileft >= 1) then
    !            if (ileft <= nchorizo-1 ) then    
    !                do k=ileft,nchorizo-1
    !                   x=ristra(k)%x
	   !                call v6proplr(x,x,-1,cwtl2r(:,:,k-1),cwtl2r(:,:,k))
    !                enddo
    !                call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))	   
    !            else ! ileft = nchorizo
    !                call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))	   
    !            endif 
    !        else ! ileft == 0
    !            call v6propr(xl,-1,cwtbegin(:,:),cwtl2r(:,:,0))	   
    !                do k=1,nchorizo-1
    !                   x=ristra(k)%x
	   !                call v6proplr(x,x,-1,cwtl2r(:,:,k-1),cwtl2r(:,:,k))
    !                enddo
    !            call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))	     
    !        endif
    !    endif
    !case (1) ! normal, from left to right
    !    icnow=0
    !    icstdmove1bisect1tot=icstdmove1bisect1tot+1
    !    !icmov1bisect1tot=abs(iright-ileft+1)*icstdmove1bisect1tot       
    !    do i=ileft,iright
    !        rn1=randn(1)
    !        if (i.eq.0) then    
    !            if (.not.imp) then
    !                if (rn1(1)<mov1ratio) then 
    !                    icmov1bisect1tot=icmov1bisect1tot+1
	   !                 xold=ristra(i+1)%x 
    !                    gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                    sigmamid=sigmamid0  
	   !                 xnew=xold+sigmamid*gauss               
	   !                 cwtold(:,:)=cwtbegin(:,:)
	   !                 tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i))))) 
	   !                 call corop(xnew,iplo,jplo,cwtnew)
	   !                 call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmp(:,:))
	   !                 tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
	   !                 logratio=tmp2-tmp1  
    !                else   
    !                    icmov1tot=icmov1tot+1
    !                    xold=ristra(i)%x 
    !                    gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                    xnew(:,:)=xold(:,:)+gauss(:,:)   
	   !                 cwtold(:,:)=cwtbegin(:,:)
	   !                 tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i)))))
	   !                 call corop(xnew,iplo,jplo,cwtnew)
	   !                 call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmp(:,:))
	   !                 tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
	   !                 xd2new(:,:)=(ristra(i+1)%x(:,:)-xnew(:,:))**2
	   !                 xd2old(:,:)=(ristra(i+1)%x(:,:)-xold(:,:))**2  
	   !                 rt2=gc*(sum(xd2new)-sum(xd2old))	  
	   !                 logratio=tmp2-tmp1+rt2                             
    !                endif         
	   !             rn=randn(1)
	   !             if (log(rn(1)).lt.logratio) then ! dt ne 0
    !                    icnow=icnow+1 
    !                    if (rn1(1)<mov1ratio) then  
	   !                     icmov1bisect1=icmov1bisect1+1	
    !                    else
    !                        icmov1=icmov1+1	  
    !                    endif         
	   !                 ristra(i)%x(:,:)=xnew(:,:)
	   !                 cwtr2l(:,:,i)=cwtr2ltmp(:,:)   
    !                ! update cwtbegin
    !                    cwtbegin(:,:)=cwtnew(:,:)  
    !                ! update cwtl2r	   
	   !                 call v6propr(xnew,-1,cwtnew(:,:),cwtl2r(:,:,i)) ! just update the i.
    !                ! prepare for the next cwtl2r
    !                    xoldnxt=ristra(i+1)%x
    !                    call v6proplr(xoldnxt,xoldnxt,-1,cwtl2r(:,:,i),cwtl2r(:,:,i+1))   
    !                    pathcoreval=exp(tmp2)
    !                endif          
    !            else ! imp i=0
    !                if (rn1(1)<mov1ratio) then        
    !                    icmov1bisect1tot=icmov1bisect1tot+1
    !                    gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                    sigmamid=sigmamid0 
    !                    xold=ristra(i+1)%x  ! ref point.
    !                    xnew=xold+sigmamid*gauss     
    !                    xold1=2*xold-ristra(i)%x 
    !                    xnew1=xold-sigmamid*gauss
    !                    call corop(xnew,iplo,jplo,cwtnew)
    !                    call corop(xold1,iplo,jplo,cwtold1)
    !                    call corop(xnew1,iplo,jplo,cwtnew1)
    !                    call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew(:,:))
    !                    call v6propl(xold1,-1,cwtr2l(:,:,i+1),cwtr2ltmpold1(:,:))
    !                    call v6propl(xnew1,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew1(:,:))      
    !                    tmp1=pathcoreval      
	   !                 tmp2=abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmpnew(:,:))))                            
    !                    tmp3=abs(real(sum(conjg(cwtold1(:,:))*cwtr2ltmpold1(:,:))))                               
	   !                 tmp4=abs(real(sum(conjg(cwtnew1(:,:))*cwtr2ltmpnew1(:,:))))                                       
    !                    q1=(tmp2+tmp4)/(tmp1+tmp3)
    !                    q2=tmp2/(tmp2+tmp4)  
    !                else     
    !                    icmov1tot=icmov1tot+1
    !                    gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                    xold=ristra(i)%x
    !                    xnew(:,:)=xold(:,:)+gauss(:,:)  
    !                    xold1=2*xnew-xold
    !                    xnew1=xold-gauss
    !                    xd2new(:,:)=(ristra(i+1)%x(:,:)-xnew(:,:))**2
    !                    xd2old(:,:)=(ristra(i+1)%x(:,:)-xold(:,:))**2
    !                    xd2new1=(ristra(i+1)%x-xnew1)**2
    !                    xd2old1=(ristra(i+1)%x-xold1)**2    
    !                    gnew=exp(gc*sum(xd2new))
    !                    gold=exp(gc*sum(xd2old))
    !                    gnew1=exp(gc*sum(xd2new1))
    !                    gold1=exp(gc*sum(xd2old1)) 
    !                    call corop(xnew,iplo,jplo,cwtnew)
    !                    call corop(xold1,iplo,jplo,cwtold1)
    !                    call corop(xnew1,iplo,jplo,cwtnew1)
    !                    call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew(:,:))
    !                    call v6propl(xold1,-1,cwtr2l(:,:,i+1),cwtr2ltmpold1(:,:))
    !                    call v6propl(xnew1,-1,cwtr2l(:,:,i+1),cwtr2ltmpnew1(:,:))               
    !                    tmp1=pathcoreval         
    !                    tmp2=abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmpnew(:,:))))           
    !                    tmp3=abs(real(sum(conjg(cwtold1(:,:))*cwtr2ltmpold1(:,:))))                                         
	   !                 tmp4=abs(real(sum(conjg(cwtnew1(:,:))*cwtr2ltmpnew1(:,:))))                         
    !                    q1=(gnew*tmp2+gnew1*tmp4)/(gold*tmp1+gold1*tmp3)
    !                    q2=gnew*tmp2/(gnew*tmp2+gnew1*tmp4)                            
    !                endif                    
	   !             rn=randn(1)
	   !             if (rn(1) < q1) then ! dt ne 0
    !                    icnow=icnow+1 
    !                    if (rn1(1)<mov1ratio) then 
	   !                     icmov1bisect1=icmov1bisect1+1	
    !                    else
    !                        icmov1=icmov1+1   
    !                    endif                     
    !                    rn=randn(1)
    !                    if (rn(1) < q2) then
	   !                     ristra(i)%x(:,:)=xnew(:,:)
	   !                     cwtr2l(:,:,i)=cwtr2ltmpnew(:,:)   
    !                        cwtbegin(:,:)=cwtnew(:,:)  
    !                        pathcoreval=tmp2                        
    !                    ! update cwtl2r	   
	   !                     call v6propr(xnew,-1,cwtnew(:,:),cwtl2r(:,:,i)) ! just update the i.
    !                    ! prepare for the next cwtl2r
    !                        xoldnxt=ristra(i+1)%x
    !                        call v6proplr(xoldnxt,xoldnxt,-1,cwtl2r(:,:,i),cwtl2r(:,:,i+1))                               
    !                    else
	   !                     ristra(i)%x=xnew1
	   !                     cwtr2l(:,:,i)=cwtr2ltmpnew1(:,:)   
    !                        cwtbegin=cwtnew1	   
    !                        pathcoreval=tmp4   
    !                    ! update cwtl2r	   
	   !                     call v6propr(xnew1,-1,cwtnew1(:,:),cwtl2r(:,:,i)) ! just update the i.
    !                    ! prepare for the next cwtl2r
    !                        xoldnxt=ristra(i+1)%x
    !                        call v6proplr(xoldnxt,xoldnxt,-1,cwtl2r(:,:,i),cwtl2r(:,:,i+1))  
    !                    endif
    !                endif                        
    !            endif  
	   !     else
	   !         if (i.eq.nchorizo) then   
    !                if (.not.imp) then
    !                    if (rn1(1)<mov1ratio) then   
    !                        icmov1bisect1tot=icmov1bisect1tot+1
	   !                     xold=ristra(i-1)%x 
    !                        gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                        sigmamid=sigmamid0  
	   !                     xnew=xold+sigmamid*gauss 
	   !                     cwtold(:,:)=cwtend(:,:)
	   !                     tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtl2r(:,:,i))*cwtold(:,:)))))
	   !                     call corop(xnew,ipro,jpro,cwtnew)
	   !                     call v6propl(xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmp(:,:))
	   !                     tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtnew(:,:)))))	  
	   !                     logratio=tmp2-tmp1            
    !                    else
    !                        icmov1tot=icmov1tot+1
    !                        xold=ristra(i)%x 
    !                        gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                        xnew(:,:)=xold(:,:)+gauss(:,:)                       
	   !                     cwtold(:,:)=cwtend(:,:)
	   !                     tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtl2r(:,:,i))*cwtold(:,:)))))
	   !                     call corop(xnew,ipro,jpro,cwtnew)
	   !                     call v6propl(xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmp(:,:))
	   !                     tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtnew(:,:)))))
	   !                     xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2
	   !                     xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2
	   !                     rt2=gc*(sum(xd2new)-sum(xd2old))	  
	   !                     logratio=tmp2-tmp1+rt2                    
    !                    endif 
	   !                 rn=randn(1)
	   !                 if (log(rn(1)).lt.logratio) then
    !                        icnow=icnow+1
    !                        if (rn1(1)<mov1ratio) then  
	   !                         icmov1bisect1=icmov1bisect1+1	
    !                        else
    !                            icmov1=icmov1+1	  
    !                        endif
	   !                     ristra(i)%x(:,:)=xnew(:,:)
	   !                     cwtl2r(:,:,i)=cwtl2rtmp(:,:)
    !                    ! update cwtend
    !                        cwtend(:,:)=cwtnew(:,:)                  
    !                        pathcoreval=exp(tmp2)
    !                    endif                        
    !                else ! imp               
    !                    if (rn1(1)<mov1ratio) then     
    !                        icmov1bisect1tot=icmov1bisect1tot+1
    !                        gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                        sigmamid=sigmamid0 
    !                        xold=ristra(i-1)%x  ! ref point
    !                        xnew=xold+sigmamid*gauss     
    !                        xold1=2*xold-ristra(i)%x 
    !                        xnew1=xold-sigmamid*gauss
    !                        call corop(xnew,ipro,jpro,cwtnew)
    !                        call corop(xold1,ipro,jpro,cwtold1)
    !                        call corop(xnew1,ipro,jpro,cwtnew1)
    !                        call v6propl(xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmpnew(:,:))
    !                        call v6propl(xold1,-1,cwtl2r(:,:,i-1),cwtl2rtmpold1(:,:))
    !                        call v6propl(xnew1,-1,cwtl2r(:,:,i-1),cwtl2rtmpnew1(:,:))                          
    !                        tmp1=pathcoreval                          
	   !                     tmp2=abs(real(sum(conjg(cwtl2rtmpnew(:,:))*cwtnew(:,:))))	                                                
	   !                     tmp3=abs(real(sum(conjg(cwtl2rtmpold1(:,:))*cwtold1(:,:))))	                                                  
	   !                     tmp4=abs(real(sum(conjg(cwtl2rtmpnew1(:,:))*cwtnew1(:,:))))	                                     
    !                        q1=(tmp2+tmp4)/(tmp1+tmp3)
    !                        q2=tmp2/(tmp2+tmp4)  
    !                    else  
    !                        icmov1tot=icmov1tot+1
    !                        gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                        xold=ristra(i)%x 
    !                        xnew(:,:)=xold(:,:)+gauss(:,:)  
    !                        xold1=2*xnew-xold
    !                        xnew1=xold-gauss
    !                        xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2
    !                        xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2
    !                        xd2new1=(ristra(i-1)%x-xnew1)**2
    !                        xd2old1=(ristra(i-1)%x-xold1)**2    
    !                        gnew=exp(gc*sum(xd2new))
    !                        gold=exp(gc*sum(xd2old))
    !                        gnew1=exp(gc*sum(xd2new1))
    !                        gold1=exp(gc*sum(xd2old1)) 
    !                        call corop(xnew,ipro,jpro,cwtnew)
    !                        call corop(xold1,ipro,jpro,cwtold1)
    !                        call corop(xnew1,ipro,jpro,cwtnew1)
    !                        call v6propl(xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmpnew(:,:))
    !                        call v6propl(xold1,-1,cwtl2r(:,:,i-1),cwtl2rtmpold1(:,:))
    !                        call v6propl(xnew1,-1,cwtl2r(:,:,i-1),cwtl2rtmpnew1(:,:))                         
    !                        tmp1=pathcoreval               
	   !                     tmp2=abs(real(sum(conjg(cwtl2rtmpnew(:,:))*cwtnew(:,:))))                        
	   !                     tmp3=abs(real(sum(conjg(cwtl2rtmpold1(:,:))*cwtold1(:,:))))                                              
	   !                     tmp4=abs(real(sum(conjg(cwtl2rtmpnew1(:,:))*cwtnew1(:,:))))                                             
    !                        q1=(gnew*tmp2+gnew1*tmp4)/(gold*tmp1+gold1*tmp3)
    !                        q2=gnew*tmp2/(gnew*tmp2+gnew1*tmp4)                            
    !                    endif                             
	   !                 rn=randn(1)
	   !                 if (rn(1) < q1) then ! dt ne 0
    !                        icnow=icnow+1 
    !                        if (rn1(1)<mov1ratio) then 
	   !                         icmov1bisect1=icmov1bisect1+1	
    !                        else
    !                            icmov1=icmov1+1   
    !                        endif                     
    !                        rn=randn(1)
    !                        if (rn(1) < q2) then
	   !                         ristra(i)%x(:,:)=xnew(:,:)
	   !                         cwtl2r(:,:,i)=cwtl2rtmpnew(:,:)
    !                            cwtend(:,:)=cwtnew(:,:)	   
    !                            pathcoreval=tmp2                             
    !                        else
	   !                         ristra(i)%x=xnew1
	   !                         cwtl2r(:,:,i)=cwtl2rtmpnew1(:,:)
    !                            cwtend=cwtnew1	   
    !                            pathcoreval=tmp4                                   
    !                        endif 
    !                    endif                                
    !                endif
    !            else            
    !                if (.not.imp) then
    !                    if (rn1(1)<mov1ratio) then 
    !                        icmov1bisect1tot=icmov1bisect1tot+1
	   !                     xold=(ristra(i-1)%x+ristra(i+1)%x)/2 
    !                        gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                        sigmamid=sigmamid1  
	   !                     xnew=xold+sigmamid*gauss 
    !                        cwtold(:,:)=cwtl2r(:,:,i)
	   !                     call v6proplr(xnew,xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmp(:,:))
	   !                     tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtr2l(:,:,i+1)))))
	   !                     tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i+1)))))
	   !                     logratio=tmp2-tmp1                  
    !                    else
    !                        icmov1tot=icmov1tot+1
    !                        xold=ristra(i)%x 
    !                        gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                        xnew(:,:)=xold(:,:)+gauss(:,:)   
    !                        cwtold(:,:)=cwtl2r(:,:,i)
	   !                     call v6proplr(xnew,xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmp(:,:)) 
	   !                     tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtr2l(:,:,i+1)))))
	   !                     tmp1=log(pathcoreval) !log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i+1)))))
	   !                     xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2+(ristra(i+1)%x(:,:)-xnew(:,:))**2
	   !                     xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2+(ristra(i+1)%x(:,:)-xold(:,:))**2	   
	   !                     rt2=gc*(sum(xd2new)-sum(xd2old))	
	   !                     logratio=tmp2-tmp1+rt2	                              
    !                    endif      
	   !                 rn=randn(1)
	   !                 if (log(rn(1)).lt.logratio) then
    !                            !write(6,*) 'chorizo i updated!'
    !                        icnow=icnow+1           
    !                        if (rn1(1)<mov1ratio) then  
	   !                         icmov1bisect1=icmov1bisect1+1	
    !                        else
    !                            icmov1=icmov1+1	  
    !                        endif	
	   !                     ristra(i)%x(:,:)=xnew(:,:)
	   !                     cwtl2r(:,:,i)=cwtl2rtmp(:,:)	    
    !                        pathcoreval=exp(tmp2)
    !                    endif       
    !                else ! imp                      
    !                    if (rn1(1)<mov1ratio) then 
    !                        icmov1bisect1tot=icmov1bisect1tot+1
    !                        gauss=reshape(gaussian(3*npart),(/3,npart/))
    !                        sigmamid=sigmamid1
    !                        xold=(ristra(i-1)%x+ristra(i+1)%x)/2  ! ref point
    !                        xnew=xold+sigmamid*gauss     
    !                        xold1=2*xold-ristra(i)%x 
    !                        xnew1=xold-sigmamid*gauss
    !                        call v6proplr(xnew,xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmpnew(:,:))
    !                        call v6proplr(xold1,xold1,-1,cwtl2r(:,:,i-1),cwtl2rtmpold1(:,:))
    !                        call v6proplr(xnew1,xnew1,-1,cwtl2r(:,:,i-1),cwtl2rtmpnew1(:,:))                 
    !                        tmp1=pathcoreval             
	   !                     tmp2=abs(real(sum(conjg(cwtl2rtmpnew(:,:))*cwtr2l(:,:,i+1))))                                  
	   !                     tmp3=abs(real(sum(conjg(cwtl2rtmpold1(:,:))*cwtr2l(:,:,i+1))))                                                  
	   !                     tmp4=abs(real(sum(conjg(cwtl2rtmpnew1(:,:))*cwtr2l(:,:,i+1))))                                      
    !                        q1=(tmp2+tmp4)/(tmp1+tmp3)
    !                        q2=tmp2/(tmp2+tmp4)                                
    !                    else            
    !                        icmov1tot=icmov1tot+1
    !                        gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    !                        xold=ristra(i)%x 
    !                        xnew(:,:)=xold(:,:)+gauss(:,:)  
    !                        xold1=2*xnew-xold
    !                        xnew1=xold-gauss
    !                        xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2+(ristra(i+1)%x(:,:)-xnew(:,:))**2
    !                        xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2+(ristra(i+1)%x(:,:)-xold(:,:))**2
    !                        xd2new1=(ristra(i-1)%x-xnew1)**2+(ristra(i+1)%x-xnew1)**2
    !                        xd2old1=(ristra(i-1)%x-xold1)**2+(ristra(i+1)%x-xold1)**2    
    !                        gnew=exp(gc*sum(xd2new))
    !                        gold=exp(gc*sum(xd2old))
    !                        gnew1=exp(gc*sum(xd2new1))
    !                        gold1=exp(gc*sum(xd2old1)) 
    !                        call v6proplr(xnew,xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmpnew(:,:)) 
    !                        call v6proplr(xold1,xold1,-1,cwtl2r(:,:,i-1),cwtl2rtmpold1(:,:)) 
    !                        call v6proplr(xnew1,xnew1,-1,cwtl2r(:,:,i-1),cwtl2rtmpnew1(:,:))                             
    !                        tmp1=pathcoreval                      
	   !                     tmp2=abs(real(sum(conjg(cwtl2rtmpnew(:,:))*cwtr2l(:,:,i+1))))                      
	   !                     tmp3=abs(real(sum(conjg(cwtl2rtmpold1(:,:))*cwtr2l(:,:,i+1))))                                            
	   !                     tmp4=abs(real(sum(conjg(cwtl2rtmpnew1(:,:))*cwtr2l(:,:,i+1))))                                                   
    !                        q1=(gnew*tmp2+gnew1*tmp4)/(gold*tmp1+gold1*tmp3)
    !                        q2=gnew*tmp2/(gnew*tmp2+gnew1*tmp4)                                  
    !                    endif                         
	   !                 rn=randn(1)
	   !                 if (rn(1) < q1) then 
    !                        icnow=icnow+1 
    !                        if (rn1(1)<mov1ratio) then 
	   !                         icmov1bisect1=icmov1bisect1+1	
    !                        else
    !                            icmov1=icmov1+1   
    !                        endif                     
    !                        rn=randn(1)
    !                        if (rn(1) < q2) then
	   !                         ristra(i)%x(:,:)=xnew(:,:)
	   !                         cwtl2r(:,:,i)=cwtl2rtmpnew(:,:)   
    !                            pathcoreval=tmp2                             
    !                        else
	   !                         ristra(i)%x=xnew1
	   !                         cwtl2r(:,:,i)=cwtl2rtmpnew1(:,:)	   
    !                            pathcoreval=tmp4                                   
    !                        endif 
    !                    endif                             
    !                endif                
    !                if ( icnow /= 0 ) then
    !            ! update l2r
    !                    if (i.le.(nchorizo-2)) then
		  !              x=ristra(i+1)%x
		  !              call v6proplr(x,x,-1,cwtl2r(:,:,i),cwtl2r(:,:,i+1))
    !                    else ! i=nchorizo-1
    !                    xr=ristra(i+1)%x
    !                    call v6propl(xr,-1,cwtl2r(:,:,i),cwtl2r(:,:,i+1))	    
    !                    endif               
    !                endif   
	   !         endif 
	   !     endif 	 
    !    enddo 
    !    if (icnow /= 0) then 
    !        xl=ristra(0)%x
    !        xr=ristra(nchorizo)%x
    !        ! update cwtl2r    
    !        if (iright <= nchorizo-3) then 
    !            do k=iright,nchorizo-3
    !                x=ristra(k+2)%x
    !                call v6proplr(x,x,-1,cwtl2r(:,:,k+1),cwtl2r(:,:,k+2))   
    !            enddo
    !            call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))	      
    !        else    
    !            if (iright == nchorizo-2) then   
    !                call v6propl(xr,-1,cwtl2r(:,:,iright+1),cwtl2r(:,:,nchorizo))	          
    !            endif  
    !        endif    
    !        ! update cwtr2l from iright th bead until the 0th.
    !        if (iright == nchorizo) then    
    !            !x=ristra(iright)%x   
	   !         call v6propr(xr,-1,cwtend(:,:),cwtr2l(:,:,iright))
	   !         do k=iright-1,1,-1
	   !         x=ristra(k)%x
	   !         call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
	   !         enddo
	   !         !xl=ristra(0)%x
	   !         call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	    
    !        else      
    !            if (iright >= 1) then  
    !                do k=iright,1,-1
	   !                x=ristra(k)%x
	   !                call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
	   !             enddo
	   !             !xl=ristra(0)%x
	   !             call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	   
    !            else ! iright = 0
    !                call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	   
    !            endif  
    !        endif
    !    endif   
    !case default 
    !    write (6,*) 'arrow for stdmove1bisect1flex has not been set, abort!', arrow
	   ! call abort
    !end select
    !return
    !end subroutine stdmove1bisect1v6flexa              
    !
    !
    !subroutine cwt1test(ntab,rangev,it,i,j,cwtold,cwtnew,val,u,evdt) ! just for a particular pair i,j
    !use v6stepcalc
    !integer(kind=i4) :: i,j,it,index,ntab	
    !real(kind=r8) :: x(3,npart) ! the configuration for the system
    !real(kind=r8) :: dx(3)
    !real(kind=r8) :: r,c1,c2,c3,c4,dr,val,e1,e5! tstep=dt*2**it
    !real(kind=r8) :: u(6),evdt(6),scalep,rangev,u1(6),u2(6),evdtinv(6)
    !complex(kind=r8) :: cwtnew(0:nspin-1,nisospin),cwtold(0:nspin-1,nisospin),cwt1(0:nspin-1,nisospin)
    !real(kind=r8) :: utab(6,0:ntab,-mmax:mmax),evdttab(6,0:ntab,-mmax:mmax) 
    !real(kind=r8) :: vtab(6,0:ntab)
    !
    !
    !
    !call getutab(utab)
    !
    !
    !call getevdttab(evdttab)  
    !
    !
    !x=ristra(0)%x
    !
    !cwt1(:,:)=cwtold(:,:)
	   !
    !!write(6,*) 'x=',x
    !
    !        dx(:)=x(:,i)-x(:,j)		   
		  !  r=sqrt(dot_product(dx,dx)) !r=sqrt(sum(dx(:)**2))
		  ! 
		  !  if (r.eq.0) then
		  !      write(6,*) 'two particles have the same position, stop'
			 !   stop
		  !  endif
	   !
		  !  scalep=ntab/rangev
		  ! 
    !        if (r.lt.rangev) then
			 !   dr=scalep*r
			 !   index=dr
			 !   index=max(1,min(index,ntab-2))
			 !   dr=dr-index
    !            c1=-dr*(dr-1.0_r8)*(dr-2.0_r8)/6.0_r8
    !            c2=(dr+1.0_r8)*(dr-1.0_r8)*(dr-2.0_r8)/2.0_r8
    !            c3=-(dr+1.0_r8)*dr*(dr-2.0_r8)/2.0_r8
    !            c4=(dr+1.0_r8)*dr*(dr-1.0_r8)/6.0_r8
    !            u(:)=c1*utab(:,index-1,it)+c2*utab(:,index,it) &
			 !       +c3*utab(:,index+1,it)+c4*utab(:,index+2,it)
			 !   evdt(:)=c1*evdttab(:,index-1,it)+c2*evdttab(:,index,it) &
			 !       +c3*evdttab(:,index+1,it)+c4*evdttab(:,index+2,it)
			 !
    !          
    !            write(6,'(''index r dr='',t50, i20, 2(f15.7,1x) )') index,r,dr
    !            write(6,'(''c1-4 ='',t50,4(e27.20,1x) )') c1,c2,c3,c4
    !          
    !            call vtableout(vtab)
    !          
    !            write(6,'(''vtab='',t50,6(e27.20,1x) )') vtab(:,index)
    !            write(6,'(''evdt='',t50,6(e27.20,1x) )') evdt
    !            write(6,'(''1/evdt='',t50,6(e27.20,1x) )') 1./evdt 
			 !   write(6,'(''lagrange u='',t50,6(e27.20,1x) )') u
    !          
    !    u1(1)=(6.*evdt(1)+3.*evdt(2)+2.*evdt(3)+evdt(4)+3.*evdt(5)+evdt(6))/16.
	   ! u1(2)=(6.*evdt(1)+3.*evdt(2)+2.*evdt(3)+evdt(4)-9.*evdt(5)-3.*evdt(6))/48.
	   ! u1(3)=(3.*evdt(1)-3.*evdt(2)+evdt(3)-evdt(4))/24.
	   ! u1(4)=(2.*evdt(1)+evdt(2)-2.*evdt(3)-evdt(4)+evdt(5)-evdt(6))/16.
	   ! u1(5)=(2.*evdt(1)+evdt(2)-2.*evdt(3)-evdt(4)-3.*evdt(5)+3.*evdt(6))/48.
	   ! u1(6)=(evdt(1)-evdt(2)-evdt(3)+evdt(4))/24.
    !   
    !    evdtinv=1./evdt 
    !    u2(1)=(6.*evdtinv(1)+3.*evdtinv(2)+2.*evdtinv(3)+evdtinv(4)+3.*evdtinv(5)+evdtinv(6))/16.
	   ! u2(2)=(6.*evdtinv(1)+3.*evdtinv(2)+2.*evdtinv(3)+evdtinv(4)-9.*evdtinv(5)-3.*evdtinv(6))/48.
	   ! u2(3)=(3.*evdtinv(1)-3.*evdtinv(2)+evdtinv(3)-evdtinv(4))/24.
	   ! u2(4)=(2.*evdtinv(1)+evdtinv(2)-2.*evdtinv(3)-evdtinv(4)+evdtinv(5)-evdtinv(6))/16.
	   ! u2(5)=(2.*evdtinv(1)+evdtinv(2)-2.*evdtinv(3)-evdtinv(4)-3.*evdtinv(5)+3.*evdtinv(6))/48.
	   ! u2(6)=(evdtinv(1)-evdtinv(2)-evdtinv(3)+evdtinv(4))/24.       
    !          
    !            write(6,'(''assume evdt is correct true u='',t50,6(e27.20,1x) )') u1
    !            write(6,'(''difference btw true u and u is='',t50,6(e27.20,1x) )') u1-u
    !            write(6,'(''true u for -dt the u2 ='',t50,6(e27.20,1x) )') u2
    !          
    !           
    !           
    !            e1=evdt(1)
    !            e5=evdt(5)
    !            val=(3./e1+1./e5)*(3.*e1+e5)/16.+(3./e1-3./e5)*(e1-e5)/16.
    !          
    !            write(6,*) 'check 1= ',val
    !           
    !          
		  !  else
			 !   write(6,*) 'Damn!!!!!!!!!!!!!!!!!!' 
			 !   u=0.
    !            u1=0.
			 !   stop
		  !  endif
		  !  call v6prop1(x,i,j,u1,cwt1,cwtnew)
    !return
    !end subroutine cwt1test    
    !    

    !
    !subroutine stdmove1rangefastv6(ileft,iright) ! this is not bisection, just propose a move, and not efficient.
    !use wavefunction
    !use estimator
    !use v6stepcalc
    !use random
    !use brut
    !integer(kind=i4) :: i,k,ileft,iright,istep
    !real(kind=r8) :: rn(1),gauss(3,npart)
    !real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart),xd2old(3,npart),xd2new(3,npart) &
	   !             ,xr(3,npart),xl(3,npart)
    !real(kind=r8) :: rt2,tmp1,tmp2,logratio
    !complex(kind=r8) :: cwtold(0:nspin-1,nisospin),cwtnew(0:nspin-1,nisospin)
    !complex(kind=r8) :: cwtl2rtmp(0:nspin-1,nisospin),cwtr2ltmp(0:nspin-1,nisospin)
    !
    !ristraold(0:nchorizo)=ristra(0:nchorizo)
    !icstdmove1tot=icstdmove1tot+1
    !if ((ileft.lt.0).or.(iright.gt.nchorizo).or.(iright.lt.0).or.(ileft.gt.nchorizo)) then
	   ! write (6,*) 'stdmove1 range error', ileft,iright
	   ! stop
    !else
	   ! if ( ileft.le.iright ) then 
		  !  istep=1
	   ! else
		  !  istep=-1
	   ! endif 
    !endif
    !do i=ileft,iright,istep
    !  
    !    !write(6,*) 'i=',i
    !   
	   ! icmov1tot=icmov1tot+1  
	   ! xold=ristraold(i)%x 
    !! ran move	 
	   ! !rn3=reshape(randn(3*npart),shape(rn3)) 
	   ! !xnew=xold+mov1step*(rn3-0.5_r8) 
    !! gaussian move	 
	   !
    !    gauss(:,:)=mov1step*reshape(gaussian(3*npart),(/3,npart/))
    ! 
    !    !gauss=mov1step*reshape((randn(3*npart)-0.5),(/3,npart/)) ! this is not gauss move. for vmc(dt=0) only
    ! 
	   ! xnew(:,:)=xold(:,:)+gauss(:,:)  
	   !
    !    if (i.eq.0) then
		  !
	   ! !call psitcwt(xold,cwtold)
	   !
	   ! cwtold(:,:)=cwtbegin(:,:)
	   !
	   ! tmp1=log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i)))))
	   ! !call psitcwt(xnew,cwtnew)  
	   ! call corop(xnew,iplo,jplo,cwtnew)
	   ! call v6propl(xnew,-1,cwtr2l(:,:,i+1),cwtr2ltmp(:,:))
	   ! tmp2=log(abs(real(sum(conjg(cwtnew(:,:))*cwtr2ltmp(:,:)))))
	   ! xd2new(:,:)=(ristra(i+1)%x(:,:)-xnew(:,:))**2
	   ! xd2old(:,:)=(ristra(i+1)%x(:,:)-xold(:,:))**2
	   !
	   ! rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	  
    !    !rt2=0. ! vmc(dt=0) only.
	   !
	   ! logratio=tmp2-tmp1+rt2
	   ! rn=randn(1)
    !
	   ! if (log(rn(1)).lt.logratio) then ! dt ne 0
    !     
	   ! icmov1=icmov1+1	
	   ! ristra(i)%x(:,:)=xnew(:,:)
	   ! cwtr2l(:,:,i)=cwtr2ltmp(:,:)
    !! update cwtbegin
    !    cwtbegin(:,:)=cwtnew(:,:)  
    !! update cwtl2r	   
	   ! call v6propr(xnew,-1,cwtnew(:,:),cwtl2r(:,:,i))
	   ! do k=i+1,nchorizo-1
	   ! x=ristra(k)%x
	   ! call v6proplr(x,x,-1,cwtl2r(:,:,k-1),cwtl2r(:,:,k))
	   ! enddo
	   ! xr=ristra(nchorizo)%x
	   ! call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))	     
	   ! !write(6,*) 'move succeed 0',i
    ! 
	   ! endif 
	   !
	   ! else
		  !
	   ! if (i.eq.nchorizo) then   
	   ! !call psitcwt(xold,cwtold)
	   ! cwtold(:,:)=cwtend(:,:)
	   ! tmp1=log(abs(real(sum(conjg(cwtl2r(:,:,i))*cwtold(:,:)))))
	   ! !call psitcwt(xnew,cwtnew)
	   ! call corop(xnew,ipro,jpro,cwtnew)
	   ! call v6propl(xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmp(:,:))
	   ! tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtnew(:,:)))))
	   ! xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2
	   ! xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2
	   !
	   ! rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	  
    !    !rt2=0. ! vmc(dt=0) only.
	   !
	   ! logratio=tmp2-tmp1+rt2
	   ! rn=randn(1)
	   !
	   ! if (log(rn(1)).lt.logratio) then
	   ! icmov1=icmov1+1	
	   ! ristra(i)%x(:,:)=xnew(:,:)
	   ! cwtl2r(:,:,i)=cwtl2rtmp(:,:)
    !! update cwtend
    !    cwtend(:,:)=cwtnew(:,:)	   
    !! update cwtr2l	   
	   ! call v6propr(xnew,-1,cwtnew(:,:),cwtr2l(:,:,i))
	   ! do k=nchorizo-1,1,-1
	   ! x=ristra(k)%x
	   ! call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
	   ! enddo
	   ! xl=ristra(0)%x
	   ! call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	    
	   ! !write(6,*) 'move succeed 0',i
	   ! endif 		     
	   !
	   ! else
	   !
    !    cwtold(:,:)=cwtl2r(:,:,i)
    !
	   ! call v6proplr(xnew,xnew,-1,cwtl2r(:,:,i-1),cwtl2rtmp(:,:))
	   !
	   ! tmp2=log(abs(real(sum(conjg(cwtl2rtmp(:,:))*cwtr2l(:,:,i+1)))))
	   ! tmp1=log(abs(real(sum(conjg(cwtold(:,:))*cwtr2l(:,:,i+1)))))
	   ! 
	   ! xd2new(:,:)=(ristra(i-1)%x(:,:)-xnew(:,:))**2+(ristra(i+1)%x(:,:)-xnew(:,:))**2
	   ! xd2old(:,:)=(ristra(i-1)%x(:,:)-xold(:,:))**2+(ristra(i+1)%x(:,:)-xold(:,:))**2
	   !
	   ! rt2=-0.5/(dt*(2.*hbar))*(sum(xd2new)-sum(xd2old))	
    !    !rt2=0.  ! vmc(dt=0) only.
	   !
	   ! logratio=tmp2-tmp1+rt2	 
	   ! rn=randn(1)
    !
	   ! if (log(rn(1)).lt.logratio) then
    !       
    !        !write(6,*) 'chorizo i updated!'
    !       
	   ! icmov1=icmov1+1	
	   ! ristra(i)%x(:,:)=xnew(:,:)
	   ! cwtl2r(:,:,i)=cwtl2rtmp(:,:)	
    !
    !    xr=ristra(nchorizo)%x
	   ! xl=ristra(0)%x
    !! update l2r
    !    if (i.le.(nchorizo-2)) then
	   ! do k=i+1,nchorizo-1
		  !  x=ristra(k)%x
		  !  call v6proplr(x,x,-1,cwtl2r(:,:,k-1),cwtl2r(:,:,k))
	   ! enddo
	   ! call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))	
    !    else ! i=nchorizo-1
    !    call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))	    
    !    endif     
    !! update r2l      
	   ! do k=i,1,-1
    !        x=ristra(k)%x
		  !  call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))			
	   ! enddo
	   ! call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))	
	   ! endif    
	   ! endif 
	   ! endif 	 
    !enddo 
    !!if (mod(icstdmove1tot,nstepdecor).eq.0) then
    !!!write (6,*) 'icmov1tot=',icmov1tot
    !!ice=ice+1
    !!call hpsi(ristra(0))
    !!eleft=ristra(0)%elocal  
    !!call hpsi(ristra(nchorizo))
    !!eright=ristra(nchorizo)%elocal 
    !!select case (irep)
    !!case (1:2)
    !!   call addval(3,eleft,1.0d0)
    !!   call addval(4,eright,1.0d0)
	   ! !call funcrmsq(ristra(nchorizomid)%x,rm,rmsq)
	   ! !call addval(5,rmsq,1.0d0)
	   ! !write(11,'(2x,3(G12.5,2x))') eleft,eright,(eleft+eright)/2.,rm!11 means write to 'values.txt'
    !!case (3:4)
    !!   call addval(4,eleft,1.0d0)
    !!   call addval(5,eright,1.0d0)
	   ! !write(11,'(2x,3(G12.5,2x))') eleft,eright,(eleft+eright)/2.!11 means write to 'values.txt'
    !!end select
    !!endif 
    !return
    !end subroutine stdmove1rangefastv6         
    !    
    
    !
    !subroutine bisectv6improvew11 ! for v6 interaction with move 1 by 1.
    !use wavefunction
    !use estimator
    !use random
    !use v6stepcalc
    !integer(kind=i4) :: ileftbisect,il,ir,i,imid,j,jp,jd,jmax,l,lmax &
				!	    ,irightbisect,k,bilvl,n1,n2 &
	   !                 ,dtexpt
    !real(kind=r8) :: rn(1),tau,sigmamid,gauss(3,npart)
    !real(kind=r8) :: xmid(3,npart),xr(3,npart),xl(3,npart)
    !real(kind=r8) :: x(3,npart),xnew(3,npart),xold(3,npart)
    !real(kind=r8) :: tmp 
    !complex(kind=r8) :: oldlv(0:mmax),newlv(0:mmax)
    !complex(kind=r8) :: cwtl2rtmp1(0:nspin-1,nisospin,0:nbisect),cwtr2ltmp1(0:nspin-1,nisospin,0:nbisect) &
	   !                 ,cwtl2rtmp2(0:nspin-1,nisospin,0:nbisect),cwtr2ltmp2(0:nspin-1,nisospin,0:nbisect) &
	   !                 ,cwtbegintmp(0:nspin-1,nisospin),cwtendtmp(0:nspin-1,nisospin)
    !logical :: reject	
    !icbisecttot=icbisecttot+1
    !if (nbisect.lt.nchorizo) then
    !call bisectpicksliceswide(ileftbisect,irightbisect) ! pick the ileftbisect here
    !
    !!write(6,*) 'ileftbisect=', ileftbisect
    !!write(12,*) 'ileftbisect=', ileftbisect  
    !
    !if (ileftbisect.le.(nchorizo-nbisect)) then
    !ibisecttot=ibisecttot+1 
    !
    !!write(6,*) 'bisection is running'
    !
    !!oldristra(0:nbisect)=ristra(ileftbisect:ileftbisect+nbisect)
    !
    !do i=0,nbisect
	   ! oldristra(i)=ristra(ileftbisect+i)
    !enddo
    !
    !newristra(0)=oldristra(0)
    !newristra(nbisect)=oldristra(nbisect)
    !
    !!do i=0,nbisect
    !!xoldtmp(i,:,:)=oldristra(i)%x
    !!enddo
    ! 
    !cwtbegintmp(:,:)=cwtl2r(:,:,ileftbisect)
    !cwtendtmp(:,:)=cwtr2l(:,:,ileftbisect+nbisect)
    !
    !! tmp1 deal with oldlv
    !cwtl2rtmp1(:,:,0:nbisect)=cwtl2r(:,:,ileftbisect:ileftbisect+nbisect)
    !cwtr2ltmp1(:,:,0:nbisect)=cwtr2l(:,:,ileftbisect:ileftbisect+nbisect)
    !! tmp2 deal with newlv   
    !cwtl2rtmp2(:,:,0:nbisect)=cwtl2r(:,:,ileftbisect:ileftbisect+nbisect)
    !cwtr2ltmp2(:,:,0:nbisect)=cwtr2l(:,:,ileftbisect:ileftbisect+nbisect)
    !
    !oldlv=1.
    !newlv=1.
    !
    !tau=dt*nbisect ! total time for the bisection slice.
    !do bilvl=1,mmax   ! level 1 to N. mmax = the total level N.
	   ! sigmamid=sqrt(hbar*tau*2.0d0**(-bilvl))  
	   ! !write(6,*) 'bilvl sigmamid=',bilvl,sigmamid
	   ! do i=0,2**(bilvl-1)-1 ! store all the possible trial config in new ristra%x
		  !  imid=2**(mmax-bilvl)+i*2**(mmax-bilvl+1)
		  !  il=imid-2**(mmax-bilvl)
		  !  ir=imid+2**(mmax-bilvl)
		  !  !write(6,*) 'bilvl il ir imid sigmamid=',bilvl, il,ir,imid,sigmamid
		  !  gauss(:,:)=sigmamid*reshape(gaussian(3*npart),(/3,npart/))
		  !  xmid(:,:)=(newristra(il)%x(:,:)+newristra(ir)%x(:,:))/2
		  !  newristra(imid)%x(:,:)=xmid(:,:)+gauss(:,:)
		  ! 
		  !  !xold=oldristra(imid)%x
		  !  !xnew=newristra(imid)%x
    !
    !!       if ( sum((newristra(imid)%x-oldristra(imid)%x)**2).eq.0) then
		  !  !  write(6,*) '-------------------------'
		  !  !  write(6,*) 'bilvl=',bilvl
		  !  !  write(6,*) 'imid,il,ir=',imid,il,ir
		  !  !  write(6,*) 'sigmamid=',sigmamid
		  !  !  write(6,*) 'gauss=',gauss
		  !  !     write(6,*) 'xnew=',newristra(imid)%x
		  !  !  write(6,*) 'xold=',oldristra(imid)%x
    !!           write(6,*) 'oldristra imid x=',oldristra(imid)%x
		  !  !  write(6,*) '-------------------------'
		  !  !
		  !  !  
    !!do j=0,nbisect
    !!write(6,*) 'j=',j  
    !!write(6,*) 'xold ristra diff=', oldristra(j)%x-xoldtmp(j,:,:)
    !!enddo			  
		  !  !  
		  !  !  stop
		  !  !  endif		   	   
		  !  !call hpsi(newristra(imid))
		  !  !call hpsi(oldristra(imid)) ! put here for save. But is not really needed if do not move beads before or after bisection
	   ! enddo 
    !enddo  
    !! up to here, newristra has been loaded.
    !
    !! judge if we accept the newristra or not.
    !! in this v6 case, v is not justpotential v, it should be the value of < ... >. lv.
		  !  !call hpsi0(newristra(imid))
		  !  !tmp=(0.5d0*(newristra(il)%v+newristra(ir)%v)-newristra(imid)%v &
    !    !          -0.5d0*(oldristra(il)%v+oldristra(ir)%v)+oldristra(imid)%v)*tau/(2.0d0**bilvl) 	      
    !
    !all: do bilvl=1,mmax  ! besection.  level 1 to N. mmax = the total level N.	   
	   ! bisecttot(bilvl)=bisecttot(bilvl)+1.	   
	   ! lmax=2**bilvl-1 
	   ! jd=2**(mmax-bilvl) ! interval
	   ! dtexpt=mmax-bilvl-1 ! for the R L part, that's why the extra -1. dt=2**dtexp
	   ! do l=1,lmax
		  !  j=l*jd
		  !  jp=(l-1)*jd
		  !  xold=oldristra(j)%x
    !        xnew(:,:)=newristra(j)%x(:,:)
		  !  call v6proplr(xold,xold,dtexpt,cwtl2rtmp1(:,:,jp),cwtl2rtmp1(:,:,j))    ! should be lr, not rl I believe.
		  !  call v6proplr(xnew,xnew,dtexpt,cwtl2rtmp2(:,:,jp),cwtl2rtmp2(:,:,j))  	   
		  !  !if ( sum((xnew-xold)**2).eq.0) then
			 !   !write(6,*) 'bilvl=',bilvl	
			 !   !write(6,*) 'l=',l,j,jp
			 !   !write(6,*) 'sigmamid=',sigmamid
		  !  !   write(6,*) 'xnew-xold=',xnew-xold
		  !  !endif   
	   ! enddo
	   ! jmax=lmax*jd
	   ! newlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp2(:,:,jmax))*cwtendtmp(:,:)) )) 	   
	   ! oldlv(bilvl)=abs(real( sum(conjg(cwtl2rtmp1(:,:,jmax))*cwtendtmp(:,:)) ))	   
	   ! tmp=log(newlv(bilvl))-log(oldlv(bilvl))+log(oldlv(bilvl-1))-log(newlv(bilvl-1))   
  	 !   !tmp=(newlv(bilvl)/oldlv(bilvl))*(oldlv(bilvl-1)/newlv(bilvl-1))
	   ! rn=randn(1)	   
    !! can check if l2r and r2l give	the same lv.		    
	   !     !write(6,*) 'tmp ln rn=',tmp,log(rn(1))
	   ! if (log(rn(1)).lt.tmp) then
	   !     reject=.false. 
			 !
			 !   !if ( sum((xnew-xold)**2).eq.0) then
			 !   !
			 !   !write(6,*) 'bilvl=',bilvl	
			 !   !write(6,*) 'sigmamid=',sigmamid
		  !  !   write(6,*) 'xnew-xold=',xnew-xold
    !    !         write(6,*) 'tmp ln(rn)=',tmp,log(rn(1))		   
		  !  !   write (6,*) 'lvnew-lvold=',newlv(bilvl)-oldlv(bilvl)
			 !   !write (6,*) 'previous lvnew-lvold=',newlv(bilvl-1)-oldlv(bilvl-1)		  
				!    ! 
		  !  !
			 !   !endif  
				!
		  !  !if (mod(icstepstd,5*nstepdecor).eq.0) then   
		  !  ! !if (tmp.gt.300.) then
		  !  !  !write(6,*) 'bug bead appear! exit!'
		  !  !     write(6,*) 'ileft bilvl il ir imid sigmamid=',ileftbisect,bilvl,il,ir,imid,sigmamid	   
	   ! !        write(6,*) 'tmp ln(rn)=',tmp,log(rn(1))		   
		  !  !     write (6,*) 'lv=',oldlv(bilvl),newlv(bilvl),oldlv(bilvl-1),newlv(bilvl-1)
		  !  !  write(6,*) 'cwtbegintmp+*cwtendtmp=',sum(conjg(cwtbegintmp(:,:))*cwtendtmp(:,:))
		  !  !  !write(6,*) 'cwtl2rtmp1=',cwtl2rtmp1(:,:,jmax),jmax
		  !  !  !write(6,*) 'xold=',xold	
    !!           !write(6,*) 'xnew=',xnew	
		  !  !  !write(6,*) 'cwtl2rtmp2=',cwtl2rtmp2(:,:,jmax)	
		  !  !  !write(6,*) 'cwtendtmp=',cwtendtmp
		  !  !
		  !  !     !stop
		  !  !    endif
	   ! else
	   !     reject=.true.  
	   !     exit all 
	   ! endif 
		  !  !write(6,*) 'imid=',imid
    !
	   ! bisectcount(bilvl)=bisectcount(bilvl)+1.
    !enddo all
    !if (reject) then
	   ! call addval(2,0.0_r8,1.0_r8)
    !else
	   ! ibisect=ibisect+1 
	   ! call addval(2,1.0_r8,1.0_r8)  
    !
	   ! do i=0,nbisect
    !        ristra(ileftbisect+i)=newristra(i)
    !    enddo	 
    !!ristra(ileftbisect:ileftbisect+nbisect)=newristra(0:nbisect) ! this look like a bad way (pointer issue I think), do loop is good.	 
	 	 !
    !! update cwtl2r, cwtr2l.	 
    !! update cwtl2r for all beads. this can be optimized in the future.
	   !
    !    xr=ristra(nchorizo)%x
    !    xl=ristra(0)%x
    ! 
    !    cwtl2r(:,:,ileftbisect+1:ileftbisect+nbisect-1)=cwtl2rtmp2(:,:,1:jmax) ! here jmax should already be 2**mmax-1 
    !    if ((ileftbisect+nbisect).le.nchorizo-1) then
    !    do k=ileftbisect+nbisect,nchorizo-1
	   ! x=ristra(k)%x
	   ! call v6proplr(x,x,-1,cwtl2r(:,:,k-1),cwtl2r(:,:,k))
	   ! enddo
	   ! call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))
    !    else
    !    call v6propl(xr,-1,cwtl2r(:,:,nchorizo-1),cwtl2r(:,:,nchorizo))    
    !    endif
    ! 
    ! 
    !! update cwtr2l for all the beads
    ! 
	   ! do k=ileftbisect+nbisect-1,1,-1
	   ! x=ristra(k)%x
	   ! call v6proplr(x,x,-1,cwtr2l(:,:,k+1),cwtr2l(:,:,k))
	   ! enddo
    !    call v6propl(xl,-1,cwtr2l(:,:,1),cwtr2l(:,:,0))  	
    !endif  
    !else
    !
    !        !write(6,*) 'bisection is not running'
    !   
	   ! ! this part is somewhat like reptate perhaps.
    !! need to modify here   
	   !
    !   
    !!!!--------------------    
    !   
	   ! if (ileftbisect.lt.nchorizo) then ! nbisect should be >= 2.
	   ! n1=nchorizo-ileftbisect
	   ! n2=nbisect-1-n1
	   ! !write (6,*) 'bisect improve', ileftbisect+1,nchorizo,n1,n2,n2-1
	   ! !call stdmove1range(ileftbisect+1,nchorizo)
    !    
    !    !call stdmove1rangefastv6(ileftbisect+1,nchorizo)	
    !    call stdmove1fastv6(ileftbisect+1,nchorizo)
    !    
	   ! if (n2.gt.0) then
	   ! !call stdmove1range(0,n2-1)
    !        
	   ! !call stdmove1rangefastv6(0,n2-1)
    !    call stdmove1fastv6(0,n2-1)
    !    
    !    endif
    !    else ! ileftbisect.eq.nchorizo
		  !  !write (6,*) 'bisect improve ileftbisect=', ileftbisect+1,nchorizo,n1,n2,n2-1
	   ! !call stdmove1range(0,nbisect-2)	
    !     
    !    !call stdmove1rangefastv6(0,nbisect-2)	
    !    call stdmove1fastv6(0,nbisect-2)
    !    endif	
    !   
    !!!!--------------------     
    ! 
    !endif
    !else
    !!call bisectnchorizo  
	   ! write (6,*) 'nbisect >= nchorizo, stop!'
	   ! stop 
    !endif
    !
    !
    !! do not need to calculate elocal too often.  
    !!if (mod(icbisecttot,nstepdecor).eq.0) then
    !!ice=ice+1
    !!call hpsi(ristra(0))
    !!eleft=ristra(0)%elocal  
    !!call hpsi(ristra(nchorizo))
    !!eright=ristra(nchorizo)%elocal 
    !!select case (irep)
    !!case (1:2)
    !!   call addval(3,eleft,1.0d0)
    !!   call addval(4,eright,1.0d0)
	   ! !call addval(5,(eleft+eright)/2.,1.0d0)
	   ! !call funcrmsq(ristra(nchorizomid)%x,rm,rmsq)
	   ! !call addval(6,rmsq,1.0d0)
	   ! !call addval(7,rm,1.0d0)
	   ! !write(11,'(2x,3(G12.5,2x))') eleft,eright,(eleft+eright)/2.,rm!11 means write to 'values.txt'
    !!case (3:4)
    !!   call addval(4,eleft,1.0d0)
    !!   call addval(5,eright,1.0d0)
	   ! !write(11,'(2x,3(G12.5,2x))') eleft,eright,(eleft+eright)/2.!11 means write to 'values.txt'
    !!end select 
    !!endif
    !return  
    !end subroutine bisectv6improvew11      
!
      

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
    !integer :: ipl(npair),jpl(npair),ipr(npair),jpr(npair)
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
    !   val1=sum(conjg(cwtl2r(:,:,nchorizo))*cwt(:,:))
    !   
    !   !write(6,*) 'cwt=', cwt(:,:)
    !   !write(6,*) 'cwtl2r=', cwtl2r(:,:,nchorizo)  
    !   !write(6,*) 'xl before=', xl
    !   
    !   call coropl(xl,cwt)
    !   
    !   !write(6,*) 'xl after=', xl
    !    
    !   val2=sum(conjg(cwtr2l(:,:,0))*cwt(:,:))
    !
    !   !write(6,*) 'cwt=', cwt(:,:)
    !   !write(6,*) 'cwtr2l=', cwtr2l(:,:,0)
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
    !   call hpsitcwt(xtmp,iplo,jplo,cwtr2l(:,:,0),psi20,e0,d20,pe0,f0,cwta1)
    !   
    !  ! write(6,*) 'xl after=', xl
    !   
    !   !write(6,*) 'cwtr2l=', cwtr2l(:,:,0)
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
    !   call hpsitcwt(xtmp,ipro,jpro,cwtl2r(:,:,nchorizo),psi2n,en,d2n,pen,fn,cwta)
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
    !gl=sum(psith0cwt(:,:)*cwtr2l(:,:,0))
    !rgl=real(gl)
    !gr=sum(conjg(cwtl2r(:,:,nchorizo))*nhpsitcwt(:,:))
    !rgr=real(gr)   
    !rg=0.5*(rgl+rgr)
    !
    !f=sum(cwtold(:,:)*cwtr2l(:,:,0))
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
    !!write(6,*) 'cwtr2l(0)=',cwtr2l(:,:,0)
    !!write(6,*) 'psith0cwt=',psith0cwt(:,:)
    !return
    !end subroutine hpsi       
    
   
   
    end module step
