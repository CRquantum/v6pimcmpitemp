!$Id: he4reptate.f90,v 1.5 2003/02/25 18:11:59 nuclear Exp $
   program v6pimc
   use random
   use random2
   use brut
   use mympi
   use step
   use v6stepcalc
   use estimator
   use estimatorristra
   use wavefunction
   use math
   use matrixmod
! perhaps make those constants a type, which can be used in all the module.

   implicit none
   integer, parameter :: i4=selected_int_kind(9) ! i4 range is about 10^-10 ~ 10^10.
   integer, parameter :: i8=selected_int_kind(15)
   integer, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i8) :: irn,irnsave
   integer(kind=i8), allocatable :: irnall(:)
   real(kind=r8) :: dummy   
   integer(kind=i4) :: npart,nprot,nneut,nspin,nisospin,nbasis,nstate,nstatenew,ispin,iispin
   integer(kind=i4) :: i,imax,imaxval,it,j,k,l,m,n &
	                  ,istep,nstep,neq,nav,nstepdecor,mmax,nbisect,ntab,nchorizo,irep,nrhobin &
	                  ,nchorizomid,ixtemp,nrepmax &
                      ,day,hour,minute
   real(kind=r8) :: rangev,rangerho,mov1step,x0step,mov2step,shftstep,rhobinsize,rangex
   real(kind=r8) :: val1,val2,u1(6),u2(6),evdt1(6),evdt2(6),val3,u3(6),evdt3(6)
   complex(kind=r8), allocatable :: cwtnew(:,:),cwtold(:,:),cwt1(:,:),cwt2(:,:),cwt1tmp(:,:),cwt2tmp(:,:)
   real(kind=r8), allocatable :: x(:,:),xtot(:,:,:),xr(:,:),xl(:,:),rhodist(:),rhodisterr(:)
   real(kind=r8) :: hbar,dt
   real(kind=r8), allocatable :: valnow(:),valaverage(:),valerr(:)
   real(kind=r8) :: valn,val,err,err1,err2,error
   real(kind=r8) :: time0,time1,time2,time00,timeall,second! standard requires these to be real w/o kind
   logical :: isite,icorrchk,iem,calcmode
   character(len=120), dimension(:), allocatable :: answer
   character(len=70) :: infile,outfile,ooo
   character(len=30) :: vpropchorizo
   integer(kind=i4) :: ipl(6),jpl(6),ipr(6),jpr(6)	
   integer(kind=i4), allocatable :: iplall(:),jplall(:),iprall(:),jprall(:)	! for mpi_gather
   real(kind=r8), allocatable :: xall(:) ! for mpi_gather
   real(kind=r8) :: randomtest(10)

    call init0 !mpi initialization must be done before reading    
    if (myrank().eq.0) then
          
    read (5,*) irn     !random number seed
    read (5,*) hbar    !hbar^2/2m
    read (5,*) dt      !dt the time step    
    read (5,*) npart   !particle number
    read (5,*) nprot   ! proton number
    read (5,*) ntab    ! table points for potential
    read (5,*) mmax    ! max level of bisection, nbisect=2**mmax
    read (5,*) nrepmax ! max # of beads for reptation
    !   read (5,*) i ! i input
    !   read (5,*) j ! j input
    !   read (5,*) nstate ! name the input state. 
    read (5,*) isite   !if true take initial
    read (5,*) neq     !equilibration blocks
    read (5,*) nav     !averaging blocks
    read (5,*) nstep   !steps per block
    read (5,*) nstepdecor   !calculate energy every this number steps
    read (5,*) nchorizo !number of chorizos, better be an even number
    read (5,*) mov1step ! move1 step. move beads 1 by1 step
    read (5,*) mov2step ! move2 step move all beads at the same time. For irep=4 the VMC, it is VMC move step delta.
    read (5,*) shftstep ! shift all beads at the same time
    read (5,*) x0step !the random step move for 0th and nth chorizo
    read (5,*) irep    !type of action. 2=pimc. 4=vmc. 5= calculation mode (consume less ram). 6=perhaps big memory calculation mode.
    read (5,*) rangev    ! potential range (fm)
    read (5,*) rangerho    ! one particle density rho range (fm)   
    read (5,*) nrhobin    ! the number of bins for one particle density rho.
    read (5,*) icorrchk    ! true=check correlation, false means not. 
    read (5,*) iem    ! true=EM force on.  
    read (5,'(a70)') infile
    read (5,'(a70)') outfile
    ! formpi, need to stike an 'Enter' to make one more line in the unit 5 file.
 
    infile=adjustl(infile)
    infile=infile(1:index(infile,' ')-1)    
    outfile=adjustl(outfile)
    outfile=outfile(1:index(outfile,' ')-1)   
   
    write (6,'(''MPI number of cores ='',t30,i20)') nproc()
    write (6,'(''random number seed ='',t30,i20)') irn
    write (6,'(''hbar^2/2m ='',t40,f10.5)') hbar
    write (6,'(''time step ='',t40,f15.7)') dt
    write (6,'(''npart ='',t40,i10)') npart
    write (6,'(''nproton ='',t40,i10)') nprot  
    write (6,'(''table size ='',t40,i10)') ntab   
    write (6,'(''bisection order ='',t40,i10)') mmax
    write (6,'(''reptation beads number ='',t40,i10)') nrepmax   
    write (6,'(''start from sites ='',t40,l10)') isite
    write (6,'(''equilibration blocks ='',t40,i10)') neq
    write (6,'(''averaging blocks ='',t40,i10)') nav
    write (6,'(''nsteps/blocks ='',t40,i10)') nstep
    write (6,'(''stepdecor ='',t40,i10)') nstepdecor
    write (6,'(''nchorizo ='',t40,i10)') nchorizo
    write (6,'(''stepmove1 ='',t40,f10.5)') mov1step
    write (6,'(''stepmove2 ='',t40,f10.5)') mov2step 
    write (6,'(''stepshift ='',t40,f10.5)') shftstep
    write (6,'(''x0step ='',t40,f10.5)') x0step
    write (6,'(''potential range ='',t40,f10.5)') rangev
    write (6,'(''one particle density rho range ='',t40,f10.5)') rangerho
    write (6,'(''bin # of rho ='',t40,i10)') nrhobin
    write (6,'(''correlation check ='',t40,l10)') icorrchk
    write (6,'(''EM force ='',t40,l10)') iem
    select case (irep)
    case (1)
        write (6,'(''Primitive Action'')')
    case (2)
        write (6,'(''Local Energy Action: Path integegral QMC'')')
    case (3)
        write (6,'(''Local Energy Action with Metropolis'')')
    case (4)
        write (6,'(''Variational Calculation'')')
    case (5)
        write (6,'(''==============calculation mode=============='')')
    case (6)
        write (6,'(''==============calculation super mode=============='')')
    case default
        write (6,'(''Illegal irep value: '')') irep
        call abort
    end select

    !-----------------------------------------------------   
    endif
   
    call bcast(irn)     !random number seed
    call bcast(hbar)    !hbar^2/2m
    call bcast(dt)      !dt the time step    
    call bcast(npart)   !particle number
    call bcast(nprot)   ! proton number
    call bcast(ntab)    ! table points for potential
    call bcast(mmax)    ! max level of bisection, nbisect=2**mmax
    call bcast(nrepmax) ! max # of beads for reptation
    call bcast(isite)   !if true take initial
    call bcast(neq)     !equilibration blocks
    call bcast(nav)     !averaging blocks
    call bcast(nstep)   !steps per block
    call bcast(nstepdecor)   !calculate energy every this number steps
    call bcast(nchorizo) !number of chorizos, better be an even number
    call bcast(mov1step) ! move1 step
    call bcast(mov2step) ! move2 step move all beads at the same time
    call bcast(shftstep) ! shift all beads at the same time
    call bcast(x0step) !the random step move for 0th and nth chorizo
    call bcast(irep)    !type of action
    call bcast(rangev)    ! potential range (fm)
    call bcast(rangerho)    ! one particle density rho range (fm)   
    call bcast(nrhobin)    ! the number of bins for one particle density rho.
    call bcast(icorrchk)    ! true=check correlation, false means not. 
    call bcast(iem)    ! true=EM force on.   
    call bcast(infile)
    call bcast(outfile)
   
    select case (irep)
    case (1:2,5:6)
        call setestnum(11,irep)
        call addest(1,'reptation %')
        call addest(2,'bisection %')
        call addest(3,'numerator')
        call addest(4,'denominator|energy')
        call addest(5,'real f')	  
	    call addest(6,'<sqrt(<rm^2>)>') ! this one is not important, it is like < sqrt( < rm^2 > ) >.
	    call addest(7,'<rm>') ! this is the one really need this is the rm.
	    call addest(8,'rho (0th bin)')	
        call addest(9,'bisection half %')
        call addest(10,'bisection L end %')
        call addest(11,'bisection R end %')
            
        allocate(answer(11))

    case (3)
        call setestnum(5,irep)
        call addest(1,'reptation %')
        call addest(2,'bisection %')
        call addest(3,'diffusion %')
        call addest(4,'energyleft')
        call addest(5,'energyright')
        
        allocate(answer(5))
	  
    case (4)	! VMC 
        call setestnum(8,irep)
        call addest(1,'reptation %')
        call addest(2,'VMC %')
        call addest(3,'numerator')
        call addest(4,'denominator|energy')
        call addest(5,'real f')	  
	    call addest(6,'rm')
	    call addest(7,'<rm>')
	    call addest(8,'rho (0th bin)')	
      
        allocate(answer(8))	
      
    !case (5:6) ! pimc calculation mode
    !    call setestnum(7,irep)
    !    call addest(1,'numerator')
    !    call addest(2,'denominator')
    !    call addest(3,'energy')
    !    call addest(4,'real f')	  
	   ! call addest(5,'rm')
	   ! call addest(6,'<rm>')
	   ! call addest(7,'rho (0th bin)')	
    !    
    !    allocate(answer(7))
    case default
        write (6,'(''Illegal irep value: '')') irep
        call abort 
    end select
   
! select irep case   
    if (irep == 2) then ! pimc

    nbisect=2**mmax   
    nchorizomid=nint(dble(nchorizo)/2.)
    rhobinsize=rangerho/dble(nrhobin)       
      
    if (myrank().eq.0) then   
    if (nbisect.gt.nchorizo) then
	    write(6,*) '# of Bisection beads must be smaller than nchorizo! Stop!',nbisect,nchorizo
	    call abort
    endif
    if (nrepmax.gt.nchorizo) then  
	    write(6,*) '# of Reptation beads must be smaller than nchorizo! Stop!',nrepmax,nchorizo
	    call abort
    endif        
    call cleanfolder(irep)

    open (12,FILE='answer.txt',FORM='FORMATTED')
    write (12,'(''MPI # of cores ='',t30,i20)') nproc()
    write (12,'(''random number seed ='',t30,i20)') irn
    write (12,'(''hbar^2/2m ='',t40,f10.5)') hbar
    write (12,'(''time step ='',t40,f15.7)') dt
    write (12,'(''npart ='',t40,i10)') npart
    write (12,'(''nproton ='',t40,i10)') nprot  
    write (12,'(''table size ='',t40,i10)') ntab   
    write (12,'(''bisection order ='',t40,i10)') mmax
    write (12,'(''# of reptation beads ='',t40,i10)') nrepmax  
    write (12,'(''start from sites ='',t40,l10)') isite
    write (12,'(''equilibration blocks ='',t40,i10)') neq
    write (12,'(''averaging blocks ='',t40,i10)') nav
    write (12,'(''nsteps/blocks ='',t40,i10)') nstep
    write (12,'(''stepdecor ='',t40,i10)') nstepdecor
    write (12,'(''nchorizo ='',t40,i10)') nchorizo
    write (12,'(''stepmove1 ='',t40,f10.5)') mov1step
    write (12,'(''x0step ='',t40,f10.5)') x0step
    write (12,'(''potential range ='',t40,f10.5)') rangev
    write (12,'(''one particle density rho range ='',t40,f10.5)') rangerho
    write (12,'(''bin # of rho ='',t40,i10)') nrhobin
    write (12,'(''correlation check ='',t40,l10)') icorrchk
    write (12,'(''EM force ='',t40,l10)') iem

    if (isite) then
        write (12,'(''input from sites'')')
    else
        write (12,'(''input from file'',t40,a20)') infile
    endif   

    endif
    
    do i=0,nproc()-1
        call ran2(dummy,irn)
        if (myrank().eq.i) then
            irnsave=irn
            call setrn(irn)    ! different cores have different seed.  
            exit
        endif 
    enddo
   
    !allocate(irnall(0:nproc()-1))
    !call gather(irnsave,irnall)
    !if (myrank().eq.0) then
    !    open(unit=9,form='formatted',file='irnallcheck.txt')
    !    do i=0,nproc()-1
    !       write(9,*) i,irnall(i)
    !    enddo
    !    close(9)
    !endif
   
   
    call v6stepcalcinit(npart,nprot,dt,hbar,ntab,rangev,mmax,nrepmax,iem,irep) 
    call wavefunctioninit(nchorizo,npart,nprot,ntab,hbar,rangev,rangerho,nrhobin,irep)  
    call stepinit(npart,nchorizo,hbar,dt,mmax,nrepmax,irep,nstepdecor,x0step,mov1step,mov2step,shftstep,nrhobin &
        ,icorrchk,nprot,nav,nstep)
    call brutinit(npart,nprot,hbar,iem) ! 3 in av18pot.f is v6' . this gives the initial L,R order for psi_T.
    call stepchorizoinit(isite,infile,outfile)   ! call brutinit first.	
	  	  
    call zerest
    call setestnumristra(1,nchorizo,irep)
    call addestristra(1,'potential.txt')
    call zerestristra
    call zerrhodist
    it=0
    call counterinit 
    call zerepstepmon
          
    if (myrank().eq.0) time00=mpi_wtime()
    do i=1,nav+neq  !do blocks
        if (i.eq.neq+1) then 
        if (myrank().eq.0) then  
            write (6,'(/,''equilibration done'')')
		    write (6,'(/,''=================='')')
		    write (12,'(/,''equilibration done'')')
		    write (12,'(/,''=================='')')
        endif
            call zerest
            call zerestristra
		    call zerrhodist
            call zerepstepmon
            it=0
        endif
        it=it+1 
	    if ((i.eq.1).and.(myrank().eq.0)) time0=mpi_wtime()
	    call counterinit
        do j=1,nstep
		    call stepstd
        enddo
        if (myrank().eq.0) then
        write (6,'(/,''iteration ='',t30,i14)') it
        write (12,'(/,''iteration ='',t30,i14)') it
        call showstat 
        endif
        call update   !collect block averages
        call updateristra
        answer=resstring()
      
        if (myrank().eq.0) then
        write (6,'(a120)') (answer(k),k=1,size(answer))
	    write (12,'(a120)') (answer(k),k=1,size(answer))  
        endif
      
    ! write out x.      
    allocate(xall(3*npart*(nchorizo+1)*nproc()))
    allocate(iplall(6*nproc()),jplall(6*nproc()),iprall(6*nproc()),jprall(6*nproc()))   
    allocate(xtot(3,npart,0:nchorizo))
    call chorizoallout(xtot,ipl,jpl,ipr,jpr)  
    call gather(reshape(xtot,(/3*npart*(nchorizo+1)/)),xall)
    call gather(ipl,iplall)
    call gather(jpl,jplall)
    call gather(ipr,iprall)
    call gather(jpr,jprall)      
      
    if (myrank().eq.0) then
        open(unit=9,form='formatted',file=trim(outfile),position='rewind')
        write(9,'(i10)') nproc()
        ixtemp=3*npart*(nchorizo+1)
        do l=0,nproc()-1
            write (9,'(6i10)') iplall((6*l+1):(6*(l+1))) 
            write (9,'(6i10)') jplall((6*l+1):(6*(l+1))) 
            write (9,'(6i10)') iprall((6*l+1):(6*(l+1))) 
            write (9,'(6i10)') jprall((6*l+1):(6*(l+1))) 
            write (9,'(3e15.7)') xall((ixtemp*l+1):(ixtemp*(l+1)))    
        enddo 
        close(9)
    endif
    deallocate(xall,iplall,jplall,iprall,jprall,xtot)
     
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
    call result(4,valn,val2,err2) ! #4 is the denominator 
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
     
    if ((i.eq.1).and.(myrank().eq.0)) then
        time1=mpi_wtime()
	    write (6,'(/,''Time for one block ='',f10.3,'' minutes'')') (time1-time0)/60.
        write (6,'(/,''Total time estimation ='',f10.3,'' hours'')') dble(nav+neq)*(time1-time0)/3600.
        write (12,'(/,''Time for one block ='',f10.3,'' minutes'')') (time1-time0)/60.
        write (12,'(/,''Total time estimation ='',f10.3,'' hours'')') dble(nav+neq)*(time1-time0)/3600.
    endif

    enddo
   
   
    call barrier
    if (myrank().eq.0) then
    time2=mpi_wtime()
    write (6,'(/,''Total time ='',f10.3,'' hours'')') (time2-time00)/3600.
    write (12,'(/,''Total time ='',f10.3,'' hours'')') (time2-time00)/3600.
    close (12)
    endif
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    else if ((irep == 5).or.(irep == 6)) then ! pure pimc calculation mode
           
    nbisect=2**mmax   
    nchorizomid=nint(dble(nchorizo)/2.0_r8)
    rhobinsize=rangerho/dble(nrhobin)       

    if (myrank().eq.0) then 
        time00=mpi_wtime()
        open (12,FILE='answer_CalcMode.txt',FORM='FORMATTED')
        write (12,'(''MPI # of cores ='',t30,i20)') nproc()
        write (12,'(''random number seed ='',t30,i20)') irn
        write (12,'(''hbar^2/2m ='',t40,f10.5)') hbar
        write (12,'(''time step ='',t40,f15.7)') dt
        write (12,'(''npart ='',t40,i10)') npart
        write (12,'(''nproton ='',t40,i10)') nprot  
        write (12,'(''table size ='',t40,i10)') ntab   
        write (12,'(''bisection order ='',t40,i10)') mmax
        write (12,'(''# of reptation beads ='',t40,i10)') nrepmax  
        write (12,'(''start from sites ='',t40,l10)') isite
        write (12,'(''equilibration blocks ='',t40,i10)') neq
        write (12,'(''averaging blocks ='',t40,i10)') nav
        write (12,'(''nsteps/blocks ='',t40,i10)') nstep
        write (12,'(''stepdecor ='',t40,i10)') nstepdecor
        write (12,'(''nchorizo ='',t40,i10)') nchorizo
        write (12,'(''stepmove1 ='',t40,f10.5)') mov1step
        write (12,'(''x0step ='',t40,f10.5)') x0step
        write (12,'(''potential range ='',t40,f10.5)') rangev
        write (12,'(''one particle density rho range ='',t40,f10.5)') rangerho
        write (12,'(''bin # of rho ='',t40,i10)') nrhobin
        write (12,'(''correlation check ='',t40,l10)') icorrchk
        write (12,'(''EM force ='',t40,l10)') iem
        write (12,'(''input from file'',t40,a20)') infile
    endif

    if (myrank().eq.0) call cleanfolder(irep)
    
    call v6stepcalcinit(npart,nprot,dt,hbar,ntab,rangev,mmax,nrepmax,iem,irep) 
    call wavefunctioninit(nchorizo,npart,nprot,ntab,hbar,rangev,rangerho,nrhobin,irep)  
    call stepinit(npart,nchorizo,hbar,dt,mmax,nrepmax,irep,nstepdecor,x0step,mov1step,mov2step,shftstep &
        ,nrhobin,icorrchk,nprot,nav,nstep)
    call brutinit(npart,nprot,hbar,iem) ! 3 in av18pot.f is v6' . this gives the initial L,R order for psi_T.
    !call stepchorizoinit(isite,infile,outfile)   ! call brutinit first.	
    
    call zerest
    call setestnumristra(1,nchorizo,irep)
    call addestristra(1,'potential.txt')
    call zerestristra
    call zerrhodist
    call counterinit 
    call zerepstepmon
 
    call computeinit(answer,rhobinsize)
    call compute
    
    
    call barrier
    if (myrank().eq.0) then
    time2=mpi_wtime() ! for mpi code, mpi_wtime() is better than cpu_time.
    timeall=(time2-time00)
    call timedhms(timeall,day,hour,minute,second)
    write (6,'(/,''Total time ='',i10,'' days'',i10,'' hours'',i10,'' minutes'',f10.3,'' seconds'')') &
           day,hour,minute,second
    write (12,'(/,''Total time ='',i10,'' days'',i10,'' hours'',i10,'' minutes'',f10.3,'' seconds'')') &
           day,hour,minute,second
    close (12)     
    endif
    
       
       
       
       
       
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   
       
    else if (irep == 4) then    ! VMC calculations.	irep=4
	
    dt=0.
    nchorizo=max(2,nchorizo)
    mmax=max(1,mmax)
    nchorizomid=nint(dble(nchorizo)/2.)
    rhobinsize=rangerho/dble(nrhobin) 	
	  	  
    call setrn(irn+myrank())  ! different thread have different seed. 
    call v6stepcalcinit(npart,nprot,dt,hbar,ntab,rangev,mmax,nrepmax,iem,irep) 
    call wavefunctioninit(nchorizo,npart,nprot,ntab,hbar,rangev,rangerho,nrhobin,irep)  
    call stepinit(npart,nchorizo,hbar,dt,mmax,nrepmax,irep,nstepdecor,x0step,mov1step,mov2step,shftstep &
        ,nrhobin,icorrchk,nprot,nav,nstep)
    call brutinit(npart,nprot,hbar,iem) ! 3 in av18pot.f is v6' . this gives the initial L,R order for psi_T.
    call stepchorizoinit(isite,infile,outfile)   ! call brutinit first.	call stepinit first.   
   
    call zerest
    call setestnumristra(1,nchorizo,irep)
    call addestristra(1,'potential')
    call zerestristra
    it=0

    if (myrank().eq.0) then
    open(11,FILE='values.txt',FORM='FORMATTED')
    open(12,FILE='answer.txt',FORM='FORMATTED')
    write (12,'(''random number seed ='',t30,i20)') irn
    write (12,'(''hbar^2/2m ='',t40,f10.5)') hbar
    write (12,'(''time step ='',t40,f15.7)') dt
    write (12,'(''npart ='',t40,i10)') npart
    write (12,'(''nproton ='',t40,i10)') nprot  
    write (12,'(''table size ='',t40,i10)') ntab   
    write (12,'(''bisection order ='',t40,i10)') mmax
    write (12,'(''start from sites ='',t40,l10)') isite
    write (12,'(''equilibration blocks ='',t40,i10)') neq
    write (12,'(''averaging blocks ='',t40,i10)') nav
    write (12,'(''nsteps/blocks ='',t40,i10)') nstep
    write (12,'(''stepdecor ='',t40,i10)') nstepdecor
    write (12,'(''nchorizo ='',t40,i10)') nchorizo
    write (12,'(''stepmove1 ='',t40,f10.5)') mov1step
    write (12,'(''x0step ='',t40,f10.5)') x0step
    write (12,'(''potential range ='',t40,f10.5)') rangev
    write (12,'(''one particle density rho range ='',t40,f10.5)') rangerho
    write (12,'(''bin # of rho ='',t40,i10)') nrhobin
    write (12,'(''correlation check ='',t40,l10)') icorrchk
    write (12,'(''EM force ='',t40,l10)') iem
   
    endif
   
   
   
    if (isite) then
        if (myrank().eq.0) write (12,'(''input from sites'')')
    else
        if (myrank().eq.0) write (12,'(''input from file'',t40,a20)') infile
    endif
   
   
    if (myrank().eq.0) open(13,FILE='energy_evolution.txt',FORM='FORMATTED')
   
    call CPU_TIME(time00)
    do i=1,nav+neq  !do blocks
        if (i.eq.neq+1) then 
        if (myrank().eq.0) then  
            write (6,'(/,''equilibration done'')')
		    write (6,'(/,''=================='')')
		    write (12,'(/,''equilibration done'')')
		    write (12,'(/,''=================='')')
        endif
            call zerest
            call zerestristra
		    call zerrhodist
            it=0
        endif
        it=it+1 
	    if (i.eq.1) call CPU_TIME(time0)
	    call counterinit
        do j=1,nstep
            !call step1
            !write(6,*) 'j=',j
            !write(12,*) 'j=',j
		    call stepstd
	    enddo
	    if (myrank().eq.0) call showstat 
        call update   !collect block averages
        call updateristra
        answer=resstring()
      
        if (myrank().eq.0) then
            write (6,'(/,''iteration = '',t30,i14)') it    
            write (6,'(a120)') (answer(k),k=1,size(answer))
	        write (12,'(/,''iteration = '',t30,i14)') it
	        write (12,'(a120)') (answer(k),k=1,size(answer))  
        endif
      
    ! write out x.      
    allocate(xall(3*npart*(nchorizo+1)*nproc()))
    allocate(iplall(6*nproc()),jplall(6*nproc()),iprall(6*nproc()),jprall(6*nproc()))   
    allocate(xtot(3,npart,0:nchorizo))
    call chorizoallout(xtot,ipl,jpl,ipr,jpr)  
    call gather(reshape(xtot,(/3*npart*(nchorizo+1)/)),xall)
    call gather(ipl,iplall)
    call gather(jpl,jplall)
    call gather(ipr,iprall)
    call gather(jpr,jprall)      
      
    if (myrank().eq.0) then
        open(unit=9,form='formatted',file=trim(outfile),position='rewind')
        write(9,'(i10)') nproc()
        ixtemp=3*npart*(nchorizo+1)
        do l=0,nproc()-1
            write (9,'(6i10)') iplall((6*l+1):(6*(l+1))) 
            write (9,'(6i10)') jplall((6*l+1):(6*(l+1))) 
            write (9,'(6i10)') iprall((6*l+1):(6*(l+1))) 
            write (9,'(6i10)') jprall((6*l+1):(6*(l+1))) 
            write (9,'(3e15.7)') xall((ixtemp*l+1):(ixtemp*(l+1)))    
        enddo 
        close(9)
    endif
    deallocate(xall,iplall,jplall,iprall,jprall,xtot)
     
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
    call result(4,valn,val2,err2) ! #4 is the denominator 
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
     
    if ((i.eq.1).and.(myrank().eq.0)) then
        call CPU_TIME(time1)  
	    write (6,'(/,''Time for one block ='',f10.3,'' minutes'')') (time1-time0)/60.
        write (6,'(/,''Total time estimation ='',f10.3,'' hours'')') dble(nav+neq)*(time1-time0)/3600.
        write (12,'(/,''Time for one block ='',f10.3,'' minutes'')') (time1-time0)/60.
        write (12,'(/,''Total time estimation ='',f10.3,'' hours'')') dble(nav+neq)*(time1-time0)/3600.
    endif

    enddo
    
    call barrier
    if (myrank().eq.0) then
    close (13)
    call CPU_TIME(time2)
    write (6,'(/,''Total time ='',f10.3,'' hours'')') (time2-time00)/3600.
    write (12,'(/,''Total time ='',f10.3,'' hours'')') (time2-time00)/3600.
    close (12)
    close(11)
    endif	  
	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
    else
	write(6,*) 'irep does not make sense. nothing is done.'  
    endif
   
    call barrier   
    if (myrank().eq.0) write(6,*) 'The program end normally!'
    call done
    end program v6pimc
	

    

    
    
    
	
