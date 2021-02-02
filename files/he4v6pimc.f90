! Rong Chen, based on $Id: he4reptate.f90,v 1.5 2003/02/25 18:11:59 nuclear Exp $
   program v6pimc
   use mympi
   use random
   use random2
   use brut
   use step
   use v6stepcalc
   use estimator
   use estimatorristra
   use wavefunction
   use math
   !use matrixmod
! perhaps make those constants a type, which can be used in all the module.
   implicit none
   integer, parameter :: i4=selected_int_kind(9) ! i4 range is about 10^-10 ~ 10^10.
   integer, parameter :: i8=selected_int_kind(15)
   integer, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i8) :: irn 
   integer(kind=i4) :: npart,nprot,lpot
   integer(kind=i4) :: i,it,ishift,j,k,l &
	                  ,nstep,neq,nav,nstepdecor,mmax,nbisect,ntab,nchorizo,irep,nrhobin &
	                  ,nchorizomid,nrepmax &
                      ,day,hour,minute &
                      ,nest,nestristra,nrespest,ne0respbeads &
                      ,idummy
   real(kind=r8) :: rangev,rangerho,mov1step,x0step,mov2step,shftstep,rhobinsize
   real(kind=r8) :: hbar,dt
   real(kind=r8) :: time0,time1,time2,time3,time4,time5,time6,time7,time8,time9,time10,time11 &
                   ,time00,timeall,t1,t2,t3,t4,t5,t6,t7,t8,t3a,t3p,t6a,t6p,tcalc1,tmpiio1 &
                   ,second! standard requires these to be real w/o kind
   logical :: isite,icorrchk,iem,spinorbit,vtune,resume
   character(len=120), dimension(:), allocatable :: answer
   character(len=70) :: infile,outfile
   character(len=120) :: filename,file_id
   real(kind=r8) :: rbisect,rmov1,rmov2,rmov3,rmov23 
   
   call init0 !mpi initialization must be done before reading  
    
   vtune=.false.
    
    if (myrank().eq.0) then
        if (.not.vtune) then   
            read (5,*) irn     !random number seed
            read (5,*) hbar    !hbar^2/2m
            read (5,*) dt      !dt the time step    
            read (5,*) npart   !particle number
            read (5,*) nprot   ! proton number
            read (5,*) lpot    ! type of the potential, 3 is v6'. 112-114,122-124 are cheft.
            read (5,*) ntab    ! table points for potential
            read (5,*) mmax    ! max level of bisection, nbisect=2**mmax. min better >=3
            read (5,*) nrepmax ! max # of beads for reptation, 0 means off.
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
            read (5,*) resume  ! set resume mode. when resume, settings shoud be the same, including the MPI core numbers.
            read (5,*) rangev    ! potential range (fm)
            read (5,*) rangerho    ! one particle density rho range (fm)   
            read (5,*) nrhobin    ! the number of bins for one particle density rho.
            read (5,*) icorrchk    ! true=check correlation, false means not. 
            read (5,*) iem    ! true=EM force on.  
            read (5,*) rbisect ! 0.9 bisection upper bar (reptation if enabled is in it.)
            read (5,*) rmov1 ! move 1 by 1 upper bar
            read (5,*) rmov2 ! move all beads  upper bar 
            read (5,*) rmov3 ! shift beads rate upper bar
            read (5,*) rmov23 ! move all then shift beads upper bar
            read (5,'(a70)') infile
            read (5,'(a70)') outfile
 
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
                    write (6,'(''============== Calculation Mode (w angle ave response) ================'')')
                case (6)
                    write (6,'(''============== Calculation Fast Mode (wo angle ave response) ================'')')
                case default
                    write (6,'(''Illegal irep value: '')') irep
                    call abort
            end select 
        
            select case (lpot) ! check lpot
                case (1)
                    write (6,'('' Argonne v18` '')')
                case (2)
                    write (6,'('' Argonne v8` '')')
                case (3)
                    write (6,'('' Argonne v6` '')')
                case (4)
                    write (6,'('' Argonne v4` '')')
                case (5)
                    write (6,'('' Argonne vx` '')')
                case (6)
                    write (6,'('' Argonne v2` '')')
                case (7)
                    write (6,'('' Argonne v1` '')')
                case (112)
                    write (6,'('' ChEFT LO, R0=1.0 '') ')
                case (113)
                    write (6,'('' ChEFT NLO, R0=1.0 '') ')
                case (114)
                    write (6,'('' ChEFT N2LO, R0=1.0 '') ')
                case (122)
                    write (6,'('' ChEFT LO, R0=1.2 '')')
                case (123)
                    write (6,'('' ChEFT NLO, R0=1.2 '') ')
                case (124)
                    write (6,'('' ChEFT N2LO, R0=1.2 '') ') 
                case default
                    write (6,'(''Illegal lpot value: '')') lpot
                    call abort          
            end select  
                     
            select case (lpot)
                case (1:2,113:114,123:124)
                    spinorbit=.true.       
                case default
                    spinorbit=.false.
                end select  ! can just add switch in the input file.
                
            if (nprot<=1) then
               iem=.false.
               write (6,'('' # of proton <=1, iem (pp) off '')') 
            endif
            
            if (irep.eq.2) then ! PIMC special case
                if ((nchorizo <= 20).or.(mmax <= 2)) then
	                write (6,*) 'nchorizo < 8 or mmax <= 2, bisecion is not needed and will not be used ',nchorizo,mmax
                    rbisect=-1 ! - means turned off. 0.9 bisection upper bar (reptation if enabled is in it.)
                    rmov1=0.7 ! move 1 by 1 upper bar
                    rmov2=0.8 ! move all beads  upper bar 
                    rmov3=0.9 ! shift beads rate upper bar
                    rmov23=1.0 ! move all then shift beads upper bar
                    mmax=2 ! set it to 2 to make sure the arrays bounds are correct, like utab.
                endif  
            endif
     
        else
            call setreadin ! for VTune Amplifier only.
        endif
        ! for mpi, need to stike an 'Enter' to make one more line in the unit 5 file. I guess this is to add an end of the file mark.
        
        infile=adjustl(infile)
        infile=infile(1:index(infile,' ')-1)    
        outfile=adjustl(outfile)
        outfile=outfile(1:index(outfile,' ')-1) 
        
        ! parepare file names.
        call stepfilesinit(irep)
        call wavefunctionfilesinit(irep)
        call estimatorfilesinit(irep)
        
        ! consistence check. check list.
        if (nproc()>=20) then
            if (icorrchk) then
                write (6,*) 'number of core >= 20, correlation check is turned off! ', nproc()
                icorrchk=.false.
            endif
        endif 
        call checkresume(irep,resume,outfile,isite) 
        ! consistence check end.   
          
        write (6,'(''MPI number of cores ='',t30,i20)') nproc()
        write (6,'(''random number seed ='',t30,i20)') irn
        write (6,'(''hbar^2/2m ='',t40,f10.5)') hbar
        write (6,'(''time step ='',t40,g15.7)') dt
        write (6,'(''npart ='',t40,i10)') npart
        write (6,'(''nproton ='',t40,i10)') nprot 
        write (6,'(''potential number ='',t40,i10)') lpot
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
        write (6,'(''Spin-Orbit force ='',t40,l10)') spinorbit
        write (6,'(''bisection (reptation) upper bar='',t40,f10.5)') rbisect 
        write (6,'(''move 1 by 1 upper bar='',t40,f10.5)') rmov1 
        write (6,'(''move all beads upper bar='',t40,f10.5)') rmov2 
        write (6,'(''shift beads upper bar='',t40,f10.5)') rmov3 
        write (6,'(''move all + shift beads upper bar='',t40,f10.5)') rmov23         
        
    endif
   
    call bcast(irn)     !random number seed
    call bcast(hbar)    !hbar^2/2m
    call bcast(dt)      !dt the time step    
    call bcast(npart)   !particle number
    call bcast(nprot)   ! proton number
    call bcast(lpot)
    call bcast(ntab)    ! table points for potential
    call bcast(mmax)    ! max level of bisection, nbisect=2**mmax
    call bcast(nrepmax) ! max # of beads for reptation
    call bcast(isite)   !if true take initial
    call bcast(neq)     !equilibration blocks
    call bcast(nav)     !averaging blocks
    call bcast(nstep)   !steps per block
    call bcast(nstepdecor)   !calculate energy or x every this number steps to check correlation.
    call bcast(nchorizo) !number of chorizos, better be an even number
    call bcast(mov1step) ! move1 step
    call bcast(mov2step) ! move2 step move all beads at the same time
    call bcast(shftstep) ! shift all beads at the same time
    call bcast(x0step) !the random step move for 0th and nth chorizo
    call bcast(irep)    !type of action
    call bcast(resume)  ! resume mode or not
    call bcast(rangev)    ! potential range (fm)
    call bcast(rangerho)    ! one particle density rho range (fm)   
    call bcast(nrhobin)    ! the number of bins for one particle density rho.
    call bcast(icorrchk)    ! true=check correlation, false means not. 
    call bcast(iem)    ! true=EM force on.   
    call bcast(spinorbit)    ! true=spin-orbit force on.   
    call bcast(rbisect)
    call bcast(rmov1)
    call bcast(rmov2)
    call bcast(rmov3)
    call bcast(rmov23)
    call bcast(infile)
    call bcast(outfile)
   
    select case (irep)
        case (1:2,5:6)
            nest=21
            allocate(answer(nest)) 
            call setestnum(nest,irep)
            call addest(1,'reptation %')
            call addest(2,'bisection %')
        
            call addest(3,'numerator')
            call addest(4,'denominator|energy')
            call addest(5,'real f')	  
        
	        call addest(6,'sqrt(<rm^2>)') ! root-mean-square radii
	        call addest(7,'<rm>') ! this is the rm, mean-root-square, <sqrt(rm^2)>
            call addest(8,'rho (0th bin)')
      
            call addest(9,'bisection half %')
            call addest(10,'bisection L end %')
            call addest(11,'bisection R end %')
            call addest(12,'bisection extra %')  
        
	        call addest(13,'N sqrt(<rm^2>)') ! root-mean-square radii
	        call addest(14,'N <rm>') ! this is the rm, mean-root-square, <sqrt(rm^2)>
	        call addest(15,'P sqrt(<rm^2>)') ! root-mean-square radii
	        call addest(16,'P <rm>') ! this is the rm, mean-root-square, <sqrt(rm^2)>
        
            call addest(17,'KE') ! kinetic energy
            call addest(18,'PE') ! total PE
            call addest(19,'V6') ! v6' local potential
            call addest(20,'Vem') ! em potential
            call addest(21,'Vls') ! LS spin-orbit potential
        
            nestristra=3
            call setestnumristra(nestristra,nchorizo,irep)
        
        
            if ((irep==5).or.(irep==6)) then
                write(file_id,'(i5)') irep
                filename='potential_CalcMode' // trim(adjustl(file_id)) // '.txt'             
            else
                write(filename,'("potential.txt")')
            endif
            call addestristra(1,trim(filename))
	  
        case (4)	! VMC 
            nest=21 ! nest the same as pimc, match pimc
            allocate(answer(nest))
            call setestnum(nest,irep)
            call addest(1,'reptation %')
            call addest(2,'VMC %')
            call addest(3,'numerator')
            call addest(4,'denominator|energy')
            call addest(5,'real f')	  
	        call addest(6,'rm')
	        call addest(7,'<rm>')
	        call addest(8,'rho (0th bin)')	

            call addest(13,'N sqrt(<rm^2>)') ! root-mean-square radii
	        call addest(14,'N <rm>') ! this is the rm, mean-root-square, <sqrt(rm^2)>
	        call addest(15,'P sqrt(<rm^2>)') ! root-mean-square radii
	        call addest(16,'P <rm>') ! this is the rm, mean-root-square, <sqrt(rm^2)>
        
            call addest(17,'KE') ! kinetic energy
            call addest(18,'PE') ! total PE
            call addest(19,'V6') ! v6' local potential
            call addest(20,'Vem') ! em potential
            call addest(21,'Vls') ! LS spin-orbit potential

            nestristra=3
            nchorizo=1  ! just set as 1, there are 2 beads, bead 0 and 1.
            call setestnumristra(nestristra,nchorizo,irep) 
            write(filename,'("potential_VMC.txt")')
            call addestristra(1,trim(filename))
      
        case default
            write (6,'(''irep value has not been defined. '')') irep
            call abort 
    end select
   
! select irep case  

    select case (irep)
        
    case (2) ! PIMC
        
        nbisect=2**mmax   
        nchorizomid=nint(dble(nchorizo)/2)
        rhobinsize=rangerho/dble(nrhobin) 
    ! protection check begin:
        if (myrank().eq.0) then   
            if (nbisect.gt.nchorizo) then
	            write(6,*) '# of Bisection beads must be smaller than nchorizo! Stop!',nbisect,nchorizo
	            call abort
            endif
            if (nrepmax.gt.nchorizo) then  
	            write(6,*) '# of Reptation beads must be smaller than nchorizo! Stop!',nrepmax,nchorizo
	            call abort
            endif           
            
    ! protection check finished.            
              
            call cleanfolder(irep,resume)      
            if (.not.resume) then    
                open (12,FILE='answer.txt',FORM='FORMATTED',position='rewind')     
                call write1(12) 
            else
                open (12,FILE='answer.txt',FORM='FORMATTED',position='append')  
                write (12,*) 'RESUMED'  
            endif
            
            if (isite) then
                write (12,'(''input from sites'')')
            else
                if (.not.resume) then
                    write (12,'(''input from file'',t40,a20)') infile
                else
                    write (12,'(''input from file'',t40,a20)') outfile
                endif
            endif   
        endif
          
        !do i=0,611
        !    call ran2(dummy,irn)
        !    if (i.eq.611) then
        !    open(unit=9,form='formatted',file='buggyirn.txt')
        !       write(9,*) i,irn
        !    close(9)
        !    endif 
        !enddo        
        !call done
           
        call setrnall(irn)
        call v6stepcalcinit(npart,nprot,lpot,dt,hbar,ntab,rangev,mmax,nrepmax,iem,irep) 
        call wavefunctioninit(nchorizo,npart,nprot,ntab,hbar,rangev,rangerho,nrhobin,irep)  
        call stepinit(npart,nchorizo,hbar,dt,mmax,nrepmax,iem,spinorbit,irep,nstepdecor,x0step,mov1step,mov2step,shftstep,nrhobin &
            ,icorrchk,nprot,nav,nstep,neq,rbisect,rmov1,rmov2,rmov3,rmov23,nest,nestristra,resume)
        call brutinit(npart,nprot,lpot,hbar,iem,spinorbit,ntab,rangev) ! 3 in av18pot.f is v6' . this gives the initial L,R order for psi_T.
        call stepchorizoinit(isite,infile,outfile)   ! call brutinit first, also call stepinit first.       
        
        call zerestristra
        call zerest
        call zerrhodist   
           
        if (.not.resume) then
            it=0 
            ishift=1
        else
            ishift=neq+1 
            call readcheckpoint(irep,it) ! this is only for rank 0.
            call bcast(it)
        endif        
                
        call counterinit 
        call zerepstepmon
 
        if (myrank().eq.0) time00=mpi_wtime()
        do i=it+ishift,nav+neq  !do blocks
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
            if ((i.eq.2).and.isite) then 
                if (myrank().eq.0) then  
                    write (6,'(/,'' garbage isite config will be dropped '')')
		            write (6,'(/,''=================='')')
		            write (12,'(/,'' garbage isite config will be dropped '')')
		            write (12,'(/,''=================='')')
                endif
                call zerest
                call zerestristra
		        call zerrhodist 
                call zerepstepmon
                !it=0
            endif
            it=it+1 ! now it = i.
	        if ((i.eq.1).and.(myrank().eq.0)) time0=mpi_wtime()
	        call counterinit
            do j=1,nstep
                if ((i.eq.1).and.(j.eq.1).and.(myrank().eq.0)) time3=mpi_wtime()         
		        call stepstd(i) ! quantum Monte Carlo sampling step. 
                if ((i.eq.1).and.(j.eq.min(50,max(1,int(nstep*0.05_r8)))).and.(myrank().eq.0)) then
                    time4=mpi_wtime()
                    t1=nproc()*dble(nav+neq)*dble(nstep)*(time4-time3)/j/3600
                    call timedhms(dble(nstep)*(time4-time3)/j,day,hour,minute,second) 
            ! dble(nstep) to make protection, convert things to r8 at early stage to prevent i4 range explosion.
        write (6,'(/,''1 block time may be: '',i5,'' days'',i5,'' hours'',i5,'' minutes'',f10.3,'' seconds'')') &
                day,hour,minute,second
        write (12,'(/,''1 block time may be: '',i5,'' days'',i5,'' hours'',i5,'' minutes'',f10.3,'' seconds'')') &
                day,hour,minute,second
                    call timedhms(dble(nav+neq)*dble(nstep)*(time4-time3)/j,day,hour,minute,second)
        write (6,'(/,''Total time may be: '',i10,'' days'',i10,'' hours'',i10,'' minutes'',f10.3,'' seconds'')') &
                day,hour,minute,second
        write (12,'(/,''Total time may be: '',i10,'' days'',i10,'' hours'',i10,'' minutes'',f10.3,'' seconds'')') &
                day,hour,minute,second                
                    write (6,'(/,''Total core hour est: '',i10,'' hours'')') nint(t1)
                    write (12,'(/,''Total core hour est: '',i10,'' hours'')') nint(t1)
                endif             
            enddo
                   
            !if (myrank().eq.0) then
            !    write (6,'(/, ''iteration='',t20,i5,t30,''('',i5,'')'' )') it,i
            !    write (12,'(/, ''iteration='',t20,i5,t30,''('',i5,'')'' )') it,i
            !    call showstat 
            !endif
            !do k=0,nproc()-1 ! optional
            !    if (myrank().eq.k) then
            !        call printstat(6)                  
            !    endif
            !call barrier   
            !enddo
            
            if ((i.eq.1).and.(myrank().eq.0)) time6=mpi_wtime()
            call checkblkstat
            
            call update   !collect block averages
            call updateristra 
            call updaterhodistrbt(nchorizomid) ! choose the middle bead
            
            if ((i.eq.1).and.(myrank().eq.0)) time7=mpi_wtime()
            if (myrank().eq.0) then
                if (i.eq.1) time5=mpi_wtime()        
                write (6,'(/, ''iteration='',t20,i5,t30,''('',i5,'')'' )') it,i
                write (12,'(/, ''iteration='',t20,i5,t30,''('',i5,'')'' )') it,i
                call showstat                 
                answer=resstring()
                call writeanswer(answer)
                !   write (6,'(a120)') (answer(k),k=1,size(answer))
	            !write (12,'(a120)') (answer(k),k=1,size(answer))  
            endif 

            ! write out x.      
            call writeoutx(irep,it,i) ! write out x in stepstd. This is important. when call writex, write checkpoint is also called.
            !call writeoutxgs(irep,it,i)
            ! show rho distribution on the fly
            call writerhodist(nchorizomid)
            ! write out a file with the Vtot values at different chorizos 
            call writevtotbeads(dt)  


            if ((i.eq.1).and.(myrank().eq.0)) then
                time1=mpi_wtime()  
                t2=time1-time0
                call getstepstdcalctime(tcalc1)
                call getstepstdMPIpathtime(tmpiio1)
                write (6,'(/,''Time 1 block ='',f15.7,'' seconds'')') t2
                write (12,'(/,''Time 1 block ='',f15.7,'' seconds'')') t2  
                write (6,'(/,''Calculation time 1 block ='',g15.7,'' seconds '' '' ('',g15.7,'' % )'')') tcalc1,tcalc1/t2*100
                write (12,'(/,''Calculation time 1 block ='',g15.7,'' seconds '' '' ('',g15.7,'' % )'')') tcalc1,tcalc1/t2*100
                write (6,'(/,''MPI IO path time 1 block ='',g15.7,'' seconds '' '' ('',g15.7,'' % )'')') tmpiio1,tmpiio1/t2*100
                write (12,'(/,''MPI IO path time 1 block ='',g15.7,'' seconds '' '' ('',g15.7,'' % )'')') tmpiio1,tmpiio1/t2*100
                   
                !t3=time1-time5               
                !t6=time7-time6
                !t6p=t6/t2*100
                !t7=time9-time8
                !t8=time11-time10
                !t3a=t3-t7-t8
                !t3p=t3a/t2*100
                !t6a=t6+t7+t8                
                !write (6,'(/,''I/O time 1 block ='',g15.7,'' seconds '' '' ('',g15.7,'' % )''  )') t3a,t3p
                !write (12,'(/,''I/O time 1 block ='',g15.7,'' seconds '' '' ('',g15.7,'' % )''  )') t3a,t3p
                !write (6,'(/,''MPI communication time 1 block ='',g15.7,'' seconds '' '' ('',g15.7,'' % )''  )') t6a,t6p
                !write (12,'(/,''MPI communication time 1 block ='',g15.7,'' seconds '' '' ('',g15.7,'' % )''  )') t6a,t6p 
                
                call timedhms((nav+neq)*(time1-time0),day,hour,minute,second)
                write (6,'(/,''Total time estimation: '',i10,'' days'',i5,'' hours'',i5,'' minutes'',f10.3,'' seconds'')') &
                    day,hour,minute,second
                write (12,'(/,''Total time estimation: '',i10,'' days'',i5,'' hours'',i5,'' minutes'',f10.3,'' seconds'')') &
                    day,hour,minute,second  
                
                    t4=nproc()*dble(nav+neq)*(time1-time0)/3600
                write (6,'(/,''Total core hours estimation: '',i10,'' hours'')') nint(t4)
                write (12,'(/,''Total core hours estimation: '',i10,'' hours'')') nint(t4) 
            endif
            call barrier
        enddo
   
        call barrier
        if (myrank().eq.0) then
        time2=mpi_wtime()
        t5=nproc()*(time2-time00)/3600
        call timedhms((time2-time00),day,hour,minute,second)
        write (6,'(/,''Total time: '',i10,'' days'',i5,'' hours'',i5,'' minutes'',f10.3,'' seconds'')') &
                day,hour,minute,second
        write (12,'(/,''Total time: '',i10,'' days'',i5,'' hours'',i5,'' minutes'',f10.3,'' seconds'')') &
                day,hour,minute,second   
        write (6,'(/,''Total core hours: '',f15.7,'' hours'')') t5
        write (12,'(/,''Total core hours: '',f15.7,'' hours'')') t5
        close (12)
        endif        
        
    case (4) ! VMC 
	
        dt=0
        mmax=max(1,mmax)
        nchorizomid=nint(dble(nchorizo)/2)
        rhobinsize=rangerho/dble(nrhobin) 	
        isite=.true. ! do not have to be true but set it as true for simplicity.
	    if (nprot<=1) iem=.false. ! so no pp em force.
        nrepmax=0 ! no reptation
        
        
        call setrnall(irn) ! different thread have different seed. 
    
        if (myrank().eq.0) call cleanfolder(irep,resume)
    
        call v6stepcalcinit(npart,nprot,lpot,dt,hbar,ntab,rangev,mmax,nrepmax,iem,irep) 
        call wavefunctioninit(nchorizo,npart,nprot,ntab,hbar,rangev,rangerho,nrhobin,irep)  
        call stepinit(npart,nchorizo,hbar,dt,mmax,nrepmax,iem,spinorbit,irep,nstepdecor,x0step,mov1step,mov2step,shftstep &
            ,nrhobin,icorrchk,nprot,nav,nstep,neq,rbisect,rmov1,rmov2,rmov3,rmov23,nest,nestristra,resume)
        call brutinit(npart,nprot,lpot,hbar,iem,spinorbit,ntab,rangev) ! 3 in av18pot.f is v6' . this gives the initial L,R order for psi_T.
        call stepchorizoinit(isite,infile,outfile)   ! call brutinit first.	call stepinit first.   
   
        call zerest
        call zerestristra
        call zerrhodist
        it=0

        if (myrank().eq.0) then
            open(12,FILE='answer_VMC.txt',FORM='FORMATTED') 
            call write1(12)
        endif

        if (isite) then
            if (myrank().eq.0) write (12,'(''input from sites'')')
        else
            if (myrank().eq.0) write (12,'(''input from file'',t40,a20)') infile
        endif
   
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
		        call stepstd(j)
	        enddo
	        if (myrank().eq.0) call showstat 
            call update   !collect block averages
            call updateristra
        
      
            if (myrank().eq.0) then
                answer=resstring()
                write (6,'(/,''iteration = '',t30,i14)') it 
                write (12,'(/,''iteration = '',t30,i14)') it
                call writeanswer(answer)
                !   write (6,'(a120)') (answer(k),k=1,size(answer))
	            !write (12,'(a120)') (answer(k),k=1,size(answer))  
            endif
      
            ! write out x.      
            call writeoutxgs(irep,it,i) 
            ! show rho distribution on the fly
	        call writerhodist(nchorizomid) 
            ! write out a file with the Vtot values at different chorizos
            call writevtotbeads(dt)   
     
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
        call CPU_TIME(time2)
        write (6,'(/,''Total time ='',f10.3,'' hours'')') (time2-time00)/3600.
        write (12,'(/,''Total time ='',f10.3,'' hours'')') (time2-time00)/3600.
        close (12)
        close(11)
        endif	 
        
    case (5:6) ! Calc Mode  
        
        nbisect=2**mmax   
        nchorizomid=nint(dble(nchorizo)/2.0_r8)
        rhobinsize=rangerho/dble(nrhobin)       

        if (myrank().eq.0) then 
            time00=mpi_wtime()
            write(file_id,'(i5)') irep
            filename='answer_CalcMode' // trim(adjustl(file_id)) // '.txt'        
            open (12,FILE=trim(filename),FORM='FORMATTED')  
            call write1(12) 
        endif

        if (myrank().eq.0) call cleanfolder(irep,resume)
    
        call v6stepcalcinit(npart,nprot,lpot,dt,hbar,ntab,rangev,mmax,nrepmax,iem,irep) 
        call wavefunctioninit(nchorizo,npart,nprot,ntab,hbar,rangev,rangerho,nrhobin,irep)  
        call stepinit(npart,nchorizo,hbar,dt,mmax,nrepmax,iem,spinorbit,irep,nstepdecor,x0step,mov1step,mov2step,shftstep &
            ,nrhobin,icorrchk,nprot,nav,nstep,neq,rbisect,rmov1,rmov2,rmov3,rmov23,nest,nestristra,resume)
        call brutinit(npart,nprot,lpot,hbar,iem,spinorbit,ntab,rangev) ! 3 in av18pot.f is v6' . this gives the initial L,R order for psi_T.
        !call stepchorizoinit(isite,infile,outfile)   ! call brutinit first.	
    

        call zerest
        call zerestristra
        call zerrhodist  

        call counterinit 
        call zerepstepmon
        
        
        call computeinit(answer,rhobinsize)
        call responseinfo(nrespest,ne0respbeads)

        call compute  ! readcheckpoint is in compute.    
       
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
    case default
        write (6,'(''irep value has not been defined. So nothing has been done '')') irep        
    end select
      
! finalize  
    if (myrank().eq.0) write(6,*) 'The program end normally!'
    call barrier 
    call done
! end
    
contains
    
    subroutine write1(myunit)
    integer(kind=i4) :: myunit
    write (myunit,'(''Resume Mode ='',t30,l10)') resume 
    write (myunit,'(''MPI # of cores ='',t30,i20)') nproc()
    write (myunit,'(''random number seed ='',t30,i20)') irn
    write (myunit,'(''hbar^2/2m ='',t40,f10.5)') hbar
    write (myunit,'(''time step ='',t40,g15.7)') dt
    write (myunit,'(''npart ='',t40,i10)') npart
    write (myunit,'(''nproton ='',t40,i10)') nprot  
    write (myunit,'(''lpot ='',t40,i10)') lpot
    write (myunit,'(''table size ='',t40,i10)') ntab   
    write (myunit,'(''bisection order ='',t40,i10)') mmax
    write (myunit,'(''# of reptation beads ='',t40,i10)') nrepmax  
    write (myunit,'(''start from sites ='',t40,l10)') isite
    write (myunit,'(''equilibration blocks ='',t40,i10)') neq
    write (myunit,'(''averaging blocks ='',t40,i10)') nav
    write (myunit,'(''nsteps/blocks ='',t40,i10)') nstep
    write (myunit,'(''stepdecor ='',t40,i10)') nstepdecor
    write (myunit,'(''nchorizo ='',t40,i10)') nchorizo
    write (myunit,'(''stepmove1 ='',t40,f10.5)') mov1step
    write (myunit,'(''stepmove2 ='',t40,f10.5)') mov2step
    write (myunit,'(''stepmoveshift ='',t40,f10.5)') shftstep
    write (myunit,'(''x0step ='',t40,f10.5)') x0step
    write (myunit,'(''potential range ='',t40,f10.5)') rangev
    write (myunit,'(''one particle density rho range ='',t40,f10.5)') rangerho
    write (myunit,'(''bin # of rho ='',t40,i10)') nrhobin
    write (myunit,'(''correlation check ='',t40,l10)') icorrchk
    write (myunit,'(''EM force ='',t40,l10)') iem
    write (myunit,'(''Spin-Orbit force ='',t40,l10)') spinorbit
    write (myunit,'(''bisection (reptation) upper bar='',t40,f10.5)') rbisect 
    write (myunit,'(''move 1 by 1 upper bar='',t40,f10.5)') rmov1 
    write (myunit,'(''move all beads upper bar='',t40,f10.5)') rmov2 
    write (myunit,'(''shift beads upper bar='',t40,f10.5)') rmov3 
    write (myunit,'(''move all + shift beads upper bar='',t40,f10.5)') rmov23     
    return
    end subroutine write1 
    
    subroutine setrnall(irn)  
    integer(kind=i8) :: irn,irnsave
    integer(kind=i8), allocatable :: irnall(:) 
    integer(kind=i4) :: i
    real(kind=r8) :: dummy
    do i=0,nproc()-1
        if (myrank() /= 0) call ran2(dummy,irn)
        if (myrank().eq.i) then
            irnsave=irn
            call setrn(irn)    ! different cores have different seed. 
            exit       
        endif 
    enddo
    allocate(irnall(0:nproc()-1))
    call gather(irnsave,irnall)
    if (myrank().eq.0) then
        open(unit=9,form='formatted',file='irnlist.txt')
        do i=0,nproc()-1
            write(9,*) i,irnall(i)
        enddo
        close(9)
    endif
    deallocate(irnall)    
    return
    end subroutine setrnall

    subroutine setreadin 
! for VTune Amplifier only.
    irn=88888888     !random number seed
    hbar=20.735    !hbar^2/2m
    dt=0.00075      !dt the time step    
    npart=4   !particle number
    nprot=2   ! proton number
    lpot=11 ! potential number
    ntab=50000   ! table points for potential
    mmax=6    ! max level of bisection, nbisect=2**mmax
    nrepmax=0 ! max # of beads for reptation, 0 means off.
    isite=.true.   !if true take initial
    neq=1     !equilibration blocks
    nav=5000     !averaging blocks
    nstep=2000   !steps per block
    nstepdecor=142   !calculate energy every this number steps
    nchorizo=200 !number of chorizos, better be an even number
    mov1step=0.08 ! move1 step. move beads 1 by1 step
    mov2step=0.01 ! move2 step move all beads at the same time. For irep=4 the VMC, it is VMC move step delta.
    shftstep=0.1 ! shift all beads at the same time
    x0step=1.0 !the random step move for 0th and nth chorizo
    irep=2    !type of action. 2=pimc. 4=vmc. 5= calculation mode (consume less ram). 6=perhaps big memory calculation mode.
    rangev=20.0    ! potential range (fm)
    rangerho=5.0    ! one particle density rho range (fm)   
    nrhobin=50    ! the number of bins for one particle density rho.
    icorrchk=.true.    ! true=check correlation, false means not. 
    iem=.true.    ! true=EM force on.  
    rbisect=0.9 ! 0.9 bisection upper bar (reptation if enabled is in it.)
    rmov1=0.95 ! move 1 by 1 upper bar
    rmov2=0.967 ! move all beads  upper bar 
    rmov3=0.984 ! shift beads rate upper bar
    rmov23=1.0 ! move all then shift beads upper bar
    infile='he4ristra.in'
    outfile='he4ristra.out'
    return
    end subroutine setreadin
    
   end program v6pimc
	

    

    
    
    
	
