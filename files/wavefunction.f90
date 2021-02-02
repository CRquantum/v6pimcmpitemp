module wavefunction
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), private, parameter :: pi=4.0_r8*atan(1.0_r8) 
   integer(kind=i4), private,save :: npart,ntab,irep,nrhobin,nchorizo
   integer(kind=i4), private, save :: nprot,nneut,nspin,nisospin,nbasis
   real(kind=r8), private, save, allocatable :: utab(:),dutab(:),d2utab(:)
   real(kind=r8), private, save, allocatable :: vtab(:)
   real(kind=r8), private, save :: el,scale,range,tail,hbar,scalep,rangeu,rangev,rangevsafe,rangerho
   real(kind=r8), private, save :: azabo,bzabo,dzabo
   real(kind=r8), private, save :: rhobinsize
   real(kind=r8), private, save, allocatable :: rhodist(:,:),rhodisterr(:,:),rhodisttot(:,:),rhodist2tot(:,:) &
                                               ,wtrhodistblk(:),wtrhodisttot(:),rhodistblk(:,:),rhodist2blk(:,:) &
                                               ,icrhodist(:)  &
                                               ,rhondist(:,:),rhondisterr(:,:),rhondisttot(:,:),rhondist2tot(:,:) &
                                               ,wtrhondistblk(:),wtrhondisttot(:),rhondistblk(:,:),rhondist2blk(:,:) &
                                               ,rhopdist(:,:),rhopdisterr(:,:),rhopdisttot(:,:),rhopdist2tot(:,:) &
                                               ,wtrhopdistblk(:),wtrhopdisttot(:),rhopdistblk(:,:),rhopdist2blk(:,:)
    character(len=120), private, save :: file_rhodist,file_id   
                                               
!
! npart = total number of particles
! ntab = number of jastrow and potential table values
! utab, dutab, d2utab = jastrow tables
! vtab = potential tables
! el = side of simulation cell
! scale = factor to convert distance to a table index
! range = jastrow range
!
contains
    subroutine wavefunctioninit(nchorizoin,npartin,nprotin,ntabin &
                ,hbarin,rangevin,rangerhoin,nrhobinin,irepin)
    use math
    integer(kind=i4) :: npartin,ntabin,nrhobinin,irepin,nchorizoin,nprotin
    real(kind=r8) :: hbarin,rangevin,rangerhoin
    
    nchorizo=nchorizoin

    hbar=hbarin
    npart=npartin
    ntab=ntabin
    irep=irepin
 
    rangev=rangevin
    rangerho=rangerhoin ! I suggest 10 fm.
    nrhobin=nrhobinin ! I suggest 100 bins.
    rhobinsize=rangerho/nrhobin
    allocate(rhodist(0:nrhobin,0:nchorizo),rhodisterr(0:nrhobin,0:nchorizo) &
    ,rhodisttot(0:nrhobin,0:nchorizo),rhodist2tot(0:nrhobin,0:nchorizo) &
    ,rhodistblk(0:nrhobin,0:nchorizo),rhodist2blk(0:nrhobin,0:nchorizo))
    allocate(icrhodist(0:nchorizo),wtrhodistblk(0:nchorizo),wtrhodisttot(0:nchorizo))
    
    allocate(rhondist(0:nrhobin,0:nchorizo),rhondisterr(0:nrhobin,0:nchorizo) &
    ,rhondisttot(0:nrhobin,0:nchorizo),rhondist2tot(0:nrhobin,0:nchorizo) &
    ,rhondistblk(0:nrhobin,0:nchorizo),rhondist2blk(0:nrhobin,0:nchorizo))
    allocate(wtrhondistblk(0:nchorizo),wtrhondisttot(0:nchorizo))
    
    allocate(rhopdist(0:nrhobin,0:nchorizo),rhopdisterr(0:nrhobin,0:nchorizo) &
    ,rhopdisttot(0:nrhobin,0:nchorizo),rhopdist2tot(0:nrhobin,0:nchorizo) &
    ,rhopdistblk(0:nrhobin,0:nchorizo),rhopdist2blk(0:nrhobin,0:nchorizo))
    allocate(wtrhopdistblk(0:nchorizo),wtrhopdisttot(0:nchorizo))
   
    nprot=nprotin
    nneut=npart-nprot
    nspin=2**npart ! n=4, it is 16.
    call combin(npart,nprot,nisospin) ! give value for nisospin, n=4, it is 4C2= 6.
    nbasis=nspin*nisospin ! n=4, it is 96.      

    return
    end subroutine wavefunctioninit
                
    subroutine wavefunctionfilesinit(irepin)
    integer(kind=i4) :: irepin
    
    select case (irepin)
    case (4)   
        write(file_rhodist,'("rhodistribution_VMC.txt")')  
    case (5,6)    
        write(file_id,'(i5)') irepin
        file_rhodist='rhodistribution_CalcMode' // trim(adjustl(file_id)) // '.txt'         
    case default
        write(file_rhodist,'("rhodistribution.txt")')  
    end select      
    return
    end subroutine wavefunctionfilesinit

    subroutine zerrhodist ! zero density distribution
    icrhodist=0
    rhodisttot=0
    rhodist2tot=0
    rhodist=0
    rhodisterr=0
    wtrhodistblk=0
    wtrhodisttot=0
    rhodistblk=0
    rhodist2blk=0
    
    rhondisttot=0
    rhondist2tot=0
    rhondist=0
    rhondisterr=0
    wtrhondistblk=0
    wtrhondisttot=0
    rhondistblk=0
    rhondist2blk=0 
    
    rhopdisttot=0
    rhopdist2tot=0
    rhopdist=0
    rhopdisterr=0
    wtrhopdistblk=0
    wtrhopdisttot=0
    rhopdistblk=0
    rhopdist2blk=0
    
    return
    end subroutine zerrhodist
	  
    subroutine funcrmsq(x,rm,rmsq) 
    real(kind=r8), dimension(3,npart) :: x
    real(kind=r8) :: rm,rmsq,rm1d(npart),rmsq1d(npart)
    integer(kind=i4) :: l,i
    real(kind=r8), dimension(3) :: rc 
    do l=1,3   
	    rc(l)=sum(x(l,:))/npart
    enddo

    do i=1,npart
        rmsq1d(i)=sum((x(:,i)-rc(:))**2)
        rm1d(i)=sqrt(rmsq1d(i))
    enddo 
    rmsq=sum(rmsq1d)/npart
    rm=sum(rm1d)/npart

 !   i=1 ! select particle
	!rmsq=sum((x(:,i)-rc(:))**2)
	!rm=sqrt(rmsq)

    return
    end subroutine funcrmsq

    subroutine samplerhodistrbt(x,i,cwtl2rm,cwtr2lm,rf,rhodistnowout,rm,rmsq,rmn,rmnsq,rmp,rmpsq) 
!,rhodistout,rhodisterrout) ! i is the picked the bead, combine rm,rmsq sampling.
    use mympi
    use v6stepcalc
    real(kind=r8), dimension(3,npart) :: x
    real(kind=r8) :: r,rsq,xnk,xkh1,xkh2,basis_k,signrf,rf,rfn,rfp
    integer(kind=i4) :: i,j,l,iknum,ip1,im1
    real(kind=r8), dimension(3) :: rc 
    real(kind=r8), dimension(0:nrhobin) :: rhodistnow,rhodistnowout,rhondistnow,rhopdistnow
    complex(kind=r8) :: cwtl2rm(0:nspin-1,nisospin),cwtr2lm(0:nspin-1,nisospin) &
                       ,cwtr2lnnew(0:nspin-1,nisospin,npart),cwtr2lpnew(0:nspin-1,nisospin,npart)
    real(kind=r8) :: rm,rmsq,rm1d(npart),rmsq1d(npart) &
                    ,rmn,rmnsq,rmn1d(npart),rmnsq1d(npart),rmp,rmpsq,rmp1d(npart),rmpsq1d(npart)
    
    
    signrf=sign(1.0_r8,rf)
    
    rhodistnow(:)=0
    rhodistnowout(:)=0
    
    rm=0
    rmsq=0
    rmn=0
    rmnsq=0
    rmp=0
    rmpsq=0
    
    do l=1,3   
	    rc(l)=sum(x(l,:))/npart
    enddo  

    call projectornp(cwtr2lm,cwtr2lnnew,cwtr2lpnew)    
    
    do j=1,npart ! select particle
        !j=1   
        
        rfn=real(sum(conjg(cwtl2rm(:,:))*cwtr2lnnew(:,:,j)))/abs(rf)
        rfp=real(sum(conjg(cwtl2rm(:,:))*cwtr2lpnew(:,:,j)))/abs(rf)
        
        rsq=sum((x(:,j)-rc(:))**2)
         
      ! sample radii, rms values.  
        rmsq1d(j)=rsq
        rm1d(j)=sqrt(rmsq1d(j)) 
        
        rmnsq1d(j)=rmsq1d(j)*rfn
        rmn1d(j)=rm1d(j)*rfn    
        
        rmpsq1d(j)=rmsq1d(j)*rfp
        rmp1d(j)=rm1d(j)*rfp         
        
        r=sqrt(rsq)
        iknum=int(r/rhobinsize)    
        xnk=4.0_r8/3*pi*rhobinsize**3*((iknum+1)**3-iknum**3)   
        if (iknum.le.nrhobin) then  	  
            !xkh1=iknum*rhobinsize
            !xkh2=(iknum+1)*rhobinsize
            !IF(  (r.GE.xkh1).and.(r.LE.xkh2)  ) THEN
            !        basis_k=signrf ! bc the sturcture of our PIMC, we need sign of rf instead of just 1.
            !ELSE
            !        basis_k=0
            !END IF	
            basis_k=signrf ! bc the sturcture of our PIMC, we need sign of rf instead of just 1.
            rhodistnow(iknum)=basis_k/xnk 
            rhodistnowout(iknum)=rhodistnowout(iknum)+rhodistnow(iknum)
            rhodistblk(iknum,i)=rhodistblk(iknum,i)+rhodistnow(iknum)  
            rhodist2blk(iknum,i)=rhodist2blk(iknum,i)+rhodistnow(iknum)**2 
            
            rhondistnow(iknum)=rfn/xnk
            rhopdistnow(iknum)=rfp/xnk
                               
            rhondistblk(iknum,i)=rhondistblk(iknum,i)+rhondistnow(iknum)  
            rhondist2blk(iknum,i)=rhondist2blk(iknum,i)+rhondistnow(iknum)**2             
            
            rhopdistblk(iknum,i)=rhopdistblk(iknum,i)+rhopdistnow(iknum)  
            rhopdist2blk(iknum,i)=rhopdist2blk(iknum,i)+rhopdistnow(iknum)**2             
            
            
        endif 
        wtrhodistblk(i)=wtrhodistblk(i)+1
     
    enddo
    
    if (npart/=0) then
        rhodistnowout(:)=rhodistnowout(:)/npart
        rmsq=signrf*sum(rmsq1d)/npart
        rm=signrf*sum(rm1d)/npart
        if (nneut/=0) then
            rmnsq=sum(rmnsq1d)/nneut
            rmn=sum(rmn1d)/nneut
        endif
        if (nprot/=0) then
            rmpsq=sum(rmpsq1d)/nprot
            rmp=sum(rmp1d)/nprot 
        endif
    endif
 
    return
    end subroutine samplerhodistrbt    
    
    
    
    subroutine samplerhodistrbt0(x,i,cwtl2rm,cwtr2lm,rf,rhodistnowout) !,rhodistout,rhodisterrout) ! i is the picked the bead
    use mympi
    use v6stepcalc
    real(kind=r8), dimension(3,npart) :: x
    real(kind=r8) :: r,xnk,xkh1,xkh2,basis_k,signrf,rf
    integer(kind=i4) :: i,j,l,iknum,ip1,im1
    real(kind=r8), dimension(3) :: rc 
    real(kind=r8), dimension(0:nrhobin) :: rhodistnow,rhodistnowout,rhondistnow,rhopdistnow
    complex(kind=r8) :: cwtl2rm(0:nspin-1,nisospin),cwtr2lm(0:nspin-1,nisospin) &
                       ,cwtr2lnnew(0:nspin-1,nisospin,npart),cwtr2lpnew(0:nspin-1,nisospin,npart)
    
    signrf=sign(1.0_r8,rf)
    
    rhodistnow(:)=0
    rhodistnowout(:)=0
    do l=1,3   
	    rc(l)=sum(x(l,:))/npart
    enddo  

    call projectornp(cwtr2lm,cwtr2lnnew,cwtr2lpnew)
 
    do j=1,npart ! select particle
        !j=1
        r=sqrt(sum((x(:,j)-rc(:))**2))
        iknum=int(r/rhobinsize)    
        xnk=4.0_r8/3*pi*rhobinsize**3*((iknum+1)**3-iknum**3)   
        if (iknum.le.nrhobin) then  	  
            !xkh1=iknum*rhobinsize
            !xkh2=(iknum+1)*rhobinsize
            !IF(  (r.GE.xkh1).and.(r.LE.xkh2)  ) THEN
            !        basis_k=signrf ! bc the sturcture of our PIMC, we need sign of rf instead of just 1.
            !ELSE
            !        basis_k=0
            !END IF	
            basis_k=signrf ! bc the sturcture of our PIMC, we need sign of rf instead of just 1.
            rhodistnow(iknum)=basis_k/xnk 
            rhodistnowout(iknum)=rhodistnowout(iknum)+rhodistnow(iknum)
            rhodistblk(iknum,i)=rhodistblk(iknum,i)+rhodistnow(iknum)  
            rhodist2blk(iknum,i)=rhodist2blk(iknum,i)+rhodistnow(iknum)**2 
            
      
            rhondistnow(iknum)=real(sum(conjg(cwtl2rm(:,:))*cwtr2lnnew(:,:,j)))/abs(rf)/xnk
            rhopdistnow(iknum)=real(sum(conjg(cwtl2rm(:,:))*cwtr2lpnew(:,:,j)))/abs(rf)/xnk
            
                    
            rhondistblk(iknum,i)=rhondistblk(iknum,i)+rhondistnow(iknum)  
            rhondist2blk(iknum,i)=rhondist2blk(iknum,i)+rhondistnow(iknum)**2             
            
            rhopdistblk(iknum,i)=rhopdistblk(iknum,i)+rhopdistnow(iknum)  
            rhopdist2blk(iknum,i)=rhopdist2blk(iknum,i)+rhopdistnow(iknum)**2             
            
            
        endif 
        wtrhodistblk(i)=wtrhodistblk(i)+1
    enddo
    rhodistnowout(:)=rhodistnowout(:)/npart
    return
    end subroutine samplerhodistrbt0
    
   
    subroutine updaterhodistrbt(i) 
    use mympi
    integer(kind=i4) :: i
    real(kind=r8), dimension(0:nchorizo) :: wtrhodistsum
    real(kind=r8), dimension(0:nrhobin) :: av,av2,avn,avn2,avp,avp2
    real(kind=r8), dimension(0:nrhobin,0:nchorizo) :: rhodistsum,rhodist2sum &
                                                     ,rhondistsum,rhondist2sum,rhopdistsum,rhopdist2sum
   
    call addall(wtrhodistblk(i),wtrhodistsum(i))
    
    call addall(rhodistblk(:,i),rhodistsum(:,i))
    call addall(rhodist2blk(:,i),rhodist2sum(:,i))
    
    call addall(rhondistblk(:,i),rhondistsum(:,i))
    call addall(rhondist2blk(:,i),rhondist2sum(:,i))  
    
    call addall(rhopdistblk(:,i),rhopdistsum(:,i))
    call addall(rhopdist2blk(:,i),rhopdist2sum(:,i))
    
    if (myrank().eq.0) then  
        
        wtrhodistblk(i)=wtrhodistsum(i) 
        wtrhodisttot(i)=wtrhodisttot(i)+wtrhodistblk(i)
        
        rhodistblk(:,i)=rhodistsum(:,i) 
        rhodist2blk(:,i)=rhodist2sum(:,i)                
        rhodisttot(:,i)=rhodisttot(:,i)+rhodistblk(:,i)
        rhodist2tot(:,i)=rhodist2tot(:,i)+rhodist2blk(:,i)
        
        rhondistblk(:,i)=rhondistsum(:,i) 
        rhondist2blk(:,i)=rhondist2sum(:,i)                
        rhondisttot(:,i)=rhondisttot(:,i)+rhondistblk(:,i)
        rhondist2tot(:,i)=rhondist2tot(:,i)+rhondist2blk(:,i)        
        
        rhopdistblk(:,i)=rhopdistsum(:,i) 
        rhopdist2blk(:,i)=rhopdist2sum(:,i)                
        rhopdisttot(:,i)=rhopdisttot(:,i)+rhopdistblk(:,i)
        rhopdist2tot(:,i)=rhopdist2tot(:,i)+rhopdist2blk(:,i)        
        
        
        
        if (wtrhodisttot(i)/=0) then
            
            rhodist(:,i)=rhodisttot(:,i)/wtrhodisttot(i)
            
            rhondist(:,i)=rhondisttot(:,i)/wtrhodisttot(i)
            
            rhopdist(:,i)=rhopdisttot(:,i)/wtrhodisttot(i)
            
            
            
        else
            
            rhodist(:,i)=0
            rhondist(:,i)=0
            rhopdist(:,i)=0
            
            
        endif
        
        av(:)=rhodist(:,i)
        
        avn(:)=rhondist(:,i)
        avp(:)=rhopdist(:,i)
        
        if (wtrhodisttot(i) /= 0 ) then
            av2(:)=rhodist2tot(:,i)/wtrhodisttot(i) 
            
            avn2(:)=rhondist2tot(:,i)/wtrhodisttot(i) 
            avp2(:)=rhopdist2tot(:,i)/wtrhodisttot(i) 
            
        else
            av2(:)=0
            avn2(:)=0
            avp2(:)=0
        endif   
        rhodisterr(:,i)=sqrt(abs(av2(:)-av(:)**2)/max(1.0_r8,wtrhodisttot(i)-1)) 
        rhondisterr(:,i)=sqrt(abs(avn2(:)-avn(:)**2)/max(1.0_r8,wtrhodisttot(i)-1)) 
        rhopdisterr(:,i)=sqrt(abs(avp2(:)-avp(:)**2)/max(1.0_r8,wtrhodisttot(i)-1)) 
        
    endif
    wtrhodistblk(i)=0 
    rhodistblk(:,i)=0
    rhodist2blk(:,i)=0
    
    rhondistblk(:,i)=0
    rhondist2blk(:,i)=0
    
    rhopdistblk(:,i)=0
    rhopdist2blk(:,i)=0    


    return
    end subroutine updaterhodistrbt
   
    
    subroutine writerhodist(n)
    use mympi
    use estimator
    integer(kind=i4) :: k,n,myunit
    real(kind=r8) :: val1,val2,err1,err2,val,error &
                    ,val1n,err1n,err2n,valn,errorn  &
                    ,val1p,err1p,err2p,valp,errorp  &
                    ,r,factrad
    if (myrank().eq.0) then
        call resultest(4,val2,err2) ! #4 is the denominator  
	    open(newunit=myunit,FILE=trim(file_rhodist),FORM='FORMATTED')	 
	    do k=0,nrhobin 
            
            val1=rhodist(k,n)
            err1=rhodisterr(k,n)
            val=val1/val2
            error=sqrt((1/val2*err1)**2+(val1/(val2**2)*err2)**2)
            
            val1n=rhondist(k,n)
            err1n=rhondisterr(k,n)
            valn=val1n/val2
            errorn=sqrt((1/val2*err1n)**2+(val1n/(val2**2)*err2)**2)
            
            val1p=rhopdist(k,n)
            err1p=rhopdisterr(k,n)
            valp=val1p/val2
            errorp=sqrt((1/val2*err1p)**2+(val1p/(val2**2)*err2)**2)            
                       
            r=(k+0.5)*rhobinsize
            factrad=4*pi*r**2 ! for radial density distribution, g(r)=4*pi*r^2*rho(r), int g(r) dr=1.
	        write (myunit,'(i10,1x,15(f10.5,1x))') k,r,val,error,valn,errorn,valp,errorp &
                                                  ,factrad*val,factrad*error & 
                                                  ,factrad*valn,factrad*errorn,factrad*valp,factrad*errorp
	    enddo
	    close(myunit) 
    endif 
    return
    end subroutine writerhodist    
    
    subroutine checkpointrhodist(input,wtrhodisttotio,rhodisttotio,rhodist2totio &
                                ,rhondisttotio,rhondist2totio,rhopdisttotio,rhopdist2totio)
    real(kind=r8) :: wtrhodisttotio(0:nchorizo),rhodisttotio(0:nrhobin,0:nchorizo),rhodist2totio(0:nrhobin,0:nchorizo) &
                    ,rhondisttotio(0:nrhobin,0:nchorizo),rhondist2totio(0:nrhobin,0:nchorizo) &
                    ,rhopdisttotio(0:nrhobin,0:nchorizo),rhopdist2totio(0:nrhobin,0:nchorizo)
    logical :: input 
    if (input) then   
	    wtrhodisttot(:)=wtrhodisttotio(:)
        
        rhodisttot(:,:)=rhodisttotio(:,:)
        rhodist2tot(:,:)=rhodist2totio(:,:)
        
        rhondisttot(:,:)=rhondisttotio(:,:)
        rhondist2tot(:,:)=rhondist2totio(:,:) 
        
        rhopdisttot(:,:)=rhopdisttotio(:,:)
        rhopdist2tot(:,:)=rhopdist2totio(:,:)        
        
    else
	    wtrhodisttotio(:)=wtrhodisttot(:)
        
        rhodisttotio(:,:)=rhodisttot(:,:)
        rhodist2totio(:,:)=rhodist2tot(:,:) 
        
        rhondisttotio(:,:)=rhondisttot(:,:)
        rhondist2totio(:,:)=rhondist2tot(:,:)         
        
        rhopdisttotio(:,:)=rhopdisttot(:,:)
        rhopdist2totio(:,:)=rhopdist2tot(:,:)         
        
        
    endif
    return
    end subroutine checkpointrhodist     
    
    
    

    !subroutine psitcwt(x,cwt) ! input x, give the 96 amplituds. for now, make it 0,+-1.
    !use v6stepcalc
    !use brut
    !real(kind=r8), dimension(3,npart) :: x
    !complex(kind=r8) :: cwt(0:nspin-1,nisospin) ! need to input nspin, nisospin first. 
    !integer(kind=i4) :: i,j,k,l
    !integer(kind=i4) :: invspintmp(nbasis),invispintmp(nbasis)
    !
    !
    !
    !return
    !end subroutine psitcwt  
   

   

end module wavefunction
