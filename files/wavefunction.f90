module wavefunction
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private,save :: npart,ntab,irep,nrhobin,nchorizo
   integer(kind=i4), private, save :: nprot,nspin,nisospin,nbasis
   real(kind=r8), private, save, allocatable :: utab(:),dutab(:),d2utab(:)
   real(kind=r8), private, save, allocatable :: vtab(:)
   real(kind=r8), private, save :: el,scale,range,tail,hbar,scalep,rangeu,rangev,rangevsafe,rangerho
   real(kind=r8), private, save :: azabo,bzabo,dzabo
   real(kind=r8), private, save :: rhobinsize
   real(kind=r8), private, save, allocatable :: rhodist(:,:),rhodisterr(:,:),rhodisttot(:,:),rhodist2tot(:,:) &
                                               ,wtrhodistblk(:),wtrhodisttot(:),rhodistblk(:,:),rhodist2blk(:,:)
    integer(kind=i4), private, save, allocatable :: icrhodist(:)  
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
    ! need to add a,b,d 
    use math
    integer(kind=i4) :: npartin,ntabin,nrhobinin,irepin,nchorizoin,nprotin
    real(kind=r8) :: elin,rangein,hbarin,ain,bin,din,rangevin,rangeuin,rangerhoin
    ! external jastrow,vpot
    integer(kind=i4) :: i,j,k,ic
    real(kind=r8) dr,r
    real(kind=r8) vintp,uintp,vexct,verr
    integer(kind=i4) :: ntabmore
    nchorizo=nchorizoin

    hbar=hbarin
    npart=npartin
    ntab=ntabin
    irep=irepin
 
    rangev=rangevin
    rangerho=rangerhoin
    nrhobin=nrhobinin
    rhobinsize=rangerho/nrhobin
    allocate(rhodist(0:nrhobin,0:nchorizo),rhodisterr(0:nrhobin,0:nchorizo) &
    ,rhodisttot(0:nrhobin,0:nchorizo),rhodist2tot(0:nrhobin,0:nchorizo) &
    ,rhodistblk(0:nrhobin,0:nchorizo),rhodist2blk(0:nrhobin,0:nchorizo))
    allocate(icrhodist(0:nchorizo),wtrhodistblk(0:nchorizo),wtrhodisttot(0:nchorizo))
   
    nprot=nprotin
    nspin=2**npart ! n=4, it is 16.
    call combin(npart,nprot,nisospin) ! give value for nisospin, n=4, it is 4C2= 6.
    nbasis=nspin*nisospin ! n=4, it is 96.      
   
    return
    end subroutine wavefunctioninit

    subroutine zerrhodist
    ! zero density distribution
    icrhodist=0
    rhodisttot=0
    rhodist2tot=0
    rhodist=0
    wtrhodistblk=0
    wtrhodisttot=0
    rhodistblk=0
    rhodist2blk=0
    return
    end subroutine zerrhodist
	  
    subroutine funcrmsq(x,rm,rmsq) 
    real(kind=r8), dimension(3,npart) :: x,f
    real(kind=r8) :: u,df,d2psi,v,val,rm,rmsq
    integer(kind=i4) :: i,j,k,l,dim,ib,jb,is,js,jbmax,index
    real(kind=r8), dimension(3) :: dx,sumcor,rc
    real(kind=r8) :: r,dr,duu,c1,c2,c3,c4   
        !  funcRMSQ=(x1-(x1+x2+x3+x4)/4.)**2
        !& +(y1-(y1+y2+y3+y4)/4.)**2
        !& +(z1-(z1+z2+z3+z4)/4.)**2 
    do l=1,3   
	    rc(l)=sum(x(l,:))/4
    enddo
	    rmsq=sum((x(:,1)-rc(:))**2)
	    rm=sqrt(rmsq)
    return
    end subroutine funcrmsq
   
    subroutine samplerhodistrbt(x,i,rhodistnowout) 
    !,rhodistout,rhodisterrout) ! i is the picked the bead
    use mympi
    real(kind=r8), PARAMETER :: pi=4.0_r8*atan(1.0_r8) 
    real(kind=r8), dimension(3,npart) :: x,f
    real(kind=r8) :: r,xnk,xkh1,xkh2,basis_k
    integer(kind=i4) :: i,j,k,l,iknum
    real(kind=r8), dimension(3) :: dx,sumcor,rc 
    real(kind=r8), dimension(0:nrhobin) :: rhodistout,rhodistnow,rhodisterrout,rhodistnowout
   
    rhodistnow=0
    rhodistnowout=0
    do l=1,3   
	    rc(l)=sum(x(l,:))/4
    enddo
    r=sqrt(sum((x(:,1)-rc(:))**2))
    iknum=int(r/rhobinsize)    
    xnk=4.0_r8/3*pi*rhobinsize**3*((iknum+1)**3-iknum**3)   
    if (iknum.le.nrhobin) then  	  
        xkh1=iknum*rhobinsize
        xkh2=(iknum+1)*rhobinsize
        IF(  (r.GE.xkh1).and.(r.LE.xkh2)  ) THEN
                basis_k=1
        ELSE
                basis_k=0
        END IF			  	  
        rhodistnow(iknum)=basis_k/xnk  	  
	    rhodistnowout=rhodistnow
            
        if ((irep.eq.5).or.(irep.eq.6)) then
            rhodistblk(iknum,i)=rhodistblk(iknum,i)+rhodistnow(iknum)  
            rhodist2blk(iknum,i)=rhodist2blk(iknum,i)+rhodistnow(iknum)**2  
            wtrhodistblk(i)=wtrhodistblk(i)+1          
        else   
            icrhodist(i)=icrhodist(i)+1
            rhodisttot(iknum,i)=rhodisttot(iknum,i)+rhodistnow(iknum)  
	        rhodist2tot(iknum,i)=rhodist2tot(iknum,i)+rhodistnow(iknum)**2   
	        rhodist(:,i)=rhodisttot(:,i)/icrhodist(i)
            !rhodistout(:)=rhodist(:,i)
	        rhodisterr(:,i)=sqrt((rhodist2tot(:,i)/icrhodist(i)-rhodist(:,i)**2)/max(1,icrhodist(i)-1))
	        !rhodisterrout=rhodisterr(:,i)                  
        endif
    
    endif  
    return
    end subroutine samplerhodistrbt
   
    subroutine updaterhodistrbt(rhodistout,rhodisterrout,i) 
    use mympi
    integer(kind=i4) :: i
    real(kind=r8), dimension(0:nrhobin) :: rhodistout,rhodistall,rhodisterrout,rhodisterrall &
                                            ,rhodistottall,rhodist2totall
    real(kind=r8), dimension(0:nchorizo) :: wtrhodistsum,av,av2
    real(kind=r8), dimension(0:nrhobin,0:nchorizo) :: rhodistsum,rhodist2sum
   
    if ((irep.eq.5).or.(irep.eq.6)) then
       
        call addall(wtrhodistblk(i),wtrhodistsum(i))
        call addall(rhodistblk(:,i),rhodistsum(:,i))
        call addall(rhodist2blk(:,i),rhodist2sum(:,i))
     
        if (myrank().eq.0) then
         
           wtrhodistblk(i)=wtrhodistsum(i) 
           rhodistblk(:,i)=rhodistsum(:,i) 
           rhodist2blk(:,i)=rhodist2sum(:,i) 
                       
           wtrhodisttot(i)=wtrhodisttot(i)+wtrhodistblk(i)
           rhodisttot(:,i)=rhodisttot(:,i)+rhodistblk(:,i)
           rhodist2tot(:,i)=rhodist2tot(:,i)+rhodist2blk(:,i)
   
        rhodistout(:)=rhodisttot(:,i)/wtrhodisttot(i)
        av=rhodistout
        if (wtrhodisttot(i) /= 0 ) then
            av2=rhodist2tot(:,i)/wtrhodisttot(i) 
        else
            av2=0
        endif   
        rhodisterrout(:)=sqrt(abs(av2(:)-av(:)**2)/max(1.0_r8,wtrhodisttot(i)-1))        
      endif
        wtrhodistblk(i)=0 
        rhodistblk(:,i)=0
        rhodist2blk(:,i)=0
 
    else ! this is for the modes in which all the core are always involved in calculations and samplings.
        call addall(rhodist(:,i),rhodistall(:))
        call addall(rhodisterr(:,i)**2,rhodisterrall(:))  
        if (myrank().eq.0) then
            rhodistout(:)=rhodistall(:)/nproc()
            rhodisterrout(:)=sqrt(rhodisterrall(:))/nproc()  ! this does not seem right. fix       
        endif        
  
    endif
   
    return
    end subroutine updaterhodistrbt
   

   
 
   
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
