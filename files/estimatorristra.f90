!$Id: estimatorristra.f90,v 1.1 2003/02/24 16:15:38 nuclear Exp $
module estimatorristra
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   character(len=30), private, dimension(:), allocatable, save :: label
   real(kind=r8), private, dimension(:,:), allocatable, save :: &
      valtot,valblk,valnow,avbad,av2bad,val2blk,val2tot
   real(kind=r8), private, dimension(:), allocatable, save :: wttot,wtblk
   integer(kind=i4), private, save :: nblock,nest,nchorizo,irep
! response functions
   character(len=30), private, dimension(:), allocatable, save :: labele0resp,labeleresp
   integer(kind=i4), private, save :: nblocke0resp,ne0respbeads,nrespest,irlintervalmax,irlintervalmax1,irlstep1,irlstep2
   real(kind=r8), private, dimension(:,:), allocatable, save :: &
       vale0resptot,vale0respnow,vale0respblk,vale0resp2blk,vale0resp2tot,ave0respbad,ave0resp2bad
   real(kind=r8), private, dimension(:), allocatable, save :: wte0respblk,wte0resptot
!
! label = label for estimator
! valtot = sum of the valblk values for blocks
! valblk = current sum of value for block
! valnow = current average for block
! avbad = sum of (valblk/wtblk)
! av2bad = sum of (valblk/wtblk)^2
! wttot = sum of wtblk
! wtblk = current sum of weights for block
! nblock = number of blocks averaged
! nest = number of estimators
!
contains
    subroutine setestnumristra(n,nchorizoin,irepin)
    ! set up arrays for n estimators
    integer(kind=i4) :: n,nchorizoin,irepin
    nest=n
    nchorizo=nchorizoin
    irep=irepin
    allocate(valtot(0:nchorizo,n),valblk(0:nchorizo,n))
    allocate(valnow(0:nchorizo,n),avbad(0:nchorizo,n),av2bad(0:nchorizo,n))
    allocate(val2blk(0:nchorizo,n),val2tot(0:nchorizo,n))
    allocate(wttot(n),wtblk(n))
    allocate(label(n))
    return
    end subroutine setestnumristra

    subroutine zerestristra
    ! zero all estimators
    valtot=0.0_r8
    val2tot=0.0_r8
    val2blk=0.0_r8
    valblk=0.0_r8
    valnow=0.0_r8
    avbad=0.0_r8
    av2bad=0.0_r8
    wtblk=0.0_r8
    wttot=0.0_r8
    nblock=0
    return
    end subroutine zerestristra

    subroutine addestristra(i,l)
    ! set label l for estimator i
    integer(kind=i4) :: i
    character(len=*) :: l
    integer(kind=i4) :: ln,j
    ln=min(len(l),len(label(i)))
    label(i)(1:ln)=l(1:ln)
    do j=ln+1,len(label(i))
        label(i)(j:j)=" "
    enddo
    return
    end subroutine addestristra

    subroutine addvalristra(i,val)
    ! add another value, val, with weight, to estimator i
    integer(kind=i4) :: i
    real(kind=r8) :: val(0:nchorizo)
    valblk(:,i)=valblk(:,i)+val(:)
    wtblk(i)=wtblk(i)+1
    val2blk(:,i)=val2blk(:,i)+val(:)**2
    return
    end subroutine addvalristra
 
    subroutine updateristra
    use mympi
    ! block all estimators
    integer(kind=i4) :: i
    real(kind=r8) :: valsum(0:nchorizo,nest),wtsum(nest),val2sum(0:nchorizo,nest)
    call addall(valblk,valsum)
    call addall(wtblk,wtsum)
    call addall(val2blk,val2sum)
    !valsum=valblk
    !wtsum=wtblk
    if (myrank().eq.0) then
        valblk=valsum
        wtblk=wtsum
        val2blk=val2sum
        
        do i=1,nest 
            if (wtblk(i)/=0) then 
                valnow(:,i)=valblk(:,i)/wtblk(i)
            else
                valnow(:,i)=0
            endif 
        enddo        
 
        wttot=wttot+wtblk
        valtot=valtot+valblk
        val2tot=val2tot+val2blk
        
    endif
    valblk=0.0_r8
    wtblk=0.0_r8
    val2blk=0.0_r8
    
    return
    end subroutine updateristra
   
    subroutine checkpointestristra(input,wttotio,valtotio,val2totio)
    real(kind=r8) :: wttotio(nest),valtotio(0:nchorizo,nest),val2totio(0:nchorizo,nest)
    logical :: input 
    if (input) then
        wttot(:)=wttotio(:)
        valtot(:,:)=valtotio(:,:)
        val2tot(:,:)=val2totio(:,:)          
    else
        wttotio(:)=wttot(:)
        valtotio(:,:)=valtot(:,:)
        val2totio(:,:)=val2tot(:,:)        
    endif
    return
    end subroutine checkpointestristra   
   
    subroutine resultristra(i,vnow,val,err,labelristra)
    ! return current result, value for current block in vnow, and
    ! current average and error in val and err
    integer(kind=i4) :: i
    real(kind=r8), dimension(0:nchorizo) :: vnow,val,err
    real(kind=r8), dimension(0:nchorizo) :: av,av2
    character(len=*) :: labelristra
    labelristra=label(i)
    vnow(:)=valnow(:,i) 
    if (wttot(i) /= 0) then
        val(:)=valtot(:,i)/wttot(i)
    else
        val=0
    endif
    av=val  
    if (wttot(i) /= 0 ) then
        av2(:)=val2tot(:,i)/wttot(i)  ! check.
    else
        av2=0
    endif   
    err=sqrt(abs(av2-av**2)/max(1.0_r8,wttot(i)-1))           
    return
    end subroutine resultristra
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   
    subroutine writevtotbeads(dt)
    use mympi
    use estimator
    integer(kind=i4) :: j,myunit
    real(kind=r8) :: dt
    real(kind=r8) :: val1,val2,valnn1,valem1
    real(kind=r8) :: vnow(0:nchorizo),valaverage(0:nchorizo),valerr(0:nchorizo) &
                ,vnnnow(0:nchorizo),valnnaverage(0:nchorizo),valnnerr(0:nchorizo) &
                ,vemnow(0:nchorizo),valemaverage(0:nchorizo),valemerr(0:nchorizo)
    real(kind=r8) :: err1,err2,error,errnn1,errem1,errornn,errorem
    character(len=30) :: vpropchorizo,vpropchorizotmp
    if (myrank().eq.0) then  

        call resultristra(1,vnow,valaverage,valerr,vpropchorizo)
        call resultristra(2,vnnnow,valnnaverage,valnnerr,vpropchorizotmp)
        call resultristra(3,vemnow,valemaverage,valemerr,vpropchorizotmp)
        open(newunit=myunit,form='formatted',file=vpropchorizo)
        call resultest(4,val2,err2) ! #4 is the denominator 
        do j=0,nchorizo
            val1=valaverage(j)
            err1=valerr(j)  
            valnn1=valnnaverage(j)
            errnn1=valnnerr(j) 
            valem1=valemaverage(j)
            errem1=valemerr(j)        
            error=sqrt((1/val2*err1)**2+(val1/(val2**2)*err2)**2)
            errornn=sqrt((1/val2*errnn1)**2+(valnn1/(val2**2)*err2)**2)
            errorem=sqrt((1/val2*errem1)**2+(valem1/(val2**2)*err2)**2)
            ! val1/val2 and error are the Vtot and error we want.
            write (myunit,'(i10,f10.5,14g15.6)') j,j*dt,val1/val2,error,valaverage(j),valerr(j) &
                                            ,valnn1/val2,errornn,valnnaverage(j),valnnerr(j) &
                                            ,valem1/val2,errorem,valemaverage(j),valemerr(j) &
                                            ,val2,err2
        enddo
        close(myunit) 
     
    endif      
    return 
    end subroutine writevtotbeads      
    
    
    
    
    
    
    !!!!!!!!!!!! response function !!!!!!!!!!!!!!!!!!!! 
   
    subroutine setestnumristraresp(nrespestin,irlintervalmaxin,irlintervalmax1in,irlstep1in,irlstep2in)
    ! set up arrays for n estimators
    integer(kind=i4) :: nrespestin,irlintervalmaxin,irlintervalmax1in,irlstep1in,irlstep2in
   
 
    nrespest=nrespestin ! set it to 5.
    irlintervalmax=irlintervalmaxin
    irlintervalmax1=irlintervalmax1in
    irlstep1=irlstep1in
    irlstep2=irlstep2in
   
   
    ne0respbeads=irlintervalmax/irlstep2  ! effective beads in the middle, which*2dt=tau
   
    allocate(vale0resptot(0:ne0respbeads,nrespest),vale0respnow(0:ne0respbeads,nrespest) &
            ,vale0resp2blk(0:ne0respbeads,nrespest),vale0resp2tot(0:ne0respbeads,nrespest) &
            ,vale0respblk(0:ne0respbeads,nrespest))
    allocate(ave0respbad(0:ne0respbeads,nrespest),ave0resp2bad(0:ne0respbeads,nrespest))
    allocate(wte0resptot(nrespest),wte0respblk(nrespest))
    allocate(labele0resp(nrespest))
   

    return
    end subroutine setestnumristraresp
                                
    subroutine zerestristrae0resp
    ! zero all estimators
    vale0resptot=0
    vale0resp2tot=0
    vale0respnow=0
    vale0resp2blk=0
    vale0respblk=0
    wte0resptot=0
    wte0respblk=0
    ave0respbad=0
    ave0resp2bad=0
   
   
    return
    end subroutine zerestristrae0resp

    subroutine addestristrae0resp(i,l)
    ! set label l for estimator i
    integer(kind=i4) :: i
    character(len=*) :: l
    integer(kind=i4) :: ln,j
    ln=min(len(l),len(labele0resp(i)))
    labele0resp(i)(1:ln)=l(1:ln)
    do j=ln+1,len(labele0resp(i))
        labele0resp(i)(j:j)=" "
    enddo
    return
    end subroutine addestristrae0resp

    subroutine addvalristrae0resp(i,vale0)
    ! add another value, val, with weight, to estimator i
    integer(kind=i4) :: i
    real(kind=r8) :: vale0(0:ne0respbeads)
    vale0respblk(:,i)=vale0respblk(:,i)+vale0(:)
    wte0respblk(i)=wte0respblk(i)+1
   
    !if (i == 1) write(6,*) 'result0:', i,vale0resp2blk(:,i)
   
    vale0resp2blk(:,i)=vale0resp2blk(:,i)+vale0(:)**2
   
    !if (i == 1) write(6,*) 'result1:', vale0(:)
    !if (i == 1) write(6,*) 'result2:', vale0(:)**2
   
    return
    end subroutine addvalristrae0resp
 
    subroutine updateristrae0resp
    use mympi
    ! block all estimators
    real(kind=r8) :: vale0respsum(0:ne0respbeads,nrespest),wte0respsum(nrespest),vale0resp2sum(0:ne0respbeads,nrespest)
    call addall(vale0respblk,vale0respsum)
    call addall(wte0respblk,wte0respsum)
    call addall(vale0resp2blk,vale0resp2sum)
   
    !write(6,*) 'result003:', vale0resp2sum(:,1)
   
    !write(6,*) 'result02:', vale0resp2tot(:,1)
   
    if (myrank().eq.0) then
        vale0respblk=vale0respsum
        wte0respblk=wte0respsum
        vale0resp2blk=vale0resp2sum
      
        !write(6,*) 'result004:', vale0resp2blk(:,1)
        wte0resptot=wte0resptot+wte0respblk
        vale0resptot=vale0resptot+vale0respblk
        !write(6,*) 'result02:', vale0resp2tot(:,1)
        vale0resp2tot=vale0resp2tot+vale0resp2blk
      
        !write(6,*) 'result00:', vale0resp2tot(:,1)
        !write(6,*) 'result01:', vale0resp2blk(:,1)
      
    endif
    vale0respblk=0
    wte0respblk=0
    vale0resp2blk=0
    return
    end subroutine updateristrae0resp
    
    subroutine checkpointresponse(input,wte0resptotio,vale0resptotio,vale0resp2totio)
    real(kind=r8) :: wte0resptotio(nrespest),vale0resptotio(0:ne0respbeads,nrespest),vale0resp2totio(0:ne0respbeads,nrespest)
    logical :: input 
    if (input) then
        wte0resptot(:)=wte0resptotio(:)
        vale0resptot(:,:)=vale0resptotio(:,:)
        vale0resp2tot(:,:)=vale0resp2totio(:,:)          
    else
        wte0resptotio(:)=wte0resptot(:)
        vale0resptotio(:,:)=vale0resptot(:,:)
        vale0resp2totio(:,:)=vale0resp2tot(:,:)        
    endif
    return
    end subroutine checkpointresponse   
    
    subroutine resultristrae0resp(i,vnow,val,err,labelristra)
    ! return current result, value for current block in vnow, and
    ! current average and error in val and err
    integer(kind=i4) :: i
    real(kind=r8), dimension(0:ne0respbeads) :: vnow,val,err
    real(kind=r8), dimension(0:ne0respbeads) :: av,av2
    character(len=*) :: labelristra
    labelristra=labele0resp(i)
    vnow(:)=vale0respnow(:,i) 
    if (wte0resptot(i) /= 0) then
        val(:)=vale0resptot(:,i)/wte0resptot(i)
    else
        val=0
    endif
  
    av=val  
    if (wte0resptot(i) /= 0 ) then
        av2(:)=vale0resp2tot(:,i)/wte0resptot(i)  ! check.
    else
        av2=0
    endif   
    err=sqrt(abs(av2-av**2)/max(1.0_r8,wte0resptot(i)-1))      
       
    !write(6,*) 'result1:', av2
    !write(6,*) 'result111:',i,vale0resp2tot(:,i)
    !write(6,*) 'result2:', av
    !write(6,*) 'result3:', wte0resptot(i)

    return
    end subroutine resultristrae0resp   
   
      
   

   

end module estimatorristra
