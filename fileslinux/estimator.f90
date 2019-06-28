module estimator
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   character(len=30), private, dimension(:), allocatable, save :: label
   real(kind=r8), private, dimension(:), allocatable, save :: &
      valtot,val2tot,val2,valblk,valnow,avbad,av2bad,wttot,wtblk
   integer(kind=i4), private, save :: nblock,nest,irep
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
    subroutine setestnum(n,irepin)
    !
    ! set up arrays for n estimators
    !
    integer(kind=i4) :: n,irepin
    nest=n
    irep=irepin
    allocate(valtot(n),val2tot(n),valblk(n),val2(n),valnow(n),avbad(n),av2bad(n))
    allocate(wttot(n),wtblk(n))
    allocate(label(n))
    return
    end subroutine setestnum

    subroutine zerest
    !
    ! zero all estimators
    !
    valtot=0
    val2tot=0
    valblk=0
    val2=0
    valnow=0
    avbad=0
    av2bad=0
    wtblk=0
    wttot=0
    nblock=0
    return
    end subroutine zerest

    subroutine addest(i,l)
    !
    ! set label l for estimator i
    !
    integer(kind=i4) :: i
    character(len=*) :: l
    integer(kind=i4) :: ln,j
    ln=min(len(l),len(label(i)))
    label(i)(1:ln)=l(1:ln)
    do j=ln+1,len(label(i))
        label(i)(j:j)=" "
    enddo
    return
    end subroutine addest

    subroutine addval(i,val,wt)
    !
    ! add another value, val, with weight, wt to estimator i
    !
    integer(kind=i4) :: i
    real(kind=r8) :: val,wt
    valblk(i)=valblk(i)+val*wt
    val2(i)=val2(i)+(val*wt)**2 ! calc mode.
    wtblk(i)=wtblk(i)+wt
    return
    end subroutine addval

    subroutine update
    use mympi
    !
    ! block all estimators
    !
    integer(kind=i4) :: i
    real(kind=r8) :: valsum(nest),wtsum(nest),val2sum(nest)
   
    call addall(valblk,valsum)
    call addall(val2,val2sum)
    call addall(wtblk,wtsum)

    !valsum=valblk
    !wtsum=wtblk
     
    if (myrank().eq.0) then
        valblk=valsum
        wtblk=wtsum 
        val2=val2sum
      
        do i=1,nest
            if (wtblk(i).eq.0) then
            valnow(i)=0
            else		   
   	        if (((irep.eq.2).or.(irep.eq.4).or.(irep.eq.5).or.(irep.eq.6)).and.(i.eq.6)) then
                valnow(i)=sqrt(valblk(i)/wtblk(i)) ! calculate rm=sqrt(ave rm^2) 
		    !write (6,*) 'SQRT!!!!',valnow(i)
	        else
		        valnow(i)=valblk(i)/wtblk(i) 
		    endif
		    endif
		    avbad(i)=avbad(i)+valnow(i)
            av2bad(i)=av2bad(i)+valnow(i)**2
            wttot(i)=wttot(i)+wtblk(i)
            valtot(i)=valtot(i)+valblk(i)	
            val2tot(i)=val2tot(i)+val2(i)
        enddo
    endif
   
    valblk=0
    wtblk=0
    val2=0
   
    nblock=nblock+1
    return
    end subroutine update
 
    subroutine resultest(i,av,err) ! return current average and err in val and err
    use mympi
    integer(kind=i4) :: i
    real(kind=r8) :: val,err
    real(kind=r8) :: av,av2
    if (wttot(i).eq.0) then
        val=0
    else
        val=valtot(i)/wttot(i)
    endif
    if ( (irep.eq.5).or.(irep.eq.6)  ) then 
    ! need to use valtot/wttot, block average does not apply
    ! bc the way I arrange the cores is different from pimc. unless cores are perfect match.  
        if (i.eq.6) val=sqrt(val)
        av=val
        if (wttot(i) /= 0 ) then
            av2=val2tot(i)/wttot(i)
        else
            av2=0
        endif
        if (wttot(i).eq.1) then
            err=0
        else    
            err=sqrt(abs(av2-av**2)/max(1.0_r8,wttot(i)-1))
        endif   
    else      
        av=avbad(i)/nblock
        av2=av2bad(i)/nblock
        if (nblock.eq.1) then
            err=0 
        else
            err=sqrt(abs(av2-av**2)/max(1,nblock-1))    
        endif              
    endif
    return
    end subroutine resultest

    function resstring()
    ! return all results in a string for printing
    use mympi
    character(len=120), dimension(nest) :: resstring
    integer(kind=i4) :: i,iounit
    real(kind=r8) :: val,err,val1,val2,err1,err2,error
    do i=1,nest
        call resultest(i,val,err) 
        if (i.eq.3) then 
            val1=val
            err1=err
        endif
        if (i.eq.4) then
            val2=val
            err2=err
        endif   
        if (i.ne.4) then
	        ! write (6,*) nest,i,label(i),val,err  
            write (resstring(i),'(a20,t20,1x,g15.7,'' +-'',g15.7)') label(i),val,err
        else
    ! error for function like A/B, given delta A anad delta B.   val1=A,val2=B       
            error=sqrt((1./val2*err1)**2+(val1/(val2**2)*err2)**2) 
	    !write (6,*) nest,i,label(i),val,err,val1/val2,error  
            write (resstring(i),'(a20,t20,1x,g15.7,'' +-'',g15.7,1x,''|'',g15.7,'' +-'',g15.7)') &
                label(i),val,err,val1/val2,error  
            if (myrank().eq.0)  then      
                open(newunit=iounit,FILE='energy_evolution.txt',FORM='FORMATTED',position='append') 
                write (iounit,'(g15.7,'' +-'',g15.7)') val1/val2,error  
                close(iounit)
            endif
        endif      
    enddo
    return
    end function resstring
   
end module estimator
