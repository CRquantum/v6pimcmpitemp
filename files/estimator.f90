module estimator
    implicit none
    integer, private, parameter :: i4=selected_int_kind(9)
    integer, private, parameter :: r8=selected_real_kind(15,9)
    character(len=30), private, dimension(:), allocatable, save :: label
    real(kind=r8), private, dimension(:), allocatable, save :: &
        valtot,val2tot,val2,valblk,valnow,avbad,av2bad,wttot,wtblk
    integer(kind=i4), private, save :: nblock,nest,irep
    character(len=120), private, save :: file_energy,file_radii
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

    subroutine estimatorfilesinit(irepin)
    integer(kind=i4) :: irepin
    select case (irepin)
    case (4)   
        write(file_energy,'("energy_evolution_VMC.txt")')  
        write(file_radii,'("radii_VMC.txt")')
    case (5)    
        write(file_energy,'("energy_evolution_CalcMode5.txt")')  
        write(file_radii,'("radii_CalcMode5.txt")')    
    case (6)    
        write(file_energy,'("energy_evolution_CalcMode6.txt")')  
        write(file_radii,'("radii_CalcMode6.txt")')
    case default
        write(file_energy,'("energy_evolution.txt")')  
        write(file_radii,'("radii.txt")') 
    end select      
    return
    end subroutine estimatorfilesinit
      
    subroutine estimatorfilesclean ! call estimatorfilesinit first
    open(unit=13,form='formatted',file=trim(file_energy),status='replace',position='rewind')
    close(13)
    open(unit=13,form='formatted',file=trim(file_radii),status='replace',position='rewind')
    close(13) 
    end subroutine estimatorfilesclean   

    subroutine zerest
    !
    ! zero all estimators
    !
    valtot=0
    val2tot=0
    wttot=0
    valblk=0
    val2=0
    valnow=0
    avbad=0
    av2bad=0
    wtblk=0 
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
    use mympi
    ! add another value, val, with weight, wt to estimator i
    integer(kind=i4) :: i,k
    real(kind=r8) :: val,wt
    valblk(i)=valblk(i)+val*wt
    val2(i)=val2(i)+(val*wt)**2 ! calc mode.
    wtblk(i)=wtblk(i)+wt
     
    !do k=1,nproc()-1
    !    if (myrank().eq.k) then
    !        write(6,*) 'myrank() i = ', myrank(), i
    !        write(6,'(''valblk val2 wtblk: '',t40, 3g15.7)') ,valblk(i),val2(i),wtblk(i)
    !        write(6,'(''val wt: '',t40, 2g15.7)') val,wt
    !         write(6,*)
    !    endif
    !    call barrier
    !    
    !enddo
 
    return
    end subroutine addval

    subroutine update
    use mympi
    ! block all estimators
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
            wttot(i)=wttot(i)+wtblk(i)
            valtot(i)=valtot(i)+valblk(i)	
            val2tot(i)=val2tot(i)+val2(i)
        enddo
    endif
   
    valblk=0
    wtblk=0
    val2=0

    return
    end subroutine update
     
    subroutine checkpointest(input,wttotio,valtotio,val2totio)
    real(kind=r8) :: wttotio(nest),valtotio(nest),val2totio(nest)
    logical :: input 
    if (input) then
        wttot(:)=wttotio(:)
        valtot(:)=valtotio(:)
        val2tot(:)=val2totio(:)          
    else
        wttotio(:)=wttot(:)
        valtotio(:)=valtot(:)
        val2totio(:)=val2tot(:)        
    endif
    return
    end subroutine checkpointest

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
    ! need to use valtot/wttot, block average does not apply
    ! bc the way I arrange the cores is different from pimc. unless cores are perfect match.  
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
    return
    end subroutine resultest

    function resstring() ! return all results in a string for printing
    use mympi
    character(len=120), dimension(nest) :: resstring,resstring0
    integer(kind=i4) :: i,iounit
    real(kind=r8) :: val,err,val1(nest),val2,err1(nest),err2,error(nest),factor   
    call calcresstring(resstring,val,err,val1,val2,err1,err2,error)   
    if (myrank().eq.0)  then      
        open(newunit=iounit,FILE=trim(file_energy),FORM='FORMATTED',position='append') 
        write (iounit,'(8(g15.7,'' +-'',g15.7,1x))') val1(3)/val2,error(4),val1(3),err1(3),val2,err2 &  ! energy 
                                                    ,val1(17)/val2,error(17),val1(18)/val2,error(18) &
                                                    ,val1(19)/val2,error(19),val1(20)/val2,error(20),val1(21)/val2,error(21)    
        close(iounit)
        open(newunit=iounit,FILE=trim(file_radii),FORM='FORMATTED',position='append')         
        write (iounit,'(6(g15.7,'' +-'',g15.7,1x))') sqrt(val1(6)/val2),error(6),val1(7)/val2,error(7) & ! general,rms, rmave
                                                            ,sqrt(val1(13)/val2),error(13),val1(14)/val2,error(14) & ! neutron
                                                            ,sqrt(val1(15)/val2),error(15),val1(16)/val2,error(16) ! proton  
        close(iounit)        
    endif       
    return     
    end 
     
    function resstring0() ! return all results in a string for printing
    use mympi
    character(len=120), dimension(nest) :: resstring,resstring0
    integer(kind=i4) :: i,iounit
    real(kind=r8) :: val,err,val1(nest),val2,err1(nest),err2,error(nest),factor
    call calcresstring(resstring0,val,err,val1,val2,err1,err2,error)       
    return
    end     
    
    subroutine calcresstring(resstring,val,err,val1,val2,err1,err2,error)
    character(len=120), dimension(nest) :: resstring
    integer(kind=i4) :: i,iounit
    real(kind=r8) :: val,err,val1(nest),val2,err1(nest),err2,error(nest),factor
    do i=1,nest
        call resultest(i,val,err) 
        val1(i)=val
        err1(i)=err      
        if (i==4) then
            val2=val
            err2=err  
        endif
        select case (i)        
        case (1:3,5,9:12) ! statistics, numerator
	        ! write (6,*) nest,i,label(i),val,err  
            write (resstring(i),'(a20,t20,1x,g15.7,'' +-'',g15.7)') label(i),val1(i),err1(i)        
        case (4) ! energy, denominator
	    !write (6,*) nest,i,label(i),val,err,val1/val2,error  
            error(i)=sqrt((1.0_r8/val2*err1(3))**2+(val1(3)/(val2**2)*err2)**2) 
            write (resstring(i),'(a20,t20,1x,g15.7,'' +-'',g15.7,1x,''|'',g15.7,'' +-'',g15.7)') &
                label(i),val2,err2,val1(3)/val2,error(i)  
    
        case (6,13,15) ! rms
            factor=0.5*sqrt(val2/val1(i))
            error(i)=factor*sqrt((1.0_r8/val2*err1(i))**2+(val1(i)/(val2**2)*err2)**2) 
            write (resstring(i),'(a20,t20,1x,g15.7,'' +-'',g15.7)') label(i),sqrt(val1(i)/val2),error(i)  
    
        case (7:8,14,16,17:21)   ! other physics quantities. 
            error(i)=sqrt((1.0_r8/val2*err1(i))**2+(val1(i)/(val2**2)*err2)**2) 
            write (resstring(i),'(a20,t20,1x,g15.7,'' +-'',g15.7)') label(i),val1(i)/val2,error(i)             
        case default
             ! do nothing.                
        end select        
    enddo     
    return
    end subroutine calcresstring

    !function resstring() ! return all results in a string for printing
    !use mympi
    !character(len=120), dimension(nest) :: resstring,resstring0
    !integer(kind=i4) :: i,iounit
    !real(kind=r8) :: val,err,val1(nest),val2,err1(nest),err2,error(nest),factor
    !do i=1,nest
    !    call resultest(i,val,err) 
    !    val1(i)=val
    !    err1(i)=err      
    !    if (i==4) then
    !        val2=val
    !        err2=err  
    !    endif
    !
    !    select case (i)        
    !    case (1:3,5,9:12) ! statistics, numerator
	   !     ! write (6,*) nest,i,label(i),val,err  
    !        write (resstring(i),'(a20,t20,1x,g15.7,'' +-'',g15.7)') label(i),val1(i),err1(i)        
    !    case (4) ! energy, denominator
	   ! !write (6,*) nest,i,label(i),val,err,val1/val2,error  
    !        error(i)=sqrt((1.0_r8/val2*err1(3))**2+(val1(3)/(val2**2)*err2)**2) 
    !        write (resstring(i),'(a20,t20,1x,g15.7,'' +-'',g15.7,1x,''|'',g15.7,'' +-'',g15.7)') &
    !            label(i),val2,err2,val1(3)/val2,error(i)  
    !
    !    case (6,13,15) ! rms
    !        factor=0.5*sqrt(val2/val1(i))
    !        error(i)=factor*sqrt((1.0_r8/val2*err1(i))**2+(val1(i)/(val2**2)*err2)**2) 
    !        write (resstring(i),'(a20,t20,1x,g15.7,'' +-'',g15.7)') label(i),sqrt(val1(i)/val2),error(i)  
    !
    !    case (7:8,14,16,17:20)   ! other physics quantities. 
    !        error(i)=sqrt((1.0_r8/val2*err1(i))**2+(val1(i)/(val2**2)*err2)**2) 
    !        write (resstring(i),'(a20,t20,1x,g15.7,'' +-'',g15.7)') label(i),val1(i)/val2,error(i)             
    !    case default
    !         ! do nothing.                
    !    end select        
    !
    !enddo
    !   
    !if (myrank().eq.0)  then      
    !    open(newunit=iounit,FILE=trim(file_energy),FORM='FORMATTED',position='append') 
    !    write (iounit,'(8(g15.7,'' +-'',g15.7,1x))') val1(3)/val2,error(4),val1(3),err1(3),val2,err2 &  ! energy 
    !                                                ,val1(17)/val2,error(17),val1(18)/val2,error(18) &
    !                                                ,val1(19)/val2,error(19),val1(20)/val2,error(20)  
    !    close(iounit)
    !    
    !    open(newunit=iounit,FILE=trim(file_radii),FORM='FORMATTED',position='append') 
    !    write (iounit,'(6(g15.7,'' +-'',g15.7,1x))') sqrt(val1(6)/val2),error(6),val1(7)/val2,error(7) & ! general
    !                                                        ,sqrt(val1(13)/val2),error(13),val1(14)/val2,error(14) & ! proton
    !                                                        ,sqrt(val1(15)/val2),error(15),val1(16)/val2,error(16) ! neutron  
    !    close(iounit)        
    !endif    
    !   
    !return
    !end function resstring    
    
    subroutine writeanswer(answerin)
    character(len=120) :: answerin(:) 
    integer(kind=i4) :: k
    select case (irep)
    case (4) ! VMC
        write (6,'(a120)') (answerin(k),k=2,8)
        write (6,'(a120)') (answerin(k),k=13,size(answerin))
        write (12,'(a120)') (answerin(k),k=2,8)
        write (12,'(a120)') (answerin(k),k=13,size(answerin))          
    case (5:6) ! Calc Mode
        write (6,'(a120)') (answerin(k),k=3,8)
        write (6,'(a120)') (answerin(k),k=13,size(answerin))
        write (12,'(a120)') (answerin(k),k=3,8)
        write (12,'(a120)') (answerin(k),k=13,size(answerin))    
    case default ! PIMC
        write (6,'(a120)') (answerin(k),k=1,size(answerin))
	    write (12,'(a120)') (answerin(k),k=1,size(answerin))     
    end select 
    return
    end subroutine writeanswer    
    
    
end module estimator
