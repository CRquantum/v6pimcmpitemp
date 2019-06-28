module math
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
contains	
    subroutine combin(n,m,nout) ! C_n^m = n!/(m!*(n-m)!), m <= n
    integer(kind=i4) :: n,m,i,nout
   
    if ((n<0).or.(m<0)) then
    write (6,*) 'm and/or n <0, stop!'
    stop    
    else if (m>n) then
    write (6,*) 'm bigger than n, stop!'
    stop
    else if (m==n) then
    nout=1
    else
    !nout=gamma(real(n+1))/(gamma(real(m+1))*gamma(real(n-m+1)))   
    nout=product((/(i,i=m+1,n)/))/product((/(i,i=1,n-m)/))
    endif
    return
    end subroutine combin

    subroutine permut(n,nout) ! nout=n!
    !  can use gamma function, for integer x, gamma(x)=(x-1)!
    integer(kind=i4) :: n,i,nout   
    if (n<0) then
    write (6,*) 'n < 0, stop!'
    stop    
    else
    !nout=gamma(real(n+1))
    nout=product((/(i,i=1,n)/))
    endif
   
    return
    end subroutine permut  
   
    subroutine timedhms(timeall,d,h,m,s) ! timeall in seconds
    real(kind=r8) :: timeall,s
    integer(kind=i4) :: d,h,m 
    d=int(timeall/3600/24)
    h=int((timeall-int(timeall/3600/24)*3600*24)/3600) 
    m=int((timeall-int(timeall/3600/24)*3600*24-int((timeall-int(timeall/3600/24)*3600*24)/3600)*3600)/60)
    s=timeall-int(timeall/3600/24)*3600*24-int((timeall-int(timeall/3600/24)*3600*24)/3600)*3600 &
        -int((timeall-int(timeall/3600/24)*3600*24-int((timeall-int(timeall/3600/24)*3600*24)/3600)*3600)/60)*60    
    return
    endsubroutine timedhms   
   
 
end module math
