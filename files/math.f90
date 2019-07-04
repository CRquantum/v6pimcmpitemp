module math
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, dimension(:,:), allocatable, save :: p3map,p4map,p5map,p6map,p7map ! permutation map. 
contains
    
    subroutine mathinit ! make the table. for stepinit.
    call allpermut(3,p3map)
    call allpermut(4,p4map)
    call allpermut(5,p5map)
    call allpermut(6,p6map)
    call allpermut(7,p7map) ! 7!, should be enough.
    return
    end subroutine mathinit   
    
    subroutine getpermutmap(n,pmap)
    integer(kind=i4) :: n,pmap(:,:)
    select case (n)
    case (3)
        pmap=p3map
    case (4)
        pmap=p4map
    case (5)
        pmap=p5map
    case (6)
        pmap=p6map
    case (7)
        pmap=p7map
    case default 
        write (6,'(''stop because permutation map has not been(or no need to be) calculated: '')') n
        stop     
    end select
    return
    end subroutine getpermutmap
    
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
    
    subroutine allpermut(n,a) 
    ! every time call arranger, should give an a(n) out of the n! possibilities.
    integer(kind=i4) :: n,p,ida(n),ic,k,i
    integer(kind=i4), allocatable :: a(:,:) ! make sure the input a is allocatable.
    call permut(n,p)
    if (.not.allocated(a)) allocate(a(n,p)) 
    ida=(/(k,k=1,n)/)
    ic=0
    do 
        ic=ic+1
        do i=1,n
        !write (6,"(i3,',')",advance = "no") ida(i)
        a(i,ic)=ida(i) 
        enddo
    !     repeat if not being finished yet, otherwise exit.
        if ( nextp(n,ida) ) then 
            cycle 
        else 
            exit 
        endif 
    enddo    
    return
    end subroutine allpermut
    
    subroutine allcombin(n,m,comball) ! n>m, pick m from n.
    integer, allocatable :: comball (:,:) 
    integer :: nc,n,m
    call combin(n,m,nc)
    allocate(comball(m,nc))
    call gen (1,n,m,comball)
    ! eg:  
    !  integer, allocatable :: comball (:,:) 
    !  n=10
    !  m=4 
    !  call allcombin(n,m,comball) 
    !  call combin(n,m,nc)  
    !  write (6,*) 'nc=', nc
    !  do i=1,nc     
    !     write (6, '( ''i='',i0, t20' // repeat (', 1x, i0', m) // ')') i, comball(:,i)
    !  enddo    
    return
    end subroutine allcombin    
 
    recursive subroutine gen (m,n_max,m_max,combout) ! use allcombin to call it.
    implicit none
    integer, intent (in) :: m
    integer :: n,mm1
    integer :: n_max,m_max
    integer, allocatable, save :: comb(:)
    integer :: combout(:,:)
    integer, save :: i=1
    if (.not.allocated(comb)) then
        allocate(comb(m_max))
        comb=0
    endif
    if (m > m_max) then
        combout(:,i)=comb
        i=i+1
    else
        do n = 1, n_max
            mm1=m-1
            if (mm1 == 0) mm1=1
        if ((m == 1) .or. (n > comb (mm1))) then
            comb (m) = n
            call gen (m + 1,n_max,m_max,combout)
        end if
        end do
    end if
    return
    end subroutine gen
    
    function nextp ( n, a ) ! for permutation.  http://rosettacode.org/wiki/Permutations#Fortran
    logical :: nextp
    integer,intent(in) :: n
    integer,dimension(n),intent(inout) :: a
!
!     local variables:
!
    integer i,j,k,t
!
    i = n-1
10 if ( a(i) .lt. a(i+1) ) goto 20
    i = i-1
    if ( i .eq. 0 ) goto 20
    goto 10
20 j = i+1
    k = n
30 t = a(j)
    a(j) = a(k)
    a(k) = t
    j = j+1
    k = k-1
    if ( j .lt. k ) goto 30
    j = i
    if (j .ne. 0 ) goto 40
!      
    nextp = .false.
!      
    return
!
40 j = j+1
    if ( a(j) .lt. a(i) ) goto 40
    t = a(i)
    a(i) = a(j)
    a(j) = t
!      
    nextp = .true.
!      
    return 
    end function    
    
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
