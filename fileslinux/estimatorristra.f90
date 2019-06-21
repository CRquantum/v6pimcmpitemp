!$Id: estimatorristra.f90,v 1.1 2003/02/24 16:15:38 nuclear Exp $
module estimatorristra
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   character(len=30), private, dimension(:), allocatable, save :: label
   real(kind=r8), private, dimension(:,:), allocatable, save :: &
      valtot,valblk,valnow,avbad,av2bad,val2,val2tot
   real(kind=r8), private, dimension(:), allocatable, save :: wttot,wtblk
   integer(kind=i4), private, save :: nblock,nest,nchorizo,irep
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
!
! set up arrays for n estimators
!
   integer(kind=i4) :: n,nchorizoin,irepin
   nest=n
   nchorizo=nchorizoin
   irep=irepin
   allocate(valtot(0:nchorizo,n),valblk(0:nchorizo,n))
   allocate(valnow(0:nchorizo,n),avbad(0:nchorizo,n),av2bad(0:nchorizo,n))
   allocate(val2(0:nchorizo,n),val2tot(0:nchorizo,n))
   allocate(wttot(n),wtblk(n))
   allocate(label(n))
   return
   end subroutine setestnumristra

   subroutine zerestristra
!
! zero all estimators
!
   valtot=0.0_r8
   val2tot=0.0_r8
   val2=0.0_r8
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
   end subroutine addestristra

   subroutine addvalristra(i,val)
! add another value, val, with weight, to estimator i
   integer(kind=i4) :: i
   real(kind=r8) :: val(0:nchorizo)
   valblk(:,i)=valblk(:,i)+val(:)
   wtblk(i)=wtblk(i)+1
   val2(:,i)=val2(:,i)+val(:)**2
   return
   end subroutine addvalristra

   subroutine updateristra
   use mympi
! block all estimators
   integer(kind=i4) :: i
   real(kind=r8) :: valsum(0:nchorizo,nest),wtsum(nest),val2sum(0:nchorizo,nest)
   call addall(valblk,valsum)
   call addall(wtblk,wtsum)
   call addall(val2,val2sum)
    !valsum=valblk
    !wtsum=wtblk
   if (myrank().eq.0) then
      valblk=valsum
      wtblk=wtsum
      val2=val2sum
      do i=1,nest
         valnow(:,i)=valblk(:,i)/wtblk(i)
      enddo
      avbad=avbad+valnow
      av2bad=av2bad+valnow**2
      wttot=wttot+wtblk
      valtot=valtot+valblk
      val2tot=val2tot+val2
   endif
   valblk=0.0_r8
   wtblk=0.0_r8
   val2=0.0_r8
   nblock=nblock+1
   return
   end subroutine updateristra
 
   subroutine resultristra(i,vnow,val,err,labelristra)
! return current result, value for current block in vnow, and
! current average and error in val and err
   integer(kind=i4) :: i
   real(kind=r8), dimension(0:nchorizo) :: vnow,val,err
   real(kind=r8), dimension(0:nchorizo) :: av,av2
   character(len=*) :: labelristra
   labelristra=label(i)
   vnow(:)=valnow(:,i) 
   val(:)=valtot(:,i)/wttot(i)
   if ( (irep.eq.5).or.(irep.eq.6)  ) then   
       av=val  
       if (wttot(i) /= 0 ) then
           av2=val2tot(:,i)/wttot(i)  ! check.
       else
           av2=0
       endif   
       err=sqrt(abs(av2-av**2)/max(1.0_r8,wttot(i)-1))           
   else 
       av=avbad(:,i)/nblock
       av2=av2bad(:,i)/nblock
       err=sqrt(abs(av2-av**2)/max(1,nblock-1))
   endif
   return
   end subroutine resultristra

end module estimatorristra
