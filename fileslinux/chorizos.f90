!$Id: chorizos.f90,v 1.1.1.1 2003/02/23 07:47:06 nuclear Exp $
module chorizos
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)

   type :: chorizo
      real(kind=r8), allocatable :: x(:,:)
      real(kind=r8), allocatable :: dpsi(:,:)
      real(kind=r8) :: psiln,d2psi,v,elocal
   end type

   integer(kind=i4), private, save :: npart

   type (chorizo), public, save :: chorizo1,chorizo2

   interface assignment (=)
      module procedure copychorizo
   end interface

!  +,-,*,/ are not defined, so do not use them for now.   

contains
   subroutine copychorizo(cl,cr)
   type (chorizo), intent(inout) :: cl !left side of assignment. here need to use intent(inout)
   type (chorizo), intent(in) :: cr !right side of assignment
   cl%x(:,:)=cr%x(:,:)
   cl%dpsi(:,:)=cr%dpsi(:,:)
   cl%psiln=cr%psiln
   cl%d2psi=cr%d2psi
   cl%v=cr%v
   cl%elocal=cr%elocal
   return
   end subroutine copychorizo

   subroutine chorizoinit(npartin)
   integer(kind=i4) :: npartin
   npart=npartin
   allocate(chorizo1%x(3,npart),chorizo1%dpsi(3,npart))
   allocate(chorizo2%x(3,npart),chorizo2%dpsi(3,npart))
   return
   end subroutine chorizoinit

   function npartchorizo()
   integer(kind=i4) :: npartchorizo
   npartchorizo=npart
   return
   end function npartchorizo

end module chorizos
