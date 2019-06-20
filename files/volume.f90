!$Id: volume.f90,v 1.2 2004/06/22 15:03:44 nuclear Exp $
module volume
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), private, parameter :: fuzz = 1.0d-10

contains
   function sphere(r,elx,ely,elz)
!
! calculate the intersection volume of a sphere of radius r and a parallepiped
! of sides 2*elx, 2*ely, 2*elz
!
! elx, ely, elz are changed to a,b, and c with a < b < c.
! The area of a circle of radius r inside a rectangle with
! sides 2*a,2*b is
!
!  A(r)  = ( pi*r**2;       r < a
!          (
!          ( 2*a*sqrt(r**2-a**2)+2*r**2*arcsin(a/r);   a < r < b
!          (
!          ( 2*a*sqrt(r**2-a**2)+2*r**2*arcsin(a/r)
!          (+2*b*sqrt(r**2-b**2)+2*r**2*arcsin(b/r)
!          (-pi*r**2;      b < r < sqrt(a**2+b**2)
!          (
!          ( 4*a*b    r > sqrt(a**2+b**2)
!
!
! The volume of a sphere in 3-d is
!
!  V(r) = 2*int from 0 to min(rin,c) A(sqrt(r**2-z**2)) dz
!
! the four integrals are given in routines vint1 through vint4
!
   real(kind=r8) :: r,elx,ely,elz,sphere
   real(kind=r8) :: a,b,c,a2,b2,c2,r2,el0,el1,el2,el3,el4,vol,rr
   a=min(elx,ely,elz)
   c=max(elx,ely,elz)
   b=elx+ely+elz-a-c
   a2=a*a
   b2=b*b
   c2=c*c
   r2=min(r*r,a2+b2+c2)
   rr=sqrt(r2)
   el0=0.0_r8
   el1=max(0.0_r8,min(r2-a2-b2,c2))
   el2=max(0.0_r8,min(r2-b2,c2))
   el3=max(0.0_r8,min(r2-a2,c2))
   el4=max(0.0_r8,min(r2,c2))
   el1=sqrt(el1)
   el2=sqrt(el2)
   el3=sqrt(el3)
   el4=sqrt(el4)
   vol=0.0_r8
   if (abs(el0-el1).gt.fuzz) vol=vol+vint4(rr,el1,a,b)-vint4(rr,el0,a,b)
   if (abs(el1-el2).gt.fuzz) vol=vol+vint3(rr,el2,a,b)-vint3(rr,el1,a,b)
   if (abs(el2-el3).gt.fuzz) vol=vol+vint2(rr,el3,a)-vint2(rr,el2,a)
   if (abs(el3-el4).gt.fuzz) vol=vol+vint1(rr,el4)-vint1(rr,el3)
   sphere=vol*2.0_r8
   return
   end function sphere

   function vint1(r,z)
!
! evaluate the integral A1 = int dz pi (r^2-z^2)
!
   real(kind=r8) :: r,z,pi,vint1
   pi=4.0_r8*atan(1.0_r8)
   vint1=pi*z*(r*r-z*z/3.0_r8)
   return
   end function vint1

   function vint2(r,z,a)
!
! evaluate the integral A2 = 2 int dz a*sqrt(r^2-a^2-z^2)
!                            +(r^2-z^2) arcsin(a/sqrt(r^2-z^2))

   real(kind=r8) :: r,z,a,vint2,r2,r3,z2,a2,rt,r2p,r2m
   r2=r*r
   r3=r2*r
   z2=z*z
   a2=a*a
   rt=sqrt(max(r2-a2-z2,0.0_r8)) !max protects against negative from roundoff
   r2p=r2-a2+r*z
   r2m=r2-a2-r*z
   vint2=a*z*rt*2.0_r8/3.0_r8 &
         +z*(r2-z2/3.0_r8)*atan2(a,rt) &
         +r3*(atan2(a*rt,r2p)-atan2(a*rt,r2m))/3.0_r8 &
         +a*(a2/3.0_r8-r2)*atan2(z,-rt)
   vint2 = 2.0_r8*vint2;
   return
   end function vint2

   function vint3(r,z,a,b)
!
! evaluate the integral A3 = 2 int dz a*sqrt(r^2-a^2-z^2)
!                            +(r^2-z^2) arcsin(a/sqrt(r^2-z^2))
!                           +2 int dz a*sqrt(r^2-a^2-z^2)
!                            +(r^2-z^2) arcsin(a/sqrt(r^2-z^2))
!                            -pi*(r*r-z*z)
!

   real(kind=r8) :: r,z,a,b,vint3,a1,a2,a3
   a1=vint1(r,z)
   a2=vint2(r,z,a)
   a3=vint2(r,z,b)
   vint3=a2+a3-a1
   return
   end function vint3

   function vint4(r,z,a,b)
!
! evaluate the integral A4 = 4 int dz ab
!
   real(kind=r8) :: r,z,a,b,vint4
   vint4=4.0_r8*a*b*z
   return
   end function vint4
end module volume
