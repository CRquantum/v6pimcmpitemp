!$Id: lj.f90,v 1.2 2003/02/25 18:11:59 nuclear Exp $
module lj
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), private, parameter :: epsilon = 10.22d0
   real(kind=r8), private, parameter :: sigma = 2.556d0
   real(kind=r8), private, save :: bmcmillan,el,rfit,alpha,u0,el2
   real(kind=r8), private, save :: azabo,bzabo,dzabo   
   integer(kind=i4), private, save :: mmcmillan

contains
   subroutine ljinit(bin,elin,min)
   real(kind=r8) :: bin,elin
   real(kind=r8) :: em1,em3,ufit
   integer(kind=i4) :: min
   bmcmillan=bin
   mmcmillan=min
   el=elin
   el2=0.5d0*el
   em1=mmcmillan+1
   em3=mmcmillan+3
   rfit=em1*el2/em3
   ufit=0.5d0*(bmcmillan/rfit)**mmcmillan
   alpha=-mmcmillan*ufit/(3.0_r8*rfit*(rfit-el2)**2)
   u0=alpha*(rfit-el2)**3-ufit
   return
   end subroutine ljinit

   subroutine mcmillan(r,u,up,upp)
   real(kind=r8) :: r,u,up,upp
   if (r.gt.el2) then
      u=0.0_r8
      up=0.0_r8
      upp=0.0_r8
   else if (r.lt.rfit) then
      u=-0.5d0*(bmcmillan/r)**mmcmillan
      up=-mmcmillan*u/r**2
      upp=mmcmillan*(mmcmillan-1)*u/r**2
      u=u-u0
   else
      u=-alpha*(r-el2)**3
      up=-alpha*3.0_r8*(r-el2)**2/r
      upp=-alpha*6.0_r8*((r-el2)+(r-el2)**2/r)
   endif
   return
   end subroutine mcmillan
    
   subroutine ljzaboinit(bin)
   real(kind=r8) :: bin,elin
   real(kind=r8) :: em1,em3,ufit
   integer(kind=i4) :: irepin
   bzabo=bin
   return
   end subroutine ljzaboinit
     
   subroutine mcmillanzabo(r,u) ! it is only used in wavefunctioninit
   real(kind=r8) :: r,u,up,upp
   !if (r.gt.22_r8/bzabo) then
   !   u=0.0_r8
   !else
	  !u=exp(-bzabo*r)
   !endif
   u=exp(-bzabo*r)
   return
   end subroutine mcmillanzabo
   
   subroutine vlj(r,v)
   real(kind=r8) :: r,v
   real(kind=r8), parameter :: hbar=197.327
   !if (r.eq.0.0_r8) r=1.0d-50 ! I set it as 1d-50 instead of 1d-20.
   !if (r.gt.20.) r=20.
! need to check if need to multiply 2m/hbar^2 or something, I guess no need.
! This is Malfliet-Tjon interaction.
    v=hbar*(7.39_r8*exp(-3.11_r8*r)/r-2.93_r8*exp(-1.55_r8*r)/r)
	!v=0.
   return
   end subroutine vlj
   
   subroutine ljtail(r,tail)
   use volume
!
! returns the integral of the potential from r to infinity
! multiply by rho/2 for energy tail
!
   real(kind=r8) :: r,tail,pi
   real(kind=r8) :: rmax,volp,volm,rp,rm,vp,vm,dr
   integer(kind=i4) :: nmax,i
   nmax=10000
   rmax=100.0_r8
   dr=(rmax-r)/nmax
   tail=0.0_r8
   pi=4.0_r8*atan(1.0_r8)
   do i=1,nmax
      rm=(i-1)*dr+r
      rp=i*dr+r
      volp=4.0_r8*pi*rp**3/3.0_r8-sphere(rp,r,r,r)
      volm=4.0_r8*pi*rm**3/3.0_r8-sphere(rm,r,r,r)
      call vlj(rp,vp)
      call vlj(rm,vm)
      tail=tail+(volp-volm)*0.5d0*(vp+vm)
   enddo
   return
   end subroutine ljtail

end module lj
