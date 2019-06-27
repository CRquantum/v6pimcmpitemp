module v6pot
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: npart,ntab,ibox
   real(kind=r8), private, allocatable, save, dimension(:,:,:) :: vtab,vtabls
   real(kind=r8), private, allocatable, save, dimension(:) :: vemtab
   real(kind=r8), private, save :: scale,range,el,eli,vcbox
   real(kind=r8), private, parameter :: amu=0.69953054_r8
   real(kind=r8), private, save :: c3b,a2p3b,a2s3b,a3p3b,ar3b
contains
   subroutine v6potinit(npartin,iboxin,elin,lpot,ntabin,vfact,lpotpr,iem)
   use cheft
   !use mympi
   real(kind=r8) :: vfact(8),elin
   integer(kind=i4) :: ntabin,lpot,i,npartin,iboxin,lpotpr
   real(kind=r8) :: dr,vv(18),vem(14)
   real(kind=r8) :: rr
   real(kind=r8) :: v0r,v0s,v0t,vr,vs,vt,kr,ks,kt
   real(kind=r8), parameter :: verytiny=1.0e-20_r8
   logical :: iem
   ibox=iboxin
   npart=npartin
   el=elin
   if (el.gt.0.0_r8) then
      eli=1.0_r8/el
      range=el*sqrt(3.0_r8)*0.5_r8*(2.0_r8*ibox+1.0_r8)
   else
      range=abs(el)
      ibox=0
      el=0.0_r8
      eli=0.0_r8
   endif
   !if (myrank().eq.0) write(6,'(''Range of the potential = '',f10.5)') range
   ntab=ntabin
   scale=ntab/range
   dr=1.0_r8/scale
   if (.not.allocated(vtab)) allocate(vtab(2,6,0:ntab))
   if (.not.allocated(vtabls)) allocate(vtabls(2,2,0:ntab))
   if (iem.and.(.not.allocated(vemtab))) allocate(vemtab(0:ntab))
     
! setup table for the second potential to calculate
   if (lpotpr.lt.26) then
      !call setpot(lpotpr,1.0_r8,1.0_r8,1.0_r8,1.0_r8,1.0_r8,1.0_r8,h2m,h2mcsb) ! if use call pot then need this.
      do i=0,ntab
         rr=i*dr
         !call pot(0,rr,vv,vp,ww)
		 call av18op(3,rr,vv)
         vtab(2,:,i)=vfact(1:6)*vv(1:6:1)
         vtabls(2,:,i)=vfact(7:8)*vv(7:8:1)
!        if (myrank().eq.0) then
!           write(80,'(9f25.15)') rr,(vtab(2,j,i),j=1,6),(vtabls(2,j,i),j=1,2)
!        endif
      enddo
   elseif (lpotpr.gt.100) then
      do i=0,ntab
         rr=i*dr
         call cheft_pot(lpotpr,rr,vv)
         vtab(2,:,i)=vfact(1:6:1)*vv(1:6:1)
         vtabls(2,:,i)=vfact(7:8:1)*vv(7:8:1)
      enddo
   endif
!
! vcbox is the contribution from our own image. This can be interpreted two
! ways. If this is toroidal boundary conditions, then the particle is its
! own image and sigma.sigma = 3 etc. It seems more physical to assume periodic
! boundary conditions so the images are other particles that happen to have
! the same spin/isospin. In that case sigma.sigma = 1 etc. We use that below.
!
   if (abs(lpot).lt.26) then
      ! call setpot(lpot,1.0_r8,1.0_r8,1.0_r8,1.0_r8,1.0_r8,1.0_r8,h2m,h2mcsb) ! if use call pot then need this.
      do i=0,ntab
         rr=i*dr
         !call pot(0,rr,vv,vp,ww)
		 call av18op(3,rr,vv) ! 3 means v6' here.
		 !vv3=abs(vv-vv2)
		 !if (abs(sum(vv3(1:6))).ge.1.d-5) then
		 !write(6,*) 'v6 seems different', i,rr
		 !write(6,'(''my v6='',t15,6f15.7,1x)') vv(1:6)
		 !write(6,'(''av18op v6='',t15,6f15.7,1x)') vv2(1:6)
		 !stop
		 !endif	
		 

         vtab(1,:,i)=vfact(1:6:1)*vv(1:6:1)
         vtabls(1,:,i)=vfact(7:8:1)*vv(7:8:1)
!        if (myrank().eq.0) then
!           write(79,'(9f25.15)') rr,(vtab(1,j,i),j=1,6),(vtabls(1,j,i),j=1,2)
!        endif 		 
		 !do j=1,6
		 ! write (12,*) i,vtab(1,j,i)
		 !enddo	
		 
		 
	  enddo
   elseif (abs(lpot).eq.26) then
! Minnesota potential:
      v0r=200.0_r8  ! MeV
      v0t=178.0_r8  ! MeV
      v0s=91.85_r8  ! MeV
      kr=1.487_r8  ! fm**-2
      kt=0.639_r8  ! fm**-2
      ks=0.465_r8  ! fm**-2
      do i=0,ntab
         rr=i*dr
         vr=v0r*exp(-kr*rr**2)
         vt=-v0t*exp(-kt*rr**2)
         vs=-v0s*exp(-ks*rr**2)
         vtab(1,1,i)=3.0_r8/8.0_r8*(vr+0.5_r8*vt+0.5_r8*vs)
         vtab(1,2,i)=1.0_r8/8.0_r8*(-vr-1.5_r8*vt+0.5_r8*vs)
         vtab(1,3,i)=1.0_r8/8.0_r8*(-vr+0.5_r8*vt-1.5_r8*vs)
         vtab(1,4,i)=1.0_r8/8.0_r8*(-vr-0.5_r8*vt-0.5_r8*vs)
         vtab(1,5,i)=0.0_r8
         vtab(1,6,i)=0.0_r8
         vtabls(1,1,i)=0.0_r8
         vtabls(1,2,i)=0.0_r8
         vtab(1,1:6,i)=vfact(1:6)*vtab(1,:,i)
      enddo
   elseif (abs(lpot).gt.100) then
      do i=0,ntab
         rr=i*dr
         call cheft_pot(abs(lpot),rr,vv)
         vtab(1,:,i)=vfact(1:6:1)*vv(1:6:1)
         if (abs(vv(7)).lt.verytiny) vv(7)=0.0_r8
         vtabls(1,:,i)=vfact(7:8:1)*vv(7:8:1)
!        if (myrank().eq.0) then
!           write(79,'(9e25.10)') rr,(vtab(1,j,i),j=1,6),(vtabls(1,j,i),j=1,2)
!        endif
      enddo
   endif 
     
   if (iem) then
      do i=0,ntab
         rr=i*dr
		 call empot(2,rr,vem) ! 3 means v6' here.
		 vemtab(i)=vem(1)	 
	  enddo	     
   endif
   
   return
   end subroutine v6potinit

   subroutine getvtab(vout) ! call v6potinit first
   use mympi
   integer(kind=i4) :: i
   real(kind=r8) :: vout(6,0:ntab)
   
   if (myrank().eq.0) then
   open(unit=9,form='formatted',file='vtable.dat')
   rewind 9
    do i=0,ntab
	  write(9,'(6(e15.7,2x))') vtab(1,:,i)
	  vout(:,i)=vtab(1,:,i)      
    enddo
   close(9) 
   else
    do i=0,ntab
	  vout(:,i)=vtab(1,:,i)      
    enddo          
   endif
 
   return
   end subroutine getvtab	
   
   subroutine getvtabschmidt(vout) ! call v6potinit first, Schmidt order. use this to linearize e^(v6).
   use mympi
   integer(kind=i4) :: i
   real(kind=r8) :: vout1(6,0:ntab),vout(6,0:ntab)
   if (myrank().eq.0) then
   open(unit=9,form='formatted',file='vtableSchmidt.dat')
   rewind 9
    do i=0,ntab
        vout1(1,i)=vtab(1,1,i)
        vout1(2,i)=vtab(1,3,i)
        vout1(3,i)=vtab(1,5,i)
        vout1(4,i)=vtab(1,2,i)
        vout1(5,i)=vtab(1,4,i)
        vout1(6,i)=vtab(1,6,i)    
	  write(9,'(6(e15.7,2x))') vout1(:,i)
    enddo
   close(9)  
   vout=vout1
   else
    do i=0,ntab
        vout1(1,i)=vtab(1,1,i)
        vout1(2,i)=vtab(1,3,i)
        vout1(3,i)=vtab(1,5,i)
        vout1(4,i)=vtab(1,2,i)
        vout1(5,i)=vtab(1,4,i)
        vout1(6,i)=vtab(1,6,i)    
    enddo 
    vout=vout1       
   endif

   return
   end subroutine getvtabschmidt	   
   
   
   subroutine getvemtab(vout) ! call v6potinit first
   use mympi
   integer(kind=i4) :: i
   real(kind=r8) :: vout(0:ntab)
   if (myrank().eq.0) then
   open(unit=9,form='formatted',file='vemtable.dat')
   rewind 9
    do i=0,ntab
	  write(9,'(6(e15.7,2x))') vemtab(i)
	  vout(i)=vemtab(i)      
    enddo
   close(9) 
   else
	  vout(0:ntab)=vemtab(0:ntab)      
   endif
   
   
   return
   end subroutine getvemtab   
   
   
   
   subroutine hspot(x,vc,v2,v3,v4,v5,v6,doall,itab)
   real(kind=r8), dimension (3,npart) :: x
   real(kind=r8), dimension(3,npart,3,npart) :: v5,v6
   real(kind=r8), dimension(npart,npart) :: v2,v3,v4
   real(kind=r8) :: vc
   real(kind=r8) :: dx0(3),dx(3),vv(6),vs(3,3),vst(3,3),r,c1,c2,c3,c4,dr
   integer(kind=i4) :: i,j,index,ix,iy,iz,itab
   logical :: doall
   real(kind=r8) :: gls(3,npart,npart),glsum(3,npart)
   v2=0.0_r8
   v3=0.0_r8
   v4=0.0_r8
   v5=0.0_r8
   v6=0.0_r8
   vc=0.0_r8
   gls=0.0_r8
   glsum=0.0_r8
   do j=2,npart
      do i=1,j-1
         dx0(:)=x(:,i)-x(:,j)
         dx0=dx0-el*nint(dx0*eli)
         vs=0.0_r8
         vst=0.0_r8
         do ix=-ibox,ibox
            do iy=-ibox,ibox
               do iz=-ibox,ibox
                  dx(1)=dx0(1)+el*ix
                  dx(2)=dx0(2)+el*iy
                  dx(3)=dx0(3)+el*iz
                  r=sqrt(dot_product(dx,dx))
                  dx=dx/r
                  if (r.lt.range) then
                     dr=scale*r
                     index=dr
                     index=max(1,min(index,ntab-2))
                     dr=dr-index
                     c1=-dr*(dr-1.0_r8)*(dr-2.0_r8)/6.0_r8
                     c2=(dr+1.0_r8)*(dr-1.0_r8)*(dr-2.0_r8)/2.0_r8
                     c3=-(dr+1.0_r8)*dr*(dr-2.0_r8)/2.0_r8
                     c4=(dr+1.0_r8)*dr*(dr-1.0_r8)/6.0_r8
                     vv(:)=c1*vtab(itab,:,index-1)+c2*vtab(itab,:,index) &
                        +c3*vtab(itab,:,index+1)+c4*vtab(itab,:,index+2)
                     vc=vc+vv(1)
                     if (doall) then
                        vs(1,1)=vs(1,1)+vv(5)*(3.0_r8*dx(1)*dx(1)-1.0_r8)
                        vs(1,2)=vs(1,2)+vv(5)*3.0_r8*dx(1)*dx(2)
                        vs(1,3)=vs(1,3)+vv(5)*3.0_r8*dx(1)*dx(3)
                        vs(2,2)=vs(2,2)+vv(5)*(3.0_r8*dx(2)*dx(2)-1.0_r8)
                        vs(2,3)=vs(2,3)+vv(5)*3.0_r8*dx(2)*dx(3)
                        vs(3,3)=vs(3,3)+vv(5)*(3.0_r8*dx(3)*dx(3)-1.0_r8)
                        vst(1,1)=vst(1,1)+vv(6)*(3.0_r8*dx(1)*dx(1)-1.0_r8)
                        vst(1,2)=vst(1,2)+vv(6)*3.0_r8*dx(1)*dx(2)
                        vst(1,3)=vst(1,3)+vv(6)*3.0_r8*dx(1)*dx(3)
                        vst(2,2)=vst(2,2)+vv(6)*(3.0_r8*dx(2)*dx(2)-1.0_r8)
                        vst(2,3)=vst(2,3)+vv(6)*3.0_r8*dx(2)*dx(3)
                        vst(3,3)=vst(3,3)+vv(6)*(3.0_r8*dx(3)*dx(3)-1.0_r8)
                        v2(i,j)=v2(i,j)+vv(2)
                        v3(i,j)=v3(i,j)+vv(3)
                        v4(i,j)=v4(i,j)+vv(4)
                     endif
                  endif
               enddo
            enddo
         enddo
         if (doall) then
            vs(2,1)=vs(1,2)
            vs(3,1)=vs(1,3)
            vs(3,2)=vs(2,3)
            vst(2,1)=vst(1,2)
            vst(3,1)=vst(1,3)
            vst(3,2)=vst(2,3)
            v2(j,i)=v2(i,j)
            v3(j,i)=v3(i,j)
            v4(j,i)=v4(i,j)
            v5(:,i,:,j)=vs(:,:)
            v5(:,j,:,i)=vs(:,:)
            v6(:,i,:,j)=vst(:,:)
            v6(:,j,:,i)=vst(:,:)
         endif
      enddo
   enddo
   end subroutine hspot

   subroutine vls(r,vvls,vvlst,itab)
   real(kind=r8) :: r,vvls,vvlst,dr,c1,c2,c3,c4
   integer(kind=i4) :: index,itab
   if (r.lt.range) then
      dr=scale*r
      index=dr
      index=max(1,min(index,ntab-2))
      dr=dr-index
      c1=-dr*(dr-1.0_r8)*(dr-2.0_r8)/6.0_r8
      c2=(dr+1.0_r8)*(dr-1.0_r8)*(dr-2.0_r8)/2.0_r8
      c3=-(dr+1.0_r8)*dr*(dr-2.0_r8)/2.0_r8
      c4=(dr+1.0_r8)*dr*(dr-1.0_r8)/6.0_r8
      vvls=c1*vtabls(itab,1,index-1)+c2*vtabls(itab,1,index) &
         +c3*vtabls(itab,1,index+1)+c4*vtabls(itab,1,index+2)
      vvlst=c1*vtabls(itab,2,index-1)+c2*vtabls(itab,2,index) &
         +c3*vtabls(itab,2,index+1)+c4*vtabls(itab,2,index+2)
   else
      vvls=0.0_r8
      vvlst=0.0_r8
   endif
   end subroutine vls

   subroutine coulomb(r,vc)
   real(kind=r8) :: r,vc,rr,alpha,hc,b,br,fcoul
   rr=r
   alpha=1.0_r8/137.03599_r8
   hc=197.327053_r8
   b=4.27_r8
   br=b*rr
   fcoul=1.0_r8-(1.0_r8+11.0_r8*br/16.0_r8+3.0_r8*br**2/16.0_r8+br**3/48.0_r8)*exp(-br)
   vc=alpha*hc*fcoul/rr
   end subroutine coulomb

   subroutine getvcsb(x,amat)
   real(kind=r8) :: x(:,:)
   real(kind=r8) :: amat(:,:)
   real(kind=r8) :: dx(3),r
   integer(kind=i4) :: i,j
   amat=0.0_r8
   do i=1,npart-1
      do j=i+1,npart
         dx(:)=x(:,i)-x(:,j)
         r=sqrt(sum(dx**2))
         call coulomb(r,amat(i,j))
         amat(j,i)=amat(i,j)
      enddo
   enddo
   end subroutine getvcsb
end module v6pot
