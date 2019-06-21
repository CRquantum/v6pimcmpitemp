module matrixmod
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   complex(kind=r8), private, parameter :: cone = (1.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: czero = (0.0_r8,0.0_r8)

interface matinv
   module procedure rmatinv,cmatinv
end interface

contains
   subroutine rsqsvd(a,u,sigma,vt,n) 
! 
!  svd of a real square matrix a. a=u*sigma*vt.
!
   integer(kind=i4), intent(in) :: n
   real(kind=r8), intent(in) :: a(n,n)
   real(kind=r8), intent(out) :: sigma(n),vt(n,n),u(n,n)
   integer(kind=i4) :: i
   vt=matmul(transpose(a),a)
   call eigenrs(vt,sigma,n)
   u=matmul(a,vt)
   sigma=sqrt(abs(sigma)) !abs protects against roundoff
   do i=1,n
      if (sigma(i).gt.0.0_r8) u(:,i)=u(:,i)/sigma(i)
   enddo
   vt=transpose(vt)
   end subroutine rsqsvd

   subroutine eigenrs(eigvec,eigval,n)
!
! compute eigenvalues and eigenvectors of a real symmetric matrix
! input matrix in eigvec (only lower triangle used), output eigenvectors
! in eigvec, output eigenvalues in eigval
!
! This stupid version are updated eispack tred2 and tql2 routines.
!
   integer(kind=i4) :: n,i,j,k,ierr,m,l1,l2,mml,ii,l,msave
   real(kind=r8) :: eigvec(n,n),eigval(n),e(n),b,el1,s,s2
   real(kind=r8) :: f,g,h,hh,scale,c3,dl1,c,c2,p,r
   if (n.lt.1) return
!
! tred2
!
   do i=n,2,-1
      scale=sum(abs(eigvec(i,1:i-1)))
      if (i.eq.2.or.scale.eq.0.0_r8) then
         e(i)=eigvec(i,i-1)
         eigval(i)=0.0_r8
         cycle
      endif
      eigvec(i,1:i-1)=eigvec(i,1:i-1)/scale
      h=sum(eigvec(i,1:i-1)**2)
      f=eigvec(i,i-1)
      g=-sign(sqrt(h),f)
      e(i)=scale*g
      h=h-f*g
      eigvec(i,i-1)=f-g
      f = 0.0_r8
      do j=1,i-1
         eigvec(j,i)=eigvec(i,j)/h
         g=sum(eigvec(j,1:j)*eigvec(i,1:j))
         do k=j+1,i-1
            g=g+eigvec(k,j)*eigvec(i,k)
         enddo
         e(j)=g/h
         f=f+e(j)*eigvec(i,j)
      enddo
      hh=f/(h+h)
      do j=1,i-1
         f=eigvec(i,j)
         g=e(j)-hh*f
         e(j)=g
         eigvec(j,1:j)=eigvec(j,1:j)-f*e(1:j)-g*eigvec(i,1:j)
      enddo
      eigval(i)=h
   enddo
   eigval(1)=0.0_r8
   e(1)=0.0_r8
   do i=1,n
      if (eigval(i).ne.0.0_r8) then
         do j=1,i-1
            g=sum(eigvec(i,1:i-1)*eigvec(1:i-1,j))
            eigvec(1:i-1,j)=eigvec(1:i-1,j)-g*eigvec(1:i-1,i)
         enddo
      endif
      eigval(i)=eigvec(i,i)
      eigvec(i,i)=1.0_r8
      eigvec(i,1:i-1)=0.0_r8
      eigvec(1:i-1,i)=0.0_r8
   enddo
!
! tql2
!
   if (n.eq.1) return
   do i=2,n
      e(i-1)=e(i)
   enddo
   f=0.0_r8
   b=0.0_r8
   e(n)=0.0_r8
   do l=1,n
      h=abs(eigval(l))+abs(e(l))
      if (b.lt.h) b=h
      do m=l,n
         msave=m
         if (b+abs(e(m)).eq.b) exit
      enddo
      m=msave
      if (m.eq.l) then
         eigval(l)=eigval(l)+f
         cycle
      endif
      do j=1,30
         if (j.eq.30) then
            write (6,'('' error in eigenrs '')')
            stop
         endif
         l1=l+1
         l2=l1+1
         g=eigval(l)
         p=(eigval(l1)-g)/(2.0_r8*e(l))
         r=pythag(p,1.0_r8)
         eigval(l)=e(l)/(p+sign(r,p))
         eigval(l1)=e(l)*(p+sign(r,p))
         dl1=eigval(l1)
         h=g-eigval(l)
         do i=l2,n
            eigval(i)=eigval(i)-h
         enddo
         f=f+h
         p=eigval(m)
         c=1.0_r8
         c2=c
         el1=e(l1)
         s=0.0_r8
         mml=m-l
         do ii=1,mml
            c3=c2
            c2=c
            s2=s
            i=m-ii
            g=c*e(i)
            h=c*p
            if (abs(p).lt.abs(e(i))) then
               c=p/e(i)
               r=sqrt(c*c+1.0_r8)
               e(i+1)=s*e(i)*r
               s=1.0_r8/r
               c=c*s
            else
               c=e(i)/p
               r=sqrt(c*c+1.0_r8)
               e(i+1)=s*p*r
               s= c/r
               c=1.0_r8/r
            endif
            p=c*eigval(i)-s*g
            eigval(i+1)=h+s*(c*g+s*eigval(i))
            do k=1,n
               h=eigvec(k,i+1)
               eigvec(k,i+1)=s*eigvec(k,i)+c*h
               eigvec(k,i)=c*eigvec(k,i)-s*h
            enddo
         enddo
         p=-s*s2*c3*el1*e(l)/dl1
         e(l)=s*p
         eigval(l)=c*p
         if (b + abs(e(l)).le. b) exit
      enddo
      eigval(l)=eigval(l)+f
   enddo
   do ii=2,n
      i=ii-1
      k=i
      p=eigval(i)
      do j=ii,n
         if (eigval(j).ge.p) cycle
         k=j
         p=eigval(j)
      enddo
      if (k.ne.i) then
         eigval(k)=eigval(i)
         eigval(i)=p
         do j=1,n
            p=eigvec(j,i)
            eigvec(j,i)=eigvec(j,k)
            eigvec(j,k)=p
         enddo
      endif
   enddo
   return
   end subroutine eigenrs

   function pythag(a,b)
   real(kind=r8) ::  a,b,pythag,q,r,s,t
   pythag=max(abs(a),abs(b))
   q=min(abs(a),abs(b))
   if (q.eq.0.0_r8) return
   t=0.0_r8
   do while (t.ne.4.0_r8)
      r=(q/pythag)**2
      t=4.0_r8+r
      s=r/t
      pythag=pythag+2.0_r8*pythag*s
      q=q*s
   enddo
   return
   end function pythag

   subroutine rmatinv(a,detl,is,n)
   integer(kind=i4) :: n,iwork(n),i,j,k,itemp,idiag,is
   real(kind=r8) :: a(n,n),detl,work(n),adiag,t
   real(kind=r8) :: aabs,adiagabs
   do i=1,n
   enddo
   do i=1,n
      iwork(i)=i
   enddo
   detl=0.0_r8
   is=1
   do i=1,n
      adiag=a(iwork(i),i)
      idiag=i
      adiagabs=abs(adiag)
      do k=i,n
         aabs=abs(a(iwork(k),i))
         if (aabs.gt.adiagabs) then
            adiag=a(iwork(k),i)
            adiagabs=aabs
            idiag=k
         endif
      enddo
      if (idiag.ne.i) then
         is=-is
         itemp=iwork(i)
         iwork(i)=iwork(idiag)
         iwork(idiag)=itemp
      endif
      detl=detl+log(abs(adiag))
      is=is*sign(1.0_r8,adiag)
      a(iwork(i),i)=1.0_r8
      a(iwork(i),:)=a(iwork(i),:)/adiag
      do j=1,n
         if (j.ne.iwork(i)) then
            t=-a(j,i)
            a(j,i)=0.0_r8
            a(j,:)=a(j,:)+t*a(iwork(i),:)
         endif
      enddo
   enddo
!
! unpivot equivalent to anew(i,iwork(j))=aold(iwork(i),j)
!
   do j=1,n
      work=a(:,j)
      do k=1,n
         a(k,j)=work(iwork(k))
      enddo
   enddo
   do i=1,n
      work=a(i,:)
      do j=1,n
         a(i,iwork(j))=work(j)
      enddo
   enddo
   return
   end subroutine rmatinv


   subroutine cmatinv(a,det,n)
   integer(kind=i4) :: n,iwork(n),i,j,k,itemp,idiag
   complex (kind=r8) :: a(n,n),det,work(n),adiag,t
   real(kind=r8) :: aabs,adiagabs
   do i=1,n
   enddo
   do i=1,n
      iwork(i)=i
   enddo
   det=cone
   do i=1,n
      adiag=a(iwork(i),i)
      idiag=i
      adiagabs=abs(adiag)
      do k=i,n
         aabs=abs(a(iwork(k),i))
         if (aabs.gt.adiagabs) then
            adiag=a(iwork(k),i)
            adiagabs=aabs
            idiag=k
         endif
      enddo
      if (idiag.ne.i) then
         det=-det
         itemp=iwork(i)
         iwork(i)=iwork(idiag)
         iwork(idiag)=itemp
      endif
      det=det*adiag
      a(iwork(i),i)=cone
      a(iwork(i),:)=a(iwork(i),:)/adiag
      do j=1,n
         if (j.ne.iwork(i)) then
            t=-a(j,i)
            a(j,i)=czero
            a(j,:)=a(j,:)+t*a(iwork(i),:)
         endif
      enddo
   enddo
!
! unpivot equivalent to anew(i,iwork(j))=aold(iwork(i),j)
!
   do j=1,n
      work=a(:,j)
      do k=1,n
         a(k,j)=work(iwork(k))
      enddo
   enddo
   do i=1,n
      work=a(i,:)
      do j=1,n
         a(i,iwork(j))=work(j)
      enddo
   enddo
   return
   end subroutine cmatinv

   function expmult(a,v,n)
!
! calculate expmult = exp(a)*v
! this stupid version calls eispack routines converted without thought
! to r8. Notice that everything gets copied around a lot.
!
   integer(kind=i4) :: n,ilo,ihi,info
   complex(kind=r8) :: a(n,n),v(n),val(n),vecr(n,n),vecl(n,n),det
   complex(kind=r8) :: expmult(n)
   real(kind=r8) :: scale(n),ar(n,n),ai(n,n),ortr(n),orti(n)
   real(kind=r8) :: valr(n),vali(n),vecrr(n,n),vecri(n,n)
   if (n.eq.1) then
      expmult(1)=exp(a(1,1))*v(1)
      return
   endif
   ar=real(a)
   ai=aimag(a)
   call cbal(n,n,ar,ai,ilo,ihi,scale)
   call corth(n,n,ilo,ihi,ar,ai,ortr,orti)
   call comqr2(n,n,ilo,ihi,ortr,orti,ar,ai,valr,vali,vecrr,vecri,info)
   if (info.ne.0) then
      write (6,'(1x,''error in comqr2'',i10)') info
      stop
   endif
   call cbabk2(n,n,ilo,ihi,scale,n,vecrr,vecri)
   vecr=cmplx(vecrr,vecri)
   val=cmplx(valr,vali)
   vecl=vecr
   call cmatinv(vecl,det,n)
   val=exp(val)
   expmult=matmul(vecr,val*matmul(vecl,v))
   end function expmult

      subroutine cbal(nm,n,ar,ai,low,igh,scale)
      integer(kind=i4) :: i,j,k,l,m,n,jj,nm,igh,low,iexc
      real(kind=r8) :: ar(nm,n),ai(nm,n),scale(n)
      real(kind=r8) :: c,f,g,r,s,b2,radix
      logical :: noconv
      radix = 16
      b2 = radix * radix
      k = 1
      l = n
      go to 100
   20 scale(m) = j
      if (j .eq. m) go to 50
      do 30 i = 1, l
         f = ar(i,j)
         ar(i,j) = ar(i,m)
         ar(i,m) = f
         f = ai(i,j)
         ai(i,j) = ai(i,m)
         ai(i,m) = f
   30 continue
      do 40 i = k, n
         f = ar(j,i)
         ar(j,i) = ar(m,i)
         ar(m,i) = f
         f = ai(j,i)
         ai(j,i) = ai(m,i)
         ai(m,i) = f
   40 continue
   50 if (iexc.eq.1) go to 80
      go to 130
   80 if (l .eq. 1) go to 280
      l = l - 1
  100 do 120 jj = 1, l
         j = l + 1 - jj
         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (ar(j,i) .ne. 0.0_r8 .or. ai(j,i) .ne. 0.0_r8) go to 120
  110    continue
         m = l
         iexc = 1
         go to 20
  120 continue
      go to 140
  130 k = k + 1
  140 do 170 j = k, l
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (ar(i,j) .ne. 0.0_r8 .or. ai(i,j) .ne. 0.0_r8) go to 170
  150    continue
         m = k
         iexc = 2
         go to 20
  170 continue
      do 180 i = k, l
      scale(i) = 1.0_r8
  180 continue
  190 noconv = .false.
      do 270 i = k, l
         c = 0.0_r8
         r = 0.0_r8
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + abs(ar(j,i)) + abs(ai(j,i))
            r = r + abs(ar(i,j)) + abs(ai(i,j))
  200    continue
         if (c .eq. 0.0_r8 .or. r .eq. 0.0_r8) go to 270
         g = r / radix
         f = 1.0_r8
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
  240    if ((c + r) / f .ge. 0.95_r8 * s) go to 270
         g = 1.0_r8 / f
         scale(i) = scale(i) * f
         noconv = .true.
         do 250 j = k, n
            ar(i,j) = ar(i,j) * g
            ai(i,j) = ai(i,j) * g
  250    continue
         do 260 j = 1, l
            ar(j,i) = ar(j,i) * f
            ai(j,i) = ai(j,i) * f
  260    continue
  270 continue
      if (noconv) go to 190
  280 low = k
      igh = l
      return
      end subroutine cbal

      subroutine corth(nm,n,low,igh,ar,ai,ortr,orti)
      integer(kind=i4) :: i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      real(kind=r8) :: ar(nm,n),ai(nm,n),ortr(igh),orti(igh)
      real(kind=r8) :: f,g,h,fi,fr,scale
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
      do 180 m = kp1, la
         h = 0.0_r8
         ortr(m) = 0.0_r8
         orti(m) = 0.0_r8
         scale = 0.0_r8
         do 90 i = m, igh
         scale = scale + abs(ar(i,m-1)) + abs(ai(i,m-1))
   90 continue
         if (scale .eq. 0.0_r8) go to 180
         mp = m + igh
         do 100 ii = m, igh
            i = mp - ii
            ortr(i) = ar(i,m-1) / scale
            orti(i) = ai(i,m-1) / scale
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
  100    continue
         g = sqrt(h)
         f = pythag(ortr(m),orti(m))
         if (f .eq. 0.0_r8) go to 103
         h = h + f * g
         g = g / f
         ortr(m) = (1.0_r8 + g) * ortr(m)
         orti(m) = (1.0_r8 + g) * orti(m)
         go to 105
  103    ortr(m) = g
         ar(m,m-1) = scale
  105    do 130 j = m, n
            fr = 0.0_r8
            fi = 0.0_r8
            do 110 ii = m, igh
               i = mp - ii
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
  110       continue
            fr = fr / h
            fi = fi / h
            do 120 i = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
  120       continue
  130    continue
         do 160 i = 1, igh
            fr = 0.0_r8
            fi = 0.0_r8
            do 140 jj = m, igh
               j = mp - jj
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
  140       continue
            fr = fr / h
            fi = fi / h
            do 150 j = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
  150       continue
  160    continue
         ortr(m) = scale * ortr(m)
         orti(m) = scale * orti(m)
         ar(m,m-1) = -g * ar(m,m-1)
         ai(m,m-1) = -g * ai(m,m-1)
  180 continue
  200 return
      end subroutine corth

      subroutine comqr2(nm,n,low,igh,ortr,orti,hr,hi,wr,wi,zr,zi,ierr)
      integer(kind=i4) :: i,j,k,l,m,n,en,ii,jj,ll,nm,nn,igh,ip1
      integer(kind=i4) :: itn,its,low,lp1,enm1,iend,ierr
      real(kind=r8) :: hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n)
      real(kind=r8) :: ortr(igh),orti(igh)
      real(kind=r8) :: si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,s1,s2
      ierr = 0
      do 100 i = 1, n
         do 101 j = 1, n
            zr(i,j) = 0.0_r8
            zi(i,j) = 0.0_r8
            if (i .eq. j) zr(i,j) = 1.0_r8
  101 continue
  100 continue
      iend = igh - low - 1
      if (iend.lt.0) then
         go to 180
      else if (iend.eq.0) then
         go to 150
      endif
      do 140 ii = 1, iend
         i = igh - ii
         if (ortr(i) .eq. 0.0_r8 .and. orti(i) .eq. 0.0_r8) go to 140
         if (hr(i,i-1) .eq. 0.0_r8 .and. hi(i,i-1) .eq. 0.0_r8) go to 140
         norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)
         ip1 = i + 1
         do 110 k = ip1, igh
            ortr(k) = hr(k,i-1)
            orti(k) = hi(k,i-1)
  110    continue
         do 130 j = i, igh
            sr = 0.0_r8
            si = 0.0_r8
            do 115 k = i, igh
               sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
               si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
  115       continue
            sr = sr / norm
            si = si / norm
            do 120 k = i, igh
               zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
               zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
  120       continue
  130    continue
  140 continue
  150 l = low + 1
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0_r8) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0_r8
         do 155 j = i, n
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
         do 160 j = 1, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
         do 165 j = low, igh
            si = yr * zi(j,i) + yi * zr(j,i)
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
            zi(j,i) = si
  165    continue
  170 continue
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
      en = igh
      tr = 0.0_r8
      ti = 0.0_r8
      itn = 30*n
  220 if (en .lt. low) go to 680
      its = 0
      enm1 = en - 1
  240 continue
      do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         s1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1)) &
                   + abs(hr(l,l)) +abs(hi(l,l))
         s2 = s1 + abs(hr(l,l-1))
         if (s2 .eq. s1) go to 300
  260 continue
  300 continue
      if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0_r8 .and. xi .eq. 0.0_r8) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0_r8
      yi = (hi(enm1,enm1) - si) / 2.0_r8
      call csroot(yr**2-yi**2+xr,2.0_r8*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0_r8) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
      si = 0.0_r8
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
      lp1 = l + 1
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0_r8
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0_r8
         hi(i,i-1) = sr / norm
         do 490 j = i, n
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
  500 continue
      si = hi(en,en)
      if (si .eq. 0.0_r8) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0_r8
      if (en .eq. n) go to 540
      ip1 = en + 1
      do 520 j = ip1, n
         yr = hr(en,j)
         yi = hi(en,j)
         hr(en,j) = sr * yr + si * yi
         hi(en,j) = sr * yi - si * yr
  520 continue
  540 continue
      do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
         do 580 i = 1, j
            yr = hr(i,j-1)
            yi = 0.0_r8
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
         do 590 i = low, igh
            yr = zr(i,j-1)
            yi = zi(i,j-1)
            zzr = zr(i,j)
            zzi = zi(i,j)
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  590    continue
  600 continue
      if (si .eq. 0.0_r8) go to 240
      do 630 i = 1, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
      do 640 i = low, igh
         yr = zr(i,en)
         yi = zi(i,en)
         zr(i,en) = sr * yr - si * yi
         zi(i,en) = sr * yi + si * yr
  640 continue
      go to 240
  660 hr(en,en) = hr(en,en) + tr
      wr(en) = hr(en,en)
      hi(en,en) = hi(en,en) + ti
      wi(en) = hi(en,en)
      en = enm1
      go to 220
  680 norm = 0.0_r8
      do 720 i = 1, n
         do 721 j = i, n
            norm = norm + abs(hr(i,j)) + abs(hi(i,j))
  721 continue
  720 continue
      if (n .eq. 1 .or. norm .eq. 0.0_r8) go to 1001
      do 800 nn = 2, n
         en = n + 2 - nn
         xr = wr(en)
         xi = wi(en)
         enm1 = en - 1
         do 780 ii = 1, enm1
            i = en - ii
            zzr = hr(i,en)
            zzi = hi(i,en)
            if (i .eq. enm1) go to 760
            ip1 = i + 1
            do 740 j = ip1, enm1
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
  740       continue
  760       yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.0_r8 .or. yi .ne. 0.0_r8) go to 775
            yr = norm
  770       yr = 0.5_r8*yr
            if (norm + yr .gt. norm) go to 770
            yr = 2.0_r8*yr
  775       call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
  780    continue
  800 continue
      enm1 = n - 1
      do  840 i = 1, enm1
         if (i .ge. low .and. i .le. igh) go to 840
         ip1 = i + 1
         do 820 j = ip1, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
  820    continue
  840 continue
      do 880 jj = low, enm1
         j = n + low - jj
         m = min0(j-1,igh)
         do 881 i = low, igh
            zzr = zr(i,j)
            zzi = zi(i,j)
            do 860 k = low, m
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
  860       continue
            zr(i,j) = zzr
            zi(i,j) = zzi
  881 continue
  880 continue
      go to 1001
 1000 ierr = en
 1001 return
      end subroutine comqr2

      subroutine cbabk2(nm,n,low,igh,scale,m,zr,zi)
      integer(kind=i4) :: i,j,k,m,n,ii,nm,igh,low
      real(kind=r8) :: scale(n),zr(nm,m),zi(nm,m)
      real(kind=r8) :: s
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
      do 110 i = low, igh
         s = scale(i)
         do 100 j = 1, m
            zr(i,j) = zr(i,j) * s
            zi(i,j) = zi(i,j) * s
  100    continue
  110 continue
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
         do 130 j = 1, m
            s = zr(i,j)
            zr(i,j) = zr(k,j)
            zr(k,j) = s
            s = zi(i,j)
            zi(i,j) = zi(k,j)
            zi(k,j) = s
  130    continue
  140 continue
  200 return
      end subroutine cbabk2

      subroutine csroot(xr,xi,yr,yi)
      real(kind=r8) xr,xi,yr,yi,s,tr,ti
      tr = xr
      ti = xi
      s = sqrt(0.5_r8*(pythag(tr,ti) + abs(tr)))
      if (tr .ge. 0.0_r8) yr = s
      if (ti .lt. 0.0_r8) s = -s
      if (tr .le. 0.0_r8) yi = s
      if (tr .lt. 0.0_r8) yr = 0.5_r8*(ti/yi)
      if (tr .gt. 0.0_r8) yi = 0.5_r8*(ti/yr)
      return
      end subroutine csroot

      subroutine cdiv(ar,ai,br,bi,cr,ci)
      real(kind=r8) :: ar,ai,br,bi,cr,ci
      real(kind=r8) :: s,ars,ais,brs,bis
      s = abs(br) + abs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end subroutine cdiv

end module matrixmod
