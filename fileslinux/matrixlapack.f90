module matrixmod
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   complex(kind=r8), private, parameter :: cone = (1.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: czero = (0.0_r8,0.0_r8)

interface matinv
   module procedure rmatinv,cmatinv
end interface

interface eigen ! eigen vectors has been normalized to 1.
   module procedure eigenh,eigenrs
   module procedure eigengc,eigengr
end interface eigen

contains  
    subroutine eigengr(a,eigvecl,eigvecr,eigval,n)
!*  The routine computes for an n-by-n real general/nonsymmetric matrix A
    integer(kind=i4) :: n,info
    real(kind=r8) :: rwork(50*n)
    real(kind=r8) :: a(n,n)
    complex(kind=r8) :: atmp(n,n),eigval(n),eigvecl(n,n),eigvecr(n,n),work(50*n)
    integer(kind=i4) :: i,lwork
    integer(kind=i4), parameter :: lwmax=1000
    if (n.lt.1) return
    lwork=size(work)
    atmp=cmplx(a)
    call zgeev('Vectors','Vectors',n,atmp,n,eigval,eigvecl,n,eigvecr,n,work,lwork,rwork,info)   
    if (info.ne.0) then
        write (6,'(1x,''error in eigengr -> zgeev'',i10)') info
        open(unit=9,form='formatted',file='zgeev.error',position='rewind')
        write (9,'(i10)') info
        write (9,'(12g15.7)') a
        close(9)
        call abort()
    endif
    return
    end subroutine eigengr
    
    subroutine eigengc(a,eigvecl,eigvecr,eigval,n)
!*  The routine computes for an n-by-n real/complex general/nonsymmetric matrix A, the
!*  eigenvalues and, optionally, the left and/or right eigenvectors. The right
!*  eigenvector v(j) of A satisfies
!*  A*v(j)= lambda(j)*v(j)
!*  where lambda(j) is its eigenvalue. The left eigenvector u(j) of A satisfies
!*  u(j)H*A = lambda(j)*u(j)H
!*  where u(j)H denotes the conjugate transpose of u(j). The computed
!*  eigenvectors are normalized to have Euclidean norm equal to 1 and
!*  largest component real.
   integer(kind=i4) :: n,info
   real(kind=r8) :: rwork(50*n)
   complex(kind=r8) :: a(n,n),atmp(n,n),eigval(n),eigvecl(n,n),eigvecr(n,n),work(50*n)
   integer(kind=i4) :: i,lwork
   integer(kind=i4), parameter :: lwmax=1000
   if (n.lt.1) return
   lwork=size(work)
   atmp=a
   call zgeev('Vectors','Vectors',n,atmp,n,eigval,eigvecl,n,eigvecr,n,work,lwork,rwork,info)   
   if (info.ne.0) then
      write (6,'(1x,''error in eigengc -> zgeev'',i10)') info
      open(unit=9,form='formatted',file='zgeev.error',position='rewind')
      write (9,'(i10)') info
      write (9,'(12g15.7)') a
      close(9)
      call abort()
   endif
   return
   end subroutine eigengc    
    
   subroutine eigenh(a,eigvec,eigval,n)
!
! compute eigenvalues and eigenvectors of a Hermitian matrix A
! input matrix in a (only lower triangle used), output eigenvectors
! in eigvec, output eigenvalues in eigval
!
   integer(kind=i4) :: n,info
!
! lwork >= max(1, 2n-1).
!
   real(kind=r8) :: eigval(n),rwork(50*n)
   complex(kind=r8) :: a(n,n),eigvec(n,n),work(50*n)
   integer(kind=i4) :: i,lwork
   if (n.lt.1) return
   lwork=size(work)
   eigvec=a
   call zheev('Vectors', 'Lower', n, eigvec, n, eigval, work, lwork, rwork, info)
   if (info.ne.0) then
      write (6,'(1x,''error in eigenh -> zheev'',i10)') info
      open(unit=9,form='formatted',file='zheev.error',position='rewind')
      write (9,'(i10)') info
      write (9,'(12g15.7)') a
      close(9)
      call abort()
   endif
   !do i=1,n
   !   if (eigvec(1,i).lt.0.0_r8) eigvec(:,i)=-eigvec(:,i)
   !enddo
   return
   end subroutine eigenh    
    
   subroutine rsqsvd(a,u,sigma,vt,n) 
! 
!  svd of a real square matrix a. a=u*sigma*vt.
!
   integer(kind=i4), intent(in) :: n
   real(kind=r8), intent(in) :: a(n,n)
   real(kind=r8), intent(out) :: sigma(n),vt(n,n),u(n,n)
   real(kind=r8) :: work(50*n),atmp(n,n)
   integer(kind=i4) :: info,lwork
   lwork=size(work) !must be >= 5n
   if (n.lt.1) return
   atmp=a !dgesvd destroys a
   call dgesvd('All','All',n,n,atmp,n,sigma,u,n,vt,n,work,lwork,info)
   if (info.ne.0) then
      write (6,'(''error in rsqsvd -> dgesvd '',i10)') info
      open(unit=9,form='formatted',file='dgesvd.error',position='rewind')
      write (9,'(i10)') info
      write (9,'(12g15.7)') a
      close(9)
      call abort
   endif
   return
   end subroutine rsqsvd
    
   subroutine eigenrs(a,eigvec,eigval,n)
!
! compute eigenvalues and eigenvectors of a real symmetric matrix
! input matrix in eigvec (only lower triangle used), output eigenvectors
! in eigvec, output eigenvalues in eigval
!
   integer(kind=i4) :: n,info
!
! lwork >= (nb+2)*n so assume block size 128 is more than enough
!
   real(kind=r8) :: a(n,n),eigvec(n,n),eigval(n),work(130*n)
   integer(kind=i4) :: i
   if (n.lt.1) return
   eigvec=a
   call dsyev('v','l',n,eigvec,n,eigval,work,130*n,info)
   if (info.ne.0) then
      write (6,'(1x,''error in eigenrs -> dsyev'',i10)') info
      open(unit=9,form='formatted',file='dsyev.error',position='rewind')
      write (9,'(i10)') info
      write (9,'(12g15.7)') a
      close(9)
      call abort()
   endif
   do i=1,n
      if (eigvec(1,i).lt.0.0_r8) eigvec(:,i)=-eigvec(:,i)
   enddo
   return
   end subroutine eigenrs

   subroutine cmatinv(a,det,n)
   integer(kind=i4) :: n,ipiv(n),info,i
   complex (kind=r8) :: a(n,n),det,cwork(n,n)
!
! lapack routine for lu factorization
!
   call zgetrf(n,n,a,n,ipiv,info)
   if (info.lt.0) then
      write (6,'(1x,''error in zgetrf'',i10)') info
      open(unit=9,form='formatted',file='zgetrf.error',position='rewind')
      write (9,'(i10)') info
      write (9,'(12g15.7)') a
      close(9)
      call abort()
   endif
   if (info.gt.0) then !determinant is zero -- a is undefined -- set to zero
      det=czero
      a=czero
      return
   endif
!
! calculate determinant
!
   det=cone
   do i=1,n
      det=det*a(i,i)
      if (ipiv(i).ne.i) det=-det
   enddo
!
! lapack routine to calculate inverse from factorization
!
   call zgetri(n,a,n,ipiv,cwork,n*n,info)
   if (info.ne.0) then
      write (6,'(1x,''error in zgetri'',i10)') info
      open(unit=9,form='formatted',file='zgetri.error',position='rewind')
      write (9,'(i10)') info
      write (9,'(12g15.7)') a
      close(9)
      call abort()
   endif
   end subroutine cmatinv

   subroutine rmatinv(a,detl,is,n)
! calculate inverse, log(abs(det)), and sign.
   integer(kind=i4) :: n,is,ipiv(n),i,info
   real(kind=r8) :: a(n,n),detl,work(n*n)
   real(kind=r8), parameter :: minushuge=-1e30_r8
   call dgetrf(n,n,a,n,ipiv,info)
   if (info.ne.0) then
      write (6,'(1x,''error in dgetrf'',i10)') info
      open(unit=9,form='formatted',file='dgetrf.error',position='rewind')
      write (9,'(i10)') info
      write (9,'(12g15.7)') a
      close(9)
      stop
   endif
   is=1
   detl=0.0_r8
   do i=1,n
      detl=detl+log(abs(a(i,i)))
      is=is*sign(1.0_r8,a(i,i))
      if (ipiv(i).ne.i) is=-is
   enddo
   call dgetri(n,a,n,ipiv,work,n*n,info)
   if (info.ne.0) then
      write (6,'(1x,''error in dgetri'',i10)') info
      open(unit=9,form='formatted',file='dgetri.error',position='rewind')
      write (9,'(i10)') info
      write (9,'(12g15.7)') a
      close(9)
      stop
   endif
   end subroutine rmatinv

   function expmult(a,vec,n)
!
! routine to perform vec=exp(a)*vec for a complex matrix a and vector vec
!
   integer(kind=i4) :: n,info
   complex(kind=r8) :: a(n,n),vec(n),expmult(n)
   complex(kind=r8) :: det,val(n),vecr(n,n),vecl(n,n),cwork(n,n)
   real(kind=r8) :: work(2*n)
   call zgeev('n','v',n,a,n,val,vecl,n,vecr,n,cwork,n*n,work,info)
   if (info.ne.0) then
      write (6,'(1x,''problem in expmult'',i5)') info
      open(unit=9,form='formatted',file='expmult.error',position='rewind')
      write (9,'(i10)') info
      write (9,'(12g15.7)') a
      close(9)
      stop
   endif
   vecl=vecr
   call matinv(vecl,det,n)
   val=exp(val)
   expmult=matmul(vecr,val*matmul(vecl,vec))
   end function expmult

!=======================================================================================
   
    subroutine checkeigenh
      INTEGER          N
      PARAMETER        ( N = 4 )
      INTEGER          LDA
      PARAMETER        ( LDA = N )   
      real(kind=r8) :: eigval(n)
      complex(kind=r8) :: a(n,n),eigvec(n,n)  
       DATA             a/ &
      ( 9.14, 0.00),(-4.37, 9.22),(-1.98, 1.72),(-8.96, 9.50), &
      ( 0.00, 0.00),(-3.35, 0.00),( 2.25, 9.51),( 2.57,-2.40), &
      ( 0.00, 0.00),( 0.00, 0.00),(-4.82, 0.00),(-3.24,-2.04), &
      ( 0.00, 0.00),( 0.00, 0.00),( 0.00, 0.00),( 8.44, 0.00) &
                       /     
      call eigen(a,eigvec,eigval,n)  
!     Print eigenvalues.
      CALL PRINT_RMATRIX( 'Eigenvalues', 1, N, eigval, 1 )
!     Print eigenvectors.
      CALL PRINT_MATRIX( 'Eigenvectors (stored columnwise)', N, N, eigvec, LDA )    
      return
      end
    
      subroutine checkeigengc
      INTEGER          N
      PARAMETER        ( N = 4 )
      INTEGER          LDA
      PARAMETER        ( LDA = N )   
      complex(kind=r8) :: a(n,n),eigvec(n,n),eigval(n),eigvecl(n,n),eigvecr(n,n)
      DATA             a/ &
      (-3.84, 2.25),(-0.66, 0.83),(-3.99,-4.73),( 7.74, 4.18), &
      (-8.94,-4.75),(-4.40,-3.82),(-5.88,-6.60),( 3.66,-7.53), &
      ( 8.95,-6.53),(-3.50,-4.26),(-3.36,-0.40),( 2.58, 3.60), &
      (-9.87, 4.82),(-3.15, 7.36),(-0.75, 5.23),( 4.59, 5.41) /   
    
      call eigen(a,eigvecl,eigvecr,eigval,n)  
!     Print eigenvalues.
      CALL PRINT_MATRIX( 'Eigenvalues', 1, n, eigval, 1 )
!     Print left eigenvectors.
      CALL PRINT_MATRIX( 'Left eigenvectors', n, n, eigvecl, LDA )
!     Print right eigenvectors.
      CALL PRINT_MATRIX( 'Right eigenvectors', n, n, eigvecr, LDA )  
      return
      end    
      
      subroutine checkeigengr
      INTEGER          N
      PARAMETER        ( N = 5 )
      INTEGER          LDA
      PARAMETER        ( LDA = N )  
      real(kind=r8) :: a(n,n)    
      complex(kind=r8) :: eigvec(n,n),eigval(n),eigvecl(n,n),eigvecr(n,n)
      DATA             A/ &
      -1.01, 3.98, 3.30, 4.43, 7.31, &
       0.86, 0.53, 8.26, 4.96,-6.43, &
      -4.60,-7.04,-3.89,-7.66,-6.16, &
       3.31, 5.29, 8.20,-7.33, 2.47, &
      -4.81, 3.55,-1.51, 6.18, 5.58 / 
    
      call eigen(a,eigvecl,eigvecr,eigval,n)  
!     Print eigenvalues.
      CALL PRINT_MATRIX( 'Eigenvalues', 1, n, eigval, 1 )
!     Print left eigenvectors.
      CALL PRINT_MATRIX( 'Left eigenvectors', n, n, eigvecl, LDA )
!     Print right eigenvectors.
      CALL PRINT_MATRIX( 'Right eigenvectors', n, n, eigvecr, LDA )  
      return
      end          
!  =============================================================================
!  from Intel Math Kernel Library LAPACK Examples
!
!     Auxiliary routine: printing a matrix.
!
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      COMPLEX*16       A( LDA, * )
!
      INTEGER          I, J
!
      !WRITE(*,*)
      !WRITE(*,*) DESC
      open(unit=19,form='formatted',file=DESC)       
      DO I = 1, M
         !WRITE(*,9998) ( A( I, J ), J = 1, N )
         WRITE(19,9998) ( A( I, J ), J = 1, N )
      END DO
!
9998 FORMAT( 11(:,1X,'(',g15.7,',',g15.7,')') )
!     Print eigenvectors.
!      CALL PRINT_MATRIX( 'Eigenvectors (stored columnwise)', N, N, A, LDA ) ! LDA is usually n.
      close(19)
      
      
      !if (LDA/=1) then      
      !  do j=1,n
      !  write(6,*) 'j, normalization factor:', j,sum(abs(A(:,j))**2)
      !  enddo
      !endif
      RETURN
      END
!
!     Auxiliary routine: printing a real matrix.
!
      SUBROUTINE PRINT_RMATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      DOUBLE PRECISION A( LDA, * )
!
      INTEGER          I, J
!
      !WRITE(*,*)
      !WRITE(*,*) DESC
      
      open(unit=19,form='formatted',file=DESC)         
      DO I = 1, M
         !WRITE(*,9998) ( A( I, J ), J = 1, N )
         WRITE(19,9998) ( A( I, J ), J = 1, N )
         
      END DO
!
9998 FORMAT( 11(:,1X,g15.7) )
!     Print eigenvalues.
!      CALL PRINT_RMATRIX( 'Eigenvalues', 1, N, W, 1 )
     
      close(19)
      RETURN
      END   
      
      SUBROUTINE PRINT_CMATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      COMPLEX*16 A( LDA, * )
!
      INTEGER          I, J
!
      !WRITE(*,*)
      !WRITE(*,*) DESC
      
      open(unit=19,form='formatted',file=DESC)         
      DO I = 1, M
         !WRITE(*,9998) ( A( I, J ), J = 1, N )
         WRITE(19,9998) ( A( I, J ), J = 1, N )
         
      END DO
!
9998 FORMAT( 11(:,1X,'(',g15.7,',',g15.7,')') )
!     Print eigenvalues.
!      CALL PRINT_RMATRIX( 'Eigenvalues', 1, N, W, 1 )
     
      close(19)
      RETURN
      END       
   
end module matrixmod
