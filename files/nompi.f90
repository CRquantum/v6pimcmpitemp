!$Id: mympi.f90,v 1.2 2011/04/19 15:45:39 schmidt Exp $
module mympi
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: irank,iproc

interface bcast ! broadcast from process 0
   module procedure bcastirn
   module procedure bcasti1,bcasti1d,bcasti2d,bcasti3d,bcasti4d
   module procedure bcastr1,bcastr1d,bcastr2d,bcastr3d,bcastr4d
   module procedure bcastlogi
   module procedure bcastchar
end interface bcast

interface addall ! return sum to process 0
   module procedure addalli1,addalli1d
   module procedure addallr1,addallr1d,addallr2d
   module procedure addallc1,addallc1d
end interface addall

interface gather ! gather to process 0
   module procedure gatheri1,gatheri1d
   module procedure gatherr1,gatherr1d
end interface gather

interface scatter ! scatter from process 0 to all process, evenly scatter. inverse of gather
   module procedure scatteri1,scatteri1d
   module procedure scatterr1,scatterr1d
end interface scatter

interface send ! send to someone else
   module procedure sendi1,sendi1d
   module procedure sendr1,sendr1d
   module procedure sendc1,sendc1d
end interface send

interface recv ! recv from someone else
   module procedure recvi1,recvi1d
   module procedure recvr1,recvr1d
   module procedure recvc1,recvc1d
end interface recv

contains
   subroutine init0 ! call this before anything else
   integer :: ierror,isize,ir,ip
   integer(kind=i4) :: itest4
   integer(kind=i8) :: itest8
   real(kind=r8) :: rtest8
   irank=0
   iproc=1
   end subroutine init0

   subroutine done ! wrapper for finalize routine
   integer :: ierror
   stop 'stop, done!'
   end subroutine done

   subroutine bcastlogi(i)
   logical :: i
   integer :: ierror
   return
   end subroutine bcastlogi   
   
   subroutine bcastirn(i)
   integer(kind=i8) :: i
   integer :: ierror
   return
   end subroutine bcastirn

   subroutine bcasti1(i)
   integer(kind=i4) :: i
   integer :: ierror
   return
   end subroutine bcasti1

   subroutine bcasti1d(i)
   integer(kind=i4) :: i(:)
   integer :: ierror
   return
   end subroutine bcasti1d

   subroutine bcasti2d(i)
   integer(kind=i4) :: i(:,:)
   integer :: ierror
   return
   end subroutine bcasti2d

   subroutine bcasti3d(i)
   integer(kind=i4) :: i(:,:,:)
   integer :: ierror
   return
   end subroutine bcasti3d

   subroutine bcasti4d(i)
   integer(kind=i4) :: i(:,:,:,:)
   integer :: ierror
   return
   end subroutine bcasti4d

   subroutine bcastr1d(r)
   real(kind=r8) :: r(:)
   integer :: ierror
   return
   end subroutine bcastr1d

   subroutine bcastr2d(r)
   real(kind=r8) :: r(:,:)
   integer :: ierror
   return
   end subroutine bcastr2d

   subroutine bcastr3d(r)
   real(kind=r8) :: r(:,:,:)
   integer :: ierror
   return
   end subroutine bcastr3d

   subroutine bcastr4d(r)
   real(kind=r8) :: r(:,:,:,:)
   integer :: ierror
   return
   end subroutine bcastr4d

   subroutine bcastr1(r)
   real(kind=r8) :: r
   integer :: ierror
   return
   end subroutine bcastr1

   subroutine bcastchar(w)
   integer(kind=i4) :: ierror
   character(len=*) :: w
   return
   end subroutine bcastchar   
   
   
   function myrank() ! which process am I?
   integer(kind=i4) :: myrank
   myrank=irank
   end function myrank

   function nproc() ! How many of use are there anyway?
   integer(kind=i4) :: nproc
   nproc=iproc
   end function nproc

   subroutine addalli1(i,isum)
   integer(kind=i4) :: ierror,i,isum
   isum=i
   return
   end subroutine addalli1

   subroutine addalli1d(i,isum)
   integer(kind=i4) :: ierror,i(:),isum(:)
   isum=i
   return
   end subroutine addalli1d

   subroutine addallr1(r,rsum)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r,rsum
   rsum=r
   return
   end subroutine addallr1

   subroutine addallr1d(r,rsum)
   real(kind=r8) :: r(:),rsum(:)
   integer(kind=i4) :: ierror
   rsum=r
   return
   end subroutine addallr1d

   subroutine addallr2d(r,rsum)
   real(kind=r8) :: r(:,:),rsum(:,:)
   integer(kind=i4) :: ierror
   rsum=r
   return
   end subroutine addallr2d

   subroutine addallc1(c,csum)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: c,csum
   csum=c
   return
   end subroutine addallc1

   subroutine addallc1d(c,csum)
   complex(kind=r8) :: c(:),csum(:)
   integer(kind=i4) :: ierror
   csum=c
   return
   end subroutine addallc1d

   subroutine gatheri1(i,igather)
   integer(kind=i4) :: i,igather(:)
   integer :: ierror
   igather(1)=i
   return
   end subroutine gatheri1

   subroutine gatheri1d(i,igather)
   integer(kind=i4) :: i(:),igather(:)
   integer :: ierror
   igather(:)=i
   return
   end subroutine gatheri1d

   subroutine gatherr1(r,rgather)
   real(kind=r8) :: r,rgather(:)
   integer :: ierror
   rgather(1)=r
   return
   end subroutine gatherr1

   subroutine gatherr1d(r,rgather)
   real(kind=r8) :: r(:),rgather(:)
   integer :: ierror
   rgather(:)=r
   return
   end subroutine gatherr1d
  
   subroutine scatteri1(i,iscatter)
   integer(kind=i4) :: i(:),iscatter
   integer :: ierror   
   iscatter=i(1)
   return
   end subroutine scatteri1     
   
   subroutine scatteri1d(i,iscatter)
   integer(kind=i4) :: i(:),iscatter(:)
   integer :: ierror,sizenow
   iscatter(:)=i 
   return
   end subroutine scatteri1d
   
   subroutine scatterr1(r,rscatter)
   real(kind=r8) :: r(:),rscatter
   integer :: ierror
   rscatter=r(1)
   return
   end subroutine scatterr1   
   
   subroutine scatterr1d(r,rscatter)
   real(kind=r8) :: r(:),rscatter(:)
   integer :: ierror,sizenow
   rscatter(:)=r
   return
   end subroutine scatterr1d   
   
   subroutine sendi1(i,idto,itag)
   integer(kind=i4) :: i,idto,itag,ierror
   end subroutine sendi1

   subroutine sendi1d(i,idto,itag)
   integer(kind=i4) :: i(:),idto,itag,ierror
   end subroutine sendi1d

   subroutine sendr1(r,idto,itag)
   integer(kind=i4) :: idto,itag,ierror
   real(kind=r8) :: r
   end subroutine sendr1

   subroutine sendr1d(r,idto,itag)
   integer(kind=i4) :: idto,itag,ierror
   real(kind=r8) :: r(:)
   end subroutine sendr1d

   subroutine sendc1(c,idto,itag)
   integer(kind=i4) :: idto,itag,ierror
   complex (kind=r8) :: c
   end subroutine sendc1

   subroutine sendc1d(c,idto,itag)
   integer(kind=i4) :: idto,itag,ierror
   complex(kind=r8) :: c(:)
   end subroutine sendc1d

   subroutine recvi1(i,idfrom,itag)
   integer(kind=i4) :: i,idfrom,itag,ierror
   end subroutine recvi1

   subroutine recvi1d(i,idfrom,itag)
   integer(kind=i4) :: i(:),idfrom,itag,ierror
   end subroutine recvi1d

   subroutine recvr1(r,idfrom,itag)
   integer(kind=i4) :: idfrom,itag,ierror
   real(kind=r8) :: r
   end subroutine recvr1

   subroutine recvr1d(r,idfrom,itag)
   integer(kind=i4) :: idfrom,itag,ierror
   real(kind=r8) :: r(:)
   end subroutine recvr1d

   subroutine recvc1(c,idfrom,itag)
   integer(kind=i4) :: idfrom,itag,ierror
   complex (kind=r8) :: c
   end subroutine recvc1

   subroutine recvc1d(c,idfrom,itag)
   integer(kind=i4) :: idfrom,itag,ierror
   complex(kind=r8) :: c(:)
   end subroutine recvc1d   
   
   subroutine barrier ! wrapper for mpi_barrier
   integer(kind=i4) :: ierror
   end subroutine barrier   
   
   subroutine abort
   integer :: ierror
   stop 'abort, stop!'
   end subroutine abort

   subroutine mpiwait
   return
   end subroutine mpiwait   
   
   
 
end module mympi
