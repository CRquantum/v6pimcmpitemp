module mympi
   implicit none
   include 'mpif.h'
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer, private, save :: mpii4,mpii8,mpir8
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
   module procedure gatheri1,gatheri1d,gatheri81
   module procedure gatherr1,gatherr1d
end interface gather

interface scatter ! scatter from process 0 to all process, evenly scatter. inverse of gather
   module procedure scatteri1,scatteri1d,scatteri81
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
   use mpi
   integer :: ierror,isize,ir,ip
   integer(kind=i4) :: itest4
   integer(kind=i8) :: itest8
   real(kind=r8) :: rtest8
   call mpi_init(ierror)
   call mpi_comm_rank(mpi_comm_world,ir,ierror)
   irank=ir
   call mpi_comm_size(mpi_comm_world,ip,ierror)
   iproc=ip
   call mpi_sizeof(itest4,isize,ierror)
   call mpi_type_match_size(mpi_typeclass_integer,isize,mpii4,ierror)
   call mpi_sizeof(itest8,isize,ierror)
   call mpi_type_match_size(mpi_typeclass_integer,isize,mpii8,ierror)
   call mpi_sizeof(rtest8,isize,ierror)
   call mpi_type_match_size(mpi_typeclass_real,isize,mpir8,ierror)
   return
   end subroutine init0

   subroutine done ! wrapper for finalize routine
   integer :: ierror
   call mpi_finalize(ierror)
   return
   end subroutine done

   subroutine bcastlogi(i)
   logical :: i
   integer :: ierror
   call mpi_bcast(i,1,mpi_logical,0,mpi_comm_world,ierror)
   return
   end subroutine bcastlogi   
   
   subroutine bcastirn(i)
   integer(kind=i8) :: i
   integer :: ierror
   call mpi_bcast(i,1,mpii8,0,mpi_comm_world,ierror)
   return
   end subroutine bcastirn

   subroutine bcasti1(i)
   integer(kind=i4) :: i
   integer :: ierror
   call mpi_bcast(i,1,mpii4,0,mpi_comm_world,ierror)
   return
   end subroutine bcasti1

   subroutine bcasti1d(i)
   integer(kind=i4) :: i(:)
   integer :: ierror
   call mpi_bcast(i,size(i),mpii4,0,mpi_comm_world,ierror)
   return
   end subroutine bcasti1d

   subroutine bcasti2d(i)
   integer(kind=i4) :: i(:,:)
   integer :: ierror
   call mpi_bcast(i,size(i),mpii4,0,mpi_comm_world,ierror)
   return
   end subroutine bcasti2d

   subroutine bcasti3d(i)
   integer(kind=i4) :: i(:,:,:)
   integer :: ierror
   call mpi_bcast(i,size(i),mpii4,0,mpi_comm_world,ierror)
   return
   end subroutine bcasti3d

   subroutine bcasti4d(i)
   integer(kind=i4) :: i(:,:,:,:)
   integer :: ierror
   call mpi_bcast(i,size(i),mpii4,0,mpi_comm_world,ierror)
   return
   end subroutine bcasti4d

   subroutine bcastr1d(r)
   real(kind=r8) :: r(:)
   integer :: ierror
   call mpi_bcast(r,size(r),mpir8,0,mpi_comm_world,ierror)
   return
   end subroutine bcastr1d

   subroutine bcastr2d(r)
   real(kind=r8) :: r(:,:)
   integer :: ierror
   call mpi_bcast(r,size(r),mpir8,0,mpi_comm_world,ierror)
   return
   end subroutine bcastr2d

   subroutine bcastr3d(r)
   real(kind=r8) :: r(:,:,:)
   integer :: ierror
   call mpi_bcast(r,size(r),mpir8,0,mpi_comm_world,ierror)
   end subroutine bcastr3d

   subroutine bcastr4d(r)
   real(kind=r8) :: r(:,:,:,:)
   integer :: ierror
   call mpi_bcast(r,size(r),mpir8,0,mpi_comm_world,ierror)
   end subroutine bcastr4d

   subroutine bcastr1(r)
   real(kind=r8) :: r
   integer :: ierror
   call mpi_bcast(r,1,mpir8,0,mpi_comm_world,ierror)
   end subroutine bcastr1

   subroutine bcastchar(w)
   integer(kind=i4) :: ierror
   character(len=*) :: w
   call mpi_bcast(w,len(w),mpi_character,0,mpi_comm_world,ierror)
   return
   end subroutine bcastchar   
   
   
   function myrank() ! which process am I?
   integer(kind=i4) :: myrank
   myrank=irank
   return
   end function myrank

   function nproc() ! How many of use are there anyway?
   integer(kind=i4) :: nproc
   nproc=iproc
   return
   end function nproc

   subroutine addalli1(i,isum)
   integer(kind=i4) :: ierror,i,isum
   call mpi_reduce(i,isum,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierror)
   return
   end subroutine addalli1

   subroutine addalli1d(i,isum)
   integer(kind=i4) :: ierror,i(:),isum(:)
   call mpi_reduce(i,isum,size(i),mpi_integer,mpi_sum,0,mpi_comm_world,ierror)
   return
   end subroutine addalli1d

   subroutine addallr1(r,rsum)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r,rsum
   call mpi_reduce(r,rsum,1,mpi_double_precision,mpi_sum,0, &
      mpi_comm_world,ierror)
   return
   end subroutine addallr1

   subroutine addallr1d(r,rsum)
   real(kind=r8) :: r(:),rsum(:)
   integer(kind=i4) :: ierror
   call mpi_reduce(r,rsum,size(r),mpi_double_precision,mpi_sum,0, &
      mpi_comm_world,ierror)
   return
   end subroutine addallr1d

   subroutine addallr2d(r,rsum)
   real(kind=r8) :: r(:,:),rsum(:,:)
   integer(kind=i4) :: ierror
   call mpi_reduce(r,rsum,size(r),mpi_double_precision,mpi_sum,0, &
      mpi_comm_world,ierror)
   return
   end subroutine addallr2d

   subroutine addallc1(c,csum)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: c,csum
   call mpi_reduce(c,csum,1,mpi_double_complex,mpi_sum,0, &
      mpi_comm_world,ierror)
   return
   end subroutine addallc1

   subroutine addallc1d(c,csum)
   complex(kind=r8) :: c(:),csum(:)
   integer(kind=i4) :: ierror
   call mpi_reduce(c,csum,size(c),mpi_double_complex,mpi_sum,0, &
      mpi_comm_world,ierror)
   return
   end subroutine addallc1d

   subroutine gatheri1(i,igather)
   integer(kind=i4) :: i,igather(:)
   integer :: ierror
   call mpi_gather(i,1,mpii4,igather,1,mpii4,0,mpi_comm_world,ierror)
   return
   end subroutine gatheri1

   subroutine gatheri1d(i,igather)
   integer(kind=i4) :: i(:),igather(:)
   integer :: ierror
   call mpi_gather(i,size(i),mpii4,igather,size(i),mpii4,0, &
      mpi_comm_world,ierror)
   return
   end subroutine gatheri1d

   subroutine gatheri81(i,igather)
   integer(kind=i8) :: i,igather(:)
   integer :: ierror
   call mpi_gather(i,1,mpii8,igather,1,mpii8,0, &
      mpi_comm_world,ierror)
   return
   end subroutine gatheri81   

   subroutine gatherr1(r,rgather)
   real(kind=r8) :: r,rgather(:)
   integer :: ierror
   call mpi_gather(r,1,mpir8,rgather,1,mpir8,0,mpi_comm_world,ierror)
   return
   end subroutine gatherr1

   subroutine gatherr1d(r,rgather)
   real(kind=r8) :: r(:),rgather(:)
   integer :: ierror
   call mpi_gather(r,size(r),mpir8,rgather,size(r),mpir8,0, &
      mpi_comm_world,ierror)
   return
   end subroutine gatherr1d
  
   subroutine scatteri1(i,iscatter)
   integer(kind=i4) :: i(:),iscatter
   integer :: ierror   
   call mpi_scatter(i,1,mpii4,iscatter,1,mpii4,0, &
      mpi_comm_world,ierror)  
   return
   end subroutine scatteri1    
   
   subroutine scatteri81(i,iscatter)
   integer(kind=i8) :: i(:),iscatter
   integer :: ierror   
   call mpi_scatter(i,1,mpii8,iscatter,1,mpii8,0, &
      mpi_comm_world,ierror)  
   return
   end subroutine scatteri81
   
   subroutine scatteri1d(i,iscatter)
   integer(kind=i4) :: i(:),iscatter(:)
   integer :: ierror
   call mpi_scatter(i,size(iscatter),mpii4,iscatter,size(iscatter),mpii4,0, &
      mpi_comm_world,ierror) 
   return
   end subroutine scatteri1d
   
   subroutine scatterr1(r,rscatter)
   real(kind=r8) :: r(:),rscatter
   integer :: ierror
   call mpi_scatter(r,1,mpir8,rscatter,1,mpir8,0, &
      mpi_comm_world,ierror)
   return
   end subroutine scatterr1   
   
   subroutine scatterr1d(r,rscatter)
   real(kind=r8) :: r(:),rscatter(:)
   integer :: ierror
   call mpi_scatter(r,size(rscatter),mpir8,rscatter,size(rscatter),mpir8,0, &
      mpi_comm_world,ierror)
   return
   end subroutine scatterr1d   
   
   subroutine sendi1(i,idto,itag)
   integer(kind=i4) :: i,idto,itag,ierror
   call mpi_send(i,1,mpi_integer,idto,itag,mpi_comm_world,ierror)
   return
   end subroutine sendi1

   subroutine sendi1d(i,idto,itag)
   integer(kind=i4) :: i(:),idto,itag,ierror
   call mpi_send(i,size(i),mpi_integer,idto,itag,mpi_comm_world,ierror)
   return
   end subroutine sendi1d

   subroutine sendr1(r,idto,itag)
   integer(kind=i4) :: idto,itag,ierror
   real(kind=r8) :: r
   call mpi_send(r,1,mpi_double_precision,idto,itag,mpi_comm_world,ierror)
   return
   end subroutine sendr1

   subroutine sendr1d(r,idto,itag)
   integer(kind=i4) :: idto,itag,ierror
   real(kind=r8) :: r(:)
   call mpi_send(r,size(r),mpi_double_precision,idto,itag,mpi_comm_world,ierror)
   return
   end subroutine sendr1d

   subroutine sendc1(c,idto,itag)
   integer(kind=i4) :: idto,itag,ierror
   complex (kind=r8) :: c
   call mpi_send(c,1,mpi_double_complex,idto,itag,mpi_comm_world,ierror)
   return
   end subroutine sendc1

   subroutine sendc1d(c,idto,itag)
   integer(kind=i4) :: idto,itag,ierror
   complex(kind=r8) :: c(:)
   call mpi_send(c,size(c),mpi_double_complex,idto,itag,mpi_comm_world,ierror)
   return
   end subroutine sendc1d

   subroutine recvi1(i,idfrom,itag)
   integer(kind=i4) :: i,idfrom,itag,ierror,STATUS(MPI_STATUS_SIZE)  
   call mpi_recv(i,1,mpi_integer,idfrom,itag,mpi_comm_world,status,ierror)
   return
   end subroutine recvi1

   subroutine recvi1d(i,idfrom,itag)
   integer(kind=i4) :: i(:),idfrom,itag,ierror,STATUS(MPI_STATUS_SIZE)  
   call mpi_recv(i,size(i),mpi_integer,idfrom,itag,mpi_comm_world,status,ierror)
   return
   end subroutine recvi1d

   subroutine recvr1(r,idfrom,itag)
   integer(kind=i4) :: idfrom,itag,ierror,STATUS(MPI_STATUS_SIZE)  
   real(kind=r8) :: r
   call mpi_recv(r,1,mpi_double_precision,idfrom,itag,mpi_comm_world,status,ierror)
   return
   end subroutine recvr1

   subroutine recvr1d(r,idfrom,itag)
   integer(kind=i4) :: idfrom,itag,ierror,STATUS(MPI_STATUS_SIZE)  
   real(kind=r8) :: r(:)
   call mpi_recv(r,size(r),mpi_double_precision,idfrom,itag,mpi_comm_world, &
      status,ierror)
   return
   end subroutine recvr1d

   subroutine recvc1(c,idfrom,itag)
   integer(kind=i4) :: idfrom,itag,ierror,STATUS(MPI_STATUS_SIZE)  
   complex (kind=r8) :: c
   call mpi_recv(c,1,mpi_double_complex,idfrom,itag,mpi_comm_world,status,ierror)
   return
   end subroutine recvc1

   subroutine recvc1d(c,idfrom,itag)
   integer(kind=i4) :: idfrom,itag,ierror,STATUS(MPI_STATUS_SIZE)  
   complex(kind=r8) :: c(:)
   call mpi_recv(c,size(c),mpi_double_complex,idfrom,itag,mpi_comm_world,status,ierror)
   return
   end subroutine recvc1d   
   
   subroutine barrier ! wrapper for mpi_barrier
   integer(kind=i4) :: ierror
   call mpi_barrier(mpi_comm_world,ierror)
   return
   end subroutine barrier   
   
   subroutine abort
   integer :: ierror
   call mpi_abort(mpi_comm_world,ierror)
   return
   end subroutine abort

   subroutine mpiwait
   integer :: REQUEST, STATUS(MPI_STATUS_SIZE), IERROR
   call MPI_WAIT(REQUEST, STATUS, IERROR)
   return
   end subroutine mpiwait   
   
   
 
end module mympi
