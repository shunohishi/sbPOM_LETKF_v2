subroutine ens_mpi_allgather(master_rank,nens,im,jm,dat,edat)

  use mpi
  implicit none

  !Common
  integer i,j,iens
  integer inum
  integer ierr
  
  real(kind = 8) send(im*jm),rec(nens*im*jm)

  !IN
  integer,intent(in) :: master_rank
  integer,intent(in) :: nens,im,jm
  real(kind = 8),intent(in) :: dat(im,jm)

  !OUT
  real(kind = 8),intent(out) :: edat(nens,im,jm)

  inum=0
  do j=1,jm
     do i=1,im
        inum=inum+1
        send(inum)=dat(i,j)
     end do
  end do

  call MPI_Gather(send(1),im*jm,MPI_DOUBLE_PRECISION, &
       & rec(1),im*jm,MPI_DOUBLE_PRECISION, &
       & master_rank,MPI_COMM_WORLD,ierr)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  call MPI_Bcast(rec(1),nens*im*jm,MPI_DOUBLE_PRECISION, &
       & master_rank,MPI_COMM_WORLD,ierr)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  inum=0
  do iens=1,nens
     do j=1,jm
        do i=1,im
           inum=inum+1
           edat(iens,i,j)=rec(inum)
        end do
     end do
  end do
  
end subroutine ens_mpi_allgather

!----------------------------------

subroutine ens_mpi_bcast(master_rank,nens,im,jm,dat)

  use mpi
  implicit none

  integer iens,i,j
  integer inum
  integer ierr
  
  real(kind = 8) rec(nens*im*jm)

  !IN
  integer,intent(in) :: master_rank
  integer,intent(in) :: nens,im,jm

  !INOUT
  real(kind = 8),intent(inout) :: dat(nens,im,jm)
  
  inum=0
  do iens=1,nens
     do j=1,jm
        do i=1,im
           inum=inum+1
           rec(inum)=dat(iens,i,j)
        end do
     end do
  end do

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_Bcast(rec(1),nens*im*jm,MPI_DOUBLE_PRECISION, &
       & master_rank,MPI_COMM_WORLD,ierr)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  inum=0
  do iens=1,nens
     do j=1,jm
        do i=1,im
           inum=inum+1
           dat(iens,i,j)=rec(inum)
        end do
     end do
  end do

end subroutine ens_mpi_bcast

!-----------------------------------------

subroutine cal_mpi_ens(nens,im,jm,km,dat,ens,factor)

  use mpi
  implicit none

  integer i,j,k
  integer inum
  integer ierr

  real(kind = 8) send(im*jm*km)
  real(kind = 8) rec(im*jm*km)
  real(kind = 8) mean(im,jm,km)
  real(kind = 8) inv_ens

  !IN
  integer,intent(in) ::nens
  integer,intent(in) :: im,jm,km
  
  real(kind = 8),intent(in) :: dat(im,jm,km)
  real(kind = 8),intent(in) :: factor

  !INOUT
  real(kind = 8),intent(inout) :: ens(im,jm,km)

  inum=0
  do k=1,km
     do j=1,jm
        do i=1,im
           inum=inum+1
           send(inum)=ens(i,j,k)
        end do
     end do
  end do

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(send(1),rec(1),im*jm*km,MPI_DOUBLE_PRECISION, &
       & MPI_SUM,MPI_COMM_WORLD,ierr)

  inum=0
  inv_ens=1.d0/dble(nens)
  do k=1,km
     do j=1,jm
        do i=1,im
           inum=inum+1
           mean(i,j,k)=rec(inum)*inv_ens
        end do
     end do
  end do

  do k=1,km
     do j=1,jm
        do i=1,im
           ens(i,j,k)=dat(i,j,k)+factor*(ens(i,j,k)-mean(i,j,k))
        end do
     end do
  end do

end subroutine cal_mpi_ens
