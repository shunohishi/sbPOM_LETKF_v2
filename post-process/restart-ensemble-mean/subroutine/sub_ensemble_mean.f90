subroutine ensemble_mean(nmem,im,jm,km,dat,mean)

  use mpi
  !$ use omp_lib
  implicit none

  !---Common
  integer i,j,k
  integer ierr
  
  !---IN
  integer,intent(in) :: nmem
  integer,intent(in) :: im,jm,km
  real(kind = 4),intent(in) :: dat(im,jm,km)

  !---OUT
  real(kind = 4),intent(out) :: mean(im,jm,km)

  call MPI_ALLREDUCE(dat,mean,im*jm*km,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

  !$omp parallel do private(i,j,k) collapse(2)
  do k=1,km
     do j=1,jm
        do i=1,im
           mean(i,j,k)=mean(i,j,k)/real(nmem)
        end do
     end do
  end do
  !$omp end parallel do

end subroutine ensemble_mean
