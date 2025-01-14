! parallel_mpi.f

! subroutines for communicating between processors using MPI

!_______________________________________________________________________
subroutine initialize_mpi

  ! set up MPI execution environment and define the POM communicator

  use mpi
  use common_pom_var
  implicit none
  
  integer ierr
  
  ! initiate MPI environment
  call MPI_Init(ierr)

  ! determine processor rank
  call MPI_Comm_rank(MPI_COMM_WORLD,my_task,ierr)
  pom_comm=MPI_COMM_WORLD
  master_task=0
  error_status=0

end subroutine initialize_mpi

!_______________________________________________________________________
subroutine finalize_mpi

  ! terminate the MPI execution environment

  use mpi
  implicit none

  integer ierr
  
  ! terminate MPI environment
  call MPI_Finalize(ierr)

end subroutine finalize_mpi

!_______________________________________________________________________
subroutine distribute_mpi
  
  ! distribute the model domain across processors

  use mpi
  use common_pom_var
  implicit none

  integer i,j
  integer ierr,nproc

  ! determine the number of processors
  call MPI_Comm_size(pom_comm,nproc,ierr)

  ! check number of processors
  if(nproc /= n_proc) then
     error_status=1
     if(my_task == master_task) &
          & write(*,'(a//a)') &
          & 'Incompatible number of processors','POM terminated with error'
     call finalize_mpi
     stop
  end if

  ! determine the number of processors in x
  if(mod(im_global-2,im_local-2) == 0) then
     nproc_x=(im_global-2)/(im_local-2)
  else
     nproc_x=(im_global-2)/(im_local-2) + 1
  end if

  ! determine the number of processors in y
  if(mod(jm_global-2,jm_local-2) == 0) then
     nproc_y=(jm_global-2)/(jm_local-2)
  else
     nproc_y=(jm_global-2)/(jm_local-2) + 1
  end if

  ! check local size
  if(nproc_x*nproc_y > n_proc)then
     error_status=1
     if(my_task == master_task) &
          & write(*,'(a//a)') &
          & 'im_local or jm_local is too low','POM terminated with error'
     call finalize_mpi
     stop
  end if

  ! detemine global and local indices
  im=im_local
  imm1=im-1
  imm2=im-2
  
  jm=jm_local
  jmm1=jm-1
  jmm2=jm-2
  
  kbm1=kb-1
  kbm2=kb-2

  i_global(1:im)=0
  j_global(1:jm)=0
  
  do i=1,im
     
     i_global(i)=i+mod(my_task,nproc_x)*(im-2)
     
     if(i_global(i) > im_global)then
        im=i-1
        i_global(i)=0
     end if
     
  end do

  do j=1,jm
     
     j_global(j)=j+(my_task/nproc_x)*(jm-2)

     if(j_global(j) > jm_global) then
        jm=j-1
        j_global(j)=0
     end if
     
  end do

  ! determine the neighbors (tasks)
  !east
  n_east=my_task+1
  if(mod(n_east,nproc_x) == 0 .and. lqglobal)then
     n_east=my_task-(nproc_x-1)
  else if(mod(n_east,nproc_x) == 0)then
     n_east=-1
  end if

  !west
  n_west=my_task-1
  if(mod(n_west+1,nproc_x) == 0 .and. lqglobal)then
     n_west=my_task+(nproc_x-1)
  else if(mod(n_west+1,nproc_x) == 0)then
     n_west=-1
  end if

  !north
  n_north=my_task+nproc_x
  if(n_north/nproc_x == nproc_y)then
     n_north=-1
  end if

  !south
  n_south=my_task-nproc_x
  if((n_south+nproc_x)/nproc_x == 0)then
     n_south=-1
  end if

  !  write(*,'(a,2i5,a,2i5)') "NS:",n_north,n_south,"EW:",n_east,n_west
  
end subroutine distribute_mpi
!_______________________________________________________________________
subroutine sum0d_mpi(work,to)
  ! send sum of WORK to node TO

  use mpi
  use common_pom_var
  implicit none

  integer itmp
  integer ierr
  
  integer,intent(in) :: to
  integer,intent(inout) :: work

  ! sum data
  call MPI_Reduce(work,itmp,1,MPI_INTEGER,MPI_SUM,to,pom_comm,ierr)
  
  work=itmp
  
end subroutine sum0d_mpi

!_______________________________________________________________________
subroutine bcast0d_mpi(work,from)
  ! send WORK to all nodes from node FROM

  use mpi
  use common_pom_var
  implicit none

  integer ierr
  
  integer,intent(in) :: from
  integer,intent(inout) :: work
  
  !broadcast data
  call MPI_Bcast(work,1,MPI_INTEGER,from,pom_comm,ierr)
  
end subroutine bcast0d_mpi

!_______________________________________________________________________
subroutine exchange2d_mpi(work,nx,ny)
  ! exchange ghost cells around 2d local grids
  ! one band at a time

  use mpi
  use common_pom_var
  implicit none

  integer ierr
  integer istatus(mpi_status_size)
  
  real(kind = r_size) send_east(ny),recv_west(ny)
  real(kind = r_size) send_west(ny),recv_east(ny)
  real(kind = r_size) send_north(nx),recv_south(nx)
  real(kind = r_size) send_south(nx),recv_north(nx)

  integer,intent(in) ::  nx,ny
  real(kind = r_size),intent(inout) :: work(nx,ny)
      
  ! send ghost cell data to the east
  if(n_east /= -1)then

     send_east(1:ny)=work(nx-1,1:ny)

     if(r_size == kind(0.0d0))then
        call MPI_Send(send_east,ny,MPI_DOUBLE_PRECISION, &
             & n_east,my_task,pom_comm,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Send(send_east,ny,MPI_REAL, &
             & n_east,my_task,pom_comm,ierr)
     end if
             
  end if

  ! recieve ghost cell data from the west
  if(n_west /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Recv(recv_west,ny,MPI_DOUBLE_PRECISION, &
             & n_west,n_west,pom_comm,istatus,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Recv(recv_west,ny,MPI_REAL, &
             & n_west,n_west,pom_comm,istatus,ierr)
     end if
     
     work(1,1:ny)=recv_west(1:ny)
        
  end if

  ! send ghost cell data to the west
  if(n_west /= -1)then

     send_west(1:ny)=work(2,1:ny)

     if(r_size == kind(0.0d0))then
        call MPI_Send(send_west,ny,MPI_DOUBLE_PRECISION, &
             & n_west,my_task,pom_comm,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Send(send_west,ny,MPI_REAL, &
             & n_west,my_task,pom_comm,ierr)
     end if
        
  end if
  
  ! recieve ghost cell data from the east
  if(n_east /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Recv(recv_east,ny,MPI_DOUBLE_PRECISION, &
             & n_east,n_east,pom_comm,istatus,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Recv(recv_east,ny,MPI_REAL, &
             & n_east,n_east,pom_comm,istatus,ierr)
     end if
     
     work(nx,1:ny)=recv_east(1:ny)
     
  end if
  
  ! send ghost cell data to the north
  if(n_north /= -1)then

     send_north(1:nx)=work(1:nx,ny-1)

     if(r_size == kind(0.0d0))then
        call MPI_Send(send_north,nx,MPI_DOUBLE_PRECISION, &
             & n_north,my_task,pom_comm,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Send(send_north,nx,MPI_REAL, &
             & n_north,my_task,pom_comm,ierr)
     end if
     
  end if

  ! recieve ghost cell data from the south
  if(n_south /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Recv(recv_south,nx,MPI_DOUBLE_PRECISION, &
             & n_south,n_south,pom_comm,istatus,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Recv(recv_south,nx,MPI_REAL, &
             & n_south,n_south,pom_comm,istatus,ierr)
     end if
     
     work(1:nx,1)=recv_south(1:nx)
        
  end if

  ! send ghost cell data to the south
  if(n_south /= -1)then

     send_south(1:nx)=work(1:nx,2)

     if(r_size == kind(0.0d0))then
        call MPI_Send(send_south,nx,MPI_DOUBLE_PRECISION, &
             & n_south,my_task,pom_comm,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Send(send_south,nx,MPI_REAL, &
             & n_south,my_task,pom_comm,ierr)
     end if
        
     
  end if
  
  ! recieve ghost cell data from the north
  if(n_north /= -1)then

     if(r_size == kind(0.0d0))then     
        call MPI_Recv(recv_north,nx,MPI_DOUBLE_PRECISION, &
             & n_north,n_north,pom_comm,istatus,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Recv(recv_north,nx,MPI_REAL, &
             & n_north,n_north,pom_comm,istatus,ierr)
     end if
        
     work(1:nx,ny)=recv_north(1:nx)
        
  end if
  
end subroutine exchange2d_mpi
!_______________________________________________________________________
subroutine exchange3d_mpi_bl(work,nx,ny,nz)
  ! exchange ghost cells around 3d local grids
  ! one band at a time

  !$use omp_lib
  use mpi
  use common_pom_var
  implicit none

  integer i,j,k
  integer ierr,tag
  integer istatus(mpi_status_size)

  real(kind = r_size) send_east(ny*nz),recv_west(ny*nz)
  real(kind = r_size) send_west(ny*nz),recv_east(ny*nz)
  real(kind = r_size) send_north(nx*nz),recv_south(nx*nz)
  real(kind = r_size) send_south(nx*nz),recv_north(nx*nz)
  
  integer,intent(in) :: nx,ny,nz
  real(kind = r_size),intent(inout) :: work(nx,ny,nz)

  tag=0
  
  ! send ghost cell data to the east
  if(n_east /= -1)then

     !$omp parallel
     !$omp do private(i,j,k)
     do k=1,nz
        do j=1,ny
           i=j+(k-1)*ny
           send_east(i)=work(nx-1,j,k)
        end do
     end do
     !$omp end do
     !$omp end parallel

     if(r_size == kind(0.0d0))then
        call MPI_Send(send_east,ny*nz,MPI_DOUBLE_PRECISION, &
             & n_east,tag,pom_comm,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Send(send_east,ny*nz,MPI_REAL, &
             & n_east,tag,pom_comm,ierr)
     end if
        
  end if
  
  ! recieve ghost cell data from the west
  if(n_west /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Recv(recv_west,ny*nz,MPI_DOUBLE_PRECISION, &
             & n_west,tag,pom_comm,istatus,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Recv(recv_west,ny*nz,MPI_REAL, &
             & n_west,tag,pom_comm,istatus,ierr)
     end if
     
     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do j=1,ny
           i=j+(k-1)*ny
           work(1,j,k)=recv_west(i)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end if

  ! send ghost cell data to the west
  if(n_west /= -1)then

     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do j=1,ny
           i=j+(k-1)*ny
           send_west(i)=work(2,j,k)
        end do
     end do
     !$omp end do
     !$omp end parallel

     if(r_size == kind(0.0d0))then
        call MPI_Send(send_west,ny*nz,MPI_DOUBLE_PRECISION, &
             & n_west,tag,pom_comm,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Send(send_west,ny*nz,MPI_REAL, &
             & n_west,tag,pom_comm,ierr)
     end if
     
  end if
  
  ! recieve ghost cell data from the east
  if(n_east /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Recv(recv_east,ny*nz,MPI_DOUBLE_PRECISION, &
             & n_east,tag,pom_comm,istatus,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Recv(recv_east,ny*nz,MPI_REAL, &
             & n_east,tag,pom_comm,istatus,ierr)
     end if
        
     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do j=1,ny
           i=j+(k-1)*ny
           work(nx,j,k)=recv_east(i)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end if
  
  ! send ghost cell data to the north
  if(n_north /= -1)then

     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do i=1,nx
           j=i+(k-1)*nx
           send_north(j)=work(i,ny-1,k)
        end do
     end do
     !$omp end do
     !$omp end parallel

     if(r_size == kind(0.0d0))then
        call MPI_Send(send_north,nx*nz,MPI_DOUBLE_PRECISION, &
             & n_north,tag,pom_comm,ierr)        
     else if(r_size == kind(0.0e0))then
        call MPI_Send(send_north,nx*nz,MPI_REAL, &
             & n_north,tag,pom_comm,ierr)        
     end if
        
  end if
  
  ! recieve ghost cell data from the south
  if(n_south /= -1)then

     if(r_size == kind(0.0d0))then     
        call MPI_Recv(recv_south,nx*nz,MPI_DOUBLE_PRECISION, &
             & n_south,tag,pom_comm,istatus,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Recv(recv_south,nx*nz,MPI_REAL, &
             & n_south,tag,pom_comm,istatus,ierr)
     end if
        
     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do i=1,nx
           j=i+(k-1)*nx
           work(i,1,k)=recv_south(j)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end if
  
  ! send ghost cell data to the south
  if(n_south /= -1)then

     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do i=1,nx
           j=i+(k-1)*nx
           send_south(j)=work(i,2,k)
        end do
     end do
     !$omp end do
     !$omp end parallel

     if(r_size == kind(0.0d0))then     
        call MPI_Send(send_south,nx*nz,MPI_DOUBLE_PRECISION, &
             & n_south,tag,pom_comm,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Send(send_south,nx*nz,MPI_REAL, &
             & n_south,tag,pom_comm,ierr)
     end if
        
  end if
  
  ! recieve ghost cell data from the north
  if(n_north /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Recv(recv_north,nx*nz,MPI_DOUBLE_PRECISION, &
             & n_north,tag,pom_comm,istatus,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Recv(recv_north,nx*nz,MPI_REAL, &
             & n_north,tag,pom_comm,istatus,ierr)
     end if
     
     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do i=1,nx
           j=i+(k-1)*nx
           work(i,ny,k)=recv_north(j)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end if
  
end subroutine exchange3d_mpi_bl

!_______________________________________________________________________
subroutine order2d_mpi(work2,work4,nx,ny)
  
  ! convert a 2nd order 2d matrix to special 4th order 2d matrix

  use mpi
  use common_pom_var
  implicit none

  integer ierr,tag
  integer istatus(mpi_status_size)
  real(kind = r_size) send_east(ny),recv_west(ny)
  real(kind = r_size) send_north(nx),recv_south(nx)

  
  integer,intent(in) :: nx,ny
  real(kind = r_size),intent(in)  :: work2(nx,ny)
  real(kind = r_size),intent(out) :: work4(0:nx,0:ny)

  work4(0:nx,0:ny)=0.d0
  work4(1:nx,1:ny)=work2(1:nx,1:ny)

  tag=0
  
  ! send ghost cell data to the east
  if(n_east /= -1)then

     send_east(1:ny)=work2(nx-2,1:ny)

     if(r_size == kind(0.0d0))then
        call MPI_Send(send_east,ny,MPI_DOUBLE_PRECISION, &
             & n_east,tag,pom_comm,ierr)        
     else if(r_size == kind(0.0e0))then
        call MPI_Send(send_east,ny,MPI_REAL, &
             & n_east,tag,pom_comm,ierr)        
     end if
     
  end if
  
  ! recieve ghost cell data from the west
  if(n_west /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Recv(recv_west,ny,MPI_DOUBLE_PRECISION, &
             & n_west,tag,pom_comm,istatus,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Recv(recv_west,ny,MPI_REAL, &
             & n_west,tag,pom_comm,istatus,ierr)
     end if
     
     work4(0,1:ny)=recv_west(1:ny)
     
  end if

  ! send ghost cell data to the north
  if(n_north /= -1) then

     send_north(1:nx)=work2(1:nx,ny-2)

     if(r_size == kind(0.0d0))then
        call MPI_Send(send_north,nx,MPI_DOUBLE_PRECISION, &
             & n_north,tag,pom_comm,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Send(send_north,nx,MPI_REAL, &
             & n_north,tag,pom_comm,ierr)
     end if
     
  end if
  
  ! recieve ghost cell data from the south
  if(n_south /= -1) then

     if(r_size == kind(0.0d0))then     
        call MPI_Recv(recv_south,nx,MPI_DOUBLE_PRECISION, &
             & n_south,tag,pom_comm,istatus,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Recv(recv_south,nx,MPI_REAL, &
             & n_south,tag,pom_comm,istatus,ierr)
     end if
     
     work4(1:nx,0)=recv_south(1:nx)
        
  end if

end subroutine order2d_mpi

!_______________________________________________________________________
subroutine order3d_mpi_bl(work2,work4,nx,ny,nz)
  ! convert a 2nd order 3d matrix to special 4th order 3d matrix

  !$use omp_lib  
  use mpi
  use common_pom_var
  implicit none
  
  integer i,j,k
  integer ierr,tag
  integer istatus(mpi_status_size)
  
  real(kind = r_size) send_east(ny*nz),recv_west(ny*nz)
  real(kind = r_size) send_north(nx*nz),recv_south(nx*nz)

  
  integer,intent(in) :: nx,ny,nz
  real(kind = r_size),intent(in)  :: work2(nx,ny,nz)
  real(kind = r_size),intent(out) :: work4(0:nx,0:ny,nz)

  work4(0:nx,0:ny,1:nz)=0.d0
  work4(1:nx,1:ny,1:nz)=work2(1:nx,1:ny,1:nz)
  
  tag=0
  
  ! send ghost cell data to the east
  if(n_east /= -1)then

     !$omp parallel
     !$omp do private(i,j,k)          
     do k=1,nz
        do j=1,ny
           i=j+(k-1)*ny
           send_east(i)=work2(nx-2,j,k)
        end do
     end do
     !$omp end do
     !$omp end parallel

     if(r_size == kind(0.0d0))then
        call MPI_Send(send_east,ny*nz,MPI_DOUBLE_PRECISION, &
             & n_east,tag,pom_comm,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Send(send_east,ny*nz,MPI_REAL, &
             & n_east,tag,pom_comm,ierr)
     end if

     
  end if
  
  ! recieve ghost cell data from the west
  if(n_west /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Recv(recv_west,ny*nz,MPI_DOUBLE_PRECISION, &
             & n_west,tag,pom_comm,istatus,ierr)        
     else if(r_size == kind(0.0e0))then
        call MPI_Recv(recv_west,ny*nz,MPI_REAL, &
             & n_west,tag,pom_comm,istatus,ierr)        
     end if
     
     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do j=1,ny
           i=j+(k-1)*ny
           work4(0,j,k)=recv_west(i)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end if

  ! send ghost cell data to the north
  if(n_north /= -1)then

     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do i=1,nx
           j=i+(k-1)*nx
           send_north(j)=work2(i,ny-2,k)
        end do
     end do
     !$omp end do
     !$omp end parallel

     if(r_size == kind(0.0d0))then
        call MPI_Send(send_north,nx*nz,MPI_DOUBLE_PRECISION, &
             & n_north,tag,pom_comm,ierr)        
     else if(r_size == kind(0.0e0))then
        call MPI_Send(send_north,nx*nz,MPI_REAL, &
             & n_north,tag,pom_comm,ierr)        
     end if

  end if
  
  ! recieve ghost cell data from the south
  if(n_south /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Recv(recv_south,nx*nz,MPI_DOUBLE_PRECISION, &
             & n_south,tag,pom_comm,istatus,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Recv(recv_south,nx*nz,MPI_REAL, &
             & n_south,tag,pom_comm,istatus,ierr)
     end if
        
     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do i=1,nx
           j=i+(k-1)*nx
           work4(i,0,k)=recv_south(j)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end if

end subroutine order3d_mpi_bl
!_______________________________________________________________________
subroutine exchange2d_mpi_nb(work,nx,ny)
  ! exchange ghost cells around 2d local grids
  ! one band at a time

  use mpi
  use common_pom_var
  implicit none

  integer ierr,tag,request
  integer istatus(mpi_status_size)
  
  real(kind = r_size) send_east(ny),recv_west(ny)
  real(kind = r_size) send_west(ny),recv_east(ny)
  real(kind = r_size) send_north(nx),recv_south(nx)
  real(kind = r_size) send_south(nx),recv_north(nx)

  integer,intent(in) ::  nx,ny
  real(kind = r_size),intent(inout) :: work(nx,ny)

  tag=0
  
  ! send ghost cell data to the east
  if(n_east /= -1)then

     send_east(1:ny)=work(nx-1,1:ny)

     if(r_size == kind(0.0d0))then
        call MPI_Isend(send_east,ny,MPI_DOUBLE_PRECISION, &
             & n_east,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Isend(send_east,ny,MPI_REAL, &
             & n_east,tag,pom_comm,request,ierr)
     end if
             
  end if

  ! recieve ghost cell data from the west
  if(n_west /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Irecv(recv_west,ny,MPI_DOUBLE_PRECISION, &
             & n_west,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Irecv(recv_west,ny,MPI_REAL, &
             & n_west,tag,pom_comm,request,ierr)
     end if

     call MPI_Wait(request,istatus,ierr)
     
     work(1,1:ny)=recv_west(1:ny)
        
  end if
  
  ! send ghost cell data to the west
  if(n_west /= -1)then

     send_west(1:ny)=work(2,1:ny)

     if(r_size == kind(0.0d0))then
        call MPI_Isend(send_west,ny,MPI_DOUBLE_PRECISION, &
             & n_west,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Isend(send_west,ny,MPI_REAL, &
             & n_west,tag,pom_comm,request,ierr)
     end if
        
  end if
  
  ! recieve ghost cell data from the east
  if(n_east /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Irecv(recv_east,ny,MPI_DOUBLE_PRECISION, &
             & n_east,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Irecv(recv_east,ny,MPI_REAL, &
             & n_east,tag,pom_comm,request,ierr)
     end if

     call MPI_Wait(request,istatus,ierr)
     
     work(nx,1:ny)=recv_east(1:ny)
     
  end if
  
  ! send ghost cell data to the north
  if(n_north /= -1)then

     send_north(1:nx)=work(1:nx,ny-1)

     if(r_size == kind(0.0d0))then
        call MPI_Isend(send_north,nx,MPI_DOUBLE_PRECISION, &
             & n_north,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Isend(send_north,nx,MPI_REAL, &
             & n_north,tag,pom_comm,request,ierr)
     end if
     
  end if

  ! recieve ghost cell data from the south
  if(n_south /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Irecv(recv_south,nx,MPI_DOUBLE_PRECISION, &
             & n_south,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Irecv(recv_south,nx,MPI_REAL, &
             & n_south,tag,pom_comm,request,ierr)
     end if

     call MPI_Wait(request,istatus,ierr)
     
     work(1:nx,1)=recv_south(1:nx)
     
  end if

  ! send ghost cell data to the south
  if(n_south /= -1)then

     send_south(1:nx)=work(1:nx,2)

     if(r_size == kind(0.0d0))then
        call MPI_Isend(send_south,nx,MPI_DOUBLE_PRECISION, &
             & n_south,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Isend(send_south,nx,MPI_REAL, &
             & n_south,tag,pom_comm,request,ierr)
     end if
     
  end if
  
  ! recieve ghost cell data from the north
  if(n_north /= -1)then

     if(r_size == kind(0.0d0))then     
        call MPI_Irecv(recv_north,nx,MPI_DOUBLE_PRECISION, &
             & n_north,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Irecv(recv_north,nx,MPI_REAL, &
             & n_north,tag,pom_comm,request,ierr)
     end if

     call MPI_Wait(request,istatus,ierr)
     
     work(1:nx,ny)=recv_north(1:nx)
        
  end if

  call MPI_Barrier(pom_comm,ierr)
  
end subroutine exchange2d_mpi_nb

!_______________________________________________________________________
subroutine exchange3d_mpi(work,nx,ny,nz)
  ! exchange ghost cells around 3d local grids
  ! one band at a time

  !$use omp_lib
  use mpi
  use common_pom_var
  implicit none

  integer i,j,k
  integer ierr,request,tag
  integer istatus(mpi_status_size)

  real(kind = r_size) send_east(ny*nz),recv_west(ny*nz)
  real(kind = r_size) send_west(ny*nz),recv_east(ny*nz)
  real(kind = r_size) send_north(nx*nz),recv_south(nx*nz)
  real(kind = r_size) send_south(nx*nz),recv_north(nx*nz)
  
  integer,intent(in) :: nx,ny,nz
  real(kind = r_size),intent(inout) :: work(nx,ny,nz)

  tag=0
  
  ! send ghost cell data to the east
  if(n_east /= -1)then

     !$omp parallel
     !$omp do private(i,j,k)
     do k=1,nz
        do j=1,ny
           i=j+(k-1)*ny
           send_east(i)=work(nx-1,j,k)
        end do
     end do
     !$omp end do
     !$omp end parallel

     if(r_size == kind(0.0d0))then
        call MPI_Isend(send_east,ny*nz,MPI_DOUBLE_PRECISION, &
             & n_east,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Isend(send_east,ny*nz,MPI_REAL, &
             & n_east,tag,pom_comm,request,ierr)
     end if
        
  end if
  
  ! recieve ghost cell data from the west
  if(n_west /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Irecv(recv_west,ny*nz,MPI_DOUBLE_PRECISION, &
             & n_west,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Irecv(recv_west,ny*nz,MPI_REAL, &
             & n_west,tag,pom_comm,request,ierr)
     end if

     call MPI_Wait(request,istatus,ierr)
     
     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do j=1,ny
           i=j+(k-1)*ny
           work(1,j,k)=recv_west(i)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end if

  ! send ghost cell data to the west
  if(n_west /= -1)then

     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do j=1,ny
           i=j+(k-1)*ny
           send_west(i)=work(2,j,k)
        end do
     end do
     !$omp end do
     !$omp end parallel

     if(r_size == kind(0.0d0))then
        call MPI_Isend(send_west,ny*nz,MPI_DOUBLE_PRECISION, &
             & n_west,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Isend(send_west,ny*nz,MPI_REAL, &
             & n_west,tag,pom_comm,request,ierr)
     end if
     
  end if
  
  ! recieve ghost cell data from the east
  if(n_east /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Irecv(recv_east,ny*nz,MPI_DOUBLE_PRECISION, &
             & n_east,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Irecv(recv_east,ny*nz,MPI_REAL, &
             & n_east,tag,pom_comm,request,ierr)
     end if

     call MPI_Wait(request,istatus,ierr)
     
     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do j=1,ny
           i=j+(k-1)*ny
           work(nx,j,k)=recv_east(i)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end if
  
  ! send ghost cell data to the north
  if(n_north /= -1)then

     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do i=1,nx
           j=i+(k-1)*nx
           send_north(j)=work(i,ny-1,k)
        end do
     end do
     !$omp end do
     !$omp end parallel

     if(r_size == kind(0.0d0))then
        call MPI_Isend(send_north,nx*nz,MPI_DOUBLE_PRECISION, &
             & n_north,tag,pom_comm,request,ierr)        
     else if(r_size == kind(0.0e0))then
        call MPI_Isend(send_north,nx*nz,MPI_REAL, &
             & n_north,tag,pom_comm,request,ierr)        
     end if
        
  end if
  
  ! recieve ghost cell data from the south
  if(n_south /= -1)then

     if(r_size == kind(0.0d0))then     
        call MPI_Irecv(recv_south,nx*nz,MPI_DOUBLE_PRECISION, &
             & n_south,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Irecv(recv_south,nx*nz,MPI_REAL, &
             & n_south,tag,pom_comm,request,ierr)
     end if

     call MPI_Wait(request,istatus,ierr)
     
     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do i=1,nx
           j=i+(k-1)*nx
           work(i,1,k)=recv_south(j)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end if
  
  ! send ghost cell data to the south
  if(n_south /= -1)then

     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do i=1,nx
           j=i+(k-1)*nx
           send_south(j)=work(i,2,k)
        end do
     end do
     !$omp end do
     !$omp end parallel

     if(r_size == kind(0.0d0))then     
        call MPI_Isend(send_south,nx*nz,MPI_DOUBLE_PRECISION, &
             & n_south,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Isend(send_south,nx*nz,MPI_REAL, &
             & n_south,tag,pom_comm,request,ierr)
     end if
        
  end if
  
  ! recieve ghost cell data from the north
  if(n_north /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Irecv(recv_north,nx*nz,MPI_DOUBLE_PRECISION, &
             & n_north,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Irecv(recv_north,nx*nz,MPI_REAL, &
             & n_north,tag,pom_comm,request,ierr)
     end if

     call MPI_Wait(request,istatus,ierr)

     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do i=1,nx
           j=i+(k-1)*nx
           work(i,ny,k)=recv_north(j)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end if

  call MPI_Barrier(pom_comm,ierr)
  
end subroutine exchange3d_mpi
!_______________________________________________________________________
subroutine order2d_mpi_nb(work2,work4,nx,ny)
  
  ! convert a 2nd order 2d matrix to special 4th order 2d matrix

  use mpi
  use common_pom_var
  implicit none

  integer ierr,request,tag
  integer istatus(mpi_status_size)
  real(kind = r_size) send_east(ny),recv_west(ny)
  real(kind = r_size) send_north(nx),recv_south(nx)
  
  integer,intent(in) :: nx,ny
  real(kind = r_size),intent(in)  :: work2(nx,ny)
  real(kind = r_size),intent(out) :: work4(0:nx,0:ny)

  tag=0
  
  work4(0:nx,0:ny)=0.d0
  work4(1:nx,1:ny)=work2(1:nx,1:ny)

  ! send ghost cell data to the east
  if(n_east /= -1)then

     send_east(1:ny)=work2(nx-2,1:ny)

     if(r_size == kind(0.0d0))then
        call MPI_Isend(send_east,ny,MPI_DOUBLE_PRECISION, &
             & n_east,tag,pom_comm,request,ierr)        
     else if(r_size == kind(0.0e0))then
        call MPI_Isend(send_east,ny,MPI_REAL, &
             & n_east,tag,pom_comm,request,ierr)        
     end if
     
  end if
  
  ! recieve ghost cell data from the west
  if(n_west /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Irecv(recv_west,ny,MPI_DOUBLE_PRECISION, &
             & n_west,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Irecv(recv_west,ny,MPI_REAL, &
             & n_west,tag,pom_comm,request,ierr)
     end if

     call MPI_Wait(request,istatus,ierr)
     
     work4(0,1:ny)=recv_west(1:ny)
     
  end if

  ! send ghost cell data to the north
  if(n_north /= -1) then

     send_north(1:nx)=work2(1:nx,ny-2)

     if(r_size == kind(0.0d0))then
        call MPI_Isend(send_north,nx,MPI_DOUBLE_PRECISION, &
             & n_north,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Isend(send_north,nx,MPI_REAL, &
             & n_north,tag,pom_comm,request,ierr)
     end if
     
  end if
  
  ! recieve ghost cell data from the south
  if(n_south /= -1) then

     if(r_size == kind(0.0d0))then     
        call MPI_Irecv(recv_south,nx,MPI_DOUBLE_PRECISION, &
             & n_south,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Irecv(recv_south,nx,MPI_REAL, &
             & n_south,tag,pom_comm,request,ierr)
     end if

     call MPI_Wait(request,istatus,ierr)
     
     work4(1:nx,0)=recv_south(1:nx)
        
  end if

  call MPI_Barrier(pom_comm,ierr)

end subroutine order2d_mpi_nb
!_______________________________________________________________________
subroutine order3d_mpi(work2,work4,nx,ny,nz)
  ! convert a 2nd order 3d matrix to special 4th order 3d matrix

  !$use omp_lib  
  use mpi
  use common_pom_var
  implicit none
  
  integer i,j,k
  integer ierr,request,tag
  integer istatus(mpi_status_size)
  
  real(kind = r_size) send_east(ny*nz),recv_west(ny*nz)
  real(kind = r_size) send_north(nx*nz),recv_south(nx*nz)

  
  integer,intent(in) :: nx,ny,nz
  real(kind = r_size),intent(in)  :: work2(nx,ny,nz)
  real(kind = r_size),intent(out) :: work4(0:nx,0:ny,nz)

  tag=0
  
  work4(0:nx,0:ny,1:nz)=0.d0
  work4(1:nx,1:ny,1:nz)=work2(1:nx,1:ny,1:nz)
  
  ! send ghost cell data to the east
  if(n_east /= -1)then

     !$omp parallel
     !$omp do private(i,j,k)          
     do k=1,nz
        do j=1,ny
           i=j+(k-1)*ny
           send_east(i)=work2(nx-2,j,k)
        end do
     end do
     !$omp end do
     !$omp end parallel

     if(r_size == kind(0.0d0))then
        call MPI_Isend(send_east,ny*nz,MPI_DOUBLE_PRECISION, &
             & n_east,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Isend(send_east,ny*nz,MPI_REAL, &
             & n_east,tag,pom_comm,request,ierr)
     end if

     
  end if
  
  ! recieve ghost cell data from the west
  if(n_west /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Irecv(recv_west,ny*nz,MPI_DOUBLE_PRECISION, &
             & n_west,tag,pom_comm,request,ierr)        
     else if(r_size == kind(0.0e0))then
        call MPI_Irecv(recv_west,ny*nz,MPI_REAL, &
             & n_west,tag,pom_comm,request,ierr)        
     end if

     call MPI_Wait(request,istatus,ierr)
     
     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do j=1,ny
           i=j+(k-1)*ny
           work4(0,j,k)=recv_west(i)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end if

  ! send ghost cell data to the north
  if(n_north /= -1)then

     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do i=1,nx
           j=i+(k-1)*nx
           send_north(j)=work2(i,ny-2,k)
        end do
     end do
     !$omp end do
     !$omp end parallel

     if(r_size == kind(0.0d0))then
        call MPI_Isend(send_north,nx*nz,MPI_DOUBLE_PRECISION, &
             & n_north,tag,pom_comm,request,ierr)        
     else if(r_size == kind(0.0e0))then
        call MPI_Isend(send_north,nx*nz,MPI_REAL, &
             & n_north,tag,pom_comm,request,ierr)        
     end if

  end if
  
  ! recieve ghost cell data from the south
  if(n_south /= -1)then

     if(r_size == kind(0.0d0))then
        call MPI_Irecv(recv_south,nx*nz,MPI_DOUBLE_PRECISION, &
             & n_south,tag,pom_comm,request,ierr)
     else if(r_size == kind(0.0e0))then
        call MPI_Irecv(recv_south,nx*nz,MPI_REAL, &
             & n_south,tag,pom_comm,request,ierr)
     end if

     call MPI_Wait(request,istatus,ierr)
     
     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,nz
        do i=1,nx
           j=i+(k-1)*nx
           work4(i,0,k)=recv_south(j)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end if

  call MPI_Barrier(pom_comm,ierr)
  
end subroutine order3d_mpi
