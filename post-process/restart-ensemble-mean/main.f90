!-------------------------------------------------------
! Ensemble mean restart file
!-------------------------------------------------------
!
! Calculate ensemble mean from ensemble restart file at end of each month
!
!-------------------------------------------------------

program main

  use mod_gridinfo
  use mod_varname
  use mpi
  implicit none

  !---Parameter
  integer,parameter :: nvar2d=18,nvar3d=19
  integer,parameter :: master_rank=0
  
  !---Common
  integer iyr,imon,iday
  integer nmem
  integer ivar
  integer k

  real(kind = 4) dat2d(im,jm),dat3d(im,jm,km)
  real(kind = 4) mean2d(im,jm),mean3d(im,jm,km)
  
  character(100) dir,letkf,varname

  !---MPI
  integer PETOT,my_rank,ierr
  integer imem  

  !---MPI(Start)
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,PETOT,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
  imem=my_rank+1
  
  if(my_rank == master_rank)then
     write(*,'(a)') "===== Start: Restart ensemble mean ====="
  end if

  !---Read information
  call read_argument(dir,letkf,iyr,imon,iday,nmem)

  if(nmem /= PETOT)then
     if(my_rank == 0)then
        write(*,*) "***Error: nmem /= PETOT ==> ",nmem,PETOT
     end if
     call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
  end if
  
  !---2D
  do ivar=1,nvar2d

     k=1
     call varname2d(ivar,varname)
     call read_restart(dir,letkf,varname,iyr,imon,iday,im,jm,k,imem,dat2d)
     call ensemble_mean(nmem,im,jm,k,dat2d,mean2d)

     if(my_rank == master_rank)then
        write(*,*) "2D:",ivar,"/",nvar2d
        call write_restart(dir,letkf,varname,iyr,imon,iday,im,jm,k,mean2d)
     end if
     call MPI_Barrier(MPI_COMM_WORLD,ierr)  
        
  end do !ivar

  !---3D
  do ivar=1,nvar3d

     call varname3d(ivar,varname)
     call read_restart(dir,letkf,varname,iyr,imon,iday,im,jm,km,imem,dat3d)
     call ensemble_mean(nmem,im,jm,km,dat3d,mean3d)

     if(my_rank == master_rank)then     
        write(*,*) "3D:",ivar,"/",nvar3d        
        call write_restart(dir,letkf,varname,iyr,imon,iday,im,jm,km,mean3d)     
     end if
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     
  end do

  call MPI_Barrier(MPI_COMM_WORLD,ierr)  
  call MPI_FINALIZE(ierr)
  
end program main
