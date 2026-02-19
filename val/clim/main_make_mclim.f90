program main

  use setting
  use mod_gridinfo, im_lora => im, jm_lora => jm, km_lora => km
  use mod_read_glorys025, im_g025 => im, jm_g025 => jm, km_g025 => km
  use mod_make_ncfile
  use mod_io
  use mpi
  implicit none

  !---Parameter
  integer,parameter :: nvar2d=1 !SSH
  integer,parameter :: nvar3d=4 !T/S/U/V

  !---Common
  integer status,system
  integer idat,ivar
  integer im,jm,km

  integer ib,nb
  integer i
  integer imon,nday
  integer,allocatable :: iyr(:),iday(:)

  character(1) varname
  character(100) filename_clim,filename_mclim

  !---Analysis
  real(kind = 8),allocatable :: lont(:),lonu(:),lonv(:)
  real(kind = 8),allocatable :: latt(:),latu(:),latv(:)
  real(kind = 8),allocatable :: dept(:,:,:),depu(:,:,:),depv(:,:,:)
  real(kind = 8),allocatable :: mask(:,:),maskt(:,:),masku(:,:),maskv(:,:)

  real(kind = 8),allocatable :: mean2d(:,:),sprd2d(:,:)
  real(kind = 8),allocatable :: mean3d(:,:,:),sprd3d(:,:,:)

  !---Climatology
  integer,allocatable :: mpass_mean2d(:,:),mmiss_mean2d(:,:)
  integer,allocatable :: mpass_sprd2d(:,:),mmiss_sprd2d(:,:)
  integer,allocatable :: pass_mean2d(:,:),miss_mean2d(:,:)
  integer,allocatable :: pass_sprd2d(:,:),miss_sprd2d(:,:)
  
  real(kind = 8),allocatable :: mclim_mean2d(:,:)
  real(kind = 8),allocatable :: mclim_sprd2d(:,:)
  real(kind = 8),allocatable :: clim_mean2d(:,:)
  real(kind = 8),allocatable :: clim_sprd2d(:,:)

  integer,allocatable :: mpass_mean3d(:,:,:),mmiss_mean3d(:,:,:)
  integer,allocatable :: mpass_sprd3d(:,:,:),mmiss_sprd3d(:,:,:)
  integer,allocatable :: pass_mean3d(:,:,:),miss_mean3d(:,:,:)
  integer,allocatable :: pass_sprd3d(:,:,:),miss_sprd3d(:,:,:)
  
  real(kind = 8),allocatable :: mclim_mean3d(:,:,:)
  real(kind = 8),allocatable :: mclim_sprd3d(:,:,:)
  real(kind = 8),allocatable :: clim_mean3d(:,:,:)
  real(kind = 8),allocatable :: clim_sprd3d(:,:,:)

  !---MPI
  integer PETOT,my_rank,ierr
  integer :: master_rank=0

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,PETOT,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

  if(my_rank == master_rank)then
     write(*,*) "=== Start monthly and annual climatologies ==="
  end if
     
  do idat=1,ndat

     if(idat == 1)then
        im=im_lora
        jm=jm_lora
        km=km_lora
     else
        im=im_g025
        jm=jm_g025
        km=km_g025        
     end if

     !===Make output netcdf file
     filename_clim="dat/"//trim(datname(idat))//"_clim.nc"
     filename_mclim="dat/"//trim(datname(idat))//"_mclim.nc"

     if(my_rank == master_rank)then
        write(*,*) "Make NetCDF file"
        status=system("rm -f "//trim(filename_clim)//" "//trim(filename_mclim))
        call make_ncfile(im,jm,km,1,filename_clim)
        call make_ncfile(im,jm,km,12,filename_mclim)
     end if
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     
     !===Grid information
     !---Allocate
     allocate(lont(im),lonu(im),lonv(im))
     allocate(latt(jm),latu(jm),latv(jm))
     allocate(dept(im,jm,km),depu(im,jm,km),depv(im,jm,km))
     allocate(mask(im,jm),maskt(im,jm),masku(im,jm),maskv(im,jm))

     !---Read grid
     !write(*,*) "Read grid"
     call read_grid(idat,im,jm,km,lont,lonu,lonv,latt,latu,latv,dept,depu,depv,maskt,masku,maskv)

     !---Write grid
     if(my_rank == master_rank)then
        call write_grid(filename_clim,im,jm,km,lont,lonu,lonv,latt,latu,latv,dept,depu,depv,maskt,masku,maskv)
        call write_grid(filename_mclim,im,jm,km,lont,lonu,lonv,latt,latu,latv,dept,depu,depv,maskt,masku,maskv)
     end if

     !---Deallocate
     deallocate(lont,lonu,lonv)
     deallocate(latt,latu,latv)
     deallocate(dept,depu,depv)

     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     
     !===2D
     !---Allocate
     allocate(mean2d(im,jm),sprd2d(im,jm))     
     allocate(mclim_mean2d(im,jm),mpass_mean2d(im,jm),mmiss_mean2d(im,jm))
     allocate(mclim_sprd2d(im,jm),mpass_sprd2d(im,jm),mmiss_sprd2d(im,jm))
     allocate(clim_mean2d(im,jm),pass_mean2d(im,jm),miss_mean2d(im,jm))
     allocate(clim_sprd2d(im,jm),pass_sprd2d(im,jm),miss_sprd2d(im,jm))

     do ivar=1,nvar2d

        if(ivar == 1)then
           varname="h"
           mask(:,:)=maskt(:,:)
        else
           write(*,*) "***Error ivar in nvar2d ==> ",ivar
           stop
        end if

        call init(im,jm,1,clim_mean2d,pass_mean2d,miss_mean2d)
        call init(im,jm,1,clim_sprd2d,pass_sprd2d,miss_sprd2d)
        
        do imon=1,12

           if(my_rank == master_rank)then
              write(*,*) "Varname:",trim(varname)," Month:",imon
           end if

           call init(im,jm,1,mclim_mean2d,mpass_mean2d,mmiss_mean2d)
           call init(im,jm,1,mclim_sprd2d,mpass_sprd2d,mmiss_sprd2d)
           
           call get_total_day(imon,syr,eyr,nday)
           allocate(iyr(nday),iday(nday))
           call get_iyr_iday(imon,syr,eyr,nday,iyr,iday)
           nb=nday/PETOT+1

           do ib=1,nb

              i=(my_rank+1)+(ib-1)*PETOT
              if(nday < i) cycle

              !---Read data
              call read_dat(varname,idat,iyr(i),imon,iday(i),im,jm,1,mask,mean2d,sprd2d)

              !---Add data
              call add(im,jm,1,mask,mean2d,mclim_mean2d,mpass_mean2d,mmiss_mean2d)
              call add(im,jm,1,mask,sprd2d,mclim_sprd2d,mpass_sprd2d,mmiss_sprd2d)
              call add(im,jm,1,mask,mean2d,clim_mean2d,pass_mean2d,miss_mean2d)
              call add(im,jm,1,mask,sprd2d,clim_sprd2d,pass_sprd2d,miss_sprd2d)
              
           end do !ib

           deallocate(iyr,iday)
           
           !---Average
           call ave_mpi(im,jm,1,mclim_mean2d,mpass_mean2d,mmiss_mean2d)
           call ave_mpi(im,jm,1,mclim_sprd2d,mpass_sprd2d,mmiss_sprd2d)

           !---Write monthly climatology
           if(my_rank == master_rank)then
              call write_clim_nc(filename_mclim,varname,im,jm,1,imon,mclim_mean2d,mclim_sprd2d)              
           end if
           call MPI_Barrier(MPI_COMM_WORLD,ierr)

        end do !imon

        !---Average
        call ave_mpi(im,jm,1,clim_mean2d,pass_mean2d,miss_mean2d)
        call ave_mpi(im,jm,1,clim_sprd2d,pass_sprd2d,miss_sprd2d)

        !---Write climatology
        if(my_rank == master_rank)then
           call write_clim_nc(filename_clim,varname,im,jm,1,0,clim_mean2d,clim_sprd2d)
        end if
        call MPI_Barrier(MPI_COMM_WORLD,ierr)

     end do !ivar
     
     !---Deallocate
     deallocate(mean2d,sprd2d)
     deallocate(mclim_mean2d,mpass_mean2d,mmiss_mean2d)
     deallocate(mclim_sprd2d,mpass_sprd2d,mmiss_sprd2d)
     deallocate(clim_mean2d,pass_mean2d,miss_mean2d)
     deallocate(clim_sprd2d,pass_sprd2d,miss_sprd2d)
     
     !===3D
     !---Allocate
     allocate(mean3d(im,jm,km),sprd3d(im,jm,km))
     allocate(mclim_mean3d(im,jm,km),mpass_mean3d(im,jm,km),mmiss_mean3d(im,jm,km))
     allocate(mclim_sprd3d(im,jm,km),mpass_sprd3d(im,jm,km),mmiss_sprd3d(im,jm,km))
     allocate(clim_mean3d(im,jm,km),pass_mean3d(im,jm,km),miss_mean3d(im,jm,km))
     allocate(clim_sprd3d(im,jm,km),pass_sprd3d(im,jm,km),miss_sprd3d(im,jm,km))

     do ivar=1,nvar3d

        if(ivar == 1)then
           varname="t"
           mask(:,:)=maskt(:,:)
        else if(ivar == 2)then
           varname="s"
           mask(:,:)=maskt(:,:)           
        else if(ivar == 3)then
           varname="u"
           mask(:,:)=masku(:,:)           
        else if(ivar == 4)then
           varname="v"
           mask(:,:)=maskv(:,:)
        else
           write(*,*) "***Error: ivar in nvar3d ==> ",ivar
           stop           
        end if
        
        call init(im,jm,km,clim_mean3d,pass_mean3d,miss_mean3d)
        call init(im,jm,km,clim_sprd3d,pass_sprd3d,miss_sprd3d)
        
        do imon=1,12

           if(my_rank == master_rank)then
              write(*,*) "Varname:",trim(varname)," Month:",imon
           end if
           
           call init(im,jm,km,mclim_mean3d,mpass_mean3d,mmiss_mean3d)
           call init(im,jm,km,mclim_sprd3d,mpass_sprd3d,mmiss_sprd3d)           

           call get_total_day(imon,syr,eyr,nday)
           allocate(iyr(nday),iday(nday))
           call get_iyr_iday(imon,syr,eyr,nday,iyr,iday)
           nb=nday/PETOT+1

           do ib=1,nb

              i=(my_rank+1)+(ib-1)*PETOT
              if(nday < i) cycle

              !---Read data
              call read_dat(varname,idat,iyr(i),imon,iday(i),im,jm,km,mask,mean3d,sprd3d)

              !--Add data
              call add(im,jm,km,mask,mean3d,mclim_mean3d,mpass_mean3d,mmiss_mean3d)
              call add(im,jm,km,mask,sprd3d,mclim_sprd3d,mpass_sprd3d,mmiss_sprd3d)
              call add(im,jm,km,mask,mean3d,clim_mean3d,pass_mean3d,miss_mean3d)
              call add(im,jm,km,mask,sprd3d,clim_sprd3d,pass_sprd3d,miss_sprd3d)
              
           end do !ib

           deallocate(iyr,iday)
           
           !---Monthly 
           call ave_mpi(im,jm,km,mclim_mean3d,mpass_mean3d,mmiss_mean3d)
           call ave_mpi(im,jm,km,mclim_sprd3d,mpass_sprd3d,mmiss_sprd3d)           
           
           !---Write data
           if(my_rank == master_rank)then
              call write_clim_nc(filename_mclim,varname,im,jm,km,imon,mclim_mean3d,mclim_sprd3d)
           end if
           call MPI_Barrier(MPI_COMM_WORLD,ierr)           

        end do !imon

        !---Annual Climatology
        call ave_mpi(im,jm,km,clim_mean3d,pass_mean3d,miss_mean3d)
        call ave_mpi(im,jm,km,clim_sprd3d,pass_sprd3d,miss_sprd3d)           
        
        !---Write data
        if(my_rank == master_rank)then
           call write_clim_nc(filename_clim,varname,im,jm,km,0,clim_mean3d,clim_sprd3d)
        end if
        call MPI_Barrier(MPI_COMM_WORLD,ierr)

     end do !ivar
        
     !---Deallocate
     deallocate(mean3d,sprd3d)
     deallocate(mclim_mean3d,mpass_mean3d,mmiss_mean3d)
     deallocate(mclim_sprd3d,mpass_sprd3d,mmiss_sprd3d)
     deallocate(clim_mean3d,pass_mean3d,miss_mean3d)
     deallocate(clim_sprd3d,pass_sprd3d,miss_sprd3d)

     !---Deallocate
     deallocate(mask,maskt,masku,maskv)

  end do !idat
  
  if(my_rank == master_rank)then
     write(*,*) "=== End monthly and annual climatologies ==="
  end if

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)
  
end program main

!-----------------------------------------------------------------------------------

subroutine get_total_day(imon,syr,eyr,nday)

  implicit none

  !---Common
  integer iyr
  integer mday
  
  !---IN
  integer,intent(in) :: imon
  integer,intent(in) :: syr,eyr

  !---OUT
  integer,intent(out) :: nday

  nday=0
  
  do iyr=syr,eyr

     call get_mday(iyr,imon,mday)
     nday=nday+mday
     
  end do

end subroutine get_total_day

!-----------------------------------------------------------------------------------

subroutine get_iyr_iday(imon,syr,eyr,nday,iyr,iday)

  implicit none

  !---Common
  integer n
  integer jyr,jday
  integer mday

  !---IN
  integer,intent(in) :: imon
  integer,intent(in) :: syr,eyr
  integer,intent(in) :: nday

  !---OUT
  integer,intent(out) :: iyr(nday),iday(nday)

  n=0
  
  do jyr=syr,eyr
     
     call get_mday(jyr,imon,mday)

     do jday=1,mday

        n=n+1
        iyr(n)=jyr
        iday(n)=jday

     end do
  end do
  
end subroutine get_iyr_iday
  
!-----------------------------------------------------------------------------------

subroutine get_mday(iyr,imon,mday)

  implicit none

  integer,intent(in) :: iyr,imon
  integer,intent(out) :: mday

  if(imon == 2 .and. mod(iyr,4) == 0)then
     mday=29
  else if(imon == 2)then
     mday=28
  else if(imon == 4 .or. imon == 6 .or. imon == 9 .or. imon == 11)then
     mday=30
  else
     mday=31
  end if

end subroutine get_mday

!-----------------------------------------------------------------------------------

subroutine init(im,jm,km,clim,pass,miss)

  implicit none

  !---IN
  integer,intent(in) :: im,jm,km

  !---INOUT
  real(kind = 8),intent(inout) :: clim(im,jm,km)
  integer,intent(inout) :: pass(im,jm,km),miss(im,jm,km)

  clim(:,:,:)=0.d0
  pass(:,:,:)=0
  miss(:,:,:)=0  

end subroutine init

!-----------------------------------------------------------------------------------

subroutine add(im,jm,km,mask,dat,clim,pass,miss)

  use mod_rmiss
  !$use omp_lib
  implicit none

  !---Common
  integer i,j,k

  !---IN
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: mask(im,jm),dat(im,jm,km)

  !---INOUT
  real(kind = 8),intent(inout) :: clim(im,jm,km)
  integer,intent(inout) :: pass(im,jm,km),miss(im,jm,km)
  
  !$omp parallel
  !$omp do private(i,j,k)  
  do k=1,km
     do j=1,jm
        do i=1,im

           if(mask(i,j) == 0.d0 .or. dat(i,j,k) == rmiss)then
              miss(i,j,k)=miss(i,j,k)+1
           else
              clim(i,j,k)=clim(i,j,k)+dat(i,j,k)
              pass(i,j,k)=pass(i,j,k)+1
           end if

        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine add

!-------------------------------------------------------------------------------------

subroutine ave_mpi(im,jm,km,clim,pass,miss)

  use mod_rmiss
  use mpi
  !$use omp_lib
  implicit none

  integer i,j,k
  integer ierr
  
  !---IN
  integer,intent(in) :: im,jm,km

  !---IN/OUT
  real(kind = 8),intent(inout) :: clim(im,jm,km)
  integer,intent(inout) :: pass(im,jm,km),miss(im,jm,km)

  call MPI_Allreduce(MPI_IN_PLACE, clim, im*jm*km, MPI_DOUBLE,  MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(MPI_IN_PLACE, pass, im*jm*km, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(MPI_IN_PLACE, miss, im*jm*km, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

  !$omp parallel
  !$omp do private(i,j,k)
  do k=1,km
     do j=1,jm
        do i=1,im
           
           if (pass(i,j,k) == 0) then
              clim(i,j,k)=rmiss
           else
              clim(i,j,k)=clim(i,j,k)/dble(pass(i,j,k))
           end if

        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine ave_mpi
