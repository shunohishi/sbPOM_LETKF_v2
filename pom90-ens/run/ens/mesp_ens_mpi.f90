!-------------------------------------------------------------------------
! Make Mean & Spread |
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! Created by S.Ohishi 2018.10
! Modified by S.Ohishi 2019.05: Add budget term
! Modified by S.Ohishi 2020.04: Add nuding term
! Modified by S.Ohishi 2020.11: Add surface ensemble
! Modified by S.Ohishi 2023.06: Add mpi, nf90
! Modiffed by S.Ohishi 2024.05: Add wr
! Modified by S.Ohishi 2025.01: Add argument part
!-------------------------------------------------------------------------
program main

  use julian_day
  use mpi
  implicit none

  !Common
  integer,parameter :: nvar2d=22,nvar=50
  integer,parameter :: ivar_ssu=23,ivar_ssv=24,ivar_sst=25,ivar_sss=26
    
  integer iyr,imon,iday,ijul
  integer syr,smon,sday
  integer im,jm,km
  integer imem,nmem
  integer it,nt
  integer iswitch_budget,iswitch_rm
  integer ivar

  character(100) dir,region
  character(4) yyyy
  character(2) mm,dd

  !MPI
  integer PETOT,my_rank,ierr
  integer :: master_rank=0

  !Variable
  real(kind = 4),allocatable :: time(:)
  real(kind = 4),allocatable :: z(:,:,:),zz(:,:,:)
  real(kind = 4),allocatable :: east_u(:,:),east_v(:,:),east_e(:,:)
  real(kind = 4),allocatable :: north_u(:,:),north_v(:,:),north_e(:,:)
  real(kind = 4),allocatable :: h(:,:),fsm(:,:),dum(:,:),dvm(:,:)
  real(kind = 4),allocatable :: dat2d(:,:),mean2d(:,:),sprd2d(:,:)
  real(kind = 4),allocatable :: dat3d(:,:,:),mean3d(:,:,:),sprd3d(:,:,:)

  character(10) var

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,PETOT,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
  imem=my_rank+1

  if(my_rank == master_rank) write(*,'(a)') "----- Start: mesp_ens_mpi -----"

  !Read ensemble run information
  call read_argument(dir,region,syr,smon,sday,iyr,imon,iday,nmem,nt,iswitch_budget,iswitch_rm)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday

  !Julian day
  call ymd_julian(iyr,imon,iday,ijul)

  !Read grid infomation (im,jm,km)
  if(my_rank == master_rank) write(*,'(a)') "Read ngrid"
  call read_ngrid(dir,region,iyr,imon,iday,im,jm,km)
  allocate(time(nt))
  allocate(z(im,jm,km),zz(im,jm,km))
  allocate(east_u(im,jm),east_v(im,jm),east_e(im,jm))
  allocate(north_u(im,jm),north_v(im,jm),north_e(im,jm))
  allocate(h(im,jm),fsm(im,jm),dum(im,jm),dvm(im,jm))
  allocate(dat2d(im,jm),mean2d(im,jm),sprd2d(im,jm))
  allocate(dat3d(im,jm,km),mean3d(im,jm,km),sprd3d(im,jm,km))
  
  !Read model information(time,...,dvm)
  if(my_rank == master_rank) write(*,'(a)') "Read info"
  call read_info(dir,region,iyr,imon,iday,nt,im,jm,km, &
       & time,z,zz,east_u,east_v,east_e,north_u,north_v,north_e, &
       & h,fsm,dum,dvm)
    
  do ivar=1,nvar
     
     !Read variable name
     call read_var(ivar,var)

     do it=1,nt
        
        call julian_ymd(ijul+(it-1),iyr,imon,iday)

        !Create NetCDF file
        if(ivar == 1  .and. my_rank == master_rank)then

           write(*,'(a)') "Create netcdf"
           call create_nc("mean",dir,region,syr,smon,sday,iyr,imon,iday,im,jm,km,imem,iswitch_budget)
           call create_nc("sprd",dir,region,syr,smon,sday,iyr,imon,iday,im,jm,km,imem,iswitch_budget)

           write(*,'(a)') "Write netcdf"
           call write_info("mean",dir,region,im,jm,km,iyr,imon,iday, &
                & time(it),z,zz,east_u,east_v,east_e,north_u,north_v,north_e, &
                & h,fsm,dum,dvm)
           call write_info("sprd",dir,region,im,jm,km,iyr,imon,iday, &
                & time(it),z,zz,east_u,east_v,east_e,north_u,north_v,north_e, &
                & h,fsm,dum,dvm)
           
        end if

        if(ivar == 1)then
           call create_nc("eens",dir,region,syr,smon,sday,iyr,imon,iday,im,jm,km,imem,iswitch_budget)
        end if
        call MPI_Barrier(MPI_COMM_WORLD,ierr)

        if(my_rank == master_rank) &
             & write(*,'(a,i3.3,a,i3.3,a,i3.3,a,i3.3)') "Variable:", ivar,"/",nvar," Day:",it,"/",nt
        
        !if(my_rank == master_rank) write(*,'(a)') "Ensemble Mean and Spread"
        !Mean & Spread
        if(1 <= ivar .and. ivar <= nvar2d)then
           !if(my_rank == master_rank) write(*,*) "Read ensemble"
           call read_ens(var,dir,region,imem,iyr,imon,iday,it,im,jm,1,dat2d)
           !if(my_rank == master_rank) write(*,*) "Calculate Ensemble mean and spread"
           call mpi_ensemble_mean_sprd(nmem,im,jm,1,dat2d,mean2d,sprd2d)           
        else if(nvar2d+1 <= ivar .and. ivar <= nvar)then
           ! write(*,*) "Read ensemble"
           call read_ens(var,dir,region,imem,iyr,imon,iday,it,im,jm,km,dat3d)
           ! write(*,*) "Calculate Ensemble mean & spread"
           call mpi_ensemble_mean_sprd(nmem,im,jm,km,dat3d,mean3d,sprd3d)
        end if

        !Write data: Ensemble mean and spread
        !write(*,*) "Write ensemble mean and spread"
        if(1 <= ivar .and. ivar <= nvar2d)then
           if(my_rank == master_rank) call write_dat("mean",var,dir,region,iyr,imon,iday,im,jm,1,mean2d)
           if(my_rank == nmem-1)      call write_dat("sprd",var,dir,region,iyr,imon,iday,im,jm,1,sprd2d)
        elseif(nvar2d+1 <= ivar .and. ivar <= nvar)then
           if(my_rank == master_rank) call write_dat("mean",var,dir,region,iyr,imon,iday,im,jm,km,mean3d)
           if(my_rank == nmem-1)      call write_dat("sprd",var,dir,region,iyr,imon,iday,im,jm,km,sprd3d)
        end if
        
        !Write data: Each ensemble (el, t, s, u, v)
        !if(my_rank == master_rank) write(*,*) "Write each ensemble"
        if(ivar == 1)then
           call write_ens("eens",var,dir,region,iyr,imon,iday,im,jm,imem,dat2d)
        else if(ivar == ivar_ssu .or. ivar == ivar_ssv .or. ivar == ivar_sst .or. ivar == ivar_sss)then
           call write_ens("eens",var,dir,region,iyr,imon,iday,im,jm,imem,dat3d(:,:,1))
        end if

        call MPI_Barrier(MPI_COMM_WORLD,ierr)

     end do !it
  end do !ivar

  deallocate(time)
  deallocate(z,zz)
  deallocate(east_u,east_v,east_e)
  deallocate(north_u,north_v,north_e)
  deallocate(h,fsm,dum,dvm)
  deallocate(dat2d,mean2d,sprd2d)
  deallocate(dat3d,mean3d,sprd3d)

  if(iswitch_rm == 1 .or. iswitch_rm == 2)then
     call remove_ens(dir,region,ijul,nt,imem)
  end if

  if(iswitch_rm == 2)then
     call remove_restart(dir,ijul,nt,imem)
  end if
  
  if(my_rank == master_rank) write(*,'(a)') "----- End: mesp_ens_mpi -----"

  call MPI_FINALIZE(ierr)

end program

!----------------------------------------------------------------------
! Read argument |
!----------------------------------------------------------------------

subroutine read_argument(dir,region,syr,smon,sday,iyr,imon,iday,nmem,nt,iswitch_budget,iswitch_rm)

  implicit none

  !---Common
  integer i,length,status

  character(:),allocatable :: arg
  
  intrinsic :: command_argument_count, get_command_argument
  
  !---Out
  character(100),intent(out) :: dir,region

  integer,intent(out) :: iyr,imon,iday
  integer,intent(out) :: syr,smon,sday
  integer,intent(out) :: nmem
  integer,intent(out) :: nt
  integer,intent(out) :: iswitch_budget,iswitch_rm


  do i=1,command_argument_count()

     call get_command_argument(i,length=length,status=status)

     if(status /= 0)then
        write(*,*) "Error: arugument ",status
     else

        allocate(character(length) :: arg)

        call get_command_argument(i,arg,status=status)

        if(i == 1)then
           dir=arg
        else if(i == 2)then
           region=arg
        else if(i == 3)then
           read(arg,'(I4)') syr
        else if(i == 4)then
           read(arg,'(I2)') smon
        else if(i == 5)then
           read(arg,'(I2)') sday
        else if(i == 6)then
           read(arg,'(I4)') iyr
        else if(i == 7)then
           read(arg,'(I2)') imon
        else if(i == 8)then
           read(arg,'(I2)') iday
        else if(i == 9)then
           read(arg,'(I5)') nmem
        else if(i == 10)then
           read(arg,'(I2)') nt
        else if(i == 11)then
           read(arg,'(I1)') iswitch_budget
        else if(i == 12)then
           read(arg,'(I1)') iswitch_rm
        end if
        
        deallocate(arg)

     end if

  end do
  
end subroutine read_argument

!----------------------------------------------------------------------
! Read variable name |
!----------------------------------------------------------------------

subroutine read_var(ivar,var)

  implicit none

  integer,intent(in) :: ivar
  character(10),intent(out) :: var

  !2D
  if(ivar == 1) var="el"
  if(ivar == 2) var="lhf"
  if(ivar == 3) var="shf"
  if(ivar == 4) var="lwr"
  if(ivar == 5) var="swr"
  if(ivar == 6) var="windu"
  if(ivar == 7) var="windv"
  if(ivar == 8) var="winds"
  if(ivar == 9) var="tauu"
  if(ivar == 10) var="tauv"
  if(ivar == 11) var="taus"
  if(ivar == 12) var="qa"
  if(ivar == 13) var="qs"
  if(ivar == 14) var="ta"
  if(ivar == 15) var="evap"
  if(ivar == 16) var="prep"
  if(ivar == 17) var="river"
  if(ivar == 18) var="eflux"
  if(ivar == 19) var="pflux"
  if(ivar == 20) var="rflux"
  if(ivar == 21) var="tsfc"
  if(ivar == 22) var="ssfc"

  !3D
  if(ivar == 23) var="u"
  if(ivar == 24) var="v"
  if(ivar == 25) var="t"
  if(ivar == 26) var="s"
  if(ivar == 27) var="w"
  if(ivar == 28) var="wr"

  if(ivar == 29) var="dtdt"
  if(ivar == 30) var="txadv"
  if(ivar == 31) var="tyadv"
  if(ivar == 32) var="tzadv"
  if(ivar == 33) var="txdif"
  if(ivar == 34) var="tydif"
  if(ivar == 35) var="tzdif"
  if(ivar == 36) var="qz"

  if(ivar == 37) var="dsdt"
  if(ivar == 38) var="sxadv"
  if(ivar == 39) var="syadv"
  if(ivar == 40) var="szadv"
  if(ivar == 41) var="sxdif"
  if(ivar == 42) var="sydif"
  if(ivar == 43) var="szdif"

  if(ivar == 44) var="aam"
  if(ivar == 45) var="kh"
  if(ivar == 46) var="km"

  if(ivar == 47) var="tnudge"
  if(ivar == 48) var="snudge"

  if(ivar == 49) var="tiau"
  if(ivar == 50) var="siau"

end subroutine read_var

!-----------------------------------------------------------------------
! Read Ensemble |
!-----------------------------------------------------------------------

subroutine read_ngrid(dir,region,iyr,imon,iday,im,jm,km)

  use netcdf
  implicit none

  !Common
  integer status,status1,status2,access
  integer ncid,dimid

  character(2) mm,dd
  character(4) yyyy
  character(200) filename,filename1,filename2

  !IN
  character(100),intent(in) :: dir,region
  integer,intent(in) :: iyr,imon,iday
  
  !OUT
  integer,intent(out) :: im,jm,km

  !Filename
  filename1=trim(dir)//"/00001/out/"//trim(region)//".nc"
  status1=access(trim(filename1)," ")
  
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday
  filename2=trim(dir)//"/00001/"//trim(region)//"."//yyyy//mm//dd//".nc"
  status2=access(trim(filename2)," ")
  
  if(status1 == 0)then
     filename=filename1
  else if(status2 == 0)then
     filename=filename2
  else
     write(*,*) "***Error: Not Found "//trim(filename1)//" and "//trim(filename2)
     stop
  end if

  !Get im,jm,km
  status=nf90_open(trim(filename),nf90_nowrite,ncid)
  call check_error(status)

  status=nf90_inq_dimid(ncid,"x",dimid)
  call check_error(status)
  status=nf90_inquire_dimension(ncid,dimid,len = im)
  call check_error(status)

  status=nf90_inq_dimid(ncid,"y",dimid)
  call check_error(status)
  status=nf90_inquire_dimension(ncid,dimid,len = jm)
  call check_error(status)

  status=nf90_inq_dimid(ncid,"z",dimid)
  call check_error(status)
  status=nf90_inquire_dimension(ncid,dimid,len = km)
  call check_error(status)

  status=nf90_close(ncid)
  call check_error(status)

end subroutine read_ngrid

!-----------------------------------------

subroutine read_info(dir,region,iyr,imon,iday,nt,im,jm,km, &
     & time,z,zz,east_u,east_v,east_e,north_u,north_v,north_e, &
     & h,fsm,dum,dvm)

  use julian_day  
  use netcdf
  implicit none

  !Common
  integer status,status1,status2,access
  integer ncid,varid
  integer it,ijul
  
  character(200) filename,filename1,filename2
  character(4) yyyy
  character(2) mm,dd

  !IN
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: nt
  integer,intent(in) :: im,jm,km

  character(100),intent(in) :: dir,region

  !OUT
  real(kind = 4),intent(out) :: time(nt)
  real(kind = 4),intent(out) :: z(im,jm,km),zz(im,jm,km)
  real(kind = 4),intent(out) :: east_u(im,jm),east_v(im,jm),east_e(im,jm)
  real(kind = 4),intent(out) :: north_u(im,jm),north_v(im,jm),north_e(im,jm)
  real(kind = 4),intent(out) :: h(im,jm),fsm(im,jm),dum(im,jm),dvm(im,jm)

  !Filename
  filename1=trim(dir)//"/00001/out/"//trim(region)//".nc"
  status1=access(trim(filename1)," ")
  
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday
  filename2=trim(dir)//"/00001/"//trim(region)//"."//yyyy//mm//dd//".nc"
  status2=access(trim(filename2)," ")

  if(status1 == 0)then
     filename=filename1
  else if(status2 == 0)then
     filename=filename2
  else
     write(*,*) "***Error: Not Found "//trim(filename1)//" and "//trim(filename2)
     stop
  end if

  !Get time
  if(status1 == 0)then
  
     status=nf90_open(trim(filename),nf90_nowrite,ncid)
     call check_error(status)

     status=nf90_inq_varid(ncid,"time",varid)
     call check_error(status)
     status=nf90_get_var(ncid,varid,time)
     call check_error(status)

     status=nf90_close(ncid)
     call check_error(status)
     
  else if(status2 == 0)then

     call ymd_julian(iyr,imon,iday,ijul)
     
     do it=1,nt

        call julian_ymd(ijul+it-1,iyr,imon,iday)

        write(yyyy,'(i4.4)') iyr
        write(mm,'(i2.2)') imon
        write(dd,'(i2.2)') iday
        filename=trim(dir)//"/00001/"//trim(region)//"."//yyyy//mm//dd//".nc"

        status=nf90_open(trim(filename),nf90_nowrite,ncid)
        call check_error(status)

        status=nf90_inq_varid(ncid,"time",varid)
        call check_error(status)
        status=nf90_get_var(ncid,varid,time(it))
        call check_error(status)

        status=nf90_close(ncid)
        call check_error(status)
        
     end do
     
  end if
     
  status=nf90_open(trim(filename),nf90_nowrite,ncid)
  call check_error(status)   

  status=nf90_inq_varid(ncid,"z_w",varid)
  call check_error(status)
  status=nf90_get_var(ncid,varid,z)
  call check_error(status)

  status=nf90_inq_varid(ncid,"z_e",varid)
  call check_error(status)
  status=nf90_get_var(ncid,varid,zz)
  call check_error(status)

  status=nf90_inq_varid(ncid,"east_u",varid)
  call check_error(status)
  status=nf90_get_var(ncid,varid,east_u)
  call check_error(status)

  status=nf90_inq_varid(ncid,"east_v",varid)
  call check_error(status)
  status=nf90_get_var(ncid,varid,east_v)
  call check_error(status)

  status=nf90_inq_varid(ncid,"east_e",varid)
  call check_error(status)
  status=nf90_get_var(ncid,varid,east_e)
  call check_error(status)

  status=nf90_inq_varid(ncid,"north_u",varid)
  call check_error(status)
  status=nf90_get_var(ncid,varid,north_u)
  call check_error(status)

  status=nf90_inq_varid(ncid,"north_v",varid)
  call check_error(status)
  status=nf90_get_var(ncid,varid,north_v)
  call check_error(status)

  status=nf90_inq_varid(ncid,"north_e",varid)
  call check_error(status)
  status=nf90_get_var(ncid,varid,north_e)
  call check_error(status)

  status=nf90_inq_varid(ncid,"h",varid)
  call check_error(status)
  status=nf90_get_var(ncid,varid,h)
  call check_error(status)

  status=nf90_inq_varid(ncid,"fsm",varid)
  call check_error(status)
  status=nf90_get_var(ncid,varid,fsm)
  call check_error(status)

  status=nf90_inq_varid(ncid,"dum",varid)
  call check_error(status)
  status=nf90_get_var(ncid,varid,dum)
  call check_error(status)

  status=nf90_inq_varid(ncid,"dvm",varid)
  call check_error(status)
  status=nf90_get_var(ncid,varid,dvm)
  call check_error(status)

  status=nf90_close(ncid)
  call check_error(status)

end subroutine read_info

!------------------------------------------------------------

subroutine read_ens(var,dir,region,imem,iyr,imon,iday,it,im,jm,km,dat)

  use netcdf
  implicit none

  !Common
  integer status,status1,status2,access
  integer ncid,varid

  character(2) mm,dd
  character(4) yyyy
  character(5) mmmmm
  character(200) filename,filename1,filename2

  !IN
  integer,intent(in) :: imem
  integer,intent(in) :: iyr,imon,iday,it
  integer,intent(in) :: im,jm,km
  
  character(10),intent(in) :: var
  character(100),intent(in) :: dir,region

  !OUT
  real(kind = 4),intent(out) :: dat(im,jm,km)

  !Filename
  write(mmmmm,'(i5.5)') imem
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday

  filename1=trim(dir)//"/"//mmmmm//"/out/"//trim(region)//".nc"
  status1=access(trim(filename1)," ")

  filename2=trim(dir)//"/"//mmmmm//"/"//trim(region)//"."//yyyy//mm//dd//".nc"
  status2=access(trim(filename2)," ")

  if(status1 == 0)then
     filename=filename1
  else if(status2 == 0)then
     filename=filename2
  else     
     write(*,*) "***Error: Not Found "//trim(filename1)//" and "//trim(filename2)
     stop
  end if

  !Get data
  status=nf90_open(trim(filename),nf90_nowrite,ncid)
  call check_error(status)

  status=nf90_inq_varid(ncid,trim(var),varid)
  if(status == nf90_noerr)then
     
     if(status1 == 0 .and. km == 1)then
        status=nf90_get_var(ncid,varid,dat,(/1,1,it/),(/im,jm,1/))
     else if(status1 == 0)then
        status=nf90_get_var(ncid,varid,dat,(/1,1,1,it/),(/im,jm,km,1/))
     else if(status2 == 0 .and. km == 1)then
        status=nf90_get_var(ncid,varid,dat,(/1,1,1/),(/im,jm,1/))
     else if(status2 == 0)then
        status=nf90_get_var(ncid,varid,dat,(/1,1,1,1/),(/im,jm,km,1/))
     end if
     
  else
     dat(:,:,:)=0.e0
  end if

  status=nf90_close(ncid)
  call check_error(status)

end subroutine read_ens

!----------------------------------------------------------------------
!  Calculate Ensemble Mean |
!----------------------------------------------------------------------

subroutine mpi_ensemble_mean_sprd(nmem,im,jm,km,dat,mean,sprd)

  use mpi
  implicit none

  integer i,j,k
  integer inum,nnum
  integer ierr

  real send(im*jm*km)
  real rec(im*jm*km)

  !IN
  integer,intent(in) :: nmem
  integer,intent(in) :: im,jm,km
  real(kind = 4),intent(in) :: dat(im,jm,km)

  !OUT
  real(kind = 4),intent(out) :: mean(im,jm,km),sprd(im,jm,km)

  !nnum
  nnum=im*jm*km

  !Ensemble mean: send --> rec --> mean
  inum=0
  do k=1,km
     do j=1,jm
        do i=1,im
           inum=inum+1
           send(inum)=dat(i,j,k)
        end do
     end do
  end do

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(send(1),rec(1),nnum,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

  do inum=1,nnum
     rec(inum)=rec(inum)/real(nmem)
  end do

  inum=0
  do k=1,km
     do j=1,jm
        do i=1,im
           inum=inum+1
           mean(i,j,k)=rec(inum)
        end do
     end do
  end do

  !Ensemble sprd: send --> rec --> sprd
  inum=0
  do k=1,km
     do j=1,jm
        do i=1,im
           inum=inum+1
           send(inum)=(dat(i,j,k)-mean(i,j,k))*(dat(i,j,k)-mean(i,j,k))
        end do
     end do
  end do

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(send(1),rec(1),nnum,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

  do inum=1,nnum
     rec(inum)=sqrt(rec(inum)/real(nmem-1))
  end do

  inum=0
  do k=1,km
     do j=1,jm
        do i=1,im
           inum=inum+1
           sprd(i,j,k)=rec(inum)
        end do
     end do
  end do
  
end subroutine mpi_ensemble_mean_sprd

!----------------------------------------------------------------
! Check netcdf error |
!----------------------------------------------------------------

subroutine check_error(status)

  use mpi
  use netcdf
  implicit none

  integer ierr
  integer,intent(in) :: status
  
  if(status /= nf90_noerr) then
    write(*,'(a)') "***Error:"//trim(nf90_strerror(status))
    call MPI_FINALIZE(ierr)
    stop
  end if

end subroutine check_error

!----------------------------------------------------------------------
! Create Netcdf |
!----------------------------------------------------------------------

subroutine define_var_netcdf(ncid,ndim,dim,varid,name,long_name,units_name)

  use netcdf
  implicit none

  integer status

  integer,intent(in) :: ncid
  integer,intent(in) :: ndim,dim(ndim)
  integer,intent(inout) :: varid

  character(*),intent(in) :: name,long_name,units_name

  status=nf90_def_var(ncid,trim(name),nf90_float,dim,varid)
  call check_error(status)

  status=nf90_def_var_deflate(ncid,varid,shuffle=1,deflate=1,deflate_level=5)
  call check_error(status)
  
  status=nf90_put_att(ncid,varid,"long_name",trim(long_name))
  call check_error(status)
  
  status=nf90_put_att(ncid,varid,"units",trim(units_name))
  call check_error(status)

end subroutine define_var_netcdf

!--------------------------------

subroutine create_nc(mesp,dir,region,syr,smon,sday,iyr,imon,iday,im,jm,km,imem,iswitch_budget)
  
  use netcdf
  implicit none

  !Common
  integer status,status1,status2,access
  integer ncid,varid
  integer x_dimid,y_dimid,z_dimid,time_dimid

  integer ndim
  integer,allocatable :: dim(:)

  integer iswitch_ff     !0: No, 1: freshwater flux
  integer iswitch_tnudge !0: No, 1: T nudging
  integer iswitch_snudge !0: No, 1: S nudging
  integer iswitch_tiau   !0: No, 1: T IAU
  integer iswitch_siau   !0: No, 1: S IAU

  character(200) dirname
  character(200) filename,filename1,filename2
  character(5) mmmmm
  character(4) mesp,yyyy
  character(2) mm,dd

  !IN
  integer,intent(in) :: syr,smon,sday
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: im,jm,km
  integer,intent(in) :: imem
  integer,intent(in) :: iswitch_budget !0:No, 1: Budget term

  character(100),intent(in) :: dir,region

  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday

  !Check freshwater flux & Nudging
  iswitch_ff=0
  iswitch_tnudge=0
  iswitch_snudge=0
  iswitch_tiau=0
  iswitch_siau=0

  !Filename
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday

  !Free ensemble simulation
  filename1=trim(dir)//"/00001/out/"//trim(region)//".nc"
  status1=access(trim(filename1)," ")

  !Ensemble simulation in LETKF
  filename2=trim(dir)//"/00001/"//trim(region)//"."//yyyy//mm//dd//".nc"
  status2=access(trim(filename2)," ")

  if(status1 == 0)then
     filename=filename1
  else if(status2 == 0)then
     filename=filename2
  else     
     write(*,*) "***Error: Not Found "//trim(filename1)//" and "//trim(filename2)
     stop
  end if
  
  status=nf90_open(trim(filename),nf90_nowrite,ncid)
  call check_error(status)

  status=nf90_inq_varid(ncid,"evap",varid)
  if(status == nf90_noerr) iswitch_ff=1
  status=nf90_inq_varid(ncid,"tnudge",varid)
  if(status == nf90_noerr) iswitch_tnudge=1
  status=nf90_inq_varid(ncid,"snudge",varid)
  if(status == nf90_noerr) iswitch_snudge=1
  status=nf90_inq_varid(ncid,"tiau",varid)
  if(status == nf90_noerr) iswitch_tiau=1
  status=nf90_inq_varid(ncid,"siau",varid)
  if(status == nf90_noerr) iswitch_siau=1

  status=nf90_close(ncid)
  call check_error(status)

  !Write filename
  dirname=trim(dir)//"/"//trim(mesp)

  status=access(trim(dirname)," ")
  if(status /= 0)then
     write(*,*) "***Error: Not Found "//trim(dirname)
     stop
  end if
  
  if(mesp == "mean" .or. mesp == "sprd")then
     filename=trim(dirname)//"/"//trim(region)//yyyy//mm//dd//".nc"
  else if(mesp == "eens")then
     write(mmmmm,'(i5.5)') imem
     filename=trim(dirname)//"/"//trim(region)//yyyy//mm//dd//"."//mmmmm//".nc"
  end if

  !NF_CREATE
  status=access(trim(filename)," ")
  if(status == 0)then
     return
  endif
  
  status=nf90_create(trim(filename),nf90_netcdf4,ncid)
  call check_error(status)

  !GLOBAL ATTRIBUTE
  status=nf90_put_att(ncid,nf90_global,"title",trim(region))
  call check_error(status)
  status=nf90_put_att(ncid,nf90_global,"description","output file")
  call check_error(status)

  !DIMID
  status=nf90_def_dim(ncid,"time",1,time_dimid)
  call check_error(status)
  status=nf90_def_dim(ncid,"x",im,x_dimid)
  call check_error(status)
  status=nf90_def_dim(ncid,"y",jm,y_dimid)
  call check_error(status)
  status=nf90_def_dim(ncid,"z",km,z_dimid)

  if(mesp == "mean" .or. mesp == "sprd")then

     !---1D
     ndim=1
     allocate(dim(ndim))
     
     !TIME
     dim(1)=time_dimid
     
     write(yyyy,'(i4.4)') syr
     write(mm,'(i2.2)') smon
     write(dd,'(i2.2)') sday
     call define_var_netcdf(ncid,ndim,dim,varid,"time","time", &
          & "days since "//yyyy//"-"//mm//"-"//dd//" 00:00:00 +00:00")

     deallocate(dim)
     
     !---2D
     ndim=2
     allocate(dim(ndim))
     dim(1)=x_dimid
     dim(2)=y_dimid
     
     call define_var_netcdf(ncid,ndim,dim,varid,"east_u","longitude of u points","degree E")
     call define_var_netcdf(ncid,ndim,dim,varid,"east_v","longitude of v points","degree E")
     call define_var_netcdf(ncid,ndim,dim,varid,"east_e","longitude of elevation points","degree E")
     
     call define_var_netcdf(ncid,ndim,dim,varid,"north_u","latitude of u points","degree N")
     call define_var_netcdf(ncid,ndim,dim,varid,"north_v","latitude of v points","degree N")
     call define_var_netcdf(ncid,ndim,dim,varid,"north_e","latitude of elevation points","degree N")
     
     call define_var_netcdf(ncid,ndim,dim,varid,"h","undisturbed water depth","meter")
     
     call define_var_netcdf(ncid,ndim,dim,varid,"dum","u surface mask [u]","-")
     call define_var_netcdf(ncid,ndim,dim,varid,"dvm","v surface mask [v]","-")
     call define_var_netcdf(ncid,ndim,dim,varid,"fsm","elevation surface mask [el]","-")
     
     deallocate(dim)

     !---3D
     ndim=3
     allocate(dim(ndim))

     dim(1)=x_dimid
     dim(2)=y_dimid
     dim(3)=z_dimid  
     !Z
     call define_var_netcdf(ncid,ndim,dim,varid,"z_w","sigma of cell face","sigma_level")
     !ZZ
     call define_var_netcdf(ncid,ndim,dim,varid,"z_e","sigma of cell center","sigma_level")
     
     !---2D+time
     dim(1)=x_dimid
     dim(2)=y_dimid
     dim(3)=time_dimid
     
     call define_var_netcdf(ncid,ndim,dim,varid,"el","surface elevation [el]","meter")
     call define_var_netcdf(ncid,ndim,dim,varid,"lhf","latent heat flux [el]","W/m^2")
     call define_var_netcdf(ncid,ndim,dim,varid,"shf","sensible heat flux [el]","W/m^2")
     call define_var_netcdf(ncid,ndim,dim,varid,"lwr","longwave radiation [el]","W/m^2")
     call define_var_netcdf(ncid,ndim,dim,varid,"swr","shortwave radiation [el]","W/m^2")
     
     call define_var_netcdf(ncid,ndim,dim,varid,"windu","zonal wind [el]","m/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"windv","meridional wind [el]","m/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"winds","wind speed [el]","m/s")
     
     call define_var_netcdf(ncid,ndim,dim,varid,"tauu","zonal wind stress [el]","N/m^2")
     call define_var_netcdf(ncid,ndim,dim,varid,"tauv","meridonal wind stress [el]","N/m^2")
     call define_var_netcdf(ncid,ndim,dim,varid,"taus","magnitude of wind stress [el]","N/m^2")
     
     call define_var_netcdf(ncid,ndim,dim,varid,"qa","air specific humidity [el]","g/kg")
     call define_var_netcdf(ncid,ndim,dim,varid,"qs","surface saturated specific humidity [el]","g/kg")
     call define_var_netcdf(ncid,ndim,dim,varid,"ta","air temperature [el]","degree C")
     
     if(iswitch_ff == 1)then
        call define_var_netcdf(ncid,ndim,dim,varid,"evap","evaporation [el]","mm/day")
        call define_var_netcdf(ncid,ndim,dim,varid,"prep","precipitation [el]","mm/day")
        call define_var_netcdf(ncid,ndim,dim,varid,"river","river discharge [el]","mm/day")
        
        call define_var_netcdf(ncid,ndim,dim,varid,"eflux","evaporation flux [el]","m/s")
        call define_var_netcdf(ncid,ndim,dim,varid,"pflux","precipitation flux [el]","m/s")
        call define_var_netcdf(ncid,ndim,dim,varid,"rflux","river discharge flux [el]","m/s")        
     end if
     
     if(iswitch_budget == 1)then
        call define_var_netcdf(ncid,ndim,dim,varid,"tsfc","surface heat flux forcing [el]","degree C/day")
        call define_var_netcdf(ncid,ndim,dim,varid,"ssfc","surface freshwater flux forcing [el]","1/day")
     end if
     
     deallocate(dim)
     
     !3D+time
     ndim=4
     allocate(dim(ndim))
     dim(1)=x_dimid
     dim(2)=y_dimid
     dim(3)=z_dimid
     dim(4)=time_dimid
     
     call define_var_netcdf(ncid,ndim,dim,varid,"u","zonal velocity [u,zz]","m/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"v","meridional velocity [v,zz]","m/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"w","vertical velocity (sigma-cordinate) [el,z]","m/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"wr","vertical velocity (z-cordinate) [el,z]","m/s")
     
     call define_var_netcdf(ncid,ndim,dim,varid,"t","potential temperature [el,zz]","degree C")
     call define_var_netcdf(ncid,ndim,dim,varid,"s","salinity [el,zz]","-")
     
     if(iswitch_budget == 1)then
        
        call define_var_netcdf(ncid,ndim,dim,varid,"dtdt","temperature tendency [el,zz]","degree C/day")
        
        call define_var_netcdf(ncid,ndim,dim,varid,"txadv","zonal temperature advection [el,zz]","degree C/day")
        call define_var_netcdf(ncid,ndim,dim,varid,"tyadv","meridional temperature advection [el,zz]","degree C/day")
        call define_var_netcdf(ncid,ndim,dim,varid,"tzadv","vertical temperature advection [el,zz]","degree C/day")
        
        call define_var_netcdf(ncid,ndim,dim,varid,"txdif","zonal temperature diffusion [el,zz]","degree C/day")
        call define_var_netcdf(ncid,ndim,dim,varid,"tydif","meridional temperature diffusion [el,zz]","degree C/day")
        call define_var_netcdf(ncid,ndim,dim,varid,"tzdif","vertical temperature diffusion [el,zz]","degree C/day")
        
        call define_var_netcdf(ncid,ndim,dim,varid,"qz","shortwave penetration [el,zz]","degree C/day")
        
        call define_var_netcdf(ncid,ndim,dim,varid,"dsdt","salinity tendency [el,zz]","1/day")
        
        call define_var_netcdf(ncid,ndim,dim,varid,"sxadv","zonal salinity advection [el,zz]","1/day")
        call define_var_netcdf(ncid,ndim,dim,varid,"syadv","meridional salinity advection [el,zz]","1/day")
        call define_var_netcdf(ncid,ndim,dim,varid,"szadv","vertical salinity advection [el,zz]","1/day")
        
        
        call define_var_netcdf(ncid,ndim,dim,varid,"sxdif","zonal salinity diffusion [el,zz]","1/day")
        call define_var_netcdf(ncid,ndim,dim,varid,"sydif","meridional salinity diffusion [el,zz]","1/day")
        call define_var_netcdf(ncid,ndim,dim,varid,"szdif","vertical salinity diffusion [el,zz]","1/day")
        
        call define_var_netcdf(ncid,ndim,dim,varid,"aam","Am (Ah=0.2*Am) [el,zz]","m^2/s")
        call define_var_netcdf(ncid,ndim,dim,varid,"km","vertical kinematic viscosity [el,z]","m^2/s")
        call define_var_netcdf(ncid,ndim,dim,varid,"kh","vertical diffusivity [el,z]","m^2/s")
        
        if(iswitch_tnudge == 1)then
           call define_var_netcdf(ncid,ndim,dim,varid,"tnudge","temperature nudging [el,zz]","degree C/day")
        end if
        
        if(iswitch_snudge == 1)then
           call define_var_netcdf(ncid,ndim,dim,varid,"snudge","salinity nudging [el,zz]","1/day")
        end if
        
        if(iswitch_tiau == 1)then
           call define_var_netcdf(ncid,ndim,dim,varid,"tiau","temperature analysis increment [el]","degree C/day")        
        end if
        
        if(iswitch_siau == 1)then
           call define_var_netcdf(ncid,ndim,dim,varid,"siau","salinity analysis increment [el]","1/day")
        end if
        
     end if !budget

     deallocate(dim)

  else if(mesp == "eens")then
     
     ndim=3
     allocate(dim(ndim))
     
     dim(1)=x_dimid
     dim(2)=y_dimid
     dim(3)=time_dimid
     
     call define_var_netcdf(ncid,ndim,dim,varid,"el","surface elevation [el]","meter")
     call define_var_netcdf(ncid,ndim,dim,varid,"u","zonal velocity [u,zz]","m/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"v","meridional velocity [v,zz]","m/s")     
     call define_var_netcdf(ncid,ndim,dim,varid,"t","potential temperature [el,zz]","degree C")
     call define_var_netcdf(ncid,ndim,dim,varid,"s","salinity [el,zz]","-")
     
     deallocate(dim)

  end if !mean/sprd
     
  !END Definition
  status=nf90_enddef(ncid)
  call check_error(status)
  status=nf90_close(ncid)
  call check_error(status)
  
end subroutine create_nc

!----------------------------------------------------------------------
! Write Netcdf |
!----------------------------------------------------------------------

subroutine write_info(mesp,dir,region,im,jm,km,iyr,imon,iday, &
     & time,z,zz,east_u,east_v,east_e,north_u,north_v,north_e, &
     & h,fsm,dum,dvm)

  use netcdf
  implicit none

  !Common
  integer status,access
  integer ncid,varid
 
  character(200) filename
  character(4) yyyy
  character(2) mm,dd

  !IN
  integer,intent(in) :: im,jm,km
  integer,intent(in) :: iyr,imon,iday

  real(kind = 4),intent(in) :: time
  real(kind = 4),intent(in) :: z(im,jm,km),zz(im,jm,km)
  real(kind = 4),intent(in) :: east_u(im,jm),east_v(im,jm),east_e(im,jm)
  real(kind = 4),intent(in) :: north_u(im,jm),north_v(im,jm),north_e(im,jm)
  real(kind = 4),intent(in) :: h(im,jm),fsm(im,jm),dum(im,jm),dvm(im,jm)

  character(4),intent(in) :: mesp
  character(100),intent(in) :: dir,region

  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday
  
  filename=trim(dir)//"/"//trim(mesp)//"/"//trim(region)//yyyy//mm//dd//".nc"
  status=access(trim(filename)," ")

  if(status /= 0)then
     write(*,*) "***Error: "//trim(filename)//" not found"
     stop
  end if
  
  status=nf90_open(trim(filename),nf90_write,ncid)
  call check_error(status)

  status=nf90_inq_varid(ncid,"time",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,time)
  call check_error(status)

  status=nf90_inq_varid(ncid,"z_w",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,z)
  call check_error(status)

  status=nf90_inq_varid(ncid,"z_e",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,zz)
  call check_error(status)

  status=nf90_inq_varid(ncid,"east_u",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,east_u)
  call check_error(status)

  status=nf90_inq_varid(ncid,"east_v",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,east_v)
  call check_error(status)

  status=nf90_inq_varid(ncid,"east_e",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,east_e)
  call check_error(status)

  status=nf90_inq_varid(ncid,"north_u",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,north_u)
  call check_error(status)

  status=nf90_inq_varid(ncid,"north_v",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,north_v)
  call check_error(status)

  status=nf90_inq_varid(ncid,"north_e",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,north_e)
  call check_error(status)

  status=nf90_inq_varid(ncid,"h",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,h)
  call check_error(status)

  status=nf90_inq_varid(ncid,"fsm",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,fsm)
  call check_error(status)

  status=nf90_inq_varid(ncid,"dum",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,dum)
  call check_error(status)

  status=nf90_inq_varid(ncid,"dvm",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,dvm)
  call check_error(status)

  status=nf90_close(ncid)
  call check_error(status)

end subroutine write_info

!-------------------------------------------
subroutine write_dat(mesp,var,dir,region,iyr,imon,iday,im,jm,km,dat)

  use netcdf
  implicit none

  !Common
  integer status,access
  integer ncid,varid
  integer ierr

  character(2) mm,dd
  character(4) yyyy
  character(200) filename

  !IN
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: im,jm,km

  real(kind = 4),intent(in) :: dat(im,jm,km)

  character(4),intent(in) :: mesp
  character(10),intent(in) :: var
  character(100),intent(in) :: dir,region

  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday

  filename=trim(dir)//"/"//trim(mesp)//"/"//trim(region)//yyyy//mm//dd//".nc"
  status=access(trim(filename)," ")

  if(status /= 0)then
     write(*,*) "***Error: "//trim(filename)//" not found"
     call MPI_FINALIZE(ierr)
     stop
  end if  
  
  status=nf90_open(trim(filename),nf90_write,ncid)
  call check_error(status)

  status=nf90_inq_varid(ncid,trim(var),varid)
  if(status == nf90_noerr)then
     if(km == 1)then
        status=nf90_put_var(ncid,varid,dat,(/1,1,1/),(/im,jm,1/))
        call check_error(status)
     else
        status=nf90_put_var(ncid,varid,dat,(/1,1,1,1/),(/im,jm,km,1/))
        call check_error(status)
     end if
  end if

  status=nf90_close(ncid)
  call check_error(status)

end subroutine write_dat

!-----------------------------------------------
subroutine write_ens(mesp,var,dir,region,iyr,imon,iday,im,jm,imem,dat)

  use mpi
  use netcdf
  implicit none

  !Common
  integer status,access
  integer ncid,varid
  integer ierr

  character(2) mm,dd
  character(4) yyyy
  character(5) mmmmm
  character(200) filename

  !IN
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: im,jm,imem

  real(kind = 4),intent(in) :: dat(im,jm)

  character(4),intent(in) :: mesp
  character(10),intent(in) :: var
  character(100),intent(in) :: dir,region

  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday

  write(mmmmm,'(i5.5)') imem

  filename=trim(dir)//"/"//trim(mesp)//"/"//trim(region)//yyyy//mm//dd//"."//mmmmm//".nc"
  status=access(trim(filename)," ")
  
  if(status /= 0)then
     write(*,*) "***Error: "//trim(filename)//" not found"
     call MPI_FINALIZE(ierr)
     stop
  end if
  
  status=nf90_open(trim(filename),nf90_write,ncid)
  call check_error(status)

  status=nf90_inq_varid(ncid,trim(var),varid)
  if(status == nf90_noerr)then
     status=nf90_put_var(ncid,varid,dat,(/1,1,1/),(/im,jm,1/))
  end if

  status=nf90_close(ncid)
  call check_error(status)

end subroutine write_ens

!---------------------------------------------

subroutine remove_ens(dir,region,ijul,nt,imem)

  use julian_day
  implicit none
  
  !Common
  integer status,status1,status2,system,access
  integer iyr,imon,iday
  integer it
  
  character(5) mmmmm
  character(200) filename1,filename2

  !IN
  integer,intent(in) :: ijul,nt
  integer,intent(in) :: imem
  
  character(4) yyyy
  character(2) mm,dd
  character(100),intent(in) :: dir,region

  call julian_ymd(ijul,iyr,imon,iday)  
  write(mmmmm,'(i5.5)') imem
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday

  !Free ensemble simulation
  filename1=trim(dir)//"/"//mmmmm//"/out/"//trim(region)//".nc"
  status1=access(trim(filename1)," ")

  !Ensemble simulation in LETKF
  filename2=trim(dir)//"/"//mmmmm//"/"//trim(region)//"."//yyyy//mm//dd//".nc"
  status2=access(trim(filename2)," ")

  if(status1 == 0)then

     status=system("rm -f "//trim(filename1))
     
  else if(status2 == 0)then
     
     do it=1,nt

        call julian_ymd(ijul+(it-1),iyr,imon,iday)
        write(yyyy,'(i4.4)') iyr
        write(mm,'(i2.2)') imon
        write(dd,'(i2.2)') iday
        filename2=trim(dir)//"/"//mmmmm//"/"//trim(region)//"."//yyyy//mm//dd//".nc"

        status=system("rm -f "//trim(filename2))
        
     end do
     
  else
     
     write(*,*) "***Error: Not found "//trim(filename1)//" and "//trim(filename2)
     stop

  end if
  
  
end subroutine remove_ens

!--------------------------------------------

subroutine remove_restart(dir,ijul,nt,imem)

  use julian_day
  implicit none
  
  !Common
  integer status,status1,status2,system,access
  integer iyr,imon,iday
  integer it
  
  character(5) mmmmm
  character(200) filename1,filename2

  !IN
  integer,intent(in) :: ijul,nt
  integer,intent(in) :: imem

  character(4) yyyy
  character(2) mm,dd
  character(100),intent(in) :: dir

  !Previous date
  call julian_ymd(ijul-1,iyr,imon,iday)
  write(mmmmm,'(i5.5)') imem
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday

  !Free ensemble simulation
  filename1=trim(dir)//"/"//mmmmm//"/out/restart.nc"
  status1=access(trim(filename1)," ")

  !Ensemble simulation in LETKF
  filename2=trim(dir)//"/"//mmmmm//"/restart."//yyyy//mm//dd//".nc"
  status2=access(trim(filename2)," ")

  if(status1 == 0)then
     !status=system("rm -f "//trim(filename1))
  else if(status2 == 0)then

     do it=1,nt

        !Current date
        call julian_ymd(ijul+it-1,iyr,imon,iday)
        if(iday == 1) cycle

        !Previous date
        call julian_ymd(ijul+it-2,iyr,imon,iday)
        write(yyyy,'(i4.4)') iyr
        write(mm,'(i2.2)') imon
        write(dd,'(i2.2)') iday
        filename2=trim(dir)//"/"//mmmmm//"/restart."//yyyy//mm//dd//".nc"

        status=system("rm -f "//trim(filename2))        

     end do
        
  else
     write(*,*) "***Error: Not found "//trim(filename1)//" and "//trim(filename2)
     stop
  end if

  
end subroutine remove_restart
