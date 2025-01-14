module setting
  
  integer,parameter :: syr_atm=1979,eyr_atm=2022 !time range
  real(kind = 8),parameter :: factor=0.2 !raw data + factor * standard deviation * normalized random number  
  
end module setting

!________________________________________________________________

program main
  
  !----------------------------------------------------------------
  ! Make Ensemble atmospheric netcdf file |
  !----------------------------------------------------------------
  !
  ! 1. Calculate the number of caluculating times
  ! 2. Read model grid
  ! 3. Read JRA55 grid & land
  ! 4. Calculate ID
  ! 5. Set ensemble year
  ! 6. Read JRA55
  ! 7. Calculate ensemble JRA55
  ! 8. Apply JRA55 land mask
  ! 9. Fill value
  ! 10. Bilinear interpolation
  ! 11. Write Data
  !
  ! *Initial hour = 0 at iday = 1
  !
  ! Created by S.Ohishi @ 2019.10
  !
  ! Added openmp by S.Ohishi @ 2020.04 (*netcdf is not thread safe)
  ! Added JRA55do by S.Ohishi @ 2021.11
  ! Removed openmp & Added mpi by S.Ohishi@ 2021.11
  ! Modified all by S.Ohishi @ 2022.04
  !
  !----------------------------------------------------------------

  use setting
  use mod_rmiss
  use mod_gridinfo
  use mod_read_jra55do, im_jra55do => im, jm_jra55do => jm, read_grid_jra55do => read_grid 
  use mpi
  implicit none

  !Common
  integer itime,ntime
  integer,allocatable :: iyr(:),imon(:),iday(:),ihour(:)
  integer syr,smon,sday,shour
  integer eyr,emon,eday,ehour
  integer iens,nens
  integer i,j
  integer status,system
  
  character(4) yyyy
  character(2) mm,dd,hh

  !MPI
  integer PETOT,my_rank,ierr
  integer :: master_rank=0

  !Ensemble
  integer,allocatable :: ens_year(:)

  !Model
  integer iqglobal
  integer idx(im),idy(jm)

  real(kind = 8) lon(im),lat(jm)
  real(kind = 8) fsm(im,jm)
  real(kind = 8) slpens(im,jm) !Sea Level Pressure [Pa]/raw, standard deviation,ensemble data
  real(kind = 8) uens(im,jm)  !Zonal wind speed [m/s]
  real(kind = 8) vens(im,jm)  !Meridional wind speed [m/s]
  real(kind = 8) taens(im,jm) !Air temperature [degree C]
  real(kind = 8) qaens(im,jm) !Specific humidity [g/g]
  real(kind = 8) tcens(im,jm) !Total cloud cover [0-1]
  real(kind = 8) swens(im,jm) !Shortwave radiation [W/m^2]
  real(kind = 8) lwens(im,jm) !Longwave radiation [W/m^2]
  real(kind = 8) null1dx(im),null1dy(jm),null2d(im,jm),null3d(im,jm,km)

  !JRA55
  integer :: ncount=10
  integer im_atm,jm_atm
  integer dt !Data time interval (JRA55: 6h, JRA55do: 3h)
  real(kind = 8),allocatable :: lon_atm(:),lat_atm(:)        !Longitude/Latitude
  real(kind = 8),allocatable :: land_atm(:,:)                !Land(1)/Sea(0)
  real(kind = 8),allocatable :: slp_atm(:,:),slpens_atm(:,:) !Sea Level Pressure [Pa]
  real(kind = 8),allocatable :: u_atm(:,:),uens_atm(:,:)     !10m Zonal wind speed [m/s]
  real(kind = 8),allocatable :: v_atm(:,:),vens_atm(:,:)     !10m Meridional wind speed [m/s]
  real(kind = 8),allocatable :: ta_atm(:,:),taens_atm(:,:)   !2m Temperature [degree C]
  real(kind = 8),allocatable :: qa_atm(:,:),qaens_atm(:,:)   !2m specific humidity [kg/kg = g/g]
  real(kind = 8),allocatable :: tc_atm(:,:),tcens_atm(:,:)   !Total cloud cover [0-1]
  real(kind = 8),allocatable :: sw_atm(:,:),swens_atm(:,:)   !Shortwave radiation [W/m^2]
  real(kind = 8),allocatable :: lw_atm(:,:),lwens_atm(:,:)   !Longwave radiation [W/m^2]

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,PETOT,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
  iens=my_rank+1
  nens=PETOT !Ensemble size = Process size 
  allocate(ens_year(nens))

  !Initial Setting
  open(1,file="atm_date.dat",status="old")
  read(1,*) syr,smon,sday,shour
  read(1,*) eyr,emon,eday,ehour
  close(1)

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(my_rank == master_rank)then
     status=system("rm -f atm_date.dat")
  end if
  
  im_atm=im_jra55do
  jm_atm=jm_jra55do
  dt=3

  call estimate_ntime(syr,smon,sday,shour,eyr,emon,eday,ehour,dt,ntime)

  if(my_rank == master_rank) write(*,'(a,i4.4,i2.2,i2.2,i2.2)') "Start time:",syr,smon,sday,shour
  if(my_rank == master_rank) write(*,'(a,i4.4,i2.2,i2.2,i2.2)') "End time:",eyr,emon,eday,ehour
  if(my_rank == master_rank) write(*,'(a,i10)') "Number of Calculating time:",ntime
  if(my_rank == master_rank) write(*,'(a,i10)') "Ensemble size:",nens

  !Read model grid data
  if(my_rank == master_rank) write(*,'(a)') "Read model grid"
  call read_grid(null1dx,null1dx,lon,null1dx, &
       & null1dy,null1dy,lat,null1dy,null3d,null3d, &
       & fsm,null2d,null2d)

  if(lon(1) == lon(im-1)-360.d0 .and. lon(2) == lon(im)-360.d0)then
     iqglobal=1
  else
     iqglobal=0
  end if
  
  !Allocate ATM
  allocate(lon_atm(im_atm),lat_atm(jm_atm))
  allocate(land_atm(im_atm,jm_atm))
  allocate(slp_atm(im_atm,jm_atm),slpens_atm(im_atm,jm_atm))
  allocate(u_atm(im_atm,jm_atm),uens_atm(im_atm,jm_atm))
  allocate(v_atm(im_atm,jm_atm),vens_atm(im_atm,jm_atm))
  allocate(ta_atm(im_atm,jm_atm),taens_atm(im_atm,jm_atm))
  allocate(qa_atm(im_atm,jm_atm),qaens_atm(im_atm,jm_atm))
  allocate(tc_atm(im_atm,jm_atm),tcens_atm(im_atm,jm_atm))
  allocate(sw_atm(im_atm,jm_atm),swens_atm(im_atm,jm_atm))
  allocate(lw_atm(im_atm,jm_atm),lwens_atm(im_atm,jm_atm))

  !Read JRA55 grid data
  if(my_rank == master_rank) write(*,'(a)') "Read JRA55do grid"
  call read_grid_jra55do(lon_atm,lat_atm,land_atm)

  !Calculate ID
  if(my_rank == master_rank) write(*,'(a)') "Calculate ID"
  call cal_idlon(im_atm,lon_atm,im,lon,idx)
  call cal_idlat(jm_atm,lat_atm,jm,lat,idy)

  !Set iyr,imon,iday,ihour
  allocate(iyr(ntime),imon(ntime),iday(ntime),ihour(ntime))
  !itime=1
  iyr(1)=syr
  imon(1)=smon
  iday(1)=sday
  ihour(1)=shour
  !itime=2-ntime
  do itime=2,ntime
     iyr(itime)=iyr(itime-1)
     imon(itime)=imon(itime-1)
     iday(itime)=iday(itime-1)
     ihour(itime)=ihour(itime-1)
     call add_time(iyr(itime),imon(itime),iday(itime),ihour(itime),dt)
  end do

  do itime=1,ntime

     write(yyyy,'(i4.4)') iyr(itime)
     write(mm,'(i2.2)') imon(itime)
     write(dd,'(i2.2)') iday(itime)
     write(hh,'(i2.2)') ihour(itime)

     if(my_rank == master_rank) write(*,'(a)') "-----Start "//yyyy//mm//dd//hh//"-----"
     
     !Make different year at ensemble member
     !*Change ens_year month by month
     !JRA:syr1979, nyr=40
     if(itime == 1 .or. (iday(itime) == 1 .and. ihour(itime) == 0))then
        call set_ens_year(syr_atm,eyr_atm-syr_atm+1,iyr(itime),imon(itime),1,0,nens,ens_year)
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ens_year(1),nens,MPI_INTEGER,master_rank,MPI_COMM_WORLD,ierr)
     end if

     !Read JRA55do at time t
     if(my_rank == master_rank) write(*,'(a)') "Read JRA55do"
     call read_jra55do(iyr(itime),imon(itime),iday(itime),ihour(itime), &
          & u_atm,v_atm,ta_atm,qa_atm,lw_atm,sw_atm,slp_atm)
     tc_atm(:,:)=rmiss
     
     !Read JRA55do at time t+dt
     if(mod(ens_year(iens),4) /= 0 .and. imon(itime) == 2 .and. iday(itime) == 29)then
        call read_jra55do(ens_year(iens),3,1,ihour(itime), &
             & uens_atm,vens_atm,taens_atm,qaens_atm,lwens_atm,swens_atm,slpens_atm)           
     else
        call read_jra55do(ens_year(iens),imon(itime),iday(itime),ihour(itime), &
             & uens_atm,vens_atm,taens_atm,qaens_atm,lwens_atm,swens_atm,slpens_atm)
     end if
     tcens_atm(:,:)=rmiss

     !Calculate Perturbed Atmospheric forcing (Kunii and Miyoshi 2012)
     if(my_rank == master_rank) write(*,'(a)') "Calculate Perturbed Atm. forcing"
     call cal_mpi_ens(nens,im_atm,jm_atm,1,qa_atm,qaens_atm,factor)
     call cal_mpi_ens(nens,im_atm,jm_atm,1,ta_atm,taens_atm,factor)
     call cal_mpi_ens(nens,im_atm,jm_atm,1,u_atm,uens_atm,factor)
     call cal_mpi_ens(nens,im_atm,jm_atm,1,v_atm,vens_atm,factor)
     call cal_mpi_ens(nens,im_atm,jm_atm,1,slp_atm,slpens_atm,factor)
     call cal_mpi_ens(nens,im_atm,jm_atm,1,tc_atm,tcens_atm,factor)
     call cal_mpi_ens(nens,im_atm,jm_atm,1,sw_atm,swens_atm,factor)
     call cal_mpi_ens(nens,im_atm,jm_atm,1,lw_atm,lwens_atm,factor)

     !Apply land for JRA55 DATA
     if(my_rank == master_rank) write(*,'(a)') "Apply JRA55 Land"
     call apply_jra55_land(im_atm,jm_atm,qaens_atm,land_atm,rmiss)
     call apply_jra55_land(im_atm,jm_atm,taens_atm,land_atm,rmiss)
     call apply_jra55_land(im_atm,jm_atm,uens_atm,land_atm,rmiss)
     call apply_jra55_land(im_atm,jm_atm,vens_atm,land_atm,rmiss)
     call apply_jra55_land(im_atm,jm_atm,slpens_atm,land_atm,rmiss)
     call apply_jra55_land(im_atm,jm_atm,tcens_atm,land_atm,rmiss)
     call apply_jra55_land(im_atm,jm_atm,swens_atm,land_atm,rmiss)
     call apply_jra55_land(im_atm,jm_atm,lwens_atm,land_atm,rmiss)

     !Fill value
     if(my_rank == master_rank) write(*,'(a)') "Fill value"
     call fillvalue_2d(ncount,im_atm,jm_atm,1,qaens_atm,rmiss)
     call fillvalue_2d(ncount,im_atm,jm_atm,1,taens_atm,rmiss)
     call fillvalue_2d(ncount,im_atm,jm_atm,1,uens_atm,rmiss)
     call fillvalue_2d(ncount,im_atm,jm_atm,1,vens_atm,rmiss)
     call fillvalue_2d(ncount,im_atm,jm_atm,1,slpens_atm,rmiss)
     call fillvalue_2d(ncount,im_atm,jm_atm,1,tcens_atm,rmiss)
     call fillvalue_2d(ncount,im_atm,jm_atm,1,swens_atm,rmiss)
     call fillvalue_2d(ncount,im_atm,jm_atm,1,lwens_atm,rmiss)

     !Bilinear interpolation
     if(my_rank == master_rank) write(*,'(a)') "Bilinear interpolation"
     call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,qaens_atm, &
          & im,jm,lon,lat,qaens,idx,idy,rmiss)
     call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,taens_atm, &
          & im,jm,lon,lat,taens,idx,idy,rmiss)
     call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,uens_atm, &
          & im,jm,lon,lat,uens,idx,idy,rmiss)
     call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,vens_atm, &
          & im,jm,lon,lat,vens,idx,idy,rmiss)
     call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,slpens_atm, &
          & im,jm,lon,lat,slpens,idx,idy,rmiss)
     call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,tcens_atm, &
          & im,jm,lon,lat,tcens,idx,idy,rmiss)
     call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,swens_atm, &
          & im,jm,lon,lat,swens,idx,idy,rmiss)
     call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,lwens_atm, &
          & im,jm,lon,lat,lwens,idx,idy,rmiss)
     
     call apply_fsm(im,jm,1,qaens,fsm)
     call apply_fsm(im,jm,1,taens,fsm)
     call apply_fsm(im,jm,1,uens,fsm)
     call apply_fsm(im,jm,1,vens,fsm)
     call apply_fsm(im,jm,1,slpens,fsm)
     call apply_fsm(im,jm,1,tcens,fsm)
     call apply_fsm(im,jm,1,swens,fsm)
     call apply_fsm(im,jm,1,lwens,fsm)

     !Check range
     !$omp parallel
     !$omp do private(i,j)
     do j=1,jm
        do i=1,im
           if(qaens(i,j) < 0.d0) qaens(i,j)=0.d0
           if(slpens(i,j) < 0.d0) slpens(i,j)=0.d0
           if(tcens(i,j) < 0.d0) tcens(i,j)=0.d0
           if(tcens(i,j) > 1.d0) tcens(i,j)=1.d0
           if(swens(i,j) < 0.d0) swens(i,j)=0.d0
        end do
     end do
     !$omp end do
     !$omp end parallel

     if(iqglobal == 1)then
        uens(1:2,:)=uens(im-1:im,:)
        vens(1:2,:)=vens(im-1:im,:)
        taens(1:2,:)=taens(im-1:im,:)
        qaens(1:2,:)=qaens(im-1:im,:)
        slpens(1:2,:)=slpens(im-1:im,:)
        lwens(1:2,:)=lwens(im-1:im,:)
        swens(1:2,:)=swens(im-1:im,:)        
     end if
     
     !Write data
     if(my_rank == master_rank) write(*,'(a)') "Write data"
     call write_data(itime,ntime,syr,smon,sday,shour,iyr,imon,iday,ihour,dt,iens,im,jm, &
          & uens,vens,taens,qaens,swens,lwens,tcens,slpens)

     if(my_rank == master_rank) write(*,'(a)') "-----End "//yyyy//mm//dd//hh//"-----"
     call MPI_Barrier(MPI_COMM_WORLD,ierr)

  end do

  deallocate(ens_year)  
  deallocate(iyr,imon,iday,ihour)
  deallocate(lon_atm,lat_atm)
  deallocate(land_atm)
  deallocate(slp_atm,slpens_atm)
  deallocate(u_atm,uens_atm)
  deallocate(v_atm,vens_atm)
  deallocate(ta_atm,taens_atm)
  deallocate(qa_atm,qaens_atm)
  deallocate(tc_atm,tcens_atm)
  deallocate(sw_atm,swens_atm)
  deallocate(lw_atm,lwens_atm)

 call MPI_FINALIZE(ierr)

  if(my_rank == master_rank) write(*,'(a)') "----- All make_atm_ens job end -----"

end program main

!---------------------------------------------------------------
! Write Data |
!---------------------------------------------------------------

subroutine write_data(itime,ntime,syr,smon,sday,shour,iyr,imon,iday,ihour,dt,iens,im,jm, &
     & u,v,ta,qa,sw,lw,tc,slp)

  use netcdf
  implicit none

  integer ijul,sjul
  integer status,system
  integer ncid,varid
  integer i

  real(kind = 4) tmp1d(1),tmp2d(im,jm)
  
  real(kind = 8) rjul
  real(kind = 8) atmtime
  
  character(100) filename
  character(5) nnnnn
  
  !IN
  integer,intent(in) :: itime,ntime
  integer,intent(in) :: syr,smon,sday,shour
  integer,intent(in) :: iyr(ntime),imon(ntime),iday(ntime),ihour(ntime)
  integer,intent(in) :: dt
  integer,intent(in) :: iens
  integer,intent(in)::  im,jm

  real(kind = 8),intent(in) :: u(im,jm),v(im,jm)
  real(kind = 8),intent(in) :: ta(im,jm),qa(im,jm)
  real(kind = 8),intent(in) :: sw(im,jm),lw(im,jm)
  real(kind = 8),intent(in) :: tc(im,jm),slp(im,jm)
    
  !Character
  write(nnnnn,'(i5.5)') iens
  
  !Makefile
  filename="../in/atm."//nnnnn//".nc"
  if(itime == 1)then
     status=system("rm -f "//trim(filename))
     call make_ncfile_atm(im,jm,syr,smon,sday,filename)
  end if
               
  !Open netcdf file
  status=nf90_open(trim(filename),nf90_write,ncid)

  !Atmtime
  if(itime == 1)then
     do i=1,ntime
        !Calculate atmtime
        call ymd_julian(syr,smon,sday,sjul)
        call ymd_julian(iyr(i),imon(i),iday(i),ijul)
        rjul=dble(ijul)+dble(ihour(i))/24.d0
        atmtime=rjul-dble(sjul)
        
        !Write atmtime
        tmp1d(1)=atmtime
        status=nf90_inq_varid(ncid,"atmtime",varid)
        status=nf90_put_var(ncid,varid,tmp1d,(/i/),(/1/))
     end do
  end if
        
  !Write windu
  tmp2d(:,:)=real(u(:,:))
  status=nf90_inq_varid(ncid,"windu",varid)
  status=nf90_put_var(ncid,varid,tmp2d,(/1,1,itime/),(/im,jm,1/))
  
  !Write windv
  tmp2d(:,:)=real(v(:,:))
  status=nf90_inq_varid(ncid,"windv",varid)
  status=nf90_put_var(ncid,varid,tmp2d,(/1,1,itime/),(/im,jm,1/))
  
  !Write airt
  tmp2d(:,:)=real(ta(:,:))
  status=nf90_inq_varid(ncid,"airt",varid)
  status=nf90_put_var(ncid,varid,tmp2d,(/1,1,itime/),(/im,jm,1/))
  
  !Write airh
  tmp2d(:,:)=real(qa(:,:))
  status=nf90_inq_varid(ncid,"airh",varid)
  status=nf90_put_var(ncid,varid,tmp2d,(/1,1,itime/),(/im,jm,1/))
  
  !Write swrad
  tmp2d(:,:)=real(sw(:,:))
  status=nf90_inq_varid(ncid,"swrad",varid)
  status=nf90_put_var(ncid,varid,tmp2d,(/1,1,itime/),(/im,jm,1/))
  
  !Write lwrad
  tmp2d(:,:)=real(lw(:,:))
  status=nf90_inq_varid(ncid,"lwrad",varid)
  status=nf90_put_var(ncid,varid,tmp2d,(/1,1,itime/),(/im,jm,1/))

  !Write slp
  tmp2d(:,:)=real(slp(:,:))
  status=nf90_inq_varid(ncid,"slp",varid)
  status=nf90_put_var(ncid,varid,tmp2d,(/1,1,itime/),(/im,jm,1/))

  status=nf90_close(ncid)
     
end subroutine write_data
