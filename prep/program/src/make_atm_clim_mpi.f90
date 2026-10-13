module setting

  integer,parameter :: syr=1979,eyr=2023 !Start and End year
  integer,parameter :: iswitch_atm=2 !1:JRA55do, 2:ERA5
  integer,parameter :: nvar_atm=8    !u,v,ta,qa,lw,sw,slp,prep
  integer,parameter :: stime_rst=1   !Restart timestep
  !logical,parameter :: lswitch_remove=.true.
  logical,parameter :: lswitch_remove=.false.
  
end module setting

!----------------------------------------------------------------
! Make Atmospheric netcdf file |
!----------------------------------------------------------------
!
! 1. Calculate the number of caluculating times
! 2. Read model grid
! 3. Read Atmospheric grid & land
! 4. Calculate ID
! 5. Read Atmospheric data
! 6. Apply Atmospheric land mask
! 7. Fill value
! 8. Bilinear interpolation
! 9. Write Data
!
! *Initial hour = 0 at iday = 1
!
! 2018.08 Created by S.Ohishi
! 2021.11 Added JRA55do by S.Ohishi
! 2023.04 Modified all by S.Ohishi
! 2023.07 Modified for serial POM & Remove JRA55 by S.Ohishi
! 2024.12 Add Rainfall by S.Ohishi
! 2026.04 Add ERA5 by S.Ohishi
! 2026.10 Add MPI by S.Ohishi
!
!----------------------------------------------------------------
! Atmospheric variable |
!----------------------------------------------------------------
!1: Zonal, 2: Meridional wind speed [m/s]
!3: Air temperature [degree C]
!4: Specific humidity [g/g]
!5: Longwave radiation [W/m^2]
!6: Shortwave radiation [W/m^2]
!7: Sea Level Pressure [Pa]
!8: Precipitation [mm/day]
!-----------------------------------------------------------------

program main

  use setting
  use mod_rmiss
  use mod_gridinfo
  use mod_read_jra55do, im_jra55do => im, jm_jra55do => jm, read_grid_jra55do => read_grid 
  use mod_read_era5, im_era5 => im, jm_era5 => jm, read_grid_era5 => read_grid
  use mpi
  implicit none

  !---Parameter
  integer,parameter :: master_rank=0
  integer,parameter :: id_u=1
  integer,parameter :: id_v=2
  integer,parameter :: id_ta=3
  integer,parameter :: id_qa=4
  integer,parameter :: id_lw=5
  integer,parameter :: id_sw=6
  integer,parameter :: id_slp=7
  integer,parameter :: id_prep=8
  
  !---Common
  integer ib,nb
  integer ivar_atm
  integer itime,ntime
  integer iyr,imon,iday,ihour
  integer smon,sday,shour
  integer emon,eday,ehour

  character(4) yyyy
  character(2) mm,dd,hh

  !---MPI
  integer PETOT,my_rank,ierr
  
  !---Model
  integer iqglobal
  integer idx(im),idy(jm)

  real(kind = 8) lon(im),lat(jm)
  real(kind = 8) fsm(im,jm)
  real(kind = 8) dat(im,jm,nvar_atm)
  real(kind = 8) null1dx(im),null1dy(jm),null2d(im,jm),null3d(im,jm,km)

  !---ATM
  integer :: ncount=10
  integer im_atm,jm_atm
  integer dt                                                 !Time interval (JRA55do: 3h, ERA5: 1h)
  real(kind = 8),allocatable :: lon_atm(:),lat_atm(:)        !Longitude/Latitude
  real(kind = 8),allocatable :: land_atm(:,:)                !Land(1)/Sea(0)
  real(kind = 8),allocatable :: dat_atm(:,:,:),sum_atm(:,:,:),clim_atm(:,:,:)

  !MPI
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,PETOT,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
  
  !Initial Setting
  if(iswitch_atm == 1)then
     im_atm=im_jra55do
     jm_atm=jm_jra55do
     dt=3
  else if(iswitch_atm == 2)then
     im_atm=im_era5
     jm_atm=jm_era5
     dt=1
  else
     write(*,*) "***Error: Incorrect iswitch_atm => ",iswitch_atm
     call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  end if

  smon=1
  sday=1
  shour=0
  emon=12
  eday=31
  ehour=24-dt
  
  call estimate_ntime(1,smon,sday,shour,1,emon,eday,ehour,dt,ntime)

  if(my_rank == master_rank)then
     write(*,'(a,i2.2,i2.2,i2.2)') "Start time:",smon,sday,shour
     write(*,'(a,i2.2,i2.2,i2.2)') "End time:",emon,eday,ehour
     write(*,'(a,i10)') "Number of Calculating time:",ntime
  end if
     
  !Read model grid data
  if(my_rank == master_rank)then
     write(*,'(a)') "Read model grid"
     call read_grid(null1dx,null1dx,lon,null1dx, &
          & null1dy,null1dy,lat,null1dy,null3d,null3d, &
          & fsm,null2d,null2d)
  end if
  call MPI_Bcast(lon,im,MPI_DOUBLE_PRECISION,master_rank,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(lat,jm,MPI_DOUBLE_PRECISION,master_rank,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(fsm,im*jm,MPI_DOUBLE_PRECISION,master_rank,MPI_COMM_WORLD,ierr)
  
  if(lon(1) == lon(im-1)-360.d0 .and. lon(2) == lon(im)-360.d0)then
     iqglobal=1
  else
     iqglobal=0
  end if

  !Allocate ATM
  allocate(lon_atm(im_atm),lat_atm(jm_atm))
  allocate(land_atm(im_atm,jm_atm))
  allocate(dat_atm(im_atm,jm_atm,nvar_atm))
  allocate(sum_atm(im_atm,jm_atm,nvar_atm))
  allocate(clim_atm(im_atm,jm_atm,nvar_atm))

  !Read ATM grid data
  if(my_rank == master_rank)then
     if(iswitch_atm == 1)then
        write(*,'(a)') "Read JRA55do grid"
        call read_grid_jra55do(lon_atm,lat_atm,land_atm)
     else if(iswitch_atm == 2)then
        write(*,'(a)') "Read ERA5 grid"
        call read_grid_era5(lon_atm,lat_atm,land_atm)
     end if
  end if
  call MPI_Bcast(lon_atm,im_atm,MPI_DOUBLE_PRECISION,master_rank,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(lat_atm,jm_atm,MPI_DOUBLE_PRECISION,master_rank,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(land_atm,im_atm*jm_atm,MPI_DOUBLE_PRECISION,master_rank,MPI_COMM_WORLD,ierr)
  
  !Calculate ID
  if(my_rank == master_rank) write(*,*) "Calculate ID"
  call cal_idlon(im_atm,lon_atm,im,lon,idx)
  call cal_idlat(jm_atm,lat_atm,jm,lat,idy)
     
  !Write time
  if(my_rank == master_rank)then
     call write_time(ntime,1,smon,sday,shour,dt,im,jm)
  end if
     
  !Set initial imon,iday,ihour
  imon=smon
  iday=sday
  ihour=shour

  do itime=1,ntime
     
     write(mm,'(i2.2)') imon
     write(dd,'(i2.2)') iday
     write(hh,'(i2.2)') ihour

     if(itime < stime_rst)then
        
        if(my_rank == master_rank)then
           write(*,'(a)') "-----Skip "//mm//dd//hh//"-----"
           write(*,'(a,i10)') "itime:",itime
        end if

        iyr=2001 !Not leap year
        call add_time(iyr,imon,iday,ihour,dt)
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        cycle
        
     end if
     
     if(my_rank == master_rank)then
        write(*,'(a)') "-----Start "//mm//dd//hh//"-----"
        write(*,'(a,i10)') "itime:",itime
     end if
     
     sum_atm(:,:,:)=0.d0        
     clim_atm(:,:,:)=0.d0
     nb=(eyr-syr+1+PETOT-1)/PETOT !Cycle size

     do ib=1,nb

        !iyr
        iyr=syr+my_rank + (ib-1)*PETOT
        write(yyyy,'(i4.4)') iyr

        if(syr <= iyr .and. iyr <= eyr)then
           !Read Atm data        
           if(iswitch_atm == 1)then
              !write(*,'(a)') "Read JRA55do at "//yyyy//mm//dd//hh
              call read_jra55do(iyr,imon,iday,ihour, &
                   & dat_atm(:,:,id_u),dat_atm(:,:,id_v),dat_atm(:,:,id_ta),dat_atm(:,:,id_qa), &
                   & dat_atm(:,:,id_lw),dat_atm(:,:,id_sw),dat_atm(:,:,id_slp),dat_atm(:,:,id_prep))
           else if(iswitch_atm == 2)then
              !write(*,'(a)') "Read ERA5 at "//yyyy//mm//dd//hh
              call read_era5(iyr,imon,iday,ihour, &
                   & dat_atm(:,:,id_u),dat_atm(:,:,id_v),dat_atm(:,:,id_ta),dat_atm(:,:,id_qa), &
                   & dat_atm(:,:,id_lw),dat_atm(:,:,id_sw),dat_atm(:,:,id_slp),dat_atm(:,:,id_prep))
           end if
           sum_atm(:,:,:)=sum_atm(:,:,:)+dat_atm(:,:,:)
        end if   
     
     end do !ib

     call clim_mpi(master_rank,im_atm,jm_atm,sum_atm,clim_atm)
     
     if(my_rank == master_rank)then

        !Apply atm land
        write(*,'(a)') "Apply Atm Land"
        do ivar_atm=1,nvar_atm
           if(ivar_atm == id_slp .or. ivar_atm == id_prep)then
              cycle
           else        
              call apply_jra55_land(im_atm,jm_atm,clim_atm(:,:,ivar_atm),land_atm,rmiss)
           end if
        end do

        write(*,'(a)') "Fill value"
        do ivar_atm=1,nvar_atm
           if(ivar_atm == id_slp .or. ivar_atm == id_prep)then
              cycle
           else        
              call fillvalue_2d(ncount,im_atm,jm_atm,1,clim_atm(:,:,ivar_atm),rmiss)
           end if
        end do
        
        !Bilinear interpolation
        write(*,'(a)') "Bilinear interpolation"
        do ivar_atm=1,nvar_atm
           call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,clim_atm(:,:,ivar_atm), &
                & im,jm,lon,lat,dat(:,:,ivar_atm),idx,idy,dble(rmiss))
        end do

        do ivar_atm=1,nvar_atm           
           call apply_fsm(im,jm,1,dat(:,:,ivar_atm),fsm)
        end do
           
        if(iqglobal == 1)then
           dat(1:2,:,:)=dat(im-1:im,:,:)
        end if
     
        !Write data
        write(*,'(a)') "Write data"
        call write_data(1,smon,sday,shour,1,imon,iday,ihour,dt,im,jm, &
             & dat(:,:,id_u),dat(:,:,id_v),dat(:,:,id_ta),dat(:,:,id_qa), &
             & dat(:,:,id_lw),dat(:,:,id_sw),dat(:,:,id_slp),dat(:,:,id_prep))

     end if

     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     
     if(my_rank == master_rank)then
        write(*,'(a)') "-----End "//mm//dd//hh//"-----"
     end if
        
     iyr=2001 !Not leap year
     call add_time(iyr,imon,iday,ihour,dt)

  end do !itime

  deallocate(lon_atm,lat_atm)
  deallocate(land_atm)
  deallocate(dat_atm,sum_atm,clim_atm)

  call MPI_FINALIZE(ierr)
  
end program main

!---------------------------------------------------------------
! Climatology |
!---------------------------------------------------------------

subroutine clim_mpi(master_rank,im,jm,sum,clim)

  use setting, only: syr,eyr,nvar_atm
  use mpi
  implicit none

  !---Common
  integer ierr
  
  !---IN
  integer,intent(in) :: master_rank
  integer,intent(in) :: im,jm
  real(kind = 8),intent(in) :: sum(im,jm,nvar_atm)

  !---OUT
  real(kind = 8),intent(inout) :: clim(im,jm,nvar_atm)
  
  call MPI_Reduce(sum,clim,im*jm*nvar_atm,MPI_DOUBLE_PRECISION,MPI_SUM,master_rank,MPI_COMM_WORLD,ierr)

  clim(:,:,:)=clim(:,:,:)/dble(eyr-syr+1)
  
end subroutine clim_mpi

!---------------------------------------------------------------
! Write Data |
!---------------------------------------------------------------

subroutine write_time(ntime,syr,smon,sday,shour,dt,im,jm)

  use setting, only: lswitch_remove, stime_rst
  use netcdf
  use mod_check_netcdf
  implicit none

  integer status,system,access
  integer ncid,varid
  integer iyr,imon,iday,ihour
  integer ijul,sjul
  integer itime

  real(kind = 4) tmp1d(1)

  real(kind = 8) rjul
  real(kind = 8) atmtime

  character(100) :: filename="../in/atm_clim.nc"

  !IN  
  integer,intent(in) :: ntime
  integer,intent(in) :: syr,smon,sday,shour
  integer,intent(in) :: dt
  integer,intent(in) :: im,jm

  !Makefile
  if(lswitch_remove .and. stime_rst == 1)then
     status=system("rm -f "//trim(filename))
  end if

  status=access(trim(filename)," ")
  if(status /= 0)then
     call make_ncfile_atm(im,jm,syr,smon,sday,filename)
  end if
  
  !Open netcdf file
  status=nf90_open(trim(filename),nf90_write,ncid)
  call check_netcdf(status)
  
  !Set iyr,imon,iday,ihour
  iyr=syr
  imon=smon
  iday=sday
  ihour=shour

  do itime=1,ntime

     !Calculate atmtime
     call ymd_julian(syr,smon,sday,sjul)
     call ymd_julian(iyr,imon,iday,ijul)
     rjul=dble(ijul)+dble(ihour)/24.d0
     atmtime=rjul-dble(sjul)

     !Write atmtime
     tmp1d(1)=atmtime
     status=nf90_inq_varid(ncid,"atmtime",varid)
     call check_netcdf(status)
     status=nf90_put_var(ncid,varid,tmp1d,(/itime/),(/1/))
     call check_netcdf(status)

     call add_time(iyr,imon,iday,ihour,dt)

  end do

  !Close netcdf file
  status=nf90_close(ncid)
  call check_netcdf(status)

end subroutine write_time

!-------------------

subroutine write_data(syr,smon,sday,shour,iyr,imon,iday,ihour,dt,im,jm, &
     & u,v,ta,qa,lw,sw,slp,prep)

  use netcdf
  use mod_check_netcdf
  implicit none

  integer itime
  integer status,ncid,varid

  real(kind = 4) tmp2d(im,jm)

  character(100) :: filename="../in/atm_clim.nc"

  !IN
  integer,intent(in) :: syr,smon,sday,shour
  integer,intent(in) :: iyr,imon,iday,ihour
  integer,intent(in) :: dt
  integer,intent(in)::  im,jm

  real(kind = 8),intent(in) :: u(im,jm),v(im,jm)
  real(kind = 8),intent(in) :: ta(im,jm),qa(im,jm)
  real(kind = 8),intent(in) :: lw(im,jm),sw(im,jm)
  real(kind = 8),intent(in) :: slp(im,jm),prep(im,jm)

  !Calculate itime
  call estimate_ntime(syr,smon,sday,shour,iyr,imon,iday,ihour,dt,itime)
  
  !Open netcdf file
  status=nf90_open(trim(filename),nf90_write,ncid)
  call check_netcdf(status)

  !Write windu
  tmp2d(:,:)=real(u(:,:))
  status=nf90_inq_varid(ncid,"windu",varid)
  call check_netcdf(status)
  status=nf90_put_var(ncid,varid,tmp2d,(/1,1,itime/),(/im,jm,1/))
  call check_netcdf(status)

  !Write windv
  tmp2d(:,:)=real(v(:,:))
  status=nf90_inq_varid(ncid,"windv",varid)
  call check_netcdf(status)
  status=nf90_put_var(ncid,varid,tmp2d,(/1,1,itime/),(/im,jm,1/))
  call check_netcdf(status)

  !Write airt
  tmp2d(:,:)=real(ta(:,:))
  status=nf90_inq_varid(ncid,"airt",varid)
  call check_netcdf(status)
  status=nf90_put_var(ncid,varid,tmp2d,(/1,1,itime/),(/im,jm,1/))
  call check_netcdf(status)

  !Write airh
  tmp2d(:,:)=real(qa(:,:))
  status=nf90_inq_varid(ncid,"airh",varid)
  call check_netcdf(status)
  status=nf90_put_var(ncid,varid,tmp2d,(/1,1,itime/),(/im,jm,1/))
  call check_netcdf(status)

  !Write lwrad
  tmp2d(:,:)=real(lw(:,:))
  status=nf90_inq_varid(ncid,"lwrad",varid)
  call check_netcdf(status)
  status=nf90_put_var(ncid,varid,tmp2d,(/1,1,itime/),(/im,jm,1/))
  call check_netcdf(status)
  
  !Write swrad
  tmp2d(:,:)=real(sw(:,:))
  status=nf90_inq_varid(ncid,"swrad",varid)
  call check_netcdf(status)
  status=nf90_put_var(ncid,varid,tmp2d,(/1,1,itime/),(/im,jm,1/))
  call check_netcdf(status)

  !Write slp
  tmp2d(:,:)=real(slp(:,:))
  status=nf90_inq_varid(ncid,"slp",varid)
  call check_netcdf(status)
  status=nf90_put_var(ncid,varid,tmp2d,(/1,1,itime/),(/im,jm,1/))
  call check_netcdf(status)

  !Write prep
  tmp2d(:,:)=real(prep(:,:))
  status=nf90_inq_varid(ncid,"prep",varid)
  call check_netcdf(status)
  status=nf90_put_var(ncid,varid,tmp2d,(/1,1,itime/),(/im,jm,1/))
  call check_netcdf(status)

  status=nf90_close(ncid)
  call check_netcdf(status)

end subroutine write_data
