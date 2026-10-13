module setting

  integer,parameter :: iswitch_atm=2 !1: JRA55do, 2: ERA5
  integer,parameter :: nvar_atm=8    !u,v,ta,qa,lw,sw,slp,prep
  !logical,parameter :: lswitch_remove=.true. !Remove file
  logical,parameter :: lswitch_remove=.false.

end module setting

!------------------------------------------------------------------

program main

  !----------------------------------------------------------------
  ! Make Atmospheric netcdf file |
  !----------------------------------------------------------------
  !
  ! 1. Calculate the number of caluculating times
  ! 2. Read model grid
  ! 3. Read atm grid & land
  ! 4. Calculate ID
  ! 5. Read atm
  ! 6. Apply atm land mask
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
  ! 2024.12 Added Rainfall by S.Ohishi
  ! 2025.08 Added ERA5 by S.Ohishi
  ! 2026.10 Modified by S.Ohishi
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

  use setting, only: iswitch_atm,nvar_atm
  use mod_rmiss
  use mod_gridinfo
  use mod_read_jra55do, im_jra55do => im, jm_jra55do => jm, read_grid_jra55do => read_grid
  use mod_read_era5, im_era5 => im, jm_era5 => jm, read_grid_era5 => read_grid
  implicit none

  !---Parameter
  integer,parameter :: id_u=1
  integer,parameter :: id_v=2
  integer,parameter :: id_ta=3
  integer,parameter :: id_qa=4
  integer,parameter :: id_lw=5
  integer,parameter :: id_sw=6
  integer,parameter :: id_slp=7
  integer,parameter :: id_prep=8
  
  !---Common
  integer itime,ntime,ntime_d
  integer iyr,imon,iday,ihour
  integer syr,smon,sday,shour
  integer eyr,emon,eday,ehour
  integer ivar_atm

  character(4) yyyy
  character(2) mm,dd,hh

  !---Model
  integer iqglobal
  integer idx(im),idy(jm)

  real(kind = 8) lon(im),lat(jm)
  real(kind = 8) fsm(im,jm)
  real(kind = 8) dat(im,jm,nvar_atm)
  real(kind = 8) null1dx(im),null1dy(jm),null2d(im,jm),null3d(im,jm,km)

  !ATM
  integer :: ncount=10
  integer im_atm,jm_atm
  integer dt                                          !Data time interval (JRA55do: 3h, ERA5: 1h)
  real(kind = 8),allocatable :: lon_atm(:),lat_atm(:) !Longitude/Latitude
  real(kind = 8),allocatable :: land_atm(:,:)         !Land(1)/Sea(0)
  real(kind = 8),allocatable :: dat_atm(:,:,:)

  !Initial Setting
  call read_argument(syr,smon,sday,shour,eyr,emon,eday,ehour)
  
  if(iswitch_atm == 1)then
     im_atm=im_jra55do
     jm_atm=jm_jra55do
     dt=3
  else if(iswitch_atm == 2)then
     im_atm=im_era5
     jm_atm=jm_era5
     dt=1
  else
     write(*,*) "***Error: Incorrect iswitch_atm => ", iswitch_atm
     stop
  end if

  call estimate_ntime(syr,smon,sday,shour,eyr,emon,eday,ehour,dt,ntime)
  ntime_d=24/dt

  write(*,'(a,i4.4,i2.2,i2.2,i2.2)') "Start time:",syr,smon,sday,shour
  write(*,'(a,i4.4,i2.2,i2.2,i2.2)') "End time:",eyr,emon,eday,ehour
  write(*,'(a,i10)') "Number of Calculating time:",ntime

  !Read model grid data
  write(*,'(a)') "Read model grid"
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
  allocate(dat_atm(im_atm,jm_atm,nvar_atm))

  !Read ATM grid data
  if(iswitch_atm == 1)then
     write(*,'(a)') "Read JRA55do grid"
     call read_grid_jra55do(lon_atm,lat_atm,land_atm)
  else if(iswitch_atm == 2)then
     write(*,'(a)') "Read ERA5 grid"
     call read_grid_era5(lon_atm,lat_atm,land_atm)
  endif
     
  !Calculate ID
  write(*,*) "Calculate ID"
  call cal_idlon(im_atm,lon_atm,im,lon,idx)
  call cal_idlat(jm_atm,lat_atm,jm,lat,idy)

  !Set initial iyr,imon,iday,ihour  
  iyr=syr
  imon=smon
  iday=sday
  ihour=shour

  do itime=1,ntime

     write(yyyy,'(i4.4)') iyr
     write(mm,'(i2.2)') imon
     write(dd,'(i2.2)') iday
     write(hh,'(i2.2)') ihour

     write(*,'(a)') "-----Start "//yyyy//mm//dd//hh//"-----"
     
     !Read Atmospheric data
     if(iswitch_atm == 1)then
        write(*,'(a)') "Read JRA55do"
        call read_jra55do(iyr,imon,iday,ihour, &
             & dat_atm(:,:,id_u),dat_atm(:,:,id_v),dat_atm(:,:,id_ta),dat_atm(:,:,id_qa), &
             & dat_atm(:,:,id_lw),dat_atm(:,:,id_sw),dat_atm(:,:,id_slp),dat_atm(:,:,id_prep))
     else if(iswitch_atm == 2)then
        write(*,'(a)') "Read ERA5"
        call read_era5(iyr,imon,iday,ihour, &
             & dat_atm(:,:,id_u),dat_atm(:,:,id_v),dat_atm(:,:,id_ta),dat_atm(:,:,id_qa), &
             & dat_atm(:,:,id_lw),dat_atm(:,:,id_sw),dat_atm(:,:,id_slp),dat_atm(:,:,id_prep))
     end if
        
     !Apply land from Atmospheric DATA
     write(*,'(a)') "Apply Land"
     do ivar_atm=1,nvar_atm
        if(ivar_atm == id_slp .or. ivar_atm == id_prep)then
           cycle
        else
           call apply_jra55_land(im_atm,jm_atm,dat_atm(:,:,ivar_atm),land_atm,rmiss)
        end if
     end do

     !Fill value
     write(*,'(a)') "Fill value"
     do ivar_atm=1,nvar_atm
        if(ivar_atm == id_slp .or. ivar_atm == id_prep)then
           cycle
        else
           call fillvalue_2d(ncount,im_atm,jm_atm,1,dat_atm(:,:,ivar_atm),rmiss)
        end if
     end do
        
     !Bilinear interpolation
     write(*,'(a)') "Bilinear interpolation"
     do ivar_atm=1,nvar_atm
        call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,dat_atm(:,:,ivar_atm), &
             & im,jm,lon,lat,dat(:,:,ivar_atm),idx,idy,dble(rmiss))
     end do

     do ivar_atm=1,nvar_atm
        call apply_fsm(im,jm,1,dat(:,:,ivar_atm),fsm)
     end do

     if(iqglobal == 1)then
        dat(1:2,:,:)=dat(im-1:im,:,:)
     end if

     !Write time
     if(ihour == 0)then
        call write_time(ntime_d,iyr,imon,iday,dt,im,jm)
     end if

     !Write data
     write(*,'(a)') "Write data"
     call write_data(iyr,imon,iday,ihour,dt,im,jm, &
          & dat(:,:,id_u),dat(:,:,id_v),dat(:,:,id_ta),dat(:,:,id_qa), &
          & dat(:,:,id_lw),dat(:,:,id_sw),dat(:,:,id_slp),dat(:,:,id_prep))

     write(*,'(a)') "-----End "//yyyy//mm//dd//hh//"-----"

     call add_time(iyr,imon,iday,ihour,dt)

  end do !itime

  deallocate(lon_atm,lat_atm)
  deallocate(land_atm)
  deallocate(dat_atm)

end program main

!--------------------------------------------------------------------------------
! Read argument |
!--------------------------------------------------------------------------------

subroutine read_argument(syr,smon,sday,shour,eyr,emon,eday,ehour)

  implicit none

  !---Common
  integer i,length,status

  character(:),allocatable :: arg
  
  intrinsic :: command_argument_count, get_command_argument
  
  !---Out
  integer,intent(out) :: syr,smon,sday,shour
  integer,intent(out) :: eyr,emon,eday,ehour

  if(command_argument_count() /= 8)then
     write(*,*) "***Error: read_argument => syr smon sday shour eyr emon eday ehour"
     stop
  end if
  
  do i=1,command_argument_count()

     call get_command_argument(i,length=length,status=status)

     if(status /= 0)then
        write(*,*) "Error: argument ",status
        stop
     else

        allocate(character(length) :: arg)

        call get_command_argument(i,arg,status=status)

        if(i == 1)then
           read(arg,'(I4)') syr
        else if(i == 2)then
           read(arg,'(I2)') smon
        else if(i == 3)then
           read(arg,'(I2)') sday
        else if(i == 4)then
           read(arg,'(I2)') shour
        else if(i == 5)then
           read(arg,'(I4)') eyr
        else if(i == 6)then
           read(arg,'(I2)') emon
        else if(i == 7)then
           read(arg,'(I2)') eday
        else if(i == 8)then
           read(arg,'(I2)') ehour
        end if
        
        deallocate(arg)

     end if

  end do
  
end subroutine read_argument

!---------------------------------------------------------------
! Write Data |
!---------------------------------------------------------------

subroutine write_time(ntime,iyr,imon,iday,dt,im,jm)

  use setting
  use netcdf
  use mod_check_netcdf
  implicit none

  integer status,system,access
  integer ncid,varid
  integer itime

  real(kind = 4) tmp1d(1)

  real(kind = 8) atmtime

  character(100) filename
  character(4) yyyy
  character(2) mm,dd

  !IN  
  integer,intent(in) :: ntime
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: dt
  integer,intent(in) :: im,jm

  !filename
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday

  filename="../in/atm."//yyyy//mm//dd//".nc"

  !Makefile
  if(lswitch_remove)then
     status=system("rm -f "//trim(filename))
  end if

  status=access(trim(filename)," ")
  if(status /= 0)then
     call make_ncfile_atm(im,jm,iyr,imon,iday,filename)
  end if
  
  !Open netcdf file
  status=nf90_open(trim(filename),nf90_write,ncid)
  call check_netcdf(status)

  do itime=1,ntime

     atmtime=dble((itime-1)*dt)/24.d0

     !Write atmtime
     tmp1d(1)=atmtime
     status=nf90_inq_varid(ncid,"atmtime",varid)
     call check_netcdf(status)
     status=nf90_put_var(ncid,varid,tmp1d,(/itime/),(/1/))
     call check_netcdf(status)

  end do

  !Close netcdf file
  status=nf90_close(ncid)
  call check_netcdf(status)

end subroutine write_time

!-------------------

subroutine write_data(iyr,imon,iday,ihour,dt,im,jm, &
     & u,v,ta,qa,lw,sw,slp,prep)

  use netcdf
  use mod_check_netcdf
  implicit none

  integer itime
  integer status,ncid,varid

  real(kind = 4) tmp2d(im,jm)

  character(100) filename  
  character(4) yyyy
  character(2) mm,dd

  !IN
  integer,intent(in) :: iyr,imon,iday,ihour
  integer,intent(in) :: dt
  integer,intent(in)::  im,jm

  real(kind = 8),intent(in) :: u(im,jm),v(im,jm)
  real(kind = 8),intent(in) :: ta(im,jm),qa(im,jm)
  real(kind = 8),intent(in) :: lw(im,jm),sw(im,jm)
  real(kind = 8),intent(in) :: slp(im,jm),prep(im,jm)

  !Calculate itime
  itime=ihour/dt+1

  !filename
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday

  filename="../in/atm."//yyyy//mm//dd//".nc"

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
