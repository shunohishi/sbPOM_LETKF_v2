module setting

  integer,parameter :: iswitch_atm=1 !1: JRA55do, 2: ERA5
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
  ! 3. Read JRA55 grid & land
  ! 4. Calculate ID
  ! 5. Read JRA55
  ! 6. Apply JRA55 land mask
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
  ! 2025.08 Added ERA5 by S. Ohishi
  !
  !----------------------------------------------------------------

  use setting, only: iswitch_atm
  use mod_rmiss
  use mod_gridinfo
  use mod_read_jra55do, im_jra55do => im, jm_jra55do => jm, read_grid_jra55do => read_grid
  use mod_read_era5, im_era5 => im, jm_era5 => jm, read_grid_era5 => read_grid
  implicit none

  !Common
  integer i,j
  integer itime,ntime,ntime_d
  integer iyr,imon,iday,ihour
  integer syr,smon,sday,shour
  integer eyr,emon,eday,ehour

  integer status,system
  character(4) yyyy
  character(2) mm,dd,hh

  !Model
  integer iqglobal
  integer idx(im),idy(jm)

  real(kind = 8) lon(im),lat(jm)
  real(kind = 8) fsm(im,jm)
  real(kind = 8) slp(im,jm)        !Sea Level Pressure [Pa]
  real(kind = 8) u(im,jm),v(im,jm) !Zonal,Meridional wind speed [m/s]
  real(kind = 8) ta(im,jm)         !Air temperature [degree C]
  real(kind = 8) qa(im,jm)         !Specific humidity [g/g]
  real(kind = 8) sw(im,jm)         !Shortwave radiation [W/m^2]
  real(kind = 8) lw(im,jm)         !Longwave radiation [W/m^2]
  real(kind = 8) prep(im,jm)       !Precipitation [mm/day]
  real(kind = 8) null1dx(im),null1dy(jm),null2d(im,jm),null3d(im,jm,km)

  !ATM
  integer :: ncount=10
  integer im_atm,jm_atm
  integer dt                                          !Data time interval (JRA55do: 3h, ERA5: 1h)
  real(kind = 8),allocatable :: lon_atm(:),lat_atm(:) !Longitude/Latitude
  real(kind = 8),allocatable :: land_atm(:,:)         !Land(1)/Sea(0)
  real(kind = 8),allocatable :: slp_atm(:,:)          !Sea Level Pressure [Pa]
  real(kind = 8),allocatable :: u_atm(:,:),v_atm(:,:) !10m Zonal,Meridional wind speed [m/s]
  real(kind = 8),allocatable :: ta_atm(:,:)           !2m Temperature [degree C]
  real(kind = 8),allocatable :: qa_atm(:,:)           !2m specific humidity [kg/kg = g/g]
  real(kind = 8),allocatable :: sw_atm(:,:)           !Shortwave radiation [W/m^2]
  real(kind = 8),allocatable :: lw_atm(:,:)           !Longwave radiation [W/m^2]
  real(kind = 8),allocatable :: prep_atm(:,:)         !Precipitation [mm/day]

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

     allocate(slp_atm(im_atm,jm_atm))
     allocate(u_atm(im_atm,jm_atm),v_atm(im_atm,jm_atm))
     allocate(ta_atm(im_atm,jm_atm),qa_atm(im_atm,jm_atm))
     allocate(sw_atm(im_atm,jm_atm),lw_atm(im_atm,jm_atm))
     allocate(prep_atm(im_atm,jm_atm))
     
     !Read Atmospheric data
     if(iswitch_atm == 1)then
        write(*,'(a)') "Read JRA55do"
        call read_jra55do(iyr,imon,iday,ihour,u_atm,v_atm,ta_atm,qa_atm,lw_atm,sw_atm,slp_atm,prep_atm)
     else if(iswitch_atm == 2)then
        write(*,'(a)') "Read ERA5"
        call read_era5(iyr,imon,iday,ihour,u_atm,v_atm,ta_atm,qa_atm,lw_atm,sw_atm,slp_atm,prep_atm)
     end if
        
     !Apply land from Atmospheri DATA
     write(*,'(a)') "Apply Land"
     call apply_jra55_land(im_atm,jm_atm,qa_atm,land_atm,rmiss)
     call apply_jra55_land(im_atm,jm_atm,ta_atm,land_atm,rmiss)
     call apply_jra55_land(im_atm,jm_atm,u_atm,land_atm,rmiss)
     call apply_jra55_land(im_atm,jm_atm,v_atm,land_atm,rmiss)
     !call apply_jra55_land(im_atm,jm_atm,slp_atm,land_atm,rmiss)
     call apply_jra55_land(im_atm,jm_atm,lw_atm,land_atm,rmiss)
     call apply_jra55_land(im_atm,jm_atm,sw_atm,land_atm,rmiss)
     !call apply_jra55_land(im_atm,jm_atm,prep_atm,land_atm,rmiss)

     !Fill value
     write(*,'(a)') "Fill value"
     call fillvalue_2d(ncount,im_atm,jm_atm,1,qa_atm,rmiss)
     call fillvalue_2d(ncount,im_atm,jm_atm,1,ta_atm,rmiss)
     call fillvalue_2d(ncount,im_atm,jm_atm,1,u_atm,rmiss)
     call fillvalue_2d(ncount,im_atm,jm_atm,1,v_atm,rmiss)
     !call fillvalue_2d(ncount,im_atm,jm_atm,1,slp_atm,rmiss)
     call fillvalue_2d(ncount,im_atm,jm_atm,1,lw_atm,rmiss)
     call fillvalue_2d(ncount,im_atm,jm_atm,1,sw_atm,rmiss)
     !call fillvalue_2d(ncount,im_atm,jm_atm,1,prep_atm,rmiss)

     !Bilinear interpolation
     write(*,'(a)') "Bilinear interpolation"
     call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,qa_atm, &
          & im,jm,lon,lat,qa,idx,idy,dble(rmiss))
     call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,ta_atm, &
          & im,jm,lon,lat,ta,idx,idy,dble(rmiss))
     call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,u_atm, &
          & im,jm,lon,lat,u,idx,idy,dble(rmiss))
     call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,v_atm, &
          & im,jm,lon,lat,v,idx,idy,dble(rmiss))
     call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,slp_atm, &
          & im,jm,lon,lat,slp,idx,idy,dble(rmiss))
     call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,lw_atm, &
          & im,jm,lon,lat,lw,idx,idy,dble(rmiss))
     call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,sw_atm, &
          & im,jm,lon,lat,sw,idx,idy,dble(rmiss))
     call bilinear_interpolation_2d(im_atm,jm_atm,lon_atm,lat_atm,prep_atm, &
          & im,jm,lon,lat,prep,idx,idy,-999.d0)

     deallocate(slp_atm)
     deallocate(u_atm,v_atm)
     deallocate(ta_atm,qa_atm)
     deallocate(sw_atm,lw_atm)
     deallocate(prep_atm)
     
     call apply_fsm(im,jm,1,u,fsm)
     call apply_fsm(im,jm,1,v,fsm)
     call apply_fsm(im,jm,1,ta,fsm)
     call apply_fsm(im,jm,1,qa,fsm)
     call apply_fsm(im,jm,1,slp,fsm)
     call apply_fsm(im,jm,1,lw,fsm)
     call apply_fsm(im,jm,1,sw,fsm)
     call apply_fsm(im,jm,1,prep,fsm)

     if(iqglobal == 1)then
        u(1:2,:)=u(im-1:im,:)
        v(1:2,:)=v(im-1:im,:)
        ta(1:2,:)=ta(im-1:im,:)
        qa(1:2,:)=qa(im-1:im,:)
        slp(1:2,:)=slp(im-1:im,:)
        lw(1:2,:)=lw(im-1:im,:)
        sw(1:2,:)=sw(im-1:im,:)        
        prep(1:2,:)=prep(im-1:im,:)        
     end if

     !Write time
     if(ihour == 0)then
        call write_time(ntime_d,iyr,imon,iday,dt,im,jm)
     end if

     !Write data
     write(*,'(a)') "Write data"
     call write_data(iyr,imon,iday,ihour,dt,im,jm, &
          & u,v,ta,qa,sw,lw,slp,prep)

     write(*,'(a)') "-----End "//yyyy//mm//dd//hh//"-----"

     call add_time(iyr,imon,iday,ihour,dt)

  end do !itime

  deallocate(lon_atm,lat_atm)
  deallocate(land_atm)

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


  do i=1,command_argument_count()

     call get_command_argument(i,length=length,status=status)

     if(status /= 0)then
        write(*,*) "Error: arugument ",status
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

  do itime=1,ntime

     atmtime=dble((itime-1)*dt)/24.d0

     !Write atmtime
     tmp1d(1)=atmtime
     status=nf90_inq_varid(ncid,"atmtime",varid)
     status=nf90_put_var(ncid,varid,tmp1d,(/itime/),(/1/))

  end do

  !Close netcdf file
  status=nf90_close(ncid)

end subroutine write_time

!-------------------

subroutine write_data(iyr,imon,iday,ihour,dt,im,jm, &
     & u,v,ta,qa,sw,lw,slp,prep)

  use netcdf
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
  real(kind = 8),intent(in) :: sw(im,jm),lw(im,jm)
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

  !Write prep
  tmp2d(:,:)=real(prep(:,:))
  status=nf90_inq_varid(ncid,"prep",varid)
  status=nf90_put_var(ncid,varid,tmp2d,(/1,1,itime/),(/im,jm,1/))

  status=nf90_close(ncid)

end subroutine write_data
