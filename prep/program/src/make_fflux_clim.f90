module type_setting

  integer,parameter :: dt=3 !Data time interval [hour]
  integer,parameter :: syr=1981,eyr=2023
!  logical,parameter :: lswitch_remove=.true.
  logical,parameter :: lswitch_remove=.false.

end module type_setting

!________________________________________________________________

program main

  !----------------------------------------------------------------
  ! Make freshwater flux netcdf file |
  !----------------------------------------------------------------
  !
  ! 1. Calculate the number of caluculating times
  ! 2. Read model grid
  ! 3. Calculate model coastline
  ! 4. Read GSMaP grid
  ! 5. Read CaMa-Flood
  ! 6. Calculate ID
  !
  !~~~ do loop
  ! 7. Read GSMaP
  ! 8. Bilinear interpolation
  ! 9. Read CaMa-Flood
  ! 10. Calculate River
  ! 11. Write Data
  ! 12. Add time
  !~~~ end do
  !
  ! *Initial hour = 0 at iday = 1
  !
  ! Created by 
  ! 2018.08 S.Ohishi
  ! 
  ! Modified by
  ! 2019.02 S.Ohishi (Modify subroutine Cama-Flood/JRA55 bin -> netcdf)
  ! 2022.07 S.Ohishi (Predeptation from JRA55 --> JRA55do)
  !         *Switch prep JRA55 --> JRA55do at 2021.01.01
  ! 2023.04 S.Ohishi (Modify all)
  ! 2024.12 S.Ohishi moved rainfall to atmospheric forcing
  ! 2026.10 S.Ohishi debugged cal_id_cama
  !
  !----------------------------------------------------------------

  use type_setting
  use mod_rmiss
  use mod_gridinfo
  use mod_read_cama, im_cama => im, jm_cama => jm, read_land_cama => read_land_netcdf

  implicit none

  !Common
  integer itime,ntime
  integer iyr,imon,iday,ihour
  integer smon,sday,shour
  integer emon,eday,ehour
  character(2) mm,dd,hh

  !Model
  integer iqglobal
  integer idx_cama(im,jm),idy_cama(im,jm),dnum_cama(im,jm)

  real(kind = 8) lon(im),lat(jm)
  real(kind = 8) fsm(im,jm)
  real(kind = 8) cline(im,jm) !1: Coast line, 0: Other
  real(kind = 8) dat(im,jm) !River discharge [mm/day] positive
  real(kind = 8) null1dx(im),null1dy(jm),null2d(im,jm),null3d(im,jm,km)
  
  !CaMa-Flood
  real(kind = 8) lon_cama(im_cama),lat_cama(jm_cama)
  real(kind = 8) land_cama(im_cama,jm_cama) !Land:1, Sea:0
  real(kind = 8) dat_cama(im_cama,jm_cama)
  real(kind = 8) clim_cama(im_cama,jm_cama)
  
  !Initial setting
  smon=1
  sday=1
  shour=0
  emon=12
  eday=31
  ehour=24-dt  
  call estimate_ntime(1,smon,sday,shour,1,emon,eday,ehour,dt,ntime)

  write(*,'(a,i2.2,i2.2,i2.2)') "Start time:",smon,sday,shour
  write(*,'(a,i2.2,i2.2,i2.2)') "End time:",emon,eday,ehour
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
  
  write(*,'(a)') "Calculate model coast line"
  call cal_coast_line(im,jm,fsm,cline)
  
  !Read CaMa-Flood 
  write(*,'(a)') "Read CaMa-Flood"
  call read_cama_jra55_netcdf(1981,1,1,1,lon_cama,lat_cama,dat_cama)
  call read_land_cama(land_cama)

  !Calculate ID
  write(*,*) "Calculate ID CaMa-Flood"
  call cal_id_cama(im,jm,lon,lat,cline, &
       & im_cama,jm_cama,lon_cama,lat_cama,land_cama, &
       & idx_cama,idy_cama,dnum_cama)

  if(iqglobal == 1)then
     idx_cama(1,:)=idx_cama(im-1,:)
     idx_cama(im,:)=idx_cama(2,:)
     idy_cama(1,:)=idy_cama(im-1,:)
     idy_cama(im,:)=idy_cama(2,:)
     dnum_cama(1,:)=dnum_cama(im-1,:)
     dnum_cama(im,:)=dnum_cama(2,:)
  end if

  !Write time
  call write_time(ntime,1,smon,sday,dt,im,jm)
  
  !Set initial iyr,imon,iday,ihour
  imon=smon
  iday=sday
  ihour=shour

  do itime=1,ntime

     write(mm,'(i2.2)') imon
     write(dd,'(i2.2)') iday
     write(hh,'(i2.2)') ihour

     write(*,'(a)') "-----Start "//mm//dd//hh//"-----"
     
     !Read CaMa-Flood & Calculate River
     write(*,*) "Read CaMa-Flood"
     do iyr=syr,eyr
        call read_cama_jra55_netcdf(iyr,imon,iday,ihour,lon_cama,lat_cama,dat_cama)
        call clim(iyr,syr,eyr,im_cama,jm_cama,dat_cama,clim_cama)
     end do

     dat_cama(:,:)=clim_cama(:,:)
     
     write(*,'(a)') "Calculate River"
     call cal_river(im_cama,jm_cama,dat_cama,im,jm,dat,idx_cama,idy_cama,dnum_cama)
     call apply_fsm(im,jm,1,dat,fsm)

     if(iqglobal == 1)then
        dat(1,:)=dat(im-1,:)
        dat(im,:)=dat(2,:)
     end if
     
     !Write data
     write(*,'(a)') "Write data"
     call write_data(1,smon,sday,shour,1,imon,iday,ihour,dt,im,jm,dat)

     write(*,*) "-----End "//mm//dd//hh//"-----"

     iyr=1 !dummy
     call add_time(iyr,imon,iday,ihour,dt)

  end do
  
end program main

!-------------------------------------------------------------
! Read Coast line |
!-------------------------------------------------------------

!Coast line:1, other: 0
subroutine cal_coast_line(im,jm,fsm,cline)

  implicit none

  integer i,j,i1,i2,j1,j2
  integer pass(im,jm),miss(im,jm)

  !IN
  integer,intent(in) :: im,jm
  real(kind = 8),intent(in) :: fsm(im,jm)

  !OUT
  real(kind = 8),intent(out) :: cline(im,jm)

  do j1=1,jm     
     do i1=1,im
                
        pass(i1,j1)=0
        miss(i1,j1)=0
        
        do j2=j1-1,j1+1
           if(j2 < 1 .or. jm < j2)cycle
           do i2=i1-1,i1+1
              if(i2 < 1 .or. im < i2)cycle
                         
              if(fsm(i2,j2) == 0.)then
                 miss(i1,j1)=miss(i1,j1)+1
              else
                 pass(i1,j1)=pass(i1,j1)+1
              end if
              
           end do
        end do
     end do
  end do

  !Make cline (miss = 0.)  
  cline(:,:)=0.d0
  
  do i=1,im
     do j=1,jm
        if(fsm(i,j) == 1.d0 .and. miss(i,j) /= 0)then
           cline(i,j)=1.d0
        endif
     end do
  end do
  
end subroutine cal_coast_line

!-----------------------------------

subroutine cal_river(im_cama,jm_cama,dat_cama,im,jm,dat,idx,idy,dnum)

  implicit none

  integer i,j

  !IN
  integer,intent(in) :: im_cama,jm_cama  
  integer,intent(in) :: im,jm
  integer,intent(in) :: idx(im,jm),idy(im,jm),dnum(im,jm)

  real(kind = 8),intent(in) :: dat_cama(im_cama,jm_cama)

  !OUT
  real(kind = 8),intent(out) ::  dat(im,jm)

  dat(:,:)=0.d0

  do j=1,jm
     do i=1,im

        if(idx(i,j) == 0 .or. idy(i,j) == 0 .or. dnum(i,j) == 0)cycle
        
        dat(i,j)=dat_cama(idx(i,j),idy(i,j))/dnum(i,j)

     end do
  end do

  do j=1,jm
     do i=1,im
        if(dat(i,j) < 0.d0)then
           dat(i,j)=0.d0
        end if
     end do
  end do

end subroutine cal_river

!---------------------------------------------------------------
! Climatology |
!---------------------------------------------------------------

subroutine clim(iyr,syr,eyr,im,jm,dat,mean)

  implicit none

  integer,intent(in) :: iyr,syr,eyr
  integer,intent(in) :: im,jm
  
  real(kind = 8),intent(in) :: dat(im,jm)
  real(kind = 8),intent(inout) :: mean(im,jm)

  if(iyr == syr)then
     mean(:,:)=0.d0
  end if

  mean(:,:)=mean(:,:)+dat(:,:)

  if(iyr == eyr)then
     mean(:,:)=mean(:,:)/dble(eyr-syr+1)
  end if

end subroutine clim

!---------------------------------------------------------------
! Write Data |
!---------------------------------------------------------------

subroutine write_time(ntime,syr,smon,sday,dt,im,jm)

  use type_setting, only: lswitch_remove
  use mod_check_netcdf
  use netcdf
  implicit none

  integer status,system,access
  integer ncid,varid
  integer itime
  
  real(kind = 4) tmp1d(ntime)

  real(kind = 8) ffluxtime(ntime)
  
  character(100) :: filename="../in/river_clim.nc"
  
  !IN  
  integer,intent(in) :: ntime
  integer,intent(in) :: syr,smon,sday
  integer,intent(in) :: dt
  integer,intent(in) :: im,jm

  !ffluxtime
  do itime=1,ntime
     ffluxtime(itime)=dble((itime-1)*dt)/24.d0
  end do
  
  !Makefile
  if(lswitch_remove)then
     status=system("rm -f "//trim(filename))
  end if

  status=access(trim(filename)," ")
  if(status /= 0)then
     call make_ncfile_fflux(im,jm,syr,smon,sday,filename)
  end if
  
  !Open netcdf file
  status=nf90_open(trim(filename),nf90_write,ncid)
  call check_netcdf(status)

  !Write atmtime  
  status=nf90_inq_varid(ncid,"rivtime",varid)
  call check_netcdf(status)
  tmp1d(:)=real(ffluxtime(:))
  status=nf90_put_var(ncid,varid,tmp1d,(/1/),(/ntime/))
  call check_netcdf(status)
  
  !Close netcdf file
  status=nf90_close(ncid)
  call check_netcdf(status)
  
end subroutine write_time

!----------------------------

subroutine write_data(syr,smon,sday,shour,iyr,imon,iday,ihour,dt,im,jm,river)

  use netcdf
  use mod_check_netcdf
  implicit none

  integer itime  
  integer status,ncid,varid

  real(kind = 4) tmp2d(im,jm)
  
  character(100) :: filename="../in/river_clim.nc"
  character(2) mm,dd
  
  !IN
  integer,intent(in) :: syr,smon,sday,shour
  integer,intent(in) :: iyr,imon,iday,ihour
  integer,intent(in) :: dt
  integer,intent(in) :: im,jm

  real(kind = 8),intent(in) :: river(im,jm)

  !Calculate itime
  call estimate_ntime(syr,smon,sday,shour,iyr,imon,iday,ihour,dt,itime)
  
  !Open netcdf file
  status=nf90_open(trim(filename),nf90_write,ncid)
  call check_netcdf(status)
  
  !Write river
  tmp2d(:,:)=real(river(:,:))
  status=nf90_inq_varid(ncid,"river",varid)
  call check_netcdf(status)
  status=nf90_put_var(ncid,varid,tmp2d,(/1,1,itime/),(/im,jm,1/))
  call check_netcdf(status)
  
  status=nf90_close(ncid)
  call check_netcdf(status)
  
end subroutine write_data
