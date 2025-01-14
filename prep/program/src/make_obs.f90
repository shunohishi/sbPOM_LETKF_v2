module setting

  !Temporal range
  integer,parameter :: syr=2002,smon=6,sday=1,shour=0 !Start date
  integer,parameter :: eyr=2002,emon=6,eday=2,ehour=0 !End date
  integer,parameter :: dt=24

  !Switch (1:On,Other:Off)
  integer,parameter :: iswitch_ssuv=0
  integer,parameter :: iswitch_sst=1
  integer,parameter :: iswitch_sss=1
  integer,parameter :: iswitch_ssh=1
  integer,parameter :: iswitch_ts=1
  
  !Observation range & Error
  real(kind = 8),parameter :: ssu_min=-5.d0,  ssu_max=5.d0,  ssu_err=0.2d0, ssu_ele=2819.d0  !SSU
  real(kind = 8),parameter :: ssv_min=-5.d0,  ssv_max=5.d0,  ssv_err=0.2d0, ssv_ele=2820.d0  !SSV
  real(kind = 8),parameter :: sst_min=-1.8d0, sst_max=50.d0, sst_err=1.d0,  sst_ele=3073.d0 !SST
  real(kind = 8),parameter :: sss_min=25.d0,  sss_max=50.d0, sss_err=0.2d0, sss_ele=3332.d0 !SSS
  real(kind = 8),parameter :: ssh_min=-5.d0,  ssh_max=5.d0,  ssh_err=0.2d0, ssh_ele=2567.d0  !SSH
  real(kind = 8),parameter :: u_min=-5.d0,  u_max=5.d0,  u_err=0.2d0, u_ele=2819.d0  !Zonal vel.
  real(kind = 8),parameter :: v_min=-5.d0,  v_max=5.d0,  v_err=0.2d0, v_ele=2820.d0  !Meridional vel.
  real(kind = 8),parameter :: t_min=-1.8d0, t_max=50.d0, t_err=1.d0,  t_ele=3073.d0 !Temperature
  real(kind = 8),parameter :: s_min=25.d0,  s_max=50.d0, s_err=0.2d0, s_ele=3332.d0 !Salinity
  
  !Grid interval
  integer,parameter :: ssu_dx=1,ssu_dy=1
  integer,parameter :: ssv_dx=1,ssv_dy=1
  integer,parameter :: sst_dx=1,sst_dy=1
  integer,parameter :: sss_dx=1,sss_dy=1
  integer,parameter :: ssh_dx=1,ssh_dy=1

  !Number of satellite SST
  integer,parameter :: nsst=4 !AMSR-E, WindSat, AMSR2, Himawari-8

  !Switch Microwave SST
  integer,parameter :: iswitch_msst=1 !1: On, 0: Off

  !Switch Himawari SST
  integer,parameter :: iswitch_hsst=0 !1: On, 0: Off
  
  !Exclude Nearshore satellite SSS
  integer,parameter :: iswitch_nearshore=1 !1:On, 0:Off
  real(kind = 8),parameter :: nearshore_range=100.d3 ![m]

  !Depth limit for SSH assimilation
  real(kind = 8),parameter :: ssh_depth=200.d0 ![m]

  !Start & End Year to calculate model mean dynamical ocean topography 
  integer,parameter :: syr_ssh=1996,eyr_ssh=2001 ![year]
  
  !Initial filename: e.g. "../obs/${ini_filename}20110101.nc"
  character(100),parameter :: ini_filename="obs"
  
  !POM filename
  character(20),parameter :: pom_dirname="pom90-ens"
  character(10),parameter :: pom_filename="test"

end module setting

!-------------------------------------------------------------
! Make observatioin data for LETKF |
!-----------------------------------
!
! 1. Estimate ntime
! 2. Read model grid (& Detect nearshore grid)
! 3. Prepare Satellite SST, if iswitch_sst=1
! 4. Prepare Satellite SSS, if iswitch_sss=1
! 5. Prepare Mean Dynamic Ocean Topography, if iswitch_ssh=1
! 6. Prepare Satellite SSH + MDOT, if iswitch_ssh=1
! 7. Prepare Temperature/Salinity from GTSPP and AQC Argo, if iswitch_ts=1
!
!--------------------------------------
!
! Created                   by S.Ohishi 2018.09
! Modified                  by S.Ohishi 2019.10
! Added nearshore scheme    by S.Ohishi 2020.03
! Added openmp              by S.Ohishi 2020.04
! Added Microwave satellite by S.Ohishi 2024.04
!
!-------------------------------------------------------------

program main

  !$ use omp_lib
  use setting
  use mod_rmiss
  use mod_gridinfo
  implicit none
  
  !Common
  integer itime,ntime
  integer,allocatable :: iyr(:),imon(:),iday(:),ihour(:)
  integer inum
  integer status,system
  character(8) yyyymmdd
  character(4) yyyy
  character(2) mm,dd,hh

  !Model
  real(kind = 8) lon(im),lat(jm),depth(im,jm,km)
  real(kind = 8) mdot(im,jm) !Mean Dynamic Ocean Topography [m]
  real(kind = 8) fsm(im,jm),nshore(im,jm)
  real(kind = 8) null1x(im),null1y(jm),null2d(im,jm),null3d(im,jm,km)

  !Estimate ntime
  call estimate_ntime(syr,smon,sday,shour,eyr,emon,eday,ehour,dt,ntime)

  write(*,'(a)') "=== Make observational data for Assimilation ==="
  write(*,'(a,i4,i2,i2,i2)') "Start date: ",syr,smon,sday,shour
  write(*,'(a,i4,i2,i2,i2)') "End date: ",eyr,emon,eday,ehour
  write(*,'(a,i10)') "Number of timestep: ",ntime

  !Read model grid data
  write(*,*) "Read model grid"
  call read_grid(&
       & null1x,null1x,lon,null1x, &
       & null1y,null1y,lat,null1y, & 
       & null3d,depth, &
       & fsm,null2d,null2d)

  !Detect nearshore
  if(iswitch_nearshore == 1)then
     call detect_nearshore(lon,lat,fsm,nshore)
  else
     nshore(:,:)=1.d0
  end if

  !Prepare mdot
  if(iswitch_ssh == 1)then
     write(*,'(a)') "===Prepare MDOT"
     call prepare_mdot(lon,lat,mdot)
  end if

  !set iyr,imon,iday,ihour
  allocate(iyr(ntime),imon(ntime),iday(ntime),ihour(ntime))
  iyr(1)=syr
  imon(1)=smon
  iday(1)=sday
  ihour(1)=shour
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
     yyyymmdd=yyyy//mm//dd
     inum=0

     write(*,'(a)') "-----Start "//yyyy//mm//dd//hh//"-----"
     status=system("rm -f ../obs/"//trim(ini_filename)//yyyymmdd//".nc")
     status=system("rm -f ../obs/ts"//yyyymmdd//".nc")
     
     if(iswitch_sst == 1)then
        write(*,'(a)') "===Prepare SST"
        call prepare_sst(iyr(itime),imon(itime),iday(itime),inum,lon,lat,fsm)
     else
        write(*,'(a)') "===Skip to prepare SST"
     end if
        
     if(iswitch_sss == 1)then
        write(*,'(a)') "===Prepare SSS"
        call prepare_sss(iyr(itime),imon(itime),iday(itime),inum,lon,lat,fsm,nshore)
     else
        write(*,'(a)') "===Skip to prepare SSS"
     end if

     if(iswitch_ssh == 1)then
        write(*,'(a)') "===Prepare SSH"
        call prepare_ssh(iyr(itime),imon(itime),iday(itime),inum,lon,lat,depth,mdot,fsm)
     else
        write(*,'(a)') "===Skip to prepare SSH"
     end if

     if(iswitch_ts == 1)then
        write(*,'(a)') "===Prepare T/S"
        call prepare_ts(iyr(itime),imon(itime),iday(itime),inum,lon,lat,depth,fsm)
     else
        write(*,'(a)') "===Skip to prepare T/S"
     end if

     if(iswitch_ssuv == 1)then
        write(*,'(a)') "===Prepare SSU/V"
        call prepare_ssuv(iyr(itime),imon(itime),iday(itime),inum,lon,lat,depth,fsm)
     else
        write(*,'(a)') "===Skip to prepare U/V"
     end if


     write(*,'(a,i10)') "===The number of observation: ",inum
     write(*,'(a)') "-----End "//yyyy//mm//dd//hh//"-----"

  end do

  deallocate(iyr,imon,iday,ihour)

end program main
