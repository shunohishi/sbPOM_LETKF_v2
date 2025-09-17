module setting

  !---Temporal range
  integer,parameter :: syr=2021,smon=1,sday=1,shour=0 !Start date
  integer,parameter :: eyr=2021,emon=1,eday=2,ehour=0 !End date
  integer,parameter :: dt=24

  !---Switch (1:On,Other:Off)
  integer,parameter :: iswitch_sst=1
  integer,parameter :: iswitch_sss=1
  integer,parameter :: iswitch_ssh=1
  integer,parameter :: iswitch_ts=1
  integer,parameter :: iswitch_ssuv=0 !*Motion vector
  
  !---Observation range & Error  
  real(kind = 8),parameter :: ssu_min=-5.d0,  ssu_max=5.d0,  ssu_err=0.2d0 !SSU
  real(kind = 8),parameter :: ssv_min=-5.d0,  ssv_max=5.d0,  ssv_err=0.2d0 !SSV
  real(kind = 8),parameter :: sst_min=-1.8d0, sst_max=50.d0, sst_err=1.d0  !SST
  real(kind = 8),parameter :: sss_min=25.d0,  sss_max=50.d0, sss_err=0.2d0 !SSS
  real(kind = 8),parameter :: ssh_min=-5.d0,  ssh_max=5.d0,  ssh_err=0.2d0 !SSH
  real(kind = 8),parameter :: u_min=-5.d0,  u_max=5.d0,  u_err=0.2d0       !Zonal velocity
  real(kind = 8),parameter :: v_min=-5.d0,  v_max=5.d0,  v_err=0.2d0       !Meridional velocity
  real(kind = 8),parameter :: t_min=-1.8d0, t_max=50.d0, t_err=1.d0        !Temperature
  real(kind = 8),parameter :: s_min=25.d0,  s_max=50.d0, s_err=0.2d0       !Salinity

  !---Element ID
  integer,parameter :: ssu_ele=2819,u_ele=2819
  integer,parameter :: ssv_ele=2820,v_ele=2820
  integer,parameter :: sst_ele=3073,t_ele=3073
  integer,parameter :: sss_ele=3332,s_ele=3332
  integer,parameter :: ssh_ele=2567

  !---Number of satellite SST
  integer,parameter :: nsst=4 !AMSR-E, WindSat, AMSR2, Himawari-8

  !---Switch Microwave SST
  integer,parameter :: iswitch_msst=1 !1: On, 0: Off

  !---Switch Himawari SST
  integer,parameter :: iswitch_hsst=0 !1: On, 0: Off
  
  !---Depth limit for SSH assimilation
  real(kind = 8),parameter :: ssh_depth=200.d0 ![m]

  !---Start & End Year to calculate model mean dynamical ocean topography 
  integer,parameter :: syr_ssh=2020,eyr_ssh=2020 ![year]

  !---Low Chl-a limit
  real(kind = 8),parameter :: chla_limit=0.1d0
  
  !---Initial filename: e.g. "../obs/${ini_filename}20110101.nc"
  character(100),parameter :: ini_filename="obs"
  
  !---POM filename
  character(20),parameter :: pom_dirname="pom90-ens"
  character(10),parameter :: pom_filename="test"

end module setting

!-------------------------------------------------------------
! Make observatioin data for LETKF |
!-----------------------------------
!
! 1. Estimate ntime
! 2. Read model grid
! 3. Make MDOT
! 4. Satellite SST if iswitch_sst=1
! 5. Satellite SSS if iswitch_sss=1 
! 6. Satellite SSH (Satellite SSHA + Simulated MDOT) if iswitch_ssh=1
! 7. In-situ Temperature/Salinity if iswitch_ts=1
!
!--------------------------------------
!
! Created                   by S.Ohishi 2018.09
! Modified                  by S.Ohishi 2019.10
! Added nearshore scheme    by S.Ohishi 2020.03
! Added openmp              by S.Ohishi 2020.04
! Added Microwave satellite by S.Ohishi 2024.04
! Modified                  by S.Ohishi 2025.07
!
!-------------------------------------------------------------

program main

  use setting
  use mod_rmiss
  use mod_gridinfo
  implicit none
  
  !---Common
  integer itime,ntime
  integer,allocatable :: iyr(:),imon(:),iday(:),ihour(:)
  integer inum
  integer status,system
  character(8) yyyymmdd
  character(4) yyyy
  character(2) mm,dd,hh

  !---Model
  real(kind = 8) lon(im),lat(jm),depth(im,jm,km) !Grid information
  real(kind = 8) mdot(im,jm)                     !Mean Dynamic Ocean Topography [m]
  real(kind = 8) fsm(im,jm)                      !Land-Sea mask (1:Sea, 0: Land)
  real(kind = 8) null1x(im),null1y(jm),null2d(im,jm),null3d(im,jm,km)

  !---Estimate ntime
  call estimate_ntime(syr,smon,sday,shour,eyr,emon,eday,ehour,dt,ntime)

  write(*,'(a)') "=== Make observational data for Assimilation ==="
  write(*,'(a,i4,i2,i2,i2)') "Start date: ",syr,smon,sday,shour
  write(*,'(a,i4,i2,i2,i2)') "End date: ",eyr,emon,eday,ehour
  write(*,'(a,i10)') "Number of timestep: ",ntime

  !---Read model grid data
  write(*,*) "Read model grid"
  call read_grid(&
       & null1x,null1x,lon,null1x, &
       & null1y,null1y,lat,null1y, & 
       & null3d,depth, &
       & fsm,null2d,null2d)

  !---Prepare mdot
  if(iswitch_ssh == 1)then
     write(*,'(a)') "===Make MDOT (Mean Dynamical Ocean Topography)"
     call prepare_mdot(lon,lat,mdot)
  end if

  !---set iyr,imon,iday,ihour
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

     !---Initialization
     write(yyyy,'(i4.4)') iyr(itime)
     write(mm,'(i2.2)') imon(itime)
     write(dd,'(i2.2)') iday(itime)
     write(hh,'(i2.2)') ihour(itime)
     yyyymmdd=yyyy//mm//dd
     inum=0

     !---Start
     write(*,'(a)') "-----Start "//yyyy//mm//dd//hh//"-----"

     !---Remove file
     status=system("rm -f ../obs/"//trim(ini_filename)//yyyymmdd//".nc")
     status=system("rm -f ../obs/ts"//yyyymmdd//".nc")

     !---SST
     if(iswitch_sst == 1)then
        write(*,'(a)') "===Start: Satellite SST observation"
        call prepare_sst(iyr(itime),imon(itime),iday(itime),inum,lon,lat,fsm)
        write(*,'(a)') "===End: Satellite SST observation"
     else
        write(*,'(a)') "===Skip: Satellite SST observation"
     end if

     !---SSS
     if(iswitch_sss == 1)then
        write(*,'(a)') "===Start: Satellite SSS observation"
        call prepare_sss(iyr(itime),imon(itime),iday(itime),inum,lon,lat,fsm)
        write(*,'(a)') "===End: Satellite SSS observation"
     else
        write(*,'(a)') "===Skip: Satellite SSS observation"
     end if

     !---SSH
     if(iswitch_ssh == 1)then
        write(*,'(a)') "===Start: Satellite SSH observation"
        call prepare_ssh(iyr(itime),imon(itime),iday(itime),inum,lon,lat,depth,mdot,fsm)
        write(*,'(a)') "===End: Satellite SSH observation"
     else
        write(*,'(a)') "===Skip: Satellite SSH observation"
     end if

     !---T/S
     if(iswitch_ts == 1)then
        write(*,'(a)') "===Start: In-situ T/S observation"
        call prepare_ts(iyr(itime),imon(itime),iday(itime),inum,lon,lat,depth,fsm)
        write(*,'(a)') "===End: In-situ T/S observation"
     else
        write(*,'(a)') "===Skip: In-situ T/S observation"
     end if

     !---SSU/V
     if(iswitch_ssuv == 1)then
        write(*,'(a)') "===Start: SSU/V based on motion vector"
        call prepare_ssuv(iyr(itime),imon(itime),iday(itime),inum,lon,lat,depth,fsm)
        write(*,'(a)') "===End: SSU/V based on motion vector"
     else
        write(*,'(a)') "===Skip: SSU/V based on motion vector"
     end if

     !---End
     write(*,'(a,i10)') "===Total number of observations: ",inum
     write(*,'(a)') "-----End "//yyyy//mm//dd//hh//"-----"

  end do !itime

  deallocate(iyr,imon,iday,ihour)

end program main
