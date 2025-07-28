!--------------------------------------------------------------
! Prepare temperature/salinity data |
!--------------------------------------------------------------
!
! 1. Read AQC Argo
! 2. Make GTSPP filename
! 3. Read GTSPP
! 4. Check duplication between GTSPP and AQC Argo
! 5. Write GTSPP 
!
! 6. Write AQC Argo not duplicating GTSPP
!
!-------------------------------------------------------------
!
! *Although sbPOM uses "potential temperature",
!  "temperature" is used for assimilation, 
!  because the difference between potential temperature and temperature is small.
!  (0.1-0.2 degree C in 1000m)
!
!-------------------------------------------------------------
!
! Created by S.Ohishi 2018.09
! Modified by S.Ohishi 2020.10
! Modified by S.Ohishi 2022.09 
!  - Add no AQC Argo case
! Modified by S.Ohishi 2023.01
! Modified by S.Ohishi 2025.07
!  - Change from GTSPP+AQC to GTSPP+EN4
!  - Add instrument information (XBT, XCTD, Argo, and Mooring buoy etc.)
!    * Prioritize EN4 over GTSPP
!--------------------------------------------------------------

subroutine prepare_ts(iyr,imon,iday,inum,lon,lat,depth,fsm)
  
  use mod_gridinfo
  use mod_read_en4
  use mod_read_gtspp
  implicit none

  !---Common
  integer idataset(1)
  integer its
  integer flag
  
  !---EN4
  integer ip
  integer np_en !Number of profile
  integer km_en !Number of level

  integer,allocatable :: iday_en(:)
  integer,allocatable :: ins_en(:)
    
  real(kind = 8),allocatable :: lon_en(:),lat_en(:),dep_en(:,:)
  real(kind = 8),allocatable :: t_en(:,:),s_en(:,:)
  
  !---GTSPP
  integer ifile
  integer nfile   !Number of file
  integer km_gts  !Number of level
  integer ins_gts !Instrument
  
  real(kind = 8),allocatable :: lon_gts(:),lat_gts(:),depth_gts(:)
  real(kind = 8),allocatable :: t_gts(:),s_gts(:)
  character(100),allocatable :: filename(:)

  !---IN
  integer,intent(in) :: iyr,imon,iday
  real(kind = 8),intent(in) :: lon(im),lat(jm),depth(im,jm,km)
  real(kind = 8),intent(in) :: fsm(im,jm)

  !---INOUT
  integer,intent(inout) :: inum

  !---Initialization
  its=0

  !---Read EN4
  write(*,*) "Read EN4"
  call read_en4(iyr,imon,iday_en,np_en,km_en,lon_en,lat_en,dep_en,t_en,s_en,ins_en)

  !---Write EN4
  write(*,*) "Write EN4"
  idataset(1)=1
  do ip=1,np_en

     !Different date
     if(iday /= iday_en(ip)) cycle     
     !Missing value due to QC
     if(lon_en(ip) == 0.d0 .and. lat_en(ip) == 0.d0 .and. ins_en(ip) == 0) cycle
     !Check position
     if(lon_en(ip) < lon(1) .or. lon(im) < lon_en(ip) .or. &
          & lat_en(ip) < lat(1) .or. lat(jm) < lat_en(ip)) cycle
     
     !Write profile information
     call write_ts_info(idataset,iyr,imon,iday_en(ip),its,km_en,lon_en(ip),lat_en(ip),dep_en(:,ip))

     !Write Temperature and Salinity
     call write_obs_ts("T",iyr,imon,iday,inum,km_en, &
          & lon(1),lon(im),lon_en(ip),lat(1),lat(jm),lat_en(ip),dep_en(:,ip), &
          & t_en(:,ip),ins_en(ip))
     call write_obs_ts("S",iyr,imon,iday,inum,km_en, &
          & lon(1),lon(im),lon_en(ip),lat(1),lat(jm),lat_en(ip),dep_en(:,ip), &
          & s_en(:,ip),ins_en(ip))
     
  end do !ip

  !---Make/Read GTSPP filename
  write(*,'(a)') "Make GTSPP filename"
  call make_filename(iyr,imon)
  write(*,'(a)') "Read GTSPP filename"
  call read_filename(iyr,imon,iday,nfile,lon_gts,lat_gts,filename)

  !Nofile
  if(nfile == 0)then
     write(*,'(a)') "Skip GTSPP because of no file"
     return
  end if

  !---Write GTSPP
  write(*,*) "Read and Write GTSPP: The number of file: ",nfile
  idataset(1)=2
  do ifile=1,nfile

     if(mod(ifile,500) == 1) write(*,*) ifile,"/",nfile

     !Check position
     if(lon_gts(ifile) < lon(1) .or. lon(im) < lon_gts(ifile) &
          & .or. lat_gts(ifile) < lat(1) .or. lat(jm) < lat_gts(ifile)) cycle

     !Match up between EN4 and GTSPP
     !(Without overlap => flag=1, Overlap => flag=0)
     call match_up_en4_gtspp(np_en,iday_en,lon_en,lat_en,iday,lon_gts(ifile),lat_gts(ifile),flag)

     if(flag == 0) cycle
     
     !Read GTSPP data
     call read_gtspp(iyr,imon,iday,filename(ifile), &
          & km_gts,lon_gts(ifile),lat_gts(ifile),depth_gts,t_gts,s_gts)
     call read_gtspp_obs_type(iyr,imon,filename(ifile),ins_gts)
     
     !Write GTSPP
     call write_ts_info(idataset,iyr,imon,iday,its,km_gts,lon_gts(ifile),lat_gts(ifile),depth_gts)
     call write_obs_ts("T",iyr,imon,iday,inum,km_gts, &
          & lon(1),lon(im),lon_gts(ifile),lat(1),lat(jm),lat_gts(ifile),depth_gts, &
          & t_gts,ins_gts)
     call write_obs_ts("S",iyr,imon,iday,inum,km_gts, &
          & lon(1),lon(im),lon_gts(ifile),lat(1),lat(jm),lat_gts(ifile),depth_gts, &
          & s_gts,ins_gts)

     call deallocate_gtspp(depth_gts,t_gts,s_gts)

  end do !ifile
     
  call deallocate_gtspp_filename(lon_gts,lat_gts,filename)
  call deallocate_en4(iday_en,lon_en,lat_en,dep_en,t_en,s_en,ins_en)
  
end subroutine prepare_ts

!----------------------------------------------------------------------
! Match up between EN4 and GTSPP |
!----------------------------------------------------------------------

subroutine match_up_en4_gtspp(np_en,iday_en,lon_en,lat_en,iday,lon_gts,lat_gts,flag)

  implicit none

  !---Parameter
  real(kind = 8),parameter :: dr=0.01d0
  
  !---Common
  integer i
  
  !---IN
  !EN4
  integer,intent(in) :: np_en
  integer,intent(in) :: iday_en(np_en)
  
  real(kind = 8),intent(in) :: lon_en(np_en),lat_en(np_en)

  !GTSPP
  integer,intent(in) :: iday

  real(kind = 8),intent(in) :: lon_gts,lat_gts

  !---OUT
  integer,intent(out) :: flag

  flag=1
  do i=1,np_en

     if(iday /= iday_en(i)) cycle
     if(abs(lon_en(i)-lon_gts) < dr .and. abs(lat_en(i)-lat_gts) < dr)then
        flag=0
        exit
     end if
     
  end do
  
end subroutine match_up_en4_gtspp

!----------------------------------------------------------------------
! Write Observation information |
!----------------------------------------------------------------------

subroutine write_ts_info(idataset,iyr,imon,iday,its,km,lon,lat,depth)

  use mod_rmiss
  use netcdf
  implicit none

  !---Common
  integer status,access,system
  integer ncid,varid
  integer k

  character(100) filename
  character(4) yyyy
  character(2) mm,dd
  
  !---IN
  integer,intent(in) :: idataset(1)
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: km

  real(kind = 8),intent(in) :: lon(1),lat(1),depth(km)

  !---INOUT
  integer,intent(inout) :: its

  !---Count up its
  its=its+1

  !---Filename
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday

  filename="../obs/ts"//yyyy//mm//dd//".nc"

  !---Access
  status=access(trim(filename)," ")
  if(status /= 0)then
     call make_ncfile_ts(filename)
  end if

  !---NetCDF
  status=nf90_open(trim(filename),nf90_write,ncid)

  status=nf90_inq_varid(ncid,"lon",varid)
  status=nf90_put_var(ncid,varid,real(lon),start=(/its/),count=(/1/))

  status=nf90_inq_varid(ncid,"lat",varid)
  status=nf90_put_var(ncid,varid,real(lat),start=(/its/),count=(/1/))

  status=nf90_inq_varid(ncid,"levt",varid)
  status=nf90_put_var(ncid,varid,real(depth(1:1)),start=(/its/),count=(/1/))

  do k=km,1,-1
     if(depth(k) == rmiss) cycle
     status=nf90_inq_varid(ncid,"levb",varid)
     status=nf90_put_var(ncid,varid,real(depth(k:k)),start=(/its/),count=(/1/))
     exit
  end do
  
  status=nf90_inq_varid(ncid,"dat",varid)
  status=nf90_put_var(ncid,varid,idataset,start=(/its/),count=(/1/))
  
  status=nf90_close(ncid)

end subroutine write_ts_info
