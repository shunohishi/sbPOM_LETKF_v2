!--------------------------------------------------------------
! Prepare temperature/salinity data |
!--------------------------------------------------------------
!
! 1. Read AQC Argo
!
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
!
!--------------------------------------------------------------

subroutine prepare_ts(iyr,imon,iday,inum,lon,lat,depth,fsm)
  
  use setting,only: t_err,s_err
  use mod_gridinfo
  use mod_read_gtspp
  use mod_read_aqc_argo
  implicit none

  !Common
  integer its

  !Model
  real(kind = 8) dx,dy

  !GTSPP
  integer ifile,nfile
  integer km_gts
  integer flag_gts
  real(kind = 8),allocatable :: lon_gts(:),lat_gts(:)
  real(kind = 8),allocatable :: depth_gts(:),t_gts(:),s_gts(:)
  real(kind = 8),allocatable :: terr_gts(:),serr_gts(:)
  character(100),allocatable :: filename(:)

  !AQC Argo
  integer,allocatable :: iday_aqc(:)
  integer,allocatable :: flag_aqc(:) 
  !0: match GTSPP(no write), 1:not match GTSPP (write)
  integer ipro_aqc,npro_aqc,km_aqc
  real(kind = 8),allocatable :: lon_aqc(:),lat_aqc(:),depth_aqc(:,:)
  real(kind = 8),allocatable :: t_aqc(:,:),s_aqc(:,:)
  real(kind = 8),allocatable :: terr_aqc(:),serr_aqc(:)

  !IN
  integer,intent(in) :: iyr,imon,iday

  !OUT
  real(kind = 8),intent(in) :: lon(im),lat(jm),depth(im,jm,km)
  real(kind = 8),intent(in) :: fsm(im,jm)

  !INOUT
  integer,intent(inout) :: inum
  
  its=0

  !Resolution
  dx=lon(2)-lon(1)
  dy=lat(2)-lat(1)

  write(*,*) "Read AQC Argo"
  call read_aqc_argo(iyr,imon,iday_aqc,npro_aqc,km_aqc,&
       & lon_aqc,lat_aqc,depth_aqc,t_aqc,s_aqc)

  if(0 < npro_aqc)then
     allocate(flag_aqc(npro_aqc))
     allocate(terr_aqc(km_aqc),serr_aqc(km_aqc))
     flag_aqc(:)=0
  endif

  !---GTSPP
  write(*,*) "---GTSPP---"
  !Make GTSPP filename.txt
  call make_filename(1,1,iyr,imon)
  call read_filename(iyr,imon,iday,nfile,lon_gts,lat_gts,filename)

  if(nfile == 0)then
     write(*,'(a)') "Skip GTSPP"
  else
     
     do ifile=1,nfile
        
        if(mod(ifile,500) == 1) write(*,*) ifile,"/",nfile

        if(lon_gts(ifile) < lon(1)-dx*0.5d0 &
             & .or. lon(im)+dx*0.5d0 < lon_gts(ifile) &
             & .or. lat_gts(ifile) < lat(1)-dy*0.5d0 &
             & .or. lat(jm)+dy*0.5d0 < lat_gts(ifile))then
           cycle
        end if

        call read_gtspp(iyr,imon,iday,filename(ifile), &
             & km_gts,lon_gts(ifile),lat_gts(ifile),depth_gts,t_gts,s_gts)

        !Check AQC_Argo vs. GTSPP
        if(0 < npro_aqc)then
           call check_aqc_argo(iyr,imon,iday,lon_gts(ifile),lat_gts(ifile),&
                & npro_aqc,iday_aqc,lon_aqc,lat_aqc,flag_aqc) !1:Not match, 0: Match
        end if

        !TS Error
        allocate(terr_gts(km_gts),serr_gts(km_gts))
        terr_gts(:)=t_err
        serr_gts(:)=s_err

        !Write GTSPP
        call write_ts_info(1,iyr,imon,iday,its,km_gts,lon_gts(ifile),lat_gts(ifile),depth_gts)
        call write_obs_ts("T",iyr,imon,iday,inum,km_gts, &
             & lon_gts(ifile),lat_gts(ifile),depth_gts,t_gts,terr_gts)
        call write_obs_ts("S",iyr,imon,iday,inum,km_gts, &
             & lon_gts(ifile),lat_gts(ifile),depth_gts,s_gts,serr_gts)
        
        call deallocate_gtspp(depth_gts,t_gts,s_gts)
        deallocate(terr_gts,serr_gts)
        
     end do !ifile
     
     call deallocate_gtspp_filename(lon_gts,lat_gts,filename)
     
  end if

  !AQC Argo
  if(0 < npro_aqc)then 
     write(*,'(a)') "---AQC Argo---"
     do ipro_aqc=1,npro_aqc

        if(flag_aqc(ipro_aqc) == 0) cycle !Match with GTSPP

        !Not Match with GTSPP
        if(lon_aqc(ipro_aqc) < lon(1)-dx*0.5d0 &
             & .or. lon(im)+dx*0.5d0 < lon_aqc(ipro_aqc) &
             & .or. lat_aqc(ipro_aqc) < lat(1)-dy*0.5d0 &
             & .or. lat(jm)+dy*0.5d0 < lat_aqc(ipro_aqc)) cycle
     
        !TSError
        terr_aqc(:)=t_err
        serr_aqc(:)=s_err
        
        call write_ts_info(2,iyr,imon,iday,its,km_aqc,&
             & lon_aqc(ipro_aqc),lat_aqc(ipro_aqc),depth_aqc)
        call write_obs_ts("T",iyr,imon,iday,inum,km_aqc,&
             & lon_aqc(ipro_aqc),lat_aqc(ipro_aqc),depth_aqc(:,ipro_aqc), &
             & t_aqc(:,ipro_aqc),terr_aqc(:))
        call write_obs_ts("S",iyr,imon,iday,inum,km_aqc,&
             & lon_aqc(ipro_aqc),lat_aqc(ipro_aqc),depth_aqc(:,ipro_aqc), &
             & s_aqc(:,ipro_aqc),serr_aqc(:))
        
     end do !ipro_aqc
  end if

  if(0 < npro_aqc)then 
     deallocate(flag_aqc)
     deallocate(terr_aqc,serr_aqc)
     call deallocate_aqc_argo(iday_aqc,lon_aqc,lat_aqc,depth_aqc,t_aqc,s_aqc)
  end if

end subroutine prepare_ts

!--------------------------------------------------------------------
! Check AQC Argo matching with GTSPP | 
!--------------------------------------------------------------------

subroutine check_aqc_argo(iyr,imon,iday,lon_gts,lat_gts,&
          & npro_aqc,iday_aqc,lon_aqc,lat_aqc,flag_aqc)

  implicit none


  !Common  
  real(kind = 8),parameter :: dr=0.1d0 !Search range

  integer ipro_aqc

  !---IN
  !GTSPP
  integer,intent(in) :: iyr,imon,iday
  real(kind = 8),intent(in) :: lon_gts,lat_gts
  
  !AQC Argo
  integer,intent(in) :: npro_aqc
  integer,intent(in) :: iday_aqc(npro_aqc)
  real(kind = 8),intent(in) :: lon_aqc(npro_aqc),lat_aqc(npro_aqc)

  !---INOUT
  integer,intent(inout) :: flag_aqc(npro_aqc)

  do ipro_aqc=1,npro_aqc

     if(iday /= iday_aqc(ipro_aqc))cycle

     if(lon_gts-dr <= lon_aqc(ipro_aqc) &
          & .and. lon_aqc(ipro_aqc) <= lon_gts+dr &
          & .and. lat_gts-dr <= lat_aqc(ipro_aqc) &
          & .and. lat_aqc(ipro_aqc) <= lat_gts+dr)then
        flag_aqc(ipro_aqc)=0 !Match with GTSPP
     else
        flag_aqc(ipro_aqc)=1 !Not Match with GTSPP
     end if

  end do !ipro_aqc

end subroutine check_aqc_argo

!----------------------------------------------------------------------
! Write Observation information |
!----------------------------------------------------------------------

subroutine write_ts_info(idataset,iyr,imon,iday,its,km,lon,lat,depth)

  use mod_rmiss
  use netcdf
  implicit none

  integer status,access,system
  integer ncid,varid
  integer k
  
  integer,intent(in) :: idataset(1)
  integer,intent(in) :: iyr,imon,iday
  integer,intent(inout) :: its
  integer,intent(in) :: km

  real(kind = 8),intent(in) :: lon(1),lat(1),depth(km)

  character(100) filename
  character(4) yyyy
  character(2) mm,dd

  its=its+1
  
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday

  filename="../obs/ts"//yyyy//mm//dd//".nc"

  status=access(trim(filename)," ")
  if(status /= 0)then
     call make_ncfile_ts(filename)
  end if

  status=nf90_open(trim(filename),nf90_write,ncid)

  status=nf90_inq_varid(ncid,"lon",varid)
  status=nf90_put_var(ncid,varid,real(lon),start=(/its/),count=(/1/))

  status=nf90_inq_varid(ncid,"lat",varid)
  status=nf90_put_var(ncid,varid,real(lat),start=(/its/),count=(/1/))

  status=nf90_inq_varid(ncid,"lev1",varid)
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
