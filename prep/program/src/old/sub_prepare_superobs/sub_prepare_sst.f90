!------------------------------------------------------
! Prepare SST
!------------------------------------------------------
!
! <GCOM-W/AMSR2>
! 1. Calculate GCOM-W/AMSR2 SST
! 2. Project GCOM-W/AMSR2 SST onto model space
!
! <Himawari-8 -9/AHI>
! 3. Calculate daily mean Himawari/AHI SST 
! 4. Project daily mean Himawari/AHI SST onto model space
! 
! <Modify SST & Write Data>
! 5. Apply fsm for both SST datasets
! 6. Merge Himawari and GCOM-W SST
! 7. Write Data
!
!------------------------------------------------------
! Created by S.Ohishi 2018.09
! Modified by S.Ohishi 2018.10
!          by S.Ohishi 2019.01
!          by S.Ohishi 2020.04
!          by S.Ohishi 2024.01
!          by S.Ohishi 2024.04
!------------------------------------------------------
subroutine prepare_sst(iyr,imon,iday,inum,lon,lat,fsm)
  use setting, only: iswitch_msst, iswitch_hsst
  use mod_rmiss
  use mod_gridinfo
  use mod_read_himawari, im_h => im, jm_h => jm
  use mod_read_amsre
  use mod_read_windsat
  use mod_read_amsr2
  implicit none

  integer ihour
  integer jyr,jmon,jday
  integer i,j

  !SST in model space
  real(kind = 8) asst(im,jm),acount(im,jm) !Aqua/AMSRE SST
  real(kind = 8) csst(im,jm),ccount(im,jm) !Coriolis/WindSat SST
  real(kind = 8) gsst(im,jm),gcount(im,jm) !GCOM-W/AMSR2 SST
  real(kind = 8) hsst(im,jm) !Himawari/AHI SST

  real(kind = 8) sst(im,jm) !SST

  !GCOM-W/AMSR2, Aqua/AMSRE, Coriolis/WindSat
  integer ifile,nfile  
  integer im_g,jm_g
  real(kind = 8),allocatable :: lon_g(:,:),lat_g(:,:),sst_g(:,:)
  character(200),allocatable :: filename(:)

  !Himawari/AHI
  real(kind = 8) lon_h(im_h),lat_h(jm_h)
  real(kind = 8) dat_h(im_h,jm_h)
  real(kind = 8) sst_h(im_h,jm_h),pass_h(im_h,jm_h),miss_h(im_h,jm_h)
  
  !IN
  integer,intent(in) :: iyr,imon,iday
  real(kind = 8),intent(in) :: lon(im),lat(jm)
  real(kind = 8),intent(in) :: fsm(im,jm)

  !IN/OUT
  integer,intent(inout) :: inum  

  if(iswitch_msst == 0)then

     asst(:,:)=rmiss
     csst(:,:)=rmiss
     gsst(:,:)=rmiss     
     
  else if(iswitch_msst == 1)then

     !---Aqua/AMSR-E SST
     write(*,*) "Calculate daily-mean Aqua/AMSR-E SST"
     asst(:,:)=0.d0
     acount(:,:)=0.d0

     !Present day
     call read_amsre_filename(iyr,imon,iday,nfile,filename)

     if(nfile == 0)then
        asst(:,:)=rmiss
     else
        do ifile=1,nfile
           call read_amsre(filename(ifile),jyr,jmon,jday,ihour,im_g,jm_g,lon_g,lat_g,sst_g)
           call project_gsst(1,im_g,jm_g,lon_g,lat_g,sst_g,&
                & im,jm,lon,lat,asst,acount)
           call deallocate_amsre(lon_g,lat_g,sst_g)
        end do
        call deallocate_amsre_filename(filename)

        call project_gsst(2,im_g,jm_g,lon_g,lat_g,sst_g,&
             & im,jm,lon,lat,asst,acount)
     end if
     
     !---Coriolis/WindSat SST
     write(*,*) "Calculate daily-mean Coriolis/WindSat SST"
     csst(:,:)=0.d0
     ccount(:,:)=0.d0

     !Present day
     call read_windsat_filename(iyr,imon,iday,nfile,filename)

     if(nfile == 0)then
        csst(:,:)=rmiss
     else
        do ifile=1,nfile
           call read_windsat(filename(ifile),jyr,jmon,jday,ihour,im_g,jm_g,lon_g,lat_g,sst_g)
           call project_gsst(1,im_g,jm_g,lon_g,lat_g,sst_g,&
                & im,jm,lon,lat,csst,ccount)
           call deallocate_windsat(lon_g,lat_g,sst_g)
        end do
        call deallocate_windsat_filename(filename)

        call project_gsst(2,im_g,jm_g,lon_g,lat_g,sst_g,&
             & im,jm,lon,lat,csst,ccount)
     end if

     !---GCOM-W/AMSR2 SST
     write(*,*) "Calculate daily-mean GCOM-W/AMSR2 SST"
     gsst(:,:)=0.d0
     gcount(:,:)=0.d0

     !Present day
     call read_amsr2_filename(iyr,imon,iday,nfile,filename)

     if(nfile == 0)then
        gsst(:,:)=rmiss
     else
        do ifile=1,nfile
           call read_amsr2(filename(ifile),jyr,jmon,jday,ihour,im_g,jm_g,lon_g,lat_g,sst_g)
           call project_gsst(1,im_g,jm_g,lon_g,lat_g,sst_g,&
                & im,jm,lon,lat,gsst,gcount)
           call deallocate_amsr2(lon_g,lat_g,sst_g)
        end do
        call deallocate_amsr2_filename(filename)

        call project_gsst(2,im_g,jm_g,lon_g,lat_g,sst_g,&
             & im,jm,lon,lat,gsst,gcount)
     end if

  end if
     
  !---Himawari daily SST
  if(iswitch_hsst == 0)then
     hsst(:,:)=rmiss
  else if(iswitch_hsst == 1)then
     write(*,'(a)') "Calculate daily-mean Himawari SST in obs. space"
     do ihour=0,23
        call read_himawari(iyr,imon,iday,ihour,lon_h,lat_h,dat_h)
        call cal_daily_mean(ihour,im_h,jm_h,dat_h,sst_h,pass_h,miss_h)
     end do
     write(*,'(a)') "Project Himawari SST in obs. space to model space"
     call project_hsst(im_h,jm_h,lon_h,lat_h,sst_h,im,jm,lon,lat,hsst)
  end if
     
  !Apply fsm
  !write(*,'(a)') "Apply fsm"
  call apply_fsm(im,jm,1,asst,fsm)
  call apply_fsm(im,jm,1,csst,fsm)
  call apply_fsm(im,jm,1,gsst,fsm)
  call apply_fsm(im,jm,1,hsst,fsm)

  !Calculate SST modified by bias
  write(*,*) "Merge SST"
  call merge_sst(im,jm,asst,csst,gsst,hsst,sst)
  call apply_fsm(im,jm,1,sst,fsm)

  !Write Data
  write(*,*) "Write Data"
  call write_sst(iyr,imon,iday,lon,lat,asst,csst,gsst,hsst,sst)
  call write_obs_surface("sst",iyr,imon,iday,inum,lon,lat,sst)

end subroutine prepare_sst

!----------------------------------------------------------------------
! Calculate Daily Mean |
!----------------------------------------------------------------------

subroutine cal_daily_mean(ihour,im,jm,dat,mean,pass,miss)

  use mod_rmiss
  implicit none

  !Common
  integer i,j

  !IN
  integer,intent(in) :: ihour  
  integer,intent(in) :: im,jm
  real(kind = 8),intent(in) :: dat(im,jm)

  !INOUT
  real(kind = 8),intent(inout) :: mean(im,jm),pass(im,jm),miss(im,jm)

  if(ihour == 0)then
     mean(:,:)=0.d0
     pass(:,:)=0.d0
     miss(:,:)=0.d0
  end if

  do j=1,jm
     do i=1,im
        if(dat(i,j) == rmiss)then
           miss(i,j)=miss(i,j)+1.d0
        else
           mean(i,j)=mean(i,j)+dat(i,j)
           pass(i,j)=pass(i,j)+1.d0
        end if
     end do
  end do

  if(ihour == 23)then
     do j=1,jm
        do i=1,im
           if(pass(i,j) == 0.d0)then
              mean(i,j)=rmiss
           else
              mean(i,j)=mean(i,j)/pass(i,j)
           end if
        end do
     end do
  end if
  
end subroutine cal_daily_mean

!-------------------------------------------------------------
! Himawari8/AHI SST in obs space --> model space
!-------------------------------------------------------------

subroutine project_hsst( &
     & im_h,jm_h,lon_h,lat_h,sst_h,&
     & im,jm,lon,lat,hsst)

  use mod_rmiss
  implicit none

  !Common
  integer iswitch
  integer ih,jh,im_h,jm_h
  integer i,j,im,jm

  real(kind = 8) dx,dy,hcount(im,jm)

  !IN
  real(kind = 8),intent(in) :: lon_h(im_h),lat_h(jm_h),sst_h(im_h,jm_h)
  real(kind = 8),intent(in) :: lon(im),lat(jm)

  !OUT
  real(kind = 8),intent(out) :: hsst(im,jm)

  dx=lon(2)-lon(1)
  dy=lat(2)-lat(1)

  hsst(:,:)=0.d0
  hcount(:,:)=0.d0

  do jh=1,jm_h
     if(lat_h(jh) < lat(1)-0.5d0*dy &
          & .or. lat(jm)+0.5d0*dy < lat_h(jh)) cycle
     do ih=1,im_h
        if(lon_h(ih) < lon(1)-0.5d0*dx &
             & .or. lon(im)+0.5d0*dx < lon_h(ih)) cycle
        if(sst_h(ih,jh) == rmiss) cycle

        do j=1,jm
           if(lat_h(jh) < lat(j)-0.5d0*dy &
                & .or. lat(j)+0.5d0*dy < lat_h(jh)) cycle
           do i=1,im
              if(lon_h(ih) < lon(i)-0.5d0*dx &
                   & .or. lon(i)+0.5d0*dx < lon_h(ih)) cycle
              hsst(i,j)=hsst(i,j)+sst_h(ih,jh)
              hcount(i,j)=hcount(i,j)+1.d0
           end do !i
        end do !j

     end do !ih
  end do !jh

  do j=1,jm
     do i=1,im
        if(hcount(i,j) == 0.d0)then
           hsst(i,j)=rmiss
        else
           hsst(i,j)=hsst(i,j)/hcount(i,j)
        end if
     end do
  end do
  
end subroutine project_hsst

!-------------------------------------------------------------
! GCOM-W/AMSR2 SST in obs space --> model space
!-------------------------------------------------------------

subroutine project_gsst(iswitch,&
     & im_g,jm_g,lon_g,lat_g,sst_g,&
     & im,jm,lon,lat,gsst,gcount)

  use mod_rmiss
  implicit none

  !Common
  integer ig,jg
  integer i,j

  real(kind = 8) dx,dy

  !IN
  integer,intent(in) :: iswitch
  integer,intent(in) :: im_g,jm_g
  integer,intent(in) :: im,jm

  real(kind = 8),intent(in) :: lon_g(im_g,jm_g),lat_g(im_g,jm_g),sst_g(im_g,jm_g)
  real(kind = 8),intent(in) :: lon(im),lat(jm)

  !INOUT
  real(kind = 8),intent(inout) ::  gsst(im,jm),gcount(im,jm)

  dx=lon(2)-lon(1)
  dy=lat(2)-lat(1)

  if(iswitch == 1)then
     do jg=1,jm_g
        do ig=1,im_g
           
           if(sst_g(ig,jg) == rmiss) cycle
           if(lon_g(ig,jg) < lon(1)-0.5d0*dx &
                & .or. lon(im)+0.5d0*dx < lon_g(ig,jg) &
                & .or. lat_g(ig,jg) < lat(1)-0.5d0*dy &
                & .or. lat(jm)+0.5d0*dy < lat_g(ig,jg)) cycle

           do j=1,jm
              if(lat_g(ig,jg) < lat(j)-0.5d0*dy &
                   & .or. lat(j)+0.5d0*dy < lat_g(ig,jg))cycle
              do i=1,im
                 if(lon_g(ig,jg) < lon(i)-0.5d0*dx &
                      & .or. lon(i)+0.5d0*dx < lon_g(ig,jg))cycle
                 gsst(i,j)=gsst(i,j)+sst_g(ig,jg)
                 gcount(i,j)=gcount(i,j)+1.d0
              end do !i
           end do !j

        end do !ig
     end do !jg
  end if
  
  if(iswitch == 2)then
     do j=1,jm
        do i=1,im
           if(gcount(i,j) == 0.d0)then
              gsst(i,j)=rmiss
           else
              gsst(i,j)=gsst(i,j)/gcount(i,j)
           end if
        end do
     end do
  end if

end subroutine project_gsst

!-------------------------------------------------------------------------
! Merge Himawari-8 and GCOM-W/AMSR2 SST |
!-------------------------------------------------------------------------

subroutine merge_sst(im,jm,asst,csst,gsst,hsst,sst)

  use mod_rmiss
  use setting, only: nsst
  implicit none
  
  !Common
  integer i,j
  integer isst
  
  real(kind = 8) tmp(im,jm,nsst)
  real(kind = 8) pass(im,jm),miss(im,jm)
  
  !IN
  integer,intent(in) :: im,jm
  real(kind = 8),intent(in) :: asst(im,jm),csst(im,jm)
  real(kind = 8),intent(in) :: gsst(im,jm),hsst(im,jm)

  !OUT
  real(kind = 8),intent(out) ::  sst(im,jm)

  tmp(:,:,1)=asst(:,:)
  tmp(:,:,2)=csst(:,:)
  tmp(:,:,3)=gsst(:,:)
  tmp(:,:,4)=hsst(:,:)
  sst(:,:)=rmiss
  pass(:,:)=0.d0
  miss(:,:)=0.d0

  !Add SST
  do isst=1,nsst
     do j=1,jm
        do i=1,im

           if(tmp(i,j,isst) == rmiss)then
              miss(i,j)=miss(i,j)+1.d0
           else
              sst(i,j)=sst(i,j)+tmp(i,j,isst)
              pass(i,j)=pass(i,j)+1.d0
           end if
                
        end do
     end do
  end do

  !Composite SST
  do j=1,jm
     do i=1,im
        if(pass(i,j) == 0.d0)then
           sst(i,j)=rmiss
        else
           sst(i,j)=sst(i,j)/pass(i,j)
        end if
     end do
  end do
     
end subroutine merge_sst

!----------------------------------------------------------------------
! Write SST |
!----------------------------------------------------------------------

subroutine write_sst(iyr,imon,iday,lon,lat,asst,csst,gsst,hsst,sst)

  use netcdf
  use mod_gridinfo
  implicit none

  !Common
  integer status,access,system
  integer ncid,varid

  character(2) mm
  character(4) yyyy
  character(100) filename
  
  !IN
  integer,intent(in) :: iyr,imon,iday
  
  real(kind = 8),intent(in) :: lon(im),lat(jm)
  real(kind = 8),intent(in) :: asst(im,jm),csst(im,jm),gsst(im,jm),hsst(im,jm)
  real(kind = 8),intent(in) :: sst(im,jm)

  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon

  filename="../obs/sst"//yyyy//mm//".nc"
  
  status=access(trim(filename)," ")
  if(status /= 0)then
     call make_ncfile_sst(im,jm,iyr,imon,filename)
  end if

  status=nf90_open(trim(filename),nf90_write,ncid)

  status=nf90_inq_varid(ncid,"time",varid)
  status=nf90_put_var(ncid,varid,real(iday-0.5e0))
  
  status=nf90_inq_varid(ncid,"lon",varid)
  status=nf90_put_var(ncid,varid,real(lon))

  status=nf90_inq_varid(ncid,"lat",varid)
  status=nf90_put_var(ncid,varid,real(lat))

  status=nf90_inq_varid(ncid,"asst",varid)
  status=nf90_put_var(ncid,varid,real(asst(1:im,1:jm)),(/1,1,iday/),(/im,jm,1/))

  status=nf90_inq_varid(ncid,"csst",varid)
  status=nf90_put_var(ncid,varid,real(csst(1:im,1:jm)),(/1,1,iday/),(/im,jm,1/))

  status=nf90_inq_varid(ncid,"gsst",varid)
  status=nf90_put_var(ncid,varid,real(gsst(1:im,1:jm)),(/1,1,iday/),(/im,jm,1/))
  
  status=nf90_inq_varid(ncid,"hsst",varid)
  status=nf90_put_var(ncid,varid,real(hsst(1:im,1:jm)),(/1,1,iday/),(/im,jm,1/))
  
  status=nf90_inq_varid(ncid,"sst",varid)
  status=nf90_put_var(ncid,varid,real(sst(1:im,1:jm)),(/1,1,iday/),(/im,jm,1/))

  status=nf90_close(ncid)
  
end subroutine write_sst
