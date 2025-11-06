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
!          by S.Ohishi 2025.07
!------------------------------------------------------
subroutine prepare_sst(iyr,imon,iday,inum,lon,lat,fsm)
  
  use setting, only: idx_msst, idy_msst, idx_hsst, idy_hsst, iswitch_msst, iswitch_hsst
  use mod_rmiss
  use mod_gridinfo
  use mod_read_himawari, im_h => im, jm_h => jm
  use mod_read_amsre
  use mod_read_windsat
  use mod_read_amsr2
  implicit none

  !---Parameter
  integer,parameter :: nsat=3

  !---Common
  integer isat
  integer ihour
  integer jyr,jmon,jday
  integer i,j
  integer ins
  
  !---GCOM-W/AMSR2, Aqua/AMSRE, Coriolis/WindSat
  integer ifile,nfile  
  integer im_g,jm_g
  real(kind = 8),allocatable :: lon_g(:,:),lat_g(:,:),sst_g(:,:)
  character(200),allocatable :: filename(:)

  !---Himawari/AHI
  real(kind = 8) lon_h(im_h),lat_h(jm_h)
  real(kind = 8) dat_h(im_h,jm_h)
  real(kind = 8) sst_h(im_h,jm_h),pass_h(im_h,jm_h),miss_h(im_h,jm_h)
  
  !---IN
  integer,intent(in) :: iyr,imon,iday
  real(kind = 8),intent(in) :: lon(im),lat(jm)
  real(kind = 8),intent(in) :: fsm(im,jm)

  !---IN/OUT
  integer,intent(inout) :: inum  

  !---Microwave Satellite (AMSR-E, WindSat, and AMSR2)
  if(iswitch_msst == 1)then

     do isat=1,nsat
        
        !---Read filename
        if(isat == 1)then
           ins=1
           call read_amsre_filename(iyr,imon,iday,nfile,filename)
           write(*,*) "Aqua/AMSR-E SST:",nfile
        else if(isat == 2)then
           ins=2
           call read_windsat_filename(iyr,imon,iday,nfile,filename)
           write(*,*) "Coriolis/WindSAT SST:",nfile
        else if(isat == 3)then
           ins=3
           call read_amsr2_filename(iyr,imon,iday,nfile,filename)
           write(*,*) "GCOM-W/AMSR-2 SST:",nfile
        end if
        if(nfile == 0) cycle

        do ifile=1,nfile

           !---Read data
           if(isat == 1)then
              call read_amsre(filename(ifile),jyr,jmon,jday,ihour,im_g,jm_g,lon_g,lat_g,sst_g)
           else if(isat == 2)then
              call read_windsat(filename(ifile),jyr,jmon,jday,ihour,im_g,jm_g,lon_g,lat_g,sst_g)
           else if(isat == 3)then
              call read_amsr2(filename(ifile),jyr,jmon,jday,ihour,im_g,jm_g,lon_g,lat_g,sst_g)
           end if

           if(iyr /= jyr .or. imon /= jmon .or. iday /= jday)then
              write(*,*) "***Error: Date is inconsistent"
              write(*,*) iyr,imon,iday
              write(*,*) jyr,jmon,jday
              stop
           end if
           
           !---Write data
           call write_obs_surface2d("sst",ins,iyr,imon,iday, &
                & idx_msst,idy_msst,im_g,jm_g,lon(1),lon(im),lon_g,lat(1),lat(jm),lat_g,sst_g,inum)

           !---Deallocate
           if(isat == 1)then
              call deallocate_amsre(lon_g,lat_g,sst_g)
           else if(isat == 2)then
              call deallocate_windsat(lon_g,lat_g,sst_g)
           else if(isat == 3)then
              call deallocate_amsr2(lon_g,lat_g,sst_g)
           end if
                
        end do !ifile

        !---Deallocate(filename)
        if(isat == 1)then
           call deallocate_amsre_filename(filename)
        else if(isat == 2)then
           call deallocate_windsat_filename(filename)
        else if(isat == 3)then
           call deallocate_amsr2_filename(filename)
        end if
        
     end do !isat

  end if

  !---Himawari-8/-9
  if(iswitch_hsst == 1)then
     ins=4
     write(*,'(a)') "Daily-mean Himawari SST"
     do ihour=0,23
        call read_himawari(iyr,imon,iday,ihour,lon_h,lat_h,dat_h)
        call cal_daily_mean(ihour,im_h,jm_h,dat_h,sst_h,pass_h,miss_h)
     end do
     call write_obs_surface2dg("sst",ins,iyr,imon,iday,im,jm,lon,lat,&
          & idx_hsst,idy_hsst,im_h,jm_h,lon(1),lon(im),lon_h,lat(1),lat(jm),lat_h,sst_h,inum)
  end if

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
