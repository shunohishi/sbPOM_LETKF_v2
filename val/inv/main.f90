module setting

  integer,parameter :: syr=2021,smon=1,sday=1
  integer,parameter :: eyr=2021,emon=1,eday=1
  
  character(10),parameter :: dir="TEST",letkf="letkf"

end module setting

!------------------------------------------------------------------------

program main

  use setting
  use mod_julian
  use mod_read_lora
  implicit none
  
  !---Common
  integer ijul,sjul,ejul
  integer iyr,imon,iday
  
  integer iobs,imem

  !---Data in obs. space
  integer nobs,nmem
  integer,allocatable :: ele(:),ins(:)

  real,allocatable :: lon(:),lat(:),lev(:)
  real,allocatable :: obs(:),err(:)
  real,allocatable :: hxfmean(:),hxfsprd(:)
  real,allocatable :: hxamean(:),hxasprd(:)
  real,allocatable :: hxf(:,:),hxa(:,:)
  
  !--Julian day
  call ymd_julian(syr,smon,sday,sjul)
  call ymd_julian(eyr,emon,eday,ejul)
  
  do ijul=sjul,ejul

     !---Date
     call julian_ymd(ijul,iyr,imon,iday)
     write(*,*) iyr,imon,iday

     !---Read data
     call read_inv(dir,letkf,iyr,imon,iday,nobs,nmem, &
          & ele,ins, &
          & lon,lat,lev,obs,err,hxfmean,hxfsprd,hxamean,hxasprd,hxf,hxa)

     !---Statistics
     !Bias & RMSD

     !Histogram
     
     !---End read
     call end_read_inv(ele,ins, &
          & lon,lat,lev,obs,err,hxfmean,hxfsprd,hxamean,hxasprd,hxf,hxa)     

  end do !ijul

end program main
