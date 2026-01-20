module setting

  integer,parameter :: ndata=4 !1: LORA, 2: GLORYS, 3: ORAS5, 4:GLORS
  integer,parameter :: syr=2003
  integer,parameter :: eyr=2023
  
end module setting

!------------------------------------------------------------------------------------------

program main

  use setting
  use mod_gridinfo, im_lora => im, jm_lora => jm, km_lora => km
  use mod_read_glorys025, im_g025 => im, jm_g025 => jm, km_g025 => km
  implicit none

  !---Common

  integer idata,iyr,imon

  !DA
  integer im,jm,km

  real(kind = 8),allocatable :: lon(:),lat(:),depth(:,:,:)
  real(kind = 8),allocatable :: mask(:,:)
  real(kind = 8),allocatable :: mclim(:,:,:)

  !Index
  real(kind = 8) nino(syr:eyr,12,ndata),nino_norm(syr:eyr,12,ndata)
  real(kind = 8) iod(syr:eyr,12,ndata),iod_norm(syr:eyr,12,ndata)
  real(kind = 8) anino(syr:eyr,12,ndata),anino_norm(syr:eyr,12,ndata)

  !---Initialization (Just in case)
  nino(:,:,:)=0.d0
  iod(:,:,:)=0.d0
  anino(:,:,:)=0.d0
  
  do idata=1,ndata

     if(idata == 1)then
        im=im_lora
        jm=jm_lora
        km=1
     else if(idata == 2 .or. idata == 3 .or. idata == 4)then
        im=im_g025
        jm=jm_g025
        km=1
     else
        write(*,*) "***Error: Incorrect idata"
        stop
     end if        

     allocate(lon(im),lat(jm),depth(im,jm,km))
     allocate(mask(im,jm))
     allocate(mclim(im,jm,km))

     !---Read grid information
     write(*,*) "Read grid"
     call read_grid(idata,im,jm,km,lon,lat,depth,mask)
     
     do imon=1,12
        
        !---Monthly climatology
        write(*,*) "Monthly Climatology"
        call estimate_mclim(idata,imon,im,jm,km,mask,mclim)

        !---Write Climatology
        call write_mclim(idata,imon,im,jm,km,lon,lat,depth,mclim)
        
        !---Monthly climate index
        write(*,*) "Monthly climate index"
        do iyr=syr,eyr
           call estimate_cindex(idata,iyr,imon,im,jm,km,lon,lat,depth,mask,mclim, &
                & nino(iyr,imon,idata),iod(iyr,imon,idata),anino(iyr,imon,idata))
        end do !iyr
        
        !---Monthly anomaly at event year
        call normalization(nino(:,imon,idata),nino_norm(:,imon,idata))
        call normalization(iod(:,imon,idata),iod_norm(:,imon,idata))
        call normalization(anino(:,imon,idata),anino_norm(:,imon,idata))
                
     end do !imon
     
     deallocate(lon,lat,depth)
     deallocate(mask)
     deallocate(mclim)
     
  end do !idata

  call write_index(nino,nino_norm,iod,iod_norm,anino,anino_norm)
  
end program main

!-----------------------------------------------------------------------

subroutine read_grid(idata,im,jm,km,lon,lat,depth,mask)

  use setting
  use mod_gridinfo, only: km_lora => km
  use mod_read_lora, only: read_grid_lora => read_grid
  use mod_read_glorys025, only: read_glorys025
  implicit none

  !---Common
  integer k

  character(10) datname
  character(1) varname
  
  !---IN
  integer,intent(in) :: idata
  integer,intent(in) :: im,jm,km

  real(kind = 8) tmp1dx(im),tmp1dy(jm),tmp1dz(km),tmp2d(im,jm),tmp3d(im,jm,km)

  !---LORA
  real(kind = 8) depth_lora(im,jm,km_lora)
  
  real(kind = 8) tmp3d_lora(im,jm,km_lora)
  
  !---OUT
  real(kind = 8),intent(out) :: lon(im),lat(jm),depth(im,jm,km)
  real(kind = 8),intent(out) :: mask(im,jm)
  
  if(idata == 1)then
     call read_grid_lora("QGLOBAL",lon,tmp1dx,tmp1dx,lat,tmp1dy,tmp1dy,depth_lora,tmp3d_lora,tmp3d_lora, &
          & mask,tmp2d,tmp2d)
     do k=1,km
        depth(1:im,1:jm,k)=depth_lora(1:im,1:jm,k)
     end do
  else if(idata == 2)then
     datname="glorys"
     varname="t"
     call read_glorys025(datname,varname,syr,1,1,km,lon,lat,tmp1dz,mask,tmp3d)
     do k=1,km
        depth(1:im,1:jm,k)=tmp1dz(k)
     end do
  else if(idata == 3)then
     datname="oras5"
     varname="t"
     call read_glorys025(datname,varname,syr,1,1,km,lon,lat,tmp1dz,mask,tmp3d)
     do k=1,km
        depth(1:im,1:jm,k)=tmp1dz(k)
     end do
  else if(idata == 4)then
     datname="cglors"
     varname="t"
     call read_glorys025(datname,varname,syr,1,1,km,lon,lat,tmp1dz,mask,tmp3d)
     do k=1,km
        depth(1:im,1:jm,k)=tmp1dz(k)
     end do
  end if
       
end subroutine read_grid

!-----------------------------------------------------------------------

subroutine read_data(idata,iyr,imon,iday,im,jm,km,mask,dat)

  use mod_read_lora
  use mod_read_glorys025, only: read_glorys025
  implicit none

  !---Common
  integer itmp
  
  real(kind = 8) tmp1dx(im),tmp1dy(jm),tmp1dz(km)
  real(kind = 8) tmp2d(im,jm)

  character(10) datname
  character(1) varname
  
  !---IN
  integer,intent(in) :: idata
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: mask(im,jm)
  
  !---OUT
  real(kind = 8),intent(out) :: dat(im,jm,km)
  
  if(idata == 1)then !LORA
     itmp=1 !Dummy
     call read_anal("QGLOBAL","letkf","qglobal","mean",itmp,"t",iyr,imon,iday,im,jm,km,mask,dat)
  else if(idata == 2)then !GLORYS
     datname="glorys"
     varname="t"
     call read_glorys025(datname,varname,iyr,imon,iday,km,tmp1dx,tmp1dy,tmp1dz,tmp2d,dat)
  else if(idata == 3)then !ORAS5
     datname="oras5"
     varname="t"
     call read_glorys025(datname,varname,iyr,imon,iday,km,tmp1dx,tmp1dy,tmp1dz,tmp2d,dat)
  else if(idata == 4)then !C-GLORS
     datname="cglors"
     varname="t"
     call read_glorys025(datname,varname,iyr,imon,iday,km,tmp1dx,tmp1dy,tmp1dz,tmp2d,dat)
  end if
  
end subroutine read_data

!-----------------------------------------------------------------------

subroutine add_data(im,jm,km,dat,mask,sum,pass,miss)

  implicit none

  !---Common
  integer i,j,k
  
  !---IN
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: dat(im,jm,km),mask(im,jm)

  !---IN/OUT
  real(kind = 8),intent(inout) :: sum(im,jm,km),pass(im,jm,km),miss(im,jm,km)

  do k=1,km
     do j=1,jm
        do i=1,im
           if(mask(i,j) == 1.d0)then
              sum(i,j,k)=sum(i,j,k)+dat(i,j,k)
              pass(i,j,k)=pass(i,j,k)+1.d0              
           else if(mask(i,j) == 0.d0)then
              miss(i,j,k)=miss(i,j,k)+1.d0
           end if           
        end do
     end do
  end do
  
end subroutine add_data

!-----------------------------------------------------------------------

subroutine ave_data(im,jm,km,mean,pass,miss)

  use mod_rmiss
  implicit none

  !---Common
  integer i,j,k

  !---IN
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: pass(im,jm,km),miss(im,jm,km)

  !---IN/OUT
  real(kind = 8),intent(inout) :: mean(im,jm,km)

  do k=1,km
     do j=1,jm
        do i=1,im
           if(pass(i,j,k) == 0.d0)then
              mean(i,j,k)=rmiss
           else
              mean(i,j,k)=mean(i,j,k)/pass(i,j,k)
           end if           
        end do
     end do
  end do
  
end subroutine ave_data

!-----------------------------------------------------------------------

subroutine estimate_anomaly(im,jm,km,mean,clim,anom)

  use mod_rmiss
  implicit none

  !---Common
  integer i,j,k
  
  !---IN
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: mean(im,jm,km),clim(im,jm,km)
  
  !---OUT
  real(kind = 8),intent(out) :: anom(im,jm,km)

  do k=1,km
     do j=1,jm
        do i=1,im
           if(mean(i,j,k) == rmiss .or. clim(i,j,k) == rmiss)then
              anom(i,j,k)=rmiss
           else
              anom(i,j,k)=mean(i,j,k)-clim(i,j,k)
           end if           
        end do
     end do
  end do
  
end subroutine estimate_anomaly
  
!-----------------------------------------------------------------------

subroutine estimate_mclim(idata,imon,im,jm,km,mask,mclim)

  use setting, only: syr,eyr
  use mod_rmiss
  implicit none

  !---Common
  integer iyr,iday,nday

  real(kind = 8) dat(im,jm,km)
  real(kind = 8) mpass(im,jm,km),mmiss(im,jm,km)
  
  !---IN
  integer,intent(in) :: idata
  integer,intent(in) :: imon
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: mask(im,jm)
  
  !---OUT
  real(kind = 8),intent(out) :: mclim(im,jm,km)

  mclim(:,:,:)=0.d0
  mpass(:,:,:)=0.d0
  mmiss(:,:,:)=0.d0
  
  do iyr=syr,eyr

     call number_of_day(iyr,imon,nday)
     
     do iday=1,nday

        write(*,*) "clim:",idata,iyr,imon,iday

        !Read data
        call read_data(idata,iyr,imon,iday,im,jm,km,mask,dat)
        
        !Add data
        call add_data(im,jm,km,dat,mask,mclim,mpass,mmiss)
        
     end do
  end do

  !Average data
  call ave_data(im,jm,km,mclim,mpass,mmiss)
          
end subroutine estimate_mclim

!------------------------------------------------------------------------

subroutine estimate_cindex(idata,iyr,imon,im,jm,km,lon,lat,depth,mask,mclim,nino,iod,anino)

  implicit none

  !---Common
  integer iday,nday

  real(kind = 8) dat(im,jm,km)
  real(kind = 8) mmean(im,jm,km),manom(im,jm,km)
  real(kind = 8) mpass(im,jm,km),mmiss(im,jm,km)
  
  !---IN
  integer,intent(in) :: idata
  integer,intent(in) :: iyr,imon
  integer,intent(in) :: im,jm,km
  
  real(kind = 8),intent(in) :: lon(im),lat(jm),depth(im,jm,km)  
  real(kind = 8),intent(in) :: mask(im,jm)
  real(kind = 8),intent(in) :: mclim(im,jm,km)

  !---OUT
  real(kind = 8) nino,iod,anino

  call number_of_day(iyr,imon,nday)

  mmean(:,:,:)=0.d0
  mpass(:,:,:)=0.d0
  mmiss(:,:,:)=0.d0

  do iday=1,nday

     write(*,*) "mmean:",idata,iyr,imon,iday
     
     !Read data
     call read_data(idata,iyr,imon,iday,im,jm,km,mask,dat)
     
     !Add data
     call add_data(im,jm,km,dat,mask,mmean,mpass,mmiss)
     
  end do
  
  !Average data
  call ave_data(im,jm,km,mmean,mpass,mmiss)
  
  !Anomaly
  call estimate_anomaly(im,jm,km,mmean,mclim,manom)

  !Write monthly mean
  call write_month(idata,iyr,imon,im,jm,km,lon,lat,depth,mmean,manom)

  !INDEX
  call nino34_index(im,jm,lon,lat,manom(:,:,1),nino)

  call iod_index(im,jm,lon,lat,manom(:,:,1),iod)

  call atlantic_nino_index(im,jm,lon,lat,manom(:,:,1),anino)
  
end subroutine estimate_cindex

