program main

  !***To be modified ==> mod_gridinfo, mod_read_glorys
  use mod_rmiss
  use setting, ndat_a => ndat
  use mod_gridinfo, im_lora => im, jm_lora => jm, km_lora => km
  use mod_read_glorys025, im_g025 => im, jm_g025 => jm, km_g025 => km
  use mod_read_duacs, im_dua => im, jm_dua => jm
  use mod_read_woa18_month, im_woa => im, jm_woa => jm, kmm_woa => km
  use mod_read_woa18_season_annual,only: kms_woa => km, read_woa18_season_annual_var 
  use mod_io  
  implicit none

  !---Parameter
  integer,parameter :: nvar3d=2 !T,S
  integer,parameter :: imon=0   !Annual
  !---Common
  integer ivar
  integer,allocatable :: idlon(:),idlat(:)
  
  real(kind = 8) rlev !Reference level (SSH)
  
  character(1) varname
  character(100) filename
  
  !---Analysis
  integer idat_a
  integer im_a,jm_a,km_a

  real(kind = 8),allocatable :: lon_a(:),lont_a(:),lonu_a(:),lonv_a(:)
  real(kind = 8),allocatable :: lat_a(:),latt_a(:),latu_a(:),latv_a(:)
  real(kind = 8),allocatable :: dep_a(:,:,:),dept_a(:,:,:),depu_a(:,:,:),depv_a(:,:,:)
  real(kind = 8),allocatable :: mask_a(:,:),maskt_a(:,:),masku_a(:,:),maskv_a(:,:)

  real(kind = 8),allocatable :: mean2d_a(:,:),sprd2d_a(:,:)
  real(kind = 8),allocatable :: mean3d_a(:,:,:),sprd3d_a(:,:,:)

  real(kind = 8) mdot_ha(im_dua,jm_dua,ndat_a)

  real(kind = 8),allocatable :: mean3d_ha(:,:,:),sprd3d_ha(:,:,:)
  real(kind = 8),allocatable :: mean2d_ha_yz(:,:)
  
  !---AVISO
  real(kind = 8) lon_dua(im_dua),lat_dua(jm_dua)
  real(kind = 8) mdot_dua(im_dua,jm_dua)
  
  !---WOA18
  integer km_woa
  
  real(kind = 8) lon_woa(im_woa),lat_woa(jm_woa)
  real(kind = 8),allocatable :: dep_woa(:)

  real(kind = 8),allocatable :: mean3d_woa(:,:,:),sprd3d_woa(:,:,:)
  real(kind = 8),allocatable :: mean2d_woa_yz(:,:)
  
  !=====SSH================================================================================================
  varname="h"     
  
  !---Read Observational MDOT
  call read_duacs_mdot(lon_dua,lat_dua,mdot_dua)

  !---Allocate
  allocate(idlon(im_dua),idlat(jm_dua))
  
  do idat_a=1,ndat_a

     !*** To be modified
     !---Read Analysis grid information
     if(idat_a == 1)then
        im_a=im_lora
        jm_a=jm_lora
        km_a=km_lora        
     else
        im_a=im_g025
        jm_a=jm_g025
        km_a=km_g025        
     end if

     !---Allocate     
     allocate(lon_a(im_a),lont_a(im_a),lonu_a(im_a),lonv_a(im_a))
     allocate(lat_a(jm_a),latt_a(jm_a),latu_a(jm_a),latv_a(jm_a))
     allocate(dep_a(im_a,jm_a,km_a),dept_a(im_a,jm_a,km_a),depu_a(im_a,jm_a,km_a),depv_a(im_a,jm_a,km_a))
     allocate(mask_a(im_a,jm_a),maskt_a(im_a,jm_a),masku_a(im_a,jm_a),maskv_a(im_a,jm_a))
     allocate(mean2d_a(im_a,jm_a),sprd2d_a(im_a,jm_a))
     
     !---Read analysis grid data
     write(*,*) "Read grid data"
     call read_grid(idat_a,im_a,jm_a,km_a, &
          & lont_a,lonu_a,lonv_a,latt_a,latu_a,latv_a,dept_a,depu_a,depv_a,maskt_a,masku_a,maskv_a)

     lon_a(:)=lont_a(:)
     lat_a(:)=latt_a(:)
     mask_a(:,:)=maskt_a(:,:)
     
     !---Read analysis SSH data
     filename="dat/"//trim(datname(idat_a))//"_clim.nc"
     call read_clim_nc(filename,varname,im_a,jm_a,1,imon,mean2d_a,sprd2d_a)

     !---ID
     call cal_idlon(im_a,lon_a,im_dua,lon_dua,idlon)
     call cal_idlat(jm_a,lat_a,jm_dua,lat_dua,idlat)
     
     !---Bilinear interpolation
     call bilinear_interpolation_2d &
          & (im_a,jm_a,lon_a,lat_a,mean2d_a,mask_a, &
          &  im_dua,jm_dua,lon_dua,lat_dua, &
          &  idlon,idlat,mdot_ha(:,:,idat_a))

     !---Deallocate
     deallocate(lon_a,lont_a,lonu_a,lonv_a)
     deallocate(lat_a,latt_a,latu_a,latv_a)
     deallocate(dep_a,dept_a,depu_a,depv_a)
     deallocate(mask_a,maskt_a,masku_a,maskv_a)
     deallocate(mean2d_a,sprd2d_a)
     
  end do !idat
        
  !---Remove reference level
  call remove_reference_level(ndat_a,im_dua,jm_dua,lat_dua,mdot_dua,mdot_ha)
              
  !---Write analysis SSH data
  call write_mdot(ndat_a,im_dua,jm_dua,lon_dua,lat_dua,mdot_dua,mdot_ha)

  !---Deallocate
  deallocate(idlon,idlat)
  
  !=====3D=============================================================================================
  !---WOA18 grid information        
  km_woa=kms_woa

  allocate(dep_woa(km_woa))
  allocate(mean3d_woa(im_woa,jm_woa,km_woa),sprd3d_woa(im_woa,jm_woa,km_woa))
  allocate(idlon(im_woa),idlat(jm_woa))
  allocate(mean3d_ha(im_woa,jm_woa,km_woa),sprd3d_ha(im_woa,jm_woa,km_woa))
  allocate(mean2d_ha_yz(jm_woa,km_woa),mean2d_woa_yz(jm_woa,km_woa))
  
  do idat_a=1,ndat_a

     !*** To be modified
     !---Read Analysis grid information
     if(idat_a == 1)then
        im_a=im_lora
        jm_a=jm_lora
        km_a=km_lora        
     else
        im_a=im_g025
        jm_a=jm_g025
        km_a=km_g025        
     end if
     
     !---Allocate     
     allocate(lon_a(im_a),lont_a(im_a),lonu_a(im_a),lonv_a(im_a))
     allocate(lat_a(jm_a),latt_a(jm_a),latu_a(jm_a),latv_a(jm_a))
     allocate(dep_a(im_a,jm_a,km_a),dept_a(im_a,jm_a,km_a),depu_a(im_a,jm_a,km_a),depv_a(im_a,jm_a,km_a))
     allocate(mask_a(im_a,jm_a),maskt_a(im_a,jm_a),masku_a(im_a,jm_a),maskv_a(im_a,jm_a))
     allocate(mean3d_a(im_a,jm_a,km_a),sprd3d_a(im_a,jm_a,km_a))

     !---Read analysis grid data
     write(*,*) "Read grid data"
     call read_grid(idat_a,im_a,jm_a,km_a, &
          & lont_a,lonu_a,lonv_a,latt_a,latu_a,latv_a,dept_a,depu_a,depv_a,maskt_a,masku_a,maskv_a)

     lon_a(:)=lont_a(:)
     lat_a(:)=latt_a(:)
     dep_a(:,:,:)=dept_a(:,:,:)
     mask_a(:,:)=maskt_a(:,:)
             
     !---Get ID        
     do ivar=1,2

        if(ivar == 1)then
           varname="t"
        else if(ivar == 2)then
           varname="s"
        end if

        write(*,*) "Data: "//trim(datname(idat_a))//" Var: "//trim(varname)
        
        !---Read WOA
        call read_woa18_season_annual_var(varname,imon,lon_woa,lat_woa,dep_woa,mean3d_woa)
        sprd3d_woa(:,:,:)=rmiss
        
        !---Read analysis data
        filename="dat/"//trim(datname(idat_a))//"_clim.nc"
        call read_clim_nc(filename,varname,im_a,jm_a,km_a,imon,mean3d_a,sprd3d_a)

        !---ID
        call cal_idlon(im_a,lon_a,im_woa,lon_woa,idlon)
        call cal_idlat(jm_a,lat_a,jm_woa,lat_woa,idlat)
        
        !---Interpolate analysis data
        call vertical_bilinear_interpolate(im_a,jm_a,km_a,lon_a,lat_a,dep_a,mean3d_a, &
             & im_woa,jm_woa,km_woa,lon_woa,lat_woa,dep_woa,idlon,idlat,mean3d_ha)

        !---Zonal average => Meridional section
        call zonal_average(im_woa,jm_woa,km_woa,mean3d_ha,mean2d_ha_yz)
        call zonal_average(im_woa,jm_woa,km_woa,mean3d_woa,mean2d_woa_yz)
        
        !---Write data
        call write_var3d_txt(datname(idat_a),varname,im_woa,jm_woa,km_woa,lon_woa,lat_woa,dep_woa,mean3d_ha,mean3d_woa)
        call write_var2d_txt(datname(idat_a),varname,jm_woa,km_woa,lat_woa,dep_woa,mean2d_ha_yz,mean2d_woa_yz)
        
     end do !ivar
     
     !---Deallocate
     deallocate(lon_a,lont_a,lonu_a,lonv_a)
     deallocate(lat_a,latt_a,latu_a,latv_a)
     deallocate(dep_a,dept_a,depu_a,depv_a)
     deallocate(mask_a,maskt_a,masku_a,maskv_a)
     deallocate(mean3d_a,sprd3d_a)
     
  end do !idat_a
  
  !---Deallocate
  deallocate(dep_woa)
  deallocate(mean3d_woa,sprd3d_woa)
  deallocate(idlon,idlat)
  deallocate(mean3d_ha,sprd3d_ha)
  deallocate(mean2d_ha_yz,mean2d_woa_yz)
  
end program main

!--------------------------------------------------------------------

subroutine remove_reference_level(ndat,im,jm,lat,obs,hdat)

  use mod_rmiss
  implicit none

  !---Parameter
  real(kind = 8),parameter :: pi=4.d0*atan(1.d0)
  
  !---Common
  integer i,j,idat

  real(kind = 8) weight,sum_weight
  real(kind = 8) rlev(0:ndat)
  real(kind = 8) tmp(im,jm,0:ndat)
  
  !---IN
  integer,intent(in) :: ndat
  integer,intent(in) :: im,jm

  real(kind = 8),intent(in) :: lat(jm)

  !---IN/OUT
  real(kind = 8),intent(inout) :: obs(im,jm),hdat(im,jm,ndat)
  
  tmp(:,:,0)=obs(:,:)
  tmp(:,:,1:ndat)=hdat(:,:,1:ndat)

  !---Fill in missing value
  do idat=0,ndat
     do j=1,jm
        do i=1,im

           if(tmp(i,j,idat) == rmiss)then
              tmp(i,j,:)=rmiss
           end if
           
        end do
     end do
  end do
        
  !---Reference level for each dataset
  do idat=0,ndat

     rlev(idat)=0.d0
     sum_weight=0.d0
     
     do j=1,jm

        weight=cos(pi*lat(j)/180.d0)

        do i=1,im

           if(tmp(i,j,idat) == rmiss) cycle
        
           rlev(idat)=rlev(idat)+tmp(i,j,idat)*weight
           sum_weight=sum_weight+weight
        
        end do
     end do

     if(sum_weight == 0.d0)then
        rlev(idat)=rmiss
     else
        rlev(idat)=rlev(idat)/sum_weight
     end if

     write(*,'(i6,a,f12.5)') idat," Reference level:",rlev(idat)
     
  end do

  !---Remove reference level
  do idat=0,ndat
     do j=1,jm
        do i=1,im
           if(tmp(i,j,idat) == rmiss) cycle
           tmp(i,j,idat)=tmp(i,j,idat)-rlev(idat)
        end do
     end do
  end do

  !---Substitute
  obs(:,:)=tmp(:,:,0)
  hdat(:,:,1:ndat)=tmp(:,:,1:ndat)
  
end subroutine remove_reference_level

!--------------------------------------------------------------------------------

subroutine zonal_average(im,jm,km,dat3d,dat2d)

  use mod_rmiss
  implicit none

  !---Common
  integer i,j,k
  integer ipass,imiss
  
  !---IN
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: dat3d(im,jm,km)

  !---OUT
  real(kind = 8),intent(out) :: dat2d(jm,km)

  do k=1,km
     do j=1,jm

        dat2d(j,k)=0.d0
        ipass=0
        imiss=0
        
        do i=1,im
           if(dat3d(i,j,k) == rmiss)then
              imiss=imiss+1
           else
              dat2d(j,k)=dat2d(j,k)+dat3d(i,j,k)
              ipass=ipass+1
           end if
        end do

        if(ipass == 0)then
           dat2d(j,k)=rmiss
        else
           dat2d(j,k)=dat2d(j,k)/dble(ipass)
        end if

     end do
  end do
          
end subroutine zonal_average
