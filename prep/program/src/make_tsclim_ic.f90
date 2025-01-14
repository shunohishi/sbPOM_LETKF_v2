module setting

  real(kind = 8),parameter :: rhoref=1025.d0,grav=9.806d0 !Consistent with POM code
  
end module setting

program main

  !------------------------------------------------------------------------
  ! Make IC & TSCLIM |
  !------------------------------------------------------------------------
  ! 
  ! WOA: z coordinate --> Model: sigma coordinate
  !
  !------------------------------------------------
  !
  ! 1. Read Data
  ! 2. Fill value in horizontal direction at some extent
  ! 3. Bilinear interpolation in horizontal direction
  ! 4. Linear interpolation in vertical direction
  ! 5. Fill value in vertical direction
  ! 6. Write Data
  !
  ! Created by 
  ! 2018.08 S.Ohishi
  ! 2023.03 S.Ohishi
  !
  !------------------------------------------------------------------------

  use mod_rmiss
  use mod_gridinfo
  use mod_read_woa18_month, im_woa => im, jm_woa => jm, kmm_woa => km
  use mod_read_woa18_season_annual,only: kms_woa => km, read_woa18_season_annual_var
  implicit none

  !Common
  integer imon,ivar
  integer k
  character(1) varname
  
  !Model
  integer idx(im),idy(jm),idz_s(im,jm,km),idz_m(im,jm,km)
  integer iqglobal

  real(kind = 8) lon(im),lat(jm),depth(im,jm,km)
  real(kind = 8) fsm(im,jm)
  real(kind = 8) dat(im,jm,km)
  real(kind = 8),allocatable :: null1x(:),null1y(:),null2d(:,:)
  real(kind = 8),allocatable :: null3d(:,:,:) !NULL

  !WOA18
  real(kind = 8) lon_woa(im_woa),lat_woa(jm_woa)
  real(kind = 8),allocatable :: depth_woa(:)
  real(kind = 8),allocatable :: t_woa(:,:,:),s_woa(:,:,:)
  real(kind = 8),allocatable :: dat_woa(:,:,:)
  
  !Intermidiate (Depth is only woa)
  real(kind = 8),allocatable :: dat_int(:,:,:)

  !Read Model information
  write(*,*) "Read Model grid information"
  allocate(null1x(im),null1y(jm),null2d(im,jm))
  allocate(null3d(im,jm,km))
  call read_grid(null1x,null1x,lon,null1x, &
       & null1y,lat,null1y,null1y, &
       & null3d,depth,fsm,null2d,null2d)
  deallocate(null1x,null1y,null2d)
  deallocate(null3d)

  write(*,'(a)') "Calculate ID"
  allocate(depth_woa(kms_woa),dat_woa(im_woa,jm_woa,kms_woa))
  call read_woa18_season_annual_var("t",0,lon_woa,lat_woa,depth_woa,dat_woa)
  call cal_idlon(im_woa,lon_woa,im,lon,idx)
  call cal_idlat(jm_woa,lat_woa,jm,lat,idy)
  call cal_idz(kms_woa,depth_woa,im,jm,km,depth,idz_s)
  deallocate(depth_woa,dat_woa)
  
  allocate(depth_woa(kmm_woa),dat_woa(im_woa,jm_woa,kmm_woa))
  call read_woa18_month_var("t",1,lon_woa,lat_woa,depth_woa,dat_woa)
  call cal_idz(kmm_woa,depth_woa,im,jm,km,depth,idz_m)
  deallocate(depth_woa,dat_woa)
  
  
  if(lon(1) == lon(im-1)-360.d0 .and. lon(2) == lon(im)-360.d0)then
     iqglobal=1
  else
     iqglobal=0
  end if

  do imon=0,12
     do ivar=1,3
        
        !---Annual/Seasonal
        allocate(depth_woa(kms_woa))
        
        if(ivar == 1) varname="t"
        if(ivar == 2) varname="s"
        if(ivar == 3) varname="r"
        
        write(*,'(a,i2,a)') "----- Start Month and variable: ",imon," "//varname//"-----"
        
        allocate(dat_woa(im_woa,jm_woa,kms_woa))

        write(*,'(a)') "Read Seasonal/Annual climatology"
        if(ivar == 3)then
           allocate(t_woa(im_woa,jm_woa,kms_woa),s_woa(im_woa,jm_woa,kms_woa))
           call read_woa18_season_annual_var("t",imon,lon_woa,lat_woa,depth_woa,t_woa)
           call read_woa18_season_annual_var("s",imon,lon_woa,lat_woa,depth_woa,s_woa)
           call density_woa(im_woa,jm_woa,kms_woa,depth_woa,t_woa,s_woa,dat_woa)
           call horizontal_ave(im_woa,jm_woa,kms_woa,dat_woa)
           deallocate(t_woa,s_woa)
        else
           call read_woa18_season_annual_var(varname,imon,lon_woa,lat_woa,depth_woa,dat_woa)
        end if
           
        write(*,'(a)') "Fill value"
        call fillvalue_2d(10,im_woa,jm_woa,kms_woa,dat_woa,rmiss)

        write(*,'(a)') "Bilinear Interpolation in horizontal direction"        
        allocate(dat_int(im,jm,kms_woa))
        do k=1,kms_woa 
           call bilinear_interpolation_2d( &
                & im_woa,jm_woa,lon_woa,lat_woa,dat_woa(:,:,k), &
                & im,jm,lon,lat,dat_int(:,:,k),idx,idy,rmiss)
        end do
        call apply_fsm(im,jm,kms_woa,dat_int,fsm)
        deallocate(dat_woa)

        write(*,'(a)') "Interpoation in vertical direction"        
        call linear_interpolation_vertical_full( &
             & kms_woa,depth_woa,dat_int,im,jm,km,depth,fsm,dat,idz_s,rmiss)
        deallocate(dat_int)

        write(*,'(a)') "Fill value in vertical direction"
        call fillvalue_vertical(im,jm,km,fsm,dat,rmiss)
        call apply_fsm(im,jm,km,dat,fsm)
        
        deallocate(depth_woa)

        !---Monthly
        if(imon /= 0)then

           allocate(depth_woa(kmm_woa))
           allocate(dat_woa(im_woa,jm_woa,kmm_woa))
           
           write(*,*) "Read Monthly climatology"
           if(ivar == 3)then
              allocate(t_woa(im_woa,jm_woa,kmm_woa),s_woa(im_woa,jm_woa,kmm_woa))              
              call read_woa18_month_var("t",imon,lon_woa,lat_woa,depth_woa,t_woa)
              call read_woa18_month_var("s",imon,lon_woa,lat_woa,depth_woa,s_woa)
              call density_woa(im_woa,jm_woa,kmm_woa,depth_woa,t_woa,s_woa,dat_woa)
              call horizontal_ave(im_woa,jm_woa,kmm_woa,dat_woa)
              deallocate(t_woa,s_woa)
           else
              call read_woa18_month_var(varname,imon,lon_woa,lat_woa,depth_woa,dat_woa)
           end if
              
           write(*,'(a)') "Fill value"
           call fillvalue_2d(10,im_woa,jm_woa,kmm_woa,dat_woa,rmiss)
           
           allocate(dat_int(im,jm,kmm_woa))
           write(*,'(a)') "Bilinear Interpolation in horizontal direction"
           do k=1,kmm_woa
              call bilinear_interpolation_2d( &
                   & im_woa,jm_woa,lon_woa,lat_woa,dat_woa(:,:,k), &
                   & im,jm,lon,lat,dat_int(:,:,k),idx,idy,rmiss)
           end do
           call apply_fsm(im,jm,kmm_woa,dat_int,fsm)

           deallocate(dat_woa)

           write(*,'(a)') "Interpoation in vertical direction"
           call linear_interpolation_vertical_part( &
                & kmm_woa,depth_woa,dat_int,im,jm,km,depth,fsm,dat,idz_m,rmiss)
           
           deallocate(dat_int)
           deallocate(depth_woa)
           
        end if !imon

        if(iqglobal == 1)then
           dat(1:2,:,:)=dat(im-1:im,:,:)
        end if

        write(*,'(a)') "Write data"        
        call write_tsclim(varname,imon,im,jm,km,dat)

        if(varname == "t" .or. varname == "s")then
           call write_ic(varname,imon,im,jm,km,dat)
        end if
           
        write(*,'(a,i2,a)') "----- End Month and variable: ",imon," "//varname//"-----"
        
     end do !ivar
  end do !imon
  
end program main

!-----------------------------------------------------------------------------
! POM Density |
!-----------------------------------------------------------------------------

subroutine density_woa(im,jm,kb,depth,t,s,rho)

  ! calculate (density-1000.)/rhoref.
  ! see: Mellor, G.L., 1991, J. Atmos. Oceanic Tech., 609-611
  ! note: if pressure is not used in dens, buoyancy term (boygr) in
  ! subroutine profq must be changed (see note in subroutine profq)

  !$use omp_lib
  use setting
  use mod_rmiss
  implicit none

  integer i,j,k
  real(kind = 8) cr,p,rhor,sr,tr,tr2,tr3,tr4

  integer,intent(in) :: im,jm,kb

  real(kind = 8),intent(in) :: depth(kb)
  real(kind = 8),intent(in)  :: t(im,jm,kb),s(im,jm,kb)
  real(kind = 8),intent(out) :: rho(im,jm,kb)

  rho(:,:,:)=0.d0

  !$omp parallel
  !$omp do private(i,j,k,tr,sr,tr2,tr3,tr4,p,rhor,cr)  
  do k=1,kb-1
     do j=1,jm
        do i=1,im

           if(t(i,j,k) == rmiss .or. s(i,j,k) == rmiss)then
              rho(i,j,k)=rmiss
              cycle
           end if
           
           tr=t(i,j,k)
           sr=s(i,j,k)
           tr2=tr*tr
           tr3=tr2*tr
           tr4=tr3*tr

           ! approximate pressure in units of bars
           !p=grav*rhoref*(-zz(k)* h(i,j))*1.d-5
           p=grav*rhoref*depth(k)*1.d-5

           rhor=-0.157406d0+6.793952d-2*tr         &
                & -9.095290d-3*tr2+1.001685d-4*tr3 &
                & -1.120083d-6*tr4+6.536332d-9*tr4*tr

           rhor=rhor &
                & +(0.824493d0-4.0899d-3*tr+7.6438d-5*tr2-8.2467d-7*tr3+5.3875d-9*tr4)*sr &
                & +(-5.72466d-3+1.0227d-4*tr-1.6546d-6*tr2)*abs(sr)**1.5d0 &
                & +4.8314d-4*sr*sr

           cr=1449.1d0+0.0821d0*p+4.55d0*tr-0.045d0*tr2+1.34d0*(sr-35.d0)
           rhor=rhor+1.d5*p/(cr*cr)*(1.d0-2.d0*p/(cr*cr))

           !rho(i,j,k)=rhor/rhoref*fsm(i,j)
           rho(i,j,k)=rhor/rhoref
           
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine density_woa

!-----------------------------------------------------------------------------
! Horizontal average |
!-----------------------------------------------------------------------------

subroutine horizontal_ave(im,jm,km,dat)

  !$use omp_lib  
  use mod_rmiss
  implicit  none

  integer i,j,k
  integer ipass,imiss
  real(kind = 8) ave(km)
  
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(inout) :: dat(im,jm,km)

  do k=1,km

     ave(k)=0.d0
     ipass=0
     imiss=0
     
     do j=1,jm
        do i=1,im
           if(dat(i,j,k) == rmiss)then
              imiss=imiss+1
           else
              ave(k)=ave(k)+dat(i,j,k)
              ipass=ipass+1
           end if
        end do
     end do

     if(ipass == 0)then
        ave(k)=rmiss
     else
        ave(k)=ave(k)/dble(ipass)
     end if
     
  end do

  !$omp parallel
  !$omp do private(i,j,k)  
  do k=1,km
     do j=1,jm
        do i=1,im
           dat(i,j,k)=ave(k)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine horizontal_ave

!-----------------------------------------------------------------------------
! Linear interpolation |
!-----------------------------------------------------------------------------

subroutine linear_interpolation_vertical_full( &
     & km1,depth1,dat1,im2,jm2,km2,depth2,fsm,dat2,id,rmiss)

  !$use omp_lib 
  implicit none
  
  integer k1
  integer i2,j2,k2
  
  integer,intent(in) :: km1  
  integer,intent(in) :: im2,jm2,km2
  integer,intent(in) :: id(im2,jm2,km2)

  real(kind = 8),intent(in) :: depth1(km1),dat1(im2,jm2,km1)
  real(kind = 8),intent(in) :: depth2(im2,jm2,km2),fsm(im2,jm2)
  real(kind = 8),intent(in) :: rmiss

  real(kind = 8),intent(out) :: dat2(im2,jm2,km2)
  
  dat2(:,:,:)=rmiss

  !$omp parallel
  !$omp do private(k1,i2,j2,k2)
  do j2=1,jm2
     do i2=1,im2

        if(fsm(i2,j2)==0.)then
           dat2(i2,j2,:)=0.
           cycle
        endif

        do k2=1,km2
           
           if(id(i2,j2,k2)==0)then
              cycle
           else
              k1=id(i2,j2,k2)
           end if
           
           if(dat1(i2,j2,k1) == rmiss .or. dat1(i2,j2,k1+1) == rmiss)then
              dat2(i2,j2,k2)=rmiss
           else
              call linear_interpolate( &
                   & depth1(k1),depth1(k1+1),dat1(i2,j2,k1),dat1(i2,j2,k1+1), &
                   & depth2(i2,j2,k2),dat2(i2,j2,k2))
           end if

        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine linear_interpolation_vertical_full

!------------------------------------------------

subroutine linear_interpolation_vertical_part( &
     & km1,depth1,dat1,im2,jm2,km2,depth2,fsm,dat2,id,rmiss)

  !$use omp_lib 
  implicit none
  
  integer k1
  integer i2,j2,k2
  
  integer,intent(in) :: km1  
  integer,intent(in) :: im2,jm2,km2
  integer,intent(in) :: id(im2,jm2,km2)

  real(kind = 8),intent(in) :: depth1(km1),dat1(im2,jm2,km1)
  real(kind = 8),intent(in) :: depth2(im2,jm2,km2),fsm(im2,jm2)
  real(kind = 8),intent(in) :: rmiss

  real(kind = 8),intent(out) :: dat2(im2,jm2,km2)
  
  !$omp parallel
  !$omp do private(k1,i2,j2,k2)
  do j2=1,jm2
     do i2=1,im2

        if(fsm(i2,j2)==0.)then
           dat2(i2,j2,:)=0.
           cycle
        endif

        do k2=1,km2
           
           if(id(i2,j2,k2)==0)then
              cycle
           else
              k1=id(i2,j2,k2)
           end if
           
           if(dat1(i2,j2,k1) == rmiss .or. dat1(i2,j2,k1+1) == rmiss)then
              cycle
           else
              call linear_interpolate( &
                   & depth1(k1),depth1(k1+1),dat1(i2,j2,k1),dat1(i2,j2,k1+1), &
                   & depth2(i2,j2,k2),dat2(i2,j2,k2))
           end if

        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine linear_interpolation_vertical_part

!----------------------------------------------------------------------
! Combine Monthly & Seasonal WOA dataset |
!-----------------------------------------
!
! Monthly: 0-1500m
! Seasonal: 0-5500m
!
!----------------------------------------------------------------------

subroutine combine_ms_woa(im,jm,km,fsm,mdat,sdat,dat,rmiss)

  !$use omp_lib  
  implicit none

  integer i,j,k

  integer,intent(in) :: im,jm,km
  
  real(kind = 8),intent(in) :: fsm(im,jm)
  real(kind = 8),intent(in) :: mdat(im,jm,km) !Monthly dataset
  real(kind = 8),intent(in) :: sdat(im,jm,km) !Seasonal dataset
  real(kind = 8),intent(in) :: rmiss

  real(kind = 8),intent(out) :: dat(im,jm,km)
  
  dat(:,:,:)=sdat(:,:,:)

  !$omp parallel
  !$omp do private(i,j,k)  
  do j=1,jm
     do i=1,im

        if(fsm(i,j) == 0.d0)then
           dat(i,j,:)=0.d0
           cycle
        end if

        do k=1,km           
           if(mdat(i,j,k) /= rmiss)then
              dat(i,j,k)=mdat(i,j,k)
           end if
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine combine_ms_woa

!-----------------------------------------------------------------------
! Write Data |
!-----------------------------------------------------------------------

!tsclim: annual/monthly climtology file
subroutine write_tsclim(varname,imon,im,jm,km,dat)

  use netcdf
  implicit none

  integer status,access,system
  integer ncid,varid

  real(kind = 4) tmp(im,jm,km)

  character(100) :: filename="../in/tsclim.nc"

  !IN
  integer,intent(in) :: imon
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: dat(im,jm,km)

  character(1),intent(in) :: varname
  
  !Check file
  status=access(trim(filename)," ")

  if(imon==0 .and. varname == "t" .and. status == 0)then
     status=system("rm -f "//trim(filename))
  endif

  if(imon==0 .and. varname == "t")then
     call make_ncfile_tsclim(im,jm,km,filename)
  !   status=system("ncgen -o "//trim(filename)//" ../hdr/tsclim.hdr")
  end if

  if(imon /= 0 .and. status /= 0)then
     write(*,'(a)') "***Error: Not Found "//trim(filename)
     stop
  end if

  status=nf90_open(trim(filename),nf90_write,ncid)

  tmp(:,:,:)=real(dat(:,:,:))
  if(imon == 0)then
     if(varname == "t") status=nf90_inq_varid(ncid,"tclim",varid)
     if(varname == "s") status=nf90_inq_varid(ncid,"sclim",varid)
     if(varname == "r") status=nf90_inq_varid(ncid,"rmean",varid)
     status=nf90_put_var(ncid,varid,tmp)
  else
     if(varname == "t") status=nf90_inq_varid(ncid,"tclimm",varid)
     if(varname == "s") status=nf90_inq_varid(ncid,"sclimm",varid)     
     if(varname == "r") status=nf90_inq_varid(ncid,"rmeanm",varid)     
     status=nf90_put_var(ncid,varid,tmp,start=(/1,1,1,imon/),count=(/im,jm,km,1/))
  end if

  status=nf90_close(ncid)

end subroutine write_tsclim

!----------------------------

!ic: Initial condition file
subroutine write_ic(varname,imon,im,jm,km,dat)

  use netcdf
  implicit none

  integer status,access,system
  integer ncid,varid

  real(kind = 8) null(im,jm,km)

  character(2) mm
  character(100) filename
  
  !IN
  integer,intent(in) :: imon
  integer,intent(in) :: im,jm,km
  
  real(kind = 8),intent(in) :: dat(im,jm,km)

  character(1),intent(in) :: varname
  
  write(mm,'(i2.2)') imon
  filename="../in/ic.woa18."//mm//".nc"

  !Because of numerical instability, horizontal velocities are set to be zero.
  null(:,:,:)=0.d0

  !Check file
  status=access(trim(filename)," ")
  if(status == 0 .and. varname == "t")then
     status=system("rm -f "//trim(filename))

  end if

  !Make file
  status=access(trim(filename)," ")
  !status=system("ncgen -o "//trim(filename)//" ../hdr/ic.hdr")
  if(status /= 0)then
     call make_ncfile_ic(im,jm,km,filename)
  end if
  
  !Write data
  status=nf90_open(trim(filename),nf90_write,ncid)

  status=nf90_inq_varid(ncid,varname,varid)
  status=nf90_put_var(ncid,varid,dat)

  if(varname == "t")then
     status=nf90_inq_varid(ncid,"u",varid)
     status=nf90_put_var(ncid,varid,null)

     status=nf90_inq_varid(ncid,"v",varid)
     status=nf90_put_var(ncid,varid,null)
  end if
     
  status=nf90_close(ncid)
  
end subroutine write_ic

!----------------------------------------------------------------------------------
