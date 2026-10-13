!---------------------------------------------------------------------
! Write Observational data at sea surface |
!---------------------------------------------------------------------

subroutine write_obs_surface(var,iyr,imon,iday,inum,lon,lat,dat)

  use mod_rmiss
  use mod_gridinfo
  use setting, only : &
       & sst_min, sst_max, sst_err, sst_ele, sst_dx, sst_dy, &
       & sss_min, sss_max, sss_err, sss_ele, sss_dx, sss_dy, &
       & ssh_min, ssh_max, ssh_err, ssh_ele, ssh_dx, ssh_dy, &
       & ssu_min, ssu_max, ssu_err, ssu_ele, ssu_dx, ssu_dy, &
       & ssv_min, ssv_max, ssv_err, ssv_ele, ssv_dx, ssv_dy, &
       & ini_filename
  use netcdf
  implicit none

  !Common
  real(kind = 8),parameter :: lev(1)=0.d0

  integer i,j,di,dj
  integer status,access,system
  integer ncid,varid

  real(kind = 8) err(im,jm)
  real(kind = 8) min,max,ele(1)

  character(100) filename
  character(8) yyyymmdd
  character(4) yyyy
  character(2) mm,dd

  !IN
  integer,intent(in) :: iyr,imon,iday
  real(kind = 8),intent(in) :: lon(im),lat(jm),dat(im,jm)
  character(3),intent(in) :: var

  !INOUT
  integer,intent(inout) :: inum
  
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday
  yyyymmdd=yyyy//mm//dd
  
  filename="../obs/"//trim(ini_filename)//yyyymmdd//".nc"

  status=access(trim(filename)," ")
  if(status /= 0)then
     call make_ncfile_obs(filename)
  end if

  if(var == "sst")then
     min=sst_min
     max=sst_max
     err(:,:)=sst_err
     ele(:)=sst_ele
     di=sst_dx
     dj=sst_dy
  else if(var == "sss")then
     min=sss_min
     max=sss_max
     err(:,:)=sss_err
     ele(:)=sss_ele     
     di=sss_dx
     dj=sss_dy
  else if(var == "ssh")then
     min=ssh_min
     max=ssh_max
     err(:,:)=ssh_err
     ele(:)=ssh_ele
     di=ssh_dx
     dj=ssh_dy
  else if(var == "ssu")then
     min=ssu_min
     max=ssu_max
     err(:,:)=ssu_err
     ele(:)=ssu_ele
     di=ssu_dx
     dj=ssu_dy
  else if(var == "ssv")then
     min=ssv_min
     max=ssv_max
     err(:,:)=ssv_err
     ele(:)=ssv_ele
     di=ssv_dx
     dj=ssv_dy
  else
     write(*,'(a)') "***Error: Choose sst/sss/ssh/ssu/ssv"
     stop
  end if

  status=nf90_open(trim(filename),nf90_write,ncid)

  do j=2,jm-1,dj
     do i=2,im-1,di
        
        if(dat(i,j) == rmiss .or. dat(i,j) /= dat(i,j) &
             & .or. dat(i,j) < min .or. max < dat(i,j)) cycle 
        
        inum=inum+1

        status=nf90_inq_varid(ncid,"ele",varid)
        status=nf90_put_var(ncid,varid,real(ele),(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"lon",varid)
        status=nf90_put_var(ncid,varid,real(lon(i:i)),(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"lat",varid)
        status=nf90_put_var(ncid,varid,real(lat(j:j)),(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"lev",varid)
        status=nf90_put_var(ncid,varid,real(lev),(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"dat",varid)
        status=nf90_put_var(ncid,varid,real(dat(i:i,j:j)),(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"err",varid)
        status=nf90_put_var(ncid,varid,real(err(i:i,j:j)),(/inum/),(/1/))

     end do
  end do
  status=nf90_close(ncid)

end subroutine write_obs_surface

!----------------------------------------------------------------------------
! Write Temperature/Salinity Observation |
!----------------------------------------------------------------------------

subroutine write_obs_ts(var,iyr,imon,iday,inum,km,lon,lat,depth,dat,err)

  use setting,only: &
       & t_min,t_max,t_err,t_ele, &
       & s_min,s_max,s_err,s_ele, &
       & ini_filename
  use mod_rmiss
  use netcdf
  implicit none

  integer k
  integer status,access,system
  integer ncid,varid

  real(kind = 8) min,max,ele(1)
  real(kind = 8) mdep(km)

  character(100) filename
  character(4) yyyy
  character(2) mm,dd
  character(1) var
  
  !IN
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: km
  real(kind = 8),intent(in) :: lon(1),lat(1),depth(km)
  real(kind = 8),intent(in) :: dat(km),err(km)

  !INOUT
  integer,intent(inout) :: inum


  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday

  !Filename
  filename="../obs/"//trim(ini_filename)//yyyy//mm//dd//".nc"

  status=access(trim(filename)," ")
  if(status /= 0)then
     call make_ncfile_obs(filename)
  end if

  !Setting
  if(var == "T")then
     min=t_min
     max=t_max
     ele(:)=t_ele
  elseif(var == "S")then
     min=s_min
     max=s_max
     ele(:)=s_ele
  else
     write(*,'(a)') "***Error: Choose T or S"
     stop
  end if

  mdep(:)=-1.d0*abs(depth(:))
  
  !Write Data
  status=nf90_open(trim(filename),nf90_write,ncid)
  do k=1,km

     if(dat(k) == rmiss .or. dat(k) /= dat(k))cycle
     if(dat(k) < min .or. max < dat(k))cycle 
     
     inum=inum+1

     status=nf90_inq_varid(ncid,"ele",varid)
     status=nf90_put_var(ncid,varid,real(ele),(/inum/),(/1/))
     
     status=nf90_inq_varid(ncid,"lon",varid)
     status=nf90_put_var(ncid,varid,real(lon),(/inum/),(/1/))
     
     status=nf90_inq_varid(ncid,"lat",varid)
     status=nf90_put_var(ncid,varid,real(lat),(/inum/),(/1/))

     status=nf90_inq_varid(ncid,"lev",varid)
     status=nf90_put_var(ncid,varid,real(mdep(k:k)),(/inum/),(/1/))

     status=nf90_inq_varid(ncid,"dat",varid)
     status=nf90_put_var(ncid,varid,real(dat(k:k)),(/inum/),(/1/))

     status=nf90_inq_varid(ncid,"err",varid)
     status=nf90_put_var(ncid,varid,real(err(k:k)),(/inum/),(/1/))

  end do
  status=nf90_close(ncid)

end subroutine write_obs_ts

!---------------------------------------------------------------------
! Write Observational data at sea surface uv|
!---------------------------------------------------------------------

subroutine write_obs_surface_uv(var,iyr,imon,iday,inum,lon,lat,dat,mask)

  use mod_rmiss
  use mod_gridinfo
  use setting, only : &
       & ssu_min, ssu_max, ssu_err, ssu_ele, ssu_dx, ssu_dy, &
       & ssv_min, ssv_max, ssv_err, ssv_ele, ssv_dx, ssv_dy, &
       & ini_filename
  use netcdf
  implicit none

  real(kind = 8),parameter :: lev(1)=0.d0

  integer,intent(in) :: iyr,imon,iday
  integer,intent(inout) :: inum
  integer i,j,di,dj
  integer status,access,system
  integer ncid,varid

  real(kind = 8),intent(in) :: lon(im),lat(jm)
  real(kind = 8),intent(in) :: dat(im,jm),mask(im,jm)
  real(kind = 8) err(im,jm)
  real(kind = 8) min,max,ele(1)

  character(3),intent(in) :: var
  character(100) filename
  character(8) yyyymmdd
  character(4) yyyy
  character(2) mm,dd
  
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday
  yyyymmdd=yyyy//mm//dd

  filename="../obs/"//trim(ini_filename)//yyyymmdd//".nc"

  status=access(trim(filename)," ")
  if(status /= 0)then
     call make_ncfile_obs(filename)
  end if

  if(var == "ssu")then
     min=ssu_min
     max=ssu_max
     err(:,:)=ssu_err
     ele(:)=ssu_ele
     di=ssu_dx
     dj=ssu_dy
  elseif(var == "ssv")then
     min=ssv_min
     max=ssv_max
     err(:,:)=ssv_err
     ele(:)=ssv_ele     
     di=ssv_dx
     dj=ssv_dy
  else
     write(*,'(a)') "***Error: Choose ssu or ssv"
     stop
  end if

  status=nf90_open(trim(filename),nf90_write,ncid)

  do j=2,jm-1,dj
     do i=2,im-1,di
        
        if(mask(i,j) == 0.d0 .or. dat(i,j) /= dat(i,j))cycle
        if(dat(i,j) < min .or. max < dat(i,j))cycle 
        
        inum=inum+1

        status=nf90_inq_varid(ncid,"ele",varid)
        status=nf90_put_var(ncid,varid,real(ele),(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"lon",varid)
        status=nf90_put_var(ncid,varid,real(lon(i:i)),(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"lat",varid)
        status=nf90_put_var(ncid,varid,real(lat(j:j)),(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"lev",varid)
        status=nf90_put_var(ncid,varid,real(lev),(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"dat",varid)
        status=nf90_put_var(ncid,varid,real(dat(i:i,j:j)),(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"err",varid)
        status=nf90_put_var(ncid,varid,real(err(i:i,j:j)),(/inum/),(/1/))

     end do
  end do
  status=nf90_close(ncid)

end subroutine write_obs_surface_uv

