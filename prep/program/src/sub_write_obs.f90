!----------------------------------------------------------------------------
! Get observation information
!----------------------------------------------------------------------------

subroutine get_obs_surface_info(varname,min,max,err,ele)

  use setting, only: &
       & ssu_min, ssu_max, ssu_err, & 
       & ssv_min, ssv_max, ssv_err, & 
       & sst_min, sst_max, sst_err, & 
       & sss_min, sss_max, sss_err, & 
       & ssh_min, ssh_max, ssh_err, &
       & ssu_ele, ssv_ele, sst_ele, ssh_ele, sss_ele

  implicit none

  !---IN
  character(*),intent(in) :: varname

  !---OUT
  integer,intent(out) :: ele
  real(kind = 8),intent(out) :: min,max,err

  !---Min,Max,Err,Ele
  if(varname == "ssu")then
     min=ssu_min
     max=ssu_max
     err=ssu_err
     ele=ssu_ele     
  else if(varname == "ssv")then
     min=ssv_min
     max=ssv_max
     err=ssv_err
     ele=ssv_ele     
  else if(varname == "sst")then
     min=sst_min
     max=sst_max
     err=sst_err
     ele=sst_ele     
  else if(varname == "sss")then
     min=sss_min
     max=sss_max
     err=sss_err
     ele=sss_ele     
  else if(varname == "ssh")then
     min=ssh_min
     max=ssh_max
     err=ssh_err
     ele=ssh_ele
  else
     write(*,*) "***Error: varname => "//trim(varname)
     stop
  end if

end subroutine get_obs_surface_info

!-----------------------------------

subroutine get_obs_interior_info(varname,min,max,err,ele)

  use setting, only: &
       & t_min, t_max, t_err, & 
       & s_min, s_max, s_err, & 
       & u_min, u_max, u_err, & 
       & v_min, v_max, v_err, & 
       & t_ele, s_ele, u_ele, v_ele

  implicit none

  !---IN
  character(*),intent(in) :: varname

  !---OUT
  integer,intent(out) :: ele
  real(kind = 8),intent(out) :: min,max,err

  !---Min,Max,Err,Ele
  if(varname == "t")then
     min=t_min
     max=t_max
     err=t_err
     ele=t_ele     
  else if(varname == "s")then
     min=s_min
     max=s_max
     err=s_err
     ele=s_ele     
  else if(varname == "u")then
     min=u_min
     max=u_max
     err=u_err
     ele=u_ele     
  else if(varname == "v")then
     min=v_min
     max=v_max
     err=v_err
     ele=v_ele     
  else
     write(*,*) "***Error: varname => "//trim(varname)
     stop
  end if

end subroutine get_obs_interior_info

!-------------------------------------------------------
! Write Data
!-------------------------------------------------------
! Note:
!  surface2d: lon(im,jm),lat(im,jm)
!  surface2dg: lon(im),lat(jm)
!  surface1d: lon(ntime),lat(ntime)
!-------------------------------------------------------

subroutine write_obs_common(filename,so,no,ele,ins,lon,lat,lev,dat,err)

  use mod_check_netcdf
  use netcdf
  implicit none

  !---Common
  integer status,access
  integer ncid,varid

  !---IN
  integer,intent(in) :: so,no
  integer,intent(in) :: ele(no),ins(no)

  real(kind = 8),intent(in) :: lon(no),lat(no),lev(no)
  real(kind = 8),intent(in) :: dat(no),err(no)

  character(200),intent(in) :: filename

  !---Check access
  status=access(trim(filename)," ")
  if(status /= 0)then
     write(*,*) "***Error: Not found "//trim(filename)
     stop
  end if
  
  !---Write data
  status=nf90_open(trim(filename),nf90_write,ncid)
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"ele",varid)
  call check_netcdf(status)
  status=nf90_put_var(ncid,varid,ele,(/so/),(/no/))
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"ins",varid)
  call check_netcdf(status)
  status=nf90_put_var(ncid,varid,ins,(/so/),(/no/))        
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"lon",varid)
  call check_netcdf(status)
  status=nf90_put_var(ncid,varid,real(lon),(/so/),(/no/))
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"lat",varid)
  call check_netcdf(status)
  status=nf90_put_var(ncid,varid,real(lat),(/so/),(/no/))
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"lev",varid)
  call check_netcdf(status)
  status=nf90_put_var(ncid,varid,real(lev),(/so/),(/no/))
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"dat",varid)
  call check_netcdf(status)
  status=nf90_put_var(ncid,varid,real(dat),(/so/),(/no/))
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"err",varid)
  call check_netcdf(status)
  status=nf90_put_var(ncid,varid,real(err),(/so/),(/no/))
  call check_netcdf(status)

  status=nf90_close(ncid)
  call check_netcdf(status)

end subroutine write_obs_common

!-------------------------------------------------------

subroutine write_obs_surface2d(varname,ins,iyr,imon,iday,&
     & idx,idy,im,jm,lons,lone,lon,lats,late,lat,dat,inum)

  use setting, only: ini_filename
  use mod_rmiss
  use mod_check_netcdf
  use netcdf
  implicit none

  !---Parameter
  real(kind = 8),parameter :: lev=0.d0

  !---Common
  integer i,j
  integer ele
  integer so,no
  integer status,access

  integer,allocatable :: ele_tmp(:),ins_tmp(:)

  real(kind = 8) min,max,err

  real(kind = 8),allocatable :: lon_tmp(:),lat_tmp(:),lev_tmp(:)
  real(kind = 8),allocatable :: dat_tmp(:),err_tmp(:)

  character(200) filename
  character(8) yyyymmdd
  character(4) yyyy
  character(2) mm,dd

  !---IN
  integer,intent(in) :: ins
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: idx,idy
  integer,intent(in) :: im,jm

  real(kind = 8),intent(in) :: lons,lone,lon(im,jm)
  real(kind = 8),intent(in) :: lats,late,lat(im,jm)
  real(kind = 8),intent(in) :: dat(im,jm)

  character(3),intent(in) :: varname

  !---INOUT
  integer,intent(inout) :: inum

  !---Filename
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday
  yyyymmdd=yyyy//mm//dd

  filename="../obs/"//trim(ini_filename)//yyyymmdd//".nc"

  call get_obs_surface_info(varname,min,max,err,ele)

  !---Make file
  status=access(trim(filename)," ")
  if(status /= 0)then
     call make_ncfile_obs(filename)
  end if

  !---Count obs.
  so=inum+1
  no=0
  do j=1,jm,idy
     do i=1,im,idx

        if(dat(i,j) == rmiss .or. dat(i,j) < min .or. max < dat(i,j)) cycle
        if(lon(i,j) < lons .or. lone < lon(i,j)) cycle
        if(lat(i,j) < lats .or. late < lat(i,j)) cycle

        inum=inum+1
        no=no+1

     end do
  end do

  !---Identify no observation 
  if(no < 1)then
     return
  end if

  !---Get data
  allocate(ele_tmp(no),ins_tmp(no))
  allocate(lon_tmp(no),lat_tmp(no),lev_tmp(no))
  allocate(dat_tmp(no),err_tmp(no))

  no=0
  ele_tmp(:)=ele
  ins_tmp(:)=ins
  lev_tmp(:)=lev
  err_tmp(:)=err
  do j=1,jm,idy
     do i=1,im,idx

        if(dat(i,j) == rmiss .or. dat(i,j) < min .or. max < dat(i,j)) cycle
        if(lon(i,j) < lons .or. lone < lon(i,j)) cycle
        if(lat(i,j) < lats .or. late < lat(i,j)) cycle

        no=no+1
        lon_tmp(no)=lon(i,j)
        lat_tmp(no)=lat(i,j)
        dat_tmp(no)=dat(i,j)

     end do
  end do

  !---Write data
  call write_obs_common(filename,so,no,ele_tmp,ins_tmp,lon_tmp,lat_tmp,lev_tmp,dat_tmp,err_tmp)

  !---Deallocate
  deallocate(ele_tmp,ins_tmp)
  deallocate(lon_tmp,lat_tmp,lev_tmp)
  deallocate(dat_tmp,err_tmp)

end subroutine write_obs_surface2d

!-----------------------------------

subroutine write_obs_surface2dg(varname,ins,iyr,imon,iday, &
     & idx,idy,im,jm,lons,lone,lon,lats,late,lat,dat,inum)

  use setting, only: ini_filename
  use mod_rmiss
  use mod_check_netcdf
  use netcdf
  implicit none

  !---Parameter
  real(kind = 8),parameter :: lev=0.d0

  !---Common
  integer i,j
  integer so,no
  integer ele
  integer status,access

  integer,allocatable :: ele_tmp(:),ins_tmp(:)

  real(kind = 8) min,max,err

  real(kind = 8),allocatable :: lon_tmp(:),lat_tmp(:),lev_tmp(:)
  real(kind = 8),allocatable :: dat_tmp(:),err_tmp(:)

  character(200) filename  
  character(8) yyyymmdd
  character(4) yyyy
  character(2) mm,dd

  !---IN
  integer,intent(in) :: ins
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: idx,idy
  integer,intent(in) :: im,jm

  real(kind = 8),intent(in) :: lons,lone,lon(im)
  real(kind = 8),intent(in) :: lats,late,lat(jm)
  real(kind = 8),intent(in) :: dat(im,jm)

  character(3),intent(in) :: varname

  !---INOUT
  integer,intent(inout) :: inum

  !---Filename
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday
  yyyymmdd=yyyy//mm//dd

  filename="../obs/"//trim(ini_filename)//yyyymmdd//".nc"

  call get_obs_surface_info(varname,min,max,err,ele)

  !---Make file
  status=access(trim(filename)," ")
  if(status /= 0)then
     call make_ncfile_obs(filename)
  end if

  !---Count obs
  so=inum+1
  no=0
  do j=1,jm,idy

     if(lat(j) < lats .or. late < lat(j)) cycle        

     do i=1,im,idx

        if(dat(i,j) == rmiss .or. dat(i,j) < min .or. max < dat(i,j)) cycle
        if(lon(i) < lons .or. lone < lon(i)) cycle

        inum=inum+1
        no=no+1

     end do
  end do

  !---Identify no observation 
  if(no < 1)then
     return
  end if

  !---Get data
  allocate(ele_tmp(no),ins_tmp(no))
  allocate(lon_tmp(no),lat_tmp(no),lev_tmp(no))
  allocate(dat_tmp(no),err_tmp(no))

  no=0
  ele_tmp(:)=ele
  ins_tmp(:)=ins
  lev_tmp(:)=lev
  err_tmp(:)=err
  do j=1,jm,idy

     if(lat(j) < lats .or. late < lat(j)) cycle        

     do i=1,im,idx

        if(dat(i,j) == rmiss .or. dat(i,j) < min .or. max < dat(i,j)) cycle
        if(lon(i) < lons .or. lone < lon(i)) cycle

        no=no+1
        lon_tmp(no)=lon(i)
        lat_tmp(no)=lat(j)
        dat_tmp(no)=dat(i,j)

     end do
  end do

  !---Write data
  call write_obs_common(filename,so,no,ele_tmp,ins_tmp,lon_tmp,lat_tmp,lev_tmp,dat_tmp,err_tmp)

  !---Deallocate
  deallocate(ele_tmp,ins_tmp)
  deallocate(lon_tmp,lat_tmp,lev_tmp)
  deallocate(dat_tmp,err_tmp)

end subroutine write_obs_surface2dg

!--------------------------------------

subroutine write_obs_surface1d(varname,ins,iyr,imon,iday,&
     & idx,im,lons,lone,lon,lats,late,lat,dat,inum)

  use setting, only: ini_filename
  use mod_rmiss
  use mod_check_netcdf
  use netcdf
  implicit none

  !---Parameter
  real(kind = 8),parameter :: lev=0.d0

  !---Common
  integer i
  integer so,no
  integer ele
  integer status,access

  integer,allocatable :: ele_tmp(:),ins_tmp(:)

  real(kind = 8) min,max,err

  real(kind = 8),allocatable :: lon_tmp(:),lat_tmp(:),lev_tmp(:)
  real(kind = 8),allocatable :: dat_tmp(:),err_tmp(:)

  character(200) filename
  character(8) yyyymmdd
  character(4) yyyy
  character(2) mm,dd

  !---IN
  integer,intent(in) :: ins
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: idx
  integer,intent(in) :: im

  real(kind = 8),intent(in) :: lons,lone,lon(im)
  real(kind = 8),intent(in) :: lats,late,lat(im)
  real(kind = 8),intent(in) :: dat(im)

  character(3),intent(in) :: varname

  !---INOUT
  integer,intent(inout) :: inum

  !---Filename
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday
  yyyymmdd=yyyy//mm//dd

  filename="../obs/"//trim(ini_filename)//yyyymmdd//".nc"

  call get_obs_surface_info(varname,min,max,err,ele)

  !---Make file
  status=access(trim(filename)," ")
  if(status /= 0)then
     call make_ncfile_obs(filename)
  end if

  !---Count data
  so=inum+1
  no=0  
  do i=1,im,idx

     if(dat(i) == rmiss .or. dat(i) < min .or. max < dat(i)) cycle
     if(lon(i) < lons .or. lone < lon(i)) cycle
     if(lat(i) < lats .or. late < lat(i)) cycle

     inum=inum+1
     no=no+1

  end do

  !---Identify no observation 
  if(no < 1)then
     return
  end if

  !---Get data
  allocate(ele_tmp(no),ins_tmp(no))
  allocate(lon_tmp(no),lat_tmp(no),lev_tmp(no))
  allocate(dat_tmp(no),err_tmp(no))

  no=0
  ele_tmp(:)=ele
  ins_tmp(:)=ins
  lev_tmp(:)=lev
  err_tmp(:)=err
  do i=1,im,idx

     if(dat(i) == rmiss .or. dat(i) < min .or. max < dat(i)) cycle
     if(lon(i) < lons .or. lone < lon(i)) cycle
     if(lat(i) < lats .or. late < lat(i)) cycle

     no=no+1
     lon_tmp(no)=lon(i)
     lat_tmp(no)=lat(i)
     dat_tmp(no)=dat(i)

  end do

  !---Write data
  call write_obs_common(filename,so,no,ele_tmp,ins_tmp,lon_tmp,lat_tmp,lev_tmp,dat_tmp,err_tmp)

  !---Deallocate
  deallocate(ele_tmp,ins_tmp)
  deallocate(lon_tmp,lat_tmp,lev_tmp)
  deallocate(dat_tmp,err_tmp)

end subroutine write_obs_surface1d

!----------------------------------------------------------------------------
! Write Temperature/Salinity Observation |
!----------------------------------------------------------------------------

subroutine write_obs_ts(var,iyr,imon,iday,&
     & inum,km,lons,lone,lon,lats,late,lat,depth,dat,ins)

  use setting,only: &
       & t_min,t_max,t_err,t_ele, &
       & s_min,s_max,s_err,s_ele, &
       & ini_filename
  use mod_rmiss
  use mod_check_netcdf
  use netcdf
  implicit none

  !---Common
  integer k
  integer so,no
  integer status,access
  integer ele


  integer,allocatable :: ele_tmp(:),ins_tmp(:)

  real(kind = 8) min,max,err
  real(kind = 8) mdep(km)
  
  real(kind = 8),allocatable :: lon_tmp(:),lat_tmp(:),lev_tmp(:)
  real(kind = 8),allocatable :: dat_tmp(:),err_tmp(:)

  character(200) filename
  character(4) yyyy
  character(2) mm,dd

  !---IN
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: km
  integer,intent(in) :: ins

  real(kind = 8),intent(in) :: lons,lone,lon
  real(kind = 8),intent(in) :: lats,late,lat
  real(kind = 8),intent(in) :: depth(km)
  real(kind = 8),intent(in) :: dat(km)

  character(1),intent(in) :: var

  !---INOUT
  integer,intent(inout) :: inum

  !---Check position
  if(lon < lons .or. lone < lon) return
  if(lat < lats .or. late < lat) return

  !---Filename  
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday

  filename="../obs/"//trim(ini_filename)//yyyy//mm//dd//".nc"

  status=access(trim(filename)," ")
  if(status /= 0)then
     call make_ncfile_obs(filename)
  end if

  !---Setting
  if(var == "T")then
     call get_obs_interior_info("t",min,max,err,ele)
  elseif(var == "S")then
     call get_obs_interior_info("s",min,max,err,ele)
  elseif(var == "U")then
     call get_obs_interior_info("u",min,max,err,ele)
  elseif(var == "V")then
     call get_obs_interior_info("v",min,max,err,ele)
  else
     write(*,'(a)') "***Error: Choose T/S/U/V"
     stop
  end if

  !---Minus Depth
  mdep(:)=-1.d0*abs(depth(:))

  !---Count data
  so=inum+1
  no=0
  do k=1,km

     !Missing value
     if(dat(k) == rmiss .or. dat(k) /= dat(k))cycle
     !Out of range
     if(dat(k) < min .or. max < dat(k))cycle 

     inum=inum+1
     no=no+1

  end do

  !---Identify no observation 
  if(no < 1)then
     return
  end if

  !---Get data
  allocate(ele_tmp(no),ins_tmp(no))
  allocate(lon_tmp(no),lat_tmp(no),lev_tmp(no))
  allocate(dat_tmp(no),err_tmp(no))

  !---Get data
  no=0
  ele_tmp(:)=ele
  ins_tmp(:)=ins
  lon_tmp(:)=lon
  lat_tmp(:)=lat
  err_tmp(:)=err
  do k=1,km

     !Missing value
     if(dat(k) == rmiss .or. dat(k) /= dat(k))cycle
     !Out of range
     if(dat(k) < min .or. max < dat(k))cycle 

     no=no+1
     lev_tmp(no)=mdep(k)
     dat_tmp(no)=dat(k)

  end do

  !---Write data
  call write_obs_common(filename,so,no,ele_tmp,ins_tmp,lon_tmp,lat_tmp,lev_tmp,dat_tmp,err_tmp)

  !---Deallocate
  deallocate(ele_tmp,ins_tmp)
  deallocate(lon_tmp,lat_tmp,lev_tmp)
  deallocate(dat_tmp,err_tmp)

end subroutine write_obs_ts
