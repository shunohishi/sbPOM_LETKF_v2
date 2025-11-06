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
  character(3),intent(in) :: varname

  !---OUT
  integer,intent(out) :: ele(1)
  real(kind = 8),intent(out) :: min,max,err(1)

  !---Min,Max,Err,Ele
  if(varname == "ssu")then
     min=ssu_min
     max=ssu_max
     err(1)=ssu_err
     ele(1)=ssu_ele     
  else if(varname == "ssv")then
     min=ssv_min
     max=ssv_max
     err(1)=ssv_err
     ele(1)=ssv_ele     
  else if(varname == "sst")then
     min=sst_min
     max=sst_max
     err(1)=sst_err
     ele(1)=sst_ele     
  else if(varname == "sss")then
     min=sss_min
     max=sss_max
     err(1)=sss_err
     ele(1)=sss_ele     
  else if(varname == "ssh")then
     min=ssh_min
     max=ssh_max
     err(1)=ssh_err
     ele(1)=ssh_ele
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
  character(3),intent(in) :: varname

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
! surface2d: lon(im,jm),lat(im,jm)
! surface2dg: lon(im),lat(jm)
! surface1d: lon(ntime),lat(ntime)
!-------------------------------------------------------

subroutine write_obs_surface2d(varname,ins,iyr,imon,iday,&
     & idx,idy,im,jm,lons,lone,lon,lats,late,lat,dat,inum)

  use setting, only: ini_filename
  use mod_rmiss
  use netcdf
  implicit none

  !---Parameter
  real(kind = 8),parameter :: lev(1)=0.d0
  
  !---Common
  integer i,j
  integer ele(1)
  integer status,access
  integer ncid,varid
  
  real(kind = 8) min,max,err(1)

  character(200) filename
  character(8) yyyymmdd
  character(4) yyyy
  character(2) mm,dd
  
  !---IN
  integer,intent(in) :: ins(1)
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

  !---Write data
  status=nf90_open(trim(filename),nf90_write,ncid)

  do j=1,jm,idy
     do i=1,im,idx

        if(dat(i,j) == rmiss .or. dat(i,j) < min .or. max < dat(i,j)) cycle
        if(lon(i,j) < lons .or. lone < lon(i,j)) cycle
        if(lat(i,j) < lats .or. late < lat(i,j)) cycle

        inum=inum+1

        status=nf90_inq_varid(ncid,"ele",varid)
        status=nf90_put_var(ncid,varid,ele(1:1),(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"ins",varid)
        status=nf90_put_var(ncid,varid,ins(1:1),(/inum/),(/1/))        
        
        status=nf90_inq_varid(ncid,"lon",varid)
        status=nf90_put_var(ncid,varid,[real(lon(i,j))],(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"lat",varid)
        status=nf90_put_var(ncid,varid,[real(lat(i,j))],(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"lev",varid)
        status=nf90_put_var(ncid,varid,[real(lev(1:1))],(/inum/),(/1/))
        
        status=nf90_inq_varid(ncid,"dat",varid)
        status=nf90_put_var(ncid,varid,[real(dat(i,j))],(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"err",varid)
        status=nf90_put_var(ncid,varid,[real(err(1:1))],(/inum/),(/1/))

     end do
  end do

  status=nf90_close(ncid)
  
end subroutine write_obs_surface2d

!-----------------------------------

subroutine write_obs_surface2dg(varname,ins,iyr,imon,iday, &
     & idx,idy,im,jm,lons,lone,lon,lats,late,lat,dat,inum)

  use setting, only: ini_filename
  use mod_rmiss
  use netcdf
  implicit none

  !---Parameter
  real(kind = 8),parameter :: lev(1)=0.d0
  
  !---Common
  integer i,j
  integer ele(1)
  integer status,access
  integer ncid,varid
  
  real(kind = 8) min,max,err(1)
  
  character(200) filename  
  character(8) yyyymmdd
  character(4) yyyy
  character(2) mm,dd
  
  !---IN
  integer,intent(in) :: ins(1)
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

  !---Write data
  status=nf90_open(trim(filename),nf90_write,ncid)

  do j=1,jm,idy

     if(lat(j) < lats .or. late < lat(j)) cycle        

     do i=1,im,idx

        if(dat(i,j) == rmiss .or. dat(i,j) < min .or. max < dat(i,j)) cycle
        if(lon(i) < lons .or. lone < lon(i)) cycle

        inum=inum+1

        status=nf90_inq_varid(ncid,"ele",varid)
        status=nf90_put_var(ncid,varid,ele(1:1),(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"ins",varid)
        status=nf90_put_var(ncid,varid,ins(1:1),(/inum/),(/1/))        
        
        status=nf90_inq_varid(ncid,"lon",varid)
        status=nf90_put_var(ncid,varid,[real(lon(i))],(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"lat",varid)
        status=nf90_put_var(ncid,varid,[real(lat(j))],(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"lev",varid)
        status=nf90_put_var(ncid,varid,[real(lev(1:1))],(/inum/),(/1/))
        
        status=nf90_inq_varid(ncid,"dat",varid)
        status=nf90_put_var(ncid,varid,[real(dat(i,j))],(/inum/),(/1/))

        status=nf90_inq_varid(ncid,"err",varid)
        status=nf90_put_var(ncid,varid,[real(err(1:1))],(/inum/),(/1/))

     end do
  end do

  status=nf90_close(ncid)
  
end subroutine write_obs_surface2dg

!--------------------------------------

subroutine write_obs_surface1d(varname,ins,iyr,imon,iday,&
     & idx,im,lons,lone,lon,lats,late,lat,dat,inum)

  use setting, only: ini_filename
  use mod_rmiss
  use netcdf
  implicit none

  !---Parameter
  real(kind = 8),parameter :: lev(1)=0.d0
  
  !---Common
  integer i
  integer ele(1)
  integer status,access
  integer ncid,varid

  real(kind = 8) min,max,err(1)
  
  character(200) filename
  character(8) yyyymmdd
  character(4) yyyy
  character(2) mm,dd
  
  !---IN
  integer,intent(in) :: ins(1)
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

  !---Write data
  status=nf90_open(trim(filename),nf90_write,ncid)

  do i=1,im,idx

     if(dat(i) == rmiss .or. dat(i) < min .or. max < dat(i)) cycle
     if(lon(i) < lons .or. lone < lon(i)) cycle
     if(lat(i) < lats .or. late < lat(i)) cycle
     
     inum=inum+1

     status=nf90_inq_varid(ncid,"ele",varid)
     status=nf90_put_var(ncid,varid,ele(1:1),(/inum/),(/1/))

     status=nf90_inq_varid(ncid,"ins",varid)
     status=nf90_put_var(ncid,varid,ins(1:1),(/inum/),(/1/))        

     status=nf90_inq_varid(ncid,"lon",varid)
     status=nf90_put_var(ncid,varid,[real(lon(i:i))],(/inum/),(/1/))

     status=nf90_inq_varid(ncid,"lat",varid)
     status=nf90_put_var(ncid,varid,[real(lat(i:i))],(/inum/),(/1/))

     status=nf90_inq_varid(ncid,"lev",varid)
     status=nf90_put_var(ncid,varid,[real(lev(1:1))],(/inum/),(/1/))

     status=nf90_inq_varid(ncid,"dat",varid)
     status=nf90_put_var(ncid,varid,[real(dat(i:i))],(/inum/),(/1/))

     status=nf90_inq_varid(ncid,"err",varid)
     status=nf90_put_var(ncid,varid,[real(err(1:1))],(/inum/),(/1/))

  end do

  status=nf90_close(ncid)
  
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
  use netcdf
  implicit none

  !---Common
  integer k
  integer status,access
  integer ncid,varid

  integer ele(1)
  
  real(kind = 8) min,max,err(km)
  real(kind = 8) mdep(km)

  character(100) filename
  character(4) yyyy
  character(2) mm,dd
  
  !---IN
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: km
  integer,intent(in) :: ins(1)

  real(kind = 8),intent(in) :: lons,lone,lon(1)
  real(kind = 8),intent(in) :: lats,late,lat(1)
  real(kind = 8),intent(in) :: depth(km)
  real(kind = 8),intent(in) :: dat(km)

  character(1),intent(in) :: var
  
  !---INOUT
  integer,intent(inout) :: inum

  !---Check position
  if(lon(1) < lons .or. lone < lon(1)) return
  if(lat(1) < lats .or. late < lat(1)) return
  
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
     min=t_min
     max=t_max
     ele(:)=t_ele
     err(:)=t_err
  elseif(var == "S")then
     min=s_min
     max=s_max
     ele(:)=s_ele
     err(:)=s_err
  else
     write(*,'(a)') "***Error: Choose T or S"
     stop
  end if

  !---Minus Depth
  mdep(:)=-1.d0*abs(depth(:))
  
  !---Write Data
  status=nf90_open(trim(filename),nf90_write,ncid)
  do k=1,km

     !Missing value
     if(dat(k) == rmiss .or. dat(k) /= dat(k))cycle
     !Out of range
     if(dat(k) < min .or. max < dat(k))cycle 

     !Count up
     inum=inum+1

     status=nf90_inq_varid(ncid,"ele",varid)
     status=nf90_put_var(ncid,varid,ele,(/inum/),(/1/))

     status=nf90_inq_varid(ncid,"ins",varid)
     status=nf90_put_var(ncid,varid,ins,(/inum/),(/1/))
     
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
