!-------------------------------------------------------------------------
! Define variable |
!-------------------------------------------------------------------------

subroutine define_var_netcdf(ncid,ndim,dim,varid,type,name,long_name,units_name)

  use netcdf
  implicit none

  integer status

  integer,intent(in) :: ncid
  integer,intent(in) :: ndim,dim(ndim)
  integer,intent(inout) :: varid

  character(*),intent(in) :: type
  character(*),intent(in) :: name,long_name,units_name

  if(trim(type) == "int")then
     status=nf90_def_var(ncid,trim(name),nf90_int,dim,varid)
  else if(trim(type) == "real")then
     status=nf90_def_var(ncid,trim(name),nf90_float,dim,varid)
  else if(trim(type) == "dble")then
     status=nf90_def_var(ncid,trim(name),nf90_double,dim,varid)
  end if
  call check_error(status)

  status=nf90_def_var_deflate(ncid,varid,shuffle=1,deflate=1,deflate_level=5)
  call check_error(status)
  
  status=nf90_put_att(ncid,varid,"long_name",trim(long_name))
  call check_error(status)
  
  status=nf90_put_att(ncid,varid,"units",trim(units_name))
  call check_error(status)

end subroutine define_var_netcdf

!----------------------------------------------------------------
! Check netcdf error |
!----------------------------------------------------------------

subroutine check_error(status)

  use netcdf
  implicit none

  integer,intent(in) :: status
 
  if(status /= nf90_noerr) then
    write(*,'(a)') trim(nf90_strerror(status))
    stop "Stopped"
  end if

end subroutine check_error

!----------------------------------------------------------------
! Make netcdf file |
!----------------------------------------------------------------

subroutine make_ncfile_web(im,jm,km,nt,filename)

  use setting, only: title
  use mod_rmiss
  use netcdf
  implicit none

  !---Parameter
  integer,parameter :: ndim_2d=3
  integer,parameter :: ndim_3d=4
  real(kind = 8),parameter :: fillvalue=rmiss

  !---Common
  integer status
  integer ncid,dimid,varid
  integer dim_2d(ndim_2d)
  integer dim_3d(ndim_3d)

  character(4) yyyy
  character(2) mm
  
  !---IN
  integer,intent(in) :: im,jm,km,nt
  character(100),intent(in) :: filename

  !---NF90_CREATE
  status=nf90_create(trim(filename),nf90_netcdf4,ncid)
  call check_error(status)

  !TITLE
  status=nf90_put_att(ncid,NF90_GLOBAL,"title",trim(title))
  call check_error(status)

  !REERENCE
  status=nf90_put_att(ncid,NF90_GLOBAL,"references","Ohishi et al. (submitted)")
  call check_error(status)

  !---DIMID
  !X
  status=nf90_def_dim(ncid,"im",im,dimid)
  call check_error(status)
  dim_2d(1)=dimid
  dim_3d(1)=dimid
  
  !Y
  status=nf90_def_dim(ncid,"jm",jm,dimid)
  call check_error(status)
  dim_2d(2)=dimid
  dim_3d(2)=dimid
  
  !Z
  status=nf90_def_dim(ncid,"km",km,dimid)
  call check_error(status)
  dim_3d(3)=dimid
  
  !T
  status=nf90_def_dim(ncid,"nt",nt,dimid)
  call check_error(status)
  dim_2d(3)=dimid
  dim_3d(4)=dimid

  !---1D
  !time
  call define_var_netcdf(ncid,1,dim_3d(4),varid,"real", &
       & "time","Time (day)","day since 1950-1-1 00:00:00")

  !lon
  call define_var_netcdf(ncid,1,dim_3d(1),varid,"real", &
       & "lont","Longitude (T, S, and SSH)","degree E")

  call define_var_netcdf(ncid,1,dim_3d(1),varid,"real", &
       & "lonu","Longitude (U)","degree E")

  call define_var_netcdf(ncid,1,dim_3d(1),varid,"real", &
       & "lonv","Longitude (V)","degree E")
  
  !lat
  call define_var_netcdf(ncid,1,dim_3d(2),varid,"real", &
       & "latt","Latitude (T, S, and SSH)","degree N")

  call define_var_netcdf(ncid,1,dim_3d(2),varid,"real", &
       & "latu","Latitude (U)","degree N")

  call define_var_netcdf(ncid,1,dim_3d(2),varid,"real", &
       & "latv","Latitude (V)","degree N")
  
  !depth
  call define_var_netcdf(ncid,1,dim_3d(3),varid,"real", &
       & "depth","Depth","meter")

  !---2D
  !hmean
  call define_var_netcdf(ncid,3,dim_2d(1:3),varid,"real", &
       & "elmean","Ensemble mean of sea surface height","meter")
  status=nf90_def_var_fill(ncid,varid,0,real(fillvalue))
  call check_error(status)

  !hsprd
  call define_var_netcdf(ncid,3,dim_2d(1:3),varid,"real", &
       & "elsprd","Ensemble spread of sea surface height","meter")
  status=nf90_def_var_fill(ncid,varid,0,real(fillvalue))
  call check_error(status)
  
  !---3D
  !tmean
  call define_var_netcdf(ncid,4,dim_3d(1:4),varid,"real", &
       & "tmean","Ensemble mean of temperature","degree Celcius")
  status=nf90_def_var_fill(ncid,varid,0,real(fillvalue))
  call check_error(status)

  !tsprd
  call define_var_netcdf(ncid,4,dim_3d(1:4),varid,"real", &
       & "tsprd","Ensemble spread of temperature","degree Celcius")
  status=nf90_def_var_fill(ncid,varid,0,real(fillvalue))
  call check_error(status)

  !smean
  call define_var_netcdf(ncid,4,dim_3d(1:4),varid,"real", &
       & "smean","Ensemble mean of salinity","-")
  status=nf90_def_var_fill(ncid,varid,0,real(fillvalue))
  call check_error(status)

  !ssprd
  call define_var_netcdf(ncid,4,dim_3d(1:4),varid,"real", &
       & "ssprd","Ensemble spread of salinity","-")
  status=nf90_def_var_fill(ncid,varid,0,real(fillvalue))
  call check_error(status)

  !umean
  call define_var_netcdf(ncid,4,dim_3d(1:4),varid,"real", &
       & "umean","Ensemble mean of zonal velocity","m/s")
  status=nf90_def_var_fill(ncid,varid,0,real(fillvalue))
  call check_error(status)

  !usprd
  call define_var_netcdf(ncid,4,dim_3d(1:4),varid,"real", &
       & "usprd","Ensemble spread of zonal velocity","m/s")
  status=nf90_def_var_fill(ncid,varid,0,real(fillvalue))
  call check_error(status)

  !vmean
  call define_var_netcdf(ncid,4,dim_3d(1:4),varid,"real", &
       & "vmean","Ensemble mean of meridional velocity","m/s")
  status=nf90_def_var_fill(ncid,varid,0,real(fillvalue))
  call check_error(status)

  !vsprd
  call define_var_netcdf(ncid,4,dim_3d(1:4),varid,"real", &
       & "vsprd","Ensemble spread of meridional velocity","m/s")
  status=nf90_def_var_fill(ncid,varid,0,real(fillvalue))
  call check_error(status)

  !NF90_ENDDEF/CLOSE
  status=nf90_enddef(ncid)
  call check_error(status)

  status=nf90_close(ncid)
  call check_error(status)  
  
end subroutine make_ncfile_web

!----------------------------------------------------------------
! Write Grid |
!----------------------------------------------------------------

subroutine write_grid(iyr,imon,im,jm,km,lont,lonu,lonv,latt,latu,latv,depth)

  use netcdf
  implicit none

  !---Common
  integer status,access
  integer ncid,varid
  integer nt
  
  !---IN
  integer,intent(in) :: iyr,imon

  integer,intent(in) :: im,jm,km
  real(kind = 8),intent(in) :: lont(im),lonu(im),lonv(im)
  real(kind = 8),intent(in) :: latt(jm),latu(jm),latv(jm)
  real(kind = 8),intent(in) :: depth(km)

  character(100) filename
  character(5) var
  character(4) yyyy
  character(2) mm

  !Filename
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon

  filename="dat/"//yyyy//mm//".nc"

  !nt
  if(imon == 2 .and. mod(iyr,4) == 0)then
     nt=29
  else if(imon == 2)then
     nt=28
  else if(imon == 6 .or. imon == 9 .or. imon == 11)then
     nt=30
  else
     nt=31
  end if     

  !Check file
  status=access(trim(filename)," ")
  if(status == 0)then
     return
  else
     call make_ncfile_web(im,jm,km,nt,filename)
  end if

  !Open
  status=nf90_open(trim(filename),nf90_write,ncid)
  call check_error(status)

  !Longitude
  status=nf90_inq_varid(ncid,"lont",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,real(lont))
  call check_error(status)

  status=nf90_inq_varid(ncid,"lonu",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,real(lonu))
  call check_error(status)

  status=nf90_inq_varid(ncid,"lonv",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,real(lonv))
  call check_error(status)

  !Latitude
  status=nf90_inq_varid(ncid,"latt",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,real(latt))
  call check_error(status)

  status=nf90_inq_varid(ncid,"latu",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,real(latu))
  call check_error(status)

  status=nf90_inq_varid(ncid,"latv",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,real(latv))
  call check_error(status)
  
  !Depth
  status=nf90_inq_varid(ncid,"depth",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,real(depth))
  call check_error(status)

  !Close
  status=nf90_close(ncid)
  call check_error(status)
  
end subroutine write_grid

!----------------------------------------------------------------------
! Write data |
!----------------------------------------------------------------------

subroutine write_data(varname,iyr,imon,iday,im,jm,km,rjul,dat)

  use netcdf
  implicit none

  !---Common
  integer status,access
  integer ncid,varid

  character(100) filename
  character(4) yyyy
  character(2) mm
  
  !---IN
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: rjul(1)
  real(kind = 8),intent(in) :: dat(im,jm,km)
  
  character(*),intent(in) :: varname

  !Filename
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon

  filename="dat/"//yyyy//mm//".nc"

  status=access(trim(filename)," ")
  if(status /= 0)then
     write(*,*) "***Error: Not found "//trim(filename)
     stop
  end if
  
  status=nf90_open(trim(filename),nf90_write,ncid)
  call check_error(status)

  status=nf90_inq_varid(ncid,"time",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,real(rjul),(/iday/),(/1/))
  call check_error(status)
  
  status=nf90_inq_varid(ncid,trim(varname),varid)
  call check_error(status)
  if(km == 1)then
     status=nf90_put_var(ncid,varid,real(dat),(/1,1,iday/),(/im,jm,1/))
  else
     status=nf90_put_var(ncid,varid,real(dat),(/1,1,1,iday/),(/im,jm,km,1/))
  end if
  call check_error(status)
  
  status=nf90_close(ncid)
  call check_error(status)
  
end subroutine write_data
