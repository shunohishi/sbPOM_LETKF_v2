!-------------------------------------------------------------------------------------
! Read ensemble mean restart file |
!-------------------------------------------------------------------------------------

subroutine read_restart(dir,letkf,varname,iyr,imon,iday,im,jm,km,dat)

  use netcdf
  implicit none

  !---Common
  integer status,access
  integer ncid,varid

  real(kind = 4) tmp(im,jm,km)
  
  character(200) filename
  character(8) yyyymmdd
  character(4) yyyy
  character(2) mm,dd

  !---IN
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: im,jm,km

  character(100),intent(in) :: dir,letkf,varname

  !---OUT
  real(kind = 8),intent(out) :: dat(im,jm,km)

  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday
  yyyymmdd=yyyy//mm//dd

  filename=trim(dir)//"/"//trim(letkf)//"/output/mean/restart."//yyyymmdd//".nc"

  status=access(trim(filename)," ")
  if(status /= 0)then
     write(*,*) "***Error: Not found "//trim(filename)
     stop
  end if

  status=nf90_open(trim(filename),nf90_nowrite,ncid)
  call check_error(status)

  status=nf90_inq_varid(ncid,trim(varname),varid)
  call check_error(status)
  if(km == 1)then
     status=nf90_get_var(ncid,varid,tmp,(/1,1,1/),(/im,jm,1/))
  else
     status=nf90_get_var(ncid,varid,tmp,(/1,1,1,1/),(/im,jm,km,1/))
  end if
  call check_error(status)
     
  status=nf90_close(ncid)
  call check_error(status)

  dat(:,:,:)=dble(tmp(:,:,:))
  
end subroutine read_restart

!-------------------------------------------------------------------------------------
! Write system grid |
!-------------------------------------------------------------------------------------

subroutine write_grid_ncfile(im,jm,km, &
     & lont,lonu,lonv, &
     & latt,latu,latv, &
     & dept,depu,depv,depw, &
     & maskt,masku,maskv)

  use setting, filename => filename_grid
  use netcdf
  implicit none

  !---Common
  integer status,access
  integer ncid,varid

  !---IN
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: lont(im),lonu(im),lonv(im)
  real(kind = 8),intent(in) :: latt(jm),latu(jm),latv(jm)
  real(kind = 8),intent(in) :: dept(im,jm,km),depu(im,jm,km),depv(im,jm,km),depw(im,jm,km)
  real(kind = 8),intent(in) :: maskt(im,jm),masku(im,jm),maskv(im,jm)

  !---Check file
  status=access(trim(filename)," ")
  if(status /= 0)then
     write(*,*) "***Error: Not found "//trim(filename)
     stop
  end if

  !---Open
  status=nf90_open(trim(filename),nf90_write,ncid)
  call check_error(status)

  !---Write data
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
  status=nf90_inq_varid(ncid,"dept",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,real(dept))
  call check_error(status)

  status=nf90_inq_varid(ncid,"depu",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,real(depu))
  call check_error(status)

  status=nf90_inq_varid(ncid,"depv",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,real(depv))
  call check_error(status)

  status=nf90_inq_varid(ncid,"depw",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,real(depw))
  call check_error(status)

  !Mask
  status=nf90_inq_varid(ncid,"maskt",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,real(maskt))
  call check_error(status)

  status=nf90_inq_varid(ncid,"masku",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,real(masku))
  call check_error(status)

  status=nf90_inq_varid(ncid,"maskv",varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,real(maskv))
  call check_error(status)

  !---Close
  status=nf90_close(ncid)
  call check_error(status)  

end subroutine write_grid_ncfile

!-------------------------------------------------------------------------------------
! Write 2D data |
!-------------------------------------------------------------------------------------

subroutine write_dat2d_ncfile(ms,varname,iyr,imon,iday,im,jm,dat)

  use netcdf
  implicit none

  !---Common
  integer status,access
  integer ncid,varid

  character(100) dirname,filename
  character(8) yyyymmdd
  character(6) yyyymm
  character(4) yyyy
  character(2) mm,dd

  !---IN
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: im,jm

  real(kind = 8),intent(in) :: dat(im,jm)

  character(4),intent(in) :: ms
  character(100),intent(in) :: varname

  !---Filename
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  if(iday == 0)then
     yyyymm=yyyy//mm
     dirname="dat/"//yyyymm//"/"//trim(ms)
     filename=trim(dirname)//"/"//trim(varname)//yyyymm//".nc"
  else
     write(dd,'(i2.2)') iday
     yyyymmdd=yyyy//mm//dd
     dirname="dat/"//yyyymmdd//"/"//trim(ms)
     filename=trim(dirname)//"/"//trim(varname)//yyyymmdd//".nc"
  end if

  !---Check file
  status=access(trim(filename)," ")
  if(status /= 0)then
     write(*,*) "***Error: Not found "//trim(filename)
     stop
  end if

  !---Open
  status=nf90_open(trim(filename),nf90_write,ncid)
  call check_error(status)

  !---Write
  status=nf90_inq_varid(ncid,trim(varname),varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,real(dat))
  call check_error(status)

  !---Close
  status=nf90_close(ncid)  
  call check_error(status)

end subroutine write_dat2d_ncfile

!-------------------------------------------------------------------------------------
! Write 3D data |
!-------------------------------------------------------------------------------------

subroutine write_dat3d_ncfile(ms,varname,iyr,imon,iday,im,jm,km,dat)

  use netcdf
  implicit none

  !---Common
  integer status,access
  integer ncid,varid

  character(100) dirname,filename
  character(8) yyyymmdd
  character(6) yyyymm
  character(4) yyyy
  character(2) mm,dd

  !---IN
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: dat(im,jm,km)

  character(4),intent(in) :: ms
  character(100),intent(in) :: varname

  !---Filename
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  if(iday == 0)then
     yyyymm=yyyy//mm
     dirname="dat/"//yyyymm//"/"//trim(ms)
     filename=trim(dirname)//"/"//trim(varname)//yyyymm//".nc"
  else  
     write(dd,'(i2.2)') iday
     yyyymmdd=yyyy//mm//dd
     dirname="dat/"//yyyymmdd//"/"//trim(ms)
     filename=trim(dirname)//"/"//trim(varname)//yyyymmdd//".nc"
  end if
  
  !---Check file
  status=access(trim(filename)," ")
  if(status /= 0)then
     write(*,*) "***Error: Not found "//trim(filename)
     stop
  end if

  !---Open
  status=nf90_open(trim(filename),nf90_write,ncid)
  call check_error(status)

  !---Write
  status=nf90_inq_varid(ncid,trim(varname),varid)
  call check_error(status)
  status=nf90_put_var(ncid,varid,real(dat))
  call check_error(status)

  !---Close
  status=nf90_close(ncid)  
  call check_error(status)

end subroutine write_dat3d_ncfile
