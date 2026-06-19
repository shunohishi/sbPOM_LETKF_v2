!----------------------------------------------------------------
! Check netcdf error |
!----------------------------------------------------------------

subroutine check_error(status)

  use netcdf
  use mpi
  implicit none

  !---Common
  integer ierr
  
  !---IN
  integer,intent(in) :: status

  if(status /= nf90_noerr) then
     write(*,'(a)') trim(nf90_strerror(status))
     call MPI_ABORT(MPI_COMM_WORLD,3,ierr)
  end if

end subroutine check_error

!----------------------------------------------------------------
! Read ensemble restart file |
!----------------------------------------------------------------

subroutine read_restart(dir,letkf,varname,iyr,imon,iday,im,jm,km,imem,dat)

  use netcdf
  use mpi
  implicit none

  !---Common
  integer status,access
  integer ncid,varid
  integer ierr

  character(200) filename
  character(8) yyyymmdd
  character(5) mmmmm
  character(4) yyyy
  character(2) mm,dd

  !---IN
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: im,jm,km
  integer,intent(in) :: imem

  character(100),intent(in) :: dir,letkf,varname

  !---OUT
  real(kind = 4),intent(out) :: dat(im,jm,km)

  write(mmmmm,'(i5.5)') imem
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday
  yyyymmdd=yyyy//mm//dd

  filename=trim(dir)//"/"//trim(letkf)//"/output/"//mmmmm//"/restart."//yyyymmdd//".nc"

  status=access(trim(filename)," ")
  if(status /= 0)then
     write(*,*) "***Error: Not found "//trim(filename)
     call MPI_ABORT(MPI_COMM_WORLD,2,ierr)
  end if

  status=nf90_open(trim(filename),nf90_nowrite,ncid)
  call check_error(status)

  status=nf90_inq_varid(ncid,trim(varname),varid)
  call check_error(status)
  if(km == 1)then
     status=nf90_get_var(ncid,varid,dat,(/1,1,1/),(/im,jm,1/))
  else
     status=nf90_get_var(ncid,varid,dat,(/1,1,1,1/),(/im,jm,km,1/))
  end if
  call check_error(status)
     
  status=nf90_close(ncid)
  call check_error(status)
  
end subroutine read_restart

!----------------------------------------------------------------
! Read ensemble restart file |
!----------------------------------------------------------------

subroutine write_restart(dir,letkf,varname,iyr,imon,iday,im,jm,km,dat)

  use netcdf
  use mpi
  implicit none

  !---Common
  integer status,access
  integer ncid,varid
  integer ierr

  character(200) filename
  character(8) yyyymmdd
  character(4) yyyy
  character(2) mm,dd

  !---IN
  integer,intent(in) :: iyr,imon,iday
  integer,intent(in) :: im,jm,km
  
  character(100),intent(in) :: dir,letkf,varname

  !---OUT
  real(kind = 4),intent(out) :: dat(im,jm,km)

  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday
  yyyymmdd=yyyy//mm//dd

  filename=trim(dir)//"/"//trim(letkf)//"/output/mean/restart."//yyyymmdd//".nc"
    
  status=access(trim(filename)," ")
  if(status /= 0)then
     write(*,*) "***Error: Not found "//trim(filename)
     call MPI_ABORT(MPI_COMM_WORLD,4,ierr)
  end if

  status=nf90_open(trim(filename),nf90_write,ncid)
  call check_error(status)

  status=nf90_inq_varid(ncid,trim(varname),varid)
  call check_error(status)
  if(km == 1)then
     status=nf90_put_var(ncid,varid,dat,(/1,1,1/),(/im,jm,1/))
  else
     status=nf90_put_var(ncid,varid,dat,(/1,1,1,1/),(/im,jm,km,1/))
  end if
  call check_error(status)
     
  status=nf90_close(ncid)
  call check_error(status)
  
end subroutine write_restart
