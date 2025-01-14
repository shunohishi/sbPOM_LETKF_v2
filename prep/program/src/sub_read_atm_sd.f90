subroutine read_atm_sd(var,imon,ihour,dt,im,jm,dat)

  implicit none
  include 'netcdf.inc'

  integer imon,ihour,dt
  integer itime
  integer im,jm
  integer status,access
  integer varid,ncid

  real dat(im,jm)

  character(2) month
  character(3) var
  character(10) ncvar
  character(100) filename

  write(month,'(i2.2)') imon
  filename="../in/atm_mclim."//month//".nc "

  status=access(trim(filename)," ")
  if(status /= 0)then
     write(*,*) "***Error: Not Found "//trim(filename)
     stop
  end if

  if(trim(var) == "u")then
     ncvar="windu"
  else if(trim(var) == "v")then
     ncvar="windv"
  else if(trim(var) == "ta")then
     ncvar="airt"
  else if(trim(var) == "qa")then
     ncvar="airh"
  else if(trim(var) == "sw")then
     ncvar="swrad"
  else if(trim(var) == "tc")then
     ncvar="cloud"
  else if(trim(var) == "slp")then
     ncvar="slp"
  else
     write(*,*) "***Error not found:"//trim(var)
     stop
  end if

  status=nf_open(trim(filename),nf_nowrite,ncid)

  status=nf_inq_varid(ncid,trim(ncvar)//"_sd",varid)
  itime=ihour/dt
  status=nf_get_vara_real(ncid,varid,(/1,1,itime+1/),(/im,jm,1/),dat)

  status=nf_close(ncid)

end subroutine read_atm_sd
