subroutine read_wind_sd(imon,im,jm,usd,vsd)

  implicit none
  include 'netcdf.inc'

  integer imon
  integer im,jm
  integer status,access
  integer varid,ncid

  real usd(im,jm),vsd(im,jm)
  
  character(100) filename

  filename="../in/wind.nc"
  
  status=access(trim(filename)," ")
  if(status /= 0)then
     write(*,*) "***Error: Not Found "//trim(filename)
     stop
  end if

  status=nf_open(trim(filename),nf_nowrite,ncid)

  status=nf_inq_varid(ncid,"windu_sd",varid)
  status=nf_get_vara_real(ncid,varid,(/1,1,imon/),(/im,jm,1/),usd)

  status=nf_inq_varid(ncid,"windv_sd",varid)
  status=nf_get_vara_real(ncid,varid,(/1,1,imon/),(/im,jm,1/),vsd)

  status=nf_close(ncid)

end subroutine read_wind_sd
