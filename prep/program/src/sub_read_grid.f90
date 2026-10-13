subroutine read_grid(lonu,lonv,lone,lonc,latu,latv,late,latc,depthz,depthzz,fsm,dum,dvm)

  !$ use omp_lib
  use mod_gridinfo
  use mod_check_netcdf
  use netcdf
  implicit none

  integer i,j,k
  integer status,access
  integer ncid,varid

  real(kind = 8) z(im,jm,km),zz(im,jm,km)
  real(kind = 8) east_u(im,jm),east_v(im,jm),east_e(im,jm),east_c(im,jm)
  real(kind = 8) north_u(im,jm),north_v(im,jm),north_e(im,jm),north_c(im,jm)
  real(kind = 8) h(im,jm)
  
  real(kind = 8),intent(out) :: lonu(im),lonv(im),lone(im),lonc(im)
  real(kind = 8),intent(out) :: latu(jm),latv(jm),late(jm),latc(jm)
  real(kind = 8),intent(out) :: depthz(im,jm,km),depthzz(im,jm,km)
  real(kind = 8),intent(out) :: fsm(im,jm),dum(im,jm),dvm(im,jm)

  character(100) filename

  filename="../in/grid.nc"
  
  status=access(trim(filename)," ")

  if(status /= 0)then
     write(*,*) "***Error: Not found "//trim(filename)
     stop
  end if

  status=nf90_open(trim(filename),nf90_nowrite,ncid)
  call check_netcdf(status)
  
  status=nf90_inq_varid(ncid,"z_w",varid)
  call check_netcdf(status)
  status=nf90_get_var(ncid,varid,z)
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"z_e",varid)
  call check_netcdf(status)
  status=nf90_get_var(ncid,varid,zz)
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"east_u",varid)
  call check_netcdf(status)
  status=nf90_get_var(ncid,varid,east_u)
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"east_v",varid)
  call check_netcdf(status)
  status=nf90_get_var(ncid,varid,east_v)
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"east_e",varid)
  call check_netcdf(status)
  status=nf90_get_var(ncid,varid,east_e)
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"east_c",varid)
  call check_netcdf(status)
  status=nf90_get_var(ncid,varid,east_c)
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"north_u",varid)
  call check_netcdf(status)
  status=nf90_get_var(ncid,varid,north_u)
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"north_v",varid)
  call check_netcdf(status)
  status=nf90_get_var(ncid,varid,north_v)
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"north_e",varid)
  call check_netcdf(status)
  status=nf90_get_var(ncid,varid,north_e)
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"north_c",varid)
  call check_netcdf(status)
  status=nf90_get_var(ncid,varid,north_c)
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"h",varid)
  call check_netcdf(status)
  status=nf90_get_var(ncid,varid,h)
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"fsm",varid)
  call check_netcdf(status)
  status=nf90_get_var(ncid,varid,fsm)
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"dum",varid)
  call check_netcdf(status)
  status=nf90_get_var(ncid,varid,dum)
  call check_netcdf(status)

  status=nf90_inq_varid(ncid,"dvm",varid)
  call check_netcdf(status)
  status=nf90_get_var(ncid,varid,dvm)
  call check_netcdf(status)

  status=nf90_close(ncid)
  call check_netcdf(status)

  lonu(:)=east_u(:,1)
  lonv(:)=east_v(:,1)
  lone(:)=east_e(:,1)
  lonc(:)=east_c(:,1)

  latu(:)=north_u(1,:)
  latv(:)=north_v(1,:)
  late(:)=north_e(1,:)
  latc(:)=north_c(1,:)

  !$omp parallel
  !$omp do private(i,j,k)
  do k=1,km
     do j=1,jm
        do i=1,im
           depthz(i,j,k)=fsm(i,j)*h(i,j)*z(i,j,k)
           depthz(i,j,k)=abs(depthz(i,j,k))
           depthzz(i,j,k)=fsm(i,j)*h(i,j)*zz(i,j,k)
           depthzz(i,j,k)=abs(depthzz(i,j,k))
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine read_grid
