subroutine read_grid(lonu,lonv,lone,lonc,latu,latv,late,latc,depthz,depthzz,fsm,dum,dvm)

  !$use omp_lib
  use mod_gridinfo
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

  status=nf90_inq_varid(ncid,"z_w",varid)
  status=nf90_get_var(ncid,varid,z)

  status=nf90_inq_varid(ncid,"z_e",varid)
  status=nf90_get_var(ncid,varid,zz)

  status=nf90_inq_varid(ncid,"east_u",varid)
  status=nf90_get_var(ncid,varid,east_u)

  status=nf90_inq_varid(ncid,"east_v",varid)
  status=nf90_get_var(ncid,varid,east_v)

  status=nf90_inq_varid(ncid,"east_e",varid)
  status=nf90_get_var(ncid,varid,east_e)

  status=nf90_inq_varid(ncid,"east_c",varid)
  status=nf90_get_var(ncid,varid,east_c)

  status=nf90_inq_varid(ncid,"north_u",varid)
  status=nf90_get_var(ncid,varid,north_u)

  status=nf90_inq_varid(ncid,"north_v",varid)
  status=nf90_get_var(ncid,varid,north_v)

  status=nf90_inq_varid(ncid,"north_e",varid)
  status=nf90_get_var(ncid,varid,north_e)

  status=nf90_inq_varid(ncid,"north_c",varid)
  status=nf90_get_var(ncid,varid,north_c)

  status=nf90_inq_varid(ncid,"h",varid)
  status=nf90_get_var(ncid,varid,h)

  status=nf90_inq_varid(ncid,"fsm",varid)
  status=nf90_get_var(ncid,varid,fsm)

  status=nf90_inq_varid(ncid,"dum",varid)
  status=nf90_get_var(ncid,varid,dum)

  status=nf90_inq_varid(ncid,"dvm",varid)
  status=nf90_get_var(ncid,varid,dvm)

  status=nf90_close(ncid)

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
