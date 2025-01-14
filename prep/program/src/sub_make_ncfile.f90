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
     write(*,*) trim(nf90_strerror(status))
     stop "Stopped"
  end if

end subroutine check_error

!-----------------------------------------------------------------------
! GRID |
!-----------------------------------------------------------------------

subroutine make_ncfile_grid(im,jm,km,filename)

  use netcdf
  implicit none

  integer,parameter :: ndim=3

  integer status,system
  integer ncid,dimid,varid
  integer dim(ndim)

  integer,intent(in) :: im,jm,km
  character(100),intent(in) :: filename

  status=system("rm -f "//trim(filename))

  status=nf90_create(trim(filename),nf90_netcdf4,ncid)
  call check_error(status)

  status=nf90_put_att(ncid,NF90_GLOBAL,"title","grid")
  call check_error(status)
  status=nf90_put_att(ncid,NF90_GLOBAL,"description","sbPOM grid file")
  call check_error(status)

  !X
  status=nf90_def_dim(ncid,"x",im,dimid)
  call check_error(status)
  dim(1)=dimid

  !Y
  status=nf90_def_dim(ncid,"y",jm,dimid)
  call check_error(status)
  dim(2)=dimid

  !Z
  status=nf90_def_dim(ncid,"z",km,dimid)
  call check_error(status)
  dim(3)=dimid

  !2D
  call define_var_netcdf(ncid,2,dim(1:2),varid,"dble","dx","grid increment in x","meter")
  call define_var_netcdf(ncid,2,dim(1:2),varid,"dble","dy","grid increment in y","meter")

  call define_var_netcdf(ncid,2,dim(1:2),varid,"dble","east_u","east of u-points","degree")
  call define_var_netcdf(ncid,2,dim(1:2),varid,"dble","east_v","east of v-points","degree")
  call define_var_netcdf(ncid,2,dim(1:2),varid,"dble","east_e","east of elevation points","degree")
  call define_var_netcdf(ncid,2,dim(1:2),varid,"dble","east_c","east of cell corners","degree")

  call define_var_netcdf(ncid,2,dim(1:2),varid,"dble","north_u","north of u-points","degree")
  call define_var_netcdf(ncid,2,dim(1:2),varid,"dble","north_v","north of v-points","degree")
  call define_var_netcdf(ncid,2,dim(1:2),varid,"dble","north_e","north of elevation points","degree")
  call define_var_netcdf(ncid,2,dim(1:2),varid,"dble","north_c","north of cell corners","degree")

  call define_var_netcdf(ncid,2,dim(1:2),varid,"dble","h","undistrubed water depth","meter")
  status=nf90_put_att(ncid,varid,"coords","east_e north_e")
  call check_error(status)

  call define_var_netcdf(ncid,2,dim(1:2),varid,"dble","hdif", &
       & "Difference between original and smooth water depth","meter")
  status=nf90_put_att(ncid,varid,"coords","east_e north_e")
  call check_error(status)  

  call define_var_netcdf(ncid,2,dim(1:2),varid,"dble","ratio", &
       & "abs[h(i+1)-h(i)]/abs[h(i+1)+h(i)]","-")
  status=nf90_put_att(ncid,varid,"coords","east_e north_e")
  call check_error(status)  

  call define_var_netcdf(ncid,2,dim(1:2),varid,"dble","fsm","free surface mask","-")
  status=nf90_put_att(ncid,varid,"coords","east_e north_e")
  call check_error(status)

  call define_var_netcdf(ncid,2,dim(1:2),varid,"dble","dum","u-velocity mask","-")
  status=nf90_put_att(ncid,varid,"coords","east_u north_u")
  call check_error(status)

  call define_var_netcdf(ncid,2,dim(1:2),varid,"dble","dvm","v-velocity mask","-")
  status=nf90_put_att(ncid,varid,"coords","east_v north_v")
  call check_error(status)

  !3D
  call define_var_netcdf(ncid,3,dim(1:3),varid,"dble","z_w","sigma of w-points","sigma level")
  call define_var_netcdf(ncid,3,dim(1:3),varid,"dble","z_e","sigma of tsuv-points","sigma level")
  
  status=nf90_enddef(ncid)
  call check_error(status)
  status=nf90_close(ncid)
  call check_error(status)

end subroutine make_ncfile_grid

!-------------------------------------------------------------------------
! TSCLIM
!-------------------------------------------------------------------------

subroutine make_ncfile_tsclim(im,jm,km,filename)

  use netcdf
  implicit none

  integer,parameter :: ndim=4

  integer status
  integer ncid,dimid,varid
  integer dim(ndim)

  integer,intent(in) :: im,jm,km
  character(100),intent(in) :: filename

  status=nf90_create(trim(filename),nf90_netcdf4,ncid)
  call check_error(status)

  status=nf90_put_att(ncid,NF90_GLOBAL,"title","tsclim")
  call check_error(status)
  status=nf90_put_att(ncid,NF90_GLOBAL,"description","sbPOM climatology from WOA18")
  call check_error(status)

  !X
  status=nf90_def_dim(ncid,"x",im,dimid)
  call check_error(status)
  dim(1)=dimid

  !Y
  status=nf90_def_dim(ncid,"y",jm,dimid)
  call check_error(status)
  dim(2)=dimid

  !Z
  status=nf90_def_dim(ncid,"z",km,dimid)
  call check_error(status)
  dim(3)=dimid

  !TIME
  status=nf90_def_dim(ncid,"t",12,dimid)
  call check_error(status)
  dim(4)=dimid

  !3D
  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "tclim","Annual potential temperature climatology","degree C")
  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "sclim","Annual salinity climatology","-")
  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "rmean","Annual rmean climatology: rmean=(density-1000)/rhoref","-")

  !4D
  call define_var_netcdf(ncid,4,dim(1:4),varid,"real", &
       & "tclimm","Monthly potential temperature climatology","degree C")
  call define_var_netcdf(ncid,4,dim(1:4),varid,"real", &
       & "sclimm","Monthly salinity climatology","-")
  call define_var_netcdf(ncid,4,dim(1:4),varid,"real", &
       & "rmeanm","Monthly rmean climatology: rmean=(density-1000)/rhoref","-")

  status=nf90_enddef(ncid)
  call check_error(status)
  status=nf90_close(ncid)
  call check_error(status)

end subroutine make_ncfile_tsclim

!-------------------------------------------------------------------------
! Initial Condition |
!-------------------------------------------------------------------------

subroutine make_ncfile_ic(im,jm,km,filename)

  use netcdf
  implicit none

  integer,parameter :: ndim=3

  integer status
  integer ncid,dimid,varid
  integer dim(ndim)

  integer,intent(in) :: im,jm,km
  character(100),intent(in) :: filename

  status=nf90_create(trim(filename),nf90_netcdf4,ncid)
  call check_error(status)

  status=nf90_put_att(ncid,NF90_GLOBAL,"title","ic")
  call check_error(status)
  status=nf90_put_att(ncid,NF90_GLOBAL,"description","sbPOM initial condition from WOA18")
  call check_error(status)

  !X
  status=nf90_def_dim(ncid,"x",im,dimid)
  call check_error(status)
  dim(1)=dimid

  !Y
  status=nf90_def_dim(ncid,"y",jm,dimid)
  call check_error(status)
  dim(2)=dimid

  !Z
  status=nf90_def_dim(ncid,"z",km,dimid)
  call check_error(status)
  dim(3)=dimid

  !3D
  call define_var_netcdf(ncid,3,dim(1:3),varid,"dble", &
       & "t","Potential temperature","degree C")
  call define_var_netcdf(ncid,3,dim(1:3),varid,"dble", &
       & "s","Salinity","-")
  call define_var_netcdf(ncid,3,dim(1:3),varid,"dble", &
       & "u","Zonal velocity","meter/sec")
  call define_var_netcdf(ncid,3,dim(1:3),varid,"dble", &
       & "v","Meridional velocity","meter/sec")

  status=nf90_enddef(ncid)
  call check_error(status)
  status=nf90_close(ncid)
  call check_error(status)

end subroutine make_ncfile_ic

!----------------------------------------------------------------
! LBC |
!----------------------------------------------------------------

subroutine make_ncfile_lbc_mclim(im,jm,km,filename)

  use netcdf
  implicit none

  integer status
  integer ncid,varid
  integer dimx,dimy,dimz,dimt
  integer,allocatable :: dim(:)

  integer,intent(in) :: im,jm,km
  character(100),intent(in) :: filename

  status=nf90_create(trim(filename),nf90_netcdf4,ncid)
  call check_error(status)

  status=nf90_put_att(ncid,NF90_GLOBAL,"title","lbc")
  call check_error(status)
  status=nf90_put_att(ncid,NF90_GLOBAL,"description","sbPOM boudary condition from SODA")
  call check_error(status)

  !X
  status=nf90_def_dim(ncid,"x",im,dimx)
  call check_error(status)

  !Y
  status=nf90_def_dim(ncid,"y",jm,dimy)
  call check_error(status)

  !Z
  status=nf90_def_dim(ncid,"z",km,dimz)
  call check_error(status)

  !Time (month)
  status=nf90_def_dim(ncid,"t",12,dimt)
  call check_error(status)

  !1D-x
  allocate(dim(2))
  dim(1)=dimx
  dim(2)=dimt

  call define_var_netcdf(ncid,2,dim(1:2),varid,"real", &
       & "eln","Sea level on northern boundary","m")

  call define_var_netcdf(ncid,2,dim(1:2),varid,"real", &
       & "els","Sea level on southern boundary","m")

  call define_var_netcdf(ncid,2,dim(1:2),varid,"real", &
       & "uabn","Vertical mean of u on northern boundary","m/s")

  call define_var_netcdf(ncid,2,dim(1:2),varid,"real", &
       & "uabs","Vertical mean of u on southern boundary","m/s")

  call define_var_netcdf(ncid,2,dim(1:2),varid,"real", &
       & "vabn","Vertical mean of v on northern boundary","m/s")

  call define_var_netcdf(ncid,2,dim(1:2),varid,"real", &
       & "vabs","Vertical mean of v on southern boundary","m/s")

  deallocate(dim)

  !1D-y
  allocate(dim(2))
  dim(1)=dimy
  dim(2)=dimt

  call define_var_netcdf(ncid,2,dim(1:2),varid,"real", &
       & "ele","Sea level on eastern boundary","m")

  call define_var_netcdf(ncid,2,dim(1:2),varid,"real", &
       & "elw","Sea level on western boundary","m")

  call define_var_netcdf(ncid,2,dim(1:2),varid,"real", &
       & "uabe","Vertical mean of u on eastern boundary","m/s")

  call define_var_netcdf(ncid,2,dim(1:2),varid,"real", &
       & "uabw","Vertical mean of u on western boundary","m/s")

  call define_var_netcdf(ncid,2,dim(1:2),varid,"real", &
       & "vabe","Vertical mean of v on eastern boundary","m/s")

  call define_var_netcdf(ncid,2,dim(1:2),varid,"real", &
       & "vabw","Vertical mean of v on western boundary","m/s")

  deallocate(dim)

  !2D-xz
  allocate(dim(3))
  dim(1)=dimx
  dim(2)=dimz
  dim(3)=dimt

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "tbn","Potential temperature on northern boundary","degree C")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "tbs","Potential temperature on southern boundary","degree C")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "sbn","Salinity on northern boundary","-")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "sbs","Salinity on southern boundary","-")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "ubn","u on northern boundary","m/s")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "ubs","u on southern boundary","m/s")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "vbn","v on northern boundary","m/s")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "vbs","v on southern boundary","m/s")

  deallocate(dim)

  !2D-yz
  allocate(dim(3))
  dim(1)=dimy
  dim(2)=dimz
  dim(3)=dimt

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "tbe","Potential temperature on eastern boundary","degree C")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "tbw","Potential temperature on western boundary","degree C")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "sbe","Salinity on eastern boundary","-")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "sbw","Salinity on western boundary","-")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "ube","u on eastern boundary","m/s")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "ubw","u on western boundary","m/s")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "vbe","v on eastern boundary","m/s")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "vbw","v on western boundary","m/s")

  deallocate(dim)

  status=nf90_enddef(ncid)
  call check_error(status)
  status=nf90_close(ncid)
  call check_error(status)

end subroutine make_ncfile_lbc_mclim

!---------------------------------------------------------------------
! TSDATA |
!---------------------------------------------------------------------

subroutine make_ncfile_tsdata_mclim(im,jm,km,filename)

  use netcdf
  implicit none

  integer,parameter :: ndim2d=3
  integer,parameter :: ndim3d=4

  integer status
  integer ncid,varid
  integer xid,yid,zid,tid
  integer,allocatable :: dim(:)

  integer,intent(in) :: im,jm,km
  character(100),intent(in) :: filename

  status=nf90_create(trim(filename),nf90_netcdf4,ncid)
  call check_error(status)

  status=nf90_put_att(ncid,NF90_GLOBAL,"title","tsdata mclim")
  call check_error(status)
  status=nf90_put_att(ncid,NF90_GLOBAL,"description","sbPOM T & S monthly climatology from SODA")
  call check_error(status)

  !X
  status=nf90_def_dim(ncid,"x",im,xid)
  call check_error(status)

  !Y
  status=nf90_def_dim(ncid,"y",jm,yid)
  call check_error(status)

  !Z
  status=nf90_def_dim(ncid,"z",km,zid)
  call check_error(status)

  !T
  status=nf90_def_dim(ncid,"t",12,tid)
  call check_error(status)
  
  !2D
  allocate(dim(ndim2d))
  dim(1)=xid
  dim(2)=yid
  dim(3)=tid
  call define_var_netcdf(ncid,ndim2d,dim(1:ndim2d),varid,"real", &
       & "ssh","Sea surface height","meter")
  call define_var_netcdf(ncid,ndim2d,dim(1:ndim2d),varid,"real", &
       & "ua","Vertical averaged u","m/s")
  call define_var_netcdf(ncid,ndim2d,dim(1:ndim2d),varid,"real", &
       & "va","Vertical averaged v","m/s")
  deallocate(dim)
  
  !3D
  allocate(dim(ndim3d))
  dim(1)=xid
  dim(2)=yid
  dim(3)=zid
  dim(4)=tid
  call define_var_netcdf(ncid,ndim3d,dim(1:ndim3d),varid,"real", &
       & "temp","Potential temperature","degree C")
  call define_var_netcdf(ncid,ndim3d,dim(1:ndim3d),varid,"real", &
       & "sal","Salinity","-")
  call define_var_netcdf(ncid,ndim3d,dim(1:ndim3d),varid,"real", &
       & "u","Zonal velocity","m/s")
  call define_var_netcdf(ncid,ndim3d,dim(1:ndim3d),varid,"real", &
       & "v","Meridional velocity","m/s")
  deallocate(dim)
  
  status=nf90_enddef(ncid)
  call check_error(status)
  status=nf90_close(ncid)
  call check_error(status)

end subroutine make_ncfile_tsdata_mclim

!---------------------------------------------------------------------
! ATM |
!---------------------------------------------------------------------

subroutine make_ncfile_atm(im,jm,iyr,imon,iday,filename)

  use netcdf
  implicit none

  integer,parameter :: ndim=3

  integer status
  integer ncid,dimid,varid
  integer dim(ndim)

  character(4) yyyy
  character(2) mm,dd

  !IN
  integer,intent(in) :: im,jm
  integer,intent(in) :: iyr,imon,iday
  character(100),intent(in) :: filename

  status=nf90_create(trim(filename),nf90_netcdf4,ncid)
  call check_error(status)

  status=nf90_put_att(ncid,NF90_GLOBAL,"title","atm")
  call check_error(status)
  status=nf90_put_att(ncid,NF90_GLOBAL,"description","sbPOM atmospheric forcing data")
  call check_error(status)

  !X
  status=nf90_def_dim(ncid,"x",im,dimid)
  call check_error(status)
  dim(1)=dimid

  !Y
  status=nf90_def_dim(ncid,"y",jm,dimid)
  call check_error(status)
  dim(2)=dimid

  !T
  status=nf90_def_dim(ncid,"t",nf90_unlimited,dimid)
  call check_error(status)
  dim(3)=dimid

  !1D
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday
  call define_var_netcdf(ncid,1,dim(3),varid,"real", &
       & "atmtime","Atmospheric data time","days since "//yyyy//"-"//mm//"-"//dd//" 00:00:00")

  !3D
  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "windu","Near-surface zonal wind velocity","m/s")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "windv","Near-surface meridional wind velocity","m/s")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "airt","Near-surface air temperature","degree C")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "airh","Near-surface air humidity","-")


  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "swrad","Surface downward shortwave radiation","W/m^2")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "lwrad","Surface downward longwave radiation","W/m^2")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "slp","Sea level pressure","Pa")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "prep","Precipitation","mm/day")
  
  status=nf90_enddef(ncid)
  call check_error(status)
  status=nf90_close(ncid)
  call check_error(status)

end subroutine make_ncfile_atm

!---------------------------------------------------------------------
! FFlux |
!---------------------------------------------------------------------

subroutine make_ncfile_fflux(im,jm,iyr,imon,iday,filename)

  use netcdf
  implicit none

  integer,parameter :: ndim=3

  integer status
  integer ncid,dimid,varid
  integer dim(ndim)

  character(4) yyyy
  character(2) mm,dd

  !IN
  integer,intent(in) :: im,jm
  integer,intent(in) :: iyr,imon,iday
  character(100),intent(in) :: filename

  status=nf90_create(trim(filename),nf90_netcdf4,ncid)
  call check_error(status)

  status=nf90_put_att(ncid,NF90_GLOBAL,"title","river")
  call check_error(status)
  status=nf90_put_att(ncid,NF90_GLOBAL,"description","sbPOM River data")
  call check_error(status)

  !X
  status=nf90_def_dim(ncid,"x",im,dimid)
  call check_error(status)
  dim(1)=dimid

  !Y
  status=nf90_def_dim(ncid,"y",jm,dimid)
  call check_error(status)
  dim(2)=dimid

  !T
  status=nf90_def_dim(ncid,"t",nf90_unlimited,dimid)
  call check_error(status)
  dim(3)=dimid

  !1D
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday
  call define_var_netcdf(ncid,1,dim(3),varid,"real", &
       & "rivtime","River data time","days since "//yyyy//"-"//mm//"-"//dd//" 00:00:00")

  !3D
  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "river","River discharge","mm/day")

  status=nf90_enddef(ncid)
  call check_error(status)
  status=nf90_close(ncid)
  call check_error(status)

end subroutine make_ncfile_fflux

!---------------------------------------------------------------------
! Assimilated SST |
!---------------------------------------------------------------------

subroutine make_ncfile_sst(im,jm,iyr,imon,filename)

  use netcdf
  implicit none

  integer,parameter :: ndim=3

  integer status
  integer ncid,dimid,varid
  integer dim(ndim)

  character(4) yyyy
  character(2) mm

  !IN
  integer,intent(in) :: im,jm
  integer,intent(in) :: iyr,imon
  character(100),intent(in) :: filename

  status=nf90_create(trim(filename),nf90_netcdf4,ncid)
  call check_error(status)

  status=nf90_put_att(ncid,NF90_GLOBAL,"title","sst")
  call check_error(status)
  status=nf90_put_att(ncid,NF90_GLOBAL,"description","Assimilated satellite SST")
  call check_error(status)

  !X
  status=nf90_def_dim(ncid,"x",im,dimid)
  call check_error(status)
  dim(1)=dimid

  !Y
  status=nf90_def_dim(ncid,"y",jm,dimid)
  call check_error(status)
  dim(2)=dimid

  !T
  status=nf90_def_dim(ncid,"t",nf90_unlimited,dimid)
  call check_error(status)
  dim(3)=dimid

  !1D
  call define_var_netcdf(ncid,1,dim(1),varid,"real", &
       & "lon","Longitude","degree E")

  call define_var_netcdf(ncid,1,dim(2),varid,"real", &
       & "lat","Latitude","degree N")

  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  call define_var_netcdf(ncid,1,dim(3),varid,"real", &
       & "time","Composite SST data time","days since "//yyyy//"-"//mm//"-01 00:00:00")

  !3D
  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "asst","Aqua/AMSR-E SST","degree C")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "csst","Coriolis/WindSat SST","degree C")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "gsst","GCOM-W/AMSR2 SST","degree C")  
  
  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "hsst","Himawari/AHI SST","degree C")
  
  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "sst","Assimilated SST","degree C")

  status=nf90_enddef(ncid)
  call check_error(status)

  status=nf90_close(ncid)
  call check_error(status)

end subroutine make_ncfile_sst

!---------------------------------------------------------------------
! Assimilated SSS |
!---------------------------------------------------------------------

subroutine make_ncfile_sss(im,jm,iyr,imon,filename)

  use netcdf
  implicit none

  integer,parameter :: ndim=3

  integer status
  integer ncid,dimid,varid
  integer dim(ndim)

  character(4) yyyy
  character(2) mm

  !IN
  integer,intent(in) :: im,jm
  integer,intent(in) :: iyr,imon
  character(100),intent(in) :: filename

  status=nf90_create(trim(filename),nf90_netcdf4,ncid)
  call check_error(status)

  status=nf90_put_att(ncid,NF90_GLOBAL,"title","sss")
  call check_error(status)
  status=nf90_put_att(ncid,NF90_GLOBAL,"description","Assimilated satellite SSS")
  call check_error(status)

  !X
  status=nf90_def_dim(ncid,"x",im,dimid)
  call check_error(status)
  dim(1)=dimid

  !Y
  status=nf90_def_dim(ncid,"y",jm,dimid)
  call check_error(status)
  dim(2)=dimid

  !T
  status=nf90_def_dim(ncid,"t",nf90_unlimited,dimid)
  call check_error(status)
  dim(3)=dimid

  !1D
  call define_var_netcdf(ncid,1,dim(1),varid,"real", &
       & "lon","Longitude","degree E")

  call define_var_netcdf(ncid,1,dim(2),varid,"real", &
       & "lat","Latitude","degree N")

  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  call define_var_netcdf(ncid,1,dim(3),varid,"real", &
       & "time","Composite SSS data time","days since "//yyyy//"-"//mm//"-01 00:00:00")

  !3D
  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "sss","Assimilated SSS (SMOS +SMAP)","-")

  status=nf90_enddef(ncid)
  call check_error(status)

  status=nf90_close(ncid)
  call check_error(status)

end subroutine make_ncfile_sss

!---------------------------------------------------------------------
! Assimilated SSH |
!---------------------------------------------------------------------

subroutine make_ncfile_ssh(im,jm,iyr,imon,filename)

  use netcdf
  implicit none

  integer,parameter :: ndim=3

  integer status
  integer ncid,dimid,varid
  integer dim(ndim)

  character(4) yyyy
  character(2) mm

  !IN
  integer,intent(in) :: im,jm
  integer,intent(in) :: iyr,imon
  character(100),intent(in) :: filename

  status=nf90_create(trim(filename),nf90_netcdf4,ncid)
  call check_error(status)

  status=nf90_put_att(ncid,NF90_GLOBAL,"title","ssh")
  call check_error(status)
  status=nf90_put_att(ncid,NF90_GLOBAL,"description","Assimilated satellite SSH")
  call check_error(status)

  !X
  status=nf90_def_dim(ncid,"x",im,dimid)
  call check_error(status)
  dim(1)=dimid

  !Y
  status=nf90_def_dim(ncid,"y",jm,dimid)
  call check_error(status)
  dim(2)=dimid

  !T
  status=nf90_def_dim(ncid,"t",nf90_unlimited,dimid)
  call check_error(status)
  dim(3)=dimid

  !1D
  call define_var_netcdf(ncid,1,dim(1),varid,"real", &
       & "lon","Longitude","degree E")

  call define_var_netcdf(ncid,1,dim(2),varid,"real", &
       & "lat","Latitude","degree N")

  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  call define_var_netcdf(ncid,1,dim(3),varid,"real", &
       & "time","Composite SSH data time","days since "//yyyy//"-"//mm//"-01 00:00:00")

  !3D
  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "mdot","Mean Dynamic Ocean Topography","meter")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "ssh","Sea surface height","meter")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "ssha","Sea surface height anomaly","meter")

  status=nf90_enddef(ncid)
  call check_error(status)

  status=nf90_close(ncid)
  call check_error(status)

end subroutine make_ncfile_ssh

!---------------------------------------------------------------------
! Assimilated SSH |
!---------------------------------------------------------------------

subroutine make_ncfile_ssuv(im,jm,iyr,imon,filename)

  use netcdf
  implicit none

  integer,parameter :: ndim=3

  integer status
  integer ncid,dimid,varid
  integer dim(ndim)

  character(4) yyyy
  character(2) mm

  !IN
  integer,intent(in) :: im,jm
  integer,intent(in) :: iyr,imon
  character(100),intent(in) :: filename

  status=nf90_create(trim(filename),nf90_netcdf4,ncid)
  call check_error(status)

  status=nf90_put_att(ncid,NF90_GLOBAL,"title","ssuv")
  call check_error(status)
  status=nf90_put_att(ncid,NF90_GLOBAL,"description","Assimilated satellite sea surface velocity")
  call check_error(status)

  !X
  status=nf90_def_dim(ncid,"x",im,dimid)
  call check_error(status)
  dim(1)=dimid

  !Y
  status=nf90_def_dim(ncid,"y",jm,dimid)
  call check_error(status)
  dim(2)=dimid

  !T
  status=nf90_def_dim(ncid,"t",nf90_unlimited,dimid)
  call check_error(status)
  dim(3)=dimid

  !1D
  call define_var_netcdf(ncid,1,dim(1),varid,"real", &
       & "lon","Longitude","degree E")

  call define_var_netcdf(ncid,1,dim(2),varid,"real", &
       & "lat","Latitude","degree N")

  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  call define_var_netcdf(ncid,1,dim(3),varid,"real", &
       & "time","Composite sea surface velocity data time", &
       & "days since "//yyyy//"-"//mm//"-01 00:00:00")

  !3D
  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "ssu","Sea surface zonal velocity","m/s")

  call define_var_netcdf(ncid,3,dim(1:3),varid,"real", &
       & "ssv","Sea surface meridional velocity","m/s")

  status=nf90_enddef(ncid)
  call check_error(status)

  status=nf90_close(ncid)
  call check_error(status)

end subroutine make_ncfile_ssuv

!---------------------------------------------------------------------
! Assimilated TS |
!---------------------------------------------------------------------

subroutine make_ncfile_ts(filename)

  use netcdf
  implicit none

  integer,parameter :: ndim=1
  
  integer status
  integer ncid,dimid,varid
  integer dim(ndim)

  !IN
  character(100),intent(in) :: filename
  
  status=nf90_create(trim(filename),nf90_netcdf4,ncid)
  call check_error(status)

  status=nf90_put_att(ncid,NF90_GLOBAL,"title","TS profile information")
  call check_error(status)
  status=nf90_put_att(ncid,NF90_GLOBAL,"description","Assimilated TS profile information")
  call check_error(status)

  !#Obs.
  status=nf90_def_dim(ncid,"nobs",nf90_unlimited,dimid)
  call check_error(status)
  dim(1)=dimid

  !1D
  call define_var_netcdf(ncid,1,dim(1),varid,"real", &
       & "lon","longitude","degree E")
  
  call define_var_netcdf(ncid,1,dim(1),varid,"real", &
       & "lat","latitude","degree N")

  call define_var_netcdf(ncid,1,dim(1),varid,"real", &
       & "lev1","most upper level","meter")

  call define_var_netcdf(ncid,1,dim(1),varid,"real", &
       & "levb","most bottom level","meter")

  call define_var_netcdf(ncid,1,dim(1),varid,"int", &
       & "dat","Obs. dataset(1: GTSPP, 2: AQC Argo","")
    
  status=nf90_enddef(ncid)
  call check_error(status)
  
  status=nf90_close(ncid)
  call check_error(status)
  
end subroutine make_ncfile_ts

!---------------------------------------------------------------------
! Assimilated Obs. |
!---------------------------------------------------------------------

subroutine make_ncfile_obs(filename)

  use netcdf
  implicit none

  integer,parameter :: ndim=1

  integer status
  integer ncid,dimid,varid
  integer dim(ndim)

  !IN
  character(100),intent(in) :: filename

  status=nf90_create(trim(filename),nf90_netcdf4,ncid)
  call check_error(status)

  status=nf90_put_att(ncid,NF90_GLOBAL,"title","observation")
  call check_error(status)
  status=nf90_put_att(ncid,NF90_GLOBAL,"description","Assimilated observation")
  call check_error(status)

  !#Obs.
  status=nf90_def_dim(ncid,"nobs",nf90_unlimited,dimid)
  call check_error(status)
  dim(1)=dimid

  !1D
  call define_var_netcdf(ncid,1,dim(1),varid,"real", &
       & "ele","element","element (h:2567, u:2819, v:2820, t:3073, s:3332)")

  call define_var_netcdf(ncid,1,dim(1),varid,"real", &
       & "lon","longitude","degree E")

  call define_var_netcdf(ncid,1,dim(1),varid,"real", &
       & "lat","latitude","degree N")

  call define_var_netcdf(ncid,1,dim(1),varid,"real", &
       & "lev","level","meter")

  call define_var_netcdf(ncid,1,dim(1),varid,"real", &
       & "dat","observation","H: meter, U and V: m/s, T: degree C, S: -")

  call define_var_netcdf(ncid,1,dim(1),varid,"real", &
       & "err","observation error","H: meter, U and V: m/s, T: degree C, S: -")

  status=nf90_enddef(ncid)
  call check_error(status)

  status=nf90_close(ncid)
  call check_error(status)

end subroutine make_ncfile_obs
