module mod_make_ncfile

contains

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
       stop "***Error: Stopped"
    end if

  end subroutine check_error

  !----------------------------------------------------------------
  ! Make obsfile |
  !----------------------------------------------------------------

  subroutine make_ncfile(im,jm,km,nt,filename)

    use netcdf
    implicit none

    !---Parameter
    integer,parameter :: ndim=4
    
    !---Common
    integer status
    integer ncid,dimid,varid

    integer dim(ndim)

    !---IN
    integer,intent(in) :: im,jm,km,nt
    
    character(100),intent(in) :: filename

    status=nf90_create(trim(filename),nf90_netcdf4,ncid)
    call check_error(status)

    status=nf90_put_att(ncid,NF90_GLOBAL,"title","Climatology")
    call check_error(status)
    status=nf90_put_att(ncid,NF90_GLOBAL,"description","Climatology")
    call check_error(status)

    !---1D
    status=nf90_def_dim(ncid,"x",im,dimid)
    call check_error(status)
    dim(1)=dimid

    call define_var_netcdf(ncid,1,dim(1),varid,"dble", &
         & "lont","longitude [tsh]","degree E")

    call define_var_netcdf(ncid,1,dim(1),varid,"dble", &
         & "lonu","longitude [u]","degree E")

    call define_var_netcdf(ncid,1,dim(1),varid,"dble", &
         & "lonv","longitude [v]","degree E")
    
    status=nf90_def_dim(ncid,"y",jm,dimid)
    call check_error(status)
    dim(2)=dimid

    call define_var_netcdf(ncid,1,dim(2),varid,"dble", &
         & "latt","latitude [tsh]","degree N")  

    call define_var_netcdf(ncid,1,dim(2),varid,"dble", &
         & "latu","latitude [u]","degree N")  

    call define_var_netcdf(ncid,1,dim(2),varid,"dble", &
         & "latv","latitude [v]","degree N")  

    !---2D
    call define_var_netcdf(ncid,2,dim(1:2),varid,"dble", &
         & "maskt","land-sea mask [tsh]","-")

    call define_var_netcdf(ncid,2,dim(1:2),varid,"dble", &
         & "masku","land-sea mask [u]","-")

    call define_var_netcdf(ncid,2,dim(1:2),varid,"dble", &
         & "maskv","land-sea mask [v]","-")    
    
    !---Time
    status=nf90_def_dim(ncid,"t",nt,dimid)
    call check_error(status)
    dim(3)=dimid

    call define_var_netcdf(ncid,1,dim(3),varid,"dble", &
         & "time","time","month")

    !---2D+Time        
    call define_var_netcdf(ncid,3,dim(1:3),varid,"dble", &
         & "hmean","analysis SSH","m")

    call define_var_netcdf(ncid,3,dim(1:3),varid,"dble", &
         & "hsprd","analysis SSH spread","m")
    
    !---3D
    status=nf90_def_dim(ncid,"z",km,dimid)
    call check_error(status)
    dim(4)=dim(3)
    dim(3)=dimid

    call define_var_netcdf(ncid,3,dim(1:3),varid,"dble", &
         & "dept","depth [tsh]","m")    

    call define_var_netcdf(ncid,3,dim(1:3),varid,"dble", &
         & "depu","depth [u]","m")    

    call define_var_netcdf(ncid,3,dim(1:3),varid,"dble", &
         & "depv","depth [v]","m")    
    
    !---3D+Time    
    call define_var_netcdf(ncid,ndim,dim,varid,"dble", &
         & "tmean","analysis temperature","degree C")

    call define_var_netcdf(ncid,ndim,dim,varid,"dble", &
         & "tsprd","analysis temperature spread","degree C")

    call define_var_netcdf(ncid,ndim,dim,varid,"dble", &
         & "smean","analysis salinity","-")

    call define_var_netcdf(ncid,ndim,dim,varid,"dble", &
         & "ssprd","analysis salinity spread","-")

    call define_var_netcdf(ncid,ndim,dim,varid,"dble", &
         & "umean","analysis zonal velocity","m/s")

    call define_var_netcdf(ncid,ndim,dim,varid,"dble", &
         & "usprd","analysis zonal velocity spread","m/s")

    call define_var_netcdf(ncid,ndim,dim,varid,"dble", &
         & "vmean","analysis meridional velocity","m/s")

    call define_var_netcdf(ncid,ndim,dim,varid,"dble", &
         & "vsprd","analysis meridional velocity spread","m/s")       
    
    status=nf90_enddef(ncid)
    call check_error(status)

    status=nf90_close(ncid)
    call check_error(status)

  end subroutine make_ncfile

end module mod_make_ncfile
