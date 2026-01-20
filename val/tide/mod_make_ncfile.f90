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
       stop "Stopped"
    end if

  end subroutine check_error

  !----------------------------------------------------------------
  ! Make obsfile |
  !----------------------------------------------------------------

  subroutine make_obsfile(nst,filename)

    use netcdf
    implicit none

    !---Parameter
    integer,parameter :: ndim=2

    !---Common
    integer status
    integer ncid,dimid,varid
    integer dim(ndim)

    !---IN
    integer,intent(in) :: nst
    character(100),intent(in) :: filename

    status=nf90_create(trim(filename),nf90_netcdf4,ncid)
    call check_error(status)

    status=nf90_put_att(ncid,NF90_GLOBAL,"title","tide gauge observation")
    call check_error(status)
    status=nf90_put_att(ncid,NF90_GLOBAL,"description","tide gauge observation vs. analysis in obs. space")
    call check_error(status)

    !Dimension
    status=nf90_def_dim(ncid,"nst",nst,dimid)
    call check_error(status)
    dim(1)=dimid

    status=nf90_def_dim(ncid,"time",nf90_unlimited,dimid)
    call check_error(status)
    dim(2)=dimid
    
    !1D
    call define_var_netcdf(ncid,2,dim,varid,"dble", &
         & "lon_a","longitude at the nearest model grid","degree E")

    call define_var_netcdf(ncid,2,dim,varid,"dble", &
         & "lat_a","latitude at the nearest model grid","degree N")  

    call define_var_netcdf(ncid,2,dim,varid,"dble", &
         & "hmean_a","analysis SSH ensemble mean at the nearest model grid","meter")

    call define_var_netcdf(ncid,2,dim,varid,"dble", &
         & "hsprd_a","analysis SSH ensemble spread at the nearest model grid","meter")    
    
    call define_var_netcdf(ncid,2,dim,varid,"dble", &
         & "lon_o","longitude in observation","degree E")

    call define_var_netcdf(ncid,2,dim,varid,"dble", &
         & "lat_o","latitude in observation","degree N")  

    call define_var_netcdf(ncid,2,dim,varid,"dble", &
         & "h_o","analysis SSH in observation","meter")

    call define_var_netcdf(ncid,2,dim,varid,"dble", &
         & "dist","distance between the nearest model grid and obs. point","meter")
    
    status=nf90_enddef(ncid)
    call check_error(status)

    status=nf90_close(ncid)
    call check_error(status)

  end subroutine make_obsfile

end module mod_make_ncfile
