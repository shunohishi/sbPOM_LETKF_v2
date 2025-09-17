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

  subroutine make_obsfile(filename)

    use netcdf
    implicit none

    !---Parameter
    integer,parameter :: ndim=1

    !---Common
    integer status
    integer ncid,dimid,varid
    integer dim(ndim)

    !---IN
    character(100),intent(in) :: filename

    status=nf90_create(trim(filename),nf90_netcdf4,ncid)
    call check_error(status)

    status=nf90_put_att(ncid,NF90_GLOBAL,"title","drifter buoy observation")
    call check_error(status)
    status=nf90_put_att(ncid,NF90_GLOBAL,"description","drifter buoy observation vs. analysis in obs. space")
    call check_error(status)

    !#Obs.
    status=nf90_def_dim(ncid,"nobs",nf90_unlimited,dimid)
    call check_error(status)
    dim(1)=dimid

    !1D
    call define_var_netcdf(ncid,1,dim(1),varid,"dble", &
         & "lon_o","longitude","degree E")

    call define_var_netcdf(ncid,1,dim(1),varid,"dble", &
         & "lat_o","latitude","degree N")  

    call define_var_netcdf(ncid,1,dim(1),varid,"dble", &
         & "ht_a","analysis SST in obs. space","degree C")

    call define_var_netcdf(ncid,1,dim(1),varid,"dble", &
         & "hu_a","analysis SSU in obs. space","m/s")

    call define_var_netcdf(ncid,1,dim(1),varid,"dble", &
         & "hv_a","analysis SSV in obs. space","m/s")

    call define_var_netcdf(ncid,1,dim(1),varid,"dble", &
         & "htsprd_a","analysis SST spread in obs. space","degree C")

    call define_var_netcdf(ncid,1,dim(1),varid,"dble", &
         & "husprd_a","analysis SSU spread in obs. space","m/s")

    call define_var_netcdf(ncid,1,dim(1),varid,"dble", &
         & "hvsprd_a","analysis SSV spread in obs. space","m/s")

    call define_var_netcdf(ncid,1,dim(1),varid,"dble", &
         & "t_o","observed SST","degree C")

    call define_var_netcdf(ncid,1,dim(1),varid,"dble", &
         & "u_o","observed SSU","m/s")

    call define_var_netcdf(ncid,1,dim(1),varid,"dble", &
         & "v_o","observed SSV","m/s")

    status=nf90_enddef(ncid)
    call check_error(status)

    status=nf90_close(ncid)
    call check_error(status)

  end subroutine make_obsfile

end module mod_make_ncfile
