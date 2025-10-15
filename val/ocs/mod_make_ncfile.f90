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

  subroutine make_ncfile(km,varname,filename)

    use netcdf
    implicit none

    !---Parameter
    integer,parameter :: im=1,jm=1
    integer,parameter :: ndim=4
    
    !---Common
    integer status
    integer ncid,dimid,varid

    integer dim(ndim)

    !---IN
    integer,intent(in) :: km
    
    character(1),intent(in) :: varname
    character(100),intent(in) :: filename

    status=nf90_create(trim(filename),nf90_netcdf4,ncid)
    call check_error(status)

    status=nf90_put_att(ncid,NF90_GLOBAL,"title","ocean climate station (OCS)")
    call check_error(status)
    status=nf90_put_att(ncid,NF90_GLOBAL,"description","OCS observation vs. analysis in obs. space")
    call check_error(status)

    status=nf90_def_dim(ncid,"x",im,dimid)
    call check_error(status)
    dim(1)=dimid

    status=nf90_def_dim(ncid,"y",jm,dimid)
    call check_error(status)
    dim(2)=dimid
    
    status=nf90_def_dim(ncid,"z",km,dimid)
    call check_error(status)
    dim(3)=dimid
    
    status=nf90_def_dim(ncid,"time",nf90_unlimited,dimid)
    call check_error(status)
    dim(4)=dimid

    call define_var_netcdf(ncid,1,dim(1),varid,"dble", &
         & "lon_o","longitude","degree E")

    call define_var_netcdf(ncid,1,dim(2),varid,"dble", &
         & "lat_o","latitude","degree N")  

    call define_var_netcdf(ncid,1,dim(3),varid,"dble", &
         & "dep_o","depth","m")  

    if(varname == "t")then
    
       call define_var_netcdf(ncid,2,dim(3:4),varid,"dble", &
            & "htmean_a","analysis T in obs. space","degree C")

       call define_var_netcdf(ncid,2,dim(3:4),varid,"dble", &
            & "htsprd_a","analysis T spread in obs. space","degree C")

       call define_var_netcdf(ncid,2,dim(3:4),varid,"dble", &
            & "t_o","observed T","degree C")
       
    else if(varname == "s")then

       call define_var_netcdf(ncid,2,dim(3:4),varid,"dble", &
            & "hsmean_a","analysis S in obs. space","-")

       call define_var_netcdf(ncid,2,dim(3:4),varid,"dble", &
            & "hssprd_a","analysis S spread in obs. space","-")

       call define_var_netcdf(ncid,2,dim(3:4),varid,"dble", &
            & "s_o","observed S","-")
       
    else if(varname == "u")then

       call define_var_netcdf(ncid,2,dim(3:4),varid,"dble", &
            & "humean_a","analysis U in obs. space","m/s")       

       call define_var_netcdf(ncid,2,dim(3:4),varid,"dble", &
            & "husprd_a","analysis U spread in obs. space","m/s")

       call define_var_netcdf(ncid,2,dim(3:4),varid,"dble", &
            & "u_o","observed U","m/s")
       
    else if(varname == "v")then

       call define_var_netcdf(ncid,2,dim(3:4),varid,"dble", &
            & "hvmean_a","analysis V in obs. space","m/s")

       call define_var_netcdf(ncid,2,dim(3:4),varid,"dble", &
            & "hvsprd_a","analysis V spread in obs. space","m/s")

       call define_var_netcdf(ncid,2,dim(3:4),varid,"dble", &
            & "v_o","observed V","m/s")
       
    end if
    
    status=nf90_enddef(ncid)
    call check_error(status)

    status=nf90_close(ncid)
    call check_error(status)

  end subroutine make_ncfile

end module mod_make_ncfile
