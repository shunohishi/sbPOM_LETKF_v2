module mod_read_woa23_season_annual

  integer,parameter :: im=1440,jm=720,km=102

contains

  !----------------------------------------------------------------
  ! Read WOA23 |
  !-------------
  !
  ! 1. Read WOA23 t,s
  ! 2. Modify missing value
  ! 3. Modify long direction
  ! * 0: Annual mean
  ! * 13(1-3):Winter, 14(4-6):Spring, 15(7-9):Summer, 16(10-12):Autumn
  !
  ! 2018.07 created by S. Ohishi
  ! 2023.03 modified by S. Ohishi
  !
  !----------------------------------------------------------------

  subroutine read_woa23_season_annual(imon,long,lati,depth,t,s)

    use mod_rmiss
    use mod_check_netcdf
    use netcdf
    implicit none

    real,parameter :: dmiss=9.96921e+36

    integer jmon
    integer i,j,k,itmp
    integer status,access,ncid,varid

    real tmp1dx(im),tmp1dy(jm),tmp1dz(km)
    real ttmp(im,jm,km),stmp(im,jm,km)

    character(2) mm
    character(100) tfilename,sfilename

    integer,intent(in) :: imon
    real(kind = 8),intent(out) :: long(im),lati(jm),depth(km)
    real(kind = 8),intent(out) :: t(im,jm,km),s(im,jm,km)

    !Filename
    if(imon == 0)then
       jmon=0
    elseif(1 <= imon .and. imon <= 3)then
       jmon=13
    else if(4 <= imon .and. imon <= 6)then
       jmon=14
    else if(7 <= imon .and. imon <= 9)then
       jmon=15
    else if(10 <= imon .and. imon <= 12)then
       jmon=16
    end if

    write(mm,'(i2.2)') jmon
    tfilename="/data/R/R2402/DATA/WOA23/woa23_decav_t"//mm//"_04.nc"
    sfilename="/data/R/R2402/DATA/WOA23/woa23_decav_s"//mm//"_04.nc"

    !Check filename
    status=access(trim(tfilename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not found "//trim(tfilename)
       stop
    end if
    status=access(trim(sfilename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not found "//trim(sfilename)
       stop
    end if

    !Read Positon & Temperature
    status=nf90_open(trim(tfilename),nf90_nowrite,ncid)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"lon",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,tmp1dx)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"lat",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,tmp1dy)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"depth",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,tmp1dz)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"t_an",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,ttmp)
    call check_netcdf(status)

    status=nf90_close(ncid)
    call check_netcdf(status)

    !Read Salinity
    status=nf90_open(trim(sfilename),nf90_nowrite,ncid)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"s_an",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,stmp)
    call check_netcdf(status)

    status=nf90_close(ncid)
    call check_netcdf(status)

    !---Post process
    !Longitude

    i=0

    do itmp=1,im
       if(0.e0 <= tmp1dx(itmp) .and. tmp1dx(itmp) < 180.e0)then
          i=i+1
          long(i)=dble(tmp1dx(itmp))
       end if
    end do

    do itmp=1,im
       if(-180.e0 <= tmp1dx(itmp) .and. tmp1dx(itmp)< 0.)then
          i=i+1
          long(i)=dble(tmp1dx(itmp))+360.d0
       end if
    end do

    !Latitude
    lati(:)=dble(tmp1dy(:))
    !Depth
    depth(:)=dble(tmp1dz(:))

    !T & S
    do k=1,km
       do j=1,jm
          do i=1,im

             if(ttmp(i,j,k) == dmiss)then
                ttmp(i,j,k)=real(rmiss)
             end if

             if(stmp(i,j,k) == dmiss)then
                stmp(i,j,k)=real(rmiss)
             end if

          end do
       end do
    end do


    do k=1,km
       do j=1,jm

          i=0

          do itmp=1,im
             if(0.e0 <= tmp1dx(itmp) .and. tmp1dx(itmp) < 180.e0)then
                i=i+1
                t(i,j,k)=dble(ttmp(itmp,j,k))
                s(i,j,k)=dble(stmp(itmp,j,k))
             end if
          end do

          do itmp=1,im
             if(-180.e0 <= tmp1dx(itmp) .and. tmp1dx(itmp)< 0.e0)then
                i=i+1
                t(i,j,k)=dble(ttmp(itmp,j,k))
                s(i,j,k)=dble(stmp(itmp,j,k))
             end if
          end do

       end do
    end do

  end subroutine read_woa23_season_annual

  !-----------------------------------------------------------

  subroutine read_woa23_season_annual_var(varname,imon,long,lati,depth,dat)

    use mod_rmiss
    use mod_check_netcdf
    use netcdf
    implicit none

    real,parameter :: dmiss=9.96921e+36

    integer jmon
    integer i,j,k,itmp
    integer status,access,ncid,varid

    real tmp1dx(im),tmp1dy(jm),tmp1dz(km)
    real tmp(im,jm,km)

    character(2) mm
    character(100) filename

    !---IN
    integer,intent(in) :: imon
    character(1),intent(in) :: varname

    !---OUT
    real(kind = 8),intent(out) :: long(im),lati(jm),depth(km)
    real(kind = 8),intent(out) :: dat(im,jm,km)

    !Filename
    if(imon == 0)then
       jmon=0
    elseif(1 <= imon .and. imon <= 3)then
       jmon=13
    else if(4 <= imon .and. imon <= 6)then
       jmon=14
    else if(7 <= imon .and. imon <= 9)then
       jmon=15
    else if(10 <= imon .and. imon <= 12)then
       jmon=16
    end if

    write(mm,'(i2.2)') jmon
    if(varname == "t")then
       filename="/data/R/R2402/DATA/WOA23/woa23_decav_t"//mm//"_04.nc"
    else if(varname == "s")then
       filename="/data/R/R2402/DATA/WOA23/woa23_decav_s"//mm//"_04.nc"
    end if

    !Check filename
    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not found "//trim(filename)
       stop
    end if

    !Read Positon
    status=nf90_open(trim(filename),nf90_nowrite,ncid)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"lon",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,tmp1dx)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"lat",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,tmp1dy)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"depth",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,tmp1dz)
    call check_netcdf(status)

    if(varname == "t")then
       status=nf90_inq_varid(ncid,"t_an",varid)
       call check_netcdf(status)
    else if(varname == "s")then
       status=nf90_inq_varid(ncid,"s_an",varid)
       call check_netcdf(status)
    end if
    status=nf90_get_var(ncid,varid,tmp)
    call check_netcdf(status)

    status=nf90_close(ncid)
    call check_netcdf(status)

    !---Post process
    !Longitude

    i=0

    do itmp=1,im
       if(0.e0 <= tmp1dx(itmp) .and. tmp1dx(itmp) < 180.e0)then
          i=i+1
          long(i)=dble(tmp1dx(itmp))
       end if
    end do

    do itmp=1,im
       if(-180.e0 <= tmp1dx(itmp) .and. tmp1dx(itmp)< 0.)then
          i=i+1
          long(i)=dble(tmp1dx(itmp))+360.d0
       end if
    end do

    !Latitude
    lati(:)=dble(tmp1dy(:))
    !Depth
    depth(:)=dble(tmp1dz(:))

    !T & S
    do k=1,km
       do j=1,jm
          do i=1,im

             if(tmp(i,j,k) == dmiss)then
                tmp(i,j,k)=real(rmiss)
             end if

          end do
       end do
    end do


    do k=1,km
       do j=1,jm

          i=0

          do itmp=1,im
             if(0.e0 <= tmp1dx(itmp) .and. tmp1dx(itmp) < 180.e0)then
                i=i+1
                dat(i,j,k)=dble(tmp(itmp,j,k))
             end if
          end do

          do itmp=1,im
             if(-180.e0 <= tmp1dx(itmp) .and. tmp1dx(itmp)< 0.e0)then
                i=i+1
                dat(i,j,k)=dble(tmp(itmp,j,k))
             end if
          end do

       end do
    end do

  end subroutine read_woa23_season_annual_var

  !-----------------------------------------------------------

  subroutine read_woa23_season_annual_err(imon,long,lati,depth,tse,sse)

    use mod_rmiss
    use mod_check_netcdf
    use netcdf
    implicit none

    real,parameter :: dmiss=9.96921e+36

    integer jmon
    integer i,j,k,itmp
    integer status,access,ncid,varid

    real tmp1dx(im),tmp1dy(jm),tmp1dz(km)
    real tsetmp(im,jm,km),ssetmp(im,jm,km)
    real tsdtmp(im,jm,km),ssdtmp(im,jm,km)

    real(kind = 8) tsd(im,jm,km),ssd(im,jm,km)

    character(2) mm
    character(100) tfilename,sfilename

    integer,intent(in) :: imon
    real(kind = 8),intent(out) :: long(im),lati(jm),depth(km)
    real(kind = 8),intent(out) :: tse(im,jm,km),sse(im,jm,km)


    !Filename
    if(imon == 0)then
       jmon=0
    elseif(1 <= imon .and. imon <= 3)then
       jmon=13
    else if(4 <= imon .and. imon <= 6)then
       jmon=14
    else if(7 <= imon .and. imon <= 9)then
       jmon=15
    else if(10 <= imon .and. imon <= 12)then
       jmon=16
    end if

    write(mm,'(i2.2)') jmon
    tfilename="/data/R/R2402/DATA/WOA23/woa23_decav_t"//mm//"_04.nc"
    sfilename="/data/R/R2402/DATA/WOA23/woa23_decav_s"//mm//"_04.nc"

    !Check filename
    status=access(trim(tfilename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not found "//trim(tfilename)
       stop
    end if
    status=access(trim(sfilename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not found "//trim(sfilename)
       stop
    end if

    !Read Positon & Temperature
    status=nf90_open(trim(tfilename),nf90_nowrite,ncid)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"lon",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,tmp1dx)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"lat",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,tmp1dy)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"depth",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,tmp1dz)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"t_se",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,tsetmp)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"t_sd",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,tsdtmp)
    call check_netcdf(status)

    status=nf90_close(ncid)
    call check_netcdf(status)

    !Read Salinity
    status=nf90_open(trim(sfilename),nf90_nowrite,ncid)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"s_se",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,ssetmp)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"s_sd",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,ssdtmp)
    call check_netcdf(status)

    status=nf90_close(ncid)
    call check_netcdf(status)

    !Post process

    i=0

    do itmp=1,im
       if(0.e0 <= tmp1dx(itmp) .and. tmp1dx(itmp) < 180.e0)then
          i=i+1
          long(i)=dble(tmp1dx(itmp))
       end if
    end do

    do itmp=1,im
       if(-180.e0 <= tmp1dx(itmp) .and. tmp1dx(itmp)< 0.e0)then
          i=i+1
          long(i)=dble(tmp1dx(itmp))+360.d0
       end if
    end do

    lati(:)=dble(tmp1dy(:))
    depth(:)=dble(tmp1dz(:))

    do k=1,km
       do j=1,jm
          do i=1,im

             if(tsetmp(i,j,k) == dmiss)then
                tsetmp(i,j,k)=real(rmiss)
             end if

             if(ssetmp(i,j,k) == dmiss)then
                ssetmp(i,j,k)=real(rmiss)
             end if

             if(tsdtmp(i,j,k) == dmiss)then
                tsdtmp(i,j,k)=real(rmiss)
             end if

             if(ssdtmp(i,j,k) == dmiss)then
                ssdtmp(i,j,k)=real(rmiss)
             end if

          end do
       end do
    end do


    do k=1,km
       do j=1,jm

          i=0

          do itmp=1,im
             if(0.e0 <= tmp1dx(itmp) .and. tmp1dx(itmp) < 180.e0)then
                i=i+1
                tse(i,j,k)=dble(tsetmp(itmp,j,k))
                sse(i,j,k)=dble(ssetmp(itmp,j,k))
                tsd(i,j,k)=dble(tsdtmp(itmp,j,k))
                ssd(i,j,k)=dble(ssdtmp(itmp,j,k))
             end if
          end do

          do itmp=1,im
             if(-180.e0 <= tmp1dx(itmp) .and. tmp1dx(itmp)< 0.e0)then
                i=i+1
                tse(i,j,k)=dble(tsetmp(itmp,j,k))
                sse(i,j,k)=dble(ssetmp(itmp,j,k))
                tsd(i,j,k)=dble(tsdtmp(itmp,j,k))
                ssd(i,j,k)=dble(ssdtmp(itmp,j,k))
             end if
          end do

       end do
    end do

  end subroutine read_woa23_season_annual_err

end module mod_read_woa23_season_annual
