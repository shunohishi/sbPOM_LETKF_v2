module mod_read_woa18_month
  
  integer,parameter :: im=360,jm=180,km=57
  
contains

  !----------------------------------------------------------------
  ! Read WOA18 |
  !-------------
  !
  ! 1. Read WOA18 t,s
  ! 2. Modify missing value
  ! 3. Modify long direction
  !
  ! 2018.07 created by S. Ohishi
  ! 2023.03 modified by S. Ohishi
  !
  !----------------------------------------------------------------
  
  subroutine read_woa18_month(imon,long,lati,depth,t,s)

    !$use omp_lib    
    use mod_rmiss
    use netcdf
    implicit none

    real(kind = 4),parameter :: dmiss=9.96921e+36

    integer i,j,k,itmp
    integer status,access,ncid,varid

    real(kind = 4) tmp1dx(im),tmp1dy(jm),tmp1dz(km)
    real(kind = 4) ttmp(im,jm,km),stmp(im,jm,km)

    character(2) mm
    character(100) tfilename,sfilename

    !IN
    integer,intent(in) :: imon

    !OUT
    real(kind = 8),intent(out) :: long(im),lati(jm),depth(km)
    real(kind = 8),intent(out) :: t(im,jm,km),s(im,jm,km)
    
    !Filename
    write(mm,'(i2.2)') imon
    tfilename="/data/R/R2402/DATA/WOA18/woa18_decav_t"//mm//"_01.nc"
    sfilename="/data/R/R2402/DATA/WOA18/woa18_decav_s"//mm//"_01.nc"
      
    !---Check filename
    status=access(trim(tfilename)," ")
    if(status /= 0)then
       write(*,*) "Not found:"//trim(tfilename)
       stop
    end if
    
    status=access(trim(sfilename)," ")
    if(status /= 0)then
       write(*,*) "Not found:"//trim(sfilename)
       stop
    end if
       
    !---Read Positon & Temperature
    status=nf90_open(trim(tfilename),nf90_nowrite,ncid)
    
    status=nf90_inq_varid(ncid,"lon",varid)
    status=nf90_get_var(ncid,varid,tmp1dx)

    status=nf90_inq_varid(ncid,"lat",varid)
    status=nf90_get_var(ncid,varid,tmp1dy)

    status=nf90_inq_varid(ncid,"depth",varid)
    status=nf90_get_var(ncid,varid,tmp1dz)

    status=nf90_inq_varid(ncid,"t_an",varid)
    status=nf90_get_var(ncid,varid,ttmp)
    
    status=nf90_close(ncid)

    !---Read Salinity
    status=nf90_open(trim(sfilename),nf90_nowrite,ncid)
    
    status=nf90_inq_varid(ncid,"s_an",varid)
    status=nf90_get_var(ncid,varid,stmp)

    status=nf90_close(ncid)

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
       if(-180.e0 <= tmp1dx(itmp) .and. tmp1dx(itmp)< 0.e0)then
          i=i+1
          long(i)=dble(tmp1dx(itmp))+360.d0
       end if
    end do

    !Latitude
    lati(:)=dble(tmp1dy(:))

    !Depth
    depth(:)=dble(tmp1dz(:))

    !T & S
    !$omp parallel
    !$omp do private(i,j,k)  
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
  !$omp end do
  !$omp end parallel
    
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

  end subroutine read_woa18_month

  !---------------------------------------------------------------

  subroutine read_woa18_month_var(varname,imon,long,lati,depth,dat)

    !$use omp_lib      
    use mod_rmiss
    use netcdf
    implicit none

    real(kind = 4),parameter :: dmiss=9.96921e+36

    integer i,j,k,itmp
    integer status,access,ncid,varid

    real(kind = 4) tmp1dx(im),tmp1dy(jm),tmp1dz(km)
    real(kind = 4) tmp(im,jm,km)

    character(2) mm
    character(100) filename

    !IN
    integer,intent(in) :: imon
    character(1),intent(in) :: varname

    !OUT
    real(kind = 8),intent(out) :: long(im),lati(jm),depth(km)
    real(kind = 8),intent(out) :: dat(im,jm,km)
    
    !Filename
    write(mm,'(i2.2)') imon
    if(varname == "t")then
       filename="/data/R/R2402/DATA/WOA18/woa18_decav_t"//mm//"_01.nc"
    else if(varname == "s")then
       filename="/data/R/R2402/DATA/WOA18/woa18_decav_s"//mm//"_01.nc"
    end if
       
    !---Check filename
    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "Not found:"//trim(filename)
       stop
    end if

    !---Open NetCDF file
    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    !---Read Positon    
    status=nf90_inq_varid(ncid,"lon",varid)
    status=nf90_get_var(ncid,varid,tmp1dx)

    status=nf90_inq_varid(ncid,"lat",varid)
    status=nf90_get_var(ncid,varid,tmp1dy)

    status=nf90_inq_varid(ncid,"depth",varid)
    status=nf90_get_var(ncid,varid,tmp1dz)

    !---Read Data
    if(varname == "t")then
       status=nf90_inq_varid(ncid,"t_an",varid)
    else if(varname == "s")then
       status=nf90_inq_varid(ncid,"s_an",varid)
    end if
    status=nf90_get_var(ncid,varid,tmp)

    !---Close NetCDF file
    status=nf90_close(ncid)

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
       if(-180.e0 <= tmp1dx(itmp) .and. tmp1dx(itmp)< 0.e0)then
          i=i+1
          long(i)=dble(tmp1dx(itmp))+360.d0
       end if
    end do

    !Latitude
    lati(:)=dble(tmp1dy(:))

    !Depth
    depth(:)=dble(tmp1dz(:))
    
    !T & S
    !$omp parallel
    !$omp do private(i,j,k)  
    do k=1,km
       do j=1,jm
          do i=1,im
             
             if(tmp(i,j,k) == dmiss)then
                tmp(i,j,k)=real(rmiss)
             end if
             
          end do
       end do
    end do
    !$omp end do
    !$omp end parallel
    
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

  end subroutine read_woa18_month_var

  !---------------------------------------------------------------------
  
  subroutine read_woa18_month_err(imon,long,lati,depth,tse,sse)

    !$use omp_lib
    use mod_rmiss
    use netcdf
    implicit none

    real,parameter :: dmiss=9.96921e+36

    integer i,j,k,itmp
    integer status,access,ncid,varid
    
    real(kind = 4) tmp1dx(im),tmp1dy(jm),tmp1dz(km)
    real(kind = 4) tsetmp(im,jm,km),ssetmp(im,jm,km)
    real(kind = 4) tsdtmp(im,jm,km),ssdtmp(im,jm,km)

    real(kind = 8) tsd(im,jm,km),ssd(im,jm,km) !Standard deviation
    
    character(2) mm
    character(100) tfilename,sfilename

    !IN
    integer,intent(in) :: imon

    !OUT
    real(kind = 8),intent(out) :: long(im),lati(jm),depth(km)
    real(kind = 8),intent(out) :: tse(im,jm,km),sse(im,jm,km) !Standard error
    
    !---Filename
    write(mm,'(i2.2)') imon
    tfilename="/data/R/R2402/DATA/WOA18/woa18_decav_t"//mm//"_01.nc"
    sfilename="/data/R/R2402/DATA/WOA18/woa18_decav_s"//mm//"_01.nc"
      
    !---Check filename
    status=access(trim(tfilename)," ")
    if(status /= 0)then
       write(*,*) "Not found:"//trim(tfilename)
       stop
    end if
    
    status=access(trim(sfilename)," ")
    if(status /= 0)then
       write(*,*) "Not found:"//trim(sfilename)
       stop
    end if
       
    !---Read Positon & Temperature
    status=nf90_open(trim(tfilename),nf90_nowrite,ncid)
    
    status=nf90_inq_varid(ncid,"lon",varid)
    status=nf90_get_var(ncid,varid,tmp1dx)

    status=nf90_inq_varid(ncid,"lat",varid)
    status=nf90_get_var(ncid,varid,tmp1dy)

    status=nf90_inq_varid(ncid,"depth",varid)
    status=nf90_get_var(ncid,varid,tmp1dz)

    status=nf90_inq_varid(ncid,"t_se",varid)
    status=nf90_get_var(ncid,varid,tsetmp)

    status=nf90_inq_varid(ncid,"t_sd",varid)
    status=nf90_get_var(ncid,varid,tsdtmp)

    status=nf90_close(ncid)

    !Read Salinity
    status=nf90_open(trim(sfilename),nf90_nowrite,ncid)
    
    status=nf90_inq_varid(ncid,"s_se",varid)
    status=nf90_get_var(ncid,varid,ssetmp)

    status=nf90_inq_varid(ncid,"s_sd",varid)
    status=nf90_get_var(ncid,varid,ssdtmp)

    status=nf90_close(ncid)

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
       if(-180.e0 <= tmp1dx(itmp) .and. tmp1dx(itmp)< 0.e0)then
          i=i+1
          long(i)=dble(tmp1dx(itmp))+360.d0
       end if
    end do

    !Latitude
    lati(:)=dble(tmp1dy(:))

    !Depth
    depth(:)=dble(tmp1dz(:))

  !$omp parallel
  !$omp do private(i,j,k)    
    do k=1,km
       do j=1,jm
          do i=1,im

             if(tsetmp(i,j,k) == dmiss)then
                tsetmp(i,j,k)=real(rmiss)
             end if

             if(tsdtmp(i,j,k) == dmiss)then
                tsdtmp(i,j,k)=real(rmiss)
             end if
             
             if(ssetmp(i,j,k) == dmiss)then
                ssetmp(i,j,k)=real(rmiss)
             end if

             if(ssdtmp(i,j,k) == dmiss)then
                ssdtmp(i,j,k)=real(rmiss)
             end if

          end do
       end do
    end do
    !$omp end do
    !$omp end parallel
    
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
             if(-180.e0 <= tmp1dx(itmp) .and. tmp1dx(itmp) < 0.e0)then
                i=i+1
                tse(i,j,k)=dble(tsetmp(itmp,j,k))
                sse(i,j,k)=dble(ssetmp(itmp,j,k))
                tsd(i,j,k)=dble(tsdtmp(itmp,j,k))
                ssd(i,j,k)=dble(ssdtmp(itmp,j,k))
             end if             
          end do

       end do
    end do

  end subroutine read_woa18_month_err

end module mod_read_woa18_month
