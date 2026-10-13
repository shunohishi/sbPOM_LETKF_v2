module mod_read_etopo1

  integer,parameter :: im=21600,jm=10801

contains

  !------------------------------------------------------------
  ! Read Etopo1 |
  !--------------
  !
  ! 1. Modify missing value
  ! 2. Modify longitudinal direction
  !
  !------------------------------------------------------------

  subroutine read_etopo1(long,lati,dat)

    use mod_rmiss
    use mod_check_netcdf
    use netcdf
    implicit none

    integer(kind = 8),parameter :: dmiss=-2147483648_8

    integer i,j,itmp
    integer status,access,ncid,varid

    integer,allocatable :: idat(:,:)    
    real(kind = 8) dlon(im),dlat(jm)

    real(kind = 8),intent(out) :: long(im),lati(jm),dat(im,jm)
    
    
    character(100) filename

    filename="/data/R/R2402/DATA/ETOPO1/ETOPO1_Ice_g_gmt4.grd"
    
    status=access(trim(filename)," ")
    if(status /= 0) write(*,*) "Not Found:"//trim(filename)

    allocate(idat(im,jm))
    
    !Read position & Data
    status=nf90_open(trim(filename),nf90_nowrite,ncid)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"x",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,dlon)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"y",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,dlat)
    call check_netcdf(status)

    status=nf90_inq_varid(ncid,"z",varid)
    call check_netcdf(status)
    status=nf90_get_var(ncid,varid,idat)
    call check_netcdf(status)

    status=nf90_close(ncid)
    call check_netcdf(status)

    !---Post process
    !Longitude
    i=0

    do itmp=1,im
       if(0.d0 <= dlon(itmp) .and. dlon(itmp) < 180.d0)then
          i=i+1
          long(i)=dlon(itmp)
       end if
    end do

    do itmp=1,im
       if(-180.d0 <= dlon(itmp) .and. dlon(itmp) < 0.)then
          i=i+1
          long(i)=dlon(itmp)+360.d0
       end if
    end do

    !Latitude
    lati(:)=dlat(:)
    
    !Dat
    do j=1,jm
       do i=1,im
          if(idat(i,j) == dmiss)then
             idat(i,j)=nint(rmiss)
          end if
       end do
    end do

    do j=1,jm

       i=0

       do itmp=1,im
          if(0.d0 <= dlon(itmp) .and. dlon(itmp) < 180.d0)then
             i=i+1
             dat(i,j)=dble(idat(itmp,j))
          end if
       end do

       do itmp=1,im
          if(-180.d0 <= dlon(itmp) .and. dlon(itmp) < 0.)then
             i=i+1
             dat(i,j)=dble(idat(itmp,j))
          end if
       end do
       
    end do

    deallocate(idat)
    
  end subroutine read_etopo1

end module mod_read_etopo1
