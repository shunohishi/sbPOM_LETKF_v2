module mod_read_cama

  integer,parameter :: im=1440,jm=720
  
contains

  !-----------------------------------------------
  ! Read Cama-Flood JRA55/GSMP |
  !-----------------------------
  !
  ! little endian (no header)
  ! units: m^3/s
  ! resolution: 0.25 degree global
  !             3 hourly
  !
  !------------------------------------------------
  
  !Land:1, Sea:0
  subroutine read_land_netcdf(land)

    use netcdf
    implicit none

    real(kind = 4),parameter :: dmiss=9.999e+20

    integer status,access
    integer ncid,varid
    integer i,j

    real(kind = 4) dland(im,jm)

    character(100) filename
    
    !OUT
    real(kind = 8),intent(out) :: land(im,jm)


    filename="/data/R/R2402/DATA/CaMa-Flood/netcdf/YEE2_JRA-55_outflw_H198101_GLB025.nc"
    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not found: "//trim(filename)
       stop
    end if
    
    !Read data
    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    status=nf90_inq_varid(ncid,"outflw",varid)
    status=nf90_get_var(ncid,varid,dland,(/1,1,1,1/),(/im,jm,1,1/))

    status=nf90_close(ncid)

    do j=1,jm
       do i=1,im
          if(dland(i,j) == dmiss)then
             land(i,j)=0.d0
          else
             land(i,j)=1.d0
          end if
       end do
    end do

  end subroutine read_land_netcdf

  !---------------------------

  subroutine read_cama_jra55_netcdf(iyr,imon,iday,ihour,long,lati,dat)

    use mod_rmiss
    use netcdf
    implicit none

    real(kind = 4),parameter :: dmiss=9.999e+20

    integer status,access,ncid,varid
    integer i,j

    real(kind = 4) tmp1dx(im),tmp1dy(jm),tmp2d(im,jm)
    
    character(100) filename
    character(4) yyyy
    character(2) mm
    
    !IN
    integer,intent(in) :: iyr,imon,iday,ihour

    !OUT
    real(kind = 8),intent(out) :: long(im),lati(jm)
    real(kind = 8),intent(out) :: dat(im,jm)

    
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon

    filename="/data/R/R2402/DATA/CaMa-Flood/netcdf/YEE2_JRA-55_outflw_H" &
         & //yyyy//mm//"_GLB025.nc"

    status=access(trim(filename)," ")
    if(status /= 0)then

       !If no data, dat=0.
       write(*,'(a)') "***Error: Not Found:"//trim(filename)

       filename="/data/R/R2402/DATA/CaMa-Flood/netcdf/YEE2_JRA-55_outflw_H200101_GLB025.nc "
       status=nf90_open(trim(filename),nf90_nowrite,ncid)
       
       status=nf90_inq_varid(ncid,"lon",varid)
       status=nf90_get_var(ncid,varid,tmp1dx)
       long(:)=dble(tmp1dx(:))
       
       status=nf90_inq_varid(ncid,"lat",varid)
       status=nf90_get_var(ncid,varid,tmp1dy)
       lati(:)=dble(tmp1dy(:))
       
       status=nf90_close(ncid)

       dat(:,:)=0.d0

       return
    end if

    !Read netcdf
    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    status=nf90_inq_varid(ncid,"lon",varid)
    status=nf90_get_var(ncid,varid,tmp1dx)
    long(:)=dble(tmp1dx(:))
    
    status=nf90_inq_varid(ncid,"lat",varid)
    status=nf90_get_var(ncid,varid,tmp1dy)
    lati(:)=dble(tmp1dy(:))

    status=nf90_inq_varid(ncid,"outflw",varid)
    status=nf90_get_var(ncid,varid,tmp2d,(/1,1,1,(iday-1)*8+ihour/3+1/),(/im,jm,1,1/))
    
    status=nf90_close(ncid)
    
    do j=1,jm
       do i=1,im
          if(tmp2d(i,j) == dmiss)then
             dat(i,j)=rmiss
          else
             dat(i,j)=dble(tmp2d(i,j))
          end if
       end do
    end do

    call modify_units(long,lati,dat)
    
  end subroutine read_cama_jra55_netcdf
  
  !--------------------------------------------
  ! Modify Unit [m^3/s] --> [mm/day]
  !--------------------------------------------

  subroutine modify_units(long,lati,dat)

    use mod_rmiss
    implicit none

    real(kind = 8),parameter :: pi=4.d0*atan(1.d0),a=6371.d3
    
    integer i,j

    real(kind = 8) dx(jm),dy

    
    real(kind = 8),intent(in) :: long(im),lati(jm)
    real(kind = 8),intent(inout) :: dat(im,jm)

    dy=0.25d0*pi*a/180.d0

    do j=1,jm
       dx(j)=0.25d0*pi*a*cos(pi*lati(j)/180.d0)/180.d0
    end do
    
    do j=1,jm
       do i=1,im
          if(dat(i,j)/=rmiss)then
             ![m^3/s] --> [mm * m^2/day]/(dx*dy) --> [mm/day]
             !dat(i,j)=dat(i,j)*1.d3*24.d0*60.d0*60.d0/(dx(j)*dy)
             dat(i,j)=dat(i,j)*86400.d3/(dx(j)*dy)
          end if
       end do
    end do

  end subroutine modify_units

end module mod_read_cama
