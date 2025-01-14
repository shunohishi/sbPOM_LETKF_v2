module mod_read_windsat

contains

  !---------------------------------------------------
  ! Read WindSat |
  !-------------
  ! Created by S.Ohishi 2024.04 based on AMSR2 code
  !---------------------------------------------------

  subroutine read_windsat_filename(iyr,imon,iday,nfile,filename)

    implicit none

    !Common
    integer status,access,system
    integer ifile

    character(8) yyyymmdd
    character(6) yyyymm
    character(4) yyyy
    character(2) mm,dd

    !IN
    integer,intent(in) :: iyr,imon,iday

    !OUT
    integer,intent(out) :: nfile
    character(200),allocatable,intent(out) :: filename(:)

    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    yyyymm=yyyy//mm
    yyyymmdd=yyyy//mm//dd

    status=access("/data/R/R2402/DATA/WindSat/"//yyyymm," ")
    if(status /= 0)then
       nfile=0
       return
    end if
    
    !Make filename file
    status=system("find /data/R/R2402/DATA/WindSat/"//yyyymm// &
         & " -name "//yyyymmdd//&
         &"*-JAXA-L2P_GHRSST-SSTsubskin-WINDSAT-v8_*-v02.0-fv01.0.nc "//&
         &"> windsat"//yyyymmdd//".dat")
    
    !Count file
    nfile=0
    open(1,file="windsat"//yyyymmdd//".dat",status="old")
    do
       read(1,*,end=100)
       nfile=nfile+1
    end do
100 close(1)

    write(*,*) "The number of file:",nfile,"at "//yyyymmdd

    !Get filename
    if(nfile /= 0)then
       allocate(filename(nfile))

       open(1,file="windsat"//yyyymmdd//".dat",status="old")
       do ifile=1,nfile
          read(1,'(a)') filename(ifile)
       end do
       close(1)
       status=system("rm -f windsat"//yyyymmdd//".dat")

    end if
    status=system("rm -f windsat"//yyyymmdd//".dat")

  end subroutine read_windsat_filename

  !-------------------

  subroutine read_windsat(filename,iyr,imon,iday,ihour,im,jm,long,lati,dat)

    use netcdf
    use mod_rmiss
    implicit none

    !Common
    real(kind = 4),parameter :: dmiss=-32768.e0,dql=5.e0

    integer i,j
    integer ifile,nfile
    integer ijul
    integer itime,ntime

    integer status,system
    integer ncid,dimid,varid

    real(kind = 4),allocatable :: rlon(:,:),rlat(:,:),rdat(:,:)
    real(kind = 4),allocatable :: rql(:,:)

    !IN
    character(200),intent(in) :: filename

    !OUT
    integer,intent(out) :: im,jm
    integer,intent(out) ::  iyr,imon,iday,ihour
    real(kind = 8),allocatable,intent(out) :: long(:,:),lati(:,:),dat(:,:)

    write(*,'(a)') "Read: "//trim(filename)

    !Read ni,nj
    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    status=nf90_inq_dimid(ncid,"ni",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = im)

    status=nf90_inq_dimid(ncid,"nj",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = jm)

    !Allocate
    allocate(rlon(im,jm),rlat(im,jm),rdat(im,jm),rql(im,jm))
    allocate(long(im,jm),lati(im,jm),dat(im,jm))

    !Read data
    status=nf90_inq_varid(ncid,"time",varid)
    status=nf90_get_var(ncid,varid,itime)

    status=nf90_inq_varid(ncid,"lon",varid)
    status=nf90_get_var(ncid,varid,rlon)

    status=nf90_inq_varid(ncid,"lat",varid)
    status=nf90_get_var(ncid,varid,rlat)

    status=nf90_inq_varid(ncid,"sea_surface_temperature",varid)
    status=nf90_get_var(ncid,varid,rdat,(/1,1,1/),(/im,jm,1/))

    status=nf90_inq_varid(ncid,"quality_level",varid)
    status=nf90_get_var(ncid,varid,rql,(/1,1,1/),(/im,jm,1/))

    status=nf90_close(ncid)

    !Calculate time
    call ymd_julian(1981,1,1,ijul)
    call julian_ymd(itime/86400+ijul,iyr,imon,iday)
    ihour=nint((dble(itime)/86400.d0-itime/86400)*24.d0)

    !Post process
    long(:,:)=dble(rlon(:,:))
    lati(:,:)=dble(rlat(:,:))
    dat(:,:)=rmiss

    do j=1,jm
       do i=1,im

          if(long(i,j) < 0.d0)then
             long(i,j)=long(i,j)+360.d0
          end if

          if(rdat(i,j) == dmiss .or. rql(i,j) /= dql)then
             dat(i,j)=rmiss
          else
             dat(i,j)=dble(rdat(i,j))*0.01d0
          end if
       end do
    end do

    !deallocate
    deallocate(rlon,rlat,rdat,rql)

  end subroutine read_windsat

  !---------------------------------------

  subroutine deallocate_windsat(long,lati,dat)

    implicit none

    real(kind = 8),allocatable :: long(:,:),lati(:,:),dat(:,:)
    deallocate(long,lati,dat)

  end subroutine deallocate_windsat

  !-----

  subroutine deallocate_windsat_filename(filename)

    implicit none

    character(200),allocatable :: filename(:)
    deallocate(filename)

  end subroutine deallocate_windsat_filename

end module mod_read_windsat
