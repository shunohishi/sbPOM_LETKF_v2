module mod_read_cmems

  !Modified by S.Ohishi 2020.05
  !Modified by S.Ohishi 2022.06
  ! switch nrt/dt

contains

  !------------------------------------------------------------
  ! Read CMEMS |
  !-------------
  !
  ! ~2020.12: DT
  ! 2021.12~2022.03.28: NRT
  !
  !------------------------------------------------------------

  subroutine read_filename(iyr,imon,iday,nfile,filename)

    implicit none

    integer status,system
    integer ifile

    character(8) yyyymmdd
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
    yyyymmdd=yyyy//mm//dd
    
    if(iyr <= 2020)then
       status=system("find /data/R/R2402/DATA/CMEMS/dt/"//&
            &" -name *"//yyyymmdd//"_*.nc "//&
            &"> cmems"//yyyymmdd//".dat")
    else if(2021 <= iyr)then
       status=system("find /data/R/R2402/DATA/CMEMS/nrt/"//&
            &" -name *"//yyyymmdd//"_*.nc "//&
            &"> cmems"//yyyymmdd//".dat")
    end if

    nfile=0
    open(1,file="cmems"//yyyymmdd//".dat",status="old")
    do
       read(1,*,end=100)
       nfile=nfile+1
    end do
100 close(1)

    write(*,'(a,i10,a)') "The number of file:",nfile," at "//yyyymmdd

    if(nfile /= 0)then
       allocate(filename(nfile))
       open(1,file="cmems"//yyyymmdd//".dat",status="old")
       do ifile=1,nfile
          read(1,'(a)') filename(ifile)
       end do
       close(1)
    end if
    status=system("rm -f cmems"//yyyymmdd//".dat")

  end subroutine read_filename

  !------------------------------

  subroutine read_cmems(filename,ntime,long,lati,dat)
    
    use mod_rmiss
    use netcdf
    implicit none

    real(kind = 4),parameter :: dmiss=32767.e0

    integer itime
    integer status,system,access
    integer ncid,dimid,varid
    
    integer,allocatable :: ilon(:),ilat(:)
    real,allocatable :: rdat(:)
    
    !IN
    character(200),intent(in) :: filename

    !OUT
    integer,intent(out) :: ntime
    real(kind = 8),allocatable,intent(out) :: long(:),lati(:),dat(:)
    
    status=access(trim(filename)," ")
    if(status == 0)then
       write(*,*) "Read "//trim(filename)
    else
       write(*,*) "***Error: Not Found "//trim(filename)
       stop
    end if

    !Get ntime
    status=nf90_open(trim(filename),nf90_nowrite,ncid)
    
    status=nf90_inq_dimid(ncid,"time",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = ntime)

    !Allocate
    allocate(ilon(ntime),ilat(ntime),rdat(ntime))
    allocate(long(ntime),lati(ntime),dat(ntime))

    status=nf90_inq_varid(ncid,"longitude",varid)
    status=nf90_get_var(ncid,varid,ilon)

    status=nf90_inq_varid(ncid,"latitude",varid)
    status=nf90_get_var(ncid,varid,ilat)

    status=nf90_inq_varid(ncid,"sla_filtered",varid)
    status=nf90_get_var(ncid,varid,rdat)

    status=nf90_close(ncid)
    
    !Modify
    do itime=1,ntime
       
       long(itime)=dble(ilon(itime))*1.d-6
       if(long(itime) < 0.d0)then
          long(itime)=long(itime)+360.d0
       end if

       lati(itime)=dble(ilat(itime))*1.d-6

       if(rdat(itime) == dmiss)then
          dat(itime)=rmiss
       else
          dat(itime)=dble(rdat(itime))*0.001d0
       end if

    end do

    !deallocate
    deallocate(ilon,ilat,rdat)

  end subroutine read_cmems

  !---------------------------------------

  subroutine deallocate_cmems(long,lati,dat)

    implicit none
    
    real(kind = 8),allocatable :: long(:),lati(:),dat(:)
    deallocate(long,lati,dat)

  end subroutine deallocate_cmems

  !----------------------------------------

  subroutine deallocate_cmems_filename(filename)

    implicit none

    character(200),allocatable :: filename(:)
    deallocate(filename)

  end subroutine deallocate_cmems_filename

end module mod_read_cmems
