module mod_read_gtspp

  !Modified S.Ohishi 2020.05

contains

  !-----------------------------------------------------------
  ! Make & Read filename |
  !-----------------------------------------------------------
  
  subroutine make_filename(iunit,nunit,iyr,imon)

    implicit none

    integer jyr,jmon,jday
    integer status,access,system
    integer ifile,nfile

    real(kind = 8) long,lati

    character(100) filename
    character(8) yyyymmdd
    character(6) yyyymm
    character(4) yyyy
    character(2) mm,dd

    !IN
    integer,intent(in) :: iunit,nunit
    integer,intent(in) :: iyr,imon
    
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    yyyymm=yyyy//mm

    status=access("/data/R/R2402/DATA/GTSPP/filename/"//yyyymm//".txt"," ")
    if(status == 0)then
       write(*,'(a)') "Exist /data/R/R2402/DATA/GTSPP/filename/"//yyyymm//".txt"
       return
    end if

    write(*,'(a)') "Start: Make monthly filename list at "//yyyymm
    status=system("ls -f /data/R/R2402/DATA/GTSPP/"//yyyymm// &
         & " > /data/R/R2402/DATA/GTSPP/filename/"//yyyymm//".txt")
    write(*,'(a)') "End: Make monthly filename list"//yyyymm

    nfile=0
    open(iunit,file="/data/R/R2402/DATA/GTSPP/filename/"//yyyymm//".txt",status="old")
    read(iunit,*)
    read(iunit,*)
    do 
       read(iunit,*,end=100)
       nfile=nfile+1
    end do
100 close(iunit)

    write(*,'(a,i10)') "Total number of file:",nfile

    write(*,'(a)') "Start: Make daily filename list at "//yyyymm
    open(iunit,file="/data/R/R2402/DATA/GTSPP/filename/"//yyyymm//".txt",status="old")
    read(iunit,*)
    read(iunit,*)
    do ifile=1,nfile

       if(mod(ifile,100) == 1) write(*,'(i10,a,i10)') ifile,"/",nfile
       read(iunit,'(a)') filename
       !write(*,*) trim(filename)
       call read_info(filename,iyr,imon,jyr,jmon,jday,long,lati)

       if(iyr == jyr .and. imon == jmon)then
          write(yyyy,'(i4.4)') jyr
          write(mm,'(i2.2)') jmon
          write(dd,'(i2.2)') jday
          yyyymmdd=yyyy//mm//dd
          open(iunit+nunit, &
               & file="/data/R/R2402/DATA/GTSPP/filename/"//yyyymmdd//".txt", &
               & access="append")
          write(iunit+nunit,'(2f12.5,x,a)') long,lati,trim(filename)
          close(iunit+nunit)
       else
          write(*,'(a)') "***Error: Not match date"
          stop
       end if
    end do
    close(iunit)

    write(*,'(a)') "End: Make daily filename list at "//yyyymm
    
  end subroutine make_filename

  !-------------------------

  subroutine read_filename(iyr,imon,iday,nfile,long,lati,filename)

    implicit none

    integer ifile
    integer status,access

    character(8) yyyymmdd
    character(4) yyyy
    character(2) mm,dd

    !IN
    integer,intent(in) :: iyr,imon,iday

    !OUT
    integer,intent(out) :: nfile
    real(kind = 8),allocatable,intent(out) :: long(:),lati(:)
    character(100),allocatable,intent(out) :: filename(:)


    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    yyyymmdd=yyyy//mm//dd

    status=access("/data/R/R2402/DATA/GTSPP/filename/"//yyyymmdd//".txt"," ")
    if(status /= 0)then
       write(*,'(a)') "***Error: Not found "//"/data/R/R2402/DATA/GTSPP/filename/"&
            &//yyyymmdd//".txt"
       stop
    end if

    nfile=0
    open(1,file="/data/R/R2402/DATA/GTSPP/filename/"//yyyymmdd//".txt",status="old")
    do 
       read(1,*,end=100)
       nfile=nfile+1
    end do
100 close(1)

    write(*,'(a,i10)') "The number of file:",nfile

    if(nfile == 0)then
       return
    else
       
       allocate(long(nfile),lati(nfile),filename(nfile))

       open(1,file="/data/R/R2402/DATA/GTSPP/filename/"//yyyymmdd//".txt",status="old")
       do ifile=1,nfile
          read(1,'(2f12.5,x,a)') long(ifile),lati(ifile),filename(ifile)
       end do
       close(1)

    end if

  end subroutine read_filename

  !-----------------------------------------------------------------
  ! Read information(date, position)
  !-----------------------------------------------------------------

  subroutine read_info(filename,iyr,imon,jyr,jmon,jday,long,lati)

    use netcdf
    implicit none

    integer status,access
    integer ncid,varid
    integer ijul

    real(kind = 4) rlon,rlat
    real(kind = 8) :: dtime

    character(200) fullfilename
    character(4) yyyy
    character(2) mm

    !IN
    integer,intent(in) :: iyr,imon
    character(100),intent(in) :: filename

    !OUT
    integer,intent(out) :: jyr,jmon,jday
    real(kind = 8),intent(out) :: long,lati

    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon

    fullfilename="/data/R/R2402/DATA/GTSPP/"//yyyy//mm//"/"//trim(filename)

    status=access(trim(fullfilename)," ")
    if(status /= 0)then
       write(*,'(a)') "***Error: Not found "//trim(fullfilename)
       stop
    end if

    status=nf90_open(trim(fullfilename),nf90_nowrite,ncid)

    status=nf90_inq_varid(ncid,"time",varid)
    status=nf90_get_var(ncid,varid,dtime)

    status=nf90_inq_varid(ncid,"longitude",varid)
    status=nf90_get_var(ncid,varid,rlon)

    status=nf90_inq_varid(ncid,"latitude",varid)
    status=nf90_get_var(ncid,varid,rlat)

    status=nf90_close(ncid)

    long=dble(rlon)
    lati=dble(rlat)
    
    if(long < 0.d0)then
       long=long+360.d0
    end if

    call ymd_julian(1900,1,1,ijul)
    call julian_ymd(int(dtime)+ijul,jyr,jmon,jday)

  end subroutine read_info

  !----------------------------------------------------------------
  ! Read data |
  !------------
  !
  ! <flag>
  ! 0: No QC done, 1: Good data, 2:Probably good data
  ! 3: Probably bad data, 4: bad data
  ! 5: Changed, 6-8 Reserved, 9: Element missing
  !
  !----------------------------------------------------------------
  
  subroutine read_gtspp(iyr,imon,iday,filename,km,long,lati,depth,t,s)

    use mod_rmiss
    use netcdf
    implicit none

    integer,parameter :: iflag=1    
    real(kind = 4),parameter :: dmiss=99999.e0

    integer status,access
    integer ncid,dimid,varid
    integer jyr,jmon,jday
    integer ifile,nfile
    integer k

    integer,allocatable :: dflag(:),tflag(:),sflag(:)

    real(kind = 4),allocatable :: rdep(:),rt(:),rs(:)
    
    character(200) fullfilename
    character(4) yyyy
    character(2) mm

    !IN
    integer,intent(in) :: iyr,imon,iday
    character(100),intent(in) :: filename

    !OUT
    integer,intent(out) :: km
    real(kind = 8),allocatable,intent(out) :: depth(:),t(:),s(:)

    !INOUT
    real(kind = 8),intent(inout) :: long,lati
    
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    fullfilename="/data/R/R2402/DATA/GTSPP/"//yyyy//mm//"/"//trim(filename)

    status=access(trim(fullfilename)," ")
    if(status /= 0)then
       write(*,'(a)') "***Error: Not Found "//trim(fullfilename)
       stop
    end if

    call read_info(filename,iyr,imon,jyr,jmon,jday,long,lati)

    status=nf90_open(trim(fullfilename),nf90_nowrite,ncid)

    !Allocate
    status=nf90_inq_dimid(ncid,"z",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = km)
    allocate(rdep(km),rt(km),rs(km))
    allocate(depth(km),t(km),s(km))
    allocate(dflag(km),tflag(km),sflag(km))

    !Read data
    !depth
    status=nf90_inq_varid(ncid,"z",varid)
    if(status == nf90_noerr)then
       status=nf90_get_var(ncid,varid,rdep)
       depth(:)=dble(rdep(:))
    else
       depth(:)=rmiss
    end if

    status=nf90_inq_varid(ncid,"z_variable_quality_flag",varid)
    if(status == nf90_noerr)then
       status=nf90_get_var(ncid,varid,dflag)
    else
       dflag(:)=rmiss
    end if

    !temperature
    status=nf90_inq_varid(ncid,"temperature",varid)
    if(status == nf90_noerr)then
       status=nf90_get_var(ncid,varid,rt,(/1,1,1,1/),(/1,1,km,1/))
       t(:)=dble(rt(:))
    else
       t(:)=rmiss
    end if

    status=nf90_inq_varid(ncid,"temperature_quality_flag",varid)
    if(status == nf90_noerr)then
       status=nf90_get_var(ncid,varid,tflag)
    else
       tflag(:)=rmiss
    end if

    !salinity
    status=nf90_inq_varid(ncid,"salinity",varid)
    if(status == nf90_noerr)then
       status=nf90_get_var(ncid,varid,rs,(/1,1,1,1/),(/1,1,km,1/))
       s(:)=dble(rs(:))
    else
       s(:)=rmiss
    end if

    status=nf90_inq_varid(ncid,"salinity_quality_flag",varid)
    if(status == nf90_noerr)then
       status=nf90_get_var(ncid,varid,sflag)
    else
       sflag(:)=rmiss
    end if

    status=nf90_close(ncid)

    !QC
    do k=1,km
       if(tflag(k) /= iflag) t(k)=rmiss
       if(sflag(k) /= iflag) s(k)=rmiss
    end do

    deallocate(rdep,rt,rs)
    deallocate(dflag,tflag,sflag)

  end subroutine read_gtspp

  !-------------------------------

  subroutine deallocate_gtspp(depth,t,s)

    implicit none
    real(kind = 8),allocatable :: depth(:),t(:),s(:)

    deallocate(depth,t,s)

  end subroutine deallocate_gtspp

  !-------------------------------

  subroutine deallocate_gtspp_filename(long,lati,filename)

    implicit none

    real(kind = 8),allocatable :: long(:),lati(:)
    character(100),allocatable :: filename(:)

    deallocate(filename)

  end subroutine deallocate_gtspp_filename

end module mod_read_gtspp
