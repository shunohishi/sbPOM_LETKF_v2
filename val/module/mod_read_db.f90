module mod_read_db

  character(100),parameter :: db_dir="/data/R/R2402/DATA/DB"
  
contains

  !-----------------------------------------------------------------------------
  ! Make Drifter buoy filename|
  !-----------------------------------------------------------------------------

  subroutine make_filename

    implicit none

    integer status,access,system
    integer ifile,nfile

    character(200) command    
    character(100) filename

    integer iyr_min,imon_min,iday_min
    integer iyr_max,imon_max,iday_max

    real(kind = 8) lon_min,lon_max
    real(kind = 8) lat_min,lat_max
    
    status=access(trim(db_dir)//"/filename/filename.txt"," ")
    if(status == 0)then
       write(*,*) trim(db_dir)//"/filename/filename.txt founded"
       write(*,*) "Skip to create filename"
       return
    end if
    
    command="find "//trim(db_dir)//" -type f -name *.nc | "//&
         &"sed 's|"//trim(db_dir)//"/||' > "//trim(db_dir)//"/filename/filename.txt"
    
    status=system(trim(command))
    
    nfile=0
    open(1,file=trim(db_dir)//"/filename/filename.txt")
    do
       read(1,*,end=100)
       nfile=nfile+1
    end do
100 close(1)

    write(*,*) "The number of files:",nfile

    open(1,file=trim(db_dir)//"/filename/filename.txt",status="old")
    do ifile=1,nfile

       if(mod(ifile,100) == 1) write(*,'(i10,a,i10)') ifile,"/",nfile
       read(1,'(a)') filename

       call make_info(filename,&
            & iyr_min,imon_min,iday_min,&
            & iyr_max,imon_max,iday_max,&
            & lon_min,lon_max,lat_min,lat_max)

       if(ifile == 1)then
          open(11,file=trim(db_dir)//"/filename/info.dat",status="replace")
       else
          open(11,file=trim(db_dir)//"/filename/info.dat",access="append")
       end if
       write(11,'(6i6,4f12.5,x,a)') &
            & iyr_min,imon_min,iday_min, &
            & iyr_max,imon_max,iday_max, &
            & lon_min,lon_max,lat_min,lat_max, &
            & trim(filename)
       close(11)
       
    end do
    close(1)
    
  end subroutine make_filename

  !-----------------------------------------------------------------------------
  ! Make information of each buoy |
  !-----------------------------------------------------------------------------
  
  subroutine make_info(filename,&
            & iyr_min,imon_min,iday_min,&
            & iyr_max,imon_max,iday_max,&
            & lon_min,lon_max,lat_min,lat_max)

    use netcdf
    use mod_julian
    implicit none

    !---Parameter
    real(kind = 8),parameter :: dmiss=-999999.d0
    
    !---Common
    integer status,ncid,dimid,varid
    integer idb,ndb
    integer iobs,nobs

    integer sjul
    
    integer,allocatable :: ijul(:,:)
    
    real(kind = 8),allocatable :: time(:,:)
    real(kind = 8),allocatable :: lon(:,:),lat(:,:)

    !---IN
    character(100),intent(in) :: filename
    
    !---OUT
    integer,intent(out) :: iyr_min,imon_min,iday_min
    integer,intent(out) :: iyr_max,imon_max,iday_max
    
    real(kind = 8),intent(out) :: lon_min,lon_max,lat_min,lat_max
    
    !Open
    status=nf90_open(trim(db_dir)//"/"//trim(filename),nf90_nowrite,ncid)

    !Read dimension
    status=nf90_inq_dimid(ncid,"traj",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = ndb)

    status=nf90_inq_dimid(ncid,"obs",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = nobs)

    !Allocate
    allocate(ijul(nobs,ndb))
    allocate(time(nobs,ndb),lon(nobs,ndb),lat(nobs,ndb))

    !Read data
    status=nf90_inq_varid(ncid,"lon360",varid)
    status=nf90_get_var(ncid,varid,lon)

    status=nf90_inq_varid(ncid,"latitude",varid)
    status=nf90_get_var(ncid,varid,lat)

    status=nf90_inq_varid(ncid,"time",varid)
    status=nf90_get_var(ncid,varid,time)

    !Close
    status=nf90_close(ncid)

    !Julian day
    call ymd_julian(1970,1,1,sjul) !Reference date: 1970/01/01    
    ijul(:,:)=0
    
    do idb=1,ndb
       do iobs=1,nobs

          if(time(iobs,idb) == dmiss .or. lon(iobs,idb) == dmiss .or. lat(iobs,idb) == dmiss) cycle

          ijul(iobs,idb)=sjul+int(time(iobs,idb)/(24.d0*60.d0*60.d0)) !julian day
          
       end do
    end do

    call julian_ymd(minval(ijul(:,:), mask=ijul > 0),iyr_min,imon_min,iday_min)
    call julian_ymd(maxval(ijul(:,:), mask=ijul > 0),iyr_max,imon_max,iday_max)

    lon_min=minval(lon(:,:), mask=lon >= 0.d0)
    lon_max=maxval(lon(:,:), mask=lon <= 360.d0)
    lat_min=minval(lat(:,:), mask=lat >= -90.d0) 
    lat_max=maxval(lat(:,:), mask=lat <= 90.d0)

    deallocate(ijul)
    deallocate(time,lon,lat)
    
  end subroutine make_info

  !-------------------------------------------------------------------------------------
  ! Read info |
  !-------------------------------------------------------------------------------------

  subroutine read_info(nfile,filename, &
       & iyr_min,imon_min,iday_min, &
       & iyr_max,imon_max,iday_max, &
       & lon_min,lon_max,lat_min,lat_max)

    implicit none
    
    !---Common
    integer ifile
    
    !---OUT
    integer,intent(out) :: nfile

    integer,allocatable,intent(out) :: iyr_min(:),imon_min(:),iday_min(:)
    integer,allocatable,intent(out) :: iyr_max(:),imon_max(:),iday_max(:)

    real(kind = 8),allocatable,intent(out) :: lon_min(:),lon_max(:)
    real(kind = 8),allocatable,intent(out) :: lat_min(:),lat_max(:)

    character(100),allocatable,intent(out) :: filename(:)
    
    nfile=0
    open(1,file=trim(db_dir)//"/filename/info.dat",status="old")
    do
       read(1,*,end=100)
       nfile=nfile+1
    end do
100 close(1)

    write(*,*) "The number of file:", nfile

    allocate(iyr_min(nfile),imon_min(nfile),iday_min(nfile))
    allocate(iyr_max(nfile),imon_max(nfile),iday_max(nfile))
    allocate(lon_min(nfile),lon_max(nfile))
    allocate(lat_min(nfile),lat_max(nfile))
    allocate(filename(nfile))

    open(1,file=trim(db_dir)//"/filename/info.dat",status="old")
    do ifile=1,nfile
       read(1,'(6i6,4f12.5,x,a)') &
            & iyr_min(ifile),imon_min(ifile),iday_min(ifile), &
            & iyr_max(ifile),imon_max(ifile),iday_max(ifile), &
            & lon_min(ifile),lon_max(ifile),lat_min(ifile),lat_max(ifile), &
            & filename(ifile)
    end do
    close(1)
    
  end subroutine read_info

  !----------------------------------------
  
  subroutine deallocate_info(filename, &
       & iyr_min,imon_min,iday_min, &
       & iyr_max,imon_max,iday_max, &
       & lon_min,lon_max,lat_min,lat_max)

    implicit none

    integer,intent(inout),allocatable :: iyr_min(:),imon_min(:),iday_min(:)
    integer,intent(inout),allocatable :: iyr_max(:),imon_max(:),iday_max(:)

    real(kind = 8),intent(inout),allocatable :: lon_min(:),lon_max(:)
    real(kind = 8),intent(inout),allocatable :: lat_min(:),lat_max(:)

    character(100),intent(inout),allocatable :: filename(:)

    if(allocated(filename)) deallocate(filename)
    
    if(allocated(iyr_min)) deallocate(iyr_min)
    if(allocated(imon_min)) deallocate(imon_min)
    if(allocated(iday_min)) deallocate(iday_min)

    if(allocated(iyr_max)) deallocate(iyr_max)
    if(allocated(imon_max)) deallocate(imon_max)
    if(allocated(iday_max)) deallocate(iday_max)
    
    if(allocated(lon_min)) deallocate(lon_min)
    if(allocated(lon_max)) deallocate(lon_max)
    if(allocated(lat_min)) deallocate(lat_min)
    if(allocated(lat_max)) deallocate(lat_max)
    
  end subroutine deallocate_info
  
  !-------------------------------------------------------------------------------------
  ! Read data |
  !------------
  !
  ! Note: To compare daily-mean dataset, take 1d average in read_db_1d
  !
  !-------------------------------------------------------------------------------------

  subroutine read_db(filename,ndb,nobs,iyr,imon,iday,lon,lat,u,v,t)

    use netcdf
    use mod_rmiss
    use mod_julian
    implicit none

    !---Parameter
    real(kind = 8),parameter :: dmiss=-999999.d0
    real(kind = 8),parameter :: pi=4.d0*atan(1.d0)
    
    !---Common
    integer status,access
    integer ncid,dimid,varid
    integer idb,iobs
    integer sjul,ijul
    integer ihour
    
    real(kind = 8),allocatable :: time(:,:)
    
    !---IN
    character(100),intent(in) :: filename

    !---OUT
    integer,intent(out) :: ndb,nobs
    integer,intent(out),allocatable :: iyr(:,:),imon(:,:),iday(:,:)
    
    real(kind = 8),intent(out),allocatable :: lon(:,:),lat(:,:)
    real(kind = 8),intent(out),allocatable :: u(:,:),v(:,:),t(:,:)

    !Access
    status=access(trim(db_dir)//"/"//trim(filename)," ")
    if(status == 0)then
       !write(*,*) "Read:"//trim(db_dir)//"/"//trim(filename)
    else
       write(*,*) "***Error: "//trim(db_dir)//"/"//trim(filename)//"not found"
       stop
    end if
    
    !Open
    status=nf90_open(trim(db_dir)//"/"//trim(filename),nf90_nowrite,ncid)

    !Read dimension
    status=nf90_inq_dimid(ncid,"traj",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = ndb)

    status=nf90_inq_dimid(ncid,"obs",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = nobs)

    if(ndb == 0 .or. nobs == 0)then
       status=nf90_close(ncid)
       return
    end if       
    
    !Allocate
    allocate(time(nobs,ndb))
    allocate(iyr(nobs,ndb),imon(nobs,ndb),iday(nobs,ndb))
    allocate(lon(nobs,ndb),lat(nobs,ndb))
    allocate(u(nobs,ndb),v(nobs,ndb),t(nobs,ndb))

    !Read data
    status=nf90_inq_varid(ncid,"lon360",varid)
    status=nf90_get_var(ncid,varid,lon)

    status=nf90_inq_varid(ncid,"latitude",varid)
    status=nf90_get_var(ncid,varid,lat)

    status=nf90_inq_varid(ncid,"time",varid)
    status=nf90_get_var(ncid,varid,time)

    status=nf90_inq_varid(ncid,"temp",varid)
    status=nf90_get_var(ncid,varid,t)

    status=nf90_inq_varid(ncid,"ve",varid)
    status=nf90_get_var(ncid,varid,u)

    status=nf90_inq_varid(ncid,"vn",varid)
    status=nf90_get_var(ncid,varid,v)
    
    !Close
    status=nf90_close(ncid)

    call ymd_julian(1970,1,1,sjul) !Reference date: 1970/01/01    
    
    do idb=1,ndb       
       do iobs=1,nobs

          if(time(iobs,idb) == dmiss .or. lon(iobs,idb) == dmiss .or. lat(iobs,idb) == dmiss)then
             iyr(iobs,idb)=int(rmiss)
             imon(iobs,idb)=int(rmiss)
             iday(iobs,idb)=int(rmiss)
             lon(iobs,idb)=rmiss
             lat(iobs,idb)=rmiss
             t(iobs,idb)=rmiss
             u(iobs,idb)=rmiss
             v(iobs,idb)=rmiss
          else
             ijul=sjul+int(time(iobs,idb)/(24.d0*60.d0*60.d0)) !julian day
             call julian_ymd(ijul,iyr(iobs,idb),imon(iobs,idb),iday(iobs,idb))
             ihour=nint((time(iobs,idb)-dble(ijul-sjul)*24.d0*60.d0*60.d0)/3600.d0)
          end if          

          if(t(iobs,idb) == dmiss)then
             t(iobs,idb)=rmiss
          end if

          if(u(iobs,idb) == dmiss)then
             u(iobs,idb)=rmiss
          end if

          if(v(iobs,idb) == dmiss)then
             v(iobs,idb)=rmiss
          end if
          
       end do !iobs       
    end do !idb
    
    !---Deallocate
    deallocate(time)
    
  end subroutine read_db

  !---------------------------
  
  subroutine read_db_1d(filename,ndb,nobs_1d,iyr_1d,imon_1d,iday_1d,lon_1d,lat_1d,u_1d,v_1d,t_1d)

    use netcdf
    use mod_rmiss
    use mod_julian
    implicit none

    !---Parameter
    real(kind = 8),parameter :: dmiss=-999999.d0
    real(kind = 8),parameter :: pi=4.d0*atan(1.d0)
    
    !---Common
    integer status,access
    integer ncid,dimid,varid
    integer idb
    integer iobs,nobs
    integer iobs_1d
    integer sjul,ijul
    integer ihour
    integer i
    
    real(kind = 8),allocatable :: time(:,:)
    real(kind = 8),allocatable :: lon(:,:),lat(:,:)
    real(kind = 8),allocatable :: u(:,:),v(:,:),t(:,:)

    integer lonpass,lonmiss
    integer latpass,latmiss
    integer upass,umiss
    integer vpass,vmiss
    integer tpass,tmiss
    
    real(kind = 8) lonx,lony
    
    !---IN
    character(100),intent(in) :: filename

    !---OUT
    integer,intent(out) :: ndb,nobs_1d
    integer,intent(out),allocatable :: iyr_1d(:,:),imon_1d(:,:),iday_1d(:,:)

    real(kind = 8),intent(out),allocatable :: lon_1d(:,:),lat_1d(:,:)
    real(kind = 8),intent(out),allocatable :: u_1d(:,:),v_1d(:,:),t_1d(:,:)

    !Access
    status=access(trim(db_dir)//"/"//trim(filename)," ")
    if(status == 0)then
       !write(*,*) "Read:"//trim(db_dir)//"/"//trim(filename)
    else
       write(*,*) "***Error: "//trim(db_dir)//"/"//trim(filename)//"not found"
       stop
    end if
        
    !Open
    status=nf90_open(trim(db_dir)//"/"//trim(filename),nf90_nowrite,ncid)

    !Read dimension
    status=nf90_inq_dimid(ncid,"traj",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = ndb)

    status=nf90_inq_dimid(ncid,"obs",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = nobs)

    if(ndb == 0 .or. nobs == 0)then
       nobs_1d=0
       status=nf90_close(ncid)
       return
    end if
    
    !Allocate
    allocate(time(nobs,ndb))
    allocate(lon(nobs,ndb),lat(nobs,ndb))
    allocate(u(nobs,ndb),v(nobs,ndb),t(nobs,ndb))

    !Read data
    status=nf90_inq_varid(ncid,"lon360",varid)
    status=nf90_get_var(ncid,varid,lon)

    status=nf90_inq_varid(ncid,"latitude",varid)
    status=nf90_get_var(ncid,varid,lat)

    status=nf90_inq_varid(ncid,"time",varid)
    status=nf90_get_var(ncid,varid,time)

    status=nf90_inq_varid(ncid,"temp",varid)
    status=nf90_get_var(ncid,varid,t)

    status=nf90_inq_varid(ncid,"ve",varid)
    status=nf90_get_var(ncid,varid,u)

    status=nf90_inq_varid(ncid,"vn",varid)
    status=nf90_get_var(ncid,varid,v)
    
    !Close
    status=nf90_close(ncid)

    call ymd_julian(1970,1,1,sjul) !Reference date: 1970/01/01    
    
    !---Count nobs_1d    
    nobs_1d=0

    do idb=1,ndb

       if(nobs < 4)then
          cycle
       end if

       i=0
       iobs_1d=0
       
       do iobs=1,nobs

          if(time(iobs,idb) == dmiss)then
             cycle
          else
             ijul=sjul+int(time(iobs,idb)/(24.d0*60.d0*60.d0)) !julian day
             !call julian_ymd(ijul,iyr,imon,iday)
             ihour=nint((time(iobs,idb)-dble(ijul-sjul)*24.d0*60.d0*60.d0)/3600.d0)
          end if

          if(i*6 == ihour)then
             i=i+1
          else
             i=0
             cycle
          end if

          if(i == 4)then
             iobs_1d=iobs_1d+1
             i=0
          end if

       end do

       if(iobs_1d > nobs_1d)then
          nobs_1d=iobs_1d
       end if
       
    end do
    
    if(nobs_1d == 0)then
       deallocate(time)
       deallocate(lon,lat)
       deallocate(u,v,t)
       return
    end if

    !--- 1-day average    
    allocate(iyr_1d(nobs_1d,ndb),imon_1d(nobs_1d,ndb),iday_1d(nobs_1d,ndb))
    allocate(lon_1d(nobs_1d,ndb),lat_1d(nobs_1d,ndb))
    allocate(u_1d(nobs_1d,ndb),v_1d(nobs_1d,ndb),t_1d(nobs_1d,ndb))

    do idb=1,ndb

       i=0
       iobs_1d=1

       do iobs=1,nobs

          if(nobs_1d < iobs_1d) exit 
          
          !Initialization
          if(i == 0)then
          
             lonx=0.d0
             lony=0.d0
             lonpass=0
             lonmiss=0
          
             lat_1d(iobs_1d,idb)=0.d0
             latpass=0
             latmiss=0
          
             u_1d(iobs_1d,idb)=0.d0
             upass=0
             umiss=0
             
             v_1d(iobs_1d,idb)=0.d0
             vpass=0
             vmiss=0
             
             t_1d(iobs_1d,idb)=0.d0
             tpass=0
             tmiss=0

          end if

          !Time
          if(lon(iobs,idb) == dmiss .or. lat(iobs,idb) == dmiss .or. time(iobs,idb) == dmiss)then
             cycle
          else
             ijul=sjul+int(time(iobs,idb)/(24.d0*60.d0*60.d0)) !julian day
             call julian_ymd(ijul,iyr_1d(iobs_1d,idb),imon_1d(iobs_1d,idb),iday_1d(iobs_1d,idb))
             ihour=int((time(iobs,idb)-dble(ijul-sjul)*24.d0*60.d0*60.d0)/3600.d0)
          end if

          if(i*6 == ihour)then
             i=i+1
          else
             i=0
             cycle
          end if
          
          !Add
          if(lon(iobs,idb) == dmiss)then
             lonmiss=lonmiss+1
          else
             lonpass=lonpass+1
             lonx=lonx+cos(lon(iobs,idb)*pi/180.d0)
             lony=lony+sin(lon(iobs,idb)*pi/180.d0)
          end if
             
          if(lat(iobs,idb) == dmiss)then
             latmiss=latmiss+1
          else
             latpass=latpass+1
             lat_1d(iobs_1d,idb)=lat_1d(iobs_1d,idb)+lat(iobs,idb)
          end if
          
          if(t(iobs,idb) == dmiss)then
             tmiss=tmiss+1
          else
             tpass=tpass+1
             t_1d(iobs_1d,idb)=t_1d(iobs_1d,idb)+t(iobs,idb)
          end if
          
          if(u(iobs,idb) == dmiss)then
             umiss=umiss+1
          else
             upass=upass+1
             u_1d(iobs_1d,idb)=u_1d(iobs_1d,idb)+u(iobs,idb)
          end if
          
          if(v(iobs,idb) == dmiss)then
             vmiss=vmiss+1
          else
             vpass=vpass+1
             v_1d(iobs_1d,idb)=v_1d(iobs_1d,idb)+v(iobs,idb)
          end if
          
          !Daily mean
          if(i == 4)then

             if(lonmiss == 0)then
                lon_1d(iobs_1d,idb)=atan2(lony,lonx)*180.d0/pi
                if(lon_1d(iobs_1d,idb) < 0.d0) lon_1d(iobs_1d,idb)=lon_1d(iobs_1d,idb)+360.d0
             else
                lon_1d(iobs_1d,idb)=rmiss
             end if
             
             if(latmiss == 0)then
                lat_1d(iobs_1d,idb)=lat_1d(iobs_1d,idb)/dble(latpass)
             else
                lat_1d(iobs_1d,idb)=rmiss
             end if
             
             if(tmiss == 0)then
                t_1d(iobs_1d,idb)=t_1d(iobs_1d,idb)/dble(tpass)
             else
                t_1d(iobs_1d,idb)=rmiss
             end if
             
             if(umiss == 0)then
                u_1d(iobs_1d,idb)=u_1d(iobs_1d,idb)/dble(upass)
             else
                u_1d(iobs_1d,idb)=rmiss
             end if
             
             if(vmiss == 0)then
                v_1d(iobs_1d,idb)=v_1d(iobs_1d,idb)/dble(vpass)
             else
                v_1d(iobs_1d,idb)=rmiss
             end if
             
             !Reset
             i=0
             iobs_1d=iobs_1d+1
             
          end if
          
       end do !iobs
    end do !idb
    
    !---Deallocate
    deallocate(time)
    deallocate(lon,lat)
    deallocate(u,v,t)
    
  end subroutine read_db_1d

  !------------------------
  
  subroutine deallocate_db(iyr,imon,iday,lon,lat,u,v,t)

    implicit none
    
    integer,intent(inout),allocatable :: iyr(:,:),imon(:,:),iday(:,:)
    real(kind = 8),intent(inout),allocatable :: lon(:,:),lat(:,:)
    real(kind = 8),intent(inout),allocatable :: u(:,:),v(:,:),t(:,:)

    if(allocated(iyr)) deallocate(iyr)
    if(allocated(imon)) deallocate(imon)
    if(allocated(iday)) deallocate(iday)
    if(allocated(lon)) deallocate(lon)
    if(allocated(lat)) deallocate(lat)
    if(allocated(u)) deallocate(u)
    if(allocated(v)) deallocate(v)
    if(allocated(t)) deallocate(t)
    
  end subroutine deallocate_db
  
end module mod_read_db
