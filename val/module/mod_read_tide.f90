module mod_read_tide

  integer,parameter :: nst=915 !Number of station
  character(100),parameter :: tide_dir="/data/R/R2402/DATA/TIDE/daily"
  
contains

  !-------------------------------------------------------------
  ! Read sea level data from tide gauge |
  !-------------------------------------------------------------
  !
  ! Provider: University of Hawaii Sea Level Center (UHSLC)
  ! Website:  https://uhslc.soest.hawaii.edu/
  ! 
  ! - Sea level unit: [mm] --> [m]
  ! - Data might not have continous time series.
  ! - If data have overlap period, the prior is the primary, and the posterir bias is offset. 
  !
  !------------------------------------------------------------
  ! Created by S.Ohishi 2025.12.03
  !
  !-------------------------------------------------------------

  subroutine get_tide_data(filename,ntime,ijul,lon,lat,dat)

    use mod_julian
    use mod_rmiss
    use netcdf
    implicit none

    !---Parameter
    integer,parameter :: dmiss=-32767
    
    !---Common
    integer status,ncid,dimid,varid
    integer ijul0
    integer itime
    integer,allocatable :: tmp1d(:)

    real(kind = 4) tmp1d_x(1),tmp1d_y(1)
    real(kind = 8),allocatable :: tmp1d_t(:)
    
    !---IN
    character(200),intent(in) :: filename
    
    !---OUT
    integer,intent(out) :: ntime
    integer,allocatable,intent(out) :: ijul(:)

    real(kind = 8),allocatable,intent(out) :: lon(:),lat(:),dat(:)

    !---Reference julian date
    call ymd_julian(1800,1,1,ijul0)
    
    !---Open
    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    !---Get dimension
    status=nf90_inq_dimid(ncid,"time",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = ntime)

    !---Check ntime (Just in case)
    if(ntime == 0)then
       return
    end if
    
    !---Allocate
    allocate(tmp1d(ntime),tmp1d_t(ntime))
    
    allocate(ijul(ntime))
    allocate(lon(ntime),lat(ntime),dat(ntime))

    !---Read    
    status=nf90_inq_varid(ncid,"time",varid)
    status=nf90_get_var(ncid,varid,tmp1d_t)
    
    status=nf90_inq_varid(ncid,"lon",varid)
    status=nf90_get_var(ncid,varid,tmp1d_x)
    
    status=nf90_inq_varid(ncid,"lat",varid)
    status=nf90_get_var(ncid,varid,tmp1d_y)

    status=nf90_inq_varid(ncid,"sea_level",varid)
    status=nf90_get_var(ncid,varid,tmp1d)

    !---Close
    status=nf90_close(ncid)

    !---Post process
    ijul(:)=int(tmp1d_t(:))+ijul0
    lon(:)=dble(tmp1d_x(1))
    lat(:)=dble(tmp1d_y(1))
    do itime=1,ntime
       if(tmp1d(itime) == dmiss)then
          dat(itime)=rmiss
       else
          dat(itime)=1.d-3*dble(tmp1d(itime)) !Unit[mm] -> [m]
       end if
    end do

    deallocate(tmp1d,tmp1d_t)
    
  end subroutine get_tide_data

  !-----------------------------

  subroutine end_get_tide_data(ijul,lon,lat,dat)

    implicit none

    integer,allocatable,intent(out) :: ijul(:)
    real(kind = 8),allocatable,intent(out) :: lon(:),lat(:),dat(:)

    if(allocated(ijul)) deallocate(ijul)
    if(allocated(lon)) deallocate(lon)
    if(allocated(lat)) deallocate(lat)
    if(allocated(dat)) deallocate(dat)

  end subroutine end_get_tide_data

  !--------------------------------------------------------------

  subroutine remove_ave(ntime,dat)

    use mod_rmiss
    implicit none

    !---Common
    integer itime
    integer n

    real(kind = 8) ave
    
    !---IN
    integer,intent(in) :: ntime

    !---IN/OUT
    real(kind = 8),intent(inout) :: dat(ntime)

    !---Average
    n=0
    ave=0.d0

    do itime=1,ntime

       if(dat(itime) == rmiss) cycle
       
       ave=ave+dat(itime)
       n=n+1
       
    end do

    if(n == 0)then
       ave=rmiss
       return
    else
       ave=ave/dble(n)
    end if
    
    !---Remove average
    do itime=1,ntime

       if(dat(itime) == rmiss) cycle

       dat(itime)=dat(itime)-ave

    end do
    
  end subroutine remove_ave
  
  !--------------------------------------------------------------

  subroutine get_overlap_period(ntime1,ijul1,ntime2,ijul2, &
       & stime1,etime1,stime2,etime2)

    implicit none

    !---Common
    integer n
    integer itime1,itime2
    
    !---IN
    integer,intent(in) :: ntime1,ijul1(ntime1)
    integer,intent(in) :: ntime2,ijul2(ntime2)

    !---OUT
    integer,intent(out) :: stime1,etime1
    integer,intent(out) :: stime2,etime2

    n=0
    stime1=0
    etime1=0
    stime2=0
    etime2=0

    if(ijul1(ntime1) < ijul2(1))then
       return
    end if
    
    do itime1=1,ntime1

       if(ijul1(itime1) < ijul2(1)) cycle
       
       do itime2=1,ntime2
          
          if(ijul1(itime1) == ijul2(itime2))then

             n=n+1

             if(n == 1)then
                stime1=itime1
                stime2=itime2
             end if

             etime1=itime1
             etime2=itime2

          end if

       end do
    end do          
    
  end subroutine get_overlap_period
  
  !--------------------------------------------------------------

  subroutine get_bias(ntime,dat1,dat2,bias)

    use mod_rmiss
    implicit none

    !---Common
    integer itime
    integer n
    
    !---IN
    integer,intent(in) :: ntime
    real(kind = 8),intent(in) :: dat1(ntime),dat2(ntime)

    !---OUT
    real(kind = 8),intent(out) :: bias
    
    n=0
    bias=0.d0

    do itime=1,ntime
       
       if(dat1(itime) == rmiss .or. dat2(itime) == rmiss) cycle

       n=n+1
       bias=bias+(dat2(itime)-dat1(itime))

    end do
    
    if(n == 0)then
       bias=0.d0
    else
       bias=bias/dble(n)
    end if    
    
  end subroutine get_bias
  
  !--------------------------------------------------------------
  
  subroutine read_tide(ist,ntime,ijul,lon,lat,dat)

    use mod_rmiss
    use netcdf
    implicit none

    !---Parameter
    integer,parameter :: nfile_max=5
    character(5),parameter :: abc="abcde"
    
    !---Common
    integer ifile,nfile
    integer itime
    integer status,access

    integer stime,etime
    integer stime1,etime1,ntime1
    integer stime2,etime2,ntime2
    integer,allocatable :: ijul1(:),ijul2(:)

    real(kind = 8) bias    
    real(kind = 8),allocatable :: lon1(:),lat1(:),dat1(:)
    real(kind = 8),allocatable :: lon2(:),lat2(:),dat2(:)
    
    character(200),allocatable :: filename_list(:)
    character(200) filename
    character(3) nnn
    
    !---IN
    integer,intent(in) :: ist

    !---OUT
    integer,intent(out) :: ntime
    integer,allocatable,intent(out) :: ijul(:)

    real(kind = 8),allocatable,intent(out) :: lon(:),lat(:),dat(:)

    !---File number
    write(nnn,'(i3.3)') ist
    
    !---Count nfile
    nfile=0
    do ifile=1,nfile_max
       
       filename=trim(tide_dir)//"/d"//nnn//abc(ifile:ifile)//".nc"
       status=access(trim(filename)," ")

       if(status == 0)then
          nfile=nfile+1
       end if
       
    end do

    if(nfile == 0)then
       ntime=0
       write(*,*) "No data at d"//trim(nnn)
       return
    end if    

    !---Make filename list
    allocate(filename_list(nfile))
       
    nfile=0
    do ifile=1,nfile_max

       filename=trim(tide_dir)//"/d"//nnn//abc(ifile:ifile)//".nc"
       status=access(trim(filename)," ")

       if(status == 0)then
          nfile=nfile+1
          filename_list(nfile)=filename
       end if

    end do

    !---Count ntime
    do ifile=1,nfile

       if(ifile == 1)then
          
          call get_tide_data(filename_list(ifile),ntime1,ijul1,lon1,lat1,dat1)
          ntime=ntime1
          call end_get_tide_data(ijul1,lon1,lat1,dat1)
          
       else
          
          call get_tide_data(filename_list(ifile-1),ntime1,ijul1,lon1,lat1,dat1)
          call get_tide_data(filename_list(ifile),ntime2,ijul2,lon2,lat2,dat2)

          !Without overlap
          if(ijul1(ntime1) < ijul2(1))then          
             ntime=ntime+ntime2
          else
             !With overlap
             call get_overlap_period(ntime1,ijul1,ntime2,ijul2,stime1,etime1,stime2,etime2)
             ntime=ntime+ntime2-etime2
          end if
          
          call end_get_tide_data(ijul1,lon1,lat1,dat1)
          call end_get_tide_data(ijul2,lon2,lat2,dat2)

       end if
          
    end do
    
    !---Allocate
    allocate(ijul(ntime))
    allocate(lon(ntime),lat(ntime),dat(ntime))
    
    !---Read tide time series
    do ifile=1,nfile
       
       if(ifile == 1)then
          
          call get_tide_data(filename_list(ifile),ntime1,ijul1,lon1,lat1,dat1)
          call remove_ave(ntime1,dat1)

          ntime=ntime1
          ijul(1:ntime1)=ijul1(1:ntime1)
          lon(1:ntime1)=lon1(1:ntime1)
          lat(1:ntime1)=lat1(1:ntime1)
          dat(1:ntime1)=dat1(1:ntime1)
          call end_get_tide_data(ijul1,lon1,lat1,dat1)
          
       else
          
          call get_tide_data(filename_list(ifile),ntime1,ijul1,lon1,lat1,dat1)
          call remove_ave(ntime1,dat1)
          
          !Without overlap
          if(ijul(ntime) < ijul1(1))then
             
             !Average
             call remove_ave(ntime,dat(1:ntime))
             call remove_ave(ntime1,dat1)
             
             ijul(ntime+1:ntime+ntime1)=ijul1(1:ntime1)
             lon(ntime+1:ntime+ntime1)=lon1(1:ntime1)
             lat(ntime+1:ntime+ntime1)=lat1(1:ntime1)
             dat(ntime+1:ntime+ntime1)=dat1(1:ntime1)
             ntime=ntime+ntime1
                          
          else
             
             !With overlap             
             call get_overlap_period(ntime,ijul,ntime1,ijul1,stime,etime,stime1,etime1)
             call get_bias(etime-stime+1,dat(stime:etime),dat1(stime1:etime1),bias)
             do itime=1,ntime1-etime1
                ijul(ntime+itime)=ijul1(etime1+itime)
                lon(ntime+itime)=lon1(etime1+itime)
                lat(ntime+itime)=lat1(etime1+itime)
                if(dat1(etime1+itime) == rmiss)then
                   dat(ntime+itime)=rmiss
                else
                   dat(ntime+itime)=dat1(etime1+itime)-bias
                end if
             end do
             ntime=ntime+ntime1-etime1

          end if
          
          call end_get_tide_data(ijul1,lon1,lat1,dat1)

       end if
       
    end do

    !---Deallicate
    deallocate(filename_list)
                
  end subroutine read_tide

  !--------------------------------------------------------------

  subroutine end_read_tide(ijul,lon,lat,dat)

    implicit none

    integer,allocatable,intent(out) :: ijul(:)
    real(kind = 8),allocatable,intent(out) :: lon(:),lat(:),dat(:)

    if(allocated(ijul)) deallocate(ijul)
    if(allocated(lon)) deallocate(lon)
    if(allocated(lat)) deallocate(lat)
    if(allocated(dat)) deallocate(dat)

  end subroutine end_read_tide

end module mod_read_tide

  
