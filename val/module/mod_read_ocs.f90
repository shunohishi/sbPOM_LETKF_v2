module mod_read_ocs

  character(100),parameter :: keo_dir="/data/R/R2402/DATA/KEO"
  character(100),parameter :: papa_dir="/data/R/R2402/DATA/PAPA"

contains

  !-------------------------------------------------------------------
  ! Read Ocean Climate Station (KEO/PAPA) data |
  !-------------------------------------------------------------------
  !
  ! Website: https://www.pmel.noaa.gov/ocs/
  ! 
  !-------------------------------------------------------------------

  subroutine get_ocs_info(buoyname,varname,syr,smon,sday,filename,datname,qcname)

    implicit none

    !---IN
    character(10),intent(in) :: buoyname
    character(1),intent(in) :: varname

    !---OUT
    integer,intent(out) :: syr,smon,sday

    character(10),intent(out) :: datname,qcname
    character(100),intent(out) :: filename

    if(buoyname == "keo" .and. varname == "t")then
       syr=2004
       smon=6
       sday=16
       datname="T_20" !degree C
       qcname="QT_5020"
       filename="t32n145e_dy.cdf"
    else if(buoyname == "keo" .and. varname == "s")then
       syr=2004
       smon=6
       sday=16
       datname="S_41" !psu
       qcname="QS_5041"
       filename="s32n145e_dy.cdf"
    else if(buoyname == "keo" .and. varname == "u")then
       syr=2005
       smon=5
       sday=30
       datname="U_320" !cm/s
       qcname="QCS_5300"
       filename="cur32n145e_dy.cdf"      
    else if(buoyname == "keo" .and. varname == "v")then
       syr=2005
       smon=5
       sday=30
       datname="V_321" !cm/s
       qcname="QCS_5300"
       filename="cur32n145e_dy.cdf"         
    else if(buoyname == "papa" .and. varname == "t")then
       syr=2007
       smon=6
       sday=8
       datname="T_20" !degree C
       qcname="QT_5020"
       filename="t50n145w_dy.cdf"
    else if(buoyname == "papa" .and. varname == "s")then
       syr=2007
       smon=6
       sday=8
       datname="S_41" !psu
       qcname="QS_5041"
       filename="s50n145w_dy.cdf"
    else if(buoyname == "papa" .and. varname == "u")then
       syr=2007
       smon=6
       sday=8
       datname="U_320" !cm/s
       qcname="QCS_5300"
       filename="cur50n145w_dy.cdf"      
    else if(buoyname == "papa" .and. varname == "v")then
       syr=2007
       smon=6
       sday=8
       datname="V_321" !cm/s
       qcname="QCS_5300"
       filename="cur50n145w_dy.cdf"         
    else
       write(*,*) "***Error: Not found ==> "//trim(buoyname)//trim(varname)
       stop
    end if

  end subroutine get_ocs_info

  !------------------

  subroutine read_ocs(buoyname,varname,ntime,ndep,iyr,imon,iday,long,lati,depth,dat)

    use mod_julian
    use mod_rmiss
    use netcdf
    implicit none

    !---Parameter
    real(kind = 4),parameter :: dmiss=1.e35

    !---Common
    integer status,access
    integer ncid,dimid,varid
    integer itime
    integer idep
    integer syr,smon,sday
    integer ijul,sjul

    real(kind = 4),allocatable :: time(:),qc(:,:)
    real(kind = 4),allocatable :: tmp1dx(:),tmp1dy(:),tmp1dz(:),tmp2d(:,:)
    
    character(10) datname,qcname
    character(100) filename
    character(200) fulfilename

    !---IN
    character(10),intent(in) :: buoyname
    character(1),intent(in) :: varname
    
    !---OUT
    integer,intent(out) :: ntime,ndep
    integer,allocatable,intent(out) ::  iyr(:),imon(:),iday(:)
    real(kind = 8),intent(out) :: long,lati
    real(kind = 8),allocatable,intent(out) :: depth(:),dat(:,:)

    !---Get info
    call get_ocs_info(buoyname,varname,syr,smon,sday,filename,datname,qcname)

    !---Filename
    if(buoyname == "keo")then
       fulfilename=trim(keo_dir)//"/"//trim(filename)
    else if(buoyname == "papa")then
       fulfilename=trim(papa_dir)//"/"//trim(filename)
    else
       write(*,*) "***Error: Incorrect buoyname ==> "//trim(buoyname)
    end if       
    
    status=access(trim(fulfilename)," ")
    if(status == 0)then
       write(*,*) "Read :",trim(fulfilename)
    else
       write(*,*) "Not found: "//trim(fulfilename)
       stop
    end if

    !---Read NetCDF file
    status=nf90_open(trim(fulfilename),nf90_nowrite,ncid)

    !---Get dimension
    !ntime
    status=nf90_inq_dimid(ncid,"time",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = ntime)

    !ndep
    status=nf90_inq_dimid(ncid,"depth",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = ndep)

    !---Allocate
    allocate(time(ntime),qc(ndep,ntime))
    allocate(tmp1dx(1),tmp1dy(1),tmp1dz(ndep),tmp2d(ndep,ntime))
    
    allocate(iyr(ntime),imon(ntime),iday(ntime))
    allocate(depth(ndep),dat(ndep,ntime))

    !---Read data
    !time
    status=nf90_inq_varid(ncid,"time",varid)
    status=nf90_get_var(ncid,varid,time)

    !long
    status=nf90_inq_varid(ncid,"lon",varid)
    status=nf90_get_var(ncid,varid,tmp1dx)

    !lati
    status=nf90_inq_varid(ncid,"lat",varid)
    status=nf90_get_var(ncid,varid,tmp1dy)

    !depth
    status=nf90_inq_varid(ncid,"depth",varid)
    status=nf90_get_var(ncid,varid,tmp1dz)

    !dat
    status=nf90_inq_varid(ncid,trim(datname),varid)
    status=nf90_get_var(ncid,varid,tmp2d,(/1,1,1,1/),(/1,1,ndep,ntime/))

    !qc
    status=nf90_inq_varid(ncid,trim(qcname),varid)
    status=nf90_get_var(ncid,varid,qc,(/1,1,1,1/),(/1,1,ndep,ntime/))

    status=nf90_close(ncid)

    !---Date
    call ymd_julian(syr,smon,sday,sjul)
    do itime=1,ntime
       ijul=sjul+int(time(itime))
       call julian_ymd(ijul,iyr(itime),imon(itime),iday(itime))
    end do

    !---Position
    long=dble(tmp1dx(1))
    lati=dble(tmp1dy(1))
    depth(:)=dble(tmp1dz(:))
    
    !---QC and Unit
    do itime=1,ntime
       do idep=1,ndep

          !See https://www.pmel.noaa.gov/ocs/oceansites-data-info
          ! 1: Good data, 2: Probably good data
          if(qc(idep,itime) == 1 .or. qc(idep,itime) == 2)then
             dat(idep,itime)=dble(tmp2d(idep,itime))
             !cm/s => m/s
             if(varname == "u" .or. varname == "v" .or. varname == "speed") dat(idep,itime)=0.01d0*dat(idep,itime) 
          else if(tmp2d(idep,itime) == dmiss)then
             dat(idep,itime)=rmiss
          else
             dat(idep,itime)=rmiss
          end if
             
       end do !idep
    end do    !itime

    deallocate(time,qc)
    deallocate(tmp1dx,tmp1dy,tmp1dz,tmp2d)

  end subroutine read_ocs

  !----------------------------------------------------------------------------------
  ! End Read KEO |
  !-----------------------------------------------------------------------------------

  subroutine deallocate_ocs(iyr,imon,iday,depth,dat)

    implicit none
    
    integer,allocatable,intent(out) :: iyr(:),imon(:),iday(:)
    real(kind = 8),allocatable,intent(out) :: depth(:),dat(:,:)

    if(allocated(iyr))   deallocate(iyr)
    if(allocated(imon))  deallocate(imon)
    if(allocated(iday))  deallocate(iday)
    if(allocated(depth)) deallocate(depth)
    if(allocated(dat))   deallocate(dat)

  end subroutine deallocate_ocs

end module mod_read_ocs

