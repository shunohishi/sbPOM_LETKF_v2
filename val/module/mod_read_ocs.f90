module mod_read_ocs

  character(100),parameter :: keo_dir="/data/R/R2402/DATA/KEO"
  character(100),parameter :: papa_dir="/data/R/R2402/DATA/PAPA"

contains

  !-------------------------------------------------------------------
  ! Read Ocean Cite Station (KEO/PAPA) data |
  !-------------------------------------------------------------------
  !
  !-------------------------------------------------------------------

  subroutine get_name(var,syr,smon,sday,filename,datname,qcname)

    implicit none

    integer syr,smon,sday

    character(10) var,datname,qcname
    character(100) filename

    if(var == "t")then
       syr=2004
       smon=6
       sday=16
       datname="T_20" !degree C
       qcname="QT_5020"
       filename="/data/R/R2402/DATA/KEO/t32n145e_dy.cdf"
    else if(var == "s")then
       syr=2004
       smon=6
       sday=16
       datname="S_41" !psu
       qcname="QS_5041"
       filename="/data/R/R2402/DATA/KEO/s32n145e_dy.cdf"
    else if(var == "u")then
       syr=2005
       smon=5
       sday=30
       datname="U_320" !cm/s
       qcname="QCS_5300"
       filename="/data/R/R2402/DATA/KEO/cur32n145e_dy.cdf"      
    else if(var == "v")then
       syr=2005
       smon=5
       sday=30
       datname="V_321" !cm/s
       qcname="QCS_5300"
       filename="/data/R/R2402/DATA/KEO/cur32n145e_dy.cdf"         
    else if(var == "speed")then
       syr=2005
       smon=5
       sday=30
       datname="CS_300" !cm/s
       qcname="QCS_5300"
       filename="/data/R/R2402/DATA/KEO/cur32n145e_dy.cdf"
    else if(var == "theta")then
       syr=2005
       smon=5
       sday=30
       datname="CD_310" !degree
       qcname="QCD_5310"
       filename="/data/R/R2402/DATA/KEO/cur32n145e_dy.cdf"
    else
       write(*,*) "***Error: Not found "//var
       stop
    end if

  end subroutine get_name

  !------------------

  subroutine read_keo(var,ntime,ndep,iyr,imon,iday,long,lati,depth,dat)

    use julian_day
    use mod_rmiss
    implicit none
    include "netcdf.inc"

    real,parameter :: dmiss=1.e35

    integer status,access
    integer ncid,dimid,varid
    integer itime,ntime
    integer idep,ndep
    integer syr,smon,sday
    integer ijul,sjul

    integer,allocatable ::  iyr(:),imon(:),iday(:)
    real long,lati
    real,allocatable :: time(:),depth(:),dat(:,:),qc(:,:)

    character(10) var,datname,qcname
    character(100) filename

    call get_name(var,syr,smon,sday,filename,datname,qcname)

    status=access(trim(filename)," ")
    if(status == 0)then
       write(*,*) "Read :",trim(filename)
    else
       write(*,*) "Not found: "//trim(filename)
       stop
    end if

    !---open netcdf
    status=nf_open(trim(filename),nf_nowrite,ncid)

    !---get dimension
    !ntime
    status=nf_inq_dimid(ncid,"time",dimid)
    status=nf_inq_dimlen(ncid,dimid,ntime)

    !ndep
    status=nf_inq_dimid(ncid,"depth",dimid)
    status=nf_inq_dimlen(ncid,dimid,ndep)

    !allocate
    allocate(iyr(ntime),imon(ntime),iday(ntime))
    allocate(depth(ndep),time(ntime))
    allocate(dat(ndep,ntime),qc(ndep,ntime))

    !---read data
    !time
    status=nf_inq_varid(ncid,"time",varid)
    status=nf_get_vara_real(ncid,varid,(/1/),(/ntime/),time)

    !long
    status=nf_inq_varid(ncid,"lon",varid)
    status=nf_get_vara_real(ncid,varid,(/1/),(/1/),long)

    !lati
    status=nf_inq_varid(ncid,"lat",varid)
    status=nf_get_vara_real(ncid,varid,(/1/),(/1/),lati)

    !depth
    status=nf_inq_varid(ncid,"depth",varid)
    status=nf_get_vara_real(ncid,varid,(/1/),(/ndep/),depth)

    !dat
    status=nf_inq_varid(ncid,trim(datname),varid)
    status=nf_get_vara_real(ncid,varid,(/1,1,1,1/),(/1,1,ndep,ntime/),dat)

    !qc
    status=nf_inq_varid(ncid,trim(qcname),varid)
    status=nf_get_vara_real(ncid,varid,(/1,1,1,1/),(/1,1,ndep,ntime/),qc)

    !---close
    status=nf_close(ncid)

    !---check QC
    do itime=1,ntime
       do idep=1,ndep

          !See https://www.pmel.noaa.gov/ocs/oceansites-data-info
          ! 1: Good data
          if(qc(idep,itime) /= 1 .and. qc(idep,itime) /= 2)then
             dat(idep,itime)=rmiss
          end if

          if(dat(idep,itime) == dmiss) dat(idep,itime)=rmiss

       end do
    end do

    !Modify unit
    if(var == "u" .or. var == "v" .or. var == "speed")then
       do itime=1,ntime
          do idep=1,ndep
             if(dat(idep,itime) /= rmiss) dat(idep,itime)=0.01*dat(idep,itime) !cm/s -> m/s
          end do
       end do
    end if

    !---Estimate date
    call ymd_julian(syr,smon,sday,sjul)
    do itime=1,ntime
       ijul=sjul+int(time(itime))
       call julian_ymd(ijul,iyr(itime),imon(itime),iday(itime))
    end do

    deallocate(time,qc)

  end subroutine read_keo

  !----------------------------------------------------------------------------------
  ! End Read KEO |
  !-----------------------------------------------------------------------------------

  subroutine end_read_keo(iyr,imon,iday,depth,dat)

    implicit none
    integer,allocatable :: iyr(:),imon(:),iday(:)
    real,allocatable :: depth(:),dat(:,:)

    deallocate(iyr,imon,iday)
    deallocate(depth,dat)

  end subroutine end_read_keo

end module mod_read_ocs

