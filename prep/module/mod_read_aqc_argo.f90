module mod_read_aqc_argo

contains

  !---------------------------------------------------------------
  ! Read AQC Argo |
  !----------------
  !
  ! Modified by S.Ohishi 2022.06.20
  !
  !    ~2019.12: AQC Argo
  !    2020.01~: AQC + DQC Argo (Salinity bias corrected)
  !
  ! Modified by S. Ohishi 2022.09.16
  !
  !    ~2022.07: AQC Argo
  !    Add no AQC Argo case
  !
  ! Modified by S. Ohishi 2025.07.08
  !
  !    Add detailed comment
  !
  !----------------------------------------------------------------

  subroutine read_aqc_argo(iyr,imon,iday,npro,km,long,lati,depth,t,s)

    use mod_rmiss
    use netcdf
    implicit none

    !---Parameter
    real(kind = 4),parameter :: dmiss1=9999.9e0,dmiss2=99.999e0

    !---Common
    integer status,access
    integer ncid,dimid,varid
    integer ipro,k

    real(kind = 4),allocatable :: rlon(:),rlat(:),rdep(:,:)
    real(kind = 4),allocatable :: rt(:,:),rs(:,:)
    
    character(1),allocatable,dimension(:,:) :: dflag,tflag,sflag
    character(1) qc
    character(2) mm
    character(4) yyyy
    character(16),allocatable,dimension(:) :: yyyymmdd
    character(16) ctmp
    character(100) filename

    !---IN
    integer,intent(in) :: iyr,imon

    !---OUT
    integer,allocatable,intent(out) :: iday(:)
    integer,intent(out) :: npro
    integer,intent(out) :: km

    real(kind = 8),allocatable,intent(out) :: long(:),lati(:),depth(:,:)
    real(kind = 8),allocatable,intent(out) :: t(:,:),s(:,:)
    
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon

    !---Filename
    filename="/data/R/R2402/DATA/AQC_Argo/AQC_Profile_Data_" &
         & //yyyy//mm//".nc"
    status=access(trim(filename)," ")

    !---Access
    if(status /= 0)then
       npro=0
       km=0
       write(*,'(a)') "***Error: Not found "//trim(filename)
       return
    end if

    !---Open netcdf file
    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    !---Get km,npro
    status=nf90_inq_dimid(ncid,"N_LEVELS",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = km)

    status=nf90_inq_dimid(ncid,"N_PROF",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = npro)

    write(*,'(a,i10)') "The number of profile: ",npro
    write(*,'(a,i10)') "The number of depth: ",km

    !---Allocate
    allocate(yyyymmdd(npro),iday(npro))
    
    allocate(rlon(npro),rlat(npro),rdep(km,npro))
    allocate(rt(km,npro),rs(km,npro))
    allocate(long(npro),lati(npro),depth(km,npro))
    allocate(t(km,npro),s(km,npro))
    allocate(dflag(km,npro),tflag(km,npro),sflag(km,npro))

    !---Read Data
    status=nf90_inq_varid(ncid,"TIME",varid)
    status=nf90_get_var(ncid,varid,yyyymmdd)

    status=nf90_inq_varid(ncid,"LONGITUDE",varid)
    status=nf90_get_var(ncid,varid,rlon)

    status=nf90_inq_varid(ncid,"LATITUDE",varid)
    status=nf90_get_var(ncid,varid,rlat)

    status=nf90_inq_varid(ncid,"PRES",varid)
    status=nf90_get_var(ncid,varid,rdep)

    status=nf90_inq_varid(ncid,"TEMP",varid)
    status=nf90_get_var(ncid,varid,rt)

    status=nf90_inq_varid(ncid,"PSAL",varid)
    status=nf90_get_var(ncid,varid,rs)

    status=nf90_inq_varid(ncid,"PRES_FLAG",varid)
    status=nf90_get_var(ncid,varid,dflag)

    status=nf90_inq_varid(ncid,"TEMP_FLAG",varid)
    status=nf90_get_var(ncid,varid,tflag)

    status=nf90_inq_varid(ncid,"PSAL_FLAG",varid)
    status=nf90_get_var(ncid,varid,sflag)

    status=nf90_close(ncid)

    !---Post process
    !iday
    do ipro=1,npro

       ctmp=yyyymmdd(ipro)
       read(ctmp(7:8),*) iday(ipro)
       
       if(rlon(ipro) < 0.e0)then
          long(ipro)=dble(rlon(ipro))+360.d0
       else
          long(ipro)=dble(rlon(ipro))
       end if
              
    end do !ipro

    lati(:)=dble(rlat(:))
    depth(:,:)=dble(rdep(:,:))
    t(:,:)=dble(rt(:,:))
    s(:,:)=dble(rs(:,:))
    
    !---QC
    qc="1"
    do ipro=1,npro
       do k=1,km
          if(rdep(k,ipro)==dmiss1 .or. dflag(k,ipro) /= qc) depth(k,ipro)=rmiss
          if(rt(k,ipro)==dmiss2 .or. tflag(k,ipro) /= qc) t(k,ipro)=rmiss
          if(rs(k,ipro)==dmiss2 .or. sflag(k,ipro) /= qc) s(k,ipro)=rmiss          
       end do !k
    end do    !ipro

    deallocate(yyyymmdd)    
    deallocate(rlon,rlat,rdep)
    deallocate(rt,rs)
    deallocate(dflag,tflag,sflag)

  end subroutine read_aqc_argo

  !------------------------------

  subroutine deallocate_aqc_argo(iday,long,lati,depth,t,s)

    implicit none
    integer,allocatable :: iday(:)
    real(kind = 8),allocatable :: long(:),lati(:),depth(:,:)
    real(kind = 8),allocatable :: t(:,:),s(:,:)
    
    deallocate(iday)
    deallocate(long,lati,depth)
    deallocate(t,s)

  end subroutine deallocate_aqc_argo

end module mod_read_aqc_argo
