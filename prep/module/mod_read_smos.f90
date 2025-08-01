module mod_read_smos

  !-----------------------------------------------------------------
  ! Read SMOS from CMEMS |
  !-----------------------------------------------------------------
  !
  ! Web: https://data.marine.copernicus.eu/product/MULTIOBS_GLO_PHY_SSS_L3_MYNRT_015_014/description
  !      https://doi.org/10.1016/j.rse.2016.02.061
  !
  !-----------------------------------------------------------------
  ! Created by S.Ohishi @ 2025.07.31 
  !-----------------------------------------------------------------

  integer,parameter :: im=1388,jm=584
  character(100),parameter :: pdir="/data/R/R2402/DATA"
  
contains

  subroutine read_smos(AD,iyr,imon,iday,lon,lat,dat)

    use mod_rmiss    
    use netcdf
    implicit none

    !---Parameter
    integer,parameter :: imiss=-32767
    
    !---Common
    integer i,j
    integer status,access
    integer ncid,varid

    integer iqc(im,jm)
    integer idat(im,jm)
    
    real(kind = 4) rlon(im),rlat(jm)
    real(kind = 4) rdat(im,jm)
    
    character(200) filename
    character(8) yyyymmdd
    character(4) yyyy
    character(2) mm,dd
    character(1) v
    
    !---IN
    integer,intent(in) :: iyr,imon,iday

    character(1) AD !A or D (Acending/Decending)
    
    !---OUT
    real(kind = 8),intent(out) :: lon(im),lat(jm)
    real(kind = 8),intent(out) :: dat(im,jm)

    !---Filename
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    yyyymmdd=yyyy//mm//dd

    !Version
    do i=3,1,-1

       write(v,'(i1.1)') i
       filename=trim(pdir)//"/SMOS/cmems/CATDS_CSF2Q"//trim(AD)//"_"&
            & //yyyymmdd//"T000000_"//yyyymmdd//"T235959_334_"//v//".nc"

       status=access(trim(filename)," ")
       if(status == 0)then
          exit
       end if
       
    end do

    !Check Access
    if(status == 0)then
       write(*,*) "Read: "//trim(filename)
    else
       write(*,*) "Skip to Read: "//trim(filename)
       lon(:)=rmiss
       lat(:)=rmiss
       dat(:,:)=rmiss
       return
    end if

    !---Read data
    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    status=nf90_inq_varid(ncid,"lon",varid)
    status=nf90_get_var(ncid,varid,rlon)

    status=nf90_inq_varid(ncid,"lat",varid)
    status=nf90_get_var(ncid,varid,rlat)

    status=nf90_inq_varid(ncid,"Sea_Surface_Salinity_Rain_Corrected",varid)
    status=nf90_get_var(ncid,varid,idat)

    status=nf90_inq_varid(ncid,"Sea_Surface_Salinity_QC",varid)
    status=nf90_get_var(ncid,varid,iqc)
    
    status=nf90_close(ncid)
    
    !---Post process
    !Longitude
    do i=1,im
       if(rlon(i) < 0.d0)then
          lon(i)=dble(rlon(i))+360.d0
       else
          lon(i)=dble(rlon(i))
       end if
    end do

    !Latitude
    lat(:)=dble(rlat(:))

    !SSS
    do j=1,jm
       do i=1,im
          if(idat(i,j) == imiss .or. iqc(i,j) /= 0)then
             dat(i,j)=rmiss
          else
             dat(i,j)=dble(idat(i,j))*0.001d0+30.d0
          end if          
       end do
    end do

  end subroutine read_smos
  
end module mod_read_smos
  
