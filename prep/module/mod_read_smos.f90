module mod_read_smos

  !Modified by S.Ohishi 2020.04
  !Modified by S.Ohishi 2025.02

  character(100),parameter :: smos_dir="/data/R/R2402/DATA/SMOS/"
  
contains

!----------------------------------------------------
! Filename |
!----------------------------------------------------

subroutine read_filename(iyr,imon,iday,nfile,filename)

  implicit none

  integer ifile
  integer status,system,access

  character(256) command
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

  !Check directory
  status=access(trim(smos_dir)//yyyy//mm," ")
  if(status /= 0)then
     write(*,*) "Not found: "//trim(smos_dir)//yyyy//mm
     nfile=0
     return
  end if
  
  !Extract filename
  command="find "//trim(smos_dir)//yyyy//mm//" -name SM_*_MIR_OSUDP2_"//yyyymmdd//"T*.nc | "//&
       &"sed 's|"//trim(smos_dir)//"||' > smos"//yyyymmdd//".txt"

  status=system(trim(command))
  
  !status=system("find "//trim(smos_dir)//yyyy//mm//&
  !&" -name SM_*_MIR_OSUDP2_"//yyyymmdd//"T*.nc "//&
  !     &"> smos"//yyyymmdd//".txt")

  status=access("smos"//yyyymmdd//".txt"," ")
  if(status /= 0)then
     write(*,*) "Not Found : smos"//yyyymmdd//".txt"
     nfile=0
     return
  end if
     
  !Count nfile
  nfile=0
  open(1,file="smos"//yyyymmdd//".txt",status="old")
  do 
     read(1,*,end=100)
     nfile=nfile+1
  end do
100 close(1)

  write(*,'(a,i10,a)') "The number of file:",nfile," at "//yyyymmdd

  if(nfile == 0)then
     status=system("rm -f smos"//yyyymmdd//".txt")
     return
  else

     !Read filename
     allocate(filename(nfile))

     open(1,file="smos"//yyyymmdd//".txt",status="old")
     do ifile=1,nfile
        read(1,'(a)') filename(ifile)     
     end do
     close(1)
     status=system("rm -f smos"//yyyymmdd//".txt")

  end if

end subroutine read_filename

!----------------------------------------------------
! Read SMAP |
!----------------------------------------------------

subroutine read_smos(filename,iyr,imon,iday,ihour,ngrid,lon,lat,dat)

  use mod_rmiss
  use netcdf
  implicit none

  !Common
  real(kind = 4),parameter :: dmiss=-999.

  integer igrid
  integer status,access,system
  integer ncid,dimid,varid

  integer,allocatable :: iqc(:)
  real(kind = 4),allocatable :: rlon(:),rlat(:),rdat(:)
  
  !IN
  character(200),intent(in) :: filename

  !OUT
  integer,intent(out) :: iyr,imon,iday,ihour
  integer,intent(out) :: ngrid
  real(kind = 8),allocatable,intent(out) :: lon(:),lat(:),dat(:)
  
  status=access(trim(smos_dir)//trim(filename)," ")
  if(status == 0)then
     write(*,*) "Read "//trim(smos_dir)//trim(filename)
  else
     write(*,*) "***Error: Not found "//trim(smos_dir)//trim(filename)
     stop
  end if
  
  read(filename(27:30),*) iyr
  read(filename(31:32),*) imon
  read(filename(33:34),*) iday
  read(filename(36:37),*) ihour

  status=nf90_open(trim(smos_dir)//trim(filename),nf90_nowrite,ncid)

  !Get im,jm
  status=nf90_inq_dimid(ncid,"n_grid_points",dimid)
  status=nf90_inquire_dimension(ncid,dimid,len = ngrid)

  !Allocate  
  allocate(rlon(ngrid),rlat(ngrid),rdat(ngrid),iqc(ngrid))
  allocate(lon(ngrid),lat(ngrid),dat(ngrid))

  !Read Data
  status=nf90_inq_varid(ncid,"Longitude",varid)
  status=nf90_get_var(ncid,varid,rlon)

  status=nf90_inq_varid(ncid,"Latitude",varid)
  status=nf90_get_var(ncid,varid,rlat)

  status=nf90_inq_varid(ncid,"Dg_quality_SSS_Corr",varid)
  status=nf90_get_var(ncid,varid,iqc)

  status=nf90_inq_varid(ncid,"SSS_corr",varid)
  status=nf90_get_var(ncid,varid,rdat)
  
  status=nf90_close(ncid)

  !Post process
  do igrid=1,ngrid
     
     if(rlon(igrid) == dmiss)then
        lon(igrid)=rmiss
     elseif(rlon(igrid) < 0.e0)then
        lon(igrid)=dble(rlon(igrid))+360.d0
     else
        lon(igrid)=dble(rlon(igrid))
     end if

     if(rlat(igrid) == dmiss)then
        lat(igrid)=rmiss
     else
        lat(igrid)=dble(rlat(igrid))
     end if
     
     !https://earth.esa.int/documents/10174/1854503/SMOS-Level-2-Ocean-Salinity-v662-release-note
     if(rdat(igrid) == dmiss .or. iqc(igrid) > 150)then 
        dat(igrid)=rmiss
     else
        dat(igrid)=dble(rdat(igrid))
     end if

  end do

  deallocate(rlon,rlat,rdat,iqc)

end subroutine read_smos

!--------------------------------------------------------

subroutine deallocate_smos(lon,lat,dat)

  implicit none
  real(kind = 8),allocatable :: lon(:),lat(:),dat(:)

  deallocate(lon,lat,dat)

end subroutine deallocate_smos

!-------------------------------------------------------

subroutine deallocate_smos_filename(filename)

  implicit none
  character(200),allocatable :: filename(:)

  deallocate(filename)

end subroutine deallocate_smos_filename

end module mod_read_smos
