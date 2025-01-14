module mod_read_smap

  !Modified by S.Ohishi 2020.04
  !If you change directory, please check "to be modified".

contains

  !----------------------------------------------------
  ! Filename |
  !----------------------------------------------------

  subroutine read_filename(iyr,imon,iday,nfile,filename)

    implicit none

    !Common
    integer ifile
    integer status,system

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

    !*to be modified: 
    !V4.3 --> V5.0 @ 2021.06.16 Ohishi
    status=system("find /data/R/R2402/DATA/SMAP/"//yyyy//"/*/"//&
         &" -name *"//yyyymmdd//"*V5.0.h5 "//&
         &"> smap"//yyyymmdd//".txt")

    nfile=0
    open(1,file="smap"//yyyymmdd//".txt",status="old")
    do 
       read(1,*,end=100)
       nfile=nfile+1
    end do
100 close(1)

    write(*,'(a,i6,a)') "The number of file:",nfile," at "//yyyymmdd

    if(nfile /= 0)then

       !Read filename
       allocate(filename(nfile))

       open(1,file="smap"//yyyymmdd//".txt",status="old")
       do ifile=1,nfile
          read(1,'(a)') filename(ifile)
       end do
       close(1)

    end if
    status=system("rm -f smap"//yyyymmdd//".txt")

  end subroutine read_filename

  !----------------------------------------------------
  ! Read SMAP |
  !----------------------------------------------------

  subroutine read_smap(filename,iyr,imon,iday,ihour,im,jm,lon,lat,dat)

    use mod_rmiss
    use netcdf
    implicit none

    real(kind = 4),parameter :: dmiss=-9999.e0

    integer i,j
    integer status,access,system
    integer ncid,dimid,varid

    integer,allocatable :: iqc(:,:)
    real(kind = 4),allocatable :: rlon(:,:),rlat(:,:),rdat(:,:)

    !IN
    character(200),intent(in) :: filename

    !OUT
    integer,intent(out) :: iyr,imon,iday,ihour
    integer,intent(out) :: im,jm

    real(kind = 8),allocatable,intent(out) :: lon(:,:),lat(:,:),dat(:,:)  

    status=access(trim(filename)," ")
    if(status == 0)then
       write(*,'(a)') "Read "//trim(filename)
    else
       write(*,'(a)') "***Error: Not found "//trim(filename)
       stop
    end if

    !To be modified
    read(filename(53:56),*) iyr
    read(filename(57:58),*) imon
    read(filename(59:60),*) iday
    read(filename(62:63),*) ihour

    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    !Get im,jm
    status=nf90_inq_dimid(ncid,"phony_dim_1",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = im)

    status=nf90_inq_dimid(ncid,"phony_dim_0",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = jm)

    !Allocate
    allocate(rlon(im,jm),rlat(im,jm),rdat(im,jm),iqc(im,jm))
    allocate(lon(im,jm),lat(im,jm),dat(im,jm))

    !Read Data
    status=nf90_inq_varid(ncid,"lon",varid)
    status=nf90_get_var(ncid,varid,rlon)

    status=nf90_inq_varid(ncid,"lat",varid)
    status=nf90_get_var(ncid,varid,rlat)

    status=nf90_inq_varid(ncid,"quality_flag",varid)
    status=nf90_get_var(ncid,varid,iqc)

    status=nf90_inq_varid(ncid,"smap_sss",varid)
    status=nf90_get_var(ncid,varid,rdat)

    status=nf90_close(ncid)

    !Modify unit
    do j=1,jm
       do i=1,im

          if(rlon(i,j) == dmiss)then
             lon(i,j)=rmiss
          elseif(rlon(i,j) < 0.e0)then
             lon(i,j)=dble(rlon(i,j))+360.d0
          else
             lon(i,j)=dble(rlon(i,j))
          end if

          if(rlat(i,j) == dmiss)then
             lat(i,j)=rmiss
          else
             lat(i,j)=dble(rlat(i,j))
          end if

          if(rdat(i,j) == dmiss .or. iqc(i,j) /= 0)then
             dat(i,j)=rmiss
          else
             dat(i,j)=dble(rdat(i,j))
          end if

       end do
    end do

    deallocate(rlon,rlat,rdat,iqc)

  end subroutine read_smap

  !--------------------------------------------------------

  subroutine deallocate_smap(lon,lat,dat)

    implicit none
    real(kind = 8),allocatable :: lon(:,:),lat(:,:),dat(:,:)

    deallocate(lon,lat,dat)

  end subroutine deallocate_smap

  !---------------------------------------------------------

  subroutine deallocate_smap_filename(filename)

    implicit none
    character(200),allocatable :: filename(:)

    deallocate(filename)

  end subroutine deallocate_smap_filename

end module mod_read_smap

