module mod_read_gcomc

contains

  !---------------------------------------------------
  ! Read GCOM-C |
  !---------------------------------------------------

  subroutine sst_filename(iyr,imon,iday,nfile,filename)

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

    status=system("find /data/R/R2402/DATA/GCOM-C/SST/"//yyyymmdd// &
         & " -name *"//yyyymmdd//"*.h5 > filename.dat")

    nfile=0
    open(1,file="filename.dat",status="old")
    do
       read(1,*,end=100)
       nfile=nfile+1
    end do
100 close(1)

    allocate(filename(nfile))
    open(1,file="filename.dat",status="old")
    do ifile=1,nfile
       write(1,'(a)') filename(ifile)
    end do
    close(1)
    status=system("rm -f filename.dat")
    
    write(*,'(a,i10,a)') "The number of file: ",nfile," at "//yyyymmdd

  end subroutine sst_filename

  !------------------

  subroutine ssuv_filename(iyr,imon,iday,nfile,filename)

    implicit none

    !Common
    integer ifile
    integer status,system
    character(256) command
    character(8) yyyymmdd
    character(4) yyyy
    character(2) mm,dd

    !IN
    integer,intent(in) :: iyr,imon,iday

    !OUT
    integer,intent(out) :: nfile
    character(256),allocatable,intent(out) :: filename(:)

    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    yyyymmdd=yyyy//mm//dd

    command="find /data/R/R2402/DATA/GCOM-C/SSUV/"//yyyy//mm//&
         &" -name curvec_S3*_OL_"//yyyymmdd//"*_GC1SG1_"//yyyymmdd//"*_LCI_J0000_v01.nc"//&
         &" > filename.dat"
    status=system(trim(command))

    !Count nfile
    nfile=0
    open(1,file="filename.dat",status="old")
    do
       read(1,*,end=100)
       nfile=nfile+1
    end do
100 close(1)

    write(*,'(a,i10,a)') "The number of file: ",nfile," at "//yyyymmdd

    !Read filename
    if(nfile /= 0)then
       allocate(filename(nfile))
       open(1,file="filename.dat",status="old")
       do ifile=1,nfile
          read(1,'(a)') filename(ifile)
       end do
       close(1)
    end if
    status=system("rm -f filename.dat")

  end subroutine ssuv_filename

  !-------------------

  subroutine read_gcomc_sst(filename,ihour,im,jm,long,lati,dat)

    use netcdf
    use mod_rmiss
    implicit none

    !Common
    integer,parameter :: dmiss1=65535,dmiss2=65534,dmiss3=65533,dmiss4=65532

    integer i,j,i1,j1,i2,j2
    integer im,jm,imc,jmc
    integer itime

    integer status,system
    integer ncid,grpid,dimid(2),varid
    integer ndim,iparent
    
    integer,allocatable :: idat(:,:),iql(:,:)
    real(kind = 4),allocatable :: clon(:,:),clat(:,:)

    character(16) ql

    !---IN
    character(200),intent(in) :: filename

    !---OUT
    integer,intent(out) :: ihour
    real(kind = 8),allocatable :: long(:,:),lati(:,:),dat(:,:)

    !--- Get imc,jmc,im,jm
    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    !Geometry_data
    status=nf90_inq_grp_ncid(ncid,"Geometry_data",grpid)
    status=nf90_inq_dimids(grpid,ndim,dimid,iparent)
    status=nf90_inquire_dimension(grpid,dimid(2),len = imc)
    status=nf90_inquire_dimension(grpid,dimid(1),len = jmc)

    !Image_data
    status=nf90_inq_grp_ncid(ncid,"Image_data",grpid)
    status=nf90_inq_dimids(grpid,ndim,dimid,iparent)
    status=nf90_inquire_dimension(grpid,dimid(2),len = im)
    status=nf90_inquire_dimension(grpid,dimid(1),len = jm)

    status=nf90_close(ncid)

    write(*,'(2a,i6)') "imc:",imc,"jmc:",jmc
    write(*,'(2a,i6)') "im:",im,"jm:",jm

    !---Allocate
    allocate(clon(imc,jmc),clat(imc,jmc))
    allocate(idat(im,jm),iql(im,jm))
    allocate(long(im,jm),lati(im,jm),dat(im,jm))

    !---Read netcdf
    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    !Geometry_data
    status=nf90_inq_grp_ncid(ncid,"Geometry_data",grpid)
    status=nf90_inq_varid(grpid,"Obs_time",varid)
    status=nf90_get_var(grpid,varid,itime)
    status=nf90_inq_varid(grpid,"Longitude",varid)
    status=nf90_get_var(grpid,varid,clon)
    status=nf90_inq_varid(grpid,"Latitude",varid)
    status=nf90_get_var(grpid,varid,clat)

    !Image_data
    status=nf90_inq_grp_ncid(ncid,"Image_data",grpid)

    status=nf90_inq_varid(grpid,"SST",varid)
    status=nf90_get_var(grpid,varid,idat)
    status=nf90_inq_varid(grpid,"QA_flag",varid)
    status=nf90_get_var(grpid,varid,iql)

    status=nf90_close(ncid)

    !Calculate time

    ihour=nint(itime*0.001d0)

    !Bit-
    !0:nodata
    !1:land
    !2:Rejected by QC
    !3:Retrieval error
    !4:No data(TIR1)
    !5:No data(TIR2)
    !6-7:no
    !8:0:nighttime or no visible data,1:daytime
    !9-10:no
    !11:unknown (clear/cloudy)
    !12:cloudy
    !13:acceptable (possibly cloudy)
    !14:good
    !15:0:unreliable (inland/too close to land),1:reliable
    !*15-left___right-0

    !---Post process
    do j=1,jm
       do i=1,im

          dat(i,j)=rmiss
          write(ql,'(b16.16)') iql(i,j)

          if(idat(i,j) == dmiss1 .or. idat(i,j) == dmiss2 &
               & .or. idat(i,j) == dmiss3 .or. idat(i,j) == dmiss4)then
             dat(i,j)=rmiss
          elseif(ql(16:16) == "1" & !nodata
               & .or. ql(15:15) == "1" & !land
               & .or. ql(14:14) == "1" & !Rejected by QC
               & .or. ql(13:13) == "1" & !Retrieval error
               & .or. ql(12:12) == "1" & !No data (TIR1)
               & .or. ql(11:11) == "1" & !No data (TIR2)
               & .or. ql(4:4) == "1" & !cloudy
               & .or. ql(1:1) == "0")then !No data (unreliable)
             dat(i,j)=rmiss
          else
             dat(i,j)=dble(idat(i,j))*0.0012d0-10.d0
          end if

       end do !i
    end do !j

    !Longitude/Latitude
    do j=1,jm

       j1=(j-1)/10+1
       j2=j1+1

       do i=1,im

          i1=(i-1)/10+1
          i2=i1+1

          if(i1 < 1 .or. imc < i2 .or. j1 < 1 .or. jmc < j2)then
             long(i,j)=rmiss
             lati(i,j)=rmiss
             dat(i,j)=rmiss
             cycle
          end if

          call bi(dble(i1),dble(i2),dble(j1),dble(j2),dble((i-1)*0.1d0+1.d0),dble((j-1)*0.1d0+1.d0), &
               & dble(clon(i1,j1)),dble(clon(i2,j1)),dble(clon(i1,j2)),dble(clon(i2,j2)),long(i,j))
          call bi(dble(i1),dble(i2),dble(j1),dble(j2),dble((i-1)*0.1d0+1.d0),dble((j-1)*0.1d0+1.d0), &
               & dble(clat(i1,j1)),dble(clat(i2,j1)),dble(clat(i1,j2)),dble(clat(i2,j2)),lati(i,j))

          if(long(i,j) < 0.d0) long(i,j)=long(i,j)+360.d0

       end do !i
    end do !j


    !deallocate
    deallocate(clon,clat)
    deallocate(idat,iql)

  end subroutine read_gcomc_sst

  !-----------------------------------------------------------------------
  ! Domain: 123-150E, 24-50N 
  ! Time dif. > 15 min.
  ! Speed < 2 m/s
  ! Correlation > 0.4
  !-----------------------------------------------------------------------

  subroutine read_gcomc_ssuv(filename,n,long,lati,u,v)

    use mod_rmiss
    use netcdf
    implicit none

    real(kind = 4),parameter :: dmiss=-999.e0

    integer i
    integer status
    integer ncid,dimid,varid

    real(kind = 4),allocatable :: cor(:)
    real(kind = 4),allocatable :: t1(:),t2(:)
    real(kind = 4),allocatable :: rlon(:),rlat(:)
    real(kind = 4),allocatable :: ru(:),rv(:)    
    
    !IN
    character(256),intent(in) :: filename

    !OUT
    integer,intent(out) :: n
    real(kind = 8),allocatable,intent(out) :: long(:),lati(:)
    real(kind = 8),allocatable,intent(out) :: u(:),v(:)

    !Open ncfile
    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    !Dimension
    status=nf90_inq_dimid(ncid,"samples",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len = n)

    allocate(cor(n))
    allocate(t1(n),t2(n))
    allocate(rlon(n),rlat(n))
    allocate(ru(n),rv(n))
    
    allocate(long(n),lati(n))
    allocate(u(n),v(n))
    
    !Read DATA
    status=nf90_inq_varid(ncid,"lon",varid)
    status=nf90_get_var(ncid,varid,rlon)

    status=nf90_inq_varid(ncid,"lat",varid)
    status=nf90_get_var(ncid,varid,rlat)

    status=nf90_inq_varid(ncid,"u",varid)
    status=nf90_get_var(ncid,varid,ru)

    status=nf90_inq_varid(ncid,"v",varid)
    status=nf90_get_var(ncid,varid,rv)

    status=nf90_inq_varid(ncid,"corr",varid)
    status=nf90_get_var(ncid,varid,cor)

    status=nf90_inq_varid(ncid,"obs_time_1",varid)
    status=nf90_get_var(ncid,varid,t1)

    status=nf90_inq_varid(ncid,"obs_time_2",varid)
    status=nf90_get_var(ncid,varid,t2)

    !Close
    status=nf90_close(ncid)

    !Post process
    long(:)=dble(rlon(:))
    lati(:)=dble(rlat(:))
    u(:)=dble(ru(:))
    v(:)=dble(rv(:))
    do i=1,n

       if(rlon(i) == dmiss .or. rlat(i) == dmiss .or. ru(i) == dmiss .or. rv(i) == dmiss .or. &
            & cor(i) < 0.4e0 .or. abs(t2(i)-t1(i)) < 15.e0/60.e0)then
          long(i)=rmiss
          lati(i)=rmiss
          u(i)=rmiss
          v(i)=rmiss
       else if(2.d0 < sqrt(u(i)*u(i)+v(i)*v(i)))then
          long(i)=rmiss
          lati(i)=rmiss
          u(i)=rmiss
          v(i)=rmiss
       end if

    end do

    deallocate(cor)
    deallocate(t1,t2)
    deallocate(rlon,rlat)
    deallocate(ru,rv)

  end subroutine read_gcomc_ssuv

  !---------------------------------------

  subroutine deallocate_gcomc_sst(long,lati,dat)

    implicit none

    real(kind = 8),allocatable :: long(:,:),lati(:,:),dat(:,:)
    deallocate(long,lati,dat)

  end subroutine deallocate_gcomc_sst

  !--------------------------------------

  subroutine deallocate_gcomc_filename(filename)

    implicit none

    character(256),allocatable :: filename(:)

    deallocate(filename)

  end subroutine deallocate_gcomc_filename

  !-----

  subroutine deallocate_gcomc_ssuv(long,lati,u,v)

    implicit none

    real(kind = 8),allocatable :: long(:),lati(:)
    real(kind = 8),allocatable :: u(:),v(:)

    deallocate(long,lati,u,v)

  end subroutine deallocate_gcomc_ssuv

  !-----------------------------------------------------------
  ! Bilinear Interpolation |
  !-----------------------------------------------------------
  !
  ! f01(x0,y1) ___________ f11(x1,y1)
  !  |                       |
  ! y|-----------fxy(x,y)    |
  !  |            |          |
  !  |            |          |
  ! f00(x0,y0) ___|_______ f10(x1,y0)
  !               x
  !-----------------------------------------------------------

  subroutine bi(x0,x1,y0,y1,x,y,f00,f10,f01,f11,fxy)

    implicit none
    real(kind = 8),intent(in) :: x0,x1,y0,y1,x,y
    real(kind = 8),intent(in) :: f00,f10,f01,f11
    real(kind = 8),intent(out) :: fxy

    fxy=(y1-y)/(y1-y0)*((x1-x)/(x1-x0)*f00+(x-x0)/(x1-x0)*f10) &
         & +(y-y0)/(y1-y0)*((x1-x)/(x1-x0)*f01+(x-x0)/(x1-x0)*f11)

  end subroutine bi

end module mod_read_gcomc
