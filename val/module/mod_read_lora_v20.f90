module mod_read_lora

!  character(100),parameter :: pdir="/vol0004/ra000007/data/a04048"
  character(100),parameter :: pdir="/data/R/R2402/ohishi"
  
contains

  !------------------------------------------------------------------
  ! Read grid |
  !------------------------------------------------------------------

  subroutine read_grid(dir,&
       & lont,lonu,lonv, &
       & latt,latu,latv, &
       & dept,depu,depv, &
       & maskt,masku,maskv)

    use mod_gridinfo
    use netcdf
    implicit none

    integer i,j,k
    integer status,access
    integer ncid,varid

    real(kind = 8) tmp2d(im,jm)
    real(kind = 8) z(im,jm,km),h(im,jm)

    real(kind = 8),intent(out) :: lont(im),latt(jm)
    real(kind = 8),intent(out) :: lonu(im),latu(jm)
    real(kind = 8),intent(out) :: lonv(im),latv(jm)
    real(kind = 8),intent(out) :: dept(im,jm,km),depu(im,jm,km),depv(im,jm,km)
    real(kind = 8),intent(out) :: maskt(im,jm),masku(im,jm),maskv(im,jm)

    character(*) dir
    character(100) filename

    filename=trim(pdir)//"/"//trim(dir)//"/prep/in/grid.nc"

    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: "//trim(filename)//" not found"
       stop
    end if

    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    !status=nf90_inq_varid(ncid,"z_w",varid)    
    status=nf90_inq_varid(ncid,"z_e",varid)
    status=nf90_get_var(ncid,varid,z)

    status=nf90_inq_varid(ncid,"east_e",varid)
    status=nf90_get_var(ncid,varid,tmp2d)
    lont(:)=tmp2d(:,1)

    status=nf90_inq_varid(ncid,"east_u",varid)
    status=nf90_get_var(ncid,varid,tmp2d)
    lonu(:)=tmp2d(:,1)

    status=nf90_inq_varid(ncid,"east_v",varid)
    status=nf90_get_var(ncid,varid,tmp2d)
    lonv(:)=tmp2d(:,1)
    
    status=nf90_inq_varid(ncid,"north_e",varid)
    status=nf90_get_var(ncid,varid,tmp2d)
    latt(:)=tmp2d(1,:)
    
    status=nf90_inq_varid(ncid,"north_u",varid)
    status=nf90_get_var(ncid,varid,tmp2d)
    latu(:)=tmp2d(1,:)

    status=nf90_inq_varid(ncid,"north_v",varid)
    status=nf90_get_var(ncid,varid,tmp2d)
    latv(:)=tmp2d(1,:)
    
    status=nf90_inq_varid(ncid,"h",varid)
    status=nf90_get_var(ncid,varid,h)

    status=nf90_inq_varid(ncid,"fsm",varid)
    status=nf90_get_var(ncid,varid,maskt)

    status=nf90_inq_varid(ncid,"dum",varid)
    status=nf90_get_var(ncid,varid,masku)

    status=nf90_inq_varid(ncid,"dvm",varid)
    status=nf90_get_var(ncid,varid,maskv)
    
    status=nf90_close(ncid)

    !Depth
    do k=1,km
       do j=1,jm
          do i=1,im
             dept(i,j,k)=abs(z(i,j,k)*h(i,j)*maskt(i,j))
             depu(i,j,k)=abs(z(i,j,k)*h(i,j)*masku(i,j))
             depv(i,j,k)=abs(z(i,j,k)*h(i,j)*maskv(i,j))
          end do
       end do
    end do

  end subroutine read_grid

  !------------------------------------------------------------------
  ! Read LORA dataset |
  !------------------------------------------------------------------

  subroutine read_anal(dir,letkf,region,ms,imem,var,iyr,imon,iday,im,jm,km,mask,dat)

    use mod_rmiss
    use netcdf
    implicit none

    real(kind = 4),parameter :: dmiss=9.96921e36

    !---Common
    integer status,access
    integer ncid,varid
    integer i,j,k

    real(kind = 4) ddat(im,jm,km)
    
    character(100) filename
    character(6) yyyymm
    character(5) mmmmm    
    character(4) yyyy
    character(2) mm,dd

    !---IN
    character(*),intent(in) :: dir,letkf,region,ms,var

    integer,intent(in) :: imem    
    integer,intent(in) :: iyr,imon,iday
    integer,intent(in) :: im,jm,km

    real(kind = 8),intent(in) :: mask(im,jm)
    
    !---OUT
    real(kind = 8),intent(out) :: dat(im,jm,km)

    write(mmmmm,'(i5.5)') imem
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    yyyymm=yyyy//mm

    if(trim(ms) == "mean" .or. trim(ms) == "sprd")then
       filename=trim(pdir)//"/"//trim(dir)//"/"//trim(letkf)//&
            &"/output/"//trim(ms)//"/"//trim(region)//yyyymm//".nc "
    else if(trim(ms) == "eens")then
       filename=trim(pdir)//"/"//trim(dir)//"/"//trim(letkf)//&
            &"/output/"//trim(ms)//"/"//trim(region)//yyyymm//"."//mmmmm//".nc "
    end if
    
    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: "//trim(filename)//" not found"
       stop
    end if
    
    status=nf90_open(trim(filename),nf90_nowrite,ncid)
    
    status=nf90_inq_varid(ncid,trim(var),varid)
    if(status == nf90_noerr)then

       if(trim(ms) == "mean" .or. trim(ms) == "sprd")then
          status=nf90_get_var(ncid,varid,ddat,(/1,1,1,iday/),(/im,jm,km,1/))
       else if(trim(ms) == "eens")then
          status=nf90_get_var(ncid,varid,ddat,(/1,1,1,iday/),(/im,jm,1,1/))
       end if
       
    end if
    status=nf90_close(ncid)

    do k=1,km
       do j=1,jm
          do i=1,im
             if(ddat(i,j,k) == dmiss .or. mask(i,j) == 0.d0)then
                dat(i,j,k)=rmiss
             else
                dat(i,j,k)=dble(ddat(i,j,k))
             end if
          end do
       end do
    end do
    
  end subroutine read_anal

  !---------------------------------------------------------------
  ! Read innovation statistics |
  !---------------------------------------------------------------

  subroutine read_inv(dir,letkf,iyr,imon,iday,nobs,nmem, &
       & ele,ins, &
       & lon,lat,lev,obs,err,hxfmean,hxfsprd,hxamean,hxasprd,hxf,hxa)

    use netcdf
    implicit none

    !---Common
    integer status,access
    integer ncid,dimid,varid

    character(200) filename    
    character(4) yyyy
    character(2) mm,dd
    
    !---IN
    integer,intent(in) :: iyr,imon,iday

    character(*),intent(in) :: dir,letkf
    
    !---OUT
    integer,intent(out) :: nobs,nmem
    integer,allocatable,intent(out) :: ele(:),ins(:)

    real,allocatable,intent(out) :: lon(:),lat(:),lev(:)
    real,allocatable,intent(out) :: obs(:),err(:)
    real,allocatable,intent(out) :: hxfmean(:),hxfsprd(:)
    real,allocatable,intent(out) :: hxamean(:),hxasprd(:)
    real,allocatable,intent(out) :: hxf(:,:),hxa(:,:)

    !Filename
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    
    filename=trim(pdir)//"/"//trim(dir)//"/"//trim(letkf)//&
         & "/output/mean/inv."//yyyy//mm//dd//".nc"

    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not found "//trim(filename)
       stop
    end if
    
    !Open
    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    !Get nobs
    status=nf90_inq_dimid(ncid,"nobs",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len=nobs)

    !Get nmem
    status=nf90_inq_dimid(ncid,"nmem",dimid)
    status=nf90_inquire_dimension(ncid,dimid,len=nmem)

    !Allocate
    allocate(ele(nobs),ins(nobs))
    allocate(lon(nobs),lat(nobs),lev(nobs))
    allocate(obs(nobs),err(nobs))
    allocate(hxfmean(nobs),hxfsprd(nobs))
    allocate(hxamean(nobs),hxasprd(nobs))
    allocate(hxf(nobs,nmem),hxa(nobs,nmem))

    !Get element
    status=nf90_inq_varid(ncid,"ele",varid)
    status=nf90_get_var(ncid,varid,ele)

    !Get instrument
    status=nf90_inq_varid(ncid,"ins",varid)
    status=nf90_get_var(ncid,varid,ins)
    
    !Get longitude
    status=nf90_inq_varid(ncid,"lon",varid)
    status=nf90_get_var(ncid,varid,lon)

    !Get latitude
    status=nf90_inq_varid(ncid,"lat",varid)
    status=nf90_get_var(ncid,varid,lat)

    !Get level
    status=nf90_inq_varid(ncid,"lev",varid)
    status=nf90_get_var(ncid,varid,lev)
    
    !Get observation
    status=nf90_inq_varid(ncid,"obs",varid)
    status=nf90_get_var(ncid,varid,obs)

    !Get observation error
    status=nf90_inq_varid(ncid,"err",varid)
    status=nf90_get_var(ncid,varid,err)

    !Get H(xfmean)
    status=nf90_inq_varid(ncid,"hxfmean",varid)
    status=nf90_get_var(ncid,varid,hxfmean)

    !Get H(xfsprd)
    status=nf90_inq_varid(ncid,"hxfsprd",varid)
    status=nf90_get_var(ncid,varid,hxfsprd)

    !Get H(xamean)
    status=nf90_inq_varid(ncid,"hxamean",varid)
    status=nf90_get_var(ncid,varid,hxamean)

    !Get H(xasprd)
    status=nf90_inq_varid(ncid,"hxasprd",varid)
    status=nf90_get_var(ncid,varid,hxasprd)

    !Get H(xf)
    status=nf90_inq_varid(ncid,"hxf",varid)
    status=nf90_get_var(ncid,varid,hxf)

    !Get H(xa)
    status=nf90_inq_varid(ncid,"hxa",varid)
    status=nf90_get_var(ncid,varid,hxa)

    !Close
    status=nf90_close(ncid)

  end subroutine read_inv

  !-----------------------------------------

  subroutine end_read_inv(ele,ins, &
       & lon,lat,lev,obs,err,hxfmean,hxfsprd,hxamean,hxasprd,hxf,hxa)

    implicit none

    integer,allocatable :: ele(:),ins(:)

    real,allocatable :: lon(:),lat(:),lev(:)
    real,allocatable :: obs(:),err(:)
    real,allocatable :: hxfmean(:),hxfsprd(:)
    real,allocatable :: hxamean(:),hxasprd(:)
    real,allocatable :: hxf(:,:),hxa(:,:)
    
    deallocate(ele,ins)
    deallocate(lon,lat,lev)
    deallocate(obs,err)
    deallocate(hxfmean,hxfsprd)
    deallocate(hxamean,hxasprd)
    deallocate(hxf,hxa)
    
  end subroutine end_read_inv
  
end module mod_read_lora
