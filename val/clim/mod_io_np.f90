module mod_io

contains

  !-----------------------------------------------------------------------------------------------------

  subroutine get_grid_size(idat,im,jm,km)

    use mod_gridinfo, im_lora => im, jm_lora => jm, km_lora => km
    use mod_read_bran2020,   only: im_bran => im, jm_bran => jm, km_bran => km
    use mod_read_fora_np60,  only: im_fora => im, jm_fora => jm, km_fora => km
    use mod_read_glorys12v1, only: im_g010 => im, jm_g010 => jm, km_g010 => km
    use mod_read_jcope_fgo,  only: im_jcope => im, jm_jcope => jm, km_jcope => km
    implicit none

    !---IN
    integer,intent(in) :: idat

    !---OUT
    integer,intent(out) :: im,jm,km    
    
    if(idat == 1)then !---LORA-NP
        im=im_lora
        jm=jm_lora
        km=km_lora
     else if(idat == 2)then !---BRAN2020
        im=im_bran
        jm=jm_bran
        km=km_bran
     else if(idat == 3)then !---FORA-NP60
        im=im_fora
        jm=jm_fora
        km=km_fora
     else if(idat == 4)then !---GLORYS12V1
        im=im_g010
        jm=jm_g010
        km=km_g010
     else if(idat == 5)then !---JCOPE-FGO
        im=im_jcope
        jm=jm_jcope
        km=km_jcope
    else
       write(*,*) "***Error: idat => ",idat
       stop
    end if

  end subroutine get_grid_size
  
  !-----------------------------------------------------------------------------------------------------

  subroutine read_grid(idat,im,jm,km,lont,lonu,lonv,latt,latu,latv,dept,depu,depv,maskt,masku,maskv)

    use setting, only: datname
    use mod_read_lora, only: read_grid_lora => read_grid
    use mod_read_bran2020,   only: read_bran2020
    use mod_read_fora_np60,  only: read_fora_np60
    use mod_read_glorys12v1, only: read_glorys12v1
    use mod_read_jcope_fgo,  only: read_grid_jcope => read_grid
    implicit none

    !---Common
    integer i,j

    real(kind = 8),allocatable :: tmp1dz(:)
    real(kind = 8),allocatable :: tmp3d(:,:,:)

    !LORA
    character(10) :: dir

    !Others
    character(1) varname

    !---IN
    integer,intent(in) :: idat
    integer,intent(in) :: im,jm,km

    !---OUT
    real(kind = 8),intent(out) :: lont(im),lonu(im),lonv(im)
    real(kind = 8),intent(out) :: latt(jm),latu(jm),latv(jm)
    real(kind = 8),intent(out) :: dept(im,jm,km),depu(im,jm,km),depv(im,jm,km)
    real(kind = 8),intent(out) :: maskt(im,jm),masku(im,jm),maskv(im,jm)
    
    allocate(tmp1dz(km))
    allocate(tmp3d(im,jm,km))
    
    if(idat == 1)then !---LORA-NP
       dir="NP"
       call read_grid_lora(dir,lont,lonu,lonv, &
            & latt,latu,latv, &
            & dept,depu,depv,tmp3d, &
            & maskt,masku,maskv)
    else if(idat == 2)then !---BRAN2020
       varname="t"
       call read_bran2020(varname,2003,1,1,km,lont,latt,tmp1dz,maskt,tmp3d)
       varname="u"
       call read_bran2020(varname,2003,1,1,km,lonu,latu,tmp1dz,masku,tmp3d)
       varname="v"
       call read_bran2020(varname,2003,1,1,km,lonv,latv,tmp1dz,maskv,tmp3d)
       do j=1,jm
          do i=1,im
             dept(i,j,1:km)=tmp1dz(1:km)
             depu(i,j,1:km)=tmp1dz(1:km)
             depv(i,j,1:km)=tmp1dz(1:km)
          end do
       end do
    else if(idat == 3)then !---FORA-NP60
       varname="t"
       call read_fora_np60(varname,2003,1,1,km,lont,latt,tmp1dz,maskt,tmp3d)
       varname="u"
       call read_fora_np60(varname,2003,1,1,km,lonu,latu,tmp1dz,masku,tmp3d)
       varname="v"
       call read_fora_np60(varname,2003,1,1,km,lonv,latv,tmp1dz,maskv,tmp3d)
       do j=1,jm
          do i=1,im
             dept(i,j,1:km)=tmp1dz(1:km)
             depu(i,j,1:km)=tmp1dz(1:km)
             depv(i,j,1:km)=tmp1dz(1:km)
          end do
       end do
    else if(idat == 4)then !---GLORYS12V1 (Probably Regridded)
       varname="t"
       call read_glorys12v1(varname,2003,1,1,km,lont,latt,tmp1dz,maskt,tmp3d)
       lonu(:)=lont(:)
       lonv(:)=lont(:)
       latu(:)=latt(:)
       latv(:)=latt(:)
       masku(:,:)=maskt(:,:)
       maskv(:,:)=maskt(:,:)
       do j=1,jm
          do i=1,im
             dept(i,j,1:km)=tmp1dz(1:km)
             depu(i,j,1:km)=tmp1dz(1:km)
             depv(i,j,1:km)=tmp1dz(1:km)
          end do
       end do   
    else if(idat == 5)then !---JCOPE-FGO
       call read_grid_jcope(lont,lonu,lonv,latt,latu,latv,maskt,dept)
       masku(:,:)=maskt(:,:)
       maskv(:,:)=maskt(:,:)
       depu(:,:,:)=dept(:,:,:)
       depv(:,:,:)=dept(:,:,:)
    end if

    deallocate(tmp1dz)
    deallocate(tmp3d)
    
  end subroutine read_grid

  !-----------------------------------------------------------------------------------

  subroutine read_dat(varname,idat,iyr,imon,iday,im,jm,km,mask,mean,sprd)

    use setting, only: datname
    use mod_read_lora,       only: read_anal
    use mod_read_bran2020,   only: read_bran2020
    use mod_read_fora_np60,  only: read_fora_np60
    use mod_read_glorys12v1, only: read_glorys12v1
    use mod_read_jcope_fgo,  only: read_jcope_fgo    
    use mod_rmiss
    implicit none

    !---Common
    !integer i,j
    integer k

    real(kind = 8),allocatable :: tmp1dx(:),tmp1dy(:),tmp1dz(:)
    real(kind = 8),allocatable :: tmp2d(:,:)

    !LORA
    integer imem !Dummy
    character(10) :: dir
    character(10) :: letkf
    character(10) :: region
    character(10) :: ms

    !GLORYS
    character(1) varname

    !---IN
    integer,intent(in) :: idat
    integer,intent(in) :: iyr,imon,iday
    integer,intent(in) :: im,jm,km

    real(kind = 8),intent(in) :: mask(im,jm)

    !---OUT
    real(kind = 8),intent(out) :: mean(im,jm,km),sprd(im,jm,km)

    allocate(tmp1dx(im),tmp1dy(jm),tmp1dz(km))
    allocate(tmp2d(im,jm))

    if(idat == 1)then !---LORA-NP

       imem=0
       dir="NP        "
       letkf="letkf     "
       region="np        "

       !Ensemble Mean
       ms="mean      "
       if(varname == "h")then
          k=1
          call read_anal(dir,letkf,region,ms,imem,"el",iyr,imon,iday,im,jm,k,mask,mean)
       else
          call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,km,mask,mean)
       end if

       !Ensemble Spread
       ms="sprd      "
       if(varname == "h")then
          k=1
          call read_anal(dir,letkf,region,ms,imem,"el",iyr,imon,iday,im,jm,k,mask,sprd)
       else
          call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,km,mask,sprd)
       end if

    else if(idat == 2)then !---BRAN2020

       !Determistic analysis
       if(varname == "h")then
          k=1          
          call read_bran2020(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,mean)       
       else
          call read_bran2020(varname,iyr,imon,iday,km,tmp1dx,tmp1dy,tmp1dz,tmp2d,mean)
       end if
          
       !Ensemble Spread
       sprd(:,:,:)=rmiss

    else if(idat == 3)then !---FORA-NP60

       !Deterministic analysis
       if(varname == "h")then
          k=1
          call read_fora_np60(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,mean)          
       else
          call read_fora_np60(varname,iyr,imon,iday,km,tmp1dx,tmp1dy,tmp1dz,tmp2d,mean)          
       end if       
       
       !Ensemble Spread
       sprd(:,:,:)=rmiss
       
    else if(idat == 4)then !---GLORYS12V1

       !Deterministic analysis
       if(varname == "h")then
          k=1
          call read_glorys12v1(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,mean)          
       else
          call read_glorys12v1(varname,iyr,imon,iday,km,tmp1dx,tmp1dy,tmp1dz,tmp2d,mean)
       end if       
       
       !Ensemble Spread
       sprd(:,:,:)=rmiss
       
    else if(idat == 5)then !---JCOPE-FGO

       !Deterministic analysis
       if(varname == "h")then
          k=1
          call read_jcope_fgo(varname,iyr,imon,iday,k,mask,mean)
       else
          call read_jcope_fgo(varname,iyr,imon,iday,km,mask,mean)
       end if       
       
       !Ensemble Spread
       sprd(:,:,:)=rmiss

    end if

    deallocate(tmp1dx,tmp1dy,tmp1dz)
    deallocate(tmp2d)

  end subroutine read_dat

  !--------------------------------------------------------------------------------

  subroutine write_grid(filename,im,jm,km,lont,lonu,lonv,latt,latu,latv,dept,depu,depv,maskt,masku,maskv)

    use netcdf
    implicit none

    !---Common
    integer status,access
    integer ncid,varid

    !---IN
    integer,intent(in) :: im,jm,km

    real(kind = 8),intent(in) :: lont(im),lonu(im),lonv(im)
    real(kind = 8),intent(in) :: latt(jm),latu(jm),latv(jm)
    real(kind = 8),intent(in) :: dept(im,jm,km),depu(im,jm,km),depv(im,jm,km)
    real(kind = 8),intent(in) :: maskt(im,jm),masku(im,jm),maskv(im,jm)

    character(100),intent(in) :: filename

    !---Check access
    status=access(trim(filename)," ")
    if(status == 0)then
       write(*,*) "Write grid information to "//trim(filename)
    else
       write(*,*) "***Error: Not found "//trim(filename)
       stop
    end if

    !Open NetCDF file
    status=nf90_open(trim(filename),nf90_write,ncid)

    !Longitude
    status=nf90_inq_varid(ncid,"lont",varid)
    status=nf90_put_var(ncid,varid,lont)

    status=nf90_inq_varid(ncid,"lonu",varid)
    status=nf90_put_var(ncid,varid,lonu)

    status=nf90_inq_varid(ncid,"lonv",varid)
    status=nf90_put_var(ncid,varid,lonv)

    !Latitude
    status=nf90_inq_varid(ncid,"latt",varid)
    status=nf90_put_var(ncid,varid,latt)

    status=nf90_inq_varid(ncid,"latu",varid)
    status=nf90_put_var(ncid,varid,latu)

    status=nf90_inq_varid(ncid,"latv",varid)
    status=nf90_put_var(ncid,varid,latv)

    !Depth
    status=nf90_inq_varid(ncid,"dept",varid)
    status=nf90_put_var(ncid,varid,dept)

    status=nf90_inq_varid(ncid,"depu",varid)
    status=nf90_put_var(ncid,varid,depu)

    status=nf90_inq_varid(ncid,"depv",varid)
    status=nf90_put_var(ncid,varid,depv)

    !Mask
    status=nf90_inq_varid(ncid,"maskt",varid)
    status=nf90_put_var(ncid,varid,maskt)

    status=nf90_inq_varid(ncid,"masku",varid)
    status=nf90_put_var(ncid,varid,masku)

    status=nf90_inq_varid(ncid,"maskv",varid)
    status=nf90_put_var(ncid,varid,maskv)

    !Close
    status=nf90_close(ncid)    

  end subroutine write_grid

  !--------------------------------------------------------------------------------

  subroutine write_clim_nc(filename,varname,im,jm,km,imon,mean,sprd)

    use netcdf
    implicit none

    !---Common
    integer status,access
    integer ncid,varid

    !---IN
    integer,intent(in) :: im,jm,km
    integer,intent(in) :: imon

    real(kind = 8),intent(in) :: mean(im,jm,km),sprd(im,jm,km)

    character(100),intent(in) :: filename    
    character(1),intent(in) :: varname

    !---Check file    
    status=access(filename," ")
    if(status == 0)then
       write(*,*) "Write to "//trim(filename)
    else
       write(*,*) "***Error: Not found "//trim(filename)
       stop
    end if

    !---Open NetCDFfile
    status=nf90_open(trim(filename),nf90_write,ncid)

    !---Mean
    status=nf90_inq_varid(ncid,trim(varname)//"mean",varid)

    if(km == 1)then !2D

       if(imon == 0)then !Annual
          status=nf90_put_var(ncid,varid,mean,(/1,1,1/),(/im,jm,1/))
       else !Monthly
          status=nf90_put_var(ncid,varid,mean,(/1,1,imon/),(/im,jm,1/))
       end if

    else !3D

       if(imon == 0)then !Annual
          status=nf90_put_var(ncid,varid,mean,(/1,1,1,1/),(/im,jm,km,1/))
       else !Monthly
          status=nf90_put_var(ncid,varid,mean,(/1,1,1,imon/),(/im,jm,km,1/))
       end if

    end if

    !---Spread
    status=nf90_inq_varid(ncid,trim(varname)//"sprd",varid)

    if(km == 1)then !2D

       if(imon == 0)then !Annual
          status=nf90_put_var(ncid,varid,sprd,(/1,1,1/),(/im,jm,1/))
       else !Monthly
          status=nf90_put_var(ncid,varid,sprd,(/1,1,imon/),(/im,jm,1/))
       end if

    else !3D

       if(imon == 0)then !Annual
          status=nf90_put_var(ncid,varid,sprd,(/1,1,1,1/),(/im,jm,km,1/))
       else !Monthly
          status=nf90_put_var(ncid,varid,sprd,(/1,1,1,imon/),(/im,jm,km,1/))
       end if

    end if

    status=nf90_close(ncid)    

  end subroutine write_clim_nc

  !-------------------------------------------------------------------------

  subroutine read_clim_nc(filename,varname,im,jm,km,imon,mean,sprd)

    use netcdf
    implicit none

    !---Common
    integer status,access
    integer ncid,varid

    !---IN
    integer,intent(in) :: im,jm,km
    integer,intent(in) :: imon

    character(100),intent(in) :: filename    
    character(1),intent(in) :: varname

    !---OUT
    real(kind = 8),intent(out) :: mean(im,jm,km),sprd(im,jm,km)
    
    !---Check file    
    status=access(filename," ")
    if(status == 0)then
       write(*,*) "Read "//trim(filename)
    else
       write(*,*) "***Error: Not found "//trim(filename)
       stop
    end if

    !---Open NetCDFfile
    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    !---Mean
    status=nf90_inq_varid(ncid,trim(varname)//"mean",varid)

    if(km == 1)then !2D

       if(imon == 0)then !Annual
          status=nf90_get_var(ncid,varid,mean(:,:,1),(/1,1,1/),(/im,jm,1/))
       else !Monthly
          status=nf90_get_var(ncid,varid,mean(:,:,1),(/1,1,imon/),(/im,jm,1/))
       end if

    else !3D

       if(imon == 0)then !Annual
          status=nf90_get_var(ncid,varid,mean,(/1,1,1,1/),(/im,jm,km,1/))
       else !Monthly
          status=nf90_get_var(ncid,varid,mean,(/1,1,1,imon/),(/im,jm,km,1/))
       end if

    end if

    !---Spread
    status=nf90_inq_varid(ncid,trim(varname)//"sprd",varid)

    if(km == 1)then !2D

       if(imon == 0)then !Annual
          status=nf90_get_var(ncid,varid,sprd(:,:,1),(/1,1,1/),(/im,jm,1/))
       else !Monthly
          status=nf90_get_var(ncid,varid,sprd(:,:,1),(/1,1,imon/),(/im,jm,1/))
       end if

    else !3D

       if(imon == 0)then !Annual
          status=nf90_get_var(ncid,varid,sprd,(/1,1,1,1/),(/im,jm,km,1/))
       else !Monthly
          status=nf90_get_var(ncid,varid,sprd,(/1,1,1,imon/),(/im,jm,km,1/))
       end if

    end if

    status=nf90_close(ncid)    

  end subroutine read_clim_nc

  !-------------------------------------------------------------------------

  subroutine write_mdot(ndat,im,jm,lon,lat,obs,hdat)

    implicit none

    !---Common
    integer i,j

    character(100) format,filename
    
    !---IN
    integer,intent(in) :: ndat
    integer,intent(in) :: im,jm

    real(kind = 8),intent(in) :: lon(im),lat(jm)
    real(kind = 8),intent(in) :: obs(im,jm),hdat(im,jm,ndat)

    write(format,'(a,I0,a)') '(3f12.5,',ndat,'f12.5)'

    filename="dat/mdot.dat"
    open(1,file=trim(filename),status="replace")
    do j=1,jm
       do i=1,im
          write(1,trim(format)) lon(i),lat(j),obs(i,j),hdat(i,j,1:ndat)
       end do
    end do       
    close(1)
    
  end subroutine write_mdot

  !----------------------------------------------------------------------------

  subroutine write_var3d_txt(datname,varname, &
       & im,jm,km,lon,lat,dep,hdat,obs)

    implicit none

    !---Common
    integer i,j,k

    character(100) filename
    
    !---IN
    integer,intent(in) :: im,jm,km

    real(kind = 8),intent(in) :: lon(im),lat(jm),dep(km)

    real(kind = 8),intent(in) :: hdat(im,jm,km),obs(im,jm,km)
    
    character(10),intent(in) :: datname
    character(1),intent(in) :: varname

    filename="dat/"//trim(datname)//"_"//trim(varname)//".dat"

    open(1,file=trim(filename),status="replace")
    do k=1,km
       do j=1,jm
          do i=1,im
             write(1,'(5f12.5)') lon(i),lat(j),dep(k),hdat(i,j,k),obs(i,j,k)
          end do
       end do
    end do
    close(1)    
    
  end subroutine write_var3d_txt

  !------------------------------------------------------------------------

  subroutine write_var2d_txt(datname,varname,jm,km,lat,dep,hdat,obs)

    use mod_rmiss
    implicit none

    !---Parameter
    real(kind = 8),parameter :: dz=5.d0,dep_max=5500.d0
    
    !---Common
    integer j,k
    integer km_out
    
    real(kind = 8),allocatable :: dep_out(:)
    real(kind = 8),allocatable :: hdat_out(:,:),obs_out(:,:)
    
    character(100) filename
    
    !---IN
    integer,intent(in) :: jm,km

    real(kind = 8),intent(in) :: lat(jm),dep(km)
    real(kind = 8),intent(in) :: hdat(jm,km),obs(jm,km)
    
    character(10),intent(in) :: datname
    character(1),intent(in) :: varname

    !---Output depth
    km_out=int(5500/dz)+1
    allocate(dep_out(km_out))
    allocate(hdat_out(jm,km_out),obs_out(jm,km_out))
    
    do k=1,km_out
       dep_out(k)=(k-1)*dz
    end do

    hdat_out(:,:)=rmiss
    obs_out(:,:)=rmiss
    do k=1,km_out
       do j=1,jm
          call vertical_interpolation(km,dep,hdat(j,:),dep_out(k),hdat_out(j,k))
          call vertical_interpolation(km,dep,obs(j,:),dep_out(k),obs_out(j,k))
       end do
    end do

    !---Write
    filename="dat/"//trim(datname)//"_"//trim(varname)//"2d.dat"
    
    open(1,file=trim(filename),status="replace")
    do k=1,km_out
       do j=1,jm
          write(1,'(4f12.5)') lat(j),dep_out(k),hdat_out(j,k),obs_out(j,k)
       end do
    end do
    close(1)

    deallocate(dep_out)
    deallocate(hdat_out,obs_out)
    
  end subroutine write_var2d_txt
    
end module mod_io
