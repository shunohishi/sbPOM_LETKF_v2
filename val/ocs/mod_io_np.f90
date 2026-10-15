module mod_io

contains

  !---------------------------------------------------------------------------------
  ! Grid size (Analysis data) |
  !---------------------------------------------------------------------------------

  ! *** To be modified ***
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
    else if(idat == 4)then !---GLORYS010
       im=im_g010
       jm=jm_g010
       km=km_g010
    else if(idat == 5)then !---JCOPE-FGO
       im=im_jcope
       jm=jm_jcope
       km=km_jcope
    else
       write(*,*) "***Error: Incorecot idat_a => ",idat
       stop
    end if

  end subroutine get_grid_size

  !-----------------------------------------------------------------------
  ! Read grid data |
  !-----------------------------------------------------------------------

  !*** To be modified
  subroutine read_grid(idat, &
       & im,jm,km,lont,lonu,lonv,latt,latu,latv,dept,depu,depv,maskt,masku,maskv)

    use setting, only: datname
    use mod_read_lora, only: read_grid_lora => read_grid
    use mod_read_bran2020,   only: extract_bran2020
    use mod_read_fora_np60,  only: extract_fora_np60
    use mod_read_glorys12v1, only: extract_grid_glorys12v1 => extract_grid
    use mod_read_jcope_fgo,  only: read_grid_jcope => read_grid
    implicit none

    !---Common
    integer i,j

    real(kind = 8),allocatable :: tmp1dz(:)
    real(kind = 8),allocatable :: tmp3d(:,:,:)

    !LORA
    character(10) dir

    !GLORYS
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
    else if(idat == 2)then !---BRAN
       varname="t"
       call extract_bran2020(varname,2003,1,1,1,im,1,jm,1,km,lont,latt,tmp1dz,maskt,tmp3d)
       varname="u"
       call extract_bran2020(varname,2003,1,1,1,im,1,jm,1,km,lonu,latu,tmp1dz,masku,tmp3d)
       varname="v"
       call extract_bran2020(varname,2003,1,1,1,im,1,jm,1,km,lonv,latv,tmp1dz,maskv,tmp3d)
       do j=1,jm
          do i=1,im
             dept(i,j,1:km)=tmp1dz(1:km)
             depu(i,j,1:km)=tmp1dz(1:km)
             depv(i,j,1:km)=tmp1dz(1:km)
          end do
       end do
    else if(idat == 3)then !---FORA-NP60
       varname="t"
       call extract_fora_np60(varname,2003,1,1,1,im,1,jm,1,km,lont,latt,tmp1dz,maskt,tmp3d)
       varname="u"
       call extract_fora_np60(varname,2003,1,1,1,im,1,jm,1,km,lonu,latu,tmp1dz,masku,tmp3d)
       varname="v"
       call extract_fora_np60(varname,2003,1,1,1,im,1,jm,1,km,lonv,latv,tmp1dz,maskv,tmp3d)
       do j=1,jm
          do i=1,im
             dept(i,j,1:km)=tmp1dz(1:km)
             depu(i,j,1:km)=tmp1dz(1:km)
             depv(i,j,1:km)=tmp1dz(1:km)
          end do
       end do
    else if(idat == 4)then !---GLORYS12V1 (Regridded)
       call extract_grid_glorys12v1(1,im,1,jm,1,km,lont,latt,tmp1dz,maskt)
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
       call read_grid_jcope(lont,lonu,lonv,latt,latu,latv,maskt,masku,maskv,dept)
       depu(:,:,:)=dept(:,:,:)
       depv(:,:,:)=dept(:,:,:)
    end if

    deallocate(tmp1dz)
    deallocate(tmp3d)

  end subroutine read_grid

  !---------------------------------------------------------------------------------
  ! Extract analysis data |
  !---------------------------------------------------------------------------------

  !***To be modified
  subroutine extract_data(varname,idat,iyr,imon,iday,is,im,js,jm,ks,km,mean,sprd)

    use setting, only: datname
    use mod_read_lora, only: extract_anal
    use mod_read_bran2020,   only: extract_bran2020
    use mod_read_fora_np60,  only: extract_fora_np60
    use mod_read_glorys12v1, only: extract_glorys12v1
    use mod_read_jcope_fgo,  only: extract_jcope_fgo
    use mod_rmiss
    implicit none

    !---Common
    real(kind = 8),allocatable :: tmp1dx(:),tmp1dy(:),tmp1dz(:)
    real(kind = 8),allocatable :: tmp2d(:,:)

    !LORA
    integer imem !Dummy
    character(10) dir
    character(10) letkf
    character(10) region
    character(10) ms

    !---IN
    integer,intent(in) :: idat
    integer,intent(in) :: iyr,imon,iday
    integer,intent(in) :: is,im,js,jm,ks,km

    !GLORYS
    character(1),intent(in) :: varname

    !---OUT
    real(kind = 8),intent(out) :: mean(im,jm,km),sprd(im,jm,km)

    allocate(tmp1dx(im),tmp1dy(jm),tmp1dz(km))
    allocate(tmp2d(im,jm))

    if(idat == 1)then !---LORA-NP
       imem=0
       dir="NP"
       letkf="letkf"
       region="np"
       ms="mean"
       call extract_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,is,im,js,jm,ks,km,mean)
       ms="sprd"
       call extract_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,is,im,js,jm,ks,km,sprd)
    else if(idat == 2)then !---BRAN2020
       call extract_bran2020(varname,iyr,imon,iday,is,im,js,jm,ks,km,tmp1dx,tmp1dy,tmp1dz,tmp2d,mean)
       sprd(:,:,:)=rmiss
    else if(idat == 3)then !---FORA-NP60
       call extract_fora_np60(varname,iyr,imon,iday,is,im,js,jm,ks,km,tmp1dx,tmp1dy,tmp1dz,tmp2d,mean)
       sprd(:,:,:)=rmiss       
    else if(idat == 4)then !---GLORYS12V1
       call extract_glorys12v1(varname,iyr,imon,iday,is,im,js,jm,ks,km,tmp1dx,tmp1dy,tmp1dz,tmp2d,mean)
       sprd(:,:,:)=rmiss
    else if(idat == 5)then !---JCOPE-FGO
       call extract_jcope_fgo(varname,iyr,imon,iday,is,im,js,jm,ks,km,mean)
       sprd(:,:,:)=rmiss
    end if

    deallocate(tmp1dx,tmp1dy,tmp1dz)
    deallocate(tmp2d)

  end subroutine extract_data

end module mod_io
