module mod_io

contains

  !---------------------------------------------------------------------------------
  ! Grid size (Analysis data) |
  !---------------------------------------------------------------------------------

  ! *** To be modified ***
  subroutine get_grid_size(idat,im,jm,km)

    use mod_gridinfo, im_lora => im, jm_lora => jm, km_lora => km
    use mod_read_glorys025,  only: im_g025 => im, jm_g025 => jm, km_g025 => km
    implicit none

    !---IN
    integer,intent(in) :: idat

    !---OUT
    integer,intent(out) :: im,jm,km

    if(idat == 1)then !---LORA-QG
       im=im_lora
       jm=jm_lora
       km=km_lora
    else if(idat == 2 .or. idat == 3 .or. idat == 4)then !---GLORYS025/ORAS5/C-GLORS
       im=im_g025
       jm=jm_g025
       km=km_g025
    else
       write(*,*) "***Error: Incorecot idat_a => ",idat
       stop
    end if

  end subroutine get_grid_size

  !-----------------------------------------------------------------------
  ! Read grid data |
  !-----------------------------------------------------------------------

  !*** To be modified
  subroutine read_grid(idat,im,jm,km,lont,lonu,lonv,latt,latu,latv,dept,depu,depv,maskt,masku,maskv)

    use setting, only: datname
    use mod_read_lora, only: read_grid_lora => read_grid
    use mod_read_glorys025, only: extract_glorys025
    implicit none

    !---Common
    integer i,j

    real(kind = 8),allocatable :: tmp3d(:,:,:)

    !LORA
    character(10) :: dir="QGLOBAL"

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

    allocate(tmp3d(im,jm,km))

    if(idat == 1)then
       dir="QGLOBAL"
       call read_grid_lora(dir,lont,lonu,lonv, &
            & latt,latu,latv, &
            & dept,depu,depv,tmp3d, &
            & maskt,masku,maskv)
    else if(idat == 2 .or. idat == 3 .or. idat == 4)then
       varname="t"
       call extract_glorys025(datname(idat),varname,2003,1,1,1,im,1,jm,1,km,lont,latt,dept(1,1,:),maskt,tmp3d)
       lonu(:)=lont(:)
       lonv(:)=lont(:)
       latu(:)=latt(:)
       latv(:)=latt(:)
       do j=1,jm
          do i=1,im
             dept(i,j,:)=dept(1,1,:)
             depu(i,j,:)=dept(1,1,:)
             depv(i,j,:)=dept(1,1,:)
          end do
       end do
       masku(:,:)=maskt(:,:)
       maskv(:,:)=maskt(:,:)
    end if

    deallocate(tmp3d)

  end subroutine read_grid

  !---------------------------------------------------------------------------------
  ! Extract analysis data |
  !---------------------------------------------------------------------------------

  !***To be modified
  subroutine extract_data(varname,idat,iyr,imon,iday,is,im,js,jm,ks,km,mean,sprd)

    use setting, only: datname
    use mod_read_lora, only: extract_anal
    use mod_read_glorys025, only: extract_glorys025
    use mod_rmiss
    implicit none

    !---Common
    real(kind = 8),allocatable :: tmp1dx(:),tmp1dy(:),tmp1dz(:)
    real(kind = 8),allocatable :: tmp2d(:,:)

    !LORA
    integer imem !Dummy
    character(10) :: dir="QGLOBAL"
    character(10) :: letkf="letkf"
    character(10) :: region="qglobal"
    character(10) :: ms

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

    if(idat == 1)then
       imem=0
       ms="mean"
       call extract_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,is,im,js,jm,ks,km,mean)
       ms="sprd"
       call extract_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,is,im,js,jm,ks,km,sprd)
    else if(idat == 2 .or. idat == 3 .or. idat == 4)then
       call extract_glorys025(datname(idat),varname, &
            & iyr,imon,iday,is,im,js,jm,ks,km,tmp1dx,tmp1dy,tmp1dz,tmp2d,mean)
       sprd=rmiss
    end if

    deallocate(tmp1dx,tmp1dy,tmp1dz)
    deallocate(tmp2d)

  end subroutine extract_data

end module mod_io
