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

  !---------------------------------------------------------------------------------
  ! Read Grid data |
  !---------------------------------------------------------------------------------

  ! ***To be modified ***
  subroutine read_grid(idat,im,jm,km,lont,lonu,lonv,latt,latu,latv,maskt,masku,maskv)

    use setting, only: datname
    use mod_read_lora, only: read_grid_lora => read_grid
    use mod_read_glorys025, only: read_glorys025
    implicit none

    !---Common
    real(kind = 8),allocatable :: tmp1dx(:),tmp1dy(:),tmp1dz(:)
    real(kind = 8),allocatable :: tmp2d(:,:)
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
    real(kind = 8),intent(out) :: maskt(im,jm),masku(im,jm),maskv(im,jm)

    allocate(tmp1dx(im),tmp1dy(jm),tmp1dz(km))
    allocate(tmp2d(im,jm))
    allocate(tmp3d(im,jm,km))

    if(idat == 1)then
       dir="QGLOBAL"
       call read_grid_lora(dir,lont,lonu,lonv, &
            & latt,latu,latv, &
            & tmp3d,tmp3d,tmp3d,tmp3d, &
            & maskt,masku,maskv)
    else if(idat == 2 .or. idat == 3 .or. idat == 4)then
       varname="t"
       call read_glorys025(datname(idat),varname,2003,1,1,km,tmp1dx,tmp1dy,tmp1dz,tmp2d,tmp3d)
       lont(:)=tmp1dx(:)
       lonu(:)=tmp1dx(:)
       lonv(:)=tmp1dx(:)
       latt(:)=tmp1dy(:)
       latu(:)=tmp1dy(:)
       latv(:)=tmp1dy(:)
       maskt(:,:)=tmp2d(:,:)
       masku(:,:)=tmp2d(:,:)
       maskv(:,:)=tmp2d(:,:)
    end if

    deallocate(tmp1dx,tmp1dy,tmp1dz)
    deallocate(tmp2d)
    deallocate(tmp3d)

  end subroutine read_grid

  !---------------------------------------------------------------------------------
  ! Read sea-surface analysis data |
  !---------------------------------------------------------------------------------

  !*** To be modified ***
  subroutine read_surface_data(idat,iyr,imon,iday,im,jm,km,maskt,masku,maskv,t,u,v,tsprd,usprd,vsprd)

    use setting, only: datname
    use mod_read_lora, only: read_anal
    use mod_read_glorys025, only: read_glorys025
    use mod_rmiss
    implicit none

    !---Common
    integer k

    real(kind = 8),allocatable :: tmp1dx(:),tmp1dy(:),tmp1dz(:)
    real(kind = 8),allocatable :: tmp2d(:,:)

    !LORA
    integer imem !Dummy
    character(10) dir
    character(10) letkf
    character(10) region
    character(10) ms

    !GLORYS
    character(1) varname

    !---IN
    integer,intent(in) :: idat
    integer,intent(in) :: iyr,imon,iday
    integer,intent(in) :: im,jm,km

    real(kind = 8),intent(in) :: maskt(im,jm),masku(im,jm),maskv(im,jm)

    !---OUT
    real(kind = 8),intent(out) :: t(im,jm),u(im,jm),v(im,jm)
    real(kind = 8),intent(out) :: tsprd(im,jm),usprd(im,jm),vsprd(im,jm)

    allocate(tmp1dx(im),tmp1dy(jm),tmp1dz(km))
    allocate(tmp2d(im,jm))

    if(idat == 1)then
       imem=0
       k=1
       dir="QGLOBAL"
       letkf="letkf"
       region="qglobal"       
       ms="mean"
       call read_anal(dir,letkf,region,ms,imem,"t",iyr,imon,iday,im,jm,k,maskt,t)
       call read_anal(dir,letkf,region,ms,imem,"u",iyr,imon,iday,im,jm,k,masku,u)
       call read_anal(dir,letkf,region,ms,imem,"v",iyr,imon,iday,im,jm,k,maskv,v)
       ms="sprd"
       call read_anal(dir,letkf,region,ms,imem,"t",iyr,imon,iday,im,jm,k,maskt,tsprd)
       call read_anal(dir,letkf,region,ms,imem,"u",iyr,imon,iday,im,jm,k,masku,usprd)
       call read_anal(dir,letkf,region,ms,imem,"v",iyr,imon,iday,im,jm,k,maskv,vsprd)     
    else if(idat == 2 .or. idat == 3 .or. idat == 4)then
       k=1
       varname="t"
       call read_glorys025(datname(idat),varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,t)
       varname="u"
       call read_glorys025(datname(idat),varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,u)
       varname="v"
       call read_glorys025(datname(idat),varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,v)
       tsprd=rmiss
       usprd=rmiss
       vsprd=rmiss
    end if

    deallocate(tmp1dx,tmp1dy,tmp1dz)
    deallocate(tmp2d)

  end subroutine read_surface_data

end module mod_io
