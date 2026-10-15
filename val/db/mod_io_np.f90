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
     else if(idat == 4)then !---GLORYS12V1
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
  
  !---------------------------------------------------------------------------------
  ! Read Grid data |
  !---------------------------------------------------------------------------------

  ! ***To be modified ***
  subroutine read_grid(idat,im,jm,km,lont,lonu,lonv,latt,latu,latv,maskt,masku,maskv)

    use setting, only: datname
    use mod_read_lora,       only: read_grid_lora => read_grid
    use mod_read_bran2020,   only: read_bran2020
    use mod_read_fora_np60,  only: read_fora_np60
    use mod_read_glorys12v1, only: read_grid_glorys12v1 => read_grid
    use mod_read_jcope_fgo,  only: read_grid_jcope => read_grid
    implicit none

    !---Common
    real(kind = 8),allocatable :: tmp1dx(:),tmp1dy(:),tmp1dz(:)
    real(kind = 8),allocatable :: tmp2d(:,:)
    real(kind = 8),allocatable :: tmp3d(:,:,:)

    !LORA
    character(10) dir

    !Others
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

    if(idat == 1)then !--- LORA-NP
       dir="NP"
       call read_grid_lora(dir,lont,lonu,lonv, &
            & latt,latu,latv, &
            & tmp3d,tmp3d,tmp3d,tmp3d, &
            & maskt,masku,maskv)
    else if(idat == 2)then !--- BRAN2020
       varname="t"
       call read_bran2020(varname,2003,1,1,km,lont,latt,tmp1dz,maskt,tmp3d)
       varname="u"
       call read_bran2020(varname,2003,1,1,km,lonu,latu,tmp1dz,masku,tmp3d)
       varname="v"
       call read_bran2020(varname,2003,1,1,km,lonv,latv,tmp1dz,maskv,tmp3d)
    else if(idat == 3)then !---FORA-NP60
       varname="t"
       call read_fora_np60(varname,2003,1,1,km,lont,latt,tmp1dz,maskt,tmp3d)
       varname="u"
       call read_fora_np60(varname,2003,1,1,km,lonu,latu,tmp1dz,masku,tmp3d)
       varname="v"
       call read_fora_np60(varname,2003,1,1,km,lonv,latv,tmp1dz,maskv,tmp3d)
    else if(idat == 4)then !---GLORYS12V1 (Probably Regridded)
       call read_grid_glorys12v1(tmp1dx,tmp1dy,tmp1dz,tmp2d)
       lont(:)=tmp1dx(:)
       lonu(:)=tmp1dx(:)
       lonv(:)=tmp1dx(:)
       latt(:)=tmp1dy(:)
       latu(:)=tmp1dy(:)
       latv(:)=tmp1dy(:)
       maskt(:,:)=tmp2d(:,:)
       masku(:,:)=tmp2d(:,:)
       maskv(:,:)=tmp2d(:,:)
    else if(idat == 5)then !---JCOPE-FGO
       call read_grid_jcope(lont,lonu,lonv,latt,latu,latv,maskt,masku,maskv,tmp3d)
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
    use mod_read_lora,       only: read_anal
    use mod_read_bran2020,   only: read_bran2020
    use mod_read_fora_np60,  only: read_fora_np60
    use mod_read_glorys12v1, only: read_glorys12v1
    use mod_read_jcope_fgo,  only: read_jcope_fgo
    use mod_rmiss
    implicit none

    !---Common
    integer k

    real(kind = 8),allocatable :: tmp1dx(:),tmp1dy(:),tmp1dz(:)
    real(kind = 8),allocatable :: tmp2d(:,:)
    real(kind = 8),allocatable :: tmp3d(:,:,:)

    !LORA
    integer imem !Dummy
    character(10) dir
    character(10) letkf
    character(10) region
    character(10) ms

    !Others
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
    allocate(tmp3d(im,jm,km))

    if(idat == 1)then !---LORA-NP
       imem=0
       k=1
       dir="NP"
       letkf="letkf"
       region="np"       
       ms="mean"
       call read_anal(dir,letkf,region,ms,imem,"t",iyr,imon,iday,im,jm,k,maskt,t)
       call read_anal(dir,letkf,region,ms,imem,"u",iyr,imon,iday,im,jm,k,masku,u)
       call read_anal(dir,letkf,region,ms,imem,"v",iyr,imon,iday,im,jm,k,maskv,v)
       ms="sprd"
       call read_anal(dir,letkf,region,ms,imem,"t",iyr,imon,iday,im,jm,k,maskt,tsprd)
       call read_anal(dir,letkf,region,ms,imem,"u",iyr,imon,iday,im,jm,k,masku,usprd)
       call read_anal(dir,letkf,region,ms,imem,"v",iyr,imon,iday,im,jm,k,maskv,vsprd)
    else if(idat == 2)then !---BRAN2020
       k=1
       varname="t"       
       call read_bran2020(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,t)
       varname="u"       
       call read_bran2020(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,u)
       varname="v"       
       call read_bran2020(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,v)
       tsprd=rmiss
       usprd=rmiss
       vsprd=rmiss
    else if(idat == 3)then !---FORA-NP60
       k=1
       varname="t"
       call read_fora_np60(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,t)
       varname="u"
       call read_fora_np60(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,u)
       varname="v"
       call read_fora_np60(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,v)
       tsprd=rmiss
       usprd=rmiss
       vsprd=rmiss       
    else if(idat == 4)then !---GLORYS12V1
       k=1
       varname="t"
       call read_glorys12v1(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,t)
       varname="u"
       call read_glorys12v1(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,u)
       varname="v"
       call read_glorys12v1(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,v)
       tsprd=rmiss
       usprd=rmiss
       vsprd=rmiss
    else if(idat == 5)then !---JCOPE-FGO
       k=1
       varname="t"
       call read_jcope_fgo(varname,iyr,imon,iday,k,maskt,t)
       varname="u"
       call read_jcope_fgo(varname,iyr,imon,iday,k,masku,u)
       varname="v"
       call read_jcope_fgo(varname,iyr,imon,iday,k,maskv,v)
       tsprd=rmiss
       usprd=rmiss
       vsprd=rmiss
    end if

    deallocate(tmp1dx,tmp1dy,tmp1dz)
    deallocate(tmp2d)
    deallocate(tmp3d)

  end subroutine read_surface_data

end module mod_io
