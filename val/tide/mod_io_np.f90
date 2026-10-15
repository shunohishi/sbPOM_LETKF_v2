module mod_io

contains


  !---------------------------------------------------------------------------------
  ! Data name |
  !---------------------------------------------------------------------------------

  ! *** To be modified ***  
  subroutine get_datname(idat,datname)

    use setting, only: ndat_a
    implicit none

    !---IN
    integer,intent(in) :: idat

    !---OUT
    character(10),intent(out) :: datname

    if(idat == 1) datname="lora      "
    if(idat == 2) datname="bran      "
    if(idat == 3) datname="fora      "
    if(idat == 4) datname="glorys    "
    if(idat == 5) datname="jcope     "
    if(idat > ndat_a)then
       write(*,*) "***Error: idat_a => ",idat
       stop
    end if

  end subroutine get_datname

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

  !---------------------------------------------------------------------------------
  ! Read grid data |
  !---------------------------------------------------------------------------------

  ! *** To be modified ***  
  subroutine read_grid(idat,im,jm,km,lon,lat,mask)

    use mod_read_lora, only: read_grid_lora => read_grid
    use mod_read_bran2020,   only: read_bran2020
    use mod_read_fora_np60,  only: read_fora_np60
    use mod_read_glorys12v1, only: read_grid_glorys12v1 => read_grid
    use mod_read_jcope_fgo,  only: read_grid_jcope => read_grid    
    implicit none

    !---Common
    integer i,j

    real(kind = 8),allocatable :: tmp1dx(:),tmp1dy(:),tmp1dz(:)
    real(kind = 8),allocatable :: tmp2d(:,:)
    real(kind = 8),allocatable :: tmp3d(:,:,:)

    !LORA
    character(10) dir

    !GLORYS
    character(10) datname
    character(1) varname

    !---IN
    integer,intent(in) :: idat
    integer,intent(in) :: im,jm,km

    !---OUT
    real(kind = 8),intent(out) :: lon(im),lat(jm),mask(im,jm)

    allocate(tmp1dx(im),tmp1dy(jm),tmp1dz(km))
    allocate(tmp2d(im,jm))
    allocate(tmp3d(im,jm,km))

    call get_datname(idat,datname)

    if(idat == 1)then !---LORA-NP
       dir="NP"
       call read_grid_lora(dir,lon,tmp1dx,tmp1dx, &
            & lat,tmp1dy,tmp1dy, &
            & tmp3d,tmp3d,tmp3d,tmp3d, &
            & mask,tmp2d,tmp2d)
    else if(idat == 2)then !--- BRAN
       varname="t"
       call read_bran2020(varname,2003,1,1,km,lon,lat,tmp1dz,mask,tmp3d)
    else if(idat == 3)then !---FORA-NP60
       varname="t"
       call read_fora_np60(varname,2003,1,1,km,lon,lat,tmp1dz,mask,tmp3d)
    else if(idat == 4)then !---GLORYS12V1 (Probably Regridded)
       call read_grid_glorys12v1(lon,lat,tmp1dz,mask)
    else if(idat == 5)then !---JCOPE-FGO
       call read_grid_jcope(lon,tmp1dx,tmp1dx,lat,tmp1dy,tmp1dy,mask,tmp2d,tmp2d,tmp3d)
    end if

    deallocate(tmp1dx,tmp1dy,tmp1dz)
    deallocate(tmp2d)
    deallocate(tmp3d)

  end subroutine read_grid

  !---------------------------------------------------------------------------------
  ! Read SSH analysis data |
  !---------------------------------------------------------------------------------

  ! *** To be modified ***  
  subroutine read_ssh_anal(idat,iyr,imon,iday,im,jm,km,mask,dat,sprd)

    use mod_read_lora, only: read_anal
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

    !LORA
    integer imem !Dummy
    character(10) dir
    character(10) letkf
    character(10) region
    character(10) ms

    !GLORYS
    character(10) datname
    character(1) varname

    !---IN
    integer,intent(in) :: idat
    integer,intent(in) :: iyr,imon,iday
    integer,intent(in) :: im,jm,km

    real(kind = 8),intent(in) :: mask(im,jm)

    !---OUT
    real(kind = 8),intent(out) :: dat(im,jm)
    real(kind = 8),intent(out) :: sprd(im,jm)

    allocate(tmp1dx(im),tmp1dy(jm),tmp1dz(km))
    allocate(tmp2d(im,jm))

    call get_datname(idat,datname)    

    if(idat == 1)then !---LORA-NP
       imem=0
       k=1
       dir="NP"
       letkf="letkf"
       region="np"
       ms="mean"
       call read_anal(dir,letkf,region,ms,imem,"el",iyr,imon,iday,im,jm,k,mask,dat)
       ms="sprd"
       call read_anal(dir,letkf,region,ms,imem,"el",iyr,imon,iday,im,jm,k,mask,sprd)
    else if(idat == 2)then !---BRAN2020
       k=1
       varname="h"       
       call read_bran2020(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,dat)
       sprd=rmiss
    else if(idat == 3)then !---FORA-NP60
       k=1
       varname="h"
       call read_fora_np60(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,dat)
       sprd=rmiss
    else if(idat == 4)then !---GLORYS12V1
       k=1
       varname="h"
       call read_glorys12v1(varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,dat)       
       sprd=rmiss
    else if(idat == 5)then !---JCOPE-FGO
       k=1
       varname="h"
       call read_jcope_fgo(varname,iyr,imon,iday,k,mask,dat)       
       sprd=rmiss
    end if

    deallocate(tmp1dx,tmp1dy,tmp1dz)
    deallocate(tmp2d)

  end subroutine read_ssh_anal

end module mod_io
