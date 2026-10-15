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
    if(idat == 2) datname="glorys"
    if(idat == 3) datname="oras5"
    if(idat == 4) datname="cglors"
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
    use mod_read_glorys025, im_g025 => im, jm_g025 => jm, km_g025 => km
    implicit none

    !---IN
    integer,intent(in) :: idat

    !---OUT
    integer,intent(out) :: im,jm,km

    if(idat == 1)then !---LORA-QGLOBAL
       im=im_lora
       jm=jm_lora
       km=km_lora
    else if(idat == 2 .or. idat == 3 .or. idat == 4)then !GLORYS/CGLORS/ORAS5
       im=im_g025
       jm=jm_g025
       km=km_g025
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
    use mod_read_glorys025, only: read_glorys025
    implicit none

    !---Common
    integer i,j,k

    real(kind = 8),allocatable :: tmp1dx(:),tmp1dy(:),tmp1dz(:)
    real(kind = 8),allocatable :: tmp2d(:,:)
    real(kind = 8),allocatable :: tmp3d(:,:,:)

    !LORA
    character(10) dir

    !Others
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

    if(idat == 1)then !---LORA-QGLOBAL
       dir="QGLOBAL"
       call read_grid_lora(dir,lon,tmp1dx,tmp1dx, &
            & lat,tmp1dy,tmp1dy, &
            & tmp3d,tmp3d,tmp3d,tmp3d, &
            & mask,tmp2d,tmp2d)
    else if(idat == 2 .or. idat == 3 .or. idat == 4)then
       varname="t"
       call read_glorys025(datname,varname,2003,1,1,km,tmp1dx,tmp1dy,tmp1dz,tmp2d,tmp3d)
       lon(:)=tmp1dx(:)
       lat(:)=tmp1dy(:)
       mask(:,:)=tmp2d(:,:)
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

    if(idat == 1)then !---LORA-QGLOBAL
       imem=0
       k=1
       dir="QGLOBAL"
       letkf="letkf"
       region="qglobal"
       ms="mean"
       call read_anal(dir,letkf,region,ms,imem,"el",iyr,imon,iday,im,jm,k,mask,dat)
       ms="sprd"
       call read_anal(dir,letkf,region,ms,imem,"el",iyr,imon,iday,im,jm,k,mask,sprd)
    else if(idat == 2 .or. idat == 3 .or. idat == 4)then
       k=1
       varname="h"
       call read_glorys025(datname,varname,iyr,imon,iday,k,tmp1dx,tmp1dy,tmp1dz,tmp2d,dat)
       sprd=rmiss
    end if

    deallocate(tmp1dx,tmp1dy,tmp1dz)
    deallocate(tmp2d)

  end subroutine read_ssh_anal

end module mod_io
