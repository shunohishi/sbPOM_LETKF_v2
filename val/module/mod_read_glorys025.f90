module mod_read_glorys025

  integer,parameter :: im=1440,jm=681,km=75
  character(100),parameter :: g025_dir="/data/R/R2402/DATA/GLORYS/025"
  
contains
  
  !---------------------------------------------------------------------------
  ! Read GLORYS025 |
  !-----------------
  !
  ! Web: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_ENS_001_031/description
  ! DOI: https://doi.org/10.48670/moi-00024
  ! 
  ! - GLORYS2V4 from Mercator Ocean (Fr): DA => SEEK  (7 days)
  ! - ORAS5     from ECMWF              : DA => 3DVAR (5 days)
  ! - C-GLORSv7 from CMCC (It)          : DA => 3DVAR (7 days)
  !
  !---------------------------------------------------------------------------
  ! Select
  !
  ! datname:
  ! - glorys
  ! - oras5
  ! - cglors
  !
  ! varname:
  ! - t,s,u,v,h
  ! 
  !---------------------------------------------------------------------------

  subroutine get_glorys025_info(datname,varname,ncname)

    implicit none

    !---Common
    character(20) ncname1,ncname2
    
    !---IN
    character(10),intent(in) :: datname
    character(1),intent(in) :: varname
    
    !---OUT
    character(20) ncname

    !---Dataset name
    if(datname == "glorys")then
       ncname2="glor"
    else if(datname == "oras5")then
       ncname2="oras"
    else if(datname == "cglors")then
       ncname2="cglo"
    else
       write(*,*) "***Error: Incorrect datname => "//trim(datname)
    end if

    !---Variable name
    if(varname == "t")then
       ncname1="thetao"
    else if(varname == "s")then
       ncname1="so"
    else if(varname == "u")then
       ncname1="uo"
    else if(varname == "v")then
       ncname1="vo"
    else if(varname == "h")then
       ncname1="zos"
    else
       write(*,*) "***Error: Incorrect varname => "//trim(varname)
       stop
    end if

    ncname=trim(ncname1)//"_"//trim(ncname2)
    
  end subroutine get_glorys025_info

  !----------------------------------------------------------------------------
  
  subroutine read_glorys025(datname,varname,iyr,imon,iday,km_in,lon,lat,depth,mask,dat)

    use mod_rmiss
    use netcdf
    implicit none

    !---Parameter
    real(kind = 4),parameter :: dmiss=9.96921e36
    
    !---Common
    integer i,j,k,n
    integer status,access
    integer ncid,varid    

    real(kind = 4) tmp1dx(im),tmp1dy(jm),tmp1dz(km_in)
    real(kind = 4) tmp3d(im,jm,km_in)
    
    character(200) filename
    character(20) ncname
    character(4) yyyy
    character(2) mm,dd
    
    !---IN
    integer,intent(in) :: iyr,imon,iday
    integer,intent(in) :: km_in

    character(10),intent(in) :: datname
    character(1),intent(in)  :: varname 

    !---OUT
    real(kind = 8),intent(out) :: lon(im),lat(jm),depth(km_in)
    real(kind = 8),intent(out) :: mask(im,jm),dat(im,jm,km_in)

    !---Filename
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    filename=trim(g025_dir)//"/cmems_mod_glo_phy-all_my_0.25deg_P1D-m-"//yyyy//mm//dd//".nc" 

    status=access(trim(filename)," ")
    if(status == 0)then
       write(*,*) "Read :"//trim(filename)
    else
       write(*,*) "***Error: Not found "//trim(filename)
       stop
    end if
    
    !---Get ncname
    call get_glorys025_info(datname,varname,ncname)

    !---Read data
    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    status=nf90_inq_varid(ncid,"longitude",varid)
    status=nf90_get_var(ncid,varid,tmp1dx)

    status=nf90_inq_varid(ncid,"latitude",varid)
    status=nf90_get_var(ncid,varid,tmp1dy)

    status=nf90_inq_varid(ncid,"depth",varid)
    status=nf90_get_var(ncid,varid,tmp1dz,(/1/),(/km_in/))

    if(varname == "h")then
       status=nf90_inq_varid(ncid,trim(ncname),varid)
       status=nf90_get_var(ncid,varid,tmp3d(:,:,1),(/1,1/),(/im,jm/))
    else
       status=nf90_inq_varid(ncid,trim(ncname),varid)
       status=nf90_get_var(ncid,varid,tmp3d,(/1,1,1/),(/im,jm,km_in/))
    end if

    status=nf90_close(ncid)
    
    !---Post process
    !Longitude
    n=0
    do i=1,im
       if(0.e0 <= tmp1dx(i))then
          n=n+1
          lon(n)=dble(tmp1dx(i))
       end if
    end do


    do i=1,im
       if(tmp1dx(i) < 0.e0)then
          n=n+1
          lon(n)=dble(tmp1dx(i))+360.d0
       end if
    end do

    !Latitude
    lat(:)=dble(tmp1dy(:))

    !Depth
    depth(:)=dble(tmp1dz(:))

    !Mask
    k=1
    do j=1,jm

       n=0
       
       do i=1,im
          if(0.e0 <= tmp1dx(i))then
             n=n+1
             if(tmp3d(i,j,k) == dmiss)then
                mask(n,j)=0.d0
             else
                mask(n,j)=1.d0
             end if
          end if
       end do

       do i=1,im
          if(tmp1dx(i) < 0.e0)then
             n=n+1
             if(tmp3d(i,j,k) == dmiss)then
                mask(n,j)=0.d0
             else
                mask(n,j)=1.d0
             end if
          end if
       end do
    end do
    
    !Data
    do k=1,km_in
       do j=1,jm

          n=0
          
          do i=1,im
             if(0.e0 <= tmp1dx(i))then
                n=n+1
                dat(n,j,k)=dble(tmp3d(i,j,k))
             end if
          end do

          do i=1,im
             if(tmp1dx(i) < 0.e0)then
                n=n+1
                dat(n,j,k)=dble(tmp3d(i,j,k))
             end if
          end do !i
          
       end do    !j
    end do       !k

    !Missing value
    do j=1,jm
       do i=1,im
          if(mask(i,j) == 0.d0)then
             dat(i,j,:)=rmiss
          end if
       end do
    end do
    
  end subroutine read_glorys025
    
end module mod_read_glorys025
