module mod_read_era5

  integer,parameter :: im=1440,jm=721
  character(200),parameter :: era5_dir="/data/R/R2402/DATA"
  
contains

  !--------------------------------------------------------
  ! Read ERA5 |
  !--------------------------------------------------------
  ! NOTE:
  ! - Horizontal resolution: 0.25 x 0.25 degree
  ! - Temporal resolution: 1 hour
  ! - Variable:
  !     10m zonal wind                       [10m_u_component_of_wind]
  !     10m meridional wind                  [10m_v_component_of_wind]
  !     2m air temperature                   [2m_temperature]
  !     2m dewpoint temperature              [2m_dewpoint_temperature]
  !     Surface downward longwave radiation  [mean_surface_downward_long_wave_radiation_flux]
  !     Surface downward shortwave radiation [mean_surface_downward_short_wave_radiation_flux]
  !     Surface pressure                     [surface_pressure]
  !     Sea level pressure                   [mean_sea_level_pressure]
  !     Total precipitation                  [total_precipitation]
  !     
  !     2m specific humidity
  !     => Estimated from 2m air and dewpoint temperatures and surface pressure
  !
  !--------------------------------------------------------
  ! Created by S.Ohishi @ 2025.08
  !
  !--------------------------------------------------------

  subroutine get_info(ivar,varname,ncname,add,mult)

    implicit none

    !---IN
    integer,intent(in) :: ivar

    !---OUT
    real(kind = 8) add,mult
    
    character(100),intent(out) :: varname,ncname

    if(ivar == 1)then ![m/s]
       varname="10m_u_component_of_wind"
       ncname="u10"
       add=0.d0
       mult=1.d0
    else if(ivar == 2)then ![m/s]
       varname="10m_v_component_of_wind"
       ncname="v10"
       add=0.d0
       mult=1.d0
    else if(ivar == 3)then ![K] --> [degree C]
       varname="2m_temperature"
       ncname="t2m"
       add=-273.15d0
       mult=1.d0
    else if(ivar == 4)then !*** Keep [K] ***
       varname="2m_dewpoint_temperature"
       ncname="d2m"
       add=0.d0
       mult=1.d0
    else if(ivar == 5)then ![W/m2]
       varname="mean_surface_downward_long_wave_radiation_flux"
       ncname="avg_sdlwrf"
       add=0.d0
       mult=1.d0
    else if(ivar == 6)then ![W/m2]
       varname="mean_surface_downward_short_wave_radiation_flux"
       ncname="avg_sdswrf"
       add=0.d0
       mult=1.d0
    else if(ivar == 7)then ![Pa]
       varname="surface_pressure"
       ncname="sp"
       add=0.d0
       mult=1.d0
    else if(ivar == 8)then ![Pa]
       varname="mean_sea_level_pressure"
       ncname="msl"
       add=0.d0
       mult=1.d0
    else if(ivar == 9)then ![m/hour] --> [mm/day]
       varname="total_precipitation"
       ncname="tp"
       add=0.d0
       mult=1000.d0*24.d0
    end if
    
  end subroutine get_info
  
  !---------------------------------------------------------

  subroutine read_grid(lon,lat,land)

    use netcdf
    implicit none

    !---Common
    integer i,j
    integer status,access
    integer ncid,varid

    real(kind = 8) tmp1dy(jm)    
    real(kind = 4) tmp(im,jm)
    
    character(200) filename
    
    !---OUT
    real(kind = 8),intent(out) :: lon(im),lat(jm)
    real(kind = 8),intent(out) :: land(im,jm) !1:Land, 0:Ocean

    filename=trim(era5_dir)//"/ERA5/land_sea_mask.nc"

    status=access(trim(filename)," ")
    if(status == 0)then
       write(*,*) "Read: "//trim(filename)
    else
       write(*,*) "***Error: Not found "//trim(filename)
       stop
    end if

    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    status=nf90_inq_varid(ncid,"longitude",varid)
    status=nf90_get_var(ncid,varid,lon)

    status=nf90_inq_varid(ncid,"latitude",varid)
    status=nf90_get_var(ncid,varid,tmp1dy)
    
    status=nf90_inq_varid(ncid,"lsm",varid)
    status=nf90_get_var(ncid,varid,tmp,(/1,1,1/),(/im,jm,1/))
    
    status=nf90_close(ncid)

    do j=1,jm
       lat(jm-j+1)=tmp1dy(j)
    end do
    
    do j=1,jm
       do i=1,im
          if(tmp(i,j) == 0.e0)then !Ocean
             land(i,jm-j+1)=0.d0
          else !Intermediate or Land
             land(i,jm-j+1)=1.d0
          end if
       end do
    end do

  end subroutine read_grid

  !---------------------------------------------------------
  
  subroutine read_era5(iyr,imon,iday,ihour,u,v,ta,qa,lw,sw,slp,prep)

    use mod_rmiss
    use netcdf
    implicit none

    !---Parameter
    integer,parameter :: nvar=9
    real(kind = 4),parameter :: dmiss=3.40282346638529e38
    
    !---Common
    integer i,j
    integer ivar
    integer status,access
    integer ncid,varid

    real(kind = 4) tmp(im,jm)
    
    real(kind = 8) add,mult
    real(kind = 8) dat(im,jm)
    real(kind = 8) td(im,jm),p(im,jm)

    character(200) filename
    character(100) varname,ncname
    character(8) yyyymmdd
    character(4) yyyy
    character(2) mm,dd
    
    !---IN
    integer,intent(in) :: iyr,imon,iday,ihour

    !---OUT
    real(kind = 8),intent(out) :: u(im,jm),v(im,jm)
    real(kind = 8),intent(out) :: ta(im,jm),qa(im,jm)
    real(kind = 8),intent(out) :: lw(im,jm),sw(im,jm)
    real(kind = 8),intent(out) :: slp(im,jm),prep(im,jm)

    !---Date
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    yyyymmdd=yyyy//mm//dd
    
    do ivar=1,nvar

       !---Varname
       call get_info(ivar,varname,ncname,add,mult)
       
       !---Filename
       filename=trim(era5_dir)//"/ERA5/"//trim(varname)//"_"//yyyymmdd//".nc"

       status=access(trim(filename)," ")
       if(status == 0)then
          write(*,*) "Read: "//trim(filename)
       else
          write(*,*) "***Error: Not found "//trim(filename)
          stop
       end if
       
       !---Read data
       !Open NetCDF
       status=nf90_open(trim(filename),nf90_nowrite,ncid)
       
       status=nf90_inq_varid(ncid,trim(ncname),varid)
       status=nf90_get_var(ncid,varid,tmp,(/1,1,ihour+1/),(/im,jm,1/)) !ihour = 1 at 0000UTC
       
       status=nf90_close(ncid)

       !---Post process
       do j=1,jm
          do i=1,im
             if(tmp(i,j) == dmiss)then
                dat(i,jm-j+1)=rmiss
             else
                dat(i,jm-j+1)=mult*dble(tmp(i,j))+add
             end if
          end do
       end do
       
       if(ivar == 1) u(:,:)=dat(:,:)
       if(ivar == 2) v(:,:)=dat(:,:)
       if(ivar == 3) ta(:,:)=dat(:,:)
       if(ivar == 4) td(:,:)=dat(:,:)
       if(ivar == 5) lw(:,:)=dat(:,:)
       if(ivar == 6) sw(:,:)=dat(:,:)
       if(ivar == 7) p(:,:)=dat(:,:)
       if(ivar == 8) slp(:,:)=dat(:,:)
       if(ivar == 9) prep(:,:)=dat(:,:)
       
    end do !ivar

    call specific_humidity_2m(td,p,qa)
    
  end subroutine read_era5
  
  !--------------------------------------------------------------------------------
  ! Specific humidity at 2m |
  !--------------------------
  !
  ! IFS documentation CY41R2
  ! Chapter 7: Clouds and large-scale precipitation
  ! Eqs. (7.4) and (7.5) 
  !
  !--------------------------------------------------------------------------------
  
  subroutine specific_humidity_2m(td,p,qa)

    implicit none

    !---Parameter
    real(kind = 8),parameter :: Rdry=287.0597d0,Rvap=461.5250d0 ![J/kg/K]
    real(kind = 8),parameter :: a1=611.21d0 ![Pa]
    real(kind = 8),parameter :: a3=17.502d0 ![-]
    real(kind = 8),parameter :: a4=32.19d0  ![K]
    real(kind = 8),parameter :: T0=273.16d0 ![K]
    
    !---Common
    integer i,j

    real(kind = 8) e
    
    !---IN
    real(kind = 8),intent(in) :: td(im,jm)  ![K]
    real(kind = 8),intent(in) :: p(im,jm)   ![Pa]

    !---OUT
    real(kind = 8),intent(out) :: qa(im,jm) ![-]

    do j=1,jm
       do i=1,im
          e=a1*exp(a3*(td(i,j)-T0)/(td(i,j)-a4)) ![Pa]
          qa(i,j)=(Rdry/Rvap)*e/(p(i,j)-(1-Rdry/Rvap)*e) ![-]
       end do
    end do

  end subroutine specific_humidity_2m
  
end module mod_read_era5
