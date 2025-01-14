module mod_read_jra55

  integer,parameter :: im=288,jm=145
  
contains
  
  !----------------------------------------------------------------
  ! Make position (1.25 degree version)
  !----------------------------------------------------------------
 
  subroutine read_grid(long,lati)
    
    implicit none
    
    integer i,j
    
    real(kind = 8),intent(out) :: long(im),lati(jm)

    do i=1,im
       long(i)=0.d0+dble(i-1)*1.25d0
    end do

    do j=1,jm
       lati(j)=-90.d0+dble(j-1)*1.25d0
    end do
    
  end subroutine read_grid

  !----------------------------------------------------------------
  ! Read data (LL125/anl_surf125/fcst_surf125/fcst_phys2m125)|
  !-----------------------------------------------------------
  !
  ! 2022.10.06 S.Ohishi modified read_land
  ! 
  !----------------------------------------------------------------
  
  !1: Land, 0: Ocean
  subroutine read_land(iunit,land)

    implicit none

    integer status,system,access
    integer i,j
    
    real(kind = 4) dland(im,jm)

    character(100) filename

    !IN
    integer,intent(in) :: iunit

    !OUT
    real(kind = 8),intent(out) :: land(im,jm)
    
    filename="/data/R/R2402/DATA/JRA55/Const/LL125.grib"

    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not found "//trim(filename) 
       stop
    end if

    status=system("wgrib "//trim(filename)//" | grep :LAND: | wgrib "//trim(filename)//" -i -ieee -nh -o land.bin")

    open(unit=iunit,file="land.bin",access="direct",form="unformatted",recl=4*im*jm)
    read(iunit,rec=1) ((dland(i,j),i=1,im),j=1,jm)
    close(iunit)
    status=system("rm -f land.bin")

    do j=1,jm
       do i=1,im

          if(0.e0 < dland(i,j))then
             land(i,jm-j+1)=1.d0
          else
             land(i,jm-j+1)=0.d0
          end if

          !if(dland(i,j) == 1.)then
          !   land(i,jm-j+1)=dland(i,j)
          !else
          !   land(i,jm-j+1)=0.
          !end if

       end do
    end do

  end subroutine read_land

  !--------------------------

  subroutine read_geo(iunit,geo)

    use mod_rmiss
    implicit none

    real(kind = 4),parameter :: dmiss=9.999e20

    integer status,system,access
    integer i,j
    
    real(kind = 4) dgeo(im,jm)

    character(100) filename

    !IN
    integer,intent(in) :: iunit

    !OUT
    real(kind = 8),intent(out) :: geo(im,jm)

    filename="/data/R/R2402/DATA/JRA55/Const/LL125.grib"

    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not found "//trim(filename) 
       stop
    end if

    status=system("wgrib "//trim(filename)//" | grep :GP: | wgrib "//trim(filename)//" -i -ieee -nh -o geo.bin")

    open(unit=iunit,file="geo.bin",access="direct",form="unformatted",recl=4*im*jm)
    read(iunit,rec=1) ((dgeo(i,j),i=1,im),j=1,jm)
    close(iunit)
    status=system("rm -f geo.bin")

    do j=1,jm
       do i=1,im
          if(dgeo(i,j) == dmiss)then
             geo(im,jm-j+1)=rmiss
          else
             geo(i,jm-j+1)=dble(dgeo(i,j))
          end if
       end do
    end do

  end subroutine read_geo

  !--------------------------

  subroutine read_anl_surf(iunit,iyr,imon,iday,ihour,slp,qa,t,u,v)
    
    use mod_rmiss
    implicit none

    real(kind = 4),parameter :: dmiss=9.999e20
    
    integer status,system,access
    integer i,j

    real(kind = 4) dslp(im,jm),dqa(im,jm),dt(im,jm),du(im,jm),dv(im,jm)

    character(100) filename
    character(15) ymdhn
    character(10) ymdh
    character(5) nnnnn
    character(4) yyyy
    character(2) mm,dd,hh

    !IN
    integer,intent(in) :: iunit
    integer,intent(in) :: iyr,imon,iday,ihour

    !OUT
    real(kind = 8),intent(out) :: slp(im,jm),qa(im,jm),t(im,jm),u(im,jm),v(im,jm)

    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    write(hh,'(i2.2)') ihour
    write(nnnnn,'(i5.5)') iunit

    ymdhn=yyyy//mm//dd//hh//nnnnn
    ymdh=yyyy//mm//dd//hh
    
    filename="/data/R/R2402/DATA/JRA55/anl_surf125/anl_surf125."//ymdh

    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not found "//trim(filename) 
       stop
    end if

    status=system("wgrib "//trim(filename)//" | grep :PRMSL: | wgrib "&
         &//trim(filename)//" -i -ieee -nh -o slp"//ymdhn//".bin")
    status=system("wgrib "//trim(filename)//" | grep :TMP: | wgrib "&
         &//trim(filename)//" -i -ieee -nh -o t"//ymdhn//".bin")
    status=system("wgrib "//trim(filename)//" | grep :SPFH: | wgrib "&
         &//trim(filename)//" -i -ieee -nh -o qa"//ymdhn//".bin")    
    status=system("wgrib "//trim(filename)//" | grep :UGRD: | wgrib "&
         &//trim(filename)//" -i -ieee -nh -o u"//ymdhn//".bin")    
    status=system("wgrib "//trim(filename)//" | grep :VGRD: | wgrib "&
         &//trim(filename)//" -i -ieee -nh -o v"//ymdhn//".bin")    

    !Sea level pressure [Pa]
    open(unit=iunit,file="slp"//ymdhn//".bin", &
         & access="direct",form="unformatted",recl=4*im*jm)
    read(iunit,rec=1) ((dslp(i,j),i=1,im),j=1,jm)
    close(iunit)
    status=system("rm -f slp"//ymdhn//".bin")

    !Temperature at 2m [K]
    open(unit=iunit,file="t"//ymdhn//".bin", &
         & access="direct",form="unformatted",recl=4*im*jm)
    read(iunit,rec=1) ((dt(i,j),i=1,im),j=1,jm)
    close(iunit)
    status=system("rm -f t"//ymdhn//".bin")

    !Specific Humidity at 2m [kg/kg]
    open(unit=iunit,file="qa"//ymdhn//".bin", &
         & access="direct",form="unformatted",recl=4*im*jm)
    read(iunit,rec=1) ((dqa(i,j),i=1,im),j=1,jm)
    close(iunit)
    status=system("rm -f qa"//ymdhn//".bin")

    !Zonal wind at 10m
    open(unit=iunit,file="u"//ymdhn//".bin", &
         & access="direct",form="unformatted",recl=4*im*jm)
    read(iunit,rec=1) ((du(i,j),i=1,im),j=1,jm)
    close(iunit)
    status=system("rm -f u"//ymdhn//".bin")

    !Meridional wind at 10m
    open(unit=iunit,file="v"//ymdhn//".bin", &
         & access="direct",form="unformatted",recl=4*im*jm)
    read(iunit,rec=1) ((dv(i,j),i=1,im),j=1,jm)
    close(iunit)
    status=system("rm -f v"//ymdhn//".bin")

    !Modify unit & Calculate Wind speed (s)
    do j=1,jm
       do i=1,im

          if(dslp(i,j) == dmiss)then
             slp(i,jm-j+1)=rmiss
          else
             slp(i,jm-j+1)=dble(dslp(i,j))
          end if

          if(dt(i,j) == dmiss)then
             t(i,jm-j+1)=rmiss
          else
             t(i,jm-j+1)=dble(dt(i,j)-273.15d0) ![K] -> [degree C]
          endif

          if(dqa(i,j) == dmiss)then
             qa(i,jm-j+1)=rmiss
          else
             qa(i,jm-j+1)=dble(dqa(i,j))
          end if

          if(du(i,j) == dmiss)then
             u(i,jm-j+1)=rmiss
          else
             u(i,jm-j+1)=dble(du(i,j))
          end if

          if(dv(i,j) == dmiss)then
             v(i,jm-j+1)=rmiss
          else
             v(i,jm-j+1)=dble(dv(i,j))
          end if

       end do
    end do
    
  end subroutine read_anl_surf

  !-------------------------------

  subroutine read_anl_surf_wind(iunit,iyr,imon,iday,ihour,u,v)
    
    use mod_rmiss
    implicit none

    real(kind = 4),parameter :: dmiss=9.999e20
    
    integer status,system,access
    integer i,j
    
    real(kind = 4) du(im,jm),dv(im,jm)

    character(100) filename
    character(15) ymdhn
    character(10) ymdh
    character(5) nnnnn
    character(4) yyyy
    character(2) mm,dd,hh

    !IN
    integer,intent(in) :: iunit
    integer,intent(in) :: iyr,imon,iday,ihour

    !OUT
    real(kind = 8),intent(out) :: u(im,jm),v(im,jm)
    
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    write(hh,'(i2.2)') ihour
    write(nnnnn,'(i5.5)') iunit

    ymdhn=yyyy//mm//dd//hh//nnnnn
    ymdh=yyyy//mm//dd//hh
    
    filename="/data/R/R2402/DATA/JRA55/anl_surf125/anl_surf125."//ymdh

    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not found "//trim(filename) 
       stop
    end if

    status=system("wgrib "//trim(filename)//" | grep :UGRD: | wgrib "&
         &//trim(filename)//" -i -ieee -nh -o u"//ymdhn//".bin")    
    status=system("wgrib "//trim(filename)//" | grep :VGRD: | wgrib "&
         &//trim(filename)//" -i -ieee -nh -o v"//ymdhn//".bin")    

    !Zonal wind at 10m
    open(unit=iunit,file="u"//ymdhn//".bin", &
         & access="direct",form="unformatted",recl=4*im*jm)
    read(iunit,rec=1) ((du(i,j),i=1,im),j=1,jm)
    close(iunit)
    status=system("rm -f u"//ymdhn//".bin")

    !Meridional wind at 10m
    open(unit=iunit,file="v"//ymdhn//".bin", &
         & access="direct",form="unformatted",recl=4*im*jm)
    read(iunit,rec=1) ((dv(i,j),i=1,im),j=1,jm)
    close(iunit)
    status=system("rm -f v"//ymdhn//".bin")

    !Modify unit & Calculate Wind speed (s)
    do j=1,jm
       do i=1,im

          if(du(i,j) == dmiss)then
             u(i,jm-j+1)=rmiss
          else
             u(i,jm-j+1)=dble(du(i,j))
          end if

          if(dv(i,j) == dmiss)then
             v(i,jm-j+1)=rmiss
          else
             v(i,jm-j+1)=dble(dv(i,j))
          end if

       end do
    end do
    
  end subroutine read_anl_surf_wind

  !-------------------------------

  subroutine read_anl_surf_var(iunit,var,iyr,imon,iday,ihour,dat)
    
    use mod_rmiss
    implicit none

    real(kind = 4),parameter :: dmiss=9.999e20
    
    integer status,system,access
    integer i,j
    
    real(kind = 4) ddat(im,jm)
    real(kind = 8) add

    character(100) filename
    character(15) ymdhn
    character(10) ymdh
    character(5) nnnnn
    character(4) yyyy
    character(2) mm,dd,hh
    character(*) var

    integer,intent(in) :: iunit
    integer,intent(in) :: iyr,imon,iday,ihour

    real(kind = 8),intent(out) :: dat(im,jm)
    
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    write(hh,'(i2.2)') ihour
    write(nnnnn,'(i5.5)') iunit

    ymdhn=yyyy//mm//dd//hh//nnnnn
    ymdh=yyyy//mm//dd//hh
    
    filename="/data/R/R2402/DATA/JRA55/anl_surf125/anl_surf125."//ymdh

    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not found "//trim(filename) 
       stop
    end if

    if(trim(var) == "slp")then
       status=system("wgrib "//trim(filename)//" | grep :PRMSL: | wgrib "&
            &//trim(filename)//" -i -ieee -nh -o atm"//ymdhn//".bin")
       add=0.d0
    else if(trim(var) == "ta")then
       status=system("wgrib "//trim(filename)//" | grep :TMP: | wgrib "&
            &//trim(filename)//" -i -ieee -nh -o atm"//ymdhn//".bin")
       add=-273.15d0
    else if(trim(var) == "qa")then
       status=system("wgrib "//trim(filename)//" | grep :SPFH: | wgrib "&
            &//trim(filename)//" -i -ieee -nh -o atm"//ymdhn//".bin")
       add=0.d0
    else if(trim(var) == "u")then
       status=system("wgrib "//trim(filename)//" | grep :UGRD: | wgrib "&
            &//trim(filename)//" -i -ieee -nh -o atm"//ymdhn//".bin")
       add=0.d0
    else if(trim(var) == "v")then
       status=system("wgrib "//trim(filename)//" | grep :VGRD: | wgrib "&
            &//trim(filename)//" -i -ieee -nh -o atm"//ymdhn//".bin")    
       add=0.d0
    else
       write(*,*) "***Error Not Found variable:"//trim(var)
       stop
    end if

    !Read data
    open(unit=iunit,file="atm"//ymdhn//".bin",&
         &access="direct",form="unformatted",recl=4*im*jm)
    read(iunit,rec=1) ((ddat(i,j),i=1,im),j=1,jm)
    close(iunit)
    status=system("rm -f atm"//ymdhn//".bin")

    !Modify unit & Calculate Wind speed (s)
    do j=1,jm
       do i=1,im

          if(ddat(i,j) == dmiss)then
             dat(i,jm-j+1)=rmiss
          else
             dat(i,jm-j+1)=dble(ddat(i,j))+add
          end if

       end do
    end do
    
  end subroutine read_anl_surf_var

  !-------------------------------

  subroutine read_fcst_surf(iunit,iyr,imon,iday,ihour,cloud)

    use mod_rmiss
    implicit none
    
    real(kind = 4),parameter :: dmiss=9.999e20

    integer status,system,access
    integer i,j

    real(kind = 4) dcloud(im,jm)

    character(100) filename
    character(15) ymdhn
    character(10) ymdh
    character(5) nnnnn
    character(4) yyyy
    character(2) mm,dd,hh

    integer,intent(in) :: iunit
    integer,intent(in) :: iyr,imon,iday,ihour

    real(kind = 8),intent(out) :: cloud(im,jm)
    
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    write(hh,'(i2.2)') ihour
    write(nnnnn,'(i5.5)') iunit

    ymdhn=yyyy//mm//dd//hh//nnnnn
    ymdh=yyyy//mm//dd//hh
    
    filename="/data/R/R2402/DATA/JRA55/fcst_surf125/fcst_surf125."//ymdh

    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not found "//trim(filename) 
       stop
    end if

    status=system("wgrib "//trim(filename)//" | grep :TCDC: | wgrib "&
         &//trim(filename)//" -i -ieee -nh -o cloud"//ymdhn//".bin")

    !Total Cloud [%]
    open(unit=iunit,file="cloud"//ymdhn//".bin", &
         & access="direct",form="unformatted",recl=4*im*jm)
    read(iunit,rec=1) ((dcloud(i,j),i=1,im),j=1,jm)
    close(iunit)
    status=system("rm -f cloud"//ymdhn//".bin")

    !Modify unit
    do j=1,jm
       do i=1,im

          if(dcloud(i,j) == dmiss)then
             cloud(i,jm-j+1)=rmiss
          else
             cloud(i,jm-j+1)=0.01d0*dcloud(i,j) ![%] -> [0-1]
          end if

       end do
    end do
    
  end subroutine read_fcst_surf
  
  !--------------------------------------

  subroutine read_fcst_phy2m(iunit,iyr,imon,iday,ihour,sw)

    use mod_rmiss
    implicit none
    
    real(kind = 4),parameter :: dmiss=9.999e20

    integer status,system,access
    integer i,j
    
    real(kind = 4) usw(im,jm),dsw(im,jm)

    character(100) filename
    character(15) ymdhn
    character(10) ymdh
    character(5) nnnnn
    character(4) yyyy
    character(2) mm,dd,hh

    integer,intent(in) :: iunit
    integer,intent(in) :: iyr,imon,iday,ihour

    real(kind = 8),intent(out) :: sw(im,jm)
    
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    write(hh,'(i2.2)') ihour
    write(nnnnn,'(i5.5)') iunit

    ymdhn=yyyy//mm//dd//hh//nnnnn
    ymdh=yyyy//mm//dd//hh
    
    filename="/data/R/R2402/DATA/JRA55/fcst_phy2m125/fcst_phy2m125."//ymdh

    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not found "//trim(filename) 
       stop
    end if

    status=system("wgrib "//trim(filename)//" | grep :USWRF: | grep :sfc: | wgrib "&
         &//trim(filename)//" -i -ieee -nh -o usw"//ymdhn//".bin")
    status=system("wgrib "//trim(filename)//" | grep :DSWRF: | grep :sfc: | wgrib "&
         &//trim(filename)//" -i -ieee -nh -o dsw"//ymdhn//".bin")

    !Upward Shortwave radiation [W/m^2]
    open(unit=iunit,file="usw"//ymdhn//".bin", &
         & access="direct",form="unformatted",recl=4*im*jm)
    read(iunit,rec=1) ((usw(i,j),i=1,im),j=1,jm)
    close(iunit)
    status=system("rm -f usw"//ymdhn//".bin")

    !Downward Shortwave radiation [W/m^2]
    open(unit=iunit,file="dsw"//ymdhn//".bin", &
         & access="direct",form="unformatted",recl=4*im*jm)
    read(iunit,rec=1) ((dsw(i,j),i=1,im),j=1,jm)
    close(iunit)
    status=system("rm -f dsw"//ymdhn//".bin")

    !Calculate Net shortwave radiation
    do j=1,jm
       do i=1,im

          if(usw(i,j) == dmiss .or. dsw(i,j) == dmiss)then
             sw(i,jm-j+1)=rmiss
          else
             sw(i,jm-j+1)=-1.d0*dble(usw(i,j))+dble(dsw(i,j))
          end if
          
       end do
    end do
    
  end subroutine read_fcst_phy2m

  !------------------------------------------

  subroutine read_fcst_phy2m_prep(iunit,iyr,imon,iday,ihour,prep)

    use mod_rmiss
    implicit none
    
    real(kind = 4),parameter :: dmiss=9.999e20

    integer status,system,access
    integer i,j

    real(kind = 4) dprep(im,jm)

    character(100) filename
    character(15) ymdhn
    character(10) ymdh
    character(5) nnnnn
    character(4) yyyy
    character(2) mm,dd,hh

    integer,intent(in) :: iunit
    integer,intent(in) :: iyr,imon,iday,ihour
    
    real(kind = 8),intent(out) :: prep(im,jm)
    
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    write(hh,'(i2.2)') ihour
    write(nnnnn,'(i5.5)') iunit

    ymdhn=yyyy//mm//dd//hh//nnnnn
    ymdh=yyyy//mm//dd//hh
    
    filename="/data/R/R2402/DATA/JRA55/fcst_phy2m125/fcst_phy2m125."//ymdh

    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: Not found "//trim(filename) 
       stop
    end if

    status=system("wgrib "//trim(filename)//" | grep :APCP: | grep :sfc: | wgrib "&
         &//trim(filename)//" -i -ieee -nh -o prep"//ymdhn//".bin")

    !Precipitation [mm/day]
    open(unit=iunit,file="prep"//ymdhn//".bin",&
         &access="direct",form="unformatted",recl=4*im*jm)
    read(iunit,rec=1) ((dprep(i,j),i=1,im),j=1,jm)
    close(iunit)
    status=system("rm -f prep"//ymdhn//".bin")

    !Calculate Precipitation
    do j=1,jm
       do i=1,im

          ![mm/day]
          if(dprep(i,j) == dmiss)then
             prep(i,jm-j+1)=rmiss
          else
             prep(i,jm-j+1)=dble(dprep(i,j))
          end if
          
       end do
    end do
    
  end subroutine read_fcst_phy2m_prep
  
end module mod_read_jra55
