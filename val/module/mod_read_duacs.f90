module mod_read_duacs

  integer,parameter :: im=1440,jm=720

  character(100),parameter :: dir_duacs="/data/R/R2402/DATA/AVISO"
  
contains

  !--------------------------------------------------------------------
  ! Read DUACS |
  !--------------------------------------------------------------------

  subroutine read_duacs_var(varname,iyr,imon,iday,lon,lat,dat)

    use mod_rmiss
    use netcdf
    implicit none

    !---Parameter
    integer,parameter :: dmiss=-2147483647
    
    !---Common
    integer i,j

    integer status,system,access
    integer ncid,varid

    integer idat(im,jm)
    integer ipass
    
    real(kind = 4) rlon(im),rlat(jm)
    
    character(200) filename
    character(8) yyyymmdd
    character(4) yyyy
    character(2) mm,dd    
    
    !---IN
    integer,intent(in) :: iyr,imon,iday

    character(*),intent(in) :: varname

    !---OUT    
    real(kind = 8),intent(out) :: lon(im),lat(jm)
    real(kind = 8),intent(out) :: dat(im,jm)

    !---Date
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    
    yyyymmdd=yyyy//mm//dd

    !---Filename
    status=system("find "//trim(dir_duacs)//"/"//yyyy//" -type f -name dt_global_allsat_phy_l4_"//yyyymmdd//"_*.nc > filename.txt")

    open(1,file="filename.txt",status="old")
    read(1,'(a)') filename
    close(1)
    
    status=system("rm -f filename.txt")

    !---Access
    status=access(trim(filename)," ")

    if(status == 0)then
       write(*,*) "Read "//trim(filename)
    else
       write(*,*) "***Error: Not found "//trim(filename)
       stop
    end if

    !---NetCDF
    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    status=nf90_inq_varid(ncid,"longitude",varid)
    status=nf90_get_var(ncid,varid,rlon)

    status=nf90_inq_varid(ncid,"latitude",varid)
    status=nf90_get_var(ncid,varid,rlat)
    
    status=nf90_inq_varid(ncid,trim(varname),varid)
    status=nf90_get_var(ncid,varid,idat)
    
    status=nf90_close(ncid)

    ipass=0
    do i=1,im
       if(0.e0 <= rlon(i) .and. rlon(i) <= 180.d0)then
          ipass=ipass+1
          lon(ipass)=rlon(i)
       end if       
    end do

    do i=1,im
       if(-180.e0 <= rlon(i) .and. rlon(i) < 0.d0)then
          ipass=ipass+1
          lon(ipass)=rlon(i)+360.d0
       end if
    end do

    lat(:)=dble(rlat(:))
    
    do j=1,jm

       ipass=0
       do i=1,im

          if(0.e0 <= rlon(i) .and. rlon(i) <= 180.d0)then
             ipass=ipass+1

             if(idat(i,j) == dmiss)then
                dat(ipass,j)=rmiss
             else
                dat(ipass,j)=dble(idat(i,j))*0.0001d0
             end if

          end if
       end do

       do i=1,im
          if(-180.e0 <= rlon(i) .and. rlon(i) < 0.d0)then
             ipass=ipass+1

             if(idat(i,j) == dmiss)then
                dat(ipass,j)=rmiss
             else
                dat(ipass,j)=dble(idat(i,j))*0.0001d0
             end if
             
          end if
       end do

    end do

  end subroutine read_duacs_var

  !---------------------------------------------------------------------
  
  subroutine read_duacs_mdot(lon,lat,mdot)

    use mod_rmiss
    implicit none

    !---Parameter
    integer,parameter :: iyr=2003,imon=1,iday=1

    !---Common
    integer i,j
    
    real(kind = 8) adt(im,jm),sla(im,jm)

    !---OUT
    real(kind = 8),intent(out) :: lon(im),lat(jm)
    real(kind = 8),intent(out) :: mdot(im,jm)


    call read_duacs_var("adt",iyr,imon,iday,lon,lat,adt)
    call read_duacs_var("sla",iyr,imon,iday,lon,lat,sla)

    do j=1,jm
       do i=1,im
          if(adt(i,j) == rmiss .or. sla(i,j) == rmiss)then
             mdot(i,j)=rmiss
          else
             mdot(i,j)=adt(i,j)-sla(i,j)
          end if
       end do
    end do
    
  end subroutine read_duacs_mdot

end module mod_read_duacs
