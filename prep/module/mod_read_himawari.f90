module mod_read_himawari

  integer,parameter :: im=6001,jm=6001

contains

  !----------------------------------------------------------
  ! Read Himawari 8 |
  !----------------------------------------------------------
  !
  ! Quality level = 5: OK, otherwise: rmiss
  !
  !----------------------------------------------------------

  subroutine read_himawari(iyr,imon,iday,ihour,long,lati,dat)

    use netcdf
    use mod_rmiss
    implicit none

    !Common
    real(kind = 4),parameter :: dmiss=-32768.e0,dql=5.e0

    integer i,j
    integer status,status1,status2,status3,status4,status_all,access
    integer ncid,varid

    real(kind = 4) rlon(im),rlat(jm)
    real(kind = 4) rdat(im,jm),rql(im,jm)

    character(200) filename,filename1,filename2,filename3,filename4
    character(10) yyyymmddhh
    character(4) yyyy
    character(2) mm,dd,hh
    
    !IN
    integer,intent(in) :: iyr,imon,iday,ihour
    
    !OUT
    real(kind = 8),intent(out) :: long(im),lati(jm)
    real(kind = 8),intent(out) :: dat(im,jm)
    
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    write(hh,'(i2.2)') ihour
    yyyymmddhh=yyyy//mm//dd//hh
    
    filename1="/data/R/R2402/DATA/Himawari8/"//yyyy//"/"&
         &//yyyymmddhh//"0000-JAXA-L3C_GHRSST-SSTskin-H08_AHI-v1.2-v02.0-fv01.0.nc"
    filename2="/data/R/R2402/DATA/Himawari8/"//yyyy//"/"&
         &//yyyymmddhh//"0000-JAXA-L3C_GHRSST-SSTskin-H08_AHI-v1.2-v02.0-fv02.0.nc"
    filename3="/data/R/R2402/DATA/Himawari8/"//yyyy//"/"&
         &//yyyymmddhh//"0000-JAXA-L3C_GHRSST-SSTskin-H08_AHI-v2.0-v02.0-fv01.0.nc"
    filename4="/data/R/R2402/DATA/Himawari9/"//yyyy//"/"&
         &//yyyymmddhh//"0000-JAXA-L3C_GHRSST-SSTskin-H09_AHI_NRT-v2.1-v02.0-fv01.0.nc"

    status1=access(trim(filename1)," ")
    status2=access(trim(filename2)," ")
    status3=access(trim(filename3)," ")
    status4=access(trim(filename4)," ")

    if(status4 == 0)then
       status_all=0
       filename=filename4
    else if(status3 == 0)then
       status_all=0
       filename=filename3
    else if(status2 == 0)then
       status_all=0
       filename=filename2
    elseif(status1 == 0)then
       status_all=0
       filename=filename1
    else
       status_all=999
       filename="/data/R/R2402/DATA/Himawari8/"//&
            &"2015/20150707000000-JAXA-L3C_GHRSST-SSTskin-H08_AHI-v1.2-v02.0-fv01.0.nc"
    end if

    if(status_all == 0)then
       write(*,'(a)') "Read: "//trim(filename)
    else
       write(*,'(a)') "Read: "//trim(filename)//" as dummy"
       write(*,'(a)') "Check: Not Found Himawari-8/-9 data at"//yyyymmddhh//"if the date is available preiod"
    end if

    status=nf90_open(trim(filename),nf90_nowrite,ncid)
    
    status=nf90_inq_varid(ncid,"lon",varid)
    status=nf90_get_var(ncid,varid,rlon)

    status=nf90_inq_varid(ncid,"lat",varid)
    status=nf90_get_var(ncid,varid,rlat)

    status=nf90_inq_varid(ncid,"sea_surface_temperature",varid)
    status=nf90_get_var(ncid,varid,rdat,(/1,1,1/),(/im,jm,1/))

    status=nf90_inq_varid(ncid,"quality_level",varid)
    status=nf90_get_var(ncid,varid,rql,(/1,1,1/),(/im,jm,1/))

    status=nf90_close(ncid)

    !Modify longitude
    long(:)=dble(rlon(:))
    do i=1,im
       if(long(i) < 0.d0) long(i)=long(i)+360.d0
    end do

    do j=1,jm
       lati(j)=dble(rlat(jm-j+1))
    end do

    !Modify SST
    do j=1,jm
       do i=1,im
          
          if(rdat(i,j) == dmiss .or. rql(i,j) /= dql)then
             rdat(i,j)=real(rmiss)
          else
             rdat(i,j)=rdat(i,j)*0.01e0
          end if

       end do !i
    end do !j

    if(status_all == 0)then
       do j=1,jm
          do i=1,im
             dat(i,j)=dble(rdat(i,jm-j+1))
          end do
       end do
    else
       dat(:,:)=rmiss
    end if
    
  end subroutine read_himawari
  
end module mod_read_himawari
