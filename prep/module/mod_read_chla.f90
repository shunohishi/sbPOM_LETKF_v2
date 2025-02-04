module mod_read_chla

  integer,parameter :: im=7200,jm=3601
  
contains

  !-----------------------------------------------------------------
  ! GCOM-C/SGLI Monthly climatology |
  !----------------------------------
  !
  ! For mask, 0:Land+No obs., 1: Obs
  !
  !-----------------------------------------------------------------
  
  subroutine read_chla_clim(imon,lon,lat,mask,dat)

    use netcdf
    use mod_rmiss
    implicit none

    integer,parameter :: dmiss=-999
    
    !---Common
    integer i,j
    integer status,access
    integer varid,ncid

    real(kind = 4) dlon(im),dlat(jm)
    integer(kind = 2) idat(im,jm)
    
    
    character(2) mm
    character(200) filename
    
    !---IN
    integer,intent(in) :: imon

    !---OUT
    real(kind = 8),intent(out) :: lon(im),lat(jm)
    real(kind = 8),intent(out) :: mask(im,jm),dat(im,jm)

    write(mm,'(i2.2)') imon
    
    filename="/data/R/R2402/DATA/CHLA/clim/"//&
         &"GC1SG1_YYYY"//mm//"00D01M_D0000_3MSG_CHLAM_CLIM.nc" 

    status=access(trim(filename)," ")
    if(status /= 0)then
       write(*,*) "***Error: not found "//trim(filename)
       stop
    end if

    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    status=nf90_inq_varid(ncid,"longitude",varid)
    status=nf90_get_var(ncid,varid,dlon)

    status=nf90_inq_varid(ncid,"latitude",varid)
    status=nf90_get_var(ncid,varid,dlat)
    
    status=nf90_inq_varid(ncid,"CHLA_climatology",varid)
    status=nf90_get_var(ncid,varid,idat)
    
    lon(:)=dble(dlon(:))
    lat(:)=dble(dlat(:))
    
    do j=1,jm
       do i=1,im

          if(dlon(i) < 0.e0)then
             lon(i)=dlon(i)+360.d0
          else
             lon(i)=dlon(i)
          end if

          lat(jm-j+1)=dlat(j)
          
          if(idat(i,j) == dmiss)then
             mask(i,jm-j+1)=0.d0
             dat(i,jm-j+1)=rmiss
          else
             mask(i,jm-j+1)=1.d0
             dat(i,jm-j+1)=dble(idat(i,j))*0.003d0
          end if
          
       end do
    end do

  end subroutine read_chla_clim
  
end module mod_read_chla
