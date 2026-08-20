module mod_read_bran2020

  integer,parameter :: im=3600,jm=1500,km=51
  character(100),parameter :: bran2020_dir="/lvs0/rccs-dart/ohishi/DATA/BRAN2020"
  
contains
  
  !---------------------------------------------------------------------------
  ! Read BRAN2020 |
  !-----------------
  !
  ! Web: https://thredds.nci.org.au/thredds/catalog/gb6/BRAN/BRAN2020/catalog.html
  ! DOI: https://dx.doi.org/10.25914/6009627c7af03
  ! 
  ! - BRAN2020 from CSIRO (Australia): DA => EnOI  (3 days)
  !
  !---------------------------------------------------------------------------
  ! Select
  !
  ! varname:
  ! - t,s,u,v,h
  ! 
  !---------------------------------------------------------------------------

  subroutine get_bran2020_info(varname,prefixname,lonname,latname,depname,ncname,add,mult)

    implicit none

    !---Common
    
    !---IN
    character(1),intent(in) :: varname
    
    !---OUT
    real(kind = 8),intent(out) :: add  !add_offset
    real(kind = 8),intent(out) :: mult !scale_factor
    character(20),intent(out) :: prefixname
    character(20),intent(out) :: lonname,latname,depname,ncname
    
    !---Variable name
    if(varname == "t")then
       prefixname="ocean_temp"
       lonname="xt_ocean"
       latname="yt_ocean"
       depname="st_ocean"
       ncname="temp"
       add=245.d0
       mult=0.00778222d0      
    else if(varname == "s")then
       prefixname="ocean_salt"
       lonname="xt_ocean"
       latname="yt_ocean"
       depname="st_ocean"
       ncname="salt"
       add=45.d0
       mult=0.001678518d0      
    else if(varname == "u")then
       prefixname="ocean_u"
       lonname="xu_ocean"
       latname="yu_ocean"
       depname="st_ocean" 
       ncname="u"
       add=0.d0
       mult=0.0003051851d0      
    else if(varname == "v")then
       prefixname="ocean_v"
       lonname="xu_ocean"
       latname="yu_ocean"
       depname="st_ocean"
       ncname="v"
       add=0.d0
       mult=0.0003051851d0     
    else if(varname == "h")then
       prefixname="ocean_eta_t"
       lonname="xt_ocean"
       latname="yt_ocean"
       depname=""
       ncname="eta_t"
       add=0.d0
       mult=1.d0
    else
       write(*,*) "***Error: Incorrect varname => "//trim(varname)
       stop
    end if
    
  end subroutine get_bran2020_info

  !----------------------------------------------------------------------------
  
  subroutine read_bran2020(varname,iyr,imon,iday,km_in,lon,lat,depth,mask,dat)

    use mod_rmiss
    use netcdf
    implicit none

    !---Parameter
    integer,parameter :: imiss=-32768
    real(kind = 4),parameter :: dmiss=-1.e20
    
    !---Common
    integer i,j,k
    integer status,access
    integer ncid,varid    

    integer itmp3d(im,jm,km_in)
    real(kind = 4) tmp3d(im,jm,km_in)
    
    real(kind = 4) tmp1dx(im),tmp1dy(jm),tmp1dz(km_in)
    real(kind = 8) add,mult
    
    character(200) filename
    character(20) prefixname
    character(20) lonname,latname,depname,ncname
    character(4) yyyy
    character(2) mm
    
    !---IN
    integer,intent(in) :: iyr,imon,iday
    integer,intent(in) :: km_in

    character(1),intent(in)  :: varname 

    !---OUT
    real(kind = 8),intent(out) :: lon(im),lat(jm),depth(km_in)
    real(kind = 8),intent(out) :: mask(im,jm),dat(im,jm,km_in)

    !---Get ncname
    call get_bran2020_info(varname,prefixname,lonname,latname,depname,ncname,add,mult)
    
    !---Filename
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon

    filename=trim(bran2020_dir)//"/"//trim(prefixname)//"_"//yyyy//"_"//mm//".nc"
    
    status=access(trim(filename)," ")
    if(status == 0)then
       write(*,*) "Read :"//trim(filename)
    else
       write(*,*) "***Error: Not found "//trim(filename)
       stop
    end if
    
    !---Read data
    status=nf90_open(trim(filename),nf90_nowrite,ncid)

    status=nf90_inq_varid(ncid,trim(lonname),varid)
    status=nf90_get_var(ncid,varid,tmp1dx)

    status=nf90_inq_varid(ncid,trim(latname),varid)
    status=nf90_get_var(ncid,varid,tmp1dy)

    if(varname == "h")then
       tmp1dz(:)=0.e0
    else
       status=nf90_inq_varid(ncid,trim(depname),varid)
       status=nf90_get_var(ncid,varid,tmp1dz,(/1/),(/km_in/))
    end if
       
    if(varname == "h")then
       status=nf90_inq_varid(ncid,trim(ncname),varid)
       status=nf90_get_var(ncid,varid,tmp3d,(/1,1,iday/),(/im,jm,1/))
    else
       status=nf90_inq_varid(ncid,trim(ncname),varid)
       status=nf90_get_var(ncid,varid,itmp3d,(/1,1,1,iday/),(/im,jm,km_in,1/))
    end if

    status=nf90_close(ncid)
    
    !---Post process
    !Longitude
    lon(:)=dble(tmp1dx(:))
    
    !Latitude
    lat(:)=dble(tmp1dy(:))

    !Depth
    depth(:)=dble(tmp1dz(:))

    !Mask
    k=1
    if(varname == "h")then
       do j=1,jm
          do i=1,im
             if(tmp3d(i,j,k) == dmiss)then
                mask(i,j)=0.d0
             else
                mask(i,j)=1.d0
             end if
          end do
       end do
    else
       do j=1,jm
          do i=1,im
             if(itmp3d(i,j,k) == imiss)then
                mask(i,j)=0.d0
             else
                mask(i,j)=1.d0
             end if
          end do
       end do       
    end if
        
    !Data
    if(varname == "h")then
       k=1
       do j=1,jm
          do i=1,im
             if(tmp3d(i,j,k) == dmiss)then
                dat(i,j,k)=rmiss
             else
                dat(i,j,k)=dble(tmp3d(i,j,k))*mult+add
             end if
          end do
       end do
    else
       do k=1,km_in
          do j=1,jm
             do i=1,im
                if(itmp3d(i,j,k) == imiss)then
                   dat(i,j,k)=rmiss
                else
                   dat(i,j,k)=dble(itmp3d(i,j,k))*mult+add
                end if
             end do
          end do
       end do
    end if
           
    !Missing value
    do j=1,jm
       do i=1,im
          if(mask(i,j) == 0.d0)then
             dat(i,j,:)=rmiss
          end if
       end do
    end do
        
  end subroutine read_bran2020
  
end module mod_read_bran2020
