module mod_read_glorys12v1

  integer,parameter :: im=4320,jm=2041,km=50
  character(100),parameter :: g12v1_dir="/lvs0/rccs-dart/ohishi/DATA/GLORYS/010"
  
contains
  
  !---------------------------------------------------------------------------
  ! Read GLORYS12V1 |
  !-----------------
  !
  ! Web: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description
  ! DOI: https://doi.org/10.48670/moi-00021
  ! 
  ! - GLORYS12Vv from Mercator Ocean (Fr): DA => SEEK  (7 days)
  !
  !---------------------------------------------------------------------------
  ! Select
  !
  ! varname:
  ! - t,s,u,v,h
  ! 
  !---------------------------------------------------------------------------

  subroutine get_glorys12v1_info(varname,ncname,add,mult)

    implicit none

    !---Common
    
    !---IN
    character(1),intent(in) :: varname
    
    !---OUT
    real(kind = 8),intent(out) :: add  !add_offset
    real(kind = 8),intent(out) :: mult !scale_factor
    character(20),intent(out) :: ncname


    
    !---Variable name
    if(varname == "t")then
       ncname="thetao"
       add=21.d0
       mult=0.000732444226741791d0      
    else if(varname == "s")then
       ncname="so"
       add=-0.00152592547237873d0
       mult=0.00152592547237873d0     
    else if(varname == "u")then
       ncname="uo"
       add=0.d0
       mult=0.000610370188951492d0  
    else if(varname == "v")then
       ncname="vo"
       add=0.d0
       mult=0.000610370188951492d0  
    else if(varname == "h")then
       ncname="zos"
       add=0.d0
       mult=0.000305185094475746d0
    else
       write(*,*) "***Error: Incorrect varname => "//trim(varname)
       stop
    end if
    
  end subroutine get_glorys12v1_info

  !----------------------------------------------------------------------------
  
  subroutine read_glorys12v1(varname,iyr,imon,iday,km_in,lon,lat,depth,mask,dat)

    use mod_rmiss
    use netcdf
    implicit none

    !---Parameter
    integer,parameter :: dmiss=-32767
    
    !---Common
    integer i,j,k,n
    integer status,system,access
    integer ncid,varid    

    integer itmp3d(im,jm,km_in)
    
    real(kind = 4) tmp1dx(im),tmp1dy(jm),tmp1dz(km_in)
    real(kind = 8) add,mult
    
    character(200) filename
    character(20) ncname
    character(8) yyyymmdd
    character(4) yyyy
    character(2) mm,dd
    
    !---IN
    integer,intent(in) :: iyr,imon,iday
    integer,intent(in) :: km_in

    character(1),intent(in)  :: varname 

    !---OUT
    real(kind = 8),intent(out) :: lon(im),lat(jm),depth(km_in)
    real(kind = 8),intent(out) :: mask(im,jm),dat(im,jm,km_in)

    !---Filename
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    yyyymmdd=yyyy//mm//dd

    status=system("find "//trim(g12v1_dir)//" -type f -name "//"mercatorglorys12v1_gl12_mean_"//yyyymmdd//"_R*.nc > glorys12v1.txt")

    open(1,file="glorys12v1.txt",status="old")
    read(1,'(a)') filename
    close(1)

    status=system("rm -f glorys12v1.txt")
    
    status=access(trim(filename)," ")
    if(status == 0)then
       write(*,*) "Read :"//trim(filename)
    else
       write(*,*) "***Error: Not found "//trim(filename)
       stop
    end if
    
    !---Get ncname
    call get_glorys12v1_info(varname,ncname,add,mult)

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
       status=nf90_get_var(ncid,varid,itmp3d,(/1,1,1/),(/im,jm,1/))
    else
       status=nf90_inq_varid(ncid,trim(ncname),varid)
       status=nf90_get_var(ncid,varid,itmp3d,(/1,1,1,1/),(/im,jm,km_in,1/))
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
             if(itmp3d(i,j,k) == dmiss)then
                mask(n,j)=0.d0
             else
                mask(n,j)=1.d0
             end if
          end if
       end do

       do i=1,im
          if(tmp1dx(i) < 0.e0)then
             n=n+1
             if(itmp3d(i,j,k) == dmiss)then
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
                if(itmp3d(i,j,k) == dmiss)then
                   dat(n,j,k)=rmiss
                else
                   dat(n,j,k)=dble(itmp3d(i,j,k))*mult+add
                end if
             end if
          end do

          do i=1,im
             if(tmp1dx(i) < 0.e0)then
                n=n+1
                if(itmp3d(i,j,k) == dmiss)then
                   dat(n,j,k)=rmiss
                else
                   dat(n,j,k)=dble(itmp3d(i,j,k))*mult+add
                end if
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
        
  end subroutine read_glorys12v1
  
end module mod_read_glorys12v1
