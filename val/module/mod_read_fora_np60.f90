module mod_read_fora_np60

  integer,parameter :: im=2049,jm=784,km=60
  character(100),parameter :: fora_np60_dir="/lvs0/rccs-dart/ohishi/DATA/FORA-JPN60/NP"
  
contains
  
  !---------------------------------------------------------------------------
  ! Read FORA-JPN/WNP60 |
  !----------------------
  !
  ! Web: https://www.jamstec.go.jp/fora/j/
  ! DOI: https://doi.org/10.1007/s10872-026-00795-x
  ! 
  ! - FORA-JPN60 from MRI (Japan): DA => 4D-VAR  (10 days)
  ! - Arakawa B grid (t,s,ssh / u, v)
  !---------------------------------------------------------------------------
  ! Select
  !
  ! varname:
  ! - t,s,u,v,h
  ! 
  !---------------------------------------------------------------------------

  subroutine get_fora_np60_info(varname,prefixname,lonname,latname,depname,ncname,add,mult)

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
       prefixname="Basic-3D/t/nc_t."
       lonname="lon"
       latname="lat"
       depname="depth"
       ncname="thetao"
       add=0.d0
       mult=1.d0      
    else if(varname == "s")then
       prefixname="Basic-3D/s/nc_s."
       lonname="lon"
       latname="lat"
       depname="depth"
       ncname="so"
       add=0.d0
       mult=1.d0     
    else if(varname == "u")then
       prefixname="Basic-3D/u/nc_u."
       lonname="lon"
       latname="lat"
       depname="depth"
       ncname="uo"
       add=0.d0
       mult=1.d-2 ![cm/s] => [m/s]     
    else if(varname == "v")then
       prefixname="Basic-3D/v/nc_v."
       lonname="lon"
       latname="lat"
       depname="depth"
       ncname="vo"
       add=0.d0
       mult=1.d-2 ![cm/s] => [m/s]     
    else if(varname == "h")then
       prefixname="Basic-2D/ssh/nc_ssh."
       lonname="lon"
       latname="lat"
       depname="depth"
       ncname="zos"
       add=0.d0
       mult=1.d-2 ![cm] => [m] 
    else
       write(*,*) "***Error: Incorrect varname => "//trim(varname)
       stop
    end if
    
  end subroutine get_fora_np60_info

  !----------------------------------------------------------------------------
  
  subroutine read_fora_np60(varname,iyr,imon,iday,km_in,lon,lat,depth,mask,dat)

    !$use omp_lib    
    use mod_rmiss
    use netcdf
    implicit none

    !---Parameter
    real(kind = 4),parameter :: dmiss=-9.99e33
    
    !---Common
    integer i,j,k
    integer status,access
    integer ncid,varid    

    real(kind = 4) tmp3d(im,jm,km_in)
    
    real(kind = 8) add,mult
    
    character(200) filename
    character(20) prefixname
    character(20) lonname,latname,depname,ncname
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

    !---Get ncname
    call get_fora_np60_info(varname,prefixname,lonname,latname,depname,ncname,add,mult)
    
    !---Filename
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    yyyymmdd=yyyy//mm//dd
    
    filename=trim(fora_np60_dir)//"/"//trim(prefixname)//yyyymmdd
    
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
    status=nf90_get_var(ncid,varid,lon)

    status=nf90_inq_varid(ncid,trim(latname),varid)
    status=nf90_get_var(ncid,varid,lat)

    if(varname == "h")then
       depth(:)=0.e0
    else
       status=nf90_inq_varid(ncid,trim(depname),varid)
       status=nf90_get_var(ncid,varid,depth,(/1/),(/km_in/))
    end if
       
    if(varname == "h")then
       status=nf90_inq_varid(ncid,trim(ncname),varid)
       status=nf90_get_var(ncid,varid,tmp3d(:,:,1),(/1,1,1/),(/im,jm,1/))
    else
       status=nf90_inq_varid(ncid,trim(ncname),varid)
       status=nf90_get_var(ncid,varid,tmp3d(:,:,:),(/1,1,1,1/),(/im,jm,km_in,1/))
    end if

    status=nf90_close(ncid)
    
    !---Post process
    !Mask
    !$omp parallel
    !$omp do private(i,j) collapse(2)
    do j=1,jm
       do i=1,im
          if(tmp3d(i,j,1) == dmiss)then
             mask(i,j)=0.d0
          else
             mask(i,j)=1.d0
          end if
       end do
    end do
    !$omp end do       
    
    !Data
    !$omp do private(i,j,k) collapse(3)
    do k=1,km_in
       do j=1,jm
          do i=1,im
             if(tmp3d(i,j,k) == dmiss)then
                dat(i,j,k)=rmiss
             else
                dat(i,j,k)=dble(tmp3d(i,j,k))*mult+add
             end if
          end do
       end do
    end do
    !$omp end do       
    
    !Missing value
    !$omp do private(i,j) collapse(2)
    do j=1,jm
       do i=1,im
          if(mask(i,j) == 0.d0)then
             dat(i,j,:)=rmiss
          end if
       end do
    end do
    !$omp end do       
    !$omp end parallel
        
  end subroutine read_fora_np60

  !----------------------------------------------------------------------------
  
  subroutine extract_fora_np60(varname,iyr,imon,iday,is,im_in,js,jm_in,ks,km_in,lon,lat,depth,mask,dat)

    !$use omp_lib    
    use mod_rmiss
    use netcdf
    implicit none

    !---Parameter
    real(kind = 4),parameter :: dmiss=-9.99e33
    
    !---Common
    integer i,j,k
    integer status,access
    integer ncid,varid    

    real(kind = 4) tmp3d(im_in,jm_in,km_in)
    
    real(kind = 8) add,mult
    
    character(200) filename
    character(20) prefixname
    character(20) lonname,latname,depname,ncname
    character(8) yyyymmdd
    character(4) yyyy
    character(2) mm,dd
    
    !---IN
    integer,intent(in) :: iyr,imon,iday    
    integer,intent(in) :: is,im_in
    integer,intent(in) :: js,jm_in
    integer,intent(in) :: ks,km_in

    character(1),intent(in)  :: varname 

    !---OUT
    real(kind = 8),intent(out) :: lon(im_in),lat(jm_in),depth(km_in)
    real(kind = 8),intent(out) :: mask(im_in,jm_in),dat(im_in,jm_in,km_in)

    !---Get ncname
    call get_fora_np60_info(varname,prefixname,lonname,latname,depname,ncname,add,mult)
    
    !---Filename
    write(yyyy,'(i4.4)') iyr
    write(mm,'(i2.2)') imon
    write(dd,'(i2.2)') iday
    yyyymmdd=yyyy//mm//dd
    
    filename=trim(fora_np60_dir)//"/"//trim(prefixname)//yyyymmdd
    
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
    status=nf90_get_var(ncid,varid,lon,(/is/),(/im_in/))

    status=nf90_inq_varid(ncid,trim(latname),varid)
    status=nf90_get_var(ncid,varid,lat,(/js/),(/jm_in/))

    if(varname == "h")then
       depth(:)=0.e0
    else
       status=nf90_inq_varid(ncid,trim(depname),varid)
       status=nf90_get_var(ncid,varid,depth,(/ks/),(/km_in/))
    end if
       
    if(varname == "h")then
       status=nf90_inq_varid(ncid,trim(ncname),varid)
       status=nf90_get_var(ncid,varid,tmp3d(:,:,1),(/is,js,1/),(/im_in,jm_in,1/))
    else
       status=nf90_inq_varid(ncid,trim(ncname),varid)
       status=nf90_get_var(ncid,varid,tmp3d(:,:,:),(/is,js,ks,1/),(/im_in,jm_in,km_in,1/))
    end if

    status=nf90_close(ncid)
    
    !---Post process
    !Mask
    !$omp parallel
    !$omp do private(i,j) collapse(2)
    do j=1,jm_in
       do i=1,im_in
          if(tmp3d(i,j,1) == dmiss)then
             mask(i,j)=0.d0
          else
             mask(i,j)=1.d0
          end if
       end do
    end do
    !$omp end do       
    
    !Data
    !$omp do private(i,j,k) collapse(3)
    do k=1,km_in
       do j=1,jm_in
          do i=1,im_in
             if(tmp3d(i,j,k) == dmiss)then
                dat(i,j,k)=rmiss
             else
                dat(i,j,k)=dble(tmp3d(i,j,k))*mult+add
             end if
          end do
       end do
    end do
    !$omp end do       
    
    !Missing value
    !$omp do private(i,j) collapse(2)
    do j=1,jm_in
       do i=1,im_in
          if(mask(i,j) == 0.d0)then
             dat(i,j,:)=rmiss
          end if
       end do
    end do
    !$omp end do       
    !$omp end parallel
        
  end subroutine extract_fora_np60
  
end module mod_read_fora_np60
