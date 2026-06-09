!-------------------------------------------------------------------------------------------------
! Post process for web visualization |
!-------------------------------------------------------------------------------------------------
! Created by Shun OHISHI 2026/06/09
!-------------------------------------------------------------------------------------------------
! Update:
! Name yyyy/mm/dd Detail comment
!
!-------------------------------------------------------------------------------------------------

program main

  use setting
  use mod_gridinfo, im_in => im, jm_in => jm, km_in => km
  use mod_read_lora
  use mod_julian
  implicit none

  integer ijul0,ijul,sjul,ejul
  integer iyr,imon,iday
  integer ivar,ims
  integer i,j,k
  integer itmp
  
  !---INPUT
  real(kind = 8) lon_in(im_in),lat_in(jm_in)
  real(kind = 8) lont_in(im_in),latt_in(jm_in)
  real(kind = 8) lonu_in(im_in),latu_in(jm_in)
  real(kind = 8) lonv_in(im_in),latv_in(jm_in)

  real(kind = 8) dep_in(im_in,jm_in,km_in)

  real(kind = 8) mask_in(im_in,jm_in)
  real(kind = 8) fsm_in(im_in,jm_in),usm_in(im_in,jm_in),vsm_in(im_in,jm_in)

  real(kind = 8) dat2d_in(im_in,jm_in)
  real(kind = 8) dat3d_in(im_in,jm_in,km_in)

  real(kind = 8) null2d_in(im_in,jm_in),null3d_in(im_in,jm_in,km_in)
  
  character(4) ms
  character(10) varname
  
  !---INT
  integer idx(im_out),idy(jm_out)
  integer idxt(im_out),idyt(jm_out)
  integer idxu(im_out),idyu(jm_out)
  integer idxv(im_out),idyv(jm_out)
  integer idz(im_in,jm_in,km_out) !IN grid number

  real(kind = 8) dat3d_int(im_in,jm_in,km_out)
  
  !---OUTPUT
  real(kind = 8) rjul
  real(kind = 8) lon_out(im_out),lat_out(jm_out)
  real(kind = 8) lont_out(im_out),latt_out(jm_out)
  real(kind = 8) lonu_out(im_out),latu_out(jm_out)
  real(kind = 8) lonv_out(im_out),latv_out(jm_out)

  real(kind = 8) dep_out(km_out)

  real(kind = 8) dat2d_out(im_out,jm_out)
  real(kind = 8) dat3d_out(im_out,jm_out,km_out)
  
  write(*,'(a)') "=== Start: Create Netcdf file for Web visualization ==="
  
  !---Julian day
  call ymd_julian(1950,1,1,ijul0)
  call ymd_julian(syr,smon,sday,sjul)
  call ymd_julian(eyr,emon,eday,ejul)

  !---Read grid data
  call read_grid(dir, &
       & lont_in,lonu_in,lonv_in, &
       & latt_in,latu_in,latv_in, &
       & dep_in,null3d_in,null3d_in, &
       & fsm_in,usm_in,vsm_in)

  !---Output grid
  lont_out(:)=lont_in(:)
  lonu_out(:)=lonu_in(:)
  lonv_out(:)=lonv_in(:)
  latt_out(:)=latt_in(:)
  latu_out(:)=latu_in(:)
  latv_out(:)=latv_in(:)
  call make_depth(dep_out)

  !---ID
  call cal_idxy(im_in,jm_in,lont_in,latt_in, &
       & im_out,jm_out,lont_out,latt_out,idxt,idyt)
  call cal_idxy(im_in,jm_in,lonu_in,latu_in, &
       & im_out,jm_out,lonu_out,latu_out,idxu,idyu)
  call cal_idxy(im_in,jm_in,lonv_in,latv_in, &
       & im_out,jm_out,lonv_out,latv_out,idxv,idyv)

  call cal_idz(im_in,jm_in,km_in,dep_in,km_out,dep_out,idz)

  do ijul=sjul,ejul

     call julian_ymd(ijul,iyr,imon,iday)
     rjul=REAL(ijul-ijul0)+0.5e0 !Daily mean
     write(*,*) "Date:",iyr,imon,iday

     !---Write grid data
     call write_grid(iyr,imon,im_out,jm_out,km_out, &
          & lont_out,lonu_out,lonv_out,latt_out,latu_out,latv_out,dep_out)

     !---2D
     do ivar=1,nvar2d

        if(ivar == 1)then
           varname="el"
           lon_in(:)=lont_in(:)
           lat_in(:)=latt_in(:)
           mask_in(:,:)=fsm_in(:,:)
           lon_out(:)=lont_out(:)
           lat_out(:)=latt_out(:)
        else
           write(*,*) "***Error: varname in 2d"
           stop
        end if
        
        do ims=1,2

           if(ims == 1)then
              ms="mean"
           else if(ims == 2)then
              ms="sprd"
           endif

           !---Read data
           itmp=0
           k=1
           call read_anal(dir,letkf,region,ms,itmp,varname,iyr,imon,iday,im_in,jm_in,k,mask_in,dat2d_in)

           !---Bilinear interpolation (skip)
           
           !---Write data
           k=1
           call write_data(trim(varname)//trim(ms),iyr,imon,iday,im_out,jm_out,k,rjul,dat2d_in)
           
        end do !i
     end do !ivar

     !---3D
     do ivar=1,nvar3d

        if(ivar == 1)then
           varname="t"
           lon_in(:)=lont_in(:)
           lat_in(:)=latt_in(:)
           mask_in(:,:)=fsm_in(:,:)
           lon_out(:)=lont_out(:)
           lat_out(:)=latt_out(:)
        else if(ivar == 2)then
           varname="s"
           lon_in(:)=lont_in(:)
           lat_in(:)=latt_in(:)
           mask_in(:,:)=fsm_in(:,:)
           lon_out(:)=lont_out(:)
           lat_out(:)=latt_out(:)
        else if(ivar == 3)then
           varname="u"
           lon_in(:)=lonu_in(:)
           lat_in(:)=latu_in(:)
           mask_in(:,:)=usm_in(:,:)
           lon_out(:)=lonu_out(:)
           lat_out(:)=latu_out(:)
        else if(ivar == 4)then
           varname="v"
           lon_in(:)=lonv_in(:)
           lat_in(:)=latv_in(:)
           mask_in(:,:)=vsm_in(:,:)
           lon_out(:)=lonv_out(:)
           lat_out(:)=latv_out(:)
        else
           write(*,*) "***Error: varname in 3d"
           stop           
        end if
           
        
        do ims=1,2

           if(ims == 1)then
              ms="mean"
           else if(ims == 2)then
              ms="sprd"
           endif

           !---Read data
           itmp=0
           call read_anal(dir,letkf,region,ms,itmp,varname,iyr,imon,iday,im_in,jm_in,km_in,mask_in,dat3d_in)

           !---Bilinear interpolation (skip)
           
           !---Vertical interpolation
           call vertical_linear_interpolation(im_in,jm_in,km_in,dep_in,mask_in,dat3d_in, &
                & idz,km_out,dep_out,dat3d_out)

           !---Write data
           call write_data(trim(varname)//trim(ms),iyr,imon,iday,im_out,jm_out,km_out,rjul,dat3d_out)
           
        end do !i

     end do !ivar
             
  end do !ijul

  write(*,'(a)') "=== End ==="
  
end program main
