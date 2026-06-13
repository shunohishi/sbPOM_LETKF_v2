program main

  use setting
  use mod_varname
  use mod_julian
  use mod_gridinfo
  use mod_read_lora
  implicit none

  !---Common
  integer ijul,sjul,ejul
  integer iyr,imon,iday
  integer iyr_n1,imon_n1,iday_n1
  integer ivar,imem,ims
  integer k

  !---Grid
  real(kind = 8) lont(im),lonu(im),lonv(im)
  real(kind = 8) latt(jm),latu(jm),latv(jm)
  real(kind = 8) dept(im,jm,km),depu(im,jm,km),depv(im,jm,km),depw(im,jm,km)

  !---Land-Sea Mask
  real(kind = 8) mask(im,jm),maskt(im,jm),masku(im,jm),maskv(im,jm)

  !---Data
  integer id_mld(im,jm)
  
  real(kind = 8) dat2d(im,jm),dat3d(im,jm,km)
  real(kind = 8) ens2d(im,jm,nmem)

  real(kind = 8) t(im,jm,km),s(im,jm,km),rho(im,jm,km)
  real(kind = 8) mlt(im,jm),mls(im,jm),mld(im,jm)

  !---MLT/MLS budget (See Kim et al. 2006)
  integer id_mld_n1(im,jm) !previous time (t-1)

  real(kind = 8) mld_n1(im,jm)
  real(kind = 8) mld_ens(im,jm) !Max(mld, mld_n1)
  real(kind = 8) dhdt(im,jm)    !MLD difference btw t-1 and t
  
  real(kind = 8) t_n1(im,jm,km),s_n1(im,jm,km),rho_n1(im,jm,km)
  real(kind = 8) t1_n1(im,jm) !Depth-mean T that stays within ML btw t-1 and t
  real(kind = 8) t2_n1(im,jm) !Depth-mean T to be detrained/entrained from t-1 to t
  real(kind = 8) s1_n1(im,jm) !Depth-mean T that stays within ML btw t-1 and t
  real(kind = 8) s2_n1(im,jm) !Depth-mean T to be detrained/entrained from t-1 to t

  character(4) ms
  character(100) varname,long_name,units_name

  !---Julian day
  call ymd_julian(syr,smon,sday,sjul)
  call ymd_julian(eyr,emon,eday,ejul)  

  !---Read and Write grid data
  write(*,*) "Read grid"
  call read_grid(dir,&
       & lont,lonu,lonv, &
       & latt,latu,latv, &
       & dept,depu,depv,depw, &
       & maskt,masku,maskv)

  write(*,*) "Make grid NetCDF file"
  call make_grid_ncfile(im,jm,km)

  write(*,*) "Write grid NetCDF file"
  !call write_grid_ncfile(im,jm,km, &
  !     & lont,lonu,lonv, &
  !     & latt,latu,latv, &
  !     & dept,depu,depv,depw, &
  !     & maskt,masku,maskv)

  !---Main loop
  do ijul=sjul,ejul

     call julian_ymd(ijul,iyr,imon,iday)
     call julian_ymd(ijul-1,iyr_n1,imon_n1,iday_n1)
     write(*,*) "Date:",iyr,imon,iday

     !---Basic variables
     do ims=1,2

        if(ims == 1)then
           ms="mean"
        else if(ims == 2)then
           ms="sprd"
        end if

        !---2D
        do ivar=1,nvar2d

           k=1
           imem=0
           mask(:,:)=maskt(:,:)

           call varname2d(ivar,varname,long_name,units_name)
           write(*,*) "2D: "//trim(varname)
           
           call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,k,mask,dat2d)
           call make_dat2d_ncfile(ms,varname,long_name,units_name,iyr,imon,iday,im,jm)
           call write_dat2d_ncfile(ms,varname,iyr,imon,iday,im,jm,dat2d)

        end do !ivar(2d)

        !---3D
        do ivar=1,nvar3d

           imem=0
           call varname3d(ivar,varname,long_name,units_name)
           write(*,*) "3D: "//trim(varname)
           
           if(varname == "u")then
              mask(:,:)=masku(:,:)
           else if(varname == "v")then
              mask(:,:)=maskv(:,:)
           else
              mask(:,:)=maskt(:,:)
           end if

           call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,km,mask,dat3d)
           call make_dat3d_ncfile(ms,varname,long_name,units_name,iyr,imon,iday,im,jm,km)
           call write_dat3d_ncfile(ms,varname,iyr,imon,iday,im,jm,km,dat3d)

           !Save ensemble mean temperature and salinity
           if(ims == 1 .and. varname == "t")then
              t(:,:,:)=dat3d(:,:,:)
           else if(ims == 1 .and. varname == "s")then
              s(:,:,:)=dat3d(:,:,:)
           end if
        end do !ivar(3d)

     end do !ims

     !---T and S at t-1
     ms="mean"
     call read_anal(dir,letkf,region,ms,imem,"t",iyr_n1,imon_n1,iday_n1,im,jm,km,maskt,t_n1)
     call read_anal(dir,letkf,region,ms,imem,"s",iyr_n1,imon_n1,iday_n1,im,jm,km,maskt,s_n1)
     
     !---MLD
     call potential_density_3d(im,jm,km,maskt,t,s,rho)
     call potential_density_3d(im,jm,km,maskt,t_n1,s_n1,rho_n1)
     
     call estimate_mld(im,jm,km,mask,dept,rho,id_mld,mld)
     call estimate_mld(im,jm,km,mask,dept,rho_n1,id_mld_n1,mld_n1)

     !---Detrained/Entrained water
     ***

     !---MLT
     call mld_ave(im,jm,km,id_mld,maskt,dept,t,mlt)
     do ivar=1,nvar_mlt

        imem=0
        mask(:,:)=maskt(:,:)

        call varname_mlt(ivar,varname,long_name,units_name)

        if(ivar == 1)then !Surface flux
           k=1
           dat3d(:,:,:)=0.d0
           call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,k,mask,dat3d(:,:,k))
           call mld_ave(im,jm,km,id_mld,maskt,dept,dat3d,dat2d)
        else if(ivar == 11)then !Detrainment/Entrainment
           
        else if(ivar == 13)then !Resitual
           
        else
           call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,km,mask,dat3d)
           call mld_ave(im,jm,km,id_mld,maskt,dept,dat3d,dat2d)           
        end if

        call make_dat2d_ncfile(ms,varname,long_name,units_name,iyr,imon,iday,im,jm)
        call write_dat2d_ncfile(ms,varname,iyr,imon,iday,im,jm,dat2d)

     end do !ivar

     !---MLS
     call mld_ave(im,jm,km,id_mld,maskt,dept,s,mls)
     do ivar=1,nvar_mls

        imem=0
        mask(:,:)=maskt(:,:)

        call varname_mls(ivar,varname,long_name,units_name)

        !***

     end do !ivar

     !---ML
     do ivar=1,nvar_mld

     end do !ivar     
     
     !---Ensemble sea surface
     do ivar=1,nvar_ens

        imem=0

        call varname_ens(ivar,varname,long_name,units_name)

        if(varname == "u")then
           mask(:,:)=masku(:,:)
        else if(varname == "v")then
           mask(:,:)=maskv(:,:)
        else
           mask(:,:)=maskt(:,:)
        end if

        do imem=1,nmem
           !call read_
        end do

        !call write_

     end do !ivar
        
  end do !ijul

end program main
