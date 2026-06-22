program main

  use setting
  use mod_varname
  use mod_julian
  use mod_density
  use mod_gridinfo
  use mod_read_lora
  implicit none

  !---Common
  integer ijul,sjul,ejul
  integer iyr,imon,iday
  integer iyr_t0,imon_t0,iday_t0
  integer ivar,imem,ims
  integer k

  !---Grid
  real(kind = 8) lont(im),lonu(im),lonv(im) !Lontigude at t, u, v points
  real(kind = 8) latt(jm),latu(jm),latv(jm) !Latitude at t, u, v points
  real(kind = 8) dept(im,jm,km),depu(im,jm,km),depv(im,jm,km),depw(im,jm,km) !Depth at t, u, v, w points
  real(kind = 8) dzt(im,jm,km) !Layer thickness at T point

  !---Land-Sea Mask
  real(kind = 8) mask(im,jm) !General
  real(kind = 8) maskt(im,jm),masku(im,jm),maskv(im,jm) !t, u, v points

  !---Data
  integer id_mld(im,jm)
  
  real(kind = 8) dat2d(im,jm)      !2D
  real(kind = 8) dat3d(im,jm,km)   !3D
  real(kind = 8) ens2d(im,jm,nmem) !2D ensemble at sea surface

  real(kind = 8) t(im,jm,km),s(im,jm,km),rho(im,jm,km)
  real(kind = 8) mld(im,jm),mlt(im,jm),mls(im,jm)
  
  character(4) ms                             !mean/sprd/eens
  character(100) varname,long_name,units_name !Variable name
  
  !---MLT and MLS budget (Kim et al. 2006)
  integer id_mld_t0(im,jm),id_mld_t1(im,jm)  !MLD ID at t=t0 (e.g., 00UTC) and t=t1 (e.g., 24UTC)

  real(kind = 8) mld_t0(im,jm),mld_t1(im,jm) !MLD at t=t0 and t1
  real(kind = 8) mld_ent(im,jm)              !MLD in detrainment/entrainment estimation: Max(mld_t0, mld_t1)
  real(kind = 8) mlt_t0(im,jm),mlt_t1(im,jm) !MLT at t=t0 and t1
  real(kind = 8) mls_t0(im,jm),mls_t1(im,jm) !MLS at t=t0 and t1
  
  real(kind = 8) dhdt(im,jm)                   !MLD difference between t0 and t1 *Unit: [m]  
  real(kind = 8) delta_t(im,jm),delta_s(im,jm) !Delta T = T1-T2 and Delta S = S1-S2
  real(kind = 8) tent(im,jm),sent(im,jm)       !Temperature and Salinity detrainment/entrainment
  real(kind = 8) tres(im,jm),sres(im,jm)       !Residual
  
  !---Instanteneous value
  real(kind = 8) t_t0(im,jm,km),t_t1(im,jm,km)     !T at t0 and t1
  real(kind = 8) s_t0(im,jm,km),s_t1(im,jm,km)     !S at t0 and t1
  real(kind = 8) rho_t0(im,jm,km),rho_t1(im,jm,km) !Density at t0 and t1

  write(*,*) "=== Start: Creat FTP file ==="
  
  !---Julian day
  call ymd_julian(syr,smon,sday,sjul)
  call ymd_julian(eyr,emon,eday,ejul)  

  if(sday /= 1)then
     write(*,*) "***Error: Set sday = 1 to read ensemble mean restart file at 1st day"
     stop
  end if
  
  !---Read and Write grid data
  write(*,*) "Read grid"
  call read_grid(dir,&
       & lont,lonu,lonv, &
       & latt,latu,latv, &
       & dept,depu,depv,depw, &
       & maskt,masku,maskv)
  call make_dzt(im,jm,km,maskt,depw,dzt)

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
     call julian_ymd(ijul-1,iyr_t0,imon_t0,iday_t0)
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

           k=1    !Dummy
           imem=0 !Dummy
           mask(:,:)=maskt(:,:)

           call varname2d(ivar,varname,long_name,units_name)
           write(*,*) "2D: "//trim(varname)//" "//ms
           
           call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,k,mask,dat2d)
           call make_dat2d_ncfile(ms,varname,long_name,units_name,iyr,imon,iday,im,jm)
           call write_dat2d_ncfile(ms,varname,iyr,imon,iday,im,jm,dat2d)

        end do !ivar(2d)

        !---3D
        do ivar=1,nvar3d

           imem=0 !Dummy
           call varname3d(ivar,varname,long_name,units_name)
           write(*,*) "3D: "//trim(varname)//" "//ms
           
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

           !Save ensemble mean temperature and salinity to calculate MLD
           if(ims == 1 .and. varname == "t")then
              t(:,:,:)=dat3d(:,:,:)
           else if(ims == 1 .and. varname == "s")then
              s(:,:,:)=dat3d(:,:,:)
           end if
           
        end do !ivar(3d)

     end do !ims

     !---1-day averaged T and S
     call potential_density_3d(im,jm,km,maskt,t,s,rho)
     call estimate_mld(im,jm,km,maskt,depw,rho,id_mld,mld)
     
     !---T and S at t=t0 (instanteneous)
     if(iday == 1)then

        call read_restart(dir_restart,letkf,"t",iyr_t0,imon_t0,iday_t0,im,jm,km,t_t0)
        call read_restart(dir_restart,letkf,"s",iyr_t0,imon_t0,iday_t0,im,jm,km,s_t0)

        call potential_density_3d(im,jm,km,maskt,t_t0,s_t0,rho_t0)
        call estimate_mld(im,jm,km,maskt,depw,rho_t0,id_mld_t0,mld_t0)

        call mld_ave(im,jm,km,id_mld_t0,maskt,dzt,t_t0,mlt_t0)
        call mld_ave(im,jm,km,id_mld_t0,maskt,dzt,s_t0,mls_t0)

     end if        
     
     !---T and S at t=t1 (instanteneous)
     ms="mean"
     imem=0 !Dummy
     varname="dtdt"
     call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,km,maskt,dat3d)
     call add(im,jm,km,maskt,t_t0,dat3d,t_t1)
     
     varname="dsdt"
     call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,km,maskt,dat3d)
     call add(im,jm,km,maskt,s_t0,dat3d,s_t1)

     call potential_density_3d(im,jm,km,maskt,t_t1,s_t1,rho_t1)
     call estimate_mld(im,jm,km,maskt,depw,rho_t1,id_mld_t1,mld_t1)

     call mld_ave(im,jm,km,id_mld_t1,maskt,dzt,t_t1,mlt_t1)
     call mld_ave(im,jm,km,id_mld_t1,maskt,dzt,s_t1,mls_t1)
     
     !---Detrainment and Entrainment
     call detrainment_entrainment(im,jm,km,maskt,dzt, &
          & id_mld_t0,mld_t0,t_t0,s_t0, &
          & id_mld_t1,mld_t1, &
          & mld_ent,dhdt,delta_t,delta_s,tent,sent)
     
     !---MLT
     do ivar=1,nvar_mlt

        imem=0 !Dummy
        mask(:,:)=maskt(:,:)
        ms="mean"

        !---Varname
        call varname_mlt(ivar,varname,long_name,units_name)
        write(*,*) "MLT: "//trim(varname)

        !---Read data
        if(ivar == 1)then !dTmix/dt = MLT(t1)-MLT(t0)
           call dmltdt(im,jm,maskt,mlt_t0,mlt_t1,dat2d)
           tres(:,:)=dat2d(:,:)
        else if(ivar == 2)then !Surface flux
           k=1
           call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,k,mask,dat2d)
           dat3d(:,:,:)=0.d0
           dat3d(:,:,1)=dat2d(:,:)
           call mld_ave(im,jm,km,id_mld_t1,maskt,dzt,dat3d,dat2d)
        else if(ivar == 10)then !Detrainment/Entrainment
           dat2d(:,:)=tent(:,:)
        else if(ivar == nvar_mlt)then !Residual
           dat2d(:,:)=tres(:,:)
        else
           call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,km,mask,dat3d)
           call mld_ave(im,jm,km,id_mld_t1,maskt,dzt,dat3d,dat2d)           
        end if

        !---Residual
        if(ivar /= 1 .and. ivar /= nvar_mlt)then
           tres(:,:)=tres(:,:)-dat2d(:,:)
        end if

        !---Write NetCDF data
        call make_dat2d_ncfile(ms,varname,long_name,units_name,iyr,imon,iday,im,jm)
        call write_dat2d_ncfile(ms,varname,iyr,imon,iday,im,jm,dat2d)
        
     end do !ivar

     !---MLS
     do ivar=1,nvar_mls

        imem=0 !Dummy
        mask(:,:)=maskt(:,:)
        ms="mean"

        !---Varname
        call varname_mls(ivar,varname,long_name,units_name)
        write(*,*) "MLS: "//trim(varname)

        !---Read data
        if(ivar == 1)then !dSmix/dt = MLS(t1)-MLS(t0)
           call dmltdt(im,jm,maskt,mls_t0,mls_t1,dat2d)
           sres(:,:)=dat2d(:,:)
        else if(ivar == 2)then !Surface flux
           k=1
           call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,k,mask,dat2d)
           dat3d(:,:,:)=0.d0
           dat3d(:,:,1)=dat2d(:,:)
           call mld_ave(im,jm,km,id_mld_t1,maskt,dzt,dat3d,dat2d)
        else if(ivar == 9)then !Detrainment/Entrainment
           dat2d(:,:)=sent(:,:)
        else if(ivar == nvar_mls)then !Residual
           dat2d(:,:)=sres(:,:)
        else
           call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,km,mask,dat3d)
           call mld_ave(im,jm,km,id_mld_t1,maskt,dzt,dat3d,dat2d)           
        end if

        !---Residual
        if(ivar /= 1 .and. ivar /= nvar_mls)then
           sres(:,:)=sres(:,:)-dat2d(:,:)
        end if

        !---Write NetCDF data
        call make_dat2d_ncfile(ms,varname,long_name,units_name,iyr,imon,iday,im,jm)
        call write_dat2d_ncfile(ms,varname,iyr,imon,iday,im,jm,dat2d)
        
     end do !ivar
     
     !---ML
     do ivar=1,nvar_mld

        call varname_mld(ivar,varname,long_name,units_name)
        write(*,*) "MLD: "//trim(varname)
        
        if(ivar == 1)then
           dat2d(:,:)=mld(:,:)
        else if(ivar == 2)then
           dat2d(:,:)=mld_ent(:,:)
        else if(ivar == 3)then
           dat2d(:,:)=dhdt(:,:)
        else if(ivar == 4)then
           dat2d(:,:)=delta_t(:,:)
        else if(ivar == 5)then
           dat2d(:,:)=delta_s(:,:)
        end if

        call make_dat2d_ncfile(ms,varname,long_name,units_name,iyr,imon,iday,im,jm)
        call write_dat2d_ncfile(ms,varname,iyr,imon,iday,im,jm,dat2d)
           
     end do !ivar     
     
     !---Ensemble sea surface
     do ivar=1,nvar_ens

        imem=0
        k=1
        ms="eens"

        call varname_ens(ivar,varname,long_name,units_name)
        write(*,*) "Ensemble: "//trim(varname)

        if(varname == "u")then
           mask(:,:)=masku(:,:)
        else if(varname == "v")then
           mask(:,:)=maskv(:,:)
        else
           mask(:,:)=maskt(:,:)
        end if

        do imem=1,nmem
           call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,k,mask,ens2d(:,:,imem))
        end do

        call make_ens2d_ncfile(ms,varname,long_name,units_name,iyr,imon,iday,im,jm,nmem)
        call write_ens2d_ncfile(ms,varname,iyr,imon,iday,im,jm,nmem,ens2d)

     end do !ivar

     !---Update
     t_t0(:,:,:)=t_t1(:,:,:)
     s_t0(:,:,:)=s_t1(:,:,:)
     rho_t0(:,:,:)=rho_t1(:,:,:)

     id_mld_t0(:,:)=id_mld_t1(:,:)
     mld_t0(:,:)=mld_t1(:,:)
     mlt_t0(:,:)=mlt_t1(:,:)
     mls_t0(:,:)=mls_t1(:,:)
     
  end do !ijul

  write(*,*) "=== End: Creat FTP file ==="
  
end program main

!--------------------------------------------------------------------------

subroutine add(im,jm,km,mask,dat0,dat1,dat)

  !$ use omp_lib
  use mod_rmiss
  implicit none

  !---Common
  integer i,j,k

  !---IN
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: mask(im,jm)
  real(kind = 8),intent(in) :: dat0(im,jm,km),dat1(im,jm,km)

  !---OUT
  real(kind = 8),intent(out) :: dat(im,jm,km)

  !$omp parallel do private(i,j,k)
  do k=1,km
     do j=1,jm
        do i=1,im
           if(mask(i,j) == 0.d0)then
              dat(i,j,k)=rmiss
           else
              dat(i,j,k)=dat0(i,j,k)+dat1(i,j,k)
           end if
        end do
     end do
  end do
  !$omp end parallel do
                         
end subroutine add
