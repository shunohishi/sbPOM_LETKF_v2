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
  integer iyr,imon,iday,jday,nday
  integer iyr_t0,imon_t0,iday_t0
  integer ivar,imem,ims
  integer i,j,k

  !---Grid
  real(kind = 8) lont(im),lonu(im),lonv(im) !Lontigude at t, u, v points
  real(kind = 8) latt(jm),latu(jm),latv(jm) !Latitude at t, u, v points
  real(kind = 8) dept(im,jm,km),depu(im,jm,km),depv(im,jm,km),depw(im,jm,km) !Depth at t, u, v, w points *Plus

  real(kind = 8) dx(jm,jm),dy(im,jm) !Longitude and latitude grid interval at T point [m]
  real(kind = 8) dzt(im,jm,km)       !Layer thickness at T point [m]
  real(kind = 8) sigmat(im,jm,km)    !Sigma at T point [-]

  !---Land-Sea Mask
  real(kind = 8) mask(im,jm)                            !General
  real(kind = 8) maskt(im,jm),masku(im,jm),maskv(im,jm) !t, u, v points

  !---Data
  integer id_mld(im,jm)
  integer,allocatable :: pass2d(:,:),pass3d(:,:,:)

  real(kind = 8),allocatable :: ave2d(:,:),ave3d(:,:,:)  
  real(kind = 8),allocatable :: dat2d(:,:),dat3d(:,:,:) !2D, 3D(im,jm,km) or (im,jm,nmem)

  real(kind = 8) ssh(im,jm)
  real(kind = 8) t(im,jm,km),s(im,jm,km),rho(im,jm,km)
  real(kind = 8) u(im,jm,km),v(im,jm,km)
  
  character(4) ms                             !mean/sprd/eens
  character(100) varname,long_name,units_name !Variable name

  !---Instanteneous value
  real(kind = 8) t_t0(im,jm,km),t_t1(im,jm,km)     !T at t0 and t1
  real(kind = 8) s_t0(im,jm,km),s_t1(im,jm,km)     !S at t0 and t1
  real(kind = 8) rho_t0(im,jm,km),rho_t1(im,jm,km) !Density at t0 and t1
  
  !---MLT and MLS budget (Kim et al. 2006)
  integer id_mld_t0(im,jm),id_mld_t1(im,jm)  !MLD ID at t=t0 (e.g., 00UTC) and t=t1 (e.g., 24UTC)

  real(kind = 8) mld_t0(im,jm),mld_t1(im,jm) !MLD at t=t0 and t1
  real(kind = 8) mld_ent(im,jm)              !MLD in detrainment/entrainment estimation: Max(mld_t0, mld_t1)
  real(kind = 8) mlt_t0(im,jm),mlt_t1(im,jm) !MLT at t=t0 and t1
  real(kind = 8) mls_t0(im,jm),mls_t1(im,jm) !MLS at t=t0 and t1

  integer,allocatable :: mltpass_out(:,:,:)
  integer,allocatable :: mlspass_out(:,:,:)
  integer,allocatable :: mldpass_out(:,:,:)
  
  real(kind = 8),allocatable :: mlt_out(:,:,:),mltave_out(:,:,:) !MLT-related variables
  real(kind = 8),allocatable :: mls_out(:,:,:),mlsave_out(:,:,:) !MLS
  real(kind = 8),allocatable :: mld_out(:,:,:),mldave_out(:,:,:) !MLD

  real(kind = 8),allocatable :: xadv_cor(:,:,:) !Zonal advection correction term (Sigma to Z coordinate)
  real(kind = 8),allocatable :: yadv_cor(:,:,:) !Meridional
  real(kind = 8),allocatable :: zadv_cor(:,:,:) !Vertical
    
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

  if(llora20)then !***LORA version 2.0
     i=1
     j=1
     k=1
     call extract_grid(dir,"dx",im,jm,k,dx)
     call extract_grid(dir,"dy",im,jm,k,dy)
     call extract_grid(dir,"zz",km,i,j,sigmat(i,j,:))
     do j=1,jm
        do i=1,im
           sigmat(i,j,:)=sigmat(1,1,:)
        end do
     end do
  else !***LORA version 2.1
     call extract_grid(dir,"z_e",im,jm,km,sigmat)
  end if
     
  write(*,*) "Make grid NetCDF file"
  call make_grid_ncfile(im,jm,km)

  write(*,*) "Write grid NetCDF file"
  call write_grid_ncfile(im,jm,km, &
       & lont,lonu,lonv, &
       & latt,latu,latv, &
       & dept,depu,depv,depw, &
       & maskt,masku,maskv)

  !---Layer thickness at T point (For MLT ave)
  call make_dzt(im,jm,km,maskt,depw,dzt)
  
  !---2D basic variable------------------------------ 
  if(out_basic2d)then

     k=1    !Dummy
     imem=0 !Dummy
     allocate(dat2d(im,jm),ave2d(im,jm),pass2d(im,jm))

     do ivar=1,nvar2d

        !---Varname
        call varname2d(ivar,varname,long_name,units_name)

        !---Land-Sea Mask
        mask(:,:)=maskt(:,:)     

        do ims=1,2  

           !---Mean/Spread
           if(ims == 1)then
              ms="mean"
           else if(ims == 2)then
              ms="sprd"
           end if

           write(*,*) "2D: "//trim(varname)//" "//ms

           do ijul=sjul,ejul

              !---Date
              call julian_ymd(ijul,iyr,imon,iday)
              call estimate_nday(iyr,imon,nday)
              write(*,*) "Date:",iyr,imon,iday

              !---Read analysis ensemble mean/spread
              call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,k,mask,dat2d)

              !---Monthly average
              call month_ave(iday,nday,im,jm,k,mask,dat2d,ave2d,pass2d)

              !---Write daily analysis ensemble mean/spread
              call make_dat2d_ncfile(ms,varname,long_name,units_name,iyr,imon,iday,im,jm)
              call write_dat2d_ncfile(ms,varname,iyr,imon,iday,im,jm,dat2d)

              !---Write monthly-mean analysis ensemble mean/spread
              if(iday == nday)then
                 jday=0 !Dummy
                 call make_dat2d_ncfile(ms,varname,long_name,units_name,iyr,imon,jday,im,jm)
                 call write_dat2d_ncfile(ms,varname,iyr,imon,jday,im,jm,ave2d)
              end if

           end do !ijul
        end do !ims
     end do !ivar

     deallocate(dat2d,ave2d,pass2d)

  end if
     
  !---3D basic variable-----------------------------------
  if(out_basic3d)then
     
     imem=0 !Dummy
     allocate(dat3d(im,jm,km),ave3d(im,jm,km),pass3d(im,jm,km))

     do ivar=1,nvar3d

        !---Varname
        call varname3d(ivar,varname,long_name,units_name)

        !---Land-Sea Mask
        if(varname == "u")then
           mask(:,:)=masku(:,:)
        else if(varname == "v")then
           mask(:,:)=maskv(:,:)
        else
           mask(:,:)=maskt(:,:)
        end if

        do ims=1,2

           !---Mean/Spread
           if(ims == 1)then
              ms="mean"
           else if(ims == 2)then
              ms="sprd"
           end if

           write(*,*) "3D: "//trim(varname)//" "//ms

           do ijul=sjul,ejul

              !---Date
              call julian_ymd(ijul,iyr,imon,iday)
              call estimate_nday(iyr,imon,nday)
              write(*,*) "Date:",iyr,imon,iday

              !---Read analysis ensemble mean/spread
              call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,km,mask,dat3d)

              !---Monthly average
              call month_ave(iday,nday,im,jm,km,mask,dat3d,ave3d,pass3d)

              !---Write daily analysis ensemble mean/spread
              call make_dat3d_ncfile(ms,varname,long_name,units_name,iyr,imon,iday,im,jm,km)
              call write_dat3d_ncfile(ms,varname,iyr,imon,iday,im,jm,km,dat3d)

              !---Write monthly-mean analysis ensemble mean/spraed
              if(iday == nday)then
                 jday=0 !Dummy
                 call make_dat3d_ncfile(ms,varname,long_name,units_name,iyr,imon,jday,im,jm,km)
                 call write_dat3d_ncfile(ms,varname,iyr,imon,jday,im,jm,km,ave3d)
              end if

           end do !ijul
        end do !ims
     end do !ivar
     deallocate(dat3d,ave3d,pass3d)

  end if
     
  !---Ensemble at the sea surface--------------------------------
  if(out_ensemble2d)then

     allocate(dat3d(im,jm,nmem),ave3d(im,jm,nmem),pass3d(im,jm,nmem))
     k=1
     ms="eens"

     do ivar=1,nvar_ens

        !---Varname
        call varname_ens(ivar,varname,long_name,units_name)
        write(*,*) "Ensemble: "//trim(varname)

        !---Land-Sea mask
        if(varname == "u")then
           mask(:,:)=masku(:,:)
        else if(varname == "v")then
           mask(:,:)=maskv(:,:)
        else
           mask(:,:)=maskt(:,:)
        end if

        do ijul=sjul,ejul

           !---Date
           call julian_ymd(ijul,iyr,imon,iday)
           call estimate_nday(iyr,imon,nday)
           write(*,*) "Date:",iyr,imon,iday

           !---Read analysis ensemble
           do imem=1,nmem
              call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,k,mask,dat3d(:,:,imem))
           end do !imem

           !---Monthly average
           call month_ave(iday,nday,im,jm,nmem,mask,dat3d,ave3d,pass3d)

           !---Write analysis ensemble
           call make_dat3d_ncfile(ms,varname,long_name,units_name,iyr,imon,iday,im,jm,nmem)
           call write_dat3d_ncfile(ms,varname,iyr,imon,iday,im,jm,nmem,dat3d)

           !---Write monthly-mean analysis ensemble
           if(iday == nday)then
              jday=0 !Dummy
              call make_dat3d_ncfile(ms,varname,long_name,units_name,iyr,imon,jday,im,jm,nmem)
              call write_dat3d_ncfile(ms,varname,iyr,imon,jday,im,jm,nmem,ave3d)
           end if

        end do !ijul      
     end do !ivar

     deallocate(dat3d,ave3d,pass3d)

  end if
     
  !---MLT and MLS budget-------------------------------------------------------------------
  if(out_mlts2d)then
     
     imem=0 !Dummy
     k=1
     ms="mean"
     mask(:,:)=maskt(:,:)
  
     allocate(dat2d(im,jm),dat3d(im,jm,km))  
     allocate(mld_out(im,jm,nvar_mld),mldave_out(im,jm,nvar_mld),mldpass_out(im,jm,nvar_mld))
     allocate(mlt_out(im,jm,nvar_mlt),mltave_out(im,jm,nvar_mlt),mltpass_out(im,jm,nvar_mlt))
     allocate(mls_out(im,jm,nvar_mls),mlsave_out(im,jm,nvar_mls),mlspass_out(im,jm,nvar_mls))
     allocate(xadv_cor(im,jm,km),yadv_cor(im,jm,km),zadv_cor(im,jm,km))

     do ijul=sjul,ejul

        !---Date
        call julian_ymd(ijul,iyr,imon,iday)
        call julian_ymd(ijul-1,iyr_t0,imon_t0,iday_t0)
        call estimate_nday(iyr,imon,nday)
        write(*,*) "Date:",iyr,imon,iday

        !---Instanteneous variables at current time step (t=t0)        
        if(iday == 1)then

           !---Read Instanteneous T and S from a restart file           
           !T
           varname="t"
           call read_restart(dir_restart,letkf,varname,iyr_t0,imon_t0,iday_t0,im,jm,km,t_t0)

           !S
           varname="s"
           call read_restart(dir_restart,letkf,varname,iyr_t0,imon_t0,iday_t0,im,jm,km,s_t0)

           !Potential density
           call potential_density_3d(im,jm,km,mask,t_t0,s_t0,rho_t0)

           !MLD
           call estimate_mld(im,jm,km,mask,depw,rho_t0,id_mld_t0,mld_t0)

           !MLT
           call mld_ave(im,jm,km,id_mld_t0,mask,dzt,t_t0,mlt_t0)

           !MLS
           call mld_ave(im,jm,km,id_mld_t0,mask,dzt,s_t0,mls_t0)

        end if

        !---Instanteneous variables at next time step (t=t1)
        !T
        varname="dtdt"
        call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,km,mask,dat3d)
        call add(im,jm,km,mask,t_t0,dat3d,t_t1)

        !S
        varname="dsdt"
        call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,km,mask,dat3d)
        call add(im,jm,km,mask,s_t0,dat3d,s_t1)

        !Potential density
        call potential_density_3d(im,jm,km,mask,t_t1,s_t1,rho_t1)

        !MLD
        call estimate_mld(im,jm,km,mask,depw,rho_t1,id_mld_t1,mld_t1)

        !MLT
        call mld_ave(im,jm,km,id_mld_t1,mask,dzt,t_t1,mlt_t1)

        !MLS
        call mld_ave(im,jm,km,id_mld_t1,mask,dzt,s_t1,mls_t1)

        !---1-day averaged variables
        !T
        varname="t"
        call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,km,mask,t)

        !S
        varname="s"
        call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,km,mask,s)

        !U
        varname="u"
        call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,km,masku,u)

        !V
        varname="v"
        call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,km,maskv,v)

        !SSH
        varname="el"
        call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,k,mask,ssh)        
        
        !Potential density
        call potential_density_3d(im,jm,km,mask,t,s,rho)

        !MLD
        call estimate_mld(im,jm,km,mask,depw,rho,id_mld,mld_out(:,:,1))        

        !---Detrainment and Entrainment
        call detrainment_entrainment(im,jm,km,mask,dzt, &
             & id_mld_t0,mld_t0,t_t0,s_t0, &
             & id_mld_t1,mld_t1, &
             & mld_out(:,:,2),mld_out(:,:,3),mld_out(:,:,4),mld_out(:,:,5), &
             & mlt_out(:,:,10),mls_out(:,:,9))

        !---MLD-related variables 
        do ivar=1,nvar_mld

           !Varname
           call varname_mld(ivar,varname,long_name,units_name)
           write(*,*) "MLD: "//trim(varname)

           !Monthly avearge
           call month_ave(iday,nday,im,jm,k,mask,mld_out(:,:,ivar),mldave_out(:,:,ivar),mldpass_out(:,:,ivar))

           !Write daily data
           call make_dat2d_ncfile(ms,varname,long_name,units_name,iyr,imon,iday,im,jm)
           call write_dat2d_ncfile(ms,varname,iyr,imon,iday,im,jm,mld_out(:,:,ivar))

           !Write monthly-mean data
           if(iday == nday)then
              jday=0 !Dummy
              call make_dat2d_ncfile(ms,varname,long_name,units_name,iyr,imon,jday,im,jm)
              call write_dat2d_ncfile(ms,varname,iyr,imon,jday,im,jm,mldave_out(:,:,ivar))
           end if

        end do !ivar

        !---MLT-related variables
        !Temprature advection correction terms
        call adv_cor(im,jm,km,maskt,masku,maskv,dx,dy,sigmat,depw,ssh,u,v,t,xadv_cor,yadv_cor,zadv_cor)        
        do ivar=1,nvar_mlt

           !Varname
           call varname_mlt(ivar,varname,long_name,units_name)
           write(*,*) "MLT: "//trim(varname)

           !Read data
           if(ivar == 1)then             !dTmix/dt = MLT(t1)-MLT(t0)
              call dmltdt(im,jm,mask,mlt_t0,mlt_t1,mlt_out(:,:,ivar))
              mlt_out(:,:,nvar_mlt)=mlt_out(:,:,ivar)
           else if(ivar == 2)then        !Surface flux
              call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,k,mask,dat2d)
              dat3d(:,:,1)=dat2d(:,:)
              dat3d(:,:,2:km)=0.d0
              call mld_ave(im,jm,km,id_mld_t1,mask,dzt,dat3d,mlt_out(:,:,ivar))              
           else if(ivar == 10)then       !Detrainment/Entrainment
              !skip
           else if(ivar == nvar_mlt)then !Residual
              !skip
           else                          !Others
              
              call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,km,mask,dat3d)

              !Advection correction
              if(ivar == 4)then
                 dat3d(:,:,:)=dat3d(:,:,:)+xadv_cor(:,:,:)
              else if(ivar == 5)then
                 dat3d(:,:,:)=dat3d(:,:,:)+yadv_cor(:,:,:)
              else if(ivar == 6)then
                 dat3d(:,:,:)=dat3d(:,:,:)+zadv_cor(:,:,:)
              end if              
              
              call mld_ave(im,jm,km,id_mld_t1,mask,dzt,dat3d,mlt_out(:,:,ivar))
              
           end if

           
           !Residual
           if(ivar == 1 .or. ivar == nvar_mlt)then
              !skip
           else
              call residual(im,jm,mask,mlt_out(:,:,ivar),mlt_out(:,:,nvar_mlt))
           end if

           call month_ave(iday,nday,im,jm,k,mask,mlt_out(:,:,ivar),mltave_out(:,:,ivar),mltpass_out(:,:,ivar))

           !---Write NetCDF data
           call make_dat2d_ncfile(ms,varname,long_name,units_name,iyr,imon,iday,im,jm)
           call write_dat2d_ncfile(ms,varname,iyr,imon,iday,im,jm,mlt_out(:,:,ivar))

           if(iday == nday)then
              jday=0 !Dummy
              call make_dat2d_ncfile(ms,varname,long_name,units_name,iyr,imon,jday,im,jm)
              call write_dat2d_ncfile(ms,varname,iyr,imon,jday,im,jm,mltave_out(:,:,ivar))
           end if

        end do !ivar

        !---MLS-related variables
        !Salinity advection correction terms
        call adv_cor(im,jm,km,maskt,masku,maskv,dx,dy,sigmat,depw,ssh,u,v,s,xadv_cor,yadv_cor,zadv_cor)
        
        do ivar=1,nvar_mls

           !Varname
           call varname_mls(ivar,varname,long_name,units_name)
           write(*,*) "MLS: "//trim(varname)

           !Read data
           if(ivar == 1)then !dSmix/dt = MLS(t1)-MLS(t0)
              call dmltdt(im,jm,mask,mls_t0,mls_t1,mls_out(:,:,ivar))
              mls_out(:,:,nvar_mls)=mls_out(:,:,ivar)
           else if(ivar == 2)then !Surface flux
              call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,k,mask,dat2d)
              dat3d(:,:,1)=dat2d(:,:)
              dat3d(:,:,2:km)=0.d0
              call mld_ave(im,jm,km,id_mld_t1,mask,dzt,dat3d,mls_out(:,:,ivar))
           else if(ivar == 9)then !Detrainment/Entrainment
              !skip
           else if(ivar == nvar_mls)then !Residual
              !skip
           else
              
              call read_anal(dir,letkf,region,ms,imem,varname,iyr,imon,iday,im,jm,km,mask,dat3d)

              !Advection correction
              if(ivar == 3)then
                 dat3d(:,:,:)=dat3d(:,:,:)+xadv_cor(:,:,:)
              else if(ivar == 4)then
                 dat3d(:,:,:)=dat3d(:,:,:)+yadv_cor(:,:,:)
              else if(ivar == 5)then
                 dat3d(:,:,:)=dat3d(:,:,:)+zadv_cor(:,:,:)
              end if
              
              call mld_ave(im,jm,km,id_mld_t1,mask,dzt,dat3d,mls_out(:,:,ivar))
              
           end if

           
           !Residual
           if(ivar == 1 .or. ivar == nvar_mls)then
              !skip
           else
              call residual(im,jm,mask,mls_out(:,:,ivar),mls_out(:,:,nvar_mls))
           end if

           call month_ave(iday,nday,im,jm,k,mask,mls_out(:,:,ivar),mlsave_out(:,:,ivar),mlspass_out(:,:,ivar))

           !Write NetCDF data
           call make_dat2d_ncfile(ms,varname,long_name,units_name,iyr,imon,iday,im,jm)
           call write_dat2d_ncfile(ms,varname,iyr,imon,iday,im,jm,mls_out(:,:,ivar))

           if(iday == nday)then
              jday=0 !Dummy
              call make_dat2d_ncfile(ms,varname,long_name,units_name,iyr,imon,jday,im,jm)
              call write_dat2d_ncfile(ms,varname,iyr,imon,jday,im,jm,mlsave_out(:,:,ivar))
           end if

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
     deallocate(dat2d,dat3d)  
     deallocate(mld_out,mldave_out,mldpass_out)
     deallocate(mlt_out,mltave_out,mltpass_out)
     deallocate(mls_out,mlsave_out,mlspass_out)
     deallocate(xadv_cor,yadv_cor,zadv_cor)

  end if
     
  write(*,*) "=== End: Creat FTP file ==="
  
end program main

