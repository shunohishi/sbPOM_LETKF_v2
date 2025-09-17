program main

  use setting
  use mod_julian
  use mod_bin
  use mod_static
  use mod_gridinfo, im_lora => im, jm_lora => jm, km_lora => km
  use mod_read_glorys025, im_g025 => im, jm_g025 => jm, km_g025 => km
  use mod_read_db
  use mod_io
  use mod_rmiss
  implicit none

  !---Common
  integer ijul,sjul,ejul
  integer iyr,imon,iday
  integer syr,smon,sday
  integer eyr,emon,eday
  integer k
  integer idat_a
  integer inum
  
  !---Analysis
  integer imem
  integer im_a,jm_a,km_a

  real(kind = 8),allocatable :: lonu_a(:),latu_a(:),masku_a(:,:)
  real(kind = 8),allocatable :: lonv_a(:),latv_a(:),maskv_a(:,:)
  real(kind = 8),allocatable :: lont_a(:),latt_a(:),maskt_a(:,:)

  real(kind = 8),allocatable :: u_a(:,:),v_a(:,:),t_a(:,:)
  real(kind = 8),allocatable :: usprd_a(:,:),vsprd_a(:,:),tsprd_a(:,:)
    
  !---Drifter Buoy (DB)
  integer ifile_o,nfile_o
  integer idb,ndb
  integer iobs,nobs

  !Info
  integer,allocatable :: ijul_min_o(:),ijul_max_o(:)
  integer,allocatable :: iyr_min_o(:),imon_min_o(:),iday_min_o(:)
  integer,allocatable :: iyr_max_o(:),imon_max_o(:),iday_max_o(:)
  
  real(kind = 8),allocatable :: lon_min_o(:),lon_max_o(:)
  real(kind = 8),allocatable :: lat_min_o(:),lat_max_o(:)

  !Data
  integer,allocatable :: ijul_o(:)
  integer,allocatable :: iyr_o(:,:),imon_o(:,:),iday_o(:,:)

  real(kind = 8),allocatable :: lon_o(:,:),lat_o(:,:)
  real(kind = 8),allocatable :: u_o(:,:),v_o(:,:),t_o(:,:)
  
  character(100),allocatable :: filename_o(:)

  !---Obs, space
  real(kind = 8),allocatable :: lon_o1d(:),lat_o1d(:)
  real(kind = 8),allocatable :: u_o1d(:),v_o1d(:),t_o1d(:)
  
  !---Analysis in obs space (Bilinear Interpolation)
  integer,allocatable :: idxt(:),idyt(:)
  integer,allocatable :: idxu(:),idyu(:)
  integer,allocatable :: idxv(:),idyv(:)  

  real(kind = 8),allocatable :: hu_a(:),hv_a(:),ht_a(:)
  real(kind = 8),allocatable :: husprd_a(:),hvsprd_a(:),htsprd_a(:)

  !---Static
  !Bin (ALL Ave)
  integer im_bin,jm_bin
  integer,allocatable :: unum_bin(:,:,:),vnum_bin(:,:,:),tnum_bin(:,:,:)
  real(kind = 8),allocatable :: lon_bin(:),lat_bin(:)
  real(kind = 8),allocatable :: ubias_bin(:,:,:),vbias_bin(:,:,:),tbias_bin(:,:,:)
  real(kind = 8),allocatable :: urmsd_bin(:,:,:),vrmsd_bin(:,:,:),trmsd_bin(:,:,:)
  real(kind = 8),allocatable :: usprd_bin(:,:,:),vsprd_bin(:,:,:),tsprd_bin(:,:,:)
  
  !Monthly
  integer,allocatable :: unum_mave(:,:,:),vnum_mave(:,:,:),tnum_mave(:,:,:)
  real(kind = 8),allocatable :: ubias_mave(:,:,:),vbias_mave(:,:,:),tbias_mave(:,:,:)
  real(kind = 8),allocatable :: urmsd_mave(:,:,:),vrmsd_mave(:,:,:),trmsd_mave(:,:,:)
  real(kind = 8),allocatable :: usprd_mave(:,:,:),vsprd_mave(:,:,:),tsprd_mave(:,:,:)
  
  !Yearly
  integer,allocatable :: unum_yave(:,:),vnum_yave(:,:),tnum_yave(:,:)
  real(kind = 8),allocatable :: ubias_yave(:,:),vbias_yave(:,:),tbias_yave(:,:)
  real(kind = 8),allocatable :: urmsd_yave(:,:),vrmsd_yave(:,:),trmsd_yave(:,:)
  real(kind = 8),allocatable :: usprd_yave(:,:),vsprd_yave(:,:),tsprd_yave(:,:)

  !ALL
  integer unum_ave(ndat_a),vnum_ave(ndat_a),tnum_ave(ndat_a)
  real(kind = 8) ubias_ave(ndat_a),vbias_ave(ndat_a),tbias_ave(ndat_a)
  real(kind = 8) urmsd_ave(ndat_a),vrmsd_ave(ndat_a),trmsd_ave(ndat_a)
  real(kind = 8) usprd_ave(ndat_a),vsprd_ave(ndat_a),tsprd_ave(ndat_a)

  write(*,*) "START: Validation vs. drifter buoy"

  !---Read start and end dates
  call read_argument(syr,smon,sday,eyr,emon,eday)
  write(*,*) "Start date:",syr,smon,sday
  write(*,*) "End date:",eyr,emon,eday

  !---Allocate
  allocate(unum_mave(ndat_a,12,syr:eyr),vnum_mave(ndat_a,12,syr:eyr),tnum_mave(ndat_a,12,syr:eyr))
  allocate(ubias_mave(ndat_a,12,syr:eyr),vbias_mave(ndat_a,12,syr:eyr),tbias_mave(ndat_a,12,syr:eyr))
  allocate(urmsd_mave(ndat_a,12,syr:eyr),vrmsd_mave(ndat_a,12,syr:eyr),trmsd_mave(ndat_a,12,syr:eyr))
  allocate(usprd_mave(ndat_a,12,syr:eyr),vsprd_mave(ndat_a,12,syr:eyr),tsprd_mave(ndat_a,12,syr:eyr))

  allocate(unum_yave(ndat_a,syr:eyr),vnum_yave(ndat_a,syr:eyr),tnum_yave(ndat_a,syr:eyr))
  allocate(ubias_yave(ndat_a,syr:eyr),vbias_yave(ndat_a,syr:eyr),tbias_yave(ndat_a,syr:eyr))
  allocate(urmsd_yave(ndat_a,syr:eyr),vrmsd_yave(ndat_a,syr:eyr),trmsd_yave(ndat_a,syr:eyr))
  allocate(usprd_yave(ndat_a,syr:eyr),vsprd_yave(ndat_a,syr:eyr),tsprd_yave(ndat_a,syr:eyr))
  
  !---Make DB filename
  call make_filename

  !---Read DB information 
  call read_info(nfile_o,filename_o, &
       & iyr_min_o,imon_min_o,iday_min_o, &
       & iyr_max_o,imon_max_o,iday_max_o, &
       & lon_min_o,lon_max_o,lat_min_o,lat_max_o)
  
  !---DB JULDAY
  allocate(ijul_min_o(nfile_o),ijul_max_o(nfile_o))
  do ifile_o=1,nfile_o
     call ymd_julian(iyr_min_o(ifile_o),imon_min_o(ifile_o),iday_min_o(ifile_o),ijul_min_o(ifile_o))
     call ymd_julian(iyr_max_o(ifile_o),imon_max_o(ifile_o),iday_max_o(ifile_o),ijul_max_o(ifile_o))
  end do
  
  !---Make lon-lat bin
  call make_bin(slon_bin,elon_bin,slat_bin,elat_bin,dx_bin,dy_bin, &
       & im_bin,jm_bin,lon_bin,lat_bin)
  allocate(unum_bin(im_bin,jm_bin,ndat_a),vnum_bin(im_bin,jm_bin,ndat_a),tnum_bin(im_bin,jm_bin,ndat_a))
  allocate(ubias_bin(im_bin,jm_bin,ndat_a),vbias_bin(im_bin,jm_bin,ndat_a),tbias_bin(im_bin,jm_bin,ndat_a))
  allocate(urmsd_bin(im_bin,jm_bin,ndat_a),vrmsd_bin(im_bin,jm_bin,ndat_a),trmsd_bin(im_bin,jm_bin,ndat_a))
  allocate(usprd_bin(im_bin,jm_bin,ndat_a),vsprd_bin(im_bin,jm_bin,ndat_a),tsprd_bin(im_bin,jm_bin,ndat_a))

  !---Initialization
  do idat_a=1,ndat_a
     call static_bin_ini(im_bin,jm_bin,unum_bin(:,:,idat_a),ubias_bin(:,:,idat_a),urmsd_bin(:,:,idat_a),usprd_bin(:,:,idat_a))
     call static_bin_ini(im_bin,jm_bin,vnum_bin(:,:,idat_a),vbias_bin(:,:,idat_a),vrmsd_bin(:,:,idat_a),vsprd_bin(:,:,idat_a))
     call static_bin_ini(im_bin,jm_bin,tnum_bin(:,:,idat_a),tbias_bin(:,:,idat_a),trmsd_bin(:,:,idat_a),tsprd_bin(:,:,idat_a))
  end do!idat_a

  do iyr=syr,eyr
     
     do idat_a=1,ndat_a
        call static_ave_ini(unum_yave(idat_a,iyr),ubias_yave(idat_a,iyr),urmsd_yave(idat_a,iyr),usprd_yave(idat_a,iyr))
        call static_ave_ini(vnum_yave(idat_a,iyr),vbias_yave(idat_a,iyr),vrmsd_yave(idat_a,iyr),vsprd_yave(idat_a,iyr))
        call static_ave_ini(tnum_yave(idat_a,iyr),tbias_yave(idat_a,iyr),trmsd_yave(idat_a,iyr),tsprd_yave(idat_a,iyr))
     end do !idat_a
  
     do imon=1,12
        do idat_a=1,ndat_a
           call static_ave_ini(unum_mave(idat_a,imon,iyr),ubias_mave(idat_a,imon,iyr),urmsd_mave(idat_a,imon,iyr),usprd_mave(idat_a,imon,iyr))
           call static_ave_ini(vnum_mave(idat_a,imon,iyr),vbias_mave(idat_a,imon,iyr),vrmsd_mave(idat_a,imon,iyr),vsprd_mave(idat_a,imon,iyr))
           call static_ave_ini(tnum_mave(idat_a,imon,iyr),tbias_mave(idat_a,imon,iyr),trmsd_mave(idat_a,imon,iyr),tsprd_mave(idat_a,imon,iyr))
        end do !idat_a
     end do    !imon
     
  end do       !iyr

  do idat_a=1,ndat_a
     call static_ave_ini(unum_ave(idat_a),ubias_ave(idat_a),urmsd_ave(idat_a),usprd_ave(idat_a))
     call static_ave_ini(vnum_ave(idat_a),vbias_ave(idat_a),vrmsd_ave(idat_a),vsprd_ave(idat_a))
     call static_ave_ini(tnum_ave(idat_a),tbias_ave(idat_a),trmsd_ave(idat_a),tsprd_ave(idat_a))
  end do

  !--- Main Loop ---!
  do idat_a=1,ndat_a

     !---Grid size (To Be Modified)
     if(idat_a == 1)then
        im_a=im_lora
        jm_a=jm_lora
        km_a=km_lora
     else if(idat_a == 2 .or. idat_a == 3 .or. idat_a == 4)then
        im_a=im_g025
        jm_a=jm_g025
        km_a=km_g025
     else
        write(*,*) "***Error: Incorecot idat_a => ",idat_a
        stop
     end if        

     allocate(lont_a(im_a),lonu_a(im_a),lonv_a(im_a))
     allocate(latt_a(jm_a),latu_a(jm_a),latv_a(jm_a))
     allocate(maskt_a(im_a,jm_a),masku_a(im_a,jm_a),maskv_a(im_a,jm_a))
     allocate(t_a(im_a,jm_a),u_a(im_a,jm_a),v_a(im_a,jm_a))
     allocate(tsprd_a(im_a,jm_a),usprd_a(im_a,jm_a),vsprd_a(im_a,jm_a))
     
     !---Read DA Grid (To be modified)
     write(*,*) "Read grid"
     call read_grid(idat_a,im_a,jm_a,km_a,lont_a,lonu_a,lonv_a,latt_a,latu_a,latv_a,maskt_a,masku_a,maskv_a)
     call read_surface_data(idat_a,syr,smon,sday,im_a,jm_a,km_a,maskt_a,masku_a,maskv_a, &
          & t_a,u_a,v_a,tsprd_a,usprd_a,vsprd_a)
     
     !---Validation loop
     call ymd_julian(syr,smon,sday,sjul)
     call ymd_julian(eyr,emon,eday,ejul)  
     do ijul=sjul,ejul

        !---Initialization
        inum=0
        call julian_ymd(ijul,iyr,imon,iday)
        write(*,*) idat_a,iyr,imon,iday

        !---Read analysis surface data
        write(*,*) "Read surface data"
        call read_surface_data(idat_a,iyr,imon,iday,im_a,jm_a,km_a,maskt_a,masku_a,maskv_a, &
             & t_a,u_a,v_a,tsprd_a,usprd_a,vsprd_a)
        
        !---vs. DB
        write(*,*) "Match up with drifter buoy"
        do ifile_o=1,nfile_o
           
           !Check range
           if(ijul < ijul_min_o(ifile_o) .or. ijul_max_o(ifile_o) < ijul) cycle
           if(lon_max_o(ifile_o) < lonu_a(1) .or. lont_a(im_a) < lon_min_o(ifile_o)) cycle
           if(lat_max_o(ifile_o) < latv_a(1) .or. latt_a(jm_a) < lat_min_o(ifile_o)) cycle

           !---Read buoy data
           call read_db_1d(filename_o(ifile_o),ndb,nobs,iyr_o,imon_o,iday_o,lon_o,lat_o,u_o,v_o,t_o)
           if(ndb == 0 .or. nobs == 0) cycle
           
           allocate(ijul_o(nobs))
           allocate(idxt(nobs),idxu(nobs),idxv(nobs))
           allocate(idyt(nobs),idyu(nobs),idyv(nobs))
           allocate(hu_a(nobs),hv_a(nobs),ht_a(nobs))
           allocate(husprd_a(nobs),hvsprd_a(nobs),htsprd_a(nobs))

           !---Make obs. space data
           do idb=1,ndb
              
              !---ijul
              do iobs=1,nobs
                 call ymd_julian(iyr_o(iobs,idb),imon_o(iobs,idb),iday_o(iobs,idb),ijul_o(iobs))
              end do
           
              !---ID
              call cal_idlon(im_a,lont_a,nobs,lon_o(:,idb),idxt)
              call cal_idlon(im_a,lonu_a,nobs,lon_o(:,idb),idxu)
              call cal_idlon(im_a,lonv_a,nobs,lon_o(:,idb),idxv)
              call cal_idlat(jm_a,latt_a,nobs,lat_o(:,idb),idyt)
              call cal_idlat(jm_a,latu_a,nobs,lat_o(:,idb),idyu)
              call cal_idlat(jm_a,latv_a,nobs,lat_o(:,idb),idyv)
           
              !---Project to obs. space
              call bilinear_interpolation_2d &
                   & (im_a,jm_a,lont_a,latt_a,t_a(:,:),maskt_a, &
                   &  nobs,lon_o(:,idb),lat_o(:,idb), &
                   &  idxt(:),idyt(:),ht_a(:))                 

              call bilinear_interpolation_2d &
                   & (im_a,jm_a,lonu_a,latu_a,u_a(:,:),masku_a, &
                   &  nobs,lon_o(:,idb),lat_o(:,idb), &
                   &  idxu(:),idyu(:),hu_a(:))

              call bilinear_interpolation_2d &
                   & (im_a,jm_a,lonv_a,latv_a,v_a(:,:),maskv_a, &
                   &  nobs,lon_o(:,idb),lat_o(:,idb), &
                   &  idxv(:),idyv(:),hv_a(:))

              call bilinear_interpolation_2d &
                   & (im_a,jm_a,lont_a,latt_a,tsprd_a(:,:),maskt_a, &
                   &  nobs,lon_o(:,idb),lat_o(:,idb), &
                   &  idxt(:),idyt(:),htsprd_a(:))                 

              call bilinear_interpolation_2d &
                   & (im_a,jm_a,lonu_a,latu_a,usprd_a(:,:),masku_a, &
                   &  nobs,lon_o(:,idb),lat_o(:,idb), &
                   &  idxu(:),idyu(:),husprd_a(:))

              call bilinear_interpolation_2d &
                   & (im_a,jm_a,lonv_a,latv_a,vsprd_a(:,:),maskv_a, &
                   &  nobs,lon_o(:,idb),lat_o(:,idb), &
                   &  idxv(:),idyv(:),hvsprd_a(:))
              
              !---Write obs space data
              call write_obs(idat_a,ijul,nobs,ijul_o(:),lon_o(:,idb),lat_o(:,idb), &
                   & ht_a(:),hu_a(:),hv_a(:),htsprd_a(:),husprd_a(:),hvsprd_a(:),t_o(:,idb),u_o(:,idb),v_o(:,idb), &
                   & inum)

           end do !idb

           !Deallocate
           deallocate(ijul_o)
           deallocate(idxt,idxu,idxv)
           deallocate(idyt,idyu,idyv)
           deallocate(ht_a,hu_a,hv_a)           
           deallocate(htsprd_a,husprd_a,hvsprd_a)           
           call deallocate_db(iyr_o,imon_o,iday_o,lon_o,lat_o,u_o,v_o,t_o)

        end do !ifile_o

        !---Read
        write(*,*) "Read 1-day obs space data at ",iyr,imon,iday
        call read_obs(idat_a,iyr,imon,iday,nobs,lon_o1d,lat_o1d,ht_a,hu_a,hv_a,htsprd_a,husprd_a,hvsprd_a,t_o1d,u_o1d,v_o1d)
        write(*,*) "End: Read 1-day obs space data"
        
        if(nobs == 0) cycle

        !---Static
        write(*,*) "Conduct validation"
        !Bin
        call static_bin_add(nobs,lon_o1d,lat_o1d,hu_a,husprd_a,u_o1d, &
             & im_bin,jm_bin,lon_bin,lat_bin, &
             & unum_bin(:,:,idat_a),ubias_bin(:,:,idat_a),urmsd_bin(:,:,idat_a),usprd_bin(:,:,idat_a))
        call static_bin_add(nobs,lon_o1d,lat_o1d,hv_a,hvsprd_a,v_o1d, &
             & im_bin,jm_bin,lon_bin,lat_bin, &
             & vnum_bin(:,:,idat_a),vbias_bin(:,:,idat_a),vrmsd_bin(:,:,idat_a),vsprd_bin(:,:,idat_a))
        call static_bin_add(nobs,lon_o1d,lat_o1d,ht_a,htsprd_a,t_o1d, &
             & im_bin,jm_bin,lon_bin,lat_bin, &
             & tnum_bin(:,:,idat_a),tbias_bin(:,:,idat_a),trmsd_bin(:,:,idat_a),tsprd_bin(:,:,idat_a))

        !Monthly
        call static_ave_add(nobs,hu_a,husprd_a,u_o1d, &
             & unum_mave(idat_a,imon,iyr),ubias_mave(idat_a,imon,iyr),urmsd_mave(idat_a,imon,iyr),usprd_mave(idat_a,imon,iyr))
        call static_ave_add(nobs,hv_a,hvsprd_a,v_o1d, &
             & vnum_mave(idat_a,imon,iyr),vbias_mave(idat_a,imon,iyr),vrmsd_mave(idat_a,imon,iyr),vsprd_mave(idat_a,imon,iyr))
        call static_ave_add(nobs,ht_a,htsprd_a,t_o1d, &
             & tnum_mave(idat_a,imon,iyr),tbias_mave(idat_a,imon,iyr),trmsd_mave(idat_a,imon,iyr),tsprd_mave(idat_a,imon,iyr))

        !Yearly
        call static_ave_add(nobs,hu_a,husprd_a,u_o1d, &
             & unum_yave(idat_a,iyr),ubias_yave(idat_a,iyr),urmsd_yave(idat_a,iyr),usprd_yave(idat_a,iyr))
        call static_ave_add(nobs,hv_a,hvsprd_a,v_o1d, &
             & vnum_yave(idat_a,iyr),vbias_yave(idat_a,iyr),vrmsd_yave(idat_a,iyr),vsprd_yave(idat_a,iyr))
        call static_ave_add(nobs,ht_a,htsprd_a,t_o1d, &
             & tnum_yave(idat_a,iyr),tbias_yave(idat_a,iyr),trmsd_yave(idat_a,iyr),tsprd_yave(idat_a,iyr))

        !ALL
        call static_ave_add(nobs,hu_a,husprd_a,u_o1d, &
             & unum_ave(idat_a),ubias_ave(idat_a),urmsd_ave(idat_a),usprd_ave(idat_a))
        call static_ave_add(nobs,hv_a,hvsprd_a,v_o1d, &
             & vnum_ave(idat_a),vbias_ave(idat_a),vrmsd_ave(idat_a),vsprd_ave(idat_a))
        call static_ave_add(nobs,ht_a,htsprd_a,t_o1d, &
             & tnum_ave(idat_a),tbias_ave(idat_a),trmsd_ave(idat_a),tsprd_ave(idat_a))

        !---Deallocate
        call deallocate_obs(lon_o1d,lat_o1d, &
             & ht_a,hu_a,hv_a,htsprd_a,husprd_a,hvsprd_a,t_o1d,u_o1d,v_o1d)
        
     end do       !ijul

     deallocate(lont_a,lonu_a,lonv_a)
     deallocate(latt_a,latu_a,latv_a)
     deallocate(maskt_a,masku_a,maskv_a)
     deallocate(t_a,u_a,v_a)
     deallocate(tsprd_a,usprd_a,vsprd_a)
     
  end do          !idat_a

  !---Post process
  write(*,*) "Post process"
  !Bin
  do idat_a=1,ndat_a
     call static_bin_end(im_bin,jm_bin, &
          & unum_bin(:,:,idat_a),ubias_bin(:,:,idat_a),urmsd_bin(:,:,idat_a),usprd_bin(:,:,idat_a))
     call static_bin_end(im_bin,jm_bin, &
          & vnum_bin(:,:,idat_a),vbias_bin(:,:,idat_a),vrmsd_bin(:,:,idat_a),vsprd_bin(:,:,idat_a))
     call static_bin_end(im_bin,jm_bin, &
          & tnum_bin(:,:,idat_a),tbias_bin(:,:,idat_a),trmsd_bin(:,:,idat_a),tsprd_bin(:,:,idat_a))
  end do !idat_a

  !Ave
  do iyr=syr,eyr

     do idat_a=1,ndat_a
        call static_ave_end(unum_yave(idat_a,iyr),ubias_yave(idat_a,iyr),urmsd_yave(idat_a,iyr),usprd_yave(idat_a,iyr))
        call static_ave_end(vnum_yave(idat_a,iyr),vbias_yave(idat_a,iyr),vrmsd_yave(idat_a,iyr),vsprd_yave(idat_a,iyr))
        call static_ave_end(tnum_yave(idat_a,iyr),tbias_yave(idat_a,iyr),trmsd_yave(idat_a,iyr),tsprd_yave(idat_a,iyr))
     end do !idat_a
     
     do imon=1,12
        do idat_a=1,ndat_a
           call static_ave_end(unum_mave(idat_a,imon,iyr),ubias_mave(idat_a,imon,iyr),urmsd_mave(idat_a,imon,iyr),usprd_mave(idat_a,imon,iyr))
           call static_ave_end(vnum_mave(idat_a,imon,iyr),vbias_mave(idat_a,imon,iyr),vrmsd_mave(idat_a,imon,iyr),vsprd_mave(idat_a,imon,iyr))
           call static_ave_end(tnum_mave(idat_a,imon,iyr),tbias_mave(idat_a,imon,iyr),trmsd_mave(idat_a,imon,iyr),tsprd_mave(idat_a,imon,iyr))
        end do !idat_a
     end do    !imon

  end do !iyr 

  do idat_a=1,ndat_a
     call static_ave_end(unum_ave(idat_a),ubias_ave(idat_a),urmsd_ave(idat_a),usprd_ave(idat_a))
     call static_ave_end(vnum_ave(idat_a),vbias_ave(idat_a),vrmsd_ave(idat_a),vsprd_ave(idat_a))
     call static_ave_end(tnum_ave(idat_a),tbias_ave(idat_a),trmsd_ave(idat_a),tsprd_ave(idat_a))
  end do
  
  !---Write data
  call write_bin(im_bin,jm_bin,ndat_a,dx_bin,dy_bin,lon_bin,lat_bin, &
       unum_bin,ubias_bin,urmsd_bin,usprd_bin, &
       vnum_bin,vbias_bin,vrmsd_bin,vsprd_bin, &
       tnum_bin,tbias_bin,trmsd_bin,tsprd_bin)

  call write_ave(syr,eyr,ndat_a, &
     & unum_mave,ubias_mave,urmsd_mave,usprd_mave, &
     & vnum_mave,vbias_mave,vrmsd_mave,vsprd_mave, &
     & tnum_mave,tbias_mave,trmsd_mave,tsprd_mave, &
     & unum_yave,ubias_yave,urmsd_yave,usprd_yave, &
     & vnum_yave,vbias_yave,vrmsd_yave,vsprd_yave, &
     & tnum_yave,tbias_yave,trmsd_yave,tsprd_yave, &
     & unum_ave,ubias_ave,urmsd_ave,usprd_ave, &
     & vnum_ave,vbias_ave,vrmsd_ave,vsprd_ave, &
     & tnum_ave,tbias_ave,trmsd_ave,tsprd_ave)
  
  !---Deallocate
  call deallocate_info(filename_o, &
       & iyr_min_o,imon_min_o,iday_min_o, &
       & iyr_max_o,imon_max_o,iday_max_o, &
       & lon_min_o,lon_max_o,lat_min_o,lat_max_o)

  deallocate(unum_mave,vnum_mave,tnum_mave)
  deallocate(ubias_mave,vbias_mave,tbias_mave)
  deallocate(urmsd_mave,vrmsd_mave,trmsd_mave)
  deallocate(usprd_mave,vsprd_mave,tsprd_mave)

  deallocate(unum_yave,vnum_yave,tnum_yave)
  deallocate(ubias_yave,vbias_yave,tbias_yave)
  deallocate(urmsd_yave,vrmsd_yave,trmsd_yave)
  deallocate(usprd_yave,vsprd_yave,tsprd_yave)

  deallocate(unum_bin,vnum_bin,tnum_bin) 
  deallocate(ubias_bin,vbias_bin,tbias_bin)
  deallocate(urmsd_bin,vrmsd_bin,trmsd_bin)
  deallocate(usprd_bin,vsprd_bin,tsprd_bin)

  write(*,*) "SUCCESS: Validation vs. drifter buoy"
  
end program main

