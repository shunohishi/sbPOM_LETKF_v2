program main

  use setting
  use mod_julian
  use mod_gridinfo, im_lora => im, jm_lora => jm, km_lora => km
  use mod_read_glorys025, im_g025 => im, jm_g025 => jm, km_g025 => km
  use mod_read_db
  use mod_io
  use mod_rmiss
  use netcdf
  implicit none

  !---Common
  integer ijul,sjul,ejul
  integer iyr,imon,iday
  integer syr,smon,sday
  integer eyr,emon,eday
  integer k
  integer idat_a
  integer inum
  integer status,ncid
  
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
  
  !---Analysis in obs space (Bilinear Interpolation)
  integer,allocatable :: idxt(:),idyt(:)
  integer,allocatable :: idxu(:),idyu(:)
  integer,allocatable :: idxv(:),idyv(:)  

  real(kind = 8),allocatable :: hu_a(:),hv_a(:),ht_a(:)
  real(kind = 8),allocatable :: husprd_a(:),hvsprd_a(:),htsprd_a(:)

  write(*,*) "### START: Make data in observation space ###"

  !---Read start and end dates
  call read_argument(syr,smon,sday,eyr,emon,eday)
  write(*,*) "Start date:",syr,smon,sday
  write(*,*) "End date:",eyr,emon,eday
  
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
        ncid=0
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
                   & ncid,inum)

           end do !idb

           !Deallocate
           deallocate(ijul_o)
           deallocate(idxt,idxu,idxv)
           deallocate(idyt,idyu,idyv)
           deallocate(ht_a,hu_a,hv_a)           
           deallocate(htsprd_a,husprd_a,hvsprd_a)           
           call deallocate_db(iyr_o,imon_o,iday_o,lon_o,lat_o,u_o,v_o,t_o)

        end do !ifile_o

        !Close
        status=nf90_close(ncid)
        
     end do    !ijul

     deallocate(lont_a,lonu_a,lonv_a)
     deallocate(latt_a,latu_a,latv_a)
     deallocate(maskt_a,masku_a,maskv_a)
     deallocate(t_a,u_a,v_a)
     deallocate(tsprd_a,usprd_a,vsprd_a)
     
  end do          !idat_a
  
  !---Deallocate
  call deallocate_info(filename_o, &
       & iyr_min_o,imon_min_o,iday_min_o, &
       & iyr_max_o,imon_max_o,iday_max_o, &
       & lon_min_o,lon_max_o,lat_min_o,lat_max_o)

  write(*,*) "### SUCCESS: Make data in observation space ###"
  
end program main

