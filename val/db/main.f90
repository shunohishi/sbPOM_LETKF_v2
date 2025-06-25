module setting

  !---Date
  integer,parameter :: syr=2021,smon=1,sday=1
  integer,parameter :: eyr=2023,emon=12,eday=31

  !---Box information
  real(kind = 8),parameter :: slon_bin=0.d0,elon_bin=360.d0
  real(kind = 8),parameter :: slat_bin=-90.d0,elat_bin=90.d0
  real(kind = 8),parameter :: dx_bin=2.5d0,dy_bin=2.5d0
  
  !---Analysis information
  integer,parameter :: ndat_a=1 !The number of analyses

  character(10),parameter :: dir="CHL"
  character(10),parameter :: letkf="letkf"
  character(10),parameter :: region="chl"
  character(10),parameter :: ms="mean"
  
end module setting

!-------------------------------------------------------------------------------------

program main

  use setting
  use mod_julian
  use mod_bin
  use mod_static
  use mod_gridinfo, im_a => im, jm_a => jm, km_a => km
  use mod_read_lora
  use mod_read_db
  use mod_rmiss
  implicit none

  !---Common
  integer ijul,sjul,ejul
  integer iyr,imon,iday
  integer k
  integer idat_a
  
  !---Analysis
  integer imem
  
  real(kind = 8) lonu_a(im_a),latu_a(jm_a)
  real(kind = 8) lonv_a(im_a),latv_a(jm_a)
  real(kind = 8) lont_a(im_a),latt_a(jm_a)
  real(kind = 8) maskt_a(im_a,jm_a),masku_a(im_a,jm_a),maskv_a(im_a,jm_a)
  
  real(kind = 8) u_a(im_a,jm_a,ndat_a),v_a(im_a,jm_a,ndat_a),t_a(im_a,jm_a,ndat_a)
  
  real(kind =8) null2d_a(im_a,jm_a),null3d_a(im_a,jm_a,km_a)
  
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

  !---Static
  !Bin (ALL Ave)
  integer im_bin,jm_bin
  integer,allocatable :: unum_bin(:,:,:),vnum_bin(:,:,:),tnum_bin(:,:,:)

  real(kind = 8),allocatable :: lon_bin(:),lat_bin(:)
  real(kind = 8),allocatable :: ubias_bin(:,:,:),vbias_bin(:,:,:),tbias_bin(:,:,:)
  real(kind = 8),allocatable :: urmsd_bin(:,:,:),vrmsd_bin(:,:,:),trmsd_bin(:,:,:)
  
  !Monthly
  integer unum_mave(ndat_a,12,syr:eyr),vnum_mave(ndat_a,12,syr:eyr),tnum_mave(ndat_a,12,syr:eyr)

  real(kind = 8) ubias_mave(ndat_a,12,syr:eyr),vbias_mave(ndat_a,12,syr:eyr),tbias_mave(ndat_a,12,syr:eyr)
  real(kind = 8) urmsd_mave(ndat_a,12,syr:eyr),vrmsd_mave(ndat_a,12,syr:eyr),trmsd_mave(ndat_a,12,syr:eyr)
  
  !Yearly
  integer unum_yave(ndat_a,syr:eyr),vnum_yave(ndat_a,syr:eyr),tnum_yave(ndat_a,syr:eyr)

  real(kind = 8) ubias_yave(ndat_a,syr:eyr),vbias_yave(ndat_a,syr:eyr),tbias_yave(ndat_a,syr:eyr)
  real(kind = 8) urmsd_yave(ndat_a,syr:eyr),vrmsd_yave(ndat_a,syr:eyr),trmsd_yave(ndat_a,syr:eyr)

  !ALL
  integer unum_ave(ndat_a),vnum_ave(ndat_a),tnum_ave(ndat_a)

  real(kind = 8) ubias_ave(ndat_a),vbias_ave(ndat_a),tbias_ave(ndat_a)
  real(kind = 8) urmsd_ave(ndat_a),vrmsd_ave(ndat_a),trmsd_ave(ndat_a)
  
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

  !---Read DA system grid
  call read_grid(dir,lont_a,lonu_a,lonv_a, &
       & latt_a,latu_a,latv_a, &
       & null3d_a,null3d_a,null3d_a, &
       & maskt_a,masku_a,maskv_a)
  
  !---Make lon-lat bin
  call make_bin(slon_bin,elon_bin,slat_bin,elat_bin,dx_bin,dy_bin, &
       & im_bin,jm_bin,lon_bin,lat_bin)
  allocate(unum_bin(im_bin,jm_bin,ndat_a),vnum_bin(im_bin,jm_bin,ndat_a),tnum_bin(im_bin,jm_bin,ndat_a))
  allocate(ubias_bin(im_bin,jm_bin,ndat_a),vbias_bin(im_bin,jm_bin,ndat_a),tbias_bin(im_bin,jm_bin,ndat_a))
  allocate(urmsd_bin(im_bin,jm_bin,ndat_a),vrmsd_bin(im_bin,jm_bin,ndat_a),trmsd_bin(im_bin,jm_bin,ndat_a))

  !---Initialization
  do idat_a=1,ndat_a
     call static_bin_ini(im_bin,jm_bin,unum_bin(:,:,idat_a),ubias_bin(:,:,idat_a),urmsd_bin(:,:,idat_a))
     call static_bin_ini(im_bin,jm_bin,vnum_bin(:,:,idat_a),vbias_bin(:,:,idat_a),vrmsd_bin(:,:,idat_a))
     call static_bin_ini(im_bin,jm_bin,tnum_bin(:,:,idat_a),tbias_bin(:,:,idat_a),trmsd_bin(:,:,idat_a))
  end do!idat_a

  do iyr=syr,eyr
     
     do idat_a=1,ndat_a
        call static_ave_ini(unum_yave(idat_a,iyr),ubias_yave(idat_a,iyr),urmsd_yave(idat_a,iyr))
        call static_ave_ini(vnum_yave(idat_a,iyr),vbias_yave(idat_a,iyr),vrmsd_yave(idat_a,iyr))
        call static_ave_ini(tnum_yave(idat_a,iyr),tbias_yave(idat_a,iyr),trmsd_yave(idat_a,iyr))
     end do !idat_a
  
     do imon=1,12
        do idat_a=1,ndat_a
           call static_ave_ini(unum_mave(idat_a,imon,iyr),ubias_mave(idat_a,imon,iyr),urmsd_mave(idat_a,imon,iyr))
           call static_ave_ini(vnum_mave(idat_a,imon,iyr),vbias_mave(idat_a,imon,iyr),vrmsd_mave(idat_a,imon,iyr))
           call static_ave_ini(tnum_mave(idat_a,imon,iyr),tbias_mave(idat_a,imon,iyr),trmsd_mave(idat_a,imon,iyr))
        end do !idat_a
     end do !imon
     
  end do    !iyr

  do idat_a=1,ndat_a
     call static_ave_ini(unum_ave(idat_a),ubias_ave(idat_a),urmsd_ave(idat_a))
     call static_ave_ini(vnum_ave(idat_a),vbias_ave(idat_a),vrmsd_ave(idat_a))
     call static_ave_ini(tnum_ave(idat_a),tbias_ave(idat_a),trmsd_ave(idat_a))
  end do
  
  !---Validation loop
  call ymd_julian(syr,smon,sday,sjul)
  call ymd_julian(eyr,emon,eday,ejul)
  do ijul=sjul,ejul

     call julian_ymd(ijul,iyr,imon,iday)
     write(*,*) iyr,imon,iday

     !---Read analysis data
     do idat_a=1,ndat_a

        imem=0
        k=1
        if(idat_a == 1)then
           call read_anal(dir,letkf,region,ms,imem,"t",iyr,imon,iday,im_a,jm_a,k,maskt_a,t_a(:,:,idat_a))
           call read_anal(dir,letkf,region,ms,imem,"u",iyr,imon,iday,im_a,jm_a,k,masku_a,u_a(:,:,idat_a))
           call read_anal(dir,letkf,region,ms,imem,"v",iyr,imon,iday,im_a,jm_a,k,maskv_a,v_a(:,:,idat_a))
           !*** Space for additional data ***
        end if
           
     end do !idat
        
     !---vs. DB
     do ifile_o=1,nfile_o
        
        !Check range
        if(ijul < ijul_min_o(ifile_o) .or. ijul_max_o(ifile_o) < ijul) cycle
        if(lon_max_o(ifile_o) < lonu_a(1) .or. lont_a(im_a) < lon_min_o(ifile_o)) cycle
        if(lat_max_o(ifile_o) < latv_a(1) .or. latt_a(jm_a) < lat_min_o(ifile_o)) cycle

        !---Read buoy data
        call read_db_1d(filename_o(ifile_o),ndb,nobs,iyr_o,imon_o,iday_o,lon_o,lat_o,u_o,v_o,t_o)
        if(ndb == 0 .or. nobs == 0) cycle

        do idb=1,ndb

           allocate(ijul_o(nobs))
           allocate(idxt(nobs),idxu(nobs),idxv(nobs))
           allocate(idyt(nobs),idyu(nobs),idyv(nobs))
           allocate(hu_a(nobs),hv_a(nobs),ht_a(nobs))

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
           
           do idat_a=1,ndat_a

              !---Project to obs. space              
              call bilinear_interpolation_2d &
                   & (im_a,jm_a,lonu_a,latu_a,u_a(:,:,idat_a),masku_a, &
                   &  nobs,lon_o(:,idb),lat_o(:,idb), &
                   &  idxu(:),idyu(:),hu_a(:))

              call bilinear_interpolation_2d &
                   & (im_a,jm_a,lonv_a,latv_a,v_a(:,:,idat_a),maskv_a, &
                   &  nobs,lon_o(:,idb),lat_o(:,idb), &
                   &  idxv(:),idyv(:),hv_a(:))
              
              call bilinear_interpolation_2d &
                   & (im_a,jm_a,lont_a,latt_a,t_a(:,:,idat_a),maskt_a, &
                   &  nobs,lon_o(:,idb),lat_o(:,idb), &
                   &  idxt(:),idyt(:),ht_a(:))                 
              
              !---Static
              !Bin
              call static_bin_add(ijul,nobs,ijul_o(:),lon_o(:,idb),lat_o(:,idb),hu_a(:),u_o(:,idb), &
                   & im_bin,jm_bin,lon_bin,lat_bin,unum_bin(:,:,idat_a),ubias_bin(:,:,idat_a),urmsd_bin(:,:,idat_a))
              call static_bin_add(ijul,nobs,ijul_o(:),lon_o(:,idb),lat_o(:,idb),hv_a(:),v_o(:,idb), &
                   & im_bin,jm_bin,lon_bin,lat_bin,vnum_bin(:,:,idat_a),vbias_bin(:,:,idat_a),vrmsd_bin(:,:,idat_a))
              call static_bin_add(ijul,nobs,ijul_o(:),lon_o(:,idb),lat_o(:,idb),ht_a(:),t_o(:,idb), &
                   & im_bin,jm_bin,lon_bin,lat_bin,tnum_bin(:,:,idat_a),tbias_bin(:,:,idat_a),trmsd_bin(:,:,idat_a))

              !Monthly
              call static_ave_add(ijul,nobs,ijul_o(:),hu_a(:),u_o(:,idb), &
                   & unum_mave(idat_a,imon,iyr),ubias_mave(idat_a,imon,iyr),urmsd_mave(idat_a,imon,iyr))
              call static_ave_add(ijul,nobs,ijul_o(:),hv_a(:),v_o(:,idb), &
                   & vnum_mave(idat_a,imon,iyr),vbias_mave(idat_a,imon,iyr),vrmsd_mave(idat_a,imon,iyr))
              call static_ave_add(ijul,nobs,ijul_o(:),ht_a(:),t_o(:,idb), &
                   & tnum_mave(idat_a,imon,iyr),tbias_mave(idat_a,imon,iyr),trmsd_mave(idat_a,imon,iyr))
              
              !Yearly
              call static_ave_add(ijul,nobs,ijul_o(:),hu_a(:),u_o(:,idb), &
                   & unum_yave(idat_a,iyr),ubias_yave(idat_a,iyr),urmsd_yave(idat_a,iyr))
              call static_ave_add(ijul,nobs,ijul_o(:),hv_a(:),v_o(:,idb), &
                   & vnum_yave(idat_a,iyr),vbias_yave(idat_a,iyr),vrmsd_yave(idat_a,iyr))
              call static_ave_add(ijul,nobs,ijul_o(:),ht_a(:),t_o(:,idb), &
                   & tnum_yave(idat_a,iyr),tbias_yave(idat_a,iyr),trmsd_yave(idat_a,iyr))
                            
              !ALL
              call static_ave_add(ijul,nobs,ijul_o(:),hu_a(:),u_o(:,idb), &
                   & unum_ave(idat_a),ubias_ave(idat_a),urmsd_ave(idat_a))
              call static_ave_add(ijul,nobs,ijul_o(:),hv_a(:),v_o(:,idb), &
                   & vnum_ave(idat_a),vbias_ave(idat_a),vrmsd_ave(idat_a))
              call static_ave_add(ijul,nobs,ijul_o(:),ht_a(:),t_o(:,idb), &
                   & tnum_ave(idat_a),tbias_ave(idat_a),trmsd_ave(idat_a))
              
           end do !idat_a

           !Deallocate
           deallocate(ijul_o)
           deallocate(idxt,idxu,idxv)
           deallocate(idyt,idyu,idyv)
           deallocate(hu_a,hv_a,ht_a)
           
        end do    !idb
        
        call deallocate_db(iyr_o,imon_o,iday_o,lon_o,lat_o,u_o,v_o,t_o)
        
     end do !ifile_o
  end do    !ijul

  !---Post process
  !Bin
  do idat_a=1,ndat_a
     call static_bin_end(im_bin,jm_bin,unum_bin(:,:,idat_a),ubias_bin(:,:,idat_a),urmsd_bin(:,:,idat_a))
     call static_bin_end(im_bin,jm_bin,vnum_bin(:,:,idat_a),vbias_bin(:,:,idat_a),vrmsd_bin(:,:,idat_a))
     call static_bin_end(im_bin,jm_bin,tnum_bin(:,:,idat_a),tbias_bin(:,:,idat_a),trmsd_bin(:,:,idat_a))
  end do !idat_a

  !Ave
  do iyr=syr,eyr

     do idat_a=1,ndat_a
        call static_ave_end(unum_yave(idat_a,iyr),ubias_yave(idat_a,iyr),urmsd_yave(idat_a,iyr))
        call static_ave_end(vnum_yave(idat_a,iyr),vbias_yave(idat_a,iyr),vrmsd_yave(idat_a,iyr))
        call static_ave_end(tnum_yave(idat_a,iyr),tbias_yave(idat_a,iyr),trmsd_yave(idat_a,iyr))
     end do !idat_a
     
     do imon=1,12
        do idat_a=1,ndat_a
           call static_ave_end(unum_mave(idat_a,imon,iyr),ubias_mave(idat_a,imon,iyr),urmsd_mave(idat_a,imon,iyr))
           call static_ave_end(vnum_mave(idat_a,imon,iyr),vbias_mave(idat_a,imon,iyr),vrmsd_mave(idat_a,imon,iyr))
           call static_ave_end(tnum_mave(idat_a,imon,iyr),tbias_mave(idat_a,imon,iyr),trmsd_mave(idat_a,imon,iyr))
        end do !idat_a
     end do    !imon

  end do !iyr 

  do idat_a=1,ndat_a
     call static_ave_end(unum_ave(idat_a),ubias_ave(idat_a),urmsd_ave(idat_a))
     call static_ave_end(vnum_ave(idat_a),vbias_ave(idat_a),vrmsd_ave(idat_a))
     call static_ave_end(tnum_ave(idat_a),tbias_ave(idat_a),trmsd_ave(idat_a))
  end do
  
  !---Write data
  call write_bin(im_bin,jm_bin,ndat_a,dx_bin,dy_bin,lon_bin,lat_bin, &
       unum_bin,ubias_bin,urmsd_bin, &
       vnum_bin,vbias_bin,vrmsd_bin, &
       tnum_bin,tbias_bin,trmsd_bin)

  call write_ave(syr,eyr,ndat_a, &
     & unum_mave,ubias_mave,urmsd_mave, &
     & vnum_mave,vbias_mave,vrmsd_mave, &
     & tnum_mave,tbias_mave,trmsd_mave, &
     & unum_yave,ubias_yave,urmsd_yave, &
     & vnum_yave,vbias_yave,vrmsd_yave, &
     & tnum_yave,tbias_yave,trmsd_yave, &
     & unum_ave,ubias_ave,urmsd_ave, &
     & vnum_ave,vbias_ave,vrmsd_ave, &
     & tnum_ave,tbias_ave,trmsd_ave)

  
  !---Deallocate
  call deallocate_info(filename_o, &
       & iyr_min_o,imon_min_o,iday_min_o, &
       & iyr_max_o,imon_max_o,iday_max_o, &
       & lon_min_o,lon_max_o,lat_min_o,lat_max_o)

  deallocate(ubias_bin,vbias_bin,tbias_bin)
  deallocate(urmsd_bin,vrmsd_bin,trmsd_bin)
  deallocate(unum_bin,vnum_bin,tnum_bin)
  
end program main

