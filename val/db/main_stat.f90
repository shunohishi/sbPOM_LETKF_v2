program main

  use setting
  use mod_julian
  use mod_bin
  use mod_stat
  use mod_io
  use mod_rmiss
  implicit none

  !---Common
  integer ijul,sjul,ejul
  integer iyr,imon,iday
  integer syr,smon,sday
  integer eyr,emon,eday
  integer idat_a,jdat_a
  integer nobs
  integer i_bin,j_bin
  
  !---Obs. space
  !Analysis
  real(kind = 8),allocatable :: hu_a(:),hv_a(:),ht_a(:)
  real(kind = 8),allocatable :: husprd_a(:),hvsprd_a(:),htsprd_a(:)
  
  !DB
  real(kind = 8),allocatable :: lon_o(:),lat_o(:)
  real(kind = 8),allocatable :: u_o(:),v_o(:),t_o(:)
  
  !---Stat
  !Bin (Info)
  integer im_bin,jm_bin
  real(kind = 8),allocatable :: lon_bin(:),lat_bin(:)

  !Bin (All)
  integer,allocatable :: unum_bin(:,:,:),vnum_bin(:,:,:),tnum_bin(:,:,:)
  integer,allocatable :: ubias_dof_bin(:,:,:,:),vbias_dof_bin(:,:,:,:),tbias_dof_bin(:,:,:,:)
  integer,allocatable :: urmsd_dof_bin(:,:,:,:),vrmsd_dof_bin(:,:,:,:),trmsd_dof_bin(:,:,:,:)
  
  real(kind = 8),allocatable :: ubias_bin(:,:,:),vbias_bin(:,:,:),tbias_bin(:,:,:)
  real(kind = 8),allocatable :: ubias_tcrit_bin(:,:,:,:),vbias_tcrit_bin(:,:,:,:),tbias_tcrit_bin(:,:,:,:)
  real(kind = 8),allocatable :: ubias_tval_bin(:,:,:,:),vbias_tval_bin(:,:,:,:),tbias_tval_bin(:,:,:,:)

  real(kind = 8),allocatable :: urmsd_bin(:,:,:),vrmsd_bin(:,:,:),trmsd_bin(:,:,:)
  real(kind = 8),allocatable :: urmsd_tcrit_bin(:,:,:,:),vrmsd_tcrit_bin(:,:,:,:),trmsd_tcrit_bin(:,:,:,:)
  real(kind = 8),allocatable :: urmsd_tval_bin(:,:,:,:),vrmsd_tval_bin(:,:,:,:),trmsd_tval_bin(:,:,:,:)

  real(kind = 8),allocatable :: usprd_bin(:,:,:),vsprd_bin(:,:,:),tsprd_bin(:,:,:)

  !Bin (Monthly)
  integer,allocatable :: unum_bin_mave(:,:,:,:,:),vnum_bin_mave(:,:,:,:,:),tnum_bin_mave(:,:,:,:,:)
  real(kind = 8),allocatable :: ubias_bin_mave(:,:,:,:,:),vbias_bin_mave(:,:,:,:,:),tbias_bin_mave(:,:,:,:,:)
  real(kind = 8),allocatable :: urmsd_bin_mave(:,:,:,:,:),vrmsd_bin_mave(:,:,:,:,:),trmsd_bin_mave(:,:,:,:,:)
  real(kind = 8),allocatable :: usprd_bin_mave(:,:,:,:,:),vsprd_bin_mave(:,:,:,:,:),tsprd_bin_mave(:,:,:,:,:) 
  
  !Monthly
  integer,allocatable :: unum_mave(:,:,:),vnum_mave(:,:,:),tnum_mave(:,:,:)
  real(kind = 8),allocatable :: ubias_mave(:,:,:),vbias_mave(:,:,:),tbias_mave(:,:,:),abias_mave(:,:,:)
  real(kind = 8),allocatable :: urmsd_mave(:,:,:),vrmsd_mave(:,:,:),trmsd_mave(:,:,:)
  real(kind = 8),allocatable :: usprd_mave(:,:,:),vsprd_mave(:,:,:),tsprd_mave(:,:,:)
  
  !Yearly
  integer,allocatable :: unum_yave(:,:),vnum_yave(:,:),tnum_yave(:,:)

  real(kind = 8),allocatable :: ubias_yave(:,:),vbias_yave(:,:),tbias_yave(:,:)
  real(kind = 8),allocatable :: urmsd_yave(:,:),vrmsd_yave(:,:),trmsd_yave(:,:)
  real(kind = 8),allocatable :: usprd_yave(:,:),vsprd_yave(:,:),tsprd_yave(:,:)

  !ALL
  integer unum_ave(ndat_a),vnum_ave(ndat_a),tnum_ave(ndat_a)
  integer ubias_dof_ave(ndat_a,ndat_a),vbias_dof_ave(ndat_a,ndat_a),tbias_dof_ave(ndat_a,ndat_a)
  integer urmsd_dof_ave(ndat_a,ndat_a),vrmsd_dof_ave(ndat_a,ndat_a),trmsd_dof_ave(ndat_a,ndat_a)

  real(kind = 8) ubias_ave(ndat_a),vbias_ave(ndat_a),tbias_ave(ndat_a)
  real(kind = 8) ubias_tcrit_ave(ndat_a,ndat_a),vbias_tcrit_ave(ndat_a,ndat_a),tbias_tcrit_ave(ndat_a,ndat_a)
  real(kind = 8) ubias_tval_ave(ndat_a,ndat_a),vbias_tval_ave(ndat_a,ndat_a),tbias_tval_ave(ndat_a,ndat_a)

  real(kind = 8) urmsd_ave(ndat_a),vrmsd_ave(ndat_a),trmsd_ave(ndat_a)
  real(kind = 8) urmsd_tcrit_ave(ndat_a,ndat_a),vrmsd_tcrit_ave(ndat_a,ndat_a),trmsd_tcrit_ave(ndat_a,ndat_a)
  real(kind = 8) urmsd_tval_ave(ndat_a,ndat_a),vrmsd_tval_ave(ndat_a,ndat_a),trmsd_tval_ave(ndat_a,ndat_a)

  real(kind = 8) usprd_ave(ndat_a),vsprd_ave(ndat_a),tsprd_ave(ndat_a)

  write(*,*) "### START: Validation vs. drifter buoy ###"

  !---Read start and end dates
  call read_argument(syr,smon,sday,eyr,emon,eday)
  write(*,*) "Start date:",syr,smon,sday
  write(*,*) "End date:",eyr,emon,eday

  !---Make lon-lat bin
  call make_bin(slon_bin,elon_bin,slat_bin,elat_bin,dx_bin,dy_bin, &
       & im_bin,jm_bin,lon_bin,lat_bin)
  
  !---Allocate
  !Bin
  !All
  allocate(unum_bin(im_bin,jm_bin,ndat_a),vnum_bin(im_bin,jm_bin,ndat_a),tnum_bin(im_bin,jm_bin,ndat_a))

  allocate(ubias_bin(im_bin,jm_bin,ndat_a),vbias_bin(im_bin,jm_bin,ndat_a),tbias_bin(im_bin,jm_bin,ndat_a))
  allocate(ubias_dof_bin(im_bin,jm_bin,ndat_a,ndat_a),vbias_dof_bin(im_bin,jm_bin,ndat_a,ndat_a),tbias_dof_bin(im_bin,jm_bin,ndat_a,ndat_a))  
  allocate(ubias_tcrit_bin(im_bin,jm_bin,ndat_a,ndat_a),vbias_tcrit_bin(im_bin,jm_bin,ndat_a,ndat_a),tbias_tcrit_bin(im_bin,jm_bin,ndat_a,ndat_a))
  allocate(ubias_tval_bin(im_bin,jm_bin,ndat_a,ndat_a),vbias_tval_bin(im_bin,jm_bin,ndat_a,ndat_a),tbias_tval_bin(im_bin,jm_bin,ndat_a,ndat_a))

  allocate(urmsd_bin(im_bin,jm_bin,ndat_a),vrmsd_bin(im_bin,jm_bin,ndat_a),trmsd_bin(im_bin,jm_bin,ndat_a))
  allocate(urmsd_dof_bin(im_bin,jm_bin,ndat_a,ndat_a),vrmsd_dof_bin(im_bin,jm_bin,ndat_a,ndat_a),trmsd_dof_bin(im_bin,jm_bin,ndat_a,ndat_a))  
  allocate(urmsd_tcrit_bin(im_bin,jm_bin,ndat_a,ndat_a),vrmsd_tcrit_bin(im_bin,jm_bin,ndat_a,ndat_a),trmsd_tcrit_bin(im_bin,jm_bin,ndat_a,ndat_a))
  allocate(urmsd_tval_bin(im_bin,jm_bin,ndat_a,ndat_a),vrmsd_tval_bin(im_bin,jm_bin,ndat_a,ndat_a),trmsd_tval_bin(im_bin,jm_bin,ndat_a,ndat_a))

  allocate(usprd_bin(im_bin,jm_bin,ndat_a),vsprd_bin(im_bin,jm_bin,ndat_a),tsprd_bin(im_bin,jm_bin,ndat_a))

  !Month
  allocate(unum_bin_mave(im_bin,jm_bin,ndat_a,12,syr:eyr),vnum_bin_mave(im_bin,jm_bin,ndat_a,12,syr:eyr),tnum_bin_mave(im_bin,jm_bin,ndat_a,12,syr:eyr))
  allocate(ubias_bin_mave(im_bin,jm_bin,ndat_a,12,syr:eyr),vbias_bin_mave(im_bin,jm_bin,ndat_a,12,syr:eyr),tbias_bin_mave(im_bin,jm_bin,ndat_a,12,syr:eyr))
  allocate(urmsd_bin_mave(im_bin,jm_bin,ndat_a,12,syr:eyr),vrmsd_bin_mave(im_bin,jm_bin,ndat_a,12,syr:eyr),trmsd_bin_mave(im_bin,jm_bin,ndat_a,12,syr:eyr))
  allocate(usprd_bin_mave(im_bin,jm_bin,ndat_a,12,syr:eyr),vsprd_bin_mave(im_bin,jm_bin,ndat_a,12,syr:eyr),tsprd_bin_mave(im_bin,jm_bin,ndat_a,12,syr:eyr))
  
  !Month
  allocate(unum_mave(ndat_a,12,syr:eyr),vnum_mave(ndat_a,12,syr:eyr),tnum_mave(ndat_a,12,syr:eyr))
  allocate(ubias_mave(ndat_a,12,syr:eyr),vbias_mave(ndat_a,12,syr:eyr),tbias_mave(ndat_a,12,syr:eyr),abias_mave(ndat_a,12,syr:eyr))
  allocate(urmsd_mave(ndat_a,12,syr:eyr),vrmsd_mave(ndat_a,12,syr:eyr),trmsd_mave(ndat_a,12,syr:eyr))
  allocate(usprd_mave(ndat_a,12,syr:eyr),vsprd_mave(ndat_a,12,syr:eyr),tsprd_mave(ndat_a,12,syr:eyr))

  !Year
  allocate(unum_yave(ndat_a,syr:eyr),vnum_yave(ndat_a,syr:eyr),tnum_yave(ndat_a,syr:eyr))
  allocate(ubias_yave(ndat_a,syr:eyr),vbias_yave(ndat_a,syr:eyr),tbias_yave(ndat_a,syr:eyr))
  allocate(urmsd_yave(ndat_a,syr:eyr),vrmsd_yave(ndat_a,syr:eyr),trmsd_yave(ndat_a,syr:eyr))
  allocate(usprd_yave(ndat_a,syr:eyr),vsprd_yave(ndat_a,syr:eyr),tsprd_yave(ndat_a,syr:eyr))

  !---Initialization
  do idat_a=1,ndat_a
     call stat_bin_ini(im_bin,jm_bin,unum_bin(:,:,idat_a),ubias_bin(:,:,idat_a),urmsd_bin(:,:,idat_a),usprd_bin(:,:,idat_a))
     call stat_bin_ini(im_bin,jm_bin,vnum_bin(:,:,idat_a),vbias_bin(:,:,idat_a),vrmsd_bin(:,:,idat_a),vsprd_bin(:,:,idat_a))
     call stat_bin_ini(im_bin,jm_bin,tnum_bin(:,:,idat_a),tbias_bin(:,:,idat_a),trmsd_bin(:,:,idat_a),tsprd_bin(:,:,idat_a))     
  end do!idat_a

  do iyr=syr,eyr
     
     do idat_a=1,ndat_a
        call stat_ave_ini(unum_yave(idat_a,iyr),ubias_yave(idat_a,iyr),urmsd_yave(idat_a,iyr),usprd_yave(idat_a,iyr))
        call stat_ave_ini(vnum_yave(idat_a,iyr),vbias_yave(idat_a,iyr),vrmsd_yave(idat_a,iyr),vsprd_yave(idat_a,iyr))
        call stat_ave_ini(tnum_yave(idat_a,iyr),tbias_yave(idat_a,iyr),trmsd_yave(idat_a,iyr),tsprd_yave(idat_a,iyr))
     end do !idat_a
  
     do imon=1,12
        do idat_a=1,ndat_a

           call stat_bin_ini(im_bin,jm_bin, &
                & unum_bin_mave(:,:,idat_a,imon,iyr),ubias_bin_mave(:,:,idat_a,imon,iyr),urmsd_bin_mave(:,:,idat_a,imon,iyr),usprd_bin_mave(:,:,idat_a,imon,iyr))
           call stat_bin_ini(im_bin,jm_bin, &
                & vnum_bin_mave(:,:,idat_a,imon,iyr),vbias_bin_mave(:,:,idat_a,imon,iyr),vrmsd_bin_mave(:,:,idat_a,imon,iyr),vsprd_bin_mave(:,:,idat_a,imon,iyr))
           call stat_bin_ini(im_bin,jm_bin, &
                & tnum_bin_mave(:,:,idat_a,imon,iyr),tbias_bin_mave(:,:,idat_a,imon,iyr),trmsd_bin_mave(:,:,idat_a,imon,iyr),tsprd_bin_mave(:,:,idat_a,imon,iyr))     
           
           call stat_ave_ini(unum_mave(idat_a,imon,iyr),ubias_mave(idat_a,imon,iyr),urmsd_mave(idat_a,imon,iyr),usprd_mave(idat_a,imon,iyr))
           call stat_ave_ini(vnum_mave(idat_a,imon,iyr),vbias_mave(idat_a,imon,iyr),vrmsd_mave(idat_a,imon,iyr),vsprd_mave(idat_a,imon,iyr))
           call stat_ave_ini(tnum_mave(idat_a,imon,iyr),tbias_mave(idat_a,imon,iyr),trmsd_mave(idat_a,imon,iyr),tsprd_mave(idat_a,imon,iyr))
        end do !idat_a
     end do    !imon
     
  end do       !iyr

  do idat_a=1,ndat_a
     call stat_ave_ini(unum_ave(idat_a),ubias_ave(idat_a),urmsd_ave(idat_a),usprd_ave(idat_a))
     call stat_ave_ini(vnum_ave(idat_a),vbias_ave(idat_a),vrmsd_ave(idat_a),vsprd_ave(idat_a))
     call stat_ave_ini(tnum_ave(idat_a),tbias_ave(idat_a),trmsd_ave(idat_a),tsprd_ave(idat_a))
  end do

  !--- Main Loop ---!
  do idat_a=1,ndat_a
          
     !---Validation loop
     call ymd_julian(syr,smon,sday,sjul)
     call ymd_julian(eyr,emon,eday,ejul)  
     do ijul=sjul,ejul

        !---Initialization
        call julian_ymd(ijul,iyr,imon,iday)
        write(*,*) idat_a,iyr,imon,iday

        !---Read
        write(*,*) "Read data in observation space"
        call read_obs(idat_a,iyr,imon,iday,nobs,lon_o,lat_o,ht_a,hu_a,hv_a,htsprd_a,husprd_a,hvsprd_a,t_o,u_o,v_o)
        
        if(nobs == 0) cycle

        !---Stat
        write(*,*) "Valdation"
        !Bin(All)
        call stat_bin_add(nobs,lon_o,lat_o,hu_a,husprd_a,u_o, &
             & im_bin,jm_bin,lon_bin,lat_bin, &
             & unum_bin(:,:,idat_a),ubias_bin(:,:,idat_a),urmsd_bin(:,:,idat_a),usprd_bin(:,:,idat_a))
        call stat_bin_add(nobs,lon_o,lat_o,hv_a,hvsprd_a,v_o, &
             & im_bin,jm_bin,lon_bin,lat_bin, &
             & vnum_bin(:,:,idat_a),vbias_bin(:,:,idat_a),vrmsd_bin(:,:,idat_a),vsprd_bin(:,:,idat_a))
        call stat_bin_add(nobs,lon_o,lat_o,ht_a,htsprd_a,t_o, &
             & im_bin,jm_bin,lon_bin,lat_bin, &
             & tnum_bin(:,:,idat_a),tbias_bin(:,:,idat_a),trmsd_bin(:,:,idat_a),tsprd_bin(:,:,idat_a))

        !Bin(Monthly)
        call stat_bin_add(nobs,lon_o,lat_o,hu_a,husprd_a,u_o, &
             & im_bin,jm_bin,lon_bin,lat_bin, &
             & unum_bin_mave(:,:,idat_a,imon,iyr),ubias_bin_mave(:,:,idat_a,imon,iyr),urmsd_bin_mave(:,:,idat_a,imon,iyr),usprd_bin_mave(:,:,idat_a,imon,iyr))
        call stat_bin_add(nobs,lon_o,lat_o,hv_a,hvsprd_a,v_o, &
             & im_bin,jm_bin,lon_bin,lat_bin, &
             & vnum_bin_mave(:,:,idat_a,imon,iyr),vbias_bin_mave(:,:,idat_a,imon,iyr),vrmsd_bin_mave(:,:,idat_a,imon,iyr),vsprd_bin_mave(:,:,idat_a,imon,iyr))
        call stat_bin_add(nobs,lon_o,lat_o,ht_a,htsprd_a,t_o, &
             & im_bin,jm_bin,lon_bin,lat_bin, &
             & tnum_bin_mave(:,:,idat_a,imon,iyr),tbias_bin_mave(:,:,idat_a,imon,iyr),trmsd_bin_mave(:,:,idat_a,imon,iyr),tsprd_bin_mave(:,:,idat_a,imon,iyr))
        
        !Monthly
        call stat_ave_add(nobs,hu_a,husprd_a,u_o, &
             & unum_mave(idat_a,imon,iyr),ubias_mave(idat_a,imon,iyr),urmsd_mave(idat_a,imon,iyr),usprd_mave(idat_a,imon,iyr))
        call stat_ave_add(nobs,hv_a,hvsprd_a,v_o, &
             & vnum_mave(idat_a,imon,iyr),vbias_mave(idat_a,imon,iyr),vrmsd_mave(idat_a,imon,iyr),vsprd_mave(idat_a,imon,iyr))
        call stat_ave_add(nobs,ht_a,htsprd_a,t_o, &
             & tnum_mave(idat_a,imon,iyr),tbias_mave(idat_a,imon,iyr),trmsd_mave(idat_a,imon,iyr),tsprd_mave(idat_a,imon,iyr))

        !Yearly
        call stat_ave_add(nobs,hu_a,husprd_a,u_o, &
             & unum_yave(idat_a,iyr),ubias_yave(idat_a,iyr),urmsd_yave(idat_a,iyr),usprd_yave(idat_a,iyr))
        call stat_ave_add(nobs,hv_a,hvsprd_a,v_o, &
             & vnum_yave(idat_a,iyr),vbias_yave(idat_a,iyr),vrmsd_yave(idat_a,iyr),vsprd_yave(idat_a,iyr))
        call stat_ave_add(nobs,ht_a,htsprd_a,t_o, &
             & tnum_yave(idat_a,iyr),tbias_yave(idat_a,iyr),trmsd_yave(idat_a,iyr),tsprd_yave(idat_a,iyr))

        !ALL
        call stat_ave_add(nobs,hu_a,husprd_a,u_o, &
             & unum_ave(idat_a),ubias_ave(idat_a),urmsd_ave(idat_a),usprd_ave(idat_a))
        call stat_ave_add(nobs,hv_a,hvsprd_a,v_o, &
             & vnum_ave(idat_a),vbias_ave(idat_a),vrmsd_ave(idat_a),vsprd_ave(idat_a))
        call stat_ave_add(nobs,ht_a,htsprd_a,t_o, &
             & tnum_ave(idat_a),tbias_ave(idat_a),trmsd_ave(idat_a),tsprd_ave(idat_a))

        !---Deallocate
        call deallocate_obs(lon_o,lat_o, &
             & ht_a,hu_a,hv_a,htsprd_a,husprd_a,hvsprd_a,t_o,u_o,v_o)
        
     end do       !ijul     
  end do          !idat_a

  !---Post process
  write(*,*) "Post process"
  !Bin
  do idat_a=1,ndat_a
     call stat_bin_end(im_bin,jm_bin, &
          & unum_bin(:,:,idat_a),ubias_bin(:,:,idat_a),urmsd_bin(:,:,idat_a),usprd_bin(:,:,idat_a))
     call stat_bin_end(im_bin,jm_bin, &
          & vnum_bin(:,:,idat_a),vbias_bin(:,:,idat_a),vrmsd_bin(:,:,idat_a),vsprd_bin(:,:,idat_a))
     call stat_bin_end(im_bin,jm_bin, &
          & tnum_bin(:,:,idat_a),tbias_bin(:,:,idat_a),trmsd_bin(:,:,idat_a),tsprd_bin(:,:,idat_a))
  end do !idat_a

  !Ave
  do iyr=syr,eyr

     do idat_a=1,ndat_a
        call stat_ave_end(unum_yave(idat_a,iyr),ubias_yave(idat_a,iyr),urmsd_yave(idat_a,iyr),usprd_yave(idat_a,iyr))
        call stat_ave_end(vnum_yave(idat_a,iyr),vbias_yave(idat_a,iyr),vrmsd_yave(idat_a,iyr),vsprd_yave(idat_a,iyr))
        call stat_ave_end(tnum_yave(idat_a,iyr),tbias_yave(idat_a,iyr),trmsd_yave(idat_a,iyr),tsprd_yave(idat_a,iyr))
     end do !idat_a
     
     do imon=1,12
        do idat_a=1,ndat_a

           call stat_bin_end(im_bin,jm_bin, &
                & unum_bin_mave(:,:,idat_a,imon,iyr),ubias_bin_mave(:,:,idat_a,imon,iyr),urmsd_bin_mave(:,:,idat_a,imon,iyr),usprd_bin_mave(:,:,idat_a,imon,iyr))
           call stat_bin_end(im_bin,jm_bin, &
                & vnum_bin_mave(:,:,idat_a,imon,iyr),vbias_bin_mave(:,:,idat_a,imon,iyr),vrmsd_bin_mave(:,:,idat_a,imon,iyr),vsprd_bin_mave(:,:,idat_a,imon,iyr))
           call stat_bin_end(im_bin,jm_bin, &
                & tnum_bin_mave(:,:,idat_a,imon,iyr),tbias_bin_mave(:,:,idat_a,imon,iyr),trmsd_bin_mave(:,:,idat_a,imon,iyr),tsprd_bin_mave(:,:,idat_a,imon,iyr))

           
           call stat_ave_end(unum_mave(idat_a,imon,iyr),ubias_mave(idat_a,imon,iyr),urmsd_mave(idat_a,imon,iyr),usprd_mave(idat_a,imon,iyr))
           call stat_ave_end(vnum_mave(idat_a,imon,iyr),vbias_mave(idat_a,imon,iyr),vrmsd_mave(idat_a,imon,iyr),vsprd_mave(idat_a,imon,iyr))
           call stat_ave_end(tnum_mave(idat_a,imon,iyr),tbias_mave(idat_a,imon,iyr),trmsd_mave(idat_a,imon,iyr),tsprd_mave(idat_a,imon,iyr))
        end do !idat_a
     end do    !imon

  end do !iyr 

  do idat_a=1,ndat_a
     call stat_ave_end(unum_ave(idat_a),ubias_ave(idat_a),urmsd_ave(idat_a),usprd_ave(idat_a))
     call stat_ave_end(vnum_ave(idat_a),vbias_ave(idat_a),vrmsd_ave(idat_a),vsprd_ave(idat_a))
     call stat_ave_end(tnum_ave(idat_a),tbias_ave(idat_a),trmsd_ave(idat_a),tsprd_ave(idat_a))
  end do

  !---Paired t-test using monthly RMSD difference
  do j_bin=1,jm_bin
     do i_bin=1,im_bin
        call substitute_rmiss(ndat_a,ubias_dof_bin(i_bin,j_bin,:,:),ubias_tcrit_bin(i_bin,j_bin,:,:),ubias_tval_bin(i_bin,j_bin,:,:))
        call substitute_rmiss(ndat_a,vbias_dof_bin(i_bin,j_bin,:,:),vbias_tcrit_bin(i_bin,j_bin,:,:),vbias_tval_bin(i_bin,j_bin,:,:))
        call substitute_rmiss(ndat_a,tbias_dof_bin(i_bin,j_bin,:,:),tbias_tcrit_bin(i_bin,j_bin,:,:),tbias_tval_bin(i_bin,j_bin,:,:))

        call substitute_rmiss(ndat_a,urmsd_dof_bin(i_bin,j_bin,:,:),urmsd_tcrit_bin(i_bin,j_bin,:,:),urmsd_tval_bin(i_bin,j_bin,:,:))
        call substitute_rmiss(ndat_a,vrmsd_dof_bin(i_bin,j_bin,:,:),vrmsd_tcrit_bin(i_bin,j_bin,:,:),vrmsd_tval_bin(i_bin,j_bin,:,:))
        call substitute_rmiss(ndat_a,trmsd_dof_bin(i_bin,j_bin,:,:),trmsd_tcrit_bin(i_bin,j_bin,:,:),trmsd_tval_bin(i_bin,j_bin,:,:))
     end do
  end do

  call substitute_rmiss(ndat_a,ubias_dof_ave,ubias_tcrit_ave,ubias_tval_ave)
  call substitute_rmiss(ndat_a,vbias_dof_ave,vbias_tcrit_ave,vbias_tval_ave)
  call substitute_rmiss(ndat_a,tbias_dof_ave,tbias_tcrit_ave,tbias_tval_ave)

  call substitute_rmiss(ndat_a,urmsd_dof_ave,urmsd_tcrit_ave,urmsd_tval_ave)
  call substitute_rmiss(ndat_a,vrmsd_dof_ave,vrmsd_tcrit_ave,vrmsd_tval_ave)
  call substitute_rmiss(ndat_a,trmsd_dof_ave,trmsd_tcrit_ave,trmsd_tval_ave)
  
  do jdat_a=1,ndat_a
     do idat_a=1,ndat_a

        if(idat_a == jdat_a) cycle
        
        do j_bin=1,jm_bin
           do i_bin=1,im_bin

              call convert_absolute_bias(ndat_a*12*(eyr-syr+1),ubias_bin_mave(i_bin,j_bin,:,:,syr:eyr),abias_mave(:,:,syr:eyr))
              call paired_t_test((eyr-syr+1)*12,abias_mave(idat_a,:,syr:eyr),abias_mave(jdat_a,:,syr:eyr), &
                   & ubias_dof_bin(i_bin,j_bin,idat_a,jdat_a),ubias_tcrit_bin(i_bin,j_bin,idat_a,jdat_a),ubias_tval_bin(i_bin,j_bin,idat_a,jdat_a))
              call convert_absolute_bias(ndat_a*12*(eyr-syr+1),vbias_bin_mave(i_bin,j_bin,:,:,syr:eyr),abias_mave(:,:,syr:eyr))
              call paired_t_test((eyr-syr+1)*12,abias_mave(idat_a,:,syr:eyr),abias_mave(jdat_a,:,syr:eyr), &
                   & vbias_dof_bin(i_bin,j_bin,idat_a,jdat_a),vbias_tcrit_bin(i_bin,j_bin,idat_a,jdat_a),vbias_tval_bin(i_bin,j_bin,idat_a,jdat_a))
              call convert_absolute_bias(ndat_a*12*(eyr-syr+1),tbias_bin_mave(i_bin,j_bin,:,:,syr:eyr),abias_mave(:,:,syr:eyr))
              call paired_t_test((eyr-syr+1)*12,abias_mave(idat_a,:,syr:eyr),abias_mave(jdat_a,:,syr:eyr), &
                   & tbias_dof_bin(i_bin,j_bin,idat_a,jdat_a),tbias_tcrit_bin(i_bin,j_bin,idat_a,jdat_a),tbias_tval_bin(i_bin,j_bin,idat_a,jdat_a))
                            
              call paired_t_test((eyr-syr+1)*12,urmsd_bin_mave(i_bin,j_bin,idat_a,:,syr:eyr),urmsd_bin_mave(i_bin,j_bin,jdat_a,:,syr:eyr), &
                   & urmsd_dof_bin(i_bin,j_bin,idat_a,jdat_a),urmsd_tcrit_bin(i_bin,j_bin,idat_a,jdat_a),urmsd_tval_bin(i_bin,j_bin,idat_a,jdat_a))
              call paired_t_test((eyr-syr+1)*12,vrmsd_bin_mave(i_bin,j_bin,idat_a,:,syr:eyr),vrmsd_bin_mave(i_bin,j_bin,jdat_a,:,syr:eyr), &
                   & vrmsd_dof_bin(i_bin,j_bin,idat_a,jdat_a),vrmsd_tcrit_bin(i_bin,j_bin,idat_a,jdat_a),vrmsd_tval_bin(i_bin,j_bin,idat_a,jdat_a))
              call paired_t_test((eyr-syr+1)*12,trmsd_bin_mave(i_bin,j_bin,idat_a,:,syr:eyr),trmsd_bin_mave(i_bin,j_bin,jdat_a,:,syr:eyr), &
                   & trmsd_dof_bin(i_bin,j_bin,idat_a,jdat_a),trmsd_tcrit_bin(i_bin,j_bin,idat_a,jdat_a),trmsd_tval_bin(i_bin,j_bin,idat_a,jdat_a))

           end do
        end do

        call convert_absolute_bias(ndat_a*12*(eyr-syr+1),ubias_mave(:,:,syr:eyr),abias_mave(:,:,syr:eyr))
        call paired_t_test((eyr-syr+1)*12,abias_mave(idat_a,:,syr:eyr),abias_mave(jdat_a,:,syr:eyr), &
             & ubias_dof_ave(idat_a,jdat_a),ubias_tcrit_ave(idat_a,jdat_a),ubias_tval_ave(idat_a,jdat_a))
        call convert_absolute_bias(ndat_a*12*(eyr-syr+1),vbias_mave(:,:,syr:eyr),abias_mave(:,:,syr:eyr))
        call paired_t_test((eyr-syr+1)*12,abias_mave(idat_a,:,syr:eyr),abias_mave(jdat_a,:,syr:eyr), &
             & vbias_dof_ave(idat_a,jdat_a),vbias_tcrit_ave(idat_a,jdat_a),vbias_tval_ave(idat_a,jdat_a))
        call convert_absolute_bias(ndat_a*12*(eyr-syr+1),tbias_mave(:,:,syr:eyr),abias_mave(:,:,syr:eyr))
        call paired_t_test((eyr-syr+1)*12,abias_mave(idat_a,:,syr:eyr),abias_mave(jdat_a,:,syr:eyr), &
             & tbias_dof_ave(idat_a,jdat_a),tbias_tcrit_ave(idat_a,jdat_a),tbias_tval_ave(idat_a,jdat_a))
        
        call paired_t_test((eyr-syr+1)*12,urmsd_mave(idat_a,:,syr:eyr),urmsd_mave(jdat_a,:,syr:eyr), &
             & urmsd_dof_ave(idat_a,jdat_a),urmsd_tcrit_ave(idat_a,jdat_a),urmsd_tval_ave(idat_a,jdat_a))
        call paired_t_test((eyr-syr+1)*12,vrmsd_mave(idat_a,:,syr:eyr),vrmsd_mave(jdat_a,:,syr:eyr), &
             & vrmsd_dof_ave(idat_a,jdat_a),vrmsd_tcrit_ave(idat_a,jdat_a),vrmsd_tval_ave(idat_a,jdat_a))
        call paired_t_test((eyr-syr+1)*12,trmsd_mave(idat_a,:,syr:eyr),trmsd_mave(jdat_a,:,syr:eyr), &
             & trmsd_dof_ave(idat_a,jdat_a),trmsd_tcrit_ave(idat_a,jdat_a),trmsd_tval_ave(idat_a,jdat_a))
        
     end do
  end do
  
  !---Write data
  write(*,*) "Write data"
  call write_bin(im_bin,jm_bin,ndat_a,dx_bin,dy_bin,lon_bin,lat_bin, &
       & unum_bin,ubias_bin,urmsd_bin,usprd_bin,       &
       & ubias_dof_bin,ubias_tcrit_bin,ubias_tval_bin, &
       & urmsd_dof_bin,urmsd_tcrit_bin,urmsd_tval_bin, &
       & vnum_bin,vbias_bin,vrmsd_bin,vsprd_bin,       &
       & vbias_dof_bin,vbias_tcrit_bin,vbias_tval_bin, &
       & vrmsd_dof_bin,vrmsd_tcrit_bin,vrmsd_tval_bin, &
       & tnum_bin,tbias_bin,trmsd_bin,tsprd_bin,       &
       & tbias_dof_bin,tbias_tcrit_bin,tbias_tval_bin, &
       & trmsd_dof_bin,trmsd_tcrit_bin,trmsd_tval_bin)

  call write_ave(syr,eyr,ndat_a, &
     & unum_mave,ubias_mave,urmsd_mave,usprd_mave, &
     & vnum_mave,vbias_mave,vrmsd_mave,vsprd_mave, &
     & tnum_mave,tbias_mave,trmsd_mave,tsprd_mave, &
     & unum_yave,ubias_yave,urmsd_yave,usprd_yave, &
     & vnum_yave,vbias_yave,vrmsd_yave,vsprd_yave, &
     & tnum_yave,tbias_yave,trmsd_yave,tsprd_yave, &
     & unum_ave,ubias_ave,urmsd_ave,usprd_ave,       &
     & ubias_dof_ave,ubias_tcrit_ave,ubias_tval_ave, &
     & urmsd_dof_ave,urmsd_tcrit_ave,urmsd_tval_ave, &
     & vnum_ave,vbias_ave,vrmsd_ave,vsprd_ave,       &
     & vbias_dof_ave,vbias_tcrit_ave,vbias_tval_ave, &
     & vrmsd_dof_ave,vrmsd_tcrit_ave,vrmsd_tval_ave, &
     & tnum_ave,tbias_ave,trmsd_ave,tsprd_ave,       &
     & tbias_dof_ave,tbias_tcrit_ave,tbias_tval_ave, &
     & trmsd_dof_ave,trmsd_tcrit_ave,trmsd_tval_ave)
  
  !---Deallocate
  deallocate(unum_bin,vnum_bin,tnum_bin)
  
  deallocate(ubias_bin,vbias_bin,tbias_bin)
  deallocate(ubias_dof_bin,vbias_dof_bin,tbias_dof_bin)
  deallocate(ubias_tcrit_bin,vbias_tcrit_bin,tbias_tcrit_bin)
  deallocate(ubias_tval_bin,vbias_tval_bin,tbias_tval_bin)
  
  deallocate(urmsd_bin,vrmsd_bin,trmsd_bin)
  deallocate(urmsd_dof_bin,vrmsd_dof_bin,trmsd_dof_bin)
  deallocate(urmsd_tcrit_bin,vrmsd_tcrit_bin,trmsd_tcrit_bin)
  deallocate(urmsd_tval_bin,vrmsd_tval_bin,trmsd_tval_bin)

  deallocate(usprd_bin,vsprd_bin,tsprd_bin)

  deallocate(unum_mave,vnum_mave,tnum_mave)
  deallocate(ubias_mave,vbias_mave,tbias_mave,abias_mave)
  deallocate(urmsd_mave,vrmsd_mave,trmsd_mave)
  deallocate(usprd_mave,vsprd_mave,tsprd_mave)
  
  deallocate(unum_yave,vnum_yave,tnum_yave)
  deallocate(ubias_yave,vbias_yave,tbias_yave)
  deallocate(urmsd_yave,vrmsd_yave,trmsd_yave)
  deallocate(usprd_yave,vsprd_yave,tsprd_yave)

  write(*,*) "### SUCCESS: Validation vs. drifter buoy ###"
  
end program main
