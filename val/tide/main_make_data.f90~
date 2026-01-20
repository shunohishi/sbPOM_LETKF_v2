program main

  use setting  
  use mod_rmiss
  use mod_julian
  use mod_gridinfo, im_lora => im, jm_lora => jm, km_lora => km
  use mod_read_glorys025, im_g025 => im, jm_g025 => jm, km_g025 => km
  use mod_read_tide
  use mod_io
  implicit none

  !---Common
  integer i,j,k
  integer iyr,imon,iday
  integer syr,smon,sday
  integer eyr,emon,eday
  integer ijul,sjul,ejul

  real(kind = 8) dist_save(nst)
  
  !---Analysis data
  integer idat_a
  integer im_a,jm_a,km_a
  integer idx,idy
  
  real(kind = 8),allocatable :: lon_a(:),lat_a(:)
  real(kind = 8),allocatable :: mask_a(:,:)
  real(kind = 8),allocatable :: dat_a(:,:),sprd_a(:,:)

  !Obs. space
  real(kind = 8) lon_a_save(nst),lat_a_save(nst)
  real(kind = 8) hdat_a_save(nst),hsprd_a_save(nst)
  
  !---Observation: Sea level from tide gauge
  integer ist
  integer idt
  
  !Input
  integer ntime_o_tmp(nst)
  integer,allocatable :: ijul_o_tmp(:)
  real(kind = 8),allocatable :: lon_o_tmp(:),lat_o_tmp(:),dat_o_tmp(:)

  integer ntime_o
  integer,allocatable :: ijul_o(:,:)
  real(kind = 8),allocatable :: lon_o(:,:),lat_o(:,:),dat_o(:,:)
  
  !Info
  integer ijul_o_min(nst),ijul_o_max(nst)
  
  real(kind = 8) lon_o_min(nst),lon_o_max(nst)
  real(kind = 8) lat_o_min(nst),lat_o_max(nst)

  !Save
  real(kind = 8) lon_o_save(nst),lat_o_save(nst)
  real(kind = 8) dat_o_save(nst)
  
  !---Julian day
  call read_argument(syr,smon,sday,eyr,emon,eday)
  call ymd_julian(syr,smon,sday,sjul)
  call ymd_julian(eyr,emon,eday,ejul)

  !---Read tide data
  !Domain and Temporal coverage
  do ist=1,nst
     
     call read_tide(ist,ntime_o_tmp(ist),ijul_o_tmp,lon_o_tmp,lat_o_tmp,dat_o_tmp)
     
     if(ntime_o_tmp(ist) == 0)then
        ijul_o_min(ist)=int(rmiss)
        ijul_o_max(ist)=int(rmiss)
        lon_o_min(ist)=rmiss
        lon_o_max(ist)=rmiss
        lat_o_min(ist)=rmiss
        lat_o_max(ist)=rmiss
     else
        ijul_o_min(ist)=minval(ijul_o_tmp)
        ijul_o_max(ist)=maxval(ijul_o_tmp)
        lon_o_min(ist)=minval(lon_o_tmp)
        lon_o_max(ist)=maxval(lon_o_tmp)
        lat_o_min(ist)=minval(lat_o_tmp)
        lat_o_max(ist)=maxval(lat_o_tmp)
        call end_read_tide(ijul_o_tmp,lon_o_tmp,lat_o_tmp,dat_o_tmp)
     end if
     
  end do

  !Allocate
  ntime_o=maxval(ntime_o_tmp)
  allocate(ijul_o(ntime_o,nst))
  allocate(lon_o(ntime_o,nst),lat_o(ntime_o,nst),dat_o(ntime_o,nst))

  !Initialize
  ijul_o(:,:)=int(rmiss)
  lon_o(:,:)=rmiss
  lat_o(:,:)=rmiss
  dat_o(:,:)=rmiss
  
  !Read data
  do ist=1,nst
     
     call read_tide(ist,ntime_o_tmp(ist),ijul_o_tmp,lon_o_tmp,lat_o_tmp,dat_o_tmp)

     if(ntime_o_tmp(ist) == 0) cycle

     ijul_o(1:ntime_o_tmp(ist),ist)=ijul_o_tmp(1:ntime_o_tmp(ist))
     lon_o(1:ntime_o_tmp(ist),ist)=lon_o_tmp(1:ntime_o_tmp(ist))
     lat_o(1:ntime_o_tmp(ist),ist)=lat_o_tmp(1:ntime_o_tmp(ist))
     dat_o(1:ntime_o_tmp(ist),ist)=dat_o_tmp(1:ntime_o_tmp(ist))

     call end_read_tide(ijul_o_tmp,lon_o_tmp,lat_o_tmp,dat_o_tmp)
     
  end do
  
  !=== Main Loop ===!
  do idat_a=1,ndat_a

     !---Grid size
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

     !---Allocate
     allocate(lon_a(im_a),lat_a(jm_a))
     allocate(mask_a(im_a,jm_a))
     allocate(dat_a(im_a,jm_a),sprd_a(im_a,jm_a))
     
     !---Read DA Grid
     call read_grid(idat_a,im_a,jm_a,km_a,lon_a,lat_a,mask_a)
     
     do ijul=sjul,ejul

        !---Initialization
        lon_a_save(:)=rmiss
        lat_a_save(:)=rmiss
        hdat_a_save(:)=rmiss
        hsprd_a_save(:)=rmiss
        lon_o_save(:)=rmiss
        lat_o_save(:)=rmiss
        dat_o_save(:)=rmiss
        dist_save(:)=rmiss
        
        !---Date
        call julian_ymd(ijul,iyr,imon,iday)
        write(*,*) idat_a,iyr,imon,iday
        
        !---Read analysis sea level data
        call read_ssh_anal(idat_a,iyr,imon,iday,im_a,jm_a,km_a,mask_a,dat_a,sprd_a)
        
        !---vs. Tide data
        do ist=1,nst

           !Check data coverage
           if(ntime_o_tmp(ist) == 0) cycle
           if(ejul < ijul_o_min(ist) .or. ijul_o_max(ist) < sjul) cycle
           if(lon_o_min(ist) < lon_a(1) .or. lon_a(im_a) < lon_o_max(ist)) cycle
           if(lat_o_min(ist) < lat_a(1) .or. lat_a(jm_a) < lat_o_max(ist)) cycle

           !Date
           call get_idt(ntime_o_tmp(ist),ijul_o(1:ntime_o_tmp(ist),ist),ijul,idt)

           if(idt == 0)then
              cycle
           end if

           !ID at the closest grid point
           call get_id(im_a,jm_a,lon_a,lat_a,mask_a,lon_o(idt,ist),lat_o(idt,ist),idx,idy,dist_save(ist))
           
           !Save data
           if(idx == 0 .or. idy == 0)then
              lon_a_save(ist)=rmiss
              lat_a_save(ist)=rmiss
              hdat_a_save(ist)=rmiss
              hsprd_a_save(ist)=rmiss              
           else if(dist_crit < dist_save(ist))then
              lon_a_save(ist)=lon_a(idx)
              lat_a_save(ist)=lat_a(idy)
              hdat_a_save(ist)=rmiss
              hsprd_a_save(ist)=rmiss
           else
              lon_a_save(ist)=lon_a(idx)
              lat_a_save(ist)=lat_a(idy)
              hdat_a_save(ist)=dat_a(idx,idy)
              hsprd_a_save(ist)=sprd_a(idx,idy)
           end if
           
           lon_o_save(ist)=lon_o(idt,ist)
           lat_o_save(ist)=lat_o(idt,ist)
           dat_o_save(ist)=dat_o(idt,ist)
           
        end do !ist
        
        !Write data in obs. space
        call write_hdat(idat_a,nst,ijul,lon_a_save,lat_a_save,hdat_a_save,hsprd_a_save, &
             & lon_o_save,lat_o_save,dat_o_save,dist_save)
        
     end do !ijul
     
     !---Deallocate
     deallocate(lon_a,lat_a)
     deallocate(mask_a)
     deallocate(dat_a,sprd_a)
     
  end do !idat_a

  deallocate(ijul_o)
  deallocate(lon_o,lat_o,dat_o)
  
end program main
