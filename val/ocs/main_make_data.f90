program main

  use setting
  use mod_julian
  use mod_read_ocs
  use mod_gridinfo, im_lora => im, jm_lora => jm, km_lora => km
  use mod_read_glorys025, im_g025 => im, jm_g025 => jm, km_g025 => km
  use mod_io
  implicit none
  
  !---Common
  integer ibuoy,ivar
  integer iyr,imon,iday
  integer syr,smon,sday
  integer eyr,emon,eday
  integer ijul,sjul,ejul
  
  !---DA system
  integer idat_a
  integer im_a,jm_a,km_a
  
  real(kind = 8),allocatable :: lon_a(:),lont_a(:),lonu_a(:),lonv_a(:)
  real(kind = 8),allocatable :: lat_a(:),latt_a(:),latu_a(:),latv_a(:)
  real(kind = 8),allocatable :: dep_a(:,:,:),dept_a(:,:,:),depu_a(:,:,:),depv_a(:,:,:)
  real(kind = 8),allocatable :: mask_a(:,:),maskt_a(:,:),masku_a(:,:),maskv_a(:,:)
  real(kind = 8),allocatable :: mean_a(:,:,:),sprd_a(:,:,:)
  
  !---ID
  integer idx,idy !Southwest corner system gridpoint nearest to obs

  !---DA system in obs. space
  real(kind = 8),allocatable :: hmean_a(:),hsprd_a(:)  
  
  !---Observation
  integer itime_o,jtime_o,ntime_o
  integer km_o
  integer sjul_o,ejul_o
  integer,allocatable :: iyr_o(:),imon_o(:),iday_o(:)
  integer,allocatable :: ijul_o(:)
  
  real(kind = 8) lon_o,lat_o
  real(kind = 8),allocatable :: dep_o(:),dat_o(:,:)

  !---Read date
  call read_argument(syr,smon,sday,eyr,emon,eday)
  call ymd_julian(syr,smon,sday,sjul)
  call ymd_julian(eyr,emon,eday,ejul)
  
  do ibuoy=1,nbuoy     
     do ivar=1,nvar
        
        !---Read buoy data
        call read_ocs(buoyname(ibuoy),varname(ivar), &
             & ntime_o,km_o,iyr_o,imon_o,iday_o,lon_o,lat_o,dep_o,dat_o)

        allocate(ijul_o(ntime_o))
        allocate(hmean_a(km_o),hsprd_a(km_o))
        
        do itime_o=1,ntime_o
           call ymd_julian(iyr_o(itime_o),imon_o(itime_o),iday_o(itime_o),ijul_o(itime_o))
        end do
        sjul_o=minval(ijul_o)
        ejul_o=maxval(ijul_o)
        
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

           allocate(lon_a(im_a),lont_a(im_a),lonu_a(im_a),lonv_a(im_a))
           allocate(lat_a(jm_a),latt_a(jm_a),latu_a(jm_a),latv_a(jm_a))
           allocate(dep_a(im_a,jm_a,km_a),dept_a(im_a,jm_a,km_a),depu_a(im_a,jm_a,km_a),depv_a(im_a,jm_a,km_a))
           allocate(mask_a(im_a,jm_a),maskt_a(im_a,jm_a),masku_a(im_a,jm_a),maskv_a(im_a,jm_a))
           allocate(mean_a(2,2,km_a),sprd_a(2,2,km_a))
           
           !---Read DA Grid
           call read_grid(idat_a,im_a,jm_a,km_a, &
                & lont_a,lonu_a,lonv_a,latt_a,latu_a,latv_a,dept_a,depu_a,depv_a,maskt_a,masku_a,maskv_a)

           if(varname(ivar) == "t" .or. varname(ivar) == "s")then
              call substitute_grid(im_a,jm_a,km_a,lont_a,latt_a,dept_a,lon_a,lat_a,dep_a)
           else if(varname(ivar) == "u")then
              call substitute_grid(im_a,jm_a,km_a,lonu_a,latu_a,depu_a,lon_a,lat_a,dep_a)
           else if(varname(ivar) == "v")then
              call substitute_grid(im_a,jm_a,km_a,lonv_a,latv_a,depv_a,lon_a,lat_a,dep_a)
           end if              
           
           !---ID
           call get_id(im_a,lon_a,1,lon_o,idx)
           if(idx == 0)then
              lon_a(:)=lon_a(:)+360.d0
              call get_id(im_a,lon_a,1,lon_o,idx)
           end if
           call get_id(jm_a,lat_a,1,lat_o,idy)
           write(*,*) "ID:",idx,idy           
           
           if(idx == 0 .or. idy == 0)then
              call deallocate_ocs(iyr_o,imon_o,iday_o,dep_o,dat_o)
              deallocate(lon_a,lont_a,lonu_a,lonv_a)
              deallocate(lat_a,latt_a,latu_a,latv_a)
              deallocate(dep_a,dept_a,depu_a,depv_a)
              deallocate(mask_a,maskt_a,masku_a,maskv_a)
              deallocate(mean_a,sprd_a)
              cycle
           end if
              
           do ijul=sjul,ejul

              !---Check available date in obs.
              if(ijul < sjul_o .or. ejul_o < ijul) cycle
              
              !---Date
              call julian_ymd(ijul,iyr,imon,iday)
              write(*,*) "Buoy:"//trim(buoyname(ibuoy))//" Variable:"//trim(varname(ivar)),iyr,imon,iday              
              
              !---Extract analysis
              call extract_data(varname(ivar),idat_a,iyr,imon,iday,idx,2,idy,2,1,km_a,mean_a,sprd_a)
              
              !---Convert to obs. space
              call convert_to_obs_space(km_a,lon_a(idx:idx+1),lat_a(idy:idy+1),dep_a(idx:idx+1,idy:idy+1,1:km_a),mean_a, &
                   & km_o,lon_o,lat_o,dep_o,hmean_a)
              call convert_to_obs_space(km_a,lon_a(idx:idx+1),lat_a(idy:idy+1),dep_a(idx:idx+1,idy:idy+1,1:km_a),sprd_a, &
                   & km_o,lon_o,lat_o,dep_o,hsprd_a)

              !---Get itime_o
              jtime_o=0
              do itime_o=1,ntime_o
                 if(ijul == ijul_o(itime_o))then
                    jtime_o=itime_o
                    exit
                 end if
              end do
                            
              !---Write H(xa)
              if(jtime_o == 0) cycle
              call write_hdata(buoyname(ibuoy),varname(ivar),idat_a,iyr,imon,iday, &
                   & km_o,lon_o,lat_o,dep_o,dat_o(:,jtime_o),hmean_a,hsprd_a)
              
           end do !ijul

           deallocate(lon_a,lont_a,lonu_a,lonv_a)
           deallocate(lat_a,latt_a,latu_a,latv_a)
           deallocate(dep_a,dept_a,depu_a,depv_a)
           deallocate(mask_a,maskt_a,masku_a,maskv_a)
           deallocate(mean_a,sprd_a)

        end do    !idat_a
                
        call deallocate_ocs(iyr_o,imon_o,iday_o,dep_o,dat_o)
        deallocate(ijul_o)
        deallocate(hmean_a,hsprd_a)
        
     end do !ivar
     
  end do    !ibuoy
     
end program main

!----------------------------------------------------------------------------------------

subroutine substitute_grid(im,jm,km,lon_in,lat_in,dep_in,lon_out,lat_out,dep_out)

  implicit none

  !---IN
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: lon_in(im),lat_in(jm),dep_in(im,jm,km)

  !---OUT
  real(kind = 8),intent(out) :: lon_out(im),lat_out(jm),dep_out(im,jm,km)
  
  lon_out(:)=lon_in(:)
  lat_out(:)=lat_in(:)
  dep_out(:,:,:)=dep_in(:,:,:)
  
end subroutine substitute_grid
