program main

  use setting
  use mod_julian
  use mod_stat
  use mod_io
  use mod_rmiss
  implicit none
  
  !---Common
  integer k
  integer ibuoy
  integer ivar
  integer ijul,sjul,ejul
  integer sjul_arg,ejul_arg
  integer idat_a,jdat_a
  
  integer iyr,imon,iday  
  integer syr,smon,sday
  integer eyr,emon,eday
  
  !---Analysis in Obs. space
  real(kind = 8),allocatable :: hmean_a(:),hsprd_a(:)
  
  !---Observation
  integer sjul_o
  integer km_o
  
  real(kind = 8) lon_o,lat_o
  real(kind = 8),allocatable :: dep_o(:),dat_o(:)

  !---Static
  !Monthly
  integer,allocatable :: num_mave(:,:,:,:)
  
  real(kind = 8),allocatable :: bias_mave(:,:,:,:) !km,ndata,nmon,nyr
  real(kind = 8),allocatable :: abias_mave(:,:,:,:)
  real(kind = 8),allocatable :: rmsd_mave(:,:,:,:)
  real(kind = 8),allocatable :: sprd_mave(:,:,:,:)
  
  !ALL
  integer,allocatable :: num_ave(:,:) !km,ndata
  integer,allocatable :: bias_dof_ave(:,:,:) !km,ndata,ndata
  integer,allocatable :: rmsd_dof_ave(:,:,:)

  real(kind = 8),allocatable :: bias_ave(:,:) !km,ndata
  real(kind = 8),allocatable :: bias_tcrit_ave(:,:,:),bias_tval_ave(:,:,:) !km,ndata,ndata
  real(kind = 8),allocatable :: rmsd_ave(:,:)
  real(kind = 8),allocatable :: rmsd_tcrit_ave(:,:,:),rmsd_tval_ave(:,:,:)
  real(kind = 8),allocatable :: sprd_ave(:,:)

  write(*,*) "### START: Validation vs. OCS buoys"
  
  !---Read start and end dates
  call read_argument(syr,smon,sday,eyr,emon,eday)

  write(*,*) "Start date (validation):",syr,smon,sday
  write(*,*) "End date (validation):",eyr,emon,eday
  
  do ibuoy=1,nbuoy
     do ivar=1,nvar

        !---Start & End julian day 
        call ymd_julian(syr,smon,sday,sjul_arg)
        call ymd_julian(eyr,emon,eday,ejul_arg)

        if(ibuoy == 1 .and. (ivar == 1 .or. ivar == 2))then
           call ymd_julian(2004,6,16,sjul_o)
        else if(ibuoy == 1 .and. (ivar == 3 .or. ivar == 4))then
           call ymd_julian(2005,5,30,sjul_o)
        else if(ibuoy == 2)then
           call ymd_julian(2007,6,8,sjul_o)
        end if

        sjul=max(sjul_arg,sjul_o)
        ejul=ejul_arg
                
        !---Get grid info (*Reference date: eyr,rmon,eday)
        call read_hdata(buoyname(ibuoy),varname(ivar),ndat_a,eyr,emon,eday, &
             & km_o,lon_o,lat_o,dep_o,dat_o,hmean_a,hsprd_a)        

        if(km_o == 0)then
           write(*,*) "***Error: Buoy => "//trim(buoyname(ibuoy))//" Var => "//trim(varname(ivar))
           write(*,*) "***Error: Read km_o => ",km_o
           cycle
        end if
        
        !---Allocate
        allocate(num_mave(km_o,ndat_a,12,syr:eyr))
        allocate(bias_mave(km_o,ndat_a,12,syr:eyr))
        allocate(abias_mave(km_o,ndat_a,12,syr:eyr))
        allocate(rmsd_mave(km_o,ndat_a,12,syr:eyr))
        allocate(sprd_mave(km_o,ndat_a,12,syr:eyr))

        allocate(num_ave(km_o,ndat_a))
        allocate(bias_ave(km_o,ndat_a))
        allocate(bias_dof_ave(km_o,ndat_a,ndat_a),bias_tcrit_ave(km_o,ndat_a,ndat_a),bias_tval_ave(km_o,ndat_a,ndat_a))
        allocate(rmsd_ave(km_o,ndat_a))
        allocate(rmsd_dof_ave(km_o,ndat_a,ndat_a),rmsd_tcrit_ave(km_o,ndat_a,ndat_a),rmsd_tval_ave(km_o,ndat_a,ndat_a))
        allocate(sprd_ave(km_o,ndat_a))

        !---Deallocate
        call deallcate_hdata(dep_o,dat_o,hmean_a,hsprd_a)
        
        !---Initialize
        num_mave(:,:,:,:)=0
        bias_mave(:,:,:,:)=0.d0
        rmsd_mave(:,:,:,:)=0.d0
        sprd_mave(:,:,:,:)=0.d0

        num_ave(:,:)=0
        bias_ave(:,:)=0.d0
        rmsd_ave(:,:)=0.d0
        sprd_ave(:,:)=0.d0
        
        do idat_a=1,ndat_a           
           do ijul=sjul,ejul

              call julian_ymd(ijul,iyr,imon,iday)
              write(*,'(a,3i6)') &
                   &"Buoy: "//trim(buoyname(ibuoy))//&
                   &" Variable: "//trim(varname(ivar))//&
                   &" Data: "//trim(datname(idat_a)), &
                   & iyr,imon,iday
              
              !---Read data
              call read_hdata(buoyname(ibuoy),varname(ivar),idat_a,iyr,imon,iday, &
                   & km_o,lon_o,lat_o,dep_o,dat_o,hmean_a,hsprd_a)        

              if(km_o == 0) cycle              
              
              !---Stat
              !Monthly
              call stat_1d_add(km_o,hmean_a,hsprd_a,dat_o, &
                   & num_mave(:,idat_a,imon,iyr),  &
                   & bias_mave(:,idat_a,imon,iyr), &
                   & rmsd_mave(:,idat_a,imon,iyr), &
                   & sprd_mave(:,idat_a,imon,iyr))
              !ALL
              call stat_1d_add(km_o,hmean_a,hsprd_a,dat_o, &
                   & num_ave(:,idat_a),bias_ave(:,idat_a),rmsd_ave(:,idat_a),sprd_ave(:,idat_a))

              !---Write data
              if(lwrite_obs)then
                 call write_obs(buoyname(ibuoy),varname(ivar),datname(idat_a),iyr,imon,iday,km_o,dep_o,hmean_a,dat_o)
              end if
              
              !---Deallocate
              call deallcate_hdata(dep_o,dat_o,hmean_a,hsprd_a)
              
           end do !ijul
        end do    !idat_a

        !---Finalize
        !Monthly
        do iyr=syr,eyr
           do imon=1,12
              do idat_a=1,ndat_a
                 call stat_1d_end(km_o,num_mave(:,idat_a,imon,iyr), &
                      & bias_mave(:,idat_a,imon,iyr), &
                      & rmsd_mave(:,idat_a,imon,iyr), &
                      & sprd_mave(:,idat_a,imon,iyr))
              end do
           end do
        end do
        !ALL
        do idat_a=1,ndat_a
           call stat_1d_end(km_o,num_ave(:,idat_a),bias_ave(:,idat_a),rmsd_ave(:,idat_a),sprd_ave(:,idat_a))
        end do
           
        !---Paired t-test
        do k=1,km_o

           call substitute_rmiss(ndat_a,bias_dof_ave(k,:,:),bias_tcrit_ave(k,:,:),bias_tval_ave(k,:,:))
           call substitute_rmiss(ndat_a,rmsd_dof_ave(k,:,:),rmsd_tcrit_ave(k,:,:),rmsd_tval_ave(k,:,:))

           do idat_a=1,ndat_a
              call convert_absolute_bias((eyr-syr+1)*12,bias_mave(k,idat_a,:,:),abias_mave(k,idat_a,:,:))
           end do !idat_a
           
           do idat_a=1,ndat_a
              do jdat_a=1,ndat_a

                 if(idat_a == jdat_a) cycle

                 !Bias
                 call paired_t_test((eyr-syr+1)*12, &
                      & abias_mave(k,idat_a,:,:),abias_mave(k,jdat_a,:,:), &
                      & bias_dof_ave(k,idat_a,jdat_a), &
                      & bias_tcrit_ave(k,idat_a,jdat_a), &
                      & bias_tval_ave(k,idat_a,jdat_a))

                 !RMSD
                 call paired_t_test((eyr-syr+1)*12, &
                      & rmsd_mave(k,idat_a,:,:),rmsd_mave(k,jdat_a,:,:), &
                      & rmsd_dof_ave(k,idat_a,jdat_a), &
                      & rmsd_tcrit_ave(k,idat_a,jdat_a), &
                      & rmsd_tval_ave(k,idat_a,jdat_a))
              end do !jdat_a

           end do    !idat_a
        end do       !k

        !---Get grid info (*Reference date: eyr,rmon,eday)
        call read_hdata(buoyname(ibuoy),varname(ivar),ndat_a,eyr,emon,eday, &
             & km_o,lon_o,lat_o,dep_o,dat_o,hmean_a,hsprd_a)        
        
        !---Write data
        call write_mave(buoyname(ibuoy),varname(ivar),syr,eyr,ndat_a,km_o,dep_o,num_mave, &
             & bias_mave,rmsd_mave,sprd_mave)
        call write_ave(buoyname(ibuoy),varname(ivar),sjul,ejul,ndat_a,km_o,dep_o,num_ave, &
             & bias_ave,bias_dof_ave,bias_tcrit_ave,bias_tval_ave, &
             & rmsd_ave,rmsd_dof_ave,rmsd_tcrit_ave,rmsd_tval_ave, &
             & sprd_ave)
                     
        !---Deallocate
        deallocate(num_mave)
        deallocate(bias_mave)
        deallocate(abias_mave)
        deallocate(rmsd_mave)
        deallocate(sprd_mave)

        deallocate(num_ave)
        deallocate(bias_ave)
        deallocate(bias_dof_ave,bias_tcrit_ave,bias_tval_ave)
        deallocate(rmsd_ave)
        deallocate(rmsd_dof_ave,rmsd_tcrit_ave,rmsd_tval_ave)
        deallocate(sprd_ave)
        
        call deallcate_hdata(dep_o,dat_o,hmean_a,hsprd_a)
        
     end do       !ivar
  end do          !ibouy

end program main
