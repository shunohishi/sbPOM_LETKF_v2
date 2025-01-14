module setting

  !Setting parameter
  integer,parameter :: syr=1993,eyr=2017
  !Start/End year for making SODA monthly climatology

  integer,parameter :: nvar2d=1 !SSH
  integer,parameter :: nvar3d=4 !TSUV
  
end module setting

!____________________________________________________________________________

program main

!-------------------------------------------------------------
! Make tsdata & lbc monthly climatology |
!-------------------------------------------------------------
! 
! SODA: z coordinate --> Model: sigma coordinate
!
!-------------------------------------------------------------
!
! 1. Read model grid data
! 2. Make SODA monthly climatology
! 3. Calculate ID
! 4. Fill value in horizontal direction
! 5. Bilinear interpolation in horizontal direction
! 6. Linear interpolation in vertical direction
! 7. Fill value in vertical direction
! 8. Calculate Vertical mean
! 9. Write Data
!
! 2018.08 Created by S.Ohishi
! 2024.07 Moified by S.Ohishi  
!
!-------------------------------------------------------------

  use setting
  use mod_rmiss
  use mod_gridinfo, km_in => km
  use mod_read_soda, im_soda => im, jm_soda => jm, km_soda_in => km
  implicit none

  !Common
  integer ivar
  integer iyr,imon
  integer i,j,k
  integer km,km_soda
  integer iw,ie,js,jn
  integer :: iwt=1,iet=im,jst=1,jnt=jm !Grid number of Boudary condition: t/s/h
  integer :: iwu=2,ieu=im,jsu=2,jnu=jm !u
  integer :: iwv=2,iev=im,jsv=2,jnv=jm !v

  character(3) varname,varname2d
  
  !Model
  integer iqglobal
  integer idx(im),idy(jm),idz(im,jm,km_in)

  real(kind = 8) lont(im),lonu(im),lonv(im)
  real(kind = 8) latt(jm),latu(jm),latv(jm)
  real(kind = 8) depth(im,jm,km_in)
  real(kind = 8) fsm(im,jm),dum(im,jm),dvm(im,jm)

  real(kind = 8),allocatable :: dat(:,:,:)
  real(kind = 8),allocatable :: dat_vave(:,:)
    
  !Boundary Condition  
  real(kind = 8) sshw(jm),sshe(jm),sshs(im),sshn(im)
  real(kind = 8) uaw(jm),uae(jm),uas(im),uan(im) !Vertical averaged zonal velocity
  real(kind = 8) vaw(jm),vae(jm),vas(im),van(im) 
  real(kind = 8) tw(jm,km_in),te(jm,km_in),ts(im,km_in),tn(im,km_in)
  real(kind = 8) sw(jm,km_in),se(jm,km_in),ss(im,km_in),sn(im,km_in)
  real(kind = 8) uw(jm,km_in),ue(jm,km_in),us(im,km_in),un(im,km_in)
  real(kind = 8) vw(jm,km_in),ve(jm,km_in),vs(im,km_in),vn(im,km_in)

  !Null
  real(kind = 8) null1x(im),null1y(jm)
  real(kind = 8),allocatable ::  null(:,:,:)

  !SODA
  real(kind = 8) lont_soda(im_soda),lonu_soda(im_soda)
  real(kind = 8) latt_soda(jm_soda),latu_soda(jm_soda)
  real(kind = 8) deptht_soda(km_soda_in),depthu_soda(km_soda_in)
  
  real(kind = 8),allocatable :: dat_soda(:,:,:)
  real(kind = 8),allocatable :: clim_soda(:,:,:),pass_soda(:,:,:),miss_soda(:,:,:)

  !Intermidiate(Depth is only soda)
  real(kind = 8),allocatable :: dat_int(:,:,:)

  !---Read Model information
  allocate(null(im,jm,km_in))
  write(*,'(a)') "Read Model Grid"
  call read_grid(lonu,lonv,lont,null1x,&
       & latu,latv,latt,null1y,&
       & null,depth,fsm,dum,dvm)
  deallocate(null)
  
  !---Read SODA information
  call read_soda_info(syr,lont_soda,lonu_soda, &
       & latt_soda,latu_soda,deptht_soda,depthu_soda)
  
  if(lont(1) == lont(im-1)-360.d0 .and. lont(2) == lont(im)-360.d0)then
     iqglobal=1
  else
     iqglobal=0
  end if

  do imon=1,12
     do ivar=1,nvar2d+nvar3d

        !---Varname
        if(ivar == 1)then
           varname="ssh"
        else if(ivar == 2)then
           varname="t"
        else if(ivar == 3)then
           varname="s"
        else if(ivar == 4)then
           varname="u"
           varname2d="ua"
        else if(ivar == 5)then
           varname="v"
           varname2d="va"
        end if

        write(*,'(a,i2,a)') "----- Start "//trim(varname)//": ",imon,"-----"
        
        !---Allocate
        if(ivar <= nvar2d)then
           km_soda=1
           km=1
        else if(nvar2d+1 <= ivar .and. ivar <= nvar2d+nvar3d)then
           km_soda=km_soda_in
           km=km_in
        else
           write(*,*) "***Error: nvar2d+nvar3d"
           stop
        end if
        allocate(dat_soda(im_soda,jm_soda,km_soda))
        allocate(dat_int(im,jm,km_soda))
        allocate(dat_vave(im,jm),dat(im,jm,km))

        !---Calculate ID
        write(*,'(a)') "Calculate ID"
        if(varname == "ssh" .or. varname == "t" .or. varname == "s")then
           call cal_idlon(im_soda,lont_soda,im,lont,idx)
           call cal_idlat(jm_soda,latt_soda,jm,latt,idy)
           call cal_idz(km_soda,deptht_soda,im,jm,km,depth,idz)
           iw=iwt
           ie=iet
           js=jst
           jn=jnt
        else if(varname == "u")then
           call cal_idlon(im_soda,lonu_soda,im,lonu,idx)
           call cal_idlat(jm_soda,latu_soda,jm,latu,idy)
           call cal_idz(km_soda,depthu_soda,im,jm,km,depth,idz)
           iw=iwu
           ie=ieu
           js=jsu
           jn=jnu
        else if(varname == "v")then
           call cal_idlon(im_soda,lonu_soda,im,lonv,idx)
           call cal_idlat(jm_soda,latu_soda,jm,latv,idy)
           call cal_idz(km_soda,depthu_soda,im,jm,km,depth,idz)
           iw=iwv
           ie=iev
           js=jsv
           jn=jnv
        end if
         
        !---Make SODA Monthly climatology
        write(*,'(a)') "Make SODA Monthly Climatology"
        allocate(clim_soda(im_soda,jm_soda,km_soda))
        allocate(pass_soda(im_soda,jm_soda,km_soda))
        allocate(miss_soda(im_soda,jm_soda,km_soda))
        do iyr=syr,eyr
           write(*,'(i4,i2)') iyr,imon

           if(varname == "ssh")then
              call read_soda_2d(varname,iyr,imon,dat_soda)
           else
              call read_soda_3d(varname,iyr,imon,dat_soda)
              
           end if
           call cal_clim(iyr,syr,eyr,im_soda,jm_soda,km_soda, &
                & dat_soda,clim_soda,pass_soda,miss_soda,rmiss)
        end do !iyr
        
        dat_soda(:,:,:)=clim_soda(:,:,:)        
        deallocate(clim_soda,pass_soda,miss_soda)
        
        write(*,'(a)') "Fill value"
        call fillvalue_2d(10,im_soda,jm_soda,km_soda,dat_soda,rmiss)

        !---Bilinear interpolation in horizontal direction        
        write(*,'(a)') "Interpolation in horizontal direction"                     
        do k=1,km_soda
           if(varname == "ssh")then !SSH: SODA --> Model grid
              call bilinear_interpolation_2d( &
                   & im_soda,jm_soda,lont_soda,latt_soda,dat_soda(:,:,k), &
                   & im,jm,lont,latt,dat(:,:,k),idx,idy,rmiss)
           else if(varname == "t" .or. varname == "s")then !SODA --> Intermidiate grid
              call bilinear_interpolation_2d( &
                   & im_soda,jm_soda,lont_soda,latt_soda,dat_soda(:,:,k), &
                   & im,jm,lont,latt,dat_int(:,:,k),idx,idy,rmiss)
           else if(varname == "u")then
              call bilinear_interpolation_2d( &
                   & im_soda,jm_soda,lonu_soda,latu_soda,dat_soda(:,:,k), &
                   & im,jm,lonu,latu,dat_int(:,:,k),idx,idy,rmiss)
           else if(varname == "v")then
              call bilinear_interpolation_2d( &
                   & im_soda,jm_soda,lonu_soda,latu_soda,dat_soda(:,:,k), &
                   & im,jm,lonv,latv,dat_int(:,:,k),idx,idy,rmiss)
           end if
        end do !k
        
        if(varname == "ssh")then
           call apply_fsm(im,jm,km,dat,fsm)
        else if(varname == "t" .or. varname == "s")then
           call apply_fsm(im,jm,km_soda,dat_int,fsm)
        else if(varname == "u")then
           call apply_fsm(im,jm,km_soda,dat_int,dum)
        else if(varname == "v")then
           call apply_fsm(im,jm,km_soda,dat_int,dvm)
        end if

        write(*,'(a)') "Interpolation in vertical direction"
        write(*,'(a)') "Fill value in vertical direction"
        !Intermidiate --> Model
        if(varname == "t" .or. varname == "s")then
           call linear_interpolation_vertical( &
                & km_soda,deptht_soda,dat_int,im,jm,km,depth,fsm,dat,idz,rmiss)
           call fillvalue_vertical(im,jm,km,fsm,dat,rmiss)
           call apply_fsm(im,jm,km,dat,fsm)
        else if(varname == "u")then
           call linear_interpolation_vertical( &
                & km_soda,depthu_soda,dat_int,im,jm,km,depth,dum,dat,idz,rmiss)
           call fillvalue_vertical(im,jm,km,dum,dat,rmiss)
           call apply_fsm(im,jm,km,dat,dum)
        else if(varname == "v")then
           call linear_interpolation_vertical( &
                & km_soda,depthu_soda,dat_int,im,jm,km,depth,dvm,dat,idz,rmiss)
           call fillvalue_vertical(im,jm,km,dvm,dat,rmiss)
           call apply_fsm(im,jm,km,dat,dvm)
        end if

        if(varname == "t" .or. varname == "s" .or. varname == "u" .or. varname == "v")then
           call bottom_value(im,jm,km,dat,0.d0)
        end if
        
        !***Remove due to memory issue
        !write(*,'(a)') "Calculate dflow"
        !call cal_dflow(iwu,ieu,jsv,jnv,im,jm,km,lonu,latu,lonv,latv,depth,u,v,dum,dvm)

        write(*,'(a)') "Calculate Vertical mean"
        if(varname == "u")then
           call vertical_mean(varname,im,jm,km,depth,dum,dat,dat_vave)
        else if(varname == "v")then
           call vertical_mean(varname,im,jm,km,depth,dvm,dat,dat_vave)
        end if

        !Quasi-global
        if(iqglobal == 1)then
           
           dat(1:2,:,:)=dat(im-1:im,:,:)

           if(varname == "u" .or. varname == "v")then
              dat_vave(1:2,:)=dat_vave(im-1:im,:)
           end if   
        end if

        write(*,'(a)') "Substitute lbc"
        !ssh,ua,va
        if(varname == "ssh")then
           call substitute_lbc1d(iw,ie,js,jn,im,jm,dat(:,:,1),sshw,sshe,sshs,sshn)
        else if(varname == "t")then        
           call substitute_lbc2d(iw,ie,js,jn,im,jm,km,dat,tw,te,ts,tn)
        else if(varname == "s")then        
           call substitute_lbc2d(iw,ie,js,jn,im,jm,km,dat,sw,se,ss,sn)
        else if(varname == "u")then !ua/u
           call substitute_lbc1d(iw,ie,js,jn,im,jm,dat_vave,uaw,uae,uas,uan)
           call substitute_lbc2d(iw,ie,js,jn,im,jm,km,dat,uw,ue,us,un)
        else if(varname == "v")then !va/v
           call substitute_lbc1d(iw,ie,js,jn,im,jm,dat_vave,vaw,vae,vas,van)
           call substitute_lbc2d(iw,ie,js,jn,im,jm,km,dat,vw,ve,vs,vn)
        end if
        
        write(*,'(a)') "Write Data"
        call write_dat(varname,imon,im,jm,km,dat)
        if(varname == "u" .or. varname == "v")then
           call write_dat(varname2d,imon,im,jm,1,dat_vave)
        end if
           
        deallocate(dat_soda,dat_int,dat_vave,dat)

        write(*,'(a,i2,a)') "----- End "//trim(varname)//": ",imon,"-----"
        
     end do !ivar
        
     write(*,'(a)') "Write LBC Data"
     call write_lbc(imon,im,jm,km,&
          & sshw,sshe,sshs,sshn, &
          & uaw,uae,uas,uan,vaw,vae,vas,van, &
          & tw,te,ts,tn,sw,se,ss,sn, &
          & uw,ue,us,un,vw,ve,vs,vn)

  end do !imon

  write(*,'(a)') "===== ALL END ====="
  
end program

!------------------------------------------------------------------------
! Calculate Monthly Climatology |
!-------------------------------------------------------------------------

subroutine cal_clim(iyr,syr,eyr,im,jm,km,dat,mean,pass,miss,rmiss)

  !$use omp_lib  
  implicit none
  integer i,j,k

  integer,intent(in) :: iyr,syr,eyr
  integer,intent(in) :: im,jm,km
  
  real(kind = 8),intent(in) :: dat(im,jm,km)
  real(kind = 8),intent(inout) :: mean(im,jm,km),pass(im,jm,km),miss(im,jm,km)
  real(kind = 8),intent(in) :: rmiss

  if(iyr == syr)then
     mean(:,:,:)=0.d0
     pass(:,:,:)=0.d0
     miss(:,:,:)=0.d0
  end if

  !$omp parallel
  !$omp do private(i,j,k)    
  do k=1,km
     do j=1,jm
        do i=1,im
           if(dat(i,j,k) == rmiss)then
              miss(i,j,k)=miss(i,j,k)+1.d0
           else
              mean(i,j,k)=mean(i,j,k)+dat(i,j,k)
              pass(i,j,k)=pass(i,j,k)+1.d0
           end if
        end do
     end do
  end do
  !$omp end do
  
  if(iyr == eyr)then
     !$omp do private(i,j,k)
     do k=1,km
        do j=1,jm
           do i=1,im
              if(pass(i,j,k) == 0.d0)then
                 mean(i,j,k)=rmiss
              else
                 mean(i,j,k)=mean(i,j,k)/pass(i,j,k)
              end if
           end do
        end do
     end do
     !$omp end do
  end if
  !$omp end parallel
  
end subroutine cal_clim

!----------------------------------------------------------------------
! Linear Interpolation in vertical direction |
!----------------------------------------------------------------------
!
! As in tsclim version, 
! but for considering surface layer shallowe than depth1(1)
!
!-----------------------------------------------------------------------

subroutine linear_interpolation_vertical( &
          & km1,depth1,dat1,im2,jm2,km2,depth2,fsm,dat2,id,rmiss)

  !$use omp_lib  
  implicit none
  
  integer k1
  integer i2,j2,k2

  integer,intent(in) :: km1
  integer,intent(in) :: im2,jm2,km2
  integer,intent(in) :: id(im2,jm2,km2)
  
  real(kind = 8),intent(in) :: depth1(km1),dat1(im2,jm2,km1)
  real(kind = 8),intent(in) :: depth2(im2,jm2,km2),fsm(im2,jm2)
  real(kind = 8),intent(in) :: rmiss

  real(kind = 8),intent(inout) :: dat2(im2,jm2,km2)

  dat2(:,:,:)=rmiss

  !$omp parallel
  !$omp do private(i2,j2,k2,k1)  
  do j2=1,jm2
     do i2=1,im2

        if(fsm(i2,j2) == 0.d0)then
           dat2(i2,j2,:)=rmiss
           cycle
        end if

        do k2=1,km2

           if(depth2(i2,j2,k2) <= depth1(1))then
              dat2(i2,j2,k2)=dat1(i2,j2,1)
              cycle
           end if

           if(id(i2,j2,k2) == 0)then
              cycle
           else
              k1=id(i2,j2,k2)
           end if

           if(dat1(i2,j2,k1) == rmiss .or. dat1(i2,j2,k1+1) == rmiss)then
              dat2(i2,j2,k2)=rmiss
           else
              call linear_interpolate( &
                   & depth1(k1),depth1(k1+1),dat1(i2,j2,k1),dat1(i2,j2,k1+1), &
                   & depth2(i2,j2,k2),dat2(i2,j2,k2))
           end if
                   
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine linear_interpolation_vertical

!------------------------------------------------------------------------
! Bottom velocity |
!------------------------------------------------------------------------

subroutine bottom_value(im,jm,km,dat,value)

  implicit none

  integer i,j,k
  integer,intent(in) :: im,jm,km
  
  real(kind = 8),intent(inout) :: dat(im,jm,km)
  real(kind = 8),intent(in) :: value

  dat(:,:,km)=value

end subroutine bottom_value

!----------------------------------------------------------------------
! Mean flow at each grid for IN-OUT FLOW = 0 
!----------------------------------------------------------------------

subroutine cal_dflow(iw,ie,js,jn,im,jm,km,lonu,latu,lonv,latv,depth,u,v,dum,dvm)

  use mod_parameter
  implicit none

  integer i,j,k

  real(kind = 8) volume,area,dflow
  real(kind = 8) dlonu(jm),dlatu !Distance at each grid
  real(kind = 8) dlonv(jm),dlatv !Distance at each grid
  
  integer,intent(in) :: iw,ie !Zonal velocity boudary grid
  integer,intent(in) :: js,jn !Meridional velocity boudary grid
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: lonu(im),latu(jm)
  real(kind = 8),intent(in) :: lonv(jm),latv(jm)
  real(kind = 8),intent(in) :: depth(im,jm,km)
  real(kind = 8),intent(in) :: dum(im,jm),dvm(im,jm)

  real(kind = 8),intent(inout) :: u(im,jm,km),v(im,jm,km)
  
  do j=1,jm
     dlonu=earth*cos(pi*latu(j)/180.d0)*(lonu(2)-lonu(1))*pi/180.d0
     dlonv=earth*cos(pi*latv(j)/180.d0)*(lonv(2)-lonv(1))*pi/180.d0
  end do
  
  dlatu=(latu(2)-latu(1))*earth*pi/180.d0
  dlatv=(latv(2)-latv(1))*earth*pi/180.d0

  volume=0.d0
  area=0.d0

  do k=2,km-1

     i=iw
     do j=1,jm
        if(dum(i,j) == 0.d0)cycle
        area=area+dlatu*0.5d0*(depth(i,j,k)+depth(i-1,j,k)-depth(i,j,k-1)-depth(i-1,j,k-1))
        volume=volume+u(i,j,k)*dlatu*0.5d0*(depth(i,j,k)+depth(i-1,j,k)-depth(i,j,k-1)-depth(i-1,j,k-1))
     end do

     i=ie
     do j=1,jm
        if(dum(i,j) == 0.d0)cycle
        area=area+dlatu*0.5d0*(depth(i,j,k)+depth(i-1,j,k)-depth(i,j,k-1)-depth(i-1,j,k-1))
        volume=volume-u(i,j,k)*dlatu*0.5d0*(depth(i,j,k)+depth(i-1,j,k)-depth(i,j,k-1)-depth(i-1,j,k-1))
     end do
     
     j=js
     do i=1,im
        if(dvm(i,j) == 0.d0)cycle
        area=area+dlonv(j)*0.5d0*(depth(i,j,k)+depth(i,j-1,k)-depth(i,j,k-1)-depth(i,j-1,k-1))
        volume=volume+v(i,j,k)*dlonv(j)*0.5d0*(depth(i,j,k)+depth(i,j-1,k)-depth(i,j,k-1)-depth(i,j-1,k-1))
     end do

     j=jn
     do i=1,im
        if(dvm(i,j) == 0.d0)cycle
        area=area+dlonv(j)*0.5d0*(depth(i,j,k)+depth(i,j-1,k)-depth(i,j,k-1)-depth(i,j-1,k-1))
        volume=volume-v(i,j,k)*dlonv(j)*0.5d0*(depth(i,j,k)+depth(i,j-1,k)-depth(i,j,k-1)-depth(i,j-1,k-1))
     end do        

  end do
  
  dflow=volume/area
  write(*,*) "dflow:",dflow,"(m/s)"

  !Modify velocity
  do k=1,km-1

     i=iw
     do j=1,jm
        if(dum(i,j) == 0.d0)cycle
        u(i,j,k)=u(i,j,k)-dflow
     end do

     i=ie
     do j=1,jm
        if(dum(i,j) == 0.d0)cycle
        u(i,j,k)=u(i,j,k)+dflow
     end do

     j=js
     do i=1,im
        if(dvm(i,j) == 0.d0)cycle
        v(i,j,k)=v(i,j,k)-dflow
     end do

     j=jn
     do i=1,im
        if(dvm(i,j) == 0.d0)cycle
        v(i,j,k)=v(i,j,k)+dflow
     end do

  end do

end subroutine cal_dflow

!------------------------------------------------------------------------
! Vertical mean |
!----------------
!
! Use tranpezoidal rule
!
!------------------------------------------------------------------------

subroutine vertical_mean(var,im,jm,km,depth,fsm,dat,mean)

  implicit none

  integer i,j,k

  real(kind = 8) dep1,dep2,bottom

  !IN
  integer,intent(in) :: im,jm,km
  
  real(kind = 8),intent(in) :: depth(im,jm,km)
  real(kind = 8),intent(in) :: fsm(im,jm)
  real(kind = 8),intent(in) :: dat(im,jm,km)

  character(1),intent(in) :: var !"T"/"U"/"V"

  !OUT  
  real(kind = 8),intent(out) :: mean(im,jm)

  do j=1,jm
     do i=1,im

        mean(i,j)=0.d0

        if(fsm(i,j)==0.d0)then
           cycle
        end if

        do k=2,km-1

           !dep1: k-1, dep2:k
           if(var == "t")then
              dep1=depth(i,j,k-1)
              dep2=depth(i,j,k)
           elseif(var == "u")then
              if(i == 1)then
                 dep1=depth(i,j,k-1)
                 dep2=depth(i,j,k)
              else
                 dep1=0.5d0*(depth(i,j,k-1)+depth(i-1,j,k-1))
                 dep2=0.5d0*(depth(i,j,k)+depth(i-1,j,k))
              end if
           elseif(var == "v")then
              if(j == 1)then
                 dep1=depth(i,j,k-1)
                 dep2=depth(i,j,k)
              else
                 dep1=0.5d0*(depth(i,j,k-1)+depth(i,j-1,k-1))
                 dep2=0.5d0*(depth(i,j,k)+depth(i,j-1,k))
              end if
           end if

           mean(i,j)=mean(i,j)+0.5d0*(dep2-dep1)*(dat(i,j,k-1)+dat(i,j,k))
           
        end do

        if(var == "t")then
           bottom=depth(i,j,km-1)
        elseif(var == "u")then
           if(i == 1)then
              bottom=depth(i,j,km-1)
           else
              bottom=0.5d0*(depth(i,j,km-1)+depth(i-1,j,km-1))
           end if
        elseif(var == "v")then
           if(j == 1)then
              bottom=depth(i,j,km-1)
           else
              bottom=0.5d0*(depth(i,j,km-1)+depth(i,j-1,km-1))
           end if
        end if

        mean(i,j)=mean(i,j)/bottom

     end do
  end do
  
end subroutine vertical_mean

!------------------------------------------------------------------------
! Substitute Lateral Boundary Condition |
!------------------------------------------------------------------------

subroutine substitute_lbc1d(iw,ie,js,jn,im,jm,dat,datw,date,dats,datn)

  implicit none

  !IN
  integer,intent(in) :: iw,ie
  integer,intent(in) :: js,jn
  integer,intent(in) :: im,jm
  
  real(kind = 8),intent(in) :: dat(im,jm)

  !OUT
  real(kind = 8),intent(out) :: datw(jm),date(jm)
  real(kind = 8),intent(out) :: dats(im),datn(im)

  datw(:)=dat(iw,:)
  date(:)=dat(ie,:)

  dats(:)=dat(:,js)
  datn(:)=dat(:,jn)

end subroutine substitute_lbc1d

!-------------------------------

subroutine substitute_lbc2d(iw,ie,js,jn,im,jm,km,dat,datw,date,dats,datn)

  implicit none

  !IN
  integer,intent(in) :: iw,ie
  integer,intent(in) :: js,jn
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: dat(im,jm,km)

  !OUT
  real(kind = 8),intent(out) :: datw(jm,km),date(jm,km)
  real(kind = 8),intent(out) :: dats(im,km),datn(im,km)

  datw(:,:)=dat(iw,:,:)
  date(:,:)=dat(ie,:,:)
  dats(:,:)=dat(:,js,:)
  datn(:,:)=dat(:,jn,:)
  
end subroutine substitute_lbc2d

!------------------------------------------------------------------------
! Write Data |
!------------------------------------------------------------------------

subroutine write_lbc(imon,im,jm,km,&
     & sshw,sshe,sshs,sshn, &
     & uaw,uae,uas,uan,vaw,vae,vas,van, &
     & tw,te,ts,tn,sw,se,ss,sn, &
     & uw,ue,us,un,vw,ve,vs,vn)

  use netcdf
  implicit none

  !---Common
  integer status,access,system
  integer ncid,varid

  character(100) filename
  
  integer,intent(in) :: imon
  integer,intent(in) :: im,jm,km

  real(kind = 4) tmp1dx(im),tmp1dy(jm)
  real(kind = 4) tmp2dxz(im,km),tmp2dyz(jm,km)

  !---IN
  real(kind = 8),intent(in) :: sshw(jm),sshe(jm),sshs(im),sshn(im)
  real(kind = 8),intent(in) :: uaw(jm),uae(jm),uas(im),uan(im)
  real(kind = 8),intent(in) :: vaw(jm),vae(jm),vas(im),van(im) 
  real(kind = 8),intent(in) :: tw(jm,km),te(jm,km),ts(im,km),tn(im,km)
  real(kind = 8),intent(in) :: sw(jm,km),se(jm,km),ss(im,km),sn(im,km)
  real(kind = 8),intent(in) :: uw(jm,km),ue(jm,km),us(im,km),un(im,km)
  real(kind = 8),intent(in) :: vw(jm,km),ve(jm,km),vs(im,km),vn(im,km)
  
  filename="../in/lbc_mclim.nc"

  status=access(trim(filename)," ")

  if(imon == 1 .and. status == 0)then
     status=system("rm -f "//trim(filename))
  end if

  if(imon == 1)then
     call make_ncfile_lbc_mclim(im,jm,km,filename)
     !status=system("ncgen -o "//trim(filename)//" ../hdr/lbc_mclim.hdr")
  end if

  status=nf90_open(trim(filename),nf90_write,ncid)
  
  !ssh
  tmp1dy(:)=real(sshe(:))
  status=nf90_inq_varid(ncid,"ele",varid)
  status=nf90_put_var(ncid,varid,tmp1dy,(/1,imon/),(/jm,1/))

  tmp1dy(:)=real(sshw(:))  
  status=nf90_inq_varid(ncid,"elw",varid)
  status=nf90_put_var(ncid,varid,tmp1dy,(/1,imon/),(/jm,1/))

  tmp1dx(:)=real(sshn(:))  
  status=nf90_inq_varid(ncid,"eln",varid)
  status=nf90_put_var(ncid,varid,tmp1dx,(/1,imon/),(/im,1/))

  tmp1dx(:)=real(sshs(:))    
  status=nf90_inq_varid(ncid,"els",varid)
  status=nf90_put_var(ncid,varid,tmp1dx,(/1,imon/),(/im,1/))

  !ua
  tmp1dy(:)=real(uae(:))  
  status=nf90_inq_varid(ncid,"uabe",varid)
  status=nf90_put_var(ncid,varid,tmp1dy,(/1,imon/),(/jm,1/))

  tmp1dy(:)=real(uaw(:))    
  status=nf90_inq_varid(ncid,"uabw",varid)
  status=nf90_put_var(ncid,varid,tmp1dy,(/1,imon/),(/jm,1/))

  tmp1dx(:)=real(uan(:))    
  status=nf90_inq_varid(ncid,"uabn",varid)
  status=nf90_put_var(ncid,varid,tmp1dx,(/1,imon/),(/im,1/))

  tmp1dx(:)=real(uas(:))      
  status=nf90_inq_varid(ncid,"uabs",varid)
  status=nf90_put_var(ncid,varid,tmp1dx,(/1,imon/),(/im,1/))

  !va
  tmp1dy(:)=real(vae(:))  
  status=nf90_inq_varid(ncid,"vabe",varid)
  status=nf90_put_var(ncid,varid,tmp1dy,(/1,imon/),(/jm,1/))

  tmp1dy(:)=real(vaw(:))    
  status=nf90_inq_varid(ncid,"vabw",varid)
  status=nf90_put_var(ncid,varid,tmp1dy,(/1,imon/),(/jm,1/))

  tmp1dx(:)=real(van(:))      
  status=nf90_inq_varid(ncid,"vabn",varid)
  status=nf90_put_var(ncid,varid,tmp1dx,(/1,imon/),(/im,1/))

  tmp1dx(:)=real(vas(:))        
  status=nf90_inq_varid(ncid,"vabs",varid)
  status=nf90_put_var(ncid,varid,tmp1dx,(/1,imon/),(/im,1/))

  !temperature
  tmp2dyz(:,:)=real(te(:,:))
  status=nf90_inq_varid(ncid,"tbe",varid)
  status=nf90_put_var(ncid,varid,tmp2dyz,(/1,1,imon/),(/jm,km,1/))

  tmp2dyz(:,:)=real(tw(:,:))  
  status=nf90_inq_varid(ncid,"tbw",varid)
  status=nf90_put_var(ncid,varid,tmp2dyz,(/1,1,imon/),(/jm,km,1/))

  tmp2dxz(:,:)=real(tn(:,:))    
  status=nf90_inq_varid(ncid,"tbn",varid)
  status=nf90_put_var(ncid,varid,tmp2dxz,(/1,1,imon/),(/im,km,1/))

  tmp2dxz(:,:)=real(ts(:,:))      
  status=nf90_inq_varid(ncid,"tbs",varid)
  status=nf90_put_var(ncid,varid,tmp2dxz,(/1,1,imon/),(/im,km,1/))

  !salinity
  tmp2dyz(:,:)=real(se(:,:))
  status=nf90_inq_varid(ncid,"sbe",varid)
  status=nf90_put_var(ncid,varid,tmp2dyz,(/1,1,imon/),(/jm,km,1/))

  tmp2dyz(:,:)=real(sw(:,:))  
  status=nf90_inq_varid(ncid,"sbw",varid)
  status=nf90_put_var(ncid,varid,tmp2dyz,(/1,1,imon/),(/jm,km,1/))

  tmp2dxz(:,:)=real(sn(:,:))      
  status=nf90_inq_varid(ncid,"sbn",varid)
  status=nf90_put_var(ncid,varid,tmp2dxz,(/1,1,imon/),(/im,km,1/))

  tmp2dxz(:,:)=real(ss(:,:))      
  status=nf90_inq_varid(ncid,"sbs",varid)
  status=nf90_put_var(ncid,varid,tmp2dxz,(/1,1,imon/),(/im,km,1/))

  !zonal velocity
  tmp2dyz(:,:)=real(ue(:,:))
  status=nf90_inq_varid(ncid,"ube",varid)
  status=nf90_put_var(ncid,varid,tmp2dyz,(/1,1,imon/),(/jm,km,1/))

  tmp2dyz(:,:)=real(uw(:,:))  
  status=nf90_inq_varid(ncid,"ubw",varid)
  status=nf90_put_var(ncid,varid,tmp2dyz,(/1,1,imon/),(/jm,km,1/))

  tmp2dxz(:,:)=real(un(:,:))      
  status=nf90_inq_varid(ncid,"ubn",varid)
  status=nf90_put_var(ncid,varid,tmp2dxz,(/1,1,imon/),(/im,km,1/))

  tmp2dxz(:,:)=real(us(:,:))      
  status=nf90_inq_varid(ncid,"ubs",varid)
  status=nf90_put_var(ncid,varid,tmp2dxz,(/1,1,imon/),(/im,km,1/))

  !meridional velocity
  tmp2dyz(:,:)=real(ve(:,:))
  status=nf90_inq_varid(ncid,"vbe",varid)
  status=nf90_put_var(ncid,varid,tmp2dyz,(/1,1,imon/),(/jm,km,1/))

  tmp2dyz(:,:)=real(vw(:,:))  
  status=nf90_inq_varid(ncid,"vbw",varid)
  status=nf90_put_var(ncid,varid,tmp2dyz,(/1,1,imon/),(/jm,km,1/))

  tmp2dxz(:,:)=real(vn(:,:))      
  status=nf90_inq_varid(ncid,"vbn",varid)
  status=nf90_put_var(ncid,varid,tmp2dxz,(/1,1,imon/),(/im,km,1/))

  tmp2dxz(:,:)=real(vs(:,:))      
  status=nf90_inq_varid(ncid,"vbs",varid)
  status=nf90_put_var(ncid,varid,tmp2dxz,(/1,1,imon/),(/im,km,1/))

  status=nf90_close(ncid)

end subroutine write_lbc
  
!----------------------------------------

subroutine write_dat(varname,imon,im,jm,km,dat)

  use mod_gridinfo, only: km_in => km
  use netcdf
  implicit none

  !---Common
  integer status,access,system
  integer ncid,varid

  character(100) filename

  !---IN  
  integer,intent(in) :: imon
  integer,intent(in) :: im,jm,km
  
  real(kind = 8),intent(in) :: dat(im,jm,km)

  character(*),intent(in) :: varname
  
  filename="../in/tsdata_mclim.nc"

  !Check file
  status=access(trim(filename)," ")

  if(imon == 1 .and. varname == "ssh" .and. status == 0)then
     status=system("rm -f "//trim(filename))
  end if

  if(imon==1 .and. varname == "ssh")then
     call make_ncfile_tsdata_mclim(im,jm,km_in,filename)
     !status=system("ncgen -o "//trim(filename)//" ../hdr/tsdata_mclim.hdr")    
  end if
  
  status=nf90_open(trim(filename),nf90_write,ncid)

  if(varname == "t")then
     status=nf90_inq_varid(ncid,"temp",varid)
  else if(varname == "s")then
     status=nf90_inq_varid(ncid,"sal",varid)
  else
     status=nf90_inq_varid(ncid,trim(varname),varid)
  end if

  if(km == 1)then
     status=nf90_put_var(ncid,varid,real(dat(:,:,:)),(/1,1,imon/),(/im,jm,1/))
  else
     status=nf90_put_var(ncid,varid,real(dat(:,:,:)),(/1,1,1,imon/),(/im,jm,km,1/))
  end if
  
  status=nf90_close(ncid)
  
end subroutine write_dat
