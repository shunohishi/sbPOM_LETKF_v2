!-------------------------------------------------------------
! Make grid.nc |
!---------------
!
! 1. Calculate im,jm
! 2. Make model grid
! 3. Read ETOPO1 & Modify
! 4. Calculate ID
! 5. Bilinear interpolation
! 6. Modify topography
!
!~~~ iswitch = 1
! 7.Make fsm.dat
! !!! 8.Modify fsm.dat !!! --> Manual edit
!~~~~~~~~~~~~~~~
!
!~~~ iswitch = 2
! 9.Read fsm.dat
! 10. Calculate fsm,dum,dvm
! 11. Make sigma coordinate
! 12. Smoothing topography using gauss filter
! 13. Write Data
!
!-------------------------------------------------------------
!
! Created by 
! S.Ohishi 2018.07
! Modified by 
! S.Ohishi 2019.07 modify make_sigma, add time_step
! S.Ohishi 2020.03 modify modify_topography_fsm
! S.Ohishi 2023.03 modify all
! S.Ohishi 2023.12 modify sub_smooth_topo.f90
! S.Ohishi 2024.07 check operation
! S.Ohishi 2024.11 change sigma coorinate to generalized s-coordinate
!-------------------------------------------------------------

module setting

  !Switch
  integer,parameter :: iswitch=2 !1:Make ascii --> 2. Make grid.nc
  integer,parameter :: ivcoord=1 !1: Sigma-coordinate, 2: JAMSTEC s-coordinate, 3: ROMS s-coordinate 

  !Model information
  integer,parameter :: nproc=64                         !the number of processor
  integer,parameter :: km=75                            !the number of sigma layer
  real(kind = 8),parameter :: lons=105.4d0,lone=255.0d0   !longitude [start,end]
                                                        !*cyclic if lons=0. and lone=360.
  real(kind = 8),parameter :: lats=10.0d0,late=62.8d0  !latitude [start,end]
  real(kind = 8),parameter :: res=0.10d0                !horizontal resolution

  !Topography information
  !Topography [hmin,hmax], *WOA18: max depth=5500 m
  real(kind = 8),parameter :: hmin=-20.d0,hmax=-6500.d0
  real(kind = 8),parameter :: botsw=0.5d0
  integer,parameter :: ngrid=5
  real(kind = 8),parameter :: lon_sea(ngrid)=(/145.d0,358.d0,20.d0,5.d0,133.4d0/)
  real(kind = 8),parameter :: lat_sea(ngrid)=(/35.d0,36.d0,35.d0,70.d0,34.2d0/)
  
  !Smoothing parameter
  real(kind = 8),parameter :: dh=0.1d0     !(H2-H1)/(H2+H1) < dh
  real(kind = 8),parameter :: escale=100.d3 !Gaussian e-folding scale[m]
  real(kind = 8),parameter :: delta=1.d-3  !For rounding error
  
end module setting

!_______________________________________________________________________

program main

  use mod_rmiss
  use mod_parameter
  use mod_read_etopo1, im_etp => im, jm_etp => jm
  use setting
  implicit none

  !Model variable
  integer i,j,im,jm
  integer,allocatable :: idlon(:),idlat(:)
  real(kind = 8),allocatable :: lonc(:),lonp(:),lonu(:),lonv(:),dlon(:)
  real(kind = 8),allocatable :: latc(:),latp(:),latu(:),latv(:)
  real(kind = 8) dlat
  real(kind = 8),allocatable :: z(:,:,:),zz(:,:,:)
  real(kind = 8),allocatable :: topo(:,:),topo_dif(:,:),ratio(:,:)
  real(kind = 8),allocatable :: fsm(:,:),dum(:,:),dvm(:,:)
  real(kind = 8),allocatable :: tmp(:,:)

  !Etopo1 variable
  real(kind = 8) lon_etp(im_etp),lat_etp(jm_etp)
  real(kind = 8) topo_etp(im_etp,jm_etp)
  real(kind = 8) tmp_etp(im_etp,jm_etp)

  !Calculate im,jm & allocate varibles
  im=nint((lone-lons)/res)+2
  jm=nint((late-lats)/res)+2
  write(*,*) "the number of grid (im,jm):",im,jm
  if(mod((im-2)*(jm-2),nproc) == 0)then
     write(*,'(a)') "Clear: mod((im-2)*(jm-2),nproc) == 0"
  else
     write(*,'(a)') "***Error: mod((im-2)*(jm-2),nproc) /= 0"
     stop
  end if

  allocate(idlon(im),idlat(jm))
  allocate(lonc(im),lonp(im),lonu(im),lonv(im),dlon(jm))
  allocate(latc(jm),latp(jm),latu(jm),latv(jm))
  allocate(z(im,jm,km),zz(im,jm,km))
  allocate(topo(im,jm),topo_dif(im,jm),ratio(im,jm))
  allocate(fsm(im,jm),dum(im,jm),dvm(im,jm))
  allocate(tmp(im,jm))

  !---Make model lon,lat
  write(*,'(a)') "Make positoin"
  call make_position(im,jm,lonc,lonp,lonu,lonv,dlon,latc,latp,latu,latv,dlat)

  !---Read Topography
  write(*,'(a)') "Read Etopo1"
  call read_etopo1(lon_etp,lat_etp,topo_etp)
  call modify_topography(im_etp,jm_etp,topo_etp)
  call fill_topography(im_etp,jm_etp,topo_etp)
  
  !---Make model topography using biliniar interpolation
  write(*,'(a)') "Apply bilinear interpolation"
  call cal_idlon(im_etp,lon_etp,im,lonp,idlon)
  call cal_idlat(jm_etp,lat_etp,jm,latp,idlat)
  call bilinear_interpolation_2d &
       & (im_etp,jm_etp,lon_etp,lat_etp,topo_etp, &
       & im,jm,lonp,latp,topo,idlon,idlat,rmiss)
  
  !---Modify topography
  write(*,'(a)') "Modify topography"
  call modify_topography(im,jm,topo)
  call landsea_mask(1,im,jm,lonp,latp,topo)
  call fill_topography(im,jm,topo)

  !---Make fsm & Write the Ascii
  if(iswitch == 1)then

     write(*,'(a)') "Make fsm & Write the Ascii"
     call make_write_fsm(im,jm,topo,fsm,rmiss)
     write(*,'(a)') "Make fsm ascii file & End program"
     write(*,'(a)') "--- Please modify fsm.dat                  ---"
     write(*,'(a)') "--- using emacs(M-x toggle-truncate-lines) ---"
     write(*,'(a)') "========== End =========="

     deallocate(idlon,idlat)
     deallocate(lonc,lonp,lonu,lonv,dlon)
     deallocate(latc,latp,latu,latv)
     deallocate(z,zz)
     deallocate(topo)
     deallocate(fsm,dum,dvm)
     deallocate(tmp)

     stop

  end if

  !---Read fsm
  if(iswitch == 2)then
     write(*,'(a)') "Read fsm"
     call read_fsm(im,jm,fsm)
     call landsea_mask(2,im,jm,lonp,latp,fsm)
  end if
  
  ! modify topography using fsm
  call modify_topography_fsm(im,jm,topo,fsm)
  
  !Calculate mask
  write(*,'(a)') "Calculate Mask(fsm,dum,dvm)"
  call cal_mask(im,jm,fsm,dum,dvm)

  if(lons == lone-360.d0)then
     fsm(im-1:im,:)=fsm(1:2,:)
     dum(im-1:im,:)=dum(1:2,:)
     dvm(im-1:im,:)=dvm(1:2,:)
  end if
  
  !Apply slpmin for smoothing topography
  write(*,'(a)') "Apply slpmin for smoothing topography"
  call smooth_topo(im,jm,lonp,latp,topo,topo_dif,ratio,fsm)

  if(lons == lone-360.d0)then
     topo(im-1:im,:)=topo(1:2,:)
     topo_dif(im-1:im,:)=topo_dif(1:2,:)
     ratio(im-1:im,:)=ratio(1:2,:)
  end if

  write(*,'(a)') "Make generalized s coordinate"
  if(ivcoord == 1)then
     call make_sigma_coordinate(im,jm,topo,fsm,z,zz)
  else if(ivcoord == 2)then
     call make_jamstec_scoordinate(im,jm,topo,fsm,z,zz)
  else if(ivcoord == 3)then
     call make_roms_scoordinate(im,jm,topo,fsm,z,zz)
  else
     write(*,*) "Error: Not found suitable ivcoord"
     stop
  end if
  
  !Write grid.nc
  write(*,'(a)') "write ../in/grid.nc"
  
  call write_data(im,jm,km,z,zz,dlon,dlat,&
       & lonc,lonp,lonu,lonv,latc,latp,latu,latv,topo,topo_dif,ratio,fsm,dum,dvm)

  !Time step
  call time_step(jm,dlon,hmax)

  !Make mod_gridinfo.f90
  write(*,'(a)') "Make mod_gridinfo.f90"
  call make_mod_gridinfo(im,jm,km)

  deallocate(idlon,idlat)
  deallocate(lonc,lonp,lonu,lonv,dlon)
  deallocate(latc,latp,latu,latv)
  deallocate(z,zz)
  deallocate(topo)
  deallocate(fsm,dum,dvm)
  deallocate(tmp)

  write(*,'(a)') "========== End =========="
  
end program main

!___________________________________________________________________
!---------------------------------------------------------------
! Make position |
!---------------------------------------------------------------- 

subroutine make_position &
     & (im,jm,lonc,lonp,lonu,lonv,dlon,latc,latp,latu,latv,dlat)

  use setting
  use mod_parameter
  implicit none
  
  integer i,j,im,jm
  real(kind = 8),intent(out) :: lonc(im),lonp(im),lonu(im),lonv(im),dlon(jm)
  real(kind = 8),intent(out) :: latc(jm),latp(jm),latu(jm),latv(jm),dlat

  !Longitude
  if(lons == lone-360.d0)then !cyclic boudary condition
     do i=1,im
        lonc(i)=(dble(i)-1.d0)*res+lons
        lonp(i)=(dble(i)-0.5d0)*res+lons !0.5*(lonc(i)+lonc(i+1))
        lonu(i)=lonc(i)
        lonv(i)=lonp(i)
     end do
  else
     do i=1,im
        lonc(i)=(dble(i)-2.d0)*res+lons
        lonp(i)=(dble(i)-1.5d0)*res+lons !0.5*(lonc(i)+lonc(i+1))
        lonu(i)=lonc(i)
        lonv(i)=lonp(i)
     end do
  end if

  !Latitude
  do j=1,jm
     latc(j)=(dble(j)-2.d0)*res+lats
     latp(j)=(dble(j)-1.5d0)*res+lats !0.5*(latc(j)+latc(j+1))
     latv(j)=latc(j)
     latu(j)=latp(j)
  end do

  do j=1,jm
     dlon(j)=earth*dcos(pi*latp(j)/180.d0)*res*pi/180.d0
  end do

  dlat=earth*res*pi/180.d0

end subroutine make_position

!------------------------------------------------------------------
! Modify topography |
!--------------------
! 1. Sea floor value (minus) --> plus, Land --> 0
! 2. Apply hmin,hmax
!------------------------------------------------------------------

subroutine modify_topography(im,jm,topo)

  use setting
  implicit none

  integer i,j

  real(kind = 8) sum,pass,miss
  
  integer,intent(in) :: im,jm
  real(kind = 8),intent(inout) :: topo(im,jm)

  do j=1,jm
     do i=1,im
        
        !Land --> 0
        if(0.d0 < topo(i,j))then
           topo(i,j)=0.d0
        !Apply hmin,hmax
        else if(hmin < topo(i,j) .and. topo(i,j) < 0.d0)then
           topo(i,j)=hmin
        else if(topo(i,j) < hmax)then
           topo(i,j)=hmax
        end if

     end do
  end do
  
end subroutine modify_topography

!----------------------------------

subroutine fill_topography(im,jm,topo)

  implicit none
  
  integer i,j
  integer count,tcount

  integer,intent(in) :: im,jm
  real(kind = 8),intent(inout) :: topo(im,jm)

  tcount=1
  
  do while(tcount > 0)

     tcount=0

     do j=2,jm-1
        do i=2,im-1
           
           count=0
           
           if(topo(i,j) == 0.d0) cycle
        
           if(topo(i-1,j) == 0.d0) count=count+1
           if(topo(i+1,j) == 0.d0) count=count+1
           if(topo(i,j-1) == 0.d0) count=count+1
           if(topo(i,j+1) == 0.d0) count=count+1
           
           if(count >= 3)then
              topo(i,j)=0.d0
              tcount=tcount+1
           end if
           
        end do
     end do

     write(*,'(a,i10,a,i10)') "Total count (filling point):",tcount,"/",im*jm

  end do

end subroutine fill_topography

!-------------------------------------

subroutine landsea_mask(iswitch,im,jm,lon,lat,topo)

  use setting,only: ngrid, lon_sea, lat_sea
  implicit none

  !---Common
  integer i,j,di,dj
  integer igrid,ipass
  integer idx(ngrid),idy(ngrid)
  integer nland(2),nsea(2)
  
  integer mask(im,jm) !1: Sea, 0: Land
  
  real(kind =8) dx(ngrid),dy(ngrid)

  logical lcycle
  
  !---IN
  integer,intent(in) :: iswitch !1:topo, 2:fsm
  integer,intent(in) :: im,jm
  real(kind = 8),intent(in) :: lon(im),lat(jm)
  
  !---INOUT
  real(kind = 8),intent(inout) :: topo(im,jm)

  idx(:)=0
  idy(:)=0
  dx(:)=999.d0
  dy(:)=999.d0

  do igrid=1,ngrid
     
     do i=1,im
        if(abs(lon(i)-lon_sea(igrid)) <= dx(igrid))then
           dx(igrid)=abs(lon(i)-lon_sea(igrid))
           idx(igrid)=i
        end if
     end do

     do j=1,jm
        if(abs(lat(j)-lat_sea(igrid)) <= dy(igrid))then
           dy(igrid)=abs(lat(j)-lat_sea(igrid))
           idy(igrid)=j
        end if
     end do

  end do

  ipass=0
  do igrid=1,ngrid
     if(idx(igrid) == 0 .or. idy(igrid) == 0) cycle
     ipass=ipass+1
     write(*,'(a,i10,a,i10)') "idx:",idx(igrid),"idy:",idy(igrid)
     write(*,'(a,f12.5,a,f12.5,a,f12.5)') &
          & "Longitude:",lon(idx(igrid))," Latitude:",lat(idy(igrid))," Topo:",topo(idx(igrid),idy(igrid))
  end do

  if(ipass == 0)then
     write(*,'(a)') "***Error: No sea grid points at (lon_sea,lat_sea)"
     stop
  end if
  
  !Set Original sea point
  mask(:,:)=0
  do igrid=1,ngrid
     if(idx(igrid) == 0 .or. idy(igrid) == 0) cycle
     mask(idx(igrid),idy(igrid))=1
  end do

  lcycle=.true.
  do while(lcycle)

     call count_landsea_grid(im,jm,nland(1),nsea(1),mask)
     
     !Increase sea grid
     do j=1,jm
        do i=1,im

           if(mask(i,j) == 1)then              
              do dj=-1,1
                 do di=-1,1
                    if(i+di < 1 .or. im < i+di .or. j+dj < 1 .or. jm < j+dj) cycle
                    if(abs(di) == 1 .and. abs(dj) == 1) cycle
                    if(iswitch == 1 .and. topo(i+di,j+dj) < 0.d0)then
                       mask(i+di,j+dj)=1
                    else if(iswitch == 2 .and. topo(i+di,j+dj) == 1)then
                       mask(i+di,j+dj)=1
                    end if
                 end do
              end do
           end if
           
        end do
     end do

     call count_landsea_grid(im,jm,nland(2),nsea(2),mask)

     write(*,*) "#Land grid:",nland(2),"/",im*jm,"#Sea grid:",nsea(2),"/",im*jm
     
     if(nland(1) == nland(2) .and. nsea(1) == nsea(2))then
        lcycle=.false.
     end if
     
  end do !while

  !Apply mask file
  do j=1,jm
     do i=1,im
        topo(i,j)=topo(i,j)*mask(i,j)
     end do
  end do
  
end subroutine landsea_mask

!----------------------------------

subroutine count_landsea_grid(im,jm,nland,nsea,mask)

  implicit none

  !---Common
  integer i,j

  !---IN
  integer,intent(in) :: im,jm
  integer,intent(in) :: mask(im,jm)
  
  !---INOUT
  integer,intent(inout) :: nland,nsea

  nland=0
  nsea=0

  do j=1,jm
     do i=1,im
        if(mask(i,j) == 0)then
           nland=nland+1
        else if(mask(i,j) == 1)then
           nsea=nsea+1
        end if
     end do
  end do

end subroutine count_landsea_grid

!----------------------------------------------------------------
! Make fsm & Write Ascii |
!-------------------------
!
! fsm: 1 --> OK, 0 --> missing value/Land
!
!----------------------------------------------------------------

subroutine make_write_fsm(im,jm,topo,fsm,rmiss)

  use setting, only: lons,lone
  implicit none

  integer i,j,im,jm
  integer status,access

  real(kind = 8),intent(in) :: topo(im,jm)
  real(kind = 8),intent(in) :: rmiss

  real(kind = 8),intent(out) :: fsm(im,jm)
  
  character(5) cim
  character(20) format

  status=access("fsm.dat"," ")
  if(status == 0)then
     write(*,'(a)') "***Please remove fsm.dat & Execute again"
     stop
  end if

  !Make fsm
  do j=1,jm
     do i=1,im

        if(topo(i,j) == rmiss .or. topo(i,j) == 0.d0)then
           fsm(i,j)=0.d0
        else
           fsm(i,j)=1.d0
        end if

     end do
  end do

  !Strait of Dover
  if(lons == lone-360.d0)then
     fsm(1:2,:)=fsm(im-1:im,:)
  end if
  
  !Write fsm
  write(cim,'(i5)') im
  format="("//trim(cim)//"i1)"
  open(1,file="fsm.dat",status="replace")
  do j=jm,1,-1
     write(1,trim(format)) nint(fsm(1:im,j))
  end do
  close(1)

end subroutine make_write_fsm

!-----------------------------

subroutine read_fsm(im,jm,fsm)

  implicit none
  
  integer i,j
  integer itmp(im,jm)
  character(5) cim
  character(10) format

  integer,intent(in) :: im,jm
  real(kind = 8),intent(out) :: fsm(im,jm)
  
  write(cim,'(i5)') im
  format="("//trim(cim)//"i1)"
  open(1,file="fsm.dat",status="old")
  do j=jm,1,-1
     read(1,trim(format)) itmp(1:im,j)
  end do
  close(1)

  fsm(:,:)=dble(itmp(:,:))

end subroutine read_fsm

!-------------------------------------------------------------------
! Modify topography using fsm |
!-------------------------------------------------------------------

subroutine modify_topography_fsm(im,jm,topo,fsm)

  use setting
  implicit none

  integer i,j
  integer di,dj
  integer icount

  integer,intent(in) :: im,jm
  real(kind = 8),intent(inout) :: topo(im,jm),fsm(im,jm)

  do j=1,jm
     do i=1,im

        ! topo becomes 0 at fsm 1 -> 0
        topo(i,j)=fsm(i,j)*topo(i,j)

     end do
  end do

  do j=1,jm
     do i=1,im
        
        if(fsm(i,j) == 1.d0 .and. topo(i,j) == 0.d0)then
           ! substitute hmin into topo at fsm 0 -> 1 
           !           topo(i,j)=hmin
           ! S.Ohishi 2020.03
           ! average adjacent topo
           icount=0
           topo(i,j)=0.d0
           do dj=-1,1
              do di=-1,1
                 if(i+di < 1 .or. im < i+di)cycle
                 if(j+dj < 1 .or. jm < j+dj)cycle
                 if(topo(i+di,j+dj) == 0.d0)cycle
                 topo(i,j)=topo(i,j)+topo(i+di,j+dj)
                 icount=icount+1
              end do
           end do

           if(icount == 0)then
              topo(i,j)=hmin
           else
              topo(i,j)=topo(i,j)/dble(icount)
           end if

        end if
       
     end do
  end do

end subroutine modify_topography_fsm

!-------------------------------------------------------------------
! Calculate mask |
!-------------------------------------------------------------------

subroutine cal_mask(im,jm,fsm,dum,dvm)

  use setting, only: lons,lone
  implicit none
  integer i,j
  
  integer,intent(in) :: im,jm
  real(kind = 8),intent(in) :: fsm(im,jm)
  real(kind = 8),intent(out) :: dum(im,jm),dvm(im,jm)
  
  dum(:,:)=0.d0
  dvm(:,:)=0.d0

  if(lone-360.d0 == lons)then
     do j=1,jm
        dum(1,j)=fsm(1,j)*fsm(im-2,j)
     end do
  else
     dum(1,:)=fsm(1,:)
  end if
  
  do j=1,jm
     do i=2,im
        dum(i,j)=fsm(i,j)*fsm(i-1,j)
     end do
  end do

  dvm(:,1)=fsm(:,1)
  do j=2,jm
     do i=1,im
        dvm(i,j)=fsm(i,j)*fsm(i,j-1)
     end do
  end do
  
end subroutine cal_mask



!---------------------------------------------------------------------
! Make sigma coordinate(z,zz)
!---------------------------------------------------------------------
!
!   0[m] ----- h1 ----- h2 -----h3 -----h4 -----hmax
! k:1    ----- k1m----- k2m-----k3m-----k4m-----km
!
!---------------------------------------------------------------------

subroutine make_sigma_coordinate(im,jm,topo,fsm,z,zz)

  use setting, only: hmax,km
  implicit none
  
  integer k,k1m,k2m,k3m,k4m
  real(kind = 8) h1,h2,h3,h4
  real(kind = 8) tmp

  integer,intent(in) :: im,jm
  real(kind = 8),intent(in) :: topo(im,jm),fsm(im,jm)
  
  real(kind = 8),intent(out) :: z(im,jm,km),zz(im,jm,km)

  if(km < 75)then
     h1=-100.d0
     k1m=10
     h2=-300.d0
     k2m=20
     h3=-800.d0
     k3m=30
     h4=-2000.d0
     k4m=42     
  else if(75 <= km .and. km < 100)then
     h1=-10.d0
     k1m=10
     h2=-200.d0
     k2m=29
     h3=-500.d0
     k3m=44
     h4=-2000.d0
     k4m=59     
  else if(100 <= km)then
     h1=-10.d0
     k1m=10
     h2=-200.d0
     k2m=48
     h3=-500.d0
     k3m=67
     h4=-2000.d0
     k4m=82     
  else
     write(*,*) "***Error in make_sigma_coorcinate"
     stop
  end if
  
  !Make z
  z(:,:,1)=0.d0                   ! Surface
                              ! (dz=1m)
  z(:,:,k1m)=-1.d0*abs(h1/hmax)   ! MLD at summer: h1=10m at k1m=10
                              ! (dz=5m)
  z(:,:,k2m)=-1.d0*abs(h2/hmax)   ! MLD at winter: h2=200m at k2m=48
                              ! (dz ~ 15m)
  z(:,:,k3m)=-1.d0*abs(h3/hmax)   ! h3=500m  at k3m=67
                              ! (dz=100m)
  z(:,:,k4m)=-1.d0*abs(h4/hmax)   ! h4=2000m at k4m=82
                              ! (dz=250m)
  z(:,:,km)=-1.d0                 ! h5=hmax(6500m) at km=100
     
  !Linear interpolation     
  do k=1,k1m
     call linear_interpolate(dble(1),dble(k1m),z(1,1,1),z(1,1,k1m),dble(k),tmp)
     z(:,:,k)=tmp
  end do

  do k=k1m+1,k2m
     call linear_interpolate(dble(k1m),dble(k2m),z(1,1,k1m),z(1,1,k2m),dble(k),tmp)
     z(:,:,k)=tmp
  end do

  do k=k2m+1,k3m
     call linear_interpolate(dble(k2m),dble(k3m),z(1,1,k2m),z(1,1,k3m),dble(k),tmp)
     z(:,:,k)=tmp
  end do

  do k=k3m+1,k4m
     call linear_interpolate(dble(k3m),dble(k4m),z(1,1,k3m),z(1,1,k4m),dble(k),tmp)
     z(:,:,k)=tmp
  end do

  do k=k4m+1,km
     call linear_interpolate(dble(k4m),dble(km),z(1,1,k4m),z(1,1,km),dble(k),tmp)
     z(:,:,k)=tmp
  end do

  !Make zz
  do k=1,km
     if(k==km)then
        zz(:,:,k)=2.d0*zz(:,:,km-1)-zz(:,:,km-2)
     else
        zz(:,:,k)=0.5d0*(z(:,:,k)+z(:,:,k+1))
     end if
  end do

  write(*,'(a)') "Sigma coordinate:"
  write(*,'(a)') "----z(k)---- ----zz(k)----"
  do k=1,km
     write(*,'(4f12.5)') z(1,1,k),zz(1,1,k),hmax*z(1,1,k),hmax*zz(1,1,k)
  end do
  
end subroutine make_sigma_coordinate

  
!------------------------------------------------------------
! Generalized s-coordinate |
!------------------------------------------------------------

subroutine make_jamstec_scoordinate(im,jm,topo,fsm,z,zz)

  use setting
  implicit none

  !---Parameter
  real(kind = 8),parameter :: c1=1.d0/8.6d0,c2=4.d0 
  
  !---Common
  integer i,j,k,kk,ke

  real(kind = 8) kdz(km) !interval
  real(kind = 8) depth(im,jm,km) !depth at w point 
  real(kind = 8) tmp1dz(km)
  
  !---IN
  integer,intent(in) :: im,jm
  real(kind = 8),intent(in) :: topo(im,jm) !sign: minus
  real(kind = 8),intent(in) :: fsm(im,jm)
  
  !---OUT
  real(kind = 8),intent(out) :: z(im,jm,km),zz(im,jm,km) !Sigma coordinate z=[-1,0]

  !---kdz
  if(km <= 36)then
     kdz(1:km)=abs(hmax)/real(km-1) !Uniform interval
  else if(km < 75)then
     kdz(1:10)=10.d0    !Surface:      10m interval at 0-100 m depth
     kdz(11:30)=25.d0   !Subsurface:   25m interval at 100-600 m
     kdz(31:35)=100.d0  !Subsurface:   100m interval at 600-1000 m
     kdz(36:km-1)=(abs(hmax)-sum(kdz(1:35)))/real(km-1-36)
  else
     kdz(1:10)=1.d0     !Surface:      1m interval at 0-10 m depth
     kdz(11:30)=5.d0    !Subsurface:   5m interval at 10-100 m
     kdz(31:40)=10.d0   !Subsurface:   10m interval at 100-200 m
     kdz(41:55)=20.d0   !Subsurface:   20m interval at 200-500 m
     kdz(46:60)=100.d0  !Intermediate: 100m interval at 500-2000 m
     kdz(61:km-1)=(abs(hmax)-sum(kdz(1:60)))/real(km-1-60)
  end if
  kdz(km)=0.d0

  !Smoothing
  tmp1dz(:)=kdz(:)
  do k=2,km-2
     kdz(k)=0.5d0*(tmp1dz(k-1)+tmp1dz(k+1))
  end do
  
  !---depth
  depth(:,:,:)=0.d0
  do k=2,km
     do j=1,jm
        do i=1,im
           depth(i,j,k)=depth(i,j,k-1)-kdz(k-1)
        end do
     end do
  end do
  
  do j=1,jm
     do i=1,im

        ke=0
        do k=1,km
           if(depth(i,j,k) < botsw*topo(i,j))then
              exit
           else
              ke=k
           end if
        end do

        if(ke == 0)then

           do k=1,km
              depth(i,j,k)=-1.d0*real(k-1)/real(km-1)
           end do
           
        else
           
           do k=ke+1,km
              depth(i,j,k)=depth(i,j,k-1)-(abs(topo(i,j))+depth(i,j,ke))/real(km-ke)
           end do
           
           if(1 < ke)then
              depth(i,j,ke)=0.5*(depth(i,j,ke+1)+depth(i,j,ke-1))
           end if

        end if

     end do !i
  end do !j

  !---z
  do k=1,km
     do j=1,jm
        do i=1,im

           if(fsm(i,j) == 0.d0)then
              z(i,j,k)=-1.d0*real(k-1)/real(km-1)
           else
              z(i,j,k)=-1.d0*depth(i,j,k)/topo(i,j)
           end if
           
        end do
     end do
  end do

  !---zz
  do k=1,km
     do j=1,jm
        do i=1,im
           if(k == km)then
              zz(i,j,k)=2.d0*zz(i,j,km-1)-zz(i,j,km-2)
           else
              zz(i,j,k)=0.5d0*(z(i,j,k)+z(i,j,k+1))
           end if
        end do
     end do
  end do

end subroutine make_jamstec_scoordinate

!---------------------------------------

subroutine make_roms_scoordinate(im,jm,topo,fsm,z,zz)

  use setting
  implicit none

  !---Parameter
  real(kind = 8),parameter :: theta_s=6.5d0 !0-10
  real(kind = 8),parameter :: theta_b=2.5d0 !0-4
  
  !---Common
  integer i,j,k

  real(kind = 8) sigma(km),c(km)
  
  !---IN
  integer,intent(in) :: im,jm
  real(kind = 8),intent(in) :: topo(im,jm) !sign: minus
  real(kind = 8),intent(in) :: fsm(im,jm)
  
  !---OUT
  real(kind = 8),intent(out) :: z(im,jm,km),zz(im,jm,km) !Sigma coordinate z=[-1,0]

  !---sigma
  do k=1,km
     sigma(k)=-1.d0*(k-1)/dble(km-1)
  end do

  !---c_s
  if(0.d0 <= theta_s)then 
     do k=1,km
        c(k)=(1-cosh(theta_s*sigma(k)))/(cosh(theta_s)-1)
     end do   
  else
     do k=1,km
        c(k)=-1.d0*sigma(k)*sigma(k)
     end do   
  end if

  !---c_b
  do k=1,km
     c(k)=(exp(theta_b*c(k))-1.d0)/(1.d0-exp(-theta_b))
  end do

  !---z
  do k=1,km
     do j=1,jm
        do i=1,im
           z(i,j,k)=(abs(hmin)*sigma(k)+abs(topo(i,j))*c(k))/(abs(hmin)+abs(topo(i,j)))
        end do
     end do
  end do

  !---zz
  do k=1,km
     do j=1,jm
        do i=1,im
           if(k == km)then
              zz(i,j,k)=2.d0*zz(i,j,km-1)-zz(i,j,km-2)
           else
              zz(i,j,k)=0.5d0*(z(i,j,k)+z(i,j,k+1))
           end if
        end do
     end do
  end do  
  
end subroutine make_roms_scoordinate

!------------------------------------------------------------
! Write Data |
! -----------------------------------------------------------

subroutine write_data(im,jm,km,z,zz,dlon,dlat,&
     & lonc,lonp,lonu,lonv,latc,latp,latu,latv,topo,topo_dif,ratio,fsm,dum,dvm)

  use netcdf
  implicit none

  integer i,j,k
  integer status,system,ncid,varid

  real(kind = 8) dx(im,jm),dy(im,jm)
  real(kind = 8) east_c(im,jm),east_e(im,jm),east_u(im,jm),east_v(im,jm)
  real(kind = 8) north_c(im,jm),north_e(im,jm),north_u(im,jm),north_v(im,jm)
  real(kind = 8) h(im,jm) !h(=abs(topo))

  character(100) :: filename="../in/grid.nc"
  
  integer,intent(in) :: im,jm,km
  
  real(kind = 8),intent(in) :: z(im,jm,km),zz(im,jm,km)
  real(kind = 8),intent(in) :: dlon(jm),dlat
  real(kind = 8),intent(in) :: lonc(im),lonp(im),lonu(im),lonv(im)
  real(kind = 8),intent(in) :: latc(jm),latp(jm),latu(jm),latv(jm)
  real(kind = 8),intent(in) :: topo(im,jm),topo_dif(im,jm),ratio(im,jm)
  real(kind = 8),intent(in) :: fsm(im,jm),dum(im,jm),dvm(im,jm)

  !Substitute variables(in) into variables(out)
  do j=1,jm
     do i=1,im
        dx(i,j)=dlon(j)
        dy(i,j)=dlat
        east_c(i,j)=lonc(i)
        east_e(i,j)=lonp(i)
        east_u(i,j)=lonu(i)
        east_v(i,j)=lonv(i)
        north_c(i,j)=latc(j)
        north_e(i,j)=latp(j)
        north_u(i,j)=latu(j)
        north_v(i,j)=latv(j)

        if(fsm(i,j) == 0.d0)then
           h(i,j)=1.d0
        else
           h(i,j)=abs(topo(i,j))
        end if

     end do
  end do

  !Make grid.nc
  status=system("rm -f "//trim(filename))
  call make_ncfile_grid(im,jm,km,filename)

  !Make grid.nc  
  !status=system("ncgen -o ../in/grid.nc ../hdr/grid.hdr")

  !Write Data
  status=nf90_open("../in/grid.nc",nf90_write,ncid)

  status=nf90_inq_varid(ncid,"z_w",varid)
  status=nf90_put_var(ncid,varid,z)

  status=nf90_inq_varid(ncid,"z_e",varid)
  status=nf90_put_var(ncid,varid,zz)

  status=nf90_inq_varid(ncid,"dx",varid)
  status=nf90_put_var(ncid,varid,dx)

  status=nf90_inq_varid(ncid,"dy",varid)
  status=nf90_put_var(ncid,varid,dy)

  status=nf90_inq_varid(ncid,"east_u",varid)
  status=nf90_put_var(ncid,varid,east_u)

  status=nf90_inq_varid(ncid,"east_v",varid)
  status=nf90_put_var(ncid,varid,east_v)

  status=nf90_inq_varid(ncid,"east_e",varid)
  status=nf90_put_var(ncid,varid,east_e)

  status=nf90_inq_varid(ncid,"east_c",varid)
  status=nf90_put_var(ncid,varid,east_c)

  status=nf90_inq_varid(ncid,"north_u",varid)
  status=nf90_put_var(ncid,varid,north_u)

  status=nf90_inq_varid(ncid,"north_v",varid)
  status=nf90_put_var(ncid,varid,north_v)

  status=nf90_inq_varid(ncid,"north_e",varid)
  status=nf90_put_var(ncid,varid,north_e)

  status=nf90_inq_varid(ncid,"north_c",varid)
  status=nf90_put_var(ncid,varid,north_c)

  status=nf90_inq_varid(ncid,"h",varid)
  status=nf90_put_var(ncid,varid,h)

  status=nf90_inq_varid(ncid,"hdif",varid)
  status=nf90_put_var(ncid,varid,topo_dif)

  status=nf90_inq_varid(ncid,"ratio",varid)
  status=nf90_put_var(ncid,varid,ratio)
  
  status=nf90_inq_varid(ncid,"fsm",varid)
  status=nf90_put_var(ncid,varid,fsm)

  status=nf90_inq_varid(ncid,"dum",varid)
  status=nf90_put_var(ncid,varid,dum)

  status=nf90_inq_varid(ncid,"dvm",varid)
  status=nf90_put_var(ncid,varid,dvm)

  status=nf90_close(ncid)
  
end subroutine write_data

!------------------------------------------------------------------
! Time step
!------------------------------------------------------------------

subroutine time_step(jm,dlon,hmax)

  implicit none

  integer,intent(in) :: jm
  real(kind = 8),intent(in) :: dlon(jm),hmax

  real(kind = 8) ce,t

  ce=2.d0*dsqrt(9.81d0*dabs(hmax))
  
  t=minval(dlon)/ce

  write(*,'(a,i6)') "Mamimum external time step [s]:",nint(t)

end subroutine time_step

!------------------------------------------------------------------
! Make mod_gridinfo.f90 |
!------------------------------------------------------------------

subroutine make_mod_gridinfo(im,jm,km)

  implicit none

  integer,intent(in) :: im,jm,km
  
  open(1,file="src/mod_gridinfo.f90",status="replace")
  write(1,'(a)') "module mod_gridinfo"
  write(1,*)
  write(1,'(a,i6)') "   integer,parameter :: im=",im
  write(1,'(a,i6)') "   integer,parameter :: jm=",jm
  write(1,'(a,i6)') "   integer,parameter :: km=",km
  write(1,*) 
  write(1,'(a)') "end module mod_gridinfo"
  close(1)

end subroutine make_mod_gridinfo
