!------------------------------------------------------------------
! Detect nearshore grid |
!------------------------------------------------------------------

subroutine detect_nearshore(lon,lat,fsm,nshore)

  use mod_gridinfo
  use setting, only : range => nearshore_range
  use mod_parameter
  implicit none
  
  integer i,j,i1,j1,i2,j2
  integer idlon(jm),idlat
  integer itmp

  real(kind = 8) r
  real(kind = 8),intent(in) :: lon(im),lat(jm),fsm(im,jm)
  real(kind = 8),intent(out) :: nshore(im,jm) !1: Ocean, 0: Land or Nearshore
  
  !---Estimate calculating idlon,idlat
  idlat=nint(range*180.d0/(pi*earth*abs(lat(2)-lat(1))))

  do j1=1,jm

     idlon(j1)=0

     do j2=j1-idlat,j1+idlat

        if(j2 < 1 .or. jm < j2)cycle
           
        itmp=nint(range*180.d0/(pi*earth*abs(lon(2)-lon(1))*cos(lat(j2)*pi/180.d0)))
        if(idlon(j1) < itmp)then
           idlon(j1)=itmp
        end if

     end do
  end do

  write(*,'(a,i6)') "idlat: ",idlat
  write(*,'(a,i6)') "idlon: ",maxval(idlon)
  
  !---Detect nearshore grid
  nshore(:,:)=1.d0
  do j1=1,jm
     
     if(mod(j1,100) == 1)then
        write(*,'(i4.4,a1,i4.4,a)') j1,"/",jm," in Nearshore"
     end if
     
     do i1=1,im
                
        if(fsm(i1,j1) == 0.d0)then
           nshore(i1,j1)=0.d0
           cycle
        end if
        
        do j2=j1-idlat,j1+idlat
           
           if(j2 < 1 .or. jm < j2)cycle
           if(nshore(i1,j1) == 0.d0) cycle

           do i2=i1-idlon(j2),i1+idlon(j2)
              
              if(i2 < 1 .or. im < i2)cycle
              if(nshore(i1,j1) == 0.d0) cycle

              !Calculate distance between (lon(i1),lat(j1)) and (lon(i2),lat(j2))
              call distance(lon(i1),lat(j1),lon(i2),lat(j2),r)
              
              if(r <= range .and. fsm(i2,j2) == 0.d0)then
                 nshore(i1,j1)=0.d0
                 cycle
              end if
              
           end do !i2
        end do !j2
        
     end do !i1
  end do !j1

end subroutine detect_nearshore

!----------------------------------

subroutine detect_nearshore_smap(lon,lat,fsm, &
     & im_smap,jm_smap,lon_smap,lat_smap,nshore_smap)

  use mod_gridinfo
  use setting, only : range => nearshore_range
  use mod_parameter
  implicit none
  
  integer i,j,j1,j2
  integer i_smap,j_smap
  integer dj

  real(kind = 8) r
  real(kind = 8) lon_smap_min,lon_smap_max
  real(kind = 8) lat_smap_min,lat_smap_max

  !---IN
  !Model
  real(kind = 8),intent(in) :: lon(im),lat(jm),fsm(im,jm)

  !SMAP
  integer,intent(in) :: im_smap,jm_smap
  real(kind = 8),intent(in) :: lon_smap(im_smap,jm_smap),lat_smap(im_smap,jm_smap)
  
  !---OUT
  real(kind = 8),intent(out) :: nshore_smap(im_smap,jm_smap) !1: Ocean, 0: Land or Nearshore

  !---Min, Max longitude latitude for SMAP
  lon_smap_min=minval(lon_smap)
  lon_smap_max=maxval(lon_smap)
  lat_smap_min=minval(lat_smap)
  lat_smap_max=maxval(lat_smap)

  !---dj
  dj=nint(range*180.d0/(pi*earth*abs(lat(2)-lat(1))))

  !---j1,j2
  j1=1
  do j=1,jm
     if(lat_smap_min <= lat(j))then
        j1=j
        exit
     end if
  end do

  j2=jm
  do j=jm,1,-1
     if(lat(j) <= lat_smap_max)then
        j2=j
        exit
     end if
  end do
  
  !---Detect nearshore grid
  nshore_smap(:,:)=1.d0
  do j_smap=1,jm_smap
     do i_smap=1,im_smap

        do j=j1-dj,j2+dj
           if(j < 1 .or. jm < j) cycle
           do i=1,im

              if(fsm(i,j) == 1.d0) cycle

              !Calculate distance between (lon(i1),lat(j1)) and (lon(i2),lat(j2))
              call distance(lon(i),lat(j),lon_smap(i_smap,j_smap),lat_smap(i_smap,j_smap),r)

              if(r <= range)then
                 nshore_smap(i_smap,j_smap)=0.d0
                 exit
              end if
              
           end do !i

           if(nshore_smap(i_smap,j_smap) == 0.d0)then
              exit
           end if
           
        end do    !j
        
     end do !i_smap
  end do    !j_smap
  
end subroutine detect_nearshore_smap

!-------------------------------------

subroutine detect_nearshore_smos(lon,lat,fsm, &
     & n_smos,lon_smos,lat_smos,nshore_smos)

  use mod_gridinfo
  use setting, only : range => nearshore_range
  use mod_parameter
  implicit none
  
  integer i,j,j1,j2
  integer i_smos
  integer dj

  real(kind = 8) r
  real(kind = 8) lon_smos_min,lon_smos_max
  real(kind = 8) lat_smos_min,lat_smos_max

  !---IN
  !Model
  real(kind = 8),intent(in) :: lon(im),lat(jm),fsm(im,jm)

  !SMOS
  integer,intent(in) :: n_smos
  real(kind = 8),intent(in) :: lon_smos(n_smos),lat_smos(n_smos)
  
  !---OUT
  real(kind = 8),intent(out) :: nshore_smos(n_smos) !1: Ocean, 0: Land or Nearshore

  !---Min, Max longitude latitude for SMOS
  lon_smos_min=minval(lon_smos)
  lon_smos_max=maxval(lon_smos)
  lat_smos_min=minval(lat_smos)
  lat_smos_max=maxval(lat_smos)

  !---dj
  dj=nint(range*180.d0/(pi*earth*abs(lat(2)-lat(1))))

  !---j1,j2
  j1=1
  do j=1,jm
     if(lat_smos_min <= lat(j))then
        j1=j
        exit
     end if
  end do

  j2=jm
  do j=jm,1,-1
     if(lat(j) <= lat_smos_max)then
        j2=j
        exit
     end if
  end do
  
  !---Detect nearshore grid
  nshore_smos(:)=1.d0
  do i_smos=1,n_smos

     do j=j1-dj,j2+dj
        if(j < 1 .or. jm < j) cycle
        do i=1,im

           if(fsm(i,j) == 1.d0) cycle

           !Calculate distance between (lon(i1),lat(j1)) and (lon(i2),lat(j2))
           call distance(lon(i),lat(j),lon_smos(i_smos),lat(i_smos),r)

           if(r <= range)then
              nshore_smos(i_smos)=0.d0
              exit
           end if

        end do !i

        if(nshore_smos(i_smos) == 0.d0)then
           exit
        end if

     end do    !j

  end do !i_smos
  
end subroutine detect_nearshore_smos
