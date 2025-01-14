!---------------------------------------------------------------
! Calculate ID (longitude)|
!--------------------------
!
! 1  ----i2----2
! i1 ----i2----i1+1
! im1----i2----1
! 
! * id denotes i1, when dat2 is at i2
!
!---------------------------------------------------------------

subroutine cal_idlon(im1,lon1,im2,lon2,id)

  !$use omp_lib
  implicit none

  !---Common
  integer i1,i2
  integer n1,n2
  
  !---IN
  integer,intent(in) :: im1,im2
  real(kind = 8),intent(in) :: lon1(im1),lon2(im2)

  !---OUT
  integer,intent(out) :: id(im2)

  id(:)=0

  !$omp parallel
  !$omp do private(i1,i2,n1,n2)
  do i2=1,im2

     n1=floor(lon1(1)/360.d0)
     n2=floor(lon2(i2)/360.d0)

     !ex. -0.5 < lon2 < 0.5
     if(lon1(im1)-(n1+1)*360.d0 <= lon2(i2)-n2*360.d0 .and. lon2(i2)-n2*360.d0  <= lon1(1)-n1*360.d0)then
        id(i2)=im1
     !ex. 359.5 < lon2 < 360.5
     else if(lon1(im1)-n1*360.d0 <= lon2(i2)-n2*360.d0 .and. lon2(i2)-n2*360.d0  <= lon1(1)-(n1-1)*360.d0)then
        id(i2)=im1
     else
        
        do i1=1,im1-1
           n1=floor(lon1(i1)/360.d0)
           if(lon1(i1)-n1*360.d0 <= lon2(i2)-n2*360.d0 .and. lon2(i2)-n2*360.d0 <= lon1(i1+1)-n1*360.d0)then
              id(i2)=i1
           end if
        end do !i1
        
     end if
  end do !i2
  !$omp end do
  !$omp end parallel
  
  do i2=1,im2
     if(id(i2) == 0)then
        write(*,*) "***Error: idx is missed at ",i2
        stop
     end if
  end do

end subroutine cal_idlon

!------------------------

subroutine cal_idlat(jm1,lat1,jm2,lat2,id)

  !$use omp_lib  
  implicit none

  !---Common
  integer j1,j2
  
  !---IN
  integer,intent(in) :: jm1,jm2
  real(kind = 8),intent(in) :: lat1(jm1),lat2(jm2)

  !---OUT
  integer,intent(out) :: id(jm2)

  id(:)=0

  !$omp parallel
  !$omp do private(j1,j2)  
  do j2=1,jm2
     do j1=1,jm1-1

        if(lat1(j1) <= lat2(j2) .and. lat2(j2) <= lat1(j1+1))then
           id(j2)=j1
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  do j2=1,jm2
     if(id(j2) == 0)then
        write(*,*) "***Error: idy is missed at ",j2
        stop
     end if
  end do
  
end subroutine cal_idlat

!---------------------------------------------------------------
! Calculate ID (3D) |
!--------------------
!
! For sigma coordinate, calculate id.
! When id cannot be detected, id is 0.
!
!---------------------------------------------------------------

subroutine cal_idz(km1,dat1,im2,jm2,km2,dat2,id)

  !$use omp_lib  
  implicit none

  integer i,j,k1,k2

  integer,intent(in) :: km1
  integer,intent(in) :: im2,jm2,km2

  real(kind = 8),intent(in) :: dat1(km1)
  real(kind = 8),intent(in) :: dat2(im2,jm2,km2)

  integer,intent(out) :: id(im2,jm2,km2)
  
  id(:,:,:)=0

  !$omp parallel
  !$omp do private(i,j,k1,k2)  
  do k2=1,km2
     do k1=1,km1-1
        do j=1,jm2
           do i=1,im2
              if(dat1(k1) <= dat2(i,j,k2) .and. dat2(i,j,k2) <= dat1(k1+1))then
                 id(i,j,k2)=k1
              end if
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine cal_idz

!--------------------------------------------------------------
! Calculate ID (3D) |
!--------------------------------------------------------------
!
! Sigma coordinate vs. sigma coordinate
!
!--------------------------------------------------------------

subroutine cal_id_sigma(im,jm,km1,dat1,km2,dat2,id)

  implicit none
  
  integer i,j,k1,k2

  !IN
  integer,intent(in) :: im,jm,km1,km2
  real(kind = 8),intent(in) :: dat1(im,jm,km1),dat2(im,jm,km2)

  !OUT
  integer,intent(out) :: id(im,jm,km2)

  
  id(:,:,:)=0

  do j=1,jm
     do i=1,im

        do k2=1,km2
           do k1=1,km1-1
              if(dat1(i,j,k1) <= dat2(i,j,k2) &
                   & .and. dat2(i,j,k2) <= dat1(i,j,k1+1))then
                 id(i,j,k2)=k1
              end if
           end do !k1
        end do !k2

     end do !i
  end do !j
  
end subroutine cal_id_sigma

!--------------------------------------------------------------
! Calculate CaMa-Flood ID at model coast line |
!--------------------------------------------------------------
!
!--- INPUT
! im,jm,lon,lat: model grid
! cline: model coast line
! 
! *_cama: CaMa-Flood grid
!
!--- OUTPUT
! idx,idy: ID
! dnum: duplication number
!
!--------------------------------------------------------------

subroutine cal_id_cama(im,jm,lon,lat,cline, &
     & im_cama,jm_cama,lon_cama,lat_cama,land_cama, &
     & idx,idy,dnum)

  use mod_rmiss
  implicit none

  real(kind = 8),parameter :: r_max=50.d3

  integer i1,i2,di,j1,j2,dj
  real(kind = 8) r,r0 !distance between 2 point

  !IN
  integer,intent(in) :: im,jm  
  real(kind = 8),intent(in) :: lon(im),lat(jm),cline(im,jm)

  integer,intent(in) :: im_cama,jm_cama
  real(kind = 8),intent(in) :: lon_cama(im_cama),lat_cama(jm_cama)
  real(kind = 8),intent(in) :: land_cama(im_cama,jm_cama)

  !OUT
  integer,intent(out) :: idx(im,jm),idy(im,jm) !river position
  integer,intent(out) :: dnum(im,jm)           !duplication number
  
  idx(:,:)=0
  idy(:,:)=0

  !grid number for minimum distance 
  do j1=1,jm
     do j2=1,jm_cama
        
        if(lat(j1) < lat_cama(j2)-0.5d0  .or. lat_cama(j2)+0.5d0 < lat(j1))cycle
        
        do i1=1,im
           if(cline(i1,j1) == 0.d0)cycle
           do i2=1,im_cama
              if(lon(i1) < lon_cama(i2)-0.5d0  .or. lon_cama(i2)+0.5d0 < lon(i1))cycle
              
              r0=r_max
              
              if(land_cama(i2,j2) == 0.d0)cycle
              
              call distance(lon(i1),lat(j1),lon_cama(i2),lat_cama(j2),r)
              
              if(r < r0)then
                 r0 = r
                 idx(i1,j1)=i2
                 idy(i1,j1)=j2
              end if
                             
           end do !i2
        end do !i1
     end do !j2
  end do !j1

  di=nint(r_max*1.d-5/(lon(2)-lon(1)))
  dj=nint(r_max*1.d-5/(lat(2)-lat(1)))
  
  !Check number
  do j1=1,jm
     do i1=1,im
        
        dnum(i1,j1)=0
        
        if(cline(i1,j1)==0.d0)cycle

        do j2=j1-dj,j1+dj
           if(j2 < 1 .or. jm < j2)cycle
           do i2=i1-di,i1+di
              if(i2 < 1 .or. im < i2)cycle
              if(idx(i1,j1) == idx(i2,j2) .and. idy(i1,j1) == idy(i2,j2))then
                 dnum(i1,j1)=dnum(i1,j1)+1
              end if
           end do !i2
        end do !j2
     end do !i1
  end do !j1
  
end subroutine cal_id_cama
