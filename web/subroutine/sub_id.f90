!---------------------------------------------------------------
! Calculate ID (2D)|
!-------------------
!
! Detect idx and idy at southwestern edge
! Regular lon-lat grid version
!
!---------------------------------------------------------------

subroutine cal_idxy(im1,jm1,lon1,lat1, &
     & im2,jm2,lon2,lat2,idx,idy)

  implicit none

  !---Common
  integer i1,j1
  integer i2,j2

  integer,intent(in) :: im1,jm1
  integer,intent(in) :: im2,jm2
  
  !---IN
  real(kind = 8),intent(in) :: lon1(im1),lat1(jm1)
  real(kind = 8),intent(in) :: lon2(im2),lat2(jm2)

  !---OUT
  integer,intent(out) :: idx(im2),idy(jm2)

  !---Initalization
  idx(:)=0
  idy(:)=0

  !---Regular lon-lat grid
  do i2=1,im2
     do i1=1,im1-1
        if(lon1(i1) <= lon2(i2) .and. lon2(i2) < lon1(i1+1))then
           idx(i2)=i1
           exit
        end if
     end do
  end do

  do j2=1,jm2
     do j1=1,jm1-1
        if(lat1(j1) <= lat2(j2) .and. lat2(j2) < lat1(j1+1))then
           idy(j2)=j1
           exit
        end if
     end do
  end do
  
end subroutine cal_idxy

!---------------------------------------------------------------
! Calculate ID (3D) |
!--------------------
!
! For sigma coordinate, calculate id.
! When id cannot be detected, id is 0.
!
!---------------------------------------------------------------

subroutine cal_idz(im1,jm1,km1,dat1,km2,dat2,idz)

  implicit none

  !---Common
  integer i1,j1,k1
  integer k2

  real(kind = 8) z1,z2
  
  !---IN
  integer,intent(in) :: im1,jm1,km1
  integer,intent(in) :: km2
  real(kind = 8),intent(in) :: dat1(im1,jm1,km1),dat2(km2)

  !---OUT
  integer,intent(out) :: idz(im1,jm1,km2)

  !---Initialization
  idz(:,:,:)=0

  !---idz
  do k1=1,km1-2
     do j1=1,jm1
        do i1=1,im1

           z1=dat1(i1,j1,k1)
           z2=dat1(i1,j1,k1+1)
           
           do k2=2,km2
              if(z1 <= dat2(k2) .and. dat2(k2) < z2)then
                 idz(i1,j1,k2)=k1
              end if
           end do
           
        end do
     end do
  end do
  
end subroutine cal_idz
