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

  real(kind = 8) dx1
  
  !---IN
  integer,intent(in) :: im1,im2
  real(kind = 8),intent(in) :: lon1(im1),lon2(im2)

  !---OUT
  integer,intent(out) :: id(im2)

  id(:)=0
  dx1=lon1(2)-lon1(1)

  !$omp parallel
  !$omp do private(i1,i2,n1,n2)
  do i2=1,im2

     n1=floor(lon1(1)/360.d0)
     n2=floor(lon2(i2)/360.d0)
        
     do i1=1,im1-1
        n1=floor(lon1(i1)/360.d0)
        if(lon1(i1)-n1*360.d0 <= lon2(i2)-n2*360.d0 .and. lon2(i2)-n2*360.d0 <= lon1(i1+1)-n1*360.d0)then
           id(i2)=i1
        end if
     end do !i1
        
  end do !i2
  !$omp end do
  !$omp end parallel
  
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
    
end subroutine cal_idlat
