subroutine convert_to_obs_space(km_a,lon_a,lat_a,dep_a,dat_a, &
     & km_o,lon_o,lat_o,dep_o,hdat_a)

  use mod_rmiss
  implicit none

  !---Parameter
  integer,parameter :: im_a=2,jm_a=2

  !---Common
  integer i,j,k_o,k_a

  real(kind = 8) tmp(im_a,jm_a,km_o)
  
  !---IN
  !Analysis
  integer,intent(in) :: km_a
  real(kind = 8),intent(in) :: lon_a(im_a),lat_a(jm_a),dep_a(im_a,jm_a,km_a)
  real(kind = 8),intent(in) :: dat_a(im_a,jm_a,km_a)

  !Observation
  integer,intent(in) :: km_o
  real(kind = 8),intent(in) :: lon_o,lat_o,dep_o(km_o)

  !---OUT
  real(kind = 8),intent(out) :: hdat_a(km_o)
  
  !---Linear interpolation in vertical direction
  tmp(:,:,:)=rmiss
  
  do k_o=1,km_o
     do k_a=1,km_a-2
        do j=1,jm_a
           do i=1,im_a
              if(dep_o(k_o) < dep_a(i,j,k_a) .or. dep_a(i,j,k_a+1) < dep_o(k_o)) cycle  
              if(dat_a(i,j,k_a) == rmiss .or. dat_a(i,j,k_a+1) == rmiss) cycle
              call linear_interpolate(dep_a(i,j,k_a),dep_a(i,j,k_a+1),dat_a(i,j,k_a),dat_a(i,j,k_a+1), &
                   & dep_o(k_o),tmp(i,j,k_o))
           end do
        end do
     end do
  end do
  
  !---Bilinear interpolation in horizontal direction
  hdat_a(:)=rmiss
  
  do k_o=1,km_o
     if(tmp(1,1,k_o) == rmiss .or. tmp(2,1,k_o) == rmiss .or. tmp(1,2,k_o) == rmiss .or. tmp(2,2,k_o) == rmiss) cycle
     call bilinear_interpolation(lon_a(1),lon_a(2),lat_a(1),lat_a(2), &
          & lon_o,lat_o,tmp(1,1,k_o),tmp(2,1,k_o),tmp(1,2,k_o),tmp(2,2,k_o),hdat_a(k_o))
  end do

end subroutine convert_to_obs_space

!-----------------------------------------------------------
! Linear interpolation |
!-----------------------------------------------------------

subroutine linear_interpolate(x1,x2,y1,y2,x,y)

  implicit none

  real(kind = 8),intent(in) :: x1,x2,y1,y2
  real(kind = 8),intent(in) :: x
  real(kind = 8),intent(out) :: y

  y=(y2-y1)/(x2-x1)*(x-x1)+y1

end subroutine linear_interpolate

!-----------------------------------------------------------
! Bilinear Interpolation |
!-----------------------------------------------------------
!
! f01(x0,y1) ___________ f11(x1,y1)
!  |                       |
! y|-----------fxy(x,y)    |
!  |            |          |
!  |            |          |
! f00(x0,y0) ___|_______ f10(x1,y0)
!               x
!-----------------------------------------------------------

SUBROUTINE bilinear_interpolation(x0,x1,y0,y1,x,y,f00,f10,f01,f11,fxy)

  IMPLICIT NONE

  !---IN
  REAL(kind = 8),INTENT(IN)  :: x0,x1,y0,y1,x,y
  REAL(kind = 8),INTENT(IN)  :: f00,f10,f01,f11

  !---OUT
  REAL(kind = 8),INTENT(OUT) :: fxy

  fxy=(y1-y)/(y1-y0)*((x1-x)/(x1-x0)*f00+(x-x0)/(x1-x0)*f10) &
       & +(y-y0)/(y1-y0)*((x1-x)/(x1-x0)*f01+(x-x0)/(x1-x0)*f11)

END SUBROUTINE bilinear_interpolation

