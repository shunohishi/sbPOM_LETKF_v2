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

subroutine bilinear_interpolate(x0,x1,y0,y1,x,y,f00,f10,f01,f11,fxy)

  implicit none

  real(kind = 8),intent(in) :: x0,x1,y0,y1,x,y
  real(kind = 8),intent(in) :: f00,f10,f01,f11
  real(kind = 8),intent(out) :: fxy

  fxy=(y1-y)/(y1-y0)*((x1-x)/(x1-x0)*f00+(x-x0)/(x1-x0)*f10) &
       & +(y-y0)/(y1-y0)*((x1-x)/(x1-x0)*f01+(x-x0)/(x1-x0)*f11)

end subroutine bilinear_interpolate

!________________________________________________________________________

subroutine bilinear_interpolation_2d &
     & (im1,jm1,lon1,lat1,dat1,im2,jm2,lon2,lat2,dat2,idlon,idlat,rmiss)

  !$use omp_lib
  implicit none

  !---Common
  integer i10,i11,j10,j11
  integer i2,j2
  integer n1,n10,n11,n2
  
  !---IN
  integer,intent(in) :: im1,jm1
  integer,intent(in) :: im2,jm2
  integer,intent(in) :: idlon(im2),idlat(jm2)
  
  real(kind = 8),intent(in) :: lon1(im1),lat1(jm1),dat1(im1,jm1)
  real(kind = 8),intent(in) :: lon2(im2),lat2(jm2)
  real(kind = 8),intent(in) :: rmiss

  !---OUT
  real(kind = 8),intent(out) :: dat2(im2,jm2)

  !$omp parallel
  !$omp do private(i2,j2,n2,i10,i11,n1,n10,n11,j10,j11)    
  do j2=1,jm2
     do i2=1,im2

        n2=floor(lon2(i2)/360.d0)
        
        if(idlon(i2) == im1)then

           i10=im1
           i11=1
           n1=floor(lon1(1)/360.d0)
           
           !ex. -0.5 < lon2 < 0.5
           if(lon1(im1)-(n1+1)*360.d0 <= lon2(i2)-n2*360.d0 .and. lon2(i2)-n2*360.d0 <= lon1(1)-n1*360.d0)then
              n10=floor(lon1(im1)/360.d0)-1
              n11=floor(lon1(1)/360.d0)
           !ex.359.5 < lon2 < 360.5
           else if(lon1(im1)-n1*360.d0 <= lon2(i2)-n2*360.d0 .and. lon2(i2)-n2*360.d0 <= lon1(1)-(n1-1)*360.d0)then
              n10=floor(lon1(im1)/360.d0)
              n11=floor(lon1(1)/360.d0)+1
           end if
           
        else
           
           i10=idlon(i2)
           i11=idlon(i2)+1
           n10=floor(lon1(i10)/360.d0)
           n11=floor(lon1(i11)/360.d0)

        end if
        
        j10=idlat(j2)
        j11=idlat(j2)+1
        
        if(dat1(i10,j10)==rmiss .or. dat1(i11,j10)==rmiss &
             & .or. dat1(i10,j11)==rmiss .or. dat1(i11,j11)==rmiss)then
           dat2(i2,j2)=rmiss
        else
           call bilinear_interpolate &
                & (lon1(i10)-n10*360.d0,lon1(i11)-n11*360.d0,lat1(j10),lat1(j11),&
                & lon2(i2)-n2*360.d0,lat2(j2), &
                & dat1(i10,j10),dat1(i11,j10),dat1(i10,j11),dat1(i11,j11), &
                & dat2(i2,j2))
        end if

     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine bilinear_interpolation_2d

!-----------------------------------------------------------------------------
