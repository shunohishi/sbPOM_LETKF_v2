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
     & (im_a,jm_a,lon_a,lat_a,dat_a,mask_a, &
     &  n_o,lon_o,lat_o, &
     &  idlon,idlat,hdat_a)

  !$use omp_lib
  use mod_rmiss
  implicit none

  !---Common
  integer i_o
  integer m_o
  integer i_a0,i_a1
  integer j_a0,j_a1
  integer m_a,m_a0,m_a1
  
  !---IN
  integer,intent(in) :: im_a,jm_a
  integer,intent(in) :: n_o
  integer,intent(in) :: idlon(n_o),idlat(n_o)
  
  real(kind = 8),intent(in) :: lon_a(im_a),lat_a(jm_a),dat_a(im_a,jm_a),mask_a(im_a,jm_a)
  real(kind = 8),intent(in) :: lon_o(n_o),lat_o(n_o)

  !---OUT
  real(kind = 8),intent(out) :: hdat_a(n_o)

  !$omp parallel
  !$omp do private(i_o,m_o,i_a0,i_a1,m_a,m_a0,m_a1,j_a0,j_a1)
  do i_o=1,n_o

     if(idlon(i_o) == 0 .or. idlat(i_o) == 0)then
        hdat_a(i_o)=rmiss
        cycle
     end if
     
     m_o=floor(lon_o(i_o)/360.d0)
        
     i_a0=idlon(i_o)
     if(i_a0 == im_a)then
        i_a1=1
     else
        i_a1=idlon(i_o)+1
     end if
        
     m_a0=floor(lon_a(i_a0)/360.d0)
     m_a1=floor(lon_a(i_a1)/360.d0)
     
     j_a0=idlat(i_o)
     j_a1=idlat(i_o)+1
        
     if(mask_a(i_a0,j_a0)==0.d0 .or. mask_a(i_a1,j_a0)==0.d0 &
          & .or. mask_a(i_a0,j_a1)==0.d0 .or. mask_a(i_a1,j_a1)==0.d0)then
        hdat_a(i_o)=rmiss
     else if(i_a0 == im_a)then
        call bilinear_interpolate &
             & (lon_a(i_a0)-m_a0*360.d0,lon_a(i_a1)-(m_a1-1)*360.d0,lat_a(j_a0),lat_a(j_a1),&
             &  lon_o(i_o)-(m_o-1)*360.d0,lat_o(i_o), &
             &  dat_a(i_a0,j_a0),dat_a(i_a1,j_a0),dat_a(i_a0,j_a1),dat_a(i_a1,j_a1), &
             &  hdat_a(i_o))
        write(*,*) lon_a(i_a0)-m_a0*360.d0,lon_a(i_a1)-(m_a1-1)*360.d0
        write(*,*) lon_o(i_o)-(m_o-1)*360.d0,lat_o(i_o)
        stop        
     else
        call bilinear_interpolate &
             & (lon_a(i_a0)-m_a0*360.d0,lon_a(i_a1)-m_a1*360.d0,lat_a(j_a0),lat_a(j_a1),&
             &  lon_o(i_o)-m_o*360.d0,lat_o(i_o), &
             &  dat_a(i_a0,j_a0),dat_a(i_a1,j_a0),dat_a(i_a0,j_a1),dat_a(i_a1,j_a1), &
             &  hdat_a(i_o))        
     end if

  end do !i_o
  !$omp end do
  !$omp end parallel
  
end subroutine bilinear_interpolation_2d

!-----------------------------------------------------------------------------
