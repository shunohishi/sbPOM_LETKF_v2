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
     & (im_a,jm_a,lon_a,lat_a,dat_a,mask_a, &
     &  im_o,jm_o,lon_o,lat_o, &
     &  idlon,idlat,hdat_a)

  !$use omp_lib
  use mod_rmiss
  implicit none

  !---Common
  integer i_o,j_o
  integer m_o
  integer i_a0,i_a1
  integer j_a0,j_a1
  integer m_a,m_a0,m_a1
  
  !---IN
  integer,intent(in) :: im_a,jm_a
  integer,intent(in) :: im_o,jm_o
  integer,intent(in) :: idlon(im_o),idlat(jm_o)
  
  real(kind = 8),intent(in) :: lon_a(im_a),lat_a(jm_a),dat_a(im_a,jm_a),mask_a(im_a,jm_a)
  real(kind = 8),intent(in) :: lon_o(im_o),lat_o(jm_o)

  !---OUT
  real(kind = 8),intent(out) :: hdat_a(im_o,jm_o)

  !$omp parallel
  !$omp do private(i_o,m_o,j_o,i_a0,i_a1,m_a,m_a0,m_a1,j_a0,j_a1)
  do j_o=1,jm_o
     do i_o=1,im_o

        if(idlon(i_o) == 0 .or. idlat(j_o) == 0)then
           hdat_a(i_o,j_o)=rmiss
           cycle
        end if

        m_o=floor(lon_o(i_o)/360.d0)

        i_a0=idlon(i_o)
        i_a1=idlon(i_o)+1

        m_a0=floor(lon_a(i_a0)/360.d0)
        m_a1=floor(lon_a(i_a1)/360.d0)

        j_a0=idlat(j_o)
        j_a1=idlat(j_o)+1

        if(mask_a(i_a0,j_a0)==0.d0 .or. mask_a(i_a1,j_a0)==0.d0 &
             & .or. mask_a(i_a0,j_a1)==0.d0 .or. mask_a(i_a1,j_a1)==0.d0)then
           hdat_a(i_o,j_o)=rmiss
        else
           call bilinear_interpolate &
                & (lon_a(i_a0)-m_a0*360.d0,lon_a(i_a1)-m_a1*360.d0,lat_a(j_a0),lat_a(j_a1),&
                &  lon_o(i_o)-m_o*360.d0,lat_o(j_o), &
                &  dat_a(i_a0,j_a0),dat_a(i_a1,j_a0),dat_a(i_a0,j_a1),dat_a(i_a1,j_a1), &
                &  hdat_a(i_o,j_o))

        end if
        
     end do
  end do !i_o
  !$omp end do
  !$omp end parallel
  
end subroutine bilinear_interpolation_2d

!-----------------------------------------------------------------------------

subroutine vertical_interpolation(km_a,dep_a,dat_a,dep_o,hdat_a)

  use mod_rmiss
  implicit none

  !---Common
  integer k
  integer id
  
  !---IN
  !Analysis
  integer,intent(in) :: km_a

  real(kind = 8),intent(in) :: dep_a(km_a),dat_a(km_a)  

  !Obs
  real(kind = 8),intent(in) :: dep_o

  !---OUT
  real(kind = 8),intent(out) :: hdat_a

  !---Check surface depth
  if(dep_o <= dep_a(1))then
     hdat_a=dat_a(1)
     return
  end if
  
  !---Get ID
  id=0
  
  do k=1,km_a-1
     if(dep_o < dep_a(k) .or. dep_a(k+1) <= dep_o) cycle
     id=k
  end do

  !---Check data
  if(id == 0)then
     hdat_a=rmiss
     return
  else if(dat_a(id) == rmiss .or. dat_a(id+1) == rmiss)then
     hdat_a=rmiss
     return
  end if

  !---Linear interpolation
  call linear_interpolate(dep_a(id),dep_a(id+1),dat_a(id),dat_a(id+1),dep_o,hdat_a)
  
end subroutine vertical_interpolation

!-----------------------------------------------------------------------------

subroutine vertical_bilinear_interpolate(im_a,jm_a,km_a,lon_a,lat_a,dep_a,dat_a, &
             & im_o,jm_o,km_o,lon_o,lat_o,dep_o,idlon,idlat,hdat_a)

  use mod_rmiss
  !$use omp_lib  
  implicit none

  !---Common

  integer i_o,j_o,k_o
  integer m_o
  
  integer i_a0,i_a1
  integer j_a0,j_a1
  integer k_a
  integer m_a0,m_a1

  real(kind = 8) hdat_a00,hdat_a10,hdat_a01,hdat_a11
  
  !---IN
  !Analysis
  integer,intent(in) :: im_a,jm_a,km_a

  real(kind = 8),intent(in) :: lon_a(im_a),lat_a(jm_a),dep_a(im_a,jm_a,km_a)
  real(kind = 8),intent(in) :: dat_a(im_a,jm_a,km_a)

  !Observation
  integer,intent(in) :: im_o,jm_o,km_o

  real(kind = 8),intent(in) :: lon_o(im_o),lat_o(jm_o),dep_o(km_o)

  !ID
  integer,intent(in) :: idlon(im_o),idlat(jm_o)
  
  !---OUT
  real(kind = 8),intent(out) :: hdat_a(im_o,jm_o,km_o)
  
  !$omp parallel
  !$omp do private(i_o,j_o,k_o,i_a0,i_a1,j_a0,j_a1,k_a,m_a0,m_a1,m_o)
  do j_o=1,jm_o
     do i_o=1,im_o
        
        if(idlon(i_o) == 0 .or. idlat(j_o) == 0)then
           hdat_a(i_o,j_o,:)=rmiss
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
        
        j_a0=idlat(j_o)
        j_a1=idlat(j_o)+1
        
        k_a=1
        if(dat_a(i_a0,j_a0,k_a) == rmiss .or. dat_a(i_a1,j_a0,k_a) == rmiss &
             & .or. dat_a(i_a0,j_a1,k_a) == rmiss .or. dat_a(i_a1,j_a1,k_a) == rmiss)then
           hdat_a(i_o,j_o,:)=rmiss
           cycle
        end if

        do k_o=1,km_o     
        
           !Vertical linear interpolation
           call vertical_interpolation(km_a,dep_a(i_a0,j_a0,:),dat_a(i_a0,j_a0,:),dep_o(k_o),hdat_a00)
           call vertical_interpolation(km_a,dep_a(i_a1,j_a0,:),dat_a(i_a1,j_a0,:),dep_o(k_o),hdat_a10)
           call vertical_interpolation(km_a,dep_a(i_a0,j_a1,:),dat_a(i_a0,j_a1,:),dep_o(k_o),hdat_a01)
           call vertical_interpolation(km_a,dep_a(i_a1,j_a1,:),dat_a(i_a1,j_a1,:),dep_o(k_o),hdat_a11)
           
           !Bilinear interpolation
           if(hdat_a00 == rmiss .or. hdat_a10 == rmiss .or. hdat_a01 == rmiss .or. hdat_a11 == rmiss)then
              hdat_a(i_o,j_o,k_o)=rmiss
           else if(i_a0 == im_a)then
              call bilinear_interpolate( &
                   & lon_a(i_a0)-m_a0*360.d0,lon_a(i_a1)-(m_a1-1)*360.d0,lat_a(j_a0),lat_a(j_a1), &
                   & lon_o(i_o)-(m_o-1)*360.d0,lat_o(j_o),hdat_a00,hdat_a10,hdat_a01,hdat_a11,hdat_a(i_o,j_o,k_o))
           else
              call bilinear_interpolate( &
                   & lon_a(i_a0)-m_a0*360.d0,lon_a(i_a1)-m_a1*360.d0,lat_a(j_a0),lat_a(j_a1), &
                   & lon_o(i_o)-m_o*360.d0,lat_o(j_o),hdat_a00,hdat_a10,hdat_a01,hdat_a11,hdat_a(i_o,j_o,k_o))
           end if        
           
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine vertical_bilinear_interpolate
