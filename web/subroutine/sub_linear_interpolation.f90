!-----------------------------------------------------------
! Linear interpolation |
!-----------------------------------------------------------

subroutine linear_interpolation(x1,x2,y1,y2,x,y)

  use mod_rmiss
  implicit none

  !---IN
  real(kind = 8),intent(in) :: x1,x2,y1,y2
  real(kind = 8),intent(in) :: x

  !---OUT
  real(kind = 8),intent(out) :: y

  if(x1 == x2)then
     y=rmiss
  else
     y=(y2-y1)/(x2-x1)*(x-x1)+y1
  end if
     
end subroutine linear_interpolation

!-----------------------------------------------------------
! Vertical linear interpolation |
!-----------------------------------------------------------

subroutine vertical_linear_interpolation(im_in,jm_in,km_in,dep_in,mask_in,dat3d_in, &
                & idz,km_out,dep_out,dat3d_out)

  use mod_rmiss
  implicit none

  !---Common
  integer i,j,k_in,k_out

  real(kind = 8) z1,z2

  !---IN
  integer,intent(in) :: im_in,jm_in,km_in
  
  real(kind = 8),intent(in) :: dep_in(im_in,jm_in,km_in)
  real(kind = 8),intent(in) :: mask_in(im_in,jm_in)
  real(kind = 8),intent(in) :: dat3d_in(im_in,jm_in,km_in)

  integer,intent(in) :: idz(im_in,jm_in,km_out)
  integer,intent(in) :: km_out

  real(kind = 8),intent(in) :: dep_out(km_out)
  
  !---OUT
  real(kind = 8),intent(out) :: dat3d_out(im_in,jm_in,km_out)

  dat3d_out(:,:,:)=rmiss
  
  do k_out=1,km_out

     if(k_out == 1)then !surface

        do j=1,jm_in
           do i=1,im_in
              if(mask_in(i,j) == 0.d0)then
                 dat3d_out(i,j,k_out)=rmiss
              else
                 dat3d_out(i,j,k_out)=dat3d_in(i,j,1)
              end if
           end do
        end do

     else
     
        do j=1,jm_in
           do i=1,im_in

              if(idz(i,j,k_out) == 0 .or. mask_in(i,j) == 0.d0)then
                 dat3d_out(i,j,k_out)=rmiss
              else

                 k_in=idz(i,j,k_out)
                 z1=dep_in(i,j,k_in)
                 z2=dep_in(i,j,k_in+1)
                 
                 call linear_interpolation(z1,z2,dat3d_in(i,j,k_in),dat3d_in(i,j,k_in+1),dep_out(k_out),dat3d_out(i,j,k_out)) 
                 
              end if

           end do
        end do

     end if
     
  end do
  
end subroutine vertical_linear_interpolation
