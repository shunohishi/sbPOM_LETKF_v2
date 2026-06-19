subroutine make_dzt(im,jm,km,maskt,depw,dzt)

  !$ use omp_lib
  implicit none

  !---Common
  integer i,j,k
  
  !---IN
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: maskt(im,jm)
  real(kind = 8),intent(in) :: depw(im,jm,km)

  !---OUT
  real(kind = 8),intent(out) :: dzt(im,jm,km)

  !---Initialize
  dzt(:,:,:)=0.d0
  
  !$omp parallel do private(i,j,k) collpase(2)
  do k=1,km-1
     do j=1,jm
        do i=1,im
           dzt(i,j,k)=maskt(i,j)*(depw(i,j,k+1)-depw(i,j,k))           
        end do
     end do
  end do
  !$omp end parallel do

end subroutine make_dzt
