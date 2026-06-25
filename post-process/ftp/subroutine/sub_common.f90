subroutine add(im,jm,km,mask,dat0,dat1,dat)

  !$ use omp_lib
  use mod_rmiss
  implicit none

  !---Common
  integer i,j,k

  !---IN
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: mask(im,jm)
  real(kind = 8),intent(in) :: dat0(im,jm,km),dat1(im,jm,km)

  !---OUT
  real(kind = 8),intent(out) :: dat(im,jm,km)

  !$omp parallel do private(i,j,k)
  do k=1,km
     do j=1,jm
        do i=1,im
           if(mask(i,j) == 0.d0)then
              dat(i,j,k)=rmiss
           else
              dat(i,j,k)=dat0(i,j,k)+dat1(i,j,k)
           end if
        end do
     end do
  end do
  !$omp end parallel do
                         
end subroutine add

!----------------------------------------------------------------------

subroutine month_ave(iday,nday,im,jm,km,mask,dat,ave,pass)

  !$ use omp_lib
  use mod_rmiss
  implicit none

  !---Common
  integer i,j,k
  
  !---IN
  integer,intent(in) :: iday,nday
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: mask(im,jm)
  real(kind = 8),intent(in) :: dat(im,jm,km)

  !---IN/OUT
  integer,intent(inout) :: pass(im,jm,km)
  
  real(kind = 8),intent(inout) :: ave(im,jm,km)

  if(iday == 1)then
     ave(:,:,:)=0.d0
     pass(:,:,:)=0
  end if

  !$omp parallel do private(i,j,k) collapse(3)
  do k=1,km
     do j=1,jm
        do i=1,im
           
           if(mask(i,j) == 0.d0 .or. dat(i,j,k) == rmiss) cycle            

           ave(i,j,k)=ave(i,j,k)+dat(i,j,k)
           pass(i,j,k)=pass(i,j,k)+1

        end do
     end do
  end do
  !$omp end parallel do

  if(iday == nday)then
     !$omp parallel do private(i,j,k) collapse(3)
     do k=1,km
        do j=1,jm
           do i=1,im

              if(pass(i,j,k) == 0)then
                 ave(i,j,k)=rmiss
              else
                 ave(i,j,k)=ave(i,j,k)/dble(pass(i,j,k))
              end if

           end do
        end do
     end do
     !$omp end parallel do
  end if
  
end subroutine month_ave
