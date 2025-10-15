!---------------------------------------------------------------
! Calculate ID |
!---------------------------------------------------------------
!
! 1  ----i2----2
! i1 ----i2----i1+1
! im1----i2----1
! 
! * id denotes i1, when dat2 is at i2
!
!---------------------------------------------------------------

subroutine get_id(im1,dat1,im2,dat2,id)

  !$use omp_lib  
  implicit none

  !---Common
  integer i1,i2
  
  !---IN
  integer,intent(in) :: im1,im2
  real(kind = 8),intent(in) :: dat1(im1),dat2(im2)

  !---OUT
  integer,intent(out) :: id(im2)

  id(:)=0

  !$omp parallel
  !$omp do private(i1,i2)
  do i2=1,im2
     do i1=1,im1-1

        if(dat1(i1) <= dat2(i2) .and. dat2(i2) <= dat1(i1+1))then
           id(i2)=i1
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel
    
end subroutine get_id
