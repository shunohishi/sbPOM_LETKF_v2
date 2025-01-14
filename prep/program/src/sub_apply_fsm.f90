!----------------------------------------------------------------------------
! Apply fsm |
!----------------------------------------------------------------------------

subroutine apply_fsm(im,jm,km,dat,fsm)

  !$use omp_lib
  implicit none
  
  integer i,j,k

  integer,intent(in) :: im,jm,km
  real(kind = 8),intent(in) :: fsm(im,jm)

  real(kind = 8),intent(inout) :: dat(im,jm,km)

  !$omp parallel
  !$omp do private(i,j,k)  
  do k=1,km
     do j=1,jm
        do i=1,im
           dat(i,j,k)=dat(i,j,k)*fsm(i,j)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
           
end subroutine apply_fsm

!------------------------------------------------------------------------
! Apply JRA55 Land |
!------------------------------------------------------------------------

subroutine apply_jra55_land(im,jm,dat,land,rmiss)

  !$use omp_lib  
  implicit none

  integer i,j
  
  integer,intent(in) :: im,jm
  real(kind = 8),intent(in) :: land(im,jm) !1:Land, 0:Sea
  real(kind = 8),intent(in) :: rmiss

  real(kind = 8),intent(inout) :: dat(im,jm)

  !$omp parallel
  !$omp do private(i,j)
  do j=1,jm
     do i=1,im
        if(land(i,j) == 1.d0)then
           dat(i,j)=rmiss
        end if
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine apply_jra55_land
