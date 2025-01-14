subroutine whole_average

  !use omp_lib
  use common_pom_var
  implicit none

  !2D
  call calculate_whole_average(iint,iend,im_local,jm_local,1, &
       & el,el_ave,fsm)

  call calculate_whole_average(iint,iend,im_local,jm_local,1, &
       & ua,ua_ave,dum)

  call calculate_whole_average(iint,iend,im_local,jm_local,1, &
       & va,va_ave,dvm)

  !3D
  call calculate_whole_average(iint,iend,im_local,jm_local,kb, &
       & t,t_ave,fsm)

  call calculate_whole_average(iint,iend,im_local,jm_local,kb, &
       & s,s_ave,fsm)

  call calculate_whole_average(iint,iend,im_local,jm_local,kb, &
       & u,u_ave,dum)

  call calculate_whole_average(iint,iend,im_local,jm_local,kb, &
       & v,v_ave,dvm)

end subroutine whole_average

!-------------------------------------------------------------------------------

subroutine calculate_whole_average(iint,nint,im,jm,km,dat,ave,fsm)

  !$use omp_lib
  use common_pom_var, only: r_size
  implicit none

  integer i,j,k

  integer,intent(in) :: iint,nint !timestep
  integer,intent(in) :: im,jm,km  !grid size

  real(kind = r_size),intent(in) :: dat(im,jm,km),fsm(im,jm)
  real(kind = r_size),intent(inout) :: ave(im,jm,km)

  if(iint == 1)then
     ave(:,:,:)=0.d0
  end if

  !$omp parallel
  !$omp do private(i,j,k)
  do k=1,km
     do j=1,jm
        do i=1,im
           ave(i,j,k)=ave(i,j,k)+fsm(i,j)*dat(i,j,k)
        end do
     end do
  end do
  !$omp end do

  if(iint == nint)then
     !$omp do private(i,j,k)
     do k=1,km
        do j=1,jm
           do i=1,im
              ave(i,j,k)=fsm(i,j)*ave(i,j,k)/dble(nint)
           end do
        end do
     end do
     !$omp end do
  end if
  !$omp end parallel

end subroutine calculate_whole_average
