!--------------------------------------------------------------
! Potential Density |
!--------------------------------------------------------------

subroutine potential_density_3d(im,jm,km,mask,t,s,rho)

  !$ use omp_lib
  use mod_rmiss
  implicit none

  !---Common
  integer i,j,k

  real(kind = 8) depth

  !----IN
  integer,intent(in) :: im,jm,km
  
  real(kind = 8),intent(in) :: mask(im,jm)
  real(kind = 8),intent(in) :: t(im,jm,km),s(im,jm,km)
  
  !---OUT
  real(kind = 8),intent(out) :: rho(im,jm,km)

  !$omp parallel do private(i,j,k,depth)
  do k=1,km
     do j=1,jm
        do i=1,im
           if(mask(i,j) == 0.d0)then
              rho(i,j,k)=rmiss
           else
              depth=0.d0
              call estimate_density(t(i,j,k),s(i,j,k),depth,rho(i,j,k))
           end if
        end do
     end do
  end do
  !$omp end parallel do

end subroutine potential_density_3d

!--------------------------------------------------------------
! MLD average |
!--------------------------------------------------------------

subroutine estimate_mld(im,jm,km,mask,depth,rho,id_mld,mld)

  !$ use omp_lib
  use mod_rmiss
  use setting,only: drho
  implicit none

  !---Common
  integer i,j,k
  
  !---IN
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: mask(im,jm)
  real(kind = 8),intent(in) :: depth(im,jm,km),rho(im,jm,km) 

  !---OUT
  integer,intent(out) :: id_mld(im,jm)

  real(kind = 8),intent(out) :: mld(im,jm)

  !$omp parallel do private(i,j) collapse(2)  
  do j=1,jm
     do i=1,im

        !Land
        if(mask(i,j) == 0.d0)then
           id_mld(i,j)=0
           mld(i,j)=rmiss
           cycle
        end if

        id_mld(i,j)=km-1
        mld(i,j)=depth(i,j,km-1)
        
        do k=1,km-2
           if(drho <= rho(i,j,k+1)-rho(i,j,1))then
              id_mld(i,j)=k
              mld(i,j)=depth(i,j,k)
              exit
           end if
        end do !k

     end do !i
  end do !j
  !$omp end parallel do
  
end subroutine estimate_mld

!--------------------------------------------------------------
! MLD average |
!--------------------------------------------------------------

subroutine mld_ave(im,jm,km,id_mld,mask,depth,dat,ave)

  !$ use omp_lib
  use mod_rmiss
  implicit none

  !---Common
  integer i,j,k

  real(kind = 8) dz,dz_sum
  
  !---IN
  integer,intent(in) :: im,jm,km
  integer,intent(in) :: id_mld(im,jm)

  real(kind = 8),intent(in) :: mask(im,jm)
  real(kind = 8),intent(in) :: depth(im,jm,km),dat(im,jm,km)

  !---OUT
  real(kind = 8) ave(im,jm)

  !$omp parallel do private(i,j,dz,dz_sum) collapse(2)
  do j=1,jm
     do i=1,im

        ave(i,j)=0.d0
        dz_sum=0.d0
        
        if(mask(i,j) == 0.d0)then
           ave(i,j)=rmiss
           cycle
        end if

        !MLD = 1st layer
        if(id_mld(i,j) == 1)then
           ave(i,j)=dat(i,j,1)
           cycle
        end if

        !Vertical integration
        do k=1,id_mld(i,j)-1

           dz=depth(i,j,k+1)-depth(i,j,k)
           ave(i,j)=ave(i,j)+0.5d0*(dat(i,j,k+1)+dat(i,j,k))*dz
           dz_sum=dz_sum+dz
           
        end do !k

        if(dz_sum == 0.d0)then
           ave(i,j)=rmiss
        else
           ave(i,j)=ave(i,j)/dz_sum
        end if
        
     end do !i
  end do !j
  !$omp end parallel do
  
end subroutine mld_ave

!--------------------------------------------------------------
! Detrained/Entrained water |
!--------------------------------------------------------------

subroutine det_ent_water(im,jm,km,id_mld_n1,id_mld,mask,depth,dat,ave,mld_ent,dhdt)

  !$ use omp_lib
  use mod_rmiss
  implicit none

  !---Common
  integer i,j,k

  real(kind = 8) dz,dz_sum
  
  !---IN
  integer,intent(in) :: im,jm,km
  integer,intent(in) :: id_mld_n1(im,jm),id_mld(im,jm)

  real(kind = 8),intent(in) :: mask(im,jm)
  real(kind = 8),intent(in) :: depth(im,jm,km),dat(im,jm,km)

  !---OUT
  real(kind = 8) ave(im,jm)

  !$omp parallel do private(i,j,dz,dz_sum) collapse(2)
  do j=1,jm
     do i=1,im

        ave(i,j)=0.d0
        dz_sum=0.d0
        
        if(mask(i,j) == 0.d0)then
           ave(i,j)=rmiss
           cycle
        end if

        
        
        !MLD = 1st layer
        if(id_mld(i,j) == 1)then
           ave(i,j)=dat(i,j,1)
           cycle
        end if

        !Vertical integration
        do k=1,id_mld(i,j)-1

           dz=depth(i,j,k+1)-depth(i,j,k)
           ave(i,j)=ave(i,j)+0.5d0*(dat(i,j,k+1)+dat(i,j,k))*dz
           dz_sum=dz_sum+dz
           
        end do !k

        if(dz_sum == 0.d0)then
           ave(i,j)=rmiss
        else
           ave(i,j)=ave(i,j)/dz_sum
        end if
        
     end do !i
  end do !j
  !$omp end parallel do
  
end subroutine mld_ave
