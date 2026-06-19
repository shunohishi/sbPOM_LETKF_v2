!--------------------------------------------------------------
! Potential Density |
!--------------------------------------------------------------

subroutine potential_density_3d(im,jm,km,mask,t,s,rho)

  !$ use omp_lib
  use mod_rmiss
  use mod_density
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
  real(kind = 8),intent(in) :: depth(im,jm,km) !*W point
  real(kind = 8),intent(in) :: rho(im,jm,km) 

  !---OUT
  integer,intent(out) :: id_mld(im,jm)

  real(kind = 8),intent(out) :: mld(im,jm)

  !---Initialize
  id_mld(:,:)=0

  !---ID and MLD
  !$omp parallel do private(i,j,k) collapse(2)  
  do j=1,jm
     do i=1,im

        !Land
        if(mask(i,j) == 0.d0)then
           id_mld(i,j)=0
           mld(i,j)=rmiss
           cycle
        end if

        !Initialize (Mixed over entire column)
        id_mld(i,j)=km-1
        mld(i,j)=depth(i,j,km)

        !MLD
        do k=1,km-2
           if(drho <= rho(i,j,k+1)-rho(i,j,1))then
              id_mld(i,j)=k
              mld(i,j)=depth(i,j,k+1)
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

subroutine mld_ave(im,jm,km,id_mld,mask,dz,dat,ave)

  !$ use omp_lib
  use mod_rmiss
  implicit none

  !---Common
  integer i,j,k

  real(kind = 8) dz_sum
  
  !---IN
  integer,intent(in) :: im,jm,km
  integer,intent(in) :: id_mld(im,jm)

  real(kind = 8),intent(in) :: mask(im,jm)
  real(kind = 8),intent(in) :: dz(im,jm,km),dat(im,jm,km)

  !---OUT
  real(kind = 8),intent(out) :: ave(im,jm)

  !$omp parallel do private(i,j,k,dz_sum) collapse(2)
  do j=1,jm
     do i=1,im
        
        if(mask(i,j) == 0.d0)then
           ave(i,j)=rmiss
           cycle
        end if

        !Vertical integration
        ave(i,j)=0.d0
        dz_sum=0.d0
        do k=1,id_mld(i,j)

           dz_sum=dz_sum+dz(i,j,k)
           ave(i,j)=ave(i,j)+dat(i,j,k)*dz(i,j,k)
           
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
! Detrainment and Entrainment layer average |
!--------------------------------------------------------------

subroutine det_ent_ave(im,jm,km,id_h1,id_mld,mask,dz,dat,ave1,ave2)

  !$ use omp_lib
  use mod_rmiss
  implicit none

  !---Common
  integer i,j,k

  real(kind = 8) dz_sum
  
  !---IN
  integer,intent(in) :: im,jm,km
  integer,intent(in) :: id_h1(im,jm),id_mld(im,jm)

  real(kind = 8),intent(in) :: mask(im,jm)
  real(kind = 8),intent(in) :: dz(im,jm,km),dat(im,jm,km),ave1(im,jm)

  !---OUT
  real(kind = 8),intent(out) :: ave2(im,jm)

  !$omp parallel do private(i,j,k,dz_sum) collapse(2)
  do j=1,jm
     do i=1,im
        
        if(mask(i,j) == 0.d0)then
           ave2(i,j)=rmiss
           cycle
        end if

        if(id_h1(i,j) == id_mld(i,j))then
           ave2(i,j)=ave1(i,j)
           cycle
        end if
        
        !Vertical integration
        ave2(i,j)=0.d0
        dz_sum=0.d0
        do k=id_h1(i,j)+1,id_mld(i,j)

           dz_sum=dz_sum+dz(i,j,k)
           ave2(i,j)=ave2(i,j)+dat(i,j,k)*dz(i,j,k)
           
        end do !k

        if(dz_sum == 0.d0)then
           ave2(i,j)=rmiss
        else
           ave2(i,j)=ave2(i,j)/dz_sum
        end if
        
     end do !i
  end do !j
  !$omp end parallel do
  
end subroutine det_ent_ave

!--------------------------------------------------------------
! Detrainment and Entrainment (Kim et al. 2006) |           
!--------------------------------------------------------------

subroutine detrainment_entrainment(im,jm,km,mask,dz, &
     & id_mld_t0,mld_t0,t_t0,s_t0, &
     & id_mld_t1,mld_t1, &
     & mld,dhdt,delta_t,delta_s,tent,sent)

  !$ use omp_lib
  use mod_rmiss
  implicit none

  !---Common
  integer i,j,k
  integer id_h1(im,jm),id_mld(im,jm)

  real(kind = 8) t1_t0(im,jm),s1_t0(im,jm) !T and S averaged over h1 at t0 
  real(kind = 8) t2_t0(im,jm),s2_t0(im,jm) !T and S averaged over h2 at t0

  !---IN
  integer,intent(in) :: im,jm,km
  integer,intent(in) :: id_mld_t0(im,jm),id_mld_t1(im,jm) !Depth ID at t0 and t1

  real(kind = 8),intent(in) :: mask(im,jm)  !Land-Sea Mask
  real(kind = 8),intent(in) :: dz(im,jm,km) !Layer thickness

  real(kind = 8),intent(in) :: mld_t0(im,jm),mld_t1(im,jm)   !MLD at t0 and t1

  real(kind = 8),intent(in) :: t_t0(im,jm,km) !Temperature at t0 and t1
  real(kind = 8),intent(in) :: s_t0(im,jm,km) !Salinity at t0 and t1

  !---OUT
  real(kind = 8),intent(out) :: mld(im,jm)     !MLD=h1+h2
  real(kind = 8),intent(out) :: dhdt(im,jm)    !MLD difference *Unit: [m]
  real(kind = 8),intent(out) :: delta_t(im,jm) !Detrained/Entrained temperature
  real(kind = 8),intent(out) :: delta_s(im,jm) !Detrained/Entrained salnity
  real(kind = 8),intent(out) :: tent(im,jm)    !Temperature detrainment/entrainement
  real(kind = 8),intent(out) :: sent(im,jm)    !Salinity detrainment/entrainement

  !---ID, MLD, dhdt
  id_mld(:,:)=0
  id_h1(:,:)=0
  !$omp parallel do private(i,j) collapse(2)
  do j=1,jm
     do i=1,im

        !---Land
        if(mask(i,j) == 0.d0)then
           mld(i,j)=rmiss
           dhdt(i,j)=rmiss
           cycle
        end if

        !---MLD=h1+h2
        id_mld(i,j)=max(id_mld_t0(i,j),id_mld_t1(i,j))
        mld(i,j)=max(mld_t0(i,j),mld_t1(i,j))

        !---h1
        id_h1(i,j)=min(id_mld_t0(i,j),id_mld_t1(i,j))

        !---dhdt
        dhdt(i,j)=mld_t1(i,j)-mld_t0(i,j)        

     end do
  end do
  !$omp end parallel do

  call mld_ave(im,jm,km,id_h1,mask,dz,t_t0,t1_t0)
  call mld_ave(im,jm,km,id_h1,mask,dz,s_t0,s1_t0)
  call det_ent_ave(im,jm,km,id_h1,id_mld,mask,dz,t_t0,t1_t0,t2_t0)
  call det_ent_ave(im,jm,km,id_h1,id_mld,mask,dz,s_t0,s1_t0,s2_t0)
  
  !$omp parallel do private(i,j) collapse(2)
  do j=1,jm
     do i=1,im

        if(mask(i,j) == 0.d0)then
           delta_t(i,j)=rmiss
           delta_s(i,j)=rmiss
           tent(i,j)=rmiss
           sent(i,j)=rmiss
           cycle
        end if

        !---Delta T and S        
        if(t1_t0(i,j) == rmiss .or. t2_t0(i,j) == rmiss)then
           delta_t(i,j)=rmiss
        else
           delta_t(i,j)=t1_t0(i,j)-t2_t0(i,j)
        end if

        if(s1_t0(i,j) == rmiss .or. s2_t0(i,j) == rmiss)then
           delta_s(i,j)=rmiss
        else
           delta_s(i,j)=s1_t0(i,j)-s2_t0(i,j)
        end if

        !---Detrainment and Entrainment
        if(mld(i,j) == rmiss .or. delta_t(i,j) == rmiss)then
           tent(i,j)=rmiss
        else
           tent(i,j)=-1.d0*dhdt(i,j)*delta_t(i,j)/mld(i,j)
        end if

        if(mld(i,j) == rmiss .or. delta_s(i,j) == rmiss)then
           sent(i,j)=rmiss
        else
           sent(i,j)=-1.d0*dhdt(i,j)*delta_s(i,j)/mld(i,j)
        end if
        
     end do !i
  end do !j
  !$omp end parallel do

end subroutine detrainment_entrainment

!--------------------------------------------------------------------
! d MLT/dt |
!--------------------------------------------------------------------

subroutine dmltdt(im,jm,mask,dat0,dat1,dat)

  !$ use omp_lib
  use mod_rmiss
  implicit none

  !---Common
  integer i,j
  
  !---IN
  integer,intent(in) :: im,jm
  
  real(kind = 8),intent(in) :: mask(im,jm)
  real(kind = 8),intent(in) :: dat0(im,jm),dat1(im,jm)

  !---OUT
  real(kind = 8),intent(out) :: dat(im,jm)

  !$omp parallel do private(i,j) collapse(2)  
  do j=1,jm
     do i=1,im
        if(mask(i,j) == 0.d0)then
           dat(i,j)=rmiss
        else
           dat(i,j)=dat1(i,j)-dat0(i,j)
        end if
     end do
  end do
  !$omp end parallel do
    
end subroutine dmltdt
