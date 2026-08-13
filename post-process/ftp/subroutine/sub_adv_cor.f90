!-------------------------------------------------------------------------------
! Advection correction term (Sigma to Z)
!-------------------------------------------------------------------------------
! D = SSH + H
! xadv_cor = -u * (dSSH/dx + sigma * dD/dx)/D * dT/dsigma
! yadv_cor = -v * (dSSH/dy + sigma * dD/dy)/D * dT/dsigma
! zadv_cor = -xadv_cor - yadv_cor
!-------------------------------------------------------------------------------

subroutine adv_cor(im,jm,km,maskt,masku,maskv,dx,dy,sigma,depth,ssh,u,v,dat,xadv_cor,yadv_cor,zadv_cor)

  use setting, only: lqglobal
  implicit none

  !---Common
  integer i,j,k
  integer i1,i2,j1,j2,k1,k2

  real(kind = 8) D(im,jm)      !D = SSH + H at T point  
  real(kind = 8) Dx_u(im,jm)   !dD/dx at U point
  real(kind = 8) Dx_t(im,jm)   !dD/dx at T point
  real(kind = 8) Dy_v(im,jm)   !dD/dy at V point
  real(kind = 8) Dy_t(im,jm)   !dD/dy at T point
  
  real(kind = 8) sshx_u(im,jm) !dSSH/dx at U point
  real(kind = 8) sshx_t(im,jm) !dSSH/dx at T point
  real(kind = 8) sshy_v(im,jm) !dSSH/dy at V point
  real(kind = 8) sshy_t(im,jm) !dSSH/dy at T point

  real(kind = 8) u_t(im,jm,km) !U at T point  
  real(kind = 8) v_t(im,jm,km) !V at T point
  
  real(kind = 8) dats_w(im,jm,km) !dT/dsigma at W point
  real(kind = 8) dats_t(im,jm,km) !dT/dsigma at T point
  
  !---IN
  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: maskt(im,jm)        !Land-Sea mask at T point
  real(kind = 8),intent(in) :: masku(im,jm)        !Land-Sea mask at U point
  real(kind = 8),intent(in) :: maskv(im,jm)        !Land-Sea mask at V point
  real(kind = 8),intent(in) :: dx(im,jm),dy(im,jm) !Longitude and Latitude grid interval [m]
  real(kind = 8),intent(in) :: sigma(im,jm,km)     !Sigma at T point
  real(kind = 8),intent(in) :: depth(im,jm,km)     !Depth at W point    
  real(kind = 8),intent(in) :: ssh(im,jm)          !Sea surface height
  real(kind = 8),intent(in) :: u(im,jm,km),v(im,jm,km) !U/V
  real(kind = 8),intent(in) :: dat(im,jm,km)       !Temperature/Salinity

  !---OUT
  real(kind = 8),intent(out) :: xadv_cor(im,jm,km)
  real(kind = 8),intent(out) :: yadv_cor(im,jm,km)
  real(kind = 8),intent(out) :: zadv_cor(im,jm,km)

  !--- D = eta + H at T point
  do j=1,jm
     do i=1,im
        D(i,j)=maskt(i,j)*(ssh(i,j)+abs(depth(i,j,km)))
     end do
  end do

  !---U at T point
  u_t(:,:,:)=0.d0
  do k=1,km-1
     do j=1,jm        
        do i=1,im

           if(i == im .and. lqglobal)then
              i1=i
              i2=1
           else if(i == im)then
              i1=i
              i2=i
           else
              i1=i
              i2=i+1
           end if

           if(masku(i1,j)+masku(i2,j) == 0.d0)then
              u_t(i,j,k)=0.d0
           else              
              u_t(i,j,k)=maskt(i,j)*( masku(i1,j)*u(i1,j,k)+masku(i2,j)*u(i2,j,k) ) / ( masku(i1,j)+masku(i2,j) )          
           end if
              
        end do
     end do     
  end do
     
  !---V at T point
  v_t(:,:,:)=0.d0
  do k=1,km-1
     do j=1,jm

        if(j == jm)then
           j1=j
           j2=j
        else
           j1=j
           j2=j+1
        end if
        
        do i=1,im

           if(maskv(i,j1)+maskv(i,j2) == 0.d0)then
              v_t(i,j,k)=0.d0
           else
              v_t(i,j,k)=maskt(i,j)*( maskv(i,j1)*v(i,j1,k)+maskv(i,j2)*v(i,j2,k) ) / ( maskv(i,j1)+maskv(i,j2) )
           end if
              
        end do
     end do
  end do

  !---dSSH/dx and dD/dx at U point
  !---dSSH/dy and dD/dy at V point
  sshx_u(:,:)=0.d0
  sshy_v(:,:)=0.d0
  Dx_u(:,:)=0.d0
  Dy_v(:,:)=0.d0
  do j=1,jm     

     if(j == 1)then
        j1=j
        j2=j+1
     else
        j1=j-1
        j2=j
     end if
        
     do i=1,im

        if(i == 1 .and. lqglobal)then
           i1=im
           i2=i
        else if(i == 1)then
           i1=i
           i2=i+1
        else
           i1=i-1
           i2=i
        end if           

        !X differentiation
        sshx_u(i,j)=masku(i,j)*( ssh(i2,j)-ssh(i1,j) )/dx(i,j)
        Dx_u(i,j)=masku(i,j)*( D(i2,j)-D(i1,j) )/dx(i,j)

        !Y differentiation
        sshy_v(i,j)=maskv(i,j)*( ssh(i,j2)-ssh(i,j1) )/dy(i,j)
        Dy_v(i,j)=maskv(i,j)*( D(i,j2)-D(i,j1) )/dy(i,j)        
        
     end do
  end do

  !---dSSH/dx, dSSH/dy, dD/dx, dD/dy at T point
  sshx_t(:,:)=0.d0
  sshy_t(:,:)=0.d0
  Dx_t(:,:)=0.d0
  Dy_t(:,:)=0.d0
  
  do j=1,jm

     if(j == jm)then
        j1=j
        j2=j
     else
        j1=j
        j2=j+1
     end if
     
     do i=1,im
        
        if(i == im .and. lqglobal)then
           i1=i
           i2=1
        else if(i == im)then
           i1=i
           i2=i
        else
           i1=i
           i2=i+1
        end if
  
        if(masku(i2,j)+masku(i1,j) == 0.d0)then
           sshx_t(i,j)=0.d0
           Dx_t(i,j)=0.d0
        else
           sshx_t(i,j)= &
                & maskt(i,j)*( masku(i2,j)*sshx_u(i2,j)+masku(i1,j)*sshx_u(i1,j) ) / ( masku(i2,j)+masku(i1,j) )
           Dx_t(i,j)= &
                & maskt(i,j)*( masku(i2,j)*Dx_u(i2,j)+masku(i1,j)*Dx_u(i1,j) ) / ( masku(i2,j)+masku(i1,j) )
        end if


        if(maskv(i,j2)+maskv(i,j1) == 0.d0)then
           sshy_t(i,j)=0.d0
           Dy_t(i,j)=0.d0
        else
           sshy_t(i,j)=maskt(i,j)*( maskv(i,j1)*sshy_v(i,j1)+maskv(i,j2)*sshy_v(i,j2) ) / ( maskv(i,j1)+maskv(i,j2) )
           Dy_t(i,j)=maskt(i,j)*( maskv(i,j1)*Dy_v(i,j1)+maskv(i,j2)*Dy_v(i,j2) ) / ( maskv(i,j1)+maskv(i,j2) )
        end if

     end do
  end do
        
  !---dT/dsigma at W point
  dats_w(:,:,:)=0.d0
  do k=1,km-1

     if(k == 1)then
        k1=k+1
        k2=k
     else
        k1=k
        k2=k-1
     end if
     
     do j=1,jm
        do i=1,im
           dats_w(i,j,k)=maskt(i,j)*( dat(i,j,k2)-dat(i,j,k1) ) / ( sigma(i,j,k2)-sigma(i,j,k1) )
        end do
     end do
  end do    

  !---dT/dsigma at T point
  dats_t(:,:,:)=0.d0  
  do k=1,km-1

     if(k == km-1)then
        k1=k
        k2=k
     else
        k1=k
        k2=k+1
     end if
     
     do j=1,jm
        do i=1,im
           dats_t(i,j,k)=maskt(i,j)*0.5d0*( dats_w(i,j,k2)+dats_w(i,j,k1) )
        end do
     end do
  end do    

  !---Advection correction term
  xadv_cor(:,:,:)=0.d0
  yadv_cor(:,:,:)=0.d0
  zadv_cor(:,:,:)=0.d0
  do k=1,km-1
     do j=1,jm
        do i=1,im

           if(maskt(i,j) == 0.d0)then
              xadv_cor(i,j,k)=0.d0
              yadv_cor(i,j,k)=0.d0
              zadv_cor(i,j,k)=0.d0
           else
              
              xadv_cor(i,j,k)= -1.d0 * 86400.d0 * u_t(i,j,k) * ( sshx_t(i,j)+sigma(i,j,k)*Dx_t(i,j) ) * dats_t(i,j,k) / D(i,j)
              yadv_cor(i,j,k)= -1.d0 * 86400.d0 * v_t(i,j,k) * ( sshy_t(i,j)+sigma(i,j,k)*Dy_t(i,j) ) * dats_t(i,j,k) / D(i,j)
              zadv_cor(i,j,k)= -xadv_cor(i,j,k)-yadv_cor(i,j,k)
           end if           
           
        end do
     end do
  end do
  
end subroutine adv_cor
