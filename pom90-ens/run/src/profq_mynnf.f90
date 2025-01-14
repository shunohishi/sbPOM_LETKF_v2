!_______________________________________________________________________
subroutine profq_mynnf
  !     solve for q2 (twice the turbulent kinetic energy), q2l (q2 x turbulent
  !     length scale), km (vertical kinematic viscosity) and kh (vertical
  !     kinematic diffusivity), using a simplified version of the level 2 1/2
  !     model of Mellor and Yamada (1982)
  !     in this version, the Craig-Banner sub-model whereby breaking wave tke
  !     is injected into the surface is included. However, we use an
  !     analytical solution to the near surface tke equation to solve for q2
  !     at the surface giving the same result as C-B diffusion. The new scheme
  !     is simpler and more robust than the latter scheme
  !     
  !     including Mellor-Yamada-Nakanishi-Niino-Furuichi version
  !     see lmynnf
  !

  !$use omp_lib  
  use MYNNF_lev25_2012
  use common_pom_var
  implicit none

  integer i,j,k

  real(kind = r_size),parameter :: a1=0.92d0,a2=16.6d0,b1=0.74d0,b2=10.1d0,c1=0.08d0
  real(kind = r_size),parameter :: e1=1.8d0,e2=1.33d0
  real(kind = r_size),parameter :: sef=1.d0
  real(kind = r_size),parameter :: cbcnst=50.d0,surfl=1.d5 !surface wave breaking

  real(kind = r_size),parameter :: coef4=18.d0*a1*a1+9.d0*a1*a2
  real(kind = r_size),parameter :: coef5=9.d0*a1*a2
  
  real(kind = r_size) a(im,jm,kb),c(im,jm,kb)
  real(kind = r_size) ee(im,jm,kb),gg(im,jm,kb)
  real(kind = r_size) sm(im,jm,kb),sh(im,jm,kb)
  real(kind = r_size) cc(im,jm,kb)
  real(kind = r_size) gh(im,jm,kb),boygr(im,jm,kb),dh(im,jm),stf(im,jm,kb)
  real(kind = r_size) prod(im,jm,kb)
  real(kind = r_size) coef1,coef2,coef3
  real(kind = r_size) const1
  real(kind = r_size) p,sp,tp
  real(kind = r_size) l0(im,jm)
  real(kind = r_size) shiw
  real(kind = r_size) utau2(im,jm)
  
  if(liwbrk)then
     shiw=0.7d0
  else
     shiw=0.0d0
  end if

  l0(1:im,1:jm)=0.d0
  a(1:im,1:jm,1:kb)=0.d0
  c(1:im,1:jm,1:kb)=0.d0
  ee(1:im,1:jm,1:kb)=0.d0
  gg(1:im,1:jm,1:kb)=0.d0
  boygr(1:im,1:jm,1:kb)=0.d0
  prod(1:im,1:jm,1:kb)=0.d0
  cc(1:im,1:jm,1:kb)=0.d0

  !$omp parallel
  !$omp do private(i,j)
  do j=1,jm
     do i=1,im
        dh(i,j)=h(i,j)+etf(i,j)
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)  
  do k=2,kbm1
     do j=2,jmm1
        do i=2,imm1
           a(i,j,k)=-dti2*(kq(i,j,k+1)+kq(i,j,k)+2.d0*umol)*0.5d0 &
                & /(dzz(i,j,k-1)*dz(i,j,k)*dh(i,j)*dh(i,j))
           c(i,j,k)=-dti2*(kq(i,j,k-1)+kq(i,j,k)+2.d0*umol)*0.5d0 &
                & /(dzz(i,j,k-1)*dz(i,j,k-1)*dh(i,j)*dh(i,j))
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  !     the following section solves the equation:
  !     dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b

  !     surface and bottom boundary conditions
  if(lmynnf)then
     const1=(24.0d0**(2.d0/3.d0))*sef
  else
     const1=(16.6d0**(2.d0/3.d0))*sef
  endif
     
  !     initialize fields that are not calculated on all boundaries
  !     but are later used there

  !$omp parallel
  !$omp do private(i,j)
  do j=2,jmm1
     do i=2,imm1
        utau2(i,j)=sqrt((0.5d0*(wusurf(i,j)+wusurf(i+1,j)))**2 &
             & +(0.5d0*(wvsurf(i,j)+wvsurf(i,j+1)))**2)
        uf(i,j,kb)=sqrt((0.5d0*(wubot(i,j)+wubot(i+1,j)))**2 &
             & +(0.5d0*(wvbot(i,j)+wvbot(i,j+1)))**2)*const1
     end do
  end do
  !$omp end do
  
  if(lwbrk) then

     !$omp do private(i,j)
     do j=2,jmm1
        do i=2,imm1
           !     wave breaking energy- a variant of Craig & Banner (1994)
           !     see Mellor and Blumberg, 2004.
           ee(i,j,1)=0.d0
           gg(i,j,1)=(15.8d0*cbcnst)**(2.d0/3.d0)*utau2(i,j)
           !     surface length scale following Stacey (1999).
           l0(i,j)=surfl*utau2(i,j)/grav
        end do
     end do
     !$omp end do
     
  end if
  
  !     calculate speed of sound squared
  !$omp do private(i,j,k,tp,sp,p)
  do k=1,kbm1
     do j=2,jmm1
        do i=2,imm1
           tp=t(i,j,k)+tbias
           sp=s(i,j,k)+sbias
           !     calculate pressure in units of decibars
           p=grav*rhoref*(-zz(i,j,k)*h(i,j))*1.d-4
           cc(i,j,k)=1449.1d0+.00821d0*p+4.55d0*tp-0.045d0*tp**2+1.34d0*(sp-35.0d0)
           cc(i,j,k)=cc(i,j,k)/sqrt((1.d0-0.01642d0*p/cc(i,j,k)) &
                & *(1.d0-0.4d0*p/cc(i,j,k)**2))
        end do
     end do
  end do
  !$omp end do
  
  !     calculate buoyancy gradient
  !$omp do private(i,j,k)
  do k=2,kbm1
     do j=2,jmm1
        do i=2,imm1
           q2b(i,j,k)=abs(q2b(i,j,k))
           q2lb(i,j,k)=abs(q2lb(i,j,k))
           boygr(i,j,k)=grav*(rho(i,j,k-1)-rho(i,j,k))/(dzz(i,j,k-1)*h(i,j)) &
           !     *** note: comment out next line if dens does not include pressure
                & +(grav**2)*2.d0/(cc(i,j,k-1)**2+cc(i,j,k)**2)
        end do
     end do
  end do
  !$omp end do

  if(.not. lmynnf)then

     !$omp do private(i,j,k)
     do k=2,kbm1
        do j=2,jmm1
           do i=2,imm1
              l(i,j,k)=abs(q2lb(i,j,k)/q2b(i,j,k))
              if(z(i,j,k) > -0.5d0) l(i,j,k)=max(l(i,j,k),kappa*l0(i,j))
              gh(i,j,k)=(l(i,j,k)**2)*boygr(i,j,k)/q2b(i,j,k)
              gh(i,j,k)=min(gh(i,j,k),0.028d0)
           end do
        end do
     end do
     !$omp end do

     !$omp do private(i,j)     
     do j=2,jmm1
        do i=2,imm1
           l(i,j,1)=kappa*l0(i,j)
           l(i,j,kb)=0.d0
           gh(i,j,1)=0.d0
           gh(i,j,kb)=0.d0
        end do
     end do
     !$omp end do
     
  end if
  
  !     calculate production of turbulent kinetic energy:
  !$omp do private(i,j,k)
  do k=2,kbm1
     do j=2,jmm1
        do i=2,imm1
           prod(i,j,k)=km(i,j,k)*0.25d0*sef &
                & *((u(i,j,k)-u(i,j,k-1)+u(i+1,j,k)-u(i+1,j,k-1))**2 &
                & +(v(i,j,k)-v(i,j,k-1)+v(i,j+1,k)-v(i,j+1,k-1))**2) &
                & /(dzz(i,j,k-1)*dh(i,j))**2 &
           !     add shear due to internal wave field
                & -shiw*km(i,j,k)*boygr(i,j,k)
           prod(i,j,k)=prod(i,j,k)+kh(i,j,k)*boygr(i,j,k)
        end do
     end do
  end do
  !$omp end do

  if(lmynnf) then

     !$omp do private(i,j,k)
     do k=1,kb
        do j=2,jmm1
           do i=2,imm1
              stf(i,j,k)=1.d0
              dtef(i,j,k)=sqrt(abs(q2b(i,j,k)))*stf(i,j,k)/(24.0d0*l(i,j,k)+small)
           end do
        end do
     end do
     !$omp end do
     
  else

     !$omp do private(i,j,k)     
     do k=1,kb
        do j=2,jmm1
           do i=2,imm1
              stf(i,j,k)=1.d0
              dtef(i,j,k)=sqrt(abs(q2b(i,j,k)))*stf(i,j,k)/(b1*l(i,j,k)+small)
           end do
        end do
     end do
     !$omp end do
     
  end if
  !$omp end parallel

  do k=2,kbm1
     do j=2,jmm1
        do i=2,imm1
           gg(i,j,k)=1.d0/(a(i,j,k)+c(i,j,k)*(1.d0-ee(i,j,k-1))-(2.d0*dti2*dtef(i,j,k)+1.d0))
           ee(i,j,k)=a(i,j,k)*gg(i,j,k)
           gg(i,j,k)=(-2.d0*dti2*prod(i,j,k)+c(i,j,k)*gg(i,j,k-1)-uf(i,j,k))*gg(i,j,k)
        end do
     end do
  end do

  do k=kbm1,1,-1
     do j=2,jmm1
        do i=2,imm1
           uf(i,j,k)=ee(i,j,k)*uf(i,j,k+1)+gg(i,j,k)
        end do
     end do
  end do
  
  if(.not. lmynnf) then

     !     the following section solves the equation:
     !     dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb

     vf(2:imm1,2:jmm1,1)=0.d0
     vf(2:imm1,2:jmm1,kb)=0.d0
     ee(2:imm1,2:jmm1,2)=0.d0

     !$omp parallel
     !$omp do private(i,j)     
     do j=2,jmm1
        do i=2,imm1
           gg(i,j,2)=-kappa*z(i,j,2)*dh(i,j)*q2(i,j,2)
           vf(i,j,kb-1)=kappa*(1+z(i,j,kbm1))*dh(i,j)*q2(i,j,kbm1)
        end do
     end do
     !$omp end do

     !$omp do private(i,j,k)     
     do k=2,kbm1
        do j=2,jmm1
           do i=2,imm1
              dtef(i,j,k)=dtef(i,j,k) &
                   & *(1.d0+e2*((1.d0/abs(z(i,j,k)-z(i,j,1)) &
                   & +1.d0/abs(z(i,j,k)-z(i,j,kb)))*l(i,j,k)/(dh(i,j)*kappa))**2)
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel
     
     do k=3,kbm1
        do j=2,jmm1
           do i=2,imm1
              gg(i,j,k)=1.d0/(a(i,j,k)+c(i,j,k)*(1.d0-ee(i,j,k-1)) &
                   & -(dti2*dtef(i,j,k)+1.d0))
              ee(i,j,k)=a(i,j,k)*gg(i,j,k)
              gg(i,j,k)=(dti2*(-prod(i,j,k)*l(i,j,k)*e1) &
                   & +c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
           end do
        end do
     end do

     do k=kbm1,2,-1
        do j=2,jmm1
           do i=2,imm1
              vf(i,j,k)=ee(i,j,k)*vf(i,j,k+1)+gg(i,j,k)
           end do
        end do
     end do
     
  end if

  !     the following is to counter the problem of the ratio of two small
  !     numbers (l = q2l/q2) or one number becoming negative. Two options are
  !     included below. In this application, the second option, l was less
  !     noisy when uf or vf is small

  !$omp parallel
  !$omp do private(i,j,k)       
  do k=2,kbm1
     do j=2,jmm1
        do i=2,imm1
           uf(i,j,k)=abs(uf(i,j,k))
           vf(i,j,k)=abs(vf(i,j,k))
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  if(lmynnf) then
     
     !     gh = shear2
     gh(1:im,1:jm,1:kb)=0.d0

     !$omp parallel
     !$omp do private(i,j,k)       
     do k=2,kb
        do j=2,jmm1
           do i=2,imm1
              if(dzz(i,j,k-1)*dh(i,j) > 0.d0) then
                 gh(i,j,k)=0.25d0*((u(i,j,k)-u(i,j,k-1)+u(i+1,j,k)-u(i+1,j,k-1))**2 &
                      & +(v(i,j,k)-v(i,j,k-1)+v(i,j+1,k)-v(i,j+1,k-1))**2) &
                      & /(dzz(i,j,k-1)*dh(i,j))**2
              endif
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel
     
     call mynn_get_l(l,UF,boygr,wtsurf,wssurf,wusurf,wvsurf,l0,z,dzz,dh,im,jm,kb,ntp)
     !     use vf as storage for q2l2 - level 2 diagnostic of q2
     call mynn_get_q2l2(l,vf,boygr,gh,im,jm,kb)
     call mynn_get_ShSm(sh,sm,l,uf,vf,boygr,gh,im,jm,kb)

  else

     !     the following section solves for km and kh

     !     note that sm and sh limit to infinity when gh approaches 0.0288

     !$omp parallel
     !$omp do private(i,j,k,coef1,coef2,coef3)     
     do k=1,kb
        do j=2,jmm1
           do i=2,imm1
              coef1=a2*(1.d0-6.d0*a1/b1*stf(i,j,k))
              coef2=3.d0*a2*b2/stf(i,j,k)+18.d0*a1*a2
              coef3=a1*(1.d0-3.d0*c1-6.d0*a1/b1*stf(i,j,k))
              sh(i,j,k)=coef1/(1.d0-coef2*gh(i,j,k))
              sm(i,j,k)=coef3+sh(i,j,k)*coef4*gh(i,j,k)
              sm(i,j,k)=sm(i,j,k)/(1.d0-coef5*gh(i,j,k))
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel

  end if
  
  if(lmynnf) then
     
     !$omp parallel
     !$omp do private(i,j,k)            
     do k=1,kb
        do j=2,jmm1
           do i=2,imm1
              prod(i,j,k)=l(i,j,k)*sqrt(abs(uf(i,j,k)))
              kq(i,j,k)=(prod(i,j,k)*3.0d0*sm(i,j,k)+kq(i,j,k))*0.5d0
              km(i,j,k)=(prod(i,j,k)*sm(i,j,k)+km(i,j,k))*0.5d0
              kh(i,j,k)=(prod(i,j,k)*sh(i,j,k)+kh(i,j,k))*0.5d0
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  else

     !$omp parallel
     !$omp do private(i,j,k)                 
     do k=1,kb
        do j=2,jmm1
           do i=2,imm1
              prod(i,j,k)=l(i,j,k)*sqrt(abs(q2(i,j,k)))
              kq(i,j,k)=(prod(i,j,k)*0.41d0*sh(i,j,k)+kq(i,j,k))*0.5d0
              km(i,j,k)=(prod(i,j,k)*sm(i,j,k)+km(i,j,k))*0.5d0
              kh(i,j,k)=(prod(i,j,k)*sh(i,j,k)+kh(i,j,k))*0.5d0
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end if

  !     cosmetics: make boundr. values as interior (even if not used, printout
  !     may show strange values)

  if(n_north == -1)then
     km(1:im,jm,1:kb)=km(1:im,jmm1,1:kb)
     kh(1:im,jm,1:kb)=kh(1:im,jmm1,1:kb)
     kq(1:im,jm,1:kb)=kq(1:im,jmm1,1:kb)
  end if
  
  if(n_south == -1)then
     km(1:im,1,1:kb)=km(1:im,2,1:kb)
     kh(1:im,1,1:kb)=kh(1:im,2,1:kb)
     kq(1:im,1,1:kb)=kq(1:im,2,1:kb)
  end if

  if(n_east == -1)then
     km(im,1:jm,1:kb)=km(imm1,1:jm,1:kb)
     kh(im,1:jm,1:kb)=kh(imm1,1:jm,1:kb)
     kq(im,1:jm,1:kb)=kq(imm1,1:jm,1:kb)
  end if

  if(n_west == -1) then
     km(1,1:jm,1:kb)=km(2,1:jm,1:kb)
     kh(1,1:jm,1:kb)=kh(2,1:jm,1:kb)
     kq(1,1:jm,1:kb)=kq(2,1:jm,1:kb)
  end if

  !$omp parallel
  !$omp do private(i,j,k)                   
  do k=1,kb
     do j=1,jm
        do i=1,im
           km(i,j,k)=km(i,j,k)*fsm(i,j)
           kh(i,j,k)=kh(i,j,k)*fsm(i,j)
           kq(i,j,k)=kq(i,j,k)*fsm(i,j)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  call exchange3d_mpi(km,im,jm,kb) ! 2018.07
  call exchange3d_mpi(kh,im,jm,kb) ! 2018.07
  call exchange3d_mpi(kq,im,jm,kb) ! 2018.08
  call exchange3d_mpi(l,im,jm,kb) !2018.08
  
end subroutine profq_mynnf
