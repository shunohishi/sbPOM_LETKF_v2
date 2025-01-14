subroutine dens(si,ti,rhoo)

  ! calculate (density-1000.)/rhoref.
  ! see: Mellor, G.L., 1991, J. Atmos. Oceanic Tech., 609-611
  ! note: if pressure is not used in dens, buoyancy term (boygr) in
  ! subroutine profq must be changed (see note in subroutine profq)

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k
  real(kind = r_size) cr,p,rhor,sr,tr,tr2,tr3,tr4

  real(kind = r_size),intent(in)  :: si(im,jm,kb),ti(im,jm,kb)
  real(kind = r_size),intent(out) :: rhoo(im,jm,kb)

  rhoo(:,:,:)=0.d0

  !$omp parallel
  !$omp do private(i,j,k,tr,sr,tr2,tr3,tr4,p,rhor,cr)  
  do k=1,kbm1
     do j=1,jm
        do i=1,im

           tr=ti(i,j,k)+tbias
           sr=si(i,j,k)+sbias
           tr2=tr*tr
           tr3=tr2*tr
           tr4=tr3*tr

           ! approximate pressure in units of bars
           p=grav*rhoref*(-zz(i,j,k)* h(i,j))*1.d-5

           rhor=-0.157406d0+6.793952d-2*tr         &
                & -9.095290d-3*tr2+1.001685d-4*tr3 &
                & -1.120083d-6*tr4+6.536332d-9*tr4*tr

           rhor=rhor &
                & +(0.824493d0-4.0899d-3*tr+7.6438d-5*tr2-8.2467d-7*tr3+5.3875d-9*tr4)*sr &
                & +(-5.72466d-3+1.0227d-4*tr-1.6546d-6*tr2)*abs(sr)**1.5d0 &
                & +4.8314d-4*sr*sr

           cr=1449.1d0+0.0821d0*p+4.55d0*tr-0.045d0*tr2+1.34d0*(sr-35.d0)
           rhor=rhor+1.d5*p/(cr*cr)*(1.d0-2.d0*p/(cr*cr))

           rhoo(i,j,k)=rhor/rhoref*fsm(i,j)
           
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine dens
!_________________________________________________________________________________________

!-----------------------------------------------------------------------
! Potential Density (UNESCO, 1981) |
!-----------------------------------------------------------------------
!
!     Reference: A. E. Gill (1982), Atmosphere-Ocean Dynamics, 
!     Appendix Three, Properties of seawater.
!
!-----------------------------------------------------------------------
!
! temp (Potential Temperature): degree C
! sal:-
! pres: dbar
!-----------------------------------------------------------------------

subroutine unesco_potential_density(temp,sal,rho)

  use common_pom_var, only: im,jm,kb,r_size,tbias,sbias
  implicit none

  integer i,j,k
  
  real(kind = r_size) t,s,pres
  real(kind = r_size) ko,k0,kw
  real(kind = r_size) rho0,rhow

  real(kind = r_size),intent(in) :: temp(im,jm,kb),sal(im,jm,kb)  
  real(kind = r_size),intent(out) :: rho(im,jm,kb)

  rho(:,:,:)=0.d0

  !$omp parallel
  !$omp do private(i,j,k,t,s,pres,rhow,rho0,kw,k0,ko)
  do k=1,kb
     do j=1,jm
        do i=1,im

           t=temp(i,j,k)+tbias
           s=sal(i,j,k)+sbias
           pres=0.d0
           
           rhow=999.842594d0 &
                & +6.793952d-2*t &
                & -9.095290d-3*t*t &
                & +1.001685d-4*t*t*t &
                & -1.120083d-6*t*t*t*t &
                & +6.536332d-9*t*t*t*t*t
           
           rho0=rhow &
                & +s*(0.824493d0 &
                & -4.0899d-3*t &
                & +7.6438d-5*t*t &
                & -8.2467d-7*t*t*t &
                & +5.3875d-9*t*t*t*t) &
                & +s**1.5d0 *(-5.72466d-3 &
                & +1.0227d-4*t &
                & -1.6546d-6*t*t) &
                & +4.8314d-4*s*s 
           
           kw=19652.21d0 &
                & +148.4206d0*t &
                & -2.327105d0*t*t &
                & +1.360477d-2*t*t*t &
                & -5.155288d-5*t*t*t*t
           
           k0=kw &
                & +s*(54.6746d0 &
                & -0.603459d0*t &
                & +1.09987d-2*t*t &
                & -6.1670d-5*t*t*t) &
                & +s**1.5d0*(7.944d-2 &
                & +1.6483d-2*t -5.3009d-4*t*t)
           
           ko=k0 &
                & +pres*(3.239908d0 &
                & +1.43713d-3*t &
                & +1.16092d-4*t*t &
                & -5.77905d-7*t*t*t) &
                & +pres*s*(2.2838d-3 &
                & -1.0981d-5*t -1.6078d-6*t*t) &
                & +1.91075d-4*pres*s**1.5d0 &
                & +pres*pres*(8.50935d-5 &
                & -6.12293d-6*t &
                & +5.2787d-8*t*t) &
                & +pres*pres*s *(-9.9348d-7 &
                & +2.0816d-8*t +9.1697d-10*t*t)
           
           rho(i,j,k)=rho0*ko/(ko-pres)
           
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine unesco_potential_density
    
!_________________________________________________________________________________________
subroutine baropg
  ! calculate  baroclinic pressure gradient

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k

  !$omp parallel  
  !$omp do private(i,j,k)  
  do k=1,kbm1
     do j=1,jm
        do i=1,im
           rho(i,j,k)=rho(i,j,k)-rmean(i,j,k)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  call exchange3d_mpi(rho(:,:,1:kbm1),im,jm,kbm1) ! 2018.07.28
           
  ! calculate x-component of baroclinic pressure gradient
  
  !$omp parallel
  !$omp do private(i,j)  
  do j=2,jmm1
     do i=2,imm1
        drhox(i,j,1)=0.5d0*grav*(-zz(i,j,1))*(dt(i,j)+dt(i-1,j))*(rho(i,j,1)-rho(i-1,j,1))
     end do
  end do
  !$omp end do
  !$omp end parallel

  do k=2,kbm1
     do j=2,jmm1
        do i=2,imm1
           drhox(i,j,k)=drhox(i,j,k-1) &
                & +grav*0.25d0*(zz(i,j,k-1)-zz(i,j,k))*(dt(i,j)+dt(i-1,j)) &
                & *(rho(i,j,k)-rho(i-1,j,k)+rho(i,j,k-1)-rho(i-1,j,k-1)) &
                & +grav*.25d0*(zz(i,j,k-1)+zz(i,j,k))*(dt(i,j)-dt(i-1,j)) &
                & *(rho(i,j,k)+rho(i-1,j,k)-rho(i,j,k-1)-rho(i-1,j,k-1))
        end do
     end do
  end do

  !$omp parallel  
  !$omp do private(i,j,k)      
  do k=1,kbm1
     do j=2,jmm1
        do i=2,imm1
           drhox(i,j,k)=0.25d0*(dt(i,j)+dt(i-1,j))*drhox(i,j,k)*dum(i,j)*(dy(i,j)+dy(i-1,j))
        end do
     end do
  end do
  !$omp end do
  
  ! calculate y-component of baroclinic pressure gradient
  !$omp do private(i,j)        
  do j=2,jmm1
     do i=2,imm1
        drhoy(i,j,1)=0.5d0*grav*(-zz(i,j,1))*(dt(i,j)+dt(i,j-1))*(rho(i,j,1)-rho(i,j-1,1))
     end do
  end do
  !$omp end do
  !$omp end parallel

  do k=2,kbm1
     do j=2,jmm1
        do i=2,imm1
           drhoy(i,j,k)=drhoy(i,j,k-1) &
                & +grav*0.25d0*(zz(i,j,k-1)-zz(i,j,k))*(dt(i,j)+dt(i,j-1)) &
                & *(rho(i,j,k)-rho(i,j-1,k)+rho(i,j,k-1)-rho(i,j-1,k-1)) &
                & +grav*0.25d0*(zz(i,j,k-1)+zz(i,j,k))*(dt(i,j)-dt(i,j-1)) &
                & *(rho(i,j,k)+rho(i,j-1,k)-rho(i,j,k-1)-rho(i,j-1,k-1))
        end do
     end do
  end do

  !$omp parallel  
  !$omp do private(i,j,k)  
  do k=1,kbm1
     do j=2,jmm1
        do i=2,imm1
           drhoy(i,j,k)=0.25d0*(dt(i,j)+dt(i,j-1))*drhoy(i,j,k)*dvm(i,j)*(dx(i,j)+dx(i,j-1))
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)    
  do k=1,kb
     do j=2,jmm1
        do i=2,imm1
           drhox(i,j,k)=ramp*drhox(i,j,k)
           drhoy(i,j,k)=ramp*drhoy(i,j,k)
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)  
  do k=1,kbm1
     do j=1,jm
        do i=1,im
           rho(i,j,k)=rho(i,j,k)+rmean(i,j,k)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  call exchange3d_mpi(rho(:,:,1:kbm1),im,jm,kbm1) ! 2018.07.28

  
end subroutine baropg
!_______________________________________________________________________
subroutine baropg_thiem

  ! calculate baroclinic pressure gradient (drhox,drhoy) 
  ! following Thiem and Berntsen (2006; Ocean Modelling, 12, 140-156)
  ! 2018.05 Created by S.Ohishi

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k

  real(kind = r_size),parameter :: w1=3.d0,w2=1.d0 !Weightning coefficient(3-4:1)

  real(kind = r_size) ww2x(im,jm),ww2y(im,jm)

  !$omp parallel  
  !$omp do private(i,j,k)  
  do k=1,kbm1
     do j=1,jm
        do i=1,im
           rho(i,j,k)=rho(i,j,k)-rmean(i,j,k)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  call exchange3d_mpi(rho(:,:,1:kbm1),im,jm,kbm1) ! 2018.07.28
  
  !$omp parallel
  !$omp do private(i,j)  
  do j=2,jmm1
     do i=2,imm1

        ww2x(i,j)= &
             & w2*fsm(i,j)*fsm(i-1,j-1)*fsm(i-1,j)*fsm(i,j-1)*fsm(i,j+1)*fsm(i-1,j+1)

        ww2y(i,j)= &
             & w2*fsm(i,j)*fsm(i-1,j-1)*fsm(i-1,j)*fsm(i,j-1)*fsm(i+1,j)*fsm(i+1,j-1)

     end do
  end do
  !$omp end do

  !     k=1
  !$omp do private(i,j)
  do j=2,jmm1
     do i=2,imm1

        drhox(i,j,1)=grav*(-zz(i,j,1)) &
             & *(w1*0.5d0*(dt(i,j)+dt(i-1,j))*(rho(i,j,1)-rho(i-1,j,1)) & !(i,j)-(i-1,j)
             & +ww2x(i,j)*0.125d0 &
             & *((dt(i,j)+dt(i-1,j-1))*(rho(i,j,1)-rho(i-1,j-1,1))  & !(i,j)-(i-1,j-1)
             & -(dt(i-1,j)+dt(i,j-1))*(rho(i-1,j,1)-rho(i,j-1,1))   & !(i-1,j)-(i,j-1)
             & +(dt(i,j+1)+dt(i-1,j))*(rho(i,j+1,1)-rho(i-1,j,1))   & !(i,j+1)-(i-1,j)
             & -(dt(i-1,j+1)+dt(i,j))*(rho(i-1,j+1,1)-rho(i,j,1)))) & !(i-1,j+1)-(i,j)
             & /(w1+ww2x(i,j))

        drhoy(i,j,1)=grav*(-zz(i,j,1)) &
             & *(w1*0.5e0*(dt(i,j)+dt(i,j-1))*(rho(i,j,1)-rho(i,j-1,1)) & !(i,j)-(i,j-1)
             & +ww2y(i,j)*0.125d0 &
             & *((dt(i,j)+dt(i-1,j-1))*(rho(i,j,1)-rho(i-1,j-1,1)) & !(i,j)-(i-1,j-1) 
             & +(dt(i-1,j)+dt(i,j-1))*(rho(i-1,j,1)-rho(i,j-1,1))  & !(i-1,j)-(i,j-1)
             & +(dt(i+1,j)+dt(i,j-1))*(rho(i+1,j,1)-rho(i,j-1,1))  & !(i+1,j)-(i,j-1)
             & +(dt(i,j)+dt(i+1,j-1))*(rho(i,j,1)-rho(i+1,j-1,1))))& !(i,j)-(i+1,j-1)
             & /(w1+ww2y(i,j))
     end do
  end do
  !$omp end do
  !$omp end parallel

  !     k=2-kbm1 integration
  do k=2,kbm1
     do j=2,jmm1
        do i=2,imm1

           drhox(i,j,k)=drhox(i,j,k-1) &
                & +0.25d0*grav* &
                & (w1*((zz(i,j,k-1)-zz(i,j,k))*(dt(i,j)+dt(i-1,j)) & !(i,j)-(i-1,j)
                & *(rho(i,j,k-1)-rho(i-1,j,k-1)+rho(i,j,k)-rho(i-1,j,k))      &
                & -(zz(i,j,k-1)+zz(i,j,k))*(dt(i,j)-dt(i-1,j))     &
                & *(rho(i,j,k-1)+rho(i-1,j,k-1)-rho(i,j,k)-rho(i-1,j,k)))     &
                & +ww2x(i,j)*0.25d0 &
                & *(((zz(i,j,k-1)-zz(i,j,k))*(dt(i,j)+dt(i-1,j-1)) & !(i,j)-(i-1,j-1)
                & *(rho(i,j,k-1)-rho(i-1,j-1,k-1)+rho(i,j,k)-rho(i-1,j-1,k))  &
                & -(zz(i,j,k-1)+zz(i,j,k))*(dt(i,j)-dt(i-1,j-1))   &
                & *(rho(i,j,k-1)+rho(i-1,j-1,k-1)-rho(i,j,k)-rho(i-1,j-1,k))) &
                & -((zz(i,j,k-1)-zz(i,j,k))*(dt(i-1,j)+dt(i,j-1))  & !(i-1,j)-(i,j-1)
                & *(rho(i-1,j,k-1)-rho(i,j-1,k-1)+rho(i-1,j,k)-rho(i,j-1,k))  &
                & -(zz(i,j,k-1)+zz(i,j,k))*(dt(i-1,j)-dt(i,j-1))   &
                & *(rho(i-1,j,k-1)+rho(i,j-1,k-1)-rho(i-1,j,k)-rho(i,j-1,k))) &
                & +((zz(i,j,k-1)-zz(i,j,k))*(dt(i,j+1)+dt(i-1,j))  & !(i,j+1)-(i-1,j)
                & *(rho(i,j+1,k-1)-rho(i-1,j,k-1)+rho(i,j+1,k)-rho(i-1,j,k))  &
                & -(zz(i,j,k-1)+zz(i,j,k))*(dt(i,j+1)-dt(i-1,j))   &
                & *(rho(i,j+1,k-1)+rho(i-1,j,k-1)-rho(i,j+1,k)-rho(i-1,j,k))) &
                & -((zz(i,j,k-1)-zz(i,j,k))*(dt(i-1,j+1)+dt(i,j))  & !(i-1,j+1)-(i,j)
                & *(rho(i-1,j+1,k-1)-rho(i,j,k-1)+rho(i-1,j+1,k)-rho(i,j,k))  &
                & -(zz(i,j,k-1)+zz(i,j,k))*(dt(i-1,j+1)-dt(i,j))   &
                & *(rho(i-1,j+1,k-1)+rho(i,j,k-1)-rho(i-1,j+1,k)-rho(i,j,k))))) &
                & /(w1+ww2x(i,j))

           drhoy(i,j,k)=drhoy(i,j,k-1) &
                & +0.25d0*grav* &
                & (w1*((zz(i,j,k-1)-zz(i,j,k))*(dt(i,j)+dt(i,j-1)) & !(i,j)-(i,j-1)
                & *(rho(i,j,k-1)-rho(i,j-1,k-1)+rho(i,j,k)-rho(i,j-1,k))  &
                & -(zz(i,j,k-1)+zz(i,j,k))*(dt(i,j)-dt(i,j-1))     &
                & *(rho(i,j,k-1)+rho(i,j-1,k-1)-rho(i,j,k)-rho(i,j-1,k))) &
                & +ww2y(i,j)*0.25d0 &
                & *(((zz(i,j,k-1)-zz(i,j,k))*(dt(i,j)+dt(i-1,j-1)) & !(i,j)-(i-1,j-1)
                & *(rho(i,j,k-1)-rho(i-1,j-1,k-1)+rho(i,j,k)-rho(i-1,j-1,k))  &
                & -(zz(i,j,k-1)+zz(i,j,k))*(dt(i,j)-dt(i-1,j-1))   &
                & *(rho(i,j,k-1)+rho(i-1,j-1,k-1)-rho(i,j,k)-rho(i-1,j-1,k))) &
                & +((zz(i,j,k-1)-zz(i,j,k))*(dt(i-1,j)+dt(i,j-1))  & !(i-1,j)-(i,j-1)
                & *(rho(i-1,j,k-1)-rho(i,j-1,k-1)+rho(i-1,j,k)-rho(i,j-1,k))  &
                & -(zz(i,j,k-1)+zz(i,j,k))*(dt(i-1,j)-dt(i,j-1))   &
                & *(rho(i-1,j,k-1)+rho(i,j-1,k-1)-rho(i-1,j,k)-rho(i,j-1,k))) &
                & +((zz(i,j,k-1)-zz(i,j,k))*(dt(i+1,j)+dt(i,j-1))  & !(i+1,j)-(i,j-1)
                & *(rho(i+1,j,k-1)-rho(i,j-1,k-1)+rho(i+1,j,k)-rho(i,j-1,k))  &
                & -(zz(i,j,k-1)+zz(i,j,k))*(dt(i+1,j)-dt(i,j-1))   &
                & *(rho(i+1,j,k-1)+rho(i,j-1,k-1)-rho(i+1,j,k)-rho(i,j-1,k))) &
                & +((zz(i,j,k-1)-zz(i,j,k))*(dt(i,j)+dt(i+1,j-1))  & !(i,j)-(i+1,j-1)
                & *(rho(i,j,k-1)-rho(i+1,j-1,k-1)+rho(i,j,k)-rho(i+1,j-1,k))  &
                & -(zz(i,j,k-1)+zz(i,j,k))*(dt(i,j)-dt(i+1,j-1))   &
                & *(rho(i,j,k-1)+rho(i+1,j-1,k-1)-rho(i,j,k)-rho(i+1,j-1,k))))) &
                & /(w1+ww2y(i,j))

        end do
     end do
  end do

  !$omp parallel  
  !$omp do private(i,j,k)  
  do k=1,kbm1
     do j=2,jmm1
        do i=2,imm1

           drhox(i,j,k)=dum(i,j)*0.25d0*(dt(i,j)+dt(i-1,j)) & 
                & *drhox(i,j,k)*(dy(i,j)+dy(i-1,j))*ramp

           drhoy(i,j,k)=dvm(i,j)*0.25d0*(dt(i,j)+dt(i,j-1)) &
                & *drhoy(i,j,k)*(dx(i,j)+dx(i,j-1))*ramp

        end do
     end do
  end do
  !$omp end do
  !$omp end parallel


  !$omp parallel  
  !$omp do private(i,j,k)  
  do k=1,kbm1
     do j=1,jm
        do i=1,im
           rho(i,j,k)=rho(i,j,k)+rmean(i,j,k)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  call exchange3d_mpi(rho(:,:,1:kbm1),im,jm,kbm1) ! 2018.07.28
  
end subroutine baropg_thiem
!_______________________________________________________________________
subroutine baropg_mcc
  ! calculate  baroclinic pressure gradient
  ! 4th order correction terms, following McCalpin

  !$use omp_lib
  use common_pom_var
  implicit none

  integer i,j,k
  real(kind = r_size) d4(im,jm),ddx(im,jm),drho(im,jm,kb),rhou(im,jm,kb)
  real(kind = r_size) rho4th(0:im,0:jm,kb),d4th(0:im,0:jm)

  !$omp parallel  
  !$omp do private(i,j,k)  
  do k=1,kbm1
     do j=1,jm
        do i=1,im
           rho(i,j,k)=rho(i,j,k)-rmean(i,j,k)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
     
  call exchange3d_mpi(rho(:,:,1:kbm1),im,jm,kbm1) ! 2018.07.28
  
  ! convert a 2nd order matrices to special 4th order
  ! special 4th order case
  call order2d_mpi(d,d4th,im,jm)
  call order3d_mpi(rho,rho4th,im,jm,kb)
  
  ! compute terms correct to 4th order
  d4(1:im,1:jm)=0.d0
  ddx(1:im,1:jm)=0.d0

  rhou(1:im,1:jm,1:kb)=0.d0
  drho(1:im,1:jm,1:kb)=0.d0
  
  ! compute DRHO, RHOU, DDX and D4
  !$omp parallel
  !$omp do private(i,j,k)  
  do k=1,kbm1
     do j=1,jm
        do i=2,im
           rhou(i,j,k)=0.5d0*(rho(i,j,k)+rho(i-1,j,k))*dum(i,j)
           drho(i,j,k)=(rho(i,j,k)-rho(i-1,j,k))*dum(i,j)
        end do
     end do
  end do
  !$omp end do
  
  !$omp do private(i,j)  
  do j=1,jm
     do i=2,im
        ddx(i,j)=(d(i,j)-d(i-1,j))*dum(i,j)
        d4(i,j)=0.5d0*(d(i,j)+d(i-1,j))*dum(i,j)
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  !$omp parallel
  if(n_west == -1)then

     !$omp do private(i,j,k)  
     do k=1,kbm1
        do j=1,jm
           do i=3,imm1
              drho(i,j,k)=drho(i,j,k) &
                   & -(1.d0/24.d0)*(dum(i+1,j)*(rho(i+1,j,k)-rho(i,j,k))-2.d0*(rho(i,j,k)-rho(i-1,j,k)) &
                   & +dum(i-1,j)*(rho(i-1,j,k)-rho(i-2,j,k)))
              rhou(i,j,k)=rhou(i,j,k) &
                   & +(1.d0/16.d0)*(dum(i+1,j)*(rho(i,j,k)-rho(i+1,j,k)) &
                   & +dum(i-1,j)*(rho(i-1,j,k)-rho(i-2,j,k)))
           end do
        end do
     end do
     !$omp end do

     !$omp do private(i,j,k)     
     do j=1,jm
        do i=3,imm1
           ddx(i,j)=ddx(i,j) &
                & -(1.d0/24.d0)*(dum(i+1,j)*(d(i+1,j)-d(i,j))-2.d0*(d(i,j)-d(i-1,j)) &
                & +dum(i-1,j)*(d(i-1,j)-d(i-2,j)))
           d4(i,j)=d4(i,j) &
                & +(1.d0/16.d0)*(dum(i+1,j)*(d(i,j)-d(i+1,j))+dum(i-1,j)*(d(i-1,j)-d(i-2,j)))
        end do
     end do
     !$omp end do
     
  else

     !$omp do private(i,j,k)       
     do k=1,kbm1
        do j=1,jm
           do i=2,imm1
              drho(i,j,k)=drho(i,j,k) &
                   & -(1.d0/24.d0)*(dum(i+1,j)*(rho(i+1,j,k)-rho(i,j,k)) &
                   & -2.d0*(rho(i,j,k)-rho(i-1,j,k))+dum(i-1,j)*(rho(i-1,j,k)-rho4th(i-2,j,k)))
              rhou(i,j,k)=rhou(i,j,k) &
                   & +(1.d0/16.d0)*(dum(i+1,j)*(rho(i,j,k)-rho(i+1,j,k))+ &
                   & dum(i-1,j)*(rho(i-1,j,k)-rho4th(i-2,j,k)))
           end do
        end do
     end do
     !$omp end do

     !$omp do private(i,j)            
     do j=1,jm
        do i=2,imm1
           ddx(i,j)=ddx(i,j) &
                & -(1.d0/24.d0)*(dum(i+1,j)*(d(i+1,j)-d(i,j)) &
                & -2.d0*(d(i,j)-d(i-1,j))+dum(i-1,j)*(d(i-1,j)-d4th(i-2,j)))
           d4(i,j)=d4(i,j) &
                & +(1.d0/16.d0)*(dum(i+1,j)*(d(i,j)-d(i+1,j))+dum(i-1,j)*(d(i-1,j)-d4th(i-2,j)))
        end do
     end do
     !$omp end do
     
  end if
  !$omp end parallel
  
  ! calculate x-component of baroclinic pressure gradient
  !$omp parallel
  !$omp do private(i,j)       
  do i=2,imm1
     do j=2,jmm1
        drhox(i,j,1)=grav*(-zz(i,j,1))*d4(i,j)*drho(i,j,1)
     end do
  end do
  !$omp end do
  !$omp end parallel

  do k=2,kbm1
     do j=2,jmm1
        do i=2,imm1
           drhox(i,j,k)=drhox(i,j,k-1) &
                & +grav*0.5d0*dzz(i,j,k-1)*d4(i,j)*(drho(i,j,k-1)+drho(i,j,k)) &
                & +grav*0.5d0*(zz(i,j,k-1)+zz(i,j,k))*ddx(i,j)*(rhou(i,j,k)-rhou(i,j,k-1))
        end do
     end do
  end do

  !$omp parallel  
  !$omp do private(i,j,k)  
  do k=1,kbm1
     do j=2,jmm1
        do i=2,imm1
           drhox(i,j,k)=0.25d0*(dt(i,j)+dt(i-1,j))*drhox(i,j,k)*dum(i,j)*(dy(i,j)+dy(i-1,j))
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  ! compute terms correct to 4th order
  
  ddx(1:im,1:jm)=0.d0
  d4(1:im,1:jm)=0.d0

  rhou(1:im,1:jm,1:kb)=0.d0
  drho(1:im,1:jm,1:kb)=0.d0
  
  ! compute DRHO, RHOU, DDX and D4
  !$omp parallel
  !$omp do private(i,j,k)  
  do k=1,kbm1
     do j=2,jm
        do i=1,im
           drho(i,j,k)=(rho(i,j,k)-rho(i,j-1,k))*dvm(i,j)
           rhou(i,j,k)=0.5d0*(rho(i,j,k)+rho(i,j-1,k))*dvm(i,j)
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j)  
  do j=2,jm
     do i=1,im
        ddx(i,j)=(d(i,j)-d(i,j-1))*dvm(i,j)
        d4(i,j)=0.5d0*(d(i,j)+d(i,j-1))*dvm(i,j)
     end do
  end do
  !$omp end do
  
  if(n_south == -1) then

     !$omp do private(i,j,k)
     do k=1,kbm1
        do j=3,jmm1
           do i=1,im
              drho(i,j,k)=drho(i,j,k) &
                   & -(1.d0/24.d0)*(dvm(i,j+1)*(rho(i,j+1,k)-rho(i,j,k)) &
                   & -2.d0*(rho(i,j,k)-rho(i,j-1,k))+dvm(i,j-1)*(rho(i,j-1,k)-rho(i,j-2,k)))
              rhou(i,j,k)=rhou(i,j,k) &
                   & +(1.d0/16.d0)*(dvm(i,j+1)*(rho(i,j,k)-rho(i,j+1,k))+ &
                   & dvm(i,j-1)*(rho(i,j-1,k)-rho(i,j-2,k)))
           end do
        end do
     end do
     !$omp end do

     !$omp do private(i,j)     
     do j=3,jmm1
        do i=1,im
           ddx(i,j)=ddx(i,j) &
                & -(1.d0/24.d0)*(dvm(i,j+1)*(d(i,j+1)-d(i,j)) &
                & -2.d0*(d(i,j)-d(i,j-1))+dvm(i,j-1)*(d(i,j-1)-d(i,j-2)))
           d4(i,j)=d4(i,j) &
                & +(1.d0/16.d0)*(dvm(i,j+1)*(d(i,j)-d(i,j+1))+ &
                & dvm(i,j-1)*(d(i,j-1)-d(i,j-2)))
        end do
     end do
     !$omp end do
     
  else

     !$omp do private(i,j,k)
     do k=1,kbm1
        do j=2,jmm1
           do i=1,im
              drho(i,j,k)=drho(i,j,k) &
                   & -(1.d0/24.d0)*(dvm(i,j+1)*(rho(i,j+1,k)-rho(i,j,k)) &
                   & -2.d0*(rho(i,j,k)-rho(i,j-1,k))+dvm(i,j-1)*(rho(i,j-1,k)-rho4th(i,j-2,k)))
              rhou(i,j,k)=rhou(i,j,k) &
                   & +(1.d0/16.d0)*(dvm(i,j+1)*(rho(i,j,k)-rho(i,j+1,k))+ &
                   & dvm(i,j-1)*(rho(i,j-1,k)-rho4th(i,j-2,k)))
           end do
        end do
     end do
     !$omp end do

     !$omp do private(i,j)     
     do j=2,jmm1
        do i=1,im
           ddx(i,j)=ddx(i,j) &
                & -(1.d0/24.d0)*(dvm(i,j+1)*(d(i,j+1)-d(i,j)) &
                & -2.d0*(d(i,j)-d(i,j-1))+dvm(i,j-1)*(d(i,j-1)-d4th(i,j-2)))
           d4(i,j)=d4(i,j) &
                & +(1.d0/16.d0)*(dvm(i,j+1)*(d(i,j)-d(i,j+1)) &
                & +dvm(i,j-1)*(d(i,j-1)-d4th(i,j-2)))
        end do
     end do
     !$omp end do
     
  end if
  !$omp end parallel
  
  ! calculate y-component of baroclinic pressure gradient
  !$omp parallel
  !$omp do private(i,j)
  do j=2,jmm1
     do i=2,imm1
        drhoy(i,j,1)=grav*(-zz(i,j,1))*d4(i,j)*drho(i,j,1)
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  do k=2,kbm1
     do j=2,jmm1
        do i=2,imm1
           drhoy(i,j,k)=drhoy(i,j,k-1) &
                & +grav*0.5d0*dzz(i,j,k-1)*d4(i,j)*(drho(i,j,k-1)+drho(i,j,k)) &
                & +grav*0.5d0*(zz(i,j,k-1)+zz(i,j,k))*ddx(i,j)*(rhou(i,j,k)-rhou(i,j,k-1))
        end do
     end do
  end do

  !$omp parallel  
  !$omp do private(i,j,k)  
  do k=1,kbm1
     do j=2,jmm1
        do i=2,imm1
           drhoy(i,j,k)=0.25d0*(dt(i,j)+dt(i,j-1))*drhoy(i,j,k)*dvm(i,j)*(dx(i,j)+dx(i,j-1))
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)  
  do k=1,kb
     do j=2,jmm1
        do i=2,imm1
           drhox(i,j,k)=ramp*drhox(i,j,k)
           drhoy(i,j,k)=ramp*drhoy(i,j,k)
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)  
  do k=1,kbm1
     do j=1,jm
        do i=1,im
           rho(i,j,k)=rho(i,j,k)+rmean(i,j,k)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  call exchange3d_mpi(rho(:,:,1:kbm1),im,jm,kbm1) ! 2018.07.28
  
end subroutine baropg_mcc

!_______________________________________________________________________
subroutine advct

  ! calculate the horizontal portions of momentum advection well in
  ! advance of their use in advu and advv so that their vertical integrals
  ! (created in the main program) may be used in the external (2-D) mode
  ! calculation

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k
  
  real(kind = r_size) xflux(im,jm,kb),yflux(im,jm,kb)
  real(kind = r_size) curv(im,jm,kb)
  real(kind = r_size) dtaam

  curv(1:im,1:jm,1:kb)=0.d0
  advx(1:im,1:jm,1:kb)=0.d0
  xflux(1:im,1:jm,1:kb)=0.d0
  yflux(1:im,1:jm,1:kb)=0.d0

  !$omp parallel  
  !$omp do private(i,j,k)  
  do k=1,kbm1
     do j=2,jmm1
        do i=2,imm1
           curv(i,j,k)= &
                & 0.25d0*((v(i,j+1,k)+v(i,j,k))*(dy(i+1,j)-dy(i-1,j)) &
                & -(u(i+1,j,k)+u(i,j,k))*(dx(i,j+1)-dx(i,j-1))) &
                & /(dx(i,j)*dy(i,j))
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  call exchange3d_mpi(curv(:,:,1:kbm1),im,jm,kbm1)
  
  ! calculate x-component of velocity advection

  ! calculate horizontal advective fluxes
  !$omp parallel
  !$omp do private(i,j,k)
  do k=1,kbm1
     do j=1,jm
        do i=2,imm1
           xflux(i,j,k)=&
                & 0.125d0*((dt(i+1,j)+dt(i,j))*u(i+1,j,k)+(dt(i,j)+dt(i-1,j))*u(i,j,k)) &
                & *(u(i+1,j,k)+u(i,j,k))
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)  
  do k=1,kbm1
     do j=2,jm
        do i=2,im
           yflux(i,j,k)= &
                & 0.125d0*((dt(i,j)+dt(i,j-1))*v(i,j,k)+(dt(i-1,j)+dt(i-1,j-1))*v(i-1,j,k)) &
                & *(u(i,j,k)+u(i,j-1,k))
        end do
     end do
  end do
  !$omp end do
  
  ! add horizontal diffusive fluxes
  !$omp do private(i,j,k,dtaam) 
  do k=1,kbm1
     do j=2,jm
        do i=2,imm1
           xflux(i,j,k)= &
                & xflux(i,j,k)-dt(i,j)*aam(i,j,k)*2.d0*(ub(i+1,j,k)-ub(i,j,k))/dx(i,j)
           dtaam= &
                & 0.25d0*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1)) &
                & *(aam(i,j,k)+aam(i-1,j,k)+aam(i,j-1,k)+aam(i-1,j-1,k))
           yflux(i,j,k)= &
                & yflux(i,j,k)-dtaam*((ub(i,j,k)-ub(i,j-1,k))/(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1)) &
                & +(vb(i,j,k)-vb(i-1,j,k))/(dx(i,j)+dx(i-1,j) &
                & +dx(i,j-1)+dx(i-1,j-1)))

           xflux(i,j,k)=dy(i,j)*xflux(i,j,k)
           yflux(i,j,k)=0.25d0*(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1))*yflux(i,j,k)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  call exchange3d_mpi(xflux(:,:,1:kbm1),im,jm,kbm1)
  call exchange3d_mpi(yflux(:,:,1:kbm1),im,jm,kbm1)

  ! do horizontal advection

  !$omp parallel
  !$omp do private(i,j,k)    
  do k=1,kbm1
     do j=2,jmm1
        do i=2,imm1
           advx(i,j,k)= &
                & xflux(i,j,k)-xflux(i-1,j,k)+yflux(i,j+1,k)-yflux(i,j,k)
        end do
     end do
  end do
  !$omp end do
  
  if(n_west == -1) then

     !$omp do private(i,j,k)
     do k=1,kbm1
        do j=2,jmm1
           do i=3,imm1
              advx(i,j,k)=advx(i,j,k) &
                   & -aru(i,j)*0.25d0*(curv(i,j,k)*dt(i,j)*(v(i,j+1,k)+v(i,j,k)) &
                   & +curv(i-1,j,k)*dt(i-1,j)*(v(i-1,j+1,k)+v(i-1,j,k)))
           end do
        end do
     end do
     !$omp end do
     
  else

     !$omp do private(i,j,k)
     do k=1,kbm1
        do j=2,jmm1
           do i=2,imm1
              advx(i,j,k)=advx(i,j,k) &
                   & -aru(i,j)*0.25d0*(curv(i,j,k)*dt(i,j)*(v(i,j+1,k)+v(i,j,k)) &
                   & +curv(i-1,j,k)*dt(i-1,j)*(v(i-1,j+1,k)+v(i-1,j,k)))
           end do
        end do
     end do
     !$omp end do
     
  end if
  !$omp end parallel
  
! calculate y-component of velocity advection

  advy(1:im,1:jm,1:kb)=0.d0
  xflux(1:im,1:jm,1:kb)=0.d0
  yflux(1:im,1:jm,1:kb)=0.d0
  
  ! calculate horizontal advective fluxes
  !$omp parallel
  !$omp do private(i,j,k)
  do k=1,kbm1
     do j=2,jm
        do i=2,im
           xflux(i,j,k)= &
                & 0.125d0*((dt(i,j)+dt(i-1,j))*u(i,j,k)+(dt(i,j-1)+dt(i-1,j-1))*u(i,j-1,k)) &
                & *(v(i,j,k)+v(i-1,j,k))
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)
  do k=1,kbm1
     do j=2,jmm1
        do i=1,im
           yflux(i,j,k)= &
                & 0.125d0*((dt(i,j+1)+dt(i,j))*v(i,j+1,k)+(dt(i,j)+dt(i,j-1))*v(i,j,k)) &
                & *(v(i,j+1,k)+v(i,j,k))
        end do
     end do
  end do
  !$omp end do

  ! add horizontal diffusive fluxes
  !$omp do private(i,j,k,dtaam)
  do k=1,kbm1
     do j=2,jmm1
        do i=2,im
           dtaam= &
                & 0.25d0*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1)) &
                & *(aam(i,j,k)+aam(i-1,j,k)+aam(i,j-1,k)+aam(i-1,j-1,k))
           xflux(i,j,k)= &
                & xflux(i,j,k)-dtaam*((ub(i,j,k)-ub(i,j-1,k))/(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1)) &
                & +(vb(i,j,k)-vb(i-1,j,k))/(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1)))
           yflux(i,j,k)= &
                & yflux(i,j,k)-dt(i,j)*aam(i,j,k)*2.d0*(vb(i,j+1,k)-vb(i,j,k))/dy(i,j)

           xflux(i,j,k)=0.25d0*(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))*xflux(i,j,k)
           yflux(i,j,k)=dx(i,j)*yflux(i,j,k)
           
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  call exchange3d_mpi(xflux(:,:,1:kbm1),im,jm,kbm1)  
  call exchange3d_mpi(yflux(:,:,1:kbm1),im,jm,kbm1)

  ! do horizontal advection
  !$omp parallel
  !$omp do private(i,j,k)    
  do k=1,kbm1
     do j=2,jmm1
        do i=2,imm1
           advy(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)+yflux(i,j,k)-yflux(i,j-1,k)
        end do
     end do
  end do
  !$omp end do

  if(n_south == -1) then
     !$omp do private(i,j,k)    
     do k=1,kbm1
        do j=3,jmm1
           do i=2,imm1
              advy(i,j,k)= &
                   & advy(i,j,k)+arv(i,j)*0.25d0*(curv(i,j,k)*dt(i,j)*(u(i+1,j,k)+u(i,j,k)) &
                   & +curv(i,j-1,k)*dt(i,j-1)*(u(i+1,j-1,k)+u(i,j-1,k)))
           end do
        end do
     end do
     !$omp end do
  else
     !$omp do private(i,j,k)    
     do k=1,kbm1
        do j=2,jmm1
           do i=2,imm1              
              advy(i,j,k)= &
                   & advy(i,j,k)+arv(i,j)*.25d0*(curv(i,j,k)*dt(i,j)*(u(i+1,j,k)+u(i,j,k)) &
                   & +curv(i,j-1,k)*dt(i,j-1)*(u(i+1,j-1,k)+u(i,j-1,k)))
           end do
        end do
     end do
     !$omp end do
  end if
  !$omp end parallel
  
end subroutine advct
!_______________________________________________________________________
subroutine advave

  ! calculate horizontal advection and diffusion

  !$use omp_lib
  use common_pom_var
  implicit none

  integer i,j

  real(kind = r_size) curv2d(im,jm)
  
  ! u-advection and diffusion
  
  ! advective fluxes
  advua(1:im,1:jm)=0.d0

  !$omp parallel  
  !$omp do private(i,j)       
  do j=2,jm
     do i=2,imm1
        fluxua(i,j)= &
             & 0.125d0*((d(i+1,j)+d(i,j))*ua(i+1,j)+(d(i,j)+d(i-1,j))*ua(i,j)) &
             & *(ua(i+1,j)+ua(i,j))
     end do
  end do
  !$omp end do

  !$omp do private(i,j)         
  do j=2,jm
     do i=2,im
        fluxva(i,j)= &
             & 0.125d0*((d(i,j)+d(i,j-1))*va(i,j)+(d(i-1,j)+d(i-1,j-1))*va(i-1,j)) &
             & *(ua(i,j)+ua(i,j-1))
     end do
  end do
  !$omp end do

  ! add viscous fluxes
  !$omp do private(i,j)
  do j=2,jm
     do i=2,imm1
        fluxua(i,j)=fluxua(i,j) &
             & -d(i,j)*2.d0*aam2d(i,j)*(uab(i+1,j)-uab(i,j))/dx(i,j)
     end do
  end do
  !$omp end do

  !$omp do private(i,j)  
  do j=2,jm
     do i=2,im
        tps(i,j)= &
             & 0.25d0*(d(i,j)+d(i-1,j)+d(i,j-1)+d(i-1,j-1)) &
             & *(aam2d(i,j)+aam2d(i,j-1)+aam2d(i-1,j)+aam2d(i-1,j-1)) &
             & *((uab(i,j)-uab(i,j-1))/(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1)) &
             & +(vab(i,j)-vab(i-1,j))/(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1)))
        fluxua(i,j)=fluxua(i,j)*dy(i,j)
        fluxva(i,j)=(fluxva(i,j)-tps(i,j)) &
             & *0.25d0*(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1))
     end do
  end do
  !$omp end do
  !$omp end parallel

  call exchange2d_mpi(fluxua,im,jm)

  !$omp parallel
  !$omp do private(i,j)
  do j=2,jmm1
     do i=2,imm1
        advua(i,j)=fluxua(i,j)-fluxua(i-1,j)+fluxva(i,j+1)-fluxva(i,j)
     end do
  end do
  !$omp end do
  !$omp end parallel

  call exchange2d_mpi(advua,im,jm) !2018.08
  
  ! v-advection and diffusion
  advva(1:im,1:jm)=0.d0

  ! advective fluxes
  !$omp parallel
  !$omp do private(i,j)
  do j=2,jm
     do i=2,im
        fluxua(i,j)= &
             & 0.125d0*((d(i,j)+d(i-1,j))*ua(i,j) &
             & +(d(i,j-1)+d(i-1,j-1))*ua(i,j-1))*(va(i-1,j)+va(i,j))
     end do
  end do
  !$omp end do

  !$omp do private(i,j)  
  do j=2,jmm1
     do i=2,im
        fluxva(i,j)= &
             & 0.125d0*((d(i,j+1)+d(i,j))*va(i,j+1) &
             & +(d(i,j)+d(i,j-1))*va(i,j))*(va(i,j+1)+va(i,j))
     end do
  end do
  !$omp end do

  ! add viscous fluxes
  !$omp do private(i,j)
  do j=2,jmm1
     do i=2,im
        fluxva(i,j)=fluxva(i,j) &
             & -d(i,j)*2.d0*aam2d(i,j)*(vab(i,j+1)-vab(i,j))/dy(i,j)
     end do
  end do
  !$omp end do

  !$omp do private(i,j)  
  do j=2,jm
     do i=2,im
        fluxva(i,j)=fluxva(i,j)*dx(i,j)
        fluxua(i,j)=(fluxua(i,j)-tps(i,j)) &
             & *0.25d0*(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  call exchange2d_mpi(fluxva,im,jm)

  !$omp parallel
  !$omp do private(i,j)
  do j=2,jmm1
     do i=2,imm1
        advva(i,j)=fluxua(i+1,j)-fluxua(i,j)+fluxva(i,j)-fluxva(i,j-1)
     end do
  end do
  !$omp end do
  !$omp end parallel

  call exchange2d_mpi(advva,im,jm) !2018.08
  
  if(mode == 2) then

     !$omp parallel
     !$omp do private(i,j)
     do j=2,jmm1
        do i=2,imm1
           wubot(i,j)= &
                & -0.5d0*(cbc(i,j)+cbc(i-1,j))*sqrt(uab(i,j)**2 &
                & +(0.25d0*(vab(i,j)+vab(i,j+1)+vab(i-1,j)+vab(i-1,j+1)))**2) &
                & *uab(i,j)
        end do
     end do
     !$omp end do

     !$omp do private(i,j)     
     do j=2,jmm1
        do i=2,imm1
           wvbot(i,j)= &
                & -0.5d0*(cbc(i,j)+cbc(i,j-1))*sqrt(vab(i,j)**2 &
                & +(0.25d0*(uab(i,j)+uab(i+1,j)+uab(i,j-1)+uab(i+1,j-1)))**2) &
                & *vab(i,j)
        end do
     end do
     !$omp end do

     !$omp do private(i,j)     
     do j=2,jmm1
        do i=2,imm1
           curv2d(i,j)= &
                & 0.25d0*((va(i,j+1)+va(i,j))*(dy(i+1,j)-dy(i-1,j)) &
                & -(ua(i+1,j)+ua(i,j))*(dx(i,j+1)-dx(i,j-1)))/(dx(i,j)*dy(i,j))
        end do
     end do
     !$omp end do
     !$omp end parallel

     call exchange2d_mpi(wubot,im,jm) !2018.07.28     
     call exchange2d_mpi(wvbot,im,jm) !2018.07.28     
     call exchange2d_mpi(curv2d,im,jm)


     !$omp parallel
     if(n_west == -1)then

        !$omp do private(i,j)
        do j=2,jmm1
           do i=3,imm1
              advua(i,j)=advua(i,j) &
                   & -aru(i,j)*0.25d0*(curv2d(i,j)*d(i,j)*(va(i,j+1)+va(i,j)) &
                   & +curv2d(i-1,j)*d(i-1,j)*(va(i-1,j+1)+va(i-1,j)))
           end do
        end do
        !$omp end do
        
     else        

        !$omp do private(i,j)        
        do j=2,jmm1
           do i=2,imm1
              advua(i,j)=advua(i,j) &
                   & -aru(i,j)*0.25d0*(curv2d(i,j)*d(i,j)*(va(i,j+1)+va(i,j)) &
                   & +curv2d(i-1,j)*d(i-1,j)*(va(i-1,j+1)+va(i-1,j)))
           end do
        end do
        !$omp end do

     end if
     
     if(n_south == -1)then

        !$omp do private(i,j)        
        do j=3,jmm1
           do i=2,imm1
              advva(i,j)=advva(i,j) &
                   & +arv(i,j)*0.25d0*(curv2d(i,j)*d(i,j)*(ua(i+1,j)+ua(i,j)) &
                   & +curv2d(i,j-1)*d(i,j-1)*(ua(i+1,j-1)+ua(i,j-1)))
           end do
        end do
        !$omp end do

     else

        !$omp do private(i,j)        
        do j=2,jmm1
           do i=2,imm1
              advva(i,j)=advva(i,j) &
                   & +arv(i,j)*0.25d0*(curv2d(i,j)*d(i,j)*(ua(i+1,j)+ua(i,j)) &
                   & +curv2d(i,j-1)*d(i,j-1)*(ua(i+1,j-1)+ua(i,j-1)))
           end do
        end do
        !$omp end do
           
     end if
     !$omp end parallel
     
     call exchange2d_mpi(advua,im,jm) !2018.08     
     call exchange2d_mpi(advva,im,jm) !2018.08
     
  endif
  
end subroutine advave
!_______________________________________________________________________
subroutine advq(qb,q,qf)

  ! calculates horizontal advection and diffusion, and vertical advection
  ! for turbulent quantities

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k
  real(kind = r_size) xflux(im,jm,kb),yflux(im,jm,kb)
  
  real(kind = r_size),intent(in) :: qb(im,jm,kb),q(im,jm,kb)
  real(kind = r_size),intent(out) :: qf(im,jm,kb)

  ! do horizontal advection
  !$omp parallel
  !$omp do private(i,j,k)
  do k=2,kbm1
     do j=2,jm
        do i=2,im
           xflux(i,j,k)=0.125d0*(q(i,j,k)+q(i-1,j,k)) &
                & *(dt(i,j)+dt(i-1,j))*(u(i,j,k)+u(i,j,k-1))
           yflux(i,j,k)=0.125d0*(q(i,j,k)+q(i,j-1,k)) &
                & *(dt(i,j)+dt(i,j-1))*(v(i,j,k)+v(i,j,k-1))
        end do
     end do
  end do
  !$omp end do
  
  ! do horizontal diffusion
  !$omp do private(i,j,k)
  do k=2,kbm1
     do j=2,jm
        do i=2,im
           xflux(i,j,k)=xflux(i,j,k) &
                & -0.25d0*(aam(i,j,k)+aam(i-1,j,k)+aam(i,j,k-1)+aam(i-1,j,k-1)) &
                & *(h(i,j)+h(i-1,j))*(qb(i,j,k)-qb(i-1,j,k))*dum(i,j) &
                & /(dx(i,j)+dx(i-1,j))
           yflux(i,j,k)=yflux(i,j,k) &
                & -0.25d0*(aam(i,j,k)+aam(i,j-1,k)+aam(i,j,k-1)+aam(i,j-1,k-1)) &
                & *(h(i,j)+h(i,j-1))*(qb(i,j,k)-qb(i,j-1,k))*dvm(i,j) &
                & /(dy(i,j)+dy(i,j-1))
           xflux(i,j,k)=0.5d0*(dy(i,j)+dy(i-1,j))*xflux(i,j,k)
           yflux(i,j,k)=0.5d0*(dx(i,j)+dx(i,j-1))*yflux(i,j,k)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  ! kii2b
  call exchange3d_mpi(xflux(:,:,1:kbm1),im,jm,kbm1)
  call exchange3d_mpi(yflux(:,:,1:kbm1),im,jm,kbm1)

  ! do vertical advection, add flux terms, then step forward in time
  !$omp parallel
  !$omp do private(i,j,k)
  do k=2,kbm1
     do j=2,jmm1
        do i=2,imm1
           qf(i,j,k)= &
                & (w(i,j,k-1)*q(i,j,k-1)-w(i,j,k+1)*q(i,j,k+1))*art(i,j)/(dz(i,j,k)+dz(i,j,k-1)) &
                & +xflux(i+1,j,k)-xflux(i,j,k)+yflux(i,j+1,k)-yflux(i,j,k)
           qf(i,j,k)=((h(i,j)+etb(i,j))*art(i,j)*qb(i,j,k)-dti2*qf(i,j,k)) &
                & /((h(i,j)+etf(i,j))*art(i,j))
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine advq
!_______________________________________________________________________
subroutine vertvl
  ! calculates vertical velocity

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k
  real(kind = r_size) xflux(im,jm,kb),yflux(im,jm,kb)

  ! reestablish boundary conditions
  !$omp parallel
  !$omp do private(i,j,k)
  do k=1,kbm1
     do j=2,jm
        do i=2,im
           xflux(i,j,k)=0.25d0*(dy(i,j)+dy(i-1,j))*(dt(i,j)+dt(i-1,j))*u(i,j,k)
           yflux(i,j,k)=0.25d0*(dx(i,j)+dx(i,j-1))*(dt(i,j)+dt(i,j-1))*v(i,j,k)
        end do
     end do
  end do
  !$omp end do
  
  ! note: if one wishes to include freshwater flux, the surface velocity
  ! should be set to vflux(i,j). See also change made to 2-D volume
  ! conservation equation which calculates elf
  !$omp do private(i,j)
  do j=2,jmm1
     do i=2,imm1
        w(i,j,1)=0.5d0*(vfluxb(i,j)+vfluxf(i,j))
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  do k=1,kbm1
     do j=2,jmm1
        do i=2,imm1
           w(i,j,k+1)=w(i,j,k) &
                & +dz(i,j,k)*((xflux(i+1,j,k)-xflux(i,j,k)+yflux(i,j+1,k)-yflux(i,j,k))/(dx(i,j)*dy(i,j)) &
                & +(etf(i,j)-etb(i,j))/dti2)
        end do
     end do     
  end do
  
end subroutine vertvl
!_______________________________________________________________________
subroutine profq
  ! solve for q2 (twice the turbulent kinetic energy), q2l (q2 x turbulent
  ! length scale), km (vertical kinematic viscosity) and kh (vertical
  ! kinematic diffusivity), using a simplified version of the level 2 1/2
  ! model of Mellor and Yamada (1982)
  ! in this version, the Craig-Banner sub-model whereby breaking wave tke
  ! is injected into the surface is included. However, we use an
  ! analytical solution to the near surface tke equation to solve for q2
  ! at the surface giving the same result as C-B diffusion. The new scheme
  ! is simpler and more robust than the latter scheme

  !$use omp_lib  
  use common_pom_var
  implicit none

  real(kind = r_size),parameter :: a1=0.92d0,b1=16.6d0,a2=0.74d0,b2=10.1d0,c1=0.08d0
  real(kind = r_size),parameter :: e1=1.8d0,e2=1.33d0
  real(kind = r_size),parameter :: sef=1.d0,cbcnst=100.d0,surfl=2.d5,shiw=0.0d0

  ! surface and bottom boundary conditions
  real(kind = r_size),parameter :: const1=(16.6d0**(2.d0/3.d0))*sef

  real(kind = r_size),parameter :: ghc=-6.0d0

  real(kind = r_size),parameter :: coef4=18.d0*a1*a1+9.d0*a1*a2
  real(kind = r_size),parameter :: coef5=9.d0*a1*a2

  integer i,j,k

  real(kind = r_size) a(im,jm,kb),c(im,jm,kb)
  real(kind = r_size) ee(im,jm,kb),gg(im,jm,kb)
  real(kind = r_size) sm(im,jm,kb),sh(im,jm,kb)
  real(kind = r_size) cc(im,jm,kb)
  real(kind = r_size) gh(im,jm,kb),boygr(im,jm,kb),dh(im,jm),stf(im,jm,kb)
  real(kind = r_size) prod(im,jm,kb)
  real(kind = r_size) coef1,coef2,coef3
  real(kind = r_size) l0(im,jm)
  real(kind = r_size) utau2(im,jm)
  real(kind = r_size) p,sp,tp

  a(1:im,1:jm,1:kb)=0.d0
  c(1:im,1:jm,1:kb)=0.d0
  cc(1:im,1:jm,1:kb)=0.d0
  ee(1:im,1:jm,1:kb)=0.d0
  gg(1:im,1:jm,1:kb)=0.d0
  l0(1:im,1:jm)=0.d0
  boygr(1:im,1:jm,1:kb)=0.d0
  prod(1:im,1:jm,1:kb)=0.d0
  
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
     do j=1,jm
        do i=1,im
           a(i,j,k)=-dti2*(kq(i,j,k+1)+kq(i,j,k)+2.d0*umol)*0.5d0 &
                & /(dzz(i,j,k-1)*dz(i,j,k)*dh(i,j)*dh(i,j))
           c(i,j,k)=-dti2*(kq(i,j,k-1)+kq(i,j,k)+2.d0*umol)*0.5d0 &
                & /(dzz(i,j,k-1)*dz(i,j,k-1)*dh(i,j)*dh(i,j))
        end do
     end do
  end do
  !$omp end do

  ! the following section solves the equation:
  !     dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b

  ! initialize fields that are not calculated on all boundaries
  ! but are later used there
  !$omp do private(i,j)  
  do j=1,jmm1
     do i=1,imm1
        utau2(i,j)=sqrt((0.5d0*(wusurf(i,j)+wusurf(i+1,j)))**2 &
             & +(0.5d0*(wvsurf(i,j)+wvsurf(i,j+1)))**2)
        uf(i,j,kb)=sqrt((0.5d0*(wubot(i,j)+wubot(i+1,j)))**2 &
             & +(0.5d0*(wvbot(i,j)+wvbot(i,j+1)))**2)*const1
     end do
  end do
  !$omp end do
  !$omp end parallel

  call exchange2d_mpi(utau2,im,jm)
  call exchange2d_mpi(uf(:,:,kb),im,jm)

  !$omp parallel
  !$omp do private(i,j)  
  do j=1,jm
     do i=1,im
        ! wave breaking energy- a variant of Craig & Banner (1994)
        ! see Mellor and Blumberg, 2003.
        ee(i,j,1)=0.d0
        gg(i,j,1)=(15.8d0*cbcnst)**(2.d0/3.d0)*utau2(i,j)
        ! surface length scale following Stacey (1999).
        l0(i,j)=surfl*utau2(i,j)/grav
     end do
  end do
  !$omp end do

  ! calculate speed of sound squared
  !$omp do private(i,j,k,tp,sp,p)
  do k=1,kbm1
     do j=1,jm
        do i=1,im
           tp=t(i,j,k)+tbias
           sp=s(i,j,k)+sbias
           ! calculate pressure in units of decibars
           p=grav*rhoref*(-zz(i,j,k)*h(i,j))*1.d-4
           cc(i,j,k)=1449.1d0+0.00821d0*p+4.55d0*tp-0.045d0*tp**2+1.34d0*(sp-35.0d0)
           cc(i,j,k)=cc(i,j,k)/sqrt((1.d0-0.01642d0*p/cc(i,j,k))*(1.d0-0.4d0*p/cc(i,j,k)**2))
        end do
     end do
  end do
  !$omp end do

  ! calculate buoyancy gradient
  !$omp do private(i,j,k)
  do k=2,kbm1
     do j=1,jm
        do i=1,im
           q2b(i,j,k)=abs(q2b(i,j,k))
           q2lb(i,j,k)=abs(q2lb(i,j,k))
           boygr(i,j,k)=grav*(rho(i,j,k-1)-rho(i,j,k))/(dzz(i,j,k-1)*h(i,j)) &
                ! *** note: comment out next line if dens does not include pressure
                & +(grav**2)*2.d0/(cc(i,j,k-1)**2+cc(i,j,k)**2)
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)  
  do k=2,kbm1
     do j=1,jm
        do i=1,im
           l(i,j,k)=abs(q2lb(i,j,k)/q2b(i,j,k))
           if(z(i,j,k) > -0.5) l(i,j,k)=max(l(i,j,k),kappa*l0(i,j))
           gh(i,j,k)=(l(i,j,k)**2)*boygr(i,j,k)/q2b(i,j,k)
           gh(i,j,k)=min(gh(i,j,k),0.028d0)
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j)    
  do j=1,jm
     do i=1,im
        l(i,j,1)=kappa*l0(i,j)
        l(i,j,kb)=0.d0
        gh(i,j,1)=0.d0
        gh(i,j,kb)=0.d0
     end do
  end do
  !$omp end do

  ! calculate production of turbulent kinetic energy:
  !$omp do private(i,j,k)    
  do k=2,kbm1
     do j=2,jmm1
        do i=2,imm1
           prod(i,j,k)=km(i,j,k)*0.25d0*sef &
                & *((u(i,j,k)-u(i,j,k-1)+u(i+1,j,k)-u(i+1,j,k-1))**2 &
                & +(v(i,j,k)-v(i,j,k-1)+v(i,j+1,k)-v(i,j+1,k-1))**2) &
                & /(dzz(i,j,k-1)*dh(i,j))**2 &
                                ! add shear due to internal wave field
                & -shiw*km(i,j,k)*boygr(i,j,k)
           prod(i,j,k)=prod(i,j,k)+kh(i,j,k)*boygr(i,j,k)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  call exchange3d_mpi(prod(:,:,2:kbm1),im,jm,kbm2)

  ! note: Richardson # dep. dissipation correction (Mellor, 2001; Ezer,
  ! 2000), depends on ghc the critical number (empirical -6 to -2) to
  ! increase mixing

  !$omp parallel
  !$omp do private(i,j,k)
  do k=1,kb
     do j=1,jm
        do i=1,im
           stf(i,j,k)=1.d0
           ! It is unclear yet if diss. corr. is needed when surf. waves are included.
           !           if(gh(i,j,k).lt.0.e0)
           !    $        stf(i,j,k)=1.0e0-0.9e0*(gh(i,j,k)/ghc)**1.5e0
           !           if(gh(i,j,k).lt.ghc) stf(i,j,k)=0.1e0
           dtef(i,j,k)=sqrt(abs(q2b(i,j,k)))*stf(i,j,k)/(b1*l(i,j,k)+small)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  do k=2,kbm1
     do j=1,jm
        do i=1,im
           gg(i,j,k)=1.d0/(a(i,j,k)+c(i,j,k)*(1.d0-ee(i,j,k-1))-(2.d0*dti2*dtef(i,j,k)+1.d0))
           ee(i,j,k)=a(i,j,k)*gg(i,j,k)
           gg(i,j,k)=(-2.d0*dti2*prod(i,j,k)+c(i,j,k)*gg(i,j,k-1)-uf(i,j,k))*gg(i,j,k)
        end do
     end do
  end do

  do k=kbm1,1,-1
     do j=1,jm
        do i=1,im
           uf(i,j,k)=ee(i,j,k)*uf(i,j,k+1)+gg(i,j,k)
        end do
     end do
  end do

  ! the following section solves the equation:
  !     dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb

  !$omp parallel  
  !$omp do private(i,j)    
  do j=1,jm
     do i=1,im
        vf(i,j,1)=0.d0
        vf(i,j,kb)=0.d0
        ee(i,j,2)=0.d0
        gg(i,j,2)=-kappa*z(i,j,2)*dh(i,j)*q2(i,j,2)
        vf(i,j,kb-1)=kappa*(1+z(i,j,kbm1))*dh(i,j)*q2(i,j,kbm1)
     end do
  end do
  !$omp end do
  
  !$omp do private(i,j,k)
  do k=2,kbm1
     do j=1,jm
        do i=1,im
           dtef(i,j,k)=dtef(i,j,k) &
                & *(1.d0+e2*((1.d0/abs(z(i,j,k)-z(i,j,1))+1.d0/abs(z(i,j,k)-z(i,j,kb))) &
                & *l(i,j,k)/(dh(i,j)*kappa))**2)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  do k=3,kbm1
     do j=1,jm
        do i=1,im
           gg(i,j,k)=1.d0/(a(i,j,k)+c(i,j,k)*(1.d0-ee(i,j,k-1)) &
                & -(dti2*dtef(i,j,k)+1.d0))
           ee(i,j,k)=a(i,j,k)*gg(i,j,k)
           gg(i,j,k)=(dti2*(-prod(i,j,k)*l(i,j,k)*e1) &
                & +c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
        end do
     end do
  end do

  do k=kbm1,2,-1
     do j=1,jm
        do i=1,im
           vf(i,j,k)=ee(i,j,k)*vf(i,j,k+1)+gg(i,j,k)
        end do
     end do
  end do

  ! the following is to counter the problem of the ratio of two small
  ! numbers (l = q2l/q2) or one number becoming negative. Two options are
  ! included below. In this application, the second option, l was less
  ! noisy when uf or vf is small

  !$omp parallel  
  !$omp do private(i,j,k)    
  do k=2,kbm1
     do j=1,jm
        do i=1,im
           !           if(uf(i,j,k).le.small.or.vf(i,j,k).le.small) then
           !             uf(i,j,k)=small
           !             vf(i,j,k)=0.1*dt(i,j)*small
           !           end if
           uf(i,j,k)=abs(uf(i,j,k))
           vf(i,j,k)=abs(vf(i,j,k))
        end do
     end do
  end do
  !$omp end do

  ! the following section solves for km and kh

  ! note that sm and sh limit to infinity when gh approaches 0.0288
  !$omp do private(i,j,k,coef1,coef2,coef3)
  do k=1,kb
     do j=1,jm
        do i=1,im
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

  ! there are 2 options for kq which, unlike km and kh, was not derived by
  ! Mellor and Yamada but was purely empirical based on neutral boundary
  ! layer data. The choice is whether or not it should be subject to the
  ! stability factor, sh. Generally, there is not a great difference in
  ! output
  !$omp do private(i,j,k)
  do k=1,kb
     do j=1,jm
        do i=1,im
           prod(i,j,k)=l(i,j,k)*sqrt(abs(q2(i,j,k)))
           kq(i,j,k)=(prod(i,j,k)*0.41d0*sh(i,j,k)+kq(i,j,k))*0.5d0
           !            kq(i,j,k)=(prod(i,j,k)*.20+kq(i,j,k))*.5e0
           km(i,j,k)=(prod(i,j,k)*sm(i,j,k)+km(i,j,k))*0.5d0
           kh(i,j,k)=(prod(i,j,k)*sh(i,j,k)+kh(i,j,k))*0.5d0
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  ! cosmetics: make boundr. values as interior (even if not used, printout
  ! may show strange values)
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

end subroutine profq
!_______________________________________________________________________
subroutine advt1(fb,f,fclim,ff)
  
  ! integrate conservative scalar equations
  ! this is centred scheme, as originally provide in POM (previously
  ! called advt)
  !     2019.04 S.Ohishi real --> double precision

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k

  real(kind = r_size),intent(inout) :: fb(im,jm,kb),f(im,jm,kb)
  real(kind = r_size),intent(in) :: fclim(im,jm,kb)
  real(kind = r_size),intent(out) :: ff(im,jm,kb)

  real(kind = r_size) xflux(im,jm,kb),yflux(im,jm,kb),zflux(im,jm,kb)
  real(kind = r_size) xyzflux(im,jm,kb)

  real(kind = r_size) dfmx(im,jm,kb),dtmx(im,jm),dxm(im,jm) !f(i-1)--- fmx(i)/u(i)--- f(i) 
  real(kind = r_size) dfmy(im,jm,kb),dtmy(im,jm),dym(im,jm) !f(j-1)--- fmy(j)/v(j)--- f(j) 
  real(kind = r_size) dfmz(im,jm,kb) !f(k-1)---fmz(k)/w(k)---f(k)
  real(kind = r_size) dtmt(im,jm) !h+efb --- dtmt --- h+etf

  !     Initialization
  ff(1:im,1:jm,1:kb)=0.d0
  f(1:im,1:jm,kb)=f(1:im,1:jm,kbm1)
  fb(1:im,1:jm,kb)=fb(1:im,1:jm,kbm1)

  if(budget == 1)then
     xadvterm(1:im,1:jm,1:kb)=0.d0
     yadvterm(1:im,1:jm,1:kb)=0.d0
     zadvterm(1:im,1:jm,1:kb)=0.d0
     advterm(1:im,1:jm,1:kb)=0.d0
     xdifterm(1:im,1:jm,1:kb)=0.d0
     ydifterm(1:im,1:jm,1:kb)=0.d0
  end if

  !     Intermidiate value
  !$omp parallel
  !$omp do private(i,j,k)
  do k=1,kbm1
     do j=2,jm
        do i=2,im
           dfmx(i,j,k)=0.5d0*(f(i-1,j,k)+f(i,j,k))
           dfmy(i,j,k)=0.5d0*(f(i,j-1,k)+f(i,j,k))
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)
  do k=2,kbm1
     do j=2,jmm1
        do i=2,imm1
           dfmz(i,j,k)=0.5d0*(f(i,j,k-1)+f(i,j,k))
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  dfmz(2:imm1,2:jmm1,1)=f(2:imm1,2:jmm1,1)
  dfmz(2:imm1,2:jmm1,kb)=0.d0

  !$omp parallel  
  !$omp do private(i,j)
  do j=2,jm
     do i=2,im
        dtmx(i,j)=0.5d0*(dt(i-1,j)+dt(i,j))
        dtmy(i,j)=0.5d0*(dt(i,j-1)+dt(i,j))
        dtmt(i,j)=h(i,j)+0.5d0*(etb(i,j)+etf(i,j))
        dxm(i,j)=0.5d0*(dx(i-1,j)+dx(i,j))
        dym(i,j)=0.5d0*(dy(i,j-1)+dy(i,j))
     end do
  end do
  !$omp end do

  !---  do horizontal fluxes
  !$omp do private(i,j,k)
  do k=1,kbm1
     do j=2,jm
        do i=2,im               
           xflux(i,j,k)= &
                & dtmx(i,j)*dfmx(i,j,k)*u(i,j,k)*dym(i,j) !UDf*dy
           yflux(i,j,k)= &
                & dtmy(i,j)*dfmy(i,j,k)*v(i,j,k)*dxm(i,j) !VDf*dx
        end do
     end do
  end do
  !$omp end do
  
  !     do vertical advection
  !$omp do private(i,j,k)
  do k=1,kb
     do j=2,jmm1
        do i=2,imm1
           zflux(i,j,k)=dfmz(i,j,k)*w(i,j,k)*art(i,j) !Wf*dxdy
        end do
     end do
  end do
  !$omp end do
  
  ! add advection and then step forward in time
  !$omp do private(i,j,k)
  do k=1,kbm1
     do j=2,jmm1
        do i=2,imm1
           xyzflux(i,j,k)= &
                & (xflux(i+1,j,k)-xflux(i,j,k))+(yflux(i,j+1,k)-yflux(i,j,k))+(zflux(i,j,k)-zflux(i,j,k+1))/dble(dz(i,j,k))
           ff(i,j,k)= &
                & (fb(i,j,k)*(h(i,j)+etb(i,j))*art(i,j) & !Dn-1*Tn-1*dxdy
                & -1.d0*dti2*xyzflux(i,j,k)) & !adv
                & /((h(i,j)+etf(i,j))*art(i,j)) !/(Dn+1*dxdy)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  !     exchanged3d_mpi
  call exchange3d_mpi(ff(:,:,1:kbm1),im,jm,kbm1)
      
  !$omp parallel
  if(budget == 1)then

     !$omp do private(i,j,k)
     do k=1,kbm1
        do j=2,jmm1
           do i=2,imm1
                  
              xadvterm(i,j,k)= &
                   & -1.d0*(xflux(i+1,j,k)-xflux(i,j,k)) & !-d(UDT)*dy
                   & +0.5d0*(dfmx(i+1,j,k)+dfmx(i,j,k)) &
                   & *(dtmx(i+1,j)*u(i+1,j,k)*dym(i+1,j)-dtmx(i,j)*u(i,j,k)*dym(i,j)) !Td(UD)*dy
              xadvterm(i,j,k)= &
                   & fsm(i,j)*xadvterm(i,j,k)/((h(i,j)+etf(i,j))*art(i,j)) !-UdT/dx = (-d(UDT)*dy+Td(UD)*dy)/(D*dxdy)
                  
              yadvterm(i,j,k)= &
                   & -1.d0*(yflux(i,j+1,k)-yflux(i,j,k)) & !-d(VDT)*dx
                   & +0.5d0*(dfmy(i,j+1,k)+dfmy(i,j,k)) &
                   & *(dtmy(i,j+1)*v(i,j+1,k)*dxm(i,j+1)-dtmy(i,j)*v(i,j,k)*dxm(i,j)) !Td(VD)*dx
              yadvterm(i,j,k)= &
                   & fsm(i,j)*yadvterm(i,j,k)/((h(i,j)+etf(i,j))*art(i,j)) !-VdT/dy = (-d(VDT)*dx+Td(VD)*dx)/(D*dxdy)
                  
              zadvterm(i,j,k)= &
                   & -1.d0*(zflux(i,j,k)-zflux(i,j,k+1)) & !-d(wT)*dxdy
                   & +0.5d0*(dfmz(i,j,k+1)+dfmz(i,j,k))*(w(i,j,k)-w(i,j,k+1))*art(i,j) !Tdw*dxdy
              zadvterm(i,j,k)=fsm(i,j)*zadvterm(i,j,k) &
                   & /((h(i,j)+etf(i,j))*art(i,j)*dz(i,j,k)) !-wdT/dz=(-d(wT)*dxdy+Tdw*dxdy)/(D*dxdy*dz)

              advterm(i,j,k)= &
                   & fsm(i,j)*(ff(i,j,k)-fb(i,j,k))/dti2

           end do
        end do
     end do
     !$omp end do

  end if
      
  !--- add horizontal diffusive fluxes
  !$omp do private(i,j,k)
  do k=1,kb
     do j=1,jm
        do i=1,im
           fb(i,j,k)=fb(i,j,k)-fclim(i,j,k)
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)
  do k=1,kbm1
     do j=2,jm
        do i=2,im
           xflux(i,j,k)= &
                & -0.5d0*(aam(i,j,k)+aam(i-1,j,k))*(h(i,j)+h(i-1,j))*tprni &
                & *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j) &
                & *(dy(i,j)+dy(i-1,j))*0.5d0/(dx(i,j)+dx(i-1,j))
           yflux(i,j,k)= &
                & -0.5d0*(aam(i,j,k)+aam(i,j-1,k))*(h(i,j)+h(i,j-1))*tprni &
                & *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j) &
                & *(dx(i,j)+dx(i,j-1))*0.5d0/(dy(i,j)+dy(i,j-1))
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)
  do k=1,kbm1
     do j=2,jmm1
        do i=2,imm1
           ff(i,j,k)=ff(i,j,k) &
                & -1.d0*dble(dti2)*(xflux(i+1,j,k)-xflux(i,j,k)+yflux(i,j+1,k)-yflux(i,j,k)) &
                & /((h(i,j)+etf(i,j))*art(i,j))
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)
  do k=1,kb
     do j=1,jm
        do i=1,im
           fb(i,j,k)=fb(i,j,k)+fclim(i,j,k)
        end do
     end do
  end do
  !$omp end do
      
  !     conserve horizontal diffusion
  if(budget == 1)then
     !$omp do private(i,j,k)
     do k=1,kbm1
        do j=2,jmm1
           do i=2,imm1
              xdifterm(i,j,k)= &
                   & -(xflux(i+1,j,k)-xflux(i,j,k))/((h(i,j)+etf(i,j))*art(i,j))
              ydifterm(i,j,k)= &
                   & -(yflux(i,j+1,k)-yflux(i,j,k))/((h(i,j)+etf(i,j))*art(i,j))
           end do
        end do
     end do
     !$omp end do
  end if
  !$omp end parallel


  f(1:im,1:jm,kb)=0.d0
  fb(1:im,1:jm,kb)=0.d0
  
  !     exchanged3d_mpi
  call exchange3d_mpi(ff(:,:,1:kbm1),im,jm,kbm1)
  
  if(budget == 1)then
     call exchange3d_mpi(xadvterm(:,:,1:kbm1),im,jm,kbm1)
     call exchange3d_mpi(yadvterm(:,:,1:kbm1),im,jm,kbm1)
     call exchange3d_mpi(zadvterm(:,:,1:kbm1),im,jm,kbm1)
     call exchange3d_mpi(advterm(:,:,1:kbm1),im,jm,kbm1)
     call exchange3d_mpi(xdifterm(:,:,1:kbm1),im,jm,kbm1)
     call exchange3d_mpi(ydifterm(:,:,1:kbm1),im,jm,kbm1)
  end if
    
end subroutine advt1
!_______________________________________________________________________
subroutine advt2(fb,f,fclim,ff)
  !     integrate conservative scalar equations
  !     this is a first-order upstream scheme, which reduces implicit
  !     diffusion using the Smolarkiewicz iterative upstream scheme with an
  !     antidiffusive velocity
  !     it is based on the subroutines of Gianmaria Sannino (Inter-university
  !     Computing Consortium, Rome, Italy) and Vincenzo Artale (Italian
  !     National Agency for New Technology and Environment, Rome, Italy)
  !     S.Ohishi modified 2019.04

  !$use omp_lib  
  use common_pom_var
  implicit none
      
  integer i,j,k,itera
      
  real(kind = r_size),intent(inout) :: fb(im,jm,kb)
  real(kind = r_size),intent(in) :: f(im,jm,kb),fclim(im,jm,kb)
  real(kind = r_size),intent(out) :: ff(im,jm,kb)
  
  real(kind = r_size) xflux(im,jm,kb),yflux(im,jm,kb)
  real(kind = r_size) zflux(im,jm,kb),xyzflux(im,jm,kb)
  real(kind = r_size) fbmem(im,jm,kb),eta(im,jm)
  real(kind = r_size) xmassflux(im,jm,kb),ymassflux(im,jm,kb)
  real(kind = r_size) zwflux(im,jm,kb)
   
  real(kind = r_size) dxm(im,jm),dym(im,jm)
  real(kind = r_size) dtmx(im,jm),dtmy(im,jm)
       
  !       Initialization
  fbmem(1:im,1:jm,1:kb)=fb(1:im,1:jm,1:kb)
  xmassflux(1:im,1:jm,1:kb)=0.d0
  ymassflux(1:im,1:jm,1:kb)=0.d0
  eta(1:im,1:jm)=etb(1:im,1:jm)
  zwflux(1:im,1:jm,1:kb)=w(1:im,1:jm,1:kb)
  
  !       Preparation
  fb(1:im,1:jm,kb)=fb(1:im,1:jm,kbm1)
     
  !$omp parallel  
  !$omp do private(i,j)  
  do j=2,jm
     do i=2,im
        dxm(i,j)=0.5d0*(dx(i,j-1)+dx(i,j))
        dym(i,j)=0.5d0*(dy(i-1,j)+dy(i,j))
        dtmx(i,j)=0.5d0*(dt(i-1,j)+dt(i,j))
        dtmy(i,j)=0.5d0*(dt(i,j-1)+dt(i,j))
     end do
  end do
  !$omp end do

  !     calculate horizontal mass fluxes
  !$omp do private(i,j,k)  
  do k=1,kbm1
     
     do j=2,jmm1
        do i=2,im
           xmassflux(i,j,k)=dym(i,j)*dtmx(i,j)*u(i,j,k)
        end do
     end do
     
     do j=2,jm
        do i=2,imm1
           ymassflux(i,j,k)=dxm(i,j)*dtmy(i,j)*v(i,j,k)
        end do
     end do
     
  end do
  !$omp end do
  !$omp end parallel
  
  !     Conservation
  if(budget == 1)then
     xadvterm(1:im,1:jm,1:kb)=0.d0
     yadvterm(1:im,1:jm,1:kb)=0.d0
     zadvterm(1:im,1:jm,1:kb)=0.d0
     advterm(1:im,1:jm,1:kb)=0.d0
     xdifterm(1:im,1:jm,1:kb)=0.d0
     ydifterm(1:im,1:jm,1:kb)=0.d0     
  end if
  
  !     start Smolarkiewicz scheme
  do itera=1,nitera
         
     !     upwind advection scheme
     !$omp parallel
     !$omp do private(i,j,k)      
     do k=1,kbm1
        do j=2,jm
           do i=2,im
              xflux(i,j,k)= &
                   & 0.5d0*((xmassflux(i,j,k)+abs(xmassflux(i,j,k)))*fbmem(i-1,j,k) &
                   & +(xmassflux(i,j,k)-abs(xmassflux(i,j,k)))*fbmem(i,j,k))
                  
              yflux(i,j,k)= &
                   & 0.5d0*((ymassflux(i,j,k)+abs(ymassflux(i,j,k)))*fbmem(i,j-1,k) &
                   & +(ymassflux(i,j,k)-abs(ymassflux(i,j,k)))*fbmem(i,j,k))
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel
     
     zflux(2:imm1,2:jmm1,1)=0.d0
     zflux(2:imm1,2:jmm1,kb)=0.d0

     if(itera == 1)then

        !$omp parallel
        !$omp do private(i,j)      
        do j=2,jmm1
           do i=2,imm1
              zflux(i,j,1)=w(i,j,1)*fb(i,j,1)*art(i,j)
           end do
        end do
        !$omp end do
        !$omp end parallel

     end if

     !$omp parallel
     !$omp do private(i,j,k)
     do k=2,kbm1
        do j=2,jmm1
           do i=2,imm1
              zflux(i,j,k)= &
                   & 0.5d0*((zwflux(i,j,k)+abs(zwflux(i,j,k)))*fbmem(i,j,k) &
                   & +(zwflux(i,j,k)-abs(zwflux(i,j,k)))*fbmem(i,j,k-1))
              zflux(i,j,k)=zflux(i,j,k)*art(i,j)
           end do
        end do
     end do
     !$omp end do
     
     !     add net advective fluxes and step forward in time
     !$omp do private(i,j,k)
     do k=1,kbm1
        do j=2,jmm1
           do i=2,imm1
                                    
              xyzflux(i,j,k)= &
                   & xflux(i+1,j,k)-xflux(i,j,k)+yflux(i,j+1,k)-yflux(i,j,k) &
                   & +(zflux(i,j,k)-zflux(i,j,k+1))/dz(i,j,k)
              ff(i,j,k)= &
                   & (fbmem(i,j,k)*(h(i,j)+eta(i,j))*art(i,j)-dti2*xyzflux(i,j,k)) &
                   & /((h(i,j)+etf(i,j))*art(i,j))
                                    
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel

     ! next line added on 22-Jul-2009 by Raffaele Bernardello
     call exchange3d_mpi(ff(:,:,1:kbm1),im,jm,kbm1)
     
     !     Conservation
     if(budget == 1 .and. itera == 1)then

        !$omp parallel
        !$omp do private(i,j,k)
        do k=1,kbm1
           do j=2,jmm1
              do i=2,imm1
                                          
                 xadvterm(i,j,k)= &
                      & -1.d0*(xflux(i+1,j,k)-xflux(i,j,k)) & !-d(UDT)*dy
                      & +fb(i,j,k)*(u(i+1,j,k)*dtmx(i+1,j)*dym(i+1,j)-u(i,j,k)*dtmx(i,j)*dym(i,j)) !Td(UD)*dy
                 xadvterm(i,j,k)=fsm(i,j)*xadvterm(i,j,k)/art(i,j) !-UdT/dx*D
                 xadvterm(i,j,k)=xadvterm(i,j,k)/(h(i,j)+etf(i,j)) !-UdT/dx
                     
                 yadvterm(i,j,k)= &
                      & -1.d0*(yflux(i,j+1,k)-yflux(i,j,k)) & !-d(VDT)*dx
                      & +fb(i,j,k)*(v(i,j+1,k)*dtmy(i,j+1)*dxm(i,j+1)-v(i,j,k)*dtmy(i,j)*dxm(i,j)) !Td(VD)*dx
                 yadvterm(i,j,k)=fsm(i,j)*yadvterm(i,j,k)/art(i,j) !-VdT/dy*D
                 yadvterm(i,j,k)=yadvterm(i,j,k)/(h(i,j)+etf(i,j)) 
                 !     -VdT/dy
                     
                 zadvterm(i,j,k)= &
                      & -1.d0*(zflux(i,j,k)-zflux(i,j,k+1)) & !-d(wT)*dxdy
                      & +fb(i,j,k)*(w(i,j,k)-w(i,j,k+1))*art(i,j) !Tdw*dxdy
                 zadvterm(i,j,k)= &
                      & fsm(i,j)*zadvterm(i,j,k)/(art(i,j)*dz(i,j,k)) !-wdT/dz
                 zadvterm(i,j,k)= &
                      & zadvterm(i,j,k)/(h(i,j)+etf(i,j)) !-w/D dT/dz
                     
                 advterm(i,j,k)= &
                      & fsm(i,j)*(ff(i,j,k)-fb(i,j,k))/dti2
                     
              end do
           end do
        end do
        !$omp end do
        !$omp end parallel
        
     else if(budget == 1)then

        !$omp parallel
        !$omp do private(i,j,k)        
        do k=1,kbm1
           do j=2,jmm1
              do i=2,imm1
                 
                 xadvterm(i,j,k)=xadvterm(i,j,k) &
                      & -1.d0*(xflux(i+1,j,k)-xflux(i,j,k)) &
                      & /((h(i,j)+etf(i,j))*art(i,j))
                     
                 yadvterm(i,j,k)=yadvterm(i,j,k) &
                      & -1.d0*(yflux(i,j+1,k)-yflux(i,j,k)) &
                      & /((h(i,j)+etf(i,j))*art(i,j))
                     
                 zadvterm(i,j,k)=zadvterm(i,j,k) &
                      & -1.d0*(zflux(i,j,k)-zflux(i,j,k+1)) &
                      & /((h(i,j)+etf(i,j))*dz(i,j,k)*art(i,j))
                     
                 advterm(i,j,k)=fsm(i,j)*(ff(i,j,k)-fb(i,j,k))/dti2
                 
              end do
           end do
        end do
        !$omp end do
        !$omp end parallel

     end if
         
     if(budget == 1)then
        call exchange3d_mpi(xadvterm(:,:,1:kbm1),im,jm,kbm1)
        call exchange3d_mpi(yadvterm(:,:,1:kbm1),im,jm,kbm1)
        call exchange3d_mpi(zadvterm(:,:,1:kbm1),im,jm,kbm1)
        call exchange3d_mpi(advterm(:,:,1:kbm1),im,jm,kbm1)
     end if
     
     !     calculate antidiffusion velocity
     call smol_adif(xmassflux,ymassflux,zwflux,ff)

     eta(1:im,1:jm)=etf(1:im,1:jm)
     fbmem(1:im,1:jm,1:kb)=ff(1:im,1:jm,1:kb)
          
     !     end of Smolarkiewicz scheme
  end do !itera
  
  ! add horizontal diffusive fluxes
  !$omp parallel
  !$omp do private(i,j,k)
  do k=1,kb
     do j=1,jm
        do i=1,im
           fb(i,j,k)=fb(i,j,k)-fclim(i,j,k)
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)  
  do k=1,kbm1
     do j=2,jm
        do i=2,im
           xflux(i,j,k)= &
                & -0.5d0*(aam(i,j,k)+aam(i-1,j,k))*(h(i,j)+h(i-1,j))*tprni &
                & *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)*(dy(i,j)+dy(i-1,j))*0.5d0 &
                & /(dx(i,j)+dx(i-1,j))
           yflux(i,j,k)= &
                & -0.5d0*(aam(i,j,k)+aam(i,j-1,k))*(h(i,j)+h(i,j-1))*tprni &
                & *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)*(dx(i,j)+dx(i,j-1))*0.5d0 &
                & /(dy(i,j)+dy(i,j-1))
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)  
  do k=1,kb
     do j=1,jm
        do i=1,im
           fb(i,j,k)=fb(i,j,k)+fclim(i,j,k)
        end do
     end do
  end do
  !$omp end do
  
  !     add net horizontal fluxes and step forward in time
  !$omp do private(i,j,k)
  do k=1,kbm1
     do j=2,jmm1
        do i=2,imm1
           ff(i,j,k)=ff(i,j,k) &
                & -dti2*(xflux(i+1,j,k)-xflux(i,j,k)+yflux(i,j+1,k)-yflux(i,j,k)) &
                & /((h(i,j)+etf(i,j))*art(i,j))
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  !     Conservation
  if(budget == 1)then

     !$omp parallel
     !$omp do private(i,j,k)
     do k=1,kbm1
        do j=2,jmm1
           do i=2,imm1
              xdifterm(i,j,k)= &
                   & -(xflux(i+1,j,k)-xflux(i,j,k))/((h(i,j)+etf(i,j))*art(i,j))
              ydifterm(i,j,k)=-(yflux(i,j+1,k)-yflux(i,j,k))/((h(i,j)+etf(i,j))*art(i,j))
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel

     call exchange3d_mpi(xdifterm(:,:,1:kbm1),im,jm,kbm1)
     call exchange3d_mpi(ydifterm(:,:,1:kbm1),im,jm,kbm1)
     
  end if

  fb(1:im,1:jm,kb)=0.d0
  call exchange3d_mpi(ff(:,:,:),im,jm,kb)
  
end subroutine advt2

!_______________________________________________________________________
subroutine smol_adif(xmassflux,ymassflux,zwflux,ff)
  ! calculate the antidiffusive velocity used to reduce the numerical
  ! diffusion associated with the upstream differencing scheme
  ! this is based on a subroutine of Gianmaria Sannino (Inter-university
  ! Computing Consortium, Rome, Italy) and Vincenzo Artale (Italian
  ! National Agency for New Technology and Environment, Rome, Italy)
  !     S. Ohishi 2019.03 real --> double precision

  !$use omp_lib  
  use common_pom_var
  implicit none

  real(kind = r_size),parameter :: value_min=1.d-9,epsilon=1.d-14

  integer i,j,k
  
  real(kind = r_size),intent(inout) :: ff(im,jm,kb)
  real(kind = r_size),intent(out) :: xmassflux(im,jm,kb),ymassflux(im,jm,kb)
  real(kind = r_size),intent(out) :: zwflux(im,jm,kb)
  real(kind = r_size) mol
  real(kind = r_size) udx,u2dt,vdy,v2dt,wdz,w2dt
  
  ! apply temperature and salinity mask
  !$omp parallel
  !$omp do private(i,j,k)
  do k=1,kb
     do i=1,im
        do j=1,jm
           ff(i,j,k)=ff(i,j,k)*fsm(i,j)
        end do
     end do
  end do
  !$omp end do
  
  !     recalculate mass fluxes with antidiffusion velocity
  !$omp do private(i,j,k,udx,u2dt,mol)
  do k=1,kbm1
     do j=2,jmm1
        do i=2,im
           if(ff(i,j,k) < value_min .or. ff(i-1,j,k) < value_min)then
              xmassflux(i,j,k)=0.d0
           else
              udx=abs(xmassflux(i,j,k))
              u2dt=dti2*xmassflux(i,j,k)*xmassflux(i,j,k)*2.d0 &
                   & /(aru(i,j)*(dt(i-1,j)+dt(i,j)))
              if(abs(udx) < abs(u2dt))then
                 xmassflux(i,j,k)=0.d0
              else
                 mol=(ff(i,j,k)-ff(i-1,j,k))/(ff(i-1,j,k)+ff(i,j,k)+epsilon)
                 xmassflux(i,j,k)=(udx-u2dt)*mol*sw
              end if
           end if
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k,vdy,v2dt,mol)  
  do k=1,kbm1
     do j=2,jm
        do i=2,imm1
           if(ff(i,j,k) < value_min .or. ff(i,j-1,k) < value_min) then
              ymassflux(i,j,k)=0.d0
           else
              vdy=abs(ymassflux(i,j,k))
              v2dt=dti2*ymassflux(i,j,k)*ymassflux(i,j,k)*2.d0 &
                   & /(arv(i,j)*(dt(i,j-1)+dt(i,j)))
              if(abs(vdy) < abs(v2dt))then
                 ymassflux(i,j,k)=0.d0
              else
                 mol=(ff(i,j,k)-ff(i,j-1,k))/(ff(i,j-1,k)+ff(i,j,k)+epsilon)
                 ymassflux(i,j,k)=(vdy-v2dt)*mol*sw
              end if
           end if
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k,wdz,w2dt,mol)  
  do k=2,kbm1
     do j=2,jmm1
        do i=2,imm1
           if(ff(i,j,k) < value_min .or. ff(i,j,k-1) < value_min) then
              zwflux(i,j,k)=0.d0
           else
              wdz=abs(zwflux(i,j,k))
              w2dt=dti2*zwflux(i,j,k)*zwflux(i,j,k)/(dzz(i,j,k-1)*dt(i,j))
              if(abs(wdz) < abs(w2dt))then
                 zwflux(i,j,k)=0.d0
              else
                 mol=(ff(i,j,k-1)-ff(i,j,k))/(ff(i,j,k)+ff(i,j,k-1)+epsilon)
                 zwflux(i,j,k)=(wdz-w2dt)*mol*sw
              end if
           end if
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine smol_adif

!_______________________________________________________________________
subroutine proft(f,wfsurf,fsurf,nbc)
  ! solves for vertical diffusion of temperature and salinity using method
  ! described by Richmeyer and Morton (1967)
  ! note: wfsurf and swrad are negative values when water column is
  ! warming or salt is being added

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer,intent(in) :: nbc  
  real(kind = r_size),intent(in) :: wfsurf(im,jm),fsurf(im,jm)

  real(kind = r_size),intent(inout) :: f(im,jm,kb)
  
  integer i,j,k

  real(kind = r_size) dh(im,jm)
  real(kind = r_size) a(im,jm,kb),c(im,jm,kb)
  real(kind = r_size) ee(im,jm,kb),gg(im,jm,kb)
  real(kind = r_size) rad(im,jm,kb)

  real(kind = r_size) fb(im,jm,kb) !for conservation (budget)
  
  ! irradiance parameters after Paulson and Simpson (1977)
  !   ntp                                  1        2       3       4       5
  !   Jerlov type                          i        ia      ib      ii      iii
  real(kind = r_size),parameter :: r(5) = (/0.58d0, 0.62d0, 0.67d0, 0.77d0, 0.78d0 /)
  real(kind = r_size),parameter :: ad1(5)=(/0.35d0, 0.60d0, 1.00d0, 1.50d0, 1.40d0 /)
  real(kind = r_size),parameter :: ad2(5)=(/23.0d0, 20.0d0, 17.0d0, 14.0d0, 7.90d0 /)

  ! surface boundary condition:
  !       nbc   prescribed    prescribed   short wave
  !             temperature      flux      penetration
  !             or salinity               (temperature
  !                                           only)
  !        1        no           yes           no
  !        2        no           yes           yes
  !        3        yes          no            no
  !        4        yes          no            yes
  ! note that only 1 and 3 are allowed for salinity.

  
  !     Conservation
  if(budget == 1)then
     fb(1:im,1:jm,1:kb)=f(1:im,1:jm,1:kb)
     radterm(1:im,1:jm,1:kb)=0.d0
  end if
  
  ! the following section solves the equation
  !     dti2*(kh*f')'-f=-fb
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
     do j=1,jm
        do i=1,im
           a(i,j,k-1)=-dti2*(kh(i,j,k)+umol)/(dz(i,j,k-1)*dzz(i,j,k-1)*dh(i,j)*dh(i,j))
           c(i,j,k)=-dti2*(kh(i,j,k)+umol)/(dz(i,j,k)*dzz(i,j,k-1)*dh(i,j)*dh(i,j))
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  ! calculate penetrative radiation. At the bottom any unattenuated
  ! radiation is deposited in the bottom layer
  rad(1:im,1:jm,1:kb)=0.d0

  if(nbc == 2 .or. nbc == 4)then
     !        do k=1,kbm1
     ! rad(kb) passes into the bottom without
     ! heating lower layer or is reflected from bottom without been adsorbed.
     ! 2013.10.05 Miyazawa
     !$omp parallel
     !$omp do private(i,j,k)
     do k=1,kb
        do j=1,jm
           do i=1,im
              rad(i,j,k)=swrad(i,j) &
                   & *(r(ntp)*exp(z(i,j,k)*dh(i,j)/ad1(ntp)) &
                   & +(1.d0-r(ntp))*exp(z(i,j,k)*dh(i,j)/ad2(ntp)))
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel
  endif

  if(nbc == 1)then

     !$omp parallel     
     !$omp do private(i,j)
     do j=1,jm
        do i=1,im
           ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.d0)
           gg(i,j,1)=-dti2*wfsurf(i,j)/(-dz(i,j,1)*dh(i,j))-f(i,j,1)
           gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.d0)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  else if(nbc == 2) then

     !$omp parallel
     !$omp do private(i,j)     
     do j=1,jm
        do i=1,im
           ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.d0)
           gg(i,j,1)= &
                & dti2*(wfsurf(i,j)+rad(i,j,1)-rad(i,j,2))/(dz(i,j,1)*dh(i,j))-f(i,j,1)
           gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.d0)
        end do
     end do
     !$omp end do
     !$omp end parallel

  else if(nbc == 3 .or. nbc == 4)then

     ee(1:im,1:jm,1)=0.d0
     gg(1:im,1:jm,1)=fsurf(1:im,1:jm)
     
  endif

  do k=2,kbm2
     do j=1,jm
        do i=1,im
           gg(i,j,k)=1.d0/(a(i,j,k)+c(i,j,k)*(1.d0-ee(i,j,k-1))-1.d0)
           ee(i,j,k)=a(i,j,k)*gg(i,j,k)
           gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-f(i,j,k)+dti2*(rad(i,j,k)-rad(i,j,k+1)) &
                & /(dh(i,j)*dz(i,j,k)))*gg(i,j,k)
        end do
     end do
  end do
  
  ! bottom adiabatic boundary condition
  !$omp parallel
  !$omp do private(i,j)  
  do j=1,jm
     do i=1,im
        f(i,j,kbm1)= &
             & (c(i,j,kbm1)*gg(i,j,kbm2)-f(i,j,kbm1)+dti2*(rad(i,j,kbm1)-rad(i,j,kb)) &
             & /(dh(i,j)*dz(i,j,kbm1))) &
             & /(c(i,j,kbm1)*(1.d0-ee(i,j,kbm2))-1.d0)
     end do
  end do
  !$omp end do
  !$omp end parallel

  do k=kbm2,1,-1
     do j=1,jm
        do i=1,im
           f(i,j,k)=(ee(i,j,k)*f(i,j,k+1)+gg(i,j,k))
        end do
     end do
  end do
    
  !     Conservation
  !$omp parallel
  if(budget == 1)then

     !$omp do private(i,j)
     do j=1,jm
        do i=1,im
           sfcterm(i,j)=-wfsurf(i,j)/(dz(i,j,1)*dh(i,j))
        end do
     end do
     !$omp end do

     !$omp do private(i,j,k)     
     do k=1,kbm1
        do j=1,jm
           do i=1,im
              radterm(i,j,k)=-(rad(i,j,k)-rad(i,j,k+1))/(dh(i,j)*dz(i,j,k))
              zdifterm(i,j,k)=(f(i,j,k)-fb(i,j,k))/dti2-radterm(i,j,k)
           end do
        end do
     end do
     !$omp end do

     !$omp do private(i,j)          
     do j=1,jm
        do i=1,im
           zdifterm(i,j,1)=zdifterm(i,j,1)-sfcterm(i,j)
        end do
     end do
     !$omp end do
     
  end if
  !$omp end parallel
  
end subroutine proft
!_______________________________________________________________________
subroutine advu
  ! do horizontal and vertical advection of u-momentum, and includes
  ! coriolis, surface slope and baroclinic terms

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k

  ! do vertical advection

  uf(1:im,1:jm,1:kb)=0.d0

  !$omp parallel
  !$omp do private(i,j,k)  
  do k=2,kbm1
     do j=1,jm
        do i=2,im
           uf(i,j,k)=0.25d0*(w(i,j,k)+w(i-1,j,k))*(u(i,j,k)+u(i,j,k-1))
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  ! combine horizontal and vertical advection with coriolis, surface
  ! slope and baroclinic terms
  do k=1,kbm1
     do j=2,jmm1
        do i=2,imm1
           uf(i,j,k)=advx(i,j,k)+(uf(i,j,k)-uf(i,j,k+1))*aru(i,j)/dz(i,j,k) &
                & -aru(i,j)*0.25d0*(cor(i,j)*dt(i,j)*(v(i,j+1,k)+v(i,j,k)) &
                & +cor(i-1,j)*dt(i-1,j)*(v(i-1,j+1,k)+v(i-1,j,k))) &
                & +grav*0.125d0*(dt(i,j)+dt(i-1,j)) &
                & *(egf(i,j)-egf(i-1,j)+egb(i,j)-egb(i-1,j)+(e_atmos(i,j)-e_atmos(i-1,j))*2.d0) &
                & *(dy(i,j)+dy(i-1,j)) &
                & +drhox(i,j,k)
        end do
     end do
  end do
  
  !  step forward in time
  !$omp parallel
  !$omp do private(i,j,k)
  do k=1,kbm1
     do j=2,jmm1
        do i=2,imm1
           uf(i,j,k)=((h(i,j)+etb(i,j)+h(i-1,j)+etb(i-1,j))*aru(i,j)*ub(i,j,k) &
                & -2.d0*dti2*uf(i,j,k))/((h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))*aru(i,j))
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine advu
!_______________________________________________________________________
subroutine advv
  ! do horizontal and vertical advection of v-momentum, and includes
  ! coriolis, surface slope and baroclinic terms

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k

  ! do vertical advection
  vf(1:im,1:jm,1:kb)=0.d0

  !$omp parallel
  !$omp do private(i,j,k)  
  do k=2,kbm1
     do j=2,jm
        do i=1,im
           vf(i,j,k)=0.25d0*(w(i,j,k)+w(i,j-1,k))*(v(i,j,k)+v(i,j,k-1))
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  ! combine horizontal and vertical advection with coriolis, surface
  ! slope and baroclinic terms
  do k=1,kbm1
     do j=2,jmm1
        do i=2,imm1
           vf(i,j,k)=advy(i,j,k) &
                & +(vf(i,j,k)-vf(i,j,k+1))*arv(i,j)/dz(i,j,k) &
                & +arv(i,j)*0.25d0*(cor(i,j)*dt(i,j)*(u(i+1,j,k)+u(i,j,k))+cor(i,j-1)*dt(i,j-1)*(u(i+1,j-1,k)+u(i,j-1,k))) &
                & +grav*0.125d0*(dt(i,j)+dt(i,j-1)) &
                & *(egf(i,j)-egf(i,j-1)+egb(i,j)-egb(i,j-1)+(e_atmos(i,j)-e_atmos(i,j-1))*2.d0) &
                & *(dx(i,j)+dx(i,j-1)) &
                & +drhoy(i,j,k)
        end do
     end do
  end do
  
  ! step forward in time
  !$omp parallel
  !$omp do private(i,j,k)
  do k=1,kbm1
     do j=2,jmm1
        do i=2,imm1
           vf(i,j,k)=((h(i,j)+etb(i,j)+h(i,j-1)+etb(i,j-1))*arv(i,j)*vb(i,j,k) &
                & -2.d0*dti2*vf(i,j,k))/((h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1)) &
                & *arv(i,j))
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine advv
!_______________________________________________________________________
subroutine profu
  ! solves for vertical diffusion of x-momentum using method described by
  ! Richmeyer and Morton (1967)
  ! note: wusurf has the opposite sign to the wind speed

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k

  real(kind = r_size) a(im,jm,kb),c(im,jm,kb)
  real(kind = r_size) ee(im,jm,kb),gg(im,jm,kb)
  real(kind = r_size) dh(im,jm)

  ! the following section solves the equation
  !   dti2*(km*u')'-u=-ub

  dh(1:im,1:jm)=1.d0
  
  !$omp parallel
  !$omp do private(i,j)  
  do j=2,jm
     do i=2,im
        dh(i,j)=(h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))*0.5d0
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)  
  do k=1,kb
     do j=2,jm
        do i=2,im
           c(i,j,k)=(km(i,j,k)+km(i-1,j,k))*0.5d0
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  do k=2,kbm1
     do j=1,jm
        do i=1,im
           a(i,j,k-1)=-dti2*(c(i,j,k)+umol)/(dz(i,j,k-1)*dzz(i,j,k-1)*dh(i,j)*dh(i,j))
           c(i,j,k)=-dti2*(c(i,j,k)+umol)/(dz(i,j,k)*dzz(i,j,k-1)*dh(i,j)*dh(i,j))
        end do
     end do
  end do

  !$omp parallel  
  !$omp do private(i,j)    
  do j=1,jm
     do i=1,im
        ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.d0)
        gg(i,j,1)=(-dti2*wusurf(i,j)/(-dz(i,j,1)*dh(i,j))-uf(i,j,1)) &
             & /(a(i,j,1)-1.d0)
     end do
  end do
  !$omp end do
  !$omp end parallel

  do k=2,kbm2
     do j=1,jm
        do i=1,im
           gg(i,j,k)=1.d0/(a(i,j,k)+c(i,j,k)*(1.d0-ee(i,j,k-1))-1.d0)
           ee(i,j,k)=a(i,j,k)*gg(i,j,k)
           gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-uf(i,j,k))*gg(i,j,k)
        end do
     end do
  end do

  !$omp parallel  
  !$omp do private(i,j)     
  do j=2,jmm1
     do i=2,imm1
        tps(i,j)=0.5d0*(cbc(i,j)+cbc(i-1,j)) &
             & *sqrt(ub(i,j,kbm1)**2+(0.25d0*(vb(i,j,kbm1)+vb(i,j+1,kbm1) &
             & +vb(i-1,j,kbm1)+vb(i-1,j+1,kbm1)))**2)
        uf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-uf(i,j,kbm1)) &
             & /(tps(i,j)*dti2/(-dz(i,j,kbm1)*dh(i,j))-1.d0 &
             & -(ee(i,j,kbm2)-1.d0)*c(i,j,kbm1))
        uf(i,j,kbm1)=uf(i,j,kbm1)*dum(i,j)
     end do
  end do
  !$omp end do
  !$omp end parallel

  do k=kbm2,1,-1
     do j=2,jmm1
        do i=2,imm1
           uf(i,j,k)=(ee(i,j,k)*uf(i,j,k+1)+gg(i,j,k))*dum(i,j)
        end do
     end do
  end do

  !$omp parallel  
  !$omp do private(i,j)     
  do j=2,jmm1
     do i=2,imm1
        wubot(i,j)=-tps(i,j)*uf(i,j,kbm1)
     end do
  end do
  !$omp end do
  !$omp end parallel

  call exchange2d_mpi(wubot,im,jm)

end subroutine profu
!_______________________________________________________________________
subroutine profv
  ! solves for vertical diffusion of x-momentum using method described by
  ! Richmeyer and Morton (1967)
  ! note: wvsurf has the opposite sign to the wind speed

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k

  real(kind = r_size) a(im,jm,kb),c(im,jm,kb)
  real(kind = r_size) ee(im,jm,kb),gg(im,jm,kb)
  real(kind = r_size) dh(im,jm)

  ! the following section solves the equation
  !     dti2*(km*u')'-u=-ub

  dh(1:im,1:jm)=1.d0
  
  !$omp parallel
  !$omp do private(i,j)   
  do j=2,jm
     do i=2,im
        dh(i,j)=0.5d0*(h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1))
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)   
  do k=1,kb
     do j=2,jm
        do i=2,im
           c(i,j,k)=(km(i,j,k)+km(i,j-1,k))*0.5d0
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  do k=2,kbm1
     do j=1,jm
        do i=1,im
           a(i,j,k-1)=-dti2*(c(i,j,k)+umol)/(dz(i,j,k-1)*dzz(i,j,k-1)*dh(i,j)*dh(i,j))
           c(i,j,k)=-dti2*(c(i,j,k)+umol)/(dz(i,j,k)*dzz(i,j,k-1)*dh(i,j)*dh(i,j))
        end do
     end do
  end do

  !$omp parallel  
  !$omp do private(i,j)   
  do j=1,jm
     do i=1,im
        ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.d0)
        gg(i,j,1)=(-dti2*wvsurf(i,j)/(-dz(i,j,1)*dh(i,j))-vf(i,j,1))/(a(i,j,1)-1.d0)
     end do
  end do
  !$omp end do
  !$omp end parallel

  do k=2,kbm2
     do j=1,jm
        do i=1,im
           gg(i,j,k)=1.d0/(a(i,j,k)+c(i,j,k)*(1.d0-ee(i,j,k-1))-1.d0)
           ee(i,j,k)=a(i,j,k)*gg(i,j,k)
           gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
        end do
     end do
  end do

  !$omp parallel    
  !$omp do private(i,j)   
  do j=2,jmm1
     do i=2,imm1
        tps(i,j)=0.5d0*(cbc(i,j)+cbc(i,j-1)) &
             & *sqrt((0.25d0*(ub(i,j,kbm1)+ub(i+1,j,kbm1)+ub(i,j-1,kbm1)+ub(i+1,j-1,kbm1)))**2 &
             & +vb(i,j,kbm1)**2)
        vf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-vf(i,j,kbm1)) &
             & /(tps(i,j)*dti2 &
             & /(-dz(i,j,kbm1)*dh(i,j))-1.d0-(ee(i,j,kbm2)-1.e0)*c(i,j,kbm1))
        vf(i,j,kbm1)=vf(i,j,kbm1)*dvm(i,j)
     end do
  end do
  !$omp end do
  !$omp end parallel

  do k=kbm2,1,-1
     do j=2,jmm1
        do i=2,imm1
           vf(i,j,k)=(ee(i,j,k)*vf(i,j,k+1)+gg(i,j,k))*dvm(i,j)
        end do
     end do
  end do

  !$omp parallel  
  !$omp do private(i,j)   
  do j=2,jmm1
     do i=2,imm1
        wvbot(i,j)=-tps(i,j)*vf(i,j,kbm1)
     end do
  end do
  !$omp end do
  !$omp end parallel

  call exchange2d_mpi(wvbot,im,jm)

end subroutine profv
!_______________________________________________________________________
subroutine realvertvl
  ! calculates real vertical velocity (wr)

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k
  real(kind = r_size) dxr,dxl,dyt,dyb

  wr(1:im,1:jm,1:kb)=0.d0

  !$omp parallel
  !$omp do private(i,j,k,dxr,dxl,dyt,dyb)  
  do k=1,kbm1

     do j=1,jm
        do i=1,im
           tps(i,j)=zz(i,j,k)*dt(i,j)+et(i,j)
        end do
     end do
     
     do j=2,jmm1
        do i=2,imm1
           dxr=2.d0/(dx(i+1,j)+dx(i,j))
           dxl=2.d0/(dx(i,j)+dx(i-1,j))
           dyt=2.d0/(dy(i,j+1)+dy(i,j))
           dyb=2.d0/(dy(i,j)+dy(i,j-1))
           wr(i,j,k)=0.5d0*(w(i,j,k)+w(i,j,k+1)) &
                & +0.5d0*(u(i+1,j,k)*(tps(i+1,j)-tps(i,j))*dxr &
                & +u(i,j,k)*(tps(i,j)-tps(i-1,j))*dxl &
                & +v(i,j+1,k)*(tps(i,j+1)-tps(i,j))*dyt &
                & +v(i,j,k)*(tps(i,j)-tps(i,j-1))*dyb) &
                & +(1.d0+zz(i,j,k))*(etf(i,j)-etb(i,j))/dti2
        end do
     end do

  end do
  !$omp end do
  !$omp end parallel
  
  call exchange3d_mpi(wr(:,:,1:kbm1),im,jm,kbm1)

  
  if(n_south == -1)then
     wr(1:im,1,1:kb)=wr(1:im,2,1:kb)
  end if

  if(n_north == -1)then
     wr(1:im,jm,1:kb)=wr(1:im,jmm1,1:kb)
  end if

  if(n_west == -1)then
     wr(1,1:jm,1:kb)=wr(2,1:jm,1:kb)
  end if
     
  if(n_east == -1)then
     wr(im,1:jm,1:kb)=wr(imm1,1:jm,1:kb)   
  end if

  !$omp parallel  
  !$omp do private(i,j,k)
  do k=1,kbm1
     do j=1,jm
        do i=1,im
           wr(i,j,k)=fsm(i,j)*wr(i,j,k)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine realvertvl

