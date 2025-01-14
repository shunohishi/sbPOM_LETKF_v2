!
! Module for definition of constants and estimation of parameters
! for Mellor-Yamada-Nakanishi-Nino-Furuichi versions of level 2.5 turbulence
! clousure models. Coded by vsm, November 2012
!

! Define GEN_SIGMA variable for the generalized sigma coordinates model,
! else regular sigma coordinates model is considered.
! Different is the vertical model layers layout description used.

!#define GEN_SIGMA

module MYNNF_lev25_2012

  !$use omp_lib
  use common_pom_var, only: r_size
  implicit none

  private
  !
  !     Fixed parameters for MYNNF-2.5 scheme 
  !
  real(kind = r_size),parameter :: g1 = 0.235d0
  real(kind = r_size),parameter :: b1 = 24.d0
  real(kind = r_size),parameter :: b2 = 15.d0
  real(kind = r_size),parameter :: c2 = 0.75d0
  real(kind = r_size),parameter :: c3 = 0.352d0
  real(kind = r_size),parameter :: c4 = 0.d0
  real(kind = r_size),parameter :: c5 = 0.2d0
  real(kind = r_size),parameter :: pr = 0.74d0
  real(kind = r_size),parameter :: vk = 0.4d0   !Karman const
  real(kind = r_size),parameter :: a1 = b1*( 1.d0-3.d0*g1 )/6.d0
  ! c1 = g1 -1.0/( 3.0*a1*b1**(1.0/3.0) ),  !on SX only integer power is allowed
  real(kind = r_size),parameter :: c1 = 0.13706763d0
  real(kind = r_size),parameter :: a2 = a1*(g1-c1)/(g1*pr)
  real(kind = r_size),parameter :: g2 = b2/b1*(1.d0-c3)+2.d0*a1/b1*(3.d0-2.d0*c2)

  real(kind = r_size),parameter :: alp1 = 0.23d0      !Multiplier for L as planetary boundary layer z-scale
  real(kind = r_size),parameter :: alp2 = 0.53d0      !0.53 in Furuichi et al, 2012 for LES; 1 in NN
  real(kind = r_size),parameter :: alp3 = 1.d0/3.7d0  !1/3.7 in NN; 1 in Furuichi et al, 2012, but no 
  real(kind = r_size),parameter :: almost_zero = 1d-12 !small value for nonsingular numerics
  !
  !     Constants for surfase boundary layer Monin-Obukhov length scale "Lmo" estimation
  !
  real(kind = r_size),parameter :: grav = 9.808d0
  real(kind = r_size),parameter :: alpha_sw = 1.7d-4  !1/K, sea water thermal expansion coefficient
  real(kind = r_size),parameter :: beta_sw = 7.5d-4   !1/PSU, sea water thermal expansion coefficient
  !
  !     Traditional values for Mellor-Yamada level 2.5 model parameters:
  !
  real(kind = r_size),parameter :: a1_my=0.92d0
  real(kind = r_size),parameter :: a2_my= 0.74d0
  real(kind = r_size),parameter :: b1_my=16.6d0
  real(kind = r_size),parameter :: b2_my=10.1d0
  real(kind = r_size),parameter :: c1_my= 0.08d0

  public mynn_get_l,mynn_get_q2l2,mynn_get_ShSm,my_get_ShSm

contains

  subroutine mynn_get_l(l,q2,boygr,wtsurf,wssurf,wusurf_t,wvsurf_t,z0,z,dzz,dh,im,jm,kb,ntp)

    !     Turbulence length scale estimation for MYNN model

    use common_pom_var, only: r_size
    implicit none
    integer,intent(in):: im,jm,kb,ntp
    real(kind = r_size),intent(out):: l(im,jm,kb)     !turbulence length scale, m 
    real(kind = r_size),intent(in)::  q2(im,jm,kb)    !square of turbulence velocity scale
    real(kind = r_size),intent(in)::  boygr(im,jm,kb) !in POM, N**2 ~= -boygr/1.025
    real(kind = r_size),intent(in)::  wtsurf(im,jm)   !Total_Heat_Flux/(pho*cp), K*(m/s), positive if ocean losses heat
    real(kind = r_size),intent(in)::  wssurf(im,jm)   !Total_Salt_Flux/Rho_fresh_water, [PSU*(m/s)], -(s+sbias)*WQ*RoFWR where
    ! WQ [kg/m**2/s] is positive mass flux for evaporation (salination),
    ! wssurf is positive for precipitation case (desalination)
    real(kind = r_size),intent(in):: wusurf_t(im,jm),wvsurf_t(im,jm) !surface_stress/Rho, defined in T-points (not in U and V points as wusrf and wvsurf of POM).
    real(kind = r_size),intent(in):: z0(im,jm)               !surface roughness length; bottom roughness is fixed as 0.01 m
    real(kind = r_size),intent(in):: z(im,jm,kb),dzz(im,jm,kb),dh(im,jm) !sea depth dh is required only for SIGMA layers model

    !     Working arrays used for PBL scale,
    !     inverse Monin-Obukhov scale and buoyancy flux at the surface

    real(kind = r_size) lt(im,jm),LmoR(im,jm),Bf(im,jm)

    integer i,j,k
    real(kind = r_size) qdz,zk,lb,ls,lh,lr,u_star,zn_MO,N,q,qc,rr,hf

    !     Irradiance parameters after Paulson and Simpson, JPO, 1977, 952-956.
    !     Same as in the PROFT subroutine. Theoretically for the surface buoyancy flux
    !     estimation have to use only part of SWR adsorbed by ocean surface mixed layer.
    !     Use ocean depth as mixed layer depth. Part of radiation reaching the bottom
    !     is considered to be removed from system. Could be important to prevent over-heating
    !     (proft) and over-stabilization near the surface in the shallow coastal zones.
    !
    !     Irradiance parameters after Paulson and Simpson (1977)
    !                 NTP         =     1      2      3      4     5
    !             JERLOV TYPE     =     I      IA     IB     II    III
    real(kind = r_size),parameter:: R(5)   = (/ 0.58d0, 0.62d0, 0.67d0, 0.77d0, 0.78d0/)
    real(kind = r_size),parameter:: AD1(5) = (/ 0.35d0, 0.60d0, 1.00d0, 1.50d0, 1.40d0/)
    real(kind = r_size),parameter:: AD2(5) = (/ 23.0d0, 20.0d0, 17.0d0, 14.0d0, 7.90d0/)

    !     Use lt and bf as working arrays for vertical integrals
    lt(1:im,1:jm) = 0.d0
    bf(1:im,1:jm) = 0.d0

    !$omp parallel
    !$omp do private(i,j,k,qdz,zk)
    do k=2,kb-1
       do j=1,jm
          do i=1,im
             qdz=sqrt(abs(q2(i,j,k)))*dzz(i,j,k-1)*dh(i,j)
             zk=abs(z(i,j,k)*dh(i,j))
             lt(i,j)=lt(i,j)+qdz*zk
             bf(i,j)=bf(i,j)+qdz
          end do
       end do
    end do
    !$omp end do

    !$omp do private(i,j)    
    do j=1,jm
       do i=1,im
          ! ** Length scale depending on the PBL depth
          lt(i,j)=alp1*max(lt(i,j)/max(bf(i,j),almost_zero),almost_zero)
          l(i,j,1)=vk*abs(z0(i,j)) !surface value, L -> ZERO: define it in main cycle
          l(i,j,kb)=vk*0.01d0      !bottom value, L -> ZERO (almost_zero)
       end do
    end do
    !$omp end do

    !
    !     Estimate surface buoyancy flux, momentum flux parameter u_star and
    !     inverse of Monin-Obukhov length scale LmoR.
    !

    !$omp do private(i,j,zk,rr,hf,u_star)    
    do j=1,jm
       do i=1,im
          !         negative sea depth
          zk=z(i,j,kb)*dh(i,j)

          !         Remove radiation that reaches bottom. It could decrease
          !         surface heating stabilization impact in shallow waters

          rr = exp(zk/ad1(ntp))*r(ntp)+exp(zk/ad2(ntp))*(1.d0-r(ntp))
          hf = wtsurf(i,j)
          Bf(i,j) = grav*(alpha_sw*hf-beta_sw*wssurf(i,j))  !Bf positive for the case of convection
          !          u_star = sqrt(0.5*sqrt(
          !     *        (wusurf(i,j)+wusurf(i+1,j))**2+
          !     *        (wvsurf(i,j)+wvsurf(i,j+1))**2))
          u_star = sqrt(sqrt(wusurf_t(i,j)**2+wvsurf_t(i,j)**2))
          LmoR(i,j)= -vk*Bf(i,j)/max(u_star**3,almost_zero)
       end do
    end do
    !$omp end do

    !
    !     Turbulence length scale estimation
    !
    !$omp do private(i,j,k,zk,zn_MO,ls,lh,N,q,lb,qc,lr)
    do k=1,kb-1
       do j=1,jm
          do i=1,im
             !
             !           Length scale in the surface ML, possibly - convective ML
             !
             zk = vk*(abs(z(i,j,k)*dh(i,j))+z0(i,j))
             zn_MO = zk/vk*LmoR(i,j)

             if(abs(z(i,j,k)*dh(i,j)) < lt(i,j))then

                if(zn_MO >= 1.d0)then  !stabilizing surface fluxes
                   ls = alp3*zk
                elseif(zn_MO >= 0.d0)then  !weekly stable to neutral stratification
                   ls = zk/(1.d0+2.7d0*zn_MO)
                else                      !convective instability could develope
                   ls = zk*(1.d0-100.d0*zn_MO)**0.2d0
                endif

             else
                ls = zk
             endif
             ls = max(ls,almost_zero)
             ! test case: no MO impact
             !            ls = max(zk,almost_zero)
             !
             !           Bottom mixed layer, distance from bottom. Assume zn_MO_bottom = 0
             !
             lh = vk*(abs((z(i,j,kb)-z(i,j,k))*dh(i,j))+0.01d0)
             !            lh = max(lh,almost_zero)
             !
             !           Length scale limited by the buoyancy effect
             !
             if(k > 1 .and. boygr(i,j,k) < 0.d0) then !statically stable stratification
                N = sqrt(-boygr(i,j,k)/1.025d0)
                q = sqrt(abs(q2(i,j,k)))
                lb=alp2*q/N
                ! test case: comment out MO impact terms
                if(abs(dh(i,j)*z(i,j,k)) < lt(i,j))then
                   !v20130313, consider surface conditions only in PBL
                   if( zn_MO < 0.d0 )then !count for convection impact
                      !                 Convective velocity scale
                      qc = (Bf(i,j)*Lt(i,j))**(1.d0/3.d0)
                      lb = lb*sqrt(1.d0+40.d0*qc/(Lt(i,j)*N))
                   endif
                endif

                lb=max(lb,almost_zero)
             else !statically unstable vertical stratification
                lb=1.d0/almost_zero
             endif

             !           Length scale controlled by the smallest length scale
             !           among the three length scales: lt, lb, and (ls-lh)

             lr = 1.d0/lt(i,j) + 1.d0/lb + 1.d0/min(ls,lh)
             l(i,j,k) = max(1.d0/lr,almost_zero)
          end do
       end do
       !        write(*,*)"MYNN a1,a2,c1,g2=",a1,a2,c1,g2
       !        write(*,*)'k,z,L,ls,lh,lb,lt,zn_MO,1/l_MO=',
       !     *             k,zk/vk,l(1,1,k),ls,lh,lb,lt(1,1),zn_MO,LmoR(1,1)
    end do
    !$omp end do
    !$omp end parallel

    
  end subroutine mynn_get_l
  
  !******************************************************************************
  ! Example of shear square (shear2) estimation in calling code
  !      do k=2,kb-1
  !        do j=1,jmm1
  !          do i=1,imm1
  !            shear2=0.25*((u(i,  j,k)-u(i,  j,k-1)
  !     $                   +u(i+1,j,k)-u(i+1,j,k-1))**2
  !     $                  +(v(i,  j,k)-v(i,j,  k-1)
  !     $                   +v(i,j+1,k)-v(i,j+1,k-1))**2)/
  !#ifdef GEN_SIGMA
  !     $                   dzz(i,j,k-1)**2
  !#else
  !     $                   (dzz(i,j,k-1)*dh(i,j))**2
  !#endif
  !          enddo
  !        enddo
  !      enddo
  !******************************************************************************
  
  subroutine mynn_get_q2l2(l,q2l2,boygr,shear2,im,jm,kb)

    !     Level 2 turbulent q square (q2l2) estimation for MYNN model
    !     Values are defined on the walls of the grid boxes.
    !     Boundaries (surface,bottom,lateral) must be treated in the calling code

    !$use omp_lib    
    use common_pom_var, only: r_size
    implicit none

    integer,intent(in):: im,jm,kb
    real(kind = r_size),intent(in) :: l(im,jm,kb), shear2(im,jm,kb)
    real(kind = r_size),intent(in) :: boygr(im,jm,kb)        !in POM N**2 ~= -boygr/1.025
    real(kind = r_size),intent(out)::   q2l2(im,jm,kb)


    !     some constants

    !     Critical Richardson number rfc:
    real(kind = r_size),parameter :: rfc = g1/(g1+g2)
    real(kind = r_size),parameter :: f1 = b1*(g1-c1)+3.d0*a2*(1.d0-c2)*(1.d0-c5)+2.d0*a1*(3.d0-2.d0*c2)
    real(kind = r_size),parameter :: f2 = b1*(g1+g2)-3.d0*a1*(1.d0-c2)
    real(kind = r_size),parameter :: rf1 = b1*(g1-c1)/f1
    real(kind = r_size),parameter :: rf2 = b1*g1/f2
    real(kind = r_size),parameter :: smc = a1/a2*f1/f2
    real(kind = r_size),parameter :: shc = 3.d0*a2*(g1+g2)

    real(kind = r_size),parameter :: ri1=0.5d0/smc
    real(kind = r_size),parameter:: ri2=rf1*smc
    real(kind = r_size),parameter:: ri3=4.d0*rf2*smc-2.d0*ri2
    real(kind = r_size),parameter:: ri4=ri2**2

    integer i,j,k
    real(kind = r_size) ri,rf,sh2,sm2

    !$omp parallel
    !$omp do private(i,j,k,Ri,Rf,sh2,sm2)            
    do k=2,kb-1
       do j=1,jm
          do i=1,im

             !Gradient Richardson number

             Ri = -boygr(i,j,k)/(1.025d0*max(shear2(i,j,k),almost_zero))

             !Flux Richardson number rf, here
             ! rfc=0.2984, rf1=0.374, rf2=0.313
             ! ri1=0.806 ri2=0.232 ri3=0.313 ri4=0.0537
             ! D < 0, sqrt is OK

             Rf = min( ri1*( Ri+ri2-sqrt(Ri**2-ri3*Ri+ri4) ), rfc )
             sh2 = shc*( rfc-Rf )/( 1.d0-Rf )
             sm2 = smc*( rf1-Rf )/( rf2-Rf ) * sh2
             q2l2(i,j,k) = b1*sm2*(1.d0-Rf) &
                  & *max(shear2(i,j,k),almost_zero)*l(i,j,k)**2
             !            write(*,*)"Ri,shc,smc,rf,rfc,rf1,rf2,sh2,sm2,l(i,j,k),"//
             !     *                "ri1,ri2,ri3,ri4="
             !            write(*,*) Ri,shc,smc,rf,rfc,rf1,rf2,sh2,sm2,l(i,j,k),
             !     *                 ri1,ri2,ri3,ri4
          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

    q2l2(1:im,1:jm,1) = 0.d0
    q2l2(1:im,1:jm,kb) = 0.d0

  end subroutine mynn_get_q2l2

  !******************************************************************************
  
  subroutine mynn_get_ShSm(sh,sm,l,q2,q2l2,boygr,shear2,im,jm,kb)

    !     Definition of level 2.5 Nakanishi-Nino stability functions.

    !$use omp_lib    
    use common_pom_var, only: r_size
    implicit none

    integer,intent(in):: im,jm,kb
    real(kind = r_size),intent(out)::   sh(im,jm,kb),sm(im,jm,kb)

    real(kind = r_size),intent(in) :: l(im,jm,kb), q2(im,jm,kb), q2l2(im,jm,kb)
    !      real,dimension(im,jm,kb),intent(in) :: shear2
    !start tmp  diagnostic test
    
    real(kind = r_size),intent(inout):: shear2(im,jm,kb)
    
    !end tmp  diagnostic test

    real(kind = r_size),intent(in):: boygr(im,jm,kb)        !in POM N**2 ~= -boygr/1.025
    integer i,j,k
    real(kind = r_size) ac,gh,gm,rr,ac2,a2c2,f1,f2,f3,f4,f5,d25

    !     vsm: limiting values for non-singular numerics
    !
    !      real,parameter:: gh_max = 4.55e-2
    !      real,parameter:: gm_max = 1e0     !1.0e-1

    !      ac = 1e0
    !      ac2 = ac*ac

    a2c2 = a2*(1.d0-c2)

    !$omp parallel
    !$omp do private(i,j,k,ac,rr,gh,gm,ac2,f1,f2,f3,f4,f5,d25)
    do k=1,kb
       do j=1,jm
          do i=1,im
             if(q2l2(i,j,k) > 0.d0 .and. abs(q2(i,j,k)) < q2l2(i,j,k))then
                ac = sqrt(abs(q2(i,j,k))/q2l2(i,j,k))
             else
                ac = 1.d0
             endif
             rr = (l(i,j,k)**2)/max(q2(i,j,k),almost_zero)
             gh = boygr(i,j,k)/1.025d0*rr
             gm = shear2(i,j,k)*rr
             !vsm, start: introduce limitations on gh and gm to avoid singularities
             !     instead of ac that is often very small in new developing unstable layers...
             !            gh = min(gh,gh_max)  !d25 could have two roots, ~0.047 and 0.554 for gm==0
             !            gm = min(gm,gm_max)
             !end
             ac2 = ac*ac
             f1 = 1.d0-3.d0*ac2*a2*b2*(1.d0-c3)*gh
             f2 = 1.d0-9.d0*ac2*a1*a2c2*gh
             f3 = f1+9.d0*ac2*a2*a2c2*(1.d0-c5)*gh
             f4 = f1-12.d0*ac2*a1*a2c2*gh
             f5 = 6.d0*ac2*a1*a1*gm
             d25 = f2*f4+f5*f3
             !start tmp diagnostic tests
            if(abs(d25) <= 0.d0)then
               shear2(i,j,k) = -1.d0+d25
            endif
            !end tmp  diagnostic test
            d25 = max(d25,almost_zero)
            sm(i,j,k)=ac*a1*(f3-3.d0*c1*f4)/d25
            sh(i,j,k)=ac*a2*(f2+3.d0*c1*f5)/d25
            !      if(i==485 .and. j==357 .and. k==2)then
            !      write(*,*)"MYNN: l,q2,q2l2,bg,sh2,sh,sm=",
            !     *l(i,j,k),q2(i,j,k),q2l2(i,j,k),boygr(i,j,k),shear2(i,j,k),
            !     *sh(i,j,k),sm(i,j,k)
            !      endif
            !      write(*,*)"MYNN gh=",gh,gm
         end do
      end do
   end do
   !$omp end do
   !$omp end parallel
   
 end subroutine mynn_get_ShSm

 !******************************************************************************
 
 subroutine my_get_ShSm(sh,sm,l,q2,boygr,im,jm,kb)

   !  Definition of level 2.5 Mellor-Yamada stability functions.

   !  NOTE: Richardson # dep. dissipation correction (Mellor, 2001; Ezer, 2000) 
   !  disabled here; stf=1.0 as initialized above.
   !  It is unclear yet if diss. corr. is needed when surf. waves are included.

   !$use omp_lib   
   use common_pom_var, only: r_size   
   implicit none
   
   integer,intent(in) :: im,jm,kb
   real(kind = r_size),intent(in) :: l(im,jm,kb),q2(im,jm,kb),boygr(im,jm,kb)  !in POM: N**2 ~= -boygr/1.025

   real(kind = r_size),intent(out) :: sh(im,jm,kb),sm(im,jm,kb)
   
   integer i,j,k
   real(kind = r_size) gh

   real(kind = r_size),parameter :: stf=1.d0 !stf(i,j,k), see above for dissipation correction

   real(kind = r_size),parameter :: coef1=a2_my*(1.d0-6.d0*a1_my/b1_my*stf)
   real(kind = r_size),parameter :: coef2=3.d0*a2_my*(b2_my/stf+6.d0*a1_my)
   real(kind = r_size),parameter :: coef3=a1_my*(1.d0-3.d0*c1_my-6.d0*a1_my/b1_my*stf)
   real(kind = r_size),parameter :: coef4=9.d0*a1_my*(2.d0*a1_my+a2_my)
   real(kind = r_size),parameter :: coef5=9.d0*a1_my*a2_my

   !$omp parallel
   !$omp do private(i,j,k,gh)           
   do k=1,kb
      do j=1,jm
         do i=1,im
            
            if(k == 1 .or. k == kb)then
               gh = 0.d0
            else
               gh=(l(i,j,k)**2)*boygr(i,j,k)/(1.025d0*q2(i,j,k))
               !              write(*,*)"MY gh=",gh,1./coef2,1./coef5

               !     sm and sh limit to infinity when gh approaches 0.0288, :

               gh=min(gh,0.028d0)
            endif

            !     vsm: for large L it is possible that sm becomes negative:
            !          it was found for week stable stratification and L~100m

            sh(i,j,k)=coef1/(1.d0-coef2*gh)
            sm(i,j,k)=(coef3+sh(i,j,k)*coef4*gh)/(1.d0-coef5*gh)
         end do
      end do
   end do
   !$omp end do
   !$omp end parallel
   
 end subroutine my_get_ShSm
 
end module MYNNF_lev25_2012
