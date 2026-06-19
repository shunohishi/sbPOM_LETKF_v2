module mod_density

contains

  !-----------------------------------------------------------------------
  !     Density (UNESCO, 1981) |
  !-----------------------------------------------------------------------
  !
  !     Reference: A. E. Gill (1982), Atmosphere-Ocean Dynamics, 
  !     Appendix Three, Properties of seawater.
  !
  !-----------------------------------------------------------------------
  !
  ! temp: degree C
  ! sal:-
  ! depth: meter = dbar
  ! pres: bar
  !-----------------------------------------------------------------------

  elemental subroutine estimate_density(temp,sal,depth,rho)

    implicit none

    !---Common
    real(kind = 8) pres ![bar]
    real(kind = 8) k,k0,kw
    real(kind = 8) rho0,rhow
    real(kind = 8) t2,t3,t4,t5
    real(kind = 8) s15

    !---IN
    real(kind = 8),intent(in) :: temp  ![degree C]
    real(kind = 8),intent(in) :: sal   ![-]
    real(kind = 8),intent(in) :: depth ![m] = [dbar]

    !---OUT
    real(kind = 8),intent(out) :: rho

    pres=0.1d0*depth

    t2=temp*temp
    t3=t2*temp
    t4=t3*temp
    t5=t4*temp

    s15=sal*sqrt(sal)

    rhow=999.842594d0 &
         & +6.793952d-2*temp -9.095290d-3*t2 &
         & +1.001685d-4*t3 -1.120083d-6*t4 &
         & +6.536332d-9*t5

    rho0=rhow &
         & +sal*(0.824493d0 -4.0899d-3*temp +7.6438d-5*t2 &
         & -8.2467d-7*t3 +5.3875d-9*t4) &
         & +s15 *(-5.72466d-3 +1.0227d-4*temp &
         & -1.6546d-6*t2) +4.8314d-4*sal*sal 

    kw=19652.21d0 &
         & +148.4206d0*temp -2.327105d0*t2 &
         & +1.360477d-2*t3 -5.155288d-5*t4

    k0=kw &
         & +sal*(54.6746d0 -0.603459d0*temp &
         & +1.09987d-2*t2 -6.1670d-5*t3) &
         & +s15*(7.944d-2 +1.6483d-2*temp -5.3009d-4*t2)

    k=k0 &
         & +pres*(3.239908d0 +1.43713d-3*temp &
         & +1.16092d-4*t2 -5.77905d-7*t3) &
         & +pres*sal*(2.2838d-3 -1.0981d-5*temp -1.6078d-6*t2) &
         & +1.91075d-4*pres*s15 &
         & +pres*pres*(8.50935d-5 &
         & -6.12293d-6*temp +5.2787d-8*t2) &
         & +pres*pres*sal *(-9.9348d-7 &
         & +2.0816d-8*temp +9.1697d-10*t2)

    rho=rho0*k/(k-pres)

  end subroutine estimate_density

  !---------------------------------------------------------------
  ! Potential Temperature |
  !---------------------------------------------------------------

  elemental subroutine estimate_pt(temp,sal,depth,pt)

    implicit none

    !---Common
    real(kind = 8) pres
    real(kind = 8) p2,p3
    real(kind = 8) t2,t3
    real(kind = 8) s35

    !---IN
    real(kind = 8),intent(in) :: temp,sal,depth

    !---OUT
    real(kind = 8),intent(out) :: pt

    pres=0.1d0*depth
    p2=pres*pres
    p3=p2*pres

    t2=temp*temp
    t3=t2*temp

    s35=sal-35.d0

    pt=temp &
         & -pres*(3.6504d-4 +8.3198d-5*temp -5.4065d-7*t2 +4.0274d-9*t3) &
         & -pres*s35*(1.7439d-5 -2.9778d-7*temp) &
         & -p2*(8.9309d-7 -3.1628d-8*temp +2.1987d-10*t2) &
         & +4.1057d-9*s35*p2 &
         & -p3*(-1.6056d-10 +5.0484d-12*temp)

  end subroutine estimate_pt

  !----------------------------------------------------------------
  ! Specific Heat |
  !----------------------------------------------------------------

  elemental subroutine estimate_cp(temp,sal,cp)

    implicit none

    !---Common
    real(kind = 8) t2,t3,t4
    real(kind = 8) s15

    !---IN
    real(kind = 8),intent(in) :: temp,sal

    !---OUT
    real(kind = 8),intent(out) :: cp

    t2=temp*temp
    t3=t2*temp
    t4=t3*temp

    s15=sal*sqrt(sal)

    cp=4217.4d0 &
         & -3.720283d0*temp &
         & +0.1412855d0*t2 &
         & -2.654387d-3*t3 &
         & +2.093236d-5*t4 &
         & +sal*(-7.6444d0 +0.107276d0*temp -1.3839d-3*t2) &
         & +s15*(0.17709d0 -4.0772d-3*temp +5.3539d-5*t2)

  end subroutine estimate_cp

  !----------------------------------------------------------------
  ! Thermal expansion coefficient |
  !----------------------------------------------------------------

  elemental subroutine estimate_alpha(temp,sal,depth,alpha)

    implicit none

    !---Common
    real(kind=8) pres
    real(kind=8) kw,k0,k
    real(kind=8) rho,rho0,rhow
    real(kind=8) dkwdt,dk0dt,dkdt
    real(kind=8) drho0dt,drhowdt

    real(kind=8) p2
    real(kind=8) t2,t3,t4,t5
    real(kind=8) s15,s2

    !---IN
    real(kind=8),intent(in) :: temp,sal,depth

    !---OUT
    real(kind=8),intent(out) :: alpha

    pres=0.1d0*depth
    p2=pres*pres

    t2=temp*temp
    t3=t2*temp
    t4=t3*temp
    t5=t4*temp

    s15=sal*sqrt(sal)
    s2=sal*sal


    !=========================
    ! Density components
    !=========================

    rhow=999.842594d0 &
         +6.793952d-2*temp -9.095290d-3*t2 +1.001685d-4*t3 &
         -1.120083d-6*t4 +6.536332d-9*t5

    rho0=rhow &
         +sal*(0.824493d0 -4.0899d-3*temp +7.6438d-5*t2 &
         -8.2467d-7*t3 +5.3875d-9*t4) &
         +s15*(-5.72466d-3 +1.0227d-4*temp -1.6546d-6*t2) +4.8314d-4*s2

    !=========================
    ! Compressibility
    !=========================

    kw=19652.21d0 &
         +148.4206d0*temp -2.327105d0*t2 +1.360477d-2*t3 -5.155288d-5*t4

    k0=kw &
         +sal*(54.6746d0 -0.603459d0*temp +1.09987d-2*t2 -6.1670d-5*t3) &
         +s15*(7.944d-2 +1.6483d-2*temp -5.3009d-4*t2)

    k=k0 &
         +pres*(3.239908d0 +1.43713d-3*temp +1.16092d-4*t2 -5.77905d-7*t3) &
         +pres*sal*(2.2838d-3 -1.0981d-5*temp -1.6078d-6*t2) &
         +1.91075d-4*pres*s15 &
         +p2*(8.50935d-5 -6.12293d-6*temp +5.2787d-8*t2) &
         +p2*sal*(-9.9348d-7 +2.0816d-8*temp +9.1697d-10*t2)

    rho=rho0*k/(k-pres)

    !=========================
    ! Derivative of density
    !=========================

    drhowdt=6.793952d-2 &
         -2.d0*9.095290d-3*temp +3.d0*1.001685d-4*t2 &
         -4.d0*1.120083d-6*t3 +5.d0*6.536332d-9*t4

    drho0dt=drhowdt &
         +sal*(-4.0899d-3 +2.d0*7.6438d-5*temp &
         -3.d0*8.2467d-7*t2 +4.d0*5.3875d-9*t3) &
         +s15*(1.0227d-4 -2.d0*1.6546d-6*temp)

    dkwdt=148.4206d0 &
         -2.d0*2.327105d0*temp +3.d0*1.360477d-2*t2 -4.d0*5.155288d-5*t3

    dk0dt=dkwdt &
         +sal*(-0.603459d0 +2.d0*1.09987d-2*temp -3.d0*6.1670d-5*t2) &
         +s15*(1.6483d-2 -2.d0*5.3009d-4*temp)

    dkdt=dk0dt &
         +pres*(1.43713d-3 +2.d0*1.16092d-4*temp -3.d0*5.77905d-7*t2) &
         +pres*sal*(-1.0981d-5 -2.d0*1.6078d-6*temp) &
         +p2*(-6.12293d-6 +2.d0*5.2787d-8*temp) &
         +p2*sal*(2.0816d-8 +2.d0*9.1697d-10*temp)

    alpha=((k*drho0dt+rho0*dkdt)*(k-pres)-rho0*k*dkdt)/(k-pres)/(k-pres)

    alpha=-1.d0*alpha/rho


  end subroutine estimate_alpha

  !-------------------------------------------------------------
  ! Saline contraction coefficient |
  !-------------------------------------------------------------

  elemental subroutine estimate_beta(temp,sal,depth,beta)

    implicit none

    !---Common
    real(kind=8) pres
    real(kind=8) kw,k0,k
    real(kind=8) rho,rho0,rhow
    real(kind=8) dk0ds,dkds
    real(kind=8) drho0ds

    real(kind=8) p2
    real(kind=8) t2,t3,t4,t5
    real(kind=8) s15,sqrt_s,s2
    real(kind=8) invkp 

    !---IN
    real(kind=8),intent(in) :: temp,sal,depth

    !---OUT
    real(kind=8),intent(out) :: beta


    pres=0.1d0*depth
    p2=pres*pres

    t2=temp*temp
    t3=t2*temp
    t4=t3*temp
    t5=t4*temp

    sqrt_s=sqrt(sal)
    s15=sal*sqrt_s
    s2=sal*sal

    !=========================
    ! Density components
    !=========================

    rhow=999.842594d0 &
         +6.793952d-2*temp -9.095290d-3*t2 &
         +1.001685d-4*t3 -1.120083d-6*t4 &
         +6.536332d-9*t5

    rho0=rhow &
         +sal*(0.824493d0 -4.0899d-3*temp +7.6438d-5*t2 &
         -8.2467d-7*t3 +5.3875d-9*t4) &
         +s15*(-5.72466d-3 +1.0227d-4*temp -1.6546d-6*t2) &
         +4.8314d-4*s2

    !=========================
    ! Compressibility
    !=========================

    kw=19652.21d0 &
         +148.4206d0*temp -2.327105d0*t2 &
         +1.360477d-2*t3 -5.155288d-5*t4

    k0=kw &
         +sal*(54.6746d0 -0.603459d0*temp +1.09987d-2*t2 -6.1670d-5*t3) &
         +s15*(7.944d-2 +1.6483d-2*temp -5.3009d-4*t2)

    k=k0 &
         +pres*(3.239908d0 +1.43713d-3*temp +1.16092d-4*t2 -5.77905d-7*t3) &
         +pres*sal*(2.2838d-3 -1.0981d-5*temp -1.6078d-6*t2) &
         +1.91075d-4*pres*s15 &
         +p2*(8.50935d-5 -6.12293d-6*temp +5.2787d-8*t2) &
         +p2*sal*(-9.9348d-7 +2.0816d-8*temp +9.1697d-10*t2)

    rho=rho0*k/(k-pres)

    !=========================
    ! Derivative wrt salinity
    !=========================

    drho0ds=0.824493d0 &
         -4.0899d-3*temp +7.6438d-5*t2 -8.2467d-7*t3 +5.3875d-9*t4 &
         +1.5d0*sqrt_s*(-5.72466d-3 +1.0227d-4*temp -1.6546d-6*t2) &
         +2.d0*4.8314d-4*sal

    dk0ds=54.6746d0 &
         -0.603459d0*temp +1.09987d-2*t2 -6.1670d-5*t3 &
         +1.5d0*sqrt_s*(7.944d-2 +1.6483d-2*temp -5.3009d-4*t2)

    dkds=dk0ds &
         +pres*(2.2838d-3 -1.0981d-5*temp -1.6078d-6*t2) &
         +1.5d0*1.91075d-4*pres*sqrt_s &
         +p2*(-9.9348d-7 +2.0816d-8*temp +9.1697d-10*t2)

    !=========================
    ! beta = 1/rho * drho/dS
    !=========================

    invkp=1.d0/(k-pres)

    beta=((k*drho0ds+rho0*dkds)*(k-pres)-rho0*k*dkds)*invkp*invkp    
    beta=beta/rho

  end subroutine estimate_beta

end module mod_density

