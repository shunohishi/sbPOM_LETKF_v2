!_______________________________________________________________________
subroutine surface_airseaflux
  ! set time dependent surface airsea fluxes
  !     upward(sea->air) heat transport: positive wtsurf and swrad
  !     upward(sea->air) salt transport: positive wssurf

  !$use omp_lib
  use mod_thermo, only: q_sat
  use mod_aerobulk, only: AEROBULK_MODEL
  use common_pom_var
  implicit none

  integer,parameter :: n_iter=4 !number of interation 
  
  integer timing
  integer i,j
  integer ierr

  real(kind = r_dble) time0_data,ratio
  real(kind = r_size) qh(im,jm),qe(im,jm)
  real(kind = r_size) sst(im,jm)
  real(kind = r_size) windur(im,jm),windvr(im,jm)

  if(iint == 1) call read_atm_netcdf

  time0_data=dble(atmtime_julday)-dble(julday_start)
  call readtiming(timing,ratio,atmtime_dayint,time0_data)
  if(iint /= 1 .and. timing == 1) then
     call read_atm_netcdf
     if(my_task == master_task) &
          & write(6,'(a,i10,2f12.5)') 'atm (timing/ratio/time): ',timing,ratio,dti*float(iint-1)/86400.d0+time0
  end if
  
  !$omp parallel
  !$omp do private(i,j)
  do j=1,jm
     do i=1,im

        windu(i,j)=(1.d0-ratio)*windu0(i,j)+ratio*windu1(i,j)
        windv(i,j)=(1.d0-ratio)*windv0(i,j)+ratio*windv1(i,j)
        windur(i,j)=windu(i,j)
        windvr(i,j)=windv(i,j)
        winds(i,j)=sqrt(windu(i,j)*windu(i,j)+windv(i,j)*windv(i,j))

        airt(i,j)=(1.d0-ratio)*airt0(i,j)+ratio*airt1(i,j)
        airt(i,j)=max(airt(i,j), -273.15d0)
        
        airh(i,j)=(1.d0-ratio)*airh0(i,j)+ratio*airh1(i,j)
        airh(i,j)=max(airh(i,j), 0.d0)

        lwrad(i,j)=(1.d0-ratio)*lwrad0(i,j)+ratio*lwrad1(i,j)
        lwrad(i,j)=max(lwrad(i,j), 0.d0)

        swrad(i,j)=(1.d0-ratio)*swrad0(i,j)+ratio*swrad1(i,j)
        swrad(i,j)=max(swrad(i,j), 0.d0)
        
        slp(i,j)=(1.d0-ratio)*slp0(i,j)+ratio*slp1(i,j)
        slp(i,j)=max(slp(i,j), 0.d0)
        
        prep(i,j)=(1.d0-ratio)*prep0(i,j)+ratio*prep1(i,j)
        pflux(i,j)=fsm(i,j)*prep(i,j)*s(i,j,1)*(-1.d-3)/86400.d0
        
     end do
  end do
  !$omp end do

  ! unit change: [degree Celcius] -> [Kelvin]
  !$omp do private(i,j)
  do j=1,jm
     do i=1,im
        sst(i,j)=t(i,j,1)+tbias+273.15d0
        airt(i,j)=airt(i,j)+273.15d0
     end do
  end do
  !$omp end do

  !     Relative velocity
  if(lrtvf)then

     !$omp do private(i,j)
     do j=1,jm
        do i=1,im 

           if(dum(i,j) == 1.d0 .and. dum(i+1,j) == 1.d0)then
              windur(i,j)=windu(i,j)-0.5d0*(u(i,j,1)+u(i+1,j,1))
           elseif(dum(i,j) == 0.d0 .and. dum(i+1,j) == 1.d0)then
              windur(i,j)=windu(i,j)-u(i+1,j,1)
           elseif(dum(i,j) == 1.d0 .and. dum(i+1,j) == 0.d0)then
              windur(i,j)=windu(i,j)-u(i,j,1)
           end if

           if(dvm(i,j) == 1.d0 .and. dvm(i,j+1) == 1.d0)then
              windvr(i,j)=windv(i,j)-0.5d0*(v(i,j,1)+v(i,j+1,1))
           elseif(dvm(i,j) == 0.d0 .and. dvm(i,j+1) == 1.d0)then
              windvr(i,j)=windv(i,j)-v(i,j+1,1)
           elseif(dvm(i,j) == 1.d0 .and. dvm(i,j+1) == 0.d0)then
              windvr(i,j)=windv(i,j)-v(i,j,1)
           end if

        end do
     end do
     !$omp end do

  end if
  !$omp end parallel
  
  !--- surface heat flux & wind stress

  if(ithf_ws == 1 .and. lrtvf)then

     call latenth(qe,sst,airt,airh,qs,windur,windvr,im,jm)
     call sensibleh(qh,sst,airt,windur,windvr,im,jm)
     call wind_stress_mb(wusurf,wvsurf,windur,windvr,im,jm)

  elseif(ithf_ws == 1 .and. .not. lrtvf)then

     call latenth(qe,sst,airt,airh,qs,windu,windv,im,jm)
     call sensibleh(qh,sst,airt,windu,windv,im,jm)
     call wind_stress_mb(wusurf,wvsurf,windu,windv,im,jm)

  elseif(ithf_ws == 2)then
     call AEROBULK_MODEL &
          & ("coare",zt,zu,sst,airt,airh,windur,windvr,slp, &
          & qe,qh,wusurf,wvsurf,n_iter)
  elseif(ithf_ws == 3)then
     call AEROBULK_MODEL &
          & ("coare35",zt,zu,sst,airt,airh,windur,windvr,slp, &
          & qe,qh,wusurf,wvsurf,n_iter)
  elseif(ithf_ws == 4)then
     call AEROBULK_MODEL &
          & ("ncar",zt,zu,sst,airt,airh,windur,windvr,slp, &
          & qe,qh,wusurf,wvsurf,n_iter)
  elseif(ithf_ws == 5)then
     call AEROBULK_MODEL &
          & ("ecmwf",zt,zu,sst,airt,airh,windur,windvr,slp, &
          & qe,qh,wusurf,wvsurf,n_iter)
  else
     write(*,"(/a)") "POM terminated with error: lthf"
     call finalize_mpi
     stop
  end if
  
  !     make sign consistency of qh,qe
  !     ithf_ws=1: plus
  !     ithf_ws=2-5: minus
  if(ithf_ws == 1)then
     !$omp parallel
     !$omp do private(i,j)
     do j=1,jm
        do i=1,im
           qe(i,j)=qe(i,j)*fsm(i,j)
           qh(i,j)=qh(i,j)*fsm(i,j)
        end do
     end do
     !$omp end do
     !$omp end parallel
  else
     qs(:,:)=0.98d0*q_sat(sst,slp) ! qs[kg/kg] using q_sat in mod_thermo.f90
     !     qh,qe: minus --> plus
     !$omp parallel
     !$omp do private(i,j)
     do j=1,jm
        do i=1,im
           qe(i,j)=-1.d0*qe(i,j)*fsm(i,j)
           qh(i,j)=-1.d0*qh(i,j)*fsm(i,j)
           qs(i,j)=qs(i,j)*fsm(i,j)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end if

  !     wind stress: unit and direction changes
  !     unit: [N/m^2] = [kg/(m s^2)] --> [m^2/s^2]
  !$omp parallel
  !$omp do private(i,j)
  do j=1,jm
     do i=1,im
        wusurf(i,j)=-(wusurf(i,j))/rhoref*dum(i,j)
        wvsurf(i,j)=-(wvsurf(i,j))/rhoref*dvm(i,j)
     end do
  end do
  !$omp end do
  !$omp end parallel

  !---     Short & Longwave radiation
  call shortrad(im,jm,pi,north_e,swrad)
  call longrad2(lwrad,sst,im,jm)
  
  !---     Substitute variables in pom.h
  !$omp parallel
  !$omp do private(i,j)
  do j=1,jm
     do i=1,im
        !     lhf, shf: plus --> minus
        lhf(i,j)=-1.d0*qe(i,j)*fsm(i,j)
        shf(i,j)=-1.d0*qh(i,j)*fsm(i,j)
        lwr(i,j)=-1.d0*lwrad(i,j)*fsm(i,j)
        swr(i,j)=swrad(i,j)*fsm(i,j)
        qa(i,j) = 1.d3*airh(i,j)*fsm(i,j) ![g/kg]
        qs(i,j) = 1.d3*qs(i,j)*fsm(i,j)   ![g/kg]
        ta(i,j) = (airt(i,j)-273.15d0)*fsm(i,j) ![degree C]
        !     unit & direction changes: [m^2/s^2] --> [N/m^2]
        tauu(i,j)=-1.d0*wusurf(i,j)*rhoref*dum(i,j) ![kg/(m s^2) = N/m^2]
        tauv(i,j)=-1.d0*wvsurf(i,j)*rhoref*dvm(i,j) ![kg/(m s^2) = N/m^2]
        taus(i,j)=dum(i,j)*dvm(i,j) &
             & *sqrt(tauu(i,j)*tauu(i,j)+tauv(i,j)*tauv(i,j))
        !     lhf/(2.5*10^9) [m/s] = lhf*10^3*24*60*60/(2.5*10^9) [mm/day]
        evap(i,j)=fsm(i,j)*(-1.d0)*lhf(i,j)*86400.d0*0.4d0*1.d-6
        eflux(i,j)=fsm(i,j)*evap(i,j)*s(i,j,1)*(+1.d-3)/86400.d0
     end do
  end do
  !$omp end do
  
  ! unit change [W/m**2] -> [K m/sec]
  !$omp do private(i,j)
  do j=1,jm
     do i=1,im
        wtsurf(i,j)= &
             & (qh(i,j)+qe(i,j)+lwrad(i,j))/(4.1876d6)*fsm(i,j)
        swrad(i,j)=-swrad(i,j)/(4.1876d6)*fsm(i,j)
     end do
  end do
  !$omp end do
  !$omp end parallel

  if(issf == 2)then
     evap(:,:)=0.d0
     prep(:,:)=0.d0
     eflux(:,:)=0.d0
     pflux(:,:)=0.d0
  end if
  
end subroutine surface_airseaflux

!_______________________________________________________________________
subroutine surface_riverflux

  ! set time dependent surface fresh water fluxes
  !     upward(sea->air) heat transport: positive wtsurf and swrad
  !     upward(sea->air) salt transport: positive wssurf
  !     S.Ohishi 2018.08.20

  use common_pom_var
  implicit none

  integer timing,i,j
  
  real(kind = r_dble) time0_data,ratio

  if(issf == 1)then

     if(iint == 1) call read_riv_netcdf
     !if(iint == 1) call read_cama_netcdf            
         
     time0_data=rivtime_julday-julday_start
     call readtiming(timing,ratio,rivtime_dayint,time0_data)
     if(iint /= 1 .and. timing == 1)then
        call read_riv_netcdf
        !call read_cama_netcdf
        if(my_task == master_task)then
           write(6,'(a,i10,2f12.5)') &
                & 'river (timing/ratio/time): ',timing,ratio,dti*float(iint-1)/86400.d0+time0
        end if
     end if

     !$omp parallel
     !$omp do private(i,j)
     do j=1,jm
        do i=1,im
           
           river(i,j)=fsm(i,j)*((1.d0-ratio)*river0(i,j)+ratio*river1(i,j))
           fflux(i,j)=fsm(i,j)*(evap(i,j)-prep(i,j)-river(i,j))

           !     unit:[mm/day] --> [psu m/sec]
           !     S.Ohishi sign change 2020.02
           rflux(i,j)=fsm(i,j)*river(i,j)*s(i,j,1)*(-1.d-3)/86400.d0

           ! surface salt flux [psu m/sec]
           !  unit[mm/day] -> [m/sec]
           wssurf(i,j)=fsm(i,j)*fflux(i,j)*s(i,j,1)*(-1.d-3)/86400.d0
           
        end do
     end do
     !$omp end do
     !$omp end parallel

  else if(issf == 2)then

     river(1:im,1:jm)=0.d0
     fflux(1:im,1:jm)=0.d0
     rflux(1:im,1:jm)=0.d0
     wssurf(1:im,1:jm)=0.d0
     
  else
     
     if(my_task == master_task) &
          & write(*,'(/a)') 'POM terminated with error (issf)'
     call finalize_mpi
     stop
     
  end if

end subroutine surface_riverflux

!_______________________________________________________________________
subroutine latenth(qe,ts,ta,qa,qs,ua,va,im,jm)

  !$use omp_lib  
  use common_pom_var, only: r_size
  implicit none
  
  !     upward(sea->air): positive
  !     ts: sea surface temperature [K]
  !     ta: air temperature above 2m [k]
  !     qa: specific humidity [g/g]
  !     ua,va: wind velocity above 10m
  
  ! -- argumants
  integer,intent(in) :: im,jm
  real(kind = r_size),intent(in) :: ts(im,jm),ta(im,jm),qa(im,jm),ua(im,jm),va(im,jm)

  real(kind = r_size),intent(out) :: qe(im,jm),qs(im,jm)

  ! -- local
  integer i,j
  real(kind = r_size),parameter :: z10=10.d0,z2=1.5d0,rkalm=0.4d0
  real(kind = r_size),parameter :: rhoa=1.2d0,rL=2.501d6,Pa=1013.d0,Ce=1.1d-3
  real(kind = r_size) ua2m(im,jm),va2m(im,jm),wmag(im,jm)
  real(kind = r_size) ce10m(im,jm),ce2m(im,jm)
  real(kind = r_size) cd10m(im,jm),cd2m(im,jm)
  real(kind = r_size) devide,z0,esat_const
    
  !     bulk coefficients

  call bulkcof(cd10m,ua,va,ts,ta,wmag,im,jm,1)
  call bulkcof(ce10m,ua,va,ts,ta,wmag,im,jm,3)

  !$omp parallel
  !$omp do private(i,j,devide,z0)
  do j=1,jm
     do i=1,im

        !      transformation from Ce(10m) to Ce(2m)

        if(cd10m(i,j) > 0.d0)then
           cd2m(i,j) = rkalm*rkalm & 
                & *(rkalm*cd10m(i,j)**(-0.5)-log(z10/z2))**(-2)
        else
           cd2m(i,j) = 0.0
        end if
        
        if(ce10m(i,j) > 0.d0) then
           devide = rkalm*sqrt(cd10m(i,j))/ce10m(i,j)+log(z2/z10) 
        else
           devide = 0.d0
        end if
        
        if( devide > 0.d0) then
           ce2m(i,j) = rkalm*sqrt(cd2m(i,j))/devide
        else
           ce2m(i,j) = 0.d0
        end if

        !      transformation from U(10m) to U(2m)

        if(cd10m(i,j) > 0.d0)then
           z0 = exp( log(z10)-rkalm*( Cd10m(i,j) )**(-0.5) )
           if( z0 > 1.d-7 ) then
              devide = log(z10/z0)
           else
              z0 = 0.d0
              devide = 0.d0
           end if
        else
           z0 = 0.d0
           devide = 0.d0
        end if
        if( devide > 1.d-7 .and. z0 > 1.d-7 ) then
           ua2m(i,j) = ua(i,j)*log(z2/z0)/devide
           va2m(i,j) = va(i,j)*log(z2/z0)/devide
        else
           ce2m(i,j) = 0.d0
           ua2m(i,j) = 0.d0
           va2m(i,j) = 0.d0
        end if
        
     end do
  end do
  !$omp end do

  !$omp do private(i,j,esat_const)
  do j=1,jm
     do i=1,im

        !               qs = 0.622*esat(ts(i,j))/Pa
        !     &             /(1.-0.378*esat(ts(i,j))/Pa)
        esat_const = &
             & 6.1078d0*10.d0**(7.5d0*(ts(i,j)-273.15d0)/(237.3d0+ts(i,j)-273.15d0))        

        !qs(i,j) = 0.622*esat_const/(Pa*(1.-0.378*esat_const/Pa))
        qs(i,j) = 0.622d0*esat_const/(Pa-0.378d0*esat_const)
        qe(i,j) = rL*rhoa*ce2m(i,j)*sqrt(ua2m(i,j)*ua2m(i,j)+va2m(i,j)*va2m(i,j))*(qs(i,j)-qa(i,j))

     enddo
  enddo
  !$omp end do
  !$omp end parallel

end subroutine latenth

!_______________________________________________________________________
subroutine sensibleh(qh,ts,ta,ua,va,im,jm)

  !$use omp_lib  
  use common_pom_var, only: r_size  
  implicit none
  
  !     upward(sea->air): positive
  !     ts: sea surface temperature [K]
  !     ta: air temperature above 2m [k]

  ! -- arugumanets

  integer,intent(in) :: im,jm
  real(kind = r_size),intent(out) :: qh(im,jm)
  real(kind = r_size),intent(in) :: ts(im,jm),ta(im,jm),ua(im,jm),va(im,jm)
  
  ! -- local
  integer i,j
  real(kind = r_size),parameter :: z10=10.0d0,z2=1.5d0,rkalm=0.4d0
  real(kind = r_size),parameter :: rhoa=1.2d0,cp=1.005d3,chc=1.1d-3
  real(kind = r_size) wmag(im,jm)
  real(kind = r_size) ch10m(im,jm),ch2m(im,jm)
  real(kind = r_size) cd10m(im,jm),cd2m(im,jm)
  real(kind = r_size) ua2m(im,jm),va2m(im,jm)
  real(kind = r_size) z0

  call bulkcof(cd10m,ua,va,ts,ta,wmag,im,jm,1)
  call bulkcof(ch10m,ua,va,ts,ta,wmag,im,jm,2)

  !      transformation from Ch(10m) to Ch(2m)

  !$omp parallel
  !$omp do private(i,j,z0)  
  do j=1,jm
     do i=1,im

        !      transformation from Ch(10m) to Ch(2m)

        if(cd10m(i,j) > 0.d0)then
           cd2m(i,j)=rkalm*rkalm*(rkalm*cd10m(i,j)**(-0.5d0)-log(z10/z2))**(-2.d0)
        else
           cd2m(i,j)=0.d0
        end if
        
        if(ch10m(i,j) > 0.d0)then
           ch2m(i,j)=rkalm*sqrt(cd2m(i,j))/(rkalm*sqrt(cd10m(i,j))/ch10m(i,j)+log(z2/z10))
        else
           ch2m(i,j)=0.d0
        end if

        !      transformation from U(10m) to U(2m)

        if(cd10m(i,j) > 0.d0)then
           z0 = exp( log(z10)-rkalm*( cd10m(i,j) )**(-0.5d0) )
        else
           z0 = 0.d0
        end if
        
        if(z0 > 1.d-7)then
           ua2m(i,j) = ua(i,j)*log(z2/z0)/log(z10/z0)
           va2m(i,j) = va(i,j)*log(z2/z0)/log(z10/z0)
        else
           ua2m(i,j) = 0.d0
           va2m(i,j) = 0.d0
        end if

     enddo
  enddo
  !$omp end do
  
  !$omp do private(i,j) 
  do j=1,jm
     do i=1,im
        !  Kondo(1975)
        qh(i,j)=cp*rhoa*ch2m(i,j) &
             & *sqrt(ua2m(i,j)*ua2m(i,j)+va2m(i,j)*va2m(i,j))*(ts(i,j)-ta(i,j)) 
        
     enddo
  enddo
  !$omp end do
  !$omp end parallel

end subroutine sensibleh

!_______________________________________________________________________
subroutine wind_stress(wusurf,wvsurf,windu,windv,im,jm)
  ! set wind stress  (Large and Pond, 1981)

  !$use omp_lib  
  use common_pom_var, only: r_size
  implicit none
  
  ! -- arguments
  integer,intent(in) :: im,jm
  real(kind = r_size),intent(in) :: windu(im,jm),windv(im,jm)
  real(kind = r_size),intent(out) :: wusurf(im,jm),wvsurf(im,jm)

  ! -- local      
  integer i,j

  real(kind = r_size),parameter :: cofvmag(3)=(/4.d0,10.d0,26.d0/)
  real(kind = r_size),parameter :: rhoa=1.2d0
  real(kind = r_size) wvel,cd

  !$omp parallel
  !$omp do private(i,j,wvel,cd)      
  do j=1,jm
     do i=1,im
        wvel=sqrt(windu(i,j)*windu(i,j)+windv(i,j)*windv(i,j))
        if(wvel < cofvmag(2))then
           cd=1.14d0
        else if(wvel < cofvmag(3))then
           cd=0.49d0+0.065d0*wvel
        else
           cd=0.49d0+0.065d0*cofvmag(3)
        end if
        cd=cd*1.d-3
        wusurf(i,j)=rhoa*cd*wvel*windu(i,j)
        wvsurf(i,j)=rhoa*cd*wvel*windv(i,j)
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine wind_stress

!_______________________________________________________________________
subroutine wind_stress_mb(wusurf,wvsurf,windu,windv,im,jm)
  !     set wind stress  (Mellor and Blumberg, 2004)

  !$use omp_lib  
  use common_pom_var, only: r_size
  implicit none
  
  ! -- arguments
  integer,intent(in) :: im,jm
  real(kind = r_size),intent(in) :: windu(im,jm),windv(im,jm)

  real(kind = r_size),intent(out) :: wusurf(im,jm),wvsurf(im,jm)

  ! -- local      
  integer i,j
  real(kind = r_size),parameter :: rhoa=1.2d0
  real(kind = r_size) wvel,cd

  !$omp parallel
  !$omp do private(i,j,wvel,cd)    
  do j=1,jm
     do i=1,im
        wvel=sqrt(windu(i,j)*windu(i,j)+windv(i,j)*windv(i,j))
        if( wvel < 3.7313d0 ) then
           cd=1.5-wvel*0.134002
        else
           cd=0.75+0.067*wvel
           cd=min(2.492,cd)
        end if
        cd=cd*1.d-3
        wusurf(i,j)=rhoa*cd*wvel*windu(i,j)
        wvsurf(i,j)=rhoa*cd*wvel*windv(i,j)
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine wind_stress_mb

!_______________________________________________________________________
subroutine bulkcof(cf,u,v,ts,ta,w,id,jd,ind)
  
  !***********************************************************************
  !*          Bulk coefficient applied to turbulent fluxes               *
  !*                        Kondo(1975)                                  *
  !*
  !*          IND=1 ----> CDD   for momentum flux
  !*          IND=2 ----> CHD   for sensible heat flux
  !*          IND=3 ----> CED   for latent heat flux
  !*
  !*          coefficient        in the nutral case
  !*          CDD(ID,JD,MD)      CD(ID,JD,MD)
  !*          CHD(ID,JD,MD)      CH(ID,JD,MD)
  !*          CED(ID,JD,MD)      CE(ID,JD,MD)
  !*
  !*     Input data      Coefficient which depends on the wind speed
  !*                     in the nutral case  
  !*
  !*                     CD(ID,JD,MD)
  !*     W(ID,JD,MD)---->CH(ID,JD,MD)
  !*                     CE(ID,JD,MD) 
  !*
  !*               the stability conditions  stability parameter 
  !*
  !* TS(ID,JD,MD)    (TS-TA)<0      stable   S<0   S=S0*(|S0|/(|S0|+0.01))
  !* TA(ID,JD,MD)--->(TS-TA)>0    unstable   S>0   S0=(TS-TA)/W**2
  !*                 (TS-TA)=0      nutral   S=0   If W=0, S=0
  !*
  !************************************************************************

  !$use omp_lib  
  use common_pom_var, only: r_size
  implicit none

  integer i,j
  
  integer,intent(in) :: id,jd,ind
  real(kind = r_size),intent(in) :: ts(id,jd),ta(id,jd)
  real(kind = r_size),intent(in) :: u(id,jd),v(id,jd)

  real(kind = r_size),intent(out) :: w(id,jd)
  real(kind = r_size),intent(out) :: cf(id,jd)

  integer kkk
  real(kind = r_size) cf0,s0,s
  real(kind = r_size) A(3,5),B(3,5),C(3,5),P(3,5),RR(3)

!---------------Coefficients of (CD, CH, CE) formula 

  data A /0.d0, 0.d0, 0.d0, 0.771d0, 0.927d0, 0.969d0, 0.867d0, 1.15d0, 1.18d0, &
       & 1.2d0, 1.17d0, 1.196d0, 0.d0, 1.652d0, 1.68d0/

  data B /1.08d0, 1.185d0, 1.23d0, 0.0858d0, 0.0546d0, 0.0521d0, 0.0667d0, 0.01d0, 0.01d0, &
       & 0.025d0, 0.0075d0, 0.008d0, 0.073d0, -0.017d0, -0.016d0/
  data C /0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, & 
       & -0.00045d0, -0.0004d0, 0.d0, 0.d0, 0.d0, 0.d0/
  data P /-0.15d0, -0.157d0, -0.16d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, & 
       & 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0/
  data RR /0.47d0, 0.63d0, 0.63d0/

!--------------------Write  coefficients
!
!      WRITE(6,*) 'Coefficients in the nutral case calculations' 
!      DO 1000 I=1, 3
!      WRITE(6,100) (A(I,J),J=1,5)
!      WRITE(6,100) (B(I,J),J=1,5)
!      WRITE(6,100) (C(I,J),J=1,5)
!      WRITE(6,100) (P(I,J),J=1,5)
! 1000 CONTINUE
!  100 FORMAT(1H , 5(F8.5,1x)/)


  !$omp parallel
  !$omp do private(i,j,kkk,cf0,s0,s)  
  do i=1,id
     do j=1,jd

        w(i,j) = sqrt( u(i,j)*u(i,j) + v(i,j)*v(i,j) )

!-----------------CF0 for the nutral case

        if(w(i,j) <= 2.2d0)then
           kkk=1
        else if(w(i,j) <= 5.d0)then
           kkk=2
        else if(w(i,j) <= 8.d0)then
           kkk=3
        else if(w(i,j) <= 25.d0)then
           kkk=4
        else
           kkk=5
        endif
 
!----------------The calculation is possible for W lager than 0.3 m/s, 
!                based on Kondo(1975) formula.
!                MSK cannot be applied for all valiables.

        if(w(i,j) < 0.3d0)then
           cf(i,j)=1.1d-3
           exit
        endif

        cf0=1.d-3*( A(ind,kkk)+B(ind,kkk)*(w(i,j)**P(ind,kkk)) &
             & +C(ind,kkk)*( w(i,j)-8.d0 )**2 )

!------------------S for stability parameter

        if(w(i,j) /= 0.d0)then
           s0=( ts(i,j)-ta(i,j) )/(w(i,j)**2)
           s=s0*( abs(s0)/(abs(s0)+0.01d0) )
        else
           s=0.d0
        end if

!-----------------CF(I,J,M) for bulk coefficient
!-----------------for the stable case

        if(s < -3.3d0)then
           cf(i,j)=0.d0
        else if(s < 0.d0)then
           cf(i,j)=cf0*( 0.1d0+0.03d0*s+0.9d0*exp(4.8d0*s) )

!-----------------for the nutral case
        else if(s == 0.d0)then
           cf(i,j)=cf0
           
!-----------------for the unstable case
           
        else if(s > 0.d0) THEN
           cf(i,j)=cf0*( 1.0d0+RR(ind)*s**0.5d0 )
           
        end if

     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine bulkcof
      
!_______________________________________________________________________
subroutine longrad1(ql,ts,ta,qa,tcdc,im,jm)

  !$use omp_lib  
  use common_pom_var, only: r_size
  implicit none

  !     Eq. (13.6) in Kagimoto et al. (2008) 
  !     http://engan.cmes.ehime-u.ac.jp/xguo/paper/Kagimoto2008book.pdf

  !     upward(sea->air): positive
  !     ts: sea surface temperature [K]
  !     ta: air temperature above 2m [k]
  !     qa: specific humidity [g/g]
  !     tcdc: total clound cover [0:1]
  !     downward: positive
  !
  ! -- arguments
  integer,intent(in) :: im,jm
  real(kind = r_size),intent(in) :: ts(im,jm),ta(im,jm),qa(im,jm),tcdc(im,jm)
  real(kind = r_size),intent(out) :: ql(im,jm)
  
  ! -- local
  integer i,j
  real(kind = r_size) ea

  real(kind = r_size),parameter :: eps = 0.97d0        ! emissivity of the ocean     
  real(kind = r_size),parameter :: sigma = 5.670d-8  ! Stefan-Boltzmann constant [W m^-2 K^-4]
  real(kind = r_size),parameter :: bb = 0.8d0          ! linear correction factor
  real(kind = r_size),parameter :: Pa = 1013.d0        ! atmospheric pressure on sea surface

  !$omp parallel
  !$omp do private(i,j,ea)
  do j=1,jm
     do i=1,im
        ea = qa(i,j)*Pa/(0.622d0+0.378d0*qa(i,j))
        ql(i,j) = eps*sigma*ts(i,j)**4 &
             & *( 0.39d0 - 0.05d0*sqrt( ea ) ) &
             & *( 1.d0 - bb*tcdc(i,j) ) &
             & + 4.d0*eps*sigma*ts(i,j)**3*(ts(i,j)-ta(i,j))

     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine longrad1

!_____________________________________________________________________

subroutine longrad2(lwrad,sst,im,jm)

  !$use omp_lib  
  use common_pom_var, only: r_size
  implicit none

  !Eq.(8) in Tsujino et al. (2018) Ocean Modelling, Vol: 130, pp 79-139

  !Common
  integer i,j
  real(kind = r_size),parameter :: eps=1.d0
  real(kind = r_size),parameter :: sbc=5.67d-8 !Stefan-Boltzman constant [W m-2 K-4]

  !IN
  integer,intent(in) :: im,jm
  real(kind = r_size),intent(in) :: sst(im,jm) ![K]
      
  !INOUT
  !Positive: Upward
  real(kind = r_size),intent(inout) :: lwrad(im,jm) !Downward --> NET

  !$omp parallel
  !$omp do private(i,j)
  do j=1,jm
     do i=1,im
        !     Downward --> NET (Positive: Downward)
        lwrad(i,j) = & 
             & lwrad(i,j) - eps * sbc * sst(i,j)**4.d0
        !     Positive: Upward
        lwrad(i,j) = -1.d0*lwrad(i,j)
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine longrad2

!_______________________________________________________________________

subroutine shortrad(im,jm,pi,north_e,swrad)

  !$use omp_lib  
  use common_pom_var,only: r_size
  implicit none
      
  !Eq. (7) Tsujino et al. (2018)
      
  !Common
  integer i,j
  real(kind = r_size) albedo !Eq. (8) Largen and Yeager (2009)

  !IN
  integer,intent(in) :: im,jm
  real(kind = r_size),intent(in) :: pi
  real(kind = r_size),intent(in) :: north_e(im,jm)
      
  !INOUT
  real(kind = r_size),intent(inout) :: swrad(im,jm) !Downward --> NET

  !$omp parallel
  !$omp do private(i,j,albedo)
  do j=1,jm
     do i=1,im
        albedo=0.069d0-0.011d0*cos(2.d0*north_e(i,j)*pi/180.d0)
        swrad(i,j)=(1.d0-albedo)*swrad(i,j)
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine shortrad
