! advance.f
! advance POM
!_______________________________________________________________________
subroutine advance
  ! advance POM 1 step in time

  use common_pom_var
  implicit none
  
  if(my_task == master_task .and. mod(iint,100) == 0)then
     write(6,'(a,i10)') "start advance routine at ",iint
  end if
  
  ! get time
  call get_time
  
  ! set time dependent surface boundary conditions
  call surface_forcing
  
  ! set time dependent lateral boundary conditions
  if(ilbc == 1)then
     call lateral_forcing_mclim
  elseif(ilbc == 2)then
     call lateral_forcing
  else
     write(5,"(/a)") "POM terminated with error: ilbc"
     call finalize_mpi
     stop
  endif
  
  ! set time dependent reference temperature and salinity
  if(ilbc == 1)then
     call set_tsdata_mclim
  elseif(ilbc == 2)then
     call set_tsdata
  endif
  
  ! set monthly climatological data
  call set_tsclim

  ! set lateral viscosity
  call lateral_viscosity
  
  ! form vertical averages of 3-D fields for use in external (2-D) mode
  call mode_interaction
  
  ! external (2-D) mode calculation
  do iext=1,isplit
     call mode_external
  end do
  
  ! internal (3-D) mode calculation
  call mode_internal
  
  ! print section
  call print_section
  
  ! check CFL condition
  call check_velocity

  !check NaN condition
  !call check_nan("ua",im,jm,1,ua)
  !call check_nan("va",im,jm,1,va)
  !call check_nan("el",im,jm,1,el)
  !call check_nan("t",im,jm,kb,t)
  !call check_nan("s",im,jm,kb,s)
  !call check_nan("u",im,jm,kb,u)
  !call check_nan("v",im,jm,kb,v)
  
  ! daily average (S.Ohishi 2018.07)
  if(idave == 1)then
     call daily_average
  end if

  if(assim == 1)then
     call whole_average
  end if

  ! write netcdf
  if(netcdf_file /= 'nonetcdf' .and. mod(iint,iprint) == 0 .and. assim /= 1)then
     call write_output_netcdf
  end if
  
  ! write restart
  if(mod(iint,irestart) == 0 .and. assim /= 1)then
     call write_restart_netcdf
  end if

  ! write forecast for ensemble DA
  if(iint == iend .and. assim == 1)then
     call write_iau_netcdf
  end if
    
  if(my_task == master_task .and. mod(iint,100) == 0)then
     write(6,'(a,i10)') "finish advance routine at ",iint     
  end if
  
end subroutine advance

!_______________________________________________________________________
subroutine get_time
! return the model time

  use common_pom_var
  implicit none

  time=dti*dble(iint)/86400.d0+time0
  if(iint >= iswtch) iprint=nint(prtd2*24.d0*3600.d0/dti)

  if(lramp) then
     ramp=time/period
     if(ramp > 1.d0) ramp=1.d0
  else
     ramp=1.d0
  end if

end subroutine get_time

!_______________________________________________________________________
subroutine surface_forcing
  ! set time dependent surface boundary conditions

  use common_pom_var
  implicit none

  
  call surface_airseaflux
  
  call surface_riverflux

  ! wind stress
  ! value is negative for westerly or southerly winds. The wind stress
  ! should be tapered along the boundary to suppress numerically induced
  ! oscilations near the boundary (Jamart and Ozer, JGR, 91, 10621-10631)
  !          wusurf(i,j)=0.d0
  !          wvsurf(i,j)=0.d0
  
  e_atmos(2:imm1,2:jmm1)=0.d0
  vfluxf(2:imm1,2:jmm1)=0.d0

  ! set w(i,j,1)=vflux(i,j).ne.0 if one wishes non-zero flow across
  ! the sea surface. See calculation of elf(i,j) below and subroutines
  ! vertvl, advt1 (or advt2). If w(1,j,1)=0, and, additionally, there
  ! is no net flow across lateral boundaries, the basin volume will be
  ! constant; if also vflux(i,j).ne.0, then, for example, the average
  ! salinity will change and, unrealistically, so will total salt
  w(2:imm1,2:jmm1,1)=vfluxf(2:imm1,2:jmm1)

  ! set wtsurf to the sensible heat, the latent heat (which involves
  ! only the evaporative component of vflux) and the long wave
  ! radiation
  !          wtsurf(i,j)=0.d0
  !
  ! set swrad to the short wave radiation
  !          swrad(i,j)=0.d0

  ! to account for change in temperature of flow crossing the sea
  ! surface (generally quite small compared to latent heat effect)
  !          tatm=t(i,j,1)+tbias    ! an approximation
  !          wtsurf(i,j)=wtsurf(i,j)+vfluxf(i,j)*(tatm-t(i,j,1)-tbias)
  !          wtsurf(i,j)=50./4.1876d6*(tclim(i,j,1)-t(i,j,1))
  !
  ! set the salinity of water vapor/precipitation which enters/leaves
  ! the atmosphere (or e.g., an ice cover)
  !          satm=0.d0
  !          wssurf(i,j)=            vfluxf(i,j)*(satm-s(i,j,1)-sbias)
        
  
end subroutine surface_forcing

!_______________________________________________________________________
subroutine lateral_forcing_mclim
  !     set monthly climatological lateral boudary forcing
  !     Created by 2018.08 S.Ohishi

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k
  integer julday_time,iyr,imon,iday,isecflg
  real(kind = r_dble) timepre,sec,ratio

  timepre=dti*float(iint-1)/86400.d0+time0
  julday_time=nint(julday_start+timepre)
  call caldat(julday_time,iyr,imon,iday)

  if(iint == 1)then
     if(iday < 15)then
        imon=imon-1
        if(imon < 1) imon=12
     end if
     call read_lbc_mclim_netcdf(imon)
  else
     sec=(timepre-int(timepre))*86400.d0
     isecflg=int(abs(sec))
     if(iday == 15 .and. isecflg == 0)then
        call read_lbc_mclim_netcdf(imon)
     end if
  end if

  call ratioclim(ratio,julday_time,timepre)
  if(mod(iint-1,100) == 0 .and. my_task == master_task)then
     write(6,'(a,f12.5,a,i10)') 'lbcclim ratio:',ratio," at ",iint
  end if

  !$omp parallel
  !$omp do private(j)
  do j=2,jmm1
     ele(j)=(1.d0-ratio)*ele0(j)+ratio*ele1(j)
     uabe(j)=(1.d0-ratio)*uabe0(j)+ratio*uabe1(j)
     vabe(j)=(1.d0-ratio)*vabe0(j)+ratio*vabe1(j)
     elw(j)=(1.d0-ratio)*elw0(j)+ratio*elw1(j)
     uabw(j)=(1.d0-ratio)*uabw0(j)+ratio*uabw1(j)
     vabw(j)=(1.d0-ratio)*vabw0(j)+ratio*vabw1(j)
  end do
  !$omp end do

  !$omp do private(j,k)  
  do k=1,kbm1
     do j=2,jmm1
        tbe(j,k)=(1.d0-ratio)*tbe0(j,k)+ratio*tbe1(j,k)
        sbe(j,k)=(1.d0-ratio)*sbe0(j,k)+ratio*sbe1(j,k)
        ube(j,k)=(1.d0-ratio)*ube0(j,k)+ratio*ube1(j,k)
        vbe(j,k)=(1.d0-ratio)*vbe0(j,k)+ratio*vbe1(j,k)
        tbw(j,k)=(1.d0-ratio)*tbw0(j,k)+ratio*tbw1(j,k)
        sbw(j,k)=(1.d0-ratio)*sbw0(j,k)+ratio*sbw1(j,k)
        ubw(j,k)=(1.d0-ratio)*ubw0(j,k)+ratio*ubw1(j,k)
        vbw(j,k)=(1.d0-ratio)*vbw0(j,k)+ratio*vbw1(j,k)
     end do
  end do
  !$omp end do

  !$omp do private(i)    
  do i=2,imm1
     eln(i)=(1.d0-ratio)*eln0(i)+ratio*eln1(i)
     uabn(i)=(1.d0-ratio)*uabn0(i)+ratio*uabn1(i)
     vabn(i)=(1.d0-ratio)*vabn0(i)+ratio*vabn1(i)
     els(i)=(1.d0-ratio)*els0(i)+ratio*els1(i)
     uabs(i)=(1.d0-ratio)*uabs0(i)+ratio*uabs1(i)
     vabs(i)=(1.d0-ratio)*vabs0(i)+ratio*vabs1(i)
  end do
  !$omp end do

  !$omp do private(i,k)    
  do k=1,kbm1
     do i=2,imm1
        tbn(i,k)=(1.d0-ratio)*tbn0(i,k)+ratio*tbn1(i,k)
        sbn(i,k)=(1.d0-ratio)*sbn0(i,k)+ratio*sbn1(i,k)
        ubn(i,k)=(1.d0-ratio)*ubn0(i,k)+ratio*ubn1(i,k)
        vbn(i,k)=(1.d0-ratio)*vbn0(i,k)+ratio*vbn1(i,k)
        tbs(i,k)=(1.d0-ratio)*tbs0(i,k)+ratio*tbs1(i,k)
        sbs(i,k)=(1.d0-ratio)*sbs0(i,k)+ratio*sbs1(i,k)
        ubs(i,k)=(1.d0-ratio)*ubs0(i,k)+ratio*ubs1(i,k)
        vbs(i,k)=(1.d0-ratio)*vbs0(i,k)+ratio*vbs1(i,k)
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine lateral_forcing_mclim
      
!_______________________________________________________________________
subroutine lateral_forcing
  ! set time dependent surface boundary conditions

  !$use omp_lib    
  use common_pom_var
  implicit none

  integer timing,i,j,k
  real(kind = r_size) time0_data,ratio

  if(iint == 1) call read_lbc_netcdf
   
  time0_data=lbctime_julday-julday_start
  call readtiming(timing,ratio,lbctime_dayint,time0_data)
  if(timing == 1) then
     call read_lbc_netcdf
     if(my_task == master_task)then
        write(6,'(a,i10,2f12.5)') 'lbc timing ratio time ',timing,ratio,dti*float(iint-1)/86400.d0+time0
     end if
  end if
  
  !$omp parallel
  !$omp do private(j)
  do j=2,jmm1
     ele(j)=(1.d0-ratio)*ele0(j)+ratio*ele1(j)
     uabe(j)=(1.d0-ratio)*uabe0(j)+ratio*uabe1(j)
     vabe(j)=(1.d0-ratio)*vabe0(j)+ratio*vabe1(j)
     elw(j)=(1.d0-ratio)*elw0(j)+ratio*elw1(j)
     uabw(j)=(1.d0-ratio)*uabw0(j)+ratio*uabw1(j)
     vabw(j)=(1.d0-ratio)*vabw0(j)+ratio*vabw1(j)
  end do
  !$omp end do

  !$omp do private(j,k)  
  do k=1,kbm1
     do j=2,jmm1
        tbe(j,k)=(1.d0-ratio)*tbe0(j,k)+ratio*tbe1(j,k)
        sbe(j,k)=(1.d0-ratio)*sbe0(j,k)+ratio*sbe1(j,k)
        ube(j,k)=(1.d0-ratio)*ube0(j,k)+ratio*ube1(j,k)
        vbe(j,k)=(1.d0-ratio)*vbe0(j,k)+ratio*vbe1(j,k)
        tbw(j,k)=(1.d0-ratio)*tbw0(j,k)+ratio*tbw1(j,k)
        sbw(j,k)=(1.d0-ratio)*sbw0(j,k)+ratio*sbw1(j,k)
        ubw(j,k)=(1.d0-ratio)*ubw0(j,k)+ratio*ubw1(j,k)
        vbw(j,k)=(1.d0-ratio)*vbw0(j,k)+ratio*vbw1(j,k)
     end do
  end do
  !$omp end do

  !$omp do private(i)  
  do i=2,imm1
     eln(i)=(1.d0-ratio)*eln0(i)+ratio*eln1(i)
     uabn(i)=(1.d0-ratio)*uabn0(i)+ratio*uabn1(i)
     vabn(i)=(1.d0-ratio)*vabn0(i)+ratio*vabn1(i)
     els(i)=(1.d0-ratio)*els0(i)+ratio*els1(i)
     uabs(i)=(1.d0-ratio)*uabs0(i)+ratio*uabs1(i)
     vabs(i)=(1.d0-ratio)*vabs0(i)+ratio*vabs1(i)
  end do
  !$omp end do

  !$omp do private(i,k)
  do k=1,kbm1
     do i=2,imm1
        tbn(i,k)=(1.d0-ratio)*tbn0(i,k)+ratio*tbn1(i,k)
        sbn(i,k)=(1.d0-ratio)*sbn0(i,k)+ratio*sbn1(i,k)
        ubn(i,k)=(1.d0-ratio)*ubn0(i,k)+ratio*ubn1(i,k)
        vbn(i,k)=(1.d0-ratio)*vbn0(i,k)+ratio*vbn1(i,k)
        tbs(i,k)=(1.d0-ratio)*tbs0(i,k)+ratio*tbs1(i,k)
        sbs(i,k)=(1.d0-ratio)*sbs0(i,k)+ratio*sbs1(i,k)
        ubs(i,k)=(1.d0-ratio)*ubs0(i,k)+ratio*ubs1(i,k)
        vbs(i,k)=(1.d0-ratio)*vbs0(i,k)+ratio*vbs1(i,k)
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine lateral_forcing

!_______________________________________________________________________
subroutine set_tsdata_mclim
  !     monthly climatology temperature and salinity
  !     2018.08.20 S.Ohishi

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k
  integer julday_time,iyr,imon,iday,isecflg
  real(kind = r_dble) timepre,sec,ratio

  timepre=dti*float(iint-1)/86400.d0+time0
  julday_time=nint(julday_start+timepre)
  call caldat(julday_time,iyr,imon,iday)

  if(iint == 1)then
     if(iday < 15)then
        imon=imon-1
        if(imon < 1) imon=12
     end if
     call read_tsdata_mclim_netcdf(imon)
  else
     sec=(timepre-int(timepre))*86400.d0
     isecflg=int(abs(sec))
     if(iday == 15 .and. isecflg == 0)then
        call read_tsdata_mclim_netcdf(imon)
     end if
  end if

  call ratioclim(ratio,julday_time,timepre)
  if(mod(iint-1,100) == 0 .and. my_task == master_task)then
     write(6,'(a,f12.5,a,i10)') 'tsdataclim ratio:',ratio," at ",iint
  end if

  !$omp parallel
  !$omp do private(i,j,k)  
  do k=1,kb
     do j=1,jm
        do i=1,im
           tref(i,j,k)=(1.d0-ratio)*tref0(i,j,k)+ratio*tref1(i,j,k)
           sref(i,j,k)=(1.d0-ratio)*sref0(i,j,k)+ratio*sref1(i,j,k)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine set_tsdata_mclim

!_______________________________________________________________________
subroutine set_tsdata
  ! set time dependent reference temperature and salinity

  !$use omp_lib  
  use common_pom_var
  implicit none

  real(kind =r_dble) time0_data,ratio
  integer timing
  integer i,j,k

  if(iint == 1) call read_tsdata_netcdf
   
  time0_data=tsdatatime_julday-julday_start
  call readtiming(timing,ratio,tsdatatime_dayint,time0_data)
  if(timing == 1) then
     call read_tsdata_netcdf
     if(my_task == master_task)then
        write(6,'(a,i10,f12.5)') 'tsdata timing ratio time ',timing,ratio,dti*float(iint-1)/86400.d0+time0
     end if
  end if

  !$omp parallel
  !$omp do private(i,j,k)  
  do k=1,kb
     do j=1,jm
        do i=1,im
           tref(i,j,k)=(1.d0-ratio)*tref0(i,j,k)+ratio*tref1(i,j,k)
           sref(i,j,k)=(1.d0-ratio)*sref0(i,j,k)+ratio*sref1(i,j,k)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine set_tsdata

!_______________________________________________________________________
subroutine set_tsclim
  ! set monthly climatological temperature and salinity data

  !$use omp_lib  
  use common_pom_var
  implicit none

  real(kind = r_dble) timepre,sec,ratio
  integer i,j,k
  integer julday_time,iyr,imon,iday,isecflg

  timepre=dti*float(iint-1)/86400.d0+time0
  julday_time=nint(julday_start+timepre)
  call caldat(julday_time,iyr,imon,iday)
     
  if(iint == 1) then
     if(iday < 15) then
        imon=imon-1
        if(imon < 1) imon=12
     end if
     call read_tsclim_monthly_netcdf(imon)
  else
     sec=(timepre-int(timepre))*86400.d0
     isecflg=int(abs(sec))
     if(iday == 15 .and. isecflg == 0) then
        call read_tsclim_monthly_netcdf(imon)
     end if
  end if

  call ratioclim(ratio,julday_time,timepre)
  if(mod(iint-1,100) == 0 .and. my_task == master_task)then
     write(6,'(a,f12.5,a,i10)') 'tsclim ratio:', ratio ," at ",iint
  end if

  !$omp parallel
  !$omp do private(i,j,k)    
  do k=1,kb
     do j=1,jm
        do i=1,im
           tclimm(i,j,k)=(1.d0-ratio)*tclim0(i,j,k)+ratio*tclim1(i,j,k)
           sclimm(i,j,k)=(1.d0-ratio)*sclim0(i,j,k)+ratio*sclim1(i,j,k)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine set_tsclim

!_______________________________________________________________________
subroutine lateral_viscosity
  ! set the lateral viscosity

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k
  ! if mode=2 then initial values of aam2d are used. If one wishes
  ! to use Smagorinsky lateral viscosity and diffusion for an
  ! external (2-D) mode calculation, then appropiate code can be
  ! adapted from that below and installed just before the end of the
  ! "if(mode.eq.2)" loop in subroutine advave

  ! calculate Smagorinsky lateral viscosity:
  ! ( hor visc = horcon*dx*dy*sqrt((du/dx)**2+(dv/dy)**2
  !                                +.5*(du/dy+dv/dx)**2) )
  if(mode /= 2)then

     call advct

     if (npg == 1) then
        call baropg
     else if (npg == 2) then
        call baropg_mcc
     else if (npg == 3) then
        call baropg_thiem
     else
        error_status=1
        write(6,'(/''Error: invalid value for npg'')')
     end if
     
     !$omp parallel
     !$omp do private(i,j,k)
     do k=1,kbm1
        do j=2,jmm1
           do i=2,imm1
              aam(i,j,k)= &
                   & horcon*dx(i,j)*dy(i,j) &
                   & *sqrt( ((u(i+1,j,k)-u(i,j,k))/dx(i,j))**2+((v(i,j+1,k)-v(i,j,k))/dy(i,j))**2 &
                   & +0.5d0*(0.25d0*(u(i,j+1,k)+u(i+1,j+1,k)-u(i,j-1,k)-u(i+1,j-1,k))/dy(i,j) &
                   & +0.25d0*(v(i+1,j,k)+v(i+1,j+1,k)-v(i-1,j,k)-v(i-1,j+1,k))/dx(i,j)) **2) &
                   & +aamadd
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel

     call exchange3d_mpi(aam(:,:,1:kbm1),im,jm,kbm1)

  end if

end subroutine lateral_viscosity

!_______________________________________________________________________
subroutine mode_interaction

  ! form vertical averages of 3-D fields for use in external (2-D) mode

  !$use omp_lib  
  use common_pom_var
  implicit none
  
  integer i,j,k

  if(mode /= 2) then

     adx2d(1:im,1:jm)=0.d0
     ady2d(1:im,1:jm)=0.d0
     drx2d(1:im,1:jm)=0.d0
     dry2d(1:im,1:jm)=0.d0
     aam2d(1:im,1:jm)=0.d0

     do k=1,kbm1
        do j=1,jm
           do i=1,im
              adx2d(i,j)=adx2d(i,j)+advx(i,j,k)*dz(i,j,k)
              ady2d(i,j)=ady2d(i,j)+advy(i,j,k)*dz(i,j,k)
              drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(i,j,k)
              dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(i,j,k)
              aam2d(i,j)=aam2d(i,j)+aam(i,j,k)*dz(i,j,k)
           end do
        end do
     end do

     call exchange2d_mpi(aam2d,im,jm) !2018.08

     call advave

     !$omp parallel
     !$omp do private(i,j)      
     do j=1,jm
        do i=1,im
           adx2d(i,j)=adx2d(i,j)-advua(i,j)
           ady2d(i,j)=ady2d(i,j)-advva(i,j)
        end do
     end do
     !$omp end do
     !$omp end parallel

     call exchange2d_mpi(drx2d,im,jm) !2018.08
     call exchange2d_mpi(dry2d,im,jm) !2018.08     
     call exchange2d_mpi(adx2d,im,jm) !2018.08
     call exchange2d_mpi(ady2d,im,jm) !2018.08

  end if

  !$omp parallel
  !$omp do private(i,j)      
  do j=1,jm
     do i=1,im
        egf(i,j)=el(i,j)*ispi
     end do
  end do
  !$omp end do

  !$omp do private(i,j)        
  do j=1,jm
     do i=2,im
        utf(i,j)=ua(i,j)*(d(i,j)+d(i-1,j))*isp2i
     end do
  end do
  !$omp end do

  !$omp do private(i,j)        
  do j=2,jm
     do i=1,im
        vtf(i,j)=va(i,j)*(d(i,j)+d(i,j-1))*isp2i
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  call exchange2d_mpi(egf,im,jm) !2018.08
  call exchange2d_mpi(utf,im,jm) !2018.08
  call exchange2d_mpi(vtf,im,jm) !2018.08
  
end subroutine mode_interaction

!_______________________________________________________________________
subroutine mode_external
  ! calculate the external (2-D) mode

  !$use omp_lib
  use common_pom_var
  implicit none

  integer i,j

  !$omp parallel
  !$omp do private(i,j)  
  do j=2,jm
     do i=2,im
        fluxua(i,j)=0.25d0*(d(i,j)+d(i-1,j))*(dy(i,j)+dy(i-1,j))*ua(i,j)
        fluxva(i,j)=0.25d0*(d(i,j)+d(i,j-1))*(dx(i,j)+dx(i,j-1))*va(i,j)
     end do
  end do
  !$omp end do
  
  ! NOTE addition of surface freshwater flux, w(i,j,1)=vflux, compared
  ! with pom98.f. See also modifications to subroutine vertvl

  !$omp do private(i,j)  
  do j=2,jmm1
     do i=2,imm1
        elf(i,j)=elb(i,j) &
             & +dte2*(-(fluxua(i+1,j)-fluxua(i,j)+fluxva(i,j+1)-fluxva(i,j))/art(i,j) &
             & -vfluxf(i,j))
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  !  if(ltide) then
  !     call tide_bcond(1)
  !  else

  call bcond(1)
  
  !  end if
  
  call exchange2d_mpi(elf,im,jm)
  
  ! kii2b
  !      if(mod(iext,ispadv).eq.0) call advave(tps)
  if(mod(iext,ispadv) == 0) call advave
  
  !$omp parallel
  !$omp do private(i,j)    
  do j=2,jmm1
     do i=2,im
        uaf(i,j)= &
             & adx2d(i,j)+advua(i,j) &
             & -aru(i,j)*0.25d0*(cor(i,j)*d(i,j)*(va(i,j+1)+va(i,j)) &
             & +cor(i-1,j)*d(i-1,j)*(va(i-1,j+1)+va(i-1,j))) &
             & +0.25d0*grav*(dy(i,j)+dy(i-1,j))*(d(i,j)+d(i-1,j)) &
             & *((1.d0-2.d0*alpha)*(el(i,j)-el(i-1,j))+alpha*(elb(i,j)-elb(i-1,j)+elf(i,j)-elf(i-1,j))+e_atmos(i,j)-e_atmos(i-1,j)) &
             & +drx2d(i,j)+aru(i,j)*(wusurf(i,j)-wubot(i,j))
     end do
  end do
  !$omp end do

  !$omp do private(i,j)  
  do j=2,jmm1
     do i=2,im
        uaf(i,j)= &
             & ((h(i,j)+elb(i,j)+h(i-1,j)+elb(i-1,j))*aru(i,j)*uab(i,j)-4.d0*dte*uaf(i,j)) &
             & /((h(i,j)+elf(i,j)+h(i-1,j)+elf(i-1,j)) &
             & *aru(i,j))
     end do
  end do
  !$omp end do

  !$omp do private(i,j)  
  do j=2,jm
     do i=2,imm1
        vaf(i,j)= &
             & ady2d(i,j)+advva(i,j) &
             & +arv(i,j)*0.25d0*(cor(i,j)*d(i,j)*(ua(i+1,j)+ua(i,j)) &
             & +cor(i,j-1)*d(i,j-1)*(ua(i+1,j-1)+ua(i,j-1))) &
             & +0.25d0*grav*(dx(i,j)+dx(i,j-1))*(d(i,j)+d(i,j-1)) &
             & *((1.d0-2.d0*alpha)*(el(i,j)-el(i,j-1))+alpha*(elb(i,j)-elb(i,j-1)+elf(i,j)-elf(i,j-1))+e_atmos(i,j)-e_atmos(i,j-1)) &
             & +dry2d(i,j)+arv(i,j)*(wvsurf(i,j)-wvbot(i,j))
     end do
  end do
  !$omp end do

  !$omp do private(i,j)    
  do j=2,jm
     do i=2,imm1
        vaf(i,j)= &
             & ((h(i,j)+elb(i,j)+h(i,j-1)+elb(i,j-1))*vab(i,j)*arv(i,j)-4.d0*dte*vaf(i,j)) &
             & /((h(i,j)+elf(i,j)+h(i,j-1)+elf(i,j-1)) &
             & *arv(i,j))
     end do
  end do
  !$omp end do
  !$omp end parallel

  !  if(ltide) then
  !     call tide_bcond(2)
  !  else
  
  call bcond(2)
  
  !  end if
  
  call exchange2d_mpi(uaf,im,jm)
  call exchange2d_mpi(vaf,im,jm)
  
  !$omp parallel  
  if(iext == (isplit-2))then

     !$omp do private(i,j)
     do j=1,jm
        do i=1,im
           etf(i,j)=0.25d0*smoth*elf(i,j)
        end do
     end do
     !$omp end do
     
  else if(iext == (isplit-1)) then

     !$omp do private(i,j)     
     do j=1,jm
        do i=1,im
           etf(i,j)=etf(i,j)+0.5d0*(1.d0-0.5d0*smoth)*elf(i,j)
        end do
     end do
     !$omp end do
     
  else if(iext == isplit) then

     !$omp do private(i,j)     
     do j=1,jm
        do i=1,im
           etf(i,j)=(etf(i,j)+0.5d0*elf(i,j))*fsm(i,j)
        end do
     end do
     !$omp end do
     
  end if
    
  ! apply filter to remove time split
  !$omp do private(i,j)       
  do j=1,jm
     do i=1,im
        ua(i,j)=ua(i,j)+0.5d0*smoth*(uab(i,j)-2.d0*ua(i,j)+uaf(i,j))
        va(i,j)=va(i,j)+0.5d0*smoth*(vab(i,j)-2.d0*va(i,j)+vaf(i,j))
        el(i,j)=el(i,j)+0.5d0*smoth*(elb(i,j)-2.d0*el(i,j)+elf(i,j))
        elb(i,j)=el(i,j)
        el(i,j)=elf(i,j)
        d(i,j)=h(i,j)+el(i,j)
        uab(i,j)=ua(i,j)
        ua(i,j)=uaf(i,j)
        vab(i,j)=va(i,j)
        va(i,j)=vaf(i,j)
     end do
  end do
  !$omp end do

  if(iext /= isplit)then
     !$omp do private(i,j)     
     do j=1,jm
        do i=1,im
           egf(i,j)=egf(i,j)+el(i,j)*ispi
        end do
     end do
     !$omp end do

     !$omp do private(i,j)          
     do j=1,jm
        do i=2,im
           utf(i,j)=utf(i,j)+ua(i,j)*(d(i,j)+d(i-1,j))*isp2i
        end do
     end do
     !$omp end do

     !$omp do private(i,j)          
     do j=2,jm
        do i=1,im
           vtf(i,j)=vtf(i,j)+va(i,j)*(d(i,j)+d(i,j-1))*isp2i
        end do
     end do
     !$omp end do
  end if
  !$omp end parallel
  
  call exchange2d_mpi(utf,im,jm) !2018.08
  call exchange2d_mpi(vtf,im,jm) !2018.08
  
end subroutine mode_external

!_______________________________________________________________________
subroutine mode_internal
  !     calculate the internal (3-D) mode

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k

  if((iint /= 1 .or. time0 /= 0.d0) .and. mode /= 2) then

     !     adjust u(z) and v(z) such that depth average of (u,v) = (ua,va)
     tps(1:im,1:jm)=0.d0

     do k=1,kbm1
        do j=1,jm
           do i=1,im
              tps(i,j)=tps(i,j)+u(i,j,k)*dz(i,j,k)
           end do
        end do
     end do

     !$omp parallel
     !$omp do private(i,j,k)
     do k=1,kbm1
        do j=1,jm
           do i=2,im
              u(i,j,k)=(u(i,j,k)-tps(i,j))+(utb(i,j)+utf(i,j))/(dt(i,j)+dt(i-1,j))
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel
     
     tps(1:im,1:jm)=0.d0

     do k=1,kbm1
        do j=1,jm
           do i=1,im
              tps(i,j)=tps(i,j)+v(i,j,k)*dz(i,j,k)
           end do
        end do
     end do

     !$omp parallel
     !$omp do private(i,j,k)     
     do k=1,kbm1
        do j=2,jm
           do i=1,im
              v(i,j,k)=(v(i,j,k)-tps(i,j))+(vtb(i,j)+vtf(i,j))/(dt(i,j)+dt(i,j-1))
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel
     
     call exchange3d_mpi(u(:,:,1:kbm1),im,jm,kbm1) ! 2018.07.27
     call exchange3d_mpi(v(:,:,1:kbm1),im,jm,kbm1) ! 2018.07.27

     !     calculate w from u, v, dt (h+et), etf and etb
     call vertvl

     !     if(ltide) then
     !        call tide_bcond(5)
     !     else

     call bcond(5)

     !end if

     call exchange3d_mpi(w,im,jm,kb)
     
     !     set uf and vf to zero
     uf(1:im,1:jm,1:kb)=0.d0
     vf(1:im,1:jm,1:kb)=0.d0         
          
     !     calculate q2f and q2lf using uf, vf, a and c as temporary variables
     call advq(q2b,q2,uf)

     if(lmynnf) then
        call profq_mynnf
     else
        call advq(q2lb,q2l,vf)
        call profq
     end if
     
     !if(ltide) then
     !   call tide_bcond(6)
     !else

     call bcond(6)
     
     !end if

     call exchange3d_mpi(uf(:,:,2:kbm1),im,jm,kbm2)
     if(.not. lmynnf) call exchange3d_mpi(vf(:,:,2:kbm1),im,jm,kbm2)
     
     if(lmynnf) then
 
        !$omp parallel
        !$omp do private(i,j,k)
        do k=1,kb
           do j=1,jm
              do i=1,im
                 q2(i,j,k)=q2(i,j,k) &
                      +0.5d0*smoth*(uf(i,j,k)+q2b(i,j,k)-2.d0*q2(i,j,k))
                 q2b(i,j,k)=q2(i,j,k)
                 q2(i,j,k)=uf(i,j,k)
              end do
           end do
        end do
        !$omp end do
        !$omp end parallel

     else

        !$omp parallel
        !$omp do private(i,j,k)        
        do k=1,kb
           do j=1,jm
              do i=1,im
                 q2(i,j,k)=q2(i,j,k) &
                       & +0.5d0*smoth*(uf(i,j,k)+q2b(i,j,k)-2.d0*q2(i,j,k))
                 q2l(i,j,k)=q2l(i,j,k) &
                      & +0.5d0*smoth*(vf(i,j,k)+q2lb(i,j,k)-2.d0*q2l(i,j,k))
                 q2b(i,j,k)=q2(i,j,k)
                 q2(i,j,k)=uf(i,j,k)
                 q2lb(i,j,k)=q2l(i,j,k)
                 q2l(i,j,k)=vf(i,j,k)
              end do
           end do
        end do
        !$omp end do
        !$omp end parallel
        
     end if
     
     !     calculate tf and sf using uf, vf, a and c as temporary variables
     if(mode /= 4)then
        
        !     Temperature advection
        if(nadv == 1)then
           call advt1(tb,t,tclim,uf)
        else if(nadv == 2)then
           call advt2(tb,t,tclim,uf)
        else
           error_status=1
           write(6,'(/''Error: invalid value for nadv'')')
        end if
        
        if(budget == 1)then
           !$omp parallel
           !$omp do private(i,j,k)
           do k=1,kb
              do j=1,jm
                 do i=1,im
                    txadv(i,j,k)=xadvterm(i,j,k)*86400.d0 !*86400.:[deg C/sec] -> [deg C/day]
                    tyadv(i,j,k)=yadvterm(i,j,k)*86400.d0
                    tzadv(i,j,k)=zadvterm(i,j,k)*86400.d0
                    txdif(i,j,k)=xdifterm(i,j,k)*86400.d0
                    tydif(i,j,k)=ydifterm(i,j,k)*86400.d0
                 end do
              end do
           end do
           !$omp end do
           !$omp end parallel
        end if
        
        !     Salinity advection
        if(nadv == 1)then
           call advt1(sb,s,sclim,vf)
        else if(nadv == 2)then
           call advt2(sb,s,sclim,vf)
        else
           error_status=1
           write(6,'(/''Error: invalid value for nadv'')')
        end if
        
        if(budget == 1)then
           !$omp parallel
           !$omp do private(i,j,k)
           do k=1,kb
              do j=1,jm
                 do i=1,im
                    sxadv(i,j,k)=xadvterm(i,j,k)*86400.d0
                    syadv(i,j,k)=yadvterm(i,j,k)*86400.d0
                    szadv(i,j,k)=zadvterm(i,j,k)*86400.d0
                    sxdif(i,j,k)=xdifterm(i,j,k)*86400.d0
                    sydif(i,j,k)=ydifterm(i,j,k)*86400.d0
                 end do
              end do
           end do
           !$omp end do
           !$omp end parallel
        end if
            
        !     Temperature diffusion
        call proft(uf,wtsurf,tsurf,nbct)
           
        if(budget == 1)then
           !$omp parallel
           !$omp do private(i,j,k)
           do k=1,kb
              do j=1,jm
                 do i=1,im
                    tzdif(i,j,k)=zdifterm(i,j,k)*86400.d0
                    qz(i,j,k)=radterm(i,j,k)*86400.d0
                 end do
              end do
           end do
           !$omp end do

           !$omp do private(i,j)           
           do j=1,jm
              do i=1,im
                 tsfc(i,j)=sfcterm(i,j)*86400.d0
              end do
           end do
           !$omp end do
           !$omp end parallel

        end if
        
        !     Salinity diffusion
        call proft(vf,wssurf,ssurf,nbcs)

        if(budget == 1)then

           !$omp parallel
           !$omp do private(i,j,k)           
           do k=1,kb
              do j=1,jm
                 do i=1,im
                    szdif(i,j,k)=zdifterm(i,j,k)*86400.d0
                 end do
              end do
           end do
           !$omp end do
           !$omp do private(i,j)           
           do j=1,jm
              do i=1,im
                 ssfc(i,j)=sfcterm(i,j)*86400.d0
              end do
           end do
           !$omp end do
           !$omp end parallel

        end if
                
        !     S.Ohishi(2018.12)
        if(assim == 2)then
           call ts_iau
        end if
        
        !     S.Ohishi(2020.04)
        call ts_nudging
        
        !     S.Ohishi (2020.02)
        if(lroff)then
           call ts_roff
        end if

        !if(ltide) then
        !   call tide_bcond(4)
        !else
        
        call bcond(4)
        
        !end if
            
        !if(lfrs) call bcond_frs(4)
        
        call exchange3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
        call exchange3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)

        !     S.Ohishi(2019.04)
        if(budget == 1)then

           !$omp parallel
           !$omp do private(i,j,k)
           !     k=1,kbm1
           do k=1,kbm1
              do j=1,jm
                 do i=1,im
                    dtdt(i,j,k)=(uf(i,j,k)-tb(i,j,k))/dti2*86400.d0
                    dsdt(i,j,k)=(vf(i,j,k)-sb(i,j,k))/dti2*86400.d0
                 end do
              end do
           end do
           !$omp end do
           !$omp end parallel

           !     k=kb
           dtdt(1:im,1:jm,kb)=0.d0
           dsdt(1:im,1:jm,kb)=0.d0
           
        end if
            
        !     Asselin filter (Asselin 1972)
        !$omp parallel
        !$omp do private(i,j,k)
        do k=1,kb
           do j=1,jm
              do i=1,im
                 t(i,j,k)=t(i,j,k) &
                      & +0.5d0*smoth*(uf(i,j,k)+tb(i,j,k)-2.d0*t(i,j,k))
                 s(i,j,k)=s(i,j,k) &
                      & +0.5d0*smoth*(vf(i,j,k)+sb(i,j,k)-2.d0*s(i,j,k))
                 tb(i,j,k)=t(i,j,k)
                 t(i,j,k)=uf(i,j,k)
                 sb(i,j,k)=s(i,j,k)
                 s(i,j,k)=vf(i,j,k)
              end do
           end do
        end do
        !$omp end do
        !$omp end parallel
        
        !     Bottom value S.Ohishi 2019.04
        t(1:im,1:jm,kb)=0.d0
        tb(1:im,1:jm,kb)=0.d0
        s(1:im,1:jm,kb)=0.d0
        sb(1:im,1:jm,kb)=0.d0

        call dens(s,t,rho)
        
     end if
          
     !     calculate uf and vf
     call advu
     call advv
     
     call profu
     call profv

     !     S.Ohishi(2019.12)
     if(assim == 2)then
        call uv_iau
     end if
     
     !if(ltide) then
     !   call tide_bcond(3)
     !else

     call bcond(3)
     
     !end if
     
     call exchange3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
     call exchange3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)
     
     !after caluculate uf and vf
     tps(1:im,1:jm)=0.d0

     do k=1,kbm1
        do j=1,jm
           do i=1,im
              tps(i,j)=tps(i,j) &
                   & +(uf(i,j,k)+ub(i,j,k)-2.d0*u(i,j,k))*dz(i,j,k)
           end do
        end do
     end do

     !$omp parallel          
     !$omp do private(i,j,k)     
     do k=1,kbm1
        do j=1,jm
           do i=1,im
              u(i,j,k)=u(i,j,k) &
                   & +0.5d0*smoth*(uf(i,j,k)+ub(i,j,k)-2.d0*u(i,j,k)-tps(i,j))
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel

     tps(1:im,1:jm)=0.d0

     do k=1,kbm1
        do j=1,jm
           do i=1,im
              tps(i,j)=tps(i,j) &
                   & +(vf(i,j,k)+vb(i,j,k)-2.d0*v(i,j,k))*dz(i,j,k)
           end do
        end do
     end do

     !$omp parallel          
     !$omp do private(i,j,k)     
     do k=1,kbm1
        do j=1,jm
           do i=1,im
              v(i,j,k)=v(i,j,k) &
                   & +0.5d0*smoth*(vf(i,j,k)+vb(i,j,k)-2.d0*v(i,j,k)-tps(i,j))
           end do
        end do
     end do
     !$omp end do

     !$omp do private(i,j,k)          
     do k=1,kb
        do j=1,jm
           do i=1,im
              ub(i,j,k)=u(i,j,k)
              u(i,j,k)=uf(i,j,k)
              vb(i,j,k)=v(i,j,k)
              v(i,j,k)=vf(i,j,k)
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end if

  egb(1:im,1:jm)=egf(1:im,1:jm)
  etb(1:im,1:jm)=et(1:im,1:jm)
  et(1:im,1:jm)=etf(1:im,1:jm)
  utb(1:im,1:jm)=utf(1:im,1:jm)
  vtb(1:im,1:jm)=vtf(1:im,1:jm)
  vfluxb(1:im,1:jm)=vfluxf(1:im,1:jm)

  
  !$omp parallel
  !$omp do private(i,j)  
  do j=1,jm
     do i=1,im
        dt(i,j)=h(i,j)+et(i,j)
     end do
  end do
  !$omp end do
  !$omp end parallel
  
  call realvertvl
    
end subroutine mode_internal

!_______________________________________________________________________
subroutine print_section
  ! print output

  use common_pom_var
  implicit none

  integer i,j,k
  real(kind = r_size) atot,darea,dvol,eaver,saver,taver,vtot,tsalt

  if(mod(iint,iprint) == 0)then

     ! print time
     if(my_task == master_task) &
          & write(6,'(/ &
          & ''**********************************************************'' &
          & /''time ='',f12.4,'', iint ='',i8,'', iext ='',i8, &
          & '', iprint ='',i8)') time,iint,iext,iprint

     ! check for errors
     call sum0d_mpi(error_status,master_task)
     call bcast0d_mpi(error_status,master_task)
     if(error_status /= 0)then
        if(my_task == master_task) &
             & write(5,'(/a,i10)') 'POM terminated with error: print_section'
        call finalize_mpi
        stop
     end if

     ! local averages
     vtot=0.d0
     atot=0.d0
     taver=0.d0
     saver=0.d0
     eaver=0.d0
     do k=1,kbm1
        do j=1,jm
           do i=1,im
              darea=dx(i,j)*dy(i,j)*fsm(i,j)
              dvol=darea*dt(i,j)*dz(i,j,k)
              vtot=vtot+dvol
              taver=taver+tb(i,j,k)*dvol
              saver=saver+sb(i,j,k)*dvol
           end do
        end do
     end do

     do j=1,jm
        do i=1,im
           darea=dx(i,j)*dy(i,j)*fsm(i,j)
           atot=atot+darea
           eaver=eaver+et(i,j)*darea
        end do
     end do

     if(vtot /= vtot)then
        write(5,"(/a)") "POM terminated with error: vtot"
        call finalize_mpi
        stop
     else if(vtot == 0.d0)then
        taver=0.d0
        saver=0.d0
     else
        taver=taver/vtot
        saver=saver/vtot
     end if

     if(atot == 0.d0)then
        eaver=0.d0
     else
        eaver=eaver/atot
     end if

     tsalt=(saver+sbias)*vtot

     ! print averages
     ! global averages requiere to transfer high amounts of data between
     ! processor - therefore, only local average for master_task is printed
     if(my_task == master_task) &
          & write(6,'(/''vtot = '',e16.7, &
          & ''   atot = '',e16.7,''  eaver ='',e16.7/''taver ='',e16.7, &
          & ''   saver ='',e16.7,''  tsalt ='',e16.7)') &
          & vtot,atot,eaver,taver,saver,tsalt

  end if

end subroutine print_section

!_______________________________________________________________________
subroutine ts_nudging
  !     Temperature & Salinity nudges toward tsclim monthly climatology
  !     Created by S.Ohishi 2020.04
  !     Modified by S.Ohishi 2023.04

  !$use omp_lib  
  use common_pom_var
  implicit none

  real(kind = r_size),parameter :: mld_crit=0.125d0
  
  integer i,j,k
      
  real(kind = r_size) time_scale
  real(kind = r_size) pd(im,jm,kb)
  
  tnudge(:,:,:)=0.d0
  snudge(:,:,:)=0.d0

  if(ts_nudge == 0.d0 .and. ti_nudge == 0.d0 &
       & .and. ss_nudge == 0.d0 .and. si_nudge == 0.d0)then
     return
  else
     call unesco_potential_density(uf,vf,pd)
  end if
  
  !     Temperature
  !$omp parallel
  !$omp do private(i,j,k,time_scale)
  do k=1,kb-1
     do j=1,jm
        do i=1,im

           !Set timescale
           if(abs(pd(i,j,1)-pd(i,j,k)) < mld_crit)then
              time_scale=ts_nudge
           else
              time_scale=ti_nudge
           end if

           !T nudging
           if(time_scale > 0.d0)then
              time_scale=dti2/(time_scale*86400.d0) ![s/step]/[day*s/day] = [/step]           
              tnudge(i,j,k)=(tclimm(i,j,k)-uf(i,j,k))*time_scale ![degree C/step]
              uf(i,j,k)=uf(i,j,k)+tnudge(i,j,k)
              tnudge(i,j,k)=tnudge(i,j,k)/dti2*86400.d0 ![degree C/step] --> [dgree C/day]
           else
              tnudge(i,j,k)=0.d0
           end if
           
        end do
     end do
  end do
  !$omp end do

  !     Salinity
  !$omp do private(i,j,k,time_scale)
  do k=1,kb-1
     do j=1,jm
        do i=1,im

           !Set timescale
           if(abs(pd(i,j,1)-pd(i,j,k)) < mld_crit)then
              time_scale=ss_nudge
           else
              time_scale=si_nudge
           end if

           !S nudging
           if(time_scale > 0.d0)then
              time_scale=dti2/(time_scale*86400.d0) ![s/step]/[day*s/day] = [/step]           
              snudge(i,j,k)=(sclimm(i,j,k)-vf(i,j,k))*time_scale ![psu/step]
              vf(i,j,k)=vf(i,j,k)+snudge(i,j,k)
              snudge(i,j,k)=snudge(i,j,k)/dti2*86400.d0 ![psu/step] --> [psu/day]
           else
              snudge(i,j,k)=0.d0
           end if
           
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine ts_nudging

!_______________________________________________________________________
subroutine ts_roff
  !     created by S.Ohishi 2020.02 

  !$use omp_lib  
  use common_pom_var
  implicit none
      
  integer i,j,k
  
  !$omp parallel
  !$omp do private(i,j,k)  
  do k=1,kb-1
     do j=1,jm
        do i=1,im
           
           if(uf(i,j,k) == 0.d0) cycle
           
           if(uf(i,j,k) < tmin)then
              troff(i,j,k)=uf(i,j,k)-tmin
              uf(i,j,k)=tmin
           else if(tmax < uf(i,j,k))then
              troff(i,j,k)=tmax-uf(i,j,k)
              uf(i,j,k)=tmax
           end if
           
        end do
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)  
  do k=1,kb-1
     do j=1,jm
        do i=1,im
           
           if(vf(i,j,k) == 0.d0) cycle

           if(vf(i,j,k) < smin)then
              sroff(i,j,k)=vf(i,j,k)-smin
              vf(i,j,k)=smin
           else if(smax < vf(i,j,k))then
              sroff(i,j,k)=smax-vf(i,j,k)
              vf(i,j,k)=smax     
           endif
           
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine ts_roff
!_______________________________________________________________________
subroutine check_velocity

  ! check if velocity condition is violated

  use common_pom_var
  implicit none

  integer i,j
  integer imax,jmax
  
  real(kind = r_size) vamax

  vamax=0.d0
  
  do j=1,jm
     do i=1,im
        if(abs(vaf(i,j)) >= vamax)then
           vamax=abs(vaf(i,j))
           imax=i
           jmax=j
        end if
     end do
  end do

  if(vamax > vmaxl) then
     
     if(error_status == 0)then
        write(6,*) my_task,t(imax,jmax,1),s(imax,jmax,1),u(imax,jmax,1),v(imax,jmax,1)
        write(6,'(/ &
          & ''Error: velocity condition violated''/''time ='',f9.4, &
          & '', iint ='',i8,'', iext ='',i8,'', iprint ='',i8,/ &
          & ''vamax ='',e12.3,''   imax,jmax ='',2i5)') &
          & time,iint,iext,iprint,vamax,imax,jmax
        write(6,*) "el:",elb(imax,jmax),el(imax,jmax),elf(imax,jmax)
        write(6,*) "ua:",uab(imax,jmax),ua(imax,jmax),uaf(imax,jmax)
        write(6,*) "va:",vab(imax,jmax),va(imax,jmax),vaf(imax,jmax)
        write(6,*) "t:",tb(imax,jmax,1),t(imax,jmax,1)
        write(6,*) "s:",sb(imax,jmax,1),s(imax,jmax,1)
        write(6,*) "u:",ub(imax,jmax,1),u(imax,jmax,1)
        write(6,*) "v:",vb(imax,jmax,1),v(imax,jmax,1)
     end if
        
     error_status=1
     call finalize_mpi
     stop

  end if

end subroutine check_velocity

!_____________________________________________________________________________

subroutine ts_iau

  !     IAU method 
  !     Relaxation of TS toward analysis TS

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k

  !$omp parallel
  !$omp do private(i,j,k)           
  do k=1,kb
     do j=1,jm
        do i=1,im
           uf(i,j,k)=uf(i,j,k)+2.d0*t_iau(i,j,k)/dble(iend)
           vf(i,j,k)=vf(i,j,k)+2.d0*s_iau(i,j,k)/dble(iend)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine ts_iau
      
!________________________________________________________________________________
subroutine uv_iau

  !     IAU method 
  !     Relaxation of TS toward analysis UV

  !$use omp_lib  
  use common_pom_var
  implicit none

  integer i,j,k

  !$omp parallel
  !$omp do private(i,j,k)  
  do k=1,kb
     do j=1,jm
        do i=1,im
           uf(i,j,k)=uf(i,j,k)+2.d0*u_iau(i,j,k)/dble(iend)
           vf(i,j,k)=vf(i,j,k)+2.d0*v_iau(i,j,k)/dble(iend)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine uv_iau

!_______________________________________________________________________
subroutine check_nan(varname,im,jm,km,dat)
  
  use common_pom_var, only: iint, r_size, my_task
  use mpi
  use,intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none

  !---Common
  integer i,j,k
  integer local_err(1),global_err(1)
  integer ierr

  !---IN
  integer,intent(in) :: im,jm,km

  real(kind = r_size),intent(in) :: dat(im,jm,km)

  character(*),intent(in) :: varname

  local_err(1)=0

  do k=1,km
    do j=1,jm
      do i=1,im
         if(ieee_is_nan(dat(i,j,k)))then
            write(6,*) 'NaN:', trim(varname), ' step=',iint, &
                 ' i,j,k=',i,j,k, ' rank=',my_task
            local_err(1)=1
         end if
      end do
    end do
 end do
 
 call MPI_ALLREDUCE(local_err,global_err,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)

 if(global_err(1) == 1)then
    call finalize_mpi
    stop 'NaN detected'
 end if

end subroutine check_nan
