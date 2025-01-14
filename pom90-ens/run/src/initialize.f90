! initialize.f90

! initialize POM:
! define constant, read initial values, grid and initial conditions

!_______________________________________________________________________
subroutine initialize
  
  ! initialize POM
  use common_pom_var
  implicit none
  
  integer iyr,imon,iday
  integer julday

  character(4) yyyy
  character(2) mm,dd
    
  ! initialize the MPI execution environment and create communicator for
  ! internal POM communications
  call initialize_mpi
  
  ! distribute the model domain across processors
  call distribute_mpi

  ! read input values and define constants
  call read_input

  ! initialize arrays for safety (may be overwritten)
  call initialize_arrays
  
  ! read in grid data
  call read_grid
  
  ! read initial and lateral boundary conditions
  call initial_conditions
  
  ! update initial conditions and set the remaining conditions
  call update_initial
  
  ! calculate the bottom friction coefficient
  call bottom_friction

  ! read restart data from a previous run
  if(nread_rst /= 0) call read_restart_netcdf
  if(nread_rst == 2) call restart_tsuv

  ! read iau data from letkf analysis (S.Ohishi 2018.12)
  if(assim == 2) call read_iau_netcdf
      
  ! set start julian day
  yyyy=time_start(1:4)
  read(yyyy,*) iyr
  mm=time_start(6:7)
  read(mm,*) imon
  dd=time_start(9:10)
  read(dd,*) iday

  julday_start=julday(iyr,imon,iday)

  ! check for errors
  call sum0d_mpi(error_status,master_task)
  call bcast0d_mpi(error_status,master_task)
  if(my_task == master_task) write(6,'(a)') 'End of initialization'

end subroutine initialize

!_______________________________________________________________________
subroutine read_input
  ! read input values and defines constants

  use common_pom_var
  implicit none
  
  namelist/pom_nml/ &
       & title,netcdf_file,mode,assim, &
       & nadv,nitera,sw,npg,dte, &
       & isplit,iens,nens,time_start, &
       & nread_rst,read_rst_file,read_iau_file, &
       & write_rst,write_rst_file,write_iau_file, &
       & budget,days,prtd1,prtd2, &
       & swtch, &
       & ts_nudge,ti_nudge,ss_nudge,si_nudge

  ! read input namelist
  open(idevinout,file='pom.nml',status='old')
  read(idevinout,nml=pom_nml)
  close(idevinout)
  
  ! Input of filenames and constants

  ! Logical for inertial ramp (.true. if inertial ramp to be applied
  ! to wind stress and baroclinic forcing, otherwise .false.)
  lramp=.false.

  ! Reference density (recommended values: 1025 for seawater,
  ! 1000 for freswater; S.I. units):
  rhoref=1025.d0

  ! Temperature bias (deg. C)
  tbias=0.d0

  ! Salinity bias
  sbias=0.d0

  ! gravity constant (S.I. units)
  grav=9.806d0

  ! von Karman's constant
  kappa=0.4d0

  ! Bottom roughness (metres)
  z0b=0.01d0

  ! Minimum bottom friction coeff.
  cbcmin=0.0025d0

  ! Maximum bottom friction coeff.
  cbcmax=1.d0

  ! Smagorinsky diffusivity coeff.
  horcon=0.1d0
  !horcon=0.3d0
  aamadd=0.d0

  ! Inverse horizontal turbulent Prandtl number (ah/am; dimensionless):
  ! NOTE that tprni=0.e0 yields zero horizontal diffusivity!
  tprni=0.2d0
  !tprni=0.5d0

  ! Background viscosity used in subroutines profq, proft, profu and
  ! profv (S.I. units):
  umol=2.d-5

  ! Maximum magnitude of vaf (used in check that essentially tests
  ! for CFL violation):
  vmaxl=100.d0

  ! Maximum allowable value of:
  !   <difference of depths>/<sum of depths>
  ! for two adjacent cells (dimensionless). This is used in subroutine
  ! slpmax. If >= 1, then slpmax is not applied:
  slmax=2.d0

  ! Water type, used in subroutine proft.
  !    ntp    Jerlov water type
  !     1            i
  !     2            ia
  !     3            ib
  !     4            ii
  !     5            iii
  ntp=2

  ! Surface temperature boundary condition, used in subroutine proft:
  !    nbct   prescribed    prescribed   short wave
  !           temperature      flux      penetration
  !     1        no           yes           no
  !     2        no           yes           yes
  !     3        yes          no            no
  !     4        yes          no            yes
  nbct=2

  ! Surface salinity boundary condition, used in subroutine proft:
  !    nbcs   prescribed    prescribed
  !            salinity      flux
  !     1        no           yes
  !     3        yes          no
  ! NOTE that only 1 and 3 are allowed for salinity.
  nbcs=1
  
  !     2018.06 added by S.Ohishi
  !     2020.04 modified by S.Ohishi
  ! Surface salinity forcing:
  !     issf  prescribed
  !             E-P-R
  !       1      yes
  !       2      no
  issf=1

  ! 2018.08 added by S.Ohishi
  ! lateral boudary forcing:
  !     ilbc   monthly climatology     daily  
  !      1           yes                no
  !      2           no                 yes
  ilbc=1

  !2018.08 added by S.Ohishi
  ! output:
  !    idave   daily average      instanteneous
  !      1         yes                 no
  !      2         no                  yes
  idave=1

  ! for wind stress formula by Mellor and Blumberg (2004)
  !      lmbws=.true.

  !2018.09 S.Ohishi
  ! turbulent heat flux & wind stress:
  ! #2-5: use https://github.com/brodeau/aerobulk, Brodeau et al. (2017) 
  !  ithf_ws
  !   1   Kondo(1975) & Mellor and Blumberg (2004)
  !   2   coare 3.0
  !   3   coare 3.5
  !   4   ncar
  !   5   ecmwf
  ithf_ws=3
      
  ! Step interval during which external (2-D) mode advective terms are
  ! not updated (dimensionless):
  ispadv=5

  ! Constant in temporal filter used to prevent solution splitting
  ! (dimensionless):
  smoth=0.10d0

  ! Weight used for surface slope term in external (2-D) dynamic
  ! equation (a value of alpha = 0.e0 is perfectly acceptable, but the
  ! value, alpha=.225e0 permits a longer time step):
  alpha=0.225d0

  ! Initial value of aam:
  aam_init=500.d0

  ! Flow Relaxation Scheme
  lfrs=.false.

  ! for surface heat flux & wind stress relative to ocean current
  lrtvf=.true.

  ! for MYNNF
  lmynnf=.true.

  ! for tide
  ltide=.false.

  ! Surface wave breaking (Mellor and Blumberg, 2004)
  lwbrk=.false.
  ! for internal wave breaking parameterization
  liwbrk=.false.

  ! S.Ohishi 2020.04
  !for temperature/salinity round off
  lroff=.true.

  ! S.Ohishi 2023.12
  !Parallel netcdf
  lpnetcdf=.true.

  ! S.Ohishi 2024.12
  alpha_atm=0.2d0
  
  ! End of input of constants

  ! calculate some constants
  small=1.d-9           ! Small value
  pi=atan(1.d0)*4.d0    ! PI
  earth=6371.d3         ! Earth radius
  
  dti=dte*float(isplit)
  dte2=dte*2
  dti2=dti*2

  iend=max0(nint(days*24.d0*3600.d0/dti),2)
  iprint=nint(prtd1*24.d0*3600.d0/dti)
  iswtch=nint(swtch*24.d0*3600.d0/dti)
  irestart=nint(write_rst*24.d0*3600.d0/dti)

  ispi=1.d0/float(isplit)
  isp2i=1.d0/(2.d0*float(isplit))

! initialise time
  time0=0.d0
  time=0.d0

  ! print initial summary
  if(my_task == master_task)then
     write(6,'(/'' title      = '',a40)') title
     write(6,'(/'' mode       = '',i10)') mode
     write(6,'(/'' assim       = '',i10)') assim
     write(6,'('' nadv       = '',i10)') nadv
     write(6,'('' nitera     = '',i10)') nitera
     write(6,'('' sw         = '',f10.4)') sw
     write(6,'('' npg         = '',i10)') npg
     write(6,'('' nread_rst  = '',i10)') nread_rst
     write(6,'('' write_rst  = '',f10.4)') write_rst
     write(6,'('' irestart   = '',i10)') irestart
     write(6,'('' dte        = '',f10.2)') dte
     write(6,'('' dti        = '',f10.1)') dti
     write(6,'('' isplit     = '',i10)') isplit
     write(6,'('' time_start = '',a26)') time_start
     write(6,'('' days       = '',f10.4)') days
     write(6,'('' iend       = '',i10)') iend
     write(6,'('' prtd1      = '',f10.4)') prtd1
     write(6,'('' iprint     = '',i10)') iprint
     write(6,'('' prtd2      = '',f10.4)') prtd2
     write(6,'('' swtch      = '',f10.2)') swtch
     write(6,'('' iswtch     = '',i10)') iswtch
     write(6,'('' lramp      = '',l10)') lramp
     write(6,'('' lfrs       = '',l10)') lfrs
     write(6,'('' lrtvf      = '',l10)') lrtvf
     write(6,'('' lmynnf     = '',l10)') lmynnf
     write(6,'('' rhoref     = '',f10.3)') rhoref
     write(6,'('' tbias      = '',f10.3)') tbias
     write(6,'('' sbias      = '',f10.3)') sbias
     write(6,'('' grav       = '',f10.4)') grav
     write(6,'('' kappa      = '',f10.4)') kappa
     write(6,'('' z0b        = '',f10.6)') z0b
     write(6,'('' cbcmin     = '',f10.6)') cbcmin
     write(6,'('' cbcmax     = '',f10.6)') cbcmax
     write(6,'('' horcon     = '',f10.3)') horcon
     write(6,'('' aamadd     = '',f10.3)') aamadd
     write(6,'('' tprni      = '',f10.4)') tprni
     write(6,'('' umol       = '',f10.4)') umol
     write(6,'('' vmaxl      = '',f10.4)') vmaxl
     write(6,'('' slmax      = '',f10.4)') slmax
     write(6,'('' ntp        = '',i10)') ntp
     write(6,'('' nbct       = '',i10)') nbct
     write(6,'('' nbcs       = '',i10)') nbcs
     write(6,'('' budget     = '',i10)') budget       
     write(6,'('' idave      = '',i10)') idave       
     write(6,'('' ilbc       = '',i10)') ilbc       
     write(6,'('' issf       = '',i10)') issf        
     write(6,'('' ithf_ws    = '',i10)') ithf_ws
     write(6,'('' ispadv     = '',i10)') ispadv
     write(6,'('' smoth      = '',f10.4)') smoth
     write(6,'('' alpha      = '',f10.4)') alpha
     write(6,'('' ts_nudge   = '',f10.4)') ts_nudge
     write(6,'('' ti_nudge   = '',f10.4)') ti_nudge
     write(6,'('' ss_nudge   = '',f10.4)') ss_nudge
     write(6,'('' si_nudge   = '',f10.4)') si_nudge
     write(6,'('' iens       = '',i10)') iens
     write(6,'('' alpha_atm  = '',f10.4)') alpha_atm
     write(6,'('' ltide      = '',l10)') ltide
     write(6,'('' lwbrk      = '',l10)') lwbrk
     write(6,'('' liwbrk     = '',l10)') liwbrk
     write(6,'('' lroff      = '',l10)') lroff
     write(6,'('' lpnetcdf      = '',l10)') lpnetcdf
     write(6,'('' lqglobal      = '',l10)') lqglobal
  end if
  
end subroutine read_input

!_______________________________________________________________________
subroutine initialize_arrays
  
  ! initialize arrays for safety

  use common_pom_var
  implicit none
  
  ! boundary arrays
  vabn(1:im)=0.d0
  vabs(1:im)=0.d0
  eln(1:im)=0.d0
  els(1:im)=0.d0
  
  vbn(1:im,1:kb)=0.d0
  vbs(1:im,1:kb)=0.d0
  tbn(1:im,1:kb)=0.d0
  tbs(1:im,1:kb)=0.d0
  sbn(1:im,1:kb)=0.d0
  sbs(1:im,1:kb)=0.d0

  uabe(1:jm)=0.d0
  uabw(1:jm)=0.d0
  ele(1:jm)=0.d0
  elw(1:jm)=0.d0

  ube(1:jm,1:kb)=0.d0
  ubw(1:jm,1:kb)=0.d0
  tbe(1:jm,1:kb)=0.d0
  tbw(1:jm,1:kb)=0.d0
  sbe(1:jm,1:kb)=0.d0
  sbw(1:jm,1:kb)=0.d0
        
  ! 2-D and 3-D arrays
  uab(1:im,1:jm)=0.d0
  vab(1:im,1:jm)=0.d0
  elb(1:im,1:jm)=0.d0
  etb(1:im,1:jm)=0.d0
  e_atmos(1:im,1:jm)=0.d0
  vfluxb(1:im,1:jm)=0.d0
  vfluxf(1:im,1:jm)=0.d0
  wusurf(1:im,1:jm)=0.d0
  wvsurf(1:im,1:jm)=0.d0
  tsurf(1:im,1:jm)=0.d0
  ssurf(1:im,1:jm)=0.d0
  wtsurf(1:im,1:jm)=0.d0
  wssurf(1:im,1:jm)=0.d0
  swrad(1:im,1:jm)=0.d0
  drx2d(1:im,1:jm)=0.d0
  dry2d(1:im,1:jm)=0.d0

  ub(1:im,1:jm,1:kb)=0.d0
  vb(1:im,1:jm,1:kb)=0.d0

  ua_iau(1:im,1:jm)=0.d0
  va_iau(1:im,1:jm)=0.d0
  el_iau(1:im,1:jm)=0.d0

  t_iau(1:im,1:jm,1:kb)=0.d0
  s_iau(1:im,1:jm,1:kb)=0.d0
  u_iau(1:im,1:jm,1:kb)=0.d0
  v_iau(1:im,1:jm,1:kb)=0.d0

  tnudge_dave(1:im,1:jm,1:kb)=0.d0
  snudge_dave(1:im,1:jm,1:kb)=0.d0
  
end subroutine initialize_arrays

!_______________________________________________________________________
subroutine read_grid
  
  ! set up vertical and horizontal grid, topography, areas and masks

  !$use omp_lib    
  use common_pom_var
  implicit none
  
  integer i,j,k
  integer imax,jmax
  real(kind = r_size) deg2rad

  ! degrees to radians
  deg2rad=pi/180.d0

  ! read grid
  call read_grid_netcdf

  ! derived vertical grid variables

  !$omp parallel
  !$omp do private(i,j,k)
  do k=1,kb-1
     do j=1,jm
        do i=1,im
           dz(i,j,k)=z(i,j,k)-z(i,j,k+1)
           dzz(i,j,k)=zz(i,j,k)-zz(i,j,k+1)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  dz(:,:,kb)=0.d0
  dzz(:,:,kb)=0.d0
  
  ! print vertical grid information
  if(my_task == master_task)then
     do j=1,jm
        do i=1,im
           if(maxval(h(:,:)) == h(i,j))then
              imax=i
              jmax=j
           end if
        end do
     end do
     
     write(6,*)
     write(6,'(4x,a,4x,a,4x,a,2x,a,2x,a,7x,a,7x,a,7x,a,7x,a,7x,a,7x,a)') &
          & "i","j","k", &
          & "depth(w)","depth(t)", &
          & "  z"," zz"," dz","dzz"
     do k=1,kb
        write(6,'(3i5,2f10.2,4f10.5)') imax,jmax,k, &
             & z(imax,jmax,k)*h(imax,jmax),zz(imax,jmax,k)*h(imax,jmax), &
             & z(imax,jmax,k),zz(imax,jmax,k),dz(imax,jmax,k),dzz(imax,jmax,k)
     end do
  end if

  ! set up Coriolis parameter
  !$omp parallel
  !$omp do private(i,j)
  do j=1,jm
     do i=1,im
        cor(i,j)=2.d0*7.2921159d-5*sin(north_e(i,j)*deg2rad)
     end do
  end do
  !$omp end do
        
  ! calculate areas of "t" and "s" cells
  !$omp do private(i,j)
  do j=1,jm
     do i=1,im
        art(i,j)=dx(i,j)*dy(i,j)
     end do
  end do
  !$omp end do
  
  ! calculate areas of "u" and "v" cells
  !$omp do private(i,j)
  do j=2,jm
     do i=2,im
        aru(i,j)=0.25d0*(dx(i,j)+dx(i-1,j))*(dy(i,j)+dy(i-1,j))
        arv(i,j)=0.25d0*(dx(i,j)+dx(i,j-1))*(dy(i,j)+dy(i,j-1))
     end do
  end do
  !$omp end do
  !$omp end parallel

  ! inertial period for temporal filter
  period=(2.d0*pi)/abs(cor(im/2,jm/2))/86400.d0
  
  call exchange2d_mpi(aru,im,jm)
  call exchange2d_mpi(arv,im,jm)

  if(n_west == -1)then
     aru(1,1:jm)=aru(2,1:jm)
     arv(1,1:jm)=arv(2,1:jm)
  end if

  if(n_south == -1)then
     aru(1:im,1)=aru(1:im,2)
     arv(1:im,1)=arv(1:im,2)
  end if
  
  !$omp parallel
  !$omp do private(i,j)
  do j=1,jm
     do i=1,im
        d(i,j)=h(i,j)+el(i,j)
        dt(i,j)=h(i,j)+et(i,j)
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine read_grid

!_______________________________________________________________________
subroutine initial_conditions
  ! set up initial and lateral boundary conditions

  !$use omp_lib
  use common_pom_var
  implicit none
  
  integer i,j,k
  integer ierr
  
  real(kind = r_size) elejmid,elwjmid
  
  call read_initial_tsuv_netcdf(tb,sb,ub,vb)
  
  !$omp parallel
  !$omp do private(i,j,k)  
  do k=1,kbm1
     do j=1,jm
        do i=1,im
           tb(i,j,k)=tb(i,j,k)-tbias
           sb(i,j,k)=sb(i,j,k)-sbias
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
           
  t(1:im,1:jm,1:kbm1)=tb(1:im,1:jm,1:kbm1)
  s(1:im,1:jm,1:kbm1)=sb(1:im,1:jm,1:kbm1)
  tclim(1:im,1:jm,1:kbm1)=tb(1:im,1:jm,1:kbm1)
  sclim(1:im,1:jm,1:kbm1)=sb(1:im,1:jm,1:kbm1)
  
  ! read climatological/mean temperature and salinity 
  call read_tsclim_netcdf
  
  !$omp parallel
  !$omp do private(i,j,k)  
  do k=1,kbm1
     do j=1,jm
        do i=1,im
           tclim(i,j,k)=tclim(i,j,k)-tbias
           sclim(i,j,k)=sclim(i,j,k)-sbias
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  ! density
  call dens(sb,tb,rho)  
  
  ! lateral viscosity: add a sponge layer at the east (downstream) end
  !$omp parallel
  if(n_east == -1) then

     !$omp do private(i,j)  
     do j=1,jm
        do i=1,im
           aam2d(i,j)=100.d0*(1.d0+10.d0*exp(-(im-i)/10.d0))
        end do
     end do
     !$omp end do
     
  end if
  !$omp end parallel
  
  ! generated horizontally averaged density field (in this application,
  ! the initial condition for density is a function of z (the vertical
  ! cartesian coordinate) -- when this is not so, make sure that rmean
  ! has been area averaged before transfer to sigma coordinates)
  !      call dens(smean,tmean,rmean)
  !      do k=1,kbm1
  !        do j=1,jm
  !          do i=1,im
  !            rmean(i,j,k)=rho(i,j,k)
  !          end do
  !        end do
  !      end do
  
  ! set lateral boundary conditions, for use in subroutine bcond (in the
  ! seamount problem, the east and west boundaries are open, while the
  ! south and north boundaries are closed through the specification of the
  ! masks fsm, dum and dvm)
  rfe=1.d0
  rfw=1.d0
  rfn=1.d0
  rfs=1.d0

  !$omp parallel  
  !$omp do private(j)
  do j=2,jmm1
     ! set geostrophically conditioned elevations at the boundaries
     ele(j)=ele(j-1)-cor(imm1,j)*uabe(j)/grav*dy(imm1,j-1)
     elw(j)=elw(j-1)-cor(2,j)*uabw(j)/grav*dy(2,j-1)
  end do
  !$omp end do
  !$omp end parallel
  
  ! adjust boundary elevations so that they are zero in the middle of the
  ! channel
  elejmid=ele(jmm1/2)
  elwjmid=elw(jmm1/2)

  !$omp parallel    
  !$omp do private(j)  
  do j=2,jmm1
     ele(j)=(ele(j)-elejmid)*fsm(im,j)
     elw(j)=(elw(j)-elwjmid)*fsm(2,j)
  end do
  !$omp end do
  !$omp end parallel
  
  ! set thermodynamic boundary conditions (for the seamount problem, and
  ! other possible applications, lateral thermodynamic boundary conditions
  ! are set equal to the initial conditions and are held constant
  ! thereafter - users may create variable boundary conditions)

  tbe(1:jm,1:kbm1)=tb(im,1:jm,1:kbm1)
  tbw(1:jm,1:kbm1)=tb(1,1:jm,1:kbm1)
  sbe(1:jm,1:kbm1)=sb(im,1:jm,1:kbm1)
  sbw(1:jm,1:kbm1)=sb(1,1:jm,1:kbm1)

  tbn(1:im,1:kbm1)=tb(1:im,jm,1:kbm1)
  tbs(1:im,1:kbm1)=tb(1:im,1,1:kbm1)
  sbn(1:im,1:kbm1)=sb(1:im,jm,1:kbm1)
  sbs(1:im,1:kbm1)=sb(1:im,1,1:kbm1)
  
end subroutine initial_conditions

!_______________________________________________________________________
subroutine update_initial
  ! update the initial conditions and set the remaining initial conditions

  !$use omp_lib
  use common_pom_var
  implicit none
  
  integer i,j,k
  
  ua(1:im,1:jm)=uab(1:im,1:jm)
  va(1:im,1:jm)=vab(1:im,1:jm)
  el(1:im,1:jm)=elb(1:im,1:jm)
  et(1:im,1:jm)=etb(1:im,1:jm)
  etf(1:im,1:jm)=et(1:im,1:jm)
  w(1:im,1:jm,1)=vfluxf(1:im,1:jm)
  q2b(1:im,1:jm,1:kb)=small         
  
  !$omp parallel
  !$omp do private(i,j)  
  do j=1,jm
     do i=1,im
        d(i,j)=h(i,j)+el(i,j)
        dt(i,j)=h(i,j)+et(i,j)
     end do
  end do
  !$omp end do
  
  !$omp do private(i,j,k)
  do k=1,kb
     do j=1,jm
        do i=1,im
           l(i,j,k)=0.1d0*dt(i,j)
           q2lb(i,j,k)=l(i,j,k)*q2b(i,j,k)
           kh(i,j,k)=l(i,j,k)*sqrt(q2b(i,j,k))
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  km(1:im,1:jm,1:kb)=kh(1:im,1:jm,1:kb)
  kq(1:im,1:jm,1:kb)=kh(1:im,1:jm,1:kb)
  aam(1:im,1:jm,1:kb)=aam_init
  
  q2(1:im,1:jm,1:kbm1)=q2b(1:im,1:jm,1:kbm1)
  q2l(1:im,1:jm,1:kbm1)=q2lb(1:im,1:jm,1:kbm1)
  t(1:im,1:jm,1:kbm1)=tb(1:im,1:jm,1:kbm1)
  s(1:im,1:jm,1:kbm1)=sb(1:im,1:jm,1:kbm1)
  u(1:im,1:jm,1:kbm1)=ub(1:im,1:jm,1:kbm1)
  v(1:im,1:jm,1:kbm1)=vb(1:im,1:jm,1:kbm1)
  
  ! S.Ohishi 2018.08.20
  if(npg == 1)then
     call baropg
  else if(npg == 2)then
     call baropg_mcc
  else if(npg == 3)then
     call baropg_thiem
  else
     error_status=1
     write(6,'(/''Error: invalid value for npg'')')
  end if
  
  do k=1,kbm1
     do j=1,jm
        do i=1,im
           drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(i,j,k)
           dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(i,j,k)
        end do
     end do
  end do
  
end subroutine update_initial

!_______________________________________________________________________
subroutine restart_tsuv
  ! restart using tsuv data 

  use common_pom_var
  implicit none

  tb(1:im,1:jm,1:kbm1)=t(1:im,1:jm,1:kbm1)
  sb(1:im,1:jm,1:kbm1)=s(1:im,1:jm,1:kbm1)
  ub(1:im,1:jm,1:kbm1)=u(1:im,1:jm,1:kbm1)
  vb(1:im,1:jm,1:kbm1)=v(1:im,1:jm,1:kbm1)
  
end subroutine restart_tsuv

!_______________________________________________________________________
subroutine bottom_friction
  ! calculate the bottom friction coefficient

  !$use omp_lib
  use common_pom_var
  implicit none

  integer i,j
  
  ! calculate bottom friction
  !$omp parallel
  !$omp do private(i,j)  
  do j=1,jm
     do i=1,im
        cbc(i,j)=(kappa/log((1.d0+zz(i,j,kbm1))*h(i,j)/z0b))**2
        cbc(i,j)=max(cbcmin,cbc(i,j))
        !     if the following is invoked, then it is probable that the wrong
        !     choice of z0b or vertical spacing has been made:
        cbc(i,j)=min(cbcmax,cbc(i,j))
     end do
  end do
  !$omp end do
  !$omp end parallel
      
end subroutine bottom_friction
