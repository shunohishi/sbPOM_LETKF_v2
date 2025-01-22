module common_pom_var
  
  implicit none
  
  ! contain parameters for the model domain and for the decomposition
  ! local domain size, and common POM variables

  !_______________________________________________________________________
  !Precision setting
  
  integer,parameter :: r_sngl=kind(0.0e0)
  integer,parameter :: r_dble=kind(0.0d0)
  integer,parameter :: r_size=r_sngl
  
  !_______________________________________________________________________
  ! Grid size parameters

  ! Correct values for im_local and jm_local are found using
  !   n_proc=(im_global-2)/(im_local-2)*(jm_global-2)/(jm_local-2)
  ! Values higher than necessary will not cause the code to fail, but
  ! will allocate more memory than is necessary. Value that are too low
  ! will cause the code to exit
  
  integer,parameter :: &
       !Quasi-GLOBAL
       !&  im_global=1442 , & ! number of global grid points in x
       !&  jm_global=562  , & ! number of global grid points in y
       !&  kb=75          , & ! number of grid points in z
       !&  im_local=362   , & ! number of local grid points in x
       !&  jm_local=282   , & ! number of local grid points in y
       !&  n_proc=8          ! number of processors
       !North Pacific
       &  im_global=1498 , & ! number of global grid points in x
       &  jm_global=530  , & ! number of global grid points in y
       &  kb=75          , & ! number of grid points in z
       !&  im_local=189   , & ! number of local grid points in x
       !&  jm_local=134   ,  & ! number of local grid points in y
       !&  n_proc=32          ! number of processors
       !&  im_local=376   , & ! number of local grid points in x
       !&  jm_local=134   ,  & ! number of local grid points in y
       !&  n_proc=16          ! number of processors
       !&  im_local=376   , & ! number of local grid points in x
       !&  jm_local=178   ,  & ! number of local grid points in y
       !&  n_proc=12          ! number of processors
       &  im_local=376   , & ! number of local grid points in x
       &  jm_local=266   ,  & ! number of local grid points in y
       &  n_proc=8          ! number of processors

  ! S.Ohishi 2024.04
  !Quasi-global
  logical,parameter :: lqglobal=.false.
  !logical,parameter :: lqglobal=.true.

  integer,parameter :: idevinout=73

  !Round off parameter
  real(kind = r_size),parameter :: tmin=-1.8d0,tmax=50.d0
  real(kind = r_size),parameter :: smin=25.d0,smax=50.d0

  integer,parameter :: &
       !Ensemble year range
       & syr_atm=1979,eyr_atm=2023,      & !Start and end year
       !JRA55do
       & im_atm=620,jm_atm=320,dt_atm=3, & ![hour]
       !CaMa-Flood
       & im_riv=1440,jm_riv=720,dt_riv=3   ![hour]
  
  !_______________________________________________________________________
  ! Efective grid size

  integer,save :: &
       &  im   , & ! number of grid points used in local x domains
       &  imm1 , & ! im-1
       &  imm2 , & ! im-2
       &  jm   , & ! number of grid points used in local y domains
       &  jmm1 , & ! jm-1
       &  jmm2 , & ! jm-2
       &  kbm1 , & ! kb-1
       &  kbm2     ! kb-2
  
  ! Note that im and jm may be different between local domains
  ! im and jm are equal to or lower than im_local and jm_local, depending
  ! on the use of correct values for im_local and jm_local

  !_______________________________________________________________________
  ! Parallel variables
  
  integer,save :: &
       &  my_task            , & ! actual parallel processor ID
       &  master_task        , & ! master processor ID
       &  pom_comm           , & ! POM model MPI group communicator
       &  i_global(im_local) , & ! global i index for each point in local domain
       &  j_global(jm_local) , & ! global j index for each point in local domain
       &  n_west             , & ! western parallel processor ID
       &  n_east             , & ! eastern parallel processor ID
       &  n_south            , & ! southern parallel processor ID
       &  n_north            , & ! northern parallel processor ID
       &  nproc_x            , & ! number of processors in x
       &  nproc_y                ! number of processors in y
  
  !_______________________________________________________________________
  ! Scalars
  integer,save :: &
       &  iint           , & 
       &  iprint         , & ! interval in iint at which variables are printed
       &  mode           , & ! calculation mode
       &  assim          , & ! assimilation (0:No assimilation, 1:IAU, 2:Intermittent)
       &  ntp            , & ! water type
       &  iend           , & ! total internal mode time steps
       &  iext           , & 
       &  ispadv         , & ! step interval for updating external advective terms
       &  isplit         , & ! dti/dte
       &  iswtch         , & ! time step interval to switch from prtd1 to prtd2
       &  nadv           , & ! advection scheme
       &  nbct           , & ! surface temperature boundary condition
       &  nbcs           , & ! surface salinity boundary condition
       &  budget         , & ! conservation each term in budget equations 
       &  idave          , & ! daily average output
       &  ilbc           , & ! Lateral boudary condition
       &  issf           , & ! fresh water flux
       &  ithf_ws        , & ! turbulent heat flux & wind stress algorithm
       &  nitera         , & ! number of iterations for Smolarkiewicz scheme
       &  npg            , & ! pressure gradient scheme
       &  nread_rst      , & ! index to start from restoart file
       &  irestart       , & ! restart flag
       &  iens           , & ! Ensemble member
       &  nens           , & ! Ensemble size
       &  error_status   , & ! error flag
       &  julday_start   , & ! julian day of start time
       &  lbctime_julday , & ! start time of lateral boundary data
       &  atmtime_julday , & ! start time of atmospheric condition data
       &  rivtime_julday , & ! start time of fresh water flux condition data
       &  tsdatatime_julday    ! start time of t/s reference data

  real(kind = r_size),save :: &
       &  alpha          , & ! weight for surface slope term in external eq
       &  dte            , & ! external (2-D) time step (s)
       &  dti            , & ! internal (3-D) time step (s)
       &  dti2           , & ! 2*dti
       &  grav           , & ! gravity constant (S.I. units)
       &  kappa          , & ! von Karman's constant'
       &  pi             , & ! pi
       &  earth          , & ! Earth radius       
       &  ramp           , & ! inertial ramp
       &  rfe            , & ! flag for eastern open boundary (see bcond)
       &  rfn            , & ! flag for northern open boundary (see bcond)
       &  rfs            , & ! flag for southern open boundary (see bcond)
       &  rfw            , & ! flag for western open boundary (see bcond)
       &  rhoref         , & ! reference density
       &  sbias          , & ! salinity bias
       &  slmax          , & 
       &  small          , & ! small value
       &  tbias          , & ! temperature bias
       &  tprni          , & ! inverse horizontal turbulent Prandtl number
       &  umol           , & ! background viscosity
       &  vmaxl          , & ! max vaf used to test for model blow-up
       &  write_rst      , & 
       &  aam_init       , & ! initial value of aam
       &  cbcmax         , & ! maximum bottom friction coefficient
       &  cbcmin         , & ! minimum bottom friction coefficient
       &  days           , & ! run duration in days
       &  dte2           , & ! 2*dte
       &  horcon         , & ! smagorinsky diffusivity coefficient
       &  aamadd         , & ! Additional viscosity constant
       &  ispi           , & ! dte/dti
       &  isp2i          , & ! dte/(2*dti)
       &  period         , & ! inertial period
       &  prtd1          , & ! initial print interval (days)
       &  prtd2          , & ! final print interval (days)
       &  smoth          , & ! constant to prevent solution splitting
       &  sw             , & ! smoothing parameter for Smolarkiewicz scheme
       &  swtch          , & ! time to switch from prtd1 to prtd2
       &  z0b            , & ! bottom roughness
       &  ts_nudge       , & ! SST nuding timescale [day]
       &  ti_nudge       , & ! T nuding timescale [day]
       &  ss_nudge       , & ! SSS nuding timescale [day]
       &  si_nudge       , & ! S nuding timescale [day]
       &  alpha_atm          ! xf+alpha*xf' for ensemble perturbation

  real(kind = r_dble),save :: &
       &  time             , & ! model time (days)
       &  time0            , & ! initial time (days)
       &  atmtime_dayint   , & ! time interval of atmospheric data
       &  rivtime_dayint , & ! time interval of fresh water flux data
       &  lbctime_dayint   , & ! time interval of lateral boundary data
       &  tsdatatime_dayint    ! time interval of t/s reference data
  
  !_______________________________________________________________________
  ! 1-D arrays
  integer,save :: &
       & idx_atm(im_local)           , & !Atmosphere grid idx
       & idy_atm(jm_local)               !                idy

  real(kind = r_size),save :: &
       & lon_atm(im_atm)             , & !Atmosphere longitude
       & lat_atm(jm_atm)             , & !           latitude
       & lon_riv(im_riv)             , & !River longitude
       & lat_riv(jm_riv)             , & !      latitude
       & dlon_riv(im_riv)            , & !
       & dlat_riv(jm_riv)
         
  !_______________________________________________________________________
  ! 2-D arrays
  integer,save :: &
       & cline(im_local,jm_local)    , & ! coast line
       & idx_riv(im_local,jm_local)  , & ! River gird idx
       & idy_riv(im_local,jm_local)  , & !            idy
       & dnum_riv(im_local,jm_local) , & ! The number of Duplicated grids
       & cline_riv(im_riv,jm_riv)
  
  real(kind = r_size),save :: & 
       &  aam2d(im_local,jm_local)   , & ! vertical average of aam
       &  advua(im_local,jm_local)   , & ! sum of the 2nd, 3rd and 4th terms in eq (18)
       &  advva(im_local,jm_local)   , & ! sum of the 2nd, 3rd and 4th terms in eq (19)
       &  adx2d(im_local,jm_local)   , & ! vertical integral of advx
       &  ady2d(im_local,jm_local)   , & ! vertical integral of advy
       &  art(im_local,jm_local)     , & ! cell area centered on T grid points
       &  aru(im_local,jm_local)     , & ! cell area centered on U grid points
       &  arv(im_local,jm_local)     , & ! cell area centered on V grid points
       &  cbc(im_local,jm_local)     , & ! bottom friction coefficient
       &  cor(im_local,jm_local)     , & ! coriolis parameter
       &  d(im_local,jm_local)       , & ! h+el
       &  drx2d(im_local,jm_local)   , & ! vertical integral of drhox
       &  dry2d(im_local,jm_local)   , & ! vertical integral of drhoy
       &  dt(im_local,jm_local)      , & ! h+et
       &  dum(im_local,jm_local)     , & ! mask for u velocity
       &  dvm(im_local,jm_local)     , & ! mask for v velocity
       &  dx(im_local,jm_local)      , & ! grid spacing in x
       &  dy(im_local,jm_local)      , & ! grid spacing in y
       &  east_c(im_local,jm_local)  , & ! horizontal coordinate of cell corner points in x
       &  east_e(im_local,jm_local)  , & ! horizontal coordinate of elevation points in x
       &  east_u(im_local,jm_local)  , & ! horizontal coordinate of U points in x
       &  east_v(im_local,jm_local)  , & ! horizontal coordinate of V points in x
       &  e_atmos(im_local,jm_local) , & ! atmospheric pressure
       &  egb(im_local,jm_local)     , & ! surface elevation use for pressure gradient at time n-1
       &  egf(im_local,jm_local)     , & ! surface elevation use for pressure gradient at time n+1
       &  el(im_local,jm_local)      , & ! surface elevation used in the external mode at time n
       &  elb(im_local,jm_local)     , & ! surface elevation used in the external mode at time n-1
       &  elf(im_local,jm_local)     , & ! surface elevation used in the external mode at time n+1
       &  el_dave(im_local,jm_local) , & ! surface elevation used in the external mode averaged 1 day
       &  el_ave(im_local,jm_local)  , & ! surface elevation averaged over assimilation window
       &  el_iau(im_local,jm_local)  , & ! surface elevation increment
       &  et(im_local,jm_local)      , & ! surface elevation used in the internal mode at time n
       &  etb(im_local,jm_local)     , & ! surface elevation used in the internal mode at time n-1
       &  etf(im_local,jm_local)     , & ! surface elevation used in the internal mode at time n+1
       &  fluxua(im_local,jm_local)  , &
       &  fluxva(im_local,jm_local)  , & 
       &  fsm(im_local,jm_local)     , & ! mask for scalar variables
       &  h(im_local,jm_local)       , & ! bottom depth
       &  north_c(im_local,jm_local) , & ! horizontal coordinate of cell corner points in y
       &  north_e(im_local,jm_local) , & ! horizontal coordinate of elevation points in y
       &  north_u(im_local,jm_local) , & ! horizontal coordinate of U points in y
       &  north_v(im_local,jm_local) , & ! horizontal coordinate of V points in y
       &  psi(im_local,jm_local)     , & 
       &  ssurf(im_local,jm_local)   , & 
       &  swrad(im_local,jm_local)   , & ! short wave radiation incident on the ocean surface
       &  lwrad(im_local,jm_local)   , & ! long wave radiation
       &  lhf(im_local,jm_local)     , & !latent heat flux [W/m^2]
       &  lhf_dave(im_local,jm_local), & !latent heat flux [W/m^2] averaged within 1day
       &  shf(im_local,jm_local)     , & !sensible heat flux [W/m^2]
       &  shf_dave(im_local,jm_local), & !sensible heat flux [W/m^2] averaged within 1day
       &  swr(im_local,jm_local)     , & !shortwave radiation [W/m^2]
       &  swr_dave(im_local,jm_local), & !shortwave radiation [W/m^2] averaged within 1day
       &  lwr(im_local,jm_local)     , & !longwave radiation [W/m^2]
       &  lwr_dave(im_local,jm_local), & !longwave radiation [W/m^2] averaged within 1day
       &  vfluxb(im_local,jm_local)  , & ! volume flux through water column surface at time n-1
       &  vfluxf(im_local,jm_local)  , & ! volume flux through water column surface at time n+1
       &  tps(im_local,jm_local)     , & 
       &  tsurf(im_local,jm_local)   , & 
       &  ua(im_local,jm_local)      , & ! vertical mean of u at time n
       &  uab(im_local,jm_local)     , & ! vertical mean of u at time n-1
       &  uaf(im_local,jm_local)     , & ! vertical mean of u at time n+1
       &  ua_ave(im_local,jm_local)  , & ! vertical mean of u averaged over assimilation window
       &  ua_iau(im_local,jm_local)  , & ! vertical mean of u increment
       &  utb(im_local,jm_local)     , & ! ua time averaged over the interval dti at time n-1
       &  utf(im_local,jm_local)     , & ! ua time averaged over the interval dti at time n+1
       &  va(im_local,jm_local)      , & ! vertical mean of v at time n
       &  vab(im_local,jm_local)     , & ! vertical mean of v at time n-1
       &  vaf(im_local,jm_local)     , & ! vertical mean of v at time n+1
       &  va_ave(im_local,jm_local)  , & ! vertical mean of v averaged over assimilation window
       &  va_iau(im_local,jm_local)  , & ! difference of vertical mean of v for iau analysis
       &  vtb(im_local,jm_local)     , & ! va time averaged over the interval dti at time n-1
       &  vtf(im_local,jm_local)     , & ! va time averaged over the interval dti at time n+1
       &  wssurf(im_local,jm_local)  , & ! <ws(0)> salinity flux at the surface
       &  wtsurf(im_local,jm_local)  , & ! <wt(0)> temperature flux at the surface
       &  wubot(im_local,jm_local)   , & ! x-momentum flux at the bottom
       &  wusurf(im_local,jm_local)  , & ! <wu(0)> momentum flux at the surface
       &  wvbot(im_local,jm_local)   , & ! y-momentum flux at the bottom
       &  wvsurf(im_local,jm_local)  , & ! <wv(0)> momentum flux at the surface
       &  sfcterm(im_local,jm_local) , & ! surface forcing term
       &  tsfcf(im_local,jm_local)   , & ! temperature surface forcing term
       &  tsfc(im_local,jm_local)    , & !
       &  tsfcb(im_local,jm_local)   , & !
       &  tsfc_dave(im_local,jm_local), & ! temperature surface forcing averaged within 1day
       &  ssfcf(im_local,jm_local)    , & ! salinity surface forcing term
       &  ssfc(im_local,jm_local)     , & !
       &  ssfcb(im_local,jm_local)    , & !
       &  ssfc_dave(im_local,jm_local), & ! salinity surface forcing averaged within 1day
       !Atmosphere
       &  land_atm(im_atm,jm_atm)         ! Atmosphere land-sea mask

      
  !_______________________________________________________________________
  ! 3-D arrays
    real(kind = r_size),save :: & 
       &  dz(im_local,jm_local,kb)          , & ! z(i,j,k)-z(i,j,k+1)
       &  dzz(im_local,jm_local,kb)         , & ! zz(i,j,k)-zz(i,j,k+1)
       &  z(im_local,jm_local,kb)           , & ! sigma coordinate from z=0 (surface) to z=-1 (bottom)
       &  zz(im_local,jm_local,kb)          , & ! sigma coordinate, intermediate between z
       &  aam(im_local,jm_local,kb)         , & ! horizontal kinematic viscosity
       &  aam_dave(im_local,jm_local,kb)    , & ! horizontal kinematic viscosity averaged within 1 day
       &  advx(im_local,jm_local,kb)        , & ! x-horizontal advection and diffusion terms
       &  advy(im_local,jm_local,kb)        , & ! y-horizontal advection and diffusion terms
       &  drhox(im_local,jm_local,kb)       , & ! x-component of the internal baroclinic pressure
       &  drhoy(im_local,jm_local,kb)       , & ! y-component of the internal baroclinic pressure
       &  dtef(im_local,jm_local,kb)        , & 
       &  kh(im_local,jm_local,kb)          , & ! vertical diffusivity
       &  kh_dave(im_local,jm_local,kb)     , & ! vertical diffusivity averaged within 1 day
       &  km(im_local,jm_local,kb)          , & ! vertical kinematic viscosity
       &  km_dave(im_local,jm_local,kb)     , & ! vertical kinematic viscosity averaged within 1 day
       &  kq(im_local,jm_local,kb)          , & 
       &  l(im_local,jm_local,kb)           , & ! turbulence length scale
       &  q2b(im_local,jm_local,kb)         , & ! twice the turbulent kinetic energy at time n-1
       &  q2(im_local,jm_local,kb)          , & ! twice the turbulent kinetic energy at time n
       &  q2lb(im_local,jm_local,kb)        , & ! q2 x l at time n-1
       &  q2l(im_local,jm_local,kb)         , & ! q2 x l at time n
       &  rho(im_local,jm_local,kb)         , & ! density
       &  rmean(im_local,jm_local,kb)       , & ! Horizontal averaged rmean climatology
       &  sb(im_local,jm_local,kb)          , & ! salinity at time n-1
       &  sclim(im_local,jm_local,kb)       , & ! horizontally averaged salinity
       &  s(im_local,jm_local,kb)           , & ! salinity at time n
       &  s_dave(im_local,jm_local,kb)      , & ! salinity averaged within 1 day
       &  s_ave(im_local,jm_local,kb)       , & ! salinity averaged over assimilation window
       &  s_iau(im_local,jm_local,kb)       , & ! salinity increment
       &  dsdt(im_local,jm_local,kb)        , & ! salinity tendency in 1 day
       &  dsdt_dave(im_local,jm_local,kb)   , & ! salinity tendency in 1 day averaged within 1day
       &  tb(im_local,jm_local,kb)          , & ! temperature at time n-1
       &  tclim(im_local,jm_local,kb)       , & ! horizontally averaged temperature
       &  t(im_local,jm_local,kb)           , & ! temperature at time n
       &  t_dave(im_local,jm_local,kb)      , & ! temperature averaged within 1 day
       &  t_ave(im_local,jm_local,kb)       , & ! temperature averaged over assimilation window
       &  t_iau(im_local,jm_local,kb)       , & ! temperature increment
       &  dtdt(im_local,jm_local,kb)        , & ! temperature tendency in 1 day
       &  dtdt_dave(im_local,jm_local,kb)   , & !                               averaged within 1 day
       &  ub(im_local,jm_local,kb)          , & ! zonal velocity at time n-1
       &  uf(im_local,jm_local,kb)          , & ! zonal velocity at time n+1
       &  u(im_local,jm_local,kb)           , & ! zonal velocity at time n
       &  u_dave(im_local,jm_local,kb)      , & ! zonal velocity averaged within 1 day
       &  u_ave(im_local,jm_local,kb)       , & ! zonal velocity averaged over assimilation window
       &  u_iau(im_local,jm_local,kb)       , & ! zonal velocity increment 
       &  vb(im_local,jm_local,kb)          , & ! meridional velocity at time n-1
       &  vf(im_local,jm_local,kb)          , & ! meridional velocity at time n+1
       &  v(im_local,jm_local,kb)           , & ! meridional velocity at time n
       &  v_dave(im_local,jm_local,kb)      , & ! meridional velocity averaged within 1 day
       &  v_ave(im_local,jm_local,kb)       , & ! meridional velocity averaged over assimilation window
       &  v_iau(im_local,jm_local,kb)       , & ! meridional velocity increment 
       &  w(im_local,jm_local,kb)           , & ! sigma coordinate vertical velocity
       &  w_dave(im_local,jm_local,kb)      , & ! sigma coordinate vertical velocity averaged within 1 day
       &  wr(im_local,jm_local,kb)          , & ! real (z coordinate) vertical velocity
       &  wr_dave(im_local,jm_local,kb)         ! real (z coordinate) vertical velocity averaged within 1 day
  ! conservation each term
  real(kind = r_size),save :: &
       &  xadvterm(im_local,jm_local,kb)    , & ! advection term - x direction
       &  yadvterm(im_local,jm_local,kb)    , & !                - y direction
       &  zadvterm(im_local,jm_local,kb)    , & !                - z direction
       &  advterm(im_local,jm_local,kb)     , & ! advection term - x+y+z direction
       &  xdifterm(im_local,jm_local,kb)    , & ! diffusion term - x direction
       &  ydifterm(im_local,jm_local,kb)    , & !                - y direction
       &  zdifterm(im_local,jm_local,kb)    , & !                - z direction
       &  radterm(im_local,jm_local,kb)     , & ! radiation term
       &  txadv(im_local,jm_local,kb)       , & ! temperature advection term -x direction
       &  tyadv(im_local,jm_local,kb)       , & !                            -y direction
       &  tzadv(im_local,jm_local,kb)       , & !                            -z direction 
       &  txdif(im_local,jm_local,kb)       , & ! temperature diffusion term -x direction
       &  tydif(im_local,jm_local,kb)       , & ! 
       &  tzdif(im_local,jm_local,kb)       , & !
       &  qz(im_local,jm_local,kb)          , & ! shortwave penetration term
       &  troff(im_local,jm_local,kb)       , & ! temperature round off
       &  tnudge(im_local,jm_local,kb)      , & ! temperature nuding
       &  txadv_dave(im_local,jm_local,kb)  , & ! daily average of temperature advection term
       &  tyadv_dave(im_local,jm_local,kb)  , & ! 
       &  tzadv_dave(im_local,jm_local,kb)  , & ! 
       &  txdif_dave(im_local,jm_local,kb)  , & !                              diffusion term
       &  tydif_dave(im_local,jm_local,kb)  , & ! 
       &  tzdif_dave(im_local,jm_local,kb)  , & !
       &  qz_dave(im_local,jm_local,kb)     , & ! daily average of shortwave penetration term
       &  troff_dave(im_local,jm_local,kb)  , & ! daily temperature round off
       &  tnudge_dave(im_local,jm_local,kb) , & ! daily temperature nuding
       &  sxadv(im_local,jm_local,kb)       , & ! salinity advection term -x direction
       &  syadv(im_local,jm_local,kb)       , & !                         -y direction
       &  szadv(im_local,jm_local,kb)       , & !                         -z direction
       &  sxdif(im_local,jm_local,kb)       , & ! salinity diffusion term - xdirection
       &  sydif(im_local,jm_local,kb)       , & !                         - ydirection
       &  szdif(im_local,jm_local,kb)       , & !                         - zdirection
       &  sroff(im_local,jm_local,kb)       , & ! salinity round off
       &  snudge(im_local,jm_local,kb)      , & ! salinity nudging
       &  sxadv_dave(im_local,jm_local,kb)  , & ! daily average of salinity advection term
       &  syadv_dave(im_local,jm_local,kb)  , & ! 
       &  szadv_dave(im_local,jm_local,kb)  , & ! 
       &  sxdif_dave(im_local,jm_local,kb)  , & !                           diffusion term
       &  sydif_dave(im_local,jm_local,kb)  , & ! 
       &  szdif_dave(im_local,jm_local,kb)  , & !
       &  sroff_dave(im_local,jm_local,kb)  , & ! daily salinify round off
       &  snudge_dave(im_local,jm_local,kb)     ! daily salinity nuding

  !_______________________________________________________________________
  ! 1 and 2-D boundary value arrays

  real(kind = r_size),save :: & 
       &  ele(jm_local),ele0(jm_local),ele1(jm_local)         , & ! elevation at the eastern open boundary
       &  eln(im_local),eln0(im_local),eln1(im_local)         , & ! elevation at the northern open boundary        
       &  els(im_local),els0(im_local),els1(im_local)         , & ! elevation at the southern open boundary
       &  elw(jm_local),elw0(jm_local),elw1(jm_local)         , & ! elevation at the western open boundary
       &  sbe(jm_local,kb),sbe0(jm_local,kb),sbe1(jm_local,kb), & ! salinity at the eastern open boundary
       &  sbn(im_local,kb),sbn0(im_local,kb),sbn1(im_local,kb), & ! salinity at the northern open boundary
       &  sbs(im_local,kb),sbs0(im_local,kb),sbs1(im_local,kb), & ! salinity at the southern open boundary
       &  sbw(jm_local,kb),sbw0(jm_local,kb),sbw1(jm_local,kb), & ! salinity at the western open boundary
       &  tbe(jm_local,kb),tbe0(jm_local,kb),tbe1(jm_local,kb), & ! temperature at the eastern open boundary
       &  tbn(im_local,kb),tbn0(im_local,kb),tbn1(im_local,kb), & ! temperature at the northern open boundary
       &  tbs(im_local,kb),tbs0(im_local,kb),tbs1(im_local,kb), & ! temperature at the southern open boundary
       &  tbw(jm_local,kb),tbw0(jm_local,kb),tbw1(jm_local,kb), & ! temperature at the western open boundary
       &  uabe(jm_local),uabe0(jm_local),uabe1(jm_local)      , & ! vertical mean of u at the eastern open boundary
       &  uabn(im_local),uabn0(im_local),uabn1(im_local)      , & ! vertical mean of u at the northern open boundary
       &  uabs(im_local),uabs0(im_local),uabs1(im_local)      , & ! vertical mean of u at the southern open boundary
       &  uabw(jm_local),uabw0(jm_local),uabw1(jm_local)      , & ! vertical mean of u at the western open boundary
       &  ube(jm_local,kb),ube0(jm_local,kb),ube1(jm_local,kb), & ! u at the eastern open boundary
       &  ubn(im_local,kb),ubn0(im_local,kb),ubn1(im_local,kb), & ! u at the northern open boundary
       &  ubs(im_local,kb),ubs0(im_local,kb),ubs1(im_local,kb), & ! u at the southern open boundary
       &  ubw(jm_local,kb),ubw0(jm_local,kb),ubw1(jm_local,kb), & ! u at the western open boundary
       &  vabe(jm_local),vabe0(jm_local),vabe1(jm_local)      , & ! vertical mean of v at the eastern open boundary
       &  vabn(im_local),vabn0(im_local),vabn1(im_local)      , & ! vertical mean of v at the northern open boundary
       &  vabs(im_local),vabs0(im_local),vabs1(im_local)      , & ! vertical mean of v at the southern open boundary
       &  vabw(jm_local),vabw0(jm_local),vabw1(jm_local)      , & ! vertical mean of v at the western open boundary
       &  vbe(jm_local,kb),vbe0(jm_local,kb),vbe1(jm_local,kb), & ! v at the eastern open boundary
       &  vbn(im_local,kb),vbn0(im_local,kb),vbn1(im_local,kb), & ! v at the northern open boundary
       &  vbs(im_local,kb),vbs0(im_local,kb),vbs1(im_local,kb), & ! v at the southern open boundary
       &  vbw(jm_local,kb),vbw0(jm_local,kb),vbw1(jm_local,kb)    ! v at the western open boundary

  !_______________________________________________________________________
  ! Character variables
  character(26),save :: &
       & time_start      ! date and time of start of initial run of model
  
  character(40),save :: &
       & source, &
       & title
  
  character(120),save :: &
       &  netcdf_file    , &
       &  read_rst_file  , &
       &  read_iau_file  , &
       &  write_rst_file , &
       &  write_iau_file
  
  !_______________________________________________________________________
  ! Logical variables
  logical,save :: &
       & lramp  , & 
       & lfrs   , &
       & lmsm   , &
       & lrtvf  , &
       & lmynnf , &
       & lmbws  , &
       & liwbrk , &
       & lroff  , &
       & ltide  , &
       & lwbrk  , &
       & lpnetcdf

  !_______________________________________________________________________
  ! Additional 2-D arrays

  real(kind = r_size),save :: &
       &  windu(im_local,jm_local)     , & ! east-west component of wind
       &  windu0(im_local,jm_local)    , & 
       &  windu1(im_local,jm_local)    , & 
       &  windu_dave(im_local,jm_local), & 
       &  windv(im_local,jm_local)     , & ! north-south component of wind
       &  windv0(im_local,jm_local)    , & 
       &  windv1(im_local,jm_local)    , & 
       &  windv_dave(im_local,jm_local), & 
       &  winds(im_local,jm_local)     , & ! wind speed
       &  winds0(im_local,jm_local)    , & 
       &  winds1(im_local,jm_local)    , & 
       &  winds_dave(im_local,jm_local), & 
       &  tauu(im_local,jm_local)      , & ! east-west component of wind stress
       &  tauu_dave(im_local,jm_local) , & 
       &  tauv(im_local,jm_local)      , & ! north-south component of wind stress
       &  tauv_dave(im_local,jm_local) , &
       &  taus(im_local,jm_local)      , & ! magnitude of wind stress
       &  taus_dave(im_local,jm_local) , &
       &  airt(im_local,jm_local)      , & ! air temperature [K]
       &  airt0(im_local,jm_local)     , & 
       &  airt1(im_local,jm_local)     , & 
       &  ta(im_local,jm_local)        , & ! air temperature [degree C]
       &  ta_dave(im_local,jm_local)   , &
       &  airh(im_local,jm_local)      , & ! air humidity [g/g]
       &  airh0(im_local,jm_local)     , &
       &  airh1(im_local,jm_local)     , &
       &  qa(im_local,jm_local)        , & ! air humidity [g/kg]
       &  qa_dave(im_local,jm_local)   , & 
       &  qs(im_local,jm_local)        , & ! surface saturated specific humidity [g/kg]
       &  qs_dave(im_local,jm_local)   , &
       &  swrad0(im_local,jm_local)    , & ! short wave radiation [W/m^2]
       &  swrad1(im_local,jm_local)    , &
       &  lwrad0(im_local,jm_local)    , & ! long wave radiation [W/m^2]
       &  lwrad1(im_local,jm_local)    , & 
       &  slp(im_local,jm_local)       , & ! sea level pressure [Pa]
       &  slp0(im_local,jm_local)      , &
       &  slp1(im_local,jm_local)      , &
       &  evap(im_local,jm_local)      , & ! evaporation (E) [mm/day]
       &  evap0(im_local,jm_local)     , &
       &  evap1(im_local,jm_local)     , &
       &  evap_dave(im_local,jm_local) , &
       &  prep(im_local,jm_local)      , & ! precipitation (P) 
       &  prep0(im_local,jm_local)     , &
       &  prep1(im_local,jm_local)     , &
       &  prep_dave(im_local,jm_local) , &
       &  river(im_local,jm_local)     , & ! river discharge (R)
       &  river0(im_local,jm_local)    , &
       &  river1(im_local,jm_local)    , & 
       &  river_dave(im_local,jm_local), &
       &  fflux(im_local,jm_local)     , & ! fresh water flux (E-P-R)
       &  fflux0(im_local,jm_local)    , &
       &  fflux1(im_local,jm_local)    , &
       &  eflux(im_local,jm_local)     , & ! evaporation flux [m/s]
       &  eflux_dave(im_local,jm_local), &
       &  pflux(im_local,jm_local)     , & ! precipitation flux
       &  pflux_dave(im_local,jm_local), &
       &  rflux(im_local,jm_local)     , & ! river discharge flux
       &  rflux_dave(im_local,jm_local)

  !_______________________________________________________________________
  ! Additional 3-D arrays
  real(kind = r_size),save :: &
       &  tclimm(im_local,jm_local,kb)  , & ! monthly climatological temperature
       &  tclim0(im_local,jm_local,kb)  , &
       &  tclim1(im_local,jm_local,kb)  , &
       &  sclimm(im_local,jm_local,kb)  , & ! monthly climatological salinity
       &  sclim0(im_local,jm_local,kb)  , &
       &  sclim1(im_local,jm_local,kb)  , &
       &  tref(im_local,jm_local,kb)    , & ! reference temperature
       &  tref0(im_local,jm_local,kb)   , &
       &  tref1(im_local,jm_local,kb)   , &
       &  sref(im_local,jm_local,kb)    , & ! reference salinity
       &  sref0(im_local,jm_local,kb)   , &
       &  sref1(im_local,jm_local,kb)

end module common_pom_var
