MODULE common_setting

  IMPLICIT NONE
  PUBLIC

  !===================================================================
  ! Parameters |
  !===================================================================
  
  !-----------------------------------------------------------------------
  ! Variable size definitions
  !-----------------------------------------------------------------------
  INTEGER,PARAMETER :: r_size=kind(0.0e0)
  INTEGER,PARAMETER :: r_dble=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  !-----------------------------------------------------------------------
  ! Constant parameters
  !-----------------------------------------------------------------------
  REAL(r_size),PARAMETER :: pi=4.d0*atan(1.d0)
  REAL(r_size),PARAMETER :: re=6371.3d3
  REAL(r_size),PARAMETER :: r_omega=7.292d-5
  REAL(r_size),PARAMETER :: undef=999.0d0

  !-----------------------------------------------------------------------
  ! General parameters
  !-----------------------------------------------------------------------
  ! INTEGER,PARAMETER :: nbv=10    !Ensemble size
  INTEGER,PARAMETER :: nbv=128   !Ensemble size
  INTEGER,PARAMETER :: nmonte=100 !Number of resampling for LPF
  INTEGER,PARAMETER :: nlon=372 !longitude
  INTEGER,PARAMETER :: nlat=362  !latitude
  INTEGER,PARAMETER :: nlev=75   !depth
  INTEGER,PARAMETER :: nv3d=4    !u,v,t,s
  INTEGER,PARAMETER :: nv2d=3    !zeta,ubar,vbar
  INTEGER,PARAMETER :: iv3d_u=1
  INTEGER,PARAMETER :: iv3d_v=2
  INTEGER,PARAMETER :: iv3d_t=3
  INTEGER,PARAMETER :: iv3d_s=4
  INTEGER,PARAMETER :: iv2d_z=1
  INTEGER,PARAMETER :: iv2d_ubar=2
  INTEGER,PARAMETER :: iv2d_vbar=3
  INTEGER,PARAMETER :: nij0=nlon*nlat
  INTEGER,PARAMETER :: nlevall=nlev*nv3d+nv2d
  INTEGER,PARAMETER :: ngpv=nij0*nlevall
  REAL(r_size),PARAMETER :: hnosea=0.d0
  LOGICAL,PARAMETER :: lqglobal=.false.

  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073
  INTEGER,PARAMETER :: id_z_obs=2567
  INTEGER,PARAMETER :: id_s_obs=3332

  INTEGER,PARAMETER :: nslots=1 ! number of time slots for 4D-LETKF
  INTEGER,PARAMETER :: nbslot=1 ! basetime slot

  INTEGER,PARAMETER :: iswitch_da=1 !1: LETKF, 2: LPF, 3: LPFGM
  
  !---------------------------------------------------------------------
  ! Localization
  !---------------------------------------------------------------------
  REAL(r_size),PARAMETER :: sigma_obs=300.0d3  ! horizontal [m]
  REAL(r_size),PARAMETER :: sigma_obsv=100.0d0 ! vertical [m] !*zero or minus --> no v localization
  REAL(r_size),PARAMETER :: sigma_obst=1.0d0   ! Assimilation interval [day] * sigma_obst 
  
  !---------------------------------------------------------------------
  ! Gross error Check using innovation (y-Hx)
  !---------------------------------------------------------------------
  LOGICAL,PARAMETER :: LGE_IO=.false.
  REAL(r_size),PARAMETER :: hqc_min=-1.d0,hqc_max=1.d0     !SSH
  REAL(r_size),PARAMETER :: sstqc_min=-5.d0,sstqc_max=5.d0 !SST
  REAL(r_size),PARAMETER :: tqc_min=-5.d0,tqc_max=5.d0     !Temperature
  REAL(r_size),PARAMETER :: sssqc_min=-1.d0,sssqc_max=1.d0 !SSS
  REAL(r_size),PARAMETER :: sqc_min=-2.d0,sqc_max=2.d0     !Salinity
  REAL(r_size),PARAMETER :: ssuqc_min=-1.d0,ssuqc_max=1.d0 !SSU
  REAL(r_size),PARAMETER :: uqc_min=-2.d0,uqc_max=2.d0     !Zonal velocity
  REAL(r_size),PARAMETER :: ssvqc_min=-1.d0,ssvqc_max=1.d0 !SSV
  REAL(r_size),PARAMETER :: vqc_min=-2.d0,vqc_max=2.d0     !Meridional velocity

  !--------------------------------------------------------------------
  ! Adaptive Background & Observation Error Inflation (AOEI)
  !--------------------------------------------------------------------
  LOGICAL,PARAMETER :: AOEI=.true.

  !--------------------------------------------------------------------
  ! Parallel NetCDF
  !--------------------------------------------------------------------
  LOGICAL,PARAMETER :: pnetcdf=.false.
  !LOGICAL,PARAMETER :: pnetcdf=.true. !*** Under construction
  
  !--------------------------------------------------------------------
  ! Round off for analysis
  !--------------------------------------------------------------------
  !Lower & Upper limit for Analysis T/S
  LOGICAL,PARAMETER :: ROFF=.true.
  REAL(r_size),PARAMETER :: tmin=-1.8d0,tmax=50.d0 
  REAL(r_size),PARAMETER :: smin=25.d0,smax=50.d0

  !--------------------------------------------------------------------
  ! Covariance inflation
  !--------------------------------------------------------------------
  !Multiplicative inflation
  REAL(r_size),PARAMETER :: cov_infl_mul = 1.00d0

  !RTPP/RTPS
  REAL(r_size),PARAMETER :: ALPHA_RTPP   = 0.9d0
  REAL(r_size),PARAMETER :: ALPHA_RTPS   = 0.0d0

  !Multiplicative inflation for LPFGM
  REAL(r_size),PARAMETER :: cov_infl_gm = 1.00d0

  !===================================================================
  ! Variable |
  !===================================================================

  !-------------------------------------------------------------------
  ! MPI
  !-------------------------------------------------------------------
  INTEGER,SAVE :: nprocs
  INTEGER,SAVE :: myrank
  INTEGER,SAVE :: file_unit

  !-------------------------------------------------------------------
  ! FILE TYPE |
  !-------------------------------------------------------------------

  INTEGER,SAVE :: file_type                         !1: restart/intermittent, 2: iau/IAU
  CHARACTER(10),SAVE :: var2df(nv2d),var2da(nv2d)   !2D-variable name
  CHARACTER(10),SAVE :: var3df(nv3d),var3da(nv3d)   !3D-variable name
  
  !-------------------------------------------------------------------
  ! Grid (Global)|
  !-------------------------------------------------------------------
  REAL(r_size),SAVE :: lon(nlon,nlat)
  REAL(r_size),SAVE :: lat(nlon,nlat)
  REAL(r_size),SAVE :: dlon(nlon,nlat)
  REAL(r_size),SAVE :: dlat(nlon,nlat)
  REAL(r_size),SAVE :: fsm(nlon,nlat)
  REAL(r_size),SAVE :: phi0(nlon,nlat)       !Bottom Depth [m] (positive)
  REAL(r_size),SAVE :: depth(nlon,nlat,nlev) !Depth [m] (negative)
  CHARACTER(4),SAVE :: element(nv3d+nv2d)

  !-------------------------------------------------------------------
  ! Grid (Local) |
  !-------------------------------------------------------------------
  
  INTEGER,SAVE :: nij1                    !Grid size for each processor 
  INTEGER,SAVE :: nij1max                 !Maximum grid size for each processor
  INTEGER,ALLOCATABLE,SAVE :: nij1node(:) !To save grid size for each processor (1:nproc)
  
  REAL(r_size),ALLOCATABLE,SAVE :: phi1(:) 
  REAL(r_size),ALLOCATABLE,SAVE :: lon1(:),lat1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: dlon1(:),dlat1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: ri1(:),rj1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: depth1(:,:)

  !-------------------------------------------------------------------
  ! Observation |
  !-------------------------------------------------------------------
  
  INTEGER,SAVE :: nobs
  INTEGER,SAVE :: nobsgrd(nlon,nlat)

  INTEGER,ALLOCATABLE,SAVE :: obsidx(:)
  INTEGER,ALLOCATABLE,SAVE :: obsidy(:)
  INTEGER,ALLOCATABLE,SAVE :: obselm(:)
  INTEGER,ALLOCATABLE,SAVE :: obsins(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslon(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslat(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslev(:) !Depth [m]: Negative
  REAL(r_size),ALLOCATABLE,SAVE :: obsdat(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obserr(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsdep(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obshdxf(:,:)

  !-----------------------------------------------------------------
  ! Cutoff scale |
  !-----------------------------------------------------------------
  
  REAL(r_size),SAVE :: dist_zero
  REAL(r_size),SAVE :: dist_zerov

END MODULE common_setting
