!-------------------------------------------------------------------------
! Parameter Setting |
!-------------------------------------------------------------------------
!
! 06/04/2026 Shun OHISHI updated 
!
!-------------------------------------------------------------------------

MODULE common_setting

  IMPLICIT NONE

  !=======================================================================
  ! Parameters |
  !=======================================================================
  
  !-----------------------------------------------------------------------
  ! Variable size definitions
  !-----------------------------------------------------------------------
  INTEGER,PARAMETER :: r_size = kind(0.0e0)
  INTEGER,PARAMETER :: r_dble = kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl = kind(0.0e0)

  !-----------------------------------------------------------------------
  ! Constant parameters
  !-----------------------------------------------------------------------
  REAL(r_size),PARAMETER :: pi = REAL( 4.d0*atan(1.d0), r_size)
  REAL(r_size),PARAMETER :: re = REAL( 6371.3d3, r_size)
  REAL(r_size),PARAMETER :: r_omega = REAL( 7.292d-5, r_size)
  REAL(r_size),PARAMETER :: undef = REAL( 999.0d0, r_size)

  !-----------------------------------------------------------------------
  ! General parameters
  !-----------------------------------------------------------------------
  !INTEGER,PARAMETER :: nmem = 10    !Ensemble size
  INTEGER,PARAMETER :: nmem = 128   !Ensemble size
  INTEGER,PARAMETER :: nlon = 1442 !longitude
  INTEGER,PARAMETER :: nlat = 562  !latitude
  INTEGER,PARAMETER :: nlev = 75   !depth
  INTEGER,PARAMETER :: nv3d = 4    !u,v,t,s
  INTEGER,PARAMETER :: nv2d = 3    !zeta,ubar,vbar
  INTEGER,PARAMETER :: iv3d_u = 1
  INTEGER,PARAMETER :: iv3d_v = 2
  INTEGER,PARAMETER :: iv3d_t = 3
  INTEGER,PARAMETER :: iv3d_s = 4
  INTEGER,PARAMETER :: iv2d_z = 1
  INTEGER,PARAMETER :: iv2d_ubar = 2
  INTEGER,PARAMETER :: iv2d_vbar = 3
  INTEGER,PARAMETER :: nij0 = nlon*nlat
  INTEGER,PARAMETER :: nlevall = nlev*nv3d+nv2d
  INTEGER,PARAMETER :: ngpv = nij0*nlevall
  INTEGER,PARAMETER :: igrid_type=1 !1: Uniform, 2: Non-uniform grid
  REAL(r_size),PARAMETER :: hnosea = REAL(0.d0,r_size)
  LOGICAL,PARAMETER :: lqglobal = .true.

  INTEGER,PARAMETER :: id_u_obs = 2819
  INTEGER,PARAMETER :: id_v_obs = 2820
  INTEGER,PARAMETER :: id_t_obs = 3073
  INTEGER,PARAMETER :: id_z_obs = 2567
  INTEGER,PARAMETER :: id_s_obs = 3332

  INTEGER,PARAMETER :: nslots = 1 ! number of time slots for 4D-LETKF
  INTEGER,PARAMETER :: nbslot = 1 ! basetime slot

  INTEGER,PARAMETER :: iswitch_da = 1 !1: LETKF, 2: LPF, 3: LPFGM

  !--------------------------------------------------------------------
  ! Covariance inflation
  !--------------------------------------------------------------------
  !Multiplicative inflation
  REAL(r_size),PARAMETER :: cov_infl_mul = REAL( 1.00d0, r_size)

  !RTPP/RTPS
  REAL(r_size),PARAMETER :: ALPHA_RTPP = REAL( 0.9d0, r_size)
  REAL(r_size),PARAMETER :: ALPHA_RTPS = REAL( 0.0d0, r_size)

  !---LPF
  !Number of Resampling 
  INTEGER,PARAMETER :: nmonte = 100
  
  !Resampling method [SR:Systematic Resampling, MR: Multinominal Resampling]
  CHARACTER(2),PARAMETER :: RM = "MR"
  
  !Diagnoal Priority (Kotsuki et al. 2022)
  LOGICAL,PARAMETER :: DP = .true.
  
  !---LPFGM
  !Multiplicative inflation
  REAL(r_size),PARAMETER :: cov_infl_gm = REAL( 1.00d0, r_size)
  
  !---------------------------------------------------------------------
  ! Covariance Localization
  !---------------------------------------------------------------------
  REAL(r_size),PARAMETER :: sigma_obs  = REAL( 300.0d3, r_size)  ! horizontal [m]
  REAL(r_size),PARAMETER :: sigma_obsv = REAL( 100.0d0, r_size) ! vertical [m] !*zero or minus --> no v localization
  REAL(r_size),PARAMETER :: sigma_obst = REAL( 1.0d0, r_size)   ! Assimilation interval [day] * sigma_obst 
  REAL(r_size),PARAMETER :: dist_zero  = REAL( 2.d0*SQRT(10.0d0/3.0d0), r_size) * sigma_obs   !Cutoff scale (horizontal)
  REAL(r_size),PARAMETER :: dist_zerov = REAL( 2.d0*SQRT(10.0d0/3.0d0), r_size) * sigma_obsv !Curoff scale (vertical)
  
  !--------------------------------------------------------------------
  ! Adaptive Background & Observation Error Inflation (AOEI)
  !--------------------------------------------------------------------
  LOGICAL,PARAMETER :: AOEI = .true.
  
  !---------------------------------------------------------------------
  ! Gross error Check using innovation (y-Hx)
  !---------------------------------------------------------------------
  LOGICAL,PARAMETER :: LGE_IO=.false.
  REAL(r_size),PARAMETER :: hqc_min = REAL( -1.d0, r_size), hqc_max = REAL( 1.d0, r_size)     !SSH
  REAL(r_size),PARAMETER :: tqc_min = REAL( -5.d0, r_size), tqc_max = REAL( 5.d0, r_size)     !Temperature
  REAL(r_size),PARAMETER :: sqc_min = REAL( -2.d0, r_size), sqc_max = REAL( 2.d0, r_size)     !Salinity
  REAL(r_size),PARAMETER :: uqc_min = REAL( -2.d0, r_size), uqc_max = REAL( 2.d0, r_size)     !Zonal velocity
  REAL(r_size),PARAMETER :: vqc_min = REAL( -2.d0, r_size), vqc_max = REAL( 2.d0, r_size)     !Meridional velocity

  REAL(r_size),PARAMETER :: sstqc_min = tqc_min, sstqc_max = tqc_max !SST
  REAL(r_size),PARAMETER :: sssqc_min = sqc_min, sssqc_max = sqc_max !SSS
  REAL(r_size),PARAMETER :: ssuqc_min = uqc_min, ssuqc_max = uqc_max !SSU
  REAL(r_size),PARAMETER :: ssvqc_min = vqc_min, ssvqc_max = vqc_max !SSV

  !--------------------------------------------------------------------
  ! MPI
  !--------------------------------------------------------------------
  INTEGER,PARAMETER :: root=0
    
  !--------------------------------------------------------------------
  ! Round off for analysis
  !--------------------------------------------------------------------
  !Lower & Upper limit for Analysis T/S
  LOGICAL,PARAMETER :: ROFF = .true.
  REAL(r_size),PARAMETER :: tmin = REAL( -1.8d0, r_size), tmax= REAL( 50.d0, r_size) 
  REAL(r_size),PARAMETER :: smin = REAL( 25.d0, r_size), smax = REAL( 50.d0,r_size)

  !===================================================================
  ! Variable |
  !===================================================================

  !-------------------------------------------------------------------
  ! MPI
  !-------------------------------------------------------------------
  INTEGER,SAVE :: nprocs
  INTEGER,SAVE :: myrank
  INTEGER,SAVE :: root_out
  
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

  REAL(r_size),ALLOCATABLE,SAVE :: obshxfmean(:) !H(xfmean)
  REAL(r_size),ALLOCATABLE,SAVE :: obshxfsprd(:) !H(xfsprd)
  REAL(r_size),ALLOCATABLE,SAVE :: obsdep(:)     !Innovation: y-H(xfmean)
  REAL(r_size),ALLOCATABLE,SAVE :: obshdxf(:,:)  !dYf = H(dXf)

  !-----------------------------------------------------------------
  ! Timer |
  !-----------------------------------------------------------------
  
  INTEGER,SAVE :: rtimerxx,rtimer00,rtimer,rate
  REAL(r_dble),SAVE :: ctimer00,ctimer
  
END MODULE common_setting
