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
  REAL(r_dble),PARAMETER :: pi = REAL( 4.d0*atan(1.d0), r_dble)
  REAL(r_dble),PARAMETER :: re = REAL( 6371.3d3, r_dble)
  REAL(r_size),PARAMETER :: undef = REAL(-999.0d0, r_size)

  !-----------------------------------------------------------------------
  ! General parameters
  !-----------------------------------------------------------------------
  !INTEGER,PARAMETER :: nmem = 10    !Ensemble size
  INTEGER,PARAMETER :: nmem = 128   !Ensemble size
  INTEGER,PARAMETER :: nlon = 1442  !longitude
  INTEGER,PARAMETER :: nlat = 562   !latitude
  INTEGER,PARAMETER :: nlev = 75    !depth
  INTEGER,PARAMETER :: k_sfc=1      !surface layer
  INTEGER,PARAMETER :: k_btm=nlev-1 !deepest valid layer (k=nlev: dummy)
  INTEGER,PARAMETER :: nv3d = 4     !u,v,t,s
  INTEGER,PARAMETER :: nv2d = 3     !zeta,ubar,vbar
  INTEGER,PARAMETER :: iv3d_u = 1
  INTEGER,PARAMETER :: iv3d_v = 2
  INTEGER,PARAMETER :: iv3d_t = 3
  INTEGER,PARAMETER :: iv3d_s = 4
  INTEGER,PARAMETER :: iv2d_z = 1
  INTEGER,PARAMETER :: iv2d_ubar = 2
  INTEGER,PARAMETER :: iv2d_vbar = 3
  INTEGER,PARAMETER :: igrid_type=1 !1: Uniform, 2: Non-uniform grid
  REAL(r_size),PARAMETER :: land = REAL(0.d0,r_size)
  LOGICAL,PARAMETER :: lqglobal = .true.

  INTEGER,PARAMETER :: id_u_obs = 2819
  INTEGER,PARAMETER :: id_v_obs = 2820
  INTEGER,PARAMETER :: id_t_obs = 3073
  INTEGER,PARAMETER :: id_z_obs = 2567
  INTEGER,PARAMETER :: id_s_obs = 3332

  INTEGER,PARAMETER :: nslots = 1 ! Number of time slots for 4D-LETKF
  INTEGER,PARAMETER :: nbslot = 1 ! Baseline timeslot

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
  INTEGER,PARAMETER :: nmonte = 1
  
  !Resampling method [SR:Systematic Resampling, MR: Multinomial Resampling]
  CHARACTER(2),PARAMETER :: RM = "MR"
  
  !Diagonal Priority (Kotsuki et al. 2022)
  LOGICAL,PARAMETER :: DP = .true.

  !Resampling criteria:Resampling if meff < m_crit (>nmem: always, 0: never)
  REAL(r_dble),PARAMETER :: m_crit=REAL(nmem+1,r_dble)
  
  !Local random number
  !True : Different random number at each grid point
  !False: Same random number at all grid points (spatially smooth resampling)
  LOGICAL,PARAMETER :: LRAND=.false. 
  
  !---LPFGM
  !Multiplicative inflation
  REAL(r_size),PARAMETER :: cov_infl_gm = REAL( 1.00d0, r_size)
  
  !---------------------------------------------------------------------
  ! Covariance Localization
  !---------------------------------------------------------------------
  REAL(r_size),PARAMETER :: sigma_obs  = REAL( 300.0d3, r_size) ! horizontal localization scale [m]
  REAL(r_size),PARAMETER :: sigma_obsv = REAL( 100.0d0, r_size) ! vertical                      [m] !*zero or minus --> no v localization
  REAL(r_size),PARAMETER :: sigma_obst = REAL( 1.0d0, r_size)   ! temporal                      [slot] 
  REAL(r_size),PARAMETER :: dist_zero  = REAL( 2.d0*SQRT(10.0d0/3.0d0), r_size) * sigma_obs  !Cutoff scale (horizontal)
  REAL(r_size),PARAMETER :: dist_zerov = REAL( 2.d0*SQRT(10.0d0/3.0d0), r_size) * sigma_obsv !Cutoff scale (vertical)
  
  !--------------------------------------------------------------------
  ! Adaptive Background & Observation Error Inflation (AOEI)
  !--------------------------------------------------------------------
  LOGICAL,PARAMETER :: AOEI = .true.
  
  !--------------------------------------------------------------------
  ! MPI
  !--------------------------------------------------------------------
  INTEGER,PARAMETER :: root=0
    
  !--------------------------------------------------------------------
  ! Round off for analysis
  !--------------------------------------------------------------------
  !Lower & Upper limit for Analysis T/S
  LOGICAL,PARAMETER :: ROFF = .true.
  REAL(r_size),PARAMETER :: tmin = REAL( -1.8d0, r_size), tmax = REAL( 50.d0, r_size) 
  REAL(r_size),PARAMETER :: smin = REAL( 25.d0, r_size),  smax = REAL( 50.d0, r_size)

  !--------------------------------------------------------------------
  ! Output data in obs space
  !--------------------------------------------------------------------  

  CHARACTER(12),PARAMETER :: obs_assim_filename="obs_assim.nc" !Assimilated obs
  CHARACTER(10),PARAMETER :: hx_filename="hxNNNNN.nc"
  CHARACTER(9),PARAMETER :: hx_mean_filename="hxmean.nc"
  CHARACTER(9),PARAMETER :: hx_sprd_filename="hxsprd.nc"
    
  !===================================================================
  ! Variable |
  !===================================================================

  !-------------------------------------------------------------------
  ! MPI
  !-------------------------------------------------------------------
  INTEGER,SAVE :: nprocs
  INTEGER,SAVE :: myrank
  INTEGER,SAVE :: root_out
  INTEGER,SAVE :: MPI_R_SIZE

  !-------------------------------------------------------------------
  ! MPI shared memory (one copy of hdxfT per node)
  !-------------------------------------------------------------------
  INTEGER,SAVE :: comm_node    !Communicator within a node
  INTEGER,SAVE :: comm_leader  !Communicator of node leaders (MPI_COMM_NULL for others)
  INTEGER,SAVE :: node_rank    !Rank within a node
  INTEGER,SAVE :: win_hdxfT    !MPI window for hdxfT
  
  !-------------------------------------------------------------------
  ! FILE TYPE |
  !-------------------------------------------------------------------

  INTEGER,SAVE :: file_type                         !1: restart/intermittent, 2: iau/IAU
  CHARACTER(10),SAVE :: var2df(nv2d),var2da(nv2d)   !2D-variable name
  CHARACTER(10),SAVE :: var3df(nv3d),var3da(nv3d)   !3D-variable name
  
  !-------------------------------------------------------------------
  ! Grid (Global)|
  !-------------------------------------------------------------------
  REAL(r_size),SAVE :: lont(nlon,nlat)   !longitude at t point [degree]
  REAL(r_size),SAVE :: lonu(nlon,nlat)   !longitude at u point [degree]
  REAL(r_size),SAVE :: lonv(nlon,nlat)   !longitude at v point [degree]

  REAL(r_size),SAVE :: latt(nlon,nlat)   !latitude at t point  [degree]
  REAL(r_size),SAVE :: latu(nlon,nlat)   !latitude at u point  [degree]
  REAL(r_size),SAVE :: latv(nlon,nlat)   !latitude at v point  [degree]

  REAL(r_size),SAVE :: dlont(nlon,nlat)  !Longitude grid interval [m]
  REAL(r_size),SAVE :: dlonu(nlon,nlat)  !Longitude grid interval [m]
  REAL(r_size),SAVE :: dlonv(nlon,nlat)  !Longitude grid interval [m]

  REAL(r_size),SAVE :: dlatt(nlon,nlat)  !Latitude grid interval [m]
  REAL(r_size),SAVE :: dlatu(nlon,nlat)  !Latitude grid interval [m]
  REAL(r_size),SAVE :: dlatv(nlon,nlat)  !Latitude grid interval [m]

  REAL(r_size),SAVE :: maskt(nlon,nlat)  !mask at t point
  REAL(r_size),SAVE :: masku(nlon,nlat)  !mask at u point
  REAL(r_size),SAVE :: maskv(nlon,nlat)  !mask at v point

  REAL(r_size),SAVE :: tpg(nlon,nlat)        !Topography [m] (positive)
  REAL(r_size),SAVE :: depth(nlon,nlat,nlev) !Depth [m] (negative)

  !-------------------------------------------------------------------
  ! Grid (Local) |
  !-------------------------------------------------------------------
  
  INTEGER,SAVE :: nij1                    !Grid size for each processor 
  INTEGER,SAVE :: nij1max                 !Maximum grid size for each processor
  INTEGER,ALLOCATABLE,SAVE :: nij1node(:) !To save grid size for each processor (1:nproc)

  REAL(r_size),ALLOCATABLE,SAVE :: ri1(:),rj1(:)
  
  REAL(r_size),ALLOCATABLE,SAVE :: lont1(:),latt1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: lonu1(:),latu1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: lonv1(:),latv1(:)

  REAL(r_size),ALLOCATABLE,SAVE :: dlont1(:),dlatt1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: dlonu1(:),dlatu1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: dlonv1(:),dlatv1(:)


  REAL(r_size),ALLOCATABLE,SAVE :: maskt1(:) 
  REAL(r_size),ALLOCATABLE,SAVE :: masku1(:) 
  REAL(r_size),ALLOCATABLE,SAVE :: maskv1(:) 
  
  REAL(r_size),ALLOCATABLE,SAVE :: tpg1(:) 

  REAL(r_size),ALLOCATABLE,SAVE :: depth1(:,:)

  !-------------------------------------------------------------------
  ! Observation |
  !-------------------------------------------------------------------
  
  INTEGER,SAVE :: nobs               !Obs size
  INTEGER,SAVE :: nobsgrd(nlon,nlat) !Cumulative obs size up to a grid (ilon,ilat) 

  INTEGER,ALLOCATABLE,SAVE :: obsidt(:) !Time id (islot)
  INTEGER,ALLOCATABLE,SAVE :: obsidxt(:),obsidyt(:) !Longitude and latitude id at t point
  INTEGER,ALLOCATABLE,SAVE :: obsidxu(:),obsidyu(:) !Longitude and latitude id at u point
  INTEGER,ALLOCATABLE,SAVE :: obsidxv(:),obsidyv(:) !Longitude and latitude id at v point
  INTEGER,ALLOCATABLE,SAVE :: obselm(:) !Element
  INTEGER,ALLOCATABLE,SAVE :: obsins(:) !Instrument
  INTEGER,ALLOCATABLE,SAVE :: obsqc(:)  !QC flag (1:OK, 0:NG)
  REAL(r_size),ALLOCATABLE,SAVE :: obslon(:) !Longitude
  REAL(r_size),ALLOCATABLE,SAVE :: obslat(:) !Latitude
  REAL(r_size),ALLOCATABLE,SAVE :: obslev(:) !Depth [m]: Negative
  REAL(r_size),ALLOCATABLE,SAVE :: obsdat(:) !Data
  REAL(r_size),ALLOCATABLE,SAVE :: obserr(:) !Error

  !-------------------------------------------------------------------
  ! State vector in obs space |
  !-------------------------------------------------------------------

  REAL(r_size),ALLOCATABLE,SAVE :: hxfmean(:) !Ensemble mean   of H(Xf)
  REAL(r_size),ALLOCATABLE,SAVE :: hxamean(:) !Ensemble mean   of H(Xa)
  REAL(r_size),ALLOCATABLE,SAVE :: hxfsprd(:) !Ensemble spread of H(Xf)
  REAL(r_size),ALLOCATABLE,SAVE :: hxasprd(:) !Ensemble spread of H(Xa)
  REAL(r_size),POINTER,CONTIGUOUS,SAVE :: hdxfT(:,:) => NULL() !dYf^T (shared within a node)
  REAL(r_size),POINTER,CONTIGUOUS,SAVE :: hdxfT_flat(:) => NULL() !1-D view of hdxf^T (64-bit index)
  REAL(r_size),ALLOCATABLE,SAVE :: hdxaT_flat(:) !dYa^T: dYa(iobs,imem) = hdxaT_flat(imem+(iobs-1)*nmem) (root_out only)
  REAL(r_size),ALLOCATABLE,SAVE :: dob(:)     !Innovation: y-mean(H(Xf))
  
  !-----------------------------------------------------------------
  ! Timer |
  !-----------------------------------------------------------------
  
  INTEGER(kind = 8),SAVE :: rtimerxx,rtimer00,rtimer,rate
  REAL(r_dble),SAVE :: ctimer00,ctimer
  
END MODULE common_setting
