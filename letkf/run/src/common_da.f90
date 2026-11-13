MODULE common_da
  !=======================================================================
  !
  ! [PURPOSE:] Module for DA (LETKF/LPF/LPFGM) with POM
  !
  ! [HISTORY:]
  !   2009/01/26 Takemasa Miyoshi  created
  !   2009/02/03 Takemasa Miyoshi  modified for ROMS
  !   2011/01/26 Yasumasa Miyazawa modified for POM (check 'pom' or 'POM')
  !   2025/07/31 Shun Ohishi       added OpenMP
  !   2025/10/16 Shun Ohishi       changed precision
  !   2025/10/26 Shun Ohishi       integrated LETKF/LPF/LPFGM
  !   2026/05/29 Shun Ohishi       minor update to avoid memory issue
  !   2026/11/13 Shun Ohishi       updated by using Claude
  !=======================================================================

CONTAINS

  !-----------------------------------------------------------------------
  ! Data Assimilation System | Main loop
  !-----------------------------------------------------------------------

  SUBROUTINE da_main(ens3d,ens2d)

    !$ USE OMP_LIB
    USE common
    USE common_letkf
    USE common_lpfgm
    USE common_setting
    IMPLICIT NONE

    INTEGER i,k,k_min,k_max
    INTEGER nobsh     !#Obs. within horizontal cutoff at grid point i
    INTEGER nobsl     !#Local obs. at (i,k)
    INTEGER iseed_grd !Seed for random number in LPF/LPFGM (0: same at all grid points)

    !---Local obs. list
    INTEGER,ALLOCATABLE      :: iobs_h(:) !Global obs. index of horizontal local obs. [1:nobsh]
    REAL(r_size),ALLOCATABLE :: dist_h(:) !Horizontal distance                      [1:nobsh]
    INTEGER,ALLOCATABLE      :: iobs_l(:) !Global obs. index of local obs. at level k [1:nobsl]
    REAL(r_size),ALLOCATABLE :: dist_l(:) !Horizontal distance of local obs.           [1:nobsl]

    REAL(r_size) fcst3dm(nv3d),fcst2dm(nv2d)       !Forecast ensemble mean
    REAL(r_size) fcst3ds(nv3d),fcst2ds(nv2d)       !Forecast ensemble spread
    REAL(r_size) dxf3d(nmem,nv3d),dxf2d(nmem,nv2d) !Forecast ensemble perturbation

    REAL(r_size) anal3dm(nv3d),anal2dm(nv2d)       !Analysis ensemble mean
    REAL(r_size) anal3ds(nv3d),anal2ds(nv2d)       !Analysis ensemble spread
    REAL(r_size) dxa3d(nmem,nv3d),dxa2d(nmem,nv2d) !Analysis ensemble perturbation

    REAL(r_size),ALLOCATABLE :: hdxfl(:,:)         !dYf=H(dXf)
    REAL(r_size),ALLOCATABLE :: rdiag(:)           !Obs. error variance
    REAL(r_size),ALLOCATABLE :: rloc(:)            !Localization weight function
    REAL(r_size),ALLOCATABLE :: dobl(:)            !Innovation

    REAL(r_size) wvec(nmem)      !w vector
    REAL(r_size) Wmat(nmem,nmem) !Transform matrix W

    !---INOUT
    REAL(r_size),INTENT(INOUT)  :: ens3d(nij1,nlev,nmem,nv3d),ens2d(nij1,nmem,nv2d) !Ensemble forecast => analysis

    IF(nobs == 0)THEN
       RETURN
    END IF

    k_min=min(k_sfc,k_btm)
    k_max=max(k_sfc,k_btm)

    !---Main Analysis Loop
    !   i (horizontal): parallel, k (vertical): sequential
    !   Horizontal local obs. are searched once per i and reused for all k.
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) &
    !$OMP & PRIVATE(i,k,iseed_grd,nobsh,nobsl,iobs_h,dist_h,iobs_l,dist_l) &
    !$OMP & PRIVATE(fcst2dm,fcst2ds,dxf2d,fcst3dm,fcst3ds,dxf3d,anal2dm,anal2ds,dxa2d,anal3dm,anal3ds,dxa3d) &
    !$OMP & PRIVATE(hdxfl,rdiag,rloc,dobl,wvec,Wmat)
    DO i=1,nij1

       !---Land: Analysis = Forecast (all levels)
       IF(maskt1(i) == REAL(0.d0,r_size)) CYCLE

       !---Horizontal local obs. at grid point i (common to all levels)
       CALL get_obs_local_horizontal(i,nobsh,iobs_h,dist_h)
       ALLOCATE(iobs_l(nobsh),dist_l(nobsh))

       DO k=k_min,k_max

          !---Local obs. at (i,k): vertical cutoff on horizontal local obs.
          CALL select_obs_local_vertical(nobsh,iobs_h,dist_h,depth1(i,k),nobsl,iobs_l,dist_l)

          !---LETKF/LPF/LPFGM
          IF(nobsl == 0)THEN

             IF(iswitch_da == 1)THEN
                CALL letkf_core_noobs(wvec,Wmat)
             ELSE IF(iswitch_da == 2 .or. iswitch_da == 3)THEN
                CALL lpfgm_core_noobs(wvec,Wmat)
             END IF

          ELSE

             ALLOCATE(hdxfl(nobsl,nmem))
             ALLOCATE(rdiag(nobsl),rloc(nobsl),dobl(nobsl))
             CALL set_obs_local(nobsl,iobs_l,dist_l,depth1(i,k), &
                  & hdxfl,rdiag,rloc,dobl)

             IF(iswitch_da == 1)THEN
                CALL letkf_core(nobsl,hdxfl,rdiag,rloc,dobl,wvec,Wmat)
             ELSE IF(iswitch_da == 2 .or. iswitch_da == 3)THEN
                IF(LRAND)THEN
                   iseed_grd = NINT(ri1(i)) + (NINT(rj1(i))-1)*nlon + (k-1)*nlon*nlat
                ELSE
                   iseed_grd = 0
                END IF
                CALL lpfgm_core(iseed_grd,nobsl,hdxfl,rdiag,rloc,dobl,wvec,Wmat)
             END IF

             DEALLOCATE(hdxfl,rdiag,rloc,dobl)

          END IF

          !---3D-variable update
          call Ensemble_Mean_Spread_Perturbation(nmem,nv3d,ens3d(i,k,:,:),fcst3dm,fcst3ds,dxf3d)
          CALL Analysis_Ensemble(nmem,nv3d,fcst3dm,dxf3d,wvec,Wmat,ens3d(i,k,:,:))
          CALL Ensemble_Mean_Spread_Perturbation(nmem,nv3d,ens3d(i,k,:,:),anal3dm,anal3ds,dxa3d)

          !---2D-variable update
          IF(k == k_sfc)THEN
             call Ensemble_Mean_Spread_Perturbation(nmem,nv2d,ens2d(i,:,:),fcst2dm,fcst2ds,dxf2d)
             CALL Analysis_Ensemble(nmem,nv2d,fcst2dm,dxf2d,wvec,Wmat,ens2d(i,:,:))
             CALL Ensemble_Mean_Spread_Perturbation(nmem,nv2d,ens2d(i,:,:),anal2dm,anal2ds,dxa2d)
          END IF

          !---RTPP
          IF(0.d0 < ALPHA_RTPP)THEN
             CALL RTPP(nmem,nv3d,dxf3d,dxa3d,anal3dm,ens3d(i,k,:,:))
             IF(k == k_sfc)THEN
                CALL RTPP(nmem,nv2d,dxf2d,dxa2d,anal2dm,ens2d(i,:,:))
             END IF
          END IF

          !---RTPS
          IF(0.d0 < ALPHA_RTPS)THEN
             CALL RTPS(nmem,nv3d,fcst3ds,anal3dm,anal3ds,dxa3d,ens3d(i,k,:,:))
             IF(k == k_sfc)THEN
                CALL RTPS(nmem,nv2d,fcst2ds,anal2dm,anal2ds,dxa2d,ens2d(i,:,:))
             END IF
          END IF

          !---Mask
          CALL APPLY_MASK(i,k,ens3d(i,k,:,:),ens2d(i,:,:))

       END DO !k

       DEALLOCATE(iobs_h,dist_h,iobs_l,dist_l)

    END DO !i
    !$OMP END PARALLEL DO

    !---Round off
    IF(ROFF)THEN
       CALL round_off(ens3d)
    END IF
    
  END SUBROUTINE da_main

  !=======================================================================
  ! Local observations
  !   1. get_obs_local_horizontal : horizontal cutoff (once per grid point i)
  !   2. select_obs_local_vertical: vertical cutoff   (each level k)
  !   3. set_obs_local            : dYf, R, localization, innovation
  !
  !   Local obs. keep the order of candidates (get_obsc_id),
  !   so the result is identical to the former count_nobsl/obs_local.
  !=======================================================================

  !-----------------------------------------------------------------------
  ! 1. Horizontal local observations at grid point i |
  !-----------------------------------------------------------------------

  SUBROUTINE get_obs_local_horizontal(i,nobsh,iobs_h,dist_h)

    USE common
    USE common_setting
    IMPLICIT NONE

    !---Common
    INTEGER imin,imax  !Longitude index (might be negative)
    INTEGER jmin,jmax  !Latitude index
    INTEGER nobsc      !Candidate obs. size in the search window
    INTEGER ioc,noc    !Candidate obs. index
    INTEGER iobs       !Global obs. index

    INTEGER,ALLOCATABLE :: obsc_id(:)

    REAL(r_size),ALLOCATABLE :: dist(:)

    !---IN
    INTEGER,INTENT(IN) :: i                              !Horizontal index

    !---OUT
    INTEGER,INTENT(OUT) :: nobsh                         !#Obs. within horizontal cutoff
    INTEGER,ALLOCATABLE,INTENT(OUT)      :: iobs_h(:)    !Global obs. index   [1:nobsh is valid]
    REAL(r_size),ALLOCATABLE,INTENT(OUT) :: dist_h(:)    !Horizontal distance [1:nobsh is valid]

    !---Candidate obs. in the search window
    CALL obs_search_range(i,imin,imax,jmin,jmax)
    CALL count_candidate_obs(imin,imax,jmin,jmax,nobsc)

    ALLOCATE(obsc_id(nobsc),dist(nobsc))
    CALL get_obsc_id(imin,imax,jmin,jmax,nobsc,noc,obsc_id)

    !---Horizontal cutoff (keep candidate order)
    !Count nobsh (distance is stored to avoid recomputation)
    nobsh = 0
    DO ioc=1,noc

       iobs = obsc_id(ioc)

       call distance_2p(obslon(iobs),obslat(iobs),lont1(i),latt1(i),dist(ioc))
       IF(dist_zero < dist(ioc)) CYCLE

       nobsh = nobsh + 1

    END DO

    !Get iobs_h and dist_h (zero-size arrays if nobsh == 0)
    ALLOCATE(iobs_h(nobsh),dist_h(nobsh))
    nobsh = 0

    DO ioc=1,noc

       IF(dist_zero < dist(ioc)) CYCLE

       nobsh = nobsh + 1
       iobs_h(nobsh) = obsc_id(ioc)
       dist_h(nobsh) = dist(ioc)

    END DO

    DEALLOCATE(obsc_id,dist)
    
  END SUBROUTINE get_obs_local_horizontal

  !-----------------------------------------------------------------------
  ! 2. Local observations at level k (vertical cutoff) |
  !-----------------------------------------------------------------------

  SUBROUTINE select_obs_local_vertical(nobsh,iobs_h,dist_h,dep_in,nobsl,iobs_l,dist_l)

    USE common_setting
    IMPLICIT NONE

    INTEGER iobsh
    REAL(r_size) dlev

    !---IN
    INTEGER,INTENT(IN)      :: nobsh
    INTEGER,INTENT(IN)      :: iobs_h(nobsh)  !Global obs. index of horizontal local obs.
    REAL(r_size),INTENT(IN) :: dist_h(nobsh)  !Horizontal distance
    REAL(r_size),INTENT(IN) :: dep_in         !Depth (negative sign) [m]

    !---OUT (allocated in da_main with size nobsh)
    INTEGER,INTENT(OUT)      :: nobsl          !Local obs. size
    INTEGER,INTENT(OUT)      :: iobs_l(nobsh)  !Global obs. index of local obs. [1:nobsl is valid]
    REAL(r_size),INTENT(OUT) :: dist_l(nobsh)  !Horizontal distance             [1:nobsl is valid]

    nobsl = 0

    DO iobsh=1,nobsh

       !Level difference
       dlev = ABS(obslev(iobs_h(iobsh)) - dep_in)

       !Compare with vertical cutoff scale
       IF(0.d0 < sigma_obsv .and. dist_zerov < dlev) CYCLE

       nobsl = nobsl + 1
       iobs_l(nobsl) = iobs_h(iobsh)
       dist_l(nobsl) = dist_h(iobsh)

    END DO

  END SUBROUTINE select_obs_local_vertical
  
  !-----------------------------------------------------------------------
  ! 3. Set local observation arrays |
  !-----------------------------------------------------------------------

  SUBROUTINE set_obs_local(nobsl,iobs_l,dist_l,dep_in,hdxfl,rdiag,rloc,dobl)  

    USE common_setting
    IMPLICIT NONE

    INTEGER iobsl  !Local obs. index
    INTEGER iobs   !Global obs. index
    REAL(r_size) dist,dlev

    !---IN
    INTEGER,INTENT(IN)      :: nobsl
    INTEGER,INTENT(IN)      :: iobs_l(:)  !Global obs. index of local obs.
    REAL(r_size),INTENT(IN) :: dist_l(:)  !Horizontal distance of local obs.
    REAL(r_size),INTENT(IN) :: dep_in     !Depth (negative sign) [m]

    !---OUT
    REAL(r_size),INTENT(OUT) :: hdxfl(nobsl,nmem) !dYf
    REAL(r_size),INTENT(OUT) :: rdiag(nobsl)      !Obs. error variance
    REAL(r_size),INTENT(OUT) :: rloc(nobsl)       !Localization function
    REAL(r_size),INTENT(OUT) :: dobl(nobsl)       !Innovation

    DO iobsl=1,nobsl

       iobs = iobs_l(iobsl)
       dist = dist_l(iobsl)
       dlev = ABS(obslev(iobs) - dep_in)

       !dYf(iobs,:): contiguous access (nmem elements) of shared hdxfT(nmem,nobs)
       hdxfl(iobsl,:) = hdxfT(:,iobs)

       dobl(iobsl) = dob(iobs)

       !Obs. error variance
       rdiag(iobsl) = obserr(iobs)**2

       !Localization function
       IF(sigma_obsv <= 0.d0)THEN !No vertical localization
          rloc(iobsl) = REAL(EXP(-0.5d0 * (dist/sigma_obs)**2), r_size)
       ELSE
          rloc(iobsl) = REAL(EXP(-0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2)), r_size)
       END IF

    END DO

  END SUBROUTINE set_obs_local

  !------------------------------------------------------------------------
  ! Search longitude and latitude grid index range at grid point i |
  !------------------------------------------------------------------------
  ! NOTE:
  ! imin,imax are NOT clipped here (wrapped/clipped in obs_idx_range)
  !------------------------------------------------------------------------

  SUBROUTINE obs_search_range(i,imin,imax,jmin,jmax)

    USE common_setting
    IMPLICIT NONE

    !---Parameter
    REAL(r_dble),PARAMETER :: d2r=pi/180.d0

    !---Common
    REAL(r_dble) :: latp,rfac

    !---IN
    INTEGER,INTENT(IN) :: i

    !---OUT
    INTEGER,INTENT(OUT) :: imin,imax
    INTEGER,INTENT(OUT) :: jmin,jmax

    !Poleward latitude edge within the search window [deg]
    latp = MIN( ABS(REAL(latt1(i),r_dble)) + REAL(dist_zero,r_dble)/re/d2r, 89.d0 )

    !dx(poleward edge)/dx(grid point)  
    rfac = COS(d2r*latp) / COS(d2r*REAL(latt1(i),r_dble))   

    !Min/Max i and j index
    imin = FLOOR  ( ri1(i) - dist_zero/(dlont1(i)*rfac) ) - 1
    imax = CEILING( ri1(i) + dist_zero/(dlont1(i)*rfac) ) + 1
    jmin = MAX( FLOOR  ( rj1(i) - dist_zero/dlatt1(i) ) - 1, 1    )
    jmax = MIN( CEILING( rj1(i) + dist_zero/dlatt1(i) ) + 1, nlat )

  END SUBROUTINE obs_search_range

  !-----------------------------------------------------------------
  ! Count Candidate observation size (nobsc) in the search window |
  !-----------------------------------------------------------------
  ! NOTE:
  ! nobsc is not completely same as the local observation size. 
  !-----------------------------------------------------------------

  SUBROUTINE count_candidate_obs(imin,imax,jmin,jmax,nobsc)

    USE common_setting
    IMPLICIT NONE

    !---Common
    INTEGER j
    INTEGER icase,ncase
    INTEGER idx_s(3),idx_e(3)
    INTEGER ido_s,ido_e

    !---IN
    INTEGER,INTENT(IN) :: imin,imax
    INTEGER,INTENT(IN) :: jmin,jmax

    !---OUT
    INTEGER,INTENT(OUT) :: nobsc

    CALL obs_idx_range(imin,imax,ncase,idx_s,idx_e)

    nobsc = 0
    DO j=jmin,jmax
       DO icase=1,ncase
          CALL obsid_range(idx_s(icase),idx_e(icase),j,ido_s,ido_e)
          nobsc = nobsc + MAX(ido_e-ido_s+1,0)
       END DO
    END DO

  END SUBROUTINE count_candidate_obs

  !-------------------------------------------------------------------
  ! Longitude grid ID range in the search window |
  !-------------------------------------------------------------------
  ! imin,imax:    Search range (not clipped) from obs_search_range
  ! lqglobal:     Quasi-global cyclic domain.
  !               idx and idx+(nlon-2) are the same longitudes.
  !               nlon-1 is the same as ix=1 (halo).
  ! idx_s and _e: Start and End idx points
  ! ncase:        Number of cases with different idx_s and idx_e
  !-------------------------------------------------------------------

  SUBROUTINE obs_idx_range(imin,imax,ncase,idx_s,idx_e)

    USE common_setting, only: nlon, lqglobal
    IMPLICIT NONE

    !---Common
    INTEGER nx
    INTEGER imin_tmp,imax_tmp

    !---IN
    INTEGER,INTENT(IN) :: imin,imax

    !---OUT
    INTEGER,INTENT(OUT) :: ncase
    INTEGER,INTENT(OUT) :: idx_s(3),idx_e(3)

    ncase=0
    idx_s(:)=0
    idx_e(:)=0

    IF(lqglobal)THEN

       !Number of independent longitude grid points     
       nx = nlon-2

       IF(nx <= imax-imin+1)THEN !Search whole longitude direction
          ncase=1
          idx_s(1)=1
          idx_e(1)=nlon-1
          RETURN
       END IF

       imin_tmp = MODULO(imin-1,nx)+1
       imax_tmp = MODULO(imax-1,nx)+1

       IF(imin_tmp <= imax_tmp)THEN !General case
          ncase=1
          idx_s(1)=imin_tmp
          idx_e(1)=imax_tmp
          IF(imin_tmp == 1)THEN !(ilon=1) = (ilon=nlon-1)
             ncase=2
             idx_s(2)=nlon-1
             idx_e(2)=nlon-1
          END IF
       ELSE !Across the periodic boundary
          ncase=3

          idx_s(1)=imin_tmp
          idx_e(1)=nx

          idx_s(2)=1
          idx_e(2)=imax_tmp

          idx_s(3)=nlon-1
          idx_e(3)=nlon-1
       END IF

    ELSE

       imin_tmp = MAX(imin,1)
       imax_tmp = MIN(imax,nlon)
       IF(imin_tmp <= imax_tmp)THEN
          ncase=1
          idx_s(1)=imin_tmp
          idx_e(1)=imax_tmp
       END IF

    END IF

  END SUBROUTINE obs_idx_range

  !-----------------------------------------------------------------------------------
  ! Obs index range [ido_s,ido_e] in longitude range [idx_s,idx_e] at jth latitude |
  !-----------------------------------------------------------------------------------
  ! NOTE:
  ! nobsgrd(i,j): Cumulative obs size up to grid (i,j)
  !-----------------------------------------------------------------------------------

  SUBROUTINE obsid_range(idx_s,idx_e,j,ido_s,ido_e)

    USE common_setting, only: nlon, nobsgrd
    IMPLICIT NONE

    !---IN
    INTEGER,INTENT(IN) :: idx_s,idx_e
    INTEGER,INTENT(IN) :: j

    !---OUT
    INTEGER,INTENT(OUT) :: ido_s,ido_e

    !---ido_s
    IF(1 < idx_s)THEN
       ido_s = nobsgrd(idx_s-1,j)+1
    ELSE IF(1 < j)THEN
       ido_s = nobsgrd(nlon,j-1)+1
    ELSE
       ido_s = 1
    END IF

    !---ido_e
    ido_e = nobsgrd(idx_e,j)

  END SUBROUTINE obsid_range


  !----------------------------------------------------
  ! Get obs candidate indices in the search window |
  !----------------------------------------------------
  !
  ! noc: Candidate observation array size
  ! obsc_id: Candidate observation index 
  !
  !----------------------------------------------------

  SUBROUTINE get_obsc_id(imin,imax,jmin,jmax,nobsc,noc,obsc_id)

    USE common_setting
    IMPLICIT NONE

    !---Common
    INTEGER j
    INTEGER icase,ncase
    INTEGER idx_s(3),idx_e(3)
    INTEGER ido,ido_s,ido_e

    !---IN
    INTEGER,INTENT(IN) :: imin,imax
    INTEGER,INTENT(IN) :: jmin,jmax
    INTEGER,INTENT(IN) :: nobsc

    !---OUT
    INTEGER,INTENT(OUT) :: noc
    INTEGER,INTENT(OUT) :: obsc_id(nobsc)

    CALL obs_idx_range(imin,imax,ncase,idx_s,idx_e)

    noc = 0
    DO j=jmin,jmax
       DO icase=1,ncase

          CALL obsid_range(idx_s(icase),idx_e(icase),j,ido_s,ido_e)

          DO ido=ido_s,ido_e
             noc = noc + 1
             obsc_id(noc) = ido
          END DO
       END DO
    END DO

  END SUBROUTINE get_obsc_id

  !-----------------------------------------------------------------------
  ! Ensemble Mean/Spread/Perturbation |
  !-----------------------------------------------------------------------

  SUBROUTINE Ensemble_Mean_Spread_Perturbation(nmem,nvd,x,xmean,xsprd,dx)

    USE common_setting, only: r_size
    implicit none

    !---Common
    INTEGER imem,ivd

    !---IN
    INTEGER,INTENT(IN) ::  nmem,nvd
    REAL(r_size),INTENT(IN) :: x(nmem,nvd)

    !---OUT
    REAL(r_size),INTENT(OUT) :: xmean(nvd),xsprd(nvd),dx(nmem,nvd)

    !xmean
    DO ivd=1,nvd

       xmean(ivd) = REAL(0.d0,r_size)

       DO imem=1,nmem
          xmean(ivd)=xmean(ivd)+x(imem,ivd)
       END DO

       xmean(ivd)=xmean(ivd) / REAL(nmem,r_size)

    END DO

    !xspread
    DO ivd=1,nvd

       xsprd(ivd) = REAL(0.d0,r_size)

       DO imem=1,nmem
          xsprd(ivd)=xsprd(ivd)+(x(imem,ivd)-xmean(ivd))*(x(imem,ivd)-xmean(ivd))
       END DO

       xsprd(ivd)=sqrt(xsprd(ivd) / REAL(nmem-1, r_size))
    END DO

    !dx
    DO ivd=1,nvd
       DO imem=1,nmem
          dx(imem,ivd)=x(imem,ivd)-xmean(ivd)
       END DO
    END DO

  END SUBROUTINE Ensemble_Mean_Spread_Perturbation

  !------------------------------------------------------------------
  ! Analysis Ensemble |
  !------------------------------------------------------------------

  SUBROUTINE Analysis_Ensemble(nmem,nvd,fcstm,dxf,wvec,Wmat,anal)

    USE common_setting, only: r_size
    IMPLICIT NONE

    !---Common
    INTEGER i,j,ivd

    !---IN
    INTEGER,INTENT(IN) :: nmem,nvd
    REAL(r_size),INTENT(IN) :: fcstm(nvd),dxf(nmem,nvd)
    REAL(r_size),INTENT(IN) :: wvec(nmem),Wmat(nmem,nmem)

    !---OUT
    REAL(r_size),INTENT(OUT) :: anal(nmem,nvd)

    DO ivd=1,nvd
       DO j=1,nmem
          anal(j,ivd)=fcstm(ivd)
          DO i=1,nmem
             anal(j,ivd)=anal(j,ivd)+dxf(i,ivd)*(wvec(i)+Wmat(i,j))
          END DO
       END DO
    END DO

  END SUBROUTINE Analysis_Ensemble

  !----------------------------------------------------------------
  ! RTPP (Zhang et al. 2004; Kotsuki et al. 2017)|
  !----------------------------------------------------------------

  SUBROUTINE RTPP(nmem,nvd,dxf,dxa,analm,anal)

    USE common_setting,only :r_size, ALPHA_RTPP
    IMPLICIT NONE

    !---Common
    INTEGER imem,ivd 

    !---IN
    INTEGER,INTENT(IN) :: nmem,nvd

    REAL(r_size),INTENT(IN) :: dxf(nmem,nvd),dxa(nmem,nvd)
    REAL(r_size),INTENT(IN) :: analm(nvd)

    !---INOUT
    REAL(r_size),INTENT(INOUT) :: anal(nmem,nvd)

    DO ivd=1,nvd
       DO imem=1,nmem
          anal(imem,ivd)=analm(ivd)+ALPHA_RTPP*dxf(imem,ivd)+(REAL(1.d0,r_size)-ALPHA_RTPP)*dxa(imem,ivd)
       END DO
    END DO

  END SUBROUTINE RTPP

  !----------------------------------------------------------------
  ! RTPS (Whitaker and Hamill 2012; Kotsuki et al. 2017)|
  !----------------------------------------------------------------

  SUBROUTINE RTPS(nmem,nvd,fcsts,analm,anals,dxa,anal)

    USE common_setting,only :r_size, ALPHA_RTPS
    IMPLICIT NONE

    !---Common
    INTEGER imem,ivd 

    REAL(r_size) rtps_factor

    !---IN
    INTEGER,INTENT(IN) :: nmem,nvd

    REAL(r_size),INTENT(IN) :: fcsts(nvd)  
    REAL(r_size),INTENT(IN) :: analm(nvd),anals(nvd)  
    REAL(r_size),INTENT(IN) :: dxa(nmem,nvd)

    !---INOUT
    REAL(r_size),INTENT(INOUT) :: anal(nmem,nvd)

    DO ivd=1,nvd

       IF(anals(ivd) == REAL(0.d0, r_size))THEN
          rtps_factor=REAL(0.d0, r_size)
       ELSE
          rtps_factor=ALPHA_RTPS*fcsts(ivd)/anals(ivd)+(REAL(1.d0,r_size)-ALPHA_RTPS)
       END IF

       DO imem=1,nmem
          anal(imem,ivd)=analm(ivd)+dxa(imem,ivd)*rtps_factor
       END DO
    END DO

  END SUBROUTINE RTPS

  !----------------------------------------------------------------
  ! Mask |
  !----------------------------------------------------------------

  SUBROUTINE APPLY_MASK(i,k,dat3d,dat2d)

    USE common_setting
    IMPLICIT NONE

    !---Common
    INTEGER iv3d,iv2d
    REAL(r_size) mask

    !---IN
    INTEGER,INTENT(IN) :: i,k

    !---INOUT
    REAL(r_size),INTENT(INOUT) :: dat3d(nmem,nv3d),dat2d(nmem,nv2d)

    DO iv3d=1,nv3d

       SELECT CASE(iv3d)
       CASE(iv3d_u)
          mask=masku1(i)
       CASE(iv3d_v)
          mask=maskv1(i)
       CASE DEFAULT
          mask=maskt1(i)
       END SELECT

       dat3d(:,iv3d)=dat3d(:,iv3d)*mask

    END DO

    IF(k == k_sfc)THEN

       DO iv2d=1,nv2d

          SELECT CASE(iv2d)
          CASE(iv2d_ubar)
             mask=masku1(i)
          CASE(iv2d_vbar)
             mask=maskv1(i)
          CASE DEFAULT
             mask=maskt1(i)
          END SELECT

          dat2d(:,iv2d)=dat2d(:,iv2d)*mask

       END DO

    END IF

  END SUBROUTINE APPLY_MASK

  !----------------------------------------------------------------
  ! Round off |
  !----------------------------------------------------------------

  SUBROUTINE round_off(anal3d)

    USE common_setting
    IMPLICIT NONE

    !---Common
    INTEGER i,k,imem
    INTEGER k_min,k_max

    !---INOUT
    REAL(r_size),INTENT(INOUT) :: anal3d(nij1,nlev,nmem,nv3d)

    k_min=MIN(k_sfc,k_btm)
    k_max=MAX(k_sfc,k_btm)

    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(imem,k,i)
    DO imem=1,nmem
       DO k=k_min,k_max
          DO i=1,nij1

             IF(maskt1(i) == REAL(0.d0, r_size)) CYCLE
             anal3d(i,k,imem,iv3d_t) = MAX(anal3d(i,k,imem,iv3d_t), tmin)
             anal3d(i,k,imem,iv3d_t) = MIN(anal3d(i,k,imem,iv3d_t), tmax)
             anal3d(i,k,imem,iv3d_s) = MAX(anal3d(i,k,imem,iv3d_s), smin)
             anal3d(i,k,imem,iv3d_s) = MIN(anal3d(i,k,imem,iv3d_s), smax)           

          END DO
       END DO
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE round_off

END MODULE common_da
