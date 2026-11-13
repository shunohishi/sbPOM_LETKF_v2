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
  !   2026/10/21 Shun Ohishi       updated by using Claude
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
    USE common_mpi
    USE common_setting
    USE MPI
    IMPLICIT NONE

    INTEGER i,k,k_min,k_max
    INTEGER nobsl
    INTEGER iseed_grd !Seed for random number in LPF/LPFGM (0: same at all grid points)

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

    !---Timer (per thread, summed)
    REAL(r_dble) t0,t1,t2
    REAL(r_dble) t_search,t_core,t_update
    REAL(r_dble) time_thread(3)  !Thread-averaged time of search/core/update (myrank)
    REAL(r_dble) time_max(3)     !Max over ranks
    INTEGER ierr
    
    !---INOUT
    REAL(r_size),INTENT(INOUT)  :: ens3d(nij1,nlev,nmem,nv3d),ens2d(nij1,nmem,nv2d) !Ensemble forecast => analysis

    !---OUT
    !REAL(r_size),INTENT(OUT) :: anal3d(nij1,nlev,nmem,nv3d),anal2d(nij1,nmem,nv2d) !Ensemble analysis
    
    IF(nobs == 0)THEN
       RETURN
    END IF

    t_search=0.d0
    t_core=0.d0
    t_update=0.d0
    
    k_min=min(k_sfc,k_btm)
    k_max=max(k_sfc,k_btm)
    !---Main Analysis Loop
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) COLLAPSE(2) PRIVATE(i,k,iseed_grd,fcst2dm,fcst2ds,dxf2d,fcst3dm,fcst3ds,dxf3d,anal2dm,anal2ds,dxa2d,anal3dm,anal3ds,dxa3d,hdxfl,rdiag,rloc,dobl,nobsl,wvec,Wmat,t0,t1,t2) &
    !$OMP & REDUCTION(+:t_search,t_core,t_update)
    DO k=1,nlev
       DO i=1,nij1

          !---Land or out of [k_min,k_max]: Analysis = Forecast
          IF(maskt1(i) == REAL(0.d0,r_size) .or. k < k_min .or. k_max < k)THEN
             CYCLE
          END IF
          
          !---Count local observation
          t0 = OMP_GET_WTIME()
          CALL count_nobsl(i,k,depth1(i,k),nobsl)

          !---LETKF/LPF/LPFGM
          IF(nobsl == 0)THEN
             
             t1 = OMP_GET_WTIME()             

             IF(iswitch_da == 1)THEN
                CALL letkf_core_noobs(wvec,Wmat)
             ELSE IF(iswitch_da == 2 .or. iswitch_da == 3)THEN
                CALL lpfgm_core_noobs(wvec,Wmat)
             END IF
             
          ELSE

             ALLOCATE(hdxfl(nobsl,nmem))
             ALLOCATE(rdiag(nobsl),rloc(nobsl),dobl(nobsl))

             CALL obs_local(i,k,nobsl,depth1(i,k), &
                  & hdxfl,rdiag,rloc,dobl)

             t1 = OMP_GET_WTIME()
             
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

          t2 = OMP_GET_WTIME()
          
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

          !---Timer
          t_search = t_search + (t1-t0)
          t_core   = t_core   + (t2-t1)
          t_update = t_update + (OMP_GET_WTIME()-t2)
          
       END DO !i
    END DO !k
    !$OMP END PARALLEL DO

    !---Round off
    IF(ROFF)THEN
       CALL round_off(ens3d)
    END IF

    !---Max over ranks (thread-summed time / #threads ~ wall time)
    time_thread = (/t_search,t_core,t_update/)/REAL(OMP_GET_MAX_THREADS(),r_dble)
    CALL MPI_REDUCE(time_thread,time_max,3,MPI_DOUBLE_PRECISION,MPI_MAX,root,MPI_COMM_WORLD,ierr)
    IF(myrank == root) WRITE(6,'(A,3F10.2)') "da_main [s] search/core/update (max rank): ",time_max
    
  END SUBROUTINE da_main

END MODULE common_da

!=======================================================================
! Project global observations to local
!     (hdxf_g,dob_g,rdiag_g) -> (hdxf,dob,rdiag)
!=======================================================================

!-----------------------------------------------------------------------
! Count local observation size |
!-----------------------------------------------------------------------

SUBROUTINE count_nobsl(i,k,dep_in,nobsl)

  USE common  
  USE common_setting
  IMPLICIT NONE

  INTEGER imin,imax !Longitude index (might be negative)
  INTEGER jmin,jmax !Latitude index
  INTEGER nobsc     !Candidate observation size (~ Rough local observation) 
  INTEGER ioc,noc   !Candidate observation index
  INTEGER iobs

  INTEGER,ALLOCATABLE :: obsc_id(:)

  REAL(r_size) dist,dlev

  !---IN
  INTEGER,INTENT(IN) :: i           !Horizontal index
  INTEGER,INTENT(IN) :: k           !Vertical index
  REAL(r_size),INTENT(IN) :: dep_in !Depth (negative sign) [m]

  !---OUT
  INTEGER,INTENT(OUT) :: nobsl    !Local observation size

  nobsl=0

  !--- Search range: imin,imax,jmin,jmax
  !    (imin,imax are clipped/wrapped in obs_idx_range)
  CALL obs_search_range(i,imin,imax,jmin,jmax)

  !--- Observation size in the search area ==> Rough local observation size: nobsc
  CALL count_candidate_obs(imin,imax,jmin,jmax,nobsc)
  IF(nobsc == 0)THEN !No candidate
     RETURN
  END IF

  !--- Get candidate observation ID
  ALLOCATE(obsc_id(nobsc))
  CALL get_obsc_id(imin,imax,jmin,jmax,nobsc,noc,obsc_id)

  !---Local observation size
  DO ioc=1,noc

     !iobs
     iobs=obsc_id(ioc)

     !Level difference
     dlev = ABS(obslev(iobs) - dep_in)

     !Compare with vertical cutoff scale
     IF(0.d0 < sigma_obsv .and. dist_zerov < dlev) CYCLE

     !Horizontal distance between analysis grid and observation points
     call distance_2p(obslon(iobs),obslat(iobs),lont1(i),latt1(i),dist)

     !Compare with horizontal cutoff scale
     IF(dist_zero < dist) CYCLE

     !Count up local observation size
     nobsl = nobsl + 1

  END DO

  DEALLOCATE(obsc_id)

END SUBROUTINE count_nobsl

!--------------------------------------------------------------------
! Local observation |
!--------------------------------------------------------------------

SUBROUTINE obs_local(i,k,nobsl,dep_in,hdxfl,rdiag,rloc,dobl)

  USE common
  USE common_mpi
  USE common_setting
  USE MPI
  IMPLICIT NONE

  INTEGER imin,imax !Longitude index (might be negative)
  INTEGER jmin,jmax !Latitude index
  INTEGER nobsc     !Candidate observation size (~ Rough local observation)
  INTEGER ioc,noc   !Candidate observation index
  INTEGER iobs      !Global observation index
  INTEGER iobsl     !Local observation index
  INTEGER(8) ihdxf  !1-D index of hdxf(iobs,imem) in hdxf_flat
  INTEGER ierr

  INTEGER,ALLOCATABLE :: obsc_id(:)

  REAL(r_size) dist,dlev

  !---IN
  INTEGER,INTENT(IN) :: i         !Horizontal index
  INTEGER,INTENT(IN) :: k         !Vertical index
  INTEGER,INTENT(IN) :: nobsl     !Local observation size (from count_nobsl)
  REAL(r_size),INTENT(IN) :: dep_in !Depth (negative sign) [m]

  !---OUT
  REAL(r_size),INTENT(OUT) :: hdxfl(nobsl,nmem) !dYf
  REAL(r_size),INTENT(OUT) :: rdiag(nobsl)     !Obs. error variance
  REAL(r_size),INTENT(OUT) :: rloc(nobsl)      !Localization function
  REAL(r_size),INTENT(OUT) :: dobl(nobsl)       !Innovation

  !--- Search range: imin,imax,jmin,jmax
  !    (imin,imax are clipped/wrapped in obs_idx_range)
  CALL obs_search_range(i,imin,imax,jmin,jmax)

  !--- Observation size in the search area ==> Rough local observation size: nobsc
  CALL count_candidate_obs(imin,imax,jmin,jmax,nobsc)

  !--- Get candidate observation ID
  ALLOCATE(obsc_id(nobsc))
  CALL get_obsc_id(imin,imax,jmin,jmax,nobsc,noc,obsc_id)

  !---Local observation
  iobsl=0
  DO ioc=1,noc

     !iobs
     iobs=obsc_id(ioc)
     
     !Level difference (same as count_nobsl)
     dlev = ABS(obslev(iobs) - dep_in)

     !Compare with vertical cutoff scale
     IF(0.d0 < sigma_obsv .and. dist_zerov < dlev) CYCLE

     !Horizontal distance between analysis grid and observation points
     call distance_2p(obslon(iobs),obslat(iobs),lont1(i),latt1(i),dist)

     !Compare with horizontal cutoff scale
     IF(dist_zero < dist) CYCLE

     !Count up local observation index
     iobsl = iobsl + 1

     IF(nobsl < iobsl) EXIT

     !Global --> Local
     !dYf(iobs,:): contiguous access (nmem elements) of shared hdxfT(nmem,nobs)
     !ihdxf = INT(iobs-1,8)*INT(nmem,8)            !Offset of dYf(iobs,:) in hdxfT_flat
     !hdxfl(iobsl,:) = hdxfT_flat(ihdxf+1:ihdxf+nmem)
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

  DEALLOCATE(obsc_id)

  !---Check consistency with count_nobsl
  IF(iobsl /= nobsl)THEN
     WRITE(6,'(A,I10,A,I10)') "***FATAL ERROR: Local #OBS=",iobsl," /= #OBS from count_nobsl=",nobsl
     WRITE(6,'(A,2I10)') "IJ,K=",i,k
     CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
  END IF

END SUBROUTINE obs_local

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
