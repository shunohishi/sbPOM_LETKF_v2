MODULE common_da
  !=======================================================================
  !
  ! [PURPOSE:] Module for DA (LETKF/LPF/LPFGM) with POM
  !
  ! [HISTORY:]
  !   01/26/2009 Takemasa Miyoshi  created
  !   02/03/2009 Takemasa Miyoshi  modified for ROMS
  !   01/26/2011 Yasumasa Miyazawa modified for POM (check 'pom' or 'POM')
  !   07/31/2025 Shun Ohishi       added OpenMP
  !   10/16/2025 Shun Ohishi       changed precision
  !   10/26/2025 Shun Ohishi       integrated LETKF/LPF/LPFGM
  !=======================================================================

CONTAINS

  !-----------------------------------------------------------------------
  ! Data Assimilation System | Main loop
  !-----------------------------------------------------------------------

  SUBROUTINE das_main(fcst3d,fcst2d,anal3d,anal2d)

    !$USE OMP_LIB
    USE common_setting
    USE common
    USE common_letkf
    USE common_lpfgm
    USE common_mpi
    IMPLICIT NONE

    INTEGER i,k
    INTEGER nobsl

    REAL(r_size) fcst3dm(nv3d),fcst2dm(nv2d)     !Forecast ensemble mean
    REAL(r_size) fcst3ds(nv3d),fcst2ds(nv2d)     !Forecast ensemble spread
    REAL(r_size) dxf3d(nbv,nv3d),dxf2d(nbv,nv2d) !Forecast ensemble perturbation

    REAL(r_size) anal3dm(nv3d),anal2dm(nv2d)     !Forecast ensemble mean
    REAL(r_size) anal3ds(nv3d),anal2ds(nv2d)     !Forecast ensemble spread
    REAL(r_size) dxa3d(nbv,nv3d),dxa2d(nbv,nv2d) !Forecast ensemble perturbation
    
    REAL(r_size),ALLOCATABLE :: hdxf(:,:) !dYf=H(dXf)
    REAL(r_size),ALLOCATABLE :: rdiag(:)  !Obs. error variance
    REAL(r_size),ALLOCATABLE :: rloc(:)   !Localization weight function
    REAL(r_size),ALLOCATABLE :: dep(:)    !Innovation

    REAL(r_size) wvec(nbv)     !w vector
    REAL(r_size) Wmat(nbv,nbv) !Transform matrix W

    !---IN
    REAL(r_size),INTENT(IN) :: fcst3d(nij1,nlev,nbv,nv3d),fcst2d(nij1,nbv,nv2d) !Ensemble forecast
    
    !---OUT
    REAL(r_size),INTENT(OUT) :: anal3d(nij1,nlev,nbv,nv3d),anal2d(nij1,nbv,nv2d) !Ensemble analysis

    WRITE(file_unit,'(A,I8)') "Target observation numbers : NOBS=",nobs

    IF(nobs == 0)THEN
       anal3d(:,:,:,:)=fcst3d(:,:,:,:)
       anal2d(:,:,:)=fcst2d(:,:,:)
       RETURN
    END IF

    anal3d(:,:,:,:)=0.d0
    anal2d(:,:,:)=0.d0
    !---Main Analysis Loop
    !$OMP PARALLEL DO PRIVATE(i,k,fcst2dm,fcst2ds,dxf2d,fcst3dm,fcst3ds,dxf3d,anal2dm,anal2ds,dxa2d,anal3dm,anal3ds,dxa3d,hdxf,rdiag,rloc,dep,nobsl,wvec,Wmat)
    DO k=1,nlev
       DO i=1,nij1

          !---Count local observation
          CALL count_nobsl(i,k,depth1(i,k),nobsl)

          !---LETKF/LPF/LPFGM
          IF(nobsl == 0)THEN
             IF(iswitch_da == 1) CALL letkf_core_noobs(wvec,Wmat)
             IF(iswitch_da == 2 .or. iswitch_da == 3) CALL lpfgm_core_noobs(wvec,Wmat)
          ELSE
             ALLOCATE(hdxf(nobsl,nbv))
             ALLOCATE(rdiag(nobsl),rloc(nobsl),dep(nobsl))
             CALL obs_local(i,k,nobsl,depth1(i,k), &
                  & hdxf,rdiag,rloc,dep)
             IF(iswitch_da == 1) CALL letkf_core(nobsl,hdxf,rdiag,rloc,dep,wvec,Wmat)
             IF(iswitch_da == 2 .or. iswitch_da == 3) CALL lpfgm_core(nobsl,hdxf,rdiag,rloc,dep,wvec,Wmat)
             DEALLOCATE(hdxf,rdiag,rloc,dep)
          END IF

          !---3D-variable update
          call Ensemble_Mean_Spread_Perturbation(nbv,nv3d,fcst3d(i,k,:,:),fcst3dm,fcst3ds,dxf3d)
          CALL Analysis_Ensemble(nbv,nv3d,fcst3dm,dxf3d,wvec,Wmat,anal3d(i,k,:,:))
          CALL Ensemble_Mean_Spread_Perturbation(nbv,nv3d,anal3d(i,k,:,:),anal3dm,anal3ds,dxa3d)

          !---2D-variable update
          IF(k == nlev)THEN
             call Ensemble_Mean_Spread_Perturbation(nbv,nv2d,fcst2d(i,:,:),fcst2dm,fcst2ds,dxf2d)
             CALL Analysis_Ensemble(nbv,nv2d,fcst2dm,dxf2d,wvec,Wmat,anal2d(i,:,:))
             CALL Ensemble_Mean_Spread_Perturbation(nbv,nv2d,anal2d(i,:,:),anal2dm,anal2ds,dxa2d)
          END IF
          
          !---RTPP
          IF(0.d0 < ALPHA_RTPP)THEN
             CALL RTPP(nbv,nv3d,dxf3d,dxa3d,anal3dm,anal3d(i,k,:,:))             
             IF(k == nlev)THEN
                CALL RTPP(nbv,nv2d,dxf2d,dxa2d,anal2dm,anal2d(i,:,:))
             END IF
          END IF

          !---RTPS
          IF(0.d0 < ALPHA_RTPS)THEN
             CALL RTPS(nbv,nv3d,fcst3ds,anal3dm,anal3ds,dxa3d,anal3d(i,k,:,:))
             IF(k == nlev)THEN
                CALL RTPS(nbv,nv2d,fcst2ds,anal2dm,anal2ds,dxa2d,anal2d(i,:,:))
             END IF
          END IF
                    
       END DO !i
    END DO !k
    !$OMP END PARALLEL DO
    
    !---Round off
    IF(ROFF)THEN
       CALL round_off(anal3d)
    END IF
    
  END SUBROUTINE das_main

END MODULE common_da

!-----------------------------------------------------------------------
! Project global observations to local
!     (hdxf_g,dep_g,rdiag_g) -> (hdxf,dep,rdiag)
!-----------------------------------------------------------------------

SUBROUTINE count_nobsl(i,k,rdep,nobsl)

  USE common_setting
  USE common_mpi
  USE common
  USE MPI
  IMPLICIT NONE

  INTEGER imin,imax,jmin,jmax
  INTEGER n,nn

  INTEGER,ALLOCATABLE:: nobs_use(:)

  REAL(r_size) :: dist,dlev

  !---IN
  INTEGER,INTENT(IN) :: i,k
  REAL(r_size),INTENT(IN) :: rdep

  !---OUT
  INTEGER,INTENT(OUT) :: nobsl               !Local #OBS

  !--- INITIALIZE
  ALLOCATE(nobs_use(nobs))
  
  !--- Search range: imin,imax,jmin,jmax
  imin = MAX(NINT(ri1(i) - dist_zero/dlon1(i)),1)
  imax = MIN(NINT(ri1(i) + dist_zero/dlon1(i)),nlon)
  jmin = MAX(NINT(rj1(i) - dist_zero/dlat1(i)),1)
  jmax = MIN(NINT(rj1(i) + dist_zero/dlat1(i)),nlat)

  CALL obs_local_sub(imin,imax,jmin,jmax, &
       & nn,nobs_use)

  IF(nn == 0)THEN
     nobsl=0
     DEALLOCATE(nobs_use)
     RETURN
  END IF

  !---Local obs.
  nobsl=0
  DO n=1,nn

     !---dlev
     IF(obselm(nobs_use(n)) == id_z_obs .AND. k < nlev)THEN
        dlev = ABS(rdep)
     ELSE IF(obselm(nobs_use(n)) /= id_z_obs) THEN
        dlev = ABS(obslev(nobs_use(n)) - rdep)
     ELSE
        dlev = 0.0d0
     END IF

     IF(dlev > dist_zerov)CYCLE

     !---dist
     call distance_2p(obslon(nobs_use(n)),obslat(nobs_use(n)),lon1(i),lat1(i),dist)

     IF(dist > dist_zero)CYCLE

     !---Global --> Local
     nobsl = nobsl + 1

  END DO !n

  IF(nobsl > nobs)THEN
     WRITE(6,'(A,I5,A,I5)') "***FATAL ERROR: Local #OBS=",nobsl," > Total #OBS=",nobs
     WRITE(6,'(A,2I10)') 'IJ,NN=',i,nn
     CALL finalize_mpi
     STOP
  END IF
  
  DEALLOCATE(nobs_use)

END SUBROUTINE count_nobsl

!---------------------------------------------------------

SUBROUTINE obs_local(i,k,nobsl,rdep,hdxf,rdiag,rloc,dep)

  USE common_setting
  USE common_mpi
  USE common
  USE MPI
  IMPLICIT NONE

  INTEGER imin,imax,jmin,jmax
  INTEGER n,nn
  INTEGER iobsl
  
  INTEGER,ALLOCATABLE:: nobs_use(:)

  REAL(r_size) :: dist,dlev

  !---IN
  INTEGER,INTENT(IN) :: i,k,nobsl
  REAL(r_size),INTENT(IN) :: rdep

  !---OUT
  REAL(r_size),INTENT(OUT) :: hdxf(nobsl,nbv) !dYf
  REAL(r_size),INTENT(OUT) :: rdiag(nobsl)    !Obs. error variance
  REAL(r_size),INTENT(OUT) :: rloc(nobsl)     !Localization function
  REAL(r_size),INTENT(OUT) :: dep(nobsl)      !Innovation

  !--- INITIALIZE
  ALLOCATE(nobs_use(nobs))
  
  !--- Search range: imin,imax,jmin,jmax
  imin = MAX(NINT(ri1(i) - dist_zero/dlon1(i)),1)
  imax = MIN(NINT(ri1(i) + dist_zero/dlon1(i)),nlon)
  jmin = MAX(NINT(rj1(i) - dist_zero/dlat1(i)),1)
  jmax = MIN(NINT(rj1(i) + dist_zero/dlat1(i)),nlat)

  CALL obs_local_sub(imin,imax,jmin,jmax, &
       & nn,nobs_use)

  !---Local obs.
  iobsl=0
  DO n=1,nn

     !---dlev
     IF(obselm(nobs_use(n)) == id_z_obs .AND. k < nlev)THEN
        dlev = ABS(rdep)
     ELSE IF(obselm(nobs_use(n)) /= id_z_obs) THEN
        dlev = ABS(obslev(nobs_use(n)) - rdep)
     ELSE
        dlev = 0.0d0
     END IF

     IF(dlev > dist_zerov)CYCLE

     !---dist
     call distance_2p(obslon(nobs_use(n)),obslat(nobs_use(n)),lon1(i),lat1(i),dist)

     IF(dist > dist_zero)CYCLE

     !---Global --> Local
     iobsl = iobsl + 1

     IF(nobsl < iobsl)THEN
        WRITE(6,'(A,I5,A,I5)') "***FATAL ERROR: Local #LOBS=",iobsl," > Total #LOBS=",nobsl
        WRITE(6,'(A,2I10)') 'IJ,NN=',i,nn
        CALL finalize_mpi
        STOP
     END IF
     
     hdxf(iobsl,:) = REAL(obshdxf(nobs_use(n),:), r_size)
     dep(iobsl)    = REAL(obsdep(nobs_use(n)), r_size)

     !---Obs. error variance
     rdiag(iobsl) = REAL(obserr(nobs_use(n))**2, r_size)

     !---Localization function
     IF(sigma_obsv <= 0.d0)THEN !NO vertical localization
        rloc(iobsl) = REAL(EXP(-0.5d0 * (dist/sigma_obs)**2),r_size)
     ELSE
        rloc(iobsl) = REAL(EXP(-0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2)),r_size)
     ENDIF

  END DO !n

  DEALLOCATE(nobs_use)

END SUBROUTINE obs_local

!-------------------------------------

SUBROUTINE obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)

  USE common_setting
  IMPLICIT NONE

  INTEGER :: j,n,ib,ie,ip

  !---IN
  INTEGER,INTENT(IN) :: imin,imax,jmin,jmax

  !---OUT
  INTEGER,INTENT(OUT) :: nn,nobs_use(nobs)

  nn = 0

  DO j=jmin,jmax

     !ib: Sum #OBS at imin-1
     IF(imin > 1)THEN
        ib = nobsgrd(imin-1,j)+1
     ELSE
        IF(j > 1)THEN
           ib = nobsgrd(nlon,j-1)+1
        ELSE
           ib = 1
        END IF
     END IF

     !ie: Sum #OBS at imax
     ie = nobsgrd(imax,j)

     !n: #OBS imin-1 ~ imax
     n = ie - ib + 1 

     !No obs.
     IF(n == 0) CYCLE

     DO ip=ib,ie
        IF(nn > nobs) THEN
           WRITE(6,'(A,2I10)') "FATALERROR, NN > NOBS", NN, NOBS
        END IF
        nn = nn + 1
        nobs_use(nn) = ip
     END DO !ip
     
  END DO !j

END SUBROUTINE obs_local_sub

!-----------------------------------------------------------------------
! Ensemble Mean/Spread/Perturbation |
!-----------------------------------------------------------------------

SUBROUTINE Ensemble_Mean_Spread_Perturbation(nbv,nvd,x,xmean,xsprd,dx)

  USE common_setting, only: r_size
  implicit none

  !---Common
  INTEGER ibv,ivd

  !---IN
  INTEGER,INTENT(IN) ::  nbv,nvd
  REAL(r_size),INTENT(IN) :: x(nbv,nvd)

  !---OUT
  REAL(r_size),INTENT(OUT) :: xmean(nvd),xsprd(nvd),dx(nbv,nvd)

  !xmean
  DO ivd=1,nvd

     xmean(ivd)=0.d0
     
     DO ibv=1,nbv
        xmean(ivd)=xmean(ivd)+x(ibv,ivd)
     END DO

     xmean(ivd)=xmean(ivd) / REAL(nbv,r_size)
     
  END DO

  !xspread
  DO ivd=1,nvd

     xsprd(ivd)=0.d0

     DO ibv=1,nbv
        xsprd(ivd)=xsprd(ivd)+(x(ibv,ivd)-xmean(ivd))*(x(ibv,ivd)-xmean(ivd))
     END DO

     xsprd(ivd)=sqrt(xsprd(ivd) / REAL(nbv-1,r_size))
  END DO

  !dx
  DO ivd=1,nvd
     DO ibv=1,nbv
        dx(ibv,ivd)=x(ibv,ivd)-xmean(ivd)
     END DO
  END DO
  
END SUBROUTINE Ensemble_Mean_Spread_Perturbation

!------------------------------------------------------------------
! Analysis Ensemble |
!------------------------------------------------------------------

SUBROUTINE Analysis_Ensemble(nbv,nvd,fcstm,dxf,wvec,Wmat,anal)

  USE common_setting, only: r_size
  IMPLICIT NONE

  !---Common
  INTEGER i,j,ivd

  !---IN
  INTEGER,INTENT(IN) :: nbv,nvd
  REAL(r_size),INTENT(IN) :: fcstm(nvd),dxf(nbv,nvd)
  REAL(r_size),INTENT(IN) :: wvec(nbv),Wmat(nbv,nbv)

  !---OUT
  REAL(r_size),INTENT(OUT) :: anal(nbv,nvd)

  anal(1,:)=fcstm(:)

  DO ivd=1,nvd
     DO j=1,nbv
        anal(j,ivd)=fcstm(ivd)
        DO i=1,nbv
           anal(j,ivd)=anal(j,ivd)+dxf(i,ivd)*(wvec(i)+Wmat(i,j))
        END DO
     END DO
  END DO

END SUBROUTINE Analysis_Ensemble

!----------------------------------------------------------------
! RTPP (Zhang et al. 2004; Kotsuki et al. 2017)|
!----------------------------------------------------------------

SUBROUTINE RTPP(nbv,nvd,dxf,dxa,analm,anal)

  USE common_setting,only :r_size, ALPHA_RTPP
  IMPLICIT NONE

  !---Common
  INTEGER ibv,ivd 

  !---IN
  INTEGER,INTENT(IN) :: nbv,nvd

  REAL(r_size),INTENT(IN) :: dxf(nbv,nvd),dxa(nbv,nvd)
  REAL(r_size),INTENT(IN) :: analm(nvd)

  !---INOUT
  REAL(r_size),INTENT(INOUT) :: anal(nbv,nvd)
  
  DO ivd=1,nvd
     DO ibv=1,nbv
        anal(ibv,ivd)=analm(ivd)+ALPHA_RTPP*dxf(ibv,ivd)+(1.d0-ALPHA_RTPP)*dxa(ibv,ivd)
     END DO
  END DO
  
END SUBROUTINE RTPP

!----------------------------------------------------------------
! RTPS (Whitaker and Hamill 2012; Kotsuki et al. 2017)|
!----------------------------------------------------------------

SUBROUTINE RTPS(nbv,nvd,fcsts,analm,anals,dxa,anal)

  USE common_setting,only :r_size, ALPHA_RTPS
  IMPLICIT NONE

  !---Common
  INTEGER ibv,ivd 

  REAL(r_size) rtps_factor
  
  !---IN
  INTEGER,INTENT(IN) :: nbv,nvd

  REAL(r_size),INTENT(IN) :: fcsts(nvd)  
  REAL(r_size),INTENT(IN) :: analm(nvd),anals(nvd)  
  REAL(r_size),INTENT(IN) :: dxa(nbv,nvd)

  !---INOUT
  REAL(r_size),INTENT(INOUT) :: anal(nbv,nvd)

  DO ivd=1,nvd

     rtps_factor=ALPHA_RTPS*fcsts(ivd)/anals(ivd)+(1.d0-ALPHA_RTPS)

     DO ibv=1,nbv
        anal(ibv,ivd)=analm(ivd)+dxa(ibv,ivd)*rtps_factor
     END DO
  END DO
  
END SUBROUTINE RTPS

!----------------------------------------------------------------
! Round off |
!----------------------------------------------------------------

SUBROUTINE round_off(anal3d)

  USE common_setting
  IMPLICIT NONE

  !---Common
  INTEGER i,k,ibv
  
  !---INOUT
  REAL(r_size),INTENT(INOUT) :: anal3d(nij1,nlev,nbv,nv3d)

  !$OMP PARALLEL DO PRIVATE(ibv,k,i)
  DO ibv=1,nbv
     DO k=2,nlev !Skip: k=1
        DO i=1,nij1
           
           IF(phi1(i) == 1.d0) CYCLE
           anal3d(i,k,ibv,iv3d_t) = MAX(anal3d(i,k,ibv,iv3d_t), tmin)
           anal3d(i,k,ibv,iv3d_t) = MIN(anal3d(i,k,ibv,iv3d_t), tmax)
           anal3d(i,k,ibv,iv3d_s) = MAX(anal3d(i,k,ibv,iv3d_s), smin)
           anal3d(i,k,ibv,iv3d_s) = MIN(anal3d(i,k,ibv,iv3d_s), smax)           
           
        END DO
     END DO
  END DO
  !$OMP END PARALLEL DO
  
END SUBROUTINE round_off
