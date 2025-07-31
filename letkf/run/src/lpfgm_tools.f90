MODULE letkf_tools
  !=======================================================================
  !
  ! [PURPOSE:] Module for LETKF with POM
  !
  ! [HISTORY:]
  !   01/26/2009 Takemasa Miyoshi  created
  !   02/03/2009 Takemasa Miyoshi  modified for ROMS
  !   01/26/2011 Yasumasa Miyazawa  modified for POM (check 'pom' or 'POM')
  !
  !=======================================================================

CONTAINS

  !-----------------------------------------------------------------------
  ! Data Assimilation
  !-----------------------------------------------------------------------

  SUBROUTINE das_letkf(fcst3d,fcst2d,anal3d,anal2d)

    !$USE OMP_LIB
    USE common_setting
    USE common
    USE common_letkf
    USE common_mpi
    IMPLICIT NONE

    INTEGER i,k
    INTEGER nobsl

    REAL(r_size) fcst3dm(nv3d),fcst2dm(nv2d) !Ensemble mean
    REAL(r_size) dxf3d(nbv,nv3d),dxf2d(nbv,nv2d) !Ensemble perturvation

    REAL(r_dble),ALLOCATABLE :: hdxf(:,:) !dYf
    REAL(r_dble),ALLOCATABLE :: rdiag(:)  !Obs. error variance
    REAL(r_dble),ALLOCATABLE :: rloc(:)   !Localization function
    REAL(r_dble),ALLOCATABLE :: dep(:)    !Innovation

    REAL(r_size) :: pa(nbv,nbv)   !Pa matrix
    REAL(r_size) :: wvec(nbv)     !w vector
    REAL(r_size) :: Wmat(nbv,nbv) !Transform matrix W
    REAL(r_size) :: trans(nbv,nbv,nv3d+nv2d) !W+w

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
    !$OMP PARALLEL DO PRIVATE(i,k,fcst2dm,dxf2d,fcst3dm,dxf3d,hdxf,rdiag,rloc,dep,nobsl,pa,wvec,Wmat,trans)
    DO k=1,nlev
       DO i=1,nij1
          
          call Ensemble_Mean_Perturbation(nbv,nv2d,fcst2d(i,:,:),fcst2dm,dxf2d)
          call Ensemble_Mean_Perturbation(nbv,nv3d,fcst3d(i,k,:,:),fcst3dm,dxf3d)

          
          CALL count_nobsl(i,k,depth1(i,k),nobsl)

          IF(nobsl == 0)THEN
             CALL letkf_core_noobs(pa,wvec,Wmat)
          ELSE
             ALLOCATE(hdxf(nobsl,nbv))
             ALLOCATE(rdiag(nobsl),rloc(nobsl),dep(nobsl))
             CALL obs_local(i,k,nobsl,depth1(i,k), &
                  & hdxf,rdiag,rloc,dep)
             CALL letkf_core(nobsl,hdxf,rdiag,rloc,dep,pa,wvec,Wmat)
             DEALLOCATE(hdxf,rdiag,rloc,dep)
          END IF
          
          CALL Trans_with_RTP(dxf3d,dxf2d,pa,wvec,Wmat,trans)

          CALL Analysis_Ensemble(nbv,nv3d,fcst3dm,dxf3d,trans(:,:,1:nv3d),anal3d(i,k,:,:))
          IF(k == nlev)THEN
             CALL Analysis_Ensemble(nbv,nv2d,fcst2dm,dxf2d,trans(:,:,nv3d+1:nv3d+nv2d),anal2d(i,:,:))
          END IF
          
       END DO !i
    END DO !k
    !$OMP END PARALLEL DO
    
    !---Round off
    IF(ROFF)THEN
       CALL round_off(anal3d)
    END IF

  END SUBROUTINE das_letkf

END MODULE letkf_tools

!-----------------------------------------------------------------------
! Ensemble perturbation |
!-----------------------------------------------------------------------

SUBROUTINE Ensemble_Mean_Perturbation(nbv,nvd,x,xmean,dx)

  USE common_setting, only: r_size
  implicit none

  INTEGER ibv,ivd
  
  INTEGER,INTENT(IN) ::  nbv,nvd
  REAL(r_size),INTENT(IN) :: x(nbv,nvd)
  
  REAL(r_size),INTENT(OUT) :: xmean(nvd),dx(nbv,nvd)

  !xmean
  DO ivd=1,nvd

     xmean(ivd)=0.d0
     
     DO ibv=1,nbv
        xmean(ivd)=xmean(ivd)+x(ibv,ivd)
     END DO

     xmean(ivd)=xmean(ivd) / REAL(nbv,r_size)
     
  END DO

  !dx
  DO ivd=1,nvd
     DO ibv=1,nbv
        dx(ibv,ivd)=x(ibv,ivd)-xmean(ivd)
     END DO
  END DO
  
END SUBROUTINE Ensemble_Mean_Perturbation

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
  REAL(r_dble),INTENT(OUT) :: hdxf(nobsl,nbv) !dYf
  REAL(r_dble),INTENT(OUT) :: rdiag(nobsl)    !Obs. error variance
  REAL(r_dble),INTENT(OUT) :: rloc(nobsl)     !Localization function
  REAL(r_dble),INTENT(OUT) :: dep(nobsl)      !Innovation

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
     
     hdxf(iobsl,:) = REAL(obshdxf(nobs_use(n),:), r_dble)
     dep(iobsl)    = REAL(obsdep(nobs_use(n)), r_dble)

     !---Obs. error variance
     rdiag(iobsl) = REAL(obserr(nobs_use(n))**2, r_dble)

     !---Localization function
     IF(sigma_obsv <= 0.d0)THEN !NO vertical localization
        rloc(iobsl) = REAL(EXP(-0.5d0 * (dist/sigma_obs)**2),r_dble)
     ELSE
        rloc(iobsl) = REAL(EXP(-0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2)),r_dble)
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

!---------------------------------------------------------
! Relaxation-to-prior perturbation/spread |
!-----------------------------------------
!
! RTPP: Zhang et al. (2004)
! RTPS: Whitake and Hamill (2012)
! *Cite: Kotsuki et al. (2017)
!
!---------------------------------------------------------

SUBROUTINE Trans_with_RTP(dxf3d,dxf2d,pa,wvec,Wmat,trans)

  USE common_setting
  USE common_mpi
  IMPLICIT NONE

  !---Common
  INTEGER i,j
  INTEGER ivd

  REAL(r_size) sprdf,sprda

  !---IN
  REAL(r_size),INTENT(IN) :: dxf3d(nbv,nv3d),dxf2d(nbv,nv2d)
  REAL(r_size),INTENT(IN) :: pa(nbv,nbv)
  REAL(r_size),INTENT(IN) :: wvec(nbv)
  
  !---INOUT
  REAL(r_size),INTENT(INOUT) :: Wmat(nbv,nbv)

  !---OUT
  REAL(r_size),INTENT(OUT) :: trans(nbv,nbv,nv3d+nv2d)

  IF(0.d0 <= ALPHA_RTPP .and. ALPHA_RTPS == 0.d0)THEN !No relaxation + RTPP

     !---W = alpha * I + (1-alpha)*W
     Wmat(:,:)=(1.d0-ALPHA_RTPP)*Wmat(:,:)

     DO i=1,nbv
        Wmat(i,i)=ALPHA_RTPP+Wmat(i,i)
     END DO

     DO j=1,nbv
        DO i=1,nbv
           trans(i,j,1)=Wmat(i,j)+wvec(i)
        END DO
     END DO

     DO ivd=2,nv3d+nv2d
        trans(:,:,ivd)=trans(:,:,1)
     END DO


  ELSE IF(0.d0 < ALPHA_RTPS .and. ALPHA_RTPP == 0.d0)THEN !RTPS

     DO ivd=1,nv3d+nv2d

        IF(ivd <= nv3d)THEN !3D
           CALL ens_sprd(dxf3d(:,ivd),pa,sprdf,sprda)
        ELSE !2D
           CALL ens_sprd(dxf2d(:,ivd-nv3d),pa,sprdf,sprda)
        END IF

        DO j=1,nbv
           DO i=1,nbv
              trans(i,j,ivd)=(ALPHA_RTPS*sprdf/sprda + 1-ALPHA_RTPS)*Wmat(i,j)
              trans(i,j,ivd)=trans(i,j,ivd)+wvec(i)
           END DO
        END DO

     END DO

  ELSE

     WRITE(6,'(A,F6.2)') "***Error: RTPP parameter: ", ALPHA_RTPP
     WRITE(6,'(A,F6.2)') "***Error: RTPS parameter: ", ALPHA_RTPS
     CALL finalize_mpi
     STOP

  END IF

END SUBROUTINE TRANS_WITH_RTP

!-----------------------------------------------------------------
! Forecast and analysis ensemble spread |
!-----------------------------------------------------------------

SUBROUTINE ens_sprd(dxf,pa,sprdf,sprda)

  USE common_setting
  IMPLICIT NONE

  !---Common
  INTEGER i,j

  !---IN
  REAL(r_size),INTENT(IN) :: dxf(nbv)
  REAL(r_size),INTENT(IN) :: pa(nbv,nbv)

  !---OUT
  REAL(r_size),INTENT(OUT) :: sprdf,sprda

  !diag(Xf Xf^T)/(nbv-1)
  sprdf=0.d0
  do i=1,nbv
     sprdf=sprdf+dxf(i)*dxf(i)
  end do
  sprdf=SQRT( sprdf/REAL(nbv-1.d0,r_size) )

  !diag(Xa Xa^T)/(nbv-1) = diag(Xf W W^T Xf^T)/(nbv-1) = diag(Xf Pa Xf^T) 
  sprda=0.d0
  do j=1,nbv
     do i=1,nbv
        sprda=sprda+dxf(i)*pa(i,j)*dxf(j)
     end do
  end do
  sprda=SQRT( sprda )

END SUBROUTINE ens_sprd

!------------------------------------------------------------------
! Ensemble analysis |
!------------------------------------------------------------------

SUBROUTINE Analysis_Ensemble(nbv,nvd,fcstm,dxf,trans,anal)

  USE common_setting, only: r_size
  IMPLICIT NONE

  !---Common
  INTEGER i,j,ivd
  
  !---IN
  INTEGER,INTENT(IN) :: nbv,nvd
  REAL(r_size),INTENT(IN) :: fcstm(nvd),dxf(nbv,nvd),trans(nbv,nbv,nvd)

  !---OUT
  REAL(r_size),INTENT(OUT) :: anal(nbv,nvd)

  anal(1,:)=fcstm(:)
  
  DO ivd=1,nvd
     DO j=1,nbv
        anal(j,ivd)=fcstm(ivd)
        DO i=1,nbv
           anal(j,ivd)=anal(j,ivd)+dxf(i,ivd)*trans(i,j,ivd)
        END DO
     END DO
  END DO
  
END SUBROUTINE Analysis_Ensemble


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
