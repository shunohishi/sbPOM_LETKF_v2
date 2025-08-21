MODULE lpfgm_tools
  !=======================================================================
  !
  ! [PURPOSE:] Module for LPFGM with POM
  !
  ! [HISTORY:]
  !   01/26/2009 Takemasa Miyoshi  created
  !   02/03/2009 Takemasa Miyoshi  modified for ROMS
  !   01/26/2011 Yasumasa Miyazawa  modified for POM (check 'pom' or 'POM')
  !   08/08/2025 Shun Ohishi modified for LPFGM
  !
  !=======================================================================

CONTAINS

  !-----------------------------------------------------------------------
  ! Data Assimilation
  !-----------------------------------------------------------------------

  SUBROUTINE das_lpfgm(fcst3d,fcst2d,anal3d,anal2d)

    !$USE OMP_LIB
    USE common_setting
    USE common
    USE common_lpfgm !Sakai-san will replace with common_lpfgm
    USE common_mpi
    USE letkf_tools
    IMPLICIT NONE

    INTEGER i,k
    INTEGER ibv
    INTEGER nobsl

    REAL(r_size) fcst3dm(nv3d),fcst2dm(nv2d)     !Forecast ensemble mean
    REAL(r_size) fcst3ds(nv3d),fcst2ds(nv2d)     !Forecast ensemble spread
    REAL(r_size) dxf3d(nbv,nv3d),dxf2d(nbv,nv2d) !Forecast ensemble perturbation

    REAL(r_size) anal3dm(nv3d),anal2dm(nv2d)     !Analysis ensemble mean
    REAL(r_size) anal3ds(nv3d),anal2ds(nv2d)     !Analysis ensemble spread
    REAL(r_size) dxa3d(nbv,nv3d),dxa2d(nbv,nv2d) !Analysis ensemble perturbation

    REAL(r_dble),ALLOCATABLE :: hdxf(:,:) !dYf
    REAL(r_dble),ALLOCATABLE :: rdiag(:)  !Obs. error variance
    REAL(r_dble),ALLOCATABLE :: rloc(:)   !Localization function
    REAL(r_dble),ALLOCATABLE :: dep(:)    !Innovation

    REAL(r_size) :: Imat(nbv,nbv) !Identity matrix
    REAL(r_size) :: pa(nbv,nbv)   !Pa matrix
    REAL(r_size) :: wvec(nbv)     !w vector
    REAL(r_size) :: Wmat(nbv,nbv) !Transform matrix W

    !---IN
    REAL(r_size),INTENT(IN) :: fcst3d(nij1,nlev,nbv,nv3d),fcst2d(nij1,nbv,nv2d) !Ensemble forecast

    !---OUT
    REAL(r_size),INTENT(OUT) :: anal3d(nij1,nlev,nbv,nv3d),anal2d(nij1,nbv,nv2d) !Ensemble analysis

    WRITE(file_unit,'(A,I8)') "Target observation numbers : NOBS=",nobs

    !---No Global Observation
    IF(nobs == 0)THEN
       anal3d(:,:,:,:)=fcst3d(:,:,:,:)
       anal2d(:,:,:)=fcst2d(:,:,:)
       RETURN
    END IF

    !---Identity matrix
    Imat(:,:)=0.d0
    do ibv=1,nbv
       Imat(ibv,ibv)=1.d0
    end do
    
    !---Initialization
    anal3d(:,:,:,:)=0.d0
    anal2d(:,:,:)=0.d0
    
    !---Main Analysis Loop
    !$OMP PARALLEL DO PRIVATE(i,k,fcst2dm,fcst2ds,dxf2d,fcst3dm,fcst3ds,dxf3d,anal2dm,anal2ds,dxa2d,anal3dm,anal3ds,dxa3d,hdxf,rdiag,rloc,dep,nobsl,pa,wvec,Wmat)
    DO k=1,nlev
       DO i=1,nij1

          CALL count_nobsl(i,k,depth1(i,k),nobsl)

          IF(nobsl == 0)THEN
             !Sakai-san will modify nobsl = 0 case.
             wvec(:)=0.d0
             Wmat(:,:)=Imat(:,:)
            !  CALL letkf_core_noobs(pa,wvec,Wmat)
          ELSE
             ALLOCATE(hdxf(nobsl,nbv))
             ALLOCATE(rdiag(nobsl),rloc(nobsl),dep(nobsl))
             CALL obs_local(i,k,nobsl,depth1(i,k), &
                  & hdxf,rdiag,rloc,dep)

             ! Sakai-san will replace with lpfgm_core
             CALL lpfgm_core(nobsl,hdxf,rdiag,rloc,dep,pa,wvec,Wmat)

             DEALLOCATE(hdxf,rdiag,rloc,dep)
          END IF

          CALL Ensemble_Mean_Spread_Perturbation(nbv,nv3d,fcst3d(i,k,:,:),fcst3dm,fcst3ds,dxf3d)          
          CALL Analysis_Ensemble_LPFGM(nbv,nv3d,fcst3dm,dxf3d,wvec,Wmat,anal3d(i,k,:,:))
          CALL Ensemble_Mean_Spread_Perturbation(nbv,nv3d,anal3d(i,k,:,:),anal3dm,anal3ds,dxa3d)

          IF(k == nlev)THEN
             CALL Ensemble_Mean_Spread_Perturbation(nbv,nv2d,fcst2d(i,:,:),fcst2dm,fcst2ds,dxf2d)
             CALL Analysis_Ensemble_LPFGM(nbv,nv2d,fcst2dm,dxf2d,wvec,Wmat,anal2d(i,:,:))
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

  END SUBROUTINE das_lpfgm

END MODULE lpfgm_tools

!------------------------------------------------------------------
! Ensemble analysis |
!------------------------------------------------------------------

SUBROUTINE Analysis_Ensemble_LPFGM(nbv,nvd,fcstm,dxf,wvec,Wmat,anal)

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

END SUBROUTINE Analysis_Ensemble_LPFGM

!----------------------------------------------------------------
! RTPP |
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
! RTPP |
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
