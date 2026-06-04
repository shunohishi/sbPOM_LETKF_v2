MODULE common_letkf

  !=======================================================================
  !
  ! [PURPOSE:] Local Ensemble Transform Kalman Filtering (LETKF)
  !            Model Independent Core Module
  !
  ! [REFERENCES:]
  !  [1] Ott et al., 2004: A local ensemble Kalman filter for atmospheric
  !    data assimilation. Tellus, 56A, 415-428.
  !  [2] Hunt et al., 2007: Efficient Data Assimilation for Spatiotemporal
  !    Chaos: A Local Ensemble Transform Kalman Filter. Physica D, 230,
  !    112-126.
  !
  ! [HISTORY:]
  !  01/21/2009 Takemasa Miyoshi  Created at U. of Maryland, College Park
  !  01/04/2024 Shun OHISHI       Modified for sbPOM-LETKF v2
  !  16/10/2025 Shun OHISHI       Change precision
  !=======================================================================

CONTAINS

  !=======================================================================
  !  Main Subroutine of LETKF Core
  !   INPUT
  !     nobsl            : Number of assimilated observations at a model grid point
  !     hdxf(nobsl,nmem)  : Forecast ensemble perturbation in obs. space (=dYf)
  !     rdiag(nobsl)     : Observation error variance (=sigma_o^2)
  !     rloc(nobsl)      : Localization weighting function
  !     dep(nobsl)       : Innovation (=y-Hxfmean)
  !
  !   OUTPUT
  !     wvec(nmem)     : w vector (Update ensemble mean)
  !     Wmat(nmem,nmem) : Transform matrix
  !=======================================================================

  SUBROUTINE letkf_core_noobs(wvec,Wmat)

    USE common_setting
    IMPLICIT NONE

    INTEGER i

    !---OUT
    REAL(r_size),INTENT(OUT) :: wvec(nmem)
    REAL(r_size),INTENT(OUT) :: Wmat(nmem,nmem)

    !Wmat = diag(rho^1/2)
    Wmat(:,:)=REAL(0.0d0,r_size)
    DO i=1,nmem
       Wmat(i,i)=SQRT(cov_infl_mul)
    END DO

    !wvec
    wvec(:)=REAL(0.0d0,r_size)

  END SUBROUTINE letkf_core_noobs

  !-------------------------------

  SUBROUTINE letkf_core(nobsl,hdxf,rdiag,rloc,dep, &
       & wvec,Wmat)

    USE common_setting
    USE common
    IMPLICIT NONE

    INTEGER i,j
    INTEGER np

    REAL(r_size) :: Rlinv_dYf(nobsl,nmem) !Rloc^-1_dYf
    REAL(r_dble) :: evec(nmem,nmem)        !Eigen vector (*double precision)
    REAL(r_dble) :: eval(nmem)            !Eigen value  (*double precision)
    REAL(r_size) :: work1(nmem,nmem)
    REAL(r_size) :: work2(nmem,nobsl)
    REAL(r_size) :: work3(nmem,nmem)

    !---IN
    INTEGER,INTENT(IN) :: nobsl    
    REAL(r_size),INTENT(IN) :: hdxf(nobsl,nmem)
    REAL(r_size),INTENT(IN) :: rdiag(nobsl)
    REAL(r_size),INTENT(IN) :: rloc(nobsl)
    REAL(r_size),INTENT(IN) :: dep(nobsl)

    !---OUT
    REAL(r_size),INTENT(OUT) :: wvec(nmem)
    REAL(r_size),INTENT(OUT) :: Wmat(nmem,nmem)
    
    !-----------------------------------------------------------------------
    !  Rloc^-1 dYf: [nobsl*nmem]
    !-----------------------------------------------------------------------

    DO j=1,nmem
       DO i=1,nobsl
          Rlinv_dYf(i,j)=hdxf(i,j)*rloc(i)/rdiag(i)
       END DO
    END DO

    !-----------------------------------------------------------------------
    !  dYf^T Rloc^-1 dYf + (m-1)/rho*I  [nmem*nmem]|
    !---------------------------------------------
    !
    ! dgemm('n(A)','n(B)',m,n,k, alpha,A,m,B,k, beta,C,m)
    ! C <= alpha*(AB)+beta*C
    !
    ! A: m*k size
    ! B: k*n size
    ! C: m*n size
    !
    !-----------------------------------------------------------------------

    IF(r_size == r_dble)THEN
       work1(:,:)=0.d0
       CALL dgemm('t','n',nmem,nmem,nobsl, &
            & 1.0d0,Rlinv_dYf,nobsl,hdxf,nobsl, &
            & 0.0d0,work1,nmem)
    ELSE IF(r_size == r_sngl)THEN
       work1(:,:)=0.e0
       CALL sgemm('t','n',nmem,nmem,nobsl, &
            & 1.0e0,Rlinv_dYf,nobsl,hdxf,nobsl, &
            & 0.0e0,work1,nmem)
    END IF

    DO i=1,nmem
       work1(i,i)=work1(i,i)+REAL(nmem-1,r_size)/REAL(cov_infl_mul,r_size)
    END DO
    
    !-----------------------------------------------------------------------
    !  Eigenvalue decomposition to work1[ dYf^T Rloc^-1 dYf + (m-1)I/rho]
    !-----------------------------------------------------------------------

    CALL mtx_eigen(nmem,work1,eval,evec,np)

    !-----------------------------------------------------------------------
    !  Pa = [ dYf^T Rloc^-1 dYf + (m-1) I/rho ]^-1 [nmem*nmem]
    !-----------------------------------------------------------------------

    !work1 = Evec * Eval^-1
    DO j=1,nmem
       DO i=1,nmem
          work1(i,j)=REAL(evec(i,j)/eval(j),r_size)
       END DO
    END DO

    !Pa = Evec * Eval^-1 * Evec^T
    IF(r_size == r_dble)THEN
       work3(:,:)=0.d0
       CALL dgemm('n','t',nmem,nmem,nmem, &
            & 1.0d0,work1,nmem,evec,nmem, &
            & 0.0d0,work3,nmem)
    ELSE IF(r_size == r_sngl)THEN
       work3(:,:)=0.e0
       CALL sgemm('n','t',nmem,nmem,nmem, &
            & 1.0e0,work1,nmem,REAL(evec,r_sngl),nmem, &
            & 0.0e0,work3,nmem)
    END IF
    
    !-----------------------------------------------------------------------
    !  work2 = Pa  (Rloc^-1 dYf)^T [nmem*nobsl]
    !-----------------------------------------------------------------------

    IF(r_size == r_dble)THEN
       work2(:,:)=0.d0
       CALL dgemm('n','t',nmem,nobsl,nmem, &
            & 1.0d0,work3,nmem,Rlinv_dYf,nobsl, &
            & 0.0d0,work2,nmem)
    ELSE IF(r_size == r_sngl)THEN
       work2(:,:)=0.e0
       CALL sgemm('n','t',nmem,nobsl,nmem, &
            & 1.0e0,work3,nmem,Rlinv_dYf,nobsl, &
            & 0.0e0,work2,nmem)
    END IF
    
    !-----------------------------------------------------------------------
    !  wvec = Pa dYf^T Rloc^-1 (y-Hxfmean) [nmem]
    !-----------------------------------------------------------------------

    wvec(:)=REAL(0.d0,r_size)
    DO i=1,nmem
       wvec(i)=work2(i,1)*dep(1)
       DO j=2,nobsl
          wvec(i)=wvec(i)+work2(i,j)*dep(j)
       END DO
    END DO

    !-----------------------------------------------------------------------
    !  Wmat = SQRT[(m-1)Pa]
    !-----------------------------------------------------------------------

    !Evec * Eval^-1/2
    DO j=1,nmem
       DO i=1,nmem
          work1(i,j)=evec(i,j)*SQRT(REAL(nmem-1,r_size)/eval(j))
       END DO
    END DO

    !Evec * Eval^-1/2 * Evec
    IF(r_size == r_dble)THEN
       Wmat(:,:)=0.d0
       CALL dgemm('n','t',nmem,nmem,nmem, &
            & 1.0d0,work1,nmem,evec,nmem, &
            & 0.0d0,Wmat,nmem)
    ELSE IF(r_size == r_sngl)THEN
       Wmat(:,:)=0.e0
       CALL sgemm('n','t',nmem,nmem,nmem, &
            & 1.0e0,work1,nmem,REAL(evec,r_sngl),nmem, &
            & 0.0e0,Wmat,nmem)
    END IF
    
  END SUBROUTINE letkf_core

END MODULE common_letkf
