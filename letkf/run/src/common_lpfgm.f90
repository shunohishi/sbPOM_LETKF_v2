MODULE common_lpfgm
  
  !=======================================================================
  !
  ! [PURPOSE] Local Particle Filter Gaussian Mixture (LPFGM)
  !           Model Independent Core Module
  !
  ! [REFERENCES]
  ! [1] Kotsuki et al., 2022: A local particle filter and its Gaussian mixture
  !  extension implemented with minor modifications to the LETKF.
  !  Geosci. Model Dev., 15, 8325–8348. 
  !
  ! [HISTORY]
  !  08/21/2025 Naoki Sakai     created based on the codes from sbPOM-LETKF v2 and speedy.
  !  10/20/2025 Shun Ohishi     reviewed and modified
  !
  !=======================================================================


CONTAINS

  !-----------------------------------------------------------------------
  !  Compute normalized particle filter weights (Gaussian likelihood)
  !    INPUT
  !      nobsl          : Number of assimilated local observations
  !      dep(nobsl)     : Innovation (y - Hxfmean)
  !      hdxf(nobsl,nbv): Ensemble perturbation in obs space (dYf=Hxf_i - Hxfmean)
  !      rloc(nobsl)    : Localization weightning function
  !      rdiag(nobsl)   : Observation error variance
  !
  !    OUTPUT
  !      wgt(nbv)       : Analysis particle weight [r_dble]
  !
  !    NOTES
  !      - r_dble is used to avoid underflow
  !-----------------------------------------------------------------------

  SUBROUTINE pf_weight_gauss_likelihood(nobsl,dep,hdxf,rloc,rdiag,wgt)
    
    USE common_setting,only: r_size,r_dble,nbv
    IMPLICIT NONE

    INTEGER i,j

    REAL(r_dble) log_like(nbv)  !Log[p(y|xf(ibv))]
    REAL(r_dble) log_like_max
    REAL(r_dble) exp_dlog(nbv) !exp(log-log_max)
    
    !REAL(r_dble) sqpf, qtmp
    !REAL(r_dble) qpf(nbv),dep2_Ri(nbv)

    !---IN
    INTEGER,INTENT(IN) :: nobsl
    REAL(r_size),INTENT(IN) :: dep(nobsl)
    REAL(r_size),INTENT(IN) :: hdxf(nobsl,nbv)
    REAL(r_size),INTENT(IN) :: rloc(nobsl)
    REAL(r_size),INTENT(IN) :: rdiag(nobsl)

    !---OUT
    REAL(r_dble),INTENT(OUT) :: wgt(nbv) 

    !---Log-Sum-Exp (LSE) function
    !-Log[p(y|xf(ibv))]: 0.5*dob(j)^T R^-1 dob(j)
    log_like(:)=0.d0
    DO j=1,nbv
       DO i=1,nobsl
          log_like(j)=log_like(j)-0.5d0*(dep(i)-hdxf(i,j))*(dep(i)-hdxf(i,j))*rloc(i)/rdiag(i)
       END DO
    END DO

    log_like_max=maxval(log_like)

    exp_dlog(:)=exp(log_like(:)-log_like_max)

    wgt(:)=exp_dlog(:)/sum(exp_dlog)

  END SUBROUTINE pf_weight_gauss_likelihood

  !-----------------------------------------------------------------------
  !  Generate resampling matrix from cumulative weights
  !
  !  INPUT
  !    cwgt(nbv)      : Cumulative weights in [0,1] [REAL(r_dble)]
  !
  !  OUTPUT
  !    RSmat(nbv,nbv) : ReSampling matrix with 0/1 entries [REAL(r_dble)]
  !-----------------------------------------------------------------------

  SUBROUTINE resampling_matrix(cwgt,RSmat)

    USE common_setting, only: r_size,r_dble,nbv,RM,DP
    USE common
    USE common_mpi
    IMPLICIT NONE

    INTEGER i,j,k
    INTEGER idx(nbv)  !Sorting index for Random number
    INTEGER inum(nbv) !Resampled particle index

    REAL(r_dble) rand(nbv) !Random number [0,1)

    LOGICAL c_used(nbv) !True if kth column has been assigned
    
    !---IN
    REAL(r_dble),INTENT(in) :: cwgt(nbv)

    !---OUT
    REAL(r_dble),INTENT(out) :: RSmat(nbv,nbv)

    !---Random number
    IF(RM == "SR")THEN      !Systematic Resampling
       
       CALL com_rand_seed(nbv,0,rand)
       
       DO i=1,nbv
          rand(i)=(rand(1)+dble(i)-1.0d0)/dble(nbv)
       END DO
    
    ELSE IF(RM == "MR")THEN !Multinominal Resampling

       DO i=1,nbv
          idx(i)=i
       END DO
       
       CALL com_rand_seed(nbv,0,rand)
       CALL quick_sort_asend(nbv,rand,idx,1,nbv)       
       
    ELSE
       
       WRITE(*,*) "***Error: Incorrect RM ==> "//RM
       CALL finalize_mpi
       STOP
       
    END IF

    !---Resampling matrix (m x m)
    IF(DP)THEN !---Diagonal Priority
       
       inum(:)=0
       
       !---Select particle
       DO j=1,nbv !jth column
          DO i=1,nbv
             IF(rand(j) <= cwgt(i) .OR. i == nbv)THEN
                inum(j)=i
                EXIT 
             ENDIF
          END DO
       END DO

       !---Make Resampling matrix
       c_used(:)=.false.
       RSmat(:,:)=0.0d0
       
       !Diagonal element (priority)
       DO i=1,nbv 

          k=inum(i)
          IF(c_used(k)) CYCLE
          
          RSmat(k,k)=1.0d0
          c_used(k)=.true.
          inum(i)=0
          
       END DO

       !Off-diagonal element: Assign remaining particlues (rows) to unused columns
       DO i=1,nbv

          IF(inum(i) == 0) CYCLE          
          j=inum(i)   
          
          DO k=1,nbv
             IF(c_used(k)) CYCLE
             RSmat(j,k)=1.0d0
             c_used(k)=.true.
             inum(i)=0
             EXIT
          END DO
          
       END DO

    ELSE !Non-diagonal Priority
       
       RSmat(:,:)=0.0d0
       DO j=1,nbv !jth column
          DO i=1,nbv
             IF(rand(j) <= cwgt(i) .OR. i == nbv)THEN
                RSmat(i,j)=1.0d0
                EXIT
             END IF
          END DO          
       END DO
       
    ENDIF

  END SUBROUTINE resampling_matrix

  !=======================================================================
  !  Main Subroutine of LPF Core
  !   INPUT
  !     nobsl : Number of assimilated observations at a model grid point
  !     hdxf(nobs,nbv) : Forecast ensemble perturbation in obs. space (=dYf)
  !     rdiag(nobs)    : Observation error variance (=sigma_o^2)
  !     rloc(nobs)     : Localization weighting function
  !     dep(nobs)      : Innovation (=y-Hxfmean)
  !
  !   OUTPUT
  !     wvec(nbv)           : w vector (=zero vector)
  !     Wmat(nbv,nbv)       : Transform matrix
  !=======================================================================

  SUBROUTINE lpf_core(nobsl,hdxf,rdiag,rloc,dep,wvec,Wmat)

    USE common_setting
    IMPLICIT NONE

    INTEGER i

    REAL(r_dble) wgt(nbv)      !Analysis particle weight
    REAL(r_dble) cwgt(nbv)     !Cumulative weight
    REAL(r_dble) meff          !Effective particle size
    REAL(r_dble) RSmat(nbv,nbv) !Resampling matrix (single resampling)

    !REAL(r_dble) acc(nbv), pmat(nbv, nbv)
    !REAL(r_dble) swgh, peff

    !---IN
    INTEGER, INTENT(IN) :: nobsl

    REAL(r_size),INTENT(IN) :: hdxf(nobsl, nbv)
    REAL(r_size),INTENT(IN) :: rdiag(nobsl)
    REAL(r_size),INTENT(IN) :: rloc(nobsl)
    REAL(r_size),INTENT(IN) :: dep(nobsl)

    !---OUT
    REAL(r_size), INTENT(OUT) :: wvec(nbv)      !Zero vector for LPF
    REAL(r_size), INTENT(OUT) :: Wmat(nbv, nbv) !Weight matrix (Resampling matrix ave.)

    !---Initialization
    wvec(:)=REAL(0.d0,r_size)
    Wmat(:,:)=REAL(0.d0,r_size)

    !---Analysis partcile wegiht
    !*** Assumption: Gaussian likelihood ***
    CALL pf_weight_gauss_likelihood(nobsl,dep,hdxf,rloc,rdiag,wgt)

    !---Effective particle size
    meff=1.0d0/sum(wgt(:)**2.0d0)

    !---Cumulative weight
    cwgt(1)=wgt(1)
    DO i=2,nbv
       cwgt(i)=cwgt(i-1)+wgt(i)
    END DO

    !---Resampling
    DO i=1,nmonte
       CALL resampling_matrix(cwgt,RSmat)
       Wmat(:,:) = Wmat(:,:) + RSmat(:,:) / dble(nmonte)
    END DO

  END SUBROUTINE lpf_core

  !=======================================================================
  !  Main Subroutine of GM Core
  !   INPUT
  !     nobsl           : Number of local observations at a model grid point
  !     hdxf(nobsl,nbv) : Forecast ensemble perturbation in obs. space (=dYf)
  !     rdiag(nobsl)    : Observation error variance (=sigma_o^2)
  !     rloc(nobsl)     : Localization weighting function
  !     dep(nobsl)      : Innovation (=y-Hxfmean)
  !
  !   OUTPUT
  !     wvec(nbv)          : w vector → zero vector in LPFGM
  !     Wmat(nbv,nbv)      : Transform matrix (LPFGM shift matrix + I)
  !=======================================================================

  SUBROUTINE gm_core(nobsl,hdxf,rdiag,rloc,dep,wvec,Wmat)

    USE common_setting
    USE common
    IMPLICIT NONE

    INTEGER i,j,k
    INTEGER np

    REAL(r_size) Rlinv_dYf(nobsl,nbv)
    REAL(r_dble) evec(nbv,nbv)
    REAL(r_dble) eval(nbv)
    REAL(r_size) work1(nbv,nbv)
    REAL(r_size) work2(nbv,nobsl)
    REAL(r_size) work3(nbv,nbv)
    
    !---IN
    INTEGER,INTENT(IN) :: nobsl
    REAL(r_size),INTENT(IN) :: hdxf(nobsl,nbv)
    REAL(r_size),INTENT(IN) :: rdiag(nobsl)
    REAL(r_size),INTENT(IN) :: rloc(nobsl)
    REAL(r_size),INTENT(IN) :: dep(nobsl)

    !---OUT
    REAL(r_size),INTENT(OUT) :: wvec(nbv)
    REAL(r_size),INTENT(OUT) :: Wmat(nbv,nbv)

    !-----------------------------------------------------------------------
    ! Rloc^-1 * dYf
    !-----------------------------------------------------------------------
    DO j=1,nbv
       DO i=1,nobsl
          Rlinv_dYf(i,j)=hdxf(i,j)*rloc(i)/rdiag(i)
       END DO
    END DO

    !-----------------------------------------------------------------------
    ! dYf^T * Rloc^-1 * dYf + (m-1)/gamma * I [nbv*nbv]
    !-----------------------------------------------------------------------
    IF(r_size == r_dble)THEN
       work1(:,:)=0.d0
       CALL dgemm('t','n',nbv,nbv,nobsl, &
            & 1.0d0,Rlinv_dYf,nobsl,hdxf,nobsl, &
            & 0.0d0,work1,nbv)
    ELSE IF(r_size == r_sngl)THEN
       work1(:,:)=0.e0
       CALL sgemm('t','n',nbv,nbv,nobsl, &
            & 1.0e0,Rlinv_dYf,nobsl,hdxf,nobsl, &
            & 0.0e0,work1,nbv)
    END IF
       
    DO i=1,nbv
       work1(i,i)=work1(i,i)+REAL(nbv - 1,r_size)/REAL(cov_infl_gm,r_size)
    END DO

    !-----------------------------------------------------------------------
    ! Eigenvalue decomposition to work1
    !-----------------------------------------------------------------------

    CALL mtx_eigen(nbv,work1,eval,evec,np)

    !-----------------------------------------------------------------------
    ! Pa = Evec * Eval^-1 * Evec^T [nbv*nbv]
    !-----------------------------------------------------------------------
    
    DO j=1,nbv
       DO i=1,nbv
          work1(i,j)=REAL(evec(i,j)/eval(j),r_size)
       END DO
    END DO

    IF(r_size == r_dble)THEN
       work3(:,:)=0.d0
       CALL dgemm('n','t',nbv,nbv,nbv, &
            & 1.0d0,work1,nbv,evec,nbv, &
            & 0.0d0,work3,nbv)
    ELSE IF(r_size == r_sngl)THEN
       work3(:,:)=0.e0
       CALL sgemm('n','t',nbv,nbv,nbv, &
            & 1.0e0,work1,nbv,REAL(evec,r_sngl),nbv, &
            & 0.0e0,work3,nbv)
    END IF

    !-----------------------------------------------------------------------
    ! work2 = Pa * (Rloc^-1 * dYf)^T [nbv*nobsl]
    !-----------------------------------------------------------------------

    IF(r_size == r_dble)THEN
       work2(:,:)=0.d0
       CALL dgemm('n','t',nbv,nobsl,nbv, &
            & 1.0d0,work3,nbv,Rlinv_dYf,nobsl, &
            & 0.0d0,work2,nbv)
    ELSE IF(r_size == r_sngl)THEN
       work2(:,:)=0.e0
       CALL sgemm('n','t',nbv,nobsl,nbv, &
            & 1.0e0,work3,nbv,Rlinv_dYf,nobsl, &
            & 0.0e0,work2,nbv)
    END IF


    !***Heareafter, different with LETKF***    
    !-----------------------------------------------------------------------
    ! Wmat = work2 * (dep - dYf) + I [nbv*nbv]
    !-----------------------------------------------------------------------
    
    DO k=1,nbv
       DO i=1,nbv
          
          Wmat(i,k)=work2(i,1)*(dep(1)-hdxf(1,k))
          
          DO j=2,nobsl
             Wmat(i,k)=Wmat(i,k)+work2(i,j)*(dep(j)-hdxf(j,k))
          END DO
       END DO
    END DO
    
    DO i=1,nbv
       Wmat(i,i)=Wmat(i,i)+REAL(1.0d0,r_size)
    END DO

    !-----------------------------------------------------------------------
    ! wvec = 0
    !-----------------------------------------------------------------------

    wvec(:)=REAL(0.0d0,r_size)

  END SUBROUTINE gm_core

  !=======================================================================
  !  LPFGM main driver: compose GM and LPF transforms
  !    INPUT
  !      nobsl          : Number of assimilated observations
  !      hdxf(nobsl,nbv): Ensemble perturbation in obs. space (=dYf)
  !      rdiag(nobsl)   : Observation error variance
  !      rloc(nobsl)    : Localization weight function
  !      dep(nobsl)     : Innovation (y - Hxfmean)
  !
  !    OUTPUT
  !      wvec(nbv)      : Mean update (0 for LPF and LPFGM, r_size)
  !      Wmat(nbv,nbv)  : Transform matrix (r_size)
  !
  !=======================================================================

  SUBROUTINE lpfgm_core_noobs(wvec,Wmat)

    USE common_setting, only: r_size,nbv
    IMPLICIT NONE

    INTEGER i

    !---OUT
    REAL(r_size),INTENT(OUT) :: wvec(nbv)
    REAL(r_size),INTENT(OUT) :: Wmat(nbv,nbv)

    wvec(:)=REAL(0.d0,r_size)

    Wmat(:,:)=REAL(0.d0,r_size)
    DO i=1,nbv
       Wmat(i,i)=REAL(1.d0,r_size)
    END DO

  END SUBROUTINE lpfgm_core_noobs

  !-----------------------------------------------------------------------
  
  SUBROUTINE lpfgm_core(nobsl,hdxf,rdiag,rloc,dep,wvec,Wmat)

    USE common_setting, only: r_size,r_sngl,r_dble,nbv,iswitch_da
    IMPLICIT NONE

    !---WORK
    REAL(r_size) W_LPF(nbv,nbv)
    REAL(r_size) W_GM(nbv,nbv)
    REAL(r_size) wtmp(nbv)

    !---IN
    INTEGER,INTENT(IN)  :: nobsl

    REAL(r_size),INTENT(IN)  :: hdxf(nobsl,nbv)
    REAL(r_size),INTENT(IN)  :: rdiag(nobsl)
    REAL(r_size),INTENT(IN)  :: rloc(nobsl)
    REAL(r_size),INTENT(IN)  :: dep(nobsl)

    !---OUT
    REAL(r_size),INTENT(OUT) :: wvec(nbv)
    REAL(r_size),INTENT(OUT) :: Wmat(nbv,nbv)

    !Initialization
    wvec(:)=REAL(0.0d0,r_size) !*wvec=0 in both LPF and LPFGM
    Wmat(:,:)=REAL(0.0d0,r_size)
    
    IF(iswitch_da == 2)THEN      !LPF
       
       CALL lpf_core(nobsl, hdxf, rdiag, rloc, dep, wtmp, Wmat)
       
    ELSE IF(iswitch_da == 3)THEN !LPFGM
       
       CALL lpf_core(nobsl, hdxf, rdiag, rloc, dep, wtmp, W_LPF)
       CALL gm_core(nobsl, hdxf, rdiag, rloc, dep, wtmp, W_GM)
       
       IF(r_size == r_dble)THEN
          CALL dgemm('n','n',nbv,nbv,nbv, &
               & 1.0d0,W_GM,nbv,W_LPF,nbv, &
               & 0.0d0,Wmat,nbv)
       ELSE IF(r_size == r_sngl)THEN
          CALL sgemm('n','n',nbv,nbv,nbv, &
               & 1.0e0,W_GM,nbv,W_LPF,nbv, &
               & 0.0e0,Wmat,nbv)
       END IF
             
    ELSE
       write(*,'(a,i6)') "***Error: Incorrect iswitch_da ==> ",iswitch_da
       STOP
    END IF

  END SUBROUTINE lpfgm_core

END MODULE common_lpfgm
