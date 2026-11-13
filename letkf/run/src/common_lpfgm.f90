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
  !  2025/08/21 Naoki Sakai             Created based on the codes from sbPOM-LETKF v2 and speedy.
  !  2025/10/20 Shun Ohishi             Reviewed and modified
  !  2026/10/16 Shun Ohishi (Claude)    Modified using Claude (Note: Random number issues still remain)
  !  2026/10/28 Shun Ohishi (Claude)    Added thread-safe and reproducible random number (rand_lecuyer)
  !                                     Modified gm_core (including nobsl < nmem case)
  !
  !=======================================================================

CONTAINS

  !-----------------------------------------------------------------------
  !  Compute normalized particle filter weights (Gaussian likelihood)
  !    INPUT
  !      nobsl            : Number of assimilated local observations
  !      dobl(nobsl)      : Innovation (y - Hxfmean)
  !      hdxfl(nobsl,nmem): Ensemble perturbation in obs space (dYf=Hxf_i - Hxfmean)
  !      rloc(nobsl)      : Localization weight function
  !      rdiag(nobsl)     : Observation error variance
  !
  !    OUTPUT
  !      wgt(nmem)        : Analysis particle weight [r_dble]
  !
  !    NOTES
  !      - r_dble is used to avoid underflow
  !-----------------------------------------------------------------------

  SUBROUTINE pf_weight_gauss_likelihood(nobsl,dobl,hdxfl,rloc,rdiag,wgt)
    
    USE common_setting,only: r_size,r_dble,nmem
    IMPLICIT NONE

    !---Common
    INTEGER i,j

    REAL(r_dble) log_like(nmem)  !Log[p(y|xf(imem))]
    REAL(r_dble) log_like_max
    REAL(r_dble) exp_dlog(nmem) !exp(log-log_max)
    REAL(r_dble),ALLOCATABLE :: rinv(:) !0.5*rloc/rdiag [nobsl]
    
    !---IN
    INTEGER,INTENT(IN) :: nobsl
    REAL(r_size),INTENT(IN) :: dobl(nobsl)
    REAL(r_size),INTENT(IN) :: hdxfl(nobsl,nmem)
    REAL(r_size),INTENT(IN) :: rloc(nobsl)
    REAL(r_size),INTENT(IN) :: rdiag(nobsl)

    !---OUT
    REAL(r_dble),INTENT(OUT) :: wgt(nmem) 

    !---Log-Sum-Exp (LSE) function
    !-Log[p(y|xf(imem))]: -0.5*[d-dY(:,j)]^T Rloc^-1 [d-dY(:,j)]
    ALLOCATE(rinv(nobsl))
    DO i=1,nobsl
       rinv(i)=0.5d0*REAL(rloc(i),r_dble)/REAL(rdiag(i),r_dble)
    END DO

    log_like(:)=0.d0
    DO j=1,nmem
       DO i=1,nobsl
          log_like(j)=log_like(j)-rinv(i)*(REAL(dobl(i),r_dble)-REAL(hdxfl(i,j),r_dble))**2
       END DO
    END DO
    
    DEALLOCATE(rinv)
    
    log_like_max=maxval(log_like)

    exp_dlog(:)=exp(log_like(:)-log_like_max)

    wgt(:)=exp_dlog(:)/sum(exp_dlog)

  END SUBROUTINE pf_weight_gauss_likelihood

  !-----------------------------------------------------------------------
  !  Generate resampling matrix from cumulative weights
  !
  !  INPUT
  !    imonte           : Monte-Carlo loop index (1:nmonte)
  !    iseed_grd        : Seed for grid point (0: same random number at all grid points)
  !    cwgt(nmem)       : Cumulative weights in [0,1] [REAL(r_dble)]
  !
  !  OUTPUT
  !    RSmat(nmem,nmem) : ReSampling matrix with 0/1 entries [REAL(r_dble)]
  !-----------------------------------------------------------------------

  SUBROUTINE resampling_matrix(imonte,iseed_grd,cwgt,RSmat)

    USE common_setting, only: r_dble,nmem,RM,DP
    USE common
    USE mpi
    IMPLICIT NONE

    !---Common
    INTEGER i,j,k
    INTEGER ierr
    INTEGER idx(nmem)  !Sorting index for Random number
    INTEGER inum(nmem) !Resampled particle index

    REAL(r_dble) rand(nmem) !Random number [0,1)
    REAL(r_dble) r0(1)      !Random number for SR [0,1)
    
    LOGICAL c_used(nmem) !True if kth column has been assigned
    
    !---IN
    INTEGER,INTENT(in) :: iseed_grd
    INTEGER,INTENT(in) :: imonte

    REAL(r_dble),INTENT(in) :: cwgt(nmem)

    !---OUT
    REAL(r_dble),INTENT(out) :: RSmat(nmem,nmem)

    !---Random number
    IF(RM == "SR")THEN      !Systematic Resampling

       CALL rand_lecuyer(hash_seed(iseed_grd),hash_seed(imonte),1,r0)
       
       DO i=1,nmem
          rand(i)=(r0(1)+dble(i)-1.0d0)/dble(nmem)
       END DO
    
    ELSE IF(RM == "MR")THEN !Multinomial Resampling

       DO i=1,nmem
          idx(i)=i
       END DO
       
       CALL rand_lecuyer(hash_seed(iseed_grd),hash_seed(imonte),nmem,rand)
       
       CALL quick_sort_asend(nmem,rand,idx,1,nmem)       
       
    ELSE
       
       WRITE(*,*) "***Error: Incorrect RM ==> "//RM
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)

    END IF

    !---Resampling matrix (m x m)
    IF(DP)THEN !---Diagonal Priority
       
       inum(:)=0
       
       !---Select particle
       DO j=1,nmem !jth column
          DO i=1,nmem
             IF(rand(j) <= cwgt(i) .OR. i == nmem)THEN
                inum(j)=i
                EXIT 
             ENDIF
          END DO
       END DO

       !---Make Resampling matrix
       c_used(:)=.false.
       RSmat(:,:)=0.0d0
       
       !Diagonal element (priority)
       DO i=1,nmem 

          k=inum(i)
          IF(c_used(k)) CYCLE
          
          RSmat(k,k)=1.0d0
          c_used(k)=.true.
          inum(i)=0
          
       END DO

       !Off-diagonal element: Assign remaining particles (rows) to unused columns
       DO i=1,nmem

          IF(inum(i) == 0) CYCLE          
          j=inum(i)   
          
          DO k=1,nmem
             IF(c_used(k)) CYCLE
             RSmat(j,k)=1.0d0
             c_used(k)=.true.
             inum(i)=0
             EXIT
          END DO
          
       END DO

    ELSE !Non-diagonal Priority
       
       RSmat(:,:)=0.0d0
       DO j=1,nmem !jth column
          DO i=1,nmem
             IF(rand(j) <= cwgt(i) .OR. i == nmem)THEN
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
  !     iseed_grd        : Seed for grid point (0: same random number at all grid points)
  !     nobsl            : Number of assimilated observations at a model grid point
  !     hdxfl(nobs,nmem) : Forecast ensemble perturbation in obs. space (=dYf)
  !     rdiag(nobs)      : Observation error variance (=sigma_o^2)
  !     rloc(nobs)       : Localization weighting function
  !     dobl(nobs)       : Innovation (=y-Hxfmean)
  !
  !   OUTPUT
  !     wvec(nmem)       : w vector (=zero vector)
  !     Wmat(nmem,nmem)  : Transform matrix
  !=======================================================================

  SUBROUTINE lpf_core(iseed_grd,nobsl,hdxfl,rdiag,rloc,dobl,wvec,Wmat)

    USE common_setting
    IMPLICIT NONE

    INTEGER i

    REAL(r_dble) wgt(nmem)        !Analysis particle weight
    REAL(r_dble) cwgt(nmem)       !Cumulative weight
    REAL(r_dble) meff             !Effective particle size
    REAL(r_dble) RSmat(nmem,nmem) !Resampling matrix (single resampling)

    !---IN
    INTEGER,INTENT(IN) :: iseed_grd
    INTEGER,INTENT(IN) :: nobsl

    REAL(r_size),INTENT(IN) :: hdxfl(nobsl, nmem)
    REAL(r_size),INTENT(IN) :: rdiag(nobsl)
    REAL(r_size),INTENT(IN) :: rloc(nobsl)
    REAL(r_size),INTENT(IN) :: dobl(nobsl)

    !---OUT
    REAL(r_size), INTENT(OUT) :: wvec(nmem)       !Zero vector for LPF
    REAL(r_size), INTENT(OUT) :: Wmat(nmem, nmem) !Weight matrix (Resampling matrix ave.)

    !---Initialization
    wvec(:)=REAL(0.d0,r_size)
    Wmat(:,:)=REAL(0.d0,r_size)

    !---Analysis particle weight
    !*** Assumption: Gaussian likelihood ***
    CALL pf_weight_gauss_likelihood(nobsl,dobl,hdxfl,rloc,rdiag,wgt)

    !---Effective particle size
    meff=1.0d0/sum(wgt(:)**2.0d0)

    !---Resampling (only if meff < m_crit)
    IF(m_crit <= meff)THEN

       DO i=1,nmem
          Wmat(i,i)=REAL(1.d0,r_size)
       END DO
       
       RETURN

    END IF
    
    !---Cumulative weight
    cwgt(1)=wgt(1)
    DO i=2,nmem
       cwgt(i)=cwgt(i-1)+wgt(i)
    END DO

    !---Resampling
    DO i=1,nmonte
       CALL resampling_matrix(i,iseed_grd,cwgt,RSmat)
       Wmat(:,:) = Wmat(:,:) + RSmat(:,:) / dble(nmonte)
    END DO

  END SUBROUTINE lpf_core

  !=======================================================================
  !  GM core: Gaussian mixture (Kalman-type) update of each particle
  !=======================================================================
  !
  !  [IN]  nobsl            : Local observation size (p)
  !        rdiag(nobsl)     : Observation error variance
  !        rloc(nobsl)      : Localization weight
  !        dobl(nobsl)      : Innovation y - H(xf_mean)
  !  [OUT] wvec(nmem)       : Zero vector
  !        Wmat(nmem,nmem)  : W_GM = I + Pa dY^T Rloc^-1 [d 1^T - dY]
  !  [INOUT] hdxfl(nobsl,nmem): IN: dY = H(Xf) - H(xf_mean), OUT: dYt (overwritten; call after lpf_core)
  !
  !  [Notation] (m: nmem, p: nobsl)
  !    alpha = (m-1)/gamma                       : gamma = cov_infl_gm
  !    dYt   = Rloc^-1/2 dY          [p*m]
  !    dt    = Rloc^-1/2 [y-H(xf_mean)] [p]
  !    Pa    = [alpha I + dYt^T dYt]^-1 [m*m]
  !    wbar  = Pa dYt^T dt           [m]
  !
  !  Since Pa dYt^T dYt = Pa (Painv - alpha I) = I - alpha Pa,
  !    W_GM  = I + wbar 1^T - Pa dYt^T dYt = alpha Pa + wbar 1^T
  !
  !  (1) m <= p: Painv = alpha I + dYt^T dYt = E diag(e) E^T
  !        wbar     = E diag(1/e) E^T dYt^T dt
  !        alpha Pa = C C^T,  C = E diag( SQRT(alpha/e) )
  !
  !  (2) p < m: G = alpha I + dYt dYt^T = U diag(d) U^T  (d >= alpha > 0)
  !        dYtU     = U^T dYt
  !        wbar     = dYtU^T diag(1/d) U^T dt             (push-through identity)
  !        alpha Pa = I - dYtU^T diag(1/d) dYtU
  !
  !  Cost: (1) O(m^3 + p m^2),  (2) O(p^3 + p^2 m + p m^2)
  !=======================================================================

  SUBROUTINE gm_core(nobsl,hdxfl,rdiag,rloc,dobl,wvec,Wmat)

    USE common_setting
    USE common
    IMPLICIT NONE

    !---Common
    INTEGER i,j
    INTEGER np                               !#Positive eigenvalues (from mtx_eigen)

    REAL(r_dble) alpha                       !(m-1)/gamma
    REAL(r_dble) wbar(nmem)                  !Pa dYt^T dt [m]

    REAL(r_size),ALLOCATABLE :: sqrt_rloc(:) !Rloc^-1/2 = SQRT(rloc/rdiag) [p]
    REAL(r_size),ALLOCATABLE :: dt(:)        !Rloc^-1/2 [y-H(xf_mean)]     [p]

    !---(1) m <= p
    REAL(r_size),ALLOCATABLE :: painv(:,:)   !Painv [m*m], then C
    REAL(r_dble),ALLOCATABLE :: emat(:,:)    !E [m*m]
    REAL(r_dble),ALLOCATABLE :: eval(:)      !e [m]
    REAL(r_dble),ALLOCATABLE :: tmp_m(:)     !dYt^T dt, then diag(1/e) E^T dYt^T dt [m]

    !---(2) p < m
    REAL(r_size),ALLOCATABLE :: gmat(:,:)    !G [p*p], then U (r_size)
    REAL(r_dble),ALLOCATABLE :: umat(:,:)    !U [p*p]
    REAL(r_dble),ALLOCATABLE :: dval(:)      !d [p]
    REAL(r_size),ALLOCATABLE :: dYtU(:,:)    !dYtU [p*m], then diag(d^-1/2) dYtU
    REAL(r_dble),ALLOCATABLE :: tmp_p(:)     !diag(1/d) U^T dt [p]

    !---IN
    INTEGER,INTENT(IN) :: nobsl
    REAL(r_size),INTENT(IN) :: rdiag(nobsl)
    REAL(r_size),INTENT(IN) :: rloc(nobsl)
    REAL(r_size),INTENT(IN) :: dobl(nobsl)

    !---OUT
    REAL(r_size),INTENT(OUT) :: wvec(nmem)
    REAL(r_size),INTENT(OUT) :: Wmat(nmem,nmem)

    !---INOUT
    REAL(r_size),INTENT(INOUT) :: hdxfl(nobsl,nmem) !IN: dY, OUT: dYt (overwritten)
    
    alpha = REAL(nmem-1,r_dble)/REAL(cov_infl_gm,r_dble) !alpha = (m-1)/rho

    !=====================================================================
    ! dYt = Rloc^-1/2 dY, dt = Rloc^-1/2 [y-H(xf_mean)]
    !=====================================================================

    ALLOCATE(dt(nobsl))
    ALLOCATE(sqrt_rloc(nobsl))

    DO i=1,nobsl
       sqrt_rloc(i)=SQRT(rloc(i)/rdiag(i))
       dt(i)=sqrt_rloc(i)*dobl(i)
    END DO

    DO j=1,nmem
       DO i=1,nobsl
          hdxfl(i,j)=hdxfl(i,j)*sqrt_rloc(i)
       END DO
    END DO

    DEALLOCATE(sqrt_rloc)

    IF(nmem <= nobsl)THEN

       !==================================================================
       ! (1) m <= p
       !==================================================================

       ALLOCATE(painv(nmem,nmem),emat(nmem,nmem),eval(nmem),tmp_m(nmem))

       !---Painv = alpha I + dYt^T dYt (upper triangle)
       painv(:,:)=REAL(0.d0, r_size)
       IF(r_size == r_dble)THEN
          CALL dsyrk('U','T',nmem,nobsl,1.d0,hdxfl,nobsl,0.d0,painv,nmem)
       ELSE IF(r_size == r_sngl)THEN
          CALL ssyrk('U','T',nmem,nobsl,1.e0,hdxfl,nobsl,0.e0,painv,nmem)
       END IF
       DO i=1,nmem
          painv(i,i)=painv(i,i)+REAL(alpha,r_size)
       END DO

       !---Painv = E diag(e) E^T
       CALL mtx_eigen(nmem,painv,eval,emat,np)

       !---wbar = E diag(1/e) E^T dYt^T dt
       DO j=1,nmem
          tmp_m(j)=DOT_PRODUCT(REAL(hdxfl(:,j),r_dble),REAL(dt(:),r_dble)) ! !dYt^T dt
       END DO
       tmp_m(:)=MATMUL(tmp_m(:),emat(:,:))/eval(:) !diag(1/e) E^T dYt^T dt
       wbar(:)=MATMUL(emat(:,:),tmp_m(:))          !E diag(1/e) E^T dYt^T dt

       !---alpha Pa = C C^T, C = E diag( SQRT(alpha/e) ) (painv is reused for C)
       DO j=1,nmem
          painv(:,j)=REAL(emat(:,j)*SQRT(alpha/eval(j)),r_size) ! C = E diag( SQRT(alpha/e) )
       END DO

       Wmat(:,:)=REAL(0.d0,r_size)
       IF(r_size == r_dble)THEN
          CALL dsyrk('U','N',nmem,nmem,1.d0,painv,nmem,0.d0,Wmat,nmem) ! alpha Pa = C C^T
       ELSE IF(r_size == r_sngl)THEN
          CALL ssyrk('U','N',nmem,nmem,1.e0,painv,nmem,0.e0,Wmat,nmem)
       END IF

       DEALLOCATE(painv,emat,eval,tmp_m)

    ELSE

       !==================================================================
       ! (2) p < m
       !==================================================================

       ALLOCATE(gmat(nobsl,nobsl),umat(nobsl,nobsl),dval(nobsl))
       ALLOCATE(dYtU(nobsl,nmem),tmp_p(nobsl))

       !---G = alpha I + dYt dYt^T (upper triangle)
       gmat=REAL(0.d0,r_size)
       IF(r_size == r_dble)THEN
          CALL dsyrk('U','N',nobsl,nmem,1.d0,hdxfl,nobsl,0.d0,gmat,nobsl) !dYt dYt^T
       ELSE IF(r_size == r_sngl)THEN
          CALL ssyrk('U','N',nobsl,nmem,1.e0,hdxfl,nobsl,0.e0,gmat,nobsl)
       END IF
       
       DO i=1,nobsl
          gmat(i,i)=gmat(i,i)+REAL(alpha,r_size) !alpha I + dYt dYt^T
       END DO

       !---G = U diag(d) U^T
       CALL mtx_eigen(nobsl,gmat,dval,umat,np)

       !---dYtU = U^T dYt (gmat is reused for U in r_size)
       gmat(:,:)=REAL(umat(:,:),r_size)
       IF(r_size == r_dble)THEN
          CALL dgemm('t','n',nobsl,nmem,nobsl,1.d0,gmat,nobsl,hdxfl,nobsl,0.d0,dYtU,nobsl) !U^T dYt
       ELSE IF(r_size == r_sngl)THEN
          CALL sgemm('t','n',nobsl,nmem,nobsl,1.e0,gmat,nobsl,hdxfl,nobsl,0.e0,dYtU,nobsl)
       END IF

       !---wbar = dYtU^T diag(1/d) U^T dt
       tmp_p(:)=MATMUL(REAL(dt(:),r_dble),umat(:,:))/dval(:) !diag(1/d) U^T dt
       
       DO j=1,nmem
          wbar(j)=DOT_PRODUCT(REAL(dYtU(:,j),r_dble),tmp_p(:)) !dYtU^T diag(1/d) U^T dt
       END DO

       !---diag(d^-1/2) dYtU (dYtU is reused)
       DO j=1,nmem
          DO i=1,nobsl
             dYtU(i,j)=dYtU(i,j)*REAL(1.d0/SQRT(dval(i)),r_size) !diag(d^-1/2) dYtU
          END DO
       END DO

       !---alpha Pa = I - dYtU^T diag(1/d) dYtU (upper triangle)
       Wmat(:,:)=REAL(0.d0,r_size)
       IF(r_size == r_dble)THEN
          CALL dsyrk('U','T',nmem,nobsl,-1.d0,dYtU,nobsl,0.d0,Wmat,nmem) !- dYtU^T diag(1/d) dYtU
       ELSE IF(r_size == r_sngl)THEN
          CALL ssyrk('U','T',nmem,nobsl,-1.e0,dYtU,nobsl,0.e0,Wmat,nmem)
       END IF
       
       DO i=1,nmem
          Wmat(i,i)=Wmat(i,i)+REAL(1.d0,r_size) !I - dYtU^T diag(1/d) dYtU
       END DO

       DEALLOCATE(gmat,umat,dval,dYtU,tmp_p)

    END IF

    DEALLOCATE(dt)

    !=====================================================================
    ! W_GM = alpha Pa + wbar 1^T (alpha Pa: upper triangle --> full)
    !=====================================================================

    DO j=1,nmem
       DO i=j+1,nmem
          Wmat(i,j)=Wmat(j,i)
       END DO
    END DO

    DO j=1,nmem
       DO i=1,nmem
          Wmat(i,j)=Wmat(i,j)+REAL(wbar(i),r_size)
       END DO
    END DO

    wvec(:)=REAL(0.d0,r_size)

  END SUBROUTINE gm_core
  
  !=======================================================================
  !  LPFGM main driver: compose GM and LPF transforms
  !    INPUT
  !      nobsl            : Number of assimilated observations
  !      hdxfl(nobsl,nmem): Ensemble perturbation in obs. space (=dYf)
  !      rdiag(nobsl)     : Observation error variance
  !      rloc(nobsl)      : Localization weight function
  !      dobl(nobsl)      : Innovation (y - Hxfmean)
  !
  !    OUTPUT
  !      wvec(nmem)       : Mean update (0 for LPF and LPFGM, r_size)
  !      Wmat(nmem,nmem)  : Transform matrix (r_size)
  !
  !=======================================================================

  SUBROUTINE lpfgm_core_noobs(wvec,Wmat)

    USE common_setting, only: r_size,nmem
    IMPLICIT NONE

    INTEGER i

    !---OUT
    REAL(r_size),INTENT(OUT) :: wvec(nmem)
    REAL(r_size),INTENT(OUT) :: Wmat(nmem,nmem)

    wvec(:)=REAL(0.d0,r_size)

    Wmat(:,:)=REAL(0.d0,r_size)
    DO i=1,nmem
       Wmat(i,i)=REAL(1.d0,r_size)
    END DO

  END SUBROUTINE lpfgm_core_noobs

  !-----------------------------------------------------------------------
  
  SUBROUTINE lpfgm_core(iseed_grd,nobsl,hdxfl,rdiag,rloc,dobl,wvec,Wmat)

    USE common_setting, only: r_size,r_sngl,r_dble,nmem,iswitch_da
    USE mpi
    IMPLICIT NONE

    !---MPI
    integer ierr
    
    !---WORK
    REAL(r_size) W_LPF(nmem,nmem)
    REAL(r_size) W_GM(nmem,nmem)
    REAL(r_size) wtmp(nmem)

    !---IN
    INTEGER,INTENT(IN)  :: iseed_grd
    INTEGER,INTENT(IN)  :: nobsl

    REAL(r_size),INTENT(IN)  :: rdiag(nobsl)
    REAL(r_size),INTENT(IN)  :: rloc(nobsl)
    REAL(r_size),INTENT(IN)  :: dobl(nobsl)

    !---OUT
    REAL(r_size),INTENT(OUT) :: wvec(nmem)
    REAL(r_size),INTENT(OUT) :: Wmat(nmem,nmem)

    !---INOUT
    REAL(r_size),INTENT(INOUT) :: hdxfl(nobsl,nmem) !Overwritten by gm_core (after lpf_core)
    
    !Initialization
    wvec(:)=REAL(0.0d0,r_size) !*wvec=0 in both LPF and LPFGM
    Wmat(:,:)=REAL(0.0d0,r_size)

    IF(iswitch_da == 2)THEN      !LPF
       
       CALL lpf_core(iseed_grd, nobsl, hdxfl, rdiag, rloc, dobl, wtmp, Wmat)
       
    ELSE IF(iswitch_da == 3)THEN !LPFGM
       
       CALL lpf_core(iseed_grd, nobsl, hdxfl, rdiag, rloc, dobl, wtmp, W_LPF)
       CALL gm_core(nobsl, hdxfl, rdiag, rloc, dobl, wtmp, W_GM)
       
       IF(r_size == r_dble)THEN
          CALL dgemm('n','n',nmem,nmem,nmem, &
               & 1.d0,W_GM,nmem,W_LPF,nmem, &
               & 0.d0,Wmat,nmem)
       ELSE IF(r_size == r_sngl)THEN
          CALL sgemm('n','n',nmem,nmem,nmem, &
               & 1.e0,W_GM,nmem,W_LPF,nmem, &
               & 0.e0,Wmat,nmem)
       END IF
             
    ELSE
       write(*,'(a,i6)') "***Error: Incorrect iswitch_da ==> ",iswitch_da
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    END IF

  END SUBROUTINE lpfgm_core

  !-----------------------------------------------------------------------
  !  Integer hash for random seed (32-bit integer mixing)
  !    Nearby inputs (e.g., neighboring grid points, Monte-Carlo index)
  !    are mapped to uncorrelated seeds in [0, 2^31-1].
  !    rand_lecuyer is a linear generator, so nearby seeds give
  !    correlated sequences without this mixing.
  !-----------------------------------------------------------------------

  INTEGER FUNCTION hash_seed(n)

    IMPLICIT NONE

    !---Parameter
    INTEGER(8),PARAMETER :: mask32=4294967295_8 !2^32-1

    !---Common
    INTEGER(8) x

    !---IN
    INTEGER,INTENT(IN) :: n

    x=IAND(INT(n,8),mask32)
    x=IEOR(x,ISHFT(x,-16))
    x=IAND(x*2146121005_8,mask32)  !x < 2^32, constant < 2^31 => no overflow
    x=IEOR(x,ISHFT(x,-15))
    x=IAND(x*1540483477_8,mask32)
    x=IEOR(x,ISHFT(x,-16))

    hash_seed=INT(IAND(x,2147483647_8)) !Positive 31-bit integer

  END FUNCTION hash_seed
  
END MODULE common_lpfgm
