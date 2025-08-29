MODULE common_lpfgm
  !=======================================================================
  !
  ! [PURPOSE:] Local Particle Filter Gaussian Mixture (LPF-GM)
  !            Model Independent Core Module
  !
  ! [REFERENCES:]
  ! [1] Kotsuki et al., 2022: A local particle filter and its Gaussian mixture
  !  extension implemented with minor modifications to the LETKF.
  !  Geosci. Model Dev., 15, 8325–8348. 
  !
  ! [HISTORY:]
  !  08/21/2025 Naoki Sakai     created based on the codes from sbPOM-LETKF v2 and speedy.
  !
  !=======================================================================


CONTAINS

  !=======================================================================
  !  Eigenvalue decomposition using subroutine rs
  !    INPUT
  !      n:      Matrix dimension
  !      A(n,n): Input Matrix
  !
  !    OUTPUT
  !      eval(n) :    Eigenvalue in decending order (i.e. eval(1) is the largest)
  !      evec(n,n) : Eivenvectors
  !      np           : Number of Positive eigenvalues
  !=======================================================================

  SUBROUTINE mtx_eigen(n,A,eval,evec,np)

    USE common_setting,only: r_size, r_dble
    USE common_mpi
    IMPLICIT NONE

    !---Common
    INTEGER i
    INTEGER lwork,info

    REAL(r_dble) eval8(n),evec8(n,n)
    REAL(r_dble) work(3*n-1)

    !---IN
    INTEGER,INTENT(IN) :: n
    REAL(r_dble),INTENT(IN) :: A(n,n)

    !---OUT
    REAL(r_dble),INTENT(OUT) :: eval(n),evec(n,n)
    INTEGER,INTENT(OUT) :: np

    lwork=3*n-1
    evec8(:,:)=REAL(A(:,:),r_dble)

    CALL DSYEV('V','U',n,evec8,n,eval8,work,lwork,info)
    !eval8: Ascending order (i.e., eval8(1) is the smallest)
    !CALL rs(n,n,a8,eival8,imode,eivec8,wrk1,wrk2,ierr) !***Eigen decomposition by netlib.f

    IF(info /= 0)THEN
       WRITE(6,'(A,I10)') "***Error: DSYEV, INFO: ",info
       WRITE(6,'(a,3f12.5)') "A:",A(1,1:3)
       WRITE(6,'(a,3f12.5)') "A:",A(2,1:3)
       WRITE(6,'(a,3f12.5)') "A:",A(3,1:3)
       WRITE(6,'(a,3f12.5)') "Evec:",evec8(1,1:3)
       WRITE(6,'(a,3f12.5)') "Evec:",evec8(2,1:3)
       WRITE(6,'(a,3f12.5)') "Evec:",evec8(3,1:3)
       CALL finalize_mpi
       STOP
    END IF

    np = n
    IF(eval8(n) < 0)THEN
       WRITE(6,'(A)') "***Error (mtx_eigen): All Eigenvalues are Negative"
       CALL finalize_mpi
       STOP
    ELSE
       !Reconditioning
       DO i=1,n
          IF(eval8(i) < 0.d0)THEN
             np = np - 1
             eval8(i) = 0.0d0
             evec8(:,i) = 0.0d0
          END IF
       END DO
    END IF

    DO i=1,n
       eval(i) = REAL(eval8(n+1-i), r_dble)
       evec(:,i) = REAL(evec8(:,n+1-i), r_dble)
    END DO

  END SUBROUTINE mtx_eigen

  !=======================================================================
  !  Random number generation with optional seeding (wrapper)
  !    INPUT
  !      ndim         : Number of random numbers to generate
  !      iseed        : Seed; 0 => use DATE_AND_TIME to seed, otherwise use provided
  !
  !    OUTPUT
  !      var(ndim)    : Uniform [0,1) random numbers (r_size)
  !
  !    NOTES
  !      - Uses Mersenne Twister interface: init_gen_rand, genrand_res53 (r_dble).
  !      - Computation is r_dble; results are stored as r_size.
  !=======================================================================

  SUBROUTINE com_rand_seed(ndim,iseed,var) ! from KK 20200426

    USE common_setting, only: r_dble
    USE mt19937, only: init_genrand, genrand_res53
    IMPLICIT NONE

    !---IN
    INTEGER,INTENT(IN) :: ndim, iseed

    !---OUT
    REAL(r_dble),INTENT(OUT) :: var(1:ndim)

    !---WORK
    INTEGER :: idate(8)
    INTEGER :: i, jseed
    LOGICAL,SAVE :: first=.true.
 

    IF (first) THEN
      !!!print *, first, iseed
      IF( iseed==0 ) THEN
        CALL DATE_AND_TIME(VALUES=idate)
        jseed = idate(8) + idate(7)*1000
      ELSE
        jseed = iseed
      ENDIF
      CALL init_genrand(jseed)
      first=.false.
    END IF

    DO i=1,ndim
      var(i) = genrand_res53()
    END DO

  END SUBROUTINE com_rand_seed


  !=======================================================================
  !  Quick sort (ascending) with index tracking
  !    INPUT
  !      var(*)       : Array to sort (r_dble, will be sorted in place)
  !      init(*)      : Index array (will be permuted accordingly)
  !      first, last  : Sorting range [first,last]
  !
  !    OUTPUT
  !      var(*), init(*) are modified in place (no separate outputs)
  !
  !    NOTES
  !      - Used to sort random numbers for multinomial resampling (MR).
  !      - Calculation uses r_dble; no r_size outputs here.
  !=======================================================================
  
  RECURSIVE SUBROUTINE quick_sort_asnd(var,init,first,last)

    USE common_setting, only: r_dble
    IMPLICIT NONE

    !---IN/OUT
    INTEGER :: first, last
    INTEGER :: init(*)
    REAL(r_dble) :: var(*)

    !---WORK
    INTEGER :: i, j, it
    REAL(r_dble) :: x, t

    x = var( (first+last) / 2 )
    i = first  ;  j = last
    do
      do while (var(i) < x)  ;  i=i+1  ;  end do
      do while (x < var(j))  ;  j=j-1  ;  end do
      if (i >= j) exit
        t       = var(i)  ;  var(i)  = var(j)  ;  var(j)  = t
        it      = init(i) ;  init(i) = init(j) ;  init(j) = it
        i       = i + 1   ;  j       = j - 1
    end do
    if (first < i - 1 ) call quick_sort_asnd(var,init, first, i - 1)
    if (j + 1 < last  ) call quick_sort_asnd(var,init, j + 1, last )

  end subroutine quick_sort_asnd


  !-----------------------------------------------------------------------
  !  Generate resampling matrix from cumulative weights
  !
  !  INPUT
  !    CC            : Resampling method code [CHAR(2)]
  !                    'SU' = systematic resampling (stratified via 1 offset)
  !                    'MR' = multinomial resampling (random + sort)
  !    DG            : Fill style [CHAR(2)]
  !                    'ON' = diagonal-oriented (prefer unique diagonal hits)
  !                    'OF' = off-diagonal allowed directly
  !    nbv           : Ensemble size [INTEGER]
  !    acc(1:nbv)    : Cumulative weights in [0,1], non-decreasing [REAL(r_dble)]
  !
  !  OUTPUT
  !    pmat(nbv,nbv) : Resampling matrix with 0/1 entries [REAL(r_dble)]
  !-----------------------------------------------------------------------

  SUBROUTINE get_resampling_mtx(CC,DG,nbv,acc,pmat)
   
    USE common_setting, only: r_dble
    IMPLICIT NONE

    !---IN
    INTEGER     , INTENT(in)  :: nbv
    CHARACTER(2), INTENT(in)  :: CC, DG
    REAL(r_dble), INTENT(in)  :: acc(1:nbv)

    !---OUT
    REAL(r_dble), INTENT(out) :: pmat(nbv,nbv)

    !---WORK
    INTEGER      :: i, j, k, init(nbv), inum(nbv)
    REAL(r_dble) :: rand(nbv), temp

    DO j=1,nbv ; init(j) = j  ; END DO
    !CALL com_rand(nbv,rand) ! [0-1]

    IF( CC == 'SU' )THEN
      CALL com_rand_seed(nbv,0,rand) ! [0-1]
      temp = rand(1)
      DO i=1,nbv
        rand(i) = ( dble(i) +  temp - 1.0d0 ) / dble(nbv)
      END DO
    ELSE IF( CC == 'MR' )THEN
      CALL com_rand_seed(nbv,0,rand) ! [0-1]
      CALL quick_sort_asnd(rand,init,1,nbv)
    ELSE
      PRINT *, "  ERROR&SROP :: NOT SUCH OPTION in get_resampling_mtx for CC ", CC
      STOP
    END IF

    !==> generate resampling mxm matrix (KK's diagonal-oriented; str)
    IF( DG =='ON' )THEN
      inum(:)   = 0
      !(1) :: gen selected particle
      DO j=1,nbv ! jth column
        DO i=1,nbv-1
          IF( rand(j)<=acc(i) )THEN
            inum(j) = i
            EXIT 
          ENDIF
        END DO
        IF (inum(j) == 0) inum(j) = nbv
      END DO
      !!print *, inum(:)

      pmat(:,:) = 0.0d0
      DO i=1,nbv ! (1) diagonal component
        k = inum(i)
        IF( k/=0 .and. pmat(k,k)<0.5d0 )THEN
          pmat(k,k) = 1.0d0
          inum(i)   = 0
        END IF 
      END DO

      DO i=1,nbv ! (2) off-diagonal component
        k = inum(i)
        IF( k/=0 )THEN
          DO j=1,nbv
            IF( sum(pmat(:,j))<0.5 )THEN
              pmat(k,j) = 1.0d0
              EXIT
            END IF
          END DO
        END IF 
      END DO

    ELSE IF( DG =='OF' )THEN  
      !==> generate resampling mxm matrix (non-diagonal-oriented; str, SK's trial)
      pmat(:,:) = 0.0d0
      DO j=1,nbv ! jth column
        DO i=1,nbv-1
          IF( rand(j)<=acc(i) )THEN
            pmat(i,j) = 1.0d0
            EXIT
          END IF
        END DO
        IF (sum(pmat(:,j)) < 0.5) pmat(nbv,j) = 1.0d0
      END DO
    ENDIF

  END SUBROUTINE get_resampling_mtx

 
  !-----------------------------------------------------------------------
  !  Compute normalized particle filter weights (Gaussian likelihood)
  !    INPUT
  !      nobsl          : Number of assimilated observations
  !      nbv            : Ensemble size
  !      dep(nobsl)     : Innovation (y - H xf_mean) [r_dble]
  !      hdxf(nobsl,nbv): Ensemble perturbation in obs space (= H xf_i - H xf_mean) [r_dble]
  !      rloc(nobsl)    : Localization weights [r_dble]
  !      rdiag(nobsl)   : Observation error variance [r_dble]
  !
  !    OUTPUT
  !      pfwgh(nbv)     : Normalized PF weights (r_dble)
  !
  !    NOTES
  !      - Computations in r_dble to avoid underflow; normalization at end.
  !-----------------------------------------------------------------------

  SUBROUTINE calc_pfwgh_norml(nobsl,nbv,dep,hdxf,rloc,rdiag,pfwgh)

    USE common_setting, only: r_dble
    IMPLICIT NONE
    !---IN
    INTEGER     , INTENT(IN)  :: nobsl, nbv
    REAL(r_dble), INTENT(IN)  :: dep(1:nobsl), hdxf(1:nobsl,1:nbv), rloc(1:nobsl), rdiag(1:nobsl)

    !---OUT
    REAL(r_dble), INTENT(OUT) :: pfwgh(1:nbv) 

    !---WORK (r_dble calculations)
    INTEGER :: i, j, k
    REAL(r_dble) :: sqpf, qtmp
    REAL(r_dble) :: qpf(1:nbv), dep2_Ri(1:nbv)

    dep2_Ri(:) = 0.0d0
    DO j=1,nbv
      DO i=1,nobsl
        dep2_Ri(j) = dep2_Ri(j) - 0.5d0*((hdxf(i,j)-dep(i))**2.0d0) * rloc(i) / rdiag(i)
      END DO
    END DO

    sqpf     = 0.0d0
    qpf(:)   = 0.0d0
    DO j=1,nbv
      qtmp = 0.0d0
      DO k=1,nbv
        qtmp = qtmp + dexp( dep2_Ri(k) - dep2_Ri(j) )
      END DO
      qpf(j) = 1.0d0 / qtmp
      sqpf   = sqpf + qpf(j)
    END DO

    pfwgh(:) = qpf(:) / sqpf

  END SUBROUTINE calc_pfwgh_norml


  !=======================================================================
  !  Main Subroutine of LPF Core
  !   INPUT
  !     nobsl            : Number of assimilated observations at a model grid point
  !     hdxf(nobs,nbv)      : Forecast ensemble perturbation in obs. space (=dYf)
  !     rdiag(nobs)         : Observation error variance (=sigma_o^2)
  !     rloc(nobs)          : Localization weighting function
  !     dep(nobs)           : Innovation (=y-Hxfmean)
  !
  !   OUTPUT
  !     wvec(nbv)           : w vector (Update ensemble mean, zero vector)
  !     Wmat(nbv,nbv)       : Transform matrix
  !=======================================================================

  SUBROUTINE lpf_core(nobsl,hdxf,rdiag,rloc,dep,wvec,Wmat)

    USE common_setting
    IMPLICIT NONE

    !---IN
  !   INTEGER,INTENT(IN) :: nobs
    INTEGER, INTENT(IN) :: nobsl
    REAL(r_dble), INTENT(IN) :: hdxf(nobsl, nbv)
    REAL(r_dble), INTENT(IN) :: rdiag(nobsl)
    REAL(r_dble), INTENT(IN) :: rloc(nobsl)
    REAL(r_dble), INTENT(IN) :: dep(nobsl)

    !---OUT
    REAL(r_size), INTENT(OUT) :: wvec(nbv)
    REAL(r_size), INTENT(OUT) :: Wmat(nbv, nbv)

    !---WORK
    REAL(r_dble) :: pfwgh(nbv)
    REAL(r_dble) :: acc(nbv), pmat(nbv, nbv)
    Real(r_dble) :: swgh, peff
    INTEGER :: j
  !-----------------------------------------------------------------------    
  ! Likelihood Computation based on Gauss
  !-----------------------------------------------------------------------
    ! dep   :: yo     - Hxf(mean)
    ! hdxf  :: Hxf(i) - Hxf(mean)
    ! rdiag :: err*err (i.e., variance)
    ! rloc  :: 0-1

    CALL calc_pfwgh_norml(nobsl,nbv,dep,hdxf,rloc,rdiag,pfwgh)
    swgh     = sum( pfwgh(:) )
    pfwgh(:) = pfwgh(:) / swgh

    peff     = 1.0d0  / sum( pfwgh(:)**2.0d0 ) ! effective particle size

    acc(:)   = 0.0d0
    acc(1)   = pfwgh(1)
    DO j=2,nbv
      acc(j) = acc(j-1) + pfwgh(j)
    END DO
    
  !-----------------------------------------------------------------------    
  ! Resampling with random numbers
  !-----------------------------------------------------------------------
    Wmat(:,:) = 0.0d0
    DO j = 1, nmonte
      CALL get_resampling_mtx('MR', 'ON', nbv, acc, pmat)
      Wmat(:,:) = Wmat(:,:) + pmat(:,:) / dble(nmonte)
    END DO
  !-----------------------------------------------------------------------
  ! No ensemble mean update in LPF
  !-----------------------------------------------------------------------s
    wvec(:) = 0.0d0

  END SUBROUTINE lpf_core


  !=======================================================================
  !  Main Subroutine of GM Core
  !   INPUT
  !     nobsl              : Number of assimilated observations at a model grid point
  !     hdxf(nobsl,nbv) : Forecast ensemble perturbation in obs. space (=dYf)
  !     rdiag(nobsl)    : Observation error variance (=sigma_o^2)
  !     rloc(nobsl)     : Localization weighting function
  !     dep(nobsl)      : Innovation (=y-Hxfmean)
  !
  !   OUTPUT
  !     pa(nbv,nbv)        : Ensemble analysis covariance matrix Pa
  !     wvec(nbv)          : w vector (Update ensemble mean) → set to 0 for LPFGM
  !     Wmat(nbv,nbv)      : Transform matrix (LPFGM shift matrix + I)
  !=======================================================================

  SUBROUTINE gm_core(nobsl, hdxf, rdiag, rloc, dep, wvec, Wmat)

    USE common_setting
    IMPLICIT NONE

    INTEGER :: i, j, k, np

    !---IN
    INTEGER, INTENT(IN) :: nobsl
    REAL(r_dble), INTENT(IN) :: hdxf(nobsl, nbv)
    REAL(r_dble), INTENT(IN) :: rdiag(nobsl)
    REAL(r_dble), INTENT(IN) :: rloc(nobsl)
    REAL(r_dble), INTENT(IN) :: dep(nobsl)

    !---OUT
    ! REAL(r_size), INTENT(OUT) :: pa(nbv, nbv)
    REAL(r_size), INTENT(OUT) :: wvec(nbv)
    REAL(r_size), INTENT(OUT) :: Wmat(nbv, nbv)

    !---WORK
    REAL(r_dble) :: Rlinv_dYf(nobsl, nbv)
    REAL(r_dble) :: evec(nbv, nbv)
    REAL(r_dble) :: eval(nbv)
    REAL(r_dble) :: work1(nbv, nbv)
    REAL(r_dble) :: work2(nbv, nobsl)
    REAL(r_dble) :: work3(nbv, nbv)
    REAL(r_dble) :: work4(nbv, nbv)

    !-----------------------------------------------------------------------
    ! Rloc^-1 * dYf
    !-----------------------------------------------------------------------
    DO j = 1, nbv
      DO i = 1, nobsl
        Rlinv_dYf(i,j) = hdxf(i,j) / rdiag(i) * rloc(i)
      END DO
    END DO

    !-----------------------------------------------------------------------
    ! dYf^T * Rloc^-1 * dYf + (m-1)/rho * I
    !-----------------------------------------------------------------------
    CALL dgemm('t','n',nbv,nbv,nobsl, 1.0d0, Rlinv_dYf, nobsl, hdxf, nobsl, 0.0d0, work1, nbv)

    DO i = 1, nbv
      work1(i,i) = work1(i,i) + REAL(nbv - 1, r_dble) / REAL(cov_infl_gm, r_dble)
    END DO

    !-----------------------------------------------------------------------
    ! Eigen decomposition
    !-----------------------------------------------------------------------
    CALL mtx_eigen(nbv, work1, eval, evec, np)

    !-----------------------------------------------------------------------
    ! Pa = Evec * Eval^-1 * Evec^T
    !-----------------------------------------------------------------------
    DO j = 1, nbv
      DO i = 1, nbv
        work1(i,j) = evec(i,j) / eval(j)
      END DO
    END DO

    CALL dgemm('n','t',nbv,nbv,nbv, 1.0d0, work1, nbv, evec, nbv, 0.0d0, work3, nbv)
    ! pa(:,:) = REAL(work3(:,:), r_size)

    !-----------------------------------------------------------------------
    ! work2 = Pa * (Rloc^-1 * dYf)^T
    !-----------------------------------------------------------------------
    CALL dgemm('n','t',nbv,nobsl,nbv, 1.0d0, work3, nbv, Rlinv_dYf, nobsl, 0.0d0, work2, nbv)

    !-----------------------------------------------------------------------
    ! work4 = work2 * (dep - dYf) + I
    !-----------------------------------------------------------------------
    DO k = 1, nbv
      DO i = 1, nbv
        work4(i,k) = work2(i,1) * (dep(1) - hdxf(1,k))
        DO j = 2, nobsl
          work4(i,k) = work4(i,k) + work2(i,j) * (dep(j) - hdxf(j,k))
        END DO
      END DO
    END DO

    DO i = 1, nbv
      work4(i,i) = work4(i,i) + 1.0d0
    END DO

    Wmat(:,:) = REAL(work4(:,:), r_size)

    !-----------------------------------------------------------------------
    ! wvec = 0 (not used in LPFGM)
    !-----------------------------------------------------------------------
    wvec(:) = 0.0d0

  END SUBROUTINE gm_core


  !=======================================================================
  !  LPFGM main driver: compose GM and LPF transforms
  !    INPUT
  !      nobsl          : Number of assimilated observations
  !      hdxf(nobsl,nbv): Ensemble perturbation in obs space (=dYf) [r_dble]
  !      rdiag(nobsl)   : Observation error variance [r_dble]
  !      rloc(nobsl)    : Localization weights [r_dble]
  !      dep(nobsl)     : Innovation (y - H xf_mean) [r_dble]
  !
  !    OUTPUT
  !      wvec(nbv)      : Mean update (0 for LPFGM, r_size)
  !      Wmat(nbv,nbv)  : Transform matrix (r_size)
  !
  !    NOTES
  !      - Inputs and calculations are r_dble; outputs are r_size.
  !=======================================================================
  
  SUBROUTINE lpfgm_core(nobsl, hdxf, rdiag, rloc, dep, wvec, Wmat)
    
    USE common_setting, only: r_size, r_dble, nbv, iswitch_da
    IMPLICIT NONE
    !---IN
    INTEGER     , INTENT(IN)  :: nobsl
    REAL(r_dble), INTENT(IN)  :: hdxf(nobsl, nbv)
    REAL(r_dble), INTENT(IN)  :: rdiag(nobsl)
    REAL(r_dble), INTENT(IN)  :: rloc(nobsl)
    REAL(r_dble), INTENT(IN)  :: dep(nobsl)

    !---OUT
    REAL(r_size), INTENT(OUT) :: wvec(nbv)
    REAL(r_size), INTENT(OUT) :: Wmat(nbv, nbv)

    !---WORK
    REAL(r_size) :: W_LPF(nbv, nbv), W_GM(nbv, nbv), wtmp(nbv)

    ! LPF transform
    CALL lpf_core(nobsl, hdxf, rdiag, rloc, dep, wtmp, W_LPF)

    IF (iswitch_da == 2) THEN
      Wmat(:,:) = W_LPF(:,:)
      wvec(:) = 0.0d0
      RETURN
    ELSE IF (iswitch_da == 3) THEN
      ! GM transform
      CALL gm_core(nobsl, hdxf, rdiag, rloc, dep, wtmp, W_GM)
      Wmat(:,:) = MATMUL(W_GM, W_LPF)
    ELSE
      Wmat(:,:) = 0.0d0
    END IF

    ! No mean update for LPFGM
    wvec(:) = 0.0d0

  END SUBROUTINE lpfgm_core

END MODULE common_lpfgm
