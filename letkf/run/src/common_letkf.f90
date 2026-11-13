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
  !  2009/01/21 Takemasa Miyoshi      Created at U. of Maryland, College Park
  !  2024/01/04 Shun OHISHI           Modified for sbPOM-LETKF v2
  !  2025/10/16 Shun OHISHI           Change precision
  !  2026/10/28 Shun OHISHI (Claude)  Update all (e.g., including nobsl < nmem case)
  !=======================================================================

CONTAINS

  !=======================================================================
  !  LETKF core: Efficient ensemble update
  !=======================================================================
  !
  !  [IN]    nobsl            : Local observation size (p)
  !          rdiag(nobsl)     : Observation error variance
  !          rloc(nobsl)      : Localization weight
  !          dobl(nobsl)      : Innovation y - H(xf_mean)
  !  [INOUT] hdxfl(nobsl,nmem): IN : dY = H(Xf) - H(xf_mean)
  !                             OUT: overwritten (work array)
  !  [OUT]   wvec(nmem)       : Pa dY^T Rloc^-1 [y - H(xf_mean)]
  !          Wmat(nmem,nmem)  : W = [(m-1) Pa]^1/2
  !
  !  [Notation] (m: nmem, p: nobsl)
  !    alpha = (m-1)/rho                         : rho = multiplicative inflation factor
  !    dYt   = Rloc^-1/2 dY          [p*m]       : delta Y tilde
  !    dt    = Rloc^-1/2 [y-H(xf_mean)] [p]      : innovation weighted by Rloc^-1/2
  !    Pa    = [alpha I + dYt^T dYt]^-1 [m*m]
  !
  !  (1) m <= p: Standard LETKF (eigenvalue decomposition of Painv [m*m])
  !        Painv = alpha I + dYt^T dYt = E diag(e) E^T
  !        wvec  = E diag(1/e) E^T dYt^T dt
  !        W     = SQRT(m-1) E diag(e^-1/2) E^T = C C^T,  C = E diag( ((m-1)/e)^(1/4) )
  !
  !  (2) p < m: Efficient update (eigenvalue decomposition of G [p*p])
  !        G     = alpha I + dYt dYt^T = U diag(d) U^T        (d >= alpha > 0)
  !        dYtU  = U^T dYt               [p*m]
  !        wvec  = dYtU^T diag(1/d) U^T dt                    (push-through identity)
  !        dYhU  = diag( [SQRT(d)(SQRT(alpha)+SQRT(d))]^-1/2 ) dYtU  [p*m]
  !        W     = SQRT(rho) [ I - dYhU^T dYhU ]
  !
  !  Cost: (1) O(m^3 + p m^2),  (2) O(p^3 + p^2 m + p m^2)
  !  Only the upper triangle of Painv, G, and W is computed by xSYRK,
  !  and W is symmetrized at the end.
  !-----------------------------------------------------------------------
  !  xSYRK (BLAS level 3): Symmetric rank-k update  (DSYRK: double, SSYRK: single)
  !
  !    CALL xSYRK(uplo, trans, n, k, s1, A, lda, s2, C, ldc)
  !
  !      trans = 'N':  C <- s1 * A   A^T + s2 * C    (A: n*k)
  !      trans = 'T':  C <- s1 * A^T A   + s2 * C    (A: k*n)
  !
  !      uplo  = 'U' (or 'L'): only the upper (or lower) triangle of C [n*n]
  !                            is computed and referenced
  !      n     : order of C
  !      k     : number of columns of A (trans='N') or rows of A (trans='T')
  !      s1,s2 : scalars (s2 = 0: C need not be initialized)
  !      lda   : leading dimension of A (= n for 'N', = k for 'T')
  !      ldc   : leading dimension of C (= n)
  !
  !    Since C is symmetric, only one triangle is computed: about half the
  !    operations of xGEMM (C <- A^T A). The other triangle is not set, so
  !    it must be filled by symmetry if the full matrix is needed (W).
  !    mtx_eigen (DSYEVD with 'U') uses only the upper triangle.
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

  !-------------------------------------------------------------------------
  
  SUBROUTINE letkf_core(nobsl,hdxfl,rdiag,rloc,dobl, &
       & wvec,Wmat)

    USE common_setting
    USE common
    IMPLICIT NONE

    !---Common
    INTEGER i,j
    INTEGER np                              !#Positive eigenvalues (from mtx_eigen)

    REAL(r_dble) alpha                      !(m-1)/rho
    REAL(r_size) sqrt_rho                   !SQRT(rho)

    REAL(r_size),ALLOCATABLE :: sqrt_rloc(:) !Rloc^-1/2 = SQRT(rloc/rdiag) [p]
    REAL(r_size),ALLOCATABLE :: dt(:)        !Rloc^-1/2 [y-H(xf_mean)]     [p]

    !---(1) m <= p
    REAL(r_size),ALLOCATABLE :: painv(:,:)  !Painv [m*m], then C
    REAL(r_dble),ALLOCATABLE :: emat(:,:)   !E [m*m]
    REAL(r_dble),ALLOCATABLE :: eval(:)     !e [m]
    REAL(r_dble),ALLOCATABLE :: tmp_m(:)    !dYt^T dt, then diag(1/e) E^T dYt^T dt [m]

    !---(2) p < m
    REAL(r_size),ALLOCATABLE :: gmat(:,:)   !G [p*p], then U (r_size)
    REAL(r_dble),ALLOCATABLE :: umat(:,:)   !U [p*p]
    REAL(r_dble),ALLOCATABLE :: dval(:)     !d [p]
    REAL(r_size),ALLOCATABLE :: dYtU(:,:)   !dYtU [p*m], then dYhU
    REAL(r_dble),ALLOCATABLE :: tmp_p(:)    !diag(1/d) U^T dt, then [SQRT(d)(SQRT(alpha)+SQRT(d))]^-1/2 [p]

    !---IN
    INTEGER,INTENT(IN) :: nobsl
    REAL(r_size),INTENT(IN) :: rdiag(nobsl)
    REAL(r_size),INTENT(IN) :: rloc(nobsl)
    REAL(r_size),INTENT(IN) :: dobl(nobsl)

    !---INOUT
    REAL(r_size),INTENT(INOUT) :: hdxfl(nobsl,nmem) !IN: dY, OUT: dYt (overwritten)

    !---OUT
    REAL(r_size),INTENT(OUT) :: wvec(nmem)
    REAL(r_size),INTENT(OUT) :: Wmat(nmem,nmem)

    alpha    = REAL(nmem-1,r_dble)/REAL(cov_infl_mul,r_dble)
    sqrt_rho = SQRT(cov_infl_mul)

    !=====================================================================
    ! dYt = Rloc^-1/2 dY (overwrite hdxfl), dt = Rloc^-1/2 [y-H(xf_mean)]
    !=====================================================================

    ALLOCATE(dt(nobsl))
    ALLOCATE(sqrt_rloc(nobsl))

    DO i=1,nobsl
       sqrt_rloc(i)=SQRT(rloc(i)/rdiag(i))
       dt(i)=sqrt_rloc(i)*dobl(i)            !Rloc^(-1/2) [y-H(xf_mean)]
    END DO

    DO j=1,nmem
       DO i=1,nobsl
          hdxfl(i,j)=hdxfl(i,j)*sqrt_rloc(i) !dYt = Rloc^(-1/2) dY
       END DO
    END DO

    DEALLOCATE(sqrt_rloc)

    IF(nmem <= nobsl)THEN

       !==================================================================
       ! (1) m <= p: Standard LETKF
       !==================================================================

       ALLOCATE(painv(nmem,nmem),emat(nmem,nmem),eval(nmem),tmp_m(nmem))

       !---Painv = alpha I + dYt^T dYt (upper triangle)
       !xSYRK('U','T',n=m,k=p, 1, dYt[p*m], lda=p, 0, Painv[m*m], ldc=m): Painv <- dYt^T dYt
       painv(:,:)=REAL(0.d0,r_size)
       IF(r_size == r_dble)THEN
          CALL dsyrk('U','T',nmem,nobsl,1.d0,hdxfl,nobsl,0.d0,painv,nmem) !dYt^T dYt
       ELSE IF(r_size == r_sngl)THEN
          CALL ssyrk('U','T',nmem,nobsl,1.e0,hdxfl,nobsl,0.e0,painv,nmem)
       END IF

       DO i=1,nmem
          painv(i,i)=painv(i,i)+REAL(alpha,r_size) !alpha I + dYt^T dYt
       END DO

       !---Painv = E diag(e) E^T
       CALL mtx_eigen(nmem,painv,eval,emat,np)

       !---wvec = E diag(1/e) E^T dYt^T dt
       DO j=1,nmem
          tmp_m(j)=DOT_PRODUCT(REAL(hdxfl(:,j),r_dble),REAL(dt(:),r_dble)) !dYt^T dt
       END DO
       tmp_m(:)=MATMUL(tmp_m(:),emat(:,:))/eval(:)                         !diag(1/e) E^T dYt^T dt
       wvec(:)=REAL(MATMUL(emat(:,:),tmp_m(:)),r_size)                     !E diag(1/e) E^T dYt^T dt

       !---W = E diag( sqrt( (m-1)/e) ) E^T = C C^T, C = E diag( ((m-1)/e)^(1/4) ) (painv is reused for C)
       DO j=1,nmem
          painv(:,j)=REAL(emat(:,j)*SQRT( SQRT(REAL(nmem-1,r_dble)/eval(j)) ),r_size) !C
       END DO

       !   xSYRK('U','N',n=m,k=m, 1, C[m*m], lda=m, 0, W[m*m], ldc=m): W <- C C^T       
       Wmat(:,:)=REAL(0.d0,r_size)
       IF(r_size == r_dble)THEN
          CALL dsyrk('U','N',nmem,nmem,1.d0,painv,nmem,0.d0,Wmat,nmem) !W = C C^T
       ELSE IF(r_size == r_sngl)THEN
          CALL ssyrk('U','N',nmem,nmem,1.e0,painv,nmem,0.e0,Wmat,nmem)
       END IF

       DEALLOCATE(painv,emat,eval,tmp_m)

    ELSE

       !==================================================================
       ! (2) p < m: Efficient update
       !==================================================================

       ALLOCATE(gmat(nobsl,nobsl),umat(nobsl,nobsl),dval(nobsl))
       ALLOCATE(dYtU(nobsl,nmem),tmp_p(nobsl))

       !---G = alpha I + dYt dYt^T (upper triangle)
       !xSYRK('U','N',n=p,k=m, 1, dYt[p*m], lda=p, 0, G[p*p], ldc=p): G <- dYt dYt^T
       gmat(:,:)=REAL(0.d0,r_size)
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
          CALL dgemm('t','n',nobsl,nmem,nobsl,1.d0,gmat,nobsl,hdxfl,nobsl,0.0d0,dYtU,nobsl) !dYtU = U^T dYt
       ELSE IF(r_size == r_sngl)THEN
          CALL sgemm('t','n',nobsl,nmem,nobsl,1.e0,gmat,nobsl,hdxfl,nobsl,0.0e0,dYtU,nobsl)
       END IF

       !---wvec = dYtU^T diag(1/d) U^T dt
       tmp_p(:)=MATMUL(REAL(dt(:),r_dble),umat(:,:))/dval(:) !diag(1/d) U^T dt
       DO j=1,nmem
          wvec(j)=REAL(DOT_PRODUCT(REAL(dYtU(:,j),r_dble),tmp_p(:)),r_size) !dYtU^T diag(1/d) U^T dt
       END DO

       !---dYhU = diag( [SQRT(d)(SQRT(alpha)+SQRT(d))]^-1/2 ) * dYtU (dYtU is reused for dYhU)
       DO i=1,nobsl
          tmp_p(i)=1.d0/SQRT( SQRT(dval(i)) * ( SQRT(alpha)+SQRT(dval(i)) ) ) !1/[ sqrt(d) * (sqrt(alpha)+sqrt(d)) ]
       END DO
       DO j=1,nmem
          DO i=1,nobsl
             dYtU(i,j)=dYtU(i,j)*REAL(tmp_p(i),r_size) !diag( [SQRT(d)(SQRT(alpha)+SQRT(d))]^-1/2 ) * dYtU
          END DO
       END DO

       !---W = SQRT(rho) [ I - dYhU^T dYhU ] (upper triangle)
       Wmat(:,:)=REAL(0.d0,r_size)
       IF(r_size == r_dble)THEN
          CALL dsyrk('U','T',nmem,nobsl,-1.d0*REAL(sqrt_rho,r_dble),dYtU,nobsl,0.d0,Wmat,nmem) ! - dYhU^T dYhU
       ELSE IF(r_size == r_sngl)THEN
          CALL ssyrk('U','T',nmem,nobsl,-1.e0*REAL(sqrt_rho,r_sngl),dYtU,nobsl,0.e0,Wmat,nmem)
       END IF
       DO i=1,nmem
          Wmat(i,i)=Wmat(i,i)+sqrt_rho ! - dYhU^T dYhU + SQRT(rho) 
       END DO

       DEALLOCATE(gmat,umat,dval)
       DEALLOCATE(dYtU,tmp_p)

    END IF

    DEALLOCATE(dt)

    !=====================================================================
    ! W: upper triangle --> lower triangle (symmetric)
    !=====================================================================

    DO j=1,nmem
       DO i=j+1,nmem
          Wmat(i,j)=Wmat(j,i)
       END DO
    END DO

  END SUBROUTINE letkf_core
  
END MODULE common_letkf
