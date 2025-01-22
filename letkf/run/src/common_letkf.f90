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
  !  01/04/2024 Shun OHISHI       modified for sbPOM-LETKF v2
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

    USE common_setting,only: r_size,r_dble
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
  !  Main Subroutine of LETKF Core
  !   INPUT
  !     nobsl            : Number of assimilated observations at a model grid point
  !     hdxf_in(nobs,nbv)   : Forecast ensemble perturbation in obs. space (=dYf)
  !     rdiag_in(nobs)      : Observation error variance (=sigma_o^2)
  !     rloc_in(nobs)       : Localization weighting function
  !     dep_in(nobs)        : Innovation (=y-Hxfmean)
  !     parm_infl        : Multiplicative inflation parameter
  !
  !   OUTPUT
  !     pa(nbv,nbv)   : Ensemble analysis covariance matrix Pa
  !     wvec(nbv)     : w vector (Update ensemble mean)
  !     Wmat(nbv,nbv) : Transform matrix
  !=======================================================================

  SUBROUTINE letkf_core_noobs(pa,wvec,Wmat)

    USE common_setting
    IMPLICIT NONE

    INTEGER i

    !---OUT
    REAL(r_size),INTENT(OUT) :: pa(nbv,nbv)
    REAL(r_size),INTENT(OUT) :: wvec(nbv)
    REAL(r_size),INTENT(OUT) :: Wmat(nbv,nbv)

    !Wmat = diag(rho^1/2)
    Wmat(:,:) = 0.0d0
    DO i=1,nbv
       Wmat(i,i) = SQRT(cov_infl_mul)
    END DO

    !wvec
    wvec(:)=0.0d0

    !Pa = [(nbv-1)/rho]^-1
    pa(:,:)=0.0d0
    DO i=1,nbv
       pa(i,i)=cov_infl_mul/REAL(nbv-1,r_size)
    END DO

  END SUBROUTINE letkf_core_noobs

  !-------------------------------

  SUBROUTINE letkf_core(nobsl,hdxf,rdiag,rloc,dep, &
       & pa,wvec,Wmat)

    USE common_setting
    IMPLICIT NONE

    INTEGER i,j
    INTEGER np

    REAL(r_dble) :: Rlinv_dYf(nobsl,nbv) !Rloc^-1_dYf
    REAL(r_dble) :: evec(nbv,nbv)        !Eigen vector
    REAL(r_dble) :: eval(nbv)            !Eigen value
    REAL(r_dble) :: work1(nbv,nbv)
    REAL(r_dble) :: work2(nbv,nobsl)
    REAL(r_dble) :: work3(nbv,nbv)

    !---IN
    INTEGER,INTENT(IN) :: nobsl    
    REAL(r_dble),INTENT(IN) :: hdxf(nobsl,nbv)
    REAL(r_dble),INTENT(IN) :: rdiag(nobsl)
    REAL(r_dble),INTENT(IN) :: rloc(nobsl)
    REAL(r_dble),INTENT(IN) :: dep(nobsl)

    !---OUT
    REAL(r_size),INTENT(OUT) :: pa(nbv,nbv)
    REAL(r_size),INTENT(OUT) :: wvec(nbv)
    REAL(r_size),INTENT(OUT) :: Wmat(nbv,nbv)
    
    !-----------------------------------------------------------------------
    !  Rloc^-1 dYf: [nobsl*nbv]
    !-----------------------------------------------------------------------

    DO j=1,nbv
       DO i=1,nobsl
          Rlinv_dYf(i,j) = hdxf(i,j) / rdiag(i) * rloc(i)
       END DO
    END DO

    !-----------------------------------------------------------------------
    !  dYf^T Rloc^-1 dYf [nbv*nbv]|
    !------------------------------
    !
    ! dgemm('n(A)','n(B)',m,n,k, alpha,A,m,B,k, beta,C,m)
    ! C <= alpha*(AB)+beta*C
    !
    ! A: m*k size
    ! B: k*n size
    ! C: m*n size
    !
    !-----------------------------------------------------------------------

    work1(:,:)=0.d0
    CALL dgemm('t','n',nbv,nbv,nobsl, &
         & 1.0d0,Rlinv_dYf,nobsl,hdxf,nobsl, &
         & 0.0d0,work1,nbv)

    !  DO j=1,nbv
    !    DO i=1,nbv
    !      work1(i,j) = Rlinv_dYf(1,i) * hdxf(1,j)
    !      DO k=2,nobsl
    !        work1(i,j) = work1(i,j) + Rlinv_dYf(k,i) * hdxf(k,j)
    !      END DO
    !    END DO
    !  END DO

    !-----------------------------------------------------------------------
    !  dYf^T Rloc^-1 dYf + (m-1) I / rho [nbv*nbv]
    !-----------------------------------------------------------------------

    DO i=1,nbv
       work1(i,i) = work1(i,i) + REAL(nbv-1,r_dble)/REAL(cov_infl_mul,r_dble)
    END DO

    !       do j=1,nbv
    !          do i=1,nbv
    !             WRITE(6,*) work1(i,j),hdxf(i,j),rdiag(i),rloc(i)
    !          end do
    !       end do

    !-----------------------------------------------------------------------
    !  Eigenvalues & vectors of [ dYf^T Rloc^-1 dYf + (m-1)I/rho]
    !-----------------------------------------------------------------------

    CALL mtx_eigen(nbv,work1,eval,evec,np)
    IF(myrank == 0 .and. np /= nbv)THEN
       WRITE(6,'(A,I10,2F12.2)') "positive eigenvalues:",np,eval(1),eval(nbv)
    END IF

    !-----------------------------------------------------------------------
    !  Pa = [ dYf^T Rloc^-1 dYf + (m-1) I/rho ]^-1 [nbv*nbv]
    !-----------------------------------------------------------------------

    !work1 = Evec * Eval^-1
    DO j=1,nbv
       DO i=1,nbv
          work1(i,j) = evec(i,j) / eval(j)
       END DO
    END DO

    !Pa = Evec * Eval^-1 * Evec^T
    work3(:,:)=0.d0
    CALL dgemm('n','t',nbv,nbv,nbv, &
         & 1.0d0,work1,nbv,evec,nbv, &
         & 0.0d0,work3,nbv)
    pa(:,:)=REAL(work3(:,:),r_size)

    !  DO j=1,nbv
    !    DO i=1,nbv
    !      pa(i,j) = work1(i,1) * evec(j,1)
    !      DO k=2,nbv
    !        pa(i,j) = pa(i,j) + work1(i,k) * evec(j,k)
    !      END DO
    !    END DO
    !  END DO

    !-----------------------------------------------------------------------
    !  work2 = Pa  (Rloc^-1 dYf)^T [nbv*nobsl]
    !-----------------------------------------------------------------------

    work2(:,:)=0.d0
    CALL dgemm('n','t',nbv,nobsl,nbv, &
         & 1.0d0,work3,nbv,Rlinv_dYf,nobsl, &
         & 0.0d0,work2,nbv)

    !  DO j=1,nobsl
    !    DO i=1,nbv
    !      work2(i,j) = pa(i,1) * Rlinv_dYf(j,1)
    !      DO k=2,nbv
    !        work2(i,j) = work2(i,j) + pa(i,k) * Rlinv_dYf(j,k)
    !      END DO
    !    END DO
    !  END DO

    !-----------------------------------------------------------------------
    !  wvec = Pa dYf^T Rloc^-1 (y-Hxfmean) [nbv]
    !-----------------------------------------------------------------------

    wvec(:)=0.d0
    DO i=1,nbv
       wvec(i) = work2(i,1) * dep(1)
       DO j=2,nobsl
          wvec(i) = wvec(i) + work2(i,j) * dep(j)
       END DO
    END DO

    !-----------------------------------------------------------------------
    !  Wmat = SQRT[(m-1)Pa]
    !-----------------------------------------------------------------------

    !Evec * Eval^-1/2
    DO j=1,nbv
       DO i=1,nbv
          work1(i,j) = evec(i,j) * SQRT( REAL(nbv-1,r_dble) / eval(j) )
       END DO
    END DO

    !Evec * Eval^-1/2 * Evec
    work3(:,:)=0.d0
    CALL dgemm('n','t',nbv,nbv,nbv, &
         & 1.0d0,work1,nbv,evec,nbv, &
         & 0.0d0,work3,nbv)
    Wmat(:,:)=REAL(work3(:,:),r_size)

    !  DO j=1,nbv
    !    DO i=1,nbv
    !      Wmat(i,j) = work1(i,1) * evec(j,1)
    !      DO k=2,nbv
    !        Wmat(i,j) = Wmat(i,j) + work1(i,k) * evec(j,k)
    !      END DO
    !    END DO
    !  END DO

  END SUBROUTINE letkf_core

END MODULE common_letkf
