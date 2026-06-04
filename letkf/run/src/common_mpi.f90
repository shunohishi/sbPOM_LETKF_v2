MODULE common_mpi
  !=======================================================================
  !
  ! [PURPOSE:] General MPI procedures
  !
  ! [HISTORY:]
  !   09/06/2005 Takemasa MIYOSHI  created
  !   07/31/2025 Shun OHISHI       added bcast
  !   06/04/2026 Shun OHISHI       updated
  !
  !=======================================================================

CONTAINS

  !----------------------------------------------------------------------

  SUBROUTINE initialize_mpi

    USE MPI    
    USE common_setting
    IMPLICIT NONE

    INTEGER :: ierr

    CALL MPI_INIT(ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

    IF(myrank == 0)THEN
       WRITE(6,'(A,I5.5)') "Number of Processor:",nprocs
    END IF

  END SUBROUTINE initialize_mpi

  !----------------------------------------------------------------------

  SUBROUTINE finalize_mpi

    USE MPI    
    IMPLICIT NONE

    INTEGER :: ierr

    CALL MPI_FINALIZE(ierr)

  END SUBROUTINE finalize_mpi

  !-----------------------------------------------------------------------
  ! Scatter gridded data to processes (nrank -> all)
  !-----------------------------------------------------------------------

  SUBROUTINE scatter_grd_mpi(nrank,v3dg,v2dg,v3d,v2d)

    !$USE OMP_LIB
    USE MPI
    USE common_setting
    USE common
    IMPLICIT NONE

    !---Common
    INTEGER j,k
    INTEGER n,n0
    INTEGER ivd
    INTEGER ierr

    REAL(r_sngl) bufs(nij1max,nlevall,nprocs)
    REAL(r_sngl) bufr(nij1max,nlevall)

    !---IN
    INTEGER,INTENT(IN) :: nrank
    REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
    REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)

    !---OUT
    REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
    REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)

    bufs(:,:,:) = 0.e0
    bufr(:,:) = 0.e0

    !v3dg and v2dgs --> bufs
    IF(myrank == nrank)THEN

       DO ivd=1,nv3d
          DO k=1,nlev

             j = (ivd-1)*nlev+k          
             CALL grd_to_buf(v3dg(:,:,k,ivd),bufs(:,j,:))
          END DO
       END DO

       DO ivd=1,nv2d
          j = nv3d*nlev+ivd
          CALL grd_to_buf(v2dg(:,:,ivd),bufs(:,j,:))
       END DO
    END IF

    !bufs --> bufr
    n = nij1max * nlevall
    n0 = n
    CALL MPI_SCATTER(bufs,n,MPI_REAL,&
         & bufr,n0,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

    !bufr --> v3d and v2d
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(ivd,k,j)
    DO ivd=1,nv3d
       DO k=1,nlev

          j = (ivd-1)*nlev+k          
          v3d(:,k,ivd) = REAL(bufr(1:nij1,j),r_size)

       END DO
    END DO
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(ivd,j)
    DO ivd=1,nv2d

       j = nv3d*nlev+ivd
       v2d(:,ivd) = REAL(bufr(1:nij1,j),r_size)

    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE scatter_grd_mpi

  !-----------------------------------------------------------------------
  ! Gather gridded data (all -> nrank)
  !-----------------------------------------------------------------------
  
  SUBROUTINE gather_grd_mpi(nrank,v3d,v2d,v3dg,v2dg)

    !$USE OMP_LIB    
    USE MPI
    USE common_setting
    USE common
    IMPLICIT NONE

    !---Common
    INTEGER j,k
    INTEGER n,n0
    INTEGER ivd
    INTEGER ierr

    REAL(r_sngl) bufs(nij1max,nlevall)
    REAL(r_sngl) bufr(nij1max,nlevall,nprocs)

    !---IN
    INTEGER,INTENT(IN) :: nrank
    REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
    REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)

    !---OUT
    REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
    REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)

    bufs(:,:) = 0.e0
    bufr(:,:,:) = 0.e0

    !bufs
    !$OMP PARALLEL PRIVATE(ivd,k,j)
    !$OMP DO COLLAPSE(2)
    DO ivd=1,nv3d
       DO k=1,nlev

          j = (ivd-1)*nlev + k
          bufs(1:nij1,j) = REAL(v3d(:,k,ivd),r_sngl)

       END DO
    END DO
    !$OMP END DO

    !$OMP DO
    DO ivd=1,nv2d
       j = nv3d*nlev + ivd
       bufs(1:nij1,j) = REAL(v2d(:,ivd),r_sngl)
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    !bufs --> bufr
    n = nij1max * nlevall
    n0 = n
    CALL MPI_GATHER(bufs,n,MPI_REAL,&
         & bufr,n0,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

    !bufr --> v3dg and v2dg
    IF(myrank == nrank)THEN

       DO ivd=1,nv3d
          DO k=1,nlev
             j = (ivd-1)*nlev + k
             CALL buf_to_grd(bufr(:,j,:),v3dg(:,:,k,ivd))
          END DO
       END DO

       DO ivd=1,nv2d
          j = nv3d*nlev + ivd
          CALL buf_to_grd(bufr(:,j,:),v2dg(:,:,ivd))
       END DO

    END IF

  END SUBROUTINE gather_grd_mpi

  !-----------------------------------------------------------------------
  ! Share nrank data with all processer
  !-----------------------------------------------------------------------

  SUBROUTINE bcast_mpi_1d(nrank,im,dat)

    USE MPI
    USE common_setting
    IMPLICIT NONE

    !---Common
    INTEGER MPI_R_SIZE,ierr

    !---IN
    INTEGER,INTENT(IN) :: nrank
    INTEGER,INTENT(IN) :: im

    !---INOUT
    REAL(r_size),INTENT(INOUT) :: dat(im)

    IF(r_size == kind(0.e0))THEN
       MPI_R_SIZE=MPI_REAL
    ELSE IF(r_size == kind(0.d0))THEN
       MPI_R_SIZE=MPI_DOUBLE_PRECISION
    END IF

    CALL MPI_BCAST(dat,im,MPI_R_SIZE,nrank,MPI_COMM_WORLD,ierr)

  END SUBROUTINE bcast_mpi_1d

  !--------------------------

  SUBROUTINE bcast_mpi_2d(nrank,im,jm,dat)

    USE MPI
    USE common_setting
    IMPLICIT NONE

    !---Common
    INTEGER MPI_R_SIZE,ierr

    !---IN
    INTEGER,INTENT(IN) :: nrank
    INTEGER,INTENT(IN) :: im,jm

    !---INOUT
    REAL(r_size),INTENT(INOUT) :: dat(im,jm)

    IF(r_size == kind(0.e0))THEN
       MPI_R_SIZE=MPI_REAL
    ELSE IF(r_size == kind(0.d0))THEN
       MPI_R_SIZE=MPI_DOUBLE_PRECISION
    END IF

    CALL MPI_BCAST(dat,im*jm,MPI_R_SIZE,nrank,MPI_COMM_WORLD,ierr)

  END SUBROUTINE bcast_mpi_2d

  !----------------------------

  SUBROUTINE bcast_mpi_3d(nrank,im,jm,km,dat)

    USE MPI
    USE common_setting
    IMPLICIT NONE

    !---Common
    INTEGER MPI_R_SIZE,ierr

    !---IN
    INTEGER,INTENT(IN) :: nrank
    INTEGER,INTENT(IN) :: im,jm,km

    !---INOUT
    REAL(r_size),INTENT(INOUT) :: dat(im,jm,km)

    IF(r_size == kind(0.e0))THEN
       MPI_R_SIZE=MPI_REAL
    ELSE IF(r_size == kind(0.d0))THEN
       MPI_R_SIZE=MPI_DOUBLE_PRECISION
    END IF

    CALL MPI_BCAST(dat,im*jm*km,MPI_R_SIZE,nrank,MPI_COMM_WORLD,ierr)

  END SUBROUTINE bcast_mpi_3d

  !-----------------------------------------------------------------------
  ! gridded data -> buffer
  !-----------------------------------------------------------------------
  SUBROUTINE grd_to_buf(grd,buf)

    !$USE OMP_LIB
    USE common_setting
    IMPLICIT NONE

    INTEGER :: i,j,m,ilon,ilat

    !IN    
    REAL(r_sngl),INTENT(IN) :: grd(nlon,nlat)

    !OUT
    REAL(r_sngl),INTENT(OUT) :: buf(nij1max,nprocs)

    !$OMP PARALLEL DO PRIVATE(m,i,j,ilon,ilat)    
    DO m=1,nprocs
       DO i=1,nij1node(m)
          j = m-1 + nprocs * (i-1)
          ilon = MOD(j,nlon) + 1
          ilat = (j-ilon+1) / nlon + 1
          buf(i,m) = grd(ilon,ilat)
       END DO
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE grd_to_buf

  !-----------------------------------------------------------------------
  ! buffer -> gridded data
  !-----------------------------------------------------------------------
  SUBROUTINE buf_to_grd(buf,grd)

    !$USE OMP_LIB
    USE common_setting
    IMPLICIT NONE

    INTEGER :: i,j,m,ilon,ilat

    !IN
    REAL(r_sngl),INTENT(IN) :: buf(nij1max,nprocs)

    !OUT
    REAL(r_sngl),INTENT(OUT) :: grd(nlon,nlat)

    !$OMP PARALLEL DO PRIVATE(m,i,j,ilon,ilat)    
    DO m=1,nprocs
       DO i=1,nij1node(m)
          j = m-1 + nprocs * (i-1)
          ilon = MOD(j,nlon) + 1
          ilat = (j-ilon+1) / nlon + 1
          grd(ilon,ilat) = buf(i,m)
       END DO
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE buf_to_grd

END MODULE common_mpi
