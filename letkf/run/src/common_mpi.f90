MODULE common_mpi
  !=======================================================================
  !
  ! [PURPOSE:] General MPI procedures
  !
  ! [HISTORY:]
  !   09/06/2005 Takemasa MIYOSHI  created
  !
  !=======================================================================
  IMPLICIT NONE
  PUBLIC

CONTAINS

  !----------------------------------------------------------------------
  
  SUBROUTINE initialize_mpi

    USE common_setting
    USE MPI
    IMPLICIT NONE

    INTEGER :: ierr
    
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

    IF(myrank == 0)        WRITE(6,'(A,I5.5,A,I5.5)') 'Hello from MYRANK ',myrank,'/',nprocs-1
    !IF(myrank == nprocs-1) WRITE(6,'(A,I5.5,A,I5.5)') 'Hello from MYRANK ',myrank,'/',nprocs-1

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

    USE MPI
    USE common_setting
    IMPLICIT NONE

    INTEGER :: j,k,n,ierr,n0

    REAL(r_sngl) :: bufs(nij1max,nlevall,nprocs)
    REAL(r_sngl) :: bufr(nij1max,nlevall)    

    !IN
    INTEGER,INTENT(IN) :: nrank
    REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
    REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)

    !OUT
    REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
    REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)


    !v3dg and v2dgs --> bufs
    IF(myrank == nrank)THEN
       j=0
       DO n=1,nv3d
          DO k=1,nlev
             j = j+1
             CALL grd_to_buf(v3dg(:,:,k,n),bufs(:,j,:))
          END DO
       END DO

       DO n=1,nv2d
          j = j+1
          CALL grd_to_buf(v2dg(:,:,n),bufs(:,j,:))
       END DO
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

    !bufs --> bufr
    n = nij1max * nlevall
    n0= n
    CALL MPI_SCATTER(bufs,n ,MPI_REAL,&
         & bufr,n0,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

    !bufr --> v3d and v2d
    j=0
    DO n=1,nv3d
       DO k=1,nlev
          j = j+1
          v3d(:,k,n) = REAL(bufr(1:nij1,j),r_size)
       END DO
    END DO

    DO n=1,nv2d
       j = j+1
       v2d(:,n) = REAL(bufr(1:nij1,j),r_size)
    END DO
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  END SUBROUTINE scatter_grd_mpi

  !-----------------------------------------------------------------------
  ! Gather gridded data (all -> nrank)
  !-----------------------------------------------------------------------
  SUBROUTINE gather_grd_mpi(nrank,v3d,v2d,v3dg,v2dg)

    USE MPI
    USE common_setting
    IMPLICIT NONE

    INTEGER j,k,n,ierr,n0
    
    REAL(r_sngl) bufs(nij1max,nlevall)
    REAL(r_sngl) bufr(nij1max,nlevall,nprocs)
    
    !---IN
    INTEGER,INTENT(IN) :: nrank
    REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
    REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)

    !---OUT
    REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
    REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)

    !bufs
    j=0
    DO n=1,nv3d
       DO k=1,nlev
          j = j+1
          bufs(1:nij1,j) = REAL(v3d(:,k,n),r_sngl)
       END DO
    END DO

    DO n=1,nv2d
       j = j+1
       bufs(1:nij1,j) = REAL(v2d(:,n),r_sngl)
    END DO

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

    !bufs --> bufr
    n = nij1max * nlevall
    n0= n
    CALL MPI_GATHER(bufs,n ,MPI_REAL,&
         & bufr,n0,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

    !bufr --> v3dg and v2dg
    IF(myrank == nrank)THEN
       
       j=0
       DO n=1,nv3d
          DO k=1,nlev
             j = j+1
             CALL buf_to_grd(bufr(:,j,:),v3dg(:,:,k,n))
          END DO
       END DO

       DO n=1,nv2d
          j = j+1
          CALL buf_to_grd(bufr(:,j,:),v2dg(:,:,n))
       END DO
       
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  END SUBROUTINE gather_grd_mpi

  !-----------------------------------------------------------------------
  ! Share nrank data with all processer
  !-----------------------------------------------------------------------

  SUBROUTINE bcast_mpi_1d(nrank,n,dat)

    USE MPI
    USE common_setting
    IMPLICIT NONE

    !---Common
    INTEGER ierr
    
    REAL(r_size) work(n)
    
    !---IN
    INTEGER,INTENT(IN) :: nrank
    INTEGER,INTENT(IN) :: n

    !---INOUT
    REAL(r_size),INTENT(INOUT) :: dat(n)

    work(:)=dat(:)
    
    IF(r_size == kind(0.e0))THEN
       CALL MPI_BCAST(dat,n,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
    ELSE IF(r_size == kind(0.d0))THEN
       CALL MPI_BCAST(dat,n,MPI_DOUBLE_PRECISION,nrank,MPI_COMM_WORLD,ierr)       
    END IF
    
  END SUBROUTINE bcast_mpi_1d
  
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
