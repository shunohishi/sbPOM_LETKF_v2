MODULE common_pom
  !=======================================================================
  !
  ! [PURPOSE:] Common Information for sbPOM
  !
  ! [HISTORY:]
  !   10/15/2004 Takemasa Miyoshi  created
  !   01/23/2009 Takemasa Miyoshi  modified
  !   02/02/2009 Takemasa Miyoshi  modified for ROMS
  !   02/01/2010 Yasumasa Miyazawa  modified for POM
  !   04/01/2024 Shun Ohishi        modified for sbPOM-LETKF v2.0
  !   05/28/2026 Shun Ohishi        modified for reading by root only
  !=======================================================================

CONTAINS

  !-----------------------------------------------------------------------
  ! Set the parameters
  !-----------------------------------------------------------------------
  SUBROUTINE set_common_pom

    !$ USE OMP_LIB
    USE common_setting
    USE common_mpi
    USE common_io
    USE NETCDF
    USE MPI
    IMPLICIT NONE

    INTEGER status,access
    INTEGER ncid
    INTEGER i,j,k
    INTEGER ierr

    REAL(r_dble),ALLOCATABLE :: tmp2d(:,:)
    REAL(r_dble),ALLOCATABLE :: sigma_level(:,:,:)
       
    ! Elements
    element(iv3d_u) = "U   "
    element(iv3d_v) = "V   "
    element(iv3d_t) = "T   "
    element(iv3d_s) = "SALT"
    element(nv3d+iv2d_z) = "ZETA"
    element(nv3d+iv2d_ubar) = "UBAR"
    element(nv3d+iv2d_vbar) = "VBAR"

    status=0
    IF(myrank == root)THEN    
       status=access("grid.nc"," ")
    END IF

    CALL MPI_BCAST(status,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
    
    IF(status /= 0)THEN
       WRITE(6,'(A)') " *** Error: Not Found grid.nc"
       call finalize_mpi
       stop
    END IF

    IF(myrank == root)THEN

       ALLOCATE(tmp2d(nlon,nlat))
       
       !Open
       status = NF90_OPEN("grid.nc",NF90_NOWRITE,ncid)
       call handle_error_netcdf("NF90_OPEN: grid.nc",status)

       !2d
       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"east_e",tmp2d)
       lon(:,:)=REAL(tmp2d,r_size)

       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"north_e",tmp2d)
       lat(:,:)=REAL(tmp2d,r_size)

       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"h",tmp2d)
       phi0(:,:)=REAL(tmp2d,r_size)

       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"fsm",tmp2d)
       fsm(:,:)=REAL(tmp2d,r_size)

       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"dx",tmp2d)
       dlon(:,:)=REAL(tmp2d,r_size)

       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"dy",tmp2d)
       dlat(:,:)=REAL(tmp2d,r_size)

       DEALLOCATE(tmp2d)

       !3d
       ALLOCATE(sigma_level(nlon,nlat,nlev))
       
       !sbPOM original sigma level [1,nlev] --> [0, -1]
       call read_netcdf_var_dble_3d(ncid,nlon,nlat,nlev,"z_e",sigma_level)

       status = NF90_CLOSE(ncid)
       call handle_error_netcdf("NF90_CLOSE: grid.nc",status)

       !sbPOM inversed sigma level [1, nlev] --> [-1,0]
       !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,k)
       DO k=1,nlev
          DO j=1,nlat
             DO i=1,nlon
                depth(i,j,k) = REAL(sigma_level(i,j,nlev-k+1)*phi0(i,j), r_size) !sigma_level
             END DO
          END DO
       END DO
       !$OMP END PARALLEL DO

       DEALLOCATE(sigma_level)
       
    END IF

    CALL bcast_mpi_2d(root,nlon,nlat,lon)
    CALL bcast_mpi_2d(root,nlon,nlat,lat)
    CALL bcast_mpi_2d(root,nlon,nlat,phi0)
    CALL bcast_mpi_2d(root,nlon,nlat,fsm)
    CALL bcast_mpi_2d(root,nlon,nlat,dlon)
    CALL bcast_mpi_2d(root,nlon,nlat,dlat)
    CALL bcast_mpi_3d(root,nlon,nlat,nlev,depth)
       
  END SUBROUTINE set_common_pom

  !--------------------------------------------------------------------------------
  
  SUBROUTINE set_common_mpi_pom

    !$ USE OMP_LIB
    USE common_setting
    USE common_mpi
    IMPLICIT NONE

    INTEGER i,j
    INTEGER iprocs

    REAL(r_sngl),ALLOCATABLE :: v3dg(:,:,:,:)
    REAL(r_sngl),ALLOCATABLE :: v2dg(:,:,:)
    REAL(r_size),ALLOCATABLE :: v3d(:,:,:)
    REAL(r_size),ALLOCATABLE :: v2d(:,:)
       
    !nij1max,nij1
    i = MOD(nlon*nlat,nprocs)
    nij1max = (nlon*nlat - i)/nprocs + 1
    IF(myrank < i)THEN
       nij1 = nij1max
    ELSE
       nij1 = nij1max - 1
    END IF

    !WRITE(6,'(A,I5.5,A,I6)') "MYRANK: ",myrank," Number of grid points: nij1= ",nij1

    !nij1node
    ALLOCATE(nij1node(nprocs))
    
    !$OMP PARALLEL DO PRIVATE(iprocs)
    DO iprocs=1,nprocs
       IF(iprocs-1 < i)THEN
          nij1node(iprocs) = nij1max
       ELSE
          nij1node(iprocs) = nij1max - 1
       END IF
    END DO
    !$OMP END PARALLEL DO

    
    ALLOCATE(phi1(nij1))
    ALLOCATE(lon1(nij1), lat1(nij1))
    ALLOCATE(dlon1(nij1), dlat1(nij1))
    ALLOCATE(ri1(nij1), rj1(nij1))
    ALLOCATE(depth1(nij1,nlev))

    ALLOCATE(v3dg(nlon,nlat,nlev,nv3d))
    ALLOCATE(v2dg(nlon,nlat,nv2d))
    ALLOCATE(v3d(nij1,nlev,nv3d))
    ALLOCATE(v2d(nij1,nv2d))

    !Initialization
    v3dg(:,:,:,:) = 0.e0
    v2dg(:,:,:) = 0.e0
    v3d(:,:,:) = REAL(0.d0, r_size)
    v2d(:,:) = REAL(0.d0, r_size)

    !lon,lat, i, j, dlon, dlat --> v3dg, pjo0 --> v2dg
    IF(myrank == root)THEN
       v3dg(:,:,1,1) = REAL(lon(:,:), r_sngl)
       v3dg(:,:,1,2) = REAL(lat(:,:), r_sngl)
       
       !$OMP PARALLEL DO PRIVATE(i,j)
       DO j=1,nlat
          DO i=1,nlon
             v3dg(i,j,2,1) = REAL(i, r_sngl)
             v3dg(i,j,2,2) = REAL(j, r_sngl)
          END DO
       END DO
       !$OMP END PARALLEL DO

       v3dg(:,:,3,1) = REAL(dlon(:,:), r_sngl)
       v3dg(:,:,3,2) = REAL(dlat(:,:), r_sngl)
       v3dg(:,:,:,3) = REAL(depth(:,:,:), r_sngl)
       v2dg(:,:,1) = REAL(phi0(:,:), r_sngl)
    END IF
       
    !v3dg --> v3d, v2dg --> v2d
    CALL scatter_grd_mpi(root,v3dg,v2dg,v3d,v2d)

    DEALLOCATE(v3dg,v2dg)
    
    !v3d --> lon,lat,ri1,ij1,dlon,dlat, v2d --> phi1
    lon1(:) = REAL(v3d(:,1,1), r_size)
    lat1(:) = REAL(v3d(:,1,2), r_size)
    ri1(:) = REAL(v3d(:,2,1), r_size)
    rj1(:) = REAL(v3d(:,2,2), r_size)
    dlon1(:) = REAL(v3d(:,3,1), r_size)
    dlat1(:) = REAL(v3d(:,3,2), r_size)
    depth1(:,:) = REAL(v3d(:,:,3), r_size)
    phi1(:) = REAL(v2d(:,1), r_size)
    
    DEALLOCATE(v3d,v2d)

  END SUBROUTINE set_common_mpi_pom
  
END MODULE common_pom
