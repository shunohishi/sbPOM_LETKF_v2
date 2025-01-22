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
  !=======================================================================

  IMPLICIT NONE
  PUBLIC

CONTAINS

  !-----------------------------------------------------------------------
  ! Set the parameters
  !-----------------------------------------------------------------------
  SUBROUTINE set_common_pom

    !$USE OMP_LIB
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

    REAL(r_dble) tmp2d(nlon,nlat)
    REAL(r_dble) sigma_level(nlon,nlat,nlev)

    WRITE(file_unit,'(A)') 'Hello from set_common_pom'

    ! Elements
    element(iv3d_u) = "U   "
    element(iv3d_v) = "V   "
    element(iv3d_t) = "T   "
    element(iv3d_s) = "SALT"
    element(nv3d+iv2d_z) = "ZETA"
    element(nv3d+iv2d_ubar) = "UBAR"
    element(nv3d+iv2d_vbar) = "VBAR"

    status=access("grid.nc"," ")
    IF(status == 0)THEN
       WRITE(file_unit,'(A)') "  >> accessing to file: grid.nc"
    ELSE
       WRITE(file_unit,'(A)') " *** Error: Not Found grid.nc"
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call finalize_mpi
       stop
    END IF

    status = NF90_OPEN("grid.nc",NF90_NOWRITE,ncid)
    call handle_error_netcdf("NF90_OPEN: grid.nc",status)

    call read_netcdf_var_dble(ncid,nlon,nlat,1,"east_e",tmp2d)
    lon(:,:)=REAL(tmp2d,r_size)
    call read_netcdf_var_dble(ncid,nlon,nlat,1,"north_e",tmp2d)
    lat(:,:)=REAL(tmp2d,r_size)
    call read_netcdf_var_dble(ncid,nlon,nlat,1,"h",tmp2d)
    phi0(:,:)=REAL(tmp2d,r_size)
    call read_netcdf_var_dble(ncid,nlon,nlat,1,"fsm",tmp2d)
    fsm(:,:)=REAL(tmp2d,r_size)
    call read_netcdf_var_dble(ncid,nlon,nlat,1,"dx",tmp2d)
    dlon(:,:)=REAL(tmp2d,r_size)
    call read_netcdf_var_dble(ncid,nlon,nlat,1,"dy",tmp2d)
    dlat(:,:)=REAL(tmp2d,r_size)
    !sbPOM original sigma level [1,nlev] --> [0, -1]
    call read_netcdf_var_dble(ncid,nlon,nlat,nlev,"z_e",sigma_level)
    
    status = NF90_CLOSE(ncid)
    call handle_error_netcdf("NF90_CLOSE: grid.nc",status)

    !sbPOM inversed sigma level [1, nlev] --> [-1,0]
    !$OMP PARALLEL    
    !$OMP DO PRIVATE(i,j,k)
    DO k=1,nlev
       DO j=1,nlat
          DO i=1,nlon
             depth(i,j,k)=sigma_level(i,j,nlev-k+1)*phi0(i,j) !sigma_level
          END DO
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

  END SUBROUTINE set_common_pom

  !--------------------------------------------------------------------------------
  
  SUBROUTINE set_common_mpi_pom

    !$USE OMP_LIB
    USE common_setting
    USE common_mpi
    IMPLICIT NONE

    INTEGER :: i,j,n

    REAL(r_sngl) :: v3dg(nlon,nlat,nlev,nv3d)
    REAL(r_sngl) :: v2dg(nlon,nlat,nv2d)
    REAL(r_size),ALLOCATABLE :: v3d(:,:,:)
    REAL(r_size),ALLOCATABLE :: v2d(:,:)

    WRITE(file_unit,'(A)') "Hello from set_common_mpi_pom"

    !nij1max,nij1
    i = MOD(nlon*nlat,nprocs)
    nij1max = (nlon*nlat - i)/nprocs + 1
    IF(myrank < i)THEN
       nij1 = nij1max
    ELSE
       nij1 = nij1max - 1
    END IF
    WRITE(file_unit,'(A,I5.5,A,I6)') "MYRANK: ",myrank," Number of grid points: nij1= ",nij1

    !nij1node
    ALLOCATE(nij1node(nprocs))
    DO n=1,nprocs
       IF(n-1 < i)THEN
          nij1node(n) = nij1max
       ELSE
          nij1node(n) = nij1max - 1
       END IF
    END DO

    ALLOCATE(phi1(nij1))
    ALLOCATE(lon1(nij1),lat1(nij1))
    ALLOCATE(dlon1(nij1), dlat1(nij1))
    ALLOCATE(ri1(nij1),rj1(nij1))
    ALLOCATE(depth1(nij1,nlev))

    ALLOCATE(v3d(nij1,nlev,nv3d))
    ALLOCATE(v2d(nij1,nv2d))

    !Initialization
    v3dg(:,:,:,:)=0.e0
    v2dg(:,:,:)=0.e0
    v3d(:,:,:)=0.d0
    v2d(:,:)=0.d0
    
    !lon,lat, i, j, dlon, dlat --> v3dg, pjo0 --> v2dg
    v3dg(:,:,1,1) = REAL(lon(:,:),r_sngl)
    v3dg(:,:,1,2) = REAL(lat(:,:),r_sngl)

    !$OMP PARALLEL DO PRIVATE(i,j)
    DO j=1,nlat
       DO i=1,nlon
          v3dg(i,j,2,1) = REAL(i,r_sngl)
          v3dg(i,j,2,2) = REAL(j,r_sngl)
       END DO
    END DO
    !$OMP END PARALLEL DO

    v3dg(:,:,3,1) = REAL(dlon(:,:),r_sngl)
    v3dg(:,:,3,2) = REAL(dlat(:,:),r_sngl)
    v3dg(:,:,:,3) = REAL(depth(:,:,:),r_sngl)
    v2dg(:,:,1) = REAL(phi0(:,:),r_sngl)

    !v3dg --> v3d, v2dg --> v2d
    CALL scatter_grd_mpi(0,v3dg,v2dg,v3d,v2d)

    !v3d --> lon,lat,ri1,ij1,dlon,dlat, v2d --> phi1
    lon1(:) = v3d(:,1,1)
    lat1(:) = v3d(:,1,2)
    ri1(:) = v3d(:,2,1)
    rj1(:) = v3d(:,2,2)
    dlon1(:) = v3d(:,3,1)
    dlat1(:) = v3d(:,3,2)
    depth1(:,:) = v3d(:,:,3)
    phi1(:) = v2d(:,1)
          
    DEALLOCATE(v3d,v2d)

  END SUBROUTINE set_common_mpi_pom
  
END MODULE common_pom
