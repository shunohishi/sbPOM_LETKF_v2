MODULE common_pom
  !=======================================================================
  !
  ! [PURPOSE:] Common Information for sbPOM
  !
  ! [HISTORY:]
  !   2004/10/15 Takemasa Miyoshi  created
  !   2009/01/23 Takemasa Miyoshi  modified
  !   2009/02/02 Takemasa Miyoshi  modified for ROMS
  !   2010/02/01 Yasumasa Miyazawa  modified for POM
  !   2024/04/01 Shun Ohishi        modified for sbPOM-LETKF v2.0
  !   2026/05/28 Shun Ohishi        modified for reading by root only
  !=======================================================================

CONTAINS

  !-----------------------------------------------------------------------
  ! Set the parameters
  !-----------------------------------------------------------------------
  SUBROUTINE set_common_pom

    !$ USE OMP_LIB
    USE common_io
    USE common_mpi
    USE common_setting
    USE MPI
    USE NETCDF
    IMPLICIT NONE

    INTEGER status,access
    INTEGER ncid
    INTEGER i,j,k
    INTEGER ierr

    REAL(r_dble),ALLOCATABLE :: tmp2d(:,:)
    REAL(r_dble),ALLOCATABLE :: sigma_level(:,:,:)
       
    status=0
    IF(myrank == root)THEN    
       status=access("grid.nc"," ")
    END IF

    CALL MPI_BCAST(status,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
    
    IF(status /= 0)THEN
       WRITE(6,'(A)') " *** Error: Not Found grid.nc"
       CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
    END IF

    IF(myrank == root)THEN

       ALLOCATE(tmp2d(nlon,nlat))
       
       !Open
       status = NF90_OPEN("grid.nc",NF90_NOWRITE,ncid)
       call handle_error_netcdf("NF90_OPEN: grid.nc",status)

       !2d
       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"east_e",tmp2d)
       lont(:,:)=REAL(tmp2d,r_size)

       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"east_u",tmp2d)
       lonu(:,:)=REAL(tmp2d,r_size)

       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"east_v",tmp2d)
       lonv(:,:)=REAL(tmp2d,r_size)
       
       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"north_e",tmp2d)
       latt(:,:)=REAL(tmp2d,r_size)

       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"north_u",tmp2d)
       latu(:,:)=REAL(tmp2d,r_size)

       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"north_v",tmp2d)
       latv(:,:)=REAL(tmp2d,r_size)
       
       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"h",tmp2d)
       tpg(:,:)=REAL(tmp2d,r_size)

       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"fsm",tmp2d)
       maskt(:,:)=REAL(tmp2d,r_size)

       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"dum",tmp2d)
       masku(:,:)=REAL(tmp2d,r_size)

       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"dvm",tmp2d)
       maskv(:,:)=REAL(tmp2d,r_size)
       
       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"dx",tmp2d)
       dlont(:,:)=REAL(tmp2d,r_size)
       dlonu(:,:)=REAL(tmp2d,r_size)
       dlonv(:,:)=REAL(tmp2d,r_size)

       call read_netcdf_var_dble_2d(ncid,nlon,nlat,"dy",tmp2d)
       dlatt(:,:)=REAL(tmp2d,r_size)
       dlatu(:,:)=REAL(tmp2d,r_size)
       dlatv(:,:)=REAL(tmp2d,r_size)

       DEALLOCATE(tmp2d)

       !3d
       ALLOCATE(sigma_level(nlon,nlat,nlev))
       
       !Sigma level (k=1: surface, k=nlev: dummy below the bottom) [0 ~ -1]
       call read_netcdf_var_dble_3d(ncid,nlon,nlat,nlev,"z_e",sigma_level)

       status = NF90_CLOSE(ncid)
       call handle_error_netcdf("NF90_CLOSE: grid.nc",status)

       !Depth [m] (negative)
       !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,k)
       DO k=1,nlev
          DO j=1,nlat
             DO i=1,nlon
                depth(i,j,k) = REAL(sigma_level(i,j,k)*tpg(i,j), r_size) !sigma_level
             END DO
          END DO
       END DO
       !$OMP END PARALLEL DO

       DEALLOCATE(sigma_level)
       
    END IF

    CALL bcast_mpi_2d(root,nlon,nlat,lont)
    CALL bcast_mpi_2d(root,nlon,nlat,lonu)
    CALL bcast_mpi_2d(root,nlon,nlat,lonv)

    CALL bcast_mpi_2d(root,nlon,nlat,latt)
    CALL bcast_mpi_2d(root,nlon,nlat,latu)
    CALL bcast_mpi_2d(root,nlon,nlat,latv)

    CALL bcast_mpi_2d(root,nlon,nlat,maskt)
    CALL bcast_mpi_2d(root,nlon,nlat,masku)
    CALL bcast_mpi_2d(root,nlon,nlat,maskv)

    CALL bcast_mpi_2d(root,nlon,nlat,dlont)
    CALL bcast_mpi_2d(root,nlon,nlat,dlonu)
    CALL bcast_mpi_2d(root,nlon,nlat,dlonv)

    CALL bcast_mpi_2d(root,nlon,nlat,dlatt)
    CALL bcast_mpi_2d(root,nlon,nlat,dlatu)
    CALL bcast_mpi_2d(root,nlon,nlat,dlatv)

    CALL bcast_mpi_2d(root,nlon,nlat,tpg)
    
    CALL bcast_mpi_3d(root,nlon,nlat,nlev,depth)
       
  END SUBROUTINE set_common_pom

  !--------------------------------------------------------------------------------
  
  SUBROUTINE set_common_mpi_pom

    !$ USE OMP_LIB
    USE common_mpi
    USE common_setting
    IMPLICIT NONE

    INTEGER i,j,k
    INTEGER iprocs

    REAL(r_size),ALLOCATABLE :: v3dg(:,:,:,:)
    REAL(r_size),ALLOCATABLE :: v2dg(:,:,:)
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

    ALLOCATE(ri1(nij1), rj1(nij1))
    
    ALLOCATE(lont1(nij1), latt1(nij1))
    ALLOCATE(lonu1(nij1), latu1(nij1))
    ALLOCATE(lonv1(nij1), latv1(nij1))

    ALLOCATE(dlont1(nij1), dlatt1(nij1))
    ALLOCATE(dlonu1(nij1), dlatu1(nij1))
    ALLOCATE(dlonv1(nij1), dlatv1(nij1))

    ALLOCATE(maskt1(nij1))
    ALLOCATE(masku1(nij1))
    ALLOCATE(maskv1(nij1))
    
    ALLOCATE(tpg1(nij1))
    ALLOCATE(depth1(nij1,nlev))

    ALLOCATE(v3dg(nlon,nlat,nlev,nv3d))
    ALLOCATE(v2dg(nlon,nlat,nv2d))
    ALLOCATE(v3d(nij1,nlev,nv3d))
    ALLOCATE(v2d(nij1,nv2d))

    !Initialization
    v3dg(:,:,:,:) = REAL(0.d0, r_size)
    v2dg(:,:,:)   = REAL(0.d0, r_size)
    v3d(:,:,:) = REAL(0.d0, r_size)
    v2d(:,:)   = REAL(0.d0, r_size)

    !Substitute grid data with v3dg and v2dg
    IF(myrank == root)THEN

       !$OMP PARALLEL PRIVATE(i,j,k)
       !$OMP DO SCHEDULE(STATIC)
       DO j=1,nlat
          DO i=1,nlon
             
             v3dg(i,j,1,1) = REAL(i, r_sngl)   !3D(1,1): i index
             v3dg(i,j,1,2) = REAL(j, r_sngl)   !3D(1,2): j index

             v3dg(i,j,2,1) = REAL(lont(i,j), r_sngl) !3D(2,1): lont
             v3dg(i,j,2,2) = REAL(latt(i,j), r_sngl) !3D(2,2): latt

             v3dg(i,j,3,1) = REAL(lonu(i,j), r_sngl) !3D(3,1): lonu
             v3dg(i,j,3,2) = REAL(latu(i,j), r_sngl) !3D(3,2): latu

             v3dg(i,j,4,1) = REAL(lonv(i,j), r_sngl) !3D(4,1): lonv
             v3dg(i,j,4,2) = REAL(latv(i,j), r_sngl) !3D(4,2): latv       

             v3dg(i,j,5,1) = REAL(dlont(i,j), r_sngl) !3D(5,1): dlont
             v3dg(i,j,5,2) = REAL(dlatt(i,j), r_sngl) !3D(5,2): dlatt

             v3dg(i,j,6,1) = REAL(dlonu(i,j), r_sngl) !3D(6,1): dlonu
             v3dg(i,j,6,2) = REAL(dlatu(i,j), r_sngl) !3D(6,2): dlatu

             v3dg(i,j,7,1) = REAL(dlonv(i,j), r_sngl) !3D(7,1): dlonv
             v3dg(i,j,7,2) = REAL(dlatv(i,j), r_sngl) !3D(7,2): dlatv

             v3dg(i,j,8,1) = REAL(maskt(i,j), r_sngl) !3D(8,1): maskt
             v3dg(i,j,8,2) = REAL(masku(i,j), r_sngl) !3D(8,2): masku
             v3dg(i,j,8,3) = REAL(maskv(i,j), r_sngl) !3D(8,3): maskv

             v2dg(i,j,1) = REAL(tpg(i,j), r_sngl) !2D(1)
             
          END DO
       END DO
       !$OMP END DO

       !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)   
       DO k=1,nlev
          DO j=1,nlat
             DO i=1,nlon
                v3dg(i,j,k,4) = REAL(depth(i,j,k), r_sngl) !3D(-,4) :depth
             END DO
          END DO
       END DO
       !$OMP END DO
       !$OMP END PARALLEL
       
    END IF
       
    !v3dg --> v3d, v2dg --> v2d
    CALL scatter_grd_mpi(root,v3dg,v2dg,v3d,v2d)

    DEALLOCATE(v3dg,v2dg)
    
    !Substitute v3d and v2d with grid data
    !$OMP PARALLEL PRIVATE(i,k)
    !$OMP DO SCHEDULE(STATIC)
    DO i=1,nij1
       ri1(i) = REAL(v3d(i,1,1), r_size)
       rj1(i) = REAL(v3d(i,1,2), r_size)

       lont1(i) = REAL(v3d(i,2,1), r_size)
       latt1(i) = REAL(v3d(i,2,2), r_size)
       
       lonu1(i) = REAL(v3d(i,3,1), r_size)
       latu1(i) = REAL(v3d(i,3,2), r_size)
       
       lonv1(i) = REAL(v3d(i,4,1), r_size)
       latv1(i) = REAL(v3d(i,4,2), r_size)
       
       dlont1(i) = REAL(v3d(i,5,1), r_size)
       dlatt1(i) = REAL(v3d(i,5,2), r_size)
       
       dlonu1(i) = REAL(v3d(i,6,1), r_size)
       dlatu1(i) = REAL(v3d(i,6,2), r_size)
       
       dlonv1(i) = REAL(v3d(i,7,1), r_size)
       dlatv1(i) = REAL(v3d(i,7,2), r_size)
       
       maskt1(i) = REAL(v3d(i,8,1), r_size)
       masku1(i) = REAL(v3d(i,8,2), r_size)
       maskv1(i) = REAL(v3d(i,8,3), r_size)

       tpg1(i) = REAL(v2d(i,1), r_size)
       
    END DO
    !$OMP END DO

    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
    DO k=1,nlev
       DO i=1,nij1
          depth1(i,k) = REAL(v3d(i,k,4), r_size)
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    
    DEALLOCATE(v3d,v2d)
    
  END SUBROUTINE set_common_mpi_pom
  
END MODULE common_pom
