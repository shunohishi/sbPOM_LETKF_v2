MODULE common_io

  IMPLICIT NONE
  PUBLIC

CONTAINS

  !-----------------------------------------------------------------------
  ! Common NetCDF
  !-----------------------------------------------------------------------

  SUBROUTINE handle_error_netcdf(routine,status)

    USE common_mpi
    USE MPI
    USE NETCDF
    implicit none

    character(*),intent(in) :: routine
    integer,intent(in) :: status

    if(status /= nf90_noerr) then
       write(6,'(A)') "Error: NetCDF "//trim(routine)
       write(6,'(A)') trim(nf90_strerror(status))
       call finalize_mpi
       stop
    end if

  END SUBROUTINE handle_error_netcdf

  !---------------------------------

  SUBROUTINE read_netcdf_var_int(ncid,im,jm,km,varname,glb)

    USE common_setting, only:r_sngl
    USE NETCDF
    implicit none

    integer status,varid

    integer,intent(in) :: ncid
    integer,intent(in) :: im,jm,km

    character(*),intent(in) :: varname

    integer,intent(out) :: glb(im,jm,km)

    status=nf90_inq_varid(ncid,varname,varid)
    call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status)
    status=nf90_get_var(ncid,varid,glb,(/1,1,1/),(/im,jm,km/))
    call handle_error_netcdf('nf90_get_var:'//trim(varname),status)

  END SUBROUTINE read_netcdf_var_int
  
  !---------------------------------
  
  SUBROUTINE read_netcdf_var_sngl(ncid,im,jm,km,varname,glb)

    USE common_setting, only:r_sngl
    USE NETCDF
    implicit none

    integer status,varid

    integer,intent(in) :: ncid
    integer,intent(in) :: im,jm,km

    character(*),intent(in) :: varname

    real(kind = r_sngl),intent(out) :: glb(im,jm,km)

    status=nf90_inq_varid(ncid,varname,varid)
    call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status)
    status=nf90_get_var(ncid,varid,glb,(/1,1,1/),(/im,jm,km/))
    call handle_error_netcdf('nf90_get_var:'//trim(varname),status)

  END SUBROUTINE read_netcdf_var_sngl

  !-----------------------------------

  SUBROUTINE read_netcdf_var_dble(ncid,im,jm,km,varname,glb)

    USE common_setting, only:r_dble    
    USE NETCDF
    implicit none

    integer status,varid

    integer,intent(in) :: ncid
    integer,intent(in) :: im,jm,km

    character(*),intent(in) :: varname

    real(kind = r_dble),intent(out) :: glb(im,jm,km)

    status=nf90_inq_varid(ncid,varname,varid)
    call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status)
    status=nf90_get_var(ncid,varid,glb,(/1,1,1/),(/im,jm,km/))
    call handle_error_netcdf('nf90_get_var:'//trim(varname),status)

  END SUBROUTINE read_netcdf_var_dble

  !----------------------------------

  SUBROUTINE write_netcdf_var_int(ncid,im,jm,km,varname,glb)

    USE common_setting, only:r_sngl
    USE NETCDF
    implicit none

    integer status,varid

    integer,intent(in) :: ncid
    integer,intent(in) :: im,jm,km

    character(*),intent(in) :: varname

    integer,intent(in) :: glb(im,jm,km)

    status=nf90_inq_varid(ncid,varname,varid)
    call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status)
    status=nf90_put_var(ncid,varid,glb,(/1,1,1/),(/im,jm,km/))
    call handle_error_netcdf('nf90_put_var:'//trim(varname),status)

  END SUBROUTINE write_netcdf_var_int
  
  !----------------------------------
  
  SUBROUTINE write_netcdf_var_sngl(ncid,im,jm,km,varname,glb)

    USE common_setting, only:r_sngl
    USE NETCDF
    implicit none

    integer status,varid

    integer,intent(in) :: ncid
    integer,intent(in) :: im,jm,km

    character(*),intent(in) :: varname

    real(kind = r_sngl),intent(in) :: glb(im,jm,km)

    status=nf90_inq_varid(ncid,varname,varid)
    call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status)
    status=nf90_put_var(ncid,varid,glb,(/1,1,1/),(/im,jm,km/))
    call handle_error_netcdf('nf90_put_var:'//trim(varname),status)

  END SUBROUTINE write_netcdf_var_sngl
    
  !-------------------------------

  SUBROUTINE write_pnetcdf_var_sngl(ncid,im,jm,k,varname,glb)

    USE common_setting, only:r_sngl
    USE NETCDF
    implicit none

    integer status,varid

    integer,intent(in) :: ncid
    integer,intent(in) :: im,jm,k

    character(*),intent(in) :: varname

    real(kind = r_sngl),intent(in) :: glb(im,jm)

    status=nf90_inq_varid(ncid,varname,varid)
    call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status)
    status=nf90_var_par_access(ncid,varid,nf90_collective)
    call handle_error_netcdf('nf90_var_par_access:'//trim(varname),status)
    status=nf90_put_var(ncid,varid,glb,start=(/1,1,k/),count=(/im,jm,1/))
    call handle_error_netcdf('nf90_put_var:'//trim(varname),status)

  END SUBROUTINE write_pnetcdf_var_sngl

  !-------------------------------------------------------------------------
  ! Define variable |
  !-------------------------------------------------------------------------
  
  SUBROUTINE define_var_netcdf(ncid,ndim,dim,varid,type,name,long_name,units_name)

    use netcdf
    implicit none

    integer status

    integer,intent(in) :: ncid
    integer,intent(in) :: ndim,dim(ndim)
    integer,intent(inout) :: varid

    character(*),intent(in) :: type
    character(*),intent(in) :: name,long_name,units_name

    if(trim(type) == "int")then
       status=NF90_DEF_VAR(ncid,trim(name),nf90_int,dim,varid)
    else if(trim(type) == "real")then
       status=NF90_DEF_VAR(ncid,trim(name),nf90_float,dim,varid)
    else if(trim(type) == "dble")then
       status=NF90_DEF_VAR(ncid,trim(name),nf90_double,dim,varid)
    end if
    call handle_error_netcdf("NF90_DEF_VAR: "//trim(name),status)
    
    status=NF90_DEF_VAR_DEFLATE(ncid,varid,shuffle=1,deflate=1,deflate_level=5)
    call handle_error_netcdf("NF90_DEF_VAR_DEFLATE",status)

    status=NF90_PUT_ATT(ncid,varid,"long_name",trim(long_name))
    call handle_error_netcdf("NF90_PUT_ATT: "//trim(long_name),status)

    status=NF90_PUT_ATT(ncid,varid,"units",trim(units_name))
    call handle_error_netcdf("NF90_PUT_ATT: "//trim(units_name),status)

  END SUBROUTINE define_var_netcdf
  
  !--------------------------------------------------------------------
  ! Variable name |
  !--------------------------------------------------------------------

  SUBROUTINE detect_varname(ncid)
    
    USE common_setting
    USE common_mpi
    USE NETCDF
    USE MPI
    implicit none

    !Common
    integer status1,status2
    integer varid,ierr

    !IN
    integer,intent(in) :: ncid

    !File type
    status1 = NF90_INQ_VARID(ncid,"elb",varid)
    status2 = NF90_INQ_VARID(ncid,"el_fcst",varid)
    if(status1 == 0)then !Restart/Intermittent
       file_type=1       
    else if(status2 == 0)then !IAU
       file_type=2
    else
       write(6,'(A)') trim(nf90_strerror(status1))
       write(6,'(A)') trim(nf90_strerror(status2))
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call finalize_mpi
       stop
    end if

    !Variable name
    if(file_type == 1)then

       var2df(iv2d_z)="elb"
       var2df(iv2d_ubar)="uab"
       var2df(iv2d_vbar)="vab"
       var3df(iv3d_u)="ub"
       var3df(iv3d_v)="vb"
       var3df(iv3d_t)="tb"
       var3df(iv3d_s)="sb"

       var2da(:)=var2df(:)
       var3da(:)=var3df(:)
       
    else if(file_type == 2)then

       var2df(iv2d_z)="el_fcst"
       var2df(iv2d_ubar)="ua_fcst"
       var2df(iv2d_vbar)="va_fcst"
       var3df(iv3d_u)="u_fcst"
       var3df(iv3d_v)="v_fcst"
       var3df(iv3d_t)="t_fcst"
       var3df(iv3d_s)="s_fcst"

       var2da(iv2d_z)="el_anal"
       var2da(iv2d_ubar)="ua_anal"
       var2da(iv2d_vbar)="va_anal"
       var3da(iv3d_u)="u_anal"
       var3da(iv3d_v)="v_anal"
       var3da(iv3d_t)="t_anal"
       var3da(iv3d_s)="s_anal"
       
    end if    
    
  END SUBROUTINE detect_varname
  
  !--------------------------------------------------------------------
  ! File I/O
  !--------------------------------------------------------------------

  SUBROUTINE read_state_vector(filename,v3d,v2d)

    !$USE OMP_LIB
    USE common_setting
    USE common_mpi
    USE NETCDF
    USE MPI
    IMPLICIT NONE

    !Common
    INTEGER i,j,k
    INTEGER status,access
    INTEGER ncid
    INTEGER ierr

    REAL(r_sngl) tmp2d(nlon,nlat),tmp3d(nlon,nlat,nlev,nv3d)
    
    !IN
    CHARACTER(12),INTENT(IN) :: filename

    !OUT
    REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
    REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)

    status=access(filename," ")
    IF(status == 0)THEN
       WRITE(file_unit,'(A,I5.5,A)') "MYRANK",myrank," is reading "//trim(filename)
    ELSE
       WRITE(file_unit,'(A)') " *** Error: Not Found "//trim(filename)
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call finalize_mpi
       stop
    END IF

    !Open
    status = NF90_OPEN(trim(filename),NF90_NOWRITE,ncid)
    call handle_error_netcdf("NF90_OPEN: "//trim(filename),status)

    call detect_varname(ncid)
        
    !2D
    call read_netcdf_var_sngl(ncid,nlon,nlat,1,trim(var2df(iv2d_z)),tmp2d)
    v2d(:,:,iv2d_z)=REAL(tmp2d(:,:),r_size)

    call read_netcdf_var_sngl(ncid,nlon,nlat,1,trim(var2df(iv2d_ubar)),tmp2d)
    v2d(:,:,iv2d_ubar)=REAL(tmp2d(:,:),r_size)

    call read_netcdf_var_sngl(ncid,nlon,nlat,1,trim(var2df(iv2d_vbar)),tmp2d)
    v2d(:,:,iv2d_vbar)=REAL(tmp2d(:,:),r_size)

    !3D
    call read_netcdf_var_sngl(ncid,nlon,nlat,nlev,trim(var3df(iv3d_u)),tmp3d(:,:,:,iv3d_u))
    call read_netcdf_var_sngl(ncid,nlon,nlat,nlev,trim(var3df(iv3d_v)),tmp3d(:,:,:,iv3d_v))
    call read_netcdf_var_sngl(ncid,nlon,nlat,nlev,trim(var3df(iv3d_t)),tmp3d(:,:,:,iv3d_t))
    call read_netcdf_var_sngl(ncid,nlon,nlat,nlev,trim(var3df(iv3d_s)),tmp3d(:,:,:,iv3d_s))

    ! POM -> ROMS type vertical order for TransXtoY    
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    DO k=1,nlev
       DO j=1,nlat
          DO i=1,nlon
             v3d(i,j,nlev-k+1,iv3d_u) = REAL(tmp3d(i,j,k,iv3d_u),r_size)
             v3d(i,j,nlev-k+1,iv3d_v) = REAL(tmp3d(i,j,k,iv3d_v),r_size)
             v3d(i,j,nlev-k+1,iv3d_t) = REAL(tmp3d(i,j,k,iv3d_t),r_size)
             v3d(i,j,nlev-k+1,iv3d_s) = REAL(tmp3d(i,j,k,iv3d_s),r_size)
          END DO
       END DO
    END DO
    !$OMP END PARALLEL DO

    status = NF90_CLOSE(ncid)
    call handle_error_netcdf("NF90_CLOSE: "//trim(filename),status)

  END SUBROUTINE read_state_vector

  !-----------------------------

  SUBROUTINE write_state_vector(nrank,filename,v3d,v2d)

    !$USE OMP_LIB
    USE common_setting
    USE common_mpi
    USE NETCDF
    USE MPI
    IMPLICIT NONE

    INTEGER i,j,k
    INTEGER status,access
    INTEGER ncid
    INTEGER ierr

    REAL(r_sngl) tmp2d(nlon,nlat,nv2d)    
    REAL(r_sngl) tmp3d(nlon,nlat,nlev,nv3d)

    !--IN
    INTEGER,INTENT(IN) :: nrank
    REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
    REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)

    CHARACTER(12),INTENT(IN) :: filename
    
    !---Post process
    tmp2d(:,:,:)=v2d(:,:,:)
    !Restore vertical order
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    DO k=1,nlev
       DO j=1,nlat
          DO i=1,nlon
             tmp3d(i,j,nlev-k+1,iv3d_u)=v3d(i,j,k,iv3d_u)
             tmp3d(i,j,nlev-k+1,iv3d_v)=v3d(i,j,k,iv3d_v)
             tmp3d(i,j,nlev-k+1,iv3d_t)=v3d(i,j,k,iv3d_t)
             tmp3d(i,j,nlev-k+1,iv3d_s)=v3d(i,j,k,iv3d_s)
          END DO
       END DO
    END DO
    !$OMP END PARALLEL DO

    !Quasi-Global
    IF(lqglobal)THEN
       tmp2d(nlon-1:nlon,:,:)=tmp2d(1:2,:,:)
       tmp3d(nlon-1:nlon,:,:,:)=tmp3d(1:2,:,:,:)
    END IF
    
    !---Access
    status=access(filename," ")
    IF(status /= 0)THEN
       WRITE(file_unit,'(A)') " *** Error: Not Found "//trim(filename)
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call finalize_mpi
       stop
    ELSE IF(myrank == nrank .or. pnetcdf)THEN
       WRITE(file_unit,'(A,I5.5,A)') "MYRANK",myrank," is writing "//trim(filename)
    END IF

    !---Open
    IF(pnetcdf)THEN
       status = NF90_OPEN(filename,NF90_WRITE,ncid, &
            & comm=MPI_COMM_WORLD,info=MPI_INFO_NULL)
    ELSE IF(myrank == nrank)THEN
       status = NF90_OPEN(filename,NF90_WRITE,ncid)
    END IF
    call handle_error_netcdf("NF90_OPEN: "//trim(filename),status)

    call detect_varname(ncid)
    
    !---Write
    IF(pnetcdf)THEN

       !***** To be constructed *****

    ELSE IF(myrank == nrank)THEN

       call write_netcdf_var_sngl(ncid,nlon,nlat,1,trim(var2da(iv2d_z)),tmp2d(:,:,iv2d_z))
       call write_netcdf_var_sngl(ncid,nlon,nlat,1,trim(var2da(iv2d_ubar)),tmp2d(:,:,iv2d_ubar))
       call write_netcdf_var_sngl(ncid,nlon,nlat,1,trim(var2da(iv2d_vbar)),tmp2d(:,:,iv2d_vbar))

       call write_netcdf_var_sngl(ncid,nlon,nlat,nlev,trim(var3da(iv3d_u)),tmp3d(:,:,:,iv3d_u))
       call write_netcdf_var_sngl(ncid,nlon,nlat,nlev,trim(var3da(iv3d_v)),tmp3d(:,:,:,iv3d_v))
       call write_netcdf_var_sngl(ncid,nlon,nlat,nlev,trim(var3da(iv3d_t)),tmp3d(:,:,:,iv3d_t))
       call write_netcdf_var_sngl(ncid,nlon,nlat,nlev,trim(var3da(iv3d_s)),tmp3d(:,:,:,iv3d_s))

    END IF

    !---CLOSE
    IF(pnetcdf .or. myrank == nrank)THEN
       status = NF90_CLOSE(ncid)
       call handle_error_netcdf("NF90_CLOSE: "//trim(filename),status)
    END IF

  END SUBROUTINE write_state_vector

  !-----------------------------------------------------------------------
  ! Read ensemble data and distribute to processes
  !-----------------------------------------------------------------------

  SUBROUTINE read_ens_mpi(file,v3d,v2d)

    USE common_setting
    USE common_mpi
    IMPLICIT NONE

    !---Common
    INTEGER i,n
    INTEGER ibv,iprocs

    REAL(r_size) v3dg(nlon,nlat,nlev,nv3d)
    REAL(r_size) v2dg(nlon,nlat,nv2d)

    CHARACTER(12) :: filename="file00000.nc"

    !---IN
    CHARACTER(4),INTENT(IN) :: file

    !---OUT
    REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nbv,nv3d)
    REAL(r_size),INTENT(OUT) :: v2d(nij1,nbv,nv2d)

    n = CEILING(REAL(nbv)/REAL(nprocs))

    DO i=1,n

       ibv = myrank+1 + (i-1)*nprocs

       IF(ibv <= nbv)THEN
          WRITE(filename(1:9),'(A4,I5.5)') file,ibv
          WRITE(file_unit,'(A,I5.5,A)') "MYRANK ",myrank," is reading a file "//filename
          CALL read_state_vector(filename,v3dg,v2dg)
       END IF

       DO iprocs=0,nprocs-1
          ibv = iprocs+1 + (i-1)*nprocs
          IF(ibv <= nbv) THEN
             CALL scatter_grd_mpi(iprocs,REAL(v3dg,r_sngl),REAL(v2dg,r_sngl), &
                  & v3d(:,:,ibv,:),v2d(:,ibv,:))
          END IF
       END DO !iprocs

    END DO !i

  END SUBROUTINE read_ens_mpi

  !-----------------------------------------------------------------------
  ! Write ensemble data after collecting data from processes
  !-----------------------------------------------------------------------
  SUBROUTINE write_ens_mpi(v3df,v2df,v3da,v2da)

    USE MPI
    USE common_setting
    USE common
    USE common_mpi
    IMPLICIT NONE

    !---Common
    INTEGER,PARAMETER :: nbv_all=nbv+4 ![1,nbv]: ensemble analysis
                                       ![nbv+1]: fcst_mean, [nbv+2]: fcst_sprd
                                       ![nbv+3]: anal_mean, [nbv+4]: anal_sprd
    INTEGER i,n
    INTEGER iprocs,ibv
    INTEGER iprocs1,iprocs2,iprocs3,iprocs4
    INTEGER ierr

    REAL(r_size) :: work2d(nobs,nbv)
    
    REAL(r_size) :: v3df_m(nij1,nlev,nv3d),v2df_m(nij1,nv2d)       !Forecast ensemble mean
    REAL(r_size) :: v3df_s(nij1,nlev,nv3d),v2df_s(nij1,nv2d)       !                  spread
    REAL(r_size) :: v3da_m(nij1,nlev,nv3d),v2da_m(nij1,nv2d)       !Analysis ensemble mean
    REAL(r_size) :: v3da_s(nij1,nlev,nv3d),v2da_s(nij1,nv2d)       !                  spread

    REAL(r_sngl) :: v3dg(nlon,nlat,nlev,nv3d),v2dg(nlon,nlat,nv2d) !Global domain

    REAL(r_size) :: hxf(nobs,nbv),hxa(nobs,nbv)                  !Forecast/Analysis ensemble in obs. space
    REAL(r_size) :: hxfmean(nobs),hxamean(nobs)                  !Forecast/Analysis ensembale mean in obs. space
    REAL(r_size) :: hxfsprd(nobs),hxasprd(nobs)                  !Forecast/Analysis ensembale spread in obs. space

    CHARACTER(12) :: filename='file00000.nc'

    !---IN    
    REAL(r_size),INTENT(IN) :: v3df(nij1,nlev,nbv,nv3d),v2df(nij1,nbv,nv2d) !Forecast (local)
    REAL(r_size),INTENT(IN) :: v3da(nij1,nlev,nbv,nv3d),v2da(nij1,nbv,nv2d) !Analysis (local)

    !---Initialization
    hxf(:,:)=0.d0
    hxa(:,:)=0.d0
    
    !---Ensemble mean/sprd
    call ensemble_mesp(v3df,v2df,v3df_m,v3df_s,v2df_m,v2df_s)    
    call ensemble_mesp(v3da,v2da,v3da_m,v3da_s,v2da_m,v2da_s)

    !---Repeat index
    n = CEILING(REAL(nbv_all)/REAL(nprocs))

    !---Forecast ensemble (No output) --> H(Xf)
    DO i=1,n

       !GATHER to iprocs
       DO iprocs=0,nprocs-1

          ibv = iprocs+1 + (i-1)*nprocs

          IF(ibv <= nbv)THEN
             !Forecast ensemble
             CALL gather_grd_mpi(iprocs,v3df(:,:,ibv,:),v2df(:,ibv,:),v3dg,v2dg)
          END IF
          
       END DO !iprocs

       ibv = myrank+1 + (i-1)*nprocs

       !Observation space
       IF(ibv <= nbv)THEN
          CALL Global_Trans_XtoY(v3dg,v2dg,hxf(:,ibv))
       END IF
       
    END DO !i

    !---Analysis ensemble, forecast/analysis ensemble mean, and forecast/analysis ensemble spread (Output) --> Obs. space
    DO i=1,n

       !GATHER to iprocs
       DO iprocs=0,nprocs-1

          ibv = iprocs+1 + (i-1)*nprocs

          IF(ibv <= nbv)THEN
             !Analysis ensemble
             CALL gather_grd_mpi(iprocs,v3da(:,:,ibv,:),v2da(:,ibv,:),v3dg,v2dg)
          ELSE IF(ibv == nbv+1)THEN
             !Forecast ensemble mean
             iprocs1=iprocs
             CALL gather_grd_mpi(iprocs,v3df_m,v2df_m,v3dg,v2dg)
          ELSE IF(ibv == nbv+2)THEN
             !Forecast Ensemble spread
             iprocs2=iprocs
             CALL gather_grd_mpi(iprocs,v3df_s,v2df_s,v3dg,v2dg)
          ELSE IF(ibv == nbv+3)THEN
             !Analysis ensemble mean
             iprocs3=iprocs
             CALL gather_grd_mpi(iprocs,v3da_m,v2da_m,v3dg,v2dg)
          ELSE IF(ibv == nbv+4)THEN
             !Analysis ensemble spread
             iprocs4=iprocs
             CALL gather_grd_mpi(iprocs,v3da_s,v2da_s,v3dg,v2dg)
          END IF
          
       END DO !iprocs
       
       ibv = myrank+1 + (i-1)*nprocs

       !Observation space
       IF(ibv <= nbv)THEN
          CALL Global_Trans_XtoY(v3dg,v2dg,hxa(:,ibv))
       ELSE IF(ibv == nbv+1)THEN
          CALL Global_Trans_XtoY(v3dg,v2dg,hxfmean)
       ELSE IF(ibv == nbv+2)THEN
          CALL Global_Trans_XtoY(v3dg,v2dg,hxfsprd)
       ELSE IF(ibv == nbv+3)THEN
          CALL Global_Trans_XtoY(v3dg,v2dg,hxamean)
       ELSE IF(ibv == nbv+4)THEN
          CALL Global_Trans_XtoY(v3dg,v2dg,hxasprd)
       END IF
       
       !Filename
       IF(ibv <= nbv)THEN
          WRITE(filename(1:9),'(A4,I5.5)') "fa01",ibv
       ELSE IF(ibv == nbv+1)THEN
          WRITE(filename(1:9),'(A9)') "fcst_mean"
       ELSE IF(ibv == nbv+2)THEN
          WRITE(filename(1:9),'(A9)') "fcst_sprd"
       ELSE IF(ibv == nbv+3)THEN
          WRITE(filename(1:9),'(A9)') "anal_mean"
       ELSE IF(ibv == nbv+4)THEN
          WRITE(filename(1:9),'(A9)') "anal_sprd"          
       END IF
       
       IF(ibv <= nbv_all)THEN
          WRITE(file_unit,'(A,I5.5,A)') "MYRANK ",myrank," is writing a file "//filename
          CALL write_state_vector(ibv-1,filename,v3dg,v2dg)
       END IF
       
    END DO !i
    
    !---Share data in observation space
    !H(Xf)
    work2d(:,:)=hxf(:,:)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    IF(r_size == kind(0.d0))THEN
       CALL MPI_ALLREDUCE(work2d,hxf,nobs*nbv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    ELSE IF(r_size == kind(0.e0))THEN
       CALL MPI_ALLREDUCE(work2d,hxf,nobs*nbv,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
    END IF
    
    !H(Xa)
    work2d(:,:)=hxa(:,:)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    IF(r_size == kind(0.d0))THEN
       CALL MPI_ALLREDUCE(work2d,hxa,nobs*nbv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    ELSE IF(r_size == kind(0.e0))THEN
       CALL MPI_ALLREDUCE(work2d,hxa,nobs*nbv,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
    END IF

    CALL bcast_mpi_1d(iprocs1,nobs,hxfmean)
    CALL bcast_mpi_1d(iprocs2,nobs,hxfsprd)
    CALL bcast_mpi_1d(iprocs3,nobs,hxamean)
    CALL bcast_mpi_1d(iprocs4,nobs,hxasprd)   
        
    !---Monitor
    IF(myrank == 0)THEN
       CALL monit_mean("fcst",hxfmean)
       CALL monit_mean("anal",hxamean)
       CALL monit_sprd("fcst",hxfsprd)
       CALL monit_sprd("anal",hxasprd)
    END IF
    
    !Make & Write innovation netcdf file
    IF(myrank == nprocs-1)THEN
       CALL make_ncfile_innovation
       CALL write_innovation(hxf,hxa,hxfmean,hxfsprd,hxamean,hxasprd)
    END IF       
    
  END SUBROUTINE write_ens_mpi

  !-----------------------------------------------------------------------
  ! Get number of observations
  !-----------------------------------------------------------------------

  SUBROUTINE get_nobs(filename,no)

    USE common_setting
    USE common_mpi
    USE NETCDF
    USE MPI
    IMPLICIT NONE

    INTEGER status,access
    INTEGER ncid,dimid
    INTEGER io    
    INTEGER nu,nv,nt,ns,nz
    INTEGER ierr

    INTEGER,ALLOCATABLE :: ele(:)

    !IN
    CHARACTER(8),INTENT(IN) :: filename

    !OUT
    INTEGER,INTENT(OUT) :: no

    status=access(filename," ")
    IF(status == 0)THEN
       WRITE(file_unit,'(A)') "  >> accessing to file: "//filename
    ELSE
       WRITE(file_unit,'(A)') " *** Error: Not Found "//filename
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call finalize_mpi
       stop
    END IF

    status=NF90_OPEN(trim(filename),NF90_NOWRITE,ncid)
    call handle_error_netcdf("NF90_OPEN: "//trim(filename),status)

    status=NF90_INQ_DIMID(ncid,"nobs",dimid)
    call handle_error_netcdf("NF90_INQ_DIMID: nobs",status)

    status=NF90_INQUIRE_DIMENSION(ncid,dimid,len = no)    
    call handle_error_netcdf("NF90_INQUIRE_DIMENSION: nobs",status)

    ALLOCATE(ele(no))

    call read_netcdf_var_int(ncid,no,1,1,"ele",ele)
    
    status=NF90_CLOSE(ncid)
    call handle_error_netcdf("NF90_CLOSE: "//trim(filename),status)

    nu=0
    nv=0
    nt=0
    ns=0
    nz=0

    DO io=1,no

       SELECT CASE(ele(io))
       CASE(id_u_obs)
          nu = nu + 1
       CASE(id_v_obs)
          nv = nv + 1
       CASE(id_t_obs)
          nt = nt + 1
       CASE(id_s_obs)
          ns = ns + 1
       CASE(id_z_obs)
          nz = nz + 1
       END SELECT

    END DO

    DEALLOCATE(ele)

    WRITE(file_unit,'(A)') ""    
    WRITE(file_unit,'(A)') "== NUMBER OF OBSERVATION ======"
    WRITE(file_unit,'(A15,I10)') "NUMBER of OBS: ",no
    WRITE(file_unit,'(A15,I10)') "            U: ",nu
    WRITE(file_unit,'(A15,I10)') "            V: ",nv
    WRITE(file_unit,'(A15,I10)') "            T: ",nt
    WRITE(file_unit,'(A15,I10)') "         SALT: ",ns
    WRITE(file_unit,'(A15,I10)') "         ZETA: ",nz
    WRITE(file_unit,'(A)') "==============================="
    WRITE(file_unit,'(A)') ""    

  END SUBROUTINE get_nobs

  !--------------------------------------------------------------------
  ! Read observation
  !--------------------------------------------------------------------

  SUBROUTINE read_obs(filename,no,ele,ins,lon,lat,lev,obs,err)

    !$USE OMP_LIB
    USE common_setting, only: r_sngl,r_size, file_unit
    USE common_mpi
    USE NETCDF
    USE MPI
    IMPLICIT NONE

    !Common
    INTEGER status,access
    INTEGER ncid
    INTEGER ierr

    INTEGER :: itmp(no)
    REAL(r_sngl) :: tmp(no)

    !IN
    INTEGER,INTENT(IN) :: no
    CHARACTER(8),INTENT(IN) :: filename

    !OUT
    INTEGER,INTENT(OUT) :: ele(no) ! element ID
    INTEGER,INTENT(OUT) :: ins(no) ! instrument ID

    REAL(r_size),INTENT(OUT) :: lon(no),lat(no),lev(no) ! Grid information [degree E, degree N, m]
    REAL(r_size),INTENT(OUT) :: obs(no) ! Obs.
    REAL(r_size),INTENT(OUT) :: err(no) ! Obs. error

    status=access(filename," ")
    IF(status == 0)THEN
       WRITE(file_unit,'(A)') "  >> accessing to file: "//filename
    ELSE
       WRITE(file_unit,'(A)') " *** Error: Not Found "//filename
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call finalize_mpi
       stop
    END IF

    !open
    status=NF90_OPEN(trim(filename),NF90_NOWRITE,ncid)
    call handle_error_netcdf("NF90_OPEN: "//trim(filename),status)

    !element
    call read_netcdf_var_int(ncid,no,1,1,"ele",itmp)    
    ele(:)=itmp(:)

    !instrument
    call read_netcdf_var_int(ncid,no,1,1,"ins",itmp)    
    ins(:)=itmp(:)
    
    !longitude
    call read_netcdf_var_sngl(ncid,no,1,1,"lon",tmp)    
    lon(:)=REAL(tmp(:),r_size)

    !latitude
    call read_netcdf_var_sngl(ncid,no,1,1,"lat",tmp)    
    lat(:)=REAL(tmp(:),r_size)

    !level
    call read_netcdf_var_sngl(ncid,no,1,1,"lev",tmp)    
    lev(:)=REAL(tmp(:),r_size)

    !data
    call read_netcdf_var_sngl(ncid,no,1,1,"dat",tmp)    
    obs(:)=REAL(tmp(:),r_size)

    !Error
    call read_netcdf_var_sngl(ncid,no,1,1,"err",tmp)    
    err(:)=REAL(tmp(:),r_size)

    status=NF90_CLOSE(ncid)
    call handle_error_netcdf("NF90_CLOSE: "//trim(filename),status)

  END SUBROUTINE read_obs

  !--------------------------------------------------------------
  ! For Innovation Statistics |
  !--------------------------------------------------------------

  SUBROUTINE make_ncfile_innovation

    USE common_setting, only: nbv,nobs
    USE NETCDF
    IMPLICIT NONE

    INTEGER,PARAMETER :: ndim=2

    INTEGER status
    INTEGER ncid,dimid,varid
    INTEGER dim(ndim)

    CHARACTER(10) :: filename="inv.nc"

    status=NF90_CREATE(trim(filename),nf90_netcdf4,ncid)
    CALL handle_error_netcdf("NF90_CREATE: "//trim(filename),status)

    status=NF90_PUT_ATT(ncid,NF90_GLOBAL,"title","Innovation Statistics")
    CALL handle_error_netcdf("NF90_PUT_ATT: "//trim(filename),status)

    status=NF90_PUT_ATT(ncid,NF90_GLOBAL,"description","Innovation Statistics")
    CALL handle_error_netcdf("NF90_PUT_ATT: "//trim(filename),status)
    
    !#Obs.
    status=NF90_DEF_DIM(ncid,"nobs",nobs,dimid)
    CALL handle_error_netcdf("NF90_DEF_DIM: "//trim(filename),status)
    dim(1)=dimid

    !#Ensemble member
    status=NF90_DEF_DIM(ncid,"nmem",nbv,dimid)
    CALL handle_error_netcdf("NF90_DEF_DIM: "//trim(filename),status)
    dim(2)=dimid
    
    !1D
    call define_var_netcdf(ncid,1,dim(1),varid,"int", &
         & "ele","element","element (h:2567, u:2819, v:2820, t:3073, s:3332)")

    call define_var_netcdf(ncid,1,dim(1),varid,"int", &
         & "ins","instrument", &
         & "instrument (1: AMSR-E, 2: WindSAT, 3: AMSR-2, 4: Himawari,"// &
         & " 11: SMOS, 12: SMAP,"// &
         & " 21: SSHA,"// &
         & " 31: Motion vector,"// &
         & " 0: Others"// &
         & " -1: XBT/MBT, -2:XCTD, -3:CTD, -4: Ship, -5: Profiling float,"// &
         & " -6: Drifter buoy, -7: Mooring buoy, -8: Animal)")
    
    call define_var_netcdf(ncid,1,dim(1),varid,"real", &
         & "lon","longitude","degree E")

    call define_var_netcdf(ncid,1,dim(1),varid,"real", &
         & "lat","latitude","degree N")

    call define_var_netcdf(ncid,1,dim(1),varid,"real", &
         & "lev","level","meter")
    
    call define_var_netcdf(ncid,1,dim(1),varid,"real", &
         & "obs","observation","H: meter, U and V: m/s, T: degree C, S: -")

    call define_var_netcdf(ncid,1,dim(1),varid,"real", &
         & "err","observation error","H: meter, U and V: m/s, T: degree C, S: -")
    
    call define_var_netcdf(ncid,1,dim(1),varid,"real", &
         & "hxfmean","forecast ensemble mean in obs. space","H: meter, U and V: m/s, T: degree C, S: -")

    call define_var_netcdf(ncid,1,dim(1),varid,"real", &
         & "hxfsprd","forecast ensemble mean in obs. space","H: meter, U and V: m/s, T: degree C, S: -")
    
    call define_var_netcdf(ncid,1,dim(1),varid,"real", &
         & "hxamean","analysis ensemble spread in obs. space","H: meter, U and V: m/s, T: degree C, S: -")

    call define_var_netcdf(ncid,1,dim(1),varid,"real", &
         & "hxasprd","analysis ensemble spread in obs. space","H: meter, U and V: m/s, T: degree C, S: -")
    
    !2D
    call define_var_netcdf(ncid,2,dim(1:2),varid,"real", &
         & "hxf","forecast ensemble","H: meter, U and V: m/s, T: degree C, S: -")

    call define_var_netcdf(ncid,2,dim(1:2),varid,"real", &
         & "hxa","analysis ensemble","H: meter, U and V: m/s, T: degree C, S: -")
    
    status=NF90_ENDDEF(ncid)
    CALL handle_error_netcdf("NF90_ENDDEF: "//trim(filename),status)

    status=NF90_CLOSE(ncid)
    CALL handle_error_netcdf("NF90_CLOSE: "//trim(filename),status)
    
  END SUBROUTINE make_ncfile_innovation

  !-------------------------------------

  SUBROUTINE write_innovation(hxf,hxa,hxfmean,hxfsprd,hxamean,hxasprd)

    USE common_setting    
    USE NETCDF
    IMPLICIT NONE

    !---Common
    INTEGER status,ncid

    CHARACTER(10) :: filename="inv.nc"
    
    !---IN
    REAL(r_size),INTENT(IN) :: hxf(nobs,nbv) !H(xf)
    REAL(r_size),INTENT(IN) :: hxa(nobs,nbv) !H(xa)
    REAL(r_size),INTENT(IN) :: hxfmean(nobs) !H(xfmean)
    REAL(r_size),INTENT(IN) :: hxfsprd(nobs) !H(xfsprd)
    REAL(r_size),INTENT(IN) :: hxamean(nobs) !H(xamean)
    REAL(r_size),INTENT(IN) :: hxasprd(nobs) !H(xasprd)

    status=NF90_OPEN(trim(filename),NF90_WRITE,ncid)
    call handle_error_netcdf("NF90_OPEN: "//trim(filename),status)
    
    CALL write_netcdf_var_int(ncid,nobs,1,1,"ele",obselm(:))      
    CALL write_netcdf_var_int(ncid,nobs,1,1,"ins",obsins(:))
    CALL write_netcdf_var_sngl(ncid,nobs,1,1,"lon",REAL(obslon(:),r_sngl))
    CALL write_netcdf_var_sngl(ncid,nobs,1,1,"lat",REAL(obslat(:),r_sngl))
    CALL write_netcdf_var_sngl(ncid,nobs,1,1,"lev",REAL(obslev(:),r_sngl))
    CALL write_netcdf_var_sngl(ncid,nobs,1,1,"obs",REAL(obsdat(:),r_sngl))
    CALL write_netcdf_var_sngl(ncid,nobs,1,1,"err",REAL(obserr(:),r_sngl))
    CALL write_netcdf_var_sngl(ncid,nobs,1,1,"hxfmean",REAL(hxfmean(:),r_sngl))
    CALL write_netcdf_var_sngl(ncid,nobs,1,1,"hxfsprd",REAL(hxfsprd(:),r_sngl))
    CALL write_netcdf_var_sngl(ncid,nobs,1,1,"hxamean",REAL(hxamean(:),r_sngl))
    CALL write_netcdf_var_sngl(ncid,nobs,1,1,"hxasprd",REAL(hxasprd(:),r_sngl))
    CALL write_netcdf_var_sngl(ncid,nobs,nbv,1,"hxf",REAL(hxf(:,:),r_sngl))
    CALL write_netcdf_var_sngl(ncid,nobs,nbv,1,"hxa",REAL(hxa(:,:),r_sngl))

    status=NF90_CLOSE(ncid)
    call handle_error_netcdf("NF90_CLOSE: "//trim(filename),status)
    
  END SUBROUTINE write_innovation
  
END MODULE common_io
