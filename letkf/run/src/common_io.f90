!--------------------------------------------------------------------------
! NetCDF Input/Output code
!--------------------------------------------------------------------------
!
! 06/04/2026 Shun OHISHI updated
!
!--------------------------------------------------------------------------

MODULE common_io

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

  !----------------------------------------------------------------------
  ! Read variable 1-3D |
  !----------------------------------------------------------------------

  SUBROUTINE read_netcdf_var_int_1d(ncid,im,varname,glb)

    USE NETCDF
    implicit none

    integer status,varid

    integer,intent(in) :: ncid
    integer,intent(in) :: im

    character(*),intent(in) :: varname

    integer,intent(out) :: glb(im)

    status=nf90_inq_varid(ncid,varname,varid)
    call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status)
    status=nf90_get_var(ncid,varid,glb,(/1/),(/im/))
    call handle_error_netcdf('nf90_get_var:'//trim(varname),status)

  END SUBROUTINE read_netcdf_var_int_1d

  !-----------------------------------

  SUBROUTINE read_netcdf_var_int_2d(ncid,im,jm,varname,glb)

    USE NETCDF
    implicit none

    integer status,varid

    integer,intent(in) :: ncid
    integer,intent(in) :: im,jm

    character(*),intent(in) :: varname

    integer,intent(out) :: glb(im,jm)

    status=nf90_inq_varid(ncid,varname,varid)
    call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status)
    status=nf90_get_var(ncid,varid,glb,(/1,1/),(/im,jm/))
    call handle_error_netcdf('nf90_get_var:'//trim(varname),status)

  END SUBROUTINE read_netcdf_var_int_2d

  !---------------------------------

  SUBROUTINE read_netcdf_var_int_3d(ncid,im,jm,km,varname,glb)

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

  END SUBROUTINE read_netcdf_var_int_3d

  !---------------------------------

  SUBROUTINE read_netcdf_var_sngl_1d(ncid,im,varname,glb)

    USE common_setting, only:r_sngl
    USE NETCDF
    implicit none

    integer status,varid

    integer,intent(in) :: ncid
    integer,intent(in) :: im

    character(*),intent(in) :: varname

    real(kind = r_sngl),intent(out) :: glb(im)

    status=nf90_inq_varid(ncid,varname,varid)
    call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status)
    status=nf90_get_var(ncid,varid,glb,(/1/),(/im/))
    call handle_error_netcdf('nf90_get_var:'//trim(varname),status)

  END SUBROUTINE read_netcdf_var_sngl_1d

  !----------------------------------

  SUBROUTINE read_netcdf_var_sngl_2d(ncid,im,jm,varname,glb)

    USE common_setting, only:r_sngl
    USE NETCDF
    implicit none

    integer status,varid

    integer,intent(in) :: ncid
    integer,intent(in) :: im,jm

    character(*),intent(in) :: varname

    real(kind = r_sngl),intent(out) :: glb(im,jm)

    status=nf90_inq_varid(ncid,varname,varid)
    call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status)
    status=nf90_get_var(ncid,varid,glb,(/1,1/),(/im,jm/))
    call handle_error_netcdf('nf90_get_var:'//trim(varname),status)

  END SUBROUTINE read_netcdf_var_sngl_2d

  !---------------------------------

  SUBROUTINE read_netcdf_var_sngl_3d(ncid,im,jm,km,varname,glb)

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

  END SUBROUTINE read_netcdf_var_sngl_3d

  !-----------------------------------

  SUBROUTINE read_netcdf_var_dble_1d(ncid,im,varname,glb)

    USE common_setting, only:r_dble    
    USE NETCDF
    implicit none

    integer status,varid

    integer,intent(in) :: ncid
    integer,intent(in) :: im

    character(*),intent(in) :: varname

    real(kind = r_dble),intent(out) :: glb(im)

    status=nf90_inq_varid(ncid,varname,varid)
    call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status)
    status=nf90_get_var(ncid,varid,glb,(/1/),(/im/))
    call handle_error_netcdf('nf90_get_var:'//trim(varname),status)

  END SUBROUTINE read_netcdf_var_dble_1d

  !-----------------------------------

  SUBROUTINE read_netcdf_var_dble_2d(ncid,im,jm,varname,glb)

    USE common_setting, only:r_dble    
    USE NETCDF
    implicit none

    integer status,varid

    integer,intent(in) :: ncid
    integer,intent(in) :: im,jm

    character(*),intent(in) :: varname

    real(kind = r_dble),intent(out) :: glb(im,jm)

    status=nf90_inq_varid(ncid,varname,varid)
    call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status)
    status=nf90_get_var(ncid,varid,glb,(/1,1/),(/im,jm/))
    call handle_error_netcdf('nf90_get_var:'//trim(varname),status)

  END SUBROUTINE read_netcdf_var_dble_2d

  !-----------------------------------

  SUBROUTINE read_netcdf_var_dble_3d(ncid,im,jm,km,varname,glb)

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

  END SUBROUTINE read_netcdf_var_dble_3d

  !----------------------------------------------------------------------
  ! Write variable 1-3D |
  !----------------------------------------------------------------------

  SUBROUTINE write_netcdf_var_int_1d(ncid,im,varname,glb)

    USE NETCDF
    implicit none

    integer status,varid

    integer,intent(in) :: ncid
    integer,intent(in) :: im

    character(*),intent(in) :: varname

    integer,intent(in) :: glb(im)

    status=nf90_inq_varid(ncid,varname,varid)
    call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status)
    status=nf90_put_var(ncid,varid,glb,(/1/),(/im/))
    call handle_error_netcdf('nf90_put_var:'//trim(varname),status)

  END SUBROUTINE write_netcdf_var_int_1d

  !----------------------------------

  SUBROUTINE write_netcdf_var_int_2d(ncid,im,jm,varname,glb)

    USE NETCDF
    implicit none

    integer status,varid

    integer,intent(in) :: ncid
    integer,intent(in) :: im,jm

    character(*),intent(in) :: varname

    integer,intent(in) :: glb(im,jm)

    status=nf90_inq_varid(ncid,varname,varid)
    call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status)
    status=nf90_put_var(ncid,varid,glb,(/1,1/),(/im,jm/))
    call handle_error_netcdf('nf90_put_var:'//trim(varname),status)

  END SUBROUTINE write_netcdf_var_int_2d

  !----------------------------------

  SUBROUTINE write_netcdf_var_int_3d(ncid,im,jm,km,varname,glb)

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

  END SUBROUTINE write_netcdf_var_int_3d

  !----------------------------------

  SUBROUTINE write_netcdf_var_sngl_1d(ncid,im,varname,glb)

    USE common_setting, only:r_sngl
    USE NETCDF
    implicit none

    integer status,varid

    integer,intent(in) :: ncid
    integer,intent(in) :: im

    character(*),intent(in) :: varname

    real(kind = r_sngl),intent(in) :: glb(im)

    status=nf90_inq_varid(ncid,varname,varid)
    call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status)
    status=nf90_put_var(ncid,varid,glb,(/1/),(/im/))
    call handle_error_netcdf('nf90_put_var:'//trim(varname),status)

  END SUBROUTINE write_netcdf_var_sngl_1d

  !----------------------------------

  SUBROUTINE write_netcdf_var_sngl_2d(ncid,im,jm,varname,glb)

    USE common_setting, only:r_sngl
    USE NETCDF
    implicit none

    integer status,varid

    integer,intent(in) :: ncid
    integer,intent(in) :: im,jm

    character(*),intent(in) :: varname

    real(kind = r_sngl),intent(in) :: glb(im,jm)

    status=nf90_inq_varid(ncid,varname,varid)
    call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status)
    status=nf90_put_var(ncid,varid,glb,(/1,1/),(/im,jm/))
    call handle_error_netcdf('nf90_put_var:'//trim(varname),status)

  END SUBROUTINE write_netcdf_var_sngl_2d

  !----------------------------------

  SUBROUTINE write_netcdf_var_sngl_3d(ncid,im,jm,km,varname,glb)

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

  END SUBROUTINE write_netcdf_var_sngl_3d

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
    integer varid

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
  ! Read state vector 2-3D |
  !--------------------------------------------------------------------

  SUBROUTINE read_state_vector_2d(filename,v2d)

    !$ USE OMP_LIB
    USE common_setting
    USE common_mpi
    USE NETCDF
    USE MPI
    IMPLICIT NONE

    !Common
    INTEGER iv2d
    INTEGER status,access
    INTEGER ncid

    REAL(r_sngl),ALLOCATABLE :: tmp2d(:,:)

    !IN
    CHARACTER(12),INTENT(IN) :: filename

    !OUT
    REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)

    status=access(filename," ")
    IF(status /= 0)THEN
       WRITE(6,'(A)') " *** Error: Not Found "//trim(filename)
       call finalize_mpi
       stop
    END IF

    !---Open
    status = NF90_OPEN(trim(filename),NF90_NOWRITE,ncid)
    call handle_error_netcdf("NF90_OPEN: "//trim(filename),status)

    call detect_varname(ncid)

    !---2D
    ALLOCATE(tmp2d(nlon,nlat))
    DO iv2d=1,nv2d
       call read_netcdf_var_sngl_2d(ncid,nlon,nlat,trim(var2df(iv2d)),tmp2d)
       v2d(:,:,iv2d)=REAL(tmp2d(:,:),r_size)
    END DO
    DEALLOCATE(tmp2d)

    !---Close
    status = NF90_CLOSE(ncid)
    call handle_error_netcdf("NF90_CLOSE: "//trim(filename),status)

  END SUBROUTINE read_state_vector_2d

  !-----------------------------------

  SUBROUTINE read_state_vector_3d(filename,iv3d,v3d)

    !$ USE OMP_LIB
    USE common_setting
    USE common_mpi
    USE NETCDF
    USE MPI
    IMPLICIT NONE

    !Common
    INTEGER k
    INTEGER status,access
    INTEGER ncid

    REAL(r_sngl),ALLOCATABLE :: tmp3d(:,:,:)

    !IN
    CHARACTER(12),INTENT(IN) :: filename
    INTEGER,INTENT(IN) :: iv3d

    !OUT
    REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev)

    status=access(filename," ")
    IF(status / = 0)THEN
       WRITE(6,'(A)') " *** Error: Not Found "//trim(filename)
       call finalize_mpi
       stop
    END IF

    !---Open
    status = NF90_OPEN(trim(filename),NF90_NOWRITE,ncid)
    call handle_error_netcdf("NF90_OPEN: "//trim(filename),status)

    call detect_varname(ncid)

    !---3D
    ALLOCATE(tmp3d(nlon,nlat,nlev))
    call read_netcdf_var_sngl_3d(ncid,nlon,nlat,nlev,trim(var3df(iv3d)),tmp3d(:,:,:))

    ! POM -> ROMS type vertical order for TransXtoY
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(k)
    DO k=1,nlev
       v3d(:,:,nlev-k+1) = REAL(tmp3d(:,:,k),r_size)
    END DO
    !$OMP END PARALLEL DO
    DEALLOCATE(tmp3d)

    !---Close
    status = NF90_CLOSE(ncid)
    call handle_error_netcdf("NF90_CLOSE: "//trim(filename),status)

  END SUBROUTINE read_state_vector_3d

  !--------------------------------------------------

  SUBROUTINE write_state_vector_all(filename,v3d,v2d)

    !$ USE OMP_LIB
    USE common_setting
    USE common
    USE common_mpi
    USE NETCDF
    USE MPI
    IMPLICIT NONE

    INTEGER k
    INTEGER iv2d,iv3d
    INTEGER status,access
    INTEGER ncid

    REAL(r_sngl),ALLOCATABLE :: tmp3d(:,:,:)

    !--IN
    REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
    REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)

    CHARACTER(12),INTENT(IN) :: filename

    !---Access
    status=access(filename," ")
    IF(status /= 0)THEN
       WRITE(6,'(A)') " *** Error: Not Found "//trim(filename)
       call finalize_mpi
       stop
    END IF

    !---Open    
    status = NF90_OPEN(filename,NF90_WRITE,ncid)
    call handle_error_netcdf("NF90_OPEN: "//trim(filename),status)

    !---Variable name
    call detect_varname(ncid)

    !---2D
    DO iv2d=1,nv2d
       call write_netcdf_var_sngl_2d(ncid,nlon,nlat,trim(var2da(iv2d)),v2d(:,:,iv2d))
    END DO

    !---3D
    ALLOCATE(tmp3d(nlon,nlat,nlev))
    DO iv3d=1,nv3d

       !Restore vertical order
       !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(k)
       DO k=1,nlev
          tmp3d(:,:,nlev-k+1)=v3d(:,:,k,iv3d)
       END DO
       !$OMP END PARALLEL DO

       call write_netcdf_var_sngl_3d(ncid,nlon,nlat,nlev,trim(var3da(iv3d)),tmp3d(:,:,:))

    END DO
    DEALLOCATE(tmp3d)

    !---CLOSE
    status = NF90_CLOSE(ncid)
    call handle_error_netcdf("NF90_CLOSE: "//trim(filename),status)

  END SUBROUTINE write_state_vector_all

  !-----------------------------------------------------------------------
  ! Read ensemble data and distribute to processes
  !-----------------------------------------------------------------------

  SUBROUTINE read_ens_mpi(islot,no, &
       & tmpidx,tmpidy,tmpelm,tmplon,tmplat,tmplev, &
       & tmpqc,tmphdxf,fcst3d,fcst2d)

    !$ USE OMP_LIB
    USE common_setting
    USE common
    USE common_mpi
    USE MPI
    IMPLICIT NONE

    !---Common
    INTEGER ib,nb
    INTEGER ivd
    INTEGER imem,imem_master
    INTEGER io

    INTEGER :: irank_read_2d,irank_read_3d(nv3d)

    REAL(r_size),ALLOCATABLE :: v2dg(:,:,:)
    REAL(r_size),ALLOCATABLE :: v3dg(:,:,:,:)

    CHARACTER(12) :: fcstfile="faTTNNNNN.nc"    

    !---MPI
    INTEGER :: MPI_R_SIZE
    INTEGER :: iprocs,tag,ierr
    INTEGER,ALLOCATABLE :: req_send(:),req_recv(:)

    !---IN
    INTEGER,INTENT(IN) :: islot
    INTEGER,INTENT(IN) :: no

    INTEGER,INTENT(IN) :: tmpidx(no),tmpidy(no)
    INTEGER,INTENT(IN) :: tmpelm(no)

    REAL(r_size),INTENT(IN) :: tmplon(no),tmplat(no),tmplev(no)

    !---OUT
    INTEGER,INTENT(OUT) :: tmpqc(no)    !QC

    REAL(r_size),INTENT(OUT) :: tmphdxf(no,nmem) !dY=HdX Ensemble perturvation in obs. space
    REAL(r_size),INTENT(OUT) :: fcst3d(nij1,nlev,nmem,nv3d),fcst2d(nij1,nmem,nv2d)

    !---ALLOCATE
    ALLOCATE(v2dg(nlon,nlat,nv2d))
    ALLOCATE(v3dg(nlon,nlat,nlev,nv3d))
    ALLOCATE(req_send(nv3d))
    ALLOCATE(req_recv(nv3d))

    !---Initialization
    tmphdxf(:,:) = REAL(0.0d0, r_size)
    tmpqc(:)=1 !1:OK, 0:Undef or NG

    !---nb (number of block/loop)
    nb=CEILING(REAL(nmem)/REAL(nprocs))

    !---MPI(r_size)
    IF(r_size == kind(0.0d0))THEN
       MPI_R_SIZE=MPI_DOUBLE_PRECISION
    ELSE IF(r_size == kind(0.0e0))THEN
       MPI_R_SIZE=MPI_REAL
    END IF

    !---Main loop    
    DO ib=1,nb

       !---Initialization
       v2dg(:,:,:) = REAL(0.d0, r_size)
       v3dg(:,:,:,:) = REAL(0.d0, r_size)

       req_send(:)=MPI_REQUEST_NULL
       req_recv(:)=MPI_REQUEST_NULL

       !---imem & imem_master
       IF(ib == 1)THEN
          imem = MOD(myrank, nmem) + 1
       ELSE
          imem = myrank + (ib-1)*nprocs + 1 
       END IF

       IF(1 <= imem .and.  imem <= nmem)THEN

          imem_master = MOD(imem-1, nprocs)

          !---Forecast filename
          WRITE(fcstfile(3:9),'(I2.2,I5.5)') islot,imem

          !---Read rank definition
          irank_read_2d = imem_master

          DO ivd=1,nv3d
             irank_read_3d(ivd) = MOD( imem-1 + ivd*nmem, nprocs)
          END DO

          !---Read 2D-ensemble forecast
          IF(myrank == irank_read_2d)THEN
             CALL read_state_vector_2d(fcstfile,v2dg)
          ENDIF

          !---Read 3D-ensemble forecast
          !MPI_IRECV
          IF(myrank == imem_master)THEN
             DO ivd=1,nv3d

                IF(myrank == imem_master .and. irank_read_3d(ivd) /= imem_master)THEN
                   tag=(imem-1)*nv3d+ivd
                   CALL MPI_IRECV(v3dg(:,:,:,ivd), nlon*nlat*nlev, MPI_R_SIZE, irank_read_3d(ivd), tag, &
                        & MPI_COMM_WORLD, req_recv(ivd), ierr)
                END IF

             END DO !ivd
          END IF

          !Read data & MPI_ISEND
          DO ivd=1,nv3d

             IF(myrank == irank_read_3d(ivd))THEN

                CALL read_state_vector_3d(fcstfile,ivd,v3dg(:,:,:,ivd))

                IF(myrank /= imem_master)THEN
                   tag=(imem-1)*nv3d+ivd
                   CALL MPI_ISEND(v3dg(:,:,:,ivd), nlon*nlat*nlev, MPI_R_SIZE, imem_master, tag, &
                        & MPI_COMM_WORLD, req_send(ivd), ierr)
                END IF

             END IF

          END DO !ivd

          !MPI_WAIT
          DO ivd=1,nv3d

             IF(myrank == imem_master .and. req_recv(ivd) /= MPI_REQUEST_NULL)THEN
                CALL MPI_WAIT(req_recv(ivd),MPI_STATUS_IGNORE,ierr)
             END IF

             IF(req_send(ivd) /= MPI_REQUEST_NULL)THEN
                CALL MPI_WAIT(req_send(ivd),MPI_STATUS_IGNORE,ierr)
             END IF

          END DO !ivd

          !---Observation operator and tmoqc
          IF(myrank == imem_master)THEN

             !---Observation operator H(xf(imem))
             CALL Trans_XtoY(no, &
                  & tmpelm,tmpidx,tmpidy,tmplon,tmplat,tmplev, &
                  & v3dg,v2dg,tmphdxf(:,imem))

             !---tmpqc
             !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(io)
             DO io=1,no
                IF(tmphdxf(io,imem) == undef)THEN
                   tmpqc(io) = 0
                ENDIF
             END DO !io
             !$OMP END PARALLEL DO

          END IF

       END IF ! 1 <= imem <= nmem

       !---fcst3d,2d
       DO iprocs=0,nprocs-1
          imem=iprocs+1+(ib-1)*nprocs
          IF(islot == 1 .and. 1 <= imem .and. imem <= nmem)THEN
             CALL scatter_grd_mpi(iprocs,REAL(v3dg,r_sngl),REAL(v2dg,r_sngl), &
                  & fcst3d(:,:,imem,:),fcst2d(:,imem,:))
          END IF
       END DO

    END DO !ib

  END SUBROUTINE read_ens_mpi

  !-----------------------------------------------------------------------
  ! Write ensemble data after collecting data from processes
  !-----------------------------------------------------------------------

  SUBROUTINE write_ens_mpi(v3df,v2df,v3da,v2da)

    !$ USE OMP_LIB
    USE MPI
    USE common_setting
    USE common
    USE common_mpi
    IMPLICIT NONE

    !---Common
    INTEGER,PARAMETER :: ntask=nmem+4
    ![1,nmem]: ensemble analysis
    ![nmem+1]: anal_mean, [nmem+2]: anal_sprd
    ![nmem+3]: fcst_mean, [nmem+4]: fcst_sprd

    INTEGER i,n
    INTEGER iprocs,imem,iobs
    INTEGER iprocs1,iprocs2,iprocs3,iprocs4
    INTEGER MPI_R_SIZE,status(MPI_STATUS_SIZE),nosend,ierr

    REAL(r_size) v3df_m(nij1,nlev,nv3d),v2df_m(nij1,nv2d)       !Forecast ensemble mean
    REAL(r_size) v3df_s(nij1,nlev,nv3d),v2df_s(nij1,nv2d)       !                  spread
    REAL(r_size) v3da_m(nij1,nlev,nv3d),v2da_m(nij1,nv2d)       !Analysis ensemble mean
    REAL(r_size) v3da_s(nij1,nlev,nv3d),v2da_s(nij1,nv2d)       !                  spread

    REAL(r_sngl),ALLOCATABLE :: v3dg(:,:,:,:),v2dg(:,:,:) !Global domain

    REAL(r_size),ALLOCATABLE :: hxa(:,:)
    REAL(r_size),ALLOCATABLE :: hxa_k(:)
    REAL(r_size) hxamean(nobs),hxasprd(nobs)
    REAL(r_size) hxfmean(nobs),hxfsprd(nobs)
        
    CHARACTER(12) :: filename='file00000.nc'

    !---IN    
    REAL(r_size),INTENT(IN) :: v3df(nij1,nlev,nmem,nv3d),v2df(nij1,nmem,nv2d) !Forecast (local)
    REAL(r_size),INTENT(IN) :: v3da(nij1,nlev,nmem,nv3d),v2da(nij1,nmem,nv2d) !Analysis (local)

    !---Initialization
    ALLOCATE(v3dg(nlon,nlat,nlev,nv3d),v2dg(nlon,nlat,nv2d))
    ALLOCATE(hxa_k(nobs))
    hxa_k(:)=REAL(0.d0,r_size)

    nosend=0
    
    IF(myrank == root_out)THEN
       ALLOCATE(hxa(nobs,nmem))
       hxa(:,:) = REAL(0.d0, r_size)
    END IF
    
    IF(r_size == kind(0.d0))THEN
       MPI_R_SIZE=MPI_DOUBLE_PRECISION
    ELSE IF(r_size == kind(0.e0))THEN
       MPI_R_SIZE=MPI_REAL
    END IF

    !---Ensemble mean/sprd
    call ensemble_mesp(v3da,v2da,v3da_m,v3da_s,v2da_m,v2da_s)
    call ensemble_mesp(v3df,v2df,v3df_m,v3df_s,v2df_m,v2df_s)

    !---Repeat index
    n = CEILING(REAL(ntask)/REAL(nprocs))

    !---Analysis ensemble, and forecast/analysis ensemble mean and spread (Output) --> Obs. space
    DO i=1,n

       !GATHER to iprocs
       DO iprocs=0,nprocs-1

          imem = iprocs+1 + (i-1)*nprocs

          IF(1 <= imem .and. imem <= nmem)THEN
             !Analysis ensemble
             CALL gather_grd_mpi(iprocs,v3da(:,:,imem,:),v2da(:,imem,:),v3dg,v2dg)
          ELSE IF(imem == nmem+1)THEN
             !Analysis ensemble mean
             iprocs1=iprocs
             CALL gather_grd_mpi(iprocs,v3da_m,v2da_m,v3dg,v2dg)
          ELSE IF(imem == nmem+2)THEN
             !Analysis ensemble spread
             iprocs2=iprocs
             CALL gather_grd_mpi(iprocs,v3da_s,v2da_s,v3dg,v2dg)
          ELSE IF(imem == nmem+3)THEN
             !Forecast ensemble mean
             iprocs3=iprocs
             CALL gather_grd_mpi(iprocs,v3df_m,v2df_m,v3dg,v2dg)
          ELSE IF(imem == nmem+4)THEN
             !Forecast Ensemble spread
             iprocs4=iprocs
             CALL gather_grd_mpi(iprocs,v3df_s,v2df_s,v3dg,v2dg)
          END IF

       END DO !iprocs

       !---Post process for quasi-global
       IF(lqglobal)THEN
          v2dg(nlon-1:nlon,:,:)=v2dg(1:2,:,:)
          v3dg(nlon-1:nlon,:,:,:)=v3dg(1:2,:,:,:)
       END IF

       imem = myrank+1 + (i-1)*nprocs

       !Observation space & MPI_ISEND
       IF(1 <= imem .and. imem <= nmem)THEN !Analysis ensemble in obs. space: H(x^a(k))
          CALL Trans_XtoY(nobs,obselm,obsidx,obsidy,obslon,obslat,obslev,v3dg,v2dg,hxa_k)

          IF(myrank == root_out)THEN
             hxa(:,imem)=hxa_k(:)
             nosend=nosend+1
          ELSE
             CALL MPI_SEND(hxa_k,nobs,MPI_R_SIZE,root_out,imem,MPI_COMM_WORLD,ierr)
          END IF
          
       ELSE IF(imem == nmem+1)THEN !Analysis ensemble mean
          CALL Trans_XtoY(nobs,obselm,obsidx,obsidy,obslon,obslat,obslev,v3dg,v2dg,hxamean)
          CALL MPI_SEND(hxamean,nobs,MPI_R_SIZE,root_out,imem,MPI_COMM_WORLD,ierr)   
       ELSE IF(imem == nmem+2)THEN !Analysis ensemble spread
          CALL Trans_XtoY(nobs,obselm,obsidx,obsidy,obslon,obslat,obslev,v3dg,v2dg,hxasprd)
          CALL MPI_SEND(hxasprd,nobs,MPI_R_SIZE,root_out,imem,MPI_COMM_WORLD,ierr)   
       ELSE IF(imem == nmem+3)THEN !Forecast ensemble mean
          CALL Trans_XtoY(nobs,obselm,obsidx,obsidy,obslon,obslat,obslev,v3dg,v2dg,hxfmean)
          CALL MPI_SEND(hxfmean,nobs,MPI_R_SIZE,root_out,imem,MPI_COMM_WORLD,ierr)   
       ELSE IF(imem == nmem+4)THEN !Forecast ensemble spread
          CALL Trans_XtoY(nobs,obselm,obsidx,obsidy,obslon,obslat,obslev,v3dg,v2dg,hxfsprd)
          CALL MPI_SEND(hxfsprd,nobs,MPI_R_SIZE,root_out,imem,MPI_COMM_WORLD,ierr)   
       END IF
       
       !Filename
       IF(1 <= imem .and. imem <= nmem)THEN
          WRITE(filename(1:9),'(A4,I5.5)') "fa01",imem
       ELSE IF(imem == nmem+1)THEN
          WRITE(filename(1:9),'(A9)') "anal_mean"
       ELSE IF(imem == nmem+2)THEN
          WRITE(filename(1:9),'(A9)') "anal_sprd"
       ELSE IF(imem == nmem+3)THEN
          WRITE(filename(1:9),'(A9)') "fcst_mean"
       ELSE IF(imem == nmem+4)THEN
          WRITE(filename(1:9),'(A9)') "fcst_sprd"
       END IF

       !Write
       IF(1 <= imem .and. imem <= ntask)THEN
          CALL write_state_vector_all(filename,v3dg,v2dg)
       END IF
       
    END DO !i

    !MPI_RECV(hxa,hxamean,hxasprd,hxfmean,hxfsprd)
    IF(myrank == root_out)THEN
       DO i=1,ntask-nosend
          
          CALL MPI_RECV(hxa_k,nobs,MPI_R_SIZE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
          imem=status(MPI_TAG)

          IF(1 <= imem .and. imem <= nmem)THEN
             hxa(:,imem)=hxa_k(:)
          ELSE IF(imem == nmem+1)THEN
             hxamean(:)=hxa_k(:)
          ELSE IF(imem == nmem+2)THEN
             hxasprd(:)=hxa_k(:)
          ELSE IF(imem == nmem+3)THEN
             hxfmean(:)=hxa_k(:)
          ELSE IF(imem == nmem+4)THEN
             hxfsprd(:)=hxa_k(:)
          END IF
             
       END DO       
    END IF
    
    !---Deallocate
    DEALLOCATE(v3dg,v2dg)
    DEALLOCATE(hxa_k)
        
    !---Forecast ensemble matrix
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iobs)
    DO iobs=1,nobs
       obshdxf(iobs,:)=obshdxf(iobs,:)+obshxfmean(iobs)
    END DO
    !$OMP END PARALLEL DO
    
    !---Monitor
    IF(myrank == root_out)THEN
       CALL monit_mean("fcst",hxfmean)
       CALL monit_mean("anal",hxamean)
       CALL monit_sprd("fcst",hxfsprd)
       CALL monit_sprd("anal",hxasprd)
    END IF

    !---Make & Write innovation netcdf file
    IF(myrank == root_out)THEN
       CALL make_ncfile_innovation
       CALL write_innovation(hxfmean,hxfsprd,hxamean,hxasprd,hxa)
    END IF

    !---Deallocate
    IF(myrank == root_out)THEN
       DEALLOCATE( hxa )
    END IF

    DEALLOCATE( obsidx, obsidy, obselm, obsins )
    DEALLOCATE( obslon, obslat, obslev )
    DEALLOCATE( obsdat, obserr )
    DEALLOCATE( obshxfmean, obshxfsprd )
    DEALLOCATE( obsdep, obshdxf) 
    
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
    INTEGER ierr

    !IN
    CHARACTER(8),INTENT(IN) :: filename

    !OUT
    INTEGER,INTENT(OUT) :: no

    status=0
    IF(myrank == root)THEN
       status=access(filename," ")
    END IF

    CALL MPI_BCAST(status,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

    IF(status /= 0)THEN
       WRITE(6,'(A)') " *** Error: Not Found "//filename
       call finalize_mpi
       stop
    END IF

    IF(myrank == root)THEN

       status=NF90_OPEN(trim(filename),NF90_NOWRITE,ncid)
       call handle_error_netcdf("NF90_OPEN: "//trim(filename),status)

       status=NF90_INQ_DIMID(ncid,"nobs",dimid)
       call handle_error_netcdf("NF90_INQ_DIMID: nobs",status)

       status=NF90_INQUIRE_DIMENSION(ncid,dimid,len = no)    
       call handle_error_netcdf("NF90_INQUIRE_DIMENSION: nobs",status)

       status=NF90_CLOSE(ncid)
       call handle_error_netcdf("NF90_CLOSE: "//trim(filename),status)

    END IF

    CALL MPI_BCAST(no,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

  END SUBROUTINE get_nobs

  !--------------------------------------------------------------------
  ! Read observation
  !--------------------------------------------------------------------

  SUBROUTINE read_obs(filename,no,ele,ins,lon,lat,lev,obs,err)

    !$ USE OMP_LIB
    USE common_setting, only: r_sngl,r_size, myrank, root
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

    status=0
    IF(myrank == root)THEN
       status=access(filename," ")
    END IF

    CALL MPI_BCAST(status,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

    IF(status /= 0)THEN
       WRITE(6,'(A)') " *** Error: Not Found "//filename
       call finalize_mpi
       stop
    END IF

    IF(myrank == root)THEN

       !open
       status=NF90_OPEN(trim(filename),NF90_NOWRITE,ncid)
       call handle_error_netcdf("NF90_OPEN: "//trim(filename),status)


       !element
       call read_netcdf_var_int_1d(ncid,no,"ele",itmp)    
       ele(:)=itmp(:)

       !instrument
       call read_netcdf_var_int_1d(ncid,no,"ins",itmp)    
       ins(:)=itmp(:)

       !longitude
       call read_netcdf_var_sngl_1d(ncid,no,"lon",tmp)    
       lon(:)=REAL(tmp(:),r_size)

       !latitude
       call read_netcdf_var_sngl_1d(ncid,no,"lat",tmp)    
       lat(:)=REAL(tmp(:),r_size)

       !level
       call read_netcdf_var_sngl_1d(ncid,no,"lev",tmp)    
       lev(:)=REAL(tmp(:),r_size)

       !data
       call read_netcdf_var_sngl_1d(ncid,no,"dat",tmp)    
       obs(:)=REAL(tmp(:),r_size)

       !Error
       call read_netcdf_var_sngl_1d(ncid,no,"err",tmp)    
       err(:)=REAL(tmp(:),r_size)

       status=NF90_CLOSE(ncid)
       call handle_error_netcdf("NF90_CLOSE: "//trim(filename),status)

    END IF

    CALL MPI_BCAST(ele,no,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(ins,no,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
    CALL bcast_mpi_1d(root,no,lon)
    CALL bcast_mpi_1d(root,no,lat)
    CALL bcast_mpi_1d(root,no,lev)
    CALL bcast_mpi_1d(root,no,obs)
    CALL bcast_mpi_1d(root,no,err)

  END SUBROUTINE read_obs

  !--------------------------------------------------------------
  ! For Innovation Statistics |
  !--------------------------------------------------------------

  SUBROUTINE make_ncfile_innovation

    USE common_setting, only: nmem,nobs
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
    status=NF90_DEF_DIM(ncid,"nmem",nmem,dimid)
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
         & "hxf","forecast ensemble matrix in obs. space","H: meter, U and V: m/s, T: degree C, S: -")

    call define_var_netcdf(ncid,2,dim(1:2),varid,"real", &
         & "hxa","analysis ensemble matrix in obs. space","H: meter, U and V: m/s, T: degree C, S: -")

    status=NF90_ENDDEF(ncid)
    CALL handle_error_netcdf("NF90_ENDDEF: "//trim(filename),status)

    status=NF90_CLOSE(ncid)
    CALL handle_error_netcdf("NF90_CLOSE: "//trim(filename),status)

  END SUBROUTINE make_ncfile_innovation

  !-------------------------------------

  SUBROUTINE write_innovation(hxfmean,hxfsprd,hxamean,hxasprd,hxa)

    USE common_setting    
    USE NETCDF
    IMPLICIT NONE

    !---Common
    INTEGER status,ncid

    CHARACTER(10) :: filename="inv.nc"

    !---IN
    REAL(r_size),INTENT(IN) :: hxfmean(nobs) !H(xfmean)
    REAL(r_size),INTENT(IN) :: hxfsprd(nobs) !H(xfsprd)
    REAL(r_size),INTENT(IN) :: hxamean(nobs) !H(xamean)
    REAL(r_size),INTENT(IN) :: hxasprd(nobs) !H(xasprd)
    REAL(r_size),INTENT(IN) :: hxa(nobs,nmem) !H(xa)

    
    status=NF90_OPEN(trim(filename),NF90_WRITE,ncid)
    call handle_error_netcdf("NF90_OPEN: "//trim(filename),status)

    CALL write_netcdf_var_int_1d(ncid,nobs,"ele",obselm(:))      
    CALL write_netcdf_var_int_1d(ncid,nobs,"ins",obsins(:))
    CALL write_netcdf_var_sngl_1d(ncid,nobs,"lon",REAL(obslon(:),r_sngl))
    CALL write_netcdf_var_sngl_1d(ncid,nobs,"lat",REAL(obslat(:),r_sngl))
    CALL write_netcdf_var_sngl_1d(ncid,nobs,"lev",REAL(obslev(:),r_sngl))
    CALL write_netcdf_var_sngl_1d(ncid,nobs,"obs",REAL(obsdat(:),r_sngl))
    CALL write_netcdf_var_sngl_1d(ncid,nobs,"err",REAL(obserr(:),r_sngl))
    CALL write_netcdf_var_sngl_1d(ncid,nobs,"hxfmean",REAL(hxfmean(:),r_sngl))
    CALL write_netcdf_var_sngl_1d(ncid,nobs,"hxfsprd",REAL(hxfsprd(:),r_sngl))
    CALL write_netcdf_var_sngl_1d(ncid,nobs,"hxamean",REAL(hxamean(:),r_sngl))
    CALL write_netcdf_var_sngl_1d(ncid,nobs,"hxasprd",REAL(hxasprd(:),r_sngl))
    CALL write_netcdf_var_sngl_2d(ncid,nobs,nmem,"hxf",REAL(obshdxf(:,:),r_sngl))
    CALL write_netcdf_var_sngl_2d(ncid,nobs,nmem,"hxa",REAL(hxa(:,:),r_sngl))

    status=NF90_CLOSE(ncid)
    call handle_error_netcdf("NF90_CLOSE: "//trim(filename),status)

  END SUBROUTINE write_innovation
    
END MODULE common_io
