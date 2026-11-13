!--------------------------------------------------------------------------
! NetCDF Input/Output code
!--------------------------------------------------------------------------
!
! 2026/06/04 Shun OHISHI updated
! 2026/10/21 Shun OHISHI updated
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

    !---Common
    integer ierr

    !---IN
    integer,intent(in) :: status

    character(*),intent(in) :: routine

    if(status /= nf90_noerr) then
       write(6,'(A)') "Error: NetCDF "//trim(routine)
       write(6,'(A)') trim(nf90_strerror(status))
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if

  END SUBROUTINE handle_error_netcdf

  !----------------------------------------------------------------------
  ! Read variable 1-3D |
  !----------------------------------------------------------------------

  SUBROUTINE read_netcdf_var_int_1d(ncid,im,varname,glb)

    USE NETCDF
    implicit none

    !---Common
    integer status,varid

    !---IN
    integer,intent(in) :: ncid
    integer,intent(in) :: im

    character(*),intent(in) :: varname

    !---OUT
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

    !---Common
    integer status,varid

    !---IN    
    integer,intent(in) :: ncid
    integer,intent(in) :: im,jm

    character(*),intent(in) :: varname

    !---OUT    
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

    !---Common    
    integer status,varid

    !---IN    
    integer,intent(in) :: ncid
    integer,intent(in) :: im,jm,km

    character(*),intent(in) :: varname

    !---OUT
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

    !---Common    
    integer status,varid


    !---IN
    integer,intent(in) :: ncid
    integer,intent(in) :: im

    character(*),intent(in) :: varname

    !---OUT    
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

    !---Common    
    integer status,varid

    !---IN    
    integer,intent(in) :: ncid
    integer,intent(in) :: im,jm

    character(*),intent(in) :: varname

    !---OUT    
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

    !---Common    
    integer status,varid

    !---IN    
    integer,intent(in) :: ncid
    integer,intent(in) :: im,jm,km

    character(*),intent(in) :: varname

    !---OUT    
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

    !---Common    
    integer status,varid

    !---IN    
    integer,intent(in) :: ncid
    integer,intent(in) :: im

    character(*),intent(in) :: varname

    !---OUT    
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

    !---Common    
    integer status,varid

    !---IN
    integer,intent(in) :: ncid
    integer,intent(in) :: im,jm

    character(*),intent(in) :: varname

    !---OUT        
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

    !---Common    
    integer status,varid

    !---IN
    integer,intent(in) :: ncid
    integer,intent(in) :: im,jm,km

    character(*),intent(in) :: varname

    !---OUT    
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

    !---Common    
    integer status,varid

    !---IN
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

    !---Common
    integer status,varid

    !---IN
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

    !---Common
    integer status,varid

    !---IN
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

    !---Common
    integer status,varid

    !---IN
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

    !---Common
    integer status,varid

    !---IN
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

    !---Common
    integer status,varid

    !---IN
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
    use MPI
    implicit none

    !---Common
    integer status
    integer ierr

    !---IN
    integer,intent(in) :: ncid
    integer,intent(in) :: ndim,dim(ndim)

    character(*),intent(in) :: type
    character(*),intent(in) :: name,long_name,units_name

    !---OUT
    integer,intent(out) :: varid

    if(trim(type) == "int")then
       status=NF90_DEF_VAR(ncid,trim(name),nf90_int,dim,varid)
    else if(trim(type) == "real")then
       status=NF90_DEF_VAR(ncid,trim(name),nf90_float,dim,varid)
    else if(trim(type) == "dble")then
       status=NF90_DEF_VAR(ncid,trim(name),nf90_double,dim,varid)
    else
       status=NF90_EBADTYPE
       write(*,*) "***Error: type is not int, real, dble in define_var_netcdf ==> "//trim(type)
       CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
    call handle_error_netcdf("NF90_DEF_VAR: "//trim(name),status)

    status=NF90_DEF_VAR_DEFLATE(ncid,varid,shuffle=1,deflate=1,deflate_level=1)
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

    USE common_mpi
    USE common_setting
    USE NETCDF
    USE MPI
    implicit none

    !---Common
    integer status1,status2
    integer varid
    integer ierr

    !---IN
    integer,intent(in) :: ncid

    !---File type
    status1 = NF90_INQ_VARID(ncid,"elb",varid)
    status2 = NF90_INQ_VARID(ncid,"el_fcst",varid)
    if(status1 == 0)then !Restart/Intermittent
       file_type=1       
    else if(status2 == 0)then !IAU
       file_type=2
    else
       write(6,'(A)') trim(nf90_strerror(status1))
       write(6,'(A)') trim(nf90_strerror(status2))
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if

    !---Variable name
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
    USE common_mpi
    USE common_setting
    USE NETCDF
    USE MPI
    IMPLICIT NONE

    !---Common
    INTEGER iv2d
    INTEGER status,access
    INTEGER ncid
    integer ierr

    REAL(r_sngl),ALLOCATABLE :: tmp2d(:,:)

    !---IN
    CHARACTER(12),INTENT(IN) :: filename

    !---OUT
    REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)

    status=access(filename," ")
    IF(status /= 0)THEN
       WRITE(6,'(A)') " *** Error: Not Found "//trim(filename)
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    END IF

    !---Open
    status = NF90_OPEN(trim(filename),NF90_NOWRITE,ncid)
    call handle_error_netcdf("NF90_OPEN: "//trim(filename),status)

    call detect_varname(ncid)

    !---Read all 2D data
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
    USE common_mpi
    USE common_setting
    USE NETCDF
    USE MPI
    IMPLICIT NONE

    !---Common
    INTEGER status,access
    INTEGER ncid
    integer ierr

    REAL(r_sngl),ALLOCATABLE :: tmp3d(:,:,:)

    !---IN
    CHARACTER(12),INTENT(IN) :: filename
    INTEGER,INTENT(IN) :: iv3d

    !---OUT
    REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev)

    status=access(filename," ")
    IF(status /= 0)THEN
       WRITE(6,'(A)') " *** Error: Not Found "//trim(filename)
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    END IF

    !---Open
    status = NF90_OPEN(trim(filename),NF90_NOWRITE,ncid)
    call handle_error_netcdf("NF90_OPEN: "//trim(filename),status)

    call detect_varname(ncid)

    !---Read 3D data
    ALLOCATE(tmp3d(nlon,nlat,nlev))
    call read_netcdf_var_sngl_3d(ncid,nlon,nlat,nlev,trim(var3df(iv3d)),tmp3d)
    v3d(:,:,:)=REAL(tmp3d(:,:,:),r_size)
    DEALLOCATE(tmp3d)

    !---Close
    status = NF90_CLOSE(ncid)
    call handle_error_netcdf("NF90_CLOSE: "//trim(filename),status)

  END SUBROUTINE read_state_vector_3d

  !--------------------------------------------------

  SUBROUTINE write_state_vector_all(filename,v3d,v2d)

    !$ USE OMP_LIB
    USE common
    USE common_mpi
    USE common_setting
    USE NETCDF
    USE MPI
    IMPLICIT NONE

    !---Common
    INTEGER iv2d,iv3d
    INTEGER status,access
    INTEGER ncid
    integer ierr

    !--IN
    REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
    REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)

    CHARACTER(12),INTENT(IN) :: filename

    !---Access
    status=access(filename," ")
    IF(status /= 0)THEN
       WRITE(6,'(A)') " *** Error: Not Found "//trim(filename)
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    END IF

    !---Open    
    status = NF90_OPEN(filename,NF90_WRITE,ncid)
    call handle_error_netcdf("NF90_OPEN: "//trim(filename),status)

    !---Variable name
    call detect_varname(ncid)

    !---Write all 2D data
    DO iv2d=1,nv2d
       call write_netcdf_var_sngl_2d(ncid,nlon,nlat,trim(var2da(iv2d)),v2d(:,:,iv2d))
    END DO

    !---Write all 3D data
    DO iv3d=1,nv3d
       call write_netcdf_var_sngl_3d(ncid,nlon,nlat,nlev,trim(var3da(iv3d)),v3d(:,:,:,iv3d))
    END DO

    !---CLOSE
    status = NF90_CLOSE(ncid)
    call handle_error_netcdf("NF90_CLOSE: "//trim(filename),status)

  END SUBROUTINE write_state_vector_all

  !-----------------------------------------------------------------------
  ! Read ensemble data and distribute to processes
  !-----------------------------------------------------------------------

  SUBROUTINE read_ens_mpi(islot,fcst3d,fcst2d)

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

    INTEGER :: irank_read_2d,irank_read_3d(nv3d)

    REAL(r_size),ALLOCATABLE :: v2dg(:,:,:)
    REAL(r_size),ALLOCATABLE :: v3dg(:,:,:,:)
    REAL(r_size),ALLOCATABLE :: hxf_k(:)   !Work: H(xf) of member imem

    CHARACTER(12) :: fcstfile="faTTNNNNN.nc"    
    CHARACTER(10) :: hxfile
    
    !---MPI
    INTEGER :: tag,ierr
    INTEGER,ALLOCATABLE :: req_send(:),req_recv(:)

    !---IN
    INTEGER,INTENT(IN) :: islot

    !---INOUT
    REAL(r_size),INTENT(INOUT) :: fcst3d(nij1,nlev,nmem,nv3d),fcst2d(nij1,nmem,nv2d)

    !---ALLOCATE
    ALLOCATE(v2dg(nlon,nlat,nv2d))
    ALLOCATE(v3dg(nlon,nlat,nlev,nv3d))
    ALLOCATE(req_send(nv3d))
    ALLOCATE(req_recv(nv3d))

    !---nb (number of do loop block)
    nb=CEILING(REAL(nmem)/REAL(nprocs))

    !---Main loop    
    DO ib=1,nb

       !---Initialization
       v2dg(1:nlon,1:nlat,1:nv2d) = REAL(0.d0, r_size)
       v3dg(1:nlon,1:nlat,1:nlev,1:nv3d) = REAL(0.d0, r_size)

       req_send(1:nv3d)=MPI_REQUEST_NULL
       req_recv(1:nv3d)=MPI_REQUEST_NULL

       !---imem
       IF(ib == 1)THEN
          imem = MOD(myrank, nmem) + 1
       ELSE
          imem = myrank + (ib-1)*nprocs + 1 
       END IF

       IF(1 <= imem .and.  imem <= nmem)THEN

          !---imem_master
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
          DO ivd=1,nv3d
             IF(myrank == imem_master .and. irank_read_3d(ivd) /= imem_master)THEN
                tag=(imem-1)*nv3d+ivd
                CALL MPI_IRECV(v3dg(:,:,:,ivd), nlon*nlat*nlev, MPI_R_SIZE, irank_read_3d(ivd), tag, &
                     & MPI_COMM_WORLD, req_recv(ivd), ierr)
             END IF
          END DO

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
          END DO

          !MPI_WAIT
          DO ivd=1,nv3d
             IF(myrank == imem_master .and. req_recv(ivd) /= MPI_REQUEST_NULL)THEN
                CALL MPI_WAIT(req_recv(ivd),MPI_STATUS_IGNORE,ierr)
             END IF

             IF(req_send(ivd) /= MPI_REQUEST_NULL)THEN
                CALL MPI_WAIT(req_send(ivd),MPI_STATUS_IGNORE,ierr)
             END IF
          END DO

          !---Observation operator H[xf(imem)] --> hdxfT(imem,:)
          !   Trans_XtoY updates obs. of islot only (INTENT(INOUT)),
          !   so hdxfT(imem,:) is copied out and back.
          IF(myrank == imem_master)THEN             

             ALLOCATE(hxf_k(nobs))
             CALL get_member_dxT_flat(imem,hdxfT_flat,hxf_k)
             CALL Trans_XtoY(islot,nobs, &
                  & obselm,obsidt,obsidxt,obsidxu,obsidxv,obsidyt,obsidyu,obsidyv,obslon,obslat,obslev, &
                  & v3dg,v2dg,hxf_k)
             CALL put_member_dxT_flat(imem,hdxfT_flat,hxf_k)

             !Write individual forecast in obs space
             IF(islot == nslots)THEN
                hxfile=hx_filename
                WRITE(hxfile(3:7),'(I5.5)') imem
                CALL make_ncfile_hx(hxfile)
                CALL write_hx(hxfile,"f",hxf_k)
             END IF
                
             DEALLOCATE(hxf_k)

          END IF

       END IF ! imem = 1 ~ nmem

       !---fcst3d,2d: distribute all members of this block at once
       IF(islot == nbslot)THEN
          CALL scatter_grd_mpi_all(ib,v3dg,v2dg,fcst3d,fcst2d)
       END IF

    END DO !ib

    DEALLOCATE(v2dg,v3dg)
    DEALLOCATE(req_send,req_recv)

  END SUBROUTINE read_ens_mpi
  
  !-----------------------------------------------------------------------
  ! Write ensemble data after collecting data from processes
  !-----------------------------------------------------------------------

  SUBROUTINE write_ens_mpi(v3da,v2da,v3df_m,v3df_s,v2df_m,v2df_s)

    !$ USE OMP_LIB
    USE common
    USE common_monit
    USE common_mpi
    USE common_setting
    USE MPI
    IMPLICIT NONE

    !---Common
    INTEGER,PARAMETER :: ntask=nmem+4
    ![1,nmem]: ensemble analysis
    ![nmem+1]: anal_mean, [nmem+2]: anal_sprd
    ![nmem+3]: fcst_mean, [nmem+4]: fcst_sprd

    INTEGER i,n,k
    INTEGER imem
    INTEGER nrecv
    INTEGER status(MPI_STATUS_SIZE),ierr

    REAL(r_size),ALLOCATABLE :: v3da_m(:,:,:),v2da_m(:,:)       !Analysis ensemble mean
    REAL(r_size),ALLOCATABLE :: v3da_s(:,:,:),v2da_s(:,:)       !                  spread

    REAL(r_sngl),ALLOCATABLE :: v3dg(:,:,:,:),v2dg(:,:,:) !Global domain

    REAL(r_size),ALLOCATABLE :: hxa_k(:)                  !Work array: H(x) of one task (all ranks)
    
    CHARACTER(12) :: filename='file00000.nc'
    CHARACTER(10) hxfile !Individual file in obs. space (hxNNNNN.nc)

    !---IN    
    REAL(r_size),INTENT(IN) :: v3da(nij1,nlev,nmem,nv3d),v2da(nij1,nmem,nv2d) !Analysis ensemble (local)
    REAL(r_size),INTENT(IN) :: v3df_m(nij1,nlev,nv3d),v2df_m(nij1,nv2d) !Forecast ensemble mean(local)
    REAL(r_size),INTENT(IN) :: v3df_s(nij1,nlev,nv3d),v2df_s(nij1,nv2d) !Forecast ensemble spread(local)

    !---H(Xa)^T on root_out: dYa(iobs,imem) = hdxaT_flat(imem+(iobs-1)*nmem)
    IF(myrank == root_out)THEN
       ALLOCATE(hdxaT_flat(INT(nmem,8)*INT(nobs,8)))
       hdxaT_flat(:) = REAL(0.d0, r_size)
    END IF

    !---Initialization
    ALLOCATE(v3da_m(nij1,nlev,nv3d),v2da_m(nij1,nv2d))
    ALLOCATE(v3da_s(nij1,nlev,nv3d),v2da_s(nij1,nv2d))
    ALLOCATE(v3dg(nlon,nlat,nlev,nv3d),v2dg(nlon,nlat,nv2d))
    ALLOCATE(hxa_k(nobs))
    hxa_k(:)=REAL(0.d0,r_size)
    obsidt(:)=nbslot

    !---Ensemble mean/sprd
    call ensemble_mesp(v3da,v2da,v3da_m,v3da_s,v2da_m,v2da_s)

    !---Loop index
    n = CEILING(REAL(ntask)/REAL(nprocs))

    !---Analysis ensemble, and forecast/analysis ensemble mean and spread (Output) --> Obs. space
    DO i=1,n

       !---GATHER all tasks of this block to their owners (rank = imem-1-(i-1)*nprocs)
       CALL gather_grd_mpi_all(i,ntask,v3da,v2da,v3da_m,v2da_m,v3da_s,v2da_s, &
            &                  v3df_m,v2df_m,v3df_s,v2df_s,v3dg,v2dg)

       !---Post process for quasi-global
       IF(lqglobal)THEN
          v2dg(1,:,:)    = v2dg(nlon-1,:,:)
          v2dg(nlon,:,:) = v2dg(2,:,:)
          v3dg(1,:,:,:)    = v3dg(nlon-1,:,:,:)
          v3dg(nlon,:,:,:) = v3dg(2,:,:,:)
       END IF

       imem = myrank+1 + (i-1)*nprocs

       !---Observation space: H(x) of own task & MPI_SEND to root_out
       IF(1 <= imem .and. imem <= nmem)THEN

          CALL Trans_XtoY(nbslot,nobs, &
               & obselm,obsidt,obsidxt,obsidxu,obsidxv,obsidyt,obsidyu,obsidyv,obslon,obslat,obslev, &
               & v3dg,v2dg,hxa_k)

          IF(myrank == root_out)THEN
             !Own task: store directly (no self-send)
             CALL put_member_dxT_flat(imem,hdxaT_flat,hxa_k)
          ELSE
             CALL MPI_SEND(hxa_k,nobs,MPI_R_SIZE,root_out,imem,MPI_COMM_WORLD,ierr)
          END IF

          !H(xa) of each member --> individual file (created in read_ens_mpi)
          hxfile=hx_filename
          WRITE(hxfile(3:7),'(I5.5)') imem
          CALL write_hx(hxfile,"a",hxa_k)
          
       END IF

       !---root_out receives all tasks of this iteration before the next collective
       IF(myrank == root_out)THEN

          nrecv = MAX( MIN(nmem,i*nprocs) - (i-1)*nprocs, 0 ) !#members in this iteration          
          IF(1 <= imem .and. imem <= nmem) nrecv = nrecv - 1  !Except own task

          DO k=1,nrecv
             CALL MPI_RECV(hxa_k,nobs,MPI_R_SIZE,MPI_ANY_SOURCE,MPI_ANY_TAG, &
                  & MPI_COMM_WORLD,status,ierr)
             CALL put_member_dxT_flat(status(MPI_TAG),hdxaT_flat,hxa_k)
          END DO

       END IF

       !---Obs. space output (root_out), overlapped with state output (task owners)
       !   hdxaT is complete after the last block (i == n).
       !   root_out owns no task if nprocs > ntask, so it is idle in write_state_vector_all.
       IF(i == n .and. myrank == root_out)THEN

          ALLOCATE(hxamean(nobs),hxasprd(nobs))
          hxamean(:)    = REAL(0.d0, r_size)
          hxasprd(:)    = REAL(0.d0, r_size)
          
          !---H(Xa) => dYa=H(dXa), Hxamean, Hxasprd
          CALL ensemble_mesp_obs(hdxaT_flat,hxamean,hxasprd)
          
          !---Monitor
          CALL monit_mean("fcst",hxfmean)
          CALL monit_sprd("fcst",hxfsprd)
          CALL monit_mean("anal",hxamean)
          CALL monit_sprd("anal",hxasprd)
          
          !---Assimilated observation
          CALL make_ncfile_obs
          CALL write_obs

          !---Ensemble mean and spread in obs. space
          CALL make_ncfile_hx(hx_mean_filename)
          CALL write_hx(hx_mean_filename,"f",hxfmean)
          CALL write_hx(hx_mean_filename,"a",hxamean)

          CALL make_ncfile_hx(hx_sprd_filename)
          CALL write_hx(hx_sprd_filename,"f",hxfsprd)
          CALL write_hx(hx_sprd_filename,"a",hxasprd)
                    
          DEALLOCATE( hdxaT_flat )
          DEALLOCATE( hxamean, hxasprd )

       END IF
       
       !---Filename
       IF(1 <= imem .and. imem <= nmem)THEN
          WRITE(filename(1:9),'(A2,I2.2,I5.5)') "fa",nbslot,imem
       ELSE IF(imem == nmem+1)THEN
          WRITE(filename(1:9),'(A9)') "anal_mean"
       ELSE IF(imem == nmem+2)THEN
          WRITE(filename(1:9),'(A9)') "anal_sprd"
       ELSE IF(imem == nmem+3)THEN
          WRITE(filename(1:9),'(A9)') "fcst_mean"
       ELSE IF(imem == nmem+4)THEN
          WRITE(filename(1:9),'(A9)') "fcst_sprd"
       END IF

       !---Write
       IF(1 <= imem .and. imem <= ntask)THEN
          CALL write_state_vector_all(filename,v3dg,v2dg)
       END IF

    END DO !i

    !---Deallocate
    DEALLOCATE(v3da_m,v2da_m)
    DEALLOCATE(v3da_s,v2da_s)
    DEALLOCATE(v3dg,v2dg)
    
    DEALLOCATE( obsidt )
    DEALLOCATE( obsidxt, obsidyt )
    DEALLOCATE( obsidxu, obsidyu )
    DEALLOCATE( obsidxv, obsidyv )
    DEALLOCATE( obselm, obsins, obsqc )
    DEALLOCATE( obslon, obslat, obslev )
    DEALLOCATE( obsdat, obserr )
    DEALLOCATE( dob )
    DEALLOCATE( hxfmean, hxfsprd )
    DEALLOCATE( hxa_k )
    
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

    !---Common
    INTEGER status,access
    INTEGER ncid,dimid
    INTEGER ierr

    !---IN
    CHARACTER(8),INTENT(IN) :: filename

    !---OUT
    INTEGER,INTENT(OUT) :: no

    status=0
    IF(myrank == root)THEN
       status=access(filename," ")
    END IF

    CALL MPI_BCAST(status,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

    IF(status /= 0)THEN
       WRITE(6,'(A)') " *** Error: Not Found "//filename
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
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

    INTEGER ,ALLOCATABLE :: itmp(:)
    REAL(r_sngl),ALLOCATABLE :: tmp(:)

    !IN
    INTEGER,INTENT(IN) :: no
    CHARACTER(8),INTENT(IN) :: filename

    !OUT
    INTEGER,INTENT(OUT) :: ele(no) ! element ID
    INTEGER,INTENT(OUT) :: ins(no) ! instrument ID

    REAL(r_size),INTENT(OUT) :: lon(no),lat(no),lev(no) ! Grid information [degree E, degree N, m]
    REAL(r_size),INTENT(OUT) :: obs(no) ! Obs.
    REAL(r_size),INTENT(OUT) :: err(no) ! Obs. error

    IF(myrank == root)THEN

       !Check file
       status=access(filename," ")

       IF(status /= 0)THEN
          WRITE(6,'(A)') " *** Error: Not Found "//filename
          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
       END IF

       !open
       status=NF90_OPEN(trim(filename),NF90_NOWRITE,ncid)
       call handle_error_netcdf("NF90_OPEN: "//trim(filename),status)

       ALLOCATE(itmp(no))

       !element
       call read_netcdf_var_int_1d(ncid,no,"ele",itmp)    
       ele(:)=itmp(:)

       !instrument
       call read_netcdf_var_int_1d(ncid,no,"ins",itmp)    
       ins(:)=itmp(:)

       DEALLOCATE(itmp)

       ALLOCATE(tmp(no))

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

       DEALLOCATE(tmp)

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
  ! NetCDF file in obs space |
  !--------------------------------------------------------------

  SUBROUTINE make_ncfile_obs

    USE common_setting, only: nobs, obs_assim_filename
    USE NETCDF
    IMPLICIT NONE

    !---Parameter
    INTEGER,PARAMETER :: ndim=1

    !---Common
    INTEGER status
    INTEGER ncid,dimid,varid
    INTEGER dim(ndim)

    status=NF90_CREATE(trim(obs_assim_filename),nf90_netcdf4,ncid)
    CALL handle_error_netcdf("NF90_CREATE: "//trim(obs_assim_filename),status)

    status=NF90_PUT_ATT(ncid,NF90_GLOBAL,"title","Assimilated Observation")
    CALL handle_error_netcdf("NF90_PUT_ATT: "//trim(obs_assim_filename),status)

    status=NF90_PUT_ATT(ncid,NF90_GLOBAL,"description","Assimilated Observation")
    CALL handle_error_netcdf("NF90_PUT_ATT: "//trim(obs_assim_filename),status)
    
    !#Obs.
    status=NF90_DEF_DIM(ncid,"nobs",nobs,dimid)
    CALL handle_error_netcdf("NF90_DEF_DIM: "//trim(obs_assim_filename),status)
    dim(1)=dimid

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

    status=NF90_ENDDEF(ncid)
    CALL handle_error_netcdf("NF90_ENDDEF: "//trim(obs_assim_filename),status)

    status=NF90_CLOSE(ncid)
    CALL handle_error_netcdf("NF90_CLOSE: "//trim(obs_assim_filename),status)

  END SUBROUTINE make_ncfile_obs

  !-------------------------------------

  SUBROUTINE make_ncfile_hx(filename)

    USE common_setting, only: nobs
    USE NETCDF
    IMPLICIT NONE

    !---Parameter
    INTEGER,PARAMETER :: ndim=1

    !---Common
    INTEGER status
    INTEGER ncid,dimid,varid
    INTEGER dim(ndim)

    !---IN
    CHARACTER(*),INTENT(IN) :: filename

    status=NF90_CREATE(trim(filename),nf90_netcdf4,ncid)
    CALL handle_error_netcdf("NF90_CREATE: "//trim(filename),status)

    status=NF90_PUT_ATT(ncid,NF90_GLOBAL,"title","State variable in observation space")
    CALL handle_error_netcdf("NF90_PUT_ATT: "//trim(filename),status)

    status=NF90_PUT_ATT(ncid,NF90_GLOBAL,"description","State variable in observation space")
    CALL handle_error_netcdf("NF90_PUT_ATT: "//trim(filename),status)
    
    !#Obs.
    status=NF90_DEF_DIM(ncid,"nobs",nobs,dimid)
    CALL handle_error_netcdf("NF90_DEF_DIM: "//trim(filename),status)
    dim(1)=dimid

    !1D
    call define_var_netcdf(ncid,1,dim(1),varid,"real", &
         & "hxf","forecast in obs. space","H: meter, U and V: m/s, T: degree C, S: -")

    call define_var_netcdf(ncid,1,dim(1),varid,"real", &
         & "hxa","analysis in obs. space","H: meter, U and V: m/s, T: degree C, S: -")

    status=NF90_ENDDEF(ncid)
    CALL handle_error_netcdf("NF90_ENDDEF: "//trim(filename),status)

    status=NF90_CLOSE(ncid)
    CALL handle_error_netcdf("NF90_CLOSE: "//trim(filename),status)

  END SUBROUTINE make_ncfile_hx
  
  !-------------------------------------

  SUBROUTINE write_hx(filename,fa,hx)

    USE common_setting, only: r_sngl, r_size, nobs
    USE MPI
    USE NETCDF
    IMPLICIT NONE

    !---Common
    INTEGER status,ncid
    INTEGER ierr
    
    !---IN
    REAL(r_size),INTENT(IN) :: hx(nobs)
    
    CHARACTER(*),INTENT(IN) :: filename
    CHARACTER(1),INTENT(IN) :: fa
    
    status=NF90_OPEN(trim(filename),NF90_WRITE,ncid)
    call handle_error_netcdf("NF90_OPEN: "//trim(filename),status)

    IF(fa == "f")THEN
       CALL write_netcdf_var_sngl_1d(ncid,nobs,"hxf",REAL(hx(:),r_sngl))
    ELSE IF(fa == "a")THEN
       CALL write_netcdf_var_sngl_1d(ncid,nobs,"hxa",REAL(hx(:),r_sngl))
    ELSE
       write(6,'(A)') "Error: fa is incorrect ==> "//trim(fa)
       CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
    END IF
    
    status=NF90_CLOSE(ncid)
    CALL handle_error_netcdf("NF90_CLOSE: "//trim(filename),status)

  END SUBROUTINE write_hx

  !-------------------------------------
  ! Write assimilated observations (obs_assim_filename)
  !-------------------------------------

  SUBROUTINE write_obs

    USE common_setting
    USE NETCDF
    IMPLICIT NONE

    !---Common
    INTEGER status,ncid

    status=NF90_OPEN(trim(obs_assim_filename),NF90_WRITE,ncid)
    call handle_error_netcdf("NF90_OPEN: "//trim(obs_assim_filename),status)

    CALL write_netcdf_var_int_1d(ncid,nobs,"ele",obselm(:))
    CALL write_netcdf_var_int_1d(ncid,nobs,"ins",obsins(:))
    CALL write_netcdf_var_sngl_1d(ncid,nobs,"lon",REAL(obslon(:),r_sngl))
    CALL write_netcdf_var_sngl_1d(ncid,nobs,"lat",REAL(obslat(:),r_sngl))
    CALL write_netcdf_var_sngl_1d(ncid,nobs,"lev",REAL(obslev(:),r_sngl))
    CALL write_netcdf_var_sngl_1d(ncid,nobs,"obs",REAL(obsdat(:),r_sngl))
    CALL write_netcdf_var_sngl_1d(ncid,nobs,"err",REAL(obserr(:),r_sngl))

    status=NF90_CLOSE(ncid)
    call handle_error_netcdf("NF90_CLOSE: "//trim(obs_assim_filename),status)

  END SUBROUTINE write_obs
    
END MODULE common_io
