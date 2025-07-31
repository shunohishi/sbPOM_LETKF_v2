MODULE common_obs
  
  !=======================================================================
  !
  ! [PURPOSE:] Observational procedures
  !
  ! [HISTORY:]
  !   01/23/2009 Takemasa MIYOSHI  created
  !   01/26/2011 Yasumasa MIYAZAWA  modified for POM (check 'pom' or 'POM')
  !   07/31/2025 Shun OHISHI modofied
  !
  !=======================================================================

  IMPLICIT NONE
  PUBLIC

CONTAINS

  !-----------------------------------------------------------------------
  ! Initialize
  !-----------------------------------------------------------------------

  SUBROUTINE set_obs(fcst3d,fcst2d)

    !$USE OMP_LIB
    USE MPI
    USE common_setting
    USE common
    USE common_mpi
    USE common_io
    IMPLICIT NONE

    !---Parameter
    REAL(r_size),PARAMETER :: gross_error=10.0d0

    !---Common
    INTEGER :: n,nn,ni
    INTEGER :: i,j
    INTEGER :: islot,ibv,iprocs,ierr

    INTEGER :: nobslots(nslots)
    INTEGER :: nj(0:nlat-1), njs(1:nlat-1)

    INTEGER :: issh=0,isst=0,isss=0,issu=0,issv=0
    INTEGER :: it=0,is=0,iu=0,iv=0
    INTEGER :: numdelqc1=0, numdelqc2=0

    INTEGER,ALLOCATABLE :: iwork2d(:,:)

    REAL(r_size) :: v2dg(nlon,nlat,nv2d)    
    REAL(r_size) :: v3dg(nlon,nlat,nlev,nv3d)

    REAL(r_size) :: obserrt
    REAL(r_size) :: dob2_sprd2 !dob(obsdep)**2 - sprd**2
    REAL(r_size),ALLOCATABLE :: sprd_hdxf(:)
    REAL(r_size),ALLOCATABLE :: work2d(:,:)

    CHARACTER(8) ::  obsfile="obsTT.nc"
    CHARACTER(12) :: fcstfile="faTTNNNNN.nc"    

    !---tmp
    INTEGER,ALLOCATABLE :: tmpqc0(:,:) !Ensemble QC
    INTEGER,ALLOCATABLE :: tmpqc(:)    !QC
    INTEGER,ALLOCATABLE :: tmpidx(:),tmpidy(:)
    INTEGER,ALLOCATABLE :: tmpelm(:),tmpins(:)

    REAL(r_size),ALLOCATABLE :: tmplon(:),tmplat(:),tmplev(:)
    REAL(r_size),ALLOCATABLE :: tmpdat(:)
    REAL(r_size),ALLOCATABLE :: tmperr(:)
    REAL(r_size),ALLOCATABLE :: tmpdep(:)    !Innovation
    REAL(r_size),ALLOCATABLE :: tmphdxf(:,:) !dY=HdX Ensemble perturvation in obs. space

    !---tmp2    
    INTEGER,ALLOCATABLE :: tmp2idx(:),tmp2idy(:)
    INTEGER,ALLOCATABLE :: tmp2elm(:),tmp2ins(:)

    REAL(r_size),ALLOCATABLE :: tmp2lon(:),tmp2lat(:),tmp2lev(:)
    REAL(r_size),ALLOCATABLE :: tmp2dat(:)
    REAL(r_size),ALLOCATABLE :: tmp2err(:)
    REAL(r_size),ALLOCATABLE :: tmp2dep(:)
    REAL(r_size),ALLOCATABLE :: tmp2hdxf(:,:)

    !---OUT
    REAL(r_size),INTENT(OUT) :: fcst3d(nij1,nlev,nbv,nv3d),fcst2d(nij1,nbv,nv2d)

    WRITE(file_unit,'(A)') "Hello from set_obs"

    !---Horizontal/Vertical cutoff scale
    dist_zero = sigma_obs * SQRT(10.0d0/3.0d0) * 2.0d0
    dist_zerov = sigma_obsv * SQRT(10.0d0/3.0d0) * 2.0d0

    !---Get Nobs
    DO islot=1,nslots
       WRITE(obsfile(4:5),'(I2.2)') islot
       CALL get_nobs(obsfile,nobslots(islot))
    END DO !islot
    nobs = SUM(nobslots)
    WRITE(file_unit,'(A,I10)') "TOTAL NUMBER of OBSERVATIONS:",nobs

    !---ALLOCATE GLOBAL OBSERVATION
    ALLOCATE( tmpelm(nobs), tmpins(nobs) )
    ALLOCATE( tmplon(nobs), tmplat(nobs), tmplev(nobs) )
    ALLOCATE( tmpdat(nobs), tmperr(nobs) )
    ALLOCATE( tmpidx(nobs), tmpidy(nobs) )
    ALLOCATE( tmpdep(nobs) )
    ALLOCATE( tmphdxf(nobs,nbv) )
    ALLOCATE( tmpqc0(nobs,nbv), tmpqc(nobs) )

    !---Initialization
    tmphdxf(:,:) = 0.0d0
    tmpqc0(:,:) = 0

    !---LOOP of timeslots
    nn=0
    DO islot=1,nslots

       IF(nobslots(islot) == 0) CYCLE

       WRITE(obsfile(4:5),'(I2.2)') islot

       !---Read obs data
       CALL read_obs(obsfile,nobslots(islot), &
            & tmpelm(nn+1:nn+nobslots(islot)), tmpins(nn+1:nn+nobslots(islot)), &
            & tmplon(nn+1:nn+nobslots(islot)), tmplat(nn+1:nn+nobslots(islot)), tmplev(nn+1:nn+nobslots(islot)), &
            & tmpdat(nn+1:nn+nobslots(islot)), tmperr(nn+1:nn+nobslots(islot)))

       !---Get obs. ID
       CALL obs_id(nobslots(islot), &
            & tmplon(nn+1:nn+nobslots(islot)),tmplat(nn+1:nn+nobslots(islot)),&
            & tmpidx(nn+1:nn+nobslots(islot)),tmpidy(nn+1:nn+nobslots(islot)))

       !---Calculate tmphdxf
       ni=CEILING(REAL(nbv)/REAL(nprocs))

       DO i=1,ni

          ibv = myrank+1 + (i-1)*nprocs !ith ensemble member read by my_rank

          IF(ibv <= nbv)THEN
             !---Read forecast at islot --> fcst3d,2d
             WRITE(fcstfile(3:9),'(I2.2,I5.5)') islot,ibv
             CALL read_state_vector(fcstfile,v3dg,v2dg)

             !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n)
             DO n=1,nobslots(islot)

                !---Observation operator H(xf(i))
                CALL Trans_XtoY(&
                     & tmpelm(nn+n),tmpidx(nn+n),tmpidy(nn+n), &
                     & tmplon(nn+n),tmplat(nn+n),tmplev(nn+n), &
                     & v3dg,v2dg,tmphdxf(nn+n,ibv))

                !---tmpqc0
                IF(tmphdxf(nn+n,ibv) == undef)THEN
                   tmpqc0(nn+n,ibv) = 0                
                ELSE
                   tmpqc0(nn+n,ibv) = 1
                ENDIF

             END DO !n
             !$OMP END PARALLEL DO

          END IF

          DO iprocs=0,nprocs-1
             ibv = iprocs+1 + (i-1)*nprocs
             IF(islot == 1 .and. ibv <= nbv)THEN          
                CALL scatter_grd_mpi(iprocs,REAL(v3dg,r_sngl),REAL(v2dg,r_sngl), &
                     & fcst3d(:,:,ibv,:),fcst2d(:,ibv,:))
             END IF
          END DO !iprocs

       END DO !i

       nn = nn + nobslots(islot)

    END DO !islot

    !---Share tmphdxf (Hxf(i)) for all processors
    ALLOCATE(work2d(nobs,nbv))
    work2d(:,:) = tmphdxf(:,:)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    IF(r_size == kind(0.0d0))THEN
       CALL MPI_ALLREDUCE(work2d,tmphdxf,nobs*nbv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    ELSE IF(r_size == kind(0.0e0))THEN
       CALL MPI_ALLREDUCE(work2d,tmphdxf,nobs*nbv,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
    END IF

    DEALLOCATE(work2d)

    !---Share tmpqc0 for all processors
    ALLOCATE(iwork2d(nobs,nbv))
    iwork2d(:,:) = tmpqc0(:,:)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(iwork2d,tmpqc0,nobs*nbv,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    DEALLOCATE(iwork2d)

    !---Hxfmean, dYf, Innovation
    !$OMP PARALLEL
    !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(n,ibv)
    DO n=1,nobs

       tmpqc(n) = MINVAL(tmpqc0(n,:))

       IF(tmpqc(n) == 1)THEN

          !tmpdep: Hxfmean
          tmpdep(n) = tmphdxf(n,1)
          DO ibv=2,nbv
             tmpdep(n) = tmpdep(n) + tmphdxf(n,ibv)
          END DO
          tmpdep(n) = tmpdep(n) / REAL(nbv,r_size)

          !tmphdxf: dYf=HdXf
          DO ibv=1,nbv
             tmphdxf(n,ibv) = tmphdxf(n,ibv) - tmpdep(n)
          END DO

          !tmpdep: Innovation y-Hxfmean
          tmpdep(n) = tmpdat(n) - tmpdep(n)

       ELSE

          tmpdep(n)=undef
          tmphdxf(n,:)=undef
          tmpdep(n)=undef

       END IF

    END DO !n
    !$OMP END DO
    !$OMP END PARALLEL

    !#OBS deleted by Obs. projection
    numdelqc1=nobs-SUM(tmpqc)

    !---Gross error QC
    IF(LGE_IO .and. .not. AOEI)THEN

       !$OMP PARALLEL
       !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(n)
       DO n=1,nobs

          IF(tmpelm(n) == id_z_obs)THEN !SSH
             CALL check_gross_error_obs(hqc_min,hqc_max,tmpdep(n),tmpqc(n))
          ELSE IF(tmpelm(n) == id_t_obs .and. tmplev(n) == 0.d0)THEN !SST
             CALL check_gross_error_obs(sstqc_min,sstqc_max,tmpdep(n),tmpqc(n))
          ELSE IF(tmpelm(n) == id_t_obs)THEN !T
             CALL check_gross_error_obs(tqc_min,tqc_max,tmpdep(n),tmpqc(n))
          ELSE IF(tmpelm(n) == id_s_obs .and. tmplev(n) == 0.d0)THEN !SSS
             CALL check_gross_error_obs(sssqc_min,sssqc_max,tmpdep(n),tmpqc(n))
          ELSE IF(tmpelm(n) == id_s_obs)THEN !S
             CALL check_gross_error_obs(sqc_min,sqc_max,tmpdep(n),tmpqc(n))
          ELSE IF(tmpelm(n) == id_u_obs .and. tmplev(n) == 0.d0)THEN !SSU
             CALL check_gross_error_obs(ssuqc_min,ssuqc_max,tmpdep(n),tmpqc(n))
          ELSE IF(tmpelm(n) == id_u_obs)THEN !U
             CALL check_gross_error_obs(uqc_min,uqc_max,tmpdep(n),tmpqc(n))
          ELSE IF(tmpelm(n) == id_v_obs .and. tmplev(n) == 0.d0)THEN !SSV
             CALL check_gross_error_obs(ssvqc_min,ssvqc_max,tmpdep(n),tmpqc(n))
          ELSE IF(tmpelm(n) == id_v_obs)THEN !V
             CALL check_gross_error_obs(vqc_min,vqc_max,tmpdep(n),tmpqc(n))
          ELSE IF(ABS(tmpdep(n)) > gross_error*tmperr(n))THEN
             tmpqc(n) = 0
          END IF

       END DO !n
       !$OMP END DO
       !$OMP END PARALLEL

    END IF
    DEALLOCATE(tmpqc0)

    numdelqc2=nobs-numdelqc1-SUM(tmpqc)

    WRITE(file_unit,'(A,I10)') "#OBS. to be ASSIMILATED: ",SUM(tmpqc)
    WRITE(file_unit,'(A,I10)')   "#OBS. DELETED by Obs. Projection: ",numdelqc1
    WRITE(file_unit,'(A,I10,A)') "#OBS. DELETED by Large Deviations: ",numdelqc2

    CALL monit_dep(nobs,tmpelm,tmpdep,tmpqc)

    !--- Temporal observation localization
    nn = 0
    DO islot=1,nslots
       tmperr(nn+1:nn+nobslots(islot)) = tmperr(nn+1:nn+nobslots(islot)) &
            & * exp(0.25d0 * (REAL(islot-nbslot,r_size) / sigma_obst)**2)
       nn = nn + nobslots(islot)
    END DO

    !--- Remove undef
    nn = 0
    DO n=1,nobs
       IF(tmpqc(n) == 1)THEN
          nn = nn+1
          tmpelm(nn) = tmpelm(n)
          tmpins(nn) = tmpins(n)
          tmplon(nn) = tmplon(n)
          tmplat(nn) = tmplat(n)
          tmplev(nn) = tmplev(n)
          tmpdat(nn) = tmpdat(n)
          tmperr(nn) = tmperr(n)
          tmpidx(nn) = tmpidx(n)
          tmpidy(nn) = tmpidy(n)
          tmpdep(nn) = tmpdep(n)
          tmphdxf(nn,:) = tmphdxf(n,:)
          tmpqc(nn) = tmpqc(n)
       END IF
    END DO !n

    !--- Total #OBS without undef
    nobs = nn
    WRITE(file_unit,'(A,I5.5,A,I10)') "MYRANK: ", myrank, " #OBS: ", nobs

    !--- SORT
    ALLOCATE( tmp2elm(nobs), tmp2ins(nobs) )
    ALLOCATE( tmp2lon(nobs), tmp2lat(nobs), tmp2lev(nobs))
    ALLOCATE( tmp2dat(nobs) )
    ALLOCATE( tmp2err(nobs) )
    ALLOCATE( tmp2idx(nobs), tmp2idy(nobs) )
    ALLOCATE( tmp2dep(nobs) )
    ALLOCATE( tmp2hdxf(nobs,nbv) )

    ALLOCATE( obselm(nobs), obsins(nobs) )
    ALLOCATE( obslon(nobs), obslat(nobs), obslev(nobs) )
    ALLOCATE( obsdat(nobs) )
    ALLOCATE( obserr(nobs) )
    ALLOCATE( obsidx(nobs), obsidy(nobs) )
    ALLOCATE( obsdep(nobs) )
    ALLOCATE( obshdxf(nobs,nbv) )

    ALLOCATE( sprd_hdxf(nobs) )

    nobsgrd(:,:) = 0
    nj(:) = 0

    !nj: #OBS at a latitude band
    !$OMP PARALLEL PRIVATE(i,j,n,nn)
    !$OMP DO SCHEDULE(DYNAMIC)
    DO j=1,nlat-1
       DO n=1,nobs
          IF(tmpidy(n) < j .OR. j+1 <= tmpidy(n)) CYCLE
          nj(j) = nj(j) + 1
       END DO !n
    END DO !j
    !$OMP END DO

    !njs: #OBS for 0 ~ j-1 at j --> Start point
    !$OMP DO SCHEDULE(DYNAMIC)
    DO j=1,nlat-1
       njs(j) = SUM(nj(0:j-1))
    END DO !j
    !$OMP END DO

    !---Sort for latitude direction --> tmp2
    !$OMP DO SCHEDULE(DYNAMIC)
    DO j=1,nlat-1
       nn = 0
       DO n=1,nobs
          IF(tmpidy(n) < j .OR. j+1 <= tmpidy(n)) CYCLE
          nn = nn + 1
          tmp2elm(njs(j)+nn) = tmpelm(n)
          tmp2ins(njs(j)+nn) = tmpins(n)
          tmp2lon(njs(j)+nn) = tmplon(n)
          tmp2lat(njs(j)+nn) = tmplat(n)
          tmp2lev(njs(j)+nn) = tmplev(n)
          tmp2dat(njs(j)+nn) = tmpdat(n)
          tmp2err(njs(j)+nn) = tmperr(n)
          tmp2idx(njs(j)+nn) = tmpidx(n)
          tmp2idy(njs(j)+nn) = tmpidy(n)
          tmp2dep(njs(j)+nn) = tmpdep(n)
          tmp2hdxf(njs(j)+nn,:) = tmphdxf(n,:)
       END DO !n
    END DO !j
    !$OMP END DO

    !Sort for longitude direction --> obs
    !$OMP DO SCHEDULE(DYNAMIC)
    DO j=1,nlat-1

       IF(nj(j) == 0)THEN
          nobsgrd(:,j) = njs(j)
          CYCLE
       END IF

       nn = 0
       DO i=1,nlon
          DO n=njs(j)+1,njs(j)+nj(j)
             IF(tmp2idx(n) < i .OR. i+1 <= tmp2idx(n)) CYCLE
             nn = nn + 1
             obselm(njs(j)+nn) = tmp2elm(n)
             obsins(njs(j)+nn) = tmp2ins(n)
             obslon(njs(j)+nn) = tmp2lon(n)
             obslat(njs(j)+nn) = tmp2lat(n)
             obslev(njs(j)+nn) = tmp2lev(n)
             obsdat(njs(j)+nn) = tmp2dat(n)
             obserr(njs(j)+nn) = tmp2err(n)
             obsidx(njs(j)+nn) = tmp2idx(n)
             obsidy(njs(j)+nn) = tmp2idy(n)
             obsdep(njs(j)+nn) = tmp2dep(n)
             obshdxf(njs(j)+nn,:) = tmp2hdxf(n,:)
          END DO !n
          nobsgrd(i,j) = njs(j) + nn
       END DO !i

       !Error check
       IF(nn /= nj(j)) THEN
          !$OMP CRITICAL
          WRITE(file_unit,'(A,2I)') 'OBS DATA SORT ERROR: ',nn,nj(j)
          WRITE(file_unit,'(F6.2,A,F6.2)') j,'< J <',j+1
          WRITE(file_unit,'(F6.2,A,F6.2)') &
               & MINVAL(tmp2idy(njs(j)+1:njs(j)+nj(j))),'< OBSJ <',MAXVAL(tmp2idy(njs(j)+1:njs(j)+nj(j)))
          !$OMP END CRITICAL
       END IF

    END DO !j
    !$OMP END DO
    !$OMP END PARALLEL

    !---AOEI (Minamide and Zhang 2017; Monthly Weather Review) S.Ohishi 2019.08
    IF(AOEI)THEN

       !Forecast ensemble spread in obs. space
       !$OMP PARALLEL
       !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(n)
       DO n=1,nobs
          sprd_hdxf(n) = obshdxf(n,1)**2
          DO i=2,nbv
             sprd_hdxf(n) = sprd_hdxf(n) + obshdxf(n,i)**2
          END DO !i
          sprd_hdxf(n) = SQRT(sprd_hdxf(n) / REAL(nbv-1,r_size))
       END DO !n
       !$OMP END DO

       !Update obserr
       !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(n,obserrt,dob2_sprd2)
       DO n=1,nobs
          obserrt = obserr(n) !Prescribed obs. err
          dob2_sprd2 = obsdep(n)*obsdep(n)-sprd_hdxf(n)*sprd_hdxf(n) !dob**2 - sprd**2
          obserr(n) = sqrt( max(obserrt*obserrt,dob2_sprd2) )
       END DO !n
       !$OMP END DO
       !$OMP END PARALLEL

    END IF

    !---Deallocate
    DEALLOCATE( tmpelm, tmpins )
    DEALLOCATE( tmplon, tmplat, tmplev )
    DEALLOCATE( tmpdat )
    DEALLOCATE( tmperr )
    DEALLOCATE( tmpidx, tmpidy )
    DEALLOCATE( tmpdep )
    DEALLOCATE( tmphdxf )
    DEALLOCATE( tmpqc )

    DEALLOCATE( tmp2elm, tmp2ins )
    DEALLOCATE( tmp2lon, tmp2lat, tmp2lev )
    DEALLOCATE( tmp2dat )
    DEALLOCATE( tmp2err )
    DEALLOCATE( tmp2idx, tmp2idy )
    DEALLOCATE( tmp2dep )
    DEALLOCATE( tmp2hdxf )

    DEALLOCATE( sprd_hdxf )

  END SUBROUTINE set_obs

  !-----------------------------------------------------------------------
  ! Check gross error observation
  !-----------------------------------------------------------------------

  SUBROUTINE check_gross_error_obs(qc_min,qc_max,inv,qc)

    USE common_setting, only: r_size
    IMPLICIT NONE

    !IN
    REAL(r_size),INTENT(IN) :: qc_min,qc_max !Min. and Max. Error gross value
    REAL(r_size),INTENT(IN) :: inv           !innovation

    !INOUT
    INTEGER,INTENT(INOUT) :: qc

    IF(inv < qc_min .or. qc_max < inv)THEN
       qc=0
    END IF

  END SUBROUTINE check_gross_error_obs

END MODULE common_obs
