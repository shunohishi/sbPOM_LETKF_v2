MODULE common_obs
  
  !=======================================================================
  !
  ! [PURPOSE:] Get observation and dYf information
  !
  ! [HISTORY:]
  !   01/23/2009 Takemasa MIYOSHI  created
  !   01/26/2011 Yasumasa MIYAZAWA  modified for POM (check 'pom' or 'POM')
  !   07/31/2025 Shun OHISHI modified
  !   06/04/2026 Shun OHISHI updated (Memory and I/O issue)
  !
  !=======================================================================

  IMPLICIT NONE
  PUBLIC

CONTAINS

  !-----------------------------------------------------------------------
  ! Set information in observation space
  !-----------------------------------------------------------------------

  SUBROUTINE set_obs(fcst3d,fcst2d)

    !$ USE OMP_LIB
    USE MPI
    USE common_setting
    USE common
    USE common_mpi
    USE common_io
    IMPLICIT NONE

    !---Parameter
    REAL(r_size),PARAMETER :: gross_error=10.0d0

    !---Common
    INTEGER i,j
    INTEGER islot
    INTEGER iobs,no,no1,no2
    INTEGER imem
        
    INTEGER nobslots(nslots)
    INTEGER nj(0:nlat-1),njs(1:nlat-1)

    INTEGER :: issh=0,isst=0,isss=0,issu=0,issv=0
    INTEGER :: it=0,is=0,iu=0,iv=0
    INTEGER :: numdelqc1=0, numdelqc2=0

    REAL(r_size) obserrt
    REAL(r_size) dob2_sprd2 !dob(obsdep)**2 - sprd**2

    CHARACTER(8) :: obsfile="obsTT.nc"

    !---MPI
    INTEGER MPI_R_SIZE,ierr
    
    !---tmp
    INTEGER,ALLOCATABLE :: tmpqc(:)    !QC
    INTEGER,ALLOCATABLE :: tmpidx(:),tmpidy(:)
    INTEGER,ALLOCATABLE :: tmpelm(:),tmpins(:)

    REAL(r_size),ALLOCATABLE :: tmplon(:),tmplat(:),tmplev(:)
    REAL(r_size),ALLOCATABLE :: tmpdat(:)
    REAL(r_size),ALLOCATABLE :: tmperr(:)
    REAL(r_size),ALLOCATABLE :: tmpdep(:)    !Inovation
    REAL(r_size),ALLOCATABLE :: tmphxfmean(:),tmphxfsprd(:) !H(xfmean/sprd)
    REAL(r_size),ALLOCATABLE :: tmphdxf(:,:) !dY=HdX Ensemble perturvation in obs. space

    !---OUT
    REAL(r_size),INTENT(OUT) :: fcst3d(nij1,nlev,nmem,nv3d),fcst2d(nij1,nmem,nv2d)

    !---MPI(r_size)
    IF(r_size == kind(0.0d0))THEN
       MPI_R_SIZE=MPI_DOUBLE_PRECISION
    ELSE IF(r_size == kind(0.0e0))THEN
       MPI_R_SIZE=MPI_REAL
    END IF
    
    !---Get Nobs
    DO islot=1,nslots
       WRITE(obsfile(4:5),'(I2.2)') islot
       CALL get_nobs(obsfile,nobslots(islot))
    END DO !islot
    
    nobs = SUM(nobslots)
    IF(myrank == root)THEN
       WRITE(6,'(A,I10)') "TOTAL NUMBER of INPUT OBSERVATIONS:",nobs
    END IF
       
    !---ALLOCATE GLOBAL OBSERVATION
    ALLOCATE( tmpqc(nobs) )
    ALLOCATE( tmpidx(nobs), tmpidy(nobs) )
    ALLOCATE( tmpelm(nobs), tmpins(nobs) )
    ALLOCATE( tmplon(nobs), tmplat(nobs), tmplev(nobs) )
    ALLOCATE( tmpdat(nobs), tmperr(nobs), tmpdep(nobs) )
    ALLOCATE( tmphxfmean(nobs), tmphxfsprd(nobs) )
    ALLOCATE( tmphdxf(nobs,nmem))
    
    !---LOOP of timeslots
    no=0

    DO islot=1,nslots

       IF(nobslots(islot) == 0) CYCLE

       WRITE(obsfile(4:5),'(I2.2)') islot
       no1=no+1
       no2=no+nobslots(islot)

       !---Read obs data
       CALL read_obs(obsfile,nobslots(islot), &
            & tmpelm(no1:no2), tmpins(no1:no2), &
            & tmplon(no1:no2), tmplat(no1:no2), tmplev(no1:no2), &
            & tmpdat(no1:no2), tmperr(no1:no2))
       
       !---Get obs. ID
       IF(igrid_type == 1)THEN
          CALL obs_id_uniform(nobslots(islot), &
               & tmplon(no1:no2),tmplat(no1:no2), &
               & tmpidx(no1:no2),tmpidy(no1:no2))
       ELSE IF(igrid_type == 2)THEN
          CALL obs_id_nonuniform(nobslots(islot), &
               & tmplon(no1:no2),tmplat(no1:no2), &
               & tmpidx(no1:no2),tmpidy(no1:no2))
       END IF

       !---Read Xf and H(Xf)
       CALL read_ens_mpi(islot,nobslots(islot), &
            & tmpidx(no1:no2),tmpidy(no1:no2),tmpelm(no1:no2), &
            & tmplon(no1:no2),tmplat(no1:no2),tmplev(no1:no2), &
            & tmpqc(no1:no2),tmphdxf(no1:no2,:), &
            & fcst3d,fcst2d)                   
              
       no = no + nobslots(islot)

    END DO !islot
    
    !---Share tmphdxf (Hxf(i)) for all processors
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,tmphdxf,nobs*nmem,MPI_R_SIZE,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,tmpqc,nobs,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    
    !---Hxfmean, dYf, Inovation
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iobs)
    DO iobs=1,nobs

       IF(tmpqc(iobs) == 1)THEN

          !tmphxfmean: Hxfmean
          tmphxfmean(iobs) = sum(tmphdxf(iobs,1:nmem)) / REAL(nmem,r_size)

          !tmphxfsprd: Hxfsprd
          tmphxfsprd(iobs) = REAL(0.d0, r_size)
          DO imem=1,nmem
             tmphxfsprd(iobs) = tmphxfsprd(iobs) + (tmphdxf(iobs,imem)-tmphxfmean(iobs))**2
          END DO

          IF(1 < nmem)THEN
             tmphxfsprd(iobs)=SQRT(tmphxfsprd(iobs)/REAL(nmem-1,r_size))
          ELSE
             tmphxfsprd(iobs)=undef
          END IF
          
          !tmphdxf: dYf=HXf-xfmean*1
          tmphdxf(iobs,:) = tmphdxf(iobs,:)-tmphxfmean(iobs)

          !tmpdep: Inovation y-Hxfmean
          tmpdep(iobs) = tmpdat(iobs) - tmphxfmean(iobs)

       ELSE

          tmphxfmean(iobs)=undef
          tmphxfsprd(iobs)=undef
          tmpdep(iobs)=undef
          tmphdxf(iobs,:)=undef
          tmpqc(iobs)=0

       END IF

    END DO !iobs
    !$OMP END PARALLEL DO
    
    !#OBS deleted by Obs. projection
    numdelqc1=nobs-SUM(tmpqc)
    
    !---Gross error QC
    IF(LGE_IO .and. .not. AOEI)THEN

       !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iobs)
       DO iobs=1,nobs

          IF(tmpelm(iobs) == id_z_obs)THEN !SSH
             CALL check_gross_error_obs(hqc_min,hqc_max,tmpdep(iobs),tmpqc(iobs))
          ELSE IF(tmpelm(iobs) == id_t_obs .and. tmplev(iobs) == 0.d0)THEN !SST
             CALL check_gross_error_obs(sstqc_min,sstqc_max,tmpdep(iobs),tmpqc(iobs))
          ELSE IF(tmpelm(iobs) == id_t_obs)THEN !T
             CALL check_gross_error_obs(tqc_min,tqc_max,tmpdep(iobs),tmpqc(iobs))
          ELSE IF(tmpelm(iobs) == id_s_obs .and. tmplev(iobs) == 0.d0)THEN !SSS
             CALL check_gross_error_obs(sssqc_min,sssqc_max,tmpdep(iobs),tmpqc(iobs))
          ELSE IF(tmpelm(iobs) == id_s_obs)THEN !S
             CALL check_gross_error_obs(sqc_min,sqc_max,tmpdep(iobs),tmpqc(iobs))
          ELSE IF(tmpelm(iobs) == id_u_obs .and. tmplev(iobs) == 0.d0)THEN !SSU
             CALL check_gross_error_obs(ssuqc_min,ssuqc_max,tmpdep(iobs),tmpqc(iobs))
          ELSE IF(tmpelm(iobs) == id_u_obs)THEN !U
             CALL check_gross_error_obs(uqc_min,uqc_max,tmpdep(iobs),tmpqc(iobs))
          ELSE IF(tmpelm(iobs) == id_v_obs .and. tmplev(iobs) == 0.d0)THEN !SSV
             CALL check_gross_error_obs(ssvqc_min,ssvqc_max,tmpdep(iobs),tmpqc(iobs))
          ELSE IF(tmpelm(iobs) == id_v_obs)THEN !V
             CALL check_gross_error_obs(vqc_min,vqc_max,tmpdep(iobs),tmpqc(iobs))
          ELSE IF(ABS(tmpdep(iobs)) > gross_error*tmperr(iobs))THEN
             tmpqc(iobs) = 0
          END IF

       END DO !iobs
       !$OMP END PARALLEL DO

    END IF

    numdelqc2=nobs-numdelqc1-SUM(tmpqc) ! previous-current SUM(tmpqc)

    CALL monit_dep(nobs,tmpelm,tmpdep,tmpqc)

    !--- Temporal observation localization
    no = 0
    DO islot=1,nslots
       no1=no+1
       no2=no+nobslots(islot)
       tmperr(no+1:no2) = tmperr(no+1:no2) &
            & * exp(0.25d0 * (REAL(islot-nbslot,r_size) / sigma_obst)**2)
       no = no + nobslots(islot)
    END DO

    !--- Remove undef
    no = 0
    DO iobs=1,nobs
       IF(tmpqc(iobs) == 1)THEN
          no = no+1
          tmpelm(no) = tmpelm(iobs)
          tmpins(no) = tmpins(iobs)
          tmplon(no) = tmplon(iobs)
          tmplat(no) = tmplat(iobs)
          tmplev(no) = tmplev(iobs)
          tmpdat(no) = tmpdat(iobs)
          tmperr(no) = tmperr(iobs)
          tmpidx(no) = tmpidx(iobs)
          tmpidy(no) = tmpidy(iobs)
          tmpdep(no) = tmpdep(iobs)
          tmphxfmean(no) = tmphxfmean(iobs)
          tmphxfsprd(no) = tmphxfsprd(iobs)
          tmphdxf(no,:) = tmphdxf(iobs,:)
          tmpqc(no) = tmpqc(iobs)
       END IF
    END DO !iobs

    !--- Total #OBS without undef
    nobs = no

    IF(myrank == root)THEN
       WRITE(6,'(A,I10)') "TOTAL ASSIMILATED OBSERVATION: ",nobs
       WRITE(6,'(A,I10)')   "TOTAL OBSERVATION DELETED by Obs. Projection: ",numdelqc1
       WRITE(6,'(A,I10,A)') "TOTAL OBSERVATION. DELETED by Large Innovation: ",numdelqc2
    END IF
    
    !--- SORT
    ALLOCATE( obsidx(nobs), obsidy(nobs) )    
    ALLOCATE( obselm(nobs), obsins(nobs) )
    ALLOCATE( obslon(nobs), obslat(nobs), obslev(nobs) )
    ALLOCATE( obsdat(nobs), obserr(nobs), obsdep(nobs) )
    ALLOCATE( obshxfmean(nobs), obshxfsprd(nobs))
    ALLOCATE( obshdxf(nobs,nmem) )

    nobsgrd(:,:) = 0
    nj(:) = 0

    !nj: #OBS at a latitude band
    !$OMP PARALLEL PRIVATE(i,j,iobs,no)
    !$OMP DO SCHEDULE(STATIC)
    DO j=1,nlat-1
       DO iobs=1,nobs
          IF(tmpidy(iobs) < j .OR. j+1 <= tmpidy(iobs)) CYCLE
          nj(j) = nj(j) + 1
       END DO !iobs
    END DO !j
    !$OMP END DO

    !njs: #OBS for 0 ~ j-1 at j --> Start point
    !$OMP DO SCHEDULE(STATIC)
    DO j=1,nlat-1
       njs(j) = SUM(nj(0:j-1))
    END DO !j
    !$OMP END DO

    !---Sort for latitude direction --> tmp2
    !$OMP DO SCHEDULE(STATIC)
    DO j=1,nlat-1
       no = 0
       DO iobs=1,nobs
          IF(tmpidy(iobs) < j .OR. j+1 <= tmpidy(iobs)) CYCLE
          no = no + 1
          obsidx(njs(j)+no) = tmpidx(iobs)
          obsidy(njs(j)+no) = tmpidy(iobs)
          obselm(njs(j)+no) = tmpelm(iobs)
          obsins(njs(j)+no) = tmpins(iobs)
          obslon(njs(j)+no) = tmplon(iobs)
          obslat(njs(j)+no) = tmplat(iobs)
          obslev(njs(j)+no) = tmplev(iobs)
          obsdat(njs(j)+no) = tmpdat(iobs)
          obserr(njs(j)+no) = tmperr(iobs)
          obsdep(njs(j)+no) = tmpdep(iobs)
          obshxfmean(njs(j)+no) = tmphxfmean(iobs)
          obshxfsprd(njs(j)+no) = tmphxfsprd(iobs)
          obshdxf(njs(j)+no,:) = tmphdxf(iobs,:)
       END DO !iobs
    END DO !j
    !$OMP END DO

    !Sort for longitude direction --> obs
    !$OMP DO SCHEDULE(STATIC)
    DO j=1,nlat-1

       IF(nj(j) == 0)THEN
          nobsgrd(:,j) = njs(j)
          CYCLE
       END IF

       no = 0
       DO i=1,nlon
          DO iobs=njs(j)+1,njs(j)+nj(j)
             IF(obsidx(iobs) < i .OR. i+1 <= obsidx(iobs)) CYCLE
             no = no + 1
             tmpidx(njs(j)+no) = obsidx(iobs)
             tmpidy(njs(j)+no) = obsidy(iobs)
             tmpelm(njs(j)+no) = obselm(iobs)
             tmpins(njs(j)+no) = obsins(iobs)
             tmplon(njs(j)+no) = obslon(iobs)
             tmplat(njs(j)+no) = obslat(iobs)
             tmplev(njs(j)+no) = obslev(iobs)
             tmpdat(njs(j)+no) = obsdat(iobs)
             tmperr(njs(j)+no) = obserr(iobs)
             tmpdep(njs(j)+no) = obsdep(iobs)
             tmphxfmean(njs(j)+no) = obshxfmean(iobs)
             tmphxfsprd(njs(j)+no) = obshxfsprd(iobs)
             tmphdxf(njs(j)+no,:) = obshdxf(iobs,:)
          END DO !n
          nobsgrd(i,j) = njs(j) + no
       END DO !i

       !Error check
       IF(no /= nj(j)) THEN
          !$OMP CRITICAL
          WRITE(6,'(A,2I)') 'OBS DATA SORT ERROR: ',no,nj(j)
          WRITE(6,'(F6.2,A,F6.2)') j,'< J <',j+1
          WRITE(6,'(F6.2,A,F6.2)') &
               & MINVAL(obsidy(njs(j)+1:njs(j)+nj(j))),'< OBSJ <',MAXVAL(obsidy(njs(j)+1:njs(j)+nj(j)))
          !$OMP END CRITICAL
       END IF

    END DO !j
    !$OMP END DO
    !$OMP END PARALLEL

    !---Exchange tmp --> obs
    obsidx(1:nobs)=tmpidx(1:nobs)
    obsidy(1:nobs)=tmpidy(1:nobs)
    obselm(1:nobs)=tmpelm(1:nobs)
    obsins(1:nobs)=tmpins(1:nobs)
    obslon(1:nobs)=tmplon(1:nobs)
    obslat(1:nobs)=tmplat(1:nobs)
    obslev(1:nobs)=tmplev(1:nobs)
    obsdat(1:nobs)=tmpdat(1:nobs)
    obserr(1:nobs)=tmperr(1:nobs)
    obsdep(1:nobs)=tmpdep(1:nobs)
    obshxfmean(1:nobs)=tmphxfmean(1:nobs)
    obshxfsprd(1:nobs)=tmphxfsprd(1:nobs)
    obshdxf(1:nobs,:)=tmphdxf(1:nobs,:)

    !---Deallocate
    DEALLOCATE( tmpidx, tmpidy )
    DEALLOCATE( tmpelm, tmpins )
    DEALLOCATE( tmplon, tmplat, tmplev )
    DEALLOCATE( tmpdat, tmperr, tmpdep )
    DEALLOCATE( tmphxfmean, tmphxfsprd )
    DEALLOCATE( tmphdxf )
    DEALLOCATE( tmpqc )
    
    !---AOEI (Minamide and Zhang 2017; Monthly Weather Review) S.Ohishi 2019.08
    IF(AOEI)THEN
       
       !Update obserr
       !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iobs,obserrt,dob2_sprd2)
       DO iobs=1,nobs
          obserrt = obserr(iobs) !Prescribed obs. err
          dob2_sprd2 = obsdep(iobs)*obsdep(iobs)-obshxfsprd(iobs)*obshxfsprd(iobs) !dob**2 - sprd**2
          obserr(iobs) = sqrt( max(obserrt*obserrt,dob2_sprd2) )
       END DO !iobs
       !$OMP END PARALLEL DO
       
    END IF
    
  END SUBROUTINE set_obs

  !-----------------------------------------------------------------------
  ! Check gross error observation
  !-----------------------------------------------------------------------

  SUBROUTINE check_gross_error_obs(qc_min,qc_max,inv,qc)

    USE common_setting, only: r_size
    IMPLICIT NONE

    !IN
    REAL(r_size),INTENT(IN) :: qc_min,qc_max !Min. and Max. Error gross value
    REAL(r_size),INTENT(IN) :: inv           !Inovation

    !INOUT
    INTEGER,INTENT(INOUT) :: qc

    IF(inv < qc_min .or. qc_max < inv)THEN
       qc=0
    END IF

  END SUBROUTINE check_gross_error_obs

END MODULE common_obs
