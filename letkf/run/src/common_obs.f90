MODULE common_obs
  
  !=======================================================================
  !
  ! [PURPOSE:] Get observation and dYf information
  !
  ! [HISTORY:]
  !   2009/01/23 Takemasa MIYOSHI  created
  !   2011/01/26 Yasumasa MIYAZAWA  modified for POM (check 'pom' or 'POM')
  !   2025/07/31 Shun OHISHI modified
  !   2026/06/04 Shun OHISHI updated (Memory and I/O issue)
  !   2026/10/23 Shun OHISHI updated (by Claude)
  !
  !=======================================================================

CONTAINS

  !-----------------------------------------------------------------------
  ! Set information in observation space
  !-----------------------------------------------------------------------

  SUBROUTINE set_obs(fcst3d,fcst2d)

    !$ USE OMP_LIB
    USE common
    USE common_io
    USE common_monit
    USE common_mpi
    USE common_setting
    USE mpi
    IMPLICIT NONE

    !---Parameter
    INTEGER(8),PARAMETER :: nchunk=16777216_8 !#elements per MPI_ALLREDUCE (64 MB)

    !---Common
    INTEGER i,j,k
    INTEGER islot
    INTEGER iobs
    INTEGER io,no,no1,no2
        
    INTEGER nobslots(nslots)
    INTEGER,ALLOCATABLE :: nobs_k(:),nobs_kbsum(:)
    INTEGER,ALLOCATABLE :: ids(:)

    INTEGER numdelqc

    REAL(r_size) obserr_pre !Prescribed obs error
    REAL(r_size) obserr_est !dob**2 - hxfsprd**2

    CHARACTER(8) :: obsfile="obsTT.nc"

    !---MPI
    INTEGER(8) ntot,is8,ie8
    INTEGER ierr
    
    !---tmp
    INTEGER nobs_tmp
    
    INTEGER,ALLOCATABLE :: tmpqc(:)    !QC (1:OK, 0: NG)
    INTEGER,ALLOCATABLE :: tmpidt(:)
    INTEGER,ALLOCATABLE :: tmpidxt(:),tmpidyt(:)
    INTEGER,ALLOCATABLE :: tmpidxu(:),tmpidyu(:)
    INTEGER,ALLOCATABLE :: tmpidxv(:),tmpidyv(:)
    INTEGER,ALLOCATABLE :: tmpelm(:),tmpins(:)

    REAL(r_size),ALLOCATABLE :: tmplon(:),tmplat(:),tmplev(:)
    REAL(r_size),ALLOCATABLE :: tmpdat(:)
    REAL(r_size),ALLOCATABLE :: tmperr(:)
    
    !---OUT
    REAL(r_size),INTENT(OUT) :: fcst3d(nij1,nlev,nmem,nv3d),fcst2d(nij1,nmem,nv2d)
    
    !=====================================================================
    ! 1. Read all observations & QC (geometry only)
    !=====================================================================
    
    !---Get total observation size
    DO islot=1,nslots
       WRITE(obsfile(4:5),'(I2.2)') islot
       CALL get_nobs(obsfile,nobslots(islot))
    END DO
    
    nobs_tmp = SUM(nobslots)

    IF(myrank == root)THEN
       WRITE(6,'(A,I0)') "TOTAL NUMBER of INPUT OBSERVATIONS:",nobs_tmp
    END IF
       
    !---ALLOCATE GLOBAL OBSERVATION
    ALLOCATE( tmpqc(nobs_tmp) )
    ALLOCATE( tmpidt(nobs_tmp))
    ALLOCATE( tmpidxt(nobs_tmp), tmpidyt(nobs_tmp) )
    ALLOCATE( tmpidxu(nobs_tmp), tmpidyu(nobs_tmp) )
    ALLOCATE( tmpidxv(nobs_tmp), tmpidyv(nobs_tmp) )
    ALLOCATE( tmpelm(nobs_tmp), tmpins(nobs_tmp) )
    ALLOCATE( tmplon(nobs_tmp), tmplat(nobs_tmp), tmplev(nobs_tmp) )
    ALLOCATE( tmpdat(nobs_tmp), tmperr(nobs_tmp) )
    
    !---LOOP of timeslots to read all observations
    no=0
    
    DO islot=1,nslots

       !---Check obs size
       IF(nobslots(islot) == 0) CYCLE

       !---Obs filename
       WRITE(obsfile(4:5),'(I2.2)') islot

       !---Observation range
       no1=no+1
       no2=no+nobslots(islot)

       !---Read obs data
       CALL read_obs(obsfile,nobslots(islot), &
            & tmpelm(no1:no2), tmpins(no1:no2), &
            & tmplon(no1:no2), tmplat(no1:no2), tmplev(no1:no2), &
            & tmpdat(no1:no2), tmperr(no1:no2))
       
       !---Temporal obs ID
       tmpidt(no1:no2)=islot

       !---Temporal observation localization
       IF(islot /= nbslot)THEN

          !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(io)
          DO io=no1,no2
             tmperr(io) = tmperr(io) &
                  & * exp(0.25d0 * (REAL(islot-nbslot,r_size) / sigma_obst)**2)
          END DO
          !$OMP END PARALLEL DO

       END IF

       !---no
       no=no+nobslots(islot)
       
    END DO !islot
    
    !---Get idx,idy
    IF(igrid_type == 1)THEN
       CALL obs_id_uniform(nobs_tmp,tmplon,tmplat,lont,latt,tmpidxt,tmpidyt)
       CALL obs_id_uniform(nobs_tmp,tmplon,tmplat,lonu,latu,tmpidxu,tmpidyu)
       CALL obs_id_uniform(nobs_tmp,tmplon,tmplat,lonv,latv,tmpidxv,tmpidyv)
    ELSE IF(igrid_type == 2)THEN
       CALL obs_id_nonuniform(nobs_tmp,tmplon,tmplat,lont,latt,tmpidxt,tmpidyt)
       CALL obs_id_nonuniform(nobs_tmp,tmplon,tmplat,lonu,latu,tmpidxu,tmpidyu)
       CALL obs_id_nonuniform(nobs_tmp,tmplon,tmplat,lonv,latv,tmpidxv,tmpidyv)
    END IF
    
    !---QC (geometry only: grid, land, sea surface/floor)
    CALL obs_check(nobs_tmp,tmpelm,tmpidxt,tmpidxu,tmpidxv,tmpidyt,tmpidyu,tmpidyv,tmplev,tmpqc)
    
    !=====================================================================
    ! 2. Sort assimilated obs. (tmpqc=1) by grid cell (counting sort, stable)
    !    key = (idy-1)*nlon + idx  (i: fastest)
    !    nobsgrd(i,j): cumulative #OBS up to grid cell (i,j)
    !=====================================================================

    ALLOCATE( nobs_k(nlon*nlat), nobs_kbsum(nlon*nlat) )

    !---nobs_k
    nobs_k(:) = 0
    DO iobs=1,nobs_tmp
       IF(tmpqc(iobs) /= 1) CYCLE               !QC
       k = (tmpidyt(iobs)-1)*nlon + tmpidxt(iobs) !kth grid
       nobs_k(k) = nobs_k(k) + 1                !Count #obs at kth grid
    END DO

    !---nobs_kbsum and nobsgrd
    no = 0
    DO j=1,nlat
       DO i=1,nlon
          k = (j-1)*nlon + i  !kth grid
          nobs_kbsum(k) = no  !#Obs before kth grid
          no = no + nobs_k(k) !Sum #obs 1-kth grid
          nobsgrd(i,j) = no   !Cumulative obs total before (i,j) grid 
       END DO
    END DO

    !---Assimilated observation
    nobs = no                  !Assimilated obs size
    numdelqc = nobs_tmp - nobs !Deleted obs size

    ALLOCATE( ids(nobs) )
    DO iobs=1,nobs_tmp
       IF(tmpqc(iobs) /= 1) CYCLE               !QC
       k = (tmpidyt(iobs)-1)*nlon + tmpidxt(iobs) !kth grid
       nobs_kbsum(k) = nobs_kbsum(k) + 1        
       ids(nobs_kbsum(k)) = iobs             !Obs ID for sort
    END DO

    DEALLOCATE( nobs_k, nobs_kbsum )

    IF(myrank == root)THEN
       WRITE(6,'(A,I10)') "TOTAL ASSIMILATED OBSERVATION: ",nobs
       WRITE(6,'(A,I10)') "TOTAL OBSERVATION DELETED by Obs. Projection: ",numdelqc
    END IF

    !---Sorted observations (assimilated obs. only)
    ALLOCATE( obsidt(nobs) )
    ALLOCATE( obsidxt(nobs), obsidyt(nobs) )    
    ALLOCATE( obsidxu(nobs), obsidyu(nobs) )    
    ALLOCATE( obsidxv(nobs), obsidyv(nobs) )    
    ALLOCATE( obselm(nobs), obsins(nobs), obsqc(nobs)  )
    ALLOCATE( obslon(nobs), obslat(nobs), obslev(nobs) )
    ALLOCATE( obsdat(nobs), obserr(nobs) )

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iobs)
    DO iobs=1,nobs
       obsidt(iobs) = tmpidt(ids(iobs))
       obsidxt(iobs) = tmpidxt(ids(iobs))
       obsidxu(iobs) = tmpidxu(ids(iobs))
       obsidxv(iobs) = tmpidxv(ids(iobs))
       obsidyt(iobs) = tmpidyt(ids(iobs))
       obsidyu(iobs) = tmpidyu(ids(iobs))
       obsidyv(iobs) = tmpidyv(ids(iobs))
       obselm(iobs) = tmpelm(ids(iobs))
       obsins(iobs) = tmpins(ids(iobs))
       obsqc(iobs)  = tmpqc(ids(iobs) )
       obslon(iobs) = tmplon(ids(iobs))
       obslat(iobs) = tmplat(ids(iobs))
       obslev(iobs) = tmplev(ids(iobs))
       obsdat(iobs) = tmpdat(ids(iobs))
       obserr(iobs) = tmperr(ids(iobs))
    END DO
    !$OMP END PARALLEL DO

    DEALLOCATE( tmpqc )
    DEALLOCATE( tmpidt)
    DEALLOCATE( tmpidxt, tmpidyt )
    DEALLOCATE( tmpidxu, tmpidyu )
    DEALLOCATE( tmpidxv, tmpidyv )
    DEALLOCATE( tmpelm, tmpins )
    DEALLOCATE( tmplon, tmplat, tmplev )
    DEALLOCATE( tmpdat, tmperr )
    DEALLOCATE( ids )
    
    !=====================================================================
    ! 3. Read Xf and H(Xf) for assimilated obs. only (directly in hdxfT)
    !=====================================================================

    CALL allocate_shared_hdxfT !ALLOCATE(hdxfT(nmem,nobs))
    ALLOCATE(hxfmean(nobs),hxfsprd(nobs),dob(nobs))
        
    !---Read Xf and H(Xf)
    DO islot=1,nslots
       CALL read_ens_mpi(islot,fcst3d,fcst2d)                
    END DO !islot
    
    CALL MPI_WIN_FENCE(0,win_hdxfT,ierr)

    !---Sum over nodes and dYf (node leaders only; once per node)
    IF(node_rank == 0)THEN
       
       ntot=INT(nobs,8)*INT(nmem,8)
       DO is8=1,ntot,nchunk
          ie8=MIN(is8+nchunk-1,ntot)
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,hdxfT_flat(is8:ie8),INT(ie8-is8+1), &
               & MPI_R_SIZE,MPI_SUM,comm_leader,ierr)
       END DO
       CALL ensemble_mesp_obs(hdxfT_flat,hxfmean,hxfsprd)
       
    END IF
    CALL MPI_BCAST(hxfmean,nobs,MPI_R_SIZE,0,comm_node,ierr)
    CALL MPI_BCAST(hxfsprd,nobs,MPI_R_SIZE,0,comm_node,ierr)
    CALL MPI_WIN_FENCE(0,win_hdxfT,ierr)
    
    !---Innovation: y-Hxfmean
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iobs)
    DO iobs=1,nobs
       dob(iobs) = obsdat(iobs) - hxfmean(iobs)
    END DO
    !$OMP END PARALLEL DO
    
    !=====================================================================
    ! 4. AOEI (Adaptive Observation Error Inflation; Minamide and Zhang 2017)
    !=====================================================================    
    IF(AOEI)THEN
       
       !Update obserr
       !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iobs,obserr_pre,obserr_est)
       DO iobs=1,nobs
          obserr_pre = obserr(iobs) !Prescribed obs. err
          obserr_est = dob(iobs)**2-hxfsprd(iobs)**2 !dob**2 - sprd**2
          obserr(iobs) = sqrt( max(obserr_pre**2,obserr_est) )
       END DO
       !$OMP END PARALLEL DO
       
    END IF

    !=====================================================================
    ! 5. Monitor
    !=====================================================================    

    !IF(myrank == root)THEN
    !   CALL monit_mean("fcst",hxfmean)
    !END IF
        
  END SUBROUTINE set_obs

  !-----------------------------------------------------------------------
  ! Check observation |
  !-----------------------------------------------------------------------

  SUBROUTINE obs_check(no,elm,idxt,idxu,idxv,idyt,idyu,idyv,levo,qco)

    !$ USE OMP_LIB
    USE common_setting
    IMPLICIT NONE

    !---Common
    INTEGER io
    INTEGER i,j
    REAL(r_size) mask00,mask10,mask01,mask11
    
    !---IN
    INTEGER,INTENT(IN) :: no
    INTEGER,INTENT(IN) :: elm(no)
    INTEGER,INTENT(IN) :: idxt(no),idxu(no),idxv(no)
    INTEGER,INTENT(IN) :: idyt(no),idyu(no),idyv(no)
    
    REAL(r_size),INTENT(IN) :: levo(no)

    !---OUT
    INTEGER,INTENT(OUT) :: qco(no)

    !---Initialization
    qco(:)=1

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(io,i,j,mask00,mask10,mask01,mask11)
    DO io=1,no

       !---Check element
       SELECT CASE(elm(io))
       CASE(id_z_obs,id_u_obs,id_v_obs,id_t_obs,id_s_obs)
       CASE DEFAULT
          qco(io)=0
          CYCLE
       END SELECT

       !---Grid index at t point (used for sorting and local search of all obs.)
       IF(idxt(io) < 1 .or. nlon-1 < idxt(io) .or. idyt(io) < 1 .or. nlat-1 < idyt(io))THEN
          qco(io)=0
          CYCLE
       END IF
       
       !---Grid index       
       SELECT CASE(elm(io))
       CASE(id_z_obs,id_t_obs,id_s_obs)
          i=idxt(io)
          j=idyt(io)
       CASE(id_u_obs)
          i=idxu(io)
          j=idyu(io)
       CASE(id_v_obs)
          i=idxv(io)
          j=idyv(io)
       END SELECT
       
       IF(i == NINT(undef) .or. i < 1 .or. nlon-1 < i .or. &
            & j == NINT(undef) .or. j < 1 .or. nlat-1 < j)then
          qco(io)=0
          CYCLE
       END IF

       !---Land Check
       SELECT CASE(elm(io))
       CASE(id_z_obs,id_t_obs,id_s_obs)
          mask00=maskt(i  ,j  )
          mask10=maskt(i+1,j  )
          mask01=maskt(i  ,j+1)
          mask11=maskt(i+1,j+1)
       CASE(id_u_obs)
          mask00=masku(i  ,j  )
          mask10=masku(i+1,j  )
          mask01=masku(i  ,j+1)
          mask11=masku(i+1,j+1)
       CASE(id_v_obs)
          mask00=maskv(i  ,j  )
          mask10=maskv(i+1,j  )
          mask01=maskv(i  ,j+1)
          mask11=maskv(i+1,j+1)
       END SELECT
       
       IF(mask00 == land .or. mask10 == land &
            .or. mask01 == land .or. mask11 == land)THEN
          qco(io)=0
          CYCLE
       END IF

       !---Above sea surface
       IF(REAL(0.d0,r_size) < levo(io))THEN
          qco(io)=0
          CYCLE
       END IF

       !---Below Sea bottom
       IF(elm(io) /= id_z_obs)THEN
          IF(levo(io) < depth(i,j,k_btm) .or. levo(io) < depth(i+1,j,k_btm) &
               & .or. levo(io) < depth(i,j+1,k_btm) .or. levo(io) < depth(i+1,j+1,k_btm))THEN
             qco(io)=0
             CYCLE
          END IF
       ENDIF
       
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE obs_check
  
END MODULE common_obs
