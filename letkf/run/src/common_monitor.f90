MODULE common_monit

  !=======================================================================
  !
  ! [PURPOSE:] Monitor statistical values
  !
  ! [HISTORY:]
  !   2026/10/21 Shun OHISHI moved from common.f90
  !
  !=======================================================================

  IMPLICIT NONE

CONTAINS
  
  !--------------------------------------------------------------------------------
  ! Monitor departure of forecast and analysis mean relative to assimilated obs.|  
  !--------------------------------------------------------------------------------

  SUBROUTINE monit_mean(name,hxmean)

    !$ USE OMP_LIB
    USE common_setting
    IMPLICIT NONE

    !---Common
    INTEGER i,ivd

    REAL(r_dble) misfit  !Misfit: Hx-y

    !3D+2D variable
    INTEGER no(nv3d+nv2d)
    REAL(r_dble) bias(nv3d+nv2d),rmsd(nv3d+nv2d)

    !3D variable
    INTEGER no3ds(nv3d),no3di(nv3d)
    REAL(r_dble) bias3ds(nv3d),rmsd3ds(nv3d) !Sea surface
    REAL(r_dble) bias3di(nv3d),rmsd3di(nv3d) !Internal

    !---IN    
    CHARACTER(4),INTENT(IN) :: name

    REAL(r_size),INTENT(IN) :: hxmean(nobs) !Ensemble mean in obs space

    !---Initialization
    no(:) = 0
    no3ds(:) = 0
    no3di(:) = 0

    bias(:) = REAL(0.d0, r_dble)
    bias3ds(:) = REAL(0.d0, r_dble)
    bias3di(:) = REAL(0.d0, r_dble)

    rmsd(:) = REAL(0.d0, r_dble)   
    rmsd3ds(:) = REAL(0.d0, r_dble)
    rmsd3di(:) = REAL(0.d0, r_dble)

    !---Innovation
    !$OMP PARALLEL
    !$OMP DO PRIVATE(i,ivd,misfit) REDUCTION(+:no,no3ds,no3di,bias,bias3ds,bias3di,rmsd,rmsd3ds,rmsd3di) SCHEDULE(STATIC)
    DO i=1,nobs

       IF(hxmean(i) == undef) CYCLE

       !Innovation
       misfit = hxmean(i)-obsdat(i) 

       !Get ID
       SELECT CASE(obselm(i))
       CASE(id_u_obs)
          ivd=iv3d_u
       CASE(id_v_obs)
          ivd=iv3d_v
       CASE(id_t_obs)
          ivd=iv3d_t
       CASE(id_s_obs)
          ivd=iv3d_s
       CASE(id_z_obs)
          ivd=nv3d+iv2d_z
       CASE DEFAULT
          CYCLE
       END SELECT

       no(ivd)=no(ivd)+1
       bias(ivd)=bias(ivd)+misfit
       rmsd(ivd)=rmsd(ivd)+misfit**2

       IF(obselm(i) == id_z_obs) CYCLE

       !3D surface and internal
       IF(obslev(i) == REAL(0.d0,r_size))THEN
          no3ds(ivd)=no3ds(ivd)+1
          bias3ds(ivd)=bias3ds(ivd)+misfit
          rmsd3ds(ivd)=rmsd3ds(ivd)+misfit**2
       ELSE
          no3di(ivd)=no3di(ivd)+1
          bias3di(ivd)=bias3di(ivd)+misfit
          rmsd3di(ivd)=rmsd3di(ivd)+misfit**2
       END IF

    END DO !i
    !$OMP END DO

    !$OMP DO SCHEDULE(STATIC)
    DO ivd=1,nv3d+nv2d
       IF(no(ivd) == 0)THEN
          bias(ivd)=undef
          rmsd(ivd)=undef
       ELSE
          bias(ivd)=bias(ivd)/REAL(no(ivd),r_dble)
          rmsd(ivd)=SQRT(rmsd(ivd)/REAL(no(ivd),r_dble))
       END IF
    END DO !ivd
    !$OMP END DO

    !$OMP DO SCHEDULE(STATIC)       
    DO ivd=1,nv3d
       !Surface
       IF(no3ds(ivd) == 0)THEN
          bias3ds(ivd)=undef
          rmsd3ds(ivd)=undef
       ELSE
          bias3ds(ivd)=bias3ds(ivd)/REAL(no3ds(ivd),r_dble)
          rmsd3ds(ivd)=SQRT(rmsd3ds(ivd)/REAL(no3ds(ivd),r_dble))
       END IF

       !Internal
       IF(no3di(ivd) == 0)THEN
          bias3di(ivd)=undef
          rmsd3di(ivd)=undef
       ELSE
          bias3di(ivd)=bias3di(ivd)/REAL(no3di(ivd),r_dble)
          rmsd3di(ivd)=SQRT(rmsd3di(ivd)/REAL(no3di(ivd),r_dble))
       END IF
    END DO !ivd
    !$OMP END DO
    !$OMP END PARALLEL

    !---Write Information       
    WRITE(6,'(A,I10)') "TOTAL #OBS:",nobs
    WRITE(6,'(A)') "== STATISTICS IN OBSERVATION SPACE("//name//") =================="
    WRITE(6,'(14A11)') &
         & "Variable:", &
         & "U","SSU","Interior U", &
         & "V","SSV","Interior V", &
         & "T","SST","Interior T", &
         & "S","SSS","Interior S", &
         & "SSH"
    WRITE(6,'(A11,13F11.5)') &
         "Bias:", &
         & bias(iv3d_u),bias3ds(iv3d_u),bias3di(iv3d_u), &
         & bias(iv3d_v),bias3ds(iv3d_v),bias3di(iv3d_v), &
         & bias(iv3d_t),bias3ds(iv3d_t),bias3di(iv3d_t), &
         & bias(iv3d_s),bias3ds(iv3d_s),bias3di(iv3d_s), &
         & bias(nv3d+iv2d_z)
    WRITE(6,'(A11,13F11.5)') &
         "RMSD:", &
         & rmsd(iv3d_u),rmsd3ds(iv3d_u),rmsd3di(iv3d_u), &
         & rmsd(iv3d_v),rmsd3ds(iv3d_v),rmsd3di(iv3d_v), &
         & rmsd(iv3d_t),rmsd3ds(iv3d_t),rmsd3di(iv3d_t), &
         & rmsd(iv3d_s),rmsd3ds(iv3d_s),rmsd3di(iv3d_s), &
         & rmsd(nv3d+iv2d_z)

    WRITE(6,'(A)') "== NUMBER OF OBSERVATIONS =================================="
    WRITE(6,'(14A11)') &
         & "Variable:", &
         & "U","SSU","Interior U", &
         & "V","SSV","Interior V", &
         & "T","SST","Interior T", &
         & "S","SSS","Interior S", &
         & "SSH"
    WRITE(6,'(A11,13I11)') &
         & "#OBS:", &
         & no(iv3d_u),no3ds(iv3d_u),no3di(iv3d_u), &
         & no(iv3d_v),no3ds(iv3d_v),no3di(iv3d_v), &
         & no(iv3d_t),no3ds(iv3d_t),no3di(iv3d_t), &
         & no(iv3d_s),no3ds(iv3d_s),no3di(iv3d_s), &
         & no(nv3d+iv2d_z)
    WRITE(6,'(A)') "============================================================"

  END SUBROUTINE monit_mean

  !--------------------------------------------------------------------------------
  ! Monitor forecast and analysis sprd in obs. space |  
  !--------------------------------------------------------------------------------

  SUBROUTINE monit_sprd(name,hxsprd)

    !$ USE OMP_LIB
    USE common_setting
    IMPLICIT NONE

    !---Common
    INTEGER i,ivd

    !3D+2D variable
    INTEGER no(nv3d+nv2d)
    REAL(r_dble) sprd(nv3d+nv2d)

    !3D variable
    INTEGER no3ds(nv3d),no3di(nv3d)
    REAL(r_dble) sprd3ds(nv3d) !Sea surface
    REAL(r_dble) sprd3di(nv3d) !Internal

    !---IN    
    CHARACTER(4),INTENT(IN) :: name

    REAL(r_size),INTENT(IN) :: hxsprd(nobs) !Ensemble spread in obs space

    !---Initialization
    no(:) = 0
    no3ds(:) = 0
    no3di(:) = 0

    sprd(:) = REAL(0.d0, r_dble)
    sprd3ds(:) = REAL(0.d0, r_dble)
    sprd3di(:) = REAL(0.d0, r_dble)

    !---Innovation
    !$OMP PARALLEL
    !$OMP DO PRIVATE(i,ivd) SCHEDULE(STATIC) REDUCTION(+:no,no3ds,no3di,sprd,sprd3ds,sprd3di)
    DO i=1,nobs

       IF(hxsprd(i) == undef) CYCLE

       !Get ID
       SELECT CASE(obselm(i))
       CASE(id_u_obs)
          ivd=iv3d_u
       CASE(id_v_obs)
          ivd=iv3d_v
       CASE(id_t_obs)
          ivd=iv3d_t
       CASE(id_s_obs)
          ivd=iv3d_s
       CASE(id_z_obs)
          ivd=nv3d+iv2d_z
       CASE DEFAULT
          CYCLE
       END SELECT

       no(ivd)=no(ivd)+1
       sprd(ivd)=sprd(ivd)+hxsprd(i)**2

       IF(obselm(i) == id_z_obs) CYCLE

       !3D surface and internal
       IF(obslev(i) == REAL(0.d0,r_size))THEN
          no3ds(ivd)=no3ds(ivd)+1
          sprd3ds(ivd)=sprd3ds(ivd)+hxsprd(i)**2
       ELSE
          no3di(ivd)=no3di(ivd)+1
          sprd3di(ivd)=sprd3di(ivd)+hxsprd(i)**2
       END IF

    END DO !i
    !$OMP END DO

    !$OMP DO PRIVATE(ivd) SCHEDULE(STATIC)
    DO ivd=1,nv3d+nv2d
       IF(no(ivd) == 0)THEN
          sprd(ivd)=undef
       ELSE
          sprd(ivd)=SQRT(sprd(ivd)/REAL(no(ivd),r_dble))
       END IF
    END DO !ivd
    !$OMP END DO

    !$OMP DO PRIVATE(ivd) SCHEDULE(STATIC)      
    DO ivd=1,nv3d
       !Surface
       IF(no3ds(ivd) == 0)THEN
          sprd3ds(ivd)=undef
       ELSE
          sprd3ds(ivd)=SQRT(sprd3ds(ivd)/REAL(no3ds(ivd),r_dble))
       END IF
       !Internal
       IF(no3di(ivd) == 0)THEN
          sprd3di(ivd)=undef
       ELSE
          sprd3di(ivd)=SQRT(sprd3di(ivd)/REAL(no3di(ivd),r_dble))
       END IF
    END DO !ivd
    !$OMP END DO
    !$OMP END PARALLEL

    !---Write Information
    WRITE(6,'(A,I10)') "TOTAL #OBS:",nobs
    WRITE(6,'(A)') "== SPREAD IN OBSERVATION SPACE ("//name//") =================="
    WRITE(6,'(14A11)') &
         & "Variable:", &
         & "U","SSU","Interior U", &
         & "V","SSV","Interior V", &
         & "T","SST","Interior T", &
         & "S","SSS","Interior S", &
         & "SSH"
    WRITE(6,'(A11,13F11.5)') &
         "Spread:", &
         & sprd(iv3d_u),sprd3ds(iv3d_u),sprd3di(iv3d_u), &
         & sprd(iv3d_v),sprd3ds(iv3d_v),sprd3di(iv3d_v), &
         & sprd(iv3d_t),sprd3ds(iv3d_t),sprd3di(iv3d_t), &
         & sprd(iv3d_s),sprd3ds(iv3d_s),sprd3di(iv3d_s), &
         & sprd(nv3d+iv2d_z)
    WRITE(6,'(A)') "== NUMBER OF OBSERVATIONS =================================="
    WRITE(6,'(14A11)') &
         & "Variable:", &
         & "U","SSU","Interior U", &
         & "V","SSV","Interior V", &
         & "T","SST","Interior T", &
         & "S","SSS","Interior S", &
         & "SSH"
    WRITE(6,'(A11,13I11)') &
         & "#OBS:", &
         & no(iv3d_u),no3ds(iv3d_u),no3di(iv3d_u), &
         & no(iv3d_v),no3ds(iv3d_v),no3di(iv3d_v), &
         & no(iv3d_t),no3ds(iv3d_t),no3di(iv3d_t), &
         & no(iv3d_s),no3ds(iv3d_s),no3di(iv3d_s), &
         & no(nv3d+iv2d_z)
    WRITE(6,'(A)') "============================================================"

  END SUBROUTINE monit_sprd

END MODULE common_monit
