PROGRAM letkf

  !=======================================================================
  !
  ! [PURPOSE:] Main program
  !
  ! [HISTORY:]
  !   2009/01/16 Takemasa Miyoshi  created
  !   2009/02/03 Takemasa Miyoshi  modified for ROMS
  !   2011/01/26 Yasumasa Miyazawa  modified for POM (check 'pom' or 'POM')
  !   2024/04/01 Shun Ohishi        modified for sbPOM-LETKF v2
  !   2025/08/08 Shun Ohishi        modified for LPFGM
  !   2026/06/04 Shun Ohishi        updated
  !   2026/10/21 Shun Ohishi        updated by using Claude
  !=======================================================================

  USE common  
  USE common_da
  USE common_io
  USE common_obs
  USE common_mpi
  USE common_pom
  USE common_setting
  USE MPI
  IMPLICIT NONE

  INTEGER ierr

  REAL(r_size),ALLOCATABLE :: ens3d(:,:,:,:),ens2d(:,:,:)   !Ensemble Forecast => Ensemble analysis
  REAL(r_size),ALLOCATABLE :: fcst3d_m(:,:,:),fcst2d_m(:,:) !Forecast ensemble mean
  REAL(r_size),ALLOCATABLE :: fcst3d_s(:,:,:),fcst2d_s(:,:) !                  spread

  !-----------------------------------------------------------------------
  ! Initial settings
  !-----------------------------------------------------------------------

  CALL initialize_mpi
  CALL TIMER_INIT

  IF(myrank == root)THEN
     WRITE(6,'(A)') '============================================='
     WRITE(6,'(A)') '  LOCAL ENSEMBLE TRANSFORM KALMAN FILTER     '
     WRITE(6,'(A)') '                                             '
     WRITE(6,'(A)') '   LL      EEEEEE  TTTTTT  KK  KK  FFFFFF    '
     WRITE(6,'(A)') '   LL      EE        TT    KK KK   FF        '
     WRITE(6,'(A)') '   LL      EEEEE     TT    KKK     FFFFF     '
     WRITE(6,'(A)') '   LL      EE        TT    KK KK   FF        '
     WRITE(6,'(A)') '   LLLLLL  EEEEEE    TT    KK  KK  FF        '
     WRITE(6,'(A)') '                                             '
     WRITE(6,'(A)') '             WITHOUT LOCAL PATCH             '
     WRITE(6,'(A)') '                                             '
     WRITE(6,'(A)') '          Coded by Takemasa Miyoshi          '
     WRITE(6,'(A)') '  Based on Ott et al (2004) and Hunt (2005)  '
     WRITE(6,'(A)') '  Tested by Miyoshi and Yamane (2006)        '
     WRITE(6,'(A)') '  Developed by Ohishi et al. (2022a,b,2026)  '
     WRITE(6,'(A)') '============================================='
     WRITE(6,'(A)') '              LETKF PARAMETERS               '
     WRITE(6,'(A)') ' ------------------------------------------- '
     WRITE(6,'(A,I15)')  '   nmem       :',nmem
     WRITE(6,'(A,I15)')  '   nslots     :',nslots
     WRITE(6,'(A,I15)')  '   nbslot     :',nbslot
     WRITE(6,'(A,I15)')  '   iswitch_da :',iswitch_da
     WRITE(6,'(A,L15)')  '   AOEI       :',AOEI
     WRITE(6,'(A,F15.2)')'   sigma_obs [km]  :',sigma_obs*1.d-3
     WRITE(6,'(A,F15.2)')'   sigma_obsv [m]  :',sigma_obsv
     WRITE(6,'(A,F15.2)')'   sigma_obst [day]:',sigma_obst
     WRITE(6,'(A,F15.2)')'   MULT parameter:',cov_infl_mul
     WRITE(6,'(A,F15.2)')'   RTPP parameter:',ALPHA_RTPP
     WRITE(6,'(A,F15.2)')'   RTPS parameter:',ALPHA_RTPS
     IF(iswitch_da == 2 .or. iswitch_da == 3)THEN
        WRITE(6,'(A,I15)')  '   nmonte     :',nmonte
        WRITE(6,'(A,A15)')  '   RM         :',RM
        WRITE(6,'(A,L15)')  '   DP         :',DP
        WRITE(6,'(A,L15)')  '   LRAND      :',LRAND
        WRITE(6,'(A,F15.2)')'   m_crit     :',m_crit
        WRITE(6,'(A,F15.2)')'   GM infl.   :',cov_infl_gm
     END IF
     WRITE(6,'(A)') '============================================='
  END IF

  IF(myrank == root &
       & .and. (ALPHA_RTPP > REAL(0.d0,r_size) .and. ALPHA_RTPS > REAL(0.d0,r_size)))THEN
     write(6,'(A)') "***Error: Positive both RTPP and RTPS parameters"
     CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
  END IF

  CALL set_common_pom
  CALL set_common_mpi_pom

  ALLOCATE(ens3d(nij1,nlev,nmem,nv3d),ens2d(nij1,nmem,nv2d))

  CALL TIMER("INITIAL SETTING")

  !-----------------------------------------------------------------------
  ! Observations and Ensemble forecasts in obs. space
  !-----------------------------------------------------------------------

  CALL set_obs(ens3d,ens2d) !forecast

  CALL TIMER("READ OBS and dYf")

  !-----------------------------------------------------------------------
  ! Forecast ensemble mean and spread
  !-----------------------------------------------------------------------

  ALLOCATE(fcst3d_m(nij1,nlev,nv3d),fcst2d_m(nij1,nv2d))
  ALLOCATE(fcst3d_s(nij1,nlev,nv3d),fcst2d_s(nij1,nv2d))

  call ensemble_mesp(ens3d,ens2d,fcst3d_m,fcst3d_s,fcst2d_m,fcst2d_s)

  CALL TIMER("Forecast ensemble mean and spread")

  !-----------------------------------------------------------------------
  ! Analysis
  !-----------------------------------------------------------------------

  CALL da_main(ens3d,ens2d) !forecast -> analysis

  !---Deallocate
  CALL deallocate_shared_hdxfT !DEALLOCATE(hdxfT)  

  CALL TIMER("ANALYSIS")

  !-----------------------------------------------------------------------
  ! Write data |
  !-------------
  ! - Forecast ensemble mean and spread
  ! - Analysis ensemble mean and spread
  ! - Ensemble analyses 
  !-----------------------------------------------------------------------

  CALL write_ens_mpi(ens3d,ens2d,fcst3d_m,fcst3d_s,fcst2d_m,fcst2d_s) !analysis

  CALL TIMER("WRITE ANALYSIS")

  !-----------------------------------------------------------------------
  ! Finalize
  !-----------------------------------------------------------------------

  DEALLOCATE(fcst3d_m,fcst2d_m)
  DEALLOCATE(fcst3d_s,fcst2d_s)
  DEALLOCATE(ens3d,ens2d)

  CALL TIMER_END

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)  
  CALL MPI_FINALIZE(ierr)

END PROGRAM letkf
