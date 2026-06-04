PROGRAM letkf
  
  !=======================================================================
  !
  ! [PURPOSE:] Main program
  !
  ! [HISTORY:]
  !   01/16/2009 Takemasa Miyoshi  created
  !   02/03/2009 Takemasa Miyoshi  modified for ROMS
  !   01/26/2011 Yasumasa Miyazawa  modified for POM (check 'pom' or 'POM')
  !   04/01/2024 Shun Ohishi        modified for sbPOM-LETKF v2
  !   08/08/2025 Shun Ohishi        modified for LPFGM
  !   06/04/2026 Shun Ohishi        updated
  !=======================================================================
  
  USE MPI
  USE common_setting
  USE common
  USE common_mpi
  USE common_io
  USE common_pom
  USE common_obs
  USE common_da

  IMPLICIT NONE

  INTEGER ierr

  REAL(r_size),ALLOCATABLE :: fcst3d(:,:,:,:),fcst2d(:,:,:)
  REAL(r_size),ALLOCATABLE :: anal3d(:,:,:,:),anal2d(:,:,:)
  
  CHARACTER(4) :: fcstf='fc00'

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
     WRITE(6,'(A,I15)')  '   nmem        :',nmem
     WRITE(6,'(A,I15)')  '   nslots     :',nslots
     WRITE(6,'(A,I15)')  '   nbslot     :',nbslot
     WRITE(6,'(A,F15.2)')'   sigma_obs [km]  :',sigma_obs*1.d-3
     WRITE(6,'(A,F15.2)')'   sigma_obsv [m]  :',sigma_obsv
     WRITE(6,'(A,F15.2)')'   sigma_obst [day]:',sigma_obst
     WRITE(6,'(A,F15.2)')'   MULT parameter:',cov_infl_mul
     WRITE(6,'(A,F15.2)')'   RTPP parameter:',ALPHA_RTPP
     WRITE(6,'(A,F15.2)')'   RTPS parameter:',ALPHA_RTPS  
     WRITE(6,'(A)') '============================================='
  END IF
     
  CALL set_common_pom
  CALL set_common_mpi_pom
  
  ALLOCATE(fcst3d(nij1,nlev,nmem,nv3d),fcst2d(nij1,nmem,nv2d))
  ALLOCATE(anal3d(nij1,nlev,nmem,nv3d),anal2d(nij1,nmem,nv2d))

  CALL TIMER("INITIAL SETTING")
    
  !-----------------------------------------------------------------------
  ! Observations and Ensemble forecasts in obs. space
  !-----------------------------------------------------------------------
  
  CALL set_obs(fcst3d,fcst2d)

  CALL TIMER("READ OBS and dYf")
    
  !-----------------------------------------------------------------------
  ! Analysis
  !-----------------------------------------------------------------------

  CALL da_main(fcst3d,fcst2d,anal3d,anal2d)

  CALL TIMER("ANALYSIS")
    
  !-----------------------------------------------------------------------
  ! Write data |
  !-------------
  ! - Forecast ensemble mean and spread
  ! - Analysis ensemble mean and spread
  ! - Ensemble analyses 
  !-----------------------------------------------------------------------

  CALL write_ens_mpi(fcst3d,fcst2d,anal3d,anal2d)
  
  CALL TIMER("WRITE ANALYSIS")
      
  !-----------------------------------------------------------------------
  ! Finalize
  !-----------------------------------------------------------------------

  DEALLOCATE(fcst3d,fcst2d)
  DEALLOCATE(anal3d,anal2d)

  CALL TIMER_END
    
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi

END PROGRAM letkf
