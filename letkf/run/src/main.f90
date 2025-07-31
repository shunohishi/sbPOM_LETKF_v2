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
  !   07/31/2025 Shun Ohishi        modified for LPFGM
  !=======================================================================
  
  USE MPI
  USE common_setting
  USE common_mpi
  USE common
  USE common_io
  USE common_pom
  USE common_obs
  USE letkf_tools

  IMPLICIT NONE

  INTEGER i,ierr
  INTEGER rtimerxx,rtimer00,rtimer,rate

  REAL(r_dble) ctimer00,ctimer

  REAL(r_size),ALLOCATABLE :: fcst3d(:,:,:,:),fcst2d(:,:,:)
  REAL(r_size),ALLOCATABLE :: anal3d(:,:,:,:),anal2d(:,:,:)
  
  CHARACTER(10) :: stdoutf='NOUT-00000'
  CHARACTER(4) :: fcstf='fc00'

  !-----------------------------------------------------------------------
  ! Initial settings
  !-----------------------------------------------------------------------

  CALL SYSTEM_CLOCK(rtimerxx,rate,i)
  CALL SYSTEM_CLOCK(rtimer00)
  CALL CPU_TIME(ctimer00)
  CALL initialize_mpi

  WRITE(stdoutf(6:10), '(I5.5)') myrank
 
  OPEN(newunit=file_unit,FILE=stdoutf,status="replace")  

  WRITE(file_unit,'(A)') "===Initial settings=========================="
  WRITE(file_unit,'(A,I5.5)') "STDOUT goes to "//trim(stdoutf)//" for MYRANK ", myrank
  WRITE(file_unit,'(A,I5.5,2A)') 'MYRANK=',myrank,', STDOUTF=',stdoutf

  WRITE(file_unit,'(A)') '============================================='
  WRITE(file_unit,'(A)') '  LOCAL ENSEMBLE TRANSFORM KALMAN FILTER     '
  WRITE(file_unit,'(A)') '                                             '
  WRITE(file_unit,'(A)') '   LL      EEEEEE  TTTTTT  KK  KK  FFFFFF    '
  WRITE(file_unit,'(A)') '   LL      EE        TT    KK KK   FF        '
  WRITE(file_unit,'(A)') '   LL      EEEEE     TT    KKK     FFFFF     '
  WRITE(file_unit,'(A)') '   LL      EE        TT    KK KK   FF        '
  WRITE(file_unit,'(A)') '   LLLLLL  EEEEEE    TT    KK  KK  FF        '
  WRITE(file_unit,'(A)') '                                             '
  WRITE(file_unit,'(A)') '             WITHOUT LOCAL PATCH             '
  WRITE(file_unit,'(A)') '                                             '
  WRITE(file_unit,'(A)') '          Coded by Takemasa Miyoshi          '
  WRITE(file_unit,'(A)') '  Based on Ott et al (2004) and Hunt (2005)  '
  WRITE(file_unit,'(A)') '  Tested by Miyoshi and Yamane (2006)        '
  WRITE(file_unit,'(A)') '============================================='
  WRITE(file_unit,'(A)') '              LETKF PARAMETERS               '
  WRITE(file_unit,'(A)') ' ------------------------------------------- '
  WRITE(file_unit,'(A,I15)')  '   nbv        :',nbv
  WRITE(file_unit,'(A,I15)')  '   nslots     :',nslots
  WRITE(file_unit,'(A,I15)')  '   nbslot     :',nbslot
  WRITE(file_unit,'(A,F15.2)')'   sigma_obs [km]  :',sigma_obs*1.d-3
  WRITE(file_unit,'(A,F15.2)')'   sigma_obsv [m]  :',sigma_obsv
  WRITE(file_unit,'(A,F15.2)')'   sigma_obst [day]:',sigma_obst
  WRITE(file_unit,'(A,F15.2)')'   MULT INF parameter:',cov_infl_mul
  WRITE(file_unit,'(A,F15.2)')'   RTPP parameter:',ALPHA_RTPP
  WRITE(file_unit,'(A,F15.2)')'   RTPS parameter:',ALPHA_RTPS  
  WRITE(file_unit,'(A)') '============================================='
  
  CALL set_common_pom
  CALL set_common_mpi_pom
  
  ALLOCATE(fcst3d(nij1,nlev,nbv,nv3d),fcst2d(nij1,nbv,nv2d))
  ALLOCATE(anal3d(nij1,nlev,nbv,nv3d),anal2d(nij1,nbv,nv2d))

  CALL SYSTEM_CLOCK(rtimer)  
  CALL CPU_TIME(ctimer)
  IF(myrank == 0)THEN
     WRITE(6,'(A)')       "=== TIMER: INITIAL ====================="
     WRITE(6,'(A,F10.2)') " TOTAL CPUTIME [s]:",ctimer
     WRITE(6,'(A,F10.2)') " SECTION CPUTIME [s]:",ctimer-ctimer00
     WRITE(6,'(A,F10.2)') " TOTAL REALTIME [s]:",dble(rtimer-rtimerxx)/dble(rate)
     WRITE(6,'(A,F10.2)') " SECTION REALTIME [s]:",dble(rtimer-rtimer00)/dble(rate)
     WRITE(6,'(A)')       "========================================"
  END IF
  rtimer00=rtimer
  ctimer00=ctimer
  
  !-----------------------------------------------------------------------
  ! Observations and Ensemble forecasts in obs. space
  !-----------------------------------------------------------------------

  WRITE(file_unit,'(A)') "===Obs. and Ensemble Forecast in obs. space=="
  CALL set_obs(fcst3d,fcst2d)
  
  CALL SYSTEM_CLOCK(rtimer)
  CALL CPU_TIME(ctimer)
  IF(myrank == 0)THEN
     WRITE(6,'(A)')       "=== TIMER: OBS and dYf ================="
     WRITE(6,'(A,F10.2)') " TOTAL CPUTIME [s]:",ctimer
     WRITE(6,'(A,F10.2)') " SECTION CPUTIME [s]:",ctimer-ctimer00
     WRITE(6,'(A,F10.2)') " TOTAL REALTIME [s]:",dble(rtimer-rtimerxx)/dble(rate)
     WRITE(6,'(A,F10.2)') " SECTION REALTIME [s]:",dble(rtimer-rtimer00)/dble(rate)
     WRITE(6,'(A)')       "========================================"
  END IF
  rtimer00=rtimer
  ctimer00=ctimer
  
  !-----------------------------------------------------------------------
  ! Analysis
  !-----------------------------------------------------------------------
  
  WRITE(file_unit,'(A)') "====================Analysis=================="
  CALL das_letkf(fcst3d,fcst2d,anal3d,anal2d)

  CALL SYSTEM_CLOCK(rtimer)
  CALL CPU_TIME(ctimer)
  IF(myrank == 0)THEN
     WRITE(6,'(A)')       "=== TIMER: Analysis ===================="
     WRITE(6,'(A,F10.2)') " TOTAL CPUTIME [s]:",ctimer
     WRITE(6,'(A,F10.2)') " SECTION CPUTIME [s]:",ctimer-ctimer00
     WRITE(6,'(A,F10.2)') " TOTAL REALTIME [s]:",dble(rtimer-rtimerxx)/dble(rate)
     WRITE(6,'(A,F10.2)') " SECTION REALTIME [s]:",dble(rtimer-rtimer00)/dble(rate)
     WRITE(6,'(A)')       "========================================"
  END IF
  rtimer00=rtimer
  ctimer00=ctimer
  
  !-----------------------------------------------------------------------
  ! Write data |
  !-------------
  ! - Forecast ensemble mean and spread
  ! - Analysis ensemble mean and spread
  ! - Ensemble analyses 
  !-----------------------------------------------------------------------
  
  WRITE(file_unit,'(A)') "=================Write DATA===================="
  CALL write_ens_mpi(fcst3d,fcst2d,anal3d,anal2d)

  CALL SYSTEM_CLOCK(rtimer)
  CALL CPU_TIME(ctimer)
  IF(myrank == 0)THEN
     WRITE(6,'(A)')       "=== TIMER: Write DATA =================="
     WRITE(6,'(A,F10.2)') " TOTAL CPUTIME [s]:",ctimer
     WRITE(6,'(A,F10.2)') " SECTION CPUTIME [s]:",ctimer-ctimer00
     WRITE(6,'(A,F10.2)') " TOTAL REALTIME [s]:",dble(rtimer-rtimerxx)/dble(rate)
     WRITE(6,'(A,F10.2)') " SECTION REALTIME [s]:",dble(rtimer-rtimer00)/dble(rate)
     WRITE(6,'(A)')       "========================================"
  END IF
  rtimer00=rtimer
  ctimer00=ctimer
    
  !-----------------------------------------------------------------------
  ! Finalize
  !-----------------------------------------------------------------------

  CLOSE(file_unit)
  DEALLOCATE(fcst3d,fcst2d)
  DEALLOCATE(anal3d,anal2d)

  IF(myrank == 0)THEN
     CALL SYSTEM_CLOCK(rtimer)  
     CALL CPU_TIME(ctimer)
     WRITE(6,'(A)') "=============== SUCCESS: LETKF ===================="
     WRITE(6,'(A,F10.2)') " TOTAL CPUTIME [s]:",ctimer
     WRITE(6,'(A,F10.2)') " TOTAL REALTIME [s]:",dble(rtimer-rtimerxx)/dble(rate)
     WRITE(6,'(A)') "==================================================="
  END IF
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi

END PROGRAM letkf
