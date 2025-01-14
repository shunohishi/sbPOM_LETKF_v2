!______________________________________________________________________

INTEGER FUNCTION JULDAY(IYYY,MONTH,DD)

  ! ********** SUBROUTINE DESCRIPTION:
  !
  ! FINDS THE JULIAN DAY FROM A DATE.
  !
  ! ********** ORIGINAL AUTHOR AND DATE:
  !
  ! PRESS,FLANNERY,TEUKOLSKY,VETTERLING 1986.
  ! NUMERICAL RECIPES
  !
  ! ********** REVISION HISTORY:
  !
  !
  ! ********** ARGUMENT DEFINITIONS:

  implicit none

  INTEGER,intent(in):: IYYY,MONTH,DD


  ! NAME   IN/OUT DESCRIPTION
  !
  ! IYYY     I    YEAR
  ! MONTH    I    MONTH (1 TO 12)
  ! DD       I    DAY OF MONTH
  ! JULDAY   O    JULIAN DAY
  !
  ! ********** COMMON BLOCKS:
  !
  ! NONE
  !
  ! ********** LOCAL PARAMETER DEFINITIONS:


  INTEGER,PARAMETER :: IGREG = 15 + 31*(10 + 12*1582)
  INTEGER iYYYi

  !
  ! ********** LOCAL VARIABLE DEFINITIONS:
  !

  INTEGER JY,JM,JA

  ! NAME   DESCRIPTION
  !
  !
  ! ********** OTHER ROUTINES AND FUNCTIONS CALLED:
  !
  ! INT    - INTRINSIC TRUNCATE
  !
  !---+67--1----+----2----+----3----+----4----+----5----+----6----+----7--


  iYYYi = iYYY
  IF (IYYY < 0) IYYYi = IYYY + 1
  IF(MONTH > 2)THEN
     JY = IYYYi
     JM = MONTH + 1
  ELSE
     JY = IYYYi - 1
     JM = MONTH + 13
  ENDIF
  JULDAY = INT(365.25d0*JY) + INT(30.6001d0*JM) + DD + 1720995
  IF (DD + 31*(MONTH + 12*IYYYi) >= IGREG) THEN
     JA = INT(0.01d0*JY)
     JULDAY = JULDAY + 2 - JA + INT(0.25d0*JA)
  ENDIF

END FUNCTION JULDAY

!__________________________________________________________________________

SUBROUTINE CALDAT(JULIAN,IYYY,MONTH,DD)

  ! ********** SUBROUTINE DESCRIPTION:
  !
  ! GIVEN THE JULIAN DAY, RETURNS THE YEAR, MONTH AND DAY OF MONTH.
  !
  ! ********** ORIGINAL AUTHOR AND DATE:
  !
  ! PRESS,FLANNERY,TEUKOLSKY,VETTERLING 1986.
  ! NUMERICAL RECIPES
  !
  ! ********** REVISION HISTORY:
  !
  !
  ! ********** ARGUMENT DEFINITIONS:
  !

  implicit none

  INTEGER,intent(in ) :: JULIAN
  INTEGER,intent(out) :: IYYY,MONTH,DD


  ! NAME   IN/OUT DESCRIPTION
  !
  ! JULIAN   I    THE JULIAN DAY
  ! IYYY     O    THE YEAR
  ! MONTH    O    THE MONTH (1 TO 12)
  ! DD       O    THE DAY OF THE MONTH
  !
  ! ********** COMMON BLOCKS:
  !
  ! NONE
  !
  ! ********** LOCAL PARAMETER DEFINITIONS:


  INTEGER,PARAMETER :: IGREG=2299161

  !
  ! ********** LOCAL VARIABLE DEFINITIONS:
  !

  INTEGER JALPHA,JA,JB,JC,JD,JE

  !
  ! NAME   DESCRIPTION
  !
  !
  ! ********** OTHER ROUTINES AND FUNCTIONS CALLED:
  !
  !
  !---+67--1----+----2----+----3----+----4----+----5----+----6----+----7--


  IF (JULIAN >= IGREG) THEN
     JALPHA = INT(((JULIAN - 1867216) - 0.25d0)/36524.25d0)
     JA = JULIAN + 1 + JALPHA - INT(0.25d0*JALPHA)
  ELSE
     JA = JULIAN
  ENDIF
  JB = JA + 1524
  JC = INT(6680.d0 + ((JB - 2439870) - 122.1d0)/365.25d0)
  JD = 365*JC + INT(0.25d0*JC)
  JE = INT((JB - JD)/30.6001d0)
  DD = JB - JD - INT(30.6001d0*JE)
  MONTH = JE - 1
  IF (MONTH > 12) MONTH = MONTH - 12
  IYYY = JC - 4715
  IF (MONTH > 2) IYYY = IYYY - 1
  IF (IYYY <= 0) IYYY = IYYY - 1

END SUBROUTINE CALDAT
!_______________________________________________________________________

!-----------------------------------------------------------------
!     Year,Month,Day --> Julian Day
!     (Jean Meeus)
!-----------------------------------------------------------------

subroutine ymd_julian(iyr,imon,iday,julian)

  implicit none

  integer,parameter :: igreg=15+31*(10+12*1582)

  integer iyr,imon,iday
  integer jyr,jmon,ja
  integer julian

  if(iyr < 0)then
     iyr=iyr+1
  end if

  if(imon > 2)then
     jyr=iyr
     jmon=imon+1
  else
     jyr=iyr-1
     jmon=imon+13
  end if

  julian=int(365.25d0*jyr)+int(30.6001d0*jmon)+iday+1720995

  if(iday+31*(imon+12*iyr) >= igreg)then
     ja=int(0.01d0*jyr)
     julian=julian+2-ja+int(0.25d0*ja)
  end if

end subroutine ymd_julian

!------------------------------------------------------------------
!     Julian Day --> Year,Month,Day
!------------------------------------------------------------------

subroutine julian_ymd(julian,iyr,imon,iday)

  implicit none

  integer,parameter :: igreg=2299161

  integer julian
  integer iyr,imon,iday

  integer jalpha,ja,jb,jc,jd,je

  if(julian >= igreg)then
     jalpha=int(((julian-1867216)-0.25d0)/36524.25d0)
     ja=julian+1+jalpha-int(0.25d0*jalpha)
  else
     ja=julian
  end if

  jb=ja+1524
  jc=int(6680.d0 + ((jb-2439870)-122.1d0)/365.25d0)
  jd=365*jc+int(0.25d0*jc)
  je=int((jb-jd)/30.6001d0)

  iday=jb-jd-int(30.6001d0*je)

  imon=je-1
  if(imon > 12)then
     imon=imon-12
  endif

  iyr=jc-4715
  if(imon > 2)then
     iyr=iyr-1
  end if
  if(iyr <= 0)then
     iyr=iyr-1
  end if

end subroutine julian_ymd

!-----------------------------------------------------------

subroutine readtiming(timing,ratio,dayint_data,time0_data)
  ! determin read timing of forcing data

  use common_pom_var
  implicit none

  ! arguments
  real(kind = r_dble),intent(in) :: dayint_data,time0_data

  integer,intent(out) :: timing
  real(kind = r_dble),intent(out) :: ratio

  ! local
  real(kind = r_dble) dev,dif,timepre

  timepre=dti*dble(iint-1)/86400.d0+time0
  dev=dble(timepre-time0_data)/dble(dayint_data)
  dif=dev-nint(dev)

  timing=0
  if(abs(dif) < dti/(86400.d0*dayint_data)) then
     timing=1
  end if

  ratio=dif
  if(abs(dif) < dti/(86400.d0*dayint_data)) then
     ratio=int(dif)
  else
     if(ratio < 0.d0) ratio=1.d0+ratio
  end if

end subroutine readtiming
!_______________________________________________________________________
subroutine ratioclim(ratio,julday_time,time)

  !     julday,time: [day]
  use common_pom_var,only: r_dble
  implicit none

  ! arguments
  integer,intent(in) :: julday_time
  real(kind = r_dble),intent(in) :: time
  real(kind = r_dble),intent(out) :: ratio

  ! local
  integer iyy0,iyy1,imm0,imm1,idd0
  integer julday,julday0,julday1

  call caldat(julday_time,iyy0,imm0,idd0)

  if(idd0 < 15) then
     imm0=imm0-1
     if(imm0 < 1) then
        imm0=12
        iyy0=iyy0-1
     end if
  end if

  imm1=imm0+1
  iyy1=iyy0

  if(imm1 > 12) then
     imm1=1
     iyy1=iyy0+1
  end if

  julday0=julday(iyy0,imm0,15)
  julday1=julday(iyy1,imm1,15)

  ratio=(julday_time+time-int(time)-julday0)/(julday1-julday0)

end subroutine ratioclim

