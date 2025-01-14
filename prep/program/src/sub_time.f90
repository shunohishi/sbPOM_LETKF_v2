!-----------------------------------------------------------------
!     Year,Month,Day --> Julian Day
!     (Jean Meeus)
!-----------------------------------------------------------------

subroutine ymd_julian(iyr,imon,iday,julian)

  implicit none

  integer,parameter :: igreg=15+31*(10+12*1582)

  integer jyr,jmon,ja
  
  integer,intent(in) :: iyr,imon,iday
  integer,intent(out) :: julian

  if(imon > 2)then
     jyr=iyr
     jmon=imon+1
  else
     jyr=iyr-1
     jmon=imon+13
  end if

  julian=int(365.25*jyr)+int(30.6001*jmon)+iday+1720995

  if(iday+31*(imon+12*iyr) >= igreg)then
     ja=int(0.01d0*jyr)
     julian=julian+2-ja+int(0.25d0*ja)
  end if

end subroutine ymd_julian

!------------------------------------------------------------------
!     Julian Day --> Year,Month,Day
!------------------------------------------------------------------

subroutine julian_ymd(ijul,iyr,imon,iday)

  implicit none

  integer,parameter :: igreg=2299161

  integer jalpha,ja,jb,jc,jd,je
  
  integer,intent(in) :: ijul
  integer,intent(out) :: iyr,imon,iday


  if(ijul >= igreg)then
     jalpha=int(((ijul-1867216)-0.25)/36524.25)
     ja=ijul+1+jalpha-int(0.25*jalpha)
  else
     ja=ijul
  end if

  jb=ja+1524
  jc=int(6680. + ((jb-2439870)-122.1)/365.25)
  jd=365*jc+int(0.25*jc)
  je=int((jb-jd)/30.6001)

  iday=jb-jd-int(30.6001*je)

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

!-----------------------------------------------------------------------------
! Estimate ntime |
!-----------------------------------------------------------------------------

  subroutine estimate_ntime(syr,smon,sday,shour,eyr,emon,eday,ehour,dt,ntime)

    implicit none

    integer sjul,ejul !Julian Day
    
    integer,intent(in) :: syr,smon,sday,shour
    integer,intent(in) :: eyr,emon,eday,ehour
    integer,intent(in) :: dt
    integer,intent(out) :: ntime

    !Calculate Julian Day
    call ymd_julian(syr,smon,sday,sjul)
    call ymd_julian(eyr,emon,eday,ejul)
    
    ntime=(ejul-sjul)*24+(ehour-shour)
    ntime=ntime/dt+1

  end subroutine estimate_ntime

!---------------------------------------------------------------
! Add 1 time step
!---------------------------------------------------------------

subroutine add_time(iyr,imon,iday,ihour,dt)

  implicit none

  integer ijul
  real(kind = 8) rjul
  
  integer,intent(inout) :: iyr,imon,iday,ihour
  integer,intent(in) :: dt


  call ymd_julian(iyr,imon,iday,ijul)
  rjul=dble(ijul)+dble(ihour+dt)/24.d0

  !Estimate iyr,imon,iday
  call julian_ymd(int(rjul),iyr,imon,iday)
  !Estimate ihour
  ihour=nint((rjul-dble(int(rjul)))*24.d0)

end subroutine add_time

!---------------------------------------------------------------
! Number of day |
!---------------------------------------------------------------

subroutine number_of_day(iyr,imon,nday)

  implicit none
  
  integer,intent(in) :: iyr,imon
  integer,intent(out) :: nday

  if(imon==2 .and. mod(iyr,4)==0)then
     nday=29
  elseif(imon==2 .and. mod(iyr,4)/=0)then 
     nday=28
  elseif(imon==4 .or. imon==6 .or. imon==9 .or. imon==11)then
     nday=30
  elseif(imon==1 .or. imon==3 .or. imon==5 .or. imon==7 .or. imon==8 &
       & .or. imon==10 .or. imon == 12)then
     nday=31
  else
     write(*,'(a)') "***Error: Incorrect year or month"
     stop
  end if

end subroutine number_of_day
