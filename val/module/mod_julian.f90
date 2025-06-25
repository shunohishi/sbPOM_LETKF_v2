module mod_julian

contains

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

end module mod_julian
