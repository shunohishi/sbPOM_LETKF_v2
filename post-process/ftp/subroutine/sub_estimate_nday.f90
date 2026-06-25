subroutine estimate_nday(iyr,imon,nday)

  implicit none

  integer,intent(in)  :: iyr,imon
  integer,intent(out) :: nday

  integer,parameter :: mdays(12) = &
       (/31,28,31,30,31,30,31,31,30,31,30,31/)

  nday = mdays(imon)

  if(imon == 2) then
     if( (mod(iyr,4) == 0 .and. mod(iyr,100) /= 0) .or. &
         (mod(iyr,400) == 0) ) then
        nday = 29
     endif
  endif

end subroutine estimate_nday
