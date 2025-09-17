subroutine number_of_day(iyr,imon,nday)

  implicit none

  integer,intent(in) :: iyr,imon
  integer,intent(out) :: nday

  if(imon == 2 .and. mod(iyr,4) == 0)then
     nday=29
  else if(imon == 2)then
     nday=28
  else if(imon == 4 .or. imon == 6 .or. imon == 9 .or. imon == 11)then
     nday=30
  else
     nday=31
  end if
     
end subroutine number_of_day
