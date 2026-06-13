subroutine last_day(iyr,imon,iday)

  implicit none

  !---IN
  integer,intent(in) :: iyr,imon

  !---OUT
  integer,intent(out) :: iday

  if(imon == 2 .and. mod(iyr,imon) == 0)then
     iday=29
  else if(imon == 2)then
     iday=28
  else if(imon == 4 .or. imon == 6 .or. imon == 9 .or. imon == 11)then
     iday=30
  else
     iday=31
  end if  
  
end subroutine last_day
