subroutine check_data_location(nobs,lono,lato,hdat)

  use mod_rmiss
  use setting, only: slon, elon, slat, elat
  implicit none

  !---Common
  integer iobs
  
  !---IN
  integer,intent(in) :: nobs

  real(kind = 8),intent(in) :: lono(nobs),lato(nobs)

  !---IN/OUT
  real(kind = 8),intent(inout) :: hdat(nobs)

  do iobs=1,nobs
     if(lono(iobs) < slon .or. elon < lono(iobs) .or. &
          & lato(iobs) < slat .or. elat < lato(iobs))then
        hdat(iobs)=rmiss
     end if
  end do  
  
end subroutine check_data_location
