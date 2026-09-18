subroutine check_data_location(nst,lon_a,lat_a,hdat_a,hsprd_a, &
     & lon_o,lat_o,dat_o,dist)

  use setting, only: slon,elon,slat,elat
  use mod_rmiss
  implicit none

  !---Common
  integer ist
  
  !---IN
  integer,intent(in) :: nst

  real(kind = 8),intent(in) :: lon_a(nst),lat_a(nst)
  real(kind = 8),intent(in) :: lon_o(nst),lat_o(nst)
  
  !---INOUT
  real(kind = 8),intent(inout) :: hdat_a(nst),hsprd_a(nst)
  real(kind = 8),intent(inout) :: dat_o(nst),dist(nst)

  do ist=1,nst

     if(lon_a(ist) < slon .or. elon < lon_a(ist) &
          & .or. lon_o(ist) < slon .or. elon < lon_o(ist) &
          & .or. lat_a(ist) < slat .or. elat < lat_a(ist) &
          & .or. lat_o(ist) < slat .or. elat < lat_o(ist))then
        hdat_a(ist)=rmiss
        hsprd_a(ist)=rmiss
        dat_o(ist)=rmiss
        dist(ist)=rmiss
     end if

  end do
  
end subroutine check_data_location
