!----------------------------------------------------
! Distance between (x1,y1) and (x2,y2)
!----------------------------------------------------
!
! x1,x2: longitude (degree)
! y1,y2: latitude (degree)
!---------------------------------------------------

subroutine distance(x1,y1,x2,y2,r)
  
  use mod_parameter
  implicit none
  
  real(kind = 8),intent(in) :: x1,y1
  real(kind = 8),intent(in) :: x2,y2
  real(kind = 8),intent(out) :: r

  if(x1 == x2 .and. y1 == y2)then
     r=0.d0
  else
     r=earth*acos( &
          & cos(pi*y1/180.d0)*cos(pi*y2/180.d0)*cos(pi*x1/180.d0)*cos(pi*x2/180.d0) &
          & +cos(pi*y1/180.d0)*cos(pi*y2/180.d0)*sin(pi*x1/180.d0)*sin(pi*x2/180.d0) &
          & +sin(pi*y1/180.d0)*sin(pi*y2/180.d0) &
          & )
  end if
  
end subroutine distance
