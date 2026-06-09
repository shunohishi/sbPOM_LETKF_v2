subroutine make_depth(depth)

  use setting, km => km_out
  implicit none

  real(kind = 8),intent(out) :: depth(km)

  depth(:)=0.d0
  
  depth(1)=0.d0
  depth(2)=50.d0
  depth(3)=100.d0
  depth(4)=200.d0
  depth(5)=400.d0

end subroutine make_depth
