subroutine make_depth(depth)

  use setting, km => km_out
  implicit none

  real(kind = 8),intent(out) :: depth(km)

  depth(1)=0.e0
  depth(2)=50.e0
  depth(3)=100.e0
  depth(4)=200.e0
  depth(5)=400.e0

end subroutine make_depth
