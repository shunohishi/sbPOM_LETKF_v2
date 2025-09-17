module setting

  !---Box information
  real(kind = 8),parameter :: slon_bin=0.d0,elon_bin=360.d0
  real(kind = 8),parameter :: slat_bin=-90.d0,elat_bin=90.d0
  real(kind = 8),parameter :: dx_bin=5.0d0,dy_bin=5.0d0
  
  !---Analysis information
  integer,parameter :: ndat_a=4 !The number of analyses (1: LORA, 2: GLORYS, 3: ORAS5, 4: CGLORS) => See mod_io.f90
  
end module setting
