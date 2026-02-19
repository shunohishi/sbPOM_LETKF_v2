module setting

  !---Observation
  integer,parameter :: nbuoy=2,nvar=4

  real(kind = 8),parameter :: obs_rate=20.d0 !Observation rate for depth average
  
  character(10),dimension(nbuoy),parameter :: buoyname=(/"keo","papa"/)
  character(1),dimension(nvar),parameter :: varname=(/"t","s","u","v"/)

  logical,parameter :: lwrite_obs=.true.
  
  !---Analysis information
  !The number of analyses (1: LORA, 2: GLORYS, 3: ORAS5, 4: CGLORS) => See mod_io.f90
  integer,parameter :: ndat_a=4
  character(10),dimension(ndat_a),parameter :: datname=(/"lora","glorys","oras5","cglors"/)
  
end module setting
