module setting

  !---Analysis information (***To be modified ***)
  integer,parameter :: ndat_a=4   !The number of analysis datasets (1:LORA, 2:BRAN, 3:GLORYS, 4:JCOPE) => See mod_io.f90
  character(10),dimension(ndat_a),parameter :: datname=(/"lora      ","bran      ","glorys    ","jcope     "/) !Output filename

  !---Observation
  integer,parameter :: nbuoy=2,nvar=4

  real(kind = 8),parameter :: obs_rate=20.d0 !Observation rate for depth average [%]
  
  character(10),dimension(nbuoy),parameter :: buoyname=(/"keo       ","papa      "/)
  character(1),dimension(nvar),parameter :: varname=(/"t","s","u","v"/)

  logical,parameter :: lwrite_obs=.true.
    
end module setting
