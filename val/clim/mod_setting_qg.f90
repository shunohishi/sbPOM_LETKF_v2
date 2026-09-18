module setting

  !---Date
  integer,parameter :: syr=2003,eyr=2023 !Start/End year

  !---Zonal average (Meridional Section)
  real(kind = 8),parameter :: slon=60.d0,elon=360.d0  
  
  !---Analysis
  integer,parameter :: ndat=4 !The number of analysis dataset
  character(10),dimension(ndat),parameter :: datname=(/"lora      ","glorys    ","oras5     ","cglors    "/)
  
end module setting
