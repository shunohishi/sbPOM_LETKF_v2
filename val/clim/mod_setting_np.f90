module setting

  !---Date
  integer,parameter :: syr=2003,eyr=2020 !Start/End year
  !integer,parameter :: syr=2003,eyr=2003 !Start/End year

  !---Zonal average (Meridional Section)
  real(kind = 8),parameter :: slon=110.d0,elon=250.d0  
  
  !---Analysis
  integer,parameter :: ndat=5 !The number of analysis dataset
  character(10),dimension(ndat),parameter :: datname=(/"lora      ","bran      ","fora      ","glorys    ","jcope     "/)
  
end module setting
