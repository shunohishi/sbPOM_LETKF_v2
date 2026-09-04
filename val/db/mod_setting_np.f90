module setting

  !*** To be modified ***
  !---Analysis domain
  real(kind = 8),parameter :: slon=105.4d0,elon=255.d0 !Longitude range
  real(kind = 8),parameter :: slat=10.d0,elat=62.8d0   !Latitude range
  
  !---Analysis dataset
  integer,parameter :: ndat_a=4                     !The number of analysis datasets (1:LORA, 2: BRAN, 3:GLORYS, 4:JCOPE) => See mod_io.f90
  character(10),dimension(ndat_a),parameter :: &
       & datname=(/"lora      ","bran      ","glorys    ","jcope     "/) !Name of analysis datasets => Output filename
  
  !---Box information
  real(kind = 8),parameter :: slon_bin=110.d0,elon_bin=250.d0 !Longitude range
  real(kind = 8),parameter :: slat_bin=15.d0,elat_bin=55.d0   !Latitude range
  real(kind = 8),parameter :: dx_bin=5.0d0,dy_bin=5.0d0       !Resolution
    
end module setting
