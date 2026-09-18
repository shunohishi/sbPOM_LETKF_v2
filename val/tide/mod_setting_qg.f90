module setting
  
  integer,parameter :: ndat_a=4 !Number of analysis dataset

  real(kind = 8),parameter :: slon=0.d0,elon=360.d0   !Longitude
  real(kind = 8),parameter :: slat=-70.d0,elat=70.d0  !Latitude
  
  real(kind = 8),parameter :: dist_crit=100.d3 !Maximam distance criteria[m]  
  
end module setting

