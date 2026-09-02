module setting

  !*** To be modified ***
  !---Analysis domain
  real(kind = 8),parameter :: slon=0.d0,elon=360.d0  !Longitude range
  real(kind = 8),parameter :: slat=-70.d0,elat=70.d0 !Latitude range

  !---Analysis information
  integer,parameter :: ndat_a=4                       !The number of analyses (1:LORA, 2:GLORYS, 3:ORAS5, 4:CGLORS) => See mod_io.f90
  character(10),dimension(ndat_a),parameter :: &
       & datname=(/"lora","glorys","oras5","cglors"/) !Name of analysis datasets => Output filename

  !---Box information
  real(kind = 8),parameter :: slon_bin=0.d0,elon_bin=360.d0  !Longitude range
  real(kind = 8),parameter :: slat_bin=-70.d0,elat_bin=70.d0 !Latitude range
  real(kind = 8),parameter :: dx_bin=5.0d0,dy_bin=5.0d0      !Resolution
    
end module setting
