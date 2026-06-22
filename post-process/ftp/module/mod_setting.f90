module setting

  integer,parameter :: syr=2002,smon=6,sday=1 !*sday=1
  integer,parameter :: eyr=2002,emon=6,eday=2

  integer,parameter :: nmem=128
  integer,parameter :: nvar2d=13,nvar3d=6,nvar_ens=5
  integer,parameter :: nvar_mlt=12,nvar_mls=12,nvar_mld=5
    
  real(kind = 8),parameter :: drho=0.125d0 !MLD criteria density [kg/m^3] 

  !----File information
  !INPUT
  character(100),parameter :: dir="QGLOBAL",letkf="letkf",region="qglobal"
  character(100),parameter :: dir_restart="/data/R/R2402/ohishi/QGLOBAL"
  
  !OUTPUT
  character(100),parameter :: filename_grid="dat/grid.nc"
  
end module setting
