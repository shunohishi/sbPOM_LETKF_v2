module setting

  use mod_gridinfo, im_in => im, jm_in => jm, km_in => km
  implicit none
  
  !---Date
  integer,parameter :: syr=2003,smon=1,sday=1
  integer,parameter :: eyr=2003,emon=1,eday=31

  !---Dataset
  integer,parameter :: im_out=im_in,jm_out=jm_in,km_out=5 !Output grid size
  integer,parameter :: nvar2d=1,nvar3d=4
    
  character(10),parameter :: dir="QGLOBAL",region="qglobal"
  character(10),parameter :: letkf="letkf"
  character(100),parameter :: title="LETKF-based Ocean Reserach Analysis (LORA) version 2.0"

end module setting
