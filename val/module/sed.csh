#!/bin/csh
#---------------------------------------------------------------------
# Change "DATA" directory pass from dirname1 to dirname2
#---------------------------------------------------------------------
set dirname1="/data/R/R2402/DATA"          #JSS3
set dirname2="/lvs0/rccs-dart/ohishi/DATA" #RCCS-Cloud
#---------------------------------------------------------------------

foreach filename(mod_read_db.f90 mod_read_duacs.f90 mod_read_glorys025.f90 mod_read_ocs.f90 mod_read_tide.f90)
    if(-f ${filename})then
	echo ${filename}
	sed -e "s|${dirname1}|${dirname2}|" ${filename} > tmp.f90
	mv tmp.f90 ${filename}
    endif
end

#---------------------------------------------------------------------
# Change "LETKF" directory from dirname1 to dirname2
#---------------------------------------------------------------------
set dirname1="/data/R/R2402/ohishi/JSS3"   #JSS3
set dirname2="/lvs0/rccs-dart/ohishi/JSS3" #RCCS-Cloud
#---------------------------------------------------------------------

foreach filename(mod_read_lora.f90 mod_read_lora_v20.f90)
    if(-f ${filename})then
	echo ${filename}
	sed -e "s|${dirname1}|${dirname2}|" ${filename} > tmp.f90
	mv tmp.f90 ${filename}
    endif
end
