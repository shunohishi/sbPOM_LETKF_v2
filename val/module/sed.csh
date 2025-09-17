#!/bin/csh
#---------------------------------------------------------------------
# Change from dirname1 to dirname2
#---------------------------------------------------------------------
set dirname1="/vol0004/ra000007/data/a04048/DATA"
set dirname2="/data/R/R2402/DATA"
#---------------------------------------------------------------------

foreach filename(mod_bin.f90            mod_julian.f90         mod_read_glorys025.f90 mod_read_lora_v20.f90  mod_static.f90 mod_gridinfo.f90       mod_read_db.f90        mod_read_lora.f90      mod_rmiss.f90)
    if(-f ${filename})then
	echo ${filename}
	sed -e "s|${dirname1}|${dirname2}|" ${filename} > tmp.f90
	mv tmp.f90 ${filename}
    endif
end
