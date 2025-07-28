#!/bin/csh
#---------------------------------------------------------------------
# Change from dirname1 to dirname2
#---------------------------------------------------------------------
set dirname1="/data/R/R2402/DATA"
set dirname2="/vol0004/ra000007/data/a04048/DATA"
#---------------------------------------------------------------------

foreach filename(mod_read_amsr2.f90 mod_read_chla.f90 mod_read_en4.f90 mod_read_etopo1.f90 mod_read_jra55do.f90 mod_read_woa18_month.f90 mod_read_amsre.f90 mod_read_gcomc.f90 mod_read_smap.f90 mod_read_woa18_season_annual.f90 mod_read_aqc_argo.f90 mod_read_gtspp.f90 mod_read_smos.f90 mod_rmiss.f90 mod_read_cama.f90 mod_read_himawari.f90 mod_read_soda.f90 mod_read_cmems.f90 mod_read_jra55.f90 mod_read_windsat.f90)
    if(-f ${filename})then
	echo ${filename}
	sed -e "s|${dirname1}|${dirname2}|" ${filename} > tmp.f90
	mv tmp.f90 ${filename}
    endif
end
