#!/bin/csh
#============================================================
# ARGV |
#============================================================

set REGION=${argv[1]}
set EXE=${argv[2]}
set DIR=${argv[3]}
set NCDIR=${argv[4]}
set ATMDIR=${argv[5]}
set RIVDIR=${argv[6]}
set NMEM=${argv[7]}

set WORKDIR=${argv[8]}
set mm_s=${argv[9]}

#============================================================
cd ${WORKDIR}/
@ IMEM = 1
while(${IMEM} <= ${NMEM})

    set MMMMM=`printf "%05d" ${IMEM}`
    if(${IMEM} == 1) echo "Link Netcdf file"
    
    #---grid.nc
    if(-f ${NCDIR}/grid.nc)then
	ln -s ${NCDIR}/grid.nc ${WORKDIR}/${MMMMM}/in/${REGION}.grid.nc
    else
	echo "***Error: Not found ${NCDIR}/grid.nc"
	exit
    endif
	    
    #---tsclim.nc
    if(-f ${NCDIR}/tsclim.nc)then
	ln -s ${NCDIR}/tsclim.nc ${WORKDIR}/${MMMMM}/in/${REGION}.tsclim.nc
    else
	echo "***Error: Not found ${NCDIR}/tsclim.nc"
	exit
    endif
	    
    #---ic.nc
    if(-f ${NCDIR}/ic.woa18.${mm_s}.nc)then
	ln -s ${NCDIR}/ic.woa18.${mm_s}.nc ${WORKDIR}/${MMMMM}/in/${REGION}.ic.nc
    else
	echo "***Error: Not found ${NCDIR}/tsclim.nc"
	exit
    endif

    #---lbc/tsdata.nc
    if(-f ${NCDIR}/tsdata_mclim.${MMMMM}.nc)then
	ln -s ${NCDIR}/tsdata_mclim.${MMMMM}.nc ${WORKDIR}/${MMMMM}/in/${REGION}.tsdata.nc
    else
	echo "***Error: Not found ${NCDIR}/tsdata_mclim.${MMMMM}.nc"
	exit
    endif

    if(-f ${NCDIR}/lbc_mclim.${MMMMM}.nc)then
	ln -s ${NCDIR}/lbc_mclim.${MMMMM}.nc ${WORKDIR}/${MMMMM}/in/${REGION}.lbc.nc
    else
	echo "***Error: Not found ${NCDIR}/lbc_mclim.${MMMMM}.nc"
	exit
    endif
			
    #---atm.nc
    if(-d ${ATMDIR})then
	ln -s ${ATMDIR} ${WORKDIR}/${MMMMM}/in/atm
    else
	echo "***Error: Not found ${ATMDIR}"
	exit
    endif

    #---atm_clim.nc
    if(-f ${NCDIR}/atm_clim.nc)then
	ln -s ${NCDIR}/atm_clim.nc ${WORKDIR}/${MMMMM}/in/${REGION}.atm_clim.nc
    else
	echo "***Error: Not found ${NCDIR}/atm_clim.nc"
	exit
    endif

    
    #---river.nc
    if(-d ${RIVDIR})then
	ln -s ${RIVDIR} ${WORKDIR}/${MMMMM}/in/river
    else
	echo "***Error: Not found ${RIVDIR}"
	exit
    endif

    #---exe
    if(-f ${DIR}/run/${EXE})then
	ln -s ${DIR}/run/${EXE} ${WORKDIR}/${MMMMM}/${EXE}
    else
	echo "***Error: Not found ${DIR}/run/${EXE}"
	exit
    endif

    @ IMEM++

end
