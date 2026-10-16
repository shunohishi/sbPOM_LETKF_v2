#!/bin/csh
#---------------------------------------------------------------------------------#
# ARGV |
#---------------------------------------------------------------------------------#

set DIR=${argv[1]}
set NCDIR=${argv[2]}
set ATMDIR=${argv[3]}
set RIVDIR=${argv[4]}
set WORKDIR=${argv[5]}
set REGION=${argv[6]}
set smm=${argv[7]}
set EXE=${argv[8]}

#---------------------------------------------------------------------------------#

echo "Link Netcdf file"
#---grid.nc
if(-f ${NCDIR}/grid.nc)then
    ln -s ${NCDIR}/grid.nc ${WORKDIR}/in/${REGION}.grid.nc
else
    echo "***Error: Not found ${NCDIR}/grid.nc"
    exit 1
endif
	
#---tsclim.nc
if(-f ${NCDIR}/tsclim.nc)then
    ln -s ${NCDIR}/tsclim.nc ${WORKDIR}/in/${REGION}.tsclim.nc
else
    echo "***Error: Not found ${NCDIR}/tsclim.nc"
    exit 1
endif
	
#---ic.nc
if(-f ${NCDIR}/ic.woa.${smm}.nc)then
    ln -s ${NCDIR}/ic.woa.${smm}.nc ${WORKDIR}/in/${REGION}.ic.nc
else
    echo "***Error: Not found ${NCDIR}/ic.woa.${smm}.nc"
    exit 1
endif

#---tsdata.nc
if(-f ${NCDIR}/tsdata_mclim.nc)then
    ln -s ${NCDIR}/tsdata_mclim.nc ${WORKDIR}/in/${REGION}.tsdata.nc
else
    echo "***Error: Not found ${NCDIR}/tsdata_mclim.nc"
    exit 1
endif

#---lbc.nc
if(-f ${NCDIR}/lbc_mclim.nc)then
    ln -s ${NCDIR}/lbc_mclim.nc ${WORKDIR}/in/${REGION}.lbc.nc
else
    echo "***Error: Not found ${NCDIR}/lbc_mclim.nc"
    exit 1
endif

#---atm dir
if(-d ${ATMDIR})then
    ln -s ${ATMDIR} ${WORKDIR}/in/atm
else
    echo "***Error: Not found ${ATMDIR}"
    exit 1
endif

#---riv dir
if(-d ${RIVDIR})then
    ln -s ${RIVDIR} ${WORKDIR}/in/river
else
    echo "***Error: Not found ${RIVDIR}"
    exit 1
endif
    
#---EXE
if(-f ${DIR}/run/${EXE})then
    ln -s ${DIR}/run/${EXE} ${WORKDIR}/${EXE}
else
    echo "***Error: Not found ${DIR}/run/${EXE}"
    exit 1
endif
    
