#!/bin/csh
#=== argv ======================================================#

set ST=${argv[1]}
set NMEM=${argv[2]}
set PDIR=${argv[3]}
set ATMDIR=${argv[4]}
set RIVDIR=${argv[5]}
set MODELDATADIR=${argv[6]}

#=== Date ======================================================#

#ST
set sdate=`perl juldays.prl ${ST}`
set yyyys=`printf "%04d" ${sdate[1]}`;set mms=`printf "%02d" ${sdate[2]}`;set dds=`printf "%02d" ${sdate[3]}`

#=== Check file ================================================#

set GRID=grid               #Model Grid data
set TSCLIM=tsclim           #T/S Climatology from WOA
set IC=ic.woa18.${mms}      #Initial Condition Netcdf file
set LBC=lbc_mclim           #Lateral Boundary Condition Netcdf file
set TSDATA=tsdata_mclim     # Nudging Netcdf file
set ATMCLIM=atm_clim

foreach file(${GRID}.nc ${TSCLIM}.nc ${IC}.nc)
    if(! -f ${MODELDATADIR}/${file})then
	echo "***Error: ${file} not found"
	exit 99
    endif
end #file

@ IMEM = 1
while(${IMEM} <= ${NMEM})
    set MMMMM=`printf "%05d" ${IMEM}`
    foreach file(${LBC}.${MMMMM}.nc ${TSDATA}.${MMMMM}.nc)
	if(! -f ${MODELDATADIR}/${file})then
	    echo "***Error: ${file} not found"
	    exit 99
	endif
    end #file
    @ IMEM++
end #IMEM

foreach dir(${ATMDIR} ${RIVDIR})
    if(! -d ${dir})then
	echo "***Error: ${dir} not found"
	exit 99
    endif
end

set file=${ATMDIR}/${ATMCLIM}.nc
if(! -f ${file})then
    echo "***Error: ${file} not found" 
endif

#============================================================#

wait
exit 0
