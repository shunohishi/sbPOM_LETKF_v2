#!/bin/csh
#=== Argument ================================================================#
set switch_iau=${argv[1]} #0: Intermittent, 1: 1st simulation, 2: 2nd simulation for IAU
set WORKDIR=${argv[2]}
set OUTPUT=${argv[3]}
set INFO=${argv[4]}
set machine=${argv[5]}

set MODELDATADIR=${argv[6]}
set NMEM=${argv[7]}

set REGION=${argv[8]}
set EXE=${argv[9]}

set yyyymmdd=${argv[10]}
set yyyymmdd_n1=${argv[11]}
set dd = `echo "${yyyymmdd}" | awk '{print substr($1, length($1)-1)}'`

#=============================================================================#
echo "Move restart and output files and log"

#Restart file(mean & sprd)
if(${dd} != "01")then
    (rm -f ${OUTPUT}/mean/fcst.${yyyymmdd_n1}.nc)
    (rm -f ${OUTPUT}/sprd/fcst.${yyyymmdd_n1}.nc)
    (rm -f ${OUTPUT}/mean/anal.${yyyymmdd_n1}.nc)
    (rm -f ${OUTPUT}/sprd/anal.${yyyymmdd_n1}.nc)
endif

@ IMEM = 1
while(${IMEM} <= ${NMEM})

    set MMMMM=`printf "%05d" ${IMEM}`

    #---Restart file
    set filename1=${WORKDIR}/${MMMMM}/out/restart.nc
    set filename2=${OUTPUT}/${MMMMM}/restart.${yyyymmdd}.nc
    if(${switch_iau} != 1 && ! -f ${filename1})then
    	echo "***Error: ${filename1} not found" && exit 99
    else if(${switch_iau} != 1)then
	(mv -f ${filename1} ${filename2})
    endif

    #---IAU file
    set filename1=${WORKDIR}/${MMMMM}/out/iau.nc
    set filename2=${OUTPUT}/${MMMMM}/iau.${yyyymmdd}.nc
    if(${switch_iau} == 1 && ! -f ${filename1})then
	echo "***Error: ${filename1} not found" && exit 99
    else if(${switch_iau} == 1)then
	(mv -f ${filename1} ${filename2})
    else if(${switch_iau} == 2)then
	(rm -f ${filename2})
    endif
    
    #---Output file
    set filename1=${WORKDIR}/${MMMMM}/out/${REGION}.nc
    set filename2=${OUTPUT}/${MMMMM}/${REGION}.${yyyymmdd}.nc
    if(${switch_iau} != 1 && ! -f ${filename1})then
	echo "***Error: ${filename1} not found" && exit 99
    else if(${switch_iau} != 1)then
	(mv -f ${filename1} ${filename2})
    endif
	
    #---NML
    (rm -f ${WORKDIR}/${MMMMM}/pom.nml)
	
    #---stout
    if(${switch_iau} == 1)then
	set filename=${INFO}/stdout.${REGION}.${yyyymmdd}_iau.${MMMMM}
    else
	set filename=${INFO}/stdout.${REGION}.${yyyymmdd}.${MMMMM}
    endif
    
    if(${machine} == "jss3" && ${IMEM} == 1 && -s ${WORKDIR}/${MMMMM}/stdout.${REGION})then
	(mv -f ${WORKDIR}/${MMMMM}/stdout.${REGION} ${filename})
    else if(${machine} == "fugaku" && ${IMEM} == 1)then
	foreach file(${WORKDIR}/${MMMMM}/stdout.${REGION}.*.0)
	    if(-s ${file})then
		(mv -f ${file} ${filename})
		break
	    endif
	end
    endif

    #---stderr
    if(${switch_iau} == 1)then
	set filename=${INFO}/stderr.${REGION}.${yyyymmdd}_iau.${MMMMM}
    else
	set filename=${INFO}/stderr.${REGION}.${yyyymmdd}.${MMMMM}
    endif

    if(${machine} == "jss3" && ${IMEM} == 1 && -s ${WORKDIR}/${MMMMM}/stderr.${REGION})then
	(mv -f ${WORKDIR}/${MMMMM}/stderr.${REGION} ${filename})
    else if(${machine} == "fugaku" && ${IMEM} == 1)then
	#foreach file(${WORKDIR}/${MMMMM}/stderr.${REGION}.*.0)
	#    if(-s ${file})then
	#	(mv -f ${file} ${filename})
	#	break
	#    endif
	#end
    endif

    #---Remove link
    (unlink ${WORKDIR}/${MMMMM}/pom.exe)

    foreach filename(grid tsclim ic tsdata lbc)
	unlink ${WORKDIR}/${MMMMM}/in/${REGION}.${filename}.nc
    end
	
    foreach dirname(atm river)
	unlink ${WORKDIR}/${MMMMM}/in/${dirname}
    end

    unlink ${WORKDIR}/${MMMMM}/in/restart.nc
    
    @ IMEM++

end
wait
exit 0

