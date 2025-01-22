#!/bin/csh
#-----------------------------------------------------------------
# Prepare Ensemble member 
#-----------------------------------------------------------------

set IT=${argv[1]}
set NMEM=${argv[2]}
set MODELOUTPUTDIR=${argv[3]}
set OUTPUT=${argv[4]}

#-----------------------------------------------------------------

@ IT0 = ${IT} - 1
set date0=`perl juldays.prl ${IT0}`
set yyyy0=`printf "%04d" ${date0[1]}`; set mm0=`printf "%02d" ${date0[2]}`; set dd0=`printf "%02d" ${date0[3]}`

#---Copy restart file
@ IMEM = 1
while(${IMEM} <= ${NMEM})

    set MMMMM=`printf "%05d" ${IMEM}`

    set filename1=${OUTPUT}/${MMMMM}/restart.${yyyy0}${mm0}${dd0}.nc
    set filename2=${MODELOUTPUTDIR}/${yyyy0}${mm0}/${MMMMM}/out/restart.nc
    #set filename2=${MODELOUTPUTDIR}/${yyyy0}${mm0}/${MMMMM}/out/ens.${dd0}.nc
    set filename=${OUTPUT}/${MMMMM}/restart.${yyyy0}${mm0}${dd0}.nc
    if(! -f ${filename1}  && ! -f ${filename2})then
	echo "***Error: Not found ${filename1} and ${filename2}"
	exit 99
    else if(! -f ${filename1} && -f ${filename2})then
	(cp ${filename2} ${filename} &)
    endif

    @ IMEM++

end
wait

#---Check restart file
@ IMEM = 1
while(${IMEM} <= ${NMEM})

    set MMMMM=`printf "%05d" ${IMEM}`

    set filename=${OUTPUT}/${MMMMM}/restart.${yyyy0}${mm0}${dd0}.nc
    if(! -f ${filename})then
	echo "***Error: ${filename} not found"
	exit 99
    endif

    @ IMEM++
    
end

wait
exit 0
