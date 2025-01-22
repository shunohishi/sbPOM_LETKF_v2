#!/bin/csh
#========= ARGUMENT ==================================================================
set switch_iau=${argv[1]} #0: Intermittent, 1: IAU 1st simulation, 2: IAU 2nd simulation 
set CDIR=${argv[2]}
set WORKDIR=${argv[3]}
set OUTPUT=${argv[4]}
set INFO=${argv[5]}
set NMEM=${argv[6]}
set MODELDIR=${argv[7]}
set MODELDATADIR=${argv[8]}
set ATMDIR=${argv[9]}
set RIVDIR=${argv[10]}
set REGION=${argv[11]}
set EPROC=${argv[12]}
set ENODE=${argv[13]}
set THREAD=${argv[14]}
set EXE=${argv[15]}
set BUDGET=${argv[16]}
set TS_NUDGE=${argv[17]}
set TI_NUDGE=${argv[18]}
set SS_NUDGE=${argv[19]}
set SI_NUDGE=${argv[20]}
set ST=${argv[21]}
set IT=${argv[22]}
set DT=${argv[23]}

#========= Date ======================================================================
set sdate=`perl juldays.prl ${ST}`
set yyyys=`printf "%04d" ${sdate[1]}`; set mms=`printf "%02d" ${sdate[2]}`; set dds=`printf "%02d" ${sdate[3]}`

@ BT = ${IT} - 1
set date_n1=`perl juldays.prl ${BT}`
set yyyy_n1=`printf "%04d" ${date_n1[1]}`; set mm_n1=`printf "%02d" ${date_n1[2]}`; set dd_n1=`printf "%02d" ${date_n1[3]}`
set yyyymmdd_n1=${yyyy_n1}${mm_n1}${dd_n1}

set date=`perl juldays.prl ${IT}`
set yyyy=`printf "%04d" ${date[1]}`; set mm=`printf "%02d" ${date[2]}`; set dd=`printf "%02d" ${date[3]}`
set yyyymmdd=${yyyy}${mm}${dd}

#=== MODEL Setting ===================================================================
echo "Start Model setting ..."
#---echo "Make in out DIR"
@ IMEM = 1
while(${IMEM} <= ${NMEM})

    set MMMMM=`printf "%05d" ${IMEM}`
    
    foreach dir(in out)
	(rm -rf ${WORKDIR}/${MMMMM}/${dir} && mkdir -p ${WORKDIR}/${MMMMM}/${dir} &)
    end #dir

    @ IMEM++

end #IMEM
wait

#---echo "Check restart & iau files"
@ IMEM = 1
while(${IMEM} <= ${NMEM})

    set MMMMM=`printf "%05d" ${IMEM}`

    #Restart file
    set filename=${OUTPUT}/${MMMMM}/restart.${yyyymmdd_n1}.nc
    if(! -f ${filename})then
	echo "***Error: ${filename} not found"; exit 99
    endif

    #IAU file
    set filename=${OUTPUT}/${MMMMM}/iau.${yyyymmdd}.nc
    if(${switch_iau} == 2 && ! -f ${filename})then
	echo "***Error: ${filename} not found"; exit 99
    endif

    @ IMEM++

end
wait

#---echo "Make Netcdf file & TIDE dir & EXE file link
csh make_sbPOM_link.csh \
${switch_iau} ${NMEM} \
${WORKDIR} ${OUTPUT} ${MODELDIR} ${MODELDATADIR} ${ATMDIR} ${RIVDIR} \
${REGION} ${EXE} ${BUDGET} \
${TS_NUDGE} ${TI_NUDGE} ${SS_NUDGE} ${SI_NUDGE} \
${ST} ${IT} ${DT}

#---echo "Make vcoord file"
cd ${CDIR}
csh make_vcoord.csh ${WORKDIR} ${NMEM} ${ENODE} ${EPROC}
if($? == 99)then
    echo "***Error: make vcoordfile"
    exit 99
endif

echo "End Model setting"
#===== Submit job ==============================================
echo "Submit ensemble forecast job: ${yyyymmdd} ..."
cd ${WORKDIR}

csh ${CDIR}/submit_job.csh ${EPROC} ${NMEM} ${THREAD} ${REGION} ${EXE} ${WORKDIR} ${PJM_JOBID}

#=== Post process ======================================================
echo "Submit post process"
cd ${CDIR}
csh sbPOM_post.csh ${switch_iau} ${WORKDIR} ${OUTPUT} ${INFO} ${MODELDATADIR} ${NMEM} ${REGION} ${EXE} ${yyyymmdd} ${yyyymmdd_n1}

wait
exit 0
