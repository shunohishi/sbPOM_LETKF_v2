#!/bin/csh
#===== ARGUMENT ============================================
set switch_iau=${argv[1]} #0: Intermittent, #1: IAU
set CDIR=${argv[2]}
set WORKDIR=${argv[3]}
set OUTPUT=${argv[4]}
set INFO=${argv[5]}
set NMEM=${argv[6]}
set machine=${argv[7]}

set OBSDIR=${argv[8]}
set OBSFILE=${argv[9]}

set LETKFDIR=${argv[10]}
set EXE=${argv[11]}
set PROC=${argv[12]}
set THREAD=${argv[13]}

set MODELDATADIR=${argv[14]}
set IT=${argv[15]}

@ BT = ${IT} - 1
set date_n1=`perl juldays.prl ${BT}`
set yyyy_n1=`printf "%04d" ${date_n1[1]}`; set mm_n1=`printf "%02d" ${date_n1[2]}`; set dd_n1=`printf "%02d" ${date_n1[3]}`
set yyyymmdd_n1=${yyyy_n1}${mm_n1}${dd_n1}

set date=`perl juldays.prl ${IT}`
set yyyy=`printf "%04d" ${date[1]}`; set mm=`printf "%02d" ${date[2]}`; set dd=`printf "%02d" ${date[3]}`
set yyyymmdd=${yyyy}${mm}${dd}

#=== Setting ===============================================

#---echo "Link LETKF EXE file"
rm -f ${WORKDIR}/${EXE}
ln -s ${LETKFDIR}/${EXE} ${WORKDIR}/${EXE}

#---echo "Link Observational data: obsTT.nc"
cd ${WORKDIR}/
set filename=${OBSDIR}/${OBSFILE}${yyyymmdd}.nc
if(-f ${filename})then
    rm -f obs01.nc
    ln -s ${filename} obs01.nc
else
    echo "***Error: ${filename} not found"; exit 99
endif
	    
#---echo "Link Grid data"
ln -s ${MODELDATADIR}/grid.nc grid.nc

#---echo "Clear & restart --> fcTT***.nc, anal***.nc"
@ IMEM = 1
while(${IMEM} <= ${NMEM})

    set MMMMM=`printf "%05d" ${IMEM}`
    if(${switch_iau} == 0) set filename="${OUTPUT}/${MMMMM}/restart.${yyyymmdd_n1}.nc"
    if(${switch_iau} == 1) set filename="${OUTPUT}/${MMMMM}/iau.${yyyymmdd}.nc"

    if(! -f ${filename})then
	echo "***Error: ${filename} not found"; exit 99
    else
	(rm -f fa01${MMMMM}.nc && ln -s ${filename} fa01${MMMMM}.nc &)
    endif

    @ IMEM++

end
wait

#---echo "Copy anal001.nc --> fcst_me,sp.nc anal_me,sp.nc"
if(${switch_iau} == 0) set filename="${OUTPUT}/00001/restart.${yyyymmdd_n1}.nc"
if(${switch_iau} == 1) set filename="${OUTPUT}/00001/iau.${yyyymmdd}.nc"
(cp ${filename} fcst_mean.nc &)
(cp ${filename} fcst_sprd.nc &)
(cp ${filename} anal_mean.nc &)
(cp ${filename} anal_sprd.nc &)
wait

#=== Submit LETKF Job ==================================================================
echo "Submit LETKF Job"

setenv OMP_NUM_THREADS ${THREAD}
setenv PARALLEL ${THREAD}

echo "PROCESSOR: ${PROC}", "THREAD: ${THREAD}"

echo ${PJM_JOBID} > JOBID

if(${machine} == "jss3")then
    (mpiexec -n ${PROC} -stdout stdout.letkf -stderr stderr.letkf ./${EXE} && echo finished > FINISHED_LETKF &)
else if(${machine} == "fugaku")then
    (mpiexec -n ${PROC} -stdout-proc stdout.letkf -stderr-proc stderr.letkf ./${EXE} && echo finished > FINISHED_LETKF &)
endif
    
echo "Wait for the LETKF job to complete"
@ isec = 0
@ imin = 0
@ iint = 10
while(${isec} >= 0)

    #TIMER
    sleep ${iint}s
    @ isec = ${isec} + ${iint}
    @ imin = ${isec} / 60
    @ imod = ${isec} % 60

    set FIN_NUM=`find ${WORKDIR} -name FINISHED_LETKF | wc -l`
    echo "LETKF FINISHED [${FIN_NUM}/1] ; ${imin}:${imod} elapsed"

    #BREAK
    if(${FIN_NUM} == 1)then
	echo "End LETKF job"
	break
    endif

end #isec

#=== Post process ===================================================================
echo "Submit post process"
cd ${CDIR}
csh LETKF_post.csh ${WORKDIR} ${OUTPUT} ${INFO} ${NMEM} ${machine} ${EXE} ${yyyymmdd_n1} ${yyyymmdd}
#============================================================================#

wait
exit 0
