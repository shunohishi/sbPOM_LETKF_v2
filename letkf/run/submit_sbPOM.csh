#!/bin/csh
#--------------------------------------------------------------------------
# Submit Job cshfile with vcoordfile
#-------------------------------------------------------------------------

#---Argument---------------------------------------------------------------
set switch_iau=${argv[1]}
set EPROC=${argv[2]}
set NMEM=${argv[3]}
set THREAD=${argv[4]}
set REGION=${argv[5]}
set EXE=${argv[6]}
set machine=${argv[7]}
set WORKDIR=${argv[8]}
set JOBID=${argv[9]}

#---Submit Job------------------------------------------------------------
@ IMEM = 1
while(${IMEM} <= ${NMEM})

    set MMMMM=`printf "%05d" ${IMEM}`

    cd ${WORKDIR}/${MMMMM}
    echo "Submit Job: "${WORKDIR}/${MMMMM}
    echo ${JOBID} > JOBID

    setenv OMP_NUM_THREADS ${THREAD}
    setenv PARALLEL ${THREAD}

    if(${machine} == "jss3")then
	(mpiexec -n ${EPROC} -stdout stdout.${REGION} -stderr stderr.${REGION} --vcoordfile vcoordfile ./${EXE} && echo finished > ${WORKDIR}/${MMMMM}/FINISHED_sbPOM &)
    else if(${machine} == "fugaku")then
	(mpiexec -n ${EPROC} -stdout-proc stdout.${REGION} -stderr-proc stderr.${REGION} --vcoordfile vcoordfile ./${EXE} && echo finished > ${WORKDIR}/${MMMMM}/FINISHED_sbPOM &)
    endif
	
    @ IMEM++

end
cd ${WORKDIR}
#-------------------------------------------------------------------------

echo "Start ensemble forecast"#-----------------------------------------
@ iint = 10 #Interval [sec.]
@ isec = 0
while(${isec} >= 0)

    sleep ${iint}s
    @ isec = ${isec} + ${iint}
    @ min = ${isec} / 60
    @ sec = ${isec} % 60

    set FIN_NUM=`find ${WORKDIR} -name FINISHED_sbPOM | wc -l`
    set FIN_OUTPUT=`find ${WORKDIR} -name ${REGION}.nc | wc -l`
    set FIN_IAU=`find ${WORKDIR} -name iau.nc | wc -l`

    if(${switch_iau} == 0 || ${switch_iau} == 2)then
	echo "sbPOM FINISHED JOB[ ${FIN_NUM} / ${NMEM} ] [${FIN_OUTPUT} / ${NMEM} ] ${min}:${sec} elapsed"
    else if(${switch_iau} == 1)then
	echo "sbPOM FINISHED JOB[ ${FIN_NUM} / ${NMEM} ] [${FIN_IAU} / ${NMEM} ] ${min}:${sec} elapsed"
    endif
    
    if((${switch_iau} == 0 || ${switch_iau} == 2) && ${FIN_NUM} == ${NMEM} && ${FIN_OUTPUT} == ${NMEM})then
	echo "End ensemble forecast"
        break
    else if(${switch_iau} == 1 && ${FIN_NUM} == ${NMEM} && ${FIN_IAU} == ${NMEM})then
	echo "End ensemble forecast"
        break    
    endif

end
#-------------------------------------------------------------------------
wait
exit 0
