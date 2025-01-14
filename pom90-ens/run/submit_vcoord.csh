#!/bin/csh
#--------------------------------------------------------------------------
# Submit Job cshfile with vcoordfile
#-------------------------------------------------------------------------

#---Argument---------------------------------------------------------------

set REGION=${argv[1]}
set EXE=${argv[2]}
set EPROC=${argv[3]}
set NMEM=${argv[4]}
set THREAD=${argv[5]}
set machine=${argv[6]}
set WORKDIR=${argv[7]}
set JOBID=${argv[8]}

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
    echo "sbPOM FINISHED JOB[ ${FIN_NUM} / ${NMEM} ] ${min}:${sec} passed"

    if(${FIN_NUM} == ${NMEM})then
	echo "End ensemble forecast"
        break
    endif

end
#-------------------------------------------------------------------------
wait
exit 0
