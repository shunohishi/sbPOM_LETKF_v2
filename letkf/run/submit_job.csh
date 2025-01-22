#!/bin/csh
#--------------------------------------------------------------------------
# Submit Job cshfile with vcoordfile
#-------------------------------------------------------------------------

#---Argument---------------------------------------------------------------
set EPROC=${argv[1]}; set NMEM=${argv[2]}; set THREAD=${argv[3]}
set REGION=${argv[4]}; set EXE=${argv[5]}; set WORKDIR=${argv[6]}
set JOBID=${argv[7]}

#---Submit Job------------------------------------------------------------
@ IMEM = 1
while(${IMEM} <= ${NMEM})

    set MMMMM=`printf "%05d" ${IMEM}`

    cd ${WORKDIR}/${MMMMM}
    echo "Submit Job: "${WORKDIR}/${MMMMM}
    echo ${JOBID} > JOBID

    setenv OMP_NUM_THREADS ${THREAD}
    setenv PARALLEL ${THREAD}
    (mpiexec -n ${EPROC} -stdout stdout.${REGION} -stderr stderr.${REGION} --vcoordfile vcoordfile ./${EXE} && echo finished > ${WORKDIR}/${MMMMM}/FINISHED_sbPOM &)

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
    echo "sbPOM FINISHED JOB[ ${FIN_NUM} / ${NMEM} ] ${min}:${sec} elapsed"

    if(${FIN_NUM} == ${NMEM})then
	echo "End ensemble forecast"
        break
    endif

end
#-------------------------------------------------------------------------
wait
exit 0
