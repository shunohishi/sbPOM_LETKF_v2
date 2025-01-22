#!/bin/csh
#=== Argument ====================================================================

set DIR=${argv[1]}
set WORKDIR=${argv[2]}
set INFO=${argv[3]}
set REGION=${argv[4]}
set nday=${argv[5]}

set yyyy_s=${argv[6]}
set mm_s=${argv[7]}
set dd_s=${argv[8]}
set yyyy=${argv[9]}
set mm=${argv[10]}
set dd=${argv[11]}

set NPROC=${argv[12]}
set NMEM=${argv[13]}
set BUDGET=${argv[14]}
set RM_ENS=${argv[15]}

set machine=${argv[16]}

@ THREAD = 12
if(${NMEM} * ${THREAD} % ${NPROC} == 0)then
    @ NODE = ${NMEM} * ${THREAD} / ${NPROC}
else
    @ NODE = ${NMEM} * ${THREAD} / ${NPROC} + 1
endif

#=================================================================================

if(${nday} == 1)then
    set elapse_time="00:30:00"
else
    set elapse_time="02:00:00"
endif

#=================================================================================

echo "Make ens_info.txt" #----------------------------
cd ${DIR}/run/
if(! -f mesp_ens_mpi.out && ! -f ${WORKDIR}/mesp_ens_mpi.out)then
    echo "***Error: Not found ${DIR}/run/mesp_ens_mpi.out"
    exit
else if(! -f ${WORKDIR}/mesp_ens_mpi.out)then
    ln -s ${DIR}/run/mesp_ens_mpi.out ${WORKDIR}/mesp_ens_mpi.out
endif

echo "Submit Job: Post process"#----------------------
cd ${WORKDIR}
echo "Ensemble size: ${NMEM}, Total Node: ${NODE}"

if(${machine} == "jss3")then

jxsub <<EOF
#JX --bizcode R2402
#JX -L rscunit=SORA
#JX -L node=${NODE}
#JX --mpi proc=${NMEM}
#JX -L node-mem=29184Mi
#JX -L elapse=${elapse_time}
#JX -L proc-crproc=1024
#JX -N mesp_ens_mpi_${yyyy}${mm}${dd}
#JX -S

export OMP_NUM_THREADS=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${NMEM} -stdout stdout.mesp_ens_mpi_${yyyy}${mm}${dd} -stderr stderr.mesp_ens_mpi_${yyyy}${mm}${dd} ./mesp_ens_mpi.out ${WORKDIR} ${REGION} ${yyyy_s} ${mm_s} ${dd_s} ${yyyy} ${mm} ${dd} ${NMEM} ${nday} ${BUDGET} ${RM_ENS} &
wait

if [ -s "${WORKDIR}/stdout.mesp_ens_mpi_${yyyy}${mm}${dd}" ]; then
    mv "${WORKDIR}/stdout.mesp_ens_mpi_${yyyy}${mm}${dd}" "${INFO}/stdout.mesp_ens_mpi_${yyyy}${mm}${dd}"
else
    rm -f ${WORKDIR}/stdout.mesp_ens_mpi_${yyyy}${mm}${dd}
fi

if [ -s "${WORKDIR}/stderr.mesp_ens_mpi_${yyyy}${mm}${dd}" ]; then
    mv "${WORKDIR}/stderr.mesp_ens_mpi_${yyyy}${mm}${dd}" "${INFO}/stderr.mesp_ens_mpi_${yyyy}${mm}${dd}"
else
    rm -f ${WORKDIR}/stderr.mesp_ens_mpi_${yyyy}${mm}${dd}
fi

EOF

else if(${machine} == "fugaku")then

pjsub <<EOF
#PJM -L rscunit=rscunit_ft01
#PJM -L rscgrp=small
#PJM -L node=${NODE}
#PJM --mpi proc=${NMEM}
#PJM -L elapse=${elapse_time}
#PJM -L proc-crproc=1024
#PJM --name mesp_ens_mpi_${yyyy}${mm}${dd}
#PJM -S

export OMP_NUM_THREAD=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${NMEM} -stdout-proc stdout.mesp_ens_mpi_${yyyy}${mm}${dd} -stderr-proc stderr.mesp_ens_mpi_${yyyy}${mm}${dd} ./mesp_ens_mpi.out ${WORKDIR} ${REGION} ${yyyy_s} ${mm_s} ${dd_s} ${yyyy} ${mm} ${dd} ${NMEM} ${nday} ${BUDGET} ${RM_ENS} &
wait

if [ -s "${WORKDIR}/stdout.mesp_ens_mpi_${yyyy}${mm}${dd}.1.0" ]; then
    mv "${WORKDIR}/stdout.mesp_ens_mpi_${yyyy}${mm}${dd}.1.0" "${INFO}/stdout.mesp_ens_mpi_${yyyy}${mm}${dd}"
else
    rm -f ${WORKDIR}/stdout.mesp_ens_mpi_${yyyy}${mm}${dd}
fi

if [ -s "${WORKDIR}/stderr.mesp_ens_mpi_${yyyy}${mm}${dd}.1.0" ]; then
    mv "${WORKDIR}/stderr.mesp_ens_mpi_${yyyy}${mm}${dd}.1.0" "${INFO}/stderr.mesp_ens_mpi_${yyyy}${mm}${dd}"
else
    rm -f ${WORKDIR}/stderr.mesp_ens_mpi_${yyyy}${mm}${dd}
fi

EOF

endif

#---------------------------------------------------

wait
exit 0
