#!/bin/csh
#=== Argument ====================================================================

set REGION=${argv[1]}
set DIR=${argv[2]}
set NPROC=${argv[3]}
set NMEM=${argv[4]}
set BUDGET=${argv[5]}
set RM_ENS=${argv[6]}
set machine=${argv[7]}

set yyyy_s=${argv[8]}
set mm_s=${argv[9]}
set dd_s=${argv[10]}
set yyyy=${argv[11]}
set mm=${argv[12]}
set dd=${argv[13]}
set nday=${argv[14]}
set WORKDIR=${argv[15]}

@ THREAD = 12
if(${NMEM} * ${THREAD} % ${NPROC} == 0)then
    @ NODE = ${NMEM} * ${THREAD} / ${NPROC}
else
    @ NODE = ${NMEM} * ${THREAD} / ${NPROC} + 1
endif

#=================================================================================
# ELAPSE TIME
#=================================================================================

if(${nday} == 1)then
    set elapse_time="00:20:00"
else
    set elapse_time="03:00:00"
endif

#=================================================================================
echo "Make ens_info.txt"
#=================================================================================

#Check mesp_ens_mpi.out
cd ${DIR}/run/
if(! -f mesp_ens_mpi.out && ! -f ${WORKDIR}/mesp_ens_mpi.out)then
    echo "***Error: Not found ${DIR}/run/mesp_ens_mpi.out"
    exit
else if(! -f ${WORKDIR}/mesp_ens_mpi.out)then
    ln -s ${DIR}/run/mesp_ens_mpi.out ${WORKDIR}/mesp_ens_mpi.out
endif

#=================================================================
echo "Submit Job: Post process"
#=================================================================

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

EOF

endif

#---------------------------------------------------
wait
exit 0
