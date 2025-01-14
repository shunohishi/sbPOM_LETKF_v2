#!/bin/csh
#--------------------------------------------------------------------------------
# ARGV |
#--------------------------------------------------------------------------------

set REGION=${argv[1]}
set DIR=${argv[2]}
set EXE=${argv[3]}

set EPROC=${argv[4]}
set NMEM=${argv[5]}
set THREAD=${argv[6]}
set NODE=${argv[7]}
set PROC=${argv[8]}
set machine=${argv[9]}

set yyyy=${argv[10]}
set mm=${argv[11]}
set WORKDIR=${argv[12]}
set elapse_time=${argv[13]}

#--------------------------------------------------------------------------------
# GROUP for Fugaku |
#--------------------------------------------------------------------------------

if(1 <= ${NODE} && ${NODE} <= 384)then
    set group="small"
else if(385 <= ${NODE} && ${NODE} <= 12288)then
    set group="large"
else
    echo "***Error: Check NODE SIZE ${NODE}"
    exit
endif
 
#--------------------------------------------------------------------------------

cd ${WORKDIR}
echo "Submit Job"

if(${machine} == "jss3")then

    
    jxsub <<EOF
#JX --bizcode R2402
#JX -L rscunit=SORA
#JX -L node=${NODE}
#JX --mpi proc=${PROC}
#JX -L node-mem=29184Mi
#JX -L elapse=${elapse_time}
#JX -L proc-crproc=1024
#JX -N sbPOM_${REGION}_${yyyy}${mm}
#JX -S

export OMP_NUM_THREADS=${THREAD}
export PARALLEL=${THREAD}

csh ${DIR}/run/submit_vcoord.csh ${REGION} ${EXE} ${EPROC} ${NMEM} ${THREAD} ${machine} ${WORKDIR} \${PJM_JOBID}

EOF


else if(${machine} == "fugaku")then

    pjsub <<EOF
#PJM -L rscunit=rscunit_ft01
#PJM -L rscgrp=${group}
#PJM -L node=${NODE}
#PJM --mpi proc=${PROC}
#PJM -L elapse=${elapse_time}
#PJM -L proc-crproc=1024
#PJM --name sbPOM_${REGION}_${yyyy}${mm}
#PJM -S

export OMP_NUM_THREAD=${THREAD}
export PARALLEL=${THREAD}

csh ${DIR}/run/submit_vcoord.csh ${REGION} ${EXE} ${EPROC} ${NMEM} ${THREAD} ${machine} ${WORKDIR} \${PJM_JOBID}

EOF

endif
