#!/bin/csh
#--------------------------------------------------------------------------------
# ARGV |
#--------------------------------------------------------------------------------

set REGION=${argv[1]}
set WORKDIR=${argv[2]}
set NODE=${argv[3]}
set PROC=${argv[4]}
set THREAD=${argv[5]}
set machine=${argv[6]}
set EXE=${argv[7]}

set yyyy=${argv[8]}
set mm=${argv[9]}
set elapse_time=${argv[10]}

#-------------------------------------------------------------------------------

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
#JX -N sbPOM_${REGION}_${yyyy}${mm}
#JX -S

export OMP_NUM_THREAD=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${PROC} -stdout stdout.${REGION} -stderr stderr.${REGION} ./${EXE}

echo finished > FINISHED

EOF

else if(${machine} == "fugaku")then

    pjsub <<EOF
#PJM -L rscunit=rscunit_ft01
#PJM -L rscgrp=small
#PJM -L node=${NODE}
#PJM --mpi proc=${PROC}
#PJM -L elapse=${elapse_time}
#PJM --name sbPOM_${REGION}_${yyyy}${mm}
#PJM -S

export OMP_NUM_THREAD=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${PROC} -stdout-proc stdout.${REGION} -stderr-proc stderr.${REGION} ./${EXE}

echo finished > FINISHED

EOF

endif
