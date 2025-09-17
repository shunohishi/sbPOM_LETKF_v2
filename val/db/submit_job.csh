#!/bin/csh

#---Argument
set machine=${argv[1]}
set syr=${argv[2]}; set smon=${argv[3]}; set sday=${argv[4]}
set eyr=${argv[5]}; set emon=${argv[6]}; set eday=${argv[7]}
set yyyy=${argv[8]}; set mm=${argv[9]}

#---Submit job

if(${machine} == "jss3")then

jxsub <<EOF
#JX --bizcode R2402
#JX -L rscunit=SORA
#JX -L node=1
#JX --mpi proc=1
#JX -L node-mem=29184Mi
#JX -L elapse=04:00:00
#JX -N make_data_${yyyy}${mm}
#JX -S

export OMP_NUM_THREAD=48
export PARALLEL=48

export FLIB_BARRIER=SOFT

numactl --interleave=all ./make_data.out ${syr} ${smon} ${sday} ${eyr} ${emon} ${eday} > ${yyyy}${mm}.log

EOF

else if(${machine} == "fugaku")then

pjsub <<EOF
#PJM -L  "node=1"
#PJM -L  "rscgrp=small"
#PJM -L  "elapse=03:00:00"
#PJM -L rscunit=rscunit_ft01
#PJM -g ra000007
#PJM -x PJM_LLIO_GFSCACHE=/vol0004
#PJM -N make_data_${yyyy}${mm}
#PJM -S

./make_data.out ${syr} ${smon} ${sday} ${eyr} ${emon} ${eday} > ${yyyy}${mm}.log

EOF

endif
