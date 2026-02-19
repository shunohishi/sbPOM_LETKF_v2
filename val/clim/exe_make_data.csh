#!/bin/csh
#---------------------------------------------------------------
# Make monthly and annual climatologies |
#---------------------------------------------------------------
# Created by S.Ohishi 2026.01
#---------------------------------------------------------------

#---------------------------------------------------------------
# Select machine |
#---------------------------------------------------------------

set machine="jss3"
#set machine="fugaku"

set RSCUNIT=SORA         #JSS3
#set RSCUNIT=RURI         #JSS3
#set RSCUNIT=rscunit_ft01 #Fugaku

set TOTAL_PROC=1920  #Total processor

set elapse_time="06:00:00"

#--------------------------------------------------------
# Machine environment |
#--------------------------------------------------------

if(${RSCUNIT} == "SORA")then
    set NODE_PROC=48 #Processor per node
    set THREAD=12    #Number of thread per mpi
else if(${RSCUNIT} == "RURI")then
    set NODE_PROC=36 #Processor per node
    set THREAD=1     #Number of thread per mpi
else if(${RSCUNIT} == "rscunit_ft01")then
    set NODE_PROC=48 #Processor per node
    set THREAD=12    #Number of thread per mpi
endif
    
if(${TOTAL_PROC} % ${NODE_PROC} == 0)then
    @ NODE = ${TOTAL_PROC} / ${NODE_PROC} * ${THREAD}
    @ NODE_MPI = ${NODE_PROC} / ${THREAD}
else
    @ NODE = (${TOTAL_PROC} / ${NODE_PROC} + 1) * ${THREAD}
    @ NODE_MPI = ${NODE_PROC} / ${THREAD}
endif

#----------------------------------------------------------
# Compiler Option |
#----------------------------------------------------------

if(${machine} == "jss3" && ${RSCUNIT} == "SORA")then

    #Parallel NetCDF on SORA
    set option="-Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign ${fflag_SORA} ${cflag_SORA} ${flib_SORA} ${clib_SORA} ${static_SORA}"

else if(${machine} == "jss3" && ${RSCUNIT} == "RURI")then
    
    #NetCDF on RURI
    set debug="-CB -traceback -g"
    #set debug=""
    set option="-assume byterecl -convert big_endian -mcmodel=medium -shared-intel ${fflag_RURI} ${cflag_RURI} ${flib_RURI} ${clib_RURI} ${static_RURI}"
    
else if(${machine} == "fugaku" && ${RSCUNIT} == "rscunit_ft01")then

    setenv SPACK_ROOT /vol0004/apps/oss/spack
    source /vol0004/apps/oss/spack/share/spack/setup-env.csh
    spack load ${netcdf_fj}

    set option="-Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign ${fflag_fj} ${cflag_fj} ${flib_fj} ${clib_fj} ${static_fj}"

else

    echo "***Error: machine or RSCUNIT"
    exit
    
endif

#---------------------------------------------------------------
# Subroutine & Module |
#---------------------------------------------------------------

set module="../module/mod_rmiss.f90 ../module/mod_gridinfo.f90 ../module/mod_read_lora_v20.f90 ../module/mod_read_glorys025.f90 mod_setting.f90 mod_make_ncfile.f90 mod_io.f90"
set subroutine="sub_cal_id.f90 sub_bilinear_interpolation.f90"

#---------------------------------------------------------------
# Remove standard output and error |
# Remove execution file            |
#---------------------------------------------------------------

rm -f make_mclim.*.out make_mclim.*.err make_mclim.*.stats stdout.make_mclim stderr.make_mclim
rm -f make_mclim.out

#---------------------------------------------------------------
# Compile |
#---------------------------------------------------------------

if(${RSCUNIT} == "SORA" || ${RSCUNIT} == "rscunit_ft01")then
    mpifrtpx ${module} main_make_mclim.f90 -o make_mclim.out ${subroutine} ${option}
else if(${RSCUNIT} == "RURI")then
    mpiifort ${module} main_make_mclim.f90 -o make_mclim.out ${subroutine} -qopenmp ${option} ${debug}
endif
    
if(! -f make_mclim.out)then
    echo "***Error: Not found make_mclim.out"
    exit
endif

#---------------------------------------------------------------
# Execution |
#---------------------------------------------------------------

#---Make dir
if(! -d dat) mkdir dat

echo "NODE: ${NODE}"
echo "MPI total processor: ${TOTAL_PROC}"
echo "MPI processor per node: ${NODE_MPI}, THREAD: ${THREAD}"

#---Execulte
if(${RSCUNIT} == "SORA")then

jxsub <<EOF
#JX --bizcode R2402
#JX -L rscunit=SORA
#JX -L node=${NODE}
#JX --mpi proc=${TOTAL_PROC}
#JX -L node-mem=29184Mi
#JX -L elapse=${elapse_time}
#JX -N make_mclim
#JX -S

export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export OMP_NUM_THREADS=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${TOTAL_PROC} -stdout stdout.make_mclim -stderr stderr.make_mclim ./make_mclim.out

rm -f make_mclim.out

EOF

else if(${RSCUNIT} == "RURI")then

jxsub <<EOF
#JX --bizcode R2402
#JX -L rscunit=RURI
#JX -L vnode=${NODE}
#JX -L vnode-core=${NODE_MPI}
#JX -L vnode-mem=172080Mi
#JX -L elapse=${elapse_time}
#JX -N make_mclim
#JX -S

export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export OMP_NUM_THREADS=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${TOTAL_PROC} -ppn ${NODE_MPI} ./make_mclim.out

rm -f make_mclim.out

EOF

else if(${RSCUNIT} == "rscunit_ft01")then

pjsub <<EOF
#PJM -L rscunit=rscunit_ft01
#PJM -L rscgrp=${rscgrp}
#PJM -L node=${NODE}
#PJM --mpi proc=${TOTAL_PROC}
#PJM -L elapse=${elapse_time}
#PJM -L retention_state=0
#PJM -g ra000007
#PJM -x PJM_LLIO_GFSCACHE=/vol0004
#PJM --name make_mclim
#PJM -S

export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export OMP_NUM_THREADS=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${TOTAL_PROC} -stdout-proc stdout.make_mclim -stderr-proc stderr.make_mclim ./make_mclim.out

rm -f make_mclim.out

EOF

endif

rm -f *.mod
