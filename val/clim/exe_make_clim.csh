#!/bin/csh
#---------------------------------------------------------------
# Make monthly and annual climatologies |
#---------------------------------------------------------------
# Created by S.Ohishi 2026.01
#---------------------------------------------------------------

#---------------------------------------------------------------
# Select machine |
#---------------------------------------------------------------

#---Machine
#set machine="jss3"
#set machine="fugaku"
set machine="rc"

#---RSC Unit (JSS3 or Fugaku)
set RSCUNIT=SORA         #JSS3
#set RSCUNIT=RURI         #JSS3
#set RSCUNIT=rscunit_ft01 #Fugaku

#---Partition (only for R-CCS Cloud)
#set partition="r340"  #Execute on r340
set partition="genoa"  #Execute on r340/genoa
#set partition="fx700" #Execute on fx700 ***Issue in MPI+NetCDF (As of 17 Sep 2026)***

#---Processor size
#set TOTAL_PROC=1920  #Total processor
set TOTAL_PROC=576  #Total processor (18 yr at genoa)

#---Elapse time
set elapse_time="24:00:00"

#--------------------------------------------------------
# Machine environment |
#--------------------------------------------------------

if(${machine} == "jss3" && ${RSCUNIT} == "SORA")then
    set NODE_PROC=48 #Processor per node
    set THREAD=12    #Number of thread per mpi
else if(${machine} == "jss3" && ${RSCUNIT} == "RURI")then
    set NODE_PROC=36 #Processor per node
    set THREAD=1     #Number of thread per mpi
else if(${machine} == "fugaku" && ${RSCUNIT} == "rscunit_ft01")then
    set NODE_PROC=48 #Processor per node
    set THREAD=12    #Number of thread per mpi
else if(${machine} == "rc" && ${partition} == "r340")then
    set NODE_PROC=4 #Processor per node
    set THREAD=1    #Number of thread per mpi
else if(${machine} == "rc" && ${partition} == "genoa")then
    set NODE_PROC=96 #Processor per node
    set THREAD=1    #Number of thread per mpi    
else if(${machine} == "rc" && ${partition} == "fx700")then
    set NODE_PROC=48 #Processor per node
    set THREAD=1    #Number of thread per mpi
endif
    
if((${machine} == "jss3" || ${machine} == "fugaku" || ${machine} == "rc") && ${TOTAL_PROC} % ${NODE_PROC} == 0)then
    @ NODE = ${TOTAL_PROC} / ${NODE_PROC} * ${THREAD}
    @ NODE_MPI = ${NODE_PROC} / ${THREAD}
else if(${machine} == "jss3" || ${machine} == "fugaku" || ${machine} == "rc")then
    @ NODE = (${TOTAL_PROC} / ${NODE_PROC} + 1) * ${THREAD}
    @ NODE_MPI = ${NODE_PROC} / ${THREAD}
endif

#----------------------------------------------------------
# Compiler Option |
#----------------------------------------------------------

if(${machine} == "jss3" && ${RSCUNIT} == "SORA")then

    #Parallel NetCDF on SORA
    set debug=""
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

    set debug=""
    set option="-Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign ${fflag_fj} ${cflag_fj} ${flib_fj} ${clib_fj} ${static_fj}"

else if(${machine} == "rc")then

    #set debug="-g -fbacktrace -fcheck=all"
    set debug=""
    set fflag=`nf-config --fflags`
    set flib=`nf-config --flibs`
    set clib=`nc-config --libs`
    set option="${fflag} ${flib} ${clib} -ffree-line-length-none"

else

    echo "***Error: machine or RSCUNIT"
    exit
    
endif

#---------------------------------------------------------------
# Subroutine & Module |
#---------------------------------------------------------------

set module="../module/mod_rmiss.f90 ../module/mod_gridinfo.f90 ../module/mod_read_lora_v20.f90 ../module/mod_read_bran2020.f90 ../module/mod_read_fora_np60.f90 ../module/mod_read_glorys010.f90 ../module/mod_read_glorys025.f90 ../module/mod_read_jcope_fgo.f90 mod_setting.f90 mod_make_ncfile.f90 mod_io.f90"
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

if((${machine} == "jss3" || ${machine} == "fugaku") && (${RSCUNIT} == "SORA" || ${RSCUNIT} == "rscunit_ft01"))then
    mpifrtpx ${module} main_make_mclim.f90 -o make_mclim.out ${subroutine} ${option}
else if(${machine} == "jss3" && ${RSCUNIT} == "RURI")then
    mpiifort ${module} main_make_mclim.f90 -o make_mclim.out ${subroutine} -qopenmp ${option} ${debug}
else if(${machine} == "rc")then
    mpifort ${module} main_make_mclim.f90 -o make_mclim.out ${subroutine} -fopenmp ${option} ${debug}
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
if(${machine} == "jss3" && ${RSCUNIT} == "SORA")then

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

else if(${machine} == "jss3" && ${RSCUNIT} == "RURI")then

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

else if(${machine} == "fugaku" && ${RSCUNIT} == "rscunit_ft01")then

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

else if(${machine} == "rc" && (${partition} == "r340" || ${partition} == "genoa"))then

sbatch -N ${NODE} -t ${elapse_time} -p ${partition} --ntasks=${TOTAL_PROC} --cpus-per-task=${THREAD} --job-name=make_clim submit_job_est_bias.sh ${THREAD} ${TOTAL_PROC}

endif

rm -f *.mod
