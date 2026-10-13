#!/bin/csh
#---------------------------------------------------------
# Make monthly climatology tsdata & lbc using SODA |
#---------------------------------------------------------
#
# tsdata: for nuding near boudary conditions
# lbc: Lateral boudary condition
# 
# Monthly climatology in syr-eyr in make_tsdata_lbc_mclim.f90
# *syr:Start year, eyr:End year
# 
# SODA version 3.12.2
# 1980.01-2017.12
#
#--------------------------------------------------------
# Created by S.Ohishi 2018.08
# S.Ohishi 2023.03 Add MPI
# S.Ohishi 2024.12 Add Fugaku
#--------------------------------------------------------

#---Machine
set machine="jss3"
#set machine="fugaku"

#---Ensemble size
set ENS_SIZE=128      #Ensemble size
set TOTAL_PROC=${ENS_SIZE}  #Total processor

#---Elapse time
set elapse_time="01:00:00"

#--------------------------------------------------------
# Machine environment |
#--------------------------------------------------------

if(${machine} == "jss3")then
    set PROC_PER_NODE=48 #Total processor for 1 node
    set ENS_PER_NODE=4  #Ensemble member for 1 node
else if(${machine} == "fugaku")then
    set PROC_PER_NODE=48 #Total processor for 1 node
    set ENS_PER_NODE=4  #Ensemble member for 1 node
endif
    
@ THREAD = ${PROC_PER_NODE} / ${ENS_PER_NODE}

if(${TOTAL_PROC} % ${ENS_PER_NODE} == 0)then
    @ TOTAL_NODE = ${TOTAL_PROC} / ${ENS_PER_NODE}
else
    @ TOTAL_NODE = ${TOTAL_PROC} / ${ENS_PER_NODE} + 1
endif

#----------------------------------------------------------
# Module & Subroutine |
#----------------------------------------------------------

set module="src/mod_rmiss.f90 src/mod_parameter.f90 src/mod_gridinfo.f90 ../module/mod_check_netcdf.f90 ../module/mod_read_soda.f90"
set subroutine="src/sub_read_grid.f90 src/sub_ensemble.f90 src/sub_distance.f90 src/sub_cal_id.f90 src/sub_fillvalue.f90 src/sub_bilinear_interpolation.f90 src/sub_apply_fsm.f90 src/sub_make_ncfile.f90"

#----------------------------------------------------------
# Compiler Option |
#----------------------------------------------------------

if(${machine} == "jss3")then

    set option="-Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign ${fflag_SORA} ${cflag_SORA} ${flib_SORA} ${clib_SORA} ${static_SORA}"
    
else if(${machine} == "fugaku")then

    setenv SPACK_ROOT /vol0004/apps/oss/spack
    source /vol0004/apps/oss/spack/share/spack/setup-env.csh
    spack load ${netcdf_fj}

    set option="-Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign ${fflag_fj} ${cflag_fj} ${flib_fj} ${clib_fj} ${static_fj}"

else

    echo "***Error: machine"
    exit
    
endif

#------------------------------------------------------------
# Remove standard output and error |
#------------------------------------------------------------

rm -f make_tsdata_lbc_mclim_ens_mpi.*.out make_tsdata_lbc_mclim_ens_mpi.*.err make_tsdata_lbc_mclim_ens_mpi.*.stats
rm -f stdout.make_tsdata_lbc_mclim_ens_mpi* stderr.make_tsdata_lbc_mclim_ens_mpi*

#-------------------------------------------------------------
# Compile |
#-------------------------------------------------------------

mpifrtpx ${module} src/make_tsdata_lbc_mclim_ens_mpi.f90 -o make_tsdata_lbc_mclim_ens_mpi.out ${subroutine} ${option}
    
if(! -f make_tsdata_lbc_mclim_ens_mpi.out)then
    echo "***Error: Nout found make_tsdata_lbc_mclim_ens_mpi.out"
    exit
endif

#---------------------------------------------------------------
# Resource Group |
#---------------------------------------------------------------

if(${machine} == "fugaku" && ${TOTAL_NODE} <= 384)then
    set rscgrp="small"
else if(${machine} == "fugaku" && 385 <= ${TOTAL_NODE} && ${TOTAL_NODE} <= 12288)then
    set rscgrp="large"
else if(${machine} == "fugaku" && 12288 < ${TOTAL_NODE})then 
    echo "***Error: Too large total node size ==> ${TOTAL_NODE}"
    exit 99
endif

#---------------------------------------------------------------
# Submit job |
#---------------------------------------------------------------

echo "Total processor: ${TOTAL_PROC}"
echo "Total Node: ${TOTAL_NODE}"

if(${machine} == "jss3")then

jxsub <<EOF
#JX --bizcode R2402
#JX -L rscunit=SORA
#JX -L node=${TOTAL_NODE}
#JX --mpi proc=${TOTAL_PROC}
#JX -L node-mem=29184Mi
#JX -L elapse=${elapse_time}
#JX -N make_tsdata_lbc_mclim_ens_mpi
#JX -S

ulimit -s unlimited
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export OMP_NUM_THREADS=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${TOTAL_PROC} -stdout stdout.make_tsdata_lbc_mclim_ens_mpi -stderr stderr.make_tsdata_lbc_mclim_ens_mpi ./make_tsdata_lbc_mclim_ens_mpi.out

rm -f make_tsdata_lbc_mclim_ens_mpi.out

EOF

else if(${machine} == "fugaku")then

pjsub <<EOF
#PJM -L rscunit=rscunit_ft01
#PJM -L rscgrp=${rscgrp}
#PJM -L node=${TOTAL_NODE}
#PJM --mpi proc=${TOTAL_PROC}
#PJM -L elapse=${elapse_time}
#PJM -L retention_state=0
#PJM -g ra000007
#PJM -x PJM_LLIO_GFSCACHE=/vol0004
#PJM --name make_tsdata_lbc_mclim_ens_mpi
#PJM -S

ulimit -s unlimited
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export OMP_NUM_THREADS=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${TOTAL_PROC} -stdout-proc stdout.make_tsdata_lbc_mclim_ens_mpi -stderr-proc stderr.make_tsdata_lbc_mclim_ens_mpi ./make_tsdata_lbc_mclim_ens_mpi.out

rm -f make_tsdata_lbc_mclim_ens_mpi.out

EOF

endif

rm -f *.mod
