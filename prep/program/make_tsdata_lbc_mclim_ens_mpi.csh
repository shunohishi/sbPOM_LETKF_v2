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


#set machine="jss3"
set machine="fugaku"

#set RSCUNIT=SORA         #JSS3
#set RSCUNIT=RURI         #JSS3
set RSCUNIT=rscunit_ft01 #Fugaku

set ENS=128      #Ensemble size
set PROC=${ENS}  #Total processor

#--------------------------------------------------------
# Machine environment |
#--------------------------------------------------------

if(${RSCUNIT} == "SORA")then
    set NPROC=48 #Total processor for 1 node
    set EPROC=4  #Processor for 1 member
else if(${RSCUNIT} == "RURI")then
    set NPROC=36 #Total processor for 1 node
    set EPROC=9  #Processor for 1 member
else if(${RSCUNIT} == "rscunit_ft01")then
    set NPROC=48 #Total processor for 1 node
    set EPROC=4  #Processor for 1 member
endif
    
@ THREAD = ${NPROC} / ${EPROC}

if(${PROC} % ${EPROC} == 0)then
    @ NODE = ${PROC} / ${EPROC}
else
    @ NODE = ${PROC} / ${EPROC} + 1
endif

#----------------------------------------------------------
# Module & Subroutine |
#----------------------------------------------------------

set DAM="../module"
set module="src/mod_rmiss.f90 src/mod_parameter.f90 src/mod_gridinfo.f90 ${DAM}/mod_read_soda.f90"
set subroutine="src/sub_read_grid.f90 src/sub_ensemble.f90 src/sub_distance.f90 src/sub_cal_id.f90 src/sub_fillvalue.f90 src/sub_bilinear_interpolation.f90 src/sub_apply_fsm.f90 src/sub_make_ncfile.f90"

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

#------------------------------------------------------------
# Remove standard output and error |
#------------------------------------------------------------

rm -f make_tsdata_lbc_mclim_ens_mpi.*.out make_tsdata_lbc_mclim_ens_mpi.*.err make_tsdata_lbc_mclim_ens_mpi.*.stats
rm -f stdout.make_tsdata_lbc_mclim_ens_mpi* stderr.make_tsdata_lbc_mclim_ens_mpi*

#-------------------------------------------------------------
# Compile |
#-------------------------------------------------------------

if(${RSCUNIT} == "SORA" || ${RSCUNIT} == "rscunit_ft01")then
    mpifrtpx ${module} src/make_tsdata_lbc_mclim_ens_mpi.f90 -o make_tsdata_lbc_mclim_ens_mpi.out ${subroutine} ${option}
else if(${RSCUNIT} == "RURI")then
    mpiifort ${module} src/make_tsdata_lbc_mclim_ens_mpi.f90 -o make_tsdata_lbc_mclim_ens_mpi.out ${subroutine} -qopenmp ${option} ${debug}
endif
    
if(! -f make_tsdata_lbc_mclim_ens_mpi.out)then
    echo "***Error: Nout found make_tsdata_lbc_mclim_ens_mpi.out"
    exit
endif

#---------------------------------------------------------------
# Resource Group |
#---------------------------------------------------------------

if(${RSCUNIT} == "rscunit_ft01" && ${NODE} <= 384)then
    set rscgrp="small"
else if(${RSCUNIT} == "rscunit_ft01" && 385 <= ${NODE} && ${NODE} <= 12288)then
    set rscgrp="large"
endif

#---------------------------------------------------------------
# Submit job |
#---------------------------------------------------------------

echo "Processor:${PROC}, Node:${NODE}"

if(${RSCUNIT} == "SORA")then

jxsub <<EOF
#JX --bizcode R2402
#JX -L rscunit=SORA
#JX -L node=${NODE}
#JX --mpi proc=${PROC}
#JX -L node-mem=29184Mi
#JX -L elapse=01:00:00
#JX -N make_tsdata_lbc_mclim_ens_mpi
#JX -S

export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export OMP_NUM_THREADS=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${ENS} -stdout stdout.make_tsdata_lbc_mclim_ens_mpi -stderr stderr.make_tsdata_lbc_mclim_ens_mpi ./make_tsdata_lbc_mclim_ens_mpi.out

rm -f make_tsdata_lbc_mclim_ens.out

EOF

else if(${RSCUNIT} == "RURI")then

jxsub <<EOF
#JX --bizcode R2402
#JX -L rscunit=RURI
#JX -L vnode=${NODE}
#JX -L vnode-core=${EPROC}
#JX -L vnode-mem=172080Mi
#JX -L elapse=01:00:00
#JX -N make_tsdata_lbc_mclim_ens_mpi
#JX -S

export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export OMP_NUM_THREADS=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${ENS} -ppn ${EPROC} ./make_tsdata_lbc_mclim_ens_mpi.out

rm -f make_tsdata_lbc_mclim_ens.out

EOF

else if(${RSCUNIT} == "rscunit_ft01")then

pjsub <<EOF
#PJM -L rscunit=rscunit_ft01
#PJM -L rscgrp=${rscgrp}
#PJM -L node=${NODE}
#PJM --mpi proc=${PROC}
#PJM -L elapse=01:00:00
#PJM -L retention_state=0
#PJM -g ra000007
#PJM -x PJM_LLIO_GFSCACHE=/vol0004
#PJM --name make_tsdata_lbc_mclim_ens_mpi
#PJM -S

export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export OMP_NUM_THREADS=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${ENS} -stdout-proc stdout.make_tsdata_lbc_mclim_ens_mpi -stderr-proc stderr.make_tsdata_lbc_mclim_ens_mpi ./make_tsdata_lbc_mclim_ens_mpi.out

rm -f make_tsdata_lbc_mclim_ens.out

EOF

endif

rm -f *.mod
