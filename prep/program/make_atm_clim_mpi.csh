#!/bin/csh
#----------------------------------------------------------------
# Make Climatological atmospheric boundary condition |
#----------------------------------------------------------------
#
# Modify "module setting" in make_atm_clim_mpi.f90
#
#-------------------------------------------------------------
# Created by S.Ohishi 2018.08
# S.Ohishi 2021.12 add JRA55do
# S.Ohishi 2024.12 add Fugaku
# S.Ohishi 2024.12 add Climatology
# S.Ohishi 2026.10 add MPI
#-------------------------------------------------------------

#---Machine
set machine="jss3"
#set machine="fugaku"

#---RSCUNIT
set RSCUNIT=SORA         #JSS3
#set RSCUNIT=rscunit_ft01 #Fugaku

#---Setting
set TOTAL_PROC=48 #Total processor
set TOTAL_NODE=4  #Total node
set elapse_time="24:00:00"

#-------------------------------------------------------------
# Machine environment |
#-------------------------------------------------------------

if(${RSCUNIT} == "SORA")then
    set PROC_PER_NODE=48 #Processor per node
else if(${RSCUNIT} == "rscunit_ft01")then
    set PROC_PER_NODE=48 #Processor per nod
endif

@ THREAD = ${PROC_PER_NODE} * ${TOTAL_NODE} / ${TOTAL_PROC}

echo "TOTAL_NODE: ${TOTAL_NODE}"
echo "TOTAL_PROC: ${TOTAL_PROC}"
echo "THREAD: ${THREAD}"

#-------------------------------------------------------------
# Module & Subroutine |
#-------------------------------------------------------------

set module="src/mod_rmiss.f90 src/mod_parameter.f90 src/mod_gridinfo.f90 ../module/mod_check_netcdf.f90 ../module/mod_read_jra55do.f90 ../module/mod_read_era5.f90"
set subroutine="src/sub_read_grid.f90 src/sub_distance.f90 src/sub_cal_id.f90 src/sub_apply_fsm.f90 src/sub_fillvalue.f90 src/sub_bilinear_interpolation.f90 src/sub_time.f90 src/sub_make_ncfile.f90"

#-------------------------------------------------------------
# Compiler Option |
#-------------------------------------------------------------

if(${machine} == "jss3")then

    set option="-Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign ${fflag_SORA} ${cflag_SORA} ${flib_SORA} ${clib_SORA} ${static_SORA}"

else if(${machine} == "fugaku")then

    setenv SPACK_ROOT /vol0004/apps/oss/spack
    source /vol0004/apps/oss/spack/share/spack/setup-env.csh
    spack load ${netcdf_fj}

    set option="-Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign ${fflag_fj} ${cflag_fj} ${flib_fj} ${clib_fj} ${static_fj}"
    
endif

#--------------------------------------------------------------
# Compile |
#--------------------------------------------------------------

mpifrtpx ${module} src/make_atm_clim_mpi.f90 -o make_atm_clim_mpi.out ${subroutine} ${option}

if(! -f make_atm_clim_mpi.out)then
    echo "***Error: Not found make_atm_clim_mpi.out"
    exit
endif

#--------------------------------------------------------------
# Fugaku resource group |
#--------------------------------------------------------------

if(${RSCUNIT} == "rscunit_ft01" && ${TOTAL_NODE} <= 384)then
    set rscgrp="small"
else if(${RSCUNIT} == "rscunit_ft01" && 385 <= ${TOTAL_NODE} && ${TOTAL_NODE} <= 12288)then
    set rscgrp="large"
endif

#---------------------------------------------------------------
# Submit job |
#---------------------------------------------------------------

rm -f make_atm_clim_mpi.*.out make_atm_clim_mpi.*.err make_atm_clim_mpi.*.stats
rm -f stdout.make_atm_clim_mpi stderr.make_atm_clim_mpi

if(${RSCUNIT} == "SORA")then

jxsub <<EOF
#JX --bizcode R2402
#JX -L rscunit=SORA
#JX -L node=${TOTAL_NODE}
#JX --mpi proc=${TOTAL_PROC}
#JX -L node-mem=29184Mi
#JX -L elapse=${elapse_time}
#JX -N make_atm_clim_mpi
#JX -S

ulimit -s unlimited
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export OMP_NUM_THREADS=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${TOTAL_PROC} -stdout stdout.make_atm_clim_mpi -stderr stderr.make_atm_clim_mpi ./make_atm_clim_mpi.out

rm -f make_atm_clim_mpi.out

EOF

else if(${RSCUNIT} == "rscunit_ft01")then

pjsub <<EOF
#PJM -L rscunit=rscunit_ft01
#PJM -L rscgrp=${rscgrp}
#PJM -L node=${TOTAL_NODE}
#PJM --mpi proc=${TOTAL_PROC}
#PJM -L elapse=${elapse_time}
#PJM -L retention_state=0
#PJM -g ra000007
#PJM -x PJM_LLIO_GFSCACHE=/vol0004
#PJM --name make_atm_clim_mpi
#PJM -S

ulimit -s unlimited
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export OMP_NUM_THREADS=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${TOTAL_PROC} -stdout-proc stdout.make_atm_clim_mpi -stderr-proc stderr.make_atm_clim_mpi ./make_atm_clim_mpi.out

rm -f make_atm_clim_mpi.out

EOF

endif
