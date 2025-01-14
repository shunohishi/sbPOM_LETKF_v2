#!/bin/csh
#-------------------------------------------------------------------
# Make T/S Monthly climatology using WOA18
#-------------------------------------------------------------------
#
# Created by S.Ohishi 2018.08
#
# S.Ohishi 2024.12 Add Fugaku
#
#-------------------------------------------------------------------

#set machine="jss3"
set machine="fugaku"

#set RSCUNIT=""              #No HPC
#set RSCUNIT="RURI"          #On JSS3 if necessary
set RSCUNIT="rscunit_ft01" #On Fugaku if necessary

#-------------------------------------------------------------------
# Module & Subroutine |
#-------------------------------------------------------------------

set DAM="../module"
set module="src/mod_rmiss.f90 src/mod_parameter.f90 src/mod_gridinfo.f90 ${DAM}/mod_read_woa18_month.f90 ${DAM}/mod_read_woa18_season_annual.f90"
set subroutine="src/sub_read_grid.f90 src/sub_fillvalue.f90 src/sub_apply_fsm.f90 src/sub_distance.f90 src/sub_cal_id.f90 src/sub_bilinear_interpolation.f90 src/sub_make_ncfile.f90"

#------------------------------------------------------------------
# Compiler Option |
#------------------------------------------------------------------

if(${machine} == "jss3")then

    #set debug="-CB -traceback -g"
    set debug=""
    set option="-assume byterecl -convert big_endian -mcmodel=large -shared-intel ${fflag_RURI} ${cflag_RURI} ${flib_RURI} ${clib_RURI} ${static_RURI}"

else if(${machine} == "fugaku" && ${RSCUNIT} == "")then

    setenv SPACK_ROOT /vol0004/apps/oss/spack
    source /vol0004/apps/oss/spack/share/spack/setup-env.csh
    spack load ${netcdf_gcc}

    #set debug="-g -fcheck=bounds -fbacktrace"
    set debug=""
    set option="${fflag_gcc} ${cflag_gcc} ${flib_gcc} ${clib_gcc} ${static_gcc} -fno-range-check -mcmodel=medium -fPIC"

else if(${machine} == "fugaku" && ${RSCUNIT} == "rscunit_ft01")then

    setenv SPACK_ROOT /vol0004/apps/oss/spack
    source /vol0004/apps/oss/spack/share/spack/setup-env.csh
    spack load ${netcdf_fj}

    set option="${fflag_fj} ${cflag_fj} ${flib_fj} ${clib_fj} ${static_fj}"
    
endif

#--------------------------------------------------------------------
# Execute |
#--------------------------------------------------------------------

if(${machine} == "jss3" && ${RSCUNIT} == "")then

    ${FC} ${module} src/make_tsclim_ic.f90 -o make_tsclim_ic.out ${subroutine} ${debug} ${option}
    ./make_tsclim_ic.out > make_tsclim_ic.log

else if(${machine} == "jss3" && ${RSCUNIT} == "RURI")then

    ${FC} ${module} src/make_tsclim_ic.f90 -o make_tsclim_ic.out ${subroutine} ${debug} ${option} -parallel -qopenmp
    #Processor: 36 [/node]
    set THREAD=36
    jxsub <<EOF
#JX --bizcode R2402
#JX -L rscunit=RURI
#JX -L vnode=1
#JX -L vnode-core=${THREAD}
#JX -L vnode-mem=172080Mi
#JX -L elapse=01:00:00
#JX -N make_tsclim_ic
#JX -S

export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export OMP_NUM_THREAD=${THREAD}
export PARALLEL=${THREAD}

./make_tsclim_ic.out

rm -f make_tsclim_ic.out

EOF

else if(${machine} == "fugaku" && ${RSCUNIT} == "")then

    ${FC} ${module} src/make_tsclim_ic.f90 -o make_tsclim_ic.out ${subroutine} ${debug} ${option}
    ./make_tsclim_ic.out > make_tsclim_ic.log

else if(${machine} == "fugaku" && ${RSCUNIT} == "rscunit_ft01")then

    mpifrtpx ${module} src/make_tsclim_ic.f90 -o make_tsclim_ic.out ${subroutine} ${option}
    set THREAD=48
    
    pjsub <<EOF
#PJM -L "node=10"
#PJM -L "rscgrp=small" 
#PJM -L "elapse=01:00:00"
#PJM --mpi "proc=1"
#PJM --name "make_tsclim_ic"
#PJM -s
#

export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export OMP_NUM_THREAD=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n 1 ./make_tsclim_ic.out

rm -f make_tsclim_ic.out

EOF

endif

rm -f *.mod
