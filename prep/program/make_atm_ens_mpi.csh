#!/bin/csh
#-------------------------------------------------------------
# Make Ensemble atmospheric boudary condition using JRA55 |
#-------------------------------------------------------------
#
# Modify "module setting" in make_atm_ens.f90
#
#-------------------------------------------------------------
# Created by S.Ohishi 2021.11
# S.Ohishi 2023.07 modify for input start and end date
# S.Ohishi 2024.12 Add Fugaku
#-------------------------------------------------------------

#set machine="jss3"
set machine="fugaku"

#-------------------------------------------------------------

set syr=${argv[1]};set smon=${argv[2]};set sday=${argv[3]}
set eyr=${argv[4]};set emon=${argv[5]};set eday=${argv[6]}

echo ${syr} ${smon} ${sday} 0 > atm_date.dat
echo ${eyr} ${emon} ${eday} 0 >> atm_date.dat

#--------------------------------------------------------------
# Elapse time |
#--------------------------------------------------------------

if(${sday} == ${eday})then
    set elapse_time="00:20:00"
else
    set elapse_time="00:05:00"
endif

#-------------------------------------------------------------
# Ensemble size, Processor, Node, Thread |
#-------------------------------------------------------------

set NENS=${argv[7]}  #Ensemble size => #Processor
set NPROC=${argv[8]} #Total processor for 1 node
@ THREAD = 12

if((${NENS} * ${THREAD}) % ${NPROC} == 0)then
    @ NODE = (${NENS} * ${THREAD}) / ${NPROC}
else
    @ NODE = (${NENS} * ${THREAD}) / ${NPROC} + 1
endif

#-------------------------------------------------------------
# Module & Subroutine |
#-------------------------------------------------------------

set DAM="../module"
set module="src/mod_rmiss.f90 src/mod_parameter.f90 src/mod_gridinfo.f90 ${DAM}/mod_read_jra55.f90 ${DAM}/mod_read_jra55do.f90"

set subroutine="src/sub_read_grid.f90 src/sub_ensemble.f90 src/sub_distance.f90 src/sub_cal_id.f90 src/sub_apply_fsm.f90 src/sub_fillvalue.f90 src/sub_bilinear_interpolation.f90 src/sub_time.f90 src/sub_mpi.f90 src/sub_make_ncfile.f90"

#---------------------------------------------------------------
# Compier option |
#---------------------------------------------------------------

if(${machine} == "jss3")then

    set option="-Kfast -Kopenmp -Kparallel -Nalloc_assign ${fflag_SORA} ${cflag_SORA} ${flib_SORA} ${clib_SORA} ${static_SORA}"
    
else if(${machine} == "fugaku")then

    setenv SPACK_ROOT /vol0004/apps/oss/spack
    source /vol0004/apps/oss/spack/share/spack/setup-env.csh
    spack load ${netcdf_fj}

    set option="-Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign ${fflag_fj} ${cflag_fj} ${flib_fj} ${clib_fj} ${static_fj}"

endif

#------------------------------------------------------------------
# Remove standard output & error etc. |
#------------------------------------------------------------------
    
rm -f FINISHED_ATM
rm -f make_atm_ens_mpi.*.out make_atm_ens_mpi.*.err make_atm_ens_mpi.*.stats
rm -f stdout.make_atm_ens_mpi* stderr.make_atm_ens_mpi*

#-------------------------------------------------------------------
# Compile |
#-------------------------------------------------------------------

mpifrtpx ${module} src/make_atm_ens_mpi.f90 -o make_atm_ens_mpi.out ${subroutine} ${option} && rm -f *.mod
if(! -f make_atm_ens_mpi.out)then
    echo "***Error: Not found make_atm_ens_mpi.out"
    exit
endif

#-------------------------------------------------------------------
# Execute |
#-------------------------------------------------------------------

echo "Processor:${NENS}, Node:${NODE}"

if(${machine} == "jss3")then

jxsub <<EOF
#JX --bizcode R2402
#JX -L rscunit=SORA
#JX -L node=${NODE}
#JX --mpi proc=${NENS}
#JX -L node-mem=29184Mi
#JX -L elapse=${elapse_time}
#JX -N make_atm_ens_mpi
#JX -S

export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export OMP_NUM_THREAD=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${NENS} -stdout stdout.make_atm_ens_mpi -stderr stderr.make_atm_ens_mpi ./make_atm_ens_mpi.out &
echo \${PJM_JOBID} > JOBID_ATM &
wait

echo FINISHED > FINISHED_ATM
rm -f make_atm_ens_mpi.out

EOF

else if(${machine} == "fugaku")then

pjsub <<EOF
#PJM -L rscunit=rscunit_ft01
#PJM -L rscgrp=small
#PJM -L node=${NODE}
#PJM --mpi proc=${NENS}
#PJM -L elapse=${elapse_time}
#PJM --name make_atm_ens_mpi
#PJM -S

export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export OMP_NUM_THREAD=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${NENS} -stdout-proc stdout.make_atm_ens_mpi -stderr-proc stderr.make_atm_ens_mpi ./make_atm_ens_mpi.out

rm -f make_atm_ens_mpi.out

EOF

endif

#----------------------------------------------------------------------
# Timer |
#-----------------------------------------------------------------------

echo "Wait for starting ensemble ATM"
@ int = 10 #Interval [sec.]
@ isec = 0
while(${isec} >= 0)
    sleep ${int}s
    @ isec = ${isec} + ${int}
    @ sec = ${isec} % 60
    @ min = ${isec} / 60
    echo "Ensemble ATM: ${min}:${sec} passed"
    if(-f FINISHED_ATM)then
	echo "End Ensemble ATM"
	rm -f JOBID_ATM FINISHED_ATM
	break
    endif
end
