#!/bin/csh
#--------------------------------------------------------
# Argument |
#--------------------------------------------------------

set dir   = ${argv[1]}
set letkf = ${argv[2]}
set iyr  = ${argv[3]}
set imon = ${argv[4]}
set iday = ${argv[5]}
set nmem = ${argv[6]}
set machine = ${argv[7]}
set RSCUNIT = ${argv[8]}

set yyyy=`printf "%04d" ${iyr}`
set mm=`printf "%02d" ${imon}`
set dd=`printf "%02d" ${iday}`
set mmmmm=`printf "%05d" ${nmem}`
set elapse_time="00:20:00"

#--------------------------------------------------------
# PROCESSOR |
#--------------------------------------------------------

set PROC=${nmem}  #Total processor

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

set module="../../prep/program/src/mod_gridinfo.f90 module/mod_varname.f90"
set subroutine="subroutine/sub_argument.f90 subroutine/sub_io.f90 subroutine/sub_ensemble_mean.f90"

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

rm -f make_restart_ensemble_mean.*.out make_restart_ensemble_mean.*.err make_restart_ensemble_mean.*.stats
rm -f stdout.make_restart_ensemble_mean* stderr.make_restart_ensemble_mean*
rm -f make_restart_ensemble_mean.out FINISHED

#-------------------------------------------------------------
# Compile |
#-------------------------------------------------------------

if(${RSCUNIT} == "SORA" || ${RSCUNIT} == "rscunit_ft01")then
    mpifrtpx ${module} main.f90 -o make_restart_ensemble_mean.out ${subroutine} ${option}
else if(${RSCUNIT} == "RURI")then
    mpiifort ${module} main.f90 -o make_restart_ensemble_mean.out ${subroutine} -qopenmp ${option} ${debug}
endif
    
if(! -f make_restart_ensemble_mean.out)then
    echo "***Error: Nout found make_restart_ensemble_mean.out"
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
#JX -L elapse=${elapse_time}
#JX -N make_restart_ensemble_mean
#JX -S

export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export PLE_MPI_STD_EMPTYFILE="off"
export OMP_NUM_THREADS=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${PROC} -stdout stdout.make_restart_ensemble_mean -stderr stderr.make_restart_ensemble_mean ./make_restart_ensemble_mean.out ${dir} ${letkf} ${yyyy} ${mm} ${dd} ${mmmmm} && echo \${PJM_JOBID} > FINISHED

rm -f make_restart_ensemble_mean.out

EOF

else if(${RSCUNIT} == "RURI")then

jxsub <<EOF
#JX --bizcode R2402
#JX -L rscunit=RURI
#JX -L vnode=${NODE}
#JX -L vnode-core=${EPROC}
#JX -L vnode-mem=172080Mi
#JX -L elapse=${elapse_time}
#JX -N make_restart_ensemble_mean
#JX -S

export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export PLE_MPI_STD_EMPTYFILE="off"
export OMP_NUM_THREADS=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${PROC} -ppn ${EPROC} ./make_restart_ensemble_mean.out ${dir} ${letkf} ${yyyy} ${mm} ${dd} ${mmmmm} && echo \${PJM_JOBID} > FINISHED

rm -f make_restart_ensemble_mean.out

EOF

else if(${RSCUNIT} == "rscunit_ft01")then

pjsub <<EOF
#PJM -L rscunit=rscunit_ft01
#PJM -L rscgrp=${rscgrp}
#PJM -L node=${NODE}
#PJM --mpi proc=${PROC}
#PJM -L elapse=${elapse_time}
#PJM -L retention_state=0
#PJM -g ra000007
#PJM -x PJM_LLIO_GFSCACHE=/vol0004
#PJM --name make_restart_ensemble_mean
#PJM -S

export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export PLE_MPI_STD_EMPTYFILE="off"
export OMP_NUM_THREADS=${THREAD}
export PARALLEL=${THREAD}

mpiexec -n ${PROC} -stdout-proc stdout.make_restart_ensemble_mean -stderr-proc stderr.make_restart_ensemble_mean ./make_restart_ensemble_mean.out ${dir} ${letkf} ${yyyy} ${mm} ${dd} ${mmmmm} && echo \${PJM_JOBID} > FINISHED

rm -f make_restart_ensemble_mean.out

EOF

endif

rm -f *.mod

#----------------------------------------------------------------------------------
# Timer
#----------------------------------------------------------------------------------

@ isec = 0
@ imin = 0
@ iint = 10
set JOBID="00000000"
while(${isec} >= 0)

    #TIMER
    sleep ${iint}s
    @ isec = ${isec} + ${iint}
    @ imin = ${isec} / 60
    @ imod = ${isec} % 60

    set FIN_NUM=`find . -name FINISHED | wc -l`
    set STATS_NUM=`find . -name make_restart_ensemble_mean.${JOBID}.stats | wc -l`
    echo "[${STATS_NUM}/1] ; ${imin}:${imod} elapsed"

    #BREAK
    if(${FIN_NUM} == 1)then
	set JOBID=`head -n 1 FINISHED`
	rm -f FINISHED
    endif
    
    if(${STATS_NUM} == 1)then

    	foreach ext(out err stats)
	    if(-f make_restart_ensemble_mean.${JOBID}.${ext})then
		mv make_restart_ensemble_mean.${JOBID}.${ext} info/make_restart_ensemble_mean.${yyyy}${mm}${dd}.${ext}
    	    endif
	end

	foreach std(stdout stderr)
	    if(-f ${std}.make_restart_ensemble_mean)then
		mv ${std}.make_restart_ensemble_mean info/${std}.make_restart_ensemble_mean.${yyyy}${mm}${dd}
	    endif
	end
	
        echo "End job"
        break

    endif

end #isec
