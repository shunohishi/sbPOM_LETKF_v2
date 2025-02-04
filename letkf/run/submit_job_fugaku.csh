#!/bin/csh
#=== Argument==============================================================================================#

#Switch
set switch_iau=${argv[1]}
set switch_anal=${argv[2]}
set switch_fcst=${argv[3]}

#General
set CDIR=${argv[4]}
set WORKDIR=${argv[5]}
set OUTPUT=${argv[6]}
set INFO=${argv[7]}
set NMEM=${argv[8]}
set DT=${argv[9]}
set machine=${argv[10]}
set elapse_time=${argv[11]}

#OBS
set OBSDIR=${argv[12]}
set OBSFILE=${argv[13]}

#LETKF
set LETKFDIR=${argv[14]}
set LEXE=${argv[15]}

#sbPOM
set MODELDIR=${argv[16]}
set MODELDATADIR=${argv[17]}
set ATMDIR=${argv[18]}
set RIVDIR=${argv[19]}
set REGION=${argv[20]}
set PPROC=${argv[21]}
set PNODE=${argv[22]}
set PEXE=${argv[23]}
set TS_NUDGE=${argv[24]}
set TI_NUDGE=${argv[25]}
set SS_NUDGE=${argv[26]}
set SI_NUDGE=${argv[27]}

#NODE/PROC/THREAD
set NODE_TOTAL=${argv[28]}
set PPROC_TOTAL=${argv[29]}
set PTHREAD=${argv[30]}
set LTHREAD=${argv[31]}
set LPROC=${argv[32]}

#Date
set yyyymmdd=${argv[33]}
set ST=${argv[34]}
set IT=${argv[35]}
set LT=${argv[36]}

#==========================================================================================================#
# Group on Fugaku
#==========================================================================================================#

if(1 <= ${NODE_TOTAL} && ${NODE_TOTAL} <= 384)then
    set group="small"
else if(385 <= ${NODE_TOTAL} && ${NODE_TOTAL} <= 12288)then
    set group="large"
else
    echo "***Error: Large NODE SIZE = ${NODE_TOTAL}"
    exit
endif

#==========================================================================================================#

pjsub <<EOF
#PJM -L rscunit=rscunit_ft01
#PJM -L rscgrp=${group}
#PJM -L node=${NODE_TOTAL}
#PJM --mpi proc=${PPROC_TOTAL}
#PJM -L elapse=${elapse_time}
#PJM -L retention_state=0
#PJM -g ra000007
#PJM -x PJM_LLIO_GFSCACHE=/vol0004
#PJM --name sbPOM_LETKF_${REGION}_${yyyymmdd}
#PJM -S

IT=${IT}
LT=${LT}
DT=${DT}

export PLE_MPI_STD_EMPTYFILE="off"

while [ "\${IT}" -le "\${LT}" ]
do

date=\$(perl juldays.prl "\${IT}")

#=== ENSEMBLE FORECAST (run sbPOM) for IAU ==========#
if [ ${switch_iau} -eq 1 ]; then
    
    echo "==============================================================="
    echo "=== Start Ensemble Forecast for IAU at \${date}"
    echo "=== Ensemble size: ${NMEM}"
    echo "==============================================================="

    csh run.ensfcst_vcoord.csh 1 \
    ${CDIR} ${WORKDIR} ${OUTPUT} ${INFO} ${NMEM} ${machine} \
    ${MODELDIR} ${MODELDATADIR} ${ATMDIR} ${RIVDIR} ${REGION} ${PPROC} ${PNODE} ${PTHREAD} ${PEXE} 0 \
    ${TS_NUDGE} ${TI_NUDGE} ${SS_NUDGE} ${SI_NUDGE} \
    ${ST} "\${IT}" ${DT}
    #0: BUDGET
        
    rm -rf ${WORKDIR}/*
    cd ${CDIR}
        
    echo " ====================================================="
    echo " === End Ensemble forecast fOR IAU at \${date}"
    echo " === Ensemble size: ${NMEM}"
    echo " ====================================================="

fi #switch_iau 1
#=====================================================#

#=== LETKF ===========================================#
if [ ${switch_anal} -eq 1 ]; then

    echo "========================================================="
    echo "=== Start LETKF at \${date}"
    echo "=== Ensemble size: ${NMEM}"
    echo "========================================================="

    if [ ${switch_iau} -eq 0 ]; then
	csh run.letkf.csh 0 \
	${CDIR} ${WORKDIR} ${OUTPUT} ${INFO} ${NMEM} ${machine} ${OBSDIR} ${OBSFILE} \
	${LETKFDIR} ${LEXE} ${LPROC} ${LTHREAD} \
	${MODELDATADIR} "\${IT}"
    else
	csh run.letkf.csh 1 \
	${CDIR} ${WORKDIR} ${OUTPUT} ${INFO} ${NMEM} ${machine} ${OBSDIR} ${OBSFILE} \
	${LETKFDIR} ${LEXE} ${LPROC} ${LTHREAD} \
	${MODELDATADIR} "\${IT}"
    fi

    rm -rf ${WORKDIR}/*
    cd ${CDIR}
        
    echo "======================================================"
    echo "=== END LETKF DATA ASSIMILATION at \${date}        ===="
    echo "======================================================"

fi #switch_anal
#================================================================

#=== ENSEMBLE FORECAST (run sbPOM) ==============================
if [ ${switch_fcst} -eq 1 ]; then

    echo " ====================================================="
    echo " === START ENSEMBLE PREDICTION at \${date}          ==="
    echo " === NUMBER of ENSEMBLE: ${NMEM}                   ==="
    echo " ====================================================="

    echo "Submit sbPOM Shell" #------------------------------------------

    if [ ${switch_iau} -eq 0 ]; then
	csh run.ensfcst_vcoord.csh 0 \
	${CDIR} ${WORKDIR} ${OUTPUT} ${INFO} ${NMEM} ${machine} \
	${MODELDIR} ${MODELDATADIR} ${ATMDIR} ${RIVDIR} ${REGION} ${PPROC} ${PNODE} ${PTHREAD} ${PEXE} 1 \
	${TS_NUDGE} ${TI_NUDGE} ${SS_NUDGE} ${SI_NUDGE} \
	${ST} "\${IT}" ${DT}
	#1: BUDGET
    else
	csh run.ensfcst_vcoord.csh 2 \
	${CDIR} ${WORKDIR} ${OUTPUT} ${INFO} ${NMEM} ${machine} \
	${MODELDIR} ${MODELDATADIR} ${ATMDIR} ${RIVDIR} ${REGION} ${PPROC} ${PNODE} ${PTHREAD} ${PEXE} 1 \
	${TS_NUDGE} ${TI_NUDGE} ${SS_NUDGE} ${SI_NUDGE} \
	${ST} "\${IT}" ${DT}
	#1: BUDGET
    fi

    rm -rf ${WORKDIR}/*
    cd ${CDIR}

    echo " ==================================================="
    echo " === END ENSEMBLE PREDICTION at \${date}          ==="
    echo " === NUMBER of ENSEMBLE: ${NMEM}                 ==="
    echo " ==================================================="
        
fi #switch_fcst
#===================================================================

IT=\$((IT+DT))

done

echo \${PJM_JOBID} > ${WORKDIR}/FINISHED_sbPOM_LETKF

EOF

