#!bin/csh
#========================================================#
# DATE |
#========================================================#

set sdate=(2001 1 1) #Spin-up start date
set idate=(2003 1 1) #Assimilation start date
set edate=(2003 1 1) #Assimilation end date

#========================================================#
# Switch(1:execute, other: skip) |
#========================================================#

set switch_prepare_ens=1 #Prepare initial ensemble member
set switch_iau=1         #sbPOM forecast to conduct the IAU (0:Intermittent , 1: IAU)
set switch_anal=1        #LETKF analysis
set switch_fcst=1        #sbPOM forecast 

#========================================================#
# GENERAL |
#========================================================#

set DIR=/data/R/R2402/ohishi/TEST
set PDIR=${DIR}/prep       #Pre-process DIR
set LDIR=${DIR}/letkf      #LETKF DIR
set CDIR=`pwd`             #Current DIR
set WORKDIR=${LDIR}/work   #Work DIR
set OUTPUT=${LDIR}/output  #Output DIR
set INFO=${LDIR}/run/info  #LOG INFORMATION
set NMEM=10                #Ensembel member
#set NMEM=128               #Ensembel member
set NPROC=48               #Proceccor at 1 node
set DT=1                   #Delta T [unit: day]
set CT=1                   #Computation window day [unit: day]
set machine="jss3"
#set machine="fugaku"
#set elapse_time="00:10:00"
set elapse_time="01:00:00"

#=========================================================#
# OBSERVATION |
#=========================================================#

set OBSDIR=${PDIR}/obs #Observational DIR for Assimilation
set OBSFILE=obs        #Observation filename (ex. ${OBSFILE}${yyyy}${mm}${dd}.nc)

#========================================================#
# LETKF |
#========================================================#

set LETKFDIR=${LDIR}/run  #LETKF DIR
set LEXE=letkf.exe        #LETKF EXE File

#========================================================#
# sbPOM |
#========================================================#

set MODELDIR=${DIR}/pom90-ens/run          #sbPOM DIR
set MODELDATADIR=${PDIR}/in                #grid,tsclim,ic,lbc, tsdata DIR
set ATMDIR=${PDIR}/in                      #ATM DIR
set RIVDIR=${PDIR}/in                      #RIVER DIR
set MODELOUTPUTDIR=${DIR}/pom90-ens/output #sbPOM output DIR
set REGION=test                            #Forecast filename
set PPROC=8                                #Number of Processor at 1 simulation
set PNODE=2                                #Number of Node at 1 simulation
set PEXE=pom.exe                           #EXE file
set BUDGET=1                               #T and S budget term (1:on, 0:off)
set TS_NUDGE=0.                            #SST nuding [day] (*0 --> not execute)
set TI_NUDGE=0.                            #T nuding [day]
set SS_NUDGE=30.                           #SSS nuding [day]
set SI_NUDGE=0.                            #S nuding [day]

#=======================================================#
# NODE/PROC/THREAD |
#=======================================================#

#ALL
@ NODE_TOTAL = ${PNODE} * ${NMEM}
@ PROC_TOTAL = ${NODE_TOTAL} * ${NPROC}
echo "Total Node: ${NODE_TOTAL}, Total Processor: ${PROC_TOTAL}"

#Simulation
@ PPROC_TOTAL = ${PPROC} * ${NMEM}
@ PTHREAD = ${PNODE} * ${NPROC} / ${PPROC}
echo "Simulation (EACH):  ${PNODE} node (${PPROC} processor and ${PTHREAD} thread)"
echo "Simulation (TOTAL): ${NODE_TOTAL} node"

#LETKF
@ LTHREAD = ${PTHREAD}
@ LPROC = ${PROC_TOTAL} / ${LTHREAD}

echo "Ensemble size: ${NMEM}"
echo "Total process: ${PROC_TOTAL}"
echo "LETKF: ${LPROC} processor and ${LTHREAD} thread"

if(1 <= ${NODE_TOTAL} && ${NODE_TOTAL} <= 384)then
    set group="small"
else if(385 <= ${NODE_TOTAL} && ${NODE_TOTAL} <= 12288)then
    set group="large"
else
    echo "***Error: NODE_TOTAL: ${NODE_TOTAL}"
    exit
endif
    
#=========================================================#
# Compile Option |
#=========================================================#

#NetCDF
if(${machine} == "jss3")then

    set option="-Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign ${fflag_SORA} ${cflag_SORA} ${flib_SORA} ${clib_SORA} ${static_SORA}"

else if(${machine} == "fugaku")then

    setenv SPACK_ROOT /vol0004/apps/oss/spack
    source /vol0004/apps/oss/spack/share/spack/setup-env.csh
    spack load ${netcdf_fj}

    set option="-Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign ${fflag_fj} ${cflag_fj} ${flib_fj} ${clib_fj} ${static_fj}"

endif

#========================================================#
# Make dir
#========================================================#

###OUTPUT
echo "Check OUTPUT Directory"
if(! -d ${OUTPUT}) mkdir -p ${OUTPUT}

@ IMEM = 1
while(${IMEM} <= ${NMEM})
    set MMMMM=`printf "%05d" ${IMEM}`
    if(! -d ${OUTPUT}/${MMMMM}) mkdir -p ${OUTPUT}/${MMMMM}
    @ IMEM++
end

foreach dir("mean" "sprd" "eens")
    if(! -d ${OUTPUT}/${dir}) mkdir -p ${OUTPUT}/${dir}
end

###WORKDIR
echo "Make Work directory"
rm -rf ${WORKDIR} && mkdir -p ${WORKDIR}

###LETKF
echo "Make LETKF EXE & Log directory"
cd ${LETKFDIR}
make
if(! -f ${LEXE})then
    echo "***Error: ${LEXE} not found"
    exit 99
endif

###sbPOM
echo "Make sbPOM EXE"
cd ${MODELDIR}
make
if(! -f ${PEXE})then
    echo "***Error: ${PEXE} not found"
    exit 99
endif

cd ${CDIR}
if(! -d ${INFO}) mkdir -p ${INFO}

#===========================================#
# DATE
#===========================================#

#Julian day: Spin-up start date
@ ST = `perl caldays.prl ${sdate[1]} ${sdate[2]} ${sdate[3]}` 
set yyyys=`printf "%04d" ${sdate[1]}`;set mms=`printf "%02d" ${sdate[2]}`;set dds=`printf "%02d" ${sdate[3]}`

#Julian day: Assimilation start date
@ IT = `perl caldays.prl ${idate[1]} ${idate[2]} ${idate[3]}` 

#Julian day: Assimilation end date
@ ET = `perl caldays.prl ${edate[1]} ${edate[2]} ${edate[3]}` 
set yyyye=`printf "%04d" ${edate[1]}`;set mme=`printf "%02d" ${edate[2]}`;set dde=`printf "%02d" ${edate[3]}`

#===========================================#
#Prepare ENSENBLE MEMBER |
#===========================================#

if(${switch_prepare_ens} == 1)then
    echo "Prepare ensemble member"
    csh prepare_ens.csh ${IT} ${NMEM} ${MODELOUTPUTDIR} ${OUTPUT}
    if($? == 99)then
       echo "***Error: Prepare ensemble member"
       exit 99
    endif
endif

#============================================#
# Make mesp_ens.out |
#============================================#

rm -f mesp_ens_mpi.out
mpifrtpx ens/julian.f90 ens/mesp_ens_mpi.f90 -o mesp_ens_mpi.out ${option}
if(! -f mesp_ens_mpi.out)then
    echo "***Error: mesp_ens_mpi.out not found"
    exit 99
endif

#============================================#
# Check Pre-process file |
#============================================#

if(${switch_iau} == 1 || ${switch_fcst} == 1)then

    cd ${LETKFDIR}
    csh check_prep.csh ${ST} ${NMEM} ${PDIR} ${ATMDIR} ${RIVDIR} ${MODELDATADIR}
    if($? == 99)then
	echo "***Error: Check pre-process file"
	exit 99
    endif
    
endif

#============================================#
# Main Loop |
#============================================#

echo "======================================================================================="
echo " Start: ${idate} - ${edate}"
echo "======================================================================================="
echo "ssssss  b       PPPPPP          MM     MM       L      EEEEEE TTTTTTT  K     K  FFFFFF "
echo "s       b       P    P          M  M  M M       L      E         T     K  K     F      "
echo "ssssss  bbbbbb  PPPPPP  oooooo  M   M   M  ===  L      EEEEEE    T     K        FFFF   "
echo "     s  b    b  P       o    o  M       M       L      E         T     K  K     F      "
echo "ssssss  bbbbbb  P       oooooo  M       M       LLLLLL EEEEEE    T     K     K  F      "
echo "======================================================================================="

while(${IT} <= ${ET} - ${CT} + 1)

    #IT Date
    set date=`perl juldays.prl ${IT}`
    set yyyy=`printf "%04d" ${date[1]}`;set mm=`printf "%02d" ${date[2]}`; set dd=`printf "%02d" ${date[3]}`
    set yyyymmdd=${yyyy}${mm}${dd}

    #LT
    @ LT = ${IT} + ${CT} - 1

    if(${machine} == "jss3")then
	csh submit_job_jss3.csh ${switch_iau} ${switch_anal} ${switch_fcst} \
	${CDIR} ${WORKDIR} ${OUTPUT} ${INFO} ${NMEM} ${DT} ${machine} ${elapse_time} \
	${OBSDIR} ${OBSFILE} ${LETKFDIR} ${LEXE} \
	${MODELDIR} ${MODELDATADIR} ${ATMDIR} ${RIVDIR} ${REGION} ${PPROC} ${PNODE} ${PEXE} \
	${TS_NUDGE} ${TI_NUDGE} ${SS_NUDGE} ${SI_NUDGE} \
	${NODE_TOTAL} ${PPROC_TOTAL} ${PTHREAD} ${LTHREAD} ${LPROC} \
	${yyyymmdd} ${ST} ${IT} ${LT}
    else if(${machine} == "fugaku")then
	csh submit_job_fugaku.csh ${switch_iau} ${switch_anal} ${switch_fcst} \
	${CDIR} ${WORKDIR} ${OUTPUT} ${INFO} ${NMEM} ${DT} ${machine} ${elapse_time} \
	${OBSDIR} ${OBSFILE} ${LETKFDIR} ${LEXE} \
	${MODELDIR} ${MODELDATADIR} ${ATMDIR} ${RIVDIR} ${REGION} ${PPROC} ${PNODE} ${PEXE} \
	${TS_NUDGE} ${TI_NUDGE} ${SS_NUDGE} ${SI_NUDGE} \
	${NODE_TOTAL} ${PPROC_TOTAL} ${PTHREAD} ${LTHREAD} ${LPROC} \
	${yyyymmdd} ${ST} ${IT} ${LT}
    endif
    
    echo "Wait for the sbPOM-LETKF ${yyyymmdd} job to complete"
    @ isec = 0
    @ imin = 0
    @ iint = 10
    set JOBID="00000000"
    while(${isec} >= 0)

	sleep ${iint}s
	@ isec = ${isec} + ${iint}
	@ imin = ${isec} / 60
	@ imod = ${isec} % 60

	if(-f ${WORKDIR}/FINISHED_sbPOM_LETKF)then
	    set JOBID=`head -n 1 ${WORKDIR}/FINISHED_sbPOM_LETKF`
	endif
	
	set FIN1=`find ${WORKDIR} -name FINISHED_sbPOM_LETKF | wc -l`
	set FIN2=`find ${CDIR} -name  sbPOM_LETKF_${REGION}_${yyyymmdd}.${JOBID}.stats | wc -l`
	echo "sbPOM_LETKF FINISHED [${FIN1}/1] [${FIN2}/1] ; ${imin}:${imod} elapsed"

	#BREAK
	if(${FIN1} == 1 && ${FIN2} == 1)then
	    echo "End sbPOM-LETKF job"
	    rm -f ${WORKDIR}/FINISHED_sbPOM_LETKF
	    break
	endif

    end
    
    #--- Ensebmle mean/spread--------------------------------------------
    if(${switch_fcst} == 1)then
	echo "Post process: Ensemble mean/spread"
        csh submit_post.csh ${LDIR} ${OUTPUT} ${INFO} ${REGION} ${CT} ${yyyys} ${mms} ${dds} ${yyyy} ${mm} ${dd} \
        ${NPROC} ${NMEM} ${BUDGET} 2 ${machine}
        # 2: Remove restart and daily-mean data
    endif
        
    rm -rf ${WORKDIR}/*
    cd ${CDIR}

    #--- INFO -----------------------------------------------------------
    foreach ext(out err stats)
	if(-s sbPOM_LETKF_${REGION}_${yyyymmdd}.${JOBID}.${ext})then #No empty
	    (mv -f sbPOM_LETKF_${REGION}_${yyyymmdd}.${JOBID}.${ext} ${INFO}/sbPOM_LETKF_${REGION}_${yyyymmdd}.${ext} &)
	else
	    (rm -f sbPOM_LETKF_${REGION}_${yyyymmdd}.${JOBID}.${ext} &)
	endif
    end
    #---------------------------------------------------------------------
    
    @ IT = ${IT} + ${CT}

end #IT

rm -f mesp_ens.out
wait

echo "======================================================================================="
echo "End: ${idate} - ${edate}"
echo "======================================================================================="
echo "ssssss  b       PPPPPP          MM     MM       L      EEEEEE TTTTTTT  K     K  FFFFFF "
echo "s       b       P    P          M  M  M M       L      E         T     K  K     F      "
echo "ssssss  bbbbbb  PPPPPP  oooooo  M   M   M  ===  L      EEEEEE    T     K        FFFF   "
echo "     s  b    b  P       o    o  M       M       L      E         T     K  K     F      "
echo "ssssss  bbbbbb  P       oooooo  M       M       LLLLLL EEEEEE    T     K     K  F      "
echo "======================================================================================="

exit 0
