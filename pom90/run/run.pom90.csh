#!/bin/csh
#=======================================================================
# Setting for single simulation|
#=======================================================================
set REGION=test  #Region name
set DIR=/data/R/R2402/ohishi/TEST/pom90         #POM Directory
set PDIR=/data/R/R2402/ohishi/TEST/prep/program #Pre-process source file directory
set NCDIR=/data/R/R2402/ohishi/TEST/prep/in     #Pre-process netcdf directory
set ATMDIR=${NCDIR}
set RIVDIR=${NCDIR}
set CURDIR=`pwd`
set NPROC=48 #Total Proceccor for 1 node (JSS3 & Fugaku: 48)
set NODE=2   #Node used 
set PROC=8   #Processor used
@ THREAD = ${NPROC} * ${NODE} / ${PROC} #THREAD used
set BUDGET=0      #Heat Salinity budget convervation (1:On, 0: Off)
set TS_NUDGE=0.  #T surface nudging [day]
set TI_NUDGE=0.  #T internal
set SS_NUDGE=30.  #S surface
set SI_NUDGE=0.  #S internal
set machine="jss3"
#set machine="fugaku"

set sdate=(1993 1) #start time
set idate=(1993 1) #initial time
set edate=(1993 1) #end time

#========================================================================

#-----------------------------------
# Execution file |
#-----------------------------------

#---POM
set EXE=pom.exe  #Execute file
echo "Compile"
make
if(! -f ${EXE})then
    echo "***Error: Not found ${EXE}"
    exit
endif

#-----------------------------------
# Start: Loop |
#-----------------------------------

@ iyr = ${idate[1]}
@ imon = ${idate[2]}

while(${iyr} <= ${edate[1]})

if(${iyr} == ${edate[1]})then
    @ emon = ${edate[2]}
else
    @ emon = 12
endif

    while(${imon} <= ${emon})

    cd ${DIR}/run/
    
    #-----------------------------------------
    # Setting date information |
    #-----------------------------------------

    if(${imon} == 12)then
	@ imon_p1 = 1
	@ iyr_p1 = ${iyr} + 1
    else
	@ imon_p1 = ${imon} + 1
	@ iyr_p1 = ${iyr}
    endif

    if(${imon} == 1)then
	@ imon_n1 = 12
	@ iyr_n1 = ${iyr} - 1
    else
	@ imon_n1 = ${imon} - 1
	@ iyr_n1 = ${iyr}
    endif

    set syyyy=`printf "%04d" ${sdate[1]}`;set smm=`printf "%02d" ${sdate[2]}`     
    set yyyy_n1=`printf "%04d" ${iyr_n1}`;set mm_n1=`printf "%02d" ${imon_n1}`
    set yyyy=`printf "%04d" ${iyr}`;set mm=`printf "%02d" ${imon}`
    set yyyy_p1=`printf "%04d" ${iyr_p1}`;set mm_p1=`printf "%02d" ${imon_p1}` 

    if(${imon} == 2 && ${iyr} % 4 == 0)then
	set nday=29.
    else if(${imon} == 2 && ${iyr} % 4 != 0)then
	set nday=28.
    else if(${imon} == 4 || ${imon} == 6 || ${imon} == 9 || ${imon} == 11)then
	set nday=30.
    else
	set nday=31.
    endif
    set nday=1. #DEBUG
    
    #-----------------------------------
    echo "Start ${yyyy}${mm}"
    #-----------------------------------

    echo "Make Workdir"
    set WORKDIR=${DIR}/output/${yyyy}${mm}
    set WORKDIR_n1=${DIR}/output/${yyyy_n1}${mm_n1}
    rm -rf ${WORKDIR}
    mkdir -p ${WORKDIR}
    mkdir -p ${WORKDIR}/in
    mkdir -p ${WORKDIR}/out
    
    #-------------------------------------
    echo "Make link"
    #-------------------------------------

    cd ${DIR}/run/
    csh make_link.csh ${DIR} ${NCDIR} ${ATMDIR} ${RIVDIR} ${WORKDIR} ${REGION} ${smm} ${EXE}

    #-------------------------------------
    # Restart file 
    #-------------------------------------

    if(${iyr} == ${sdate[1]} && ${imon} == ${sdate[2]})then
	@ switch_rst = 0
    else if(-f ${WORKDIR_n1}/out/restart.nc)then
	@ switch_rst = 1
	ln -s ${WORKDIR_n1}/out/restart.nc ${WORKDIR}/in/restart.nc
    else
	echo "***Error: Not found restart.nc"
	exit
    endif

    #--------------------------------------
    # Namelist
    #--------------------------------------

    cd ${DIR}/run/    
    csh make_namelist.csh ${REGION} ${WORKDIR} \
    ${BUDGET} ${TS_NUDGE} ${TI_NUDGE} ${SS_NUDGE} ${SI_NUDGE} \
    ${syyyy} ${smm} ${switch_rst} ${nday}

    #---------------------------------------
    # Submite job 
    #---------------------------------------

    if(${nday} == 1.)then
	set elapse_time="00:10:00"
    else
	set elapse_time="03:00:00"
    endif
    
    cd ${DIR}/run/    
    csh submit_job.csh ${REGION} ${WORKDIR} ${NODE} ${PROC} ${THREAD} ${machine} ${EXE} \
    ${yyyy} ${mm} ${elapse_time}

    #---------------------------------------
    # TIMER
    #---------------------------------------

    cd ${WORKDIR}
    echo "Wait for Job finish"
    @ int = 10 #Interval [sec.]
    @ isec = 0
    while(${isec} >= 0)
	sleep ${int}s
	@ isec = ${isec} + ${int}
	@ sec = ${isec} % 60
	@ min = ${isec} / 60
	echo "sbPOM: ${min}:${sec} elapsed" 
	if(-f FINISHED)then
	    echo "End ${yyyy}${mm}"
	    break
	endif
    end

    #---Post process
    foreach filename(grid ic lbc tsclim tsdata)
	if(-f ${WORKDIR}/in/${REGION}.${filename}.nc)then
	    unlink ${WORKDIR}/in/${REGION}.${filename}.nc
	endif
    end
    foreach dirname(atm river)
	if(-d ${WORKDIR}/in/${dirname})then
	    unlink ${WORKDIR}/in/${dirname}
	endif
    end
    if(-f ${WORKDIR_n1}/out/restart.nc) unlink ${WORKDIR}/in/restart.nc
    cd ${CURDIR}

    @ imon++

    end #imon

@ iyr++
@ imon=1

end #iyr

#-----------------------------------
# End: Loop |
#-----------------------------------
echo "#####ALL END######"
