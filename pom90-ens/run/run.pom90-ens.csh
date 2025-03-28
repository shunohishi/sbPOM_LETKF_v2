#!/bin/csh
#=======================================================================
# Setting for ensemble simulation |
#=======================================================================

set REGION=test
set EXE=pom.exe
set DIR=/data/R/R2402/ohishi/TEST/pom90-ens
set PDIR=/data/R/R2402/ohishi/TEST/prep/
set NCDIR=/data/R/R2402/ohishi/TEST/prep/in
set ATMDIR=${NCDIR}
set RIVDIR=${NCDIR}
set CURDIR=`pwd`

set NPROC=48   #Total processor for one node (JSS3: 48)
set EPROC=8    #Processor for each simulation
#set NMEM=10    #Ensemble size
set NMEM=32    #Ensemble size
#set NMEM=128   #Ensemble size
set ENODE=2    #Node for each simulation
@ THREAD = ${NPROC} * ${ENODE} / ${EPROC}
@ NODE = ${ENODE} * ${NMEM}
@ PROC = ${EPROC} * ${NMEM}

set BUDGET=0     #Heat/Salinity budget convervation [1:On, 0:Off]
set TS_NUDGE=30. #T surface nudging [day]
set TI_NUDGE=30. #T internal
set SS_NUDGE=30. #S surface
set SI_NUDGE=30. #S internal
set RM_ENS=1     #Remove ensemble member 1: On, 0: Off
set machine="jss3"
#set machine="fugaku"

set sdate=(2001 1) #start time
set idate=(2001 1) #initial time
set edate=(2002 12) #end time

#=========================================================================
# Compile Option |
#=========================================================================

#NetCDF
if(${machine} == "jss3")then

    set option="-Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign ${fflag_SORA} ${cflag_SORA} ${flib_SORA} ${clib_SORA} ${static_SORA}"

else if(${machine} == "fugaku")then

    setenv SPACK_ROOT /vol0004/apps/oss/spack
    source /vol0004/apps/oss/spack/share/spack/setup-env.csh
    spack load ${netcdf_fj}

    set option="-Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign ${fflag_fj} ${cflag_fj} ${flib_fj} ${clib_fj} ${static_fj}"

endif

#===========================================================================
# Compile |
#===========================================================================

cd ${DIR}/run/
echo "Compile"
make
if(! -f ${EXE})then
    echo "***Error: Not found pom.exe" 
    exit
endif

mpifrtpx ens/julian.f90 ens/mesp_ens_mpi.f90 -o mesp_ens_mpi.out ${option}
if(! -f mesp_ens_mpi.out)then
    echo "***Error: Not found mesp_ens_mpi.out"
    exit
endif

#============================================================================
# Main loop |
#============================================================================

@ iyr = ${idate[1]}
@ imon = ${idate[2]}

while(${iyr} <= ${edate[1]})

    if(${iyr} == ${edate[1]})then
	@ emon = ${edate[2]}
    else
	@ emon = 12
    endif

    while(${imon} <= ${emon})

    #-----------------------------------------
    # Setting dat information
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

    set yyyy_s=`printf "%04d" ${sdate[1]}`;set mm_s=`printf "%02d" ${sdate[2]}`;set dd_s="01"
    set yyyy_n1=`printf "%04d" ${iyr_n1}`;set mm_n1=`printf "%02d" ${imon_n1}`
    set yyyy=`printf "%04d" ${iyr}`;set mm=`printf "%02d" ${imon}`;set dd="01"
    set yyyy_p1=`printf "%04d" ${iyr_p1}`;set mm_p1=`printf "%02d" ${imon_p1}` 

    if(${imon} == 2 && ${iyr} % 4 == 0)then
	set nday=29
    else if(${imon} == 2 && ${iyr} % 4 != 0)then
	set nday=28
    else if(${imon} == 4 || ${imon} == 6 || ${imon} == 9 || ${imon} == 11)then
	set nday=30
    else
	set nday=31
    endif
    #set nday=1
    
    #-----------------------------------
    echo "Start ${yyyy}${mm}"
    #-----------------------------------

    echo "Make Workdir"
    set WORKDIR=${DIR}/output/${yyyy}${mm}
    set WORKDIR_n1=${DIR}/output/${yyyy_n1}${mm_n1}
    rm -rf ${WORKDIR}

    if(! -d ${WORKDIR}/mean) mkdir -p ${WORKDIR}/mean
    if(! -d ${WORKDIR}/sprd) mkdir -p ${WORKDIR}/sprd
    if(! -d ${WORKDIR}/eens) mkdir -p ${WORKDIR}/eens
    
    @ IMEM = 1
    while(${IMEM} <= ${NMEM})

	set MMMMM=`printf "%05d" ${IMEM}`

	#--------------------------------------------------------
	if(${IMEM} == 1) echo "Make in/out DIR"
	#--------------------------------------------------------
	foreach dir(in out)
	    if( ! -d ${WORKDIR}/${MMMMM}/${dir})then
		mkdir -p ${WORKDIR}/${MMMMM}/${dir}
	    endif
	end

	#---------------------------------------------------------
	if(${IMEM} == 1) echo "Make restart file link"
	#---------------------------------------------------------
	if(${iyr} == ${sdate[1]} && ${imon} == ${sdate[2]})then
	    @ switch_rst = 0
	else if(-f ${WORKDIR_n1}/${MMMMM}/out/restart.nc)then
	    @ switch_rst = 1
	    ln -s ${WORKDIR_n1}/${MMMMM}/out/restart.nc ${WORKDIR}/${MMMMM}/in/restart.nc
	else
	    echo "***Error: Not found restart.nc"
	    exit
	endif

	@ IMEM++

    end #IMEM

    #-------------------------------------------------
    echo "Make link"
    #-------------------------------------------------

    csh ${DIR}/run/make_link.csh ${REGION} ${EXE} ${DIR} ${NCDIR} ${ATMDIR} ${RIVDIR} ${NMEM} ${WORKDIR} ${mm_s}

    #--------------------------------------------------
    echo "Make namelist"
    #--------------------------------------------------
    
    csh ${DIR}/run/make_namelist.csh ${REGION} ${DIR} ${NMEM} ${BUDGET} \
    ${TS_NUDGE} ${TI_NUDGE} ${SS_NUDGE} ${SI_NUDGE} \
    ${WORKDIR} ${yyyy_s} ${mm_s} ${nday} ${switch_rst}
    
    #---------------------------------------------------
    echo "Make vcoordfile"
    #---------------------------------------------------
    
    csh ${DIR}/run/make_vcoord.csh ${WORKDIR} ${NMEM} ${ENODE} ${EPROC}
    if($? == 99)then
	echo "***Error: make_vcoord.csh"
    	exit
    endif

    #-------------------------------------------------------
    echo "Submit Job: Ensemble simulation"
    #-------------------------------------------------------

    if(${nday} == 1)then
	set elapse_time="00:10:00"
    else
	set elapse_time="03:00:00"
    endif
    
    cd ${WORKDIR}
    echo "Ensemble size: ${NMEM}"
    echo "Each Node: ${ENODE}, Each Processor: ${EPROC}, Each Thread: ${THREAD}"
    echo "Total Node: ${NODE}, Total Processor: ${PROC}"
    
    csh ${DIR}/run/submit_job.csh ${REGION} ${DIR} ${EXE} \
    ${EPROC} ${NMEM} ${THREAD} ${NODE} ${PROC} ${machine} \
    ${yyyy} ${mm} ${WORKDIR} ${elapse_time}
    
    #------------------------------------------------------
    # Timer |
    #------------------------------------------------------
    
    echo "Wait for ensemble forecast to finish"
    @ int = 10 #Interval [sec.]
    @ isec = 0
    while(${isec} >= 0)
	sleep ${int}s
	@ isec = ${isec} + ${int}
	@ sec = ${isec} % 60
	@ min = ${isec} / 60
	set FIN_NUM=`find ${WORKDIR} -name FINISHED_sbPOM | wc -l`
	echo "sbPOM FINISHED [ ${FIN_NUM} / ${NMEM} ], ${min}:${sec} passed"
	if(${FIN_NUM} == ${NMEM})then
	    echo "End ${yyyy}${mm}"
	    break
	endif
    end

    #------------------------------------------------------
    echo "Remove link"
    #------------------------------------------------------
    
    @ IMEM = 1
    while(${IMEM} <= ${NMEM})
    
    	set MMMMM=`printf "%05d" ${IMEM}`
	if(-f ${WORKDIR}/${MMMMM}/pom.exe) unlink ${WORKDIR}/${MMMMM}/pom.exe
	rm -f ${WORKDIR}/${MMMMM}/pom.nml
	
	foreach filename(grid tsclim ic tsdata lbc)
	    if(-f ${WORKDIR}/${MMMMM}/in/${REGION}.${filename}.nc) unlink ${WORKDIR}/${MMMMM}/in/${REGION}.${filename}.nc
	end
	
	foreach dirname(atm river)
	    if(-f ${WORKDIR}/${MMMMM}/in/${dirname}) unlink ${WORKDIR}/${MMMMM}/in/${dirname}
	end

	if(-f ${WORKDIR}/${MMMMM}/in/restart.nc) unlink ${WORKDIR}/${MMMMM}/in/restart.nc
	
	@ IMEM++
	
    end

    #-----------------------------------------------------
    echo "Post process"
    #-----------------------------------------------------
    cd ${DIR}/run/
    csh submit_post.csh ${REGION} ${DIR} ${NPROC} ${NMEM} ${BUDGET} ${RM_ENS} ${machine} \
    ${yyyy_s} ${mm_s} ${dd_s} ${yyyy} ${mm} ${dd} ${nday} ${WORKDIR}
    #-----------------------------------------------------
    
    @ imon++

    end #imon

@ iyr++
@ imon=1

end #iyr

#=========================================================================#
echo "#####ALL END######"
#=========================================================================#
