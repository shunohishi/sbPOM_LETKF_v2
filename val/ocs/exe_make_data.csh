#!/bin/csh
#---------------------------------------------------------------
# Date |
#---------------------------------------------------------------
#
###### KEO
# T & S: 2004.06.16-
# U & V: 2005.05.30-
###### Papa
# T & S: 2007.06.08-
# U & V: 2007.06.08-
#---------------------------------------------------------------

set sdate=(2004 6)
set edate=(2020 12)

#---------------------------------------------------------------
# Validation using KEO and Papa buoys |
#---------------------------------------------------------------

#---Machine
#set machine="jss3"
#set machine="fugaku"
set machine="rc"

#---Partition (only for R-CCS Cloud)
set partition="genoa" #Execute on r340/genoa
#set partition="fx700"  #Execute on fx700

#---------------------------------------------------------------
# Option |
#---------------------------------------------------------------

if(${machine} == "jss3")then

    #set debug="-CB -traceback -g"
    set debug=""
    set option="-assume byterecl -convert big_endian -mcmodel=medium -shared-intel -qopenmp ${fflag_RURI} ${cflag_RURI} ${flib_RURI} ${clib_RURI} ${static_RURI}"
    set ncpu="72" #Login node
    set thread="1" 
    
else if(${machine} == "fugaku")then

    setenv SPACK_ROOT /vol0004/apps/oss/spack
    source /vol0004/apps/oss/spack/share/spack/setup-env.csh
    spack load ${netcdf_gcc}

    #set debug="-g -fcheck=bounds -fbacktrace"
    set debug=""
    set option="${fflag_gcc} ${cflag_gcc} ${flib_gcc} ${clib_gcc} ${static_gcc} -fopenmp -fno-range-check"
    set ncpu="64" #Login node
    set thread="1"
    
else if(${machine} == "rc")then

    #set debug="-g -fcheck=bounds -fbacktrace"
    set debug=""
    set fflag=`nf-config --fflags`
    set flib=`nf-config --flibs`
    set clib=`nc-config --libs`
    set option="${fflag} ${flib} ${clib} -fopenmp -ffree-line-length-none"

    if(${partition} == "genoa")then
	@ ncpu = 96
	@ nproc = 24
	@ thread = ${ncpu} / ${nproc}
    else if(${partition} == "fx700")then
	@ ncpu = 48
	@ nproc = 4
	@ thread = ${ncpu} / ${nproc}
    else
	@ ncpu = 1
        @ nproc = 1
	@ thread = 1	
    endif
    
endif

#---------------------------------------------------------------
# Subroutine & Module |
#---------------------------------------------------------------

set module="../module/mod_rmiss.f90  ../module/mod_julian.f90  ../module/mod_read_ocs.f90 ../module/mod_gridinfo.f90 ../module/mod_read_lora_v20.f90 ../module/mod_read_bran2020.f90 ../module/mod_read_fora_np60.f90 ../module/mod_read_glorys010.f90 ../module/mod_read_glorys025.f90 ../module/mod_read_jcope_fgo.f90 mod_setting.f90 mod_make_ncfile.f90 mod_io.f90"
set subroutine="sub_get_id.f90 sub_convert.f90"

#---------------------------------------------------------------
# Compile |
#---------------------------------------------------------------

#---Compile
rm -f make_data.out
${FC} ${module} main_make_data.f90 ${subroutine} ${option} ${debug} -o make_data.out

#---Check
if(! -f make_data.out)then
    echo "***Error: Compile main_make_data.f90"
    exit
endif

#---------------------------------------------------------------
# Execution |
#---------------------------------------------------------------

#---Make dir
foreach var(t s u v)
    if(! -d dat/keo/${var}) mkdir -p dat/keo/${var}
    if(! -d dat/papa/${var}) mkdir -p dat/papa/${var}
end

#---Execulte
@ iyr=${sdate[1]}
@ imon=${sdate[2]}
@ ijob=0
set args="${thread}"

while($iyr <= ${edate[1]})

    if(${iyr} == ${edate[1]})then
	@ emon = ${edate[2]}
    else
	@ emon = 12
    endif

    while(${imon} <= ${emon})
    
	if(${imon} == 2 && ${iyr} % 4 == 0)then
            set nday=29
        else if(${imon} == 2)then
            set nday=28
        else if(${imon} == 4 || ${imon} == 6 || ${imon} == 9 || ${imon} == 11)then
            set nday=30
        else
            set nday=31
        endif

	set yyyy=`printf "%04d" ${iyr}`
        set mm=`printf "%02d" ${imon}`

	if(${machine} == "rc")then

	    set args="${args} ${iyr} ${imon} ${nday} ${yyyy} ${mm}"
	    @ ijob++

	    if(${ijob} == ${nproc})then
		sbatch -p ${partition} --cpus-per-task=${ncpu} --exclusive --job-name=make_ocs_data submit_job_make_data.sh ${args}
		set args="${thread}"
		@ ijob=0
	    endif

	endif
	    
	@ imon++

    end

    if(${machine} == "jss3" || ${machine} == "fugaku")then
	./make_data.out ${iyr} 1 1 ${iyr} 12 31 > ${yyyy}.log &
    endif
    
    @ iyr++
    @ imon = 1
    
end

if(${machine} == "rc" && ${ijob} > 0)then
    sbatch -p ${partition} --cpus-per-task=${ncpu} --exclusive --job-name=make_ocs_data submit_job_make_data.sh ${args}
endif
        
rm -f *.mod
