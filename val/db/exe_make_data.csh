#!/bin/csh
#---------------------------------------------------------------
# Make data in observation space |
#---------------------------------------------------------------

###Machine
#set machine="jss3"
#set machine="fugaku"
set machine="rc"

###Partition (only for R-CCS Cloud)
#set partition="r340"  #Execute on r340
#set partition="genoa" #Execute on r340/genoa
set partition="fx700"  #Execute on fx700

###Period
set sdate=(2003 1)
set edate=(2023 12)

#---------------------------------------------------------------
# Option |
#---------------------------------------------------------------

if(${machine} == "jss3")then

    set option="-Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign ${fflag_SORA} ${cflag_SORA} ${flib_SORA} ${clib_SORA} ${static_SORA}"

else if(${machine} == "fugaku")then

    setenv SPACK_ROOT /vol0004/apps/oss/spack
    source /vol0004/apps/oss/spack/share/spack/setup-env.csh
    
    spack load ${netcdf_fj}
    set option="-Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign ${fflag_fj} ${cflag_fj} ${flib_fj} ${clib_fj} ${static_fj}"

else if(${machine} == "rc")then

    set fflag=`nf-config --fflags`
    set flib=`nf-config --flibs`
    set clib=`nc-config --libs`
    set option="${fflag} ${flib} ${clib} -ffree-line-length-none"
    
endif

#---------------------------------------------------------------
# Subroutine & Module |
#---------------------------------------------------------------

set module="../module/mod_julian.f90 ../module/mod_rmiss.f90 ../module/mod_read_db.f90 ../module/mod_gridinfo.f90 ../module/mod_read_lora_v20.f90 ../module/mod_read_bran2020.f90 ../module/mod_read_jcope_fgo.f90 ../module/mod_read_glorys010.f90 ../module/mod_read_glorys025.f90 mod_setting.f90 mod_make_ncfile.f90 mod_io.f90"
set subroutine="sub_bilinear_interpolation.f90 sub_cal_id.f90 sub_check_data_location.f90"

#---------------------------------------------------------------
# Compile |
#---------------------------------------------------------------

rm -f make_data.out

if(${machine} == "rc")then
    gfortran ${module} main_make_data.f90 ${subroutine} ${option} -o make_data.out
else
    mpifrtpx ${module} main_make_data.f90 ${subroutine} ${option} -o make_data.out
endif

#---------------------------------------------------------------
# Execution |
#---------------------------------------------------------------

#---Check
if(! -f make_data.out)then
    echo "***Error: Compile main_make_data.f90"
    exit
endif

#---Annual
@ iyr=${sdate[1]}
@ imon=${sdate[2]}

while($iyr <= ${edate[1]})

    if($iyr == $edate[1])then
	@ emon = $edate[2]
    else
	@ emon = 12
    endif

    while($imon <= $emon) 

	if($imon == 2 && $iyr % 4 == 0)then
	    set nday=29
	else if($imon == 2)then
	    set nday=28
	else if($imon == 4 || $imon == 6 || $imon == 9 || $imon == 11)then
	    set nday=30
	else
	    set nday=31
	endif
    
	set yyyy=`printf "%04d" ${iyr}`
	set mm=`printf "%02d" ${imon}`
	if(! -d dat/${yyyy}${mm})then
	    mkdir -p dat/${yyyy}${mm}
	endif

	if(${machine} == "rc")then
	    sbatch -p ${partition} --job-name=make_data_${yyyy}${mm} submit_job_rc.sh ${iyr} ${imon} 1 ${iyr} ${imon} ${nday} ${yyyy} ${mm}
	else
	    csh submit_job.csh ${machine} ${iyr} ${imon} 1 ${iyr} ${imon} ${nday} ${yyyy} ${mm}
	endif
	    
	@ imon++
	
    end

    @ iyr++
    @ imon = 1
	
end    

rm -f *.mod
