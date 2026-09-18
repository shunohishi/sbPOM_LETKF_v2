#!/bin/csh
#---------------------------------------------------------------
# Machine |
#---------------------------------------------------------------

#---Machine
set machine="jss3"
#set machine="fugaku"
set machine="rc"

#---Partition (only for R-CCS Cloud)
#set partition="r340"  #Execute on r340
set partition="genoa"  #Execute on r340/genoa
#set partition="fx700" #Execute on fx700

#---------------------------------------------------------------
# Option |
#---------------------------------------------------------------

if(${machine} == "jss3")then

    set debug="-CB -traceback -g"
    #set debug=""
    set option="-assume byterecl -convert big_endian -mcmodel=medium -shared-intel ${fflag_RURI} ${cflag_RURI} ${flib_RURI} ${clib_RURI} ${static_RURI}"

else if(${machine} == "fugaku")then

    setenv SPACK_ROOT /vol0004/apps/oss/spack
    source /vol0004/apps/oss/spack/share/spack/setup-env.csh
    spack load ${netcdf_gcc}

    #set debug="-g -fcheck=bounds -fbacktrace"
    set debug=""
    set option="${fflag_gcc} ${cflag_gcc} ${flib_gcc} ${clib_gcc} ${static_gcc} -fno-range-check"

else if(${machine} == "rc")then

    #set debug="-g -fbacktrace -fcheck=all"
    set debug=""
    set fflag=`nf-config --fflags`
    set flib=`nf-config --flibs`
    set clib=`nc-config --libs`
    set option="${fflag} ${flib} ${clib} -ffree-line-length-none"

else

    echo "***Error: machine or RSCUNIT"
    exit
    
endif

#---------------------------------------------------------------
# Subroutine & Module |
#---------------------------------------------------------------

set module="../module/mod_rmiss.f90 ../module/mod_read_duacs.f90 ../../prep/module/mod_read_woa18_month.f90 ../../prep/module/mod_read_woa18_season_annual.f90 ../module/mod_gridinfo.f90 ../module/mod_read_lora_v20.f90 ../module/mod_read_glorys025.f90 mod_setting.f90 mod_make_ncfile.f90 mod_io.f90"
set subroutine="sub_cal_id.f90 sub_bilinear_interpolation.f90"

#---------------------------------------------------------------
# Compile |
#---------------------------------------------------------------

rm -f est_bias.out
${FC} ${module} main_est_bias.f90 ${subroutine} ${option} ${debug} -o est_bias.out

#---------------------------------------------------------------
# Execution |
#---------------------------------------------------------------

#---Check
if(! -f est_bias.out)then
    echo "***Error: Compile main_est_bias.f90"
    exit
endif

#---Execute
if(${machine} == "rc")then
    sbatch -p ${partition} --job-name=est_bias submit_job_est_bias.sh
else
    ./est_bias.out > est_bias.log &
endif

rm -f *.mod
