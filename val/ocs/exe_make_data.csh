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

#set sdate=(2004 6 16)
#set edate=(2018 12 31)
set sdate=(2023 6 1)
set edate=(2023 6 30)

#---------------------------------------------------------------
# Validation using KEO and Papa buoys |
#---------------------------------------------------------------

set machine="jss3"
#set machine="fugaku"

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
    
endif

#---------------------------------------------------------------
# Subroutine & Module |
#---------------------------------------------------------------

set module="../module/mod_rmiss.f90  ../module/mod_julian.f90  ../module/mod_read_ocs.f90 ../module/mod_gridinfo.f90 ../module/mod_read_lora_v20.f90 ../module/mod_read_glorys025.f90 mod_setting.f90 mod_make_ncfile.f90 mod_io.f90"
set subroutine="sub_get_id.f90 sub_convert.f90"

#---------------------------------------------------------------
# Compile |
#---------------------------------------------------------------

rm -f make_data.out
${FC} ${module} main_make_data.f90 ${subroutine} ${option} ${debug} -o make_data.out

#---------------------------------------------------------------
# Execution |
#---------------------------------------------------------------

#---Make dir
foreach var(t s u v)
    if(! -d dat/keo/${var}) mkdir -p dat/keo/${var}
    if(! -d dat/papa/${var}) mkdir -p dat/papa/${var}
end

#---Check
if(! -f make_data.out)then
    echo "***Error: Compile main_make_data.f90"
    exit
endif

#---Execulte
./make_data.out ${sdate} ${edate}

rm -f *.mod

