#!/bin/csh
#---------------------------------------------------------------
# Validation using surface current from drifter buoys |
#---------------------------------------------------------------

set machine="jss3"
#set machine="fugaku"

#---------------------------------------------------------------
# Option |
#---------------------------------------------------------------

if(${machine} == "jss3")then

    #set debug="-CB -traceback -g"
    set debug=""
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

set module="../module/mod_julian.f90 ../module/mod_rmiss.f90 ../module/mod_gridinfo.f90 ../module/mod_read_lora.f90"
set subroutine=""

#---------------------------------------------------------------
# Compile |
#---------------------------------------------------------------

${FC} ${module} main.f90 ${subroutine} ${option}

#---------------------------------------------------------------
# Output directory |
#---------------------------------------------------------------

if(! -d dat)then
    mkdir dat
endif

#---------------------------------------------------------------
# Execution |
#---------------------------------------------------------------

./a.out
rm -f *.mod

