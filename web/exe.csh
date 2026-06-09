#!/bin/csh

set machine="jss3"
#set machine="fugaku"

#-----------------------------------------------------------------
# Option |
#-----------------------------------------------------------------

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

#-----------------------------------------------------------------
# Module and subroutine |
#-----------------------------------------------------------------

set module="../prep/program/src/mod_gridinfo.f90 module/mod_julian.f90 module/mod_setting.f90 module/mod_rmiss.f90 module/mod_read_lora_v20.f90"
set subroutine="subroutine/sub_id.f90 subroutine/sub_io.f90 subroutine/sub_linear_interpolation.f90 subroutine/sub_make_depth.f90"

#-----------------------------------------------------------------
# Compile | 
#-----------------------------------------------------------------

rm -f a.out
${FC} ${module} main.f90 ${subroutine} ${option}

#------------------------------------------------------------------
# Execute |
#------------------------------------------------------------------

if(! -f a.out)then
    echo "***Error: Compile main.f90"
    exit
endif

if(! -d dat) mkdir dat

./a.out

rm -f *.mod


