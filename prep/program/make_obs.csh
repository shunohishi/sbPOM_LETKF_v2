#!/bin/csh
#
# 1. Modify "module setting" in make_obs.f90
#
#-----------------------------------------------------------------
#
# Created by
# S.Ohishi 2018.09
#_________________________________________________________________


set machine="jss3"
#set machine="fugaku"

#=================================================================
# Make dir |
#=================================================================

if(! -d ../obs)then
    mkdir ../obs
endif

#=================================================================
# Compile option |
#=================================================================

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
    set option="${fflag_gcc} ${cflag_gcc} ${flib_gcc} ${clib_gcc} ${static_gcc} -fno-range-check -mcmodel=medium -fPIC"

endif

#==================================================================
# Module & Subroutine |
#==================================================================

set DAM="../module"
set module="src/mod_rmiss.f90 src/mod_parameter.f90 src/mod_gridinfo.f90 ${DAM}/mod_read_himawari.f90 ${DAM}/mod_read_amsre.f90 ${DAM}/mod_read_windsat.f90 ${DAM}/mod_read_amsr2.f90 ${DAM}/mod_read_smap.f90 ${DAM}/mod_read_smos.f90 ${DAM}/mod_read_cmems.f90 ${DAM}/mod_read_gtspp.f90 ${DAM}/mod_read_en4.f90 ${DAM}/mod_read_gcomc.f90"

set subroutine="src/sub_read_grid.f90 src/sub_detect_nearshore.f90 src/sub_time.f90 src/sub_distance.f90 src/sub_cal_id.f90 src/sub_bilinear_interpolation.f90 src/sub_read_ssh.f90 src/sub_fillvalue.f90 src/sub_prepare_sst.f90 src/sub_prepare_sss.f90 src/sub_prepare_mdot.f90 src/sub_prepare_ssh.f90 src/sub_prepare_ts.f90 src/sub_prepare_ssuv.f90 src/sub_make_ncfile.f90 src/sub_write_obs.f90"

rm -f make_obs.out

${FC} ${module} src/make_obs.f90 -o make_obs.out ${subroutine} ${debug} ${option}
if(-f make_obs.out)then
    ./make_obs.out > make_obs.log
else
    echo "***Error: Not found make_obs.out"
endif

rm -f make_obs.out
rm -f *.mod
