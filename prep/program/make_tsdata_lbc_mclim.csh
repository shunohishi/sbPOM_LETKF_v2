#!/bin/csh
#---------------------------------------------------------
# Make monthly climatology tsdata & lbc using SODA |
#---------------------------------------------------------
#
# tsdata: for nuding near boudary conditions
# lbc: Lateral boudary condition
# 
# Monthly climatology in syr-eyr in make_tsdata_lbc_mclim.f90
# *syr:Start year, eyr:End year
# 
# SODA version 3.7.2
# 1980.01-2015.12
#
#--------------------------------------------------------
# Created by S.Ohishi 2018.08
#
# S.Ohishi 2024.12 Add Fugaku 
#--------------------------------------------------------

#set machine="jss3"
set machine="fugaku"

#---------------------------------------------------------
# Module & Subroutine |
#---------------------------------------------------------

set DAM="../module"
set module="src/mod_rmiss.f90 src/mod_parameter.f90 src/mod_gridinfo.f90 ${DAM}/mod_read_soda.f90"
set subroutine="src/sub_read_grid.f90 src/sub_distance.f90 src/sub_cal_id.f90 src/sub_fillvalue.f90 src/sub_bilinear_interpolation.f90 src/sub_apply_fsm.f90 src/sub_make_ncfile.f90"

#---------------------------------------------------------
# Compiler Option |
#---------------------------------------------------------

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
    set option="${fflag_gcc} ${cflag_gcc} ${flib_gcc} ${clib_gcc} ${static_gcc} -fno-range-check -mcmodel=medium -fPIC"

endif

#--------------------------------------------------------
# Execute |
#---------------------------------------------------------

${FC} ${module} src/make_tsdata_lbc_mclim.f90 -o make_tsdata_lbc_mclim.out ${subroutine} ${debug} ${option}
./make_tsdata_lbc_mclim.out > make_tsdata_lbc_mclim.log
rm -f make_tsdata_lbc_mclim.out
rm -f *.mod
