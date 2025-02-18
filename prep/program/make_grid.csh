#!/bin/csh
#
# 1. Modify "module setting" in make_grid.f90
# 2. run with iswitch=1
# 3. Modify fsm.dat
# 4. run with iswitch =2
# 
#---------------------------------------------------------------------------
#
# Created by 
# S.Ohishi 2018.07
# S.Ohishi 2024.12 Add Fugaku version
#___________________________________________________________________________

#set machine="jss3"
set machine="fugaku"

#---------------------------------------------------------------------------
# Make dir |
#---------------------------------------------------------------------------

if(! -d ../in)then
    mkdir ../in
endif

#---------------------------------------------------------------------------
# Compile Option |
#---------------------------------------------------------------------------

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

#---------------------------------------------------------------------------
# Subroutine & Module |
#---------------------------------------------------------------------------

set DAM="../module"
set module="src/mod_rmiss.f90 src/mod_parameter.f90 ${DAM}/mod_read_etopo1.f90"
set subroutine="src/sub_cal_id.f90 src/sub_bilinear_interpolation.f90 src/sub_smooth_topo.f90 src/sub_distance.f90 src/sub_gauss_filter.f90 src/sub_make_ncfile.f90"

#---------------------------------------------------------------------------
# Compile |
#---------------------------------------------------------------------------

rm -f make_grid.out

${FC} $module src/make_grid.f90 -o make_grid.out $subroutine $debug $option

#---------------------------------------------------------------------------
# Execute |
#---------------------------------------------------------------------------

if(-f make_grid.out)then
    ./make_grid.out > make_grid.log
else
    echo "***Error: Not found make_grid.out"
endif

rm -f make_grid.out
rm -f *.mod
