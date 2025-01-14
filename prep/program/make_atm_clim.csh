#!/bin/csh
#----------------------------------------------------------------
# Make Climatological atmospheric boudary condition using JRA55 |
#----------------------------------------------------------------
#
# Modify "module time_setting" in make_atm.f90
#
#-------------------------------------------------------------
# Created by S.Ohishi 2018.08
# S.Ohishi 2021.12 add JRA55do
# S.Ohishi 2024.12 add Fugaku
# S.Ohishi 2024.12 add Climatology
#-------------------------------------------------------------

set machine="jss3"
#set machine="fugaku"

#-------------------------------------------------------------
# Module & Subroutine |
#-------------------------------------------------------------

set DAM="../module"
set module="src/mod_rmiss.f90 src/mod_parameter.f90 src/mod_gridinfo.f90 ${DAM}/mod_read_jra55.f90 ${DAM}/mod_read_jra55do.f90"
set subroutine="src/sub_read_grid.f90 src/sub_distance.f90 src/sub_cal_id.f90 src/sub_apply_fsm.f90 src/sub_fillvalue.f90 src/sub_bilinear_interpolation.f90 src/sub_time.f90 src/sub_make_ncfile.f90"

#-------------------------------------------------------------
# Compiler Option |
#-------------------------------------------------------------

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

#--------------------------------------------------------------
# Compile & Execute |
#--------------------------------------------------------------

${FC} ${module} src/make_atm_clim.f90 -o make_atm_clim.out ${subroutine} ${debug} ${option}

if(! -f make_atm_clim.out)then
    echo "***Error: Not found make_atm_clim.out"
    exit
else
    ./make_atm_clim.out > make_atm_clim.log && rm -f make_atm_clim.out *.mod &
endif
