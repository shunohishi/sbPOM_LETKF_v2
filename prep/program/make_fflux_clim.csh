#!/bin/csh
#--------------------------------------------------------------------
# Make freshwater flux using GSMaP and CaMa-Flood |
#--------------------------------------------------------------------
#
# Modify time_setting
#
#--------------------------------------------------------------------
# Created by S.Ohishi 2018.08
# S.Ohishi 2019.02 Modify
# S.Ohishi 2023.07 Input date
# S.Ohishi 2024.12 Add Fugaku
#--------------------------------------------------------------------

set machine="jss3"
#set machine="fugaku"

#---------------------------------------------------------------------
# Module & Subroutine |
#---------------------------------------------------------------------

set module="src/mod_rmiss.f90 src/mod_gridinfo.f90 src/mod_parameter.f90 ../module/mod_check_netcdf.f90 ../module/mod_read_cama.f90"
set subroutine="src/sub_read_grid.f90 src/sub_cal_id.f90 src/sub_bilinear_interpolation.f90 src/sub_apply_fsm.f90 src/sub_distance.f90 src/sub_time.f90 src/sub_make_ncfile.f90"

#--------------------------------------------------------------------
# Compiler Option |
#--------------------------------------------------------------------

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

#-----------------------------------------------------------------------
# Compile |
#------------------------------------------------------------------------

${FC} ${module} src/make_fflux_clim.f90 -o make_fflux_clim.out ${subroutine} ${debug} ${option}

#-------------------------------------------------------------------------
# Execute |
#-------------------------------------------------------------------------

if(! -f make_fflux_clim.out)then
    echo "***Error: Not found make_fflux_clim.out"
    exit
else
    ./make_fflux_clim.out > make_fflux_clim.log &
endif
