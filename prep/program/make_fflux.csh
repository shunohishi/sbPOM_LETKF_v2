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

#set machine="jss3"
set machine="fugaku"

#--------------------------------------------------------------------
# DATE |
#--------------------------------------------------------------------

set syr=${argv[1]};set smon=${argv[2]};set sday=${argv[3]}
set eyr=${argv[4]};set emon=${argv[5]};set eday=${argv[6]}

set shour=0
set ehour=0

set yyyys=`printf "%04d" ${syr}`
set mms=`printf "%02d" ${smon}`
set dds=`printf "%02d" ${sday}`

set yyyye=`printf "%04d" ${eyr}`
set mme=`printf "%02d" ${emon}`
set dde=`printf "%02d" ${eday}`

#---------------------------------------------------------------------
# Module & Subroutine |
#---------------------------------------------------------------------

set DAM="../module"
set module="src/mod_rmiss.f90 src/mod_gridinfo.f90 src/mod_parameter.f90 ${DAM}/mod_read_jra55do.f90 ${DAM}/mod_read_cama.f90"
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

${FC} ${module} src/make_fflux.f90 -o make_fflux.${yyyys}${mms}${dds}-${yyyye}${mme}${dde}.out ${subroutine} ${debug} ${option}

#-------------------------------------------------------------------------
# Execute |
#-------------------------------------------------------------------------

if(! -f make_fflux.${yyyys}${mms}${dds}-${yyyye}${mme}${dde}.out)then
    echo "***Error: Not found make_fflux.out"
    exit
else
    ./make_fflux.${yyyys}${mms}${dds}-${yyyye}${mme}${dde}.out ${syr} ${smon} ${sday} ${shour} ${eyr} ${emon} ${eday} ${ehour} > make_fflux.${yyyys}${mms}${dds}-${yyyye}${mme}${dde}.log &
endif
