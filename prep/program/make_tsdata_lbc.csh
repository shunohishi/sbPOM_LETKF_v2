#!/bin/csh
#--------------------------------------------------------------
# Make daily tsdata & lbc from Parent Model |
#--------------------------------------------------------------
#
# tsdata: for nuding near boudary conditions
# lbc: Lateral boudary condition
#
# Check filename in 
# subroutine read_parent_grid
# subroutine read_parent      in mod_read_parent.f90
#
#--------------------------------------------------------------
# Created by
# S.Ohishi 2018.09
#--------------------------------------------------------------
# Elapsed time
# wpac --> scs: 15min./month
#--------------------------------------------------------------

#set debug="-CB -traceback -g"
set debug=""
set option="-assume byterecl -convert big_endian -mcmodel=medium -shared-intel -I/opt/JX/oss/x86_64/intel/2022.3.1/netcdf-fortran/4.5.2/include -I/opt/JX/oss/x86_64/intel/2022.3.1/hdf5/1.12.0/include -I/opt/JX/oss/x86_64/intel/2022.3.1/netcdf/4.7.3/include -L/opt/JX/oss/x86_64/intel/2022.3.1/netcdf-fortran/4.5.2/lib -lnetcdff -L/opt/JX/oss/x86_64/intel/2022.3.1/hdf5/1.12.0/lib -L/opt/JX/oss/x86_64/intel/2022.3.1/netcdf/4.7.3/lib -lnetcdf -lnetcdf -lm"

set module="src/mod_rmiss.f90 src/mod_parameter.f90 src/mod_gridinfo.f90 src/mod_read_parent.f90"
set subroutine="src/sub_read_grid.f90 src/sub_distance.f90 src/sub_cal_id.f90 src/sub_fillvalue.f90 src/sub_bilinear_interpolation.f90 src/sub_apply_fsm.f90 src/sub_time.f90"

ifort $module src/make_tsdata_lbc.f90 -o make_tsdata_lbc.out $subroutine $debug $option
./make_tsdata_lbc.out > make_tsdata_lbc.log
rm -f make_tsdata_lbc.out
rm -f *.mod
