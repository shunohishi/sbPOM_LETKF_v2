#!/bin/csh
#-----------------------------------------------------------
# Replace from ***1 to ***2
#----------------------------------------------------------

#debug
set debug1="-CB -traceback -g"
set debug2="-CB -traceback -g"

#option
set option1="-assume byterecl -convert big_endian -mcmodel=medium -shared-intel -I/opt/JX/oss/x86_64/netcdf-fortran/4.5.3/include -I/opt/JX/oss/x86_64/netcdf/4.7.3/include -I/opt/JX/oss/x86_64/hdf5/1.12.0/include -L/opt/JX/oss/x86_64/netcdf-fortran/4.5.3/lib -lnetcdff -L/opt/JX/oss/x86_64/hdf5/1.12.0/lib -L/opt/JX/oss/x86_64/netcdf/4.7.3/lib -lnetcdf -lnetcdf -lhdf5_hl -lhdf5 -lm -lz -lcurl"
set option2="-assume byterecl -convert big_endian -mcmodel=medium -shared-intel -I/opt/JX/oss/x86_64/netcdf-fortran/4.5.3/include -I/opt/JX/oss/x86_64/netcdf/4.7.3/include -I/opt/JX/oss/x86_64/hdf5/1.12.0/include -L/opt/JX/oss/x86_64/netcdf-fortran/4.5.3/lib -lnetcdff -L/opt/JX/oss/x86_64/hdf5/1.12.0/lib -L/opt/JX/oss/x86_64/netcdf/4.7.3/lib -lnetcdf -lnetcdf -lhdf5_hl -lhdf5 -lm -lz -lcurl"

#DATA
set data1="/data/R/R2402/DATA"
set data2="/data/R/R2402/DATA"

#Compiler
set compiler1="ifort"
set compiler2="ifort"

#-----------------------------------------------------------

sed -i -e "s|${debug1}|${debug2}|g" *.csh
sed -i -e "s|${option1}|${option2}|g" *.csh
sed -i -e "s|${data1}|${data2}|g" *.csh
sed -i -e "s|${compiler1}|${compiler2}|g" *.csh

