#!/bin/csh
#---------------------------------------------------------------
# Date |
#---------------------------------------------------------------
#
###### KEO
# T & S: 2004.06.16-
# U & V: 2005.05.30-
###### Papa
# T & S: 2007.06.08-
# U & V: 2007.06.08-
#---------------------------------------------------------------

set sdate=(2004 6 16)
set edate=(2020 12 31)

#---------------------------------------------------------------
# Validation using KEO and Papa buoys |
#---------------------------------------------------------------

#---Machine
#set machine="jss3"
#set machine="fugaku"
set machine="rc"

#---Partition (only for R-CCS Cloud)
set partition="genoa" #Execute on r340/genoa
#set partition="fx700"  #Execute on fx700

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

else if(${machine} == "rc")then

    #set debug="-g -fcheck=bounds -fbacktrace"
    set debug=""
    set fflag=`nf-config --fflags`
    set flib=`nf-config --flibs`
    set clib=`nc-config --libs`
    set option="${fflag} ${flib} ${clib} -ffree-line-length-none"
    
endif

#---------------------------------------------------------------
# Subroutine & Module |
#---------------------------------------------------------------

set module="../module/mod_rmiss.f90  ../module/mod_julian.f90  ../module/mod_read_ocs.f90 ../module/mod_gridinfo.f90 ../module/mod_read_lora_v20.f90 ../module/mod_read_bran2020.f90 ../module/mod_read_fora_np60.f90 ../module/mod_read_glorys010.f90 ../module/mod_read_glorys025.f90 ../module/mod_read_jcope_fgo.f90 mod_setting.f90 mod_make_ncfile.f90 mod_io.f90"
set subroutine="sub_get_id.f90 sub_convert.f90"

#---------------------------------------------------------------
# Compile |
#---------------------------------------------------------------

rm -f make_data.out
${FC} ${module} main_make_data.f90 ${subroutine} ${option} ${debug} -o make_data.out

#---------------------------------------------------------------
# Execution |
#---------------------------------------------------------------

#---Make dir
foreach var(t s u v)
    if(! -d dat/keo/${var}) mkdir -p dat/keo/${var}
    if(! -d dat/papa/${var}) mkdir -p dat/papa/${var}
end

#---Check
if(! -f make_data.out)then
    echo "***Error: Compile main_make_data.f90"
    exit
endif

#---Execulte
@ nb = 2 #1: KEO, 2: Papa
@ nvar = 4 #1:T, 2:S, 3:U, 4:V

if(${machine} == "rc")then

    sbatch -p ${partition} --job-name=make_ocs_data submit_job_make_data.sh ${sdate} ${edate}

else

    @ ib = 1

    while($ib <= $nb)

	@ ivar = 1

	if(${ib} == 1)then
	    set bname="keo"
	else if(${ib} == 2)then
	    set bname="papa"
	endif

	while($ivar <= $nvar)

	    if(${ivar} == 1)then
		set varname="t"
	    else if(${ivar} == 2)then
		set varname="s"
	    else if(${ivar} == 3)then
		set varname="u"
	    else if(${ivar} == 4)then
		set varname="v"
	    endif
    
	    ./make_data.out ${sdate} ${edate} ${ib} ${ivar} > make_data_${bname}_${varname}.log

	    @ ivar++
	    
	end
	@ ib++
    end

endif
        
rm -f *.mod
