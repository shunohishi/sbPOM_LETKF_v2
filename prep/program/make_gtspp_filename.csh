#!/bin/csh
#
# Created by
# S.Ohishi 2025.01
#
#_______________________________________________________________________________

#set machine="jss3"
set machine="fugaku"

set idate=(2003 1)
set edate=(2003 1)

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

#------------------------------------------------------------------------------
# Module & Subroutine
#------------------------------------------------------------------------------

set DAM="../module"
set module="${DAM}/mod_rmiss.f90 ${DAM}/mod_read_gtspp.f90"
set subroutine="src/sub_time.f90"
${FC} ${module} src/make_gtspp_filename.f90 -o make_gtspp_filename.out ${subroutine} ${option} ${debug}

#------------------------------------------------------------------------------

rm -f make_gtspp_filename.log

@ iyr=${idate[1]}
@ imon=${idate[2]}

while($iyr <= $edate[1])

    if($iyr == $edate[1])then
	@ emon = ${edate[2]}
    else
	@ emon = 12
    endif

    while($imon <= $emon)

	set yyyy=`printf "%04d" ${iyr}`
	set mm=`printf "%02d" ${imon}`
	
	./make_gtspp_filename.out ${iyr} ${imon} >> make_gtspp_filename${yyyy}${mm}.log &

	@ imon++
    
    end

    @ iyr++
    @ imon = 1

end
wait

rm -f *.mod
