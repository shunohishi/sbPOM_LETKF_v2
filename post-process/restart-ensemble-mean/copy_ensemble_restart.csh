#!/bin/csh
#========================================
# Copy data from JSPACE to local |
#========================================

set JSPACE=${argv[1]}
set dir=${argv[2]}
set letkf=${argv[3]}
set iyr=${argv[4]}
set imon=${argv[5]}
set iday=${argv[6]}
set nmem=${argv[7]}

#=========================================
# Date
#=========================================

set yyyy=`printf "%04d" ${iyr}`
set mm=`printf "%02d" ${imon}`
set dd=`printf "%02d" ${iday}`

#=========================================

#---Ensemble
@ imem = 1
while($imem <= $nmem)

    set mmmmm=`printf "%05d" ${imem}`

    set from_file="${JSPACE}/${letkf}/output/${mmmmm}/restart.${yyyy}${mm}${dd}.nc"
    set to_dir="${dir}/${letkf}/output/${mmmmm}/"

    if(! -d ${to_dir})then
	mkdir -p ${to_dir}
    endif

    if(-f ${from_file})then
	rsync -auv ${from_file} ${to_dir}
    else
	echo "***Error: Not found ${from_file}"
    	exit 99
    endif
        
    @ imem++

end

#---Ensemble mean
if(! -d ${dir}/${letkf}/output/mean)then
    mkdir -p ${dir}/${letkf}/output/mean
endif

cp -v ${dir}/${letkf}/output/00001/restart.${yyyy}${mm}${dd}.nc ${dir}/${letkf}/output/mean/restart.${yyyy}${mm}${dd}.nc
