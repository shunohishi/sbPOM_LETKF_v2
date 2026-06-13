#!/bin/csh
#========================================
# Copy data from JSPACE to local |
#========================================

set JSPACE=${argv[1]}
set dir=${argv[2]}
set letkf=${argv[3]}
set iyr=${argv[4]}
set imon=${argv[5]}
set nmem=${argv[6]}

#=========================================
# Date
#=========================================

if($imon == 2 && $iyr % 4 == 0)then
    @ iday = 29
else if($imon == 2)then
    @ iday = 28
else if($imon == 4 || $imon ==  6 || $imon == 9 || $imon == 11)then
    @ iday = 30
else
    @ iday = 31
endif

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

    if(! -f ${from_file})then
	echo "***Error: Not found ${from_file}"
    	exit 99
    endif
    
    cp -v ${from_file} ${to_dir}
    
    @ imem++

end

#---Ensemble mean
if(! -d ${dir}/${letkf}/output/mean)then
    mkdir -p ${dir}/${letkf}/output/mean
endif

cp -v ${dir}/${letkf}/output/${mmmmm}/restart.${yyyy}${mm}${dd}.nc ${dir}/${letkf}/output/mean/restart.${yyyy}${mm}${dd}.nc
