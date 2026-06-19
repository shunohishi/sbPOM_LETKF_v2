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
set iswitch_remove_restart=${argv[8]}

#=========================================
# Date
#=========================================

set yyyy=`printf "%04d" ${iyr}`
set mm=`printf "%02d" ${imon}`
set dd=`printf "%02d" ${iday}`

#=========================================
# Move ensemble mean restart file |
#=========================================

set from_file="${dir}/${letkf}/output/mean/restart.${yyyy}${mm}${dd}.nc"
set to_dir="${JSPACE}/${letkf}/output/mean/"

if(-f ${from_file})then
    mv ${from_file} ${to_dir}
else
    echo "***Error: Not found ${from_file}"
    exit 99
endif

#==========================================
# Remove ensemble restart file |
#==========================================

if(${iswitch_remove_restart} == "on")then

    @ imem = 1
    while($imem <= $nmem)

	set mmmmm=`printf "%05d" ${imem}`

	set file="${dir}/${letkf}/output/${mmmmm}/restart.${yyyy}${mm}${dd}.nc"
	rm -f ${file}
        
	@ imem++

    end

endif
