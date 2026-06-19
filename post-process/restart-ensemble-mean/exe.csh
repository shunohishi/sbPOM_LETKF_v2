#!/bin/csh
#============================================================
# Make mean restart file (Post process) |
#============================================================
# Created by S.Ohishi @ 2026.06.15
#
#============================================================

set sdate=(2002 5)
set edate=(2024 1)
set JSPACE="/hpss_chofu/home/D/DR29900/LORA/QGLOBAL"
set dir="/data/R/R2402/ohishi/QGLOBAL/post-process/restart-ensemble-mean"
set letkf="letkf"
set nmem=128

#JSPACE
set iswitch_copy_jspace="on"    #Copy restart file from JSPACE

#REMOVE restart
#set iswitch_remove_restart="off" #Remove restart file
set iswitch_remove_restart="on" #Remove restart file

#MACHINE
set machine="jss3"
#set machine="fugaku"

#RESOURCE UNIT
set RSCUNIT=SORA         #JSS3
#set RSCUNIT=RURI         #JSS3
#set RSCUNIT=rscunit_ft01 #Fugaku

#============================================================

if(! -d info) mkdir info

#============================================================

@ iyr  = $sdate[1]
@ imon = $sdate[2]

while ( $iyr <= $edate[1] )

    if ( $iyr == $edate[1] ) then
        @ emon = $edate[2]
    else
        @ emon = 12
    endif

    while ( $imon <= $emon )

	#---Set iday
	if($imon == 2 && $iyr % 4 == 0)then
	    @ iday = 29
	else if($imon == 2)then
	    @ iday = 28
	else if($imon == 4 || $imon ==  6 || $imon == 9 || $imon == 11)then
	    @ iday = 30
	else
	    @ iday = 31
	endif
	
	echo "Date: ${iyr} ${imon} ${iday}"
	
	#---Copy ensemble restart file from JSPACE
	if(${iswitch_copy_jspace} == "on")then
	    csh copy_ensemble_restart.csh ${JSPACE} ${dir} ${letkf} ${iyr} ${imon} ${iday} ${nmem}
	    if($? == 99)then
		exit
	    endif	    
	endif
	    
        #---Execute
	csh submit_job.csh ${dir} ${letkf} ${iyr} ${imon} ${iday} ${nmem} ${machine} ${RSCUNIT}

	#---Move ensemble-mean restart file to JSPACE
	if(${iswitch_copy_jspace} == "on")then
	    csh move_ensemble_mean_restart.csh ${JSPACE} ${dir} ${letkf} ${iyr} ${imon} ${iday} ${nmem} ${iswitch_remove_restart}
	endif
	
        @ imon = $imon + 1

    end #imon

    @ iyr = $iyr + 1
    @ imon = 1

end #iyr
