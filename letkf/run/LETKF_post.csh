#!/bin/csh
#=== Argument ===============================================================#
set WORKDIR=${argv[1]}
set OUTPUT=${argv[2]}
set INFO=${argv[3]}
set NMEM=${argv[4]}
set machine=${argv[5]}

set EXE=${argv[6]}
set yyyymmdd_n1=${argv[7]}
set yyyymmdd=${argv[8]}
#============================================================================#

#---Forecast/Analysis ensemble mean/sprd
(mv ${WORKDIR}/fcst_mean.nc ${OUTPUT}/mean/fcst.${yyyymmdd}.nc &)
(mv ${WORKDIR}/fcst_sprd.nc ${OUTPUT}/sprd/fcst.${yyyymmdd}.nc &)
(mv ${WORKDIR}/anal_mean.nc ${OUTPUT}/mean/anal.${yyyymmdd}.nc &)
(mv ${WORKDIR}/anal_sprd.nc ${OUTPUT}/sprd/anal.${yyyymmdd}.nc &)

#---Innovation statistics
(mv ${WORKDIR}/inv.nc ${OUTPUT}/mean/inv.${yyyymmdd}.nc &)

#---NOUT
if(-f ${WORKDIR}/NOUT-00000)then
    (mv ${WORKDIR}/NOUT-00000 ${INFO}/letkf.${yyyymmdd}.log &)
else
    echo "***Error: NOUT-00000 not found"
    exit 99
endif

#---stdout
if(${machine} == "jss3" && -s ${WORKDIR}/stdout.letkf)then
    (mv -f ${WORKDIR}/stdout.letkf ${INFO}/stdout.letkf.${yyyymmdd} &)
else if(${machine} == "jss3")then
    (rm -f ${WORKDIR}/stdout.letkf &)
else if(${machine} == "fugaku")then
    foreach file(${WORKDIR}/stdout.letkf.*.0)
	if(-s ${file})then
	    (mv -f ${file} ${INFO}/stdout.letkf.${yyyymmdd} &)
	    break
	endif
    end
endif

#---stderr
if(${machine} == "jss3" && -s ${WORKDIR}/stderr.letkf)then
    (mv -f ${WORKDIR}/stderr.letkf ${INFO}/stderr.letkf.${yyyymmdd} &)
else if(${machine} == "jss3")then
    (rm -f ${WORKDIR}/stderr.letkf &)
else if(${machine} == "fugaku")then

#    foreach file(${WORKDIR}/stderr.letkf.*.0)
#	if(-s ${file})then
#	    (mv -f ${file} ${INFO}/stderr.letkf.${yyyymmdd} &)
#	    break
#	endif
#    end

endif

#---LETKF_yyyymmdd.*.stats/out/err (*Removed)
set JOBID=`head -n 1 ${WORKDIR}/JOBID`

#foreach ext(stats out err)
#    if(-s ${WORKDIR}/LETKF_${yyyymmdd}.${JOBID}.${ext})then
#	(mv -f ${WORKDIR}/LETKF_${yyyymmdd}.${JOBID}.${ext} ${INFO}/ &)
#    else
#	(rm -f ${WORKDIR}/LETKF_${yyyymmdd}.${JOBID}.${ext} &)
#    endif
#end #ext

#---fa01MMMMM.nc
@ IMEM = 1
while(${IMEM} <= ${NMEM})
    set MMMMM=`printf "%05d" ${IMEM}`
    (unlink ${WORKDIR}/fa01${MMMMM}.nc &)
    @ IMEM++
end

#---letkf.exe obs01.nc grid.nc
foreach filename(${EXE} obs01.nc grid.nc)
    (unlink ${WORKDIR}/${filename} &)
end

#---Others
(rm -f ${WORKDIR}/NOUT-* ${WORKDIR}/JOBID ${WORKDIR}/FINISHED_LETKF &)

wait
exit 0
