#!/bin/csh
#=== Argument ===============================================================#
set WORKDIR=${argv[1]}; set OUTPUT=${argv[2]}; set INFO=${argv[3]}
set NMEM=${argv[4]}; set EXE=${argv[5]}
set yyyymmdd_n1=${argv[6]}; set yyyymmdd=${argv[7]}
#============================================================================#

#---Forecast/Analysis ensemble mean/sprd
(mv ${WORKDIR}/fcst_mean.nc ${OUTPUT}/mean/fcst.${yyyymmdd}.nc &)
(mv ${WORKDIR}/fcst_sprd.nc ${OUTPUT}/sprd/fcst.${yyyymmdd}.nc &)
(mv ${WORKDIR}/anal_mean.nc ${OUTPUT}/mean/anal.${yyyymmdd}.nc &)
(mv ${WORKDIR}/anal_sprd.nc ${OUTPUT}/sprd/anal.${yyyymmdd}.nc &)

#---NOUT
if(-f ${WORKDIR}/NOUT-00000)then
    (mv ${WORKDIR}/NOUT-00000 ${INFO}/letkf.${yyyymmdd}.log &)
else
    echo "***Error: NOUT-00000 not found"
    exit 99
endif

#---stdout/stderr
foreach filename(stdout stderr)
    if(-s ${WORKDIR}/${filename}.letkf)then
	(mv -f ${WORKDIR}/${filename}.letkf ${INFO}/${filename}.letkf.${yyyymmdd} &)
    else
	(rm -f ${WORKDIR}/${filename}.letkf &)
    endif
end #filename

#---LETKF_yyyymmdd.*.stats/out/err
set JOBID=`head -n 1 ${WORKDIR}/JOBID`

foreach ext(stats out err)
    if(-s ${WORKDIR}/LETKF_${yyyymmdd}.${JOBID}.${ext})then
	(mv -f ${WORKDIR}/LETKF_${yyyymmdd}.${JOBID}.${ext} ${INFO}/ &)
    else
	(rm -f ${WORKDIR}/LETKF_${yyyymmdd}.${JOBID}.${ext} &)
    endif
end #ext

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
