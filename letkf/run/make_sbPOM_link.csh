#!/bin/csh
#=== ARGUMENT =================================#
set switch_iau=${argv[1]}
set NMEM=${argv[2]}
set WORKDIR=${argv[3]}
set OUTPUT=${argv[4]}
set MODELDIR=${argv[5]}
set MODELDATADIR=${argv[6]}
set ATMDIR=${argv[7]}
set RIVDIR=${argv[8]}
set REGION=${argv[9]}
set EXE=${argv[10]}
set BUDGET=${argv[11]}
set TS_NUDGE=${argv[12]}
set TI_NUDGE=${argv[13]}
set SS_NUDGE=${argv[14]}
set SI_NUDGE=${argv[15]}
set ST=${argv[16]}
set IT=${argv[17]}
set DT=${argv[18]}

#=== Date =================================#
set sdate=`perl juldays.prl ${ST}`
set yyyys=`printf "%04d" ${sdate[1]}`; set mms=`printf "%02d" ${sdate[2]}`; set dds=`printf "%02d" ${sdate[3]}`

@ BT = ${IT} - 1
set date_n1=`perl juldays.prl ${BT}`
set yyyy_n1=`printf "%04d" ${date_n1[1]}`; set mm_n1=`printf "%02d" ${date_n1[2]}`; set dd_n1=`printf "%02d" ${date_n1[3]}`
set yyyymmdd_n1=${yyyy_n1}${mm_n1}${dd_n1}

set date=`perl juldays.prl ${IT}`
set yyyy=`printf "%04d" ${date[1]}`; set mm=`printf "%02d" ${date[2]}`; set dd=`printf "%02d" ${date[3]}`
set yyyymmdd=${yyyy}${mm}${dd}

#Output timing for restart file
if(${switch_iau} == 0 || ${switch_iau} == 2)then #Intermittent/2nd simulation in IAU
    @ RT = ${DT}
else if(${switch_iau} == 1)then #1st simulation in IAU
    @ RT = 999
endif

#==========================================#

set GRID=grid               #Model Grid data
set TSCLIM=tsclim           #T/S Climatology from WOA
set IC=ic.woa18.${mms}      #Initial Condition Netcdf file
set LBC=lbc_mclim           #Lateral Boundary Condition Netcdf file
set TSDATA=tsdata_mclim     # Nudging Netcdf file

@ IMEM = 1
while(${IMEM} <= ${NMEM})

    set MMMMM=`printf "%05d" ${IMEM}`

    #---echo "Make link"
    #EXE
    (ln -s ${MODELDIR}/${EXE} ${WORKDIR}/${MMMMM}/${EXE} &)
    #GRID
    (ln -s ${MODELDATADIR}/${GRID}.nc ${WORKDIR}/${MMMMM}/in/${REGION}.grid.nc &)
    #TSCLIM
    (ln -s ${MODELDATADIR}/${TSCLIM}.nc ${WORKDIR}/${MMMMM}/in/${REGION}.tsclim.nc &)
    #IC
    (ln -s ${MODELDATADIR}/${IC}.nc ${WORKDIR}/${MMMMM}/in/${REGION}.ic.nc &)
    #TSDATA
    (ln -s ${MODELDATADIR}/${TSDATA}.${MMMMM}.nc ${WORKDIR}/${MMMMM}/in/${REGION}.tsdata.nc &)
    #LBC
    (ln -s ${MODELDATADIR}/${LBC}.${MMMMM}.nc ${WORKDIR}/${MMMMM}/in/${REGION}.lbc.nc &)
    #ATM
    (ln -s ${ATMDIR} ${WORKDIR}/${MMMMM}/in/atm &)
    (ln -s ${ATMDIR}/atm_clim.nc ${WORKDIR}/${MMMMM}/in/${REGION}.atm_clim.nc &)
    #RIVER
    (ln -s ${RIVDIR} ${WORKDIR}/${MMMMM}/in/river &)
    
    #Restart
    (ln -s ${OUTPUT}/${MMMMM}/restart.${yyyymmdd_n1}.nc ${WORKDIR}/${MMMMM}/in/restart.nc &)

    #IAU
    if(${switch_iau} == 2)then
	(ln -s ${OUTPUT}/${MMMMM}/iau.${yyyymmdd}.nc ${WORKDIR}/${MMMMM}/in/iau.nc &)
    endif
	
    #---echo "Make namelist"    
    cd ${WORKDIR}/${MMMMM}/

    if(${switch_iau} == 0) set assim=3 #Intermittent method (Initial condition: Analysis)
    if(${switch_iau} == 1) set assim=1 #IAU1 (Free simulation but save only IAU file)
    if(${switch_iau} == 2) set assim=2 #IAU2 (Adding increment)

    cat <<EOF > pom.nml
&pom_nml
    title = "${REGION}"
    netcdf_file = "${REGION}"
    mode = 3
    assim = ${assim}
    nadv = 2
    nitera = 2
    sw = 1.0
    npg = 2
    dte = 10.0
    isplit = 30
    iens = ${IMEM}
    nens = ${NMEM}
    time_start = "${yyyys}-${mms}-${dds} 00:00:00 +00:00"
    nread_rst = 1
    read_rst_file = "restart.nc"
    read_iau_file = "iau.nc"
    write_rst = ${RT}.
    write_rst_file = "restart"
    write_iau_file = "iau"
    budget = ${BUDGET}
    days =  ${DT}.
    prtd1 = 1.0
    prtd2 = 1.0
    swtch = 9999.
    ts_nudge = ${TS_NUDGE}
    ti_nudge = ${TI_NUDGE}
    ss_nudge = ${SS_NUDGE}
    si_nudge = ${SI_NUDGE}
/
EOF

    @ IMEM++

end #IMEM

wait
exit 0
