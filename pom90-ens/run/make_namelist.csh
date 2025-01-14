#!/bin/csh
#====================================================
# ARGV |
#====================================================

set REGION=${argv[1]}
set DIR=${argv[2]}
set NMEM=${argv[3]}
set BUDGET=${argv[4]}
set TS_NUDGE=${argv[5]}
set TI_NUDGE=${argv[6]}
set SS_NUDGE=${argv[7]}
set SI_NUDGE=${argv[8]}

set WORKDIR=${argv[9]}
set yyyy_s=${argv[10]}
set mm_s=${argv[11]}
set nday=${argv[12]}
set switch_rst=${argv[13]}

#====================================================
# Make Namelist |
#====================================================

@ IMEM = 1
while(${IMEM} <= ${NMEM})

    set MMMMM=`printf "%05d" ${IMEM}`

    cat <<EOF > ${WORKDIR}/${MMMMM}/pom.nml
&pom_nml
  title = "${REGION}"
  netcdf_file = "${REGION}"
  mode = 3
  assim = 0
  nadv = 2
  nitera = 2
  sw = 1.0
  npg = 2
  dte = 10.0
  isplit = 30
  iens = ${IMEM}
  nens = ${NMEM}
  time_start = "${yyyy_s}-${mm_s}-01 00:00:00 +00:00"
  nread_rst = ${switch_rst}
  read_rst_file = "restart.nc"
  read_iau_file = "iau.nc"
  write_rst = ${nday}
  write_rst_file = "restart"
  write_iau_file = "iau"
  budget = ${BUDGET}
  days = ${nday}
  prtd1 = 1.0
  prtd2 = 1.0
  swtch = 9999.
  ts_nudge = ${TS_NUDGE}
  ti_nudge = ${TI_NUDGE}
  ss_nudge = ${SS_NUDGE}
  si_nudge = ${SI_NUDGE}
/
EOF

    if(! -f ${WORKDIR}/${MMMMM}/pom.nml)then
	echo "***Error: Not found ${WORKDIR}/${MMMMM}/pom.nml"
	exit
    endif

    @ IMEM++

end #IMEM
