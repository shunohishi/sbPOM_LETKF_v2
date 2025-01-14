#!/bin/csh
#---------------------------------------------------------#
# ARGV |
#---------------------------------------------------------#
set REGION=${argv[1]}
set WORKDIR=${argv[2]}
set BUDGET=${argv[3]}
set TS_NUDGE=${argv[4]}
set TI_NUDGE=${argv[5]}
set SS_NUDGE=${argv[6]}
set SI_NUDGE=${argv[7]}

set syyyy=${argv[8]}
set smm=${argv[9]}
set switch_rst=${argv[10]}
set nday=${argv[11]}

#---------------------------------------------------------#

echo "Make namelist"
cd ${WORKDIR}    
cat <<EOF > pom.nml
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
  iens=0
  nens=1
  time_start = "${syyyy}-${smm}-01 00:00:00 +00:00"
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
