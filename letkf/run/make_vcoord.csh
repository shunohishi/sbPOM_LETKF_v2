#!/bin/csh
#--------------------------------------------------------------------------
# Make vcoordfile |
#------------------
# 2020.03 Created by S.Ohishi
#--------------------------------------------------------------------------
set WORKDIR=$argv[1];set NMEM=$argv[2];set ENODE=$argv[3]; set EPROC=$argv[4]
#--------------------------------------------------------------------------

#--- echo "Make vcoordfile" #----------------------------------------------
echo ${WORKDIR} > vcoordfile_info.txt
echo ${NMEM} >> vcoordfile_info.txt
echo ${ENODE} >> vcoordfile_info.txt
echo ${EPROC} >> vcoordfile_info.txt
${FC} ens/make_vcoordfile.f90 -o make_vcoordfile.out
if(! -f make_vcoordfile.out)then
    echo "***Error: make_vcoordfile.out not found"
    exit 99
else
    ./make_vcoordfile.out
    rm -f make_vcoordfile.out
endif
#--------------------------------------------------------------------------

#--- echo "Check for vcoordfile" #-----------------------------------------
@ IMEM = 1
while(${IMEM} <= ${NMEM})

    set MMMMM=`printf "%05d" ${IMEM}`
    set vcoordfile=${WORKDIR}/${MMMMM}/vcoordfile
    set nline=`wc -l ${vcoordfile} | awk '{print $1}'`

    if(${EPROC} == ${nline})then
	#echo "vcoordfile FINISHED [${IMEM}/${NMEM}]"
	@ IMEM++
    else
	echo "***Error @ ${MMMMM} vcoordfile"
	exit 99
    endif

end    
#--------------------------------------------------------------------------

wait
exit 0
