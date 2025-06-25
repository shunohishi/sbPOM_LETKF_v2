#!/bin/csh
#===============================================#
# Machine
#===============================================#

set machine="fugaku"
#set machine=""

#===============================================#
# Spack load GMT6
#===============================================#

if(${machine} == "fugaku")then
    setenv SPACK_ROOT /vol0004/apps/oss/spack
    source /vol0004/apps/oss/spack/share/spack/setup-env.csh
    spack load /mnrvuuq
endif

#===============================================#
# Option 
#===============================================#

rm -f gmt.conf

gmt set FONT=14p,Helvetica,black
gmt set FORMAT_DATE_MAP=o FORMAT_TIME_PRIMARY_MAP=c FORMAT_TIME_SECONDARY_MAP=f

#================================================#
# Figure setting
#================================================#

set size=12/6
set BApx=a1O
set BAsx=a1Y
set BAl=WSne

set sdate=2021-01-01
set edate=2021-12-31
set label=("(a) Zonal velocity (surface)" "(b) Meridional velocity (surface)" "(c) SST") 

#================================================#
# Make Figure
#================================================#

rm -f time_series.png

gmt begin "time_series" png

    @ i = 1

    foreach var(u v t)

	echo ${var}

	#---Label
	if(${var} == "u" || ${var} == "v")then
	    set sx=0
	    set ex=0.4
	    set dxa=a0.1
	    set dxf=f0.05
	    set ylabel="RMSD\040(m/s)"
	else if(${var} == "t")then
	    set sx=0
	    set ex=1.5
	    set dxa=a0.5
	    set dxf=f0.1
	    set ylabel="RMSD\040(\237C)"
	endif
	
	set range=${sdate}/${edate}/${sx}/${ex}
	set BAy="${dxa}${dxf}+l${ylabel}"

	#---Data
	set input=dat/${var}_mave.dat
	gawk '{if($4 != -999.) print $1,$4 > "dat.20"}' ${input}

	#---Figure
	if($i == 1)then
	    gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -X3 -Y23
	else
	    gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -Y-14.5
	endif

	gmt psxy dat.20 -W2,black
	gmt psxy dat.20 -Sc0.15 -Gblack

	gmt text -F+f14p,0,black+jLB -N -Y6.2 <<EOF
${sdate} 0 ${label[$i]}
EOF
		
	rm -f dat.20

	@ i++
	
    end
	
gmt end
