#!/bin/csh
#===============================================#
# Machine
#===============================================#

#set machine="fugaku"
set machine="jss3"

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
gmt set GMT_AUTO_DOWNLOAD=off

if(${machine} == "jss3")then
    gmt set PS_CONVERT=I+m0.4/0.4/0.4/0.4 #WESN
endif

#================================================#
# Figure setting
#================================================#

set size=12/6
set BApx=f1Y
set BAsx=a5Y
set BAl=WSne

set sdate=2003-01-01
set edate=2018-12-31
set label=("(a) Nino 3.4" "(b) IOD" "(c) Atlantic Nino") 

#================================================#
# Make Figure
#================================================#

rm -f time_series.png

gmt begin "time_series" png

    @ i = 1

    foreach index(nino34 iod anino)

	echo ${index}

	#---Label
	if(${index} == "nino34")then
	    set sx=-3.0
	    set ex=3.0
	    set dxa=a1.0
	    set dxf=f0.1
	else
	    set sx=-1.5
	    set ex=1.5
	    set dxa=a0.5
	    set dxf=f0.1
	endif
	set ylabel="Index\040(\260C)"	
	set range=${sdate}/${edate}/${sx}/${ex}
	set BAy="${dxa}${dxf}g99+l${ylabel}"

	#---Data
	set input=dat/${index}.dat
	foreach data(lora glorys oras5 cglors)
	    if(${data} == "lora")   gawk -v out=${data}.20 '{print $1,$2 > out}' ${input}
	    if(${data} == "glorys") gawk -v out=${data}.20 '{print $1,$3 > out}' ${input}
	    if(${data} == "oras5")  gawk -v out=${data}.20 '{print $1,$4 > out}' ${input}
	    if(${data} == "cglors") gawk -v out=${data}.20 '{print $1,$5 > out}' ${input}
	    if(${data} == "lora")   gawk -v out=${data}-event.20 '{if($6 <= -1 || 1 <= $6) print $1,$2 > out}' ${input}
	    if(${data} == "glorys") gawk -v out=${data}-event.20 '{if($7 <= -1 || 1 <= $7) print $1,$3 > out}' ${input}
	    if(${data} == "oras5")  gawk -v out=${data}-event.20 '{if($8 <= -1 || 1 <= $8) print $1,$4 > out}' ${input}
	    if(${data} == "cglors") gawk -v out=${data}-event.20 '{if($9 <= -1 || 1 <= $9) print $1,$5 > out}' ${input}

	end
	
	#---Figure
	if($i == 1)then
	    gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -X3 -Y23
	else
	    gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -Y-14.5
	endif

	foreach data(lora glorys oras5 cglors)
	
	    if(${data} == "lora")then
		set color="black"
		set legend="LORA"
	    else if(${data} == "glorys")then
		set color="green"
		set legend="GLORYS"
	    else if(${data} == "oras5")then
		set color="cyan"
		set legend="ORAS5"
	    else if(${data} == "cglors")then
		set color="orange"
		set legend="CGLORS"
	    endif
		
	    if(${data} == "lora")then
		set thick=4.0
	    else
		set thick=1.0
	    endif

	    if(${index} == "nino34")then
		gmt psxy ${data}.20 -W${thick},${color} -l${legend}
	    else
		gmt psxy ${data}.20 -W${thick},${color}
	    endif
	    
	end

	if(${index} == "nino34")then
	    gmt legend -DjLT+jLT+o0.2/0.2 -F+gwhite+pblack --FONT_ANNOT_PRIMARY=10p
	endif
	
	gmt text -F+f14p,0,black+jLB -Wwhite -N -Y6.2 <<EOF
${sdate} ${sx} ${label[$i]}
EOF

	rm -f *.20
	

	@ i++
	
    end
	
gmt end
