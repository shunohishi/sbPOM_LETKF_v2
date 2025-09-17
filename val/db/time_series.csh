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

#=======================================================
# Option
#=======================================================

rm -f gmt.conf

gmt set MAP_FRAME_TYPE=PLAIN FORMAT_GEO_MAP=dddF
gmt set FONT=10p,Helvetica,black
gmt set GMT_AUTO_DOWNLOAD=off

if(${machine} == "jss3")then
    gmt set PS_CONVERT=I+m0.6/0.6/0.6/0.6 #WESN
endif
    
#=======================================================
# Figure setting
#=======================================================

gmt set FONT=14p,Helvetica,black
gmt set FORMAT_DATE_MAP=o FORMAT_TIME_PRIMARY_MAP=c FORMAT_TIME_SECONDARY_MAP=f

set start_date=2003-01-01
set end_date=2019-01-01

set size=12/6
set BApx=f1Y
set BAsx=a5Y+l"Year"
set BAl=WSne
set label=("(a) Sea-surface zonal velocity" "(b) Sea-surface meridional velocity" "(c) Sea-surface temperature")
set color=("black" "cyan" "orange" "green")
set legend=("LORA" "GLORYS2V4" "ORAS5" "C-GLORS")

gmt begin time_series png

    @ i = 1
    @ n = 3
    @ ndat = 4
    
    foreach var(u v t)

	echo ${var}
    
	#---DATA
	set input=dat/${var}_mave.dat
	gawk '{if($10 != -999) print $1,$10 > "rmsd1.20"}' ${input}
	gawk '{if($11 != -999) print $1,$11 > "rmsd2.20"}' ${input}
	gawk '{if($12 != -999) print $1,$12 > "rmsd3.20"}' ${input}
	gawk '{if($13 != -999) print $1,$13 > "rmsd4.20"}' ${input}
	gawk '{if($14 != -999) print $1,$14 > "sprd1.20"}' ${input}
	gawk '{if($15 != -999) print $1,$15 > "sprd2.20"}' ${input}
	gawk '{if($16 != -999) print $1,$16 > "sprd3.20"}' ${input}
	gawk '{if($17 != -999) print $1,$17 > "sprd4.20"}' ${input}
    
	if(${var} == "u" || ${var} == "v")then
	    set ye=0.3
	    set BAy=a0.1f0.02+l"RMSD\040(m/s),\040Spread\040(m/s)"
	else if(${var} == "t")then
	    set ye=1.0
	    set BAy=a0.5f0.1+l"RMSD\040(\260C),\040Spread\040(\260C)"
	endif
	set range=${start_date}/${end_date}/0/${ye}

	if($i == 1)then
	    gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -X3 -Y22
	else
	    gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -Y-8
	endif

	@ idat=1
	
	foreach dat(rmsd1 rmsd2 rmsd3 rmsd4)
		
	    if($i == 1)then
		if(-f ${dat}.20) gmt psxy ${dat}.20 -W1,${color[$idat]} -l${legend[$idat]}
	    else
		if(-f ${dat}.20) gmt psxy ${dat}.20 -W1,${color[$idat]}
	    endif
	    
	    @ idat++
	    
	end # dat
	    
	@ idat=1
	
	foreach dat(sprd1 sprd2 sprd3 sprd4)

	    if(-f ${dat}.20)then
		gmt psxy ${dat}.20 -W1,${color[$idat]},-
	    endif
	    
	end # dat	

	if($i == 1)then
	    gmt legend -DjRB+jRB+o0.2/0.2 -F+gwhite+pblack --FONT=12p
	endif
	
	echo "2003-01-01 ${ye} ${label[$i]}" | gmt text -F+f14p,0,black+jLT -Dj0.2c/0.2c -N
	
	@ i++
	
    end #var

gmt end
