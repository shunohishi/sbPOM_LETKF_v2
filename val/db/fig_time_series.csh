#!/bin/csh

if(! -d fig) mkdir fig

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
set end_date=2024-01-01

set size=10/5
set BApx=f1Y
set BAsx=a5Y+l"Year"
set BAl=WSne
set label=("(a) Sea-surface zonal velocity" "(b) Sea-surface meridional velocity" "(c) Sea-surface temperature" "(d) Surface drifter buoy")
set color=("black" "cyan" "orange" "green")
set legend1=("LORA" "GLORYS2V4" "ORAS5" "C-GLORS")
set legend2=("U" "V" "T")

gmt begin fig/time_series png

    @ i = 1
    @ n = 3
    @ ndat = 4
    
    foreach var(u v t)

	echo ${var}
    
	#---DATA
	foreach index(rmsd sprd)
	    set input=dat/${var}${index}_mave.dat
	    gawk -v out=${index}1.20 '{if($6 != -999) print $1,$6 > out}' ${input}
	    gawk -v out=${index}2.20 '{if($7 != -999) print $1,$7 > out}' ${input}
	    gawk -v out=${index}3.20 '{if($8 != -999) print $1,$8 > out}' ${input}
	    gawk -v out=${index}4.20 '{if($9 != -999) print $1,$9 > out}' ${input}
	end
	gawk -v out=${var}nobs.20 '{if($2 != -999) print $1,$2 > out}' ${input}
	    
	if(${var} == "u" || ${var} == "v")then
	    set ye=0.3
	    set BAy=a0.1f0.02+l"RMSD\040(m/s)\040\046\040Spread\040(m/s)"
	else if(${var} == "t")then
	    set ye=1.0
	    set BAy=a0.5f0.1+l"RMSD\040(\260C)\040\046\040Spread\040(\260C)"
	endif
	set range=${start_date}/${end_date}/0/${ye}

	if($i == 1)then
	    gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -X3 -Y23
	else
	    gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -Y-7
	endif

	@ idat=1
	
	foreach dat(rmsd1 rmsd2 rmsd3 rmsd4)
		
	    if($i == 1)then
		if(-f ${dat}.20) gmt psxy ${dat}.20 -W2,${color[$idat]} -l${legend1[$idat]}
	    else
		if(-f ${dat}.20) gmt psxy ${dat}.20 -W2,${color[$idat]}
	    endif
	    
	    @ idat++
	    
	end # dat
	    
	@ idat=1
	
	foreach dat(sprd1 sprd2 sprd3 sprd4)

	    if(-f ${dat}.20)then
		gmt psxy ${dat}.20 -W0.5,${color[$idat]},-
	    endif
	    
	end # dat	

	if($i == 1)then
	    gmt legend -DjRB+jRB+o0.2/0.2 -F+gwhite+pblack --FONT=10p
	endif
	
	echo "2003-01-01 ${ye} ${label[$i]}" | gmt text -F+f14p,0,black+jLT -Dj0.2c/0.2c -N
	
	@ i++

	rm -f rmsd*.20 sprd*.20
	
    end #var

    #---Number of observation
    set ye=50000
    set BAy=a10000f1000+l"Number\040of\040observation"
    set range=${start_date}/${end_date}/0/${ye}

    gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -Y-7

    @ i=1
    foreach var(u v t)
	if(${var} == "u")then
	    gmt psxy ${var}nobs.20 -W4,${color[$i]} -l${legend2[$i]}
	else
	    gmt psxy ${var}nobs.20 -W2,${color[$i]} -l${legend2[$i]}
	endif
	@ i++
    end

    gmt legend -DjRB+jRB+o0.2/0.2 -F+gwhite+pblack --FONT=10p

    echo "2003-01-01 ${ye} ${label[$i]}" | gmt text -F+f14p,0,black+jLT -Dj0.2c/0.2c -N
    
gmt end

rm -f *.20
