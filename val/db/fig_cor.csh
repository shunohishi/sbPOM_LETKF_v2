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
    gmt set PS_CONVERT=I+m0.4/0.4/0.4/0.8 #WESN
endif
    
#=======================================================
# Figure setting
#=======================================================

set size=8d/4d
set range1=0/360/-90/90
set range2=2.5/357.5/-87.5/87.5 #Describe exact edge points
set BAx=a60f10
set BAy=a30f10
set BAl=WSne
set int1=5/5
set int2=10/10
set label=("(a) Surface zonal velocity" "(b) Surface meridional velocity" "(c) Sea surface temperature")

#=======================================================
# Color
#=======================================================

gmt makecpt -T-0.5/0.5/0.1 -Cvik -D > color.cpt

set drange=0.5/-1+w7/0.25+h+e0.5
set dBA=0.5f0.1+l"Correlation"

#========================================================
# Figure
#========================================================

gmt begin fig/cor png

@ i = 1

foreach var(u v t)

    set input=dat/${var}cor_bin.dat

    #---Data
    gawk '{if($3 != -999. && $3 != 1. && $3 != -1.) print $1,$2,$3 > "dat.20"}' ${input}
    gmt xyz2grd dat.20 -Gdat.grd -R${range2} -I${int1}

    if($i == 1)then
	gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl} -X3 -Y20
    else
	gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl} -Y-5.5
    endif
	
    gmt psmask dat.20 -R${range2} -I${int1}
    gmt grdimage dat.grd -Ccolor.cpt
    gmt psmask -C

    gmt coast -R${range1} -Dl -W0.2,black -Gwhite

    echo "0 95 ${label[$i]}" | gmt text -F+f14p,0,black+jLB -N  
    
    set input=dat/${var}cor.dat
    set ave=`gawk '{printf "%.3f", $1}' ${input}`
    set ave="Correlation: ${ave}"
    echo "360 -90 ${ave}" | gmt text -F+f12p,black+jRB -N
    
    @ i++
    
    rm -f dat.20 dat.grd

end #var

gmt colorbar -Dx${drange} -Bx${dBA} -Ccolor.cpt --FONT=20p

gmt end
    
rm -f color.cpt
