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

gmt set FONT=10p,Helvetica,black
gmt set MAP_FRAME_TYPE=PLAIN FORMAT_GEO_MAP=dddF
gmt set GMT_AUTO_DOWNLOAD=off

if(${machine} == "jss3")then
    gmt set PS_CONVERT=I+m0.6/0.6/0.6/0.6 #WESN
endif

#=======================================================
# Make directory
#=======================================================

if(! -d fig) mkdir fig

#=======================================================
# Figure setting
#=======================================================

set size=8d/4d
set range=0/360/-90/90
set BAx=a60f10
set BAy=a30f10
set BAl=WSne

#=======================================================
# Color
#=======================================================

gmt makecpt -T-0.5/0.5/0.1 -Cvik -D > color.cpt
set drange=0.5/-1+w7/0.25+h+e0.5
set dBA=a0.5f0.1+l"Correlation"

#=======================================================
# DATA
#=======================================================

@ nobs=365

#---RMSD
set input=dat/cor_ave.dat
gawk -v nobs=${nobs} '{if(nobs < $4 && $8 != -999 && $8 <= 0.) print $2,$3,$8 > "n.20"}' ${input}
gawk -v nobs=${nobs} '{if(nobs < $4 && $8 != -999 && 0. < $8) print $2,$3,$8 > "p.20"}' ${input}

#=======================================================
# Figure
#=======================================================

gmt begin fig/tide_cor png

gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X3 -Y20

gmt coast -Dl -W0.2,black -Ggray
gmt psxy n.20 -Sc0.2 -W0.2 -Ccolor.cpt
gmt psxy p.20 -Sc0.2 -W0.2 -Ccolor.cpt

#Ave
set input=dat/cor_ave_all.dat
set ave=`gawk '{printf "%.4f", $5}' ${input}`
set ave="Correlation (station-mean): ${ave}"
echo "360 -90 ${ave}" | gmt text -F+f12p,black+jRB -N 
    
gmt colorbar -Dx${drange} -Bx${dBA} -Ccolor.cpt --FONT=20p

gmt end

rm -f *.20 *.cpt


