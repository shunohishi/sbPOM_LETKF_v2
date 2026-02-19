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
    gmt set PS_CONVERT=I+m0.4/0.4/0.4/0.8 #WESN
endif
    
#=======================================================
# Figure setting
#=======================================================

@ ndat=5 #CMEMS + Analysis data

set size=8d/4d
set range=0/360/-90/90
set BAx=a60f10
set BAy=a30f10
set BAl=WSne

#For Data conversion
set range1=0/360/-89.75/89.75
set range2=0.25/359.75/-89.75/89.75
set int1=0.25/0.25
set int2=0.25/0.25

set label=("(a) CMEMS" "(b) LORA" "(c) GLORYS2V4" "(d) ORAS5" "(e) C-GLORS")

#=======================================================
# Color
#=======================================================

#MDOT
gmt makecpt -T-2.0/2.0/0.2 -Croma -D -I > color.cpt
gmt makecpt -T-10/10/0.2 -Croma -D -I > contour1.cpt
gmt makecpt -T-10/10/1.0 -Croma -D -I > contour2.cpt

#MDOT difference
gmt makecpt -T-0.20/0.20/0.025 -Cvik -D > dif.cpt

#=======================================================
# Data
#=======================================================

echo "Start: Data conversion..."

echo "Make txt..."

set input=dat/mdot.dat

gawk '{if($3 != -999) print $1,$2,$3 > "dat1.20"}' ${input} &
gawk '{if($4 != -999) print $1,$2,$4 > "dat2.20"}' ${input} &
gawk '{if($5 != -999) print $1,$2,$5 > "dat3.20"}' ${input} &
gawk '{if($6 != -999) print $1,$2,$6 > "dat4.20"}' ${input} &
gawk '{if($7 != -999) print $1,$2,$7 > "dat5.20"}' ${input} &

gawk '{if($3 != -999 || $4 != -999) print $1,$2,$4-$3 > "dif2.20"}' ${input} &
gawk '{if($3 != -999 || $5 != -999) print $1,$2,$5-$3 > "dif3.20"}' ${input} &
gawk '{if($3 != -999 || $6 != -999) print $1,$2,$6-$3 > "dif4.20"}' ${input} &
gawk '{if($3 != -999 || $7 != -999) print $1,$2,$7-$3 > "dif5.20"}' ${input} &
wait

echo "Convert txt => grd..."

@ idat = 1
while($idat <= $ndat)

    if(-f dat${idat}.20) gmt surface dat${idat}.20 -Gdat${idat}.grd -R${range1} -I${int1} -fg -T0.50 &
    if(-f dif${idat}.20) gmt surface dif${idat}.20 -Gdif${idat}.grd -R${range1} -I${int1} -fg -T0.50 &

    @ idat++
    
end
wait

echo "End: Data conversion"

#=======================================================
# Figure
#=======================================================

echo "Start: Make figure..."

if(! -d fig) mkdir fig

gmt begin fig/mdot png

@ idat = 1

while($idat <= $ndat)

    echo "[${idat}/${ndat}]"

    if($idat == 1)then
	gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X3 -Y20
    else if($idat == 2)then
	gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -Y-5.5
    else if($idat % 2 == 1)then
	gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X9.5
    else
	gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X-9.5 -Y-5.5    
    endif

    if($idat == 1)then

	set color="color"
	set drange=8.5/0.5+w3/0.25+e0.5
	set dBA=a1.0f0.2+l"SSH\040(m)"
    
	gmt psmask dat${idat}.20 -R${range2} -I${int1}
	#gmt grdview dat${idat}.grd -C${color}.cpt -Qs
	gmt grdimage dat${idat}.grd -C${color}.cpt
	gmt grdcontour dat${idat}.grd -Ccontour1.cpt -W0.2 -A-
	gmt grdcontour dat${idat}.grd -Ccontour2.cpt -W1.0 -A-
	gmt psmask -C

    else

	set color="dif"
	set drange=-4.5/-1+w7/0.25+e0.5+h
	set dBA=a0.10f0.025+l"SSH\040bias\040(m)"
    
	gmt psmask dif${idat}.20 -R${range2} -I${int1}
	#gmt grdview dif${idat}.grd -C${color}.cpt -Qs
	gmt grdimage dif${idat}.grd -C${color}.cpt
	gmt psmask -C

	gmt psmask dat${idat}.20 -R${range2} -I${int1}
	gmt grdcontour dat${idat}.grd -Ccontour1.cpt -W0.2 -A-
	gmt grdcontour dat${idat}.grd -Ccontour2.cpt -W1.0 -A-
	gmt psmask -C
    
    endif

    gmt coast -R${range} -Dl -W0.2,black -Gwhite

    gmt text -F+f14p,0,black+jLB -N <<EOF
0 95 ${label[$idat]}
EOF
    
    if($idat == 1 || $idat == $ndat)then
	gmt colorbar -Dx${drange} -Bx${dBA} -C${color}.cpt --FONT_ANNOT=20p --FONT_LABEL=20p
    endif
    
    @ idat++
    
end

gmt end

echo "End: Make figure"

rm -f *.20 *.grd *.cpt
