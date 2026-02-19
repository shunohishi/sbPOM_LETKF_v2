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

set size=8/-4
set range=-90/90/0/5000
set BAx=a30f10+l"Latitude\040(\260N)"
set BAy=a1000f200+l"Depth\040(m)"
set BAl=WSne
set int1=1/5
set int2=2/10

@ ndat=5 #WOA+Analysis
set label=("(a) WOA" "(b) LORA" "(c) GLORYS2V4" "(d) ORAS5" "(e) C-GLORS")
set datname=("lora" "lora" "glorys" "oras5" "cglors")

#=======================================================
# Color
#=======================================================

#Temperature
gmt makecpt -T0/30/2 -Croma -D -I > t_color.cpt
gmt makecpt -T-2/40/2 -Croma -D -I > t_contour1.cpt
gmt makecpt -T0/40/10 -Croma -D -I > t_contour2.cpt
gmt makecpt -T-1/1/0.25 -Cvik -D > t_dif.cpt

#Salinity
gmt makecpt -T34/36/0.25 -CbatlowK -D > s_color.cpt
gmt makecpt -T20/40/0.25 -Croma -D -I > s_contour1.cpt
gmt makecpt -T20/40/1 -Croma -D -I > s_contour2.cpt
gmt makecpt -T-0.1/0.1/0.025 -Cbam -D -I > s_dif.cpt

#=======================================================
# Figure
#=======================================================

if(! -d fig) mkdir fig

foreach var(t s)

echo "Make figure..."
gmt begin fig/${var}_ysec png

    @ idat = 1
    while($idat <= $ndat)

	echo "Data conversion..."
	set input=dat/${datname[$idat]}_${var}2d.dat
	if($idat == 1)then
	    gawk '{if($4 != -999) print $1,$2,$4 > "obs.20"}' ${input}
	else
	    gawk '{if($3 != -999) print $1,$2,$3 > "dat.20"}' ${input}
	    gawk '{if($3 != -999 && $4 != -999) print $1,$2,$3-$4 > "dif.20"}' ${input}
	endif

	foreach dat(obs dat dif)
	    if(-f ${dat}.20)then
		gmt blockmean ${dat}.20 -R${range} -I${int1} | \
		gmt surface -G${dat}.grd -R${range} -I${int1} -T0.35
	    endif
	end
	
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

	    set color="${var}_color"
	    set drange=8.5/0.5+w3/0.25+e0.5
	    if(${var} == "t") set dBA=a10f2+l"Temperature\040(\260C)"
	    if(${var} == "s") set dBA=a1.0f0.2+l"Salinity"
	    	

	    gmt psmask obs.20 -R${range} -I${int2}
	    gmt grdview obs.grd -C${color}.cpt -Qs
	    gmt grdcontour obs.grd -C${var}_contour1.cpt -W0.2 -A-
	    gmt grdcontour obs.grd -C${var}_contour2.cpt -W1.0 -A-
	    gmt psmask -C

	else

	    set color="${var}_dif"
	    set drange=-4.5/-2+w7/0.25+e0.5+h
	    if(${var} == "t") set dBA=a0.5f0.25+l"Temperature\040bias\040(\260C)"
	    if(${var} == "s") set dBA=a0.05f0.025+l"Salinity\040bias"
	    
	    gmt psmask dif.20 -R${range} -I${int2}
	    gmt grdview dif.grd -C${color}.cpt -Qs
	    gmt psmask -C

	    gmt psmask dat.20 -R${range} -I${int2}
	    gmt grdcontour dat.grd -C${var}_contour1.cpt -W0.2 -A-
	    gmt grdcontour dat.grd -C${var}_contour2.cpt -W1.0 -A-
	    gmt psmask -C	

	endif

	gmt text -F+f14p,0,black+jLB -N <<EOF
-90 -100 ${label[$idat]}
EOF
	        
	if($idat == 1 || $idat == $ndat)then
	    gmt colorbar -Dx${drange} -Bx${dBA} -C${color}.cpt --FONT_ANNOT=20p --FONT_LABEL=20p
	endif

	rm -f *.20 *.grd
		
	@ idat++

    end #idat

    gmt end
    
end #var

rm -f *.cpt
