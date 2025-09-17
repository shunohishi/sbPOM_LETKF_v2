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
gmt set MAP_FRAME_TYPE=PLAIN FORMAT_GEO_MAP=dddF
gmt set GMT_AUTO_DOWNLOAD=off

if(${machine} == "jss3")then
    gmt set PS_CONVERT=I+m0.5/0.5/0.5/1.0 #WESN
endif

#================================================#
# Common figure setting
#================================================#

set int1=0.25/0.25
set int2=0.50/0.50
set drange=6.0/-1.5+w6/0.25+e0.5+h

set label=("(a) LORA" "(b) LORA" "(c) GLORYS" "(d) GLORYS" \
"(e) ORAS5" "(f) ORAS5" "(g) CGLORS" "(h) CGLORS")

gmt makecpt -T-2/2/0.25 -Cvik -D > enso.cpt
gmt makecpt -T-1.5/1.5/0.25 -Cvik -D > iod.cpt
gmt makecpt -T-1.5/1.5/0.25 -Cvik -D > anino.cpt
gmt makecpt -T-5/40/1 -Cno_green > contour1.cpt
gmt makecpt -T0/30/5 -Cno_green > contour2.cpt

#================================================#
# Make Figure
#================================================#

foreach var(enso iod anino)

    gmt begin "${var}" png

    #---Domain
    if(${var} == "enso")then
	set size=8d/4d
	set sx=120; set ex=300
	set sy=-30; set ey=30
	set range=${sx}/${ex}/${sy}/${ey}
	set BAx=a60f5; set BAy=a15f5
	set dx=10; set dy=5
	set yyyy_p=2015; set mm_p=11
	set yyyy_n=2010; set mm_n=12
	set event_p="El Nino at ${yyyy_p}.${mm_p}"
	set event_n="La Nina at ${yyyy_n}.${mm_n}"
	set dBA=a1f0.25+l"SST\040anomaly\040(\260C)"
    else if(${var} == "iod")then
	set size=8d/4d
	set sx=30; set ex=120
	set sy=-30; set ey=30
	set range=${sx}/${ex}/${sy}/${ey}
	set BAx=a60f5; set BAy=a15f5
	set dx=10; set dy=5
	set yyyy_p=2006; set mm_p=10
	set yyyy_n=2010; set mm_n=09
	set event_p="Positive IOD at ${yyyy_p}.${mm_p}"
	set event_n="Negative IOD at ${yyyy_n}.${mm_n}"
	set dBA=a0.75f0.25+l"SST\040anomaly\040(\260C)"
    else if(${var} == "anino")then
    	set size=8d/4d
	set sx=300; set ex=380
	set sy=-30; set ey=30
	set range=${sx}/${ex}/${sy}/${ey}
	set BAx=a20f5; set BAy=a15f5
	set dx=10; set dy=5
	set yyyy_p=2009; set mm_p=05
	set yyyy_n=2012; set mm_n=01
	set event_p="Atlantic Nino at ${yyyy_p}.${mm_p}"
	set event_n="Atlantic Nina at ${yyyy_n}.${mm_n}"
	set dBA=a0.75f0.25+l"SST\040anomaly\040(\260C)"
    endif
    
    @ i = 1
    @ n = 8

    foreach data(lora glorys oras5 cglors)

	echo ${var}:${data}

	#---Data
	if(${data} == "lora")   set dd="01"
	if(${data} == "glorys") set dd="02"
	if(${data} == "oras5")  set dd="03"
	if(${data} == "cglors") set dd="04"

	#Positive event
	set input=dat/mmean${yyyy_p}${mm_p}-d${dd}.dat
	gawk '{if($4 != -999) print $1,$2,$4 > "sst_p.20"}' ${input}
	gawk '{if($5 != -999) print $1,$2,$5 > "ssta_p.20"}' ${input}
	gmt surface sst_p.20 -Gsst_p.grd -R${range} -I${int1} -fg -T0.8
	gmt surface ssta_p.20 -Gssta_p.grd -R${range} -I${int2} -fg -T0.8

	#Negative event
	set input=dat/mmean${yyyy_n}${mm_n}-d${dd}.dat
	gawk '{if($4 != -999) print $1,$2,$4 > "sst_n.20"}' ${input}
	gawk '{if($5 != -999) print $1,$2,$5 > "ssta_n.20"}' ${input}
	gmt surface sst_n.20 -Gsst_n.grd -R${range} -I${int1} -fg -T0.9
	gmt surface ssta_n.20 -Gssta_n.grd -R${range} -I${int2} -fg -T0.9
	
	#---Figure
	foreach event(p n)

	    #Axis
	    if($i == $n - 1)then
		set BAl=WSne
	    else if($i == $n)then
		set BAl=wSne
	    else if($i % 2 == 1)then
		set BAl=Wsne
	    else
		set BAl=wsne
	    endif    
	
	    if(${event} == "p" && $i == 1)then
		gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl}+t"${event_p}" -X3 -Y23
	    else if(${event} == "n" && $i == 2)then
		gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl}+t"${event_n}" -X${dx}
	    else if(${event} == "p")then
		gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X-${dx} -Y-${dy}
	    else if(${event} == "n")then
		gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X${dx}	    
	    endif
	
	    gmt psmask ssta_${event}.20 -I${int2}
	    gmt grdview ssta_${event}.grd -C${var}.cpt -Qs
	    gmt psmask -C

	    gmt psmask sst_${event}.20 -I${int2}
	    gmt grdcontour sst_${event}.grd -Ccontour1.cpt -W0.1 -A-
	    gmt grdcontour sst_${event}.grd -Ccontour2.cpt -W1.0 -A+a0+gwhite
	    gmt psmask -C
    
	    gmt coast -Dl -W0.1,black -Ggray
	    
	    echo "${sx} ${ey} ${label[$i]}" | gmt text -F+f14p,0,black+jLT -Gwhite	    

	    if(${var} == "iod")then
		gmt psxy ${var}_west_domain.dat -W1,white -A
		gmt psxy ${var}_east_domain.dat -W1,white -A
	    else
		gmt psxy ${var}_domain.dat -W1,white -A
	    endif
		
	    if($i == $n - 1)then
		gmt colorbar -Dx${drange} -Bx${dBA} -C${var}.cpt --FONT_ANNOT=20p --FONT_LABEL=20p	
	    endif

	    @ i++

	end #event

	rm -f *.20 *.grd
	
    end #data
	
    gmt end

end #$var

rm -f *.cpt
