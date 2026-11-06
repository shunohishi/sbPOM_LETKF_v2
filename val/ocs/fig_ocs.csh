#!/bin/csh
#===============================================#
# Machine
#===============================================#

#set machine="fugaku"
#set machine="jss3"
set machine="rc"

#===============================================#
# Setting
#===============================================#

#---Criteria of observation availability
set tcrit=20
set scrit=20
set ucrit=20
set vcrit=20

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

gmt set FONT=14p,Helvetica,black
gmt set MAP_TITLE_OFFSET=25p
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

gmt set MAP_FRAME_TYPE=PLAIN

gmt begin fig/ocs png

    #---Buoy position
    #Setting
    set size=12d/6d         #Box size
    set range=120/230/20/60 #Domain
    set BAx=a30f10          #X axis interval
    set BAy=a10f5           #Y axis interval
    set BAl=WSne            #Axis label

    #Box
    gmt psbasemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X3 -Y20

    #Coast
    gmt coast -Dl -W0.1,black -Gblack

    #KEO
    echo "144.6 32.3" | gmt psxy -Sa0.5 -Gorange
    echo "146.6 32.3 KEO" | gmt text -F+f14p,0,black+jLM

    #Papa
    echo "215.1 50.1" | gmt psxy -Sa0.5 -Gcyan
    echo "217.1 50.1 Papa" | gmt text -F+f14p,0,black+jLM

    #Title
    echo "120 62 (a) Position of Ocean Climate Stations" |  gmt text -F+f14p,0,black+jLB -N

    #---Observation availability
    set size=4/-6 #Box
    set range=0/100/0/600 #Range
    set BAx=a50f10+l"Obs.\040availability(\045)" #Xaxis
    set BAy=a100f20+l"Depth\040(m)"              #Yaxis
    set label=("(b) T" "(c) S" "(d) U" "(e) V" "(f) T" "(g) S" "(h) U" "(i) V") #Label

    @ i = 1
    @ n = 8
    
    foreach buoy(keo papa)

	if(${buoy} == "keo")then
	    set title="KEO"    #Title
	    set color="orange" #Line color
	else if(${buoy} == "papa")then
	    set title="Papa"
	    set color="cyan"
	endif
    
	foreach var(t s u v)

	    echo ${buoy}:${var}

	    #---Criteria
	    eval set crit = \$${var}crit
	    
	    #---Data
	    set input=dat/${buoy}/${var}rmsd_ave.dat
	    gawk '{print $2,$1 > "dat.20"}' ${input}
	    echo "${crit} 0" > line.20
	    echo "${crit} 600" >> line.20	    
	    
	    #---Figure
	    #Axis
	    if($i % 4 == 1)then
		set BAl=WSne
	    else
		set BAl=wSne
	    endif

	    #Box
	    if($i == 1)then
		gmt psbasemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl}+t"" -Y-8
	    else if($i % 4 == 1)then
		gmt psbasemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X-15 -Y-8
	    else
		gmt psbasemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X5
	    endif

	    #Obs availiability
	    gmt psxy dat.20 -Sc0.1 -G${color}
	    gmt psxy dat.20 -W1,${color}

	    #Criteria
	    gmt psxy line.20 -W1,.

	    #Label
	    echo "0 -10 ${label[$i]}" |  gmt text -F+f14p,0,black+jLB -N

	    #Title
	    if($i % 4 == 1)then
		echo "-50 300 ${title}" |  gmt text -F+f18p,0,black+jLB+a90 -N
	    endif
	    	    
	    @ i++

	    rm -f *.20
	    
	end
    end
    
gmt end

