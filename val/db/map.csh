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

#=======================================================
# Option
#=======================================================

rm -f gmt.conf

gmt set MAP_FRAME_TYPE=PLAIN FORMAT_GEO_MAP=dddF
gmt set FONT=14p,Helvetica,black
gmt set GMT_AUTO_DOWNLOAD=off

#=======================================================
# Figure setting
#=======================================================

set size=6d/6d
set range=118.75/153.75/18.75/48.75 #Describe exaxt edge points
set BAx=a10f5
set BAy=a10f5
set BAl=WSne
set int1=2.5/2.5
set int2=5/5
set drange=6.5/0.5+w5/0.25+ef0.5
set label=("(a) U: RMSD" "U: Obs." "(c) V: RMSD" "(d) V: Obs." "(e) SST: RMSD" "(f) SST: Obs.") 

#=======================================================
# Color
#=======================================================

gmt makecpt -T0/0.4/0.01 -Clajolla -D > urmsd_color.cpt
gmt makecpt -T0/0.4/0.01 -Clajolla -D > vrmsd_color.cpt
gmt makecpt -T0/1.5/0.25 -Clajolla -D > trmsd_color.cpt

gmt makecpt -T0/500/50 -Cbilbao -D > uobs_color.cpt 
gmt makecpt -T0/500/50 -Cbilbao -D > vobs_color.cpt 
gmt makecpt -T0/500/50 -Cbilbao -D > tobs_color.cpt 

#=======================================================
# Figure
#=======================================================

gmt begin map png

    @ i = 1

    foreach var(u v t)
        
	set input=dat/${var}_bin.dat
	gawk '{if($3 != 0) print $1,$2,$3 > "obs.20"}' ${input}
	gawk '{if($5 != -999.) print $1,$2,$5 > "rmsd.20"}' ${input}
	gmt xyz2grd obs.20 -Gobs.grd -R${range} -I${int1}
	gmt xyz2grd rmsd.20 -Grmsd.grd -R${range} -I${int1}

	set input=dat/${var}_ave.dat
	set rmsd=`cat ${input} | gawk '{print $3}'`
	set obs=`cat ${input} | gawk '{print $1}'`
		
	foreach type(rmsd obs)

	    echo ${var}${type}

	    if(${type} == "obs")then
		set dBA=a100f50+l"Num.\040of\040Obs."	
	    else if(${var} == "u" || ${var} == "v")then
		set dBA=a0.1f0.05+l"RMSD\040(m/s)"
	    else if(${var} == "t")then
		set dBA=a0.5f0.25+l"RMSD\040(\260C)"
	    endif
	
	    if($i == 1)then
		gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X3 -Y20
	    else if($i % 2 == 0)then
		gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X10
	    else
		gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X-10 -Y-7.5
	    endif

	    gmt psmask ${type}.20 -I${int1}
	    gmt grdimage ${type}.grd -C${var}${type}_color.cpt
	    gmt psmask -C
	    
	    gmt coast -Dl -W0.2,black -Ggray

	    gmt text -F+f14p,0,black+jLB -N <<EOF
118.75 50 ${label[$i]}
EOF

	    if(${type} == "rmsd")then
		set value="RMSD: ${rmsd}"
	    else if(${type} == "obs")then
		set value="Num: ${obs}"
	    endif

	    gmt text -F+f14p,0,black+jLT -N <<EOF
119 48 ${value}
EOF

	    gmt colorbar -Dx${drange} -Bx${dBA} -C${var}${type}_color.cpt --FONT_ANNOT=14p --FONT_LABEL=14p

	    @ i++

	    rm -f ${type}.20 ${type}.grd

	end	    
    end
    
gmt end

rm -f *.cpt
