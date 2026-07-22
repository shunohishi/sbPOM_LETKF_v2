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

set size=6/-6
set BAl=WSne
set label=("(a) Temperature" "(b) Salinity" "(c) Zonal velocity" "(d) Meridional velocity")
set buoyname=("KEO" "Papa")

@ ibuoy = 1
@ nbuoy = 2
foreach buoy(keo papa)

    set title="${buoyname[$ibuoy]} buoy"
    echo ${title}

    gmt begin fig/${buoy}-cor png

    @ ivar = 1
    @ nvar = 4
    foreach var(t s u v)
    
	#---NOBS
	if(${var} == "t" || ${var} == "s")then
	    @ ratio = 20 #[%]
	else if(${var} == "u" || ${var} == "v")then
	    @ ratio = 20 #[%]
	endif	
               
	if(${buoy} == "keo" && (${var} == "t" || ${var} == "s"))then
	    set xs=-1; set xe=1
	    set ys=0; set ye=550
	    set BAx=a0.5f0.1g2+l"Correlation"
	    set BAy=a100f10+l"Depth\040(m)"
	else if(${buoy} == "keo" && (${var} == "u" || ${var} == "v"))then
	    set xs=-1; set xe=1
	    set ys=0; set ye=40
	    set BAx=a0.5f0.1g2+l"Correlation"
	    set BAy=a10f5+l"Depth\040(m)"
	else if(${buoy} == "papa" && (${var} == "t" || ${var} == "s"))then
	    set xs=-1; set xe=1
	    set ys=0; set ye=350
	    set BAx=a0.5f0.1g2+l"Correlation"
	    set BAy=a100f10+l"Depth\040(m)"
	else if(${buoy} == "papa" && (${var} == "u" || ${var} == "v"))then
	    set xs=-1; set xe=1
	    set ys=0; set ye=40
	    set BAx=a0.5f0.1g2+l"Correlation"
	    set BAy=a10f5+l"Depth\040(m)"
	endif
	set range=${xs}/${xe}/${ys}/${ye}


	set input=dat/${buoy}/${var}cor_ave.dat
	gawk -v ratio=${ratio} '{if(ratio < $2 && $6 != -999.) print $6,$1 > "dat.20"}' ${input}    

	set input=dat/${buoy}/${var}cor_dave.dat
	gawk -v ye=${ye} '{print $1, ye > "ave.20"}' ${input}
	
	if(${ivar} == 1)then
	    gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl}+t"${title}" -X2.5 -Y20
	else if(${ivar} % 2 == 0)then
	    gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X8.5
	else
	    gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X-8.5 -Y-8.5
	endif

	gmt psxy dat.20 -W2
	gmt psxy dat.20 -Sc0.2 -Gblack
	gmt psxy ave.20 -Si0.4 -Gblack -N

	echo "${xs} 0 ${label[$ivar]}" | gmt text -F+f14p,0,black+jLB -Dj0.c/0.2c -N 
	
	rm -f dat.20 ave.20
	
	@ ivar++
	
    end

    gmt end

    @ ibuoy++
    
end
	
