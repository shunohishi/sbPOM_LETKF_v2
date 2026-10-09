#!/bin/csh

set nyr=18 #Year 
set obsa=1 #Obs. availability limit (Unit: month^-1)
set ndat=5 #Number of dataset

#===============================================#
# Machine
#===============================================#

#set machine="fugaku" #Fugaku
#set machine="jss3"   #JSS3
set machine="rc"      #R-CCS Cloud

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

set xsize=8; set ysize=4   #Figure size                             
set slon=110; set elon=250 #Longitude
set slat=15;  set elat=60  #Latitude
set res=5                  #Bin resolution
set BAx=a30f10; set BAy=a20f10; set BAl=WSne #Axis option
set label=("(a) Surface zonal velocity" "(b) Surface meridional velocity" "(c) Sea surface temperature")
set plat=62 #Latitude (Label plot)

#=======================================================
# Color setting
#=======================================================

gmt makecpt -T-0.5/0.5/0.1 -Cvik -D > color.cpt

set dBA=0.5f0.1+l"Correlation\040(RMSD\040vs.\040Spread)"

#========================================================
# Figure
#========================================================

if(! -d fig) mkdir fig

set size=${xsize}d/${ysize}d
set range1=${slon}/${elon}/${slat}/${elat}
set slon_bin = `gawk -v slon=$slon -v res=$res 'BEGIN{print slon + res/2}'`; set elon_bin = `gawk -v elon=$elon -v res=$res 'BEGIN{print elon - res/2}'`
set slat_bin = `gawk -v slat=$slat -v res=$res 'BEGIN{print slat + res/2}'`; set elat_bin = `gawk -v elat=$elat -v res=$res 'BEGIN{print elat - res/2}'`
set range2=${slon_bin}/${elon_bin}/${slat_bin}/${elat_bin}
set int=${res}/${res}

gmt begin fig/cor png

@ ifig = 1

foreach var(u v t)

    set input=dat/${var}cor_bin.dat

    #---Data
    @ obs_col = 2 + 1
    @ dat_col = 2 + ${ndat} + 1
    gawk -v obs_col=${obs_col} -v dat_col=${dat_col} -v obsa=${obsa} -v nyr=${nyr} \
    '{if(obsa <= $obs_col/(nyr*12) &&  $dat_col != -999. && $dat_col*$dat_col != 1.) print $1,$2,$dat_col > "dat.20"}' ${input}
    gmt xyz2grd dat.20 -Gdat.grd -R${range2} -I${int}


    #---Figure
    #Box
    set ix = `gawk -v xsize=$xsize 'BEGIN{print (xsize + 2.5)}'`
    set iy = `gawk -v ysize=$ysize 'BEGIN{print (ysize + 1.5)}'`
    if(${ifig} == 1)then
	gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl} -X3 -Y20
    else if(${ifig} % 2 == 0)then
	gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl} -X${ix}
    else
	gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl} -X-${ix} -Y-${iy}
    endif

    #Plot data
    gmt psmask dat.20 -R${range2} -I${int}
    gmt grdimage dat.grd -Ccolor.cpt
    gmt psmask -C

    #Coast
    gmt coast -R${range1} -Dl -W0.2,black -Gwhite

    #Label
    echo "${slon} ${plat} ${label[$ifig]}" | gmt text -F+f14p,0,black+jLB -N  
        
    @ ifig++
    
    rm -f dat.20 dat.grd

end #var

#Color scale
set cxsize = `gawk -v xsize=$xsize 'BEGIN{print (xsize + 0.5)}'`
set cysize = `gawk -v ysize=$ysize 'BEGIN{print (ysize - 1)}'`
set drange=${cxsize}/0.5+w${cysize}/0.25+e0.5
gmt colorbar -Dx${drange} -Bx${dBA} -Ccolor.cpt --FONT=20p

gmt end
    
rm -f color.cpt
