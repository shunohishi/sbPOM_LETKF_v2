#!/bin/csh
#===============================================#
# Setting
#===============================================#

#---Machine
#set machine="fugaku"
#set machine="jss3"
set machine="rc"

#--Setting
@ ndat = 5
@ nobs = 365 #Minimum number of observation
@ id_ave = 556
@ id_worst = 49 
@ id_best = 40

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
# Figure setting
#=======================================================

set xsize=8; set ysize=4  #Figure size
set slon=105; set elon=255 #Longitude range
set slat=-5; set elat=65  #Latitude range
set plat=67 #Latitude (Plot label)
set BAx=a60f10; set BAy=a30f10; set BAl=WSne #Axis option
set label1=("(a) LORA-NP" "(b) BRAN2020" "(c) FORA-JPN60" "(d) GLORYS12V1" "(e) JCOPE-FGO")
set label2=("(f) Ensemble spread of LORA-NP" "(g) BRAN2020 vs. LORA-NP" "(h) FORA-JPN60 vs. LORA-NP" "(i) GLORYS12V1 vs. LORA-NP" "(j) JCOPE-FGO vs. LORA-NP")

#=======================================================
# Color setting
#=======================================================

#---RMSD
gmt makecpt -T0/0.2/0.025 -Clajolla -D -I > rmsd.cpt
set dBA_rmsd=a0.1f0.025+l"RMSD\040(m)"

#---Spread
cp rmsd.cpt sprd.cpt
set dBA_sprd=a0.1f0.025+l"Spread\040(m)"

#---RMSD ratio
gmt makecpt -T-40/40/10 -Croma -D -I > ratio.cpt
set dBA_ratio=a20f10+l"RMSD\040ratio\040(\045)"

#=======================================================
# DATA
#=======================================================

set input=dat/rmsd_ave.dat

#---RMSD
@ idat = 1
while(${idat} <= ${ndat})

    @ obs_col = 1 + 2 + ${idat}
    @ dat_col = 1 + 2 + ${ndat} + ${idat}
    gawk -v obs_col=${obs_col} -v dat_col=${dat_col} -v nobs=${nobs} -v out="rmsd${idat}.20" \
    '{if(nobs < $obs_col && $dat_col != -999) print $2,$3,$dat_col > out}' ${input}

    @ idat++

end    
    
#---RMSD ratio
@ idat = 2
while(${idat} <= ${ndat})

    @ obs_col = 1 + 2 + ${idat}
    @ dat1_col = 1 + 2 + ${ndat} + 1
    @ dat2_col = 1 + 2 + ${ndat} + ${idat}
    @ dif1_col = 1 + 2 + 2 * ${ndat} + ${idat}
    @ dif2_col = 1 + 2 + 3 * ${ndat} + ${idat}

    gawk -v obs_col=${obs_col} -v dat1_col=${dat1_col} -v dat2_col=${dat2_col} -v dif1_col=${dif1_col} -v dif2_col=${dif2_col} \
    -v nobs=${nobs} -v out="rmsd_sig${idat}.20" \
    '{if(nobs < $obs_col && $dat1_col != -999 && $dat2_col != -999 && $dif1_col != -999 && $dif2_col != -999 && $dif1_col*$dif2_col > 0) \
    print $2,$3,100*($dat2_col-$dat1_col)/$dat1_col > out}' ${input}

    gawk -v obs_col=${obs_col} -v dat1_col=${dat1_col} -v dat2_col=${dat2_col} -v dif1_col=${dif1_col} -v dif2_col=${dif2_col} \
    -v nobs=${nobs} -v out="rmsd_nosig${idat}.20" \
    '{if(nobs < $obs_col && $dat1_col != -999 && $dat2_col != -999 && $dif1_col != -999 && $dif2_col != -999 && $dif1_col*$dif2_col <= 0) \
    print $2,$3,100*($dat2_col-$dat1_col)/$dat1_col > out}' ${input}

    @ idat++

end

#---Representative Station location
gawk -v id=${id_ave}   -v out="rmsd_ave.20"   '{if($1 == id) print $2,$3 > out}' ${input}
gawk -v id=${id_worst} -v out="rmsd_worst.20" '{if($1 == id) print $2,$3 > out}' ${input}
gawk -v id=${id_best}  -v out="rmsd_best.20"  '{if($1 == id) print $2,$3 > out}' ${input}

#---Station and time averaged RMSD
set input=dat/rmsd_ave_all.dat
@ idat = 1
while(${idat} <= ${ndat})

    @ dat_col = ${idat}
    @ dif1_col = ${ndat} + ${idat}
    @ dif2_col = 2 * ${ndat} + ${idat}
    
    gawk -v dat_col=${dat_col} -v out="rmsd${idat}_ave_all.20" '{print $dat_col > out}' ${input}
    gawk -v dif1_col=${dif1_col} -v dif2_col=${dif2_col} -v out="rmsd${idat}_ave_info.20" \
    '{if($dif1_col != -999 && $dif2_col != -999 && $dif1_col*$dif2_col > 0) print 3 > out; else print 0 > out}' ${input} 

    @ idat++
end

#---Spread
set input=dat/sprd_ave.dat
@ idat = 1
while(${idat} <= ${ndat})

    @ obs_col = 1 + 2 + ${idat}
    @ dat_col = 1 + 2 + ${ndat} + ${idat}
    gawk -v obs_col=${obs_col} -v dat_col=${dat_col} -v nobs=${nobs} -v out="sprd${idat}.20" \
    '{if(nobs < $obs_col && $dat_col != -999) print $2,$3,$dat_col > out}' ${input}

    @ idat++

end    

#=======================================================
# Figure
#=======================================================

if(! -d fig) mkdir fig

set size=${xsize}d/${ysize}d
set range=${slon}/${elon}/${slat}/${elat}


gmt begin fig/tide_map png

#---RMSD
@ idat = 1
while(${idat} <= ${ndat})

    echo "RMSD $idat"

    #Box
    if($idat == 1)then
	gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X3 -Y25
    else
	set iy = `gawk -v ysize=$ysize 'BEGIN{print ysize + 1.5}'`
	gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -Y-${iy}
    endif

    #Coast
    gmt coast -Dl -W0.2,black -Ggray

    #Plot
    if(-f rmsd${idat}.20)then
	gmt psxy rmsd${idat}.20 -Sc0.20 -W0.2 -Crmsd.cpt
    endif

    #Label
    echo "${slon} ${plat} ${label1[$idat]}" | gmt text -F+f14p,0,black+jLB -N

    #Station location
    if(${idat} == 1)then
	gmt psxy rmsd_ave.20 -Sc0.20 -W2,green
	gmt psxy rmsd_best.20 -Sc0.20 -W2,yellow
	gmt psxy rmsd_worst.20 -Sc0.20 -W2,cyan
    endif

    #Station time averaged RMSD
    set input=rmsd${idat}_ave_all.20
    set ave=`gawk '{printf "%.4f", $1}' ${input}`
    set ave="RMSD (station-time): ${ave} m"
    set input=rmsd${idat}_ave_info.20
    set font=`gawk '{print $1}' ${input}`
    echo "${elon} ${slat} ${ave}" | gmt text -F+f12p,${font},black+jRB -N 
    
    @ idat++
    
end

#Color scale
set cxsize = `gawk -v xsize=$xsize 'BEGIN{print xsize - 1}'`
set drange = 0.5/-1+w${cxsize}/0.25+ef0.5+h
gmt colorbar -Dx${drange} -Bx${dBA_rmsd} -Crmsd.cpt --FONT=20p

#---Spread
@ idat = 1

echo "LORA ensemble spread"

#Box
set ix = `gawk -v xsize=$xsize 'BEGIN{print xsize + 2}'`
set iy = `gawk -v ysize=$ysize -v ndat=${ndat} 'BEGIN{print (ndat - 1) * (ysize + 1.5) }'`
gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X${ix} -Y${iy}

#Coast
gmt coast -Dl -W0.2,black -Ggray

#Plot
gmt psxy sprd${idat}.20 -Sc0.20 -W0.2 -Csprd.cpt

#Label
echo "${slon} ${plat} ${label2[$idat]}"|gmt text -F+f14p,0,black+jLB -N

#Color scale
set cxsize = `gawk -v xsize=$xsize 'BEGIN{print xsize + 0.5}'`
set cysize = `gawk -v ysize=$ysize 'BEGIN{print ysize - 1}'`
set drange = ${cxsize}/0.5+w${cysize}/0.25+ef0.5
gmt colorbar -Dx${drange} -Bx${dBA_sprd} -Csprd.cpt --FONT=20p

#---RMSD ratio
@ idat = 2
while($idat <= $ndat)

    echo "RMSD ratio $idat"

    #Box
    set iy = `gawk -v ysize=$ysize 'BEGIN{print ysize + 1.5 }'`
    gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -Y-${iy}

    #Coast
    gmt coast -Dl -W0.2,black -Ggray

    #Plot
    if(-f rmsd_nosig${idat}.20)then
	gmt psxy rmsd_nosig${idat}.20 -Sx0.20 -W2 -Cratio.cpt    
    endif
    
    if(-f rmsd_sig${idat}.20)then
	gmt psxy rmsd_sig${idat}.20 -Sc0.20 -W0.2 -Cratio.cpt    
    endif

    #Label
    echo "${slon} ${plat} ${label2[$idat]}" | gmt text -F+f14p,0,black+jLB -N

    #Sig label
    if(${idat} == 2)then
	cat << EOF > legend.20
N 2
S 0.15c c 0.2c gray 0.2p 0.4c Significant
S 0.15c x 0.2c - 2p 0.4c Not significant
EOF
	gmt legend legend.20 -DjBR+w5.5c+o0.1c -F+gwhite+p0.5p
    endif

    @ idat++
    
end

#Color scale
set cxsize = `gawk -v xsize=$xsize 'BEGIN{print xsize - 1}'`
set drange=0.5/-1+w${cxsize}/0.25+e0.5+h
gmt colorbar -Dx${drange} -Bx${dBA_ratio} -Cratio.cpt --FONT=20p

gmt end

rm -f *.20 *.cpt

