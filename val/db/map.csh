#!/bin/csh

#===============================================#

set nyr=16

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
    gmt set PS_CONVERT=I+m0.4/0.4/0.4/0.4 #WESN
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
set label=("(a) LORA" "(b) GLORYS2V4" "(c) ORAS5" "(d) C-GLORS" "(e) LORA's ensemble spread" "(f) (b)GLORYS2V4 vs. (a)LORA" "(g) (c)ORAS5 vs. (a)LORA" "(h) (d)C-GLORS vs. (a)LORA" "(i) Observation frequency")

#=======================================================
# Color
#=======================================================

gmt makecpt -T0/0.4/0.05 -Clajolla -D -I > urmsd.cpt
gmt makecpt -T0/0.4/0.05 -Clajolla -D -I > vrmsd.cpt
gmt makecpt -T0/1.5/0.25 -Clajolla -D -I > trmsd.cpt

gmt makecpt -T-40/40/10 -Cvik -D > uratio.cpt
gmt makecpt -T-40/40/10 -Cvik -D > vratio.cpt
gmt makecpt -T-40/40/10 -Cvik -D > tratio.cpt

gmt makecpt -T0/50/5 -Cbilbao -D -I > unobs.cpt 
gmt makecpt -T0/50/5 -Cbilbao -D -I > vnobs.cpt 
gmt makecpt -T0/50/5 -Cbilbao -D -I > tnobs.cpt 

#=======================================================
# Figure
#=======================================================


foreach var(u v t)

    gmt begin ${var} png

    @ i = 1
    @ ndata=4

    #---DATA
    set input=dat/${var}_bin.dat
    set nobs=0
    #RMSD
    gawk -v nobs=${nobs} '{if(nobs < $3) print $1,$2,$11 > "rmsd1.20"}' ${input}
    gawk -v nobs=${nobs} '{if(nobs < $4) print $1,$2,$12 > "rmsd2.20"}' ${input}
    gawk -v nobs=${nobs} '{if(nobs < $5) print $1,$2,$13 > "rmsd3.20"}' ${input}
    gawk -v nobs=${nobs} '{if(nobs < $6) print $1,$2,$14 > "rmsd4.20"}' ${input}

    #Spread
    gawk -v nobs=${nobs} '{if(nobs < $3) print $1,$2,$15 > "sprd1.20"}' ${input}
	
    #RMSD ratio
    gawk -v nobs=${nobs} '{if(nobs < $3 && nobs < $4) print $1,$2,($12-$11)/$11*100 > "ratio2.20"}' ${input}
    gawk -v nobs=${nobs} '{if(nobs < $3 && nobs < $5) print $1,$2,($13-$11)/$11*100 > "ratio3.20"}' ${input}
    gawk -v nobs=${nobs} '{if(nobs < $3 && nobs < $6) print $1,$2,($14-$11)/$11*100 > "ratio4.20"}' ${input}

    #nobs
    gawk -v nyr=${nyr} '{if(0 < $4) print $1,$2,$4/(nyr*12) > "nobs.20"}' ${input} #NOTE: Monthly frequency
    
    #Significant
    set input=dat/${var}tval_bin.dat
    gawk -v idat=2 '{if($3 == idat && $8 != -999 && $12 != -999 && $8*$8 < $12*$12) print $1,$2 > "sig2.20"}' ${input}
    gawk -v idat=3 '{if($3 == idat && $8 != -999 && $12 != -999 && $8*$8 < $12*$12) print $1,$2 > "sig3.20"}' ${input}
    gawk -v idat=4 '{if($3 == idat && $8 != -999 && $12 != -999 && $8*$8 < $12*$12) print $1,$2 > "sig4.20"}' ${input}
    
    foreach dat(rmsd1 rmsd2 rmsd3 rmsd4 sprd1 ratio2 ratio3 ratio4 nobs)
    
    echo ${var} ${dat}

    #---DATA
    gmt xyz2grd ${dat}.20 -G${dat}.grd -R${range2} -I${int1}

    #---Color
    if(${dat} == "rmsd1" || ${dat} == "rmsd2" || ${dat} == "rmsd3" || ${dat} == "rmsd4")then
	set color="${var}rmsd"
	set drange=0.5/-1+w7/0.25+ef0.5+h
	if(${var} == "u" || ${var} == "v") set dBA=a0.10f0.05+l"RMSD\040(m/s)"
	if(${var} == "t") set dBA=a0.50f0.25+l"RMSD\040(\260C)"
    else if(${dat} == "sprd1")then
	set color="${var}rmsd"
	set drange=8.5/0.5+w3/0.25+ef0.5
	if(${var} == "u" || ${var} == "v") set dBA=a0.10f0.05+l"Spread\040(m/s)"
	if(${var} == "t") set dBA=a0.50f0.25+l"Spread\040(\260C)"    
    else if(${dat} == "nobs")then
	set color="${var}nobs"
	set drange=8.5/0.5+w3/0.25+ef0.5
	set dBA=a10f5+l"Obs\040frequency\040(month@+-1@+)"
    else
	set color="${var}ratio"
	set drange=0.5/-1+w7/0.25+e0.5+h
	set dBA=a20f10+l"RMSD\040ratio\040(\045)"
    endif
    
    #---Figure
    if($i == 1)then
	gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl} -X3 -Y20
    else if($i == $ndata + 1)then
	gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl} -X10 -Y16.5    
    else if($i % $ndata == 1)then
	gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl} -X11.5 -Y16.5
    else
	gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl} -Y-5.5
    endif

    gmt psmask ${dat}.20 -R${range2} -I${int1}
    gmt grdimage ${dat}.grd -C${color}.cpt
    gmt psmask -C

    if(${dat} == "ratio2" || ${dat} == "ratio3" || ${dat} == "ratio4")then
	if(${dat} == "ratio2") set sig="sig2"
	if(${dat} == "ratio3") set sig="sig3"
	if(${dat} == "ratio4") set sig="sig4"
	if(-f ${sig}.20) gmt psxy ${sig}.20 -Sx0.075 -Gblack
    endif
    
    gmt coast -R${range1} -Dl -W0.2,black -Gwhite

    echo "0 95 ${label[$i]}" | gmt text -F+f14p,0,black+jLB -N     

    set input=dat/${var}_ave.dat
    set ave=""
    if(${dat} == "rmsd1") set ave=`awk '{printf "%.3f", $9}' ${input}`
    if(${dat} == "rmsd2") set ave=`awk '{printf "%.3f", $10}' ${input}`	
    if(${dat} == "rmsd3") set ave=`awk '{printf "%.3f", $11}' ${input}`
    if(${dat} == "rmsd4") set ave=`awk '{printf "%.3f", $12}' ${input}`
     
    if(${var} == "u" || ${var} == "v") set ave="Avg.: ${ave} m/s"
    if(${var} == "t") set ave="Avg.:${ave} \260C"
    
    set input=dat/${var}tval_ave.dat
    if(${dat} == "rmsd1")then
	set font=`gawk -v idat=1 '{if ($1 == idat && $6 != -999 && $10 != -999 && $6*$6 < $10*$10) print 3; else if ($1 == idat) print 0}' ${input}`
    else if(${dat} == "rmsd2")then
	set font=`gawk -v idat=2 '{if ($1 == idat && $6 != -999 && $10 != -999 && $6*$6 < $10*$10) print 3; else if ($1 == idat) print 0}' ${input}`
    else if(${dat} == "rmsd3")then
    	set font=`gawk -v idat=3 '{if ($1 == idat && $6 != -999 && $10 != -999 && $6*$6 < $10*$10) print 3; else if ($1 == idat) print 0}' ${input}`
    else if(${dat} == "rmsd4")then
	set font=`gawk -v idat=4 '{if ($1 == idat && $6 != -999 && $10 != -999 && $6*$6 < $10*$10) print 3; else if ($1 == idat) print 0}' ${input}`
    endif        
    
    if(${dat} == "rmsd1" || ${dat} == "rmsd2" || ${dat} == "rmsd3" || ${dat} == "rmsd4")then    
	echo "360 -90 ${ave}" | gmt text -F+f14p,${font},black+jRB -N     
    endif

    echo "0 95 ${label[$i]}" | gmt text -F+f14p,0,black+jLB -N     
    
    if($i % $ndata == 0 || ($i != 1 && $i % $ndata == 1))then
	gmt colorbar -Dx${drange} -Bx${dBA} -C${color}.cpt --FONT=20p
    endif

    @ i++
    
    end

    gmt end

    rm -f *.20 *.grd
    
end    

rm -f *.cpt
