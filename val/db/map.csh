#!/bin/csh

#===============================================#

set nyr=16 #Year
set n=4 #Data

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

set size=8d/4d
set range1=0/360/-90/90
set range2=2.5/357.5/-87.5/87.5 #Describe exact edge points
set BAx=a60f10
set BAy=a30f10
set BAl=WSne
set int1=5/5
set int2=10/10
set label1=("(a) LORA" "(b) GLORYS2V4" "(c) ORAS5" "(d) C-GLORS")
set label2=("(e)" "(f) (b)GLORYS2V4 vs. (a)LORA" "(g) (c)ORAS5 vs. (a)LORA" "(h) (d)C-GLORS vs. (a)LORA")

#=======================================================
# Color
#=======================================================

gmt makecpt -T-0.1/0.1/0.01 -Croma -D -I > ubias.cpt
gmt makecpt -T-0.1/0.1/0.01 -Croma -D -I > vbias.cpt
gmt makecpt -T-0.2/0.2/0.02 -Croma -D -I > tbias.cpt

gmt makecpt -T0/0.4/0.05 -Clajolla -D -I > urmsd.cpt
gmt makecpt -T0/0.4/0.05 -Clajolla -D -I > vrmsd.cpt
gmt makecpt -T0/1.5/0.25 -Clajolla -D -I > trmsd.cpt

cp urmsd.cpt usprd1.cpt
cp vrmsd.cpt vsprd1.cpt
cp trmsd.cpt tsprd1.cpt

#gmt makecpt -T-100/100/25 -Cvik -D > ubias_ratio.cpt
#gmt makecpt -T-100/100/25 -Cvik -D > vbias_ratio.cpt
#gmt makecpt -T-100/100/25 -Cvik -D > tbias_ratio.cpt

gmt makecpt -T-0.1/0.1/0.01 -Cvik -D > ubias_dif.cpt
gmt makecpt -T-0.1/0.1/0.01 -Cvik -D > vbias_dif.cpt
gmt makecpt -T-0.2/0.2/0.02 -Cvik -D > tbias_dif.cpt

gmt makecpt -T-40/40/10 -Cvik -D > urmsd_ratio.cpt
gmt makecpt -T-40/40/10 -Cvik -D > vrmsd_ratio.cpt
gmt makecpt -T-40/40/10 -Cvik -D > trmsd_ratio.cpt

gmt makecpt -T0/50/5 -Cbilbao -D -I > unobs.cpt 
gmt makecpt -T0/50/5 -Cbilbao -D -I > vnobs.cpt 
gmt makecpt -T0/50/5 -Cbilbao -D -I > tnobs.cpt 

#=======================================================
# Figure
#=======================================================


foreach var(u v t)

    #---Unit
    if(${var} == "u")then
	set unit="m/s"
	set title="Sea-surface zonal velocity"
    else if(${var} == "v")then
	set unit="m/s"
	set title="Sea-surface meridional velocity"
    else if(${var} == "t")then
	set unit="\260C"
	set title="Sea-surface temperature"
    endif

    #---DATA
    set nobs=0

    #Bias, RMSD, Spread
    #echo "Bias, RMSD, and Spread"
    foreach index(bias rmsd sprd)
    
	set input=dat/${var}${index}_bin.dat
	
	gawk -v nobs=${nobs} -v out="${index}1.20" '{if(nobs < $3) print $1,$2,$7 > out}' ${input}
	gawk -v nobs=${nobs} -v out="${index}2.20" '{if(nobs < $4) print $1,$2,$8 > out}' ${input}
	gawk -v nobs=${nobs} -v out="${index}3.20" '{if(nobs < $5) print $1,$2,$9 > out}' ${input}
	gawk -v nobs=${nobs} -v out="${index}4.20" '{if(nobs < $6) print $1,$2,$10 > out}' ${input}
	
	@ i = 1
	while($i <= $n)
	    if(-f ${index}${i}.20) gmt xyz2grd ${index}${i}.20 -G${index}${i}.grd -R${range2} -I${int1}
	    @ i++
	end
	
	set input=dat/${var}${index}_ave.dat
	gawk -v out=${index}1_ave.20 '{printf "%.3f", $5 > out}' ${input}
	gawk -v out=${index}2_ave.20 '{printf "%.3f", $6 > out}' ${input}
	gawk -v out=${index}3_ave.20 '{printf "%.3f", $7 > out}' ${input}
	gawk -v out=${index}4_ave.20 '{printf "%.3f", $8 > out}' ${input}
	
	if(${index} == "bias" || ${index} == "rmsd")then
	    set input=dat/${var}${index}_tval_ave.dat
	    gawk -v idat=1 -v out=${index}1_ave_sig.20 \
	    '{if ($1 == idat && $6 != -999 && $10 != -999 && $6*$6 < $10*$10) print 3 > out; else if ($1 == idat) print 0 > out}' ${input}
	    gawk -v idat=2 -v out=${index}2_ave_sig.20 \
	    '{if ($1 == idat && $6 != -999 && $10 != -999 && $6*$6 < $10*$10) print 3 > out; else if ($1 == idat) print 0 > out}' ${input}
	    gawk -v idat=3 -v out=${index}3_ave_sig.20 \
	    '{if ($1 == idat && $6 != -999 && $10 != -999 && $6*$6 < $10*$10) print 3 > out; else if ($1 == idat) print 0 > out}' ${input}
	    gawk -v idat=4 -v out=${index}4_ave_sig.20 \
	    '{if ($1 == idat && $6 != -999 && $10 != -999 && $6*$6 < $10*$10) print 3 > out; else if ($1 == idat) print 0 > out}' ${input}
	endif
	    
    end #index
	    
    #Absolute bias ratio, rmsd ratio
    #echo "Absolute bias & RMSD ratio"
    foreach index(bias rmsd)

    	set input=dat/${var}${index}_bin.dat    

	if(${index} == "bias")then
	    set type="dif"
	    gawk -v nobs=${nobs} -v out="${index}_${type}2.20" \
	    '{if(nobs < $3 && nobs < $4) print $1,$2,sqrt($8*$8)-sqrt($7*$7) > out}' ${input}
	    gawk -v nobs=${nobs} -v out="${index}_${type}3.20" \
	    '{if(nobs < $3 && nobs < $5) print $1,$2,sqrt($9*$9)-sqrt($7*$7) > out}' ${input}
	    gawk -v nobs=${nobs} -v out="${index}_${type}4.20" \
	    '{if(nobs < $3 && nobs < $6) print $1,$2,sqrt($10*$10)-sqrt($7*$7) > out}' ${input}
	else if(${index} == "rmsd")then
	    set type="ratio"
	    gawk -v nobs=${nobs} -v out="${index}_${type}2.20" \
	    '{if(nobs < $3 && nobs < $4 && $7 != 0.) print $1,$2,(sqrt($8*$8)-sqrt($7*$7))/sqrt($7*$7)*100 > out}' ${input}
	    gawk -v nobs=${nobs} -v out="${index}_${type}3.20" \
	    '{if(nobs < $3 && nobs < $5 && $7 != 0.) print $1,$2,(sqrt($9*$9)-sqrt($7*$7))/sqrt($7*$7)*100 > out}' ${input}
	    gawk -v nobs=${nobs} -v out="${index}_${type}4.20" \
	    '{if(nobs < $3 && nobs < $6 && $7 != 0.) print $1,$2,(sqrt($10*$10)-sqrt($7*$7))/sqrt($7*$7)*100 > out}' ${input}
	endif
	    
	@ i = 2
	while($i <= $n)
	    if(-f ${index}_${type}${i}.20) gmt xyz2grd ${index}_${type}${i}.20 -G${index}_${type}${i}.grd -R${range2} -I${int1}
	    @ i++
	end	
    end
	    
    #Significant
    #echo "Statistical test"
    foreach index(bias rmsd)
	set input=dat/${var}${index}_tval_bin.dat
	gawk -v idat=2 -v out="${index}_sig2.20" \
	'{if($3 == idat && $8 != -999 && $12 != -999 && $8*$8 < $12*$12) print $1,$2 > out}' ${input}
	gawk -v idat=3 -v out="${index}_sig3.20" \
	'{if($3 == idat && $8 != -999 && $12 != -999 && $8*$8 < $12*$12) print $1,$2 > out}' ${input}
	gawk -v idat=4 -v out="${index}_sig4.20" \
	'{if($3 == idat && $8 != -999 && $12 != -999 && $8*$8 < $12*$12) print $1,$2 > out}' ${input}
    end

    #Number of obs.
    #echo "Number of observation"
    set input=dat/${var}rmsd_bin.dat
    gawk -v nobs=${nobs} -v nyr=${nyr} -v out="nobs.20" \
    '{if(nobs < $4) print $1,$2,$4/(nyr*12) > out}' ${input} #NOTE: Monthly frequency and GLORYS
    if(-f nobs.20) gmt xyz2grd nobs.20 -Gnobs.grd -R${range2} -I${int1}
    
    #---Figure
    foreach index(bias rmsd)

	echo "Make ${var}${index}.png"
    
	gmt begin ${var}${index} png
	
	#---Color setting
	if((${var} == "u" || ${var} == "v") && ${index} == "bias")then
	    set drange=0.5/-1+w7/0.25+e0.5+h
	    set dBA=a0.05f0.01+l"Bias\040(m/s)"
	else if(${var} == "t" && ${index} == "bias")then
	    set drange=0.5/-1+w7/0.25+e0.5+h
	    set dBA=a0.1f0.02+l"Bias\040(\040C)"
	else if((${var} == "u" || ${var} == "v") && ${index} == "rmsd")then
	    set drange=0.5/-1+w7/0.25+ef0.5+h
	    set dBA=a0.10f0.05+l"RMSD\040(m/s)"
	else if(${var} == "t" && ${index} == "rmsd")then
	    set drange=0.5/-1+w7/0.25+ef0.5+h
	    set dBA=a0.50f0.25+l"RMSD\040(\260C)"
	endif
    
	#---Bias/RMSD
	@ i = 1
	while($i <= $n)

	    if($i == 1)then
		gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl}+t"${title}" -X3 -Y20 --FONT_TITLE=16p
	    else
		gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl} -Y-5.5
	    endif
	    
	    gmt psmask ${index}${i}.20 -R${range2} -I${int1}
	    gmt grdimage ${index}${i}.grd -C${var}${index}.cpt
	    gmt psmask -C

	    gmt coast -R${range1} -Dl -W0.2,black -Gwhite

	    echo "0 95 ${label1[$i]}" | gmt text -F+f14p,0,black+jLB -N	    

	    set input=${index}${i}_ave.20
	    set ave=`gawk '{printf "%.3f", $1}' ${input}`
	    if(${var} == "u" || ${var} == "v")then
		set ave="Avg.: ${ave} m/s"
	    else if(${var} == "t")then
		set ave="Avg.:${ave} \260C"
	    endif
	    
	    set input=${index}${i}_ave_sig.20
	    set font=`gawk '{printf $1}' ${input}`
	    echo "360 -90 ${ave}" | gmt text -F+f14p,${font},black+jRB -N     

	    if($i == $n)then
		gmt colorbar -Dx${drange} -Bx${dBA} -C${var}${index}.cpt --FONT=20p
	    endif

	    @ i++
	    
	end #i

	#---Color setting
	if((${var} == "u" || ${var} == "v") && ${index} == "bias")then
	    set drange=0.5/-1+w7/0.25+e0.5+h
	    set dBA=a0.05f0.01+l"Absolute\040bias\040difference\040(m/s)"
	else if(${var} == "t" && ${index} == "bias")then
	    set drange=0.5/-1+w7/0.25+e0.5+h
	    set dBA=a0.1f0.02+l"Absolute\040bias\040difference\040(\260C)"
	else if((${var} == "u" || ${var} == "v") && ${index} == "rmsd")then
	    set drange=0.5/-1+w7/0.25+e0.5+h
	    set dBA=a20f10+l"RMSD\040ratio\040(\045)"
	else if(${var} == "t" && ${index} == "rmsd")then
	    set drange=0.5/-1+w7/0.25+e0.5+h
	    set dBA=a20f10+l"RMSD\040ratio\040(\045)"
	endif	
	

	#---Absolute bias difference/RMSD ratio
	@ i = 2
	if(${index} == "bias")then
	    set type="dif"
	else if(${index} == "rmsd")then
	    set type="ratio"
	endif
	while($i <= $n)

	    if($i == 2)then
		gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl} -X10.5 -Y11
	    else
		gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl} -Y-5.5
	    endif
	    
	    gmt psmask ${index}_${type}${i}.20 -R${range2} -I${int1}
	    gmt grdimage ${index}_${type}${i}.grd -C${var}${index}_${type}.cpt
	    gmt psmask -C

	    if(-f ${index}_sig${i}.20)then
		gmt psxy ${index}_sig${i}.20 -Sx0.075 -Gblack
	    endif

	    gmt coast -R${range1} -Dl -W0.2,black -Gwhite

	    echo "0 95 ${label2[$i]}" | gmt text -F+f14p,0,black+jLB -N	    

	    if($i == $n)then
		gmt colorbar -Dx${drange} -Bx${dBA} -C${var}${index}_${type}.cpt --FONT=20p
	    endif

	    @ i++
	    
	end #i


	#---Additional figure
	if(${index} == "bias")then
	    set dat="nobs"
	    set label="(e) Observation frequency"
	else if(${index} == "rmsd")then
	    set dat="sprd1"
	    set label="(e) LORA's ensemble spread"
	endif

	if(${index} == "bias")then
	    set dBA=a25f5+l"Observation\040frequency\040(month@+-1@+)"
	else if((${var} == "u" || ${var} == "v") && ${index} == "rmsd")then
	    set drange=0.5/-1+w7/0.25+ef0.5+h
	    set dBA=a0.10f0.05+l"Ensemble\040spread\040(m/s)"
	else if(${var} == "t" && ${index} == "rmsd")then
	    set drange=0.5/-1+w7/0.25+ef0.5+h
	    set dBA=a0.50f0.25+l"Ensemble\040spread\040(\260C)"
	endif
	set drange=8.5/0.5+w3/0.25+ef0.5	
	
	gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl} -Y16.5	

	gmt psmask ${dat}.20 -R${range2} -I${int1}
	gmt grdimage ${dat}.grd -C${var}${dat}.cpt
	gmt psmask -C
	
	gmt coast -R${range1} -Dl -W0.2,black -Gwhite

	gmt colorbar -Dx${drange} -Bx${dBA} -C${var}${dat}.cpt --FONT=20p
	
	echo "0 95 ${label}" | gmt text -F+f14p,0,black+jLB -N	    
	
	gmt end
	
    end #index

    rm -f *.20 *.grd
    
end #var

rm -f *.cpt
