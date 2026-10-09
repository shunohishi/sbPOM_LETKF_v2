#!/bin/csh

#===============================================
# Setting
#===============================================

set obsa=1 #Obs. availability limit (Unit: month^-1) 
set nyr=18 #Year --> Plot where the obs frequency > ${obsa}
set ndat=5 #Data

#===============================================
# Machine
#===============================================

#set machine="fugaku" #Fugaku
#set machine="jss3"    #JSS3
set machine="rc"      #RCCS Cloud

#===============================================
# Spack load GMT6
#===============================================

if(${machine} == "fugaku")then
    setenv SPACK_ROOT /vol0004/apps/oss/spack
    source /vol0004/apps/oss/spack/share/spack/setup-env.csh
    spack load /mnrvuuq
endif

#===============================================
# Option
#===============================================

rm -f gmt.conf

gmt set MAP_FRAME_TYPE=PLAIN FORMAT_GEO_MAP=dddF
gmt set FONT=10p,Helvetica,black
gmt set GMT_AUTO_DOWNLOAD=off

if(${machine} == "jss3")then
    gmt set PS_CONVERT=I+m0.4/0.4/0.4/0.8 #WESN
endif
    
#===============================================
# Figure setting
#===============================================

set xsize=8; set ysize=4   #Figure size                             
set slon=110; set elon=250 #Longitude
set slat=15;  set elat=60  #Latitude
set res=5                  #Bin resolution
set BAx=a30f10; set BAy=a20f10; set BAl=WSne #Axis option
set label1=("(a) LORA-NP" "(b) BRAN2020" "(c) FORA-JPN60" "(d) GLORYS12V1" "(e) JCOPE-FGO") #Label (Bias/RMSD)
set label2=("(f) xxx" "(g) BRAN2020 vs. LORA-NP" "(h) FORA-JPN60 vs. LORA-NP" "(i) GLORYS12V1 vs. LORA-NP" "(j) JCOPE-FGO vs. LORA-NP") #Label (Abs. bias dif./RMSD ratio)
set label2_bias=("(f) Observation availability")   #Label (Obs. avail)
set label2_rmsd=("(f) Ensemble spread of LORA-NP") #Label (Ens. spread)
set ypos=2                                         #Label Y Position

#=======================================================
# Color setting
#=======================================================

#---Bias
gmt makecpt -T-0.10/0.10/0.01 -Cbam -D -I > ubias.cpt
gmt makecpt -T-0.10/0.10/0.01 -Cbam -D -I > vbias.cpt
gmt makecpt -T-0.20/0.20/0.02 -Cbam -D -I > tbias.cpt

set dBA_ubias=a0.05f0.01+l"Bias\040(m/s)"
set dBA_vbias=a0.05f0.01+l"Bias\040(m/s)"
set dBA_tbias=a0.10f0.02+l"Bias\040(\260C)"

#---RMSD
gmt makecpt -T0/0.40/0.025 -Clajolla -D -I > urmsd.cpt
gmt makecpt -T0/0.40/0.025 -Clajolla -D -I > vrmsd.cpt
gmt makecpt -T0/1.50/0.125 -Clajolla -D -I > trmsd.cpt

set dBA_urmsd=a0.10f0.025+l"RMSD\040(m/s)"
set dBA_vrmsd=a0.10f0.025+l"RMSD\040(m/s)"
set dBA_trmsd=a0.50f0.125+l"RMSD\040(\260C)"

#---Ensemble spread
cp urmsd.cpt usprd.cpt
cp vrmsd.cpt vsprd.cpt
cp trmsd.cpt tsprd.cpt

set dBA_usprd=a0.10f0.05+l"Ensemble\040spread\040(m/s)"
set dBA_vsprd=a0.10f0.05+l"Ensemble\040spread\040(m/s)"
set dBA_tsprd=a0.50f0.25+l"Ensemble\040spread\040(\260C)"

#---Absolute bias difference
gmt makecpt -T-0.10/0.10/0.01 -Cvik -D > ubias_dif.cpt
gmt makecpt -T-0.10/0.10/0.01 -Cvik -D > vbias_dif.cpt
gmt makecpt -T-0.20/0.20/0.02 -Cvik -D > tbias_dif.cpt

set dBA_udif=a0.05f0.01+l"Absolute\040bias\040difference\040(m/s)"
set dBA_vdif=a0.05f0.01+l"Absolute\040bias\040difference\040(m/s)"
set dBA_tdif=a0.10f0.02+l"Absolute\040bias\040difference\040(\260C)"

#---RMSD ratio
gmt makecpt -T-40/40/5 -Cvik -D > urmsd_ratio.cpt
gmt makecpt -T-40/40/5 -Cvik -D > vrmsd_ratio.cpt
gmt makecpt -T-40/40/5 -Cvik -D > trmsd_ratio.cpt

set dBA_uratio=a20f5+l"RMSD\040ratio\040(\045)"
set dBA_vratio=a20f5+l"RMSD\040ratio\040(\045)"
set dBA_tratio=a20f5+l"RMSD\040ratio\040(\045)"


#---Observation availability
gmt makecpt -T0/50/10 -Cbilbao -D -I > uobsa.cpt 
gmt makecpt -T0/50/10 -Cbilbao -D -I > vobsa.cpt 
gmt makecpt -T0/50/10 -Cbilbao -D -I > tobsa.cpt

set dBA_uobsa=a10f10+l"Observation\040availability\040(month@+-1@+)"
set dBA_vobsa=a10f10+l"Observation\040availability\040(month@+-1@+)"
set dBA_tobsa=a10f10+l"Observation\040availability\040(month@+-1@+)"

#=======================================================
# Figure
#=======================================================

if(! -d fig) mkdir fig

set size=${xsize}d/${ysize}d
set range1=${slon}/${elon}/${slat}/${elat}
set slon_bin = `gawk -v slon=$slon -v res=$res 'BEGIN{print slon + res/2}'`; set elon_bin = `gawk -v elon=$elon -v res=$res 'BEGIN{print elon - res/2}'`
set slat_bin = `gawk -v slat=$slat -v res=$res 'BEGIN{print slat + res/2}'`; set elat_bin = `gawk -v elat=$elat -v res=$res 'BEGIN{print elat - res/2}'`
set range2=${slon_bin}/${elon_bin}/${slat_bin}/${elat_bin}
set int=${res}/${res}

foreach var(u v t)

    #---Unit
    if(${var} == "u")then
	set unit="m/s"
	set title="Surface zonal velocity"
    else if(${var} == "v")then
	set unit="m/s"
	set title="Surface meridional velocity"
    else if(${var} == "t")then
	set unit="\260C"
	set title="Sea-surface temperature"
    endif

    #---DATA----------------------------------------------------------------------------------------------------
    #---Obs. availability
    set input=dat/${var}rmsd_bin.dat
    @ idat = 1
    while(${idat} <= ${ndat})

	@ obs_col = 2 + ${idat}
	    
	gawk -v obsa=${obsa} -v nyr=${nyr} -v obs_col=${obs_col} -v out="obsa${idat}.20" \
	'{if(obsa <= $obs_col/(nyr*12)) print $1,$2,$obs_col/(nyr*12) > out}' ${input}

	if(-f obsa${idat}.20)then
	    (gmt xyz2grd obsa${idat}.20 -Gobsa${idat}.grd -R${range2} -I${int} &)
	endif

	@ idat++
    end
	
    #---Bias, RMSD, Spread
    #echo "Bias, RMSD, and Spread"
    foreach index(bias rmsd sprd)

	#Bin
	set input=dat/${var}${index}_bin.dat
	@ idat = 1
	while(${idat} <= ${ndat})
	
	    @ obs_col = 2 + ${idat}
	    @ dat_col = 2 + $ndat + $idat
	    gawk -v obsa=${obsa} -v nyr=${nyr} \
	    -v obs_col=${obs_col} -v dat_col=${dat_col} \
	    -v out="${index}${idat}.20" \
	    '{if(obsa <= $obs_col/(nyr*12)) print $1,$2,$dat_col > out}' ${input}
	    
	    if(-f ${index}${idat}.20)then
		(gmt xyz2grd ${index}${idat}.20 -G${index}${idat}.grd -R${range2} -I${int} &)
	    endif
	    
	    @ idat++

	end

        #Spatiotemporal Average				
	set input=dat/${var}${index}_ave.dat
	@ idat = 1
	while(${idat} <= ${ndat})

	    @ dat_col = ${ndat} + ${idat}
	    gawk -v dat_col=${dat_col} -v out=${index}${idat}_ave.20 \
	    '{printf "%.3f", $dat_col > out}' ${input}
	    
	    @ idat++

	end

	#Significant test
	if(${index} == "bias" || ${index} == "rmsd")then

	    set input=dat/${var}${index}_dif_ave.dat
	    @ idat = 1
	    while(${idat} <= ${ndat})
	    
		@ dif1_col = 1 + 1
		@ dif2_col = 1 + ${ndat} + 1
		gawk -v idat=${idat} -v dif1_col=${dif1_col} -v dif2_col=${dif2_col} -v out=${index}${idat}_ave_sig.20 \
		'{if ($1 == idat && $dif1_col != -999 && $dif2_col != -999 && $dif1_col * $dif2_col > 0.) print 3 > out; else if ($1 == idat) print 0 > out}' ${input}

		@ idat++
	    end
	    
	endif
	    	    
    end #index
    wait
	    
    #---Absolute bias diffrence, rmsd ratio
    #echo "Absolute bias difference & RMSD ratio"
    foreach index(bias rmsd)

	#Absolute bias difference/RMSD ratio
    	set input=dat/${var}${index}_bin.dat
	@ idat = 2
	while(${idat} <= ${ndat})

	    @ obs1_col = 2 + 1
	    @ obs2_col = 2 + ${idat}
	    @ dat1_col = 2 + ${ndat} + 1
	    @ dat2_col = 2 + ${ndat} + ${idat}

	    if(${index} == "bias")then
		set type="dif"
		gawk -v obsa=${obsa} -v nyr=${nyr} \
		     -v obs1_col=${obs1_col} -v obs2_col=${obs2_col} \
		     -v dat1_col=${dat1_col} -v dat2_col=${dat2_col} \
		     -v out="${index}_${type}${idat}.20" \
		    '{if(obsa <= $obs1_col/(nyr*12) && obsa <= $obs2_col/(nyr*12)) print $1,$2,sqrt($dat2_col*$dat2_col)-sqrt($dat1_col*$dat1_col) > out}' ${input}
	    else if(${index} == "rmsd")then
		set type="ratio"
		gawk -v obsa=${obsa} -v nyr=${nyr} \
		     -v obs1_col=${obs1_col} -v obs2_col=${obs2_col} \
		     -v dat1_col=${dat1_col} -v dat2_col=${dat2_col} \
		     -v out="${index}_${type}${idat}.20" \
		    '{if(obsa <= $obs1_col/(nyr*12) && obsa <= $obs2_col/(nyr*12) && $dat1_col != 0.) print $1,$2,(sqrt($dat2_col*$dat2_col)-sqrt($dat1_col*$dat1_col))/sqrt($dat1_col*$dat1_col)*100 > out}' ${input}
	    endif

	    if(-f ${index}_${type}${idat}.20)then
		(gmt xyz2grd ${index}_${type}${idat}.20 -G${index}_${type}${idat}.grd -R${range2} -I${int} &)
	    endif
	    
	    @ idat++
	end

	#Significant test
	set input=dat/${var}${index}_dif_bin.dat
	@ idat = 2
	while(${idat} <= ${ndat})

	    @ obs1_col = 3 + 1
	    @ obs2_col = 3 + ${idat}
	    @ dif1_col = 3 + ${ndat} + 1
	    @ dif2_col = 3 + 2 * ${ndat} + 1
	    
	    gawk -v idat=${idat} -v obsa=${obsa} -v nyr=${nyr} \
	    -v obs1_col=${obs1_col} -v obs2_col=${obs2_col} -v dif1_col=${dif1_col} -v dif2_col=${dif2_col} \
	    -v out="${index}_sig${idat}.20" \
	    '{if($3 == idat && obsa <= $obs1_col/(nyr*12) && obsa <= $obs2_col/(nyr*12) && $dif1_col != -999 && $dif2_col != -999 && $dif1_col*$dif2_col > 0.) print $1,$2 > out}' ${input}

	    @ idat++

	end
	
    end #index
    wait
    	        
    #---Figure-----------------------------------------------------------------------------------------------------
    foreach index(bias rmsd)

	echo "Make ${var}${index}.png"
    
	gmt begin fig/${var}${index} png
			    
	#---1. Bias/RMSD
	#Color setting
	set cxsize = `gawk -v xsize=$xsize 'BEGIN{print xsize - 1}'`
	if(${index} == "bias")then
	    set drange=0.5/-1+w${cxsize}/0.25+e0.5+h
	else if(${index} == "rmsd")then
	    set drange=0.5/-1+w${cxsize}/0.25+ef0.5+h
	endif
	eval set dBA = \$dBA_${var}${index}	

	#Figure
	@ idat = 1
	while(${idat} <= ${ndat})

	    #Box
	    if($idat == 1)then
		gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl}+t"${title}" -X3 -Y25 --FONT_TITLE=16p
	    else
		set iy = `gawk -v ysize=$ysize 'BEGIN{print -1 * (ysize + 1.5)}'`
		gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl} -Y${iy}
	    endif

	    #Bias/RMSD
	    gmt psmask ${index}${idat}.20 -R${range2} -I${int}
	    gmt grdimage ${index}${idat}.grd -C${var}${index}.cpt
	    gmt psmask -C

	    #Coast
	    gmt coast -R${range1} -Dl -W0.2,black -Gwhite

	    #Label1
	    @ ix = ${slon}; @ iy = ${elat} + ${ypos}
	    echo "${ix} ${iy} ${label1[$idat]}" | gmt text -F+f14p,0,black+jLB -N	    

	    #Spatiotemporal ave	    
	    if(${index} == "rmsd")then
	    	set input=${index}${idat}_ave.20
		set ave=`gawk '{printf "%.3f", $1}' ${input}`
		set ave="RMSD (space-time): ${ave} ${unit}"
		
		set input=${index}${idat}_ave_sig.20
		set font=`gawk '{printf $1}' ${input}`

		set ix = ${elon}; set iy = ${slat}
		echo "${ix} ${iy} ${ave}" | gmt text -F+f12p,${font},black+jRB -N
		
	    endif

	    #Color scale
	    if(${idat} == ${ndat})then
		gmt colorbar -Dx${drange} -Bx${dBA} -C${var}${index}.cpt --FONT=20p
	    endif

	    @ idat++
	    
	end #idat

	#---2. Obs. vail (Bias)/Ensemble spread (RMSD)
	@ idat = 1
	if(${index} == "bias")then
	    set dat="obsa"
	    set label="${label2_bias}"
	else if(${index} == "rmsd")then
	    set dat="sprd"
	    set label="${label2_rmsd}"
	endif
	eval set dBA = \$dBA_${var}${dat}
	set ix = `gawk -v xsize=$xsize 'BEGIN{print xsize + 0.5}'`
	set cysize = `gawk -v ysize=$ysize 'BEGIN{print ysize - 1}'`
	set drange=${ix}/0.5+w${cysize}/0.25+ef0.5

	#Box
	set ix = `gawk -v xsize=$xsize 'BEGIN{print xsize + 2.5}'`
	set iy = `gawk -v ysize=$ysize -v ndat=$ndat 'BEGIN{print (ysize + 1.5) * (ndat - 1)}'`	
	gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl} -X${ix} -Y${iy}

	#Obs. avail/Ens. Spead
	gmt psmask ${dat}${idat}.20 -R${range2} -I${int}
	gmt grdimage ${dat}${idat}.grd -C${var}${dat}.cpt
	gmt psmask -C

	#Coast
	gmt coast -R${range1} -Dl -W0.2,black -Gwhite

	#Colorbar
	gmt colorbar -Dx${drange} -Bx${dBA} -C${var}${dat}.cpt --FONT=20p

	@ ix = ${slon}; @ iy = ${elat} + ${ypos}
	echo "${ix} ${iy} ${label}" | gmt text -F+f14p,0,black+jLB -N	    


	#---3. Absolute bias difference/RMSD ratio	
	if(${index} == "bias")then
	    set type="dif"
	else if(${index} == "rmsd")then
	    set type="ratio"
	endif

	#Color setting	
	set drange=0.5/-1+w${cxsize}/0.25+e0.5+h
	eval set dBA = \$dBA_${var}${type}

	#Figure
	@ idat = 2
	while(${idat} <= ${ndat})

	    #Box
	    set iy = `gawk -v ysize=$ysize 'BEGIN{print -1 * (ysize + 1.5)}'`
	    gmt basemap -JX${size} -R${range1} -Bx${BAx} -By${BAy} -B${BAl} -Y${iy}

	    #Absolute bias difference/RMSD ratio
	    gmt psmask ${index}_${type}${idat}.20 -R${range2} -I${int}
	    gmt grdimage ${index}_${type}${idat}.grd -C${var}${index}_${type}.cpt
	    gmt psmask -C
	    
	    #Significant test
	    if(-f ${index}_sig${idat}.20)then
		gmt psxy ${index}_sig${idat}.20 -Sc0.05 -Gblack
	    endif

	    #Coast
	    gmt coast -R${range1} -Dl -W0.2,black -Gwhite

	    @ ix = ${slon}; @ iy = ${elat} + ${ypos}
	    echo "${ix} ${iy} ${label2[$idat]}" | gmt text -F+f14p,0,black+jLB -N	    

	    if(${idat} == ${ndat})then
		gmt colorbar -Dx${drange} -Bx${dBA} -C${var}${index}_${type}.cpt --FONT=20p
	    endif

	    @ idat++
	    
	end #idat
	
	gmt end
	
    end #index

    rm -f *.20 *.grd
    
end #var

rm -f *.cpt
