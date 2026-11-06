#!/bin/csh
#===============================================#
# Machine
#===============================================#

#set machine="fugaku"
#set machine="jss3"
set machine="rc"

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

#---Data size
set ndat=5 #Data size

#---Criteria of observation availability (%)
set tcrit=20
set scrit=20
set ucrit=20
set vcrit=20

#---Box size
set size=6/-6 #Size

#---Axis
set BAl=WSne  #Axis label

#---Label
set label=("(a) RMSD & Spread" "(b) Bias" "(c) RMSD ratio" "(d) Absolute bias difference")

#---Line color
set color=("black" "#0072B2" "#E69F00" "#009E73" "#D55E00")

#---Legend
set legend=("LORA-NP" "BRAN2020" "FORA-JPN60" "GLORYS12V1" "JCOPE-FGO")
set legend_position_keo_t="RT" #LRTB (L:Left, R:Right, T:Top, B:Bottom)
set legend_position_keo_s="LB"
set legend_position_keo_u="LT"
set legend_position_keo_v="LT"
set legend_position_papa_t="RB"
set legend_position_papa_s="RB"
set legend_position_papa_u="LT"
set legend_position_papa_v="LT"

#---Buoy name
set buoyname=("KEO" "Papa")

#---Variable name
set varname=("Temperature" "Salinity" "Zonal velocity" "Meridional velocity")

#---Horizontal & Vertical axis
#---KEO
#Y
set ys_keo_t=0; set ye_keo_t=550      #Y range
set ys_keo_s=0; set ye_keo_s=550
set ys_keo_u=0; set ye_keo_u=40
set ys_keo_v=0; set ye_keo_v=40

set BAy_keo_t=a100f10+l"Depth\040(m)" #Y label
set BAy_keo_s=a100f10+l"Depth\040(m)"
set BAy_keo_u=a10f5+l"Depth\040(m)"
set BAy_keo_v=a10f5+l"Depth\040(m)"

#KEO/T/RMSD
set xs_keo_trmsd=0; set xe_keo_trmsd=3                    #X range
set BAx_keo_trmsd=a1f0.2+l"RMSD\040\046\040Spread\040(\260C)" #X label

#KEO/T/Bias
set xs_keo_tbias=-1.5; set xe_keo_tbias=1.5
set BAx_keo_tbias=a0.5f0.1g999+l"Bias\040(\260C)"

#KEO/T/Ratio
set xs_keo_tratio=-75; set xe_keo_tratio=75
set BAx_keo_tratio=a25f5g999+l"RMSD\040ratio\040(\045)"

#KEO/T/Abias
set xs_keo_tabias=-1.5; set xe_keo_tabias=1.5
set BAx_keo_tabias=a0.5f0.1g999+l"Absolute\040bias\040difference\040(\260C)"

#KEO/S/RMSD
set xs_keo_srmsd=0; set xe_keo_srmsd=0.2
set BAx_keo_srmsd=a0.1f0.02+l"RMSD\040\046\040Spread"

#KEO/S/Bias
set xs_keo_sbias=-0.1; set xe_keo_sbias=0.1
set BAx_keo_sbias=a0.05f0.01g999+l"Bias"

#KEO/S/Ratio
set xs_keo_sratio=-75; set xe_keo_sratio=75
set BAx_keo_sratio=a25f5g999+l"RMSD\040ratio\040(\045)"

#KEO/S/Abias
set xs_keo_sabias=-0.1; set xe_keo_sabias=0.1
set BAx_keo_sabias=a0.05f0.01g999+l"Absolute\040bias\040difference"

#KEO/U/RMSD
set xs_keo_urmsd=0; set xe_keo_urmsd=0.4
set BAx_keo_urmsd=a0.1f0.05+l"RMSD\040\046\040Spread\040(m/s)"

#KEO/U/Bias
set xs_keo_ubias=-0.30; set xe_keo_ubias=0.30
set BAx_keo_ubias=a0.1f0.02g999+l"Bias\040(m/s)"

#KEO/U/Ratio
set xs_keo_uratio=-75; set xe_keo_uratio=75
set BAx_keo_uratio=a25f5g999+l"RMSD\040ratio\040(\045)"

#KEO/U/Abias
set xs_keo_uabias=-0.30; set xe_keo_uabias=0.30
set BAx_keo_uabias=a0.1f0.02g999+l"Absolute\040bias\040difference\040(m/s)"

#KEO/V/RMSD
set xs_keo_vrmsd=0; set xe_keo_vrmsd=0.4
set BAx_keo_vrmsd=a0.1f0.05+l"RMSD\040\046\040Spread\040(m/s)"

#KEO/V/Bias
set xs_keo_vbias=-0.30; set xe_keo_vbias=0.30
set BAx_keo_vbias=a0.1f0.02g999+l"Bias\040(m/s)"

#KEO/V/Ratio
set xs_keo_vratio=-75; set xe_keo_vratio=75
set BAx_keo_vratio=a25f5g999+l"RMSD\040ratio\040(\045)"

#KEO/V/Abias
set xs_keo_vabias=-0.30; set xe_keo_vabias=0.30
set BAx_keo_vabias=a0.1f0.02g999+l"Absolute\040bias\040difference\040(m/s)"

#---Papa
#Y
set ys_papa_t=0; set ye_papa_t=350
set ys_papa_s=0; set ye_papa_s=350
set ys_papa_u=0; set ye_papa_u=40
set ys_papa_v=0; set ye_papa_v=40

set BAy_papa_t=a50f10+l"Depth\040(m)"
set BAy_papa_s=a50f10+l"Depth\040(m)"
set BAy_papa_u=a10f5+l"Depth\040(m)"
set BAy_papa_v=a10f5+l"Depth\040(m)"

#Papa/T/RMSD
set xs_papa_trmsd=0; set xe_papa_trmsd=1
set BAx_papa_trmsd=a0.5f0.1+l"RMSD\040\046\040Spread\040(\260C)"

#Papa/T/Bias
set xs_papa_tbias=-0.6; set xe_papa_tbias=0.6
set BAx_papa_tbias=a0.20f0.05g999+l"Bias\040(\260C)"

#Papa/T/Ratio
set xs_papa_tratio=-100; set xe_papa_tratio=100
set BAx_papa_tratio=a50f10g999+l"RMSD\040ratio\040(\045)"

#Papa/T/Abias
set xs_papa_tabias=-0.6; set xe_papa_tabias=0.6
set BAx_papa_tabias=a0.20f0.05g999+l"Absolute\040bias\040difference\040(\260C)"

#Papa/S/RMSD
set xs_papa_srmsd=0; set xe_papa_srmsd=0.4
set BAx_papa_srmsd=a0.1f0.02+l"RMSD\040\046\040Spread"

#Papa/S/Bias
set xs_papa_sbias=-0.2; set xe_papa_sbias=0.2
set BAx_papa_sbias=a0.10f0.02g999+l"Bias"

#Papa/S/Ratio
set xs_papa_sratio=-100; set xe_papa_sratio=100
set BAx_papa_sratio=a50f10g999+l"RMSD\040ratio\040(\045)"

#Papa/S/Abias
set xs_papa_sabias=-0.2; set xe_papa_sabias=0.2
set BAx_papa_sabias=a0.1f0.02g999+l"Absolute\040bias\040difference"

#Papa/U/RMSD
set xs_papa_urmsd=0; set xe_papa_urmsd=0.1
set BAx_papa_urmsd=a0.05f0.01+l"RMSD\040\046\040Spread\040(m/s)"

#Papa/U/Bias
set xs_papa_ubias=-0.02; set xe_papa_ubias=0.02
set BAx_papa_ubias=a0.01f0.002g999+l"Bias\040(m/s)"

#Papa/U/Ratio
set xs_papa_uratio=-75; set xe_papa_uratio=75
set BAx_papa_uratio=a25f5g999+l"RMSD\040ratio\040(\045)"

#Papa/U/Abias
set xs_papa_uabias=-0.05; set xe_papa_uabias=0.05
set BAx_papa_uabias=a0.05f0.01g999+l"Absolute\040bias\040difference\040(m/s)"

#Papa/V/RMSD
set xs_papa_vrmsd=0; set xe_papa_vrmsd=0.1
set BAx_papa_vrmsd=a0.05f0.01+l"RMSD\040\046\040Spread\040(m/s)"

#Papa/V/Bias
set xs_papa_vbias=-0.02; set xe_papa_vbias=0.02
set BAx_papa_vbias=a0.01f0.002g999+l"Bias\040(m/s)"

#Papa/V/Ratio
set xs_papa_vratio=-75; set xe_papa_vratio=75
set BAx_papa_vratio=a25f5g999+l"RMSD\040ratio\040(\045)"

#Papa/V/Abias
set xs_papa_vabias=-0.02; set xe_papa_vabias=0.02
set BAx_papa_vabias=a0.01f0.002g999+l"Absolute\040bias\040difference\040(m/s)"

#========================================================

#===Data

echo "=== Start: Extract necessary data ==="

@ ibuoy = 1
@ nbuoy = 2
foreach buoy(keo papa)

    @ ivar = 1
    @ nvar = 4
    foreach var(t s u v)
    
	#---Observation availability criteria
	eval set crit = \$${var}crit

	#---Y axis
	eval set ye = \$ye_${buoy}_${var}
    
	#===Data
	#---RMSD/Spread/Bias
	foreach index(rmsd sprd bias)
	
	    set input=dat/${buoy}/${var}${index}_ave.dat
	    if(! -f ${input})then
		echo "***Error: Not found ${input}"
		exit
	    endif

	    @ idat = 1
	    while(${idat} <= ${ndat})

	    @ obs_col = 1 + ${idat}
	    @ dat_col = 1 + ${ndat} + ${idat}
	    gawk -v crit=${crit} -v obs_col=${obs_col} -v dat_col=${dat_col} -v out=${buoy}_${var}${index}${idat}.20 \
	    '{if(crit < $obs_col && $dat_col != -999.) print $dat_col,$1 > out}' ${input}

	    @ idat++

	    end #idat   
	end     #index

	#---Depth avraged RMSD
	set input=dat/${buoy}/${var}rmsd_dave.dat
	if(! -f ${input})then
	    echo "***Error: Not found ${input}"
	    exit
	endif
		
	@ idat = 1
	while(${idat} <= ${ndat})

	    @ dat_col = ${idat}
	    @ dif1_col = ${ndat} + ${idat}
	    @ dif2_col = 2 * ${ndat} + ${idat}
		    
	    gawk -v dat_col=${dat_col} -v dif1_col=${dif1_col} -v dif2_col=${dif2_col} -v ye=${ye} \
	    -v sig=${buoy}_${var}rmsd_sig${idat}.20 -v nosig=${buoy}_${var}rmsd_nosig${idat}.20 \
	    '{if($dif1_col * $dif2_col > 0) print $dat_col, ye > sig; else print $dat_col, ye > nosig}' ${input}

	    @ idat++

	end #idat
	
	#---RMSD ratio
	set input=dat/${buoy}/${var}rmsd_ave.dat
	if(! -f ${input})then
	    echo "***Error: Not found ${input}"
	    exit
	endif
	
        @ idat = 2
	while(${idat} <= ${ndat})

	    @ obs_col = 1 + 1
	    @ dat1_col = 1 + ${ndat} + 1
	    @ dat2_col = 1 + ${ndat} + ${idat}
	    @ dif1_col = 1 + 2 * ${ndat} + ${idat}
	    @ dif2_col = 1 + 3 * ${ndat} + ${idat}
	    
	    gawk -v crit=${crit} -v obs_col=${obs_col} -v dat1_col=${dat1_col} -v dat2_col=${dat2_col} \
	    -v out=${buoy}_${var}ratio${idat}.20 \
	    '{if(crit < $obs_col && $dat1_col != 0 && $dat1_col != -999 && $dat2_col != -999) print 100*($dat2_col-$dat1_col)/$dat1_col,$1 > out}' ${input}

	    gawk -v crit=${crit} -v obs_col=${obs_col} -v dat1_col=${dat1_col} -v dat2_col=${dat2_col} \
	    -v dif1_col=${dif1_col} -v dif2_col=${dif2_col} -v out=${buoy}_${var}ratio_sig${idat}.20 \
	    '{if(crit < $obs_col && $dif1_col * $dif2_col > 0. && $dat1_col != 0 && $dat1_col != -999 && $dat2_col != -999) print 100*($dat2_col-$dat1_col)/$dat1_col,$1 > out}' ${input}
	    
	    @ idat++
	    
	end #idat
	
	#---Absolute bias
	set input=dat/${buoy}/${var}bias_ave.dat
	if(! -f ${input})then
	    echo "***Error: Not found ${input}"
	    exit
	endif

	@ idat = 2
	while(${idat} <= ${ndat})

	    @ obs_col = 1 + 1
	    @ dat1_col = 1 + ${ndat} + 1
	    @ dat2_col = 1 + ${ndat} + ${idat}
	    @ dif1_col = 1 + 2 * ${ndat} + ${idat}
	    @ dif2_col = 1 + 3 * ${ndat} + ${idat}
	
	    gawk -v crit=${crit} -v obs_col=${obs_col} -v dat1_col=${dat1_col} -v dat2_col=${dat2_col} \
	    -v out=${buoy}_${var}abias${idat}.20 '{if(crit < $obs_col && $dat1_col != -999 && $dat2_col != -999) print sqrt($dat2_col*$dat2_col)-sqrt($dat1_col*$dat1_col),$1 > out}' ${input}

	    gawk -v crit=${crit} -v obs_col=${obs_col} -v dat1_col=${dat1_col} -v dat2_col=${dat2_col} \
	    -v dif1_col=${dif1_col} -v dif2_col=${dif2_col} -v out=${buoy}_${var}abias_sig${idat}.20 \
	    '{if(crit < $obs_col && $dif1_col*$dif2_col > 0. && $dat1_col != -999 && $dat2_col != -999) print sqrt($dat2_col*$dat2_col)-sqrt($dat1_col*$dat1_col),$1 > out}' ${input}

	    @ idat++

	end #idat
    end #var
end #buoy

echo "=== End: Extract necessary data ==="

#===Figure

echo "=== Start: Make figure ==="

@ ibuoy = 1
@ nbuoy = 2
foreach buoy(keo papa)

    @ ivar = 1
    @ nvar = 4
    foreach var(t s u v)

	gmt begin fig/${buoy}-${var} png

	#---Title
	set title="${varname[$ivar]} at ${buoyname[$ibuoy]} buoy"
	echo "Buoy [${ibuoy}/${nbuoy}] Variable [${ivar}/${nvar}] ${title}"

	#---Y axis
	eval set ys = \$ys_${buoy}_${var}
	eval set ye = \$ye_${buoy}_${var}
	eval set BAy = \$BAy_${buoy}_${var}
	
	#---Legend position
	eval set legend_position = \$legend_position_${buoy}_${var}
    
	@ ifig = 1
	@ nfig = 4
	foreach index(rmsd bias ratio abias)

	    #---Figure setting
	    eval set xs = \$xs_${buoy}_${var}${index}
	    eval set xe = \$xe_${buoy}_${var}${index}	    
	    set range=${xs}/${xe}/${ys}/${ye}
	    eval set BAx = \$BAx_${buoy}_${var}${index}
		
	    #---Box
	    if(${ifig} == 1)then
		gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl}+t"${title}" -X2.5 -Y20
	    else if(${ifig} % 2 == 0)then
		gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X8.5
	    else
		gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X-8.5 -Y-8.5	    
	    endif

	    #---Label
	    echo "${xs} 0 ${label[$ifig]}" | gmt text -F+f14p,0,black+jLB -Dj0.c/0.2c -N 
	    
	    #sdat
	    if(${index} == "rmsd" || ${index} == "bias")then
		@ sdat = 1
	    else if(${index} == "ratio" || ${index} == "abias")then
		@ sdat = 2
	    endif

	    #---RMSD/Bias/RMSD ratio/Absolute bias
	    @ idat = ${sdat}
	    while(${idat} <= ${ndat})
	
		set file="${buoy}_${var}${index}${idat}.20"		
		#Line
		if(${index} == "rmsd" && -f ${file})then
		    gmt psxy ${file} -W2,${color[$idat]} -l${legend[$idat]}
		else if(-f ${file})then
		    gmt psxy ${file} -W2,${color[$idat]}
		endif
		
		#Closed circle
		if(-f ${file})then
		    gmt psxy ${file} -Sc0.2 -G${color[$idat]} && rm -f ${file}
		endif
					
		#Open circle for Siginficant differenct
		set file="${buoy}_${var}${index}_sig${idat}.20"
		if((${index} == "ratio" || ${index} == "abias") && -f ${file})then
		    gmt psxy ${file} -Sc0.2 -Gwhite -W2,${color[$idat]} && rm -f ${file}
		endif		
		
		@ idat++
		
	    end #idat

	    #---Spread, Depth-averaged RMSD, legen for RMSD only
	    if(${index} == "rmsd")then

		#---Spread
		@ idat = ${sdat}
		while(${idat} <= ${ndat})

		    set file="${buoy}_${var}sprd${idat}.20"
		    if(-f ${file}) gmt psxy ${file} -W0.5,${color[${idat}]},4_1
		    if(-f ${file}) gmt psxy ${file} -Sc0.2 -G${color[$idat]} && rm -f ${file}
		    		    
		    @ idat++
		    
		end #idat

		#---Depth-averaged RMSD
		@ idat = ${sdat}
		while(${idat} <= ${ndat})
		
		    set file="${buoy}_${var}${index}_sig${idat}.20"
		    if(-f ${file}) gmt psxy ${file} -Si0.4 -Gwhite -W2,${color[${idat}]} -N && rm -f ${file}
			
		    set file="${buoy}_${var}${index}_nosig${idat}.20"
		    if(-f ${file}) gmt psxy ${file} -Si0.4 -G${color[${idat}]} -W2,${color[${idat}]} -N && rm -f ${file}

		    @ idat++
		    
		end #idat

		#---Legend
		gmt legend -Dj${legend_position}+j${legend_position}+o0.2/0.2 -F+gwhite+pblack --FONT=6p,Helvetica,black
		
	    endif	    

	    @ ifig++
	    
	end #index
            
	gmt end

	@ ivar++
	
    end #var

    @ ibuoy++

end #buoy

echo "=== End: Make figure ==="
