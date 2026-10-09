#!/bin/csh

if(! -d fig) mkdir fig

#===============================================#
# Machine
#===============================================#

#set machine="fugaku" #Fugaku
#set machine="jss3"    #JSS3
set machine="rc"      #RCCS Cloud

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
gmt set FORMAT_DATE_MAP=o FORMAT_TIME_PRIMARY_MAP=c FORMAT_TIME_SECONDARY_MAP=f
gmt set GMT_AUTO_DOWNLOAD=off

if(${machine} == "jss3")then
    gmt set PS_CONVERT=I+m0.6/0.6/0.6/0.6 #WESN
endif
    
#=======================================================
# Figure setting
#=======================================================

set ndat=5                #Data size

#---Box size
set xsize=10; set ysize=5 #Figure size

#---Horizontal axis
set start_date=2003-01-01          #Start date
set end_date=2021-01-01            #End date
set BApx=f1Y; set BAsx=a5Y+l"Year" #Axis option

#---Vertical axis
#U RMSD & Spread
set ys_urmsd=0.0; set ye_urmsd=0.3                                  #Range
set BAy_urmsd=a0.1f0.02+l"RMSD\040(m/s)\040\046\040Spread\040(m/s)" #Label

#V RMSD & Spread
set ys_vrmsd=0.0; set ye_vrmsd=0.3
set BAy_vrmsd=a0.1f0.02+l"RMSD\040(m/s)\040\046\040Spread\040(m/s)"

#T RMSD & Spread
set ys_trmsd=0.0; set ye_trmsd=1.0
set BAy_trmsd=a0.5f0.1+l"RMSD\040(\260C)\040\046\040Spread\040(\260C)"

#U Bias
set ys_ubias=-0.04; set ye_ubias=0.04
set BAy_ubias=a0.02f0.005g99+l"Bias\040(m/s)"

#V Bias
set ys_vbias=-0.04; set ye_vbias=0.04
set BAy_vbias=a0.02f0.005g99+l"Bias\040(m/s)"

#T Bias
set ys_tbias=-0.2; set ye_tbias=0.2
set BAy_tbias=a0.1f0.02g99+l"Bias\040(\260C)"	    

#Obs availability
set ys_obsa=0; set ye_obsa=15000
set BAy_obsa=a5000f1000+l"Observation\040availability\040(month@+-1@+)"

#---Axis
set BAl=WSne #Axis option

#---Line Color
set color=("black" "#0072B2" "#E69F00" "#009E73" "#D55E00")

#---Dataset legend
set legend1=("LORA-NP" "BRAN2020" "FORA-JPN60" "GLORYS12V1" "JCOPE-FGO")

#=======================================================
# Figure 
#=======================================================

set size=${xsize}/${ysize}
set label=("(a) Surface zonal velocity" "(b) Surface meridional velocity" "(c) Sea-surface temperature" "(d) Surface drifter buoy")
set legend2=("U" "V" "T")

foreach stat(rmsd bias)

gmt begin fig/${stat}_time_series png

    @ ifig = 1
    @ nfig = 3
    
    foreach var(u v t)

	echo ${stat} ${var} 
    
	#---DATA
	@ idat = 1
	while(${idat} <= ${ndat})

	    #Obs availability
	    set input=dat/${var}${stat}_mave.dat		
	    @ obs_col = 1 + ${idat}
	    gawk -v obs_col=${obs_col} -v out=${var}obsa${idat}.20 '{if($obs_col != -999) print $1,$obs_col > out}' ${input}

	    #Bias/RMSD
	    @ dat_col = 1 + ${ndat} + ${idat}
	    gawk -v dat_col=${dat_col} -v out=${stat}${idat}.20 '{if($dat_col != -999) print $1,$dat_col > out}' ${input}

	    #Spread
	    set input=dat/${var}sprd_mave.dat
	    gawk -v dat_col=${dat_col} -v out=sprd${idat}.20 '{if($dat_col != -999) print $1,$dat_col > out}' ${input}
	    
	    @ idat++
	    
	end
		
	#---Axis & Label
	eval set ys = \$ys_${var}${stat}
	eval set ye = \$ye_${var}${stat}
	eval set BAy = \$BAy_${var}${stat}	
	set range=${start_date}/${end_date}/${ys}/${ye}

	#---Draw figure
	#Box
	set ix = `gawk -v xsize=$xsize 'BEGIN{print xsize + 3}'`
	set iy = `gawk -v ysize=$ysize 'BEGIN{print ysize + 2.5}'`	
	if(${ifig} == 1)then
	    gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -X3 -Y23
	else if(${ifig} % 2 == 0)then
	    gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -X${ix}
	else
	    gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -X-${ix} -Y-${iy}
	endif

	#Bias/RMSD/Spread
	@ idat=1
	while(${idat} <= ${ndat})

	    #Bias/RMSD
	    set dat="${stat}${idat}.20"
	    			
	    if(${ifig} == 1 && -f ${dat})then
		gmt psxy ${dat} -W2,${color[$idat]} -l${legend1[$idat]}
	    else if(-f ${dat})then
		gmt psxy ${dat} -W2,${color[$idat]}
	    endif

	    #Spread
	    set dat="sprd${idat}.20"
	    if(${stat} == "rmsd" && -f ${dat})then
		gmt psxy ${dat} -W0.5,${color[$idat]},-
	    endif	    
	    
	    @ idat++
	    
	end #idat

	#Legend
	if(${ifig} == 1)then
	    gmt legend -DjRT+jRT+o0.2/0.2 -F+gwhite+pblack --FONT=6p
	endif

	#Label
	echo "${start_date} ${ye} ${label[$ifig]}" | gmt text -F+f14p,0,black+jLT -Dj0.2c/0.2c -N
	
	@ ifig++

	rm -f ${stat}*.20 sprd*.20
	
    end #var

    #---Observation availability
    set ys = ${ys_obsa}
    set ye = ${ye_obsa}
    set BAy = ${BAy_obsa}
    
    set range=${start_date}/${end_date}/${ys}/${ye}

    set ix = `gawk -v xsize=$xsize 'BEGIN{print xsize + 3}'`
    gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -X${ix}

    @ ivar=1
    @ nvar=3
    foreach var(u v t)
	if(-f ${var}obsa1.20)then
	    gmt psxy ${var}obsa1.20 -W2,${color[$ivar]} -l${legend2[$ivar]}
	endif
	@ ivar++
    end

    gmt legend -DjRB+jRB+o0.2/0.2 -F+gwhite+pblack --FONT=10p

    echo "${start_date} ${ye_obsa} ${label[$ifig]}" | gmt text -F+f14p,0,black+jLT -Dj0.2c/0.2c -N
    
gmt end

rm -f *.20

end #stat
