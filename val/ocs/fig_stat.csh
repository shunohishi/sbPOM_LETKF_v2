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
set label=("(a) RMSD & Spread" "Spread" "(b) Bias" "(c) RMSD ratio" "(d) Absolute bias difference")
set color=("black" "cyan" "orange" "green")
set legend=("LORA" "GLORYS" "ORAS5" "C-GLORS")
set buoyname=("KEO" "Papa")
set varname=("Temperature" "Salinity" "Zonal velocity" "Meridional velocity")

@ ibuoy = 1
@ nbuoy = 2
foreach buoy(keo papa)

@ ivar = 1
@ nvar = 4
foreach var(t s u v)

    set title="${varname[$ivar]} at ${buoyname[$ibuoy]} buoy"
    echo ${title}
    
    #---NOBS
    if(${var} == "t" || ${var} == "s")then
	@ ratio = 20 #[%]
    else if(${var} == "u" || ${var} == "v")then
	@ ratio = 20 #[%]
    endif	
        
    gmt begin fig/${buoy}-${var} png

	@ i = 1
	@ n = 3
	#---RMSD & Bias
	foreach index(rmsd sprd bias)

	    #---Data
	    set input=dat/${buoy}/${var}${index}_ave.dat
	    gawk -v ratio=${ratio} '{if(ratio < $2 && $6 != -999.) print $6,$1 > "dat1.20"}' ${input}
	    gawk -v ratio=${ratio} '{if(ratio < $3 && $7 != -999.) print $7,$1 > "dat2.20"}' ${input}
	    gawk -v ratio=${ratio} '{if(ratio < $4 && $8 != -999.) print $8,$1 > "dat3.20"}' ${input}
	    gawk -v ratio=${ratio} '{if(ratio < $5 && $9 != -999.) print $9,$1 > "dat4.20"}' ${input}
	    
	    #---Figure setting
	    #---KEO
	    if(${buoy} == "keo" && ${var} == "t" && ${index} == "rmsd")then
		set xs=0; set xe=3
		set ys=0; set ye=550
		set BAx=a1f0.2+l"RMSD\040\046\040Spread(\260C)"
		set BAy=a100f10+l"Depth\040(m)"
	    else if(${buoy} == "keo" && ${var} == "t" && ${index} == "bias")then
		set xs=-1.0; set xe=1.0
		set ys=0; set ye=550
		set BAx=a0.5f0.1g999+l"Bias\040(\260C)"
		set BAy=a100f10+l"Depth\040(m)"
	    else if(${buoy} == "keo" && ${var} == "s" && ${index} == "rmsd")then
		set xs=0; set xe=0.3
		set ys=0; set ye=550
		set BAx=a0.1f0.02+l"RMSD\040\046\040Spread"
		set BAy=a100f10+l"Depth\040(m)"
	    else if(${buoy} == "keo" && ${var} == "s" && ${index} == "bias")then
		set xs=-0.1; set xe=0.1
		set ys=0; set ye=550
		set BAx=a0.05f0.01g999+l"Bias"
		set BAy=a100f10+l"Depth\040(m)"
	    else if(${buoy} == "keo" && (${var} == "u" || ${var} == "v") && ${index} == "rmsd")then
		set xs=0; set xe=0.5
		set ys=0; set ye=40
		set BAx=a0.1f0.05+l"RMSD\040\046\040Spread(m/s)"
		set BAy=a10f5+l"Depth\040(m)"
	    else if(${buoy} == "keo" && (${var} == "u" || ${var} == "v") && ${index} == "bias")then
		set xs=-0.20; set xe=0.20
		set ys=0; set ye=40
		set BAx=a0.1f0.02g999+l"Bias\040(m/s)"
		set BAy=a10f5+l"Depth\040(m)"
	    #---Papa
	    else if(${buoy} == "papa" && ${var} == "t" && ${index} == "rmsd")then
		set xs=0; set xe=1
		set ys=0; set ye=350
		set BAx=a0.5f0.1+l"RMSD\040\046\040Spread(\260C)"
		set BAy=a50f10+l"Depth\040(m)"
	    else if(${buoy} == "papa" && ${var} == "t" && ${index} == "bias")then
		set xs=-0.3; set xe=0.3
		set ys=0; set ye=350
		set BAx=a0.15f0.05g999+l"Bias"
		set BAy=a50f10+l"Depth\040(m)"
	    else if(${buoy} == "papa" && ${var} == "s" && ${index} == "rmsd")then
		set xs=0; set xe=0.4
		set ys=0; set ye=350
		set BAx=a0.1f0.02+l"RMSD\040\046\040Spread"
		set BAy=a50f10+l"Depth\040(m)"
	    else if(${buoy} == "papa" && ${var} == "s" && ${index} == "bias")then
		set xs=-0.3; set xe=0.3
		set ys=0; set ye=350
		set BAx=a0.10f0.02g999+l"Bias"
		set BAy=a50f10+l"Depth\040(m)"
	    else if(${buoy} == "papa" && (${var} == "u" || ${var} == "v") && ${index} == "rmsd")then
		set xs=0; set xe=0.1
		set ys=0; set ye=40
		set BAx=a0.05f0.01+l"RMSD\040\046\040Spread(m/s)"
		set BAy=a10f5+l"Depth\040(m)"
	    else if(${buoy} == "papa" && (${var} == "u" || ${var} == "v") && ${index} == "bias")then
		set xs=-0.05; set xe=0.05
		set ys=0; set ye=40
		set BAx=a0.05f0.01g999+l"Bias\040(m/s)"
		set BAy=a10f5+l"Depth\040(m)"
	    endif
	    set range=${xs}/${xe}/${ys}/${ye}

	    #---Depth average data
	    if(${index} == "rmsd")then
		set input=dat/${buoy}/${var}${index}_dave.dat
		gawk -v ye=${ye} '{if($5*$9 > 0) print $1, ye > "sig1.20"; else print $1, ye > "nosig1.20"}' ${input}
		gawk -v ye=${ye} '{if($6*$10 > 0) print $2, ye > "sig2.20"; else print $2, ye > "nosig2.20"}' ${input}
		gawk -v ye=${ye} '{if($7*$11 > 0) print $3, ye > "sig3.20"; else print $3, ye > "nosig3.20"}' ${input}
		gawk -v ye=${ye} '{if($8*$12 > 0) print $4, ye > "sig4.20"; else print $4, ye > "nosig4.20"}' ${input}
	    endif
	    
	    #---Figure
	    if(${index} == "rmsd")then
		gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl}+t"${title}" -X2.5 -Y20
	    else if(${index} == "bias")then
		gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X8.5
	    endif

	    @ idat = 1
	    @ ndat = 4
	    while($idat <= $ndat)

		#Check data
		if(! -f dat${idat}.20)then
		    @ idat++
		    continue
		endif

		#Draw
		if(${index} == "rmsd")then
		    gmt psxy dat${idat}.20 -W2,${color[$idat]} -l${legend[$idat]}
		else if(${index} == "sprd")then 
		    gmt psxy dat${idat}.20 -W0.5,${color[$idat]},-	    
		else if(${index} == "bias")then
		    gmt psxy dat${idat}.20 -W2,${color[$idat]}
		endif
		gmt psxy dat${idat}.20 -Sc0.2 -G${color[$idat]}

		if(-f sig${idat}.20)then
		    gmt psxy sig${idat}.20 -Si0.4 -Gwhite -W2,${color[${idat}]} -N
		else if(-f nosig${idat}.20)then
		    gmt psxy nosig${idat}.20 -Si0.4 -G${color[${idat}]} -W2,${color[${idat}]} -N
		endif
		
		@ idat++
		
	    end #idat

	    if(${buoy} == "keo" && (${var} == "t" || ${var} == "s") && ${index} == "sprd")then
		gmt legend -DjRT+jRT+o0.2/0.2 -F+gwhite+pblack --FONT=8p,Helvetica,black
	    else if(${buoy} == "papa" && (${var} == "t" || ${var} == "s") && ${index} == "sprd")then
		gmt legend -DjRB+jRB+o0.2/0.2 -F+gwhite+pblack --FONT=8p,Helvetica,black
	    else if((${var} == "u" || ${var} == "v") && ${index} == "sprd")then
		gmt legend -DjLT+jLT+o0.2/0.2 -F+gwhite+pblack --FONT=8p,Helvetica,black
	    endif

	    if(${index} == "rmsd" || ${index} == "bias")then
		echo "${xs} 0 ${label[$i]}" | gmt text -F+f14p,0,black+jLB -Dj0.c/0.2c -N 
	    endif
		
	    @ i++
	    rm -f *.20
	    
	end #index

	#---RMSD ratio & Absolute Bias
	@ i = 4
	@ n = 5
	foreach index(rmsd bias)

	    #---Data
	    set input=dat/${buoy}/${var}${index}_ave.dat
	    if(${index} == "rmsd")then
		gawk -v ratio=${ratio} '{if(ratio < $3) print 100*($7-$6)/$6,$1 > "dat2.20"}' ${input}
		gawk -v ratio=${ratio} '{if(ratio < $4) print 100*($8-$6)/$6,$1 > "dat3.20"}' ${input}
		gawk -v ratio=${ratio} '{if(ratio < $5) print 100*($9-$6)/$6,$1 > "dat4.20"}' ${input}
		gawk -v ratio=${ratio} \
		'{if(ratio < $3 && $11*$15 > 0.) print 100*($7-$6)/$6,$1 > "sig2.20"}' ${input}
		gawk -v ratio=${ratio} \
		'{if(ratio < $4 && $12*$16 > 0.) print 100*($8-$6)/$6,$1 > "sig3.20"}' ${input}
		gawk -v ratio=${ratio} \
		'{if(ratio < $5 && $13*$17 > 0.) print 100*($9-$6)/$6,$1 > "sig4.20"}' ${input}
	    else if(${index} == "bias")then
		gawk -v ratio=${ratio} '{if(ratio < $3) print sqrt($7*$7)-sqrt($6*$6),$1 > "dat2.20"}' ${input}
		gawk -v ratio=${ratio} '{if(ratio < $4) print sqrt($8*$8)-sqrt($6*$6),$1 > "dat3.20"}' ${input}
		gawk -v ratio=${ratio} '{if(ratio < $5) print sqrt($9*$9)-sqrt($6*$6),$1 > "dat4.20"}' ${input}
		gawk -v ratio=${ratio} \
		'{if(ratio < $3 && $11*$15 > 0.) print sqrt($7*$7)-sqrt($6*$6),$1 > "sig2.20"}' ${input}
		gawk -v ratio=${ratio} \
		'{if(ratio < $4 && $12*$16 > 0.) print sqrt($8*$8)-sqrt($6*$6),$1 > "sig3.20"}' ${input}
		gawk -v ratio=${ratio} \
		'{if(ratio < $5 && $13*$17 > 0.) print sqrt($9*$9)-sqrt($6*$6),$1 > "sig4.20"}' ${input}
	    endif	    
	    
	    #---Figure setting
	    #---KEO
	    if(${buoy} == "keo" && ${var} == "t" && ${index} == "rmsd")then
		set xs=-50; set xe=50
		set ys=0; set ye=550
		set BAx=a25f5g999+l"RMSD\040ratio\040(\045)"
		set BAy=a100f10+l"Depth\040(m)"
	    else if(${buoy} == "keo" && ${var} == "t" && ${index} == "bias")then
		set xs=-1.0; set xe=1.0
		set ys=0; set ye=550
		set BAx=a0.5f0.1g999+l"Absolute\040bias\040difference\040(\260C)"
		set BAy=a100f10+l"Depth\040(m)"
	    else if(${buoy} == "keo" && ${var} == "s" && ${index} == "rmsd")then
		set xs=-50; set xe=50
		set ys=0; set ye=550
		set BAx=a25f5g999+l"RMSD\040ratio\040(\045)"
		set BAy=a100f10+l"Depth\040(m)"
	    else if(${buoy} == "keo" && ${var} == "s" && ${index} == "bias")then
		set xs=-0.1; set xe=0.1
		set ys=0; set ye=550
		set BAx=a0.05f0.01g999+l"Absolute\040bias\040difference"
		set BAy=a100f10+l"Depth\040(m)"
	    else if(${buoy} == "keo" && (${var} == "u" || ${var} == "v") && ${index} == "rmsd")then
		set xs=-50; set xe=50
		set ys=0; set ye=40
		set BAx=a25f5g999+l"RMSD\040ratio\040(\045)"
		set BAy=a10f5+l"Depth\040(m)"
	    else if(${buoy} == "keo" && (${var} == "u" || ${var} == "v") && ${index} == "bias")then
		set xs=-0.20; set xe=0.20
		set ys=0; set ye=40
		set BAx=a0.1f0.02g999+l"Absolute\040bias\040difference\040(m/s)"
		set BAy=a10f5+l"Depth\040(m)"
	    #---Papa
	    else if(${buoy} == "papa" && ${var} == "t" && ${index} == "rmsd")then
		set xs=-100; set xe=100
		set ys=0; set ye=350
		set BAx=a50f10g999+l"RMSD\040ratio\040(\045)"
		set BAy=a50f10+l"Depth\040(m)"
	    else if(${buoy} == "papa" && ${var} == "t" && ${index} == "bias")then
		set xs=-0.4; set xe=0.4
		set ys=0; set ye=350
		set BAx=a0.2f0.1g999+l"Absolute\040bias\040difference\040(\260C)"
		set BAy=a50f10+l"Depth\040(m)"
	    else if(${buoy} == "papa" && ${var} == "s" && ${index} == "rmsd")then
		set xs=-100; set xe=100
		set ys=0; set ye=350
		set BAx=a50f10g999+l"RMSD\040ratio\040(\045)"
		set BAy=a50f10+l"Depth\040(m)"
	    else if(${buoy} == "papa" && ${var} == "s" && ${index} == "bias")then
		set xs=-0.2; set xe=0.2
		set ys=0; set ye=350
		set BAx=a0.1f0.02g999+l"Absolute\040bias\040difference"
		set BAy=a50f10+l"Depth\040(m)"
	    else if(${buoy} == "papa" && (${var} == "u" || ${var} == "v") && ${index} == "rmsd")then
		set xs=-50; set xe=50
		set ys=0; set ye=40
		set BAx=a25f5g999+l"RMSD\040ratio\040(\045)"
		set BAy=a10f5+l"Depth\040(m)"
	    else if(${buoy} == "papa" && (${var} == "u" || ${var} == "v") && ${index} == "bias")then
		set xs=-0.05; set xe=0.05
		set ys=0; set ye=40
		set BAx=a0.05f0.01g999+l"Absolute\040bias\040difference\040(m/s)"
		set BAy=a10f5+l"Depth\040(m)"
	    endif
	    set range=${xs}/${xe}/${ys}/${ye}

	    #---Figure
	    if($i == 4)then
		gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X-8.5 -Y-8.5 
	    else
		gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X8.5
	    endif

	    @ idat = 2
	    @ ndat = 4
	    while($idat <= $ndat)

		gmt psxy dat${idat}.20 -W2,${color[$idat]}
		gmt psxy dat${idat}.20 -Sc0.2 -G${color[$idat]}

		if(-f sig${idat}.20)then
		    gmt psxy sig${idat}.20 -Sc0.2 -Gwhite -W2,${color[$idat]}
		endif
		
		@ idat++
		
	    end #idat

	    echo "${xs} 0 ${label[$i]}" | gmt text -F+f14p,0,black+jLB -Dj0.c/0.2c -N 
	    
	    @ i++
	    rm -f *.20
	    
	end #index
	
    gmt end

@ ivar++
    
end

@ ibuoy++

end
