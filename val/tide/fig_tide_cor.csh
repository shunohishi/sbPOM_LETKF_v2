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

gmt set FONT=10p,Helvetica,black
gmt set MAP_FRAME_TYPE=PLAIN FORMAT_GEO_MAP=dddF
gmt set GMT_AUTO_DOWNLOAD=off

if(${machine} == "jss3")then
    gmt set PS_CONVERT=I+m0.6/0.6/0.6/0.6 #WESN
endif

#=======================================================
# Make directory
#=======================================================

if(! -d fig) mkdir fig

#=======================================================
# Color
#=======================================================

gmt makecpt -T-0.5/0.5/0.1 -Cvik -D > color.cpt
set drange=8.5/0.5+w3/0.25+e0.5
set dBA=a0.5f0.1+l"Correlation\040(RMSD\040vs.\040Spread)"

#=======================================================
# DATA
#=======================================================

@ nobs=365
@ best=746
@ worst=292

#---Correlation
set input=dat/cor_ave.dat
gawk -v nobs=${nobs} '{if(nobs < $4 && $8 != -999 && $8 != 1. && $8 <= 0.) print $2,$3,$8 > "n.20"}' ${input}
gawk -v nobs=${nobs} '{if(nobs < $4 && $8 != -999 && $8 != -1. && 0. < $8) print $2,$3,$8 > "p.20"}' ${input}

#---Station location
gawk -v i=${best} '{if($1 == i) print $2,$3 > "best.20"}' ${input}
gawk -v i=${worst} '{if($1 == i) print $2,$3 > "worst.20"}' ${input}

#=======================================================
# Figure
#=======================================================

gmt begin fig/tide_cor png

#=======================================================
# Figure (Map)
#=======================================================

set size=8d/4d
set range=0/360/-90/90
set BAx=a60f10
set BAy=a30f10
set BAl=WSne

gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X3 -Y20

gmt coast -Dl -W0.2,black -Ggray
gmt psxy n.20 -Sc0.2 -W0.2 -Ccolor.cpt
gmt psxy p.20 -Sc0.2 -W0.2 -Ccolor.cpt

gmt psxy best.20 -Sc0.20 -W2,yellow
gmt psxy worst.20 -Sc0.20 -W2,cyan

#Ave
#set input=dat/cor_ave_all.dat
#set ave=`gawk '{printf "%.4f", $5}' ${input}`
#set ave="Correlation (station-mean): ${ave}"
#echo "360 -90 ${ave}" | gmt text -F+f12p,black+jRB -N 

gmt colorbar -Dx${drange} -Bx${dBA} -Ccolor.cpt --FONT=20p

gmt text -F+f14p,0,black+jLB -N <<EOF
0 95 (a) Correlation (RMSD vs. Spread)
EOF


#=======================================================
# Figure (Map)
#=======================================================

set size=8T/4
set range=2003-01-01/2023-12-31/0/0.25
set BApx=f1Y
set BAsx=a5Y+l"Date"
set BAy=a0.1f0.02g99+l"RMSD\046Spread"
set BAl=WSne
set label=("(b) Max. correlation (= 0.75): Cape May" "(c) Min. correlation (= -0.57): St. Helena")
set legend=("RMSD" "Spread")
set color=("black" "cyan")

@ i = 1

foreach ist(${best} ${worst})

    echo "Station $ist"

    #---Extract data
    set nnn=`printf "%03d" ${ist}`
    
    #RMSD
    set input=dat/rmsd_mave.dat
    gawk -v ist=$ist '{if($1 == ist && $9 != -999) print $4,$9 > "rmsd.20"}' ${input}
    gawk -v ist=$ist '$1 == ist {if($9 == -999){inblock=0}else{if(inblock==0){n++;outfile=sprintf("rmsd_%02d.20",n);inblock=1} print $4,$9 > outfile}}' ${input}

    #Spread
    set input=dat/sprd_mave.dat
    gawk -v ist=$ist '{if($1 == ist && $9 != -999) print $4,$9 > "sprd.20"}' ${input}
    gawk -v ist=$ist '$1 == ist {if($9 == -999){inblock=0}else{if(inblock==0){n++;outfile=sprintf("sprd_%02d.20",n);inblock=1} print $4,$9 > outfile}}' ${input}
    
    if($i == 1)then
	gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -Y-6
    else
	gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -X10
    endif

    gmt text -F+f14p,0,black+jLB -N <<EOF
2003-01-01 0.26 ${label[$i]}
EOF
    
    #---For legend
    if($i == 1)then
	@ idat = 1
	@ ndat = 2
	while($idat <= $ndat)

	if($idat == 1)then
	    set dat="rmsd"
	else if($idat == 2)then
	    set dat="sprd"
	endif
	
	if(! -f ${dat}.20)then
	    @ idat++
	    continue
	endif

	gmt psxy ${dat}.20 -W2,${color[$idat]} -l${legend[$idat]}
	
	@ idat++

	end

	gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl}+gwhite
	
    endif

    #---Main plot
    foreach idat(1 2)

	if($idat == 1)then
	    set dat="rmsd"
	else if($idat == 2)then
	    set dat="sprd"
	endif
    
	@ ifile = 1
	@ nfile = 99
	while($ifile <= $nfile)
	    
	    set nn=`printf "%02d" ${ifile}`
	    if(-f ${dat}_${nn}.20)then
		gmt psxy ${dat}_${nn}.20 -W2,${color[$idat]}
	    endif
		
	    @ ifile++
			
	end
    
    end
    
    if($i == 2)then
	gmt legend -DjRT+jRT+o0.2/0.2 -F+gwhite+pblack --FONT=10p
    endif
	
    @ i++

    rm -f *.20
    
end

gmt end

rm -f *.20 *.cpt
