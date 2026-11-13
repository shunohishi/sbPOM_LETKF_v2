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
gmt set FORMAT_DATE_MAP=o FORMAT_TIME_PRIMARY_MAP=c FORMAT_TIME_SECONDARY_MAP=f
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

#---Box
# Size
set xsize=12 #X Size
set ysize=6  #Y Size

#X axis
set start_date="2003-01-01"; set end_date="2020-12-31" #Start/End date
set BApx=f1Y; set BAsx=a5Y+l"Date"                     #X Axis

#Y axis
set ys_1=-0.4; set ye_1=0.4                  #Start/End point of Sea level [m]
set BAy_1=a0.2f0.1g99+l"Sea\040level\040(m)" #Y Axis

set ys_2=-0.4; set ye_2=0.4                  #Start/End point of Sea level [m]
set BAy_2=a0.2f0.1g99+l"Sea\040level\040(m)" #Y Axis

set ys_3=-0.4; set ye_3=0.4                  #Start/End point of Sea level [m]
set BAy_3=a0.2f0.1g99+l"Sea\040level\040(m)" #Y Axis

set BAl=WSen

#---Label
set name_ave="Crescent City"
set name_best="Kaumalapau"
set name_worst="Minamitorishima"
set label=("(a) ${name_ave} (closest to mean RMSD) " "(b) ${name_best} (smallest RMSD) " "(c) ${name_worst} (largest RMSD)")

#---Sation ID
set id_ave=556  #Crescent City
set id_best=548 #Kaumalapau
set id_worst=49 #Minamitorishima

#---Legend
set ndat=5 #Number of data 
set legend=("LORA-NP" "BRAN2020" "FORA-JPN60" "GLORYS12V1" "JCOPE-FGO")
set color=("black" "#0072B2" "#E69F00" "#009E73" "#D55E00")

#==========================================================
# Figure |
#==========================================================

set size=${xsize}T/${ysize}

gmt begin fig/tide_case png

@ ifig = 1

foreach ist(${id_ave} ${id_best} ${id_worst})

    echo "Station $ist"

    #===DATA
    set nnn=`printf "%03d" ${ist}`
    set input=dat/${nnn}_mave.dat
    
    #---Obs.
    @ obs_col = 1 + 2 * ${ndat} + 1
    gawk -v obs_col=${obs_col} '{if($obs_col != -999) print $1,$obs_col > "obs.20"}' ${input}
    gawk -v obs_col=${obs_col} 'BEGIN{n=0;inblock=0} {if($obs_col == -999){inblock=0}else{if(inblock==0){n++;outfile=sprintf("obs%02d.20",n);inblock=1} print $1,$obs_col > outfile}}' ${input}

    #---Analysis
    @ idat = 1
    while(${idat} <= ${ndat})
	@ dat_col = 1 + ${idat}
	gawk -v dat_col=${dat_col} -v out="dat${idat}.20" '{if($dat_col != -999) print $1,$dat_col > out}' ${input}
	@ idat++
    end

    #===Figure
    #---Box
    eval set ys = \$ys_${ifig} 
    eval set ye = \$ye_${ifig}
    eval set BAy = \$BAy_${ifig}
    set iy = `gawk -v ysize=$ysize 'BEGIN{print ysize + 2.5}'`

    set range=${start_date}/${end_date}/${ys}/${ye}
    if(${ifig} == 1)then
	gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -X3 -Y20
    else
	gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -Y-${iy}
    endif

    #---Label
    gmt text -F+f14p,0,black+jLB -N -D0/0.2c <<EOF
${start_date} ${ye} ${label[$ifig]}
EOF

    #---RMSD (top left)
    echo "RMSD (m):" | gmt text -F+cTL+jTL+f9p,Helvetica,black -Dj0.2c/0.2c
    set input=dat/rmsd_ave.dat
    @ idat = 1
    while(${idat} <= ${ndat})
	@ dat_col = 3 + ${ndat} + ${idat}
	set rmsd=`gawk -v ist=${ist} -v dat_col=${dat_col} '$1 == ist {printf "%.3f", $dat_col}' ${input}`
	set xoff=`gawk -v i=${idat} 'BEGIN{print 1.8 + 1.0*(i-1)}'`
	echo "${rmsd}" | gmt text -F+cTL+jTL+f9p,Helvetica-Bold,${color[$idat]} -Dj${xoff}c/0.2c
	@ idat++
    end

    #---Obs
    @ ifile = 1
    @ nfile =99
    while(${ifile} <= ${nfile})

	set nn=`printf "%02d" ${ifile}`
	set obsfile="obs${nn}.20"
	if(-f ${obsfile} && ${ifile} == 1 && ${ifig} == 1)then
	    gmt psxy ${obsfile} -W4,gray -t50 -l"Obs"+N3
	else if(-f ${obsfile})then
	    gmt psxy ${obsfile} -W4,gray -t50
	endif

	@ ifile++

    end	

    #---Analysis
    @ idat = 1
    while(${idat} <= ${ndat})

	if(${idat} == 1)then
	    set pt = 2
	else
	    set pt = 0.5
	endif
    
	set datfile="dat${idat}.20"
	if(-f ${datfile} && ${ifig} == 1)then
	    gmt psxy ${datfile} -W${pt},${color[$idat]} -l${legend[$idat]}
	else
	    gmt psxy ${datfile} -W${pt},${color[$idat]}	
	endif

	@ idat++
    end
    
    if(${ifig} == 1)then	
	gmt legend -DjBC+jBC+o0/0.2c+w11c -F+gwhite+pblack --FONT_ANNOT_PRIMARY=9p
    	#gmt legend -DjRB+jRB+o0.2/0.2 -F+gwhite+pblack --FONT=7p
    endif
	
    @ ifig++

    rm -f *.20
    
end

gmt end
