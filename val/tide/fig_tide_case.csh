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

set size=12T/6
set range=2003-01-01/2023-12-31/-1/1
set BApx=f1Y
set BAsx=a5Y+l"Date"
set BAy=a0.5f0.1g99+l"Sea\040level\040(m)"
set BAl=WSne
set label=("(a) Typical case: Tauranga" "(b) Best case for LORA: Kukup" "(C) Worst case for LORA: Booby Island")
set legend=("Obs" "LORA" "GLORYS" "ORAS5" "C-GLORS")
set color=("black" "blue" "cyan" "orange" "green")

gmt begin fig/tide_case png

@ i = 1

foreach ist(73 325 336)

    echo "Station $ist"

    #---Extract data
    set nnn=`printf "%03d" ${ist}`
    set input=dat/${nnn}.dat
    
    #Obs.
    gawk '{if($6 != -999) print $1,$6 > "dat1.20"}' ${input}
    gawk 'BEGIN{n=0;inblock=0} {if($6 == -999){inblock=0}else{if(inblock==0){n++;outfile=sprintf("dat1_%02d.20",n);inblock=1} print $1,$6 > outfile}}' ${input}

    #Analysis
    gawk '{if($2 != -999) print $1,$2 > "dat2.20"}' ${input}
    gawk '{if($3 != -999) print $1,$3 > "dat3.20"}' ${input}
    gawk '{if($4 != -999) print $1,$4 > "dat4.20"}' ${input}
    gawk '{if($5 != -999) print $1,$5 > "dat5.20"}' ${input}

    if($i == 1)then
	gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -X3 -Y20
    else
	gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -Y-8.5
    endif

    gmt text -F+f14p,0,black+jLB -N <<EOF
2003-01-01 1.05 ${label[$i]}
EOF
    
    #---For legend
    if($i == 1)then
	@ idat = 1
	@ ndat = 5
	while($idat <= $ndat)

	if(! -f dat${idat}.20)then
	    @ idat++
	    continue
	endif

	gmt psxy dat${idat}.20 -W4,${color[$idat]} -l${legend[$idat]}
	
	@ idat++

	end

	gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl}+gwhite
	
    endif

    #---Main plot
    foreach idat(1 5 4 3 2)

	if($idat == 1)then
	    @ ifile = 1
	    @ nfile = 99
	    while($ifile <= $nfile)
	    
		set nn=`printf "%02d" ${ifile}`
		if(-f dat${idat}_${nn}.20)then
		    gmt psxy dat${idat}_${nn}.20 -W4,${color[$idat]}
		endif

		@ ifile++
			
	    end
	else
	    gmt psxy dat${idat}.20 -W0.2,${color[$idat]}
	endif
    
    end
    
    if($i == 1)then
	gmt legend -DjRB+jRB+o0.2/0.2 -F+gwhite+pblack --FONT=10p
    endif
	
    @ i++

    rm -f *.20
    
end

gmt end
