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
set range=2003-01-01/2023-12-31/0/0.2
set BApx=f1Y
set BAsx=a5Y+l"Date"
set BAy=a0.1f0.02g99+l"RMSD\040\046\040Spread\040(m)"
set BAl=WSne
set label=("RMSD" "Spread")
set color=("black" "cyan")

@ i = 1
@ n = 915

while($i <= $n)

    echo $i
    set nnn=`printf "%03d" ${i}`

    #---Extract data
    set input=dat/rmsd_mave.dat
    gawk -v i=$i '{if($1 == i && $9 != -999) print $4,$9 > "rmsd.20"}' ${input}
    
    set input=dat/sprd_mave.dat
    gawk -v i=$i '{if($1 == i && $9 != -999) print $4,$9 > "sprd.20"}' ${input}

    set input=dat/cor_ave.dat
    gawk -v i=$i '{if($1 == i && $8 != -999) print $8 > "cor.20"}' ${input}
    if(-f cor.20)then
	set input=cor.20
	set cor=`gawk '{printf "%.3f", $1}' ${input}`
    endif

    if(-f rmsd.20 && -f sprd.20 && -f cor.20)then
    
	gmt begin fig/stat${nnn} png
	
	    gmt basemap -JX${size} -R${range} -Bpx${BApx} -Bsx${BAsx} -By${BAy} -B${BAl} -X3 -Y20
        
	    gmt psxy rmsd.20 -W1,${color[1]} -l${label[1]}
	    gmt psxy sprd.20 -W1,${color[2]} -l${label[2]}

	    gmt legend -DjRT+jRT+o0.2/0.2 -F+gwhite+pblack --FONT=10p

	    gmt text -F+f14p,0,black+jLT -N <<EOF
2003-01-01 0.2 Station:${i} Correlation: ${cor}
EOF
        
	gmt end

    endif

    rm -f *.20
    @ i++

end
