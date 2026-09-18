#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 #nb*nvar=8
#SBATCH -t 24:00:00

syr=$1; smon=$2; sday=$3
eyr=$4; emon=$5; eday=$6

source /etc/profile
source ~/.bashrc

nb=2
nvar=4

for ib in $(seq 1 $nb); do
    [ $ib -eq 1 ] && bname="keo" || bname="papa"
    for ivar in $(seq 1 $nvar); do
        case $ivar in
            1) varname="t" ;;
            2) varname="s" ;;
            3) varname="u" ;;
            4) varname="v" ;;
        esac
        ./make_data.out ${syr} ${smon} ${sday} ${eyr} ${emon} ${eday} ${ib} ${ivar} > make_data_${bname}_${varname}.log &
    done
done
wait
