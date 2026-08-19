#!/bin/bash
#SBATCH -N 1
#SBATCH -t 04:00:00

syr=$1
smon=$2
sday=$3
eyr=$4
emon=$5
eday=$6
yyyy=$7
mm=$8
    
source /etc/profile
./make_data.out ${syr} ${smon} ${sday} ${eyr} ${emon} ${eday} > ${yyyy}${mm}.log
