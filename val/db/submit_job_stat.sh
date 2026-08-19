#!/bin/bash
#SBATCH -N 1
#SBATCH -t 10:00:00

syr=$1
smon=$2
sday=$3
eyr=$4
emon=$5
eday=$6
    
source /etc/profile
./stat.out ${syr} ${smon} ${sday} ${eyr} ${emon} ${eday} > stat.log
