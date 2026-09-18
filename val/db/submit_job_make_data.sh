#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH -t 24:00:00

source /etc/profile
source ~/.bashrc
export OMP_NUM_THREADS=1

while [ $# -gt 0 ]; do
    
    iyr=$1; imon=$2; nday=$3; yyyy=$4; mm=$5
    shift 5

    ./make_data.out ${iyr} ${imon} 1 ${iyr} ${imon} ${nday} > ${yyyy}${mm}.log &

done
wait
    
