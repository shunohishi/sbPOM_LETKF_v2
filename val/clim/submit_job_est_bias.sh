#!/bin/bash
#SBATCH -N 1
#SBATCH -t 24:00:00
    
source /etc/profile
source ~/.bashrc

./est_bias.out > est_bias.log

