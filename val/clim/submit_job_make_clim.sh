#!/bin/bash

source /etc/profile
source ~/.bashrc

thread=$1
nproc=$2

export OMP_NUM_THREADS=${thread}
export PARALLEL=${thread}

mpirun -np ${nproc} ./make_mclim.out > make_mclim.log

