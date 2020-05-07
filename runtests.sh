#!/bin/sh

export OMP_NUM_THREADS=2
export RUNJOB="mpirun -np 2"
#export RUNJOB="mpirun -np 2 --bind-to none"
#export RUNJOB="true"

./testscript.sh
